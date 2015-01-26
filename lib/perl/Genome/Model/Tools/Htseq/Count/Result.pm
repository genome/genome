package Genome::Model::Tools::Htseq::Count::Result;

class Genome::Model::Tools::Htseq::Count::Result {
    is => ['Genome::SoftwareResult::StageableSimple'],
    has_param => [
        app_version => {
            is => 'SoftwareVersion',
            valid_values => [ Genome::Sys->sw_versions('htseq','htseq-count') ],
            is_param => 1,
            doc => 'the version of htseq-count to use',
        },
        result_version => {
            is => 'Integer',
            valid_values => [ 1, 2 ],
            doc => 'the version of results, which may iterate as this logic iterates',
        },
        mode => {
            is => 'Text',
            valid_values => [ "intersection-strict" ],
            doc => 'mode',
        },
        minaqual => {
            is => 'Integer',
            doc => 'minaqual',
        },
        whitelist_alignments_flags => {
            is => 'Text',
            is_optional => 1,
            doc => 'require alignments to match the specified flags (-f): 0x0002 limits to only properly-paired alignments',
        },
        blacklist_alignments_flags => {
            is => 'Text',
            is_optional => 1,
            doc => 'exclude alignments which match the specified flags (-F): 0x0104 excludes non-primary alignments and unaligned reads',
        },
        limit => {
            is => 'Number',
            is_optional => 1,
            doc => 'limit the number of alignments to the first N (for testing)',
        },
    ],
    has_input => [
        alignment_results => {
            is => 'Genome::InstrumentData::AlignmentResult',
            where => [ 'instrument_data.sample.extraction_type in' => [ "rna", "cdna", "total rna" ] ],
            example_values => [ { "instrument_data.sample.individual.common_name like" => "HCC%" } ],
            is_many => 1,
            doc => 'alignment results, typically from an RNA aligner',
        },
    ],
};

sub _run {
    my $self = shift;

    my @alignment_results = $self->alignment_results;
    if (@alignment_results > 1) {
        my @intermediate_results = $self->_parallelize(@alignment_results);
        return $self->_merge(@intermediate_results);
    } else {
        return $self->_run_htseq_count(@alignment_results);
    }
}

sub _parallelize {
    my $self = shift;
    my @alignment_results = @_;

    my $params = { };
    for my $param ($self->params) {
        my $name = $param->name;
        next if $name =~ /^alignment_results/; #derived params

        $params->{$name} = $self->$name;
    }

    # TODO: compose a workflow here instead of a linear run
    my $n = 0;
    for my $alignment_result (@alignment_results) {
        $n++;
        print "Breakdown $n: " . $alignment_result->__display_name__ . "\n";
    }

    my @intermediate_results;
    for my $alignment_result (@alignment_results) {
        my $class = ref $self;
        my $intermediate_result = $class->get_or_create(%$params, alignment_results => [$alignment_result]);
        unless ($intermediate_result) {
            die "failed to create partial for " . $alignment_result->__display_name__;
        }

        $intermediate_result->add_user(label => 'composes', user => $self);
        push @intermediate_results, $intermediate_result;
    }

    return @intermediate_results;
}

sub _merge {
    my $self = shift;
    my @underlying = @_;

    my $subdir = $self->temp_staging_directory . '/underlying_results';
    Genome::Sys->create_directory($subdir);

    for my $r (@underlying) {
        my $sub_subdir = $r->id;
        my $path = $subdir . '/' . $sub_subdir;
        Genome::Sys->create_symlink($r->output_dir, $path);
    }

    my $cmd1 = '';
    my $cmd2 = '';

    for my $r (@underlying) {
        my $sub_subdir = $r->id;
        my $path = $subdir . '/' . $sub_subdir;

        if ($r == $underlying[0]) {
            $cmd1 .= 'cat ' . $sub_subdir . '/gene-counts.tsv';
            $cmd2 .= 'cat ' . $sub_subdir . '/transcript-counts.tsv';
        } else {
            $cmd1 .= '| join -a 1 -a 2 - ' . $sub_subdir . '/gene-counts.tsv';
            $cmd2 .= '| join -a 1 -a 2 - ' . $sub_subdir . '/transcript-counts.tsv';
        }
    }

    my $sum = q/perl -nae '$sum = 0; for (@F[1..$#F]) { $sum += $_ }; print $F[0],"\t",$sum,"\n" '/;
    $cmd1 .= " | $sum > ../gene-counts.tsv";
    $cmd2 .= " | $sum > ../transcript-counts.tsv";

    Genome::Sys->shellcmd(cmd => "cd $subdir; $cmd1");
    Genome::Sys->shellcmd(cmd => "cd $subdir; $cmd2");

    return 1;
}

sub _htseq_stranded_param_v2 {
    my $self = shift;
    my $transcript_strand = shift;

    if ($transcript_strand eq 'unstranded') {
        return 'no';
    } elsif ($transcript_strand eq 'firststrand') {
        return 'reverse';
    } elsif ($transcript_strand eq 'secondstrand') {
        return 'yes';
    } else {
        die $self->error_message("Unknown transcript_strand $transcript_strand!  expected unstranded, firstread or secondread.");
    }
}

sub _htseq_stranded_param_v1 {
    my $self = shift;
    my $transcript_strand = shift;

    if ($transcript_strand eq 'unstranded') {
        return 'no';
    } else {
        die $self->error_message("Htseq::Count v1 only supports unstranded libraries, use v2 instead.");
    }
}

sub _htseq_stranded_param {
    my $self = shift;

    if ($self->result_version == 1) {
        $self->_htseq_stranded_param_v1(@_);
    } elsif ($self->result_version == 2) {
        $self->_htseq_stranded_param_v2(@_);
    } else {
        die "no _htseq_stranded_param implementation for version " . $self->result_version;
    }
}

sub _run_htseq_count {
    my $self = shift;

    $self->debug_message("Using HTSeq version " . $self->app_version . ', result version ' . $self->result_version . '.');

    my @alignment_result = @_;
    if (@alignment_result > 1) {
        # This should never happen with the current logic, since it will
        # dividue executions up by alignment result, but we may later support
        # batching.
        Carp::confess("Support for multiple alignment result inputs in the same execution is not implemented yet!");
    }
    my $alignment_result = $alignment_result[0];
    $self->debug_message("Using alignment result: " . $alignment_result->__display_name__);

    my $instrument_data = $alignment_result->instrument_data;
    unless ($instrument_data->sample->is_rna) {
        die $self->error_message(
            "This step can only run on alignments of RNA (cDNA), but sample "
            . $instrument_data->sample->__display_name__
            . " is " . $instrument_data->sample->extraction_type . '!'
        );
    }
    $self->debug_message("Sample extraction type: " . $instrument_data->sample->extraction_type);

    my $transcript_strand = $instrument_data->library->transcript_strand;
    unless ($transcript_strand) {
        die $self->error_message("Transcript strand is not set for instrument data " . $instrument_data->__display_name__ . "!");
    }

    my $htseq_stranded_param = $self->_htseq_stranded_param($transcript_strand);
    $self->debug_message("Strandedness: $transcript_strand (htseq-count stranded: $htseq_stranded_param)");

    my $annotation_build = $alignment_result->annotation_build;
    $self->debug_message("Annotation build: " . $annotation_build->__display_name__);

    my $gtf_file = $annotation_build->annotation_file('gtf',$annotation_build->reference_sequence->id);

    $self->debug_message("Using annotation features from: " . $gtf_file);

    my $htseq_count_path = Genome::Sys->sw_path("htseq", $self->app_version, "htseq-count");
    $self->debug_message("Executable htseq-count " . $self->app_version . " running from $htseq_count_path");

    my $output_dir = $self->temp_staging_directory;
    $self->debug_message("Output staging directory: $output_dir");
    $self->debug_message('Output destination directory: %s', $self->output_dir);

    # The samtools version is not part of the params because it is not yet required for it to vary.
    # If it does need to vary in the future a param should be addeed and backfilled
    # to represent samtools 0.1.18 (current at the time of this writing).  Because
    # we only use this do dump the bam contents and to sort it this may not be necessary
    # to ever upgrade, and is unlikely to produce different results if it is upgraded.
    my $samtools_version = '0.1.18';
    my $samtools_path = Genome::Sys->sw_path("samtools", $samtools_version);
    $self->debug_message("samtools version: $samtools_version (running from $samtools_path)");

    my $tmp = Genome::Sys->create_temp_directory();

    my $sorted_bam;
    if ($alignment_result->temp_scratch_directory) {
        # if the AR is in the middle of being built it will have a sorted bam already present in scratch space
        $self->debug_message("found temp_scratch directory ...using a presorted bam available during the alignment process");
        $sorted_bam = $alignment_result->temp_scratch_directory . '/raw_all_sequences.bam.sort.bam';
        unless (-e $sorted_bam) {
            die $self->error_message("found temp scratch dir but no sorted BAM file $sorted_bam!");
        }
    } else {
        # a completed alignment result will need to have a sorted bam created
        $self->debug_message("No temp_scratch_directory found: name sort the BAM in temp space.");
        my $unsorted_bam = $alignment_result->output_dir . '/all_sequences.bam';
        my $sorted_bam_noprefix = "$tmp/all_sequences.namesorted";
        $sorted_bam = $sorted_bam_noprefix . '.bam';

        Genome::Sys->shellcmd(
            cmd => "$samtools_path sort -n $unsorted_bam $sorted_bam_noprefix",
            input_files => [$unsorted_bam],
            output_files => [$sorted_bam],
        );
    }

    my @header = `$samtools_path view -H $sorted_bam`;
    my $header_size = scalar(@header);
    unless ($header_size > 0) {
        die "Unexpected missing header on $sorted_bam!!!";
    }

    my $limit = $self->limit;
    if ($limit) {
        $self->warning_message("SUBSAMPLING ALIGNMENTS because the limit parameter is set to $limit:");
    }

    my $whitelist_alignments_flags = $self->whitelist_alignments_flags;
    my $blacklist_alignments_flags = $self->blacklist_alignments_flags;
    my $minaqual = $self->minaqual;
    my $mode = $self->mode;

    for my $type (qw/transcript gene/) {
        $self->debug_message("Produce per-$type results...");

        # from Malachi's notes in JIRA issue TD-490
        # samtools view -h chr22_tumor_nonstranded_sorted.bam | htseq-count --mode intersection-strict --stranded no --minaqual 1 --type exon --idattr transcript_id - chr22.gff > transcript_tumor_read_counts_table.tsv

        my $cmd = $samtools_path
            . " view -h '$sorted_bam' "
            . ($whitelist_alignments_flags ? " -f $whitelist_alignments_flags " : '')
            . ($blacklist_alignments_flags ? " -F $blacklist_alignments_flags " : '')
            . ' | '
            . ($limit ? "head -n " . ($limit + $header_size) . " | " : '')
            . $htseq_count_path
            . " --mode $mode "
            . " --minaqual $minaqual "
            . " --stranded $htseq_stranded_param"
            . " --type exon " # used for both gene and transcript iterations
            . " --idattr ${type}_id"
            . " -"
            . " '$gtf_file'"
            . " 1>'$output_dir/${type}-counts.tsv'"
            . " 2>'$output_dir/${type}.err'";

        Genome::Sys->shellcmd(
            cmd => $cmd, # "touch $output_dir/${type}-counts.tsv", #$cmd,
            input_files => [$sorted_bam, $gtf_file],
            #output_files => ["$output_dir/${type}-counts.tsv"],
        );

    }

    return 1;
}

1;
