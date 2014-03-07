package Genome::Model::Tools::Htseq::Count;
use strict;
use warnings;
use Genome;

class Genome::Model::Tools::Htseq::Count {
    is => 'Genome::Command::WithSavedResults',
    parallelize_by => ['alignment_results'],
    has_input => [
        alignment_results => { 
            is => 'Genome::InstrumentData::AlignmentResult',
            is_many => 1,
            example_values => [
                { 'instrument_data.sample.patient.common_name like' => 'HCC%' }
            ],
            where => [ 'instrument_data.sample.extraction_type in' => ['rna','cdna','total rna'] ], 
            doc => 'alignment results, typically from an RNA aligner',
        },
    ],
    has_param => [
        app_version => {
            is => 'SoftwareVersion',
            default_value => '0.5.4p1',
            valid_values => [ Genome::Sys->sw_versions('htseq','htseq-count') ],
            is_param => 1,
            doc => 'the version of htseq-count to use',
        },
        result_version => {
            # Required by all ::WithSavedResults
            is => 'Integer',
            valid_values => [1],
            default_value => '1',
            doc => 'the version of results, which may iterate as this logic iterates',
        },
        mode => {
            is => 'Text',
            valid_values => ['intersection-strict'],
            default_value => 'intersection-strict',
            doc => 'mode',
        },
        minaqual => {
            is => 'Integer',
            default_value => 1,
            doc => 'minaqual'
        },
        whitelist_alignments_flags => {
            is => 'Text',
            is_optional => 1,
            doc => 'require alignments to match the specified flags (-f): 0x0002 limits to only properly-paired alignments',
        },
        blacklist_alignments_flags => {
            is => 'Text',
            is_optional => 1,
            default_value => '0x0104',
            doc => 'exclude alignments which match the specified flags (-F): 0x0104 excludes non-primary alignments and unaligned reads',
        },
        limit => {
            is => 'Number',
            is_optional => 1,
            is_param => 1,
            doc => 'limit the number of alignments to the first N (for testing)',
        },
        #sort_strategy => { 
        #    is => 'Text',
        #    is_optional => 1,
        #    doc => 'samtools $VERSION [-p1 -p2]'
        #},
        #bam_view_strategy => { 
        #    is => 'Text',
        #    is_optional => 1,
        #    doc => 'samtools $VERSION [-p1 -p2]'
        #},
    ],
    has_optional_output => [
        # outputs should NOT need to be flagged as optional, 
        # but until a bug is fixed they must be
        gene_hits_file_path => { 
            is => 'FilesystemPath',
            doc => 'the path to the file of hit counts by gene',
        },
        transcript_hits_file_path => {
            is => 'FilesystemPath',
            doc => 'the path to the file of hit counts by transcript',
        },
        output_dir => {
            is => 'FilesystemPath',
            is_optional => 1,
            is_input => 1,
            doc => 'the results directory',
        },
    ],
    doc => 'generate htseq results for an (annotation-based) alignment result',
};

sub _execute_v1 {
    my $self = shift;
    $self->debug_message("Using HTSeq version " . $self->app_version . ', result version ' . $self->result_version . '.');

    my @alignment_result = $self->alignment_results;
    if (@alignment_result > 1) {
        # This should never happen with the current logic, since it will
        # dividue executions up by alignment result, but we may later support
        # batching.
        Carp::confess("Support for multiple alignment result inputs in the same excution is not implemented yet!");
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
    
    my $output_dir = $self->output_dir;
    $self->debug_message("Output destination directory: $output_dir");

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
    }
    else {
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

    my $final_bam = $sorted_bam;
    if (0) { 
        # This logic is only needed if we are not filtering out unaligned reads,
        # and we are using alignments from tophat before version < 2.0.7.
        # TODO: determine if we should keep this
        $self->debug_message("older tophat bams require cleanup of query names and mate information");

        $self->debug_message("cleaning up bam for fixmate");
        my $sorted_cleaned_bam = "$tmp/all_sequences.namesorted.cleaned.bam";
        my $clean_cmd = "$samtools_path view -h $sorted_bam | "
            . ($whitelist_alignments_flags ? " -f $whitelist_alignments_flags " : '')
            . ($blacklist_alignments_flags ? " -F $blacklist_alignments_flags " : '')
            . ' | '
            . ($limit ? "head -n " . ($limit + $header_size) . " | " : '')
            . q| perl -ne 'if (substr($_,0,1) eq q{@}) { print } else { @F = split(qq{\\t}, $_); $F[0] =~ s/\/[12]$//; $F[6] = "*"; $F[7] = "0";  print join(qq{\\t},@F) }' |
            . ' | '
            . " samtools view -S - -b > $sorted_cleaned_bam";
        Genome::Sys->shellcmd(cmd => $clean_cmd);

        $self->debug_message("re-running fixmate with cleaned-up queries and fixmate information");
        $final_bam = "$tmp/all_sequences.namesorted.cleaned.fixed.bam";
        Genome::Sys->shellcmd(cmd => "$samtools_path fixmate $sorted_cleaned_bam $final_bam");
    }

    for my $type (qw/transcript gene/) {
        $self->debug_message("Produce per-$type results...");
        
        # from Malachi's notes in JIRA issue TD-490
        # samtools view -h chr22_tumor_nonstranded_sorted.bam | htseq-count --mode intersection-strict --stranded no --minaqual 1 --type exon --idattr transcript_id - chr22.gff > transcript_tumor_read_counts_table.tsv 
   
        my $cmd = $samtools_path 
            . " view -h '$final_bam' "
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

        my $method = $type . '_hits_file_path';
        $self->$method("$output_dir/${type}-counts.tsv");
    }

    return 1;
}

sub _htseq_stranded_param {
    my $self = shift;
    my $transcript_strand = shift;

    if ($transcript_strand eq 'unstranded') {
        return 'no';
    }
    elsif ($transcript_strand eq 'firststrand') {
        return 'yes';
    }
    elsif ($transcript_strand eq 'secondstrand') {
        return 'reverse';
    }
    else {
        die $self->error_message("Unknown transcript_strand $transcript_strand!  expected unstranded, firstread or secondread.");
    }
}

sub _merge_v1 {
    my $self = shift;
    my @underlying = @_;
   
    $self->_default_merge(@underlying);

    my $subdir = $self->output_dir . '/underlying_results';

    my $cmd1 = '';
    my $cmd2 = '';

    for my $r (@underlying) {
        my $sub_subdir = $r->id;
        my $path = $subdir . '/' . $sub_subdir;
    
        if ($r == $underlying[0]) {
            $cmd1 .= 'cat ' . $sub_subdir . '/gene-counts.tsv';
            $cmd2 .= 'cat ' . $sub_subdir . '/transcript-counts.tsv';
        }
        else {
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

sub help_synopsis {
    return <<EOS

gmt htseq count --alignment-results "instrument_data.id=2890686892" --app-version 0.5.4p1

gmt htseq count --alignment-results "instrument_data.sample.name='H_MU-752713-1209062'" --app-version 0.5.4p1

gmt htseq count --alignment-results "instrument_data.sample.patient.common_name like 'HCC%'" --app-version 0.5.4p1

# skip any data sets flagged as test data
gmt htseq count --alignment-results "instrument_data.sample.patient.common_name like 'HCC%' and test_name is null"
EOS
}

sub help_detail {
    return <<EOS
This tool runs "htseq-count" from the HTSeq package, developed by Simon Anders at EMBL Heidelberg (Genome Biology Unit).
http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html
EOS
}

sub _additional_help_sections {
    return (
       "INPUTS",
       <<EOS,
It operates on "alignment results" from RNA/cDNA instrument data, such as those produced by tophat or rna-star.  These
results have associated annotation used during alignment, and that annotation is fed into htseq.

When run on more than one alignment result, each is processed individually, and the results are also aggregated into an additinal data product which is returned.
EOS
       "OUTPUTS",
       <<EOS,
The output is a list if "hit counts" per gene, and a second list of hit counts per transcript.
EOS
        "NOTE",
        <<EOS,
This tool saves software results with each execution, and shortcuts on subsequent runs to avoid duplicating effort.
EOS
  );
}

sub _doc_manual_body {
    my $help = shift->help_detail;
    $help =~ s/\n+$/\n/g;
    return $help;
    # expect to return POD
}

sub _doc_authors {
    return <<EOS
 Scott Smith
 Malachi Griffith, Ph.D.
 Obi Griffith, Ph.D.
EOS
}

sub _doc_copyright_years {
    (2013);
}

sub _doc_license {
    my $self = shift;
    my (@y) = $self->_doc_copyright_years;  
    my $range;
    if (@y == 1) {
        $range = "$y[0]";
    }
    elsif (@y > 1) {
        $range = "$y[0]-$y[-1]";
    }
    return <<EOS
Copyright (C) $range Washington University in St. Louis.

It is released under the Lesser GNU Public License (LGPL) version 3.  See the 
associated LICENSE file in this distribution.
EOS
}

sub _doc_credits {
    return ('','Simon Anders at EMBL Heidelberg (Genome Biology Unit) is the author of HTSeq.');
}

sub _doc_see_also {
    return <<EOS
B<Genome::Model::RnaSeq>(3), B<Genome::InstrumentData::AlignmentResult::PerLaneTophat>(3)
B<genome-model-rnaseq>(1), B<genome-instrument-data-align-per-lane-tophat>(1)
EOS
}

1;

