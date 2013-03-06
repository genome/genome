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
        output_dir => { 
            is => 'FilesystemPath',
            is_optional => 1,
            doc => 'the results directory will be symlinked here',
        },
    ],
    has_param => [
        app_version => {
            is => 'SoftwareVersion',
            valid_values => [ Genome::Sys->sw_versions('htseq','htseq-count') ],
            doc => 'the version of htseq-count to use',
        },
        result_version => {
            is => 'Integer',
            default_value => '1',
            doc => 'the version of results, which may iterate as this logic iterates'
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
            default_value => '0x0104',
            doc => 'exclude alignments which match the specified flags (-F): 0x0104 excludes non-primary alignments and unaligned reads',
        },
        limit => {
            is => 'Number',
            is_optional => 1,
            doc => 'limit the number of alignments to the first N (for testing)'
        },
        #sort_strategy => { 
        #    is => 'Text',
        #    is_optional => 1,
        #    doc => 'samtools $VERSION, etc,'
        #},
        #bam_view_strategy => { 
        #    is => 'Text',
        #    is_optional => 1,
        #    doc => 'samtools $VERSION, etc,'
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
    ],
    doc => 'generate htseq results for an (annotation-based) alignment result',
};

sub _execute_v1 {
    my $self = shift;
    $self->status_message("tool version " . $self->app_version . ', result version ' . $self->result_version);

    my @alignment_result = $self->alignment_results;
    if (@alignment_result > 1) {
        die "support for multiple alignment result inputs in the same excution is not implemented yet!";
    }
    my $alignment_result = $alignment_result[0];
    $self->status_message("using alignment result " . $alignment_result->__display_name__);

    my $instrument_data = $alignment_result->instrument_data;
    unless ($instrument_data->sample->extraction_type =~ /rna|cdna/i) {
        die $self->error_message(
            "this step can only run on alignments of RNA, but sample " 
            . $instrument_data->sample->__display_name__ 
            . " is " . $instrument_data->sample->extraction_type
        );
    }
    $self->status_message("sample extraction type: " . $instrument_data->sample->extraction_type);

    my $transcript_strand = $instrument_data->transcript_strand;
    unless ($transcript_strand) {
        die $self->error_message("transcript strand is not set for instrument data " . $instrument_data->__display_name__); 
    }

    my $htseq_stranded_param;
    if ($transcript_strand eq 'unstranded') {
        $htseq_stranded_param = 'no';
    }
    elsif ($transcript_strand eq 'firstread') {
        $htseq_stranded_param = 'yes';
    }
    elsif ($transcript_strand eq 'secondread') {
        $htseq_stranded_param = 'reverse';
    }
    else {
        die "unknown transcript_strand $transcript_strand!  expected unstranded, firstread or secondread";
    }
    $self->status_message("strandedness: $transcript_strand (htseq-count strandedness: $htseq_stranded_param");
    
    my $annotation_build = $alignment_result->annotation_build;
    $self->status_message("annotation build: " . $annotation_build->__display_name__);

    my $gff_file = $annotation_build->rna_features_gff_path;
    $self->status_message("using annotation features from: " . $gff_file);
   
    my $htseq_count_path = Genome::Sys->sw_path("htseq", $self->app_version, "htseq-count");
    $self->status_message("htseq-count version: " . $self->app_version . " running from $htseq_count_path");
    
    my $output_dir = $self->output_dir;
    $self->status_message("output dir: $output_dir");

    # The samtools version is not part of the params because it is not yet required for it to vary.
    # If it does need to vary in the future a param should be addeed and backfilled
    # to represent samtools 0.18.1 (current at the time of this writing).  Because
    # we only use this do dump the bam contents and to sort it this may not be necessary
    # to ever upgrade, and is unlikely to produce different results if it is upgraded.
    my $samtools_version = '0.1.18';
    my $samtools_path = Genome::Sys->sw_path("samtools", $samtools_version);
    $self->status_message("samtools version: $samtools_version (running from $samtools_path)");
   
    my $tmp = Genome::Sys->create_temp_directory();

    my $sorted_bam;
    if ($alignment_result->temp_scratch_directory) {
        # if the AR is in the middle of being built it will have a sorted bam already present in scratch space
        $self->status_message("found temp_scratch directory ...using a presorted bam available during the alignment process");
        $sorted_bam = $alignment_result->temp_scratch_directory . '/raw_all_sequences.bam.sort.bam';
        unless (-e $sorted_bam) {
            die $self->error_message("found temp scratch dir but no sorted BAM file $sorted_bam!");
        }
    }
    else {
        # a completed alignment result will need to have a sorted bam created
        $self->status_message("no temp_scratch_directory found: name sort the BAM in temp space");
        my $unsorted_bam = $alignment_result->output_dir . '/all_sequences.bam';
        my $sorted_bam_noprefix = "$tmp/all_sequences.namesorted";
        $sorted_bam = $sorted_bam_noprefix . '.bam';

        Genome::Sys->shellcmd(
            cmd => "$samtools_path sort -n -m 402653184 $unsorted_bam $sorted_bam_noprefix",
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
        $self->status_message("older tophat bams require cleanup of query names and mate information");

        $self->status_message("cleaning up bam for fixmate");
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

        $self->status_message("re-running fixmate with cleaned-up queries and fixmate information");
        $final_bam = "$tmp/all_sequences.namesorted.cleaned.fixed.bam";
        Genome::Sys->shellcmd(cmd => "$samtools_path fixmate $sorted_cleaned_bam $final_bam");
    }

    for my $type (qw/transcript gene/) {
        $self->status_message("produce per-$type results...");
        
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
            . " '$gff_file'"
            . " 1>'$output_dir/${type}-counts.tsv'"
            . " 2>'$output_dir/${type}.err'";
        
        Genome::Sys->shellcmd(
            cmd => $cmd, # "touch $output_dir/${type}-counts.tsv", #$cmd,
            input_files => [$sorted_bam, $gff_file],
            #output_files => ["$output_dir/${type}-counts.tsv"],
       );

        my $method = $type . '_hits_file_path';
        $self->$method("$output_dir/${type}-counts.tsv");
    }

    return 1;
}

1;

__END__

    if (wantarray) {
        my @result = $self->$method(@_);
    }
    elsif (not defined wantarray) {
        $self->$method;
    }
    else {
        my $result = $self->method;
    }
