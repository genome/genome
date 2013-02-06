package Genome::Model::Tools::Htseq::Count;
use strict;
use warnings;
use Genome;

class Genome::Model::Tools::Htseq::Count {
    is => 'Command::V2',
    #shortcut_execute => 1,
    #save_results_as => 'Genome::Model::Tools::Htseq::Count::Result',
    has_input => [
        alignment_result => { 
            is => 'Genome::InstrumentData::AlignmentResult',
            doc => 'alignment results, typically from an RNA aligner',
        },
        output_dir => { 
            is => 'FilesystemPath',
            is_optional => 1,
            doc => 'the results directory',
        },
    ],
    has_param => [
        is_stranded => {
            is => 'Boolean',
            doc => 'indicate whether reads are only in the direction of transcription'
        },
        app_version => {
            is => 'SoftwareVersion',
            valid_values => [ Genome::Sys->sw_versions('htseq','htseq-count') ],
            doc => 'the version of htseq-count to use',
        },
        results_version => {
            is => 'Integer',
            default_value => '0.01',
            doc => 'the version of results, which may iterate as this logic iterates'
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
Genome::Model::Tools::Htseq::Count::Result->class;

sub execute {
    my $self = shift;
    my $results_version = $self->results_version;
    my $method = "_execute_v$results_version";
    $method =~ s/\./_/g;
    unless ($self->can($method)) {
        die "no implementation ($method) for version $results_version!";
    }
    $self->$method(@_);
}

# if output changes iterate the results_version and implement a new _execute_vX method

sub _execute_v1 {
    my $self = shift;

    my $alignment_result = $self->alignment_result;
    $self->status_message("using alignment result " . $alignment_result->__display_name__);

    my $annotation_build = $alignment_result->annotation_build;
    $self->status_message("annotation build is " . $annotation_build->__display_name__);
   
    # The samtools version is not part of the params because it is not yet required for it to vary.
    # If it does need to vary in the future a param should be addeed and backfilled
    # to represent samtools 0.18.1 (current at the time of this writing).  Because
    # we only use this do dump the bam contents and to sort it this may not be necessary
    # to ever upgrade, and is unlikely to produce different results if it is upgraded.
    my $samtools_version = '0.1.18';
    my $samtools_path = Genome::Sys->sw_path("samtools", $samtools_version);
    $self->status_message("samtools version $samtools_version running from $samtools_path");
   
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
        $self->status_message("No temp_scratch_directory found, presumably doing this post alignment.  Name sort the BAM in temp space.");
        my $unsorted_bam = $alignment_result->output_dir . '/all_sequences.bam';
        my $tmp = Genome::Sys->create_temp_directory();
        my $sorted_bam_noprefix = "$tmp/all_sequences.namesorted";
        $sorted_bam = $sorted_bam_noprefix . '.bam';
        Genome::Sys->shellcmd(
            cmd => "$samtools_path sort -n -m 402653184 $unsorted_bam $sorted_bam_noprefix",
            input_files => [$unsorted_bam],
            output_files => [$sorted_bam],
        );
    }

    my $htseq_count_path = Genome::Sys->sw_path("htseq", $self->app_version, "htseq-count");
    $self->status_message("htseq-count version " . $self->app_version . " running from $htseq_count_path");

    my $stranded = ($self->is_stranded ? 'yes' : 'no');
    $self->status_message("reads go only in the same direction as transcription: " . $stranded);

    my $gff_file = $annotation_build->rna_features_gff_path;
    $self->status_message("using annotation features from " . $gff_file);
    
    my $output_dir = $self->output_dir;
    $self->status_message("output dir: $output_dir");

    for my $type (qw/transcript gene/) {
        $self->status_message("produce per-$type results...");
        
        # from Malachi's notes in JIRA issue TD-490
        # samtools view -h chr22_tumor_nonstranded_sorted.bam | htseq-count --mode intersection-strict --stranded no --minaqual 1 --type exon --idattr transcript_id - chr22.gff > transcript_tumor_read_counts_table.tsv 
        my $cmd = "$samtools_path view -h '$sorted_bam' |"
            . $htseq_count_path
            . " --mode intersection-strict"
            . " --stranded $stranded"
            . " --minaqual 1"
            . " --type exon "
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

