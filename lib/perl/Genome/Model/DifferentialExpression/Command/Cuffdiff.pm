package Genome::Model::DifferentialExpression::Command::Cuffdiff;

use strict;
use warnings;

use Genome;

my $DEFAULT_LSF_RESOURCE = "-R 'select[mem>=64000] rusage[mem=64000] span[hosts=1]' -M 64000000 -n 4";

class Genome::Model::DifferentialExpression::Command::Cuffdiff {
    is => ['Command::V2'],
    has => [
        build => { is => 'Genome::Model::Build::DifferentialExpression', id_by => 'build_id', },
    ],
    has_input_output => {
        build_id => {},
    },
    has_param => [
        lsf_resource => { default_value => $DEFAULT_LSF_RESOURCE },
    ],
};

sub execute {
    my $self = shift;

    my $build = $self->build;
    my $model = $build->model;

    my $reference_sequence_build = $model->reference_sequence_build;
    my $annotation_build = $model->annotation_build();

    my $output_directory = $build->differential_expression_directory;
    unless (-d $output_directory) {
        Genome::Sys->create_directory($output_directory);
    }

    # This provides a sorted map of labels -> ids as an array of hashes
    my @condition_pairs = $model->condition_pairs_sorted_mapping;

    my $params = $model->differential_expression_params;
    $params .= ' --labels '. (join ',', map { $_->{label} } @condition_pairs);

    my $reference_fasta_path = $reference_sequence_build->full_consensus_path('fa');
    $params .= ' --frag-bias-correct '. $reference_fasta_path;
    if ($model->differential_expression_mask_reference_transcripts) {
        my $mask_file_method = $model->differential_expression_mask_reference_transcripts .'_file';
        my $mask_gtf_path = $annotation_build->$mask_file_method('gtf',$reference_sequence_build->id);
        $params .= ' --mask-file '. $mask_gtf_path;
    }
    my @condition_model_ids = map { $_->{model_ids} } @condition_pairs;
    my $bam_file_paths = '';
    for my $model_ids (@condition_model_ids) {
        my @condition_bam_file_paths;
        for my $model_id (@$model_ids) {
            my $rna_seq_model = Genome::Model->get($model_id);
            # Or should we just fine the latest AlignmentResult....
            my $last_succeeded_build = $rna_seq_model->last_succeeded_build;
            my $bam_file = $last_succeeded_build->alignment_result->bam_file;
            push @condition_bam_file_paths, $bam_file;
        }
        my $condition_bam_files_string = join(',',@condition_bam_file_paths);
        $bam_file_paths .= $condition_bam_files_string .' ';
    }
    my %cuffdiff_params = (
        transcript_gtf_file => $build->transcript_gtf_file_path,
        bam_file_paths => $bam_file_paths,
        params => $params,
        use_version => $model->differential_expression_version,
        output_directory => $output_directory,
    );
    my $cuffdiff_cmd = Genome::Model::Tools::Cufflinks::Cuffdiff->execute(%cuffdiff_params);
    unless ($cuffdiff_cmd and $cuffdiff_cmd->result) {
        $self->error_message('Failed to execute Cuffdiff with params: '. Data::Dumper::Dumper(%cuffdiff_params));
        return;
    }
    return 1;
}

1;

