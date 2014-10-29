package Genome::Model::DifferentialExpression::Command::Cuffmerge;

use strict;
use warnings;

use Genome;

my $DEFAULT_LSF_RESOURCE = "-R 'select[mem>=64000] rusage[mem=64000] span[hosts=1]' -M 64000000 -n 4";

class Genome::Model::DifferentialExpression::Command::Cuffmerge {
    is => ['Command::V2'],
    has => [
        build => { is => 'Genome::Model::Build', id_by => 'build_id', },
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

    my $output_directory = $build->transcript_convergence_directory;
    unless (-d $output_directory) {
        Genome::Sys->create_directory($output_directory);
    }

    my @transcript_gtf_paths;
    for my $input_model ($model->input_models) {
        my $rna_seq_build = $input_model->last_succeeded_build;
        my $transcript_gtf_path = $rna_seq_build->expression_directory .'/transcripts.gtf';
        push @transcript_gtf_paths, $transcript_gtf_path;
    }

    my $reference_fasta_path = $model->reference_sequence_build->full_consensus_path('fa');
    my $annotation_gtf_path = $model->annotation_build->annotation_file('gtf',$model->reference_sequence_build->id);
    if ($model->transcript_convergence_biotypes) {
        my $transcript_info_tsv_file = $model->annotation_build->transcript_info_file($model->reference_sequence_build->id);
        my $tmp_annotation_gtf_path = Genome::Sys->create_temp_file_path();
        my $limit_cmd = Genome::Model::Tools::Gtf::LimitByBiotype->execute(
            input_gtf_file => $annotation_gtf_path,
            output_gtf_file => $tmp_annotation_gtf_path,
            gene_biotypes => $model->transcript_convergence_biotypes,
            transcript_info_tsv_file => $transcript_info_tsv_file,
        );
        unless ($limit_cmd and $limit_cmd->result) {
            $self->error_message('Failed to limit transcripts by gene biotypes: '. $model->transcript_convergence_biotypes);
            return;
        }
        $annotation_gtf_path = $tmp_annotation_gtf_path;
    }

    my $transcript_convergence_params = eval($model->transcript_convergence_params);
    $transcript_convergence_params->{use_version} = $model->transcript_convergence_version;
    unless ($transcript_convergence_params->{input_gtf_paths}) {
        $transcript_convergence_params->{input_gtf_paths} = \@transcript_gtf_paths;
    }
    $transcript_convergence_params->{reference_fasta_path} = $reference_fasta_path;
    $transcript_convergence_params->{reference_gtf_path} = $annotation_gtf_path;

    # TODO: Setup as SoftwareResult
    $transcript_convergence_params->{output_directory} = $output_directory;
    my $merge_cmd = Genome::Model::Tools::Cufflinks::Cuffmerge->execute($transcript_convergence_params);
    unless ($merge_cmd and $merge_cmd->result) {
        $self->error_message('Failed to execute Cuffmerge with params: '. Data::Dumper::Dumper($transcript_convergence_params));
        return;
    }
    return 1;
}

1;

