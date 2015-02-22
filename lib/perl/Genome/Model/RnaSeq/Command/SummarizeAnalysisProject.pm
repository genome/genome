package Genome::Model::RnaSeq::Command::SummarizeAnalysisProject;

use strict;
use warnings;

use Genome;

class Genome::Model::RnaSeq::Command::SummarizeAnalysisProject {
    is => 'Command::V2',
    has => [
        analysis_project => {
            is => 'Genome::Config::AnalysisProject',
            shell_args_position => 1,
            doc => 'The analysis project to summarize RNA-seq results.',
        },
        output_directory => {
            doc => 'A directory to output the resulting summary files.',
        },
        model_identifier => {
            is_optional => 1,
            default_value => 'subject_name',
            valid_values => ['name','id','subject_name','individual_common_name'],
        },
    ],
};

sub help_synopsis {
    return <<"EOS"
    genome model rna-seq summarize-analysis-project
EOS
}

sub help_brief {
    return "Accumulate RNAseq results for all models in an Analysis Project.";
}

sub help_detail {
    return <<EOS
SOMETHING ELSE.
EOS
}

sub execute {
    my $self = shift;

    my @models = $self->analysis_project->models;
    my @rna_models = grep { $_->isa('Genome::Model::RnaSeq') } @models;
    unless (@rna_models) {
        die('No RNA-seq models found!');
    }
    unless (-d $self->output_directory) {
        Genome::Sys->create_directory($self->output_directory);
    }
    my $base_file_path = $self->output_directory .'/'. $self->analysis_project->id;
    my $gene_fpkm_tsv_file =  $base_file_path .'-gene_fpkm_matrix.tsv';
    unless (-e $gene_fpkm_tsv_file) {
        unless (Genome::Model::RnaSeq::Command::FpkmMatrix->execute(
            models => \@rna_models,
            gene_fpkm_tsv_file => $gene_fpkm_tsv_file,
            model_identifier => $self->model_identifier,
        )->result) {
            die('Failed to execute FpkmMatrix!');
        }
    }
    my $rna_seq_metrics_tsv_file = $base_file_path .'-rna_seq_metrics.tsv';
    my $normalized_coverage_tsv_file = $base_file_path .'-normalized_coverage.tsv';
    unless (-e $rna_seq_metrics_tsv_file && -e $normalized_coverage_tsv_file) {
        unless (Genome::Model::RnaSeq::Command::RnaSeqMetrics->execute(
            models => \@rna_models,
            metrics_tsv_file => $rna_seq_metrics_tsv_file,
            normalized_transcript_coverage_file => $normalized_coverage_tsv_file,
            model_identifier => $self->model_identifier,
        )->result) {
            die('Failed to execute RnaSeqMetrics!');
        }
    }
    my $splice_junction_metrics_tsv_file = $base_file_path .'-splice_junction_metrics.tsv';
    unless (-e $splice_junction_metrics_tsv_file ) {
        unless (Genome::Model::RnaSeq::Command::SpliceJunctionMetrics->execute(
            models => \@rna_models,
            splice_junction_metrics_tsv_file => $splice_junction_metrics_tsv_file,
            model_identifier => $self->model_identifier,
        )->result) {
            die('Failed to execute SpliceJunctionMetrics!');
        }
    }

    my $rawcount_tsv_file = $base_file_path .'-gene_rawcount_matrix.tsv';
    unless (-e $rawcount_tsv_file) {
        unless (Genome::Model::RnaSeq::Command::RawcountMatrix->execute(
            models => \@rna_models,
            output => $rawcount_tsv_file,
            model_identifier => $self->model_identifier,
        )->result) {
            die ('Failed to execute RawcountMatrix!');
        }
    }

    my @naming_headers = ($self->model_identifier,'model_id','build_id');

    # TODO: Summarize Chimerascan Results
    my $model_identifier_method = $self->model_identifier;
    my $naming_key_tsv_file = $base_file_path .'-naming_key.tsv';
    my $naming_key_writer = Genome::Utility::IO::SeparatedValueWriter->create(
        output => $naming_key_tsv_file,
        separator => "\t",
        headers => \@naming_headers,
    );
    unless ($naming_key_writer) {
        die('Failed to open naming key file: '. $naming_key_tsv_file);
    }
    for my $model (@rna_models) {
        my $label = $model->$model_identifier_method;
        my $build = $model->last_succeeded_build;

        my %data = (
            $model_identifier_method => $label,
            model_id => $model->id,
            build_id => $build->id,
        );
        $naming_key_writer->write_one(\%data);

        my $fusion_output_dir = $self->output_directory .'/fusions/'. $model->id;
        unless (-d $fusion_output_dir) {
            Genome::Sys->create_directory ($fusion_output_dir);
        }
        my $fusion_build_dir = $build->data_directory .'/fusions';
        my @fusion_file_paths = glob($fusion_build_dir .'/*.bedpe');
        for my $fusion_file_path (@fusion_file_paths) {
            my ($filename,$dirname,$suffix) = File::Basename::fileparse($fusion_file_path,qw/\.bedpe/);
            my $new_file_path = $fusion_output_dir .'/'. $filename . $suffix;
            unless (File::Copy::copy($fusion_file_path,$new_file_path)) {
                die('Failed to copy fusion bedpe file: '. $fusion_file_path .' => '. $new_file_path);
            }
        }
    }
    $naming_key_writer->output->close;
    return 1;
}


1;
