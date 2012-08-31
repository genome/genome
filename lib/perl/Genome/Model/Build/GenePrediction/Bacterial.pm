package Genome::Model::Build::GenePrediction::Bacterial;

use strict;
use warnings;

use Genome;
use YAML;
use Carp;
use File::Basename;

class Genome::Model::Build::GenePrediction::Bacterial {
    is => 'Genome::Model::Build::GenePrediction',
};

sub locus_tag {
    my $self = shift;
    my $model = $self->model;
    return $model->locus_id . $model->run_type;
}

sub assembly_name {
    my $self = shift;
    my $model = $self->model;
    return ucfirst($model->organism_name) . '_' . $self->locus_tag . '.velv.amgap';
}

sub org_dirname {
    my $self = shift;
    my $model = $self->model;
    return substr(ucfirst($model->organism_name), 0, 1) .
           substr($model->organism_name, index($model->organism_name, "_"));
}

sub organism_name {
    my $self = shift;
    my $model = $self->model;
    return ucfirst($model->organism_name);
}

sub config_file_path {
    my $self = shift;
    return $self->data_directory . '/' . $self->config_file_name;
}

sub config_file_name {
    return "config.yaml";
}

sub sequence_file_directory {
    my $self = shift;
    my $model = $self->model;
    my ($name, $path) = fileparse($model->assembly_contigs_file);
    return $path;
}

sub sequence_file_name {
    my $self = shift;
    my $model = $self->model;
    my ($name, $path) = fileparse($model->assembly_contigs_file);
    return $name;
}

sub create_config_file {
    my $self = shift;
    my $model = $self->model;
    my $config_file_path = $self->config_file_path;

    if (-e $config_file_path) {
        $self->status_message("Removing existing configuration file at $config_file_path");
        my $unlink_rv = unlink $config_file_path;
        confess "Trouble removing configuration file at $config_file_path!" unless $unlink_rv;
    }

    $self->status_message("Creating configuration file at $config_file_path");

    my %params = (
        acedb_version    => $model->acedb_version,
        assembly_name    => $self->assembly_name,
        assembly_version => $model->assembly_version,
        cell_type        => uc($model->domain),
        gram_stain       => $model->gram_stain,
        locus_id         => $model->locus_id,
        locus_tag        => $self->locus_tag,
        minimum_length   => $model->minimum_sequence_length,
        ncbi_taxonomy_id => $model->ncbi_taxonomy_id,
        nr_db            => $model->nr_database_location,
        org_dirname      => $self->org_dirname,
        organism_name    => $self->organism_name,
        path             => $self->data_directory,
        pipe_version     => $model->pipeline_version,
        project_type     => $model->project_type,
        runner_count     => $model->runner_count,
        seq_file_dir     => $self->sequence_file_directory,
        seq_file_name    => $self->sequence_file_name,
        skip_acedb_parse => $model->skip_acedb_parse,
    );

    my $rv = YAML::DumpFile($config_file_path, \%params);
    unless ($rv) {
        $self->error_message("Could not create config file at $config_file_path!");
        return;
    }

    return $config_file_path;
}

1;
