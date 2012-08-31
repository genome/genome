package Genome::Model::Command::Define::GenePrediction::Bacterial;

use strict;
use warnings;

use Genome;
use Carp;

class Genome::Model::Command::Define::GenePrediction::Bacterial {
    is => 'Genome::Model::Command::Define::GenePrediction::Helper',
    has_optional => [
        dev => {
            is => 'Boolean',
            default => 0,
            doc => 'If set, use dev databases instead of production databases',
        },
        run_type => {
            is => 'String',
            default => 'DFT',
            doc => 'A three letter identifier appended to the locus ID',
        },
        assembly_version => {
            is => 'String',
            default => 'Version_1.0',
            doc => 'Notes the assembly version',
        },
        project_type => {
            is => 'String',
            default => 'HGMI',
            doc => 'The type of project this data is being generated for',
        },
        pipeline_version => {
            is => 'String',
            default => 'Version_1.0',
            doc => 'Notes the pipeline version',
        },
        acedb_version => {
            is => 'String',
            default => 'Version_5.0',
            doc => 'Notes the version of aceDB to upload results to',
        },
        nr_database_location => {
            is => 'String',
            default => '/gscmnt/gpfstest2/analysis/blast_db/gsc_bacterial/bacterial_nr/bacterial_nr',
            doc => 'Default location of the NR database, may be overridden with a local copy if specified',
        },
    ],
};

sub help_detail {
    return <<"EOS"
Two things are needed to define a bacterial gene prediction model: a processing profile
name and a taxon ID. The taxon ID is used to find an assembly model that would correspond
to this annotation model. If such a model cannot be found, it will be created if the 
--create-assembly-model flag is set. If a model IS found and no successful build is found,
then a build is started if the --start-assembly-build flag is set.

With a successful assembly build, all the information needed to create the new gene prediction
model is available. Once the model is created, model links are set up so all further builds of 
the assembly model will kick off a build of this gene prediction model. 
EOS
}

sub execute {
    my $self = shift;

    $self->status_message("Creating bacterial gene prediction model!");

    my $rv = $self->SUPER::_execute_body();
    unless ($rv) {
        $self->error_message("Could not create new model!");
        confess;
    }

    my $model = Genome::Model->get($self->result_model_id);
    unless ($model) {
        $self->error_message("Could not get newly created gene prediction model with ID " . $self->result_model_id);
        confess;
    }

    $self->status_message("Successfully created gene prediction model!");
    return 1;
}

sub type_specific_parameters_for_create {
    my $self = shift;
    return (
        assembly_contigs_file => $self->assembly_contigs_file,
        create_assembly_model => $self->create_assembly_model, 
        start_assembly_build => $self->start_assembly_build,
        assembly_processing_profile_name => $self->assembly_processing_profile_name,
        nr_database_location => $self->nr_database_location,
        acedb_version => $self->acedb_version,
        pipeline_version => $self->pipeline_version,
        assembly_version => $self->assembly_version,
        dev => $self->dev, 
        run_type => $self->run_type,
        project_type => $self->project_type,
    );
}
1;
