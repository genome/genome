package Genome::Model::Command::Define::GenePrediction::Eukaryotic;

use strict;
use warnings;
use Genome;
use Carp 'confess';

class Genome::Model::Command::Define::GenePrediction::Eukaryotic {
    is => 'Genome::Model::Command::Define::GenePrediction::Helper',
    has => [
        repeat_library => {
            is => 'Text',
            is_optional => 1,
            doc => 'Path to repeat library used by repeat masker',
        },
        snap_models => {
            is => 'Text',
            doc => 'Paths to prediction models to be used by snap, comma delimited',
        },
        fgenesh_model => {
            is => 'Path',
            doc => 'Path to prediction model to be used by fgenesh',
        },
    ],
};

sub help_brief {
    return "Create a new eukaryotic gene prediction model";
}

sub help_synopsis {
    return "Create a new eukaryotic gene prediction model";
}

sub help_detail {
    return <<"EOS"
Two things are needed to define a eukaryotic gene annotation model: a processing profile
name and a taxon ID. The taxon ID is used to find an assembly model that would correspond
to this prediction model. If such a model cannot be found, it will be created if the 
--create-assembly-model flag is set. If a model IS found and no successful build is found,
then a build is started if the --start-assembly-build flag is set.

With a successful assembly build, all the information needed to create the new gene prediction
model is available. Once the model is created, model links are set up so all further builds of 
the assembly model will kick off a build of this gene prediction model. 
EOS
}

sub execute {
    my $self = shift;
    
    $self->status_message("Creating eukaryotic gene prediction model!");

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

    $self->status_message("Successfully create gene prediction model!");
    return 1;
}

sub type_specific_parameters_for_create {
    my $self = shift;
    my %params = (
        assembly_contigs_file => $self->assembly_contigs_file,
        create_assembly_model => $self->create_assembly_model, 
        start_assembly_build => $self->start_assembly_build,
        assembly_processing_profile_name => $self->assembly_processing_profile_name,
        snap_models => $self->snap_models,
        fgenesh_model => $self->fgenesh_model,
    );
    $params{repeat_library} = $self->repeat_library if $self->repeat_library;
    return %params;
}

1;

