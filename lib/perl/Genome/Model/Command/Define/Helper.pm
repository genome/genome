package Genome::Model::Command::Define::Helper;

# TODO: inherit from ::MinimalBase, which presumes less about instdata, profiles, subjects

use strict;
use warnings;

use Genome;
use File::Path;
use Carp 'confess';

class Genome::Model::Command::Define::Helper {
    is => 'Genome::Model::Command::Define::BaseMinimal',
    is_abstract => 1,
    has_input => [
        processing_profile => {
            is => 'Genome::ProcessingProfile',
            id_by => 'processing_profile_id',
            doc => 'Processing profile to be used by model, can provide either a name or an ID',
        },
    ],
    has_optional_input => [
        subject => {
            is => 'Genome::Subject',
            id_by => 'subject_id',
            doc => 'Subject for the model, can provide either a name or an id. If instrument data is provided and this is not, ' .
                'an attempt will be made to resolve it based on the provided instrument data'
        },
        model_name => {
            is => 'Text',
            doc => 'User meaningful name for this model, a default is used if none is provided',
        },
        result_model_id => {
            is => 'Number',
            is_transient => 1,
            doc => 'Stores the ID of the newly created model, useful when running this command from a script',
        },
        run_as => {
            is => 'Text',
            doc => 'Specify who the model should run_as if run in production'
        },
    ],
    has_optional_param => [
        auto_assign_inst_data => {
            is => 'Boolean',
            default_value => 0,
            doc => 'Assigning instrument data to the model is performed automatically',
        },
        auto_build_alignments => {
            is => 'Boolean',
            default_value => 1,
            doc => 'The building of the model is performed automatically',
        },
    ],
    has_many_optional_input => [
        instrument_data => {
            is => 'Genome::InstrumentData',
            doc => 'Instrument data to be assigned to the model, can provide a query to resolve, a list of ids, etc'
        },
        groups => {
            is => 'Genome::ModelGroup',
            doc => 'Model groups to put the newly created model into',
        },        
    ],
};

sub _resolve_subject_from_inputs {
    my $self = shift;
    my $subject = $self->_resolve_subject_from_instrument_data;
    return $subject;
}

sub _resolve_subject_from_instrument_data {
    my $self = shift;

    my @instrument_data = $self->instrument_data;
    unless (@instrument_data) {
        return; 
    }

    my @samples = map { $_->sample } @instrument_data;
    unless (grep { $_->id ne $samples[0]->id } @samples) {
        return $samples[0];
    }

    my @individuals = map { $_->source } @samples;
    unless (grep { $_->id ne $individuals[0]->id } @individuals) {
        return $individuals[0];
    }

    my @taxons = map { $_->taxon } @individuals;
    unless (grep { $_->id ne $taxons[0]->id } @taxons) {
        return $taxons[0];
    }

    return;
}

sub type_specific_parameters_for_create {
    my $self = shift;
    my @params = (
        processing_profile_id => $self->processing_profile->id,
        name => $self->model_name,
        auto_assign_inst_data => $self->auto_assign_inst_data,
        auto_build_alignments => $self->auto_build_alignments,
    );
    push @params, $self->SUPER::type_specific_parameters_for_create(@_); 
    push @params, 'instrument_data' => [$self->instrument_data] if $self->instrument_data;
    push @params, 'model_groups' => [$self->groups] if $self->groups;
    #push @params, 'projects' => [$self->projects] if $self->projects;
    return @params;
}

1;


