package Genome::InstrumentData::Command::Import::CreateEntities;

use strict;
use warnings;

use Genome;

use Genome::InstrumentData::Command::Import::Inputs::Factory;
require File::Basename;
require List::MoreUtils;
use Params::Validate qw( :types );

class Genome::InstrumentData::Command::Import::CreateEntities {
    is => 'Command::V2',
    has_input => {
        file => {
            is => 'Text',
            doc => Genome::InstrumentData::Command::Import::Inputs::Factory->csv_help,
        },
    },
    has_optional_transient => {
        _names_seen => { is => 'HASH', default_value => {}, },
    },
    doc => 'create individuals, samples, and libraries that are necessary to import instrument data',
};

sub help_detail {
    return 'Using a modified metadata spreadsheet saved as a comma/tab separated file, create the necessary entites (individual, sample, library) to import instrument data.';
}

sub execute {
    my $self = shift;
    
    my $inputs_factory = Genome::InstrumentData::Command::Import::Inputs::Factory->create(
        file => $self->file,
    );
    while ( my $inputs = $inputs_factory->next ) {
        my $entity_params = $inputs->entity_params;
        my $library = Genome::Library->get(name => $entity_params->{library}->{name});
        next if $library;

        $self->_create_individual_if_needed($entity_params);
        $self->_create_sample_if_needed($entity_params);
        $self->_create_library_if_needed($entity_params);
    }

    return 1;
}

sub _create_individual_if_needed {
    my ($self, $entity_params) = Params::Validate::validate_pos(@_, {type => HASHREF}, {type => HASHREF});
    my %params = %{$entity_params->{individual}};
    return 1 if $self->_does_enitity_exist('individual', $params{name});
    my $taxon_name = $params{taxon};
    die $self->error_message('No taxon for %s!', $params{name}) if not $taxon_name;
    my $taxon = Genome::Taxon->get(name => $taxon_name);
    die $self->error_message('Taxon does not exist for %s!', $taxon_name) if not $taxon;
    $params{taxon} = $taxon;
    return $self->_create_entity('individual', \%params);
}

sub _create_sample_if_needed {
    my ($self, $entity_params) = Params::Validate::validate_pos(@_, {type => HASHREF}, {type => HASHREF});
    my $params = $entity_params->{sample};
    return 1 if $self->_does_enitity_exist('sample', $params->{name});
    die $self->error_message('Sample extraction_type is required to create a sample!') if not $params->{extraction_type};
    my $source = Genome::Individual->get(name => $entity_params->{individual}->{name});
    die $self->error_message('No source for %s',  $entity_params->{individual}->{name}) if not $source;
    $params->{source} = $source; # add the individual name to the params
    return $self->_create_entity('sample', $params);
}

sub _create_library_if_needed {
    my ($self, $entity_params) = Params::Validate::validate_pos(@_, {type => HASHREF}, {type => HASHREF});
    my $params = $entity_params->{library};
    return 1 if $self->_does_enitity_exist('library', $params->{name});
    my $sample = Genome::Sample->get(name => $entity_params->{sample}->{name});
    die $self->error_message('No sample for %s',  $entity_params->{sample}->{name}) if not $sample;
    $params->{sample} = $sample;
    return $self->_create_entity('library', $params);
}

sub _does_enitity_exist {
    my ($self, $type, $name) = @_;
    return 1 if $self->_names_seen->{$type}->{$name};
    my $entity_class = 'Genome::'.ucfirst($type);
    return $entity_class->get(name => $name);
}

sub _create_entity {
    my ($self, $type, $params, $dependent_params) = Params::Validate::validate_pos(
        @_, {type => HASHREF}, {type => SCALAR}, {type => HASHREF},
    );

    my $entity_class = 'Genome::'.ucfirst($type);
    my $entity = $entity_class->create(%$params);
    if ( not $entity or my @errors = $entity->__errors__ ) {
        for ( @errors ) { $self->error_message($_->__display_name__); }
        die $self->error_message('Problems creating %s for %s.', $type, $params->{name});
    }

    $self->status_message('CREATE %s %s', $type, $params->{name});
    return 1;
}

1;

