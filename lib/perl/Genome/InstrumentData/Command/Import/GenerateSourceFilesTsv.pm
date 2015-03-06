package Genome::InstrumentData::Command::Import::GenerateSourceFilesTsv;

use strict;
use warnings;

use Genome;

require File::Basename;
require List::MoreUtils;
use Params::Validate qw( :types );
use Text::CSV;

class Genome::InstrumentData::Command::Import::GenerateSourceFilesTsv {
    is => 'Genome::InstrumentData::Command::Import::GenerateBase',
    doc => 'generate source files tsv to be used in the manager to import instrument data',
    has_optional_transient => {
         _instdata_property_names => { is => 'ARRAY', },
    },
};

sub help_detail {
    return;
}

sub _output_header {
    my $self = shift;
    return join("\t", 'library_name', @{$self->_instdata_property_names})."\n";
}

sub execute {
    my $self = shift;
    
    my $parser = $self->_open_file_parser;

    my @instdata_property_names = sort map { $_->{attribute} } grep { $_->{type} eq 'instdata' } @{$self->_entity_attributes};
    die 'No source_files attribute for instdata in file: '.$self->output_file if not List::MoreUtils::any { $_ eq 'source_files' } @instdata_property_names;
    $self->_instdata_property_names(\@instdata_property_names);

    while ( my $line_ref = $parser->() ) {
        my $entity_params = $self->_resolve_entity_params_for_values($line_ref);
        $self->_resolve_names_for_entities($entity_params);
        $self->_add_instdata_source_files_line($entity_params);
    }

    $self->_generate_output;

    return 1;
}

sub _add_instdata_source_files_line {
    my ($self, $entity_params) = Params::Validate::validate_pos(@_, {type => HASHREF}, {type => HASHREF});

    die $self->error_message('No instrument data source files specified for library: '.$entity_params->{library}->{name}) if not $entity_params->{instdata}->{source_files};

    push @{$self->_output}, join(
        "\t",
        $entity_params->{library}->{name}, 
        map { $entity_params->{instdata}->{$_} } @{$self->_instdata_property_names},
    )."\n";

    return 1;
}

1;

