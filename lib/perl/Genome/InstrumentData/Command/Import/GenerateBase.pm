package Genome::InstrumentData::Command::Import::GenerateBase;

use strict;
use warnings;

use Genome;

require File::Basename;
require List::MoreUtils;
use Params::Validate qw( :types );
use Text::CSV;

class Genome::InstrumentData::Command::Import::GenerateBase {
    is => 'Command::V2',
    is_abstract => 1,
    has_input => {
        file => {
            is => 'Text',
            doc => 'Comma (.csv) or tab (.tsv) separated file of entity names, attributes and other meta data. Separartor is determined by file extension.',
        },
    },
    has_optional_output => {
        output_file => {
            is => 'Text',
            default_value => '-',
            doc => 'Output file to put the commands to create the needed entities for import.',
        },
    },
    has_optional_transient => {
        _output => { is => 'ARRAY', default_value => [], },
        _instdata_property_names => { is => 'ARRAY', },
        _entity_attributes => { is => 'ARRAY', },
        _nomenclature => { is => 'Text', },
        _names_seen => { is => 'HASH', default_value => {}, },
    }
};

sub help_detail {
    return <<HELP;
The file should be a comma or tab separated values and indicated with the appropriate extension (csv and tsv). Column headers to use to generate the create commands should start with the entity (individual, sample, library, instdata) name then a period (.) and then then attribute name (Ex: sample.name_part). Here are some required and optional columns. For more, see each entity's create command (Ex: genome sample create --h).

Individual:
 Required:
  individual.name or name_part => Name or id from external source.
  individual.taxon => Species name of the taxon.

Sample:
 Required:
  sample.name or sample.name_part => Name or id from external source. If name is given, the individual and library names will be resolved from it.

 Optional, but recommended:
  sample.common_name => Usually normal or tumor to indicate disease state.

Library:
 Optional:
  library.ext   => Extension to append to the sample name. Deault is 'extlibs'.
 
Instrument Data (only for generating source-files.tsv)
 Required:
  instdata.source_files  => Local copy of the source files to import.
  instdata.run_name      => The run name or id.

HELP
}

sub entity_types {
    return (qw/ individual sample library instdata/);
}

sub _open_file_parser {
    my $self = shift;

    my $file = $self->file;
    my ($dir, $basename, $ext) = File::Basename::fileparse($file, 'csv', 'tsv');
    die $self->error_message("Cannot determine type for file: %s. It needs to end with .csv or .tsv.", $file) if not $ext;
    my $parser = Text::CSV->new({
            sep_char => ( $ext eq 'csv' ? ',' : "\t" ),
            empty_is_undef => 1,
        });
    die $self->error_message('Failed to create Text::CSV parser!') if not $parser;

    my $fh = Genome::Sys->open_file_for_reading($file);
    my $headers = $parser->getline($fh);
    die $self->error_message('File (%s) does not have any headers! Is it empty?', $file) if not $headers;
    $parser->column_names($headers);

    my $entity_attributes_ok = $self->_resolve_headers($headers);
    return if not $entity_attributes_ok;

    return sub{ return $parser->getline_hr($fh); };
}

sub _resolve_headers {
    my ($self, $headers) = Params::Validate::validate_pos(@_, {type => HASHREF}, {type => ARRAYREF});

    my @entity_attributes;
    for my $header ( @$headers ) {
        my ($type, $attribute) = split(/\./, $header, 2);
        next if not defined $attribute; # silently skip columns w/o an entity type
        die $self->error_message('Invalid entity type: %s', $type) if not List::MoreUtils::any { $type eq $_ } $self->entity_types;
        push @entity_attributes, {
            header => $header,
            type => $type,
            attribute => $attribute,
        };
    }
    $self->_entity_attributes(\@entity_attributes);

    return 1;
}

sub _resolve_entity_params_for_values {
    my ($self, $line_ref) = Params::Validate::validate_pos(@_, {type => HASHREF}, {type => HASHREF});

    my %entity_params = map { $_ => {} } $self->entity_types;
    for my $entity_attribute ( @{$self->_entity_attributes} ) {
        my $value = $line_ref->{ $entity_attribute->{header} };
        next if not defined $value;
        $entity_params{ $entity_attribute->{type} }->{ $entity_attribute->{attribute} } = $value;
    }

    return \%entity_params;
}

sub _resolve_names_for_entities {
    my ($self, $entity_params) = Params::Validate::validate_pos(@_, {type => HASHREF}, {type => HASHREF});

    my $nomenclature = delete $entity_params->{sample}->{nomenclature};
    my $individual_name_part = delete $entity_params->{individual}->{name_part};
    my $sample_name_part = delete $entity_params->{sample}->{name_part};

    # sample
    my $sample_name = $entity_params->{sample}->{name};
    if ( $sample_name ) {
        my @tokens = split(/\-/, $sample_name);
        if ( @tokens <  3 ) {
            die $self->error_message('Invalid sample name: %s. It must have at least 3 parts separated by dashes.', $sample_name);
        }
        ($nomenclature, $individual_name_part) = @tokens[0..1];
    }
    else {
        die $self->error_message('No sample.nomenclature column given! It is required to resolve entity names when no sample name is given.') if not $nomenclature;
        die $self->error_message('No sample.name_part column_given! It is required to resolve entity names when no sample name is given.') if not $sample_name_part;
        die $self->error_message('No individual.name_part column_given! It is required to resolve entity names when no sample name is given.') if not $individual_name_part;
        $sample_name = join('-', $nomenclature, $individual_name_part, $sample_name_part);
        $entity_params->{sample}->{name} = $sample_name;
    }

    # individual
    my $individual_name = $entity_params->{individual}->{name};
    if ( $individual_name ) {
        die $self->error_message('Invalid individual name: %s. It must include the first part of the sample name: %s.', $individual_name, $sample_name) if $sample_name !~ /^$individual_name/;
    }
    else {
        $entity_params->{individual}->{name} = join('-', $nomenclature, $individual_name_part);
    }
    $entity_params->{individual}->{upn} = $individual_name_part if not $entity_params->{individual}->{upn};
    
    # library - add ext or use default
    $entity_params->{library}->{name} = $sample_name.( 
        $entity_params->{library}->{ext} ? $entity_params->{library}->{ext} : '-extlibs'
    );

    # nomenclature
    $self->_nomenclature($nomenclature);

    return 1;
}

sub _generate_output {
    my $self = shift;
    my @output = @{$self->_output};
    return 1 if not @output;
    return Genome::Sys->write_file($self->output_file, $self->_output_header, @output);
}

1;

