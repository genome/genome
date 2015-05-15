package Genome::InstrumentData::Command::Import::CsvParser;

use strict;
use warnings;

use Genome;

require File::Basename;
require List::MoreUtils;
use Params::Validate qw( :types );
use Text::CSV;

class Genome::InstrumentData::Command::Import::CsvParser {
    is => 'UR::Object',
    has => {
        file => {
            is => 'Text',
            doc => 'Comma (.csv) or tab (.tsv) separated file of entity names, attributes and other metadata. Separator is determined by file extension.',
        },
    },
    has_optional_transient => {
        _entity_attributes => { is => 'ARRAY', },
        _fh => { },
        _parser => { },
    }
};

sub csv_help {
    return <<HELP;
This is a comma (.csv) or tab (.tsv) separated file of entity names, attributes and other metadata. Separator is determined by file extension. Column headers should start with the entity (individual, sample, library, instdata) name then a period (.) and then then attribute name (Ex: sample.name_part). Here are some required and optional columns. For more, see each entity's create command (Ex: genome sample create --h). Please see Confluence documentation for more information and a full example.

Individual\n
 Required
  individual.name_part      => Name or id from external source.
   OR
  individual.name           => Full individual name. Use when the name is desired to have a different value than being derived from the sample/library name.

  individual.taxon          => Species name of the taxon.
 Optional
  individual.upn            => External name/identifier. Often the second part of the new sample name.
  individual.common_name        => Usually the project name plus a number

Sample\n
 Required
  sample.name               => The full sample name. The individual and library names will be dervied from it, unless they are given. 
   OR
  sample.name_part          => Name or id from external source. If name is given, the individual and library names will be resolved from it.
  sample.extraction_type    => 'genomic dna' or 'rna'

 Optional, but recommended:
  sample.common_name        => Usually normal or tumor to indicate disease state.

Library\n
 Optional
  library.ext               => Extension to append to the sample name. Default is 'extlibs'.
 
Instrument Data (needed for generating source-files.tsv)
 Required
  instdata.source_files     => Local copy of the source files to import.
 Optional
  instdata.run_name         => The run name or id.
  instdata.flow_cell_id     => The flow cell id of the run.
  instdata.lane             => The lane of the run.
  instdata.index_sequence   => The run name or id.

HELP
}

sub entity_types {
    return (qw/ individual sample library instdata/);
}

sub resolve_sep_char_from_file_extension {
    my ($class, $file) = Params::Validate::validate_pos(@_, {isa => __PACKAGE__}, {type => SCALAR});

    my ($dir, $basename, $ext) = File::Basename::fileparse($file, 'csv', 'tsv');
    die $class->error_message("Cannot determine type for file: %s. It needs to end with .csv or .tsv.", $file) if not $ext;

    return ( $ext eq 'csv' ? ',' : "\t" ),
}

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);
    return if not $self;

    my $file = $self->file;
    my $sep_char = $self->resolve_sep_char_from_file_extension($file);
    my $parser = Text::CSV->new({
            sep_char => $sep_char,
            empty_is_undef => 1,
        });
    die $self->error_message('Failed to create Text::CSV parser!') if not $parser;
    $self->_parser($parser);

    die $self->error_message('File (%s) is empty!', $file) if not -s $file;
    my $fh = Genome::Sys->open_file_for_reading($file);
    $self->_fh($fh);
    my $headers = $parser->getline($fh);
    $parser->column_names($headers);

    my $entity_attributes_ok = $self->_resolve_headers($headers);
    return if not $entity_attributes_ok;

    return $self;
}

sub next {
    my $self = shift;

    my $line_ref = $self->_parser->getline_hr($self->_fh);
    return if not $line_ref;

    my $entity_params = $self->_resolve_entity_params_for_values($line_ref);
    $self->_resolve_names_for_entities($entity_params);

    return $entity_params;
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

    # set sample name if library name given and sample name not given
    if ( $entity_params->{library}->{name} and  not $entity_params->{sample}->{name} ) {
        my @tokens = split(/\-/, $entity_params->{library}->{name});
        pop @tokens; #rm lib ext
        $entity_params->{sample}->{name} = join('-', @tokens);
    }

    my $nomenclature = delete $entity_params->{sample}->{nomenclature};
    my $individual_name_part = delete $entity_params->{individual}->{name_part};
    my $sample_name_part = delete $entity_params->{sample}->{name_part};

    # sample - use name to fill in as needed
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
    $entity_params->{sample}->{nomenclature} = $nomenclature;

    # individual
    my $individual_name = $entity_params->{individual}->{name};
    if ( $individual_name ) {
        die $self->error_message('Invalid individual name: %s. It must include the first part of the sample name: %s.', $individual_name, $sample_name) if $sample_name !~ /^$individual_name/;
    }
    else {
        $entity_params->{individual}->{name} = join('-', $nomenclature, $individual_name_part);
    }
    $entity_params->{individual}->{nomenclature} = $nomenclature;
    $entity_params->{individual}->{upn} = $individual_name_part if not $entity_params->{individual}->{upn};
    
    # library name - add ext or use default if name not given/set
    if ( not $entity_params->{library}->{name} ) {
        $entity_params->{library}->{name} = $sample_name.( 
            $entity_params->{library}->{ext} ? $entity_params->{library}->{ext} : '-extlibs'
        );
    }

    return 1;
}

sub _generate_output {
    my $self = shift;
    my @output = @{$self->_output};
    return 1 if not @output;
    return Genome::Sys->write_file($self->output_file, $self->_output_header, @output);
}

1;

