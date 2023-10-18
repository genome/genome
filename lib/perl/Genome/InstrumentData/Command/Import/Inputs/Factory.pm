package Genome::InstrumentData::Command::Import::Inputs::Factory;

use strict;
use warnings;

use Genome;

use Genome::InstrumentData::Command::Import::Inputs;
require File::Basename;
require List::MoreUtils;
use Params::Validate qw( :types );
use Text::CSV;
use Tie::File;
use Fcntl;

class Genome::InstrumentData::Command::Import::Inputs::Factory {
    has_optional => {
        analysis_project => {
            is => 'Genome::Config::AnalysisProject',
        },
        file => { is => 'Text', },
        process => {
            is => 'Genome::InstrumentData::Command::Import::Process',
        },
    },
    has_optional_calculated => {
        _increment_line_number => {
            calculate_from => [qw/ _line_number /],
            calculate => q| return $self->_line_number( ++$_line_number ); |,
        },
    },
    has_optional_transient => {
        _entity_attributes => { is => 'ARRAY', },
        _lines => { is => 'ARRAY', },
        _line_number => { is => 'Number', default => 0, },
        _parser => { is => 'Text::CSV', },
    },
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

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);
    return if not $self;

    $self->_load_file if $self->file;

    return $self;
}

sub from_inputs_id {
    my ($self, $id) = Params::Validate::validate_pos(@_, {isa => __PACKAGE__}, {type => SCALAR});

    my $inputs = UR::Object::get('Genome::InstrumentData::Command::Import::Inputs', $id);
    return $inputs if $inputs; # in cache

    my ($process_id, $line_number) = split(/\t/, $id, 2);
    $self->fatal_message('Failed to parse inputs id! %s', $id) if not defined $line_number;
    my $process = Genome::InstrumentData::Command::Import::Process->get($process_id);
    $self->fatal_message('Failed to get instdata import process for id: %s', $process_id) if not $process;
    $self->process($process);

    my $import_file = $process->import_file;
    $self->file($import_file);
    $self->_load_file($import_file);

    return $self->from_line_number($line_number);
}

sub entity_types {
    return (qw/ individual sample library instdata/);
}

sub resolve_sep_char_from_file_extension {
    my ($class, $file) = Params::Validate::validate_pos(@_, {isa => __PACKAGE__}, {type => SCALAR});

    return "\t" if $file eq '-';

    my ($dir, $basename, $ext) = File::Basename::fileparse($file, 'csv', 'tsv');
    die $class->error_message("Cannot determine type for file: %s. It needs to end with .csv or .tsv.", $file) if not $ext;

    return ( $ext eq 'csv' ? ',' : "\t" ),
}

sub _load_file {
    my $self = shift;

    my $file = $self->file;
    my $sep_char = $self->resolve_sep_char_from_file_extension($file);
    my $parser = Text::CSV->new({
            sep_char => $sep_char,
            empty_is_undef => 1,
        });
    die $self->error_message('Failed to create Text::CSV parser!') if not $parser;
    $self->_parser($parser);

    die $self->error_message('File (%s) is empty!', $file) if not -s $file;
    $self->file($file);
    tie my @lines, 'Tie::File', $file, mode => O_RDONLY
        or $self->fatal_message('Failed to load file %s: %s', $file, $!);
    $self->_lines(\@lines);
    $parser->parse($lines[0])
        or $self->fatal_message('Failed to parse header line! %s', $lines[0]);
    my @headers = $parser->fields;

    my $entity_attributes_ok = $self->_resolve_headers(\@headers);
    return if not $entity_attributes_ok;


    return $self;
}

sub next {
    my $self = shift;
    return $self->from_line_number( $self->_increment_line_number );
}

sub from_line_number {
    my ($self, $line_number) = Params::Validate::validate_pos(@_, {isa => __PACKAGE__}, {type => SCALAR});

    my $line = $self->_lines->[$line_number];
    return if not $line;

    my $entity_params = $self->_resolve_entity_params_from_line($line);
    $self->_resolve_names_for_entities($entity_params);

    # FIXME what if this was CSV and needs to be split on spaces?
    my $source_paths = [ split(',', delete $entity_params->{instdata}->{source_files}) ];

    return $self->from_params({
            line_number => $line_number,
            source_paths => $source_paths,
            entity_params => $entity_params,
        });
}

my $line_number = 1;
sub from_params {
    my ($self, $params) = Params::Validate::validate_pos(
        @_, {isa => __PACKAGE__},
        { entity_params => { type => HASHREF  }, },
    );

    # Set input id properties process_id and line_number
    my $process_id;
    if ( $self->process ) {
        $params->{process_id} = $self->process->id;
        $params->{entity_params}->{instdata}->{process_id} = $self->process->id;
    }
    elsif ( $self->file ) { # use md5 of file name
        $params->{process_id} = Genome::Sys->md5sum_data($self->file);
    }
    else {
        $params->{process_id} = $process_id++;
    }

    if ( not $params->{line_number} ) {
        $params->{line_number} = $line_number++;
    }
    
    if ( $params->{base_working_directory} ) {
        Genome::Sys->validate_directory_for_write_access( $params->{base_working_directory} );
    }

    # Check cache - get directly with UR::Object
    my $inputs = UR::Object::get(
        'Genome::InstrumentData::Command::Import::Inputs', 
        join("\t", $params->{process_id}, $params->{line_number}),
    );
    return $inputs if $inputs;

    # If not in cache, need source paths
    if ( not $params->{source_paths} ) {
        $self->fatal_message('No source paths given to create inputs! %s', Data::Dumper::Dumper($params));
    }
    $params->{entity_params}->{instdata}->{original_data_path} = join(',', @{$params->{source_paths}});

    # Add AnP
    if ( $self->analysis_project ) {
        $params->{analysis_project_id} = $self->analysis_project->id;
    }

    # Create using UR::Object, not the inputs class create, which redirects here
    return UR::Object::create('Genome::InstrumentData::Command::Import::Inputs', %$params);
}

sub _resolve_headers {
    my ($self, $headers) = Params::Validate::validate_pos(@_, {type => HASHREF}, {type => ARRAYREF});

    my @entity_attributes;
    for my $header ( @$headers ) {
        my ($type, $attribute) = split(/\./, $header, 2);
        if ( not defined $attribute ) {
            push @entity_attributes, undef; # add to indicate unrecognized field
            next;
        }
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

sub _resolve_entity_params_from_line {
    my ($self, $line) = Params::Validate::validate_pos(@_, {type => HASHREF}, {type => SCALAR});

    $self->_parser->parse($line)
        or $self->fatal_message('Failed to parse line! %s', $line);
    my @values = $self->_parser->fields;

    my %entity_params = map { $_ => {} } $self->entity_types;
    my $entity_attributes = $self->_entity_attributes;
    for ( my $i = 0; $i <= $#$entity_attributes; $i++ ) {
        my $entity_attribute = $entity_attributes->[$i];
        next if not $entity_attribute; # skip unrecognized fields
        my $value = $values[$i];
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
    my $ext = delete $entity_params->{library}->{ext} || 'extlibs';
    if ( not $entity_params->{library}->{name} ) {
        $entity_params->{library}->{name} = join('-', $sample_name, $ext);
    }

    return 1;
}

1;

