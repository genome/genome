package Genome::InstrumentData::Command::Import::Basic;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Command::Import::Basic { 
    is => [qw/ Command::V2 Genome::Model::Tools::Picard::WithDownsampleRatio /],
    has_input => [
        analysis_project => {
            is => 'Genome::Config::AnalysisProject',
            doc => 'Analysis project to assign to the created instrument data.',
        },
        import_source_name => {
            is => 'Text',
            doc => "Organization or site name/abbreviation from where the source was generated or downloaded.",
        },
        library => {
            is => 'Genome::Library',
            doc => 'Library to use.  It must exist.',
        },
        source_files => {
            is => 'Text',
            is_many => 1,
            doc => 'Source files to import. If importing fastqs, put the file containing the forward [read 1] reads first.',
        },
    ],
    has_optional_input => {
        description  => {
            is => 'Text',
            doc => 'Description of the data.',
        },
        instrument_data_properties => {
            is => 'Text',
            is_many => 1,
            doc => 'Name and value pairs to add to the instrument data. Separate name and value with an equals (=) and name/value pairs with a comma (,).',
        },
        base_working_directory => {
            is => 'Text',
            doc => 'Base working directory to use when running the import. A temporary directory will be made inside this directory, and removed when the import fails or completes.',
        },
    },
    has_optional_transient => [
        _new_instrument_data => {},
    ],
    doc => 'import sequence files as instrument data into GMS',
};

sub help_detail {
    return <<HELP;
Import sequence files into GMS. All files will be converted to SAM format and stored as BAM.

Source Files 
 Types       Notes
  FASTQ       Can be remote, tar'd and/or gzipped.
  BAM, SAM    Will be split by read group.
  SRA         Aligned and unaligned reads will be dumped. SRAs are known to produce unreliable BAM files.

Instrument Data Properties
 Name and value pairs to add to the instrument data. Separate name and value with an equals (=)
  and name/value pairs with a comma (,).
  
  Example...set flow_cell_id to 'AXXAX' and the index sequence to 'AATTGGCC' on the created instrument
   data entities:

  flow_cell_id=AXXAX,index_sequence=AATTGGCC

HELP
}

sub execute {
    my $self = shift;

    my $work_flow_inputs = $self->_resolve_work_flow_inputs;

    my $anp = $self->analysis_project;
    my $config_dir = $anp->environment_config_dir;
    unless ($config_dir) {
        $self->fatal_message('No analysis project environment configuration defined for %s!', $anp->__display_name__);
    }
    local $ENV{XGENOME_CONFIG_PROJECT_DIR} = $config_dir;

    my $run = Genome::InstrumentData::Command::Import::WorkFlow::ImportInstData->create(
        work_flow_inputs => $work_flow_inputs,
    );
    die $self->error_message(
        'Source files (%s) have existing instrument data (%s). Cannot reimport!', 
        $work_flow_inputs->source_files->original_data_path,
        join(' ', map { $_->id } $run->instrument_data),
    ) if $run->shortcut;
    $run->execute;
    die 'Failed to run importer!' if not $run or not $run->result;

    $self->_new_instrument_data([ $run->instrument_data ]);
    return 1;
}

sub _resolve_work_flow_inputs {
    my $self = shift;

    my $factory = Genome::InstrumentData::Command::Import::Inputs::Factory->create(
        analysis_project => $self->analysis_project,
    );
    return $factory->from_params({
            base_working_directory => $self->base_working_directory,
            entity_params => {
                library => { id => $self->library->id, },
                instdata => $self->_resolve_instrument_data_properties,
            },
            source_paths => [ $self->source_files ],
        });
}

sub _resolve_instrument_data_properties {
    my $self = shift;

    my @incoming_properties = $self->instrument_data_properties;
    push @incoming_properties, 'description='.$self->description if defined $self->description;		
    push @incoming_properties, 'downsample_ratio='.$self->downsample_ratio if defined $self->downsample_ratio;

    my %properties;
    return \%properties if not @incoming_properties;

    for my $key_value_pair ( @incoming_properties ) {
        my ($label, $value) = split('=', $key_value_pair);
        if ( not defined $value or $value eq '' ) {
            die $self->error_message('Failed to parse with instrument data property label/value! '.$key_value_pair);
        }
        if ( exists $properties{$label} and $value ne $properties{$label} ) {
            die $self->error_message(
                "Multiple values for instrument data property! $label => ".join(', ', sort $value, $properties{$label})
            );
        }
        $properties{$label} = $value;
    }

    return \%properties;
}

1;

