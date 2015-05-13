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
    has_optional_input => [
        description  => {
            is => 'Text',
            doc => 'Description of the data.',
        },
        instrument_data_properties => {
            is => 'Text',
            is_many => 1,
            doc => 'Name and value pairs to add to the instrument data. Separate name and value with an equals (=) and name/value pairs with a comma (,).',
        },
    ],
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
    my $run = Genome::InstrumentData::Command::Import::WorkFlow::Run->create(
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

    my @instrument_data_properties = $self->instrument_data_properties;
    push @instrument_data_properties, 'description='.$self->description if defined $self->description;		
    push @instrument_data_properties, 'downsample_ratio='.$self->downsample_ratio if defined $self->downsample_ratio;
    return Genome::InstrumentData::Command::Import::WorkFlow::Inputs->create(
        analysis_project => $self->analysis_project,
        library => $self->library,
        source_files => [ $self->source_files ],
        instrument_data_properties => \@instrument_data_properties,
    );
}

1;

