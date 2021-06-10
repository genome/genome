package Genome::InstrumentData::Command::Import::TrustedData;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Command::Import::TrustedData {
    is => ['Command::V2'],
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
        source_directory => {
            is => 'DirectoryPath',
            doc => 'Source directory to import.',
        },
        import_format => {
            is => 'Text',
            example_values => ['bam', 'fastq', 'sanger fastq', 'genotype microarray'],
            doc => 'The format of the instrument data. Setting certain values (e.g. "bam") will enable additional features in the GMS.',
        },
        read_count => {
            is => 'Number',
            doc => 'The read count for the instrument data',
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
        remove_source_files => {
            is => 'Boolean',
            default => 0,
            doc => 'By default this tool makes a copy of the instrument data files in an allocation.  This option will delete the original source files leaving only the allocated copy.',
        },
    ],
    has_optional_transient => [
        _new_instrument_data => {
            is => 'Genome::InstrumentData::Imported',
        },
    ],
    doc => 'import an instrument data directory as-is',
};


sub help_detail {
    return <<HELP;
Import a directory of sequence file(s).  The directory will be copied as-is with no validation of the data format.

Instrument Data Properties
 Name and value pairs to add to the instrument data. Separate name and value with an equals (=)
  and name/value pairs with a comma (,).

  e.g., 'read_length=100,sequencing_platform=solexa,run_name=TEST12345-6-AAAAAA'

  If the data has already been aligned, the reference_sequence_build_id property can be used to link the correct reference.

HELP
}

sub execute {
    my $self = shift;

    my $anp = $self->analysis_project;
    my $guard = $anp->set_env;

    my $source = $self->source_directory;
    unless (-e $source) {
        $self->fatal_message('Source directory %s not found.', $source);
    }

    my %properties;
    for ($self->instrument_data_properties) {
        my ($key, $value) = split '=', $_;
        $properties{$key} = $value;
    }
    $properties{subset_name} //= 'unknown';
    $properties{sequencing_platform} //= 'solexa';

    $properties{description} = $self->description if $self->description;

    my $entity = Genome::InstrumentData::Imported->create(
        %properties,
        library => $self->library,
        import_format => $self->import_format,
        import_source_name => $self->import_source_name,
        original_data_path => $self->source_directory,
        read_count => $self->read_count,
    );
    unless ($entity) {
        $self->fatal_message('Could not instantiate new instrument data.');
    }

    my $allocation = Genome::Disk::Allocation->create(
        allocation_path => join('/', 'instrument_data', 'imported', $entity->id),
        disk_group_name => Genome::Config::get('disk_group_alignments'),
        kilobytes_requested => $entity->calculate_alignment_estimated_kb_usage,
        owner_class_name => $entity->class,
        owner_id => $entity->id,
    );
    unless ($allocation) {
        $self->fatal_message('Failed to create allocation for instrument data.');
    }

    Genome::Sys->rsync_directory(
        source_directory => $self->source_directory,
        target_directory => $allocation->absolute_path,
        chmod => 'Dug=rx,Fug=r',
        chown => ':' . Genome::Config::get('sys_group'),
        remove_source_files => $self->remove_source_files,
    );

    Genome::Config::AnalysisProject::InstrumentDataBridge->create(
        instrument_data => $entity,
        analysis_project => $anp,
    );

    $self->_new_instrument_data( $entity );


    return 1;
}

1;
