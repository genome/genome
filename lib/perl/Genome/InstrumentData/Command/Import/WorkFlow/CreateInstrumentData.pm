package Genome::InstrumentData::Command::Import::WorkFlow::CreateInstrumentData;

use strict;
use warnings;

use Genome;

require Cwd;
require File::Basename;

class Genome::InstrumentData::Command::Import::WorkFlow::CreateInstrumentData { 
    is => 'Command::V2',
    has_input => [
        source_files => {
            is => 'Text',
            is_many => 1,
            doc => 'Source files to import. If importing multiple files, put the file containing the forward reads first.',
        },
        sample => {
            is => 'Genome::Sample',
            doc => 'Sample to use. The external library for the instrument data will be gotten or created.',
        },
        instrument_data_properties => {
            is => 'Text',
            is_many => 1,
            doc => 'Name and value pairs to add to the instrument data. Separate name and value with an equals (=) and name/value pairs with a comma (,).',
        },
    ],
    has_output => [
        instrument_data => {
            is => 'Genome::InstrumentData',
            doc => 'The instrument data entity to be imported.',
        },
    ],
    has_transient_optional => [
        library => { is => 'Genome::Library', },
        original_format => { is => 'Text', },
        import_format => { is => 'Text', },
        kilobytes_requested => { is => 'Number', },
    ],
};

sub execute {
    my $self = shift;
    $self->status_message('Create instrument data...');

    my $library = $self->_resvolve_library;
    return if not $library;

    my $validate_source_files = $self->_validate_source_files;
    return if not $validate_source_files;

    my $instrument_data = $self->_create_instrument_data;
    return if not $instrument_data;

    $self->status_message('Create instrument data...done');
    return 1;
}

sub _resvolve_library {
    my $self = shift;
    $self->status_message('Resolve library...');

    my $sample = $self->sample;
    $self->status_message('Sample name: '.$sample->name);
    $self->status_message('Sample id: '.$sample->id);
    my $library_name = $sample->name.'-extlibs';
    $self->status_message('Library name: '.$library_name);
    my $library = Genome::Library->get(
        name => $library_name,
        sample => $sample,
    );
    if ( not $library ) {
        $library = Genome::Library->create(
            name => $library_name,
            sample => $sample,
        );
        if ( not $library ) {
            $self->error_message('Failed to get or create external library for sample! '.$sample->id);
            return;
        }
    }
    $self->status_message('Library id: '.$library->id);

    $self->status_message('Resolve library...done');
    return $self->library($library);
}

#< Validate Source Files >#
sub _validate_source_files {
    my $self = shift;
    $self->status_message('Validate source files...');

    my @source_files = $self->source_files;
    for my $source_file ( @source_files ) { $self->status_message("Source file(s): $source_file"); }
    $self->status_message("Source file count: ".@source_files);

    my $resolve_formats_ok = $self->_resolve_start_and_import_format(@source_files);
    return if not $resolve_formats_ok;

    my $max_source_files = ( $self->original_format =~ /fast[aq]/ ? 2 : 1 );
    if ( @source_files > $max_source_files ) {
        $self->error_message("Cannot handle more than $max_source_files source files!");
        return;
    }

    my $kilobytes_requested = $self->_resolve_kilobytes_requested(@source_files);
    return if not $kilobytes_requested;
    $kilobytes_requested += 51_200; # a little extra
    $self->kilobytes_requested($kilobytes_requested);
    $self->status_message('Kilobytes requested: '.$self->kilobytes_requested);

    $self->status_message('Validate source files...done');
    return 1;
}

sub _resolve_start_and_import_format {
    my ($self, @source_files) = @_;
    $self->status_message('Resolve start and import format...');

    my %suffixes;
    for my $source_file ( @source_files ) {
        $source_file =~ s/\.gz$//;
        my ($suffix) = $source_file =~ /\.(\w+)$/;
        if ( not $suffix ) {
            $self->error_message("Failed to get suffix from source file! $source_file");
            return;
        }
        $suffixes{$suffix}++;
    }

    my %suffixes_to_original_format = (
        txt => 'fastq',
        fastq => 'fastq',
        fq => 'fastq',
        #fasta => 'fasta',
        bam => 'bam',
        sra => 'sra',
    );
    my %formats;
    for my $suffix ( keys %suffixes ) {
        if ( not exists $suffixes_to_original_format{$suffix} ) {
            $self->error_message('Unrecognized suffix! '.$suffix);
            return;
        }
        $formats{ $suffixes_to_original_format{$suffix} } = 1;
    }

    my @formats = keys %formats;
    if ( @formats > 1 ) {
        $self->error_message('Got more than one format when trying to determine format!');
        return;
    }
    my $original_format = $formats[0];
    $self->original_format($original_format);
    $self->status_message('Start format: '.$self->original_format);

    my %original_format_to_import_format = ( # temp, as everything will soon be bam
        fastq => 'sanger fastq',
        #fastq => 'bam',
        bam => 'bam',
        sra => 'bam',
    );
    $self->import_format( $original_format_to_import_format{$original_format} );
    $self->status_message('Import format: '.$self->import_format);

    $self->status_message('Resolve start and import format...done');
    return 1;
}

sub _resolve_kilobytes_requested {
    my ($self, @source_files) = @_;

    my $kilobytes_requested;
    for my $source_file ( @source_files ) {
        my $size = -s $source_file;
        if ( not $size ) {
            $self->error_message("Source file does not exist! $source_file");
            return;
        }
        $size = int( $size / 1024 );
        $size *= 3 if $source_file =~ /\.gz$/; # assume ~30% compression rate for gzipped fasta/q
        $kilobytes_requested += $size;
    }

    my $multiplier = ( $self->original_format =~ /bam/ ? 3 : 2 );
    return $kilobytes_requested * $multiplier;
}

sub _create_instrument_data {
    my $self = shift;

    $self->status_message('Checking if source files were previously imported...');
    my %properties = (
        library => $self->library,
        original_data_path => join(',',  $self->source_files),
    );
    my $instrument_data = Genome::InstrumentData::Imported->get(%properties);
    if ( $instrument_data ) {
        $self->error_message('Found existing instrument data for library and source files. Were these previously imported? Exiting instrument data id: '.$instrument_data->id.', source files: '.$properties{original_data_path});
        return;
    }
    $self->status_message('Source files were NOT previously imported!');

    $self->status_message('Create instrument data...');
    $properties{import_format} = 'bam';
    $properties{sequencing_platform} = 'solexa';
    for my $name_value ( $self->instrument_data_properties ) {
        my ($name, $value) = split('=', $name_value);
        if ( not defined $value or $value eq '' ) {
            $self->error_message('Failed to parse with instrument data property name/value! '.$name_value);
            return;
        }
        if ( exists $properties{$name} and $value ne $properties{$name} ) {
            $self->error_message(
                "Multiple values for instrument data property! $name => ".join(', ', sort $value, $properties{$name})
            );
            return;
        }
        $properties{$name} = $value;
    }

    $instrument_data = Genome::InstrumentData::Imported->create(%properties);
    if ( not $instrument_data ) {
        $self->error_message('Failed to create instrument data!');
        return;
    }
    $self->status_message('Instrument data id: '.$instrument_data->id);
    $instrument_data->add_attribute(attribute_label => 'original_format', attribute_value => $self->original_format);

    my $allocation = Genome::Disk::Allocation->create(
        disk_group_name => 'info_alignments',
        allocation_path => 'instrument_data/imported/'.$instrument_data->id,
        kilobytes_requested => $self->kilobytes_requested,
        owner_class_name => $instrument_data->class,
        owner_id => $instrument_data->id,
    );
    if ( not $allocation ) {
        $self->error_message('Failed to create allocation for instrument data! '.$instrument_data->id);
        return;
    }
    $self->status_message('Allocation id: '.$allocation->id);
    $self->status_message('Allocation path: '.$allocation->absolute_path);
    
    my $tmp_dir = $allocation->absolute_path.'/tmp';
    Genome::Sys->create_directory($tmp_dir);
    $self->status_message('Allocation tmp path: '.$tmp_dir);

    $self->status_message('Create instrument data...done');
    return $self->instrument_data($instrument_data);
}

1;

