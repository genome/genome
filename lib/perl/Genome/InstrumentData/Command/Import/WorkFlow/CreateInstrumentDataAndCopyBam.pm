package Genome::InstrumentData::Command::Import::WorkFlow::CreateInstrumentDataAndCopyBam;

use strict;
use warnings;

use Genome;

require Cwd;
require File::Basename;

class Genome::InstrumentData::Command::Import::WorkFlow::CreateInstrumentDataAndCopyBam { 
    is => 'Command::V2',
    has_input => [
        bam_paths => {
            is => 'Text',
            is_many => 1,
            doc => 'The paths of the bams to verify and move.',
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
        source_md5_paths => {
            is => 'Text',
            is_many => 1,
            doc => 'The MD5 paths of the retrieved source paths.',
        },
    ],
    has_output => [
        instrument_data => {
            is => 'Genome::InstrumentData',
            is_many => 1,
            doc => 'The instrument data entities created.',
        },
    ],
    has_transient_optional => [
        library => { is => 'Genome::Library', },
        kilobytes_requested => { is => 'Number', },
    ],
    has_calculated => [
        helpers => { calculate => q( Genome::InstrumentData::Command::Import::WorkFlow::Helpers->get; ), },
    ],
};

sub execute {
    my $self = shift;
    $self->status_message('Create instrument data and copy bam...');

    my $library = $self->_resvolve_library;
    return if not $library;

    my $helpers = $self->helpers;
    my @instrument_data;
    my @bam_paths = $self->bam_paths;
    for my $bam_path ( @bam_paths ) {
        # Create inst data
        my $instrument_data = $self->_create_instrument_data_for_bam_path($bam_path);
    $DB::single=1;
        return if not $instrument_data;

        # Move bam
        my $final_bam_path = $instrument_data->data_directory.'/all_sequences.bam';
        my $move_ok = $helpers->move_path($bam_path, $final_bam_path);
        return if not $move_ok;

        # Flagstat
        my $flagstat_path = $final_bam_path.'.flagstat';
        $self->status_message("Flagstat path: $flagstat_path");
        my $flagstat = $helpers->validate_bam($final_bam_path, $flagstat_path);
        return if not $flagstat;

        # Attrs
        $instrument_data->add_attribute(attribute_label => 'bam_path', attribute_value => $final_bam_path);
        $instrument_data->add_attribute(attribute_label => 'read_count', attribute_value => $flagstat->{total_reads});
        $instrument_data->add_attribute(attribute_label => 'is_paired_end', attribute_value => $flagstat->{is_paired_end});

        $self->status_message('Reallocate...');
        $instrument_data->allocations->reallocate;# with move??
        $self->status_message('Reallocate...done');

        push @instrument_data, $instrument_data;
    }
    $self->instrument_data(\@instrument_data);

    $self->status_message('Create instrument data and copy bam...done');
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

sub _create_instrument_data_for_bam_path {
    my ($self, $bam_path) = @_;
    $self->status_message('Create instrument data for bam path...');

    Carp::confess('No bam path to create instrument data!') if not $bam_path;
    $self->status_message('Bam path: '.$bam_path);

    my $additional_properties = $self->helpers->key_value_pairs_to_hash( $self->instrument_data_properties );
    return if not $additional_properties;

    my $original_data_path = delete $additional_properties->{original_data_path};

    my $read_group_ids_from_bam = $self->helpers->load_read_groups_from_bam($bam_path);
    return if not $read_group_ids_from_bam; # should only be one or 0
    if ( @$read_group_ids_from_bam > 1 ) {
        $self->error_message('More than one read group in bam! '.$bam_path);
        return;
    }

    if ( @$read_group_ids_from_bam ) {
        $self->status_message("Read groups in bam: @$read_group_ids_from_bam");
        if ( @$read_group_ids_from_bam > 1 ) {
            $self->error_message('Multiple read groups in bam! '.$bam_path);
            return;
        }
        $additional_properties->{segment_id} = $read_group_ids_from_bam->[0];
    }

    my @source_md5s = $self->_load_source_md5s;
    return if not @source_md5s;

    $self->status_message('Checking if source files were previously imported...');
    my %properties = (
        library => $self->library,
        original_data_path => $original_data_path,
    );
    my @existing_instrument_data = Genome::InstrumentData::Imported->get(%properties);
    if ( defined $read_group_ids_from_bam->[0] ) {
        @existing_instrument_data = grep { 
            my $attr = $_->attributes(attribute_label => 'segment_id', attribute_value => $read_group_ids_from_bam->[0]); 
            $attr ? $_ : undef
        } @existing_instrument_data;
    }

    if ( @existing_instrument_data ) {
        $self->error_message('Found existing instrument data for library and source files. Were these previously imported? Exiting instrument data id: '.join(', ', map { $_->id } @existing_instrument_data).', source files: '.$properties{original_data_path});
        return;
    }
    $self->status_message('Source files were NOT previously imported!');

    $properties{import_format} = 'bam';
    $properties{sequencing_platform} = 'solexa';

    my $instrument_data = Genome::InstrumentData::Imported->create(%properties);
    if ( not $instrument_data ) {
        $self->error_message('Failed to create instrument data!');
        return;
    }
    $self->status_message('Instrument data id: '.$instrument_data->id);

    for my $name ( keys %$additional_properties ) {
        $self->status_message('Add attribute: '.$name.' => '.$additional_properties->{$name});
        $instrument_data->add_attribute(attribute_label => $name, attribute_value => $additional_properties->{$name});
    }

    for my $md5 ( @source_md5s ) {
        $self->status_message('Add attribute: original_data_path_md5 => '.$md5);
        $instrument_data->add_attribute(attribute_label => 'original_data_path_md5', attribute_value => $md5);
    }

    my $bam_size = -s $bam_path;
    my $kilobytes_requested = int($bam_size / 1024) + 1024;
    my $allocation = Genome::Disk::Allocation->create(
        disk_group_name => 'info_alignments',
        allocation_path => 'instrument_data/imported/'.$instrument_data->id,
        kilobytes_requested => $kilobytes_requested,
        owner_class_name => $instrument_data->class,
        owner_id => $instrument_data->id,
    );
    if ( not $allocation ) {
        $self->error_message('Failed to create allocation for instrument data! '.$instrument_data->id);
        return;
    }
    $self->status_message('Allocation id: '.$allocation->id);
    $self->status_message('Allocation path: '.$allocation->absolute_path);

    $self->status_message('Create instrument data for bam path...done');
    return $instrument_data;
}

sub _load_source_md5s {
    my $self = shift;
    $self->status_message('Load source MD5...');

    my @md5s;
    for my $source_md5_path ( $self->source_md5_paths ) {
        $self->status_message('Source MD5 path: '.$source_md5_path);
        if ( not -s $source_md5_path ) {
            $self->error_message('Source MD5 path does not exist!');
            return;
        }
        my $md5 = $self->helpers->load_md5($source_md5_path);
        return if not $md5;
        $self->status_message('MD5: '.$md5);
        push @md5s, $md5;
    }

    $self->status_message('Load source MD5...done');
    return @md5s;
}

1;

