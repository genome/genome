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
        library => {
            is => 'Genome::Library',
            doc => 'Library to use. The library must exist.',
        },
        analysis_project => {
            is => 'Genome::Config::AnalysisProject',
            doc => 'Analysis project to assign to the created instrument data.',
        },
        instrument_data_properties => {
            is => 'Hash',
            doc => 'Hash of instrument data properties to store.',
        },
        source_md5s => {
            is => 'Text',
            is_many => 1,
            doc => 'The MD5s of the retrieved source paths.',
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
    instrument_data_subset_name => { is => 'Text', },
        kilobytes_requested => { is => 'Number', },
    ],
    has_calculated => [
        helpers => { calculate => q( Genome::InstrumentData::Command::Import::WorkFlow::Helpers->get; ), },
    ],
};

sub __errors__ {
    my $self = shift;

    my @errors = $self->SUPER::__errors__;
    return @errors if @errors;

    my $instrument_data_properties = $self->instrument_data_properties;
    if ( not $instrument_data_properties->{original_data_path} ) {
        push @errors, UR::Object::Tag->create(
            type => 'invalid',
            properties => [qw/ instrument_data_properties /],
            desc => 'No original data path!',
        );
    }

    my $subset_name = delete $instrument_data_properties->{subset_name};
    $subset_name = 'unknown' if not $subset_name;
    $self->instrument_data_subset_name($subset_name);

    $instrument_data_properties->{import_format} = 'bam';
    delete $instrument_data_properties->{is_paired_end};

    $self->instrument_data_properties($instrument_data_properties);

    return @errors;
}

sub execute {
    my $self = shift;
    $self->debug_message('Create instrument data and copy bam...');

    my @source_md5s = $self->source_md5s;
    my $previously_imported = $self->helpers->were_original_path_md5s_previously_imported(md5s => \@source_md5s);
    return if $previously_imported;

    my $library = defined($self->library)
        ? $self->library
        : $self->_resvolve_library;
    return if not $library;

    my $helpers = $self->helpers;
    my @instrument_data;
    my @bam_paths = $self->bam_paths;
    for my $bam_path ( @bam_paths ) {
        # Create inst data
        my $instrument_data = $self->_create_instrument_data_for_bam_path($bam_path);
        return if not $instrument_data;

        # Create allocation
        my $allocation = $self->_create_allocation($instrument_data, $bam_path);
        return if not $allocation;

        # Move bam and flagstat
        my $final_bam_path = $instrument_data->data_directory.'/all_sequences.bam';
        my $move_ok = $helpers->move_path($bam_path, $final_bam_path);
        return if not $move_ok;

        my $flagstat_path = $bam_path.'.flagstat';
        my $final_flagstat_path = $final_bam_path.'.flagstat';
        $move_ok = $helpers->move_path($flagstat_path, $final_flagstat_path);
        return if not $move_ok;

        # Attrs
        $self->helpers->update_bam_metrics_for_instrument_data($instrument_data);

        # Analysis Project
        $instrument_data->add_analysis_project_bridge(analysis_project => $self->analysis_project);

        # Reallocate
        $self->debug_message('Reallocate...');
        $instrument_data->reallocate_disk_allocations;
        $self->debug_message('Reallocate...done');

        push @instrument_data, $instrument_data;
    }
    $self->instrument_data(\@instrument_data);

    $self->debug_message('Imported instrument data bam paths:');
    for my $instrument_data ( @instrument_data ) {
        printf("%s\n", $instrument_data->bam_path);
    }

    $self->debug_message('Create instrument data and copy bam...done');
    return 1;
}

sub _create_instrument_data_for_bam_path {
    my ($self, $bam_path) = @_;
    $self->debug_message('Create instrument data for bam path...');

    Carp::confess('No bam path to create instrument data!') if not $bam_path;
    $self->debug_message('Bam path: '.$bam_path);

    my $read_group_ids_from_bam = $self->helpers->load_read_groups_from_bam($bam_path);
    return if not $read_group_ids_from_bam;
    if ( not @$read_group_ids_from_bam ) {
        die $self->error_message('No read groups in bam!');
    }
    $self->debug_message("Read groups in bam: @$read_group_ids_from_bam");
    if ( @$read_group_ids_from_bam > 1 ) {
        $self->error_message('More than one read group in bam! '.$bam_path);
        return;
    }

    my $properties = $self->instrument_data_properties;
    $properties->{segment_id} = $read_group_ids_from_bam->[0];

    my $instrument_data = Genome::InstrumentData::Imported->create(
        id => $read_group_ids_from_bam->[0],
        library => $self->library,
        subset_name => $self->instrument_data_subset_name,
        sequencing_platform => 'solexa',
    );
    if ( not $instrument_data ) {
        $self->error_message('Failed to create instrument data!');
        return;
    }
    $self->status_message('Instrument data id: '.$instrument_data->id);

    for my $name ( keys %$properties ) {
        $self->debug_message('Add attribute: '.$name.' => '.$properties->{$name});
        $instrument_data->add_attribute(attribute_label => $name, attribute_value => $properties->{$name});
    }

    for my $md5 ( $self->source_md5s ) {
        $self->debug_message('Add attribute: original_data_path_md5 => '.$md5);
        $instrument_data->add_attribute(attribute_label => 'original_data_path_md5', attribute_value => $md5);
    }

    return $instrument_data;
}

sub _create_allocation {
    my ($self, $instrument_data, $bam_path) = @_;

    Carp::confess('No instrument data given to create allocation!') if not $instrument_data;

    my $bam_size = -s $bam_path;
    my $kilobytes_requested = int($bam_size / 1024) + 1024;
    my $allocation = Genome::Disk::Allocation->create(
        disk_group_name => Genome::Config::get('disk_group_alignments'),
        allocation_path => 'instrument_data/imported/'.$instrument_data->id,
        kilobytes_requested => $kilobytes_requested,
        owner_class_name => $instrument_data->class,
        owner_id => $instrument_data->id,
    );
    if ( not $allocation ) {
        $self->error_message('Failed to create allocation for instrument data! '.$instrument_data->id);
        return;
    }
    $self->debug_message('Allocation id: '.$allocation->id);
    $self->debug_message('Allocation path: '.$allocation->absolute_path);

    $self->debug_message('Create instrument data for bam path...done');
    return $instrument_data;
}

1;

