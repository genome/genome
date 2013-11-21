package Genome::InstrumentData::Microarray;

use strict;
use warnings;

use Genome;

require File::Path;

class Genome::InstrumentData::Microarray {
};

sub get_snpid_hash_for_variant_list {
    my $self = shift;
    my $instrument_data = shift || die 'No instrument data given';
    my $variation_list_build = shift || die 'No variation list build';
    $self->status_message('Get snp id hash for variant list...');

    $self->status_message("Instrument data: ".$instrument_data->id);
    my $platform = lc($instrument_data->sequencing_platform);
    $self->status_message("Platform: $platform");
    my $chip_attribute = $instrument_data->attributes(attribute_label => 'chip_name');
    my $chip = ( $chip_attribute ? lc($chip_attribute->attribute_value) : '' );
    $self->status_message("Chip: $chip");
    my $version_attribute = $instrument_data->attributes(attribute_label => 'version');
    my $version = ( $version_attribute ? lc($version_attribute->attribute_value) : '' );
    $self->status_message("Version: $version");

    $self->status_message('Looking for allocation for snp mapping file...');
    my $allocation = Genome::Disk::Allocation->get(allocation_path => "microarray_data/$platform-$chip-$version");
    if ( not $allocation ) {
        $self->status_message('No allocation found! This may be exepected.');
        return;
    }
    $self->status_message('Found allocation: '.$allocation->id);

    my $snp_map_path = $allocation->absolute_path . '/mapping.tsv';
    $self->status_message('Snp mapping file: '.$snp_map_path);
    Carp::confess('No snp id mapping file in '.$allocation->absolute_path) if not -s $snp_map_path;

    my $fh = eval{ Genome::Sys->open_file_for_reading($snp_map_path); };
    Carp::confess('Failed to open snp id mapping file! '.$snp_map_path) if not $fh;

    my %old_to_new_snpid_map;
    while ( my $mapping = $fh->getline ) {
        chomp $mapping;
        my ($old_id, $new_id) = split /\t/, $mapping;
        Carp::confess("Invalid line in snp mapping file! '$mapping'") if not defined $old_id and not defined $new_id;
        $old_to_new_snpid_map{$old_id} = $new_id;
    }

    $self->status_message('Get snp id hash for variant list...OK');
    return \%old_to_new_snpid_map;
}

#   Copy a genotype_file to a new disk allocation and update the instrument_data's attribute that
# points to that file.
sub update_genotype_file {
    my ($self, $instrument_data, $genotype_file) = @_;
    die $self->error_message("No instrument data was given to update_genotype_file.") if not $instrument_data;
    die $self->error_message("No genotype file was given to update_genotype_file.") if not $genotype_file;

    my $new_genotype_file = $self->_copy_genotype_file($instrument_data, $genotype_file);
    $self->_update_instrument_data_genotype_file_attribute($instrument_data, $new_genotype_file);
    return $new_genotype_file;
}

sub _update_instrument_data_genotype_file_attribute {
    my($self, $instrument_data, $new_genotype_file) = @_;

    # Remove any genotype_file attributes that already exist on this instrument_data.
    my @genotype_file_attributes = $instrument_data->attributes(attribute_label => 'genotype_file');
    for my $to_be_removed (@genotype_file_attributes) {
        $to_be_removed->delete(); 
    }

    # Create new genotype_file attribute and attach it.
    my $new_attribute = Genome::InstrumentDataAttribute->create(
        instrument_data_id => $instrument_data->id,
        attribute_label => 'genotype_file',
        attribute_value => $new_genotype_file, 
    );
    if (not $new_attribute){
        die $self->error_message("Failed to create new genotype_file attribute for file: ".$new_genotype_file);
    }
    return 1;
}

# Copy the genotype_file to a new disk allocation (if one doesn't already exist)
sub _copy_genotype_file {
    my ($self, $instrument_data, $genotype_file) = @_;

    my $size = -s $genotype_file;
    if (not $size) {
        die $self->error_message("The genotype file has no size: $genotype_file");
    }

    my ($disk_allocation) = $instrument_data->allocations;
    my $kilobytes_requested = int( $size / 1024 ) + 20; # extra for directory etc
    if (not $disk_allocation) {
        $disk_allocation = Genome::Disk::Allocation->allocate (
            disk_group_name => $ENV{GENOME_DISK_GROUP_ALIGNMENTS},
            allocation_path => 'instrument_data/imported/'.$instrument_data->id,
            kilobytes_requested => $kilobytes_requested,
            owner_class_name => $instrument_data->class,
            owner_id => $instrument_data->id,
        );
        if (not $disk_allocation) {
            die $self->error_message("Failed to create disk allocation for instrument data: ".$instrument_data->id);
        }
    }
    else {
        $disk_allocation->kilobytes_requested($kilobytes_requested);
        File::Path::mkpath($disk_allocation->absolute_path) if not -d $disk_allocation->absolute_path;
    }

    # check that the dest dir exists
    if ( not -d $disk_allocation->absolute_path ) {
        Carp::confess('Microarray ('.$self->id.') disk allocation ('.$disk_allocation->id.') absolute path does not exist! '.$disk_allocation->absolute_path);
    }

    # remove existing genotype_file (if exists)
    my $new_genotype_file = $disk_allocation->absolute_path.'/'.$instrument_data->sample->id.'.genotype';
    unlink $new_genotype_file if -e $new_genotype_file;

    # copy *this* genotype_file into the allocation
    my $copy = File::Copy::copy($genotype_file, $new_genotype_file);
    if (not $copy) {
        die $self->error_message("Failed to copy genotype file from $genotype_file to $new_genotype_file => $!");
    }

    # resize the allocation based on what is stored there.
    my $realloc = eval{ $disk_allocation->reallocate(allow_reallocate_with_move => 1); };
    if (not $realloc) {
        die $self->error_message('Failed to reallocate instrument data ('.$disk_allocation->owner_id.')');
    }

    return $new_genotype_file;
}
