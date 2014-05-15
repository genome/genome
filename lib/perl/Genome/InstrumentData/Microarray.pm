package Genome::InstrumentData::Microarray;

use strict;
use warnings;

use Genome;

require File::Basename;
require File::Path;
require Sub::Install;

class Genome::InstrumentData::Microarray {
};

BEGIN: {
    use Genome::InstrumentData::Imported;
    #Genome::InstrumentData::Imported->class;
    for my $method (qw/ 
            genotype_file update_genotype_file
            _update_instrument_data_genotype_file_name_attribute
            _copy_genotype_file
        /) {
        Sub::Install::install_sub({
           code => $method,
           into => 'Genome::InstrumentData::Imported',
           as   => $method
         });
    }
}

sub get_snpid_hash_for_variant_list {
    my $self = shift;
    my $instrument_data = shift || die 'No instrument data given';
    my $variation_list_build = shift || die 'No variation list build';
    $self->debug_message('Get snp id hash for variant list...');

    $self->debug_message("Instrument data: ".$instrument_data->id);
    my $platform = lc($instrument_data->sequencing_platform);
    $self->debug_message("Platform: $platform");
    my $chip_attribute = $instrument_data->attributes(attribute_label => 'chip_name');
    my $chip = ( $chip_attribute ? lc($chip_attribute->attribute_value) : '' );
    $self->debug_message("Chip: $chip");
    my $version_attribute = $instrument_data->attributes(attribute_label => 'version');
    my $version = ( $version_attribute ? lc($version_attribute->attribute_value) : '' );
    $self->debug_message("Version: $version");

    $self->debug_message('Looking for allocation for snp mapping file...');
    my $allocation = Genome::Disk::Allocation->get(allocation_path => "microarray_data/$platform-$chip-$version");
    if ( not $allocation ) {
        $self->status_message('No allocation found! This may be exepected.');
        return;
    }
    $self->debug_message('Found allocation: '.$allocation->id);

    my $snp_map_path = $allocation->absolute_path . '/mapping.tsv';
    $self->debug_message('Snp mapping file: '.$snp_map_path);
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

    $self->debug_message('Get snp id hash for variant list...OK');
    return \%old_to_new_snpid_map;
}

#< Genotype File >#
#
# These methods are installed into Genome::Instrument::Imported. The first parameter is $gm_inst_data,
#  and is the microarry inst data.
#
sub genotype_file {
    my ($gm_inst_data, $genotype_file) = @_;

    # Error if given a file to update
    if ( $genotype_file ) {
        $gm_inst_data->error_message("To update the genotype file, please use the 'update_genotype_file' method.");
        return;
    }

    # Check allocation exists and is not archived
    my ($allocation) = $gm_inst_data->allocations;
    if ( not $allocation ) {
        return;
    }
    if ( $allocation->is_archived ) {
        $gm_inst_data->error_message("Genotype file is archived! Returning path to archived location.");
    }

    # Get genotype file name via 
    #  genotype_file_name attr 
    #   or basename of genotype_file attr
    my $genotype_file_name;
    if ( my $genotype_file_name_attr = $gm_inst_data->attributes(attribute_label => 'genotype_file_name') ) {
        $genotype_file_name = $genotype_file_name_attr->attribute_value;
    }
    elsif ( my $genotype_file_attr = $gm_inst_data->attributes(attribute_label => 'genotype_file') ) {
        $genotype_file_name = File::Basename::basename($genotype_file_attr->attribute_value);
    }
    else {
        return;
    }
    return if not $genotype_file_name;

    # return abs path + genotype file name
    return join('/', $allocation->absolute_path, $genotype_file_name);
}

sub update_genotype_file {
    my ($gm_inst_data, $genotype_file) = @_;
    die $gm_inst_data->error_message("No genotype file was given to update_genotype_file.") if not $genotype_file;

    my $new_genotype_file = $gm_inst_data->_copy_genotype_file($genotype_file);
    $gm_inst_data->_update_instrument_data_genotype_file_name_attribute($new_genotype_file);
    return $new_genotype_file;
}

sub _update_instrument_data_genotype_file_name_attribute {
    my($gm_inst_data, $new_genotype_file) = @_;

    # Remove any genotype_file and genotype_file_name attributes that already exist
    my @old_attributes = $gm_inst_data->attributes('attribute_label in' => [ 'genotype_file', 'genotype_file_name']);
    for ( @old_attributes ) { $_->delete(); }

    # Create new genotype_file_name attribute
    my $new_attribute = Genome::InstrumentDataAttribute->create(
        instrument_data_id => $gm_inst_data->id,
        attribute_label => 'genotype_file_name',
        attribute_value => File::Basename::basename($new_genotype_file), 
    );
    if (not $new_attribute){
        die $gm_inst_data->error_message("Failed to create new genotype_file attribute for file: ".$new_genotype_file);
    }

    return 1;
}

sub _copy_genotype_file {
    my ($gm_inst_data, $genotype_file) = @_;

    my $size = -s $genotype_file;
    if (not $size) {
        die $gm_inst_data->error_message("The genotype file has no size: $genotype_file");
    }

    my ($disk_allocation) = $gm_inst_data->allocations;
    my $kilobytes_requested = int( $size / 1024 ) + 20; # extra for directory etc
    if (not $disk_allocation) {
        $disk_allocation = Genome::Disk::Allocation->allocate (
            disk_group_name => $ENV{GENOME_DISK_GROUP_ALIGNMENTS},
            allocation_path => 'instrument_data/imported/'.$gm_inst_data->id,
            kilobytes_requested => $kilobytes_requested,
            owner_class_name => $gm_inst_data->class,
            owner_id => $gm_inst_data->id,
        );
        if (not $disk_allocation) {
            die $gm_inst_data->error_message("Failed to create disk allocation for instrument data: ".$gm_inst_data->id);
        }
    }
    else {
        $disk_allocation->kilobytes_requested($kilobytes_requested);
        File::Path::mkpath($disk_allocation->absolute_path) if not -d $disk_allocation->absolute_path;
    }

    # check that the dest dir exists
    if ( not -d $disk_allocation->absolute_path ) {
        Carp::confess('Microarray ('.$gm_inst_data->id.') disk allocation ('.$disk_allocation->id.') absolute path does not exist! '.$disk_allocation->absolute_path);
    }

    # remove existing genotype_file (if exists)
    my $new_genotype_file = $disk_allocation->absolute_path.'/'.$gm_inst_data->sample->id.'.genotype';
    unlink $new_genotype_file if -e $new_genotype_file;

    # copy *this* genotype_file into the allocation
    my $copy = File::Copy::copy($genotype_file, $new_genotype_file);
    if (not $copy) {
        die $gm_inst_data->error_message("Failed to copy genotype file from $genotype_file to $new_genotype_file => $!");
    }

    # resize the allocation based on what is stored there.
    my $realloc = eval{ $disk_allocation->reallocate(allow_reallocate_with_move => 1); };
    if (not $realloc) {
        $gm_inst_data->status_message('Failed to reallocate instrument data ('.$disk_allocation->owner_id.'), continuing...');
    }

    return $new_genotype_file;
}
#<>#

1;

