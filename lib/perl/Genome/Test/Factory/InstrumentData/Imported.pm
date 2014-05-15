package Genome::Test::Factory::InstrumentData::Imported;
use base qw(Genome::Test::Factory::Base);

use strict;
use warnings;

use Genome;
use Genome::Test::Factory::Library;

require File::Basename;

our @required_params = qw(library_id sequencing_platform);

sub generate_obj {
    my ($self, %params) = @_;

    my $genotype_file = delete $params{genotype_file};
    my ($mount_path, $genotype_file_name);
    if ( $genotype_file ) {
        ($genotype_file_name, $mount_path) = File::Basename::fileparse($genotype_file);
    }

    my $instrument_data = Genome::InstrumentData::Imported->create(%params);
    die 'Failed to create instrument data' if not $instrument_data;
    return $instrument_data if not $genotype_file;

    my $volume = Genome::Disk::Volume->get(mount_path => $mount_path);
    if ( not $volume ) {
        $volume = Genome::Disk::Volume->__define__(mount_path => $mount_path, disk_status => 'active');
    }

    my $alloc = Genome::Disk::Allocation->__define__(
        owner => $instrument_data,
        mount_path => $volume->mount_path,
        group_subdirectory => '',
        allocation_path => '',
    );
    die "Failed to define allocation for genotype file!" if not $alloc;

    $instrument_data->add_attribute(attribute_label => 'genotype_file_name', attribute_value => $genotype_file_name);

    return $instrument_data;
}

sub create_library_id {
    my $lib = Genome::Test::Factory::Library->setup_object();
    return $lib->id;
}

sub create_sequencing_platform {
    return 'plink';
}

1;
