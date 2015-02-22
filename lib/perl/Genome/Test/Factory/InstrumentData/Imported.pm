package Genome::Test::Factory::InstrumentData::Imported;
use base qw(Genome::Test::Factory::Base);

use strict;
use warnings;

use Genome;

require File::Basename;
use Genome::Test::Factory::DiskAllocation;
use Genome::Test::Factory::Library;

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
    return $instrument_data if not $mount_path;

    my $alloc = Genome::Test::Factory::DiskAllocation->generate_obj(
        owner => $instrument_data,
        mount_path => $mount_path,
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
