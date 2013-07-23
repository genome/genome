package Genome::TestObjGenerator::Build;
use Genome::TestObjGenerator::Base;
@ISA = (Genome::TestObjGenerator::Base);

use strict;
use warnings;
use Genome;

my @required_params = ("model_id", "data_directory");

sub generate_obj {
    my $self = shift;
    return Genome::Model::Build->create(@_);
}

sub get_required_params {
    return \@required_params;
}

sub create_data_directory {
    return Genome::Sys->create_temp_directory();
}

1;

