package Genome::TestObjGenerator::Individual;
use base qw(Genome::TestObjGenerator::Base);

use strict;
use warnings;
use Genome;

our @required_params = qw(name common_name);

sub generate_obj {
    my $self = shift;
    return Genome::Individual->create(@_);
}

sub create_name {
    return Genome::TestObjGenerator::Util::generate_name('individual_name');
}

sub create_common_name {
    return Genome::TestObjGenerator::Util::generate_name('TEST');
}

1;
