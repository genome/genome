package Genome::Test::Factory::ProcessingProfile;
use Genome::Test::Factory::Base;
@ISA = (Genome::Test::Factory::Base);

use strict;
use warnings;
use Genome;
use Genome::Test::Factory::Util;

our @required_params = qw(name);

sub generate_obj {
    my $class = shift;
    (my $pp_class = $class) =~ s/::Test::Factory::/::/;
    return $pp_class->create(@_);
}

sub create_name {
    return Genome::Test::Factory::Util::generate_name("test_processing_profile");
}

1;
