package Genome::Test::Factory::ProcessingProfile;
use Genome::Test::Factory::Base;
@ISA = (Genome::Test::Factory::Base);

use strict;
use warnings;
use Genome;
use Genome::ProcessingProfile;
use Genome::Test::Factory::Util;
use Sub::Override;

our @required_params = qw(name);

sub generate_obj {
    my $class = shift;

    (my $pp_class = $class) =~ s/::Test::Factory::/::/;

    my $pp;
    {
        # Override this method so we don't crash because of similarities with existing production profiles
        my $override = Sub::Override->new('Genome::ProcessingProfile::_validate_no_existing_processing_profiles_with_identical_params', sub {return 1});
        $pp = $pp_class->create(@_);
    }
    return $pp;
}

sub create_name {
    return Genome::Test::Factory::Util::generate_name("test_processing_profile");
}

1;
