package Genome::Test::Factory::Config::Profile::Item;
use base qw(Genome::Test::Factory::Base);

use strict;
use warnings;

use Genome;
use Genome::Test::Factory::DiskAllocation;

our @required_params = qw(analysis_project);

sub generate_obj {
    my $self = shift;
    return Genome::Config::Profile::Item->create(@_);
}

sub create_analysis_project {
    return Genome::Test::Factory::AnalysisProject->setup_object(name => Genome::Test::Factory::Util::generate_name('anp_name'));
}
