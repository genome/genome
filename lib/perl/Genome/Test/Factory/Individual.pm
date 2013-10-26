package Genome::Test::Factory::Individual;
use base qw(Genome::Test::Factory::Base);

use strict;
use warnings;

use Genome;
use Genome::Test::Factory::Taxon;

our @required_params = qw(name common_name taxon_id);

sub generate_obj {
    my $self = shift;
    return Genome::Individual->create(@_);
}

sub create_name {
    return Genome::Test::Factory::Util::generate_name('individual_name');
}

sub create_common_name {
    return Genome::Test::Factory::Util::generate_name('TEST');
}

sub create_taxon_id {
    my $taxon = Genome::Test::Factory::Taxon->setup_object();
    return $taxon->id;
}

1;
