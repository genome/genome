package Genome::Test::Factory::Library;
use base qw(Genome::Test::Factory::Base);

use strict;
use warnings;

use Genome;
use Genome::Test::Factory::Sample;

our @required_params = qw(name sample_id);

sub generate_obj {
    my $self = shift;
    return Genome::Library->create(@_);
}

sub create_name {
    return Genome::Test::Factory::Util::generate_name('library_name');
}

sub create_sample_id {
    my $sample = Genome::Test::Factory::Sample->setup_object();
    return $sample->id;
}

1;
