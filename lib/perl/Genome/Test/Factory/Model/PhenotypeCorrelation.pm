package Genome::Test::Factory::Model::PhenotypeCorrelation;

use strict;
use warnings;

use Genome;
use base qw(Genome::Test::Factory::Model);

use Genome::Test::Factory::Sample;
use Genome::Test::Factory::ProcessingProfile::PhenotypeCorrelation;

our @required_params = qw(subject_id processing_profile_id);

sub create_processing_profile_id {
    my $p = Genome::Test::Factory::ProcessingProfile::PhenotypeCorrelation->setup_object();
    return $p->id;
}

sub create_subject_id {
    my $subject = Genome::Test::Factory::Sample->setup_object();
    return $subject->id;
}

sub generate_obj {
    my $self = shift;
    my $m = Genome::Model::PhenotypeCorrelation->create(@_);
    return $m;
}

1;
