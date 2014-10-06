package Genome::Test::Factory::ProcessingProfile::ClinSeq;
use Genome::Test::Factory::ProcessingProfile;
use Genome::Model::Tools::Sam::Readcount;
@ISA = (Genome::Test::Factory::ProcessingProfile);

use strict;
use warnings;

our @required_params = qw(bamrc_version bamrc_version_clonality);

sub create_bamrc_version {
    return Genome::Model::Tools::Sam::Readcount->default_version;
}

sub create_bamrc_version_clonality {
    return Genome::Model::Tools::Sam::Readcount->default_version;
}

1;
