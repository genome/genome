package Genome::Test::Factory::ProcessingProfile::ClinSeq;
use Genome::Test::Factory::ProcessingProfile;
use Genome::Model::Tools::Sam::Readcount;
@ISA = (Genome::Test::Factory::ProcessingProfile);

use strict;
use warnings;

our @required_params = qw(bam_readcount_version);

sub create_bam_readcount_version {
    return Genome::Model::Tools::Sam::Readcount->default_version;
}

1;
