package Genome::Test::Factory::ProcessingProfile::SingleSampleGenotype;
use Genome::Test::Factory::ProcessingProfile;
@ISA = (Genome::Test::Factory::ProcessingProfile);

use strict;
use warnings;

our @required_params = qw(aligner_version aligner_api_version qc_config haplotype_caller_version);

sub create_aligner_version { return '1' }
sub create_aligner_api_version { return 'v6' }
sub create_qc_config { return 'test' }
sub create_haplotype_caller_version { return '3.4' }

1;
