package Genome::Test::Factory::ProcessingProfile::ClinSeq;
use Genome::Test::Factory::ProcessingProfile;
use Genome::Model::Tools::Sam::Readcount;
@ISA = (Genome::Test::Factory::ProcessingProfile);

use strict;
use warnings;

our @required_params = qw(bam_readcount_version
    sireport_min_coverage
    sireport_min_tumor_vaf
    sireport_max_normal_vaf);

sub create_bam_readcount_version {
    return Genome::Model::Tools::Sam::Readcount->default_version;
}

sub create_sireport_min_coverage {
    return 20;
}

sub create_sireport_min_tumor_vaf {
    return 2.5;
}

sub create_sireport_max_normal_vaf {
    return 10;
}

sub create_sireport_min_mq{
    return "30,40";
}

sub create_sireport_min_bq {
    return "10,20";
}

1;
