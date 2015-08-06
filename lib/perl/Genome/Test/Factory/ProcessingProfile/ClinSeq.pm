package Genome::Test::Factory::ProcessingProfile::ClinSeq;
use Genome::Test::Factory::ProcessingProfile;
use Genome::Model::Tools::Sam::Readcount;
@ISA = (Genome::Test::Factory::ProcessingProfile);

use strict;
use warnings;

our @required_params = qw(bam_readcount_version
    exome_cnv
    sireport_min_coverage
    sireport_min_tumor_vaf
    sireport_max_normal_vaf
    sireport_min_mq_bq
);

#Version of bam-readcount version to use
sub create_bam_readcount_version {
    return Genome::Model::Tools::Sam::Readcount->default_version;
}

#Run exome-cnv step or not?
sub create_exome_cnv {
    return 1;
}

#Minimum coverage parameter for SNV-Indel report
sub create_sireport_min_coverage {
    return 20;
}

#Minimum tumor VAF for SNV-Indel report
sub create_sireport_min_tumor_vaf {
    return 2.5;
}

#Mamximum normal VAF for SNV-Indel report
sub create_sireport_max_normal_vaf {
    return 10;
}

#Minimum base and mapping qualities for SNV Indel report
sub create_sireport_min_mq_bq {
    return "30,40;10,20";
}


1;
