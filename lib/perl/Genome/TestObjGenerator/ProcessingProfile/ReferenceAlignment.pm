package Genome::TestObjGenerator::ProcessingProfile::ReferenceAlignment;
use Genome::TestObjGenerator::ProcessingProfile;
@ISA = (Genome::TestObjGenerator::ProcessingProfile);

use strict;
use warnings;

my @required_params = ("sequencing_platform", "dna_type", "read_aligner_name", "snv_detection_strategy");

sub create_sequencing_platform {
    return "solexa";
}

sub create_dna_type {
    return "cdna";
}

sub create_read_aligner_name {
    return "bwa";
}

sub create_snv_detection_strategy {
    return "samtools";
}

sub get_required_params {
    my $class = shift;
    my $super = $class->SUPER::get_required_params;
    my @all =  (@$super, @required_params);
    return \@all;
}

1;
