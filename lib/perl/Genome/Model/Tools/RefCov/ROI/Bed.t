#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More;

BEGIN {
    use_ok('Genome::Model::Tools::RefCov::ROI::RegionI');
    use_ok('Genome::Model::Tools::RefCov::ROI::Region');
    use_ok('Genome::Model::Tools::RefCov::ROI::FileI');
    use_ok('Genome::Model::Tools::RefCov::ROI::Bed');
}

# TODO: Subset bed file into one or two entries per chr
my $file = $ENV{GENOME_TEST_INPUTS} . '/Genome-RefCov-ROI-Bed/SANGER.bed';
if (-s $file) {
    my $region_set = Genome::Model::Tools::RefCov::ROI::Bed->create(file => $file);
    my @chromosomes = $region_set->chromosomes;
    is(scalar(@chromosomes), 25, 'got 25 chromosomes');
    while (my $region = $region_set->next_region) {
        isa_ok($region, 'HASH', 'got region as hash');
        last;
    }
}

done_testing();
