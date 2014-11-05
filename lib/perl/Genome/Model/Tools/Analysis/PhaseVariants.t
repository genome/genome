#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Test::More tests => 1;

my $plot = Genome::Model::Tools::Analysis::PhaseVariants->create(
    distance  => 700,
    command => "C",
    vcf_file => "/gscuser/cfederer/Documents/mysnvs.vcf",
    bam_file => "/gscuser/cfederer/Documents/p53.bam",
    chromosome => 17,
    sample => "none",
    relax => "false",
    snvs => "7579646T7579653C",
    print_reads => "t",
    out_file => "/gscuser/cfederer/git/genome/lib/perl/Genome/Model/phaseVariantsOutPut.txt"
);
ok($plot->execute,'Executed');

