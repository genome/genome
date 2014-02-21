#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More;

###
# Writing is tested in ReadWrite.t
###

my $class = 'Genome::Model::GenotypeMicroarray::GenotypeFile::WriteVcf';
use_ok($class) or die;

my $genotype;
ok(!eval{$class->format_string_for_genotype($genotype); }, 'failed to get format string for undef genotype');
$genotype = { reference => 'A', allele1 => 'A', allele2 => 'A', };
is($class->format_string_for_genotype($genotype), "GT\t0/0", 'homozygous ref format string');
$genotype = { reference => 'A', allele1 => 'A', allele2 => 'T', };
is($class->format_string_for_genotype($genotype), "GT\t0/1", 'heterozygous ref format string');
$genotype = { reference => 'A', allele1 => 'T', allele2 => 'T', };
is($class->format_string_for_genotype($genotype), "GT\t1/1", 'homozygous alt format string');
$genotype = { reference => 'A', allele1 => 'G', allele2 => 'T', };
is($class->format_string_for_genotype($genotype), "GT\t1/2", 'nothing matches format string');
$genotype = { reference => 'A', allele1 => '-', allele2 => '-', };
is($class->format_string_for_genotype($genotype), "GT\t./.", 'alleles are dashes');

done_testing();
