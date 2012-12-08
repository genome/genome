#!/usr/bin/env genome-perl

use strict;
use warnings;
use File::Slurp "read_file";
use Test::More;

use above 'Genome';
use Genome::Utility::Vcf qw(get_vcf_header);

my $input_filename = join( '/',  __FILE__ . ".d",  'input.clean.vcf');
my $expected_filename = join( '/',  __FILE__ . ".d",  'expected.txt');

my $header = get_vcf_header($input_filename);
my $expected_header = read_file($expected_filename);

is($header, $expected_header, "Read in the header of a vcf file correctly.");

done_testing();
1;
