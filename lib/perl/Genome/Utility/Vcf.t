#!/usr/bin/env genome-perl

use strict;
use warnings;
use File::Slurp "read_file";
use Test::More;

use above 'Genome';
use Genome::Utility::Vcf qw(get_vcf_header convert_file_with_alpha_gt_values_to_numeric diff_vcf_file_vs_file);
                    

# Test header reading
{
    my $input_filename = join( '/',  __FILE__ . ".d",  'input.clean.vcf');
    my $expected_filename = join( '/',  __FILE__ . ".d",  'expected.txt');

    my $header = get_vcf_header($input_filename);
    my $expected_header = read_file($expected_filename);

    is($header, $expected_header, "Read in the header of a vcf file correctly.");
}

# Test alpha GT value -> numeric alt value fixing from polymutt denovo files
{
    my $input_filename = join( '/',  __FILE__ . ".d",  'input_denovo.vcf.gz');
    my $expected_filename = join( '/',  __FILE__ . ".d",  'expected_denovo.vcf.gz');

    my $output_filename = Genome::Sys->create_temp_file_path;
    convert_file_with_alpha_gt_values_to_numeric($input_filename, $output_filename);

    my $diff = diff_vcf_file_vs_file($expected_filename, $output_filename);
    is($diff, "", "convert_file_with_alpha_gt_values_to_numeric outputs as expected");
}

done_testing();
1;
