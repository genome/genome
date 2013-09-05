#!/usr/bin/env genome-perl

use strict;
use warnings;
use above 'Genome';
use Test::More tests => 5;

use_ok('Genome::Model::Tools::Vcf::FilterVcfInfo');

# To rebuild results, run like:
# perl FilterVcfInfo.t REBUILD v2
# then update the version on the next line
my $version = "v1";

my $rebuild = 0;
if(defined($ARGV[0]) && ($ARGV[0] eq "REBUILD")){
    $rebuild = 1;
}

my $test_base = $ENV{GENOME_TEST_INPUTS} . "/Genome-Model-Tools-Vcf-FilterVcfInfo";
my $test_dir = $test_base . "/" . $version;

if($rebuild){
    mkdir "$test_base/$ARGV[1]";
    `cp $test_dir/small.vcf.gz $test_base/$ARGV[1]`;
    $test_dir = "$test_base/$ARGV[1]";
}

print STDERR "using data in $test_dir \n";

my $output_file = Genome::Sys->create_temp_file_path;
my $input_dir = "$test_dir";
my $input_vcf = "$input_dir/small.vcf.gz";
my $expected_result = "$input_dir/small.filtered.vcf.gz";

my $command= Genome::Model::Tools::Vcf::FilterVcfInfo->create(
    output_file => $output_file,
    vcf_file => $input_vcf,
    filters => "GMAFFILTER:GMAF>0.05,TIER2FILTER:TIER=2",
    filter_descriptions => "has GMAF>0.05,tier is 2",
);

ok($command, 'Command created');
my $rv = $command->execute;
ok($rv, 'Command completed successfully');
ok(-s $output_file, "output file created");

if($rebuild){
    `cp $output_file $test_dir/small.filtered.vcf.gz`;
}

my $diff = Genome::Sys->diff_file_vs_file($output_file, $expected_result);
ok(!$diff, 'output matched expected result');
