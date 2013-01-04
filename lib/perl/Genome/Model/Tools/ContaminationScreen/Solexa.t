#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More tests => 3;
use File::Copy;

my $datadir = __FILE__ . '.d';

BEGIN {use_ok('Genome::Model::Tools::ContaminationScreen::Solexa');}

my $temp_dir = Genome::Sys->base_temp_directory();

# crossmatch leaves a bunch of files lingering around in the same directory as the 
# input reads file, so copy this out to temp to let it do its thing.

copy($datadir . "/test3.fna", $temp_dir . "/input_reads.fna");


my $summary_file = $temp_dir . "/solexa_contamination_screen_summary";
my $output_file = $temp_dir . "/solexa_contamination_screen_summary";

my %params;
$params{input_file} = $temp_dir . "/input_reads.fna";
$params{database} = $datadir . '/test2.fna';
$params{minscore} = 42;
$params{output_file} = $output_file;
$params{summary_file} = $summary_file;

my $solexa = Genome::Model::Tools::ContaminationScreen::Solexa->create(%params);

isa_ok($solexa, 'Genome::Model::Tools::ContaminationScreen::Solexa');

ok($solexa->execute, "solexa executing");
