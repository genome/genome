#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Test::More tests => 4;

use_ok('Genome::Model::Tools::Sam::MarkDuplicates');

my $input = $ENV{GENOME_TEST_INPUTS} . '/Genome-Tools-Sam-MarkDuplicates/sample.bam';

# step 1: test 1 

my $tmp_dir = Genome::Sys->create_temp_directory();

my $output_file = Genome::Sys->create_temp_file_path("out.bam");
my $metrics_file = Genome::Sys->create_temp_file_path("out.metrics");
my $log_file = Genome::Sys->create_temp_file_path("out.log");

#uncomment to inspect output 
#$log_file->unlink_on_destroy(0);
#$output_file->unlink_on_destroy(0);
#$metrics_file->unlink_on_destroy(0);
#$log_file->unlink_on_destroy(0);

my $cmd_1 = Genome::Model::Tools::Sam::MarkDuplicates->create(file_to_mark=>$input,
                                                              marked_file=>$output_file,
                                                              metrics_file=>$metrics_file,
                                                              log_file=>$log_file,
                                                              tmp_dir=>$tmp_dir,
                                                              remove_duplicates=>1,
                                                              max_jvm_heap_size=>2,
                                                              dedup_params => "read_name_regex null",
                                                            );


ok($cmd_1, "created command");
ok($cmd_1->execute, "executed");
ok(-s $output_file, "output file is nonzero");
