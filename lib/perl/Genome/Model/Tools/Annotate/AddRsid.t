#!/usr/bin/env genome-perl

use strict;
use warnings;
use File::Temp;
use above "Genome";
use Test::More tests => 6;

use_ok('Genome::Model::Tools::Annotate::AddRsid');
my $input_dir = $ENV{GENOME_TEST_INPUTS} . "/Genome-Model-Tools-Annotate-AddRsid";
ok (-d $input_dir, "Input directory exists");


my $tmp_file = Genome::Sys->create_temp_file_path();
my $cmd = Genome::Model::Tools::Annotate::AddRsid->create(anno_file => "$input_dir/test.annotate.top2", vcf_file => "$input_dir/test.vcf.gz", output_file => $tmp_file );
ok ($cmd);

my $rv = $cmd->execute;

ok ($rv);
ok (-e $tmp_file);
my $knowndiff = `diff $input_dir/test.out $tmp_file`;
ok($knowndiff eq '', "known output as expected");


