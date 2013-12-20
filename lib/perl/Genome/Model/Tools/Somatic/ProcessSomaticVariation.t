#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;

use Genome::Utility::Test qw/compare_ok/;
my $pkg = 'Genome::Model::Tools::Somatic::ProcessSomaticVariation';

use_ok($pkg);

my $data_dir = Genome::Utility::Test->data_dir_ok($pkg, "data");

subtest "bedFileToAnnoFile" => sub {
    my $input_file  = "$data_dir/H_NS-POET0092.var.SNV.bed";
    my $output_file = Genome::Sys->create_temp_file_path;

    my $annotation_file = Genome::Model::Tools::Somatic::ProcessSomaticVariation::bedFileToAnnoFile("$input_file", "$output_file");

    is($annotation_file, $output_file, "Annotation file written to the desired location");
    compare_ok("$data_dir/H_NS-POET0092.var.SNV", "$annotation_file", "Content of output file as expected");
};

#This produces warning - refactor to avoid
subtest "bedFileToAnnoFile_slashed_input" => sub {
    my $slashed_input_file  = "$data_dir/H_NS-POET0092.var.SNV.slashed.bed";
    my $output_file = Genome::Sys->create_temp_file_path;

    my $annotation_file = Genome::Model::Tools::Somatic::ProcessSomaticVariation::bedFileToAnnoFile("$slashed_input_file", "$output_file");

    is($annotation_file, $output_file, "Annotation file written to the desired location");
    compare_ok("$data_dir/H_NS-POET0092.var.SNV", "$annotation_file", "Content of output file as expected");
};

subtest "annoFileToBedFile" => sub {
    my $input_file = "$data_dir/H_NS-POET0092.var.SNV";
    my $output_file = Genome::Sys->create_temp_file_path;

    my $bed_file = Genome::Model::Tools::Somatic::ProcessSomaticVariation::annoFileToBedFile("$input_file", $output_file);

    is($bed_file, $output_file, "BED file written to the desired location");
    compare_ok("$data_dir/H_NS-POET0092.var.SNV.bed", "$bed_file", "Content of output file as expected");
};

subtest "annoFileToSlashedBedFile" => sub {
    my $input_file = "$data_dir/H_NS-POET0092.var.SNV";
    my $output_file = Genome::Sys->create_temp_file_path;

    my $bed_file = Genome::Model::Tools::Somatic::ProcessSomaticVariation::annoFileToSlashedBedFile("$input_file", $output_file);

    is($bed_file, $output_file, "BED file written to the desired location");
    compare_ok("$data_dir/H_NS-POET0092.var.SNV.slashed.bed", "$bed_file", "Content of output file as expected");

};

done_testing();
