#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;
use Genome::Test::Factory::Model::SomaticVariation;
use Genome::Test::Factory::Model::ImportedAnnotation;
use Genome::Test::Factory::Build;
use Genome::Utility::Test qw/compare_ok/;

my $pkg = 'Genome::Model::Tools::Somatic::ProcessSomaticVariation';

use_ok($pkg);

my $data_dir = Genome::Utility::Test->data_dir_ok($pkg, "data");
my $tiering_files_dir = Genome::Utility::Test->data_dir_ok($pkg, "tiering_bed_files_v3");

subtest "annoFileToBedFile" => sub {
    my $input_file = "$data_dir/snvs_before_tiering.anno";
    my $output_file = Genome::Sys->create_temp_file_path;

    my $bed_file = Genome::Model::Tools::Somatic::ProcessSomaticVariation::annoFileToBedFile($input_file, $output_file);

    is($bed_file, $output_file, "BED file written to the desired location");
    compare_ok("$data_dir/snvs_before_tiering.bed", "$bed_file", "Content of output file as expected");
};

subtest "bedFileToAnnoFile" => sub {
    my $input_file  = "$data_dir/snvs_after_tiering.bed";
    my $output_file = Genome::Sys->create_temp_file_path;

    my $annotation_file = Genome::Model::Tools::Somatic::ProcessSomaticVariation::bedFileToAnnoFile($input_file, $output_file);

    is($annotation_file, $output_file, "Annotation file written to the desired location");
    compare_ok("$data_dir/snvs_after_tiering.anno", "$annotation_file", "Content of output file as expected");
};

# subtest "annoFileToSlashedBedFile" => sub {
    # my $input_file = "$data_dir/H_NS-POET0092.var.SNV";
    # my $output_file = Genome::Sys->create_temp_file_path;

    # my $bed_file = Genome::Model::Tools::Somatic::ProcessSomaticVariation::annoFileToSlashedBedFile($input_file, $output_file);

    # is($bed_file, $output_file, "BED file written to the desired location");
    # compare_ok("$data_dir/H_NS-POET0092.var.SNV.slashed.bed", "$bed_file", "Content of output file as expected");
# };

subtest "addTiering_annotation_input" => sub {
    my $input_file = "$data_dir/snvs_before_tiering.anno";
    my $output_file = Genome::Sys->create_temp_file_path;

    my $tiering_file = Genome::Model::Tools::Somatic::ProcessSomaticVariation::addTiering(
        $input_file,
        $tiering_files_dir,
        $output_file
    );

    is($tiering_file, $output_file, "BED file written to the expected location");
    compare_ok("$data_dir/snvs_after_tiering.bed", "$tiering_file", "Content of output file as expected");
};

subtest "addTiering_bed_input" => sub {
    my $input_file = "$data_dir/snvs_before_tiering.bed";
    my $output_file = Genome::Sys->create_temp_file_path;

    my $tiering_file = Genome::Model::Tools::Somatic::ProcessSomaticVariation::addTiering(
        $input_file,
        $tiering_files_dir,
        $output_file
    );

    is($tiering_file, $output_file, "BED file written to the expected location");
    compare_ok("$data_dir/snvs_after_tiering.bed", "$tiering_file", "Content of output file as expected");
};
done_testing();
