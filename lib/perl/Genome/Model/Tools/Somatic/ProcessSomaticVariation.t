#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;
use Genome::Test::Factory::Model::SomaticVariation;
use Genome::Test::Factory::Model::ImportedAnnotation;
use Genome::Test::Factory::Build;
use File::Spec;
use Genome::Utility::Test qw/compare_ok/;

my $pkg = 'Genome::Model::Tools::Somatic::ProcessSomaticVariation';

use_ok($pkg);
use_ok("Genome::Test::Factory::Model::SomaticVariation");

my $TEST_DATA_VERSION = 4;

my $data_dir = Genome::Utility::Test->data_dir_ok($pkg, File::Spec->join($TEST_DATA_VERSION, "data"));
my $input_dir = Genome::Utility::Test->data_dir_ok($pkg, File::Spec->join($TEST_DATA_VERSION, "input"));
my $somatic_variation_build_data_dir = Genome::Utility::Test->data_dir_ok($pkg, File::Spec->join($TEST_DATA_VERSION, "somatic_variation_build_data"));

my $normal_model = Genome::Test::Factory::Model::ReferenceAlignment->setup_object();
ok($normal_model->isa("Genome::Model::ReferenceAlignment"), "Generated a reference alignment model for normal");
my $tumor_model  = Genome::Test::Factory::Model::ReferenceAlignment->setup_object(
    subject_id            => $normal_model->subject_id,
    processing_profile_id => $normal_model->processing_profile->id,
);
ok($tumor_model->isa("Genome::Model::ReferenceAlignment"), "Generated a reference alignment model for tumor");

my $annotation_model = Genome::Test::Factory::Model::ImportedAnnotation->setup_object();
ok($annotation_model->isa("Genome::Model::ImportedAnnotation"), "Generated an annotation model");
my $annotation_build = Genome::Test::Factory::Build->setup_object(
    model_id   => $annotation_model->id,
    name       => "NCBI-human.ensembl/67_37l_v2",
    status     => "Succeeded",
);
ok($annotation_build->isa("Genome::Model::Build::ImportedAnnotation"), "Generated an annotation build");
$annotation_build->name("NCBI-human.ensembl/67_37l_v2");
is($annotation_build->name, "NCBI-human.ensembl/67_37l_v2", "Annotation build name is 'NCBI-human.ensembl/67_37l_v2'");

my $somatic_variation_model = Genome::Test::Factory::Model::SomaticVariation->setup_object(
    normal_model     => $normal_model,
    tumor_model      => $tumor_model,
    annotation_build => $annotation_build,
);
ok($somatic_variation_model->isa("Genome::Model::SomaticVariation"), "Generated a somatic variation model");

my $somatic_variation_build = Genome::Test::Factory::Build->setup_object(
    model_id         => $somatic_variation_model->id,
    data_directory   => $somatic_variation_build_data_dir,
    status           => "Succeeded",
);
ok($somatic_variation_build->isa("Genome::Model::Build::SomaticVariation"), "Generated a somatic variation build");

my $output_dir = Genome::Sys->create_temp_directory;
my $process_somatic_variation = $pkg->create(
    somatic_variation_model => $somatic_variation_model,
    output_dir              => $output_dir,
    target_regions          => "$input_dir/target_regions.bed",
    required_snv_callers    => 2,
);
ok($process_somatic_variation->isa("Genome::Model::Tools::Somatic::ProcessSomaticVariation"), "Generated a process somatic variation object");



subtest "create_directories" => sub {
    my $result = $process_somatic_variation->create_directories();

    ok($result, "Directories successfully created");
    is($process_somatic_variation->_sample_name_dir, $somatic_variation_model->subject->name, "Sample name directory as expected");
    is($process_somatic_variation->_full_output_dir, $output_dir . "/" . $somatic_variation_model->subject->name, "Full output dir as expected");
};

my $full_output_dir = $process_somatic_variation->_full_output_dir;

subtest "stage_snv_file" => sub {
    my $snv_file = $process_somatic_variation->stage_snv_file();

    is($snv_file, $full_output_dir . "/snvs/snvs.hq.bed", "snv file name as expected");
    compare_ok($snv_file, "$data_dir/snvs/snvs.hq.bed", "Contents of snv file as expected");
};

subtest "stage_indel_file" => sub {
    my $indel_file = $process_somatic_variation->stage_indel_file();

    is($indel_file, $full_output_dir . "/indels/indels.hq.bed", "indels file name as expected");
    compare_ok($indel_file, "$data_dir/indels/indels.hq.bed", "Contents of indels file as expected");
};

subtest "cleanFile" => sub {
    my $snv_file = $data_dir . "/snvs/snvs.hq.bed";
    my $cleaned_file = $process_somatic_variation->cleanFile($snv_file, 'snvs');

    is($cleaned_file, $full_output_dir . "/snvs/snvs.hq.clean.bed", "clean snv file name as expected");
    compare_ok($cleaned_file, "$data_dir/snvs/snvs.hq.clean.bed", "Contents of feature list file as expected");
};


#Tests for filtering out the off-target regions, if target regions are available
subtest "get_or_create_featurelist_file" => sub {
    my $featurelist_file = $process_somatic_variation->get_or_create_featurelist_file();

    is($featurelist_file, $full_output_dir . "/featurelist", "Feature list file path as expected");
    compare_ok($featurelist_file, "$data_dir/featurelist", "Contents of feature list file as expected");
};

subtest "filter_off_target_regions" => sub {
    my $file_to_be_filtered = $data_dir . "/snvs/snvs.hq.clean.bed";

    my $filtered_file = $process_somatic_variation->_filter_off_target_regions($file_to_be_filtered, "snvs");
    is($filtered_file, $full_output_dir . "/snvs/snvs.hq.clean.ontarget.bed", "Filtered file path as expected");
    compare_ok($filtered_file, "$data_dir/snvs/snvs.hq.clean.ontarget.bed", "Contents of filtered file as expected");
};

subtest "removeUnsupportedSites" => sub {
    my $file_to_be_filtered = $data_dir . "/snvs/snvs.hq.clean.ontarget.filtered.bed.filteredReg";

    my $filtered_file = $process_somatic_variation->removeUnsupportedSites($file_to_be_filtered);
    is($filtered_file, $full_output_dir . "/snvs/snvs.hq.clean.ontarget.filtered.bed.filteredReg.gt2callers", "Filtered file path as expected");
    compare_ok($filtered_file, "$data_dir/snvs/snvs.hq.clean.ontarget.filtered.bed.filteredReg.gt2callers", "Contents of filtered file as expected");
#file should be empty
};

subtest "doAnnotation" => sub {
    my $file_to_annotate = $data_dir . "/snvs/snvs.hq.clean.ontarget.filtered.bed.filteredReg";

    my $annotated_file = $process_somatic_variation->doAnnotation($file_to_annotate, 'snvs');
    is($annotated_file, $full_output_dir . "/snvs/snvs.hq.clean.ontarget.filtered.bed.filteredReg.anno", "Filtered file path as expected");
};

#subtest "annoFileToBedFile" => sub {
#    my $input_file = "$data_dir/snvs_before_tiering.anno";
#    my $output_file = Genome::Sys->create_temp_file_path;
#
#    my $bed_file = Genome::Model::Tools::Somatic::ProcessSomaticVariation::annoFileToBedFile($input_file, $output_file);
#
#    is($bed_file, $output_file, "BED file written to the desired location");
#    compare_ok("$data_dir/snvs_before_tiering.bed", "$bed_file", "Content of output file as expected");
#};

# subtest "bedFileToAnnoFile" => sub {
    # my $input_file  = "$data_dir/snvs_after_tiering.bed";
    # my $output_file = Genome::Sys->create_temp_file_path;

    # my $annotation_file = Genome::Model::Tools::Somatic::ProcessSomaticVariation::bedFileToAnnoFile($input_file, $output_file);

    # is($annotation_file, $output_file, "Annotation file written to the desired location");
    # compare_ok("$data_dir/snvs_after_tiering.anno", "$annotation_file", "Content of output file as expected");
# };

# subtest "annoFileToSlashedBedFile" => sub {
    # my $input_file = "$data_dir/review_anno_before_annoFileToSlashedBedFile.anno";
    # my $output_file = Genome::Sys->create_temp_file_path;

    # my $bed_file = Genome::Model::Tools::Somatic::ProcessSomaticVariation::annoFileToSlashedBedFile($input_file, $output_file);

    # is($bed_file, $output_file, "BED file written to the desired location");
    # compare_ok("$data_dir/review_bed_after_annoFileToSlashedBedFile.bed", $bed_file, "Content of output file as expected");
# };

# subtest "doAnnotation_bed_input" => sub {
    # my $input_file = "$data_dir/snvs_before_annotation.bed";
    # my $output_file = Genome::Sys->create_temp_file_path;

    # my $annotation_file = Genome::Model::Tools::Somatic::ProcessSomaticVariation::doAnnotation(
        # $input_file,
        # "NCBI-human.ensembl/67_37l_v2",
        # $output_file
    # );

    # is($annotation_file, $output_file, "BED file written to the expected location");
    # compare_ok("$data_dir/snvs_before_tiering.anno", "$annotation_file", "Content of output file as expected");
# };

# subtest "addTiering_annotation_input" => sub {
    # my $input_file = "$data_dir/snvs_before_tiering.anno";
    # my $output_file = Genome::Sys->create_temp_file_path;

    # my $tiering_file = Genome::Model::Tools::Somatic::ProcessSomaticVariation::addTiering(
        # $input_file,
        # $tiering_files_dir,
        # $output_file
    # );

    # is($tiering_file, $output_file, "BED file written to the expected location");
    # compare_ok("$data_dir/snvs_after_tiering.bed", "$tiering_file", "Content of output file as expected");
# };

# subtest "addTiering_bed_input" => sub {
    # my $input_file = "$data_dir/snvs_before_tiering.bed";
    # my $output_file = Genome::Sys->create_temp_file_path;

    # my $tiering_file = Genome::Model::Tools::Somatic::ProcessSomaticVariation::addTiering(
        # $input_file,
        # $tiering_files_dir,
        # $output_file
    # );

    # is($tiering_file, $output_file, "BED file written to the expected location");
    # compare_ok("$data_dir/snvs_after_tiering.bed", "$tiering_file", "Content of output file as expected");
# };


done_testing();
