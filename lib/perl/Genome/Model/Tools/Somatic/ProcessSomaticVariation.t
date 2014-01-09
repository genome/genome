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
use_ok("Genome::Test::Factory::Model::SomaticVariation");

my $data_dir = Genome::Utility::Test->data_dir_ok($pkg, "data");
my $input_dir = Genome::Utility::Test->data_dir_ok($pkg, "input");
my $tiering_files_dir = Genome::Utility::Test->data_dir_ok($pkg, "tiering_bed_files_v3");
my $somatic_variation_build_data_dir = Genome::Utility::Test->data_dir_ok($pkg, "somatic_variation_build_data");

my $normal_model = Genome::Test::Factory::Model::ReferenceAlignment->setup_object();
ok($normal_model->isa("Genome::Model::ReferenceAlignment"), "Generated a reference alignment model for normal");
my $tumor_model  = Genome::Test::Factory::Model::ReferenceAlignment->setup_object(
    subject_id            => $normal_model->subject_id,
    processing_profile_id => $normal_model->processing_profile->id
);
ok($tumor_model->isa("Genome::Model::ReferenceAlignment"), "Generated a reference alignment model for tumor");

my $somatic_variation_model = Genome::Test::Factory::Model::SomaticVariation->setup_object(
    normal_model => $normal_model,
    tumor_model  => $tumor_model
);
ok($somatic_variation_model->isa("Genome::Model::SomaticVariation"), "Generated a somatic variation model");

my $somatic_variation_build = Genome::Test::Factory::Build->setup_object(
    model_id        => $somatic_variation_model->id,
    data_directory  => $somatic_variation_build_data_dir,
    status          => "Succeeded",
);
ok($somatic_variation_build->isa("Genome::Model::Build::SomaticVariation"), "Generated a somatic variation build");

my $output_dir = Genome::Sys->create_temp_directory;
my $process_somatic_variation = $pkg->create(
    somatic_variation_model => $somatic_variation_model,
    output_dir              => $output_dir,
    target_regions          => "$input_dir/target_regions.bed",
    filter_regions          => "$input_dir/filter.bed",
    filter_sites            => "$input_dir/filter_sites.bed",
);
ok($process_somatic_variation->isa("Genome::Model::Tools::Somatic::ProcessSomaticVariation"), "Generated a process somatic variation object");

subtest "create_directories" => sub {
    my $result = $process_somatic_variation->create_directories();

    ok($result, "Directories successfully created");
    is($process_somatic_variation->sample_name_dir, $somatic_variation_model->subject->name, "Sample name directory as expected");
    is($process_somatic_variation->full_output_dir, $output_dir . "/" . $somatic_variation_model->subject->name, "Full output dir as expected");
};

subtest "stage_snv_file" => sub {
    my $snv_file = $process_somatic_variation->stage_snv_file();

    is($snv_file, $process_somatic_variation->full_output_dir . "/snvs/snvs.hq.bed", "snv file name as expected");
    compare_ok($snv_file, "$data_dir/snvs/snvs.hq.bed", "Contents of snv file as expected");
};

subtest "stage_indel_file" => sub {
    my $indel_file = $process_somatic_variation->stage_indel_file();

    is($indel_file, $process_somatic_variation->full_output_dir . "/indels/indels.hq.bed", "indels file name as expected");
    compare_ok($indel_file, "$data_dir/indels/indels.hq.bed", "Contents of indels file as expected");
};

subtest "cleanFile" => sub {
    my $snv_file = $process_somatic_variation->full_output_dir . "/snvs/snvs.hq.bed";
#Test that the file exists;
    $DB::single=1;
    my $cleaned_file = Genome::Model::Tools::Somatic::ProcessSomaticVariation::cleanFile($snv_file);

    is($cleaned_file, $process_somatic_variation->full_output_dir . "/snvs/snvs.hq.clean.bed", "clean snv file name as expected");    
    compare_ok($cleaned_file, "$data_dir/snvs/snvs.hq.clean.bed", "Contents of feature list file as expected");
};

subtest "get_or_create_featurelist_file" => sub {
    my $featurelist_file = $process_somatic_variation->get_or_create_featurelist_file();

    is($featurelist_file, $process_somatic_variation->full_output_dir . "/featurelist", "Feature list file path as expected");
    compare_ok($featurelist_file, "$data_dir/featurelist", "Contents of feature list file as expected");
};

subtest "_filter_off_target_regions" => sub {
    my $file_to_be_filtered = $process_somatic_variation->full_output_dir . "/snvs/snvs.hq.clean.bed";
#Test that the file exists;
    compare_ok($file_to_be_filtered, "$data_dir/snvs/snvs.hq.clean.bed", "Contents of file to be filtered as expected");

    my $filtered_file = $process_somatic_variation->_filter_off_target_regions($file_to_be_filtered);
    is($filtered_file, $process_somatic_variation->full_output_dir . "/snvs/snvs.hq.clean.ontarget.bed", "Filtered file path as expected");
    compare_ok($filtered_file, "$data_dir/snvs/snvs.hq.clean.ontarget.bed", "Contents of filtered file as expected");
};


#Test for removing filter sites specified by the user

subtest "getFilterSites" => sub {
    my $filter_sites = $process_somatic_variation->getFilterSites($process_somatic_variation->filter_sites);
    my $temp_file = Genome::Sys->create_temp_file_path;
    my $temp_fh   = Genome::Sys->open_file_for_writing($temp_file);
    print $temp_fh Data::Dumper::Dumper($filter_sites);
    close ($temp_fh);
    compare_ok($temp_file, "$data_dir/filter_sites", "Filter sites as expected");
};

subtest "get_or_create_filter_file" => sub {
    my $filter_file = $process_somatic_variation->get_or_create_filter_file();

    is($filter_file, $process_somatic_variation->full_output_dir . "/filter", "Filter file path as expected");
    compare_ok($filter_file, "$data_dir/filter", "Contents of filter file as expected");
};

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

subtest "annoFileToSlashedBedFile" => sub {
    my $input_file = "$data_dir/review_anno_before_annoFileToSlashedBedFile.anno";
    my $output_file = Genome::Sys->create_temp_file_path;

    my $bed_file = Genome::Model::Tools::Somatic::ProcessSomaticVariation::annoFileToSlashedBedFile($input_file, $output_file);

    is($bed_file, $output_file, "BED file written to the desired location");
    compare_ok("$data_dir/review_bed_after_annoFileToSlashedBedFile.bed", $bed_file, "Content of output file as expected");
};

subtest "doAnnotation_bed_input" => sub {
    my $input_file = "$data_dir/snvs_before_annotation.bed";
    my $output_file = Genome::Sys->create_temp_file_path;

    my $annotation_file = Genome::Model::Tools::Somatic::ProcessSomaticVariation::doAnnotation(
        $input_file,
        "NCBI-human.ensembl/67_37l_v2",
        $output_file
    );

    is($annotation_file, $output_file, "BED file written to the expected location");
    compare_ok("$data_dir/snvs_before_tiering.anno", "$annotation_file", "Content of output file as expected");

};

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
