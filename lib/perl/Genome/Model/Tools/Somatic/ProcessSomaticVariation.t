#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;
use Genome::Test::Factory::Model::SomaticVariation;
use Genome::Test::Factory::Model::ImportedAnnotation;
use Genome::Test::Factory::Model::ReferenceSequence;
use Genome::Test::Factory::InstrumentData::Solexa;
use Genome::Test::Factory::Build;
use File::Spec;
use Genome::Utility::Test qw/compare_ok/;
use Sub::Install qw();

my $pkg = 'Genome::Model::Tools::Somatic::ProcessSomaticVariation';

use_ok($pkg);
use_ok("Genome::Test::Factory::Model::SomaticVariation");

my $TEST_DATA_VERSION = 7;

my $main_dir = Genome::Utility::Test->data_dir_ok($pkg, $TEST_DATA_VERSION);
my $data_dir = Genome::Utility::Test->data_dir_ok($pkg, File::Spec->join($TEST_DATA_VERSION, "data"));
my $input_dir = Genome::Utility::Test->data_dir_ok($pkg, File::Spec->join($TEST_DATA_VERSION, "input"));
my $somatic_variation_build_data_dir = Genome::Utility::Test->data_dir_ok($pkg, File::Spec->join($TEST_DATA_VERSION, "somatic_variation_build_data"));
my $annotation_build_data_dir = Genome::Utility::Test->data_dir_ok($pkg, File::Spec->join($TEST_DATA_VERSION, "annotation_build_data"));
my $instrument_data_dir = Genome::Utility::Test->data_dir_ok($pkg, File::Spec->join($TEST_DATA_VERSION, "instrument_data"));

my $normal_instrument_data = Genome::Test::Factory::InstrumentData::Solexa->setup_object(
    bam_path => File::Spec->join($instrument_data_dir, "138572668_normal.bam"),
);

my $normal_model = Genome::Test::Factory::Model::ReferenceAlignment->setup_object();
ok($normal_model->isa("Genome::Model::ReferenceAlignment"), "Generated a reference alignment model for normal");
my $normal_build = Genome::Test::Factory::Build->setup_object(
    model_id         => $normal_model->id,
    status           => "Succeeded",
);

Sub::Install::reinstall_sub({
    into => "Genome::Model::Build::ReferenceAlignment",
    as   => 'merged_alignment_result',
    code => sub {
        my $self = shift;
        my $value = shift;

        if (defined($value)) {
            $self->{__merged_alignment_result} = $value;
        }
        return $self->{__merged_alignment_result};
    },
});
$normal_build->merged_alignment_result($normal_instrument_data);
ok(-s $normal_build->merged_alignment_result->bam_path, "Normal bam path correct");
ok($normal_build->isa("Genome::Model::Build::ReferenceAlignment"), "Generated a normal build");

my $tumor_reference_model = Genome::Test::Factory::Model::ReferenceSequence->setup_object();
my $tumor_reference_build = Genome::Test::Factory::Build->setup_object(
    model_id => $tumor_reference_model->id,
    data_directory => $main_dir,
);

my $tumor_instrument_data = Genome::Test::Factory::InstrumentData::Solexa->setup_object(
    bam_path => File::Spec->join($instrument_data_dir, "139691307_tumor.bam"),
);
my $tumor_model  = Genome::Test::Factory::Model::ReferenceAlignment->setup_object(
    subject_id            => $normal_model->subject_id,
    processing_profile_id => $normal_model->processing_profile->id,
    reference_sequence_build => $tumor_reference_build,
);
ok($tumor_model->isa("Genome::Model::ReferenceAlignment"), "Generated a reference alignment model for tumor");
my $tumor_build = Genome::Test::Factory::Build->setup_object(
    model_id         => $tumor_model->id,
    status           => "Succeeded",
);
$tumor_build->merged_alignment_result($tumor_instrument_data);
ok(-s $tumor_build->merged_alignment_result->bam_path, "Tumor bam path correct");
ok($tumor_build->isa("Genome::Model::Build::ReferenceAlignment"), "Generated a turmor build");

my $annotation_model = Genome::Test::Factory::Model::ImportedAnnotation->setup_object();
ok($annotation_model->isa("Genome::Model::ImportedAnnotation"), "Generated an annotation model");
my $annotation_build = Genome::Test::Factory::Build->setup_object(
    model_id       => $annotation_model->id,
    name           => "NCBI-human.ensembl/67_37l_v2",
    status         => "Succeeded",
    data_directory => $annotation_build_data_dir,
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
my $cmd = $pkg->create(
    somatic_variation_model => $somatic_variation_model,
    output_dir              => $output_dir,
    target_regions          => "$input_dir/target_regions.bed",
);
ok($cmd->isa("Genome::Model::Tools::Somatic::ProcessSomaticVariation"), "Generated a process somatic variation object");

ok($cmd->execute(), 'Command executed');

ok(-s $cmd->report, 'Found "report" output: ' . $cmd->report);
ok(-s $cmd->report_xls, 'Found "report.xls" output'  .  $cmd->report_xls);

compare_ok($cmd->report, File::Spec->join($data_dir, 'snvs.indels.annotated'), 'report is as expected');

done_testing();
