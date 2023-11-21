package Genome::Model::SomaticVariation::Command::TestHelpers;

use strict;
use warnings;

use above 'Genome';
use Test::More;
use Genome::Test::Factory::Model::SomaticVariation;
use Genome::Test::Factory::Model::ImportedAnnotation;
use Genome::Test::Factory::Model::ImportedReferenceSequence;
use Genome::Test::Factory::InstrumentData::MergedAlignmentResult;
use Genome::Test::Factory::Build;
use File::Spec;
use Genome::Utility::Test qw/compare_ok/;
use Sub::Install qw();
use IPC::System::Simple qw(capture);
use Sub::Override;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(
    create_test_objects
    run_test
);

sub create_test_objects {
    my $main_dir = shift;

    my $annotation_build_data_dir = File::Spec->join($main_dir, "annotation_build_data");
    my $instrument_data_dir = File::Spec->join($main_dir, "instrument_data");

    my $reference_model = Genome::Test::Factory::Model::ImportedReferenceSequence->setup_object();
    my $reference_build = Genome::Test::Factory::Build->setup_object(
        model_id => $reference_model->id,
        data_directory => $main_dir,
    );
    $reference_build->name('GRCh37-lite-build37');

    my $normal_merged_alignment_result = Genome::Test::Factory::InstrumentData::MergedAlignmentResult->setup_object(
        output_dir => $instrument_data_dir,
        id => '138572668_normal'
    );

    my $normal_model = Genome::Test::Factory::Model::ReferenceAlignment->setup_object(
        reference_sequence_build => $reference_build,
    );
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
    $normal_build->merged_alignment_result($normal_merged_alignment_result);
    ok(-s $normal_build->whole_rmdup_bam_file, "Normal bam path correct");
    ok($normal_build->isa("Genome::Model::Build::ReferenceAlignment"), "Generated a normal build");



    my $tumor_merged_alignment_result = Genome::Test::Factory::InstrumentData::MergedAlignmentResult->setup_object(
        output_dir => $instrument_data_dir,
        id => '139691307_tumor'
    );
    my $tumor_model  = Genome::Test::Factory::Model::ReferenceAlignment->setup_object(
        subject_id            => $normal_model->subject_id,
        processing_profile_id => $normal_model->processing_profile->id,
        reference_sequence_build => $reference_build,
    );
    ok($tumor_model->isa("Genome::Model::ReferenceAlignment"), "Generated a reference alignment model for tumor");
    my $tumor_build = Genome::Test::Factory::Build->setup_object(
        model_id         => $tumor_model->id,
        status           => "Succeeded",
    );
    $tumor_build->merged_alignment_result($tumor_merged_alignment_result);
    ok(-s $tumor_build->whole_rmdup_bam_file, "Tumor bam path correct");
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
    $somatic_variation_model->subject->name("H_NS-POET0092-4");

    my $somatic_variation_build = Genome::Test::Factory::Build->setup_object(
        model_id         => $somatic_variation_model->id,
        data_directory   => '',
        status           => "Succeeded",
    );
    ok($somatic_variation_build->isa("Genome::Model::Build::SomaticVariation"), "Generated a somatic variation build");

    my $override = Sub::Override->new('Genome::Model::Build::SomaticVariation::calculate_estimated_kb_usage', sub { return 102400; });
    my $somatic_variation_build_data_dir = $somatic_variation_build->get_or_create_data_directory();
    $override->restore();

    Genome::Sys->rsync_directory(
        source_directory => File::Spec->join($main_dir, "somatic_variation_build_data"),
        target_directory => $somatic_variation_build_data_dir,
    );

    return $somatic_variation_build;
}

sub run_test {
    my $pkg = shift;
    my $main_dir = shift;
    my $output_exists = shift;
    my %params = @_;

    my $cmd = $pkg->create(
        %params,
    );
    ok($cmd->isa("Genome::Model::SomaticVariation::Command::CreateReport"), "Generated a somatic variation create report object");

    ok($cmd->execute(), 'Command executed');

    if ($output_exists) {
        test_output_exists($cmd, $main_dir);
    }
    else {
        test_no_output_exists($cmd);
    }
}

sub test_output_exists {
    my $cmd = shift;
    my $main_dir = shift;

    ok(-s $cmd->report, 'Found "report" output: ' . $cmd->report);
    ok(-s $cmd->report_xls, 'Found "report.xls" output'  .  $cmd->report_xls);

    my $data_dir = File::Spec->join($main_dir, "data");

    compare_ok($cmd->report, File::Spec->join($data_dir, 'snvs.indels.annotated'), 'report is as expected');

    ok(-s $cmd->review_bed, 'Found "review.bed" output: ' . $cmd->review_bed);
    ok(-s $cmd->review_xml, 'Found "review.xml" output: ' . $cmd->review_xml);

    compare_ok($cmd->review_bed, File::Spec->join($data_dir, 'review.bed'), 'review.bed is as expected');

    ensure_symlinks_in_build_directory($cmd);
}

sub test_no_output_exists {
    my $cmd = shift;

    is($cmd->report, 'skipped', "Report file is 'skipped'");
    is($cmd->report_xls, 'skipped', "Report xls is 'skipped'");
    is($cmd->review_dir, 'skipped', "Review directory is 'skipped'");
    is($cmd->review_bed, 'skipped', "Review bed file is 'skipped'");
    is($cmd->review_xml, 'skipped', "Review cml file is 'skipped'");
}

sub ensure_symlinks_in_build_directory {
    my $cmd = shift;

    my $build = $cmd->somatic_variation_build;
    my $reports_dir = File::Spec->join($build->data_directory, 'reports');

    my %SYMLINKS = (
        'snvs.indels.annotated' => $cmd->report,
        'snvs.indels.annotated.xls' => $cmd->report_xls,
        'review.bed' => $cmd->review_bed,
        'review.xml' => $cmd->review_xml,
    );

    while (my ($file_in_build, $path_in_allocation) = each(%SYMLINKS)) {
        my $path_in_build = File::Spec->join($reports_dir, $file_in_build);

        my $result = capture([0,1], 'readlink', '-e', $path_in_build);
        chomp($result);
        is($result, $path_in_allocation, sprintf('found symlink to allocation: %s', $result));
    }
}

1;
