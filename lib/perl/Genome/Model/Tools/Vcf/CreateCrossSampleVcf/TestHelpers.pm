package Genome::Model::Tools::Vcf::CreateCrossSampleVcf::TestHelpers;

use strict;
use warnings;

use above 'Genome';
use Test::More;
use File::Spec;
use Genome::Utility::Test 'compare_ok';
use Genome::File::Vcf::Reader;
use Genome::File::Vcf::Differ;

use Exporter 'import';

our @EXPORT_OK = qw(
    get_test_dir
    create_test_builds
    get_roi_list
    test_cmd
);

sub get_test_dir {
    my ($class, $version) = @_;

    my $test_dir = File::Spec->join(Genome::Utility::Test->data_dir($class), "v-$version");
    diag "Test data located at $test_dir\n";
    return $test_dir;
}

sub get_reference_sequence_build {
    my $reference_sequence_build = Genome::Model::Build->get(106942997);
    ok($reference_sequence_build, "Found Reference Sequence Build") || die;
    return $reference_sequence_build;
}

sub define_test_classes {
    class Genome::Model::Test {
        is => 'Genome::ModelDeprecated',
    };

    class Genome::ProcessingProfile::Test {
        is => 'Genome::ProcessingProfile',
        has => [
            snv_detection_strategy => {
                is => 'String',
                is_calculated => 1,
                calculate => q( 'samtools r963 filtered by snp-filter v1 then false-positive v1 [--max-mm-qualsum-diff 100 --bam-readcount-version 0.4 --bam-readcount-min-base-quality 15] unique union varscan 2.2.9 [--min-coverage 3 --min-var-freq 0.20 --p-value 0.10 --strand-filter 1 --map-quality 10] filtered by false-positive v1 [--max-mm-qualsum-diff 100 --bam-readcount-version 0.4 --bam-readcount-min-base-quality 15]' ),
            },
        ],
    };

    class Genome::Model::Build::Test {
        is => 'Genome::Model::Build',
        has_optional => [
            reference_sequence_build => {
                is => 'Genome::Model::Build',
            },
            whole_rmdup_bam_file => {
                is => 'File',
            },
            get_snvs_vcf => {
                is => 'File',
            },
            get_indels_vcf => {
                is => 'File',
            },
        ],
    }
}

sub create_test_sample {
    my ($sample_name) = @_;

    my $sample = Genome::Sample->create(
        name => $sample_name,
    );
    ok($sample, sprintf('created test sample with name: %s, id: %s',
            $sample->name, $sample->id)) or die;
    return $sample;
}

sub create_test_pp {
    my ($pp_name) = @_;

    my $pp = Genome::ProcessingProfile::Test->create(
        name => $pp_name,
    );
    ok($pp, sprintf('created test pp with name: %s, id: %s',
            $pp->name, $pp->id)) or die;
    return $pp;
}

sub create_test_model {
    my ($sample, $pp, $name) = @_;

    my $model = Genome::Model::Test->create(
        subject_id => $sample->id,
        subject_class_name => $sample->class,
        processing_profile_id => $pp->id,
        name => $name
    );
    ok($model, sprintf('created test model with name: %s, id: %s',
            $model->name, $model->id)) or die;
    return $model;
}

sub _create_test_build {
    my ($model, $bam_file, $snv_vcf, $indel_vcf, $reference_sequence_build) = @_;

    my $build = Genome::Model::Build::Test->create(
        model => $model,
        whole_rmdup_bam_file => $bam_file,
        get_snvs_vcf => $snv_vcf,
        get_indels_vcf => $indel_vcf,
        reference_sequence_build => $reference_sequence_build,
    );
    $build->the_master_event->event_status('Succeeded');
    $build->the_master_event->date_completed("2013-07-11 20:47:51");
    return $build;
}

sub create_test_builds {
    my ($test_dir) = @_;

    define_test_classes();
    my $rsb = get_reference_sequence_build();
    my $pp = create_test_pp('create_cross_sample_vcf_test_pp');
    my @builds;
    for my $i (1,2,3) {
        my $dir = File::Spec->join($test_dir, "build_$i");
        my $sample = create_test_sample("create_cross_sample_vcf_test_sample_$i");
        my $model = create_test_model($sample, $pp, "create_cross_sample_vcf_test_model_$i");
        my $build = _create_test_build(
            $model,
            File::Spec->join($dir, 'whole_rmdup_bam_file.bam'),
            File::Spec->join($dir, 'snvs_file.vcf.gz'),
            File::Spec->join($dir, 'indels_file.vcf.gz'),
            $rsb,
        );
        ok($build, "Created test build $i");
        ok(-s $build->whole_rmdup_bam_file, "Found bam_file for build $i at " . $build->whole_rmdup_bam_file);
        ok(-s $build->get_snvs_vcf, "Found snvs_vcf for build $i at " . $build->get_snvs_vcf);
        ok(-s $build->get_indels_vcf, "Found indel_vcf for build $i at " . $build->get_indels_vcf);
        ok($build->reference_sequence_build, "Found reference sequence build for build $i at " . $build->reference_sequence_build);
        push @builds, $build;
    }

    return @builds;
}


#construct the FeatureList needed
sub get_roi_list {
    my $test_data_directory = shift;

    my $region_file = join("/", $test_data_directory, "roi.bed");
    my $region_file_content_hash = Genome::Sys->md5sum($region_file);
    my $roi_list = Genome::FeatureList->create(
        name                => "test feature-list",
        format              => 'multi-tracked',
        content_type        => 'targeted',
        file_path           => $region_file,
        file_content_hash   => $region_file_content_hash,
        reference_id        => 108563338,
    );
    return $roi_list;
}

sub test_cmd {
    my ($variant_type, $version, $use_mg, $no_region_limiting) = @_;

    my $class = 'Genome::Model::Tools::Vcf::CreateCrossSampleVcf::CreateCrossSampleVcf' . ucfirst($variant_type);

    my $sr_class = $class . "::Result";
    use_ok($class);
    use_ok($sr_class);

    my $test_dir = get_test_dir($class, $version);

    my $roi_list;
    $roi_list = get_roi_list($test_dir, 'roi.bed') unless $no_region_limiting;
    my $joinx_version = "1.7";
    my $varscan_version = "2.3.6";
    my @input_builds = create_test_builds($test_dir);

    my %params = (
            roi_list => $roi_list,
            max_files_per_merge => 2,
            wingspan => 500,
            allow_multiple_processing_profiles => undef,
            joinx_version => $joinx_version,
            varscan_version => $varscan_version,
            builds => \@input_builds,
    );

    if ($no_region_limiting){
       delete $params{'roi_list'};
       delete $params{'wingspan'};
    }

    if ($variant_type ne 'indels') {
        delete $params{'varscan_version'};
    }

    my %sr_params = %params;

    if ($use_mg) {
        my @models = map {$_->model} @input_builds;
        my $mg = Genome::ModelGroup->create(
            name => 'Test Model Group',
            models => \@models,
        );
        $params{model_group} = $mg;
        delete $params{builds};
    }

    my $output_directory = Genome::Sys->create_temp_directory();
    $params{output_directory} = $output_directory;

    my $cmd = $class->create(%params);
    ok($cmd, "created CreateCrossSampleVcf object");

    ok($cmd->execute(), "executed CreateCrossSampleVcf");
    my $result = $cmd->final_result;
    ok(-s $result, "result of CreateCrossSampleVcf: $result exists");

    delete $sr_params{'output_directory'};
    my $sr = $sr_class->get_with_lock(%sr_params);
    ok($sr, "found software result for test1");
    is($sr, $cmd->software_result, "found software result via cmd for test1");

    my $expected_result = get_expected_result($variant_type, $test_dir);

    my $differ = Genome::File::Vcf::Differ->new($result, $expected_result);
    is($differ->diff, undef, "Found No differences between $result and (expected) $expected_result");
}

sub get_expected_result {
    my ($variant_type, $result_dir) = @_;
    my $expected_result = File::Spec->join($result_dir, "$variant_type.merged.vcf.gz");
    ok(-s $expected_result, "expected result exists: $expected_result");
    return $expected_result;
}


1;
