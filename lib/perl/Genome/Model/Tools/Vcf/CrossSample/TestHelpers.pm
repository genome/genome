package Genome::Model::Tools::Vcf::CrossSample::TestHelpers;

use strict;
use warnings;

use above 'Genome';
use Test::More;
use File::Spec;
use Genome::Utility::Test 'compare_ok';
use Genome::File::Vcf::Reader;
use Genome::File::Vcf::Differ;
use File::Copy qw(cp);
use File::Slurp qw(write_file);
use Memoize;

use Exporter 'import';

our @EXPORT_OK = qw(
    get_test_dir
    create_test_builds
    get_roi_list
    test_indel_cmd
);

sub get_test_dir {
    my ($class, $version) = @_;

    my $test_dir = File::Spec->join(Genome::Utility::Test->data_dir($class), "v-$version");
    note "Test data located at $test_dir\n";
    return $test_dir;
}

sub get_reference_sequence_build {
    my $reference_sequence_build = Genome::Model::Build->get(106942997);
    _ensure_object($reference_sequence_build);
    return $reference_sequence_build;
}

sub define_test_classes {
    class Genome::Model::Test {
        is => 'Genome::ModelDeprecated',
    };

    class Genome::ProcessingProfile::Test {
        is => 'Genome::ProcessingProfile',
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
    };
    return 1;
}
Memoize::memoize('define_test_classes');

sub create_test_sample {
    my ($sample_name) = @_;

    my $sample = Genome::Sample->create(
        name => $sample_name,
    );
    _ensure_object($sample);
    return $sample;
}

sub create_test_pp {
    my ($pp_name) = @_;

    my $pp = Genome::ProcessingProfile::Test->create(
        name => $pp_name,
    );
    _ensure_object($pp);
    return $pp;
}
Memoize::memoize('create_test_pp');

sub create_test_model {
    my ($sample, $pp, $name) = @_;

    my $model = Genome::Model::Test->create(
        subject_id => $sample->id,
        subject_class_name => $sample->class,
        processing_profile_id => $pp->id,
        name => $name
    );
    _ensure_object($model);
    return $model;
}

sub _create_test_build {
    my ($model, $bam_file, $snv_vcf, $indel_vcf, $reference_sequence_build, $date_completed) = @_;

    my $build = Genome::Model::Build::Test->create(
        model => $model,
        whole_rmdup_bam_file => $bam_file,
        get_snvs_vcf => $snv_vcf,
        get_indels_vcf => $indel_vcf,
        reference_sequence_build => $reference_sequence_build,
    );
    $build->the_master_event->event_status('Succeeded');
    $build->the_master_event->date_completed($date_completed);
    return $build;
}

sub create_test_builds {
    my ($test_dir) = @_;

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
            $i,
        );
        _ensure_object($build);

        my $bam_index = File::Spec->join($dir, "whole_rmdup_bam_file.bam.bai");
        for my $file ($build->whole_rmdup_bam_file,
                      $bam_index,
                      $build->get_snvs_vcf,
                      $build->get_indels_vcf) {
            _ensure_file($file);
        }

        _ensure_object($build->reference_sequence_build);
        push @builds, $build;
    }

    return @builds;
}
Memoize::memoize('create_test_builds');

sub _ensure_object {
    my $object = shift;
    if ($object) {
        note "object (" . ref($object) . ") exists as expected\n";
    } else {
        die "object (" . ref($object) . ") does not exist with size";
    }
}

sub _ensure_file {
    my $file = shift;
    if (-s $file) {
        note "File ($file) exists as expected\n";
    } else {
        die "File ($file) does not exist with size";
    }
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
    _ensure_object($roi_list);
    return $roi_list;
}

sub test_indel_cmd {
    my ($version, $use_mg, $region_limiting) = @_;
    my $class = 'Genome::Model::Tools::Vcf::CrossSample::Indel';
    use_ok($class);
    define_test_classes();

    my $test_dir = get_test_dir($class, $version);
    my $output_dir = Genome::Sys->create_temp_directory();
    my %params = (
            roi_list => get_roi_list($test_dir, 'roi.bed'),
            wingspan => 500,
            allow_multiple_processing_profiles => 0,
            joinx_version => '1.8',
            varscan_version => '2.3.6',
            output_directory => $output_dir,
    );

    my @input_builds = create_test_builds($test_dir);
    $params{builds} = \@input_builds;

    if ($use_mg) {
        my @models = map {$_->model} @input_builds;
        my $mg = Genome::ModelGroup->get_or_create(
            name => 'Test Model Group',
            models => \@models,
        );
        _ensure_object($mg);
        $params{model_group} = $mg;
        delete $params{builds};
    }

    my $output_tsv = File::Spec->join($output_dir, 'inputs.tsv');
    my $output_source_path = File::Spec->join($output_dir, 'source_path.txt');
    Sub::Install::reinstall_sub({
        into => 'Genome::Model::Tools::Park::Base',
        as => '_run_rex_process',
        code => sub {
            my ($self, $source_path, $inputs_filename) = @_;

            cp($inputs_filename, $output_tsv);
            write_file($output_source_path, $source_path);
            return "/v1/processes/1/";
        },
    });

    my $expected_inputs;
    my $expected_source;
    if ($region_limiting){
       $expected_inputs = File::Spec->join($test_dir, 'expected_region_limited.tsv');
       $expected_source = File::Spec->join($test_dir, 'expected-source_path_region_limited.txt');
    } else {
       delete $params{'roi_list'};
       delete $params{'wingspan'};
       $expected_inputs = File::Spec->join($test_dir, 'expected.tsv');
       $expected_source = File::Spec->join($test_dir, 'expected-source_path.txt');
    }

    my $cmd = $class->create(%params);
    ok($cmd, "created CrossSample Indel object");

    ok($cmd->execute(), "executed CrossSample Indel");

    my %compare_args = (
        replace => [ [ qr(^region_bed_file\t.*$) => "region_bed_file\tSOMEPATH"] ],
    );

    compare_ok($output_tsv, $expected_inputs, "expected inputs file ($expected_inputs) matches what we made ($output_tsv)", %compare_args);
    compare_ok($output_source_path, $expected_source, "expected source-path file ($expected_source) matches what we made ($output_source_path)");
}

1;
