#!/usr/bin/env genome-perl

use strict;
use warnings;

use feature qw(state);

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

use File::Spec qw();
use List::MoreUtils qw(uniq);

use above 'Genome';

use Genome::Test::Factory::Individual;
use Genome::Test::Factory::InstrumentData::Imported;
use Genome::Test::Factory::InstrumentData::MergedAlignmentResult;
use Genome::Test::Factory::Model::SomaticValidation;
use Genome::Test::Factory::Model::SomaticVariation;
use Genome::Test::Factory::Model::ReferenceAlignment;
use Genome::Test::Factory::Model::RnaSeq;
use Genome::Test::Factory::Model::ClinSeq;
use Genome::Test::Factory::ProcessingProfile::RnaSeq;
use Genome::Test::Factory::Sample;

use Genome::Utility::Text;
use Genome::Utility::Test qw(compare_ok);

use constant INDIVIDUAL_NAME => 'H_TEST-ClinSeqGenomeCommands';
use constant FLOW_CELL_ID => 'TEST1ABXX';

use Test::More tests => 103;

my $temp_dir = Genome::Sys->create_temp_directory();
ok(-d $temp_dir, "created temp directory: $temp_dir") or die;

my $expected_output_dir = __FILE__ .  '.d/';
ok(-e $expected_output_dir, "Found test dir: $expected_output_dir") or die;

my ($individual, $model_group) = setup_test_data();
isa_ok($individual, 'Genome::Individual', 'created test individual');
isa_ok($model_group, 'Genome::ModelGroup', 'created test model-group');

#CLIN-SEQ UPDATE-ANALYSIS
#Test clin-seq update-analysis - make sure the following command correctly obtains two expected samples (this has been broken in the past)
my @sample_ids = map { $_->id } ($individual->samples)[1,2];
my @cmd = (qw(genome model clin-seq update-analysis), '--individual=' . $individual->name, '--samples=id in ["' . join('","', @sample_ids) . '"]', '--display-samples');
run_ok(\@cmd, 'genome model clin-seq update-analysis');

#GENOME SAMPLE LIST
@cmd = (qw(genome sample list --filter), sprintf('name like "%s%s"', INDIVIDUAL_NAME, '%'), '--show', 'id,name,common_name,tissue_desc,extraction_type,extraction_label');
run_ok(\@cmd, 'genome sample list1');

#GENOME MODEL CLIN-SEQ LIST
@cmd = (qw(genome model clin-seq list --filter), sprintf('model_groups.id=%s',$model_group->id), '--show', 'wgs_model.last_succeeded_build.id,wgs_model.last_succeeded_build.data_directory');
run_ok(\@cmd, 'genome model clin-seq list1');

@cmd = (qw(genome model clin-seq list --filter), sprintf('model_groups.id=%s',$model_group->id), '--style=tsv', '--show', 'id,name,wgs_model,tumor_rnaseq_model,subject.common_name');
run_ok(\@cmd, 'genome model clin-seq list2');

@cmd = (qw(genome model clin-seq list --style csv --filter), sprintf('model_groups.id=%s',$model_group->id), '--show', 'wgs_model.last_succeeded_build.normal_build.subject.name,wgs_model.last_succeeded_build.normal_build.whole_rmdup_bam_file');
run_ok(\@cmd, 'genome model clin-seq list3');

#GENOME MODEL SOMATIC-VARIATION LIST
@cmd = (qw(genome model somatic-variation list --filter), sprintf('group_ids=%s', $model_group->id), '--show', 'subject.individual_common_name,subject.name,id');
run_ok(\@cmd, 'genome model somatic-variation list1');

@cmd = (qw(genome model somatic-variation list --filter), sprintf('group_ids=%s', $model_group->id), '--show', 'subject.name,last_succeeded_build_directory', '--noheaders');
run_ok(\@cmd, 'genome model somatic-variation list2');

@cmd = (qw(genome model somatic-variation list), sprintf('model_groups.id=%s', $model_group->id), '--show', 'tumor_model.subject.name,tumor_model.subject.common_name', '--style=csv');
run_ok(\@cmd, 'genome model somatic-variation list3');

#GENOME MODEL RNA-SEQ LIST
my ($rna_seq_model) = grep { $_->isa('Genome::Model::RnaSeq') } $model_group->models;
@cmd = (qw(genome model rna-seq list --filter), sprintf('genome_model_id=%s', $rna_seq_model->id));
run_ok(\@cmd, 'genome model rna-seq list1');

@cmd = (qw(genome model rna-seq list), sprintf('group_ids=%s', $model_group->id), '--show', 'id,name,processing_profile,last_succeeded_build.id,last_succeeded_build.alignment_result.bam_file', '--style', 'tsv');
run_ok(\@cmd, 'genome model rna-seq list2');

#GENOME INSTRUMENT-DATA LIST SOLEXA
@cmd = (qw(genome instrument-data list solexa --show), 'id,flow_cell_id,lane,index_sequence,sample_name,library_name,clusters,read_length,bam_path', '--filter', sprintf('flow_cell_id=%s', FLOW_CELL_ID));
run_ok(\@cmd, 'genome instrument-data list1');

@cmd = (qw(genome instrument-data list solexa --filter), sprintf('sample_name=%s',$rna_seq_model->subject->name));
run_ok(\@cmd, 'genome instrument-data list2');

#GENOME MODEL-GROUP MEMBER LIST
@cmd = (qw(genome model-group member list --filter), sprintf('model.subclass_name=Genome::Model::ClinSeq,model_group_id=%s', $model_group->id), '--show', 'model.wgs_model.id,model.wgs_model.subject.individual_common_name,model.last_succeeded_build,model.last_succeeded_build.data_directory');
run_ok(\@cmd, 'genome model-group member list1');

#GENOME MODEL SOMATIC-VALIDATION LIST
@cmd = (qw(genome model somatic-validation list --filter), sprintf('model_groups.id=%s', $model_group->id), '--show', 'tumor_sample.individual_common_name,tumor_sample.name,last_complete_build.tumor_bam');
run_ok(\@cmd, 'genome somatic-validation list1');

$DB::single = 1;

my @expected_files = map { (/^$expected_output_dir\/?(.*)/)[0] } glob($expected_output_dir . '/*');
my @temp_files = map { (/^$temp_dir\/?(.*)/)[0] } glob($temp_dir . '/*');
my @files = uniq (@expected_files, @temp_files);
@files = grep { $_ =~ /\.(err|out)$/ } @files;
chomp @files;
is(scalar(@files), 28, 'found expected number of files');
for my $file (@files) {
    my $temp_file = File::Spec->join($temp_dir, $file);
    ok(-f $temp_file, "file exists in temp_dir: $file");
    my $expected_file = File::Spec->join($expected_output_dir, $file);
    ok(-f $expected_file, "file exists in expected_output_dir: $file");
    compare_ok($temp_file, $expected_file,
        name => "compared $file",
    );
}


sub run_ok {
    my $cmd = shift;
    my $test = shift;

    my @argv = @$cmd;
    shift @argv if $argv[0] eq 'genome';
    my $file_name = $test;
    $file_name =~ s/ /-/g;
    my $file_path = File::Spec->join($temp_dir, Genome::Utility::Text::sanitize_string_for_filesystem($file_name));
    my $rv;
    do {
        local *STDOUT = Genome::Sys->open_file_for_writing("$file_path.out");
        local *STDERR = Genome::Sys->open_file_for_writing("$file_path.err");
        $rv = Genome::Command->_cmdline_run(@argv);
    };
    ok(!$rv, sprintf('tested %s', $test)) or diag join(' ', @$cmd);
}

#mock all of the various sorts of data needed to run this test
sub setup_test_data {
    my $individual = Genome::Test::Factory::Individual->setup_object(
        id => _next_id(),
        name => INDIVIDUAL_NAME
    );
    my $tumor_genotype = Genome::Test::Factory::InstrumentData::Imported->setup_object(
        id => _next_id(),
    );
    my $tumor_sample = Genome::Test::Factory::Sample->setup_object(
        id => _next_id(),
        name => join('-', INDIVIDUAL_NAME, 'tumor'),
        common_name => 'tumor',
        extraction_type => 'genomic dna',
        extraction_label => 'TESTY',
        tissue_desc => 'liver',
        default_genotype_data_id => $tumor_genotype->id,
        source_id => $individual->id,
    );


    my $normal_genotype = Genome::Test::Factory::InstrumentData::Imported->setup_object(
        id => _next_id(),
    );
    my $normal_sample = Genome::Test::Factory::Sample->setup_object(
        id => _next_id(),
        name => join('-', INDIVIDUAL_NAME, 'normal'),
        common_name => 'normal',
        extraction_type => 'genomic dna',
        extraction_label => 'TESTY control',
        tissue_desc => 'skin',
        default_genotype_data_id => $normal_genotype->id,
        source_id => $individual->id,
    );

    my $rna_sample = Genome::Test::Factory::Sample->setup_object(
        id => _next_id(),
        name => join('-', INDIVIDUAL_NAME, 'tumor_cdna'),
        common_name => 'tumor',
        extraction_type => 'rna',
        extraction_label => 'TESTY rna',
        tissue_desc => 'liver',
        source_id => $individual->id,
    );

    for my $sample ($tumor_sample, $normal_sample, $rna_sample) {
        for(1..3) {
            my $library = Genome::Test::Factory::Library->setup_object(
                id => _next_id(),
                name => $sample->name . "-lib$_",
                sample_id => $sample->id,
            );
            for(1..2) {
                my $index_sequence = _next_id();
                my $id = _next_id();
                Genome::Test::Factory::InstrumentData::Solexa->setup_object(
                    flow_cell_id => FLOW_CELL_ID,
                    lane => $_,
                    index_sequence => $index_sequence,
                    subset_name => join('-', $_, $index_sequence),
                    id => $id,
                    library_id => $library->id,
                    read_length => 100,
                    clusters => (12345678 + $id),
                    old_median_insert_size => '50',
                    old_sd_above_insert_size => '85',
                    bam_path => sprintf('/tmp/fake-csf-dir/%s.bam', $id),
                );
            }
        }
    }

    my $model_group = Genome::ModelGroup->create(
        name => 'TestGenomeCommands',
        id => _next_id(),
    );

    my $rnaseq_processing_profile = Genome::Test::Factory::ProcessingProfile::RnaSeq->setup_object(
        id => _next_id(),
        name => 'TestGenomeCommands RnaSeq ProcessingProfile',
        read_aligner_name => 'tophat',
        read_aligner_version => '2.0.8',
    );
    my $rnaseq_model = Genome::Test::Factory::Model::RnaSeq->setup_object(
        subject_id => $rna_sample->id,
        instrument_data => [$rna_sample->instrument_data],
        id => _next_id(),
        processing_profile_id => $rnaseq_processing_profile->id,
        name => $rna_sample->name . '.rnaseq',
    );
    my $rnaseq_build = Genome::Test::Factory::Build->setup_object(
        model_id => $rnaseq_model->id,
        id => _next_id(),
        status => 'Succeeded',
    );
    my $rnaseq_alignment = Genome::Test::Factory::InstrumentData::MergedAlignmentResult->setup_object(
        output_dir => '/tmp/fake-rna-path',
        id => _next_id(),
    );
    Genome::SoftwareResult::User->__define__(
        label => 'uses',
        user => $rnaseq_build,
        software_result => $rnaseq_alignment,
    );


    my $somval_model = Genome::Test::Factory::Model::SomaticValidation->setup_object(
        subject_id => $individual->id,
        tumor_sample_id => $tumor_sample->id,
        normal_sample_id => $normal_sample->id,
        instrument_data => [$tumor_sample->instrument_data, $normal_sample->instrument_data],
        id => _next_id(),
        name => $individual->name . '.somatic-validation',
    );
    my $somval_build = Genome::Test::Factory::Build->setup_object(
        id => _next_id(),
        model_id => $somval_model->id,
        status => 'Succeeded',
    );
    my $somval_alignment = Genome::Test::Factory::InstrumentData::MergedAlignmentResult->setup_object(
        output_dir => '/tmp/fake-somval-path',
        id => _next_id(),
    );
    Genome::SoftwareResult::User->__define__(
        label => 'merged_alignment',
        user => $somval_build,
        software_result => $somval_alignment,
    );

    my $tumor_refalign_model = Genome::Test::Factory::Model::ReferenceAlignment->setup_object(
        subject_id => $tumor_sample->id,
        instrument_data => [$tumor_sample->instrument_data],
        id => _next_id(),
        name => $tumor_sample->name . '.refalign',
    );
    Genome::Test::Factory::Build->setup_object(
        id => _next_id(),
        model_id => $tumor_refalign_model->id,
        status => 'Succeeded',
    );

    my $normal_refalign_model = Genome::Test::Factory::Model::ReferenceAlignment->setup_object(
        subject_id => $normal_sample->id,
        instrument_data => [$normal_sample->instrument_data],
        id => _next_id(),
        name => $normal_sample->name . '.refalign',
    );
    my $normal_refalign_build = Genome::Test::Factory::Build->setup_object(
        id => _next_id(),
        model_id => $normal_refalign_model->id,
        status => 'Succeeded',
    );
    my $normal_alignment = Genome::Test::Factory::InstrumentData::MergedAlignmentResult->setup_object(
        output_dir => '/tmp/fake-normal-path',
        id => _next_id(),
    );
    Genome::SoftwareResult::User->__define__(
        label => 'uses',
        user => $normal_refalign_build,
        software_result => $normal_alignment,
    );

    my $somvar_model = Genome::Test::Factory::Model::SomaticVariation->setup_object(
        subject_id => $tumor_sample->id,
        tumor_model => $tumor_refalign_model,
        normal_model => $normal_refalign_model,
        id => _next_id(),
        name => $tumor_sample->name . '.somatic-variation',
    );
    Genome::Test::Factory::Build->setup_object(
        id => _next_id(),
        model_id => $somvar_model->id,
        status => 'Succeeded',
        data_directory => '/tmp/somvar-build-dir',
    );

    my $clinseq_model = Genome::Test::Factory::Model::ClinSeq->setup_object(
        subject_id => $individual->id,
        wgs_model => $somvar_model,
        tumor_rnaseq_model => $rnaseq_model,
        id => _next_id(),
        name => $individual->name . '.clin-seq',
    );
    Genome::Test::Factory::Build->setup_object(
        id => _next_id(),
        model_id => $clinseq_model->id,
        status => 'Succeeded',
        data_directory => '/tmp/clinseq-build-dir',
    );

    $model_group->assign_models($rnaseq_model, $somval_model, $tumor_refalign_model, $normal_refalign_model, $somvar_model, $clinseq_model);

    return ($individual, $model_group);
}

sub _next_id {
    state $next_id = -1;
    return $next_id--;
}
