#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;
use Test::Exception;
use Sub::Override;
use File::Copy::Recursive qw(dircopy);
use Genome::Utility::Test qw(compare_ok);

BEGIN {
    $ENV{UR_USE_DUMMY_AUTOGENERATED_IDS} = 1;
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

use above 'Genome';
use Genome::Test::Factory::SoftwareResult::User;

my $pkg = 'Genome::InstrumentData::AlignmentResult::Merged';
use_ok($pkg);

# Override methods for testing so we don't need to worry about commit observers in the unit test
Sub::Install::install_sub({code => sub { my ($self, $alignment) = @_; $alignment->remove_bam; },
        into => 'Genome::InstrumentData::AlignmentResult::Merged', as => '_remove_per_lane_bam_post_commit'});

Sub::Install::install_sub({code => sub { return 1; },
        into => 'Genome::InstrumentData::AlignmentResult::Merged', as => '_size_up_allocation'});


#
# Gather up versions for the tools used herein
#
###############################################################################
my $aligner_name = "bwa";
my $aligner_tools_class_name = "Genome::Model::Tools::" . Genome::InstrumentData::AlignmentResult->_resolve_subclass_name_for_aligner_name($aligner_name);
my $alignment_result_class_name = "Genome::InstrumentData::AlignmentResult::" . Genome::InstrumentData::AlignmentResult->_resolve_subclass_name_for_aligner_name($aligner_name);

my $samtools_version = Genome::Model::Tools::Sam->default_samtools_version;
my $picard_version = Genome::Model::Tools::Picard->default_picard_version;

my $test_name = 'merged_unit_test';

my $aligner_version  = $aligner_tools_class_name->default_version;
my $aligner_label    = $aligner_name.$aligner_version;
$aligner_label       =~ s/\./\_/g;

my $expected_base_dir = Genome::Utility::Test->data_dir_ok($pkg, 'v1');
my $expected_shortcut_path = $expected_base_dir .'/bwa/',
my $expected_dir = $expected_base_dir .'/expected';

my $FAKE_INSTRUMENT_DATA_ID = -123456;

my $reference_model = Genome::Model::ImportedReferenceSequence->get(name => 'TEST-human');
ok($reference_model, "got reference model");

my $reference_build = $reference_model->build_by_version('1');
ok($reference_build, "got reference build");

my $result_users = Genome::Test::Factory::SoftwareResult::User->setup_user_hash(
    reference_sequence_build => $reference_build,
);

my $tmp_shortcut_path1 = Genome::Sys->create_temp_directory();
my $tmp_shortcut_path2 = Genome::Sys->create_temp_directory();
ok(dircopy($expected_shortcut_path, $tmp_shortcut_path1), 'Copied expected dir to path1');
ok(dircopy($expected_shortcut_path, $tmp_shortcut_path2), 'Copied expected dir to path2');

my @instrument_data    = generate_fake_instrument_data();
my @individual_results = generate_individual_alignment_results(@instrument_data);

my %parameters = (
     aligner_name                => $aligner_name,
     aligner_version             => $aligner_version,
     samtools_version            => $samtools_version,
     picard_version              => $picard_version,
     reference_build             => $reference_build,
     merger_name                 => 'picard',
     merger_version              => $picard_version,
     duplication_handler_name    => 'picard',
     duplication_handler_version => $picard_version,
     instrument_data_id          => [map($_->id, @instrument_data)],
     test_name                   => $test_name,
     instrument_data_segment     => [map {$_->id . ':A:2:read_group'} @instrument_data],
);

my $merged_alignment_result = $pkg->create(%parameters, _user_data_for_nested_results => $result_users);
isa_ok($merged_alignment_result, $pkg, 'produced merged alignment result');



subtest 'testing merge without filter_name' => sub {
    compare_ok($merged_alignment_result->bam_file, File::Spec->join($expected_dir, '-120573001.bam'),
        name => 'merged bam matches expected result',
        diag => 0,
    );
    compare_ok($merged_alignment_result->merged_alignment_bam_flagstat, File::Spec->join($expected_dir, '-120573001.bam.flagstat'), 'flagstat matches expected result');

    my @individual_alignments = $merged_alignment_result->collect_individual_alignments;
    is(scalar @individual_alignments, 2, 'got back expected number of alignments');

    for my $i (@individual_alignments) {
        ok(!defined($i->filter_name), 'filter_name is not defined as expected');
        my $id = $i->id;
        for my $remove_type ('', '.bai', '.md5') {
            my $basename = 'all_sequences.bam'.$remove_type;
            my $file = File::Spec->join($i->output_dir, $basename);
            ok(!-s $file, "$id $basename removed as expected");
        }
        for my $keep_type ('header', 'flagstat') {
            my $basename = 'all_sequences.bam.'.$keep_type;
            my $file = File::Spec->join($i->output_dir, $basename);
            ok(-s $file, "$id $basename kept as expected");
            if ($keep_type eq 'header') {
                my ($num) = $i->output_dir =~ /(\d)$/;
                my $expected_header = $expected_dir ."/$num".'_all_sequences.bam.header';
                compare_ok($file, $expected_header, "$id bam header file created ok as expected");
            }
        }
    }

    my $existing_alignment_result = $pkg->get_or_create(%parameters, users => $result_users);
    is($existing_alignment_result, $merged_alignment_result, 'got back the previously created result');
};

subtest 'testing merge with filter_name' => sub {
    my @filtered_params = (
        %parameters,
        filter_name => [$instrument_data[0]->id . ':forward-only', $instrument_data[1]->id . ':forward-only'],
    );

    my $existing_alignment_result = $pkg->get_or_create(%parameters, users => $result_users);
    is($existing_alignment_result, $merged_alignment_result, 'got back the previously created result');


    my $filtered_alignment_result = $pkg->get_or_create(@filtered_params, users => $result_users);
    isa_ok($filtered_alignment_result, $pkg, 'produced merged alignment result with filter applied');
    my @filtered_individual_alignments = $filtered_alignment_result->collect_individual_alignments;
    is(scalar @filtered_individual_alignments, 2, 'got back expected number of alignments');

    #same expected files since we faked the alignment results to use the same data
    compare_ok($filtered_alignment_result->bam_file, File::Spec->join($expected_dir, '-120573001.bam'),
        name => 'merged bam matches expected result',
        diag => 0,
    );
    compare_ok($filtered_alignment_result->merged_alignment_bam_flagstat, File::Spec->join($expected_dir, '-120573001.bam.flagstat'), 'flagstat matches expected result');

    for my $i (@filtered_individual_alignments) {
        is($i->filter_name, 'forward-only', 'filter_name is defined as expected');
    }

    my $existing_filtered_alignment_result = $pkg->get_or_create(@filtered_params, users => $result_users);
    isnt($filtered_alignment_result, $existing_alignment_result, 'produced a different result when filter applied');

    is($existing_filtered_alignment_result, $filtered_alignment_result, 'got back the previously created filtered result');

    my $gotten_alignment_result = $pkg->get_with_lock(%parameters, users => $result_users);
    is($gotten_alignment_result, $existing_alignment_result, 'using get returns same result as get_or_create');
};

subtest 'testing invalid merged alignment' => sub {
    my @segmented_params = (
        %parameters,
        instrument_data_segment => [$instrument_data[0]->id . ':test:read_group', $instrument_data[0]->id . ':test2:read_group'],
    );

    my $segmented_alignment_result = eval { $pkg->get_or_create(@segmented_params, users => $result_users) };
    my $error = $@;

    ok(!defined $segmented_alignment_result, 'no result returned for nonexistent segments');
    like($error, qr/Failed to find individual alignments for all instrument_data/, 'failed for expected reason');
};


subtest 'testing supersede merged alignment' => sub {
    #need recopy bam back to ar output dir since it was removed after merge
    for my $type ('', '.bai') {
        my $ar_base = 'all_sequences.bam'.$type;
        my $ar_file = File::Spec->join($tmp_shortcut_path1, 0, $ar_base);
        Genome::Sys->copy_file(File::Spec->join($expected_shortcut_path, 0, $ar_base), $ar_file);
        ok(-s $ar_file, "$ar_base copied over ok");
    }

    $parameters{instrument_data_id} = [$instrument_data[0]->id];
    my $small_merged_alignment_result = $pkg->create(%parameters, _user_data_for_nested_results => $result_users);
    isa_ok($small_merged_alignment_result, $pkg, 'produced small merged alignment result');

    # We need to override this because get_merged_alignment_results only returns objects in the database. 
    Sub::Install::install_sub({code => sub { my $self = shift; return @_; }, into => 'Genome::InstrumentData::AlignmentResult', as => 'filter_non_database_objects'});
    is_deeply([$small_merged_alignment_result->get_superseding_results(1)], [$merged_alignment_result], 'Got supersede merged alignment for small merged alignment');
};

done_testing();


sub generate_individual_alignment_results {
    my @instrument_data = @_;
    my @alignment_results;

    my %params = (
        subclass_name    => $alignment_result_class_name,
        module_version   => '12345',
        aligner_name     => $aligner_name,
        aligner_version  => $aligner_version,
        samtools_version => $samtools_version,
        picard_version   => $picard_version,
        reference_build  => $reference_build,
        instrument_data_segment_type => 'read_group',
        instrument_data_segment_id   => 'A:2',
        test_name        => $test_name,
    );

    for my $i (0, 1) {
        my $alignment_result = $alignment_result_class_name->__define__(
            %params,
            id                 => -8765432 + $i,
            output_dir         => $tmp_shortcut_path1 ."/$i",
            instrument_data_id => $instrument_data[$i]->id,
        );
        $alignment_result->recalculate_lookup_hash;

        isa_ok($alignment_result, 'Genome::InstrumentData::AlignmentResult');
        push @alignment_results, $alignment_result;
    }

    for my $i (0, 1) {
        my $alignment_result = $alignment_result_class_name->__define__(
            %params,
            id                 => -98765432 + $i,
            output_dir         => $tmp_shortcut_path2 ."/$i",
            instrument_data_id => $instrument_data[$i]->id,
            filter_name        => 'forward-only',
        );
        $alignment_result->recalculate_lookup_hash;
        isa_ok($alignment_result, 'Genome::InstrumentData::AlignmentResult');
        push @alignment_results, $alignment_result;
    }
    return @alignment_results;
}


sub generate_fake_instrument_data {
    my $fastq_directory = Genome::Config::get('test_inputs') . '/Genome-InstrumentData-Align-Maq/test_sample_name';
    my @instrument_data;

    for my $i (0, 2) {
        my $instrument_data = Genome::InstrumentData::Solexa->create(
            id                  => $FAKE_INSTRUMENT_DATA_ID + $i,
            sequencing_platform => 'solexa',
            flow_cell_id        => '12345',
            lane                => 4 + $i,
            median_insert_size  => '22',
            run_name            => 'test_run_name',
            subset_name         => 4 + $i,
            run_type            => 'Paired End Read 2',
            gerald_directory    => $fastq_directory,
            bam_path            => Genome::Config::get('test_inputs') .'/Genome-InstrumentData-AlignmentResult-Bwa/input.bam',
            library_id          => 2792100280,
        );
        isa_ok($instrument_data, 'Genome::InstrumentData::Solexa');
        push @instrument_data, $instrument_data;
    }
    return @instrument_data;
}
