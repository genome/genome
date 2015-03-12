#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

use above "Genome";

use Test::More tests => 23;

use File::Copy qw(copy);
use Genome::Utility::Test qw(compare_ok);

my $class = 'Genome::FeatureList::Command::Import';
use_ok($class);

my $reference = Genome::Model::Build::ReferenceSequence->get(name => 'NCBI-human-build36');

my @test_beds = (undef); #fill 0th position
push @test_beds, map { sprintf('%s.d/%s.bed', __FILE__, $_) } 1..6;


#1 has unknown track name
#3 has mismatched chromosomes
my @should_fail = (1,3);

#2 is normal
#4 has browser lines and alternate track names
my @should_pass = (2,4,5,6);

my @new_feature_lists;
for my $passing_bed (@test_beds[@should_pass]) {
    my $cmd = Genome::FeatureList::Command::Import->create(
        file_path => $passing_bed,
        name => 'test ' . substr($passing_bed, -5),
        reference => $reference,
    );
    isa_ok($cmd, 'Genome::FeatureList::Command::Import', 'created command');
    ok($cmd->execute, 'executed command');
    my $new_feature_list = $cmd->new_feature_list;
    isa_ok($new_feature_list, 'Genome::FeatureList', 'created new feature-list');
    push @new_feature_lists, $new_feature_list;
}
is($new_feature_lists[0]->file_content_hash, $new_feature_lists[1]->file_content_hash, 'sanitization of browser lines led to same BED file');
is($new_feature_lists[0]->file_content_hash, $new_feature_lists[2]->file_content_hash, 'sanitization of CRs led to same BED file');
is($new_feature_lists[0]->file_content_hash, $new_feature_lists[3]->file_content_hash, 'selected correct BED file out of directory');

for my $failing_bed(@test_beds[@should_fail]) {
    my $cmd = Genome::FeatureList::Command::Import->create(
        file_path => $failing_bed,
        name => 'test ' . substr($failing_bed, -5),
        reference => $reference,
    );
    isa_ok($cmd, 'Genome::FeatureList::Command::Import', 'created command');
    my $rv = eval { $cmd->execute; };
    my $error = $@;
    ok(!$rv && $error, 'failed to execute');
}

subtest _nimblegen_capture_primary_pair => sub {
    plan tests => 2;
    my @bed_files = qw(
        150203_HG19_CRC_OID42357_EZ_HX1_capture_targets.bed
        150203_HG19_CRC_OID42357_EZ_HX1_coverage_summary.txt
        150203_HG19_CRC_OID42357_EZ_HX1_primary_targets.bed
        CRC_nimbegen_design_hg19_v3.bed
    );
    my @pair = $class->_nimblegen_capture_primary_pair(@bed_files);
    is(scalar(@pair), 2, 'two items are returned');
    isnt($pair[0], $pair[1], 'the items are distinct');
};

subtest _has_nimblegen_capture_primary_pair => sub {
    plan tests => 3;

    do {
        my @bed_files = qw(
            150203_HG19_CRC_OID42357_EZ_HX1_capture_targets.bed
            150203_HG19_CRC_OID42357_EZ_HX1_coverage_summary.txt
            150203_HG19_CRC_OID42357_EZ_HX1_primary_targets.bed
            CRC_nimbegen_design_hg19_v3.bed
        );
        ok($class->_has_nimblegen_capture_primary_pair(@bed_files),
            'single match: should have Nimblegen capture/primary pair');
    };

    do {
        my @bed_files = qw(
            150203_HG19_CRC_OID42357_EZ_HX1_coverage_summary.txt
            150203_HG19_CRC_OID42357_EZ_HX1_primary_targets.bed
            CRC_nimbegen_design_hg19_v3.bed
        );
        ok(not($class->_has_nimblegen_capture_primary_pair(@bed_files)),
            'empty match: should not have Nimblegen capture/primary pair');
    };

    do {
        my @bed_files = qw(
            150203_HG19_CRC_OID42357_EZ_HX1_capture_targets.bed
            150203_HG19_CRC_OID42357_EZ_HX1_coverage_summary.txt
            150203_HG19_CRC_OID42357_EZ_HX1_primary_targets.bed
            250203_HG19_CRC_OID42357_EZ_HX1_capture_targets.bed
            250203_HG19_CRC_OID42357_EZ_HX1_primary_targets.bed
            CRC_nimbegen_design_hg19_v3.bed
        );
        ok(not($class->_has_nimblegen_capture_primary_pair(@bed_files)),
            'multiple match: should not have Nimblegen capture/primary pair');
    };
};

subtest 'Nimblegen capture/primary Import' => sub {
    plan tests => 2;

    my $reference = Genome::Model::Build::ReferenceSequence->get(name => 'NCBI-human-build36');
    my ($expected_output, @originals) = map { sprintf('%s.d/%s.bed', __FILE__, $_) }
        qw(X_multitrack X_capture_targets X_primary_targets);

    my $temp_dir = Genome::Sys->create_temp_directory();
    for (@originals) {
        copy($_, $temp_dir);
    }

    my $cmd = $class->create(
        file_path => $temp_dir,
        name => 'Nimblegen capture/primary Import Test',
        reference => $reference,
    );
    ok($cmd->execute, 'import executed successfully');

    my $fl = $cmd->new_feature_list;
    compare_ok($fl->file_path, $expected_output);
};
