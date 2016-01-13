#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::Exception;
use Test::MockObject;
use Test::More tests => 3;

use above 'Genome';
use Genome::Test::Factory::Model::ImportedAnnotation;

my $class = 'Genome::InstrumentData::AlignmentResult::Command::PicardRnaSeqMetrics';
use_ok($class) or die;

my $annotation_build = Genome::Test::Factory::Model::ImportedAnnotation->create_mock_build;
my $reference_build =  Test::MockObject->new;
$reference_build->set_always('id', '1');

subtest 'file_from_annotation_build' => sub{
    plan tests => 3;

    throws_ok(
        sub{ $class->file_from_annotation_build(); },
        qr/but 4 were expected/,
        'verify_annotation_build_has_required_files fails w/o params',
    );

    $annotation_build->rRNA_MT_file_exists;
    is(
        $class->file_from_annotation_build($annotation_build, $reference_build, 'rRNA_MT_file'),
        $annotation_build->rRNA_MT_file,
        'file_from_annotation_build rRNA_MT_file',
    );

    $annotation_build->rRNA_MT_file_exists_without_reference_build;
    is(
        $class->file_from_annotation_build($annotation_build, $reference_build, 'rRNA_MT_file'),
        $annotation_build->rRNA_MT_file,
        'file_from_annotation_build rRNA_MT_file',
    );

};

subtest 'verify_annotation_build_has_required_files' => sub{
    plan tests => 5;

    throws_ok(
        sub{ $class->verify_annotation_build_has_required_files(); },
        qr/but 3 were expected/,
        'verify_annotation_build_has_required_files fails w/o params',
    );

    $annotation_build->rRNA_MT_file_does_not_exist;
    $annotation_build->annotation_file_does_not_exist;
    throws_ok(
        sub{ $class->verify_annotation_build_has_required_files($annotation_build, $reference_build); },
        qr/Cannot proceed\! Missing required files from annotation build: rRNA_MT_file annotation_file/,
        'verify_annotation_build_has_required_files fails when annotation build does not have required files',
    );

    for my $file_method ( $class->required_files_from_annotation_build ) {
        my $exists_method = $file_method.'_exists';
        $annotation_build->$exists_method;
        throws_ok(
            sub{ $class->verify_annotation_build_has_required_files($annotation_build, $reference_build); },
            qr/Cannot proceed\! Missing required files from annotation build: /,
            "verify_annotation_build_has_required_files fails when annotation build does not have $file_method",
        );
        my $does_not_exist_method = $file_method.'_does_not_exist';
        $annotation_build->$does_not_exist_method;
    }

    $annotation_build->rRNA_MT_file_exists;
    $annotation_build->annotation_file_exists;
    lives_ok(
        sub{$class->verify_annotation_build_has_required_files($annotation_build, $reference_build);},
        'verify_annotation_build_has_required_files succeeds with annotation files',
    );

};

done_testing();
