#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::Exception;
use Test::MockObject;
use Test::More tests => 2;

use above 'Genome';

my $class = 'Genome::InstrumentData::AlignmentResult::Command::PicardRnaSeqMetrics';
use_ok($class) or die;

my $annotation_build = Test::MockObject->new;
$annotation_build->set_always('rRNA_MT_file', 'foo');
$annotation_build->set_always('annotation_file', 'foo');

my $reference_build =  Test::MockObject->new;
$reference_build->set_always('id', '1');

subtest 'does_annotation_build_have_required_files' => sub{
    plan tests => 5;

    throws_ok(
        sub{ $class->does_annotation_build_have_required_files(); },
        qr/but 3 were expected/,
        'does_annotation_build_have_required_files fails w/o params',
    );
    throws_ok(
        sub{ $class->does_annotation_build_have_required_files($annotation_build, $reference_build); },
        qr/Cannot proceed\! Missing required files from annotation build: rRNA_MT_file annotation_file/,
        'does_annotation_build_have_required_files fails when annotation build does not have required files',
    );

    for my $file_method ( $class->required_files_from_annotation_build ) {
        $annotation_build->set_always($file_method, $0);
        throws_ok(
            sub{ $class->does_annotation_build_have_required_files($annotation_build, $reference_build); },
            qr/Cannot proceed\! Missing required files from annotation build: /,
            "does_annotation_build_have_required_files fails when annotation build does not have $file_method",
        );
        $annotation_build->set_always($file_method, 'foo');
    }

    map { $annotation_build->set_always($_, $0) } $class->required_files_from_annotation_build;
    $annotation_build->set_always('annotation_file', $0);
    lives_ok(
        sub{$class->does_annotation_build_have_required_files($annotation_build, $reference_build);},
        'does_annotation_build_have_required_files succeeds with annotation files',
    );

};

done_testing();
