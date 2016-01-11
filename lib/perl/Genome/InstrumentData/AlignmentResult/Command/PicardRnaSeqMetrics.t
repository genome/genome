#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::Exception;
use Test::MockObject;
use Test::More tests => 3;

use above 'Genome';

my $class = 'Genome::InstrumentData::AlignmentResult::Command::PicardRnaSeqMetrics';
use_ok($class) or die;

my $existing_file = Genome::Sys->write_file(Genome::Sys->create_temp_file_path('existing'), 1);
my $existing_file_sub = sub{ $existing_file; };
my $existing_file_no_ref_build = Genome::Sys->write_file(Genome::Sys->create_temp_file_path('no-ref-build'), 1);
my $existing_file_with_no_ref_build_sub = sub{ return if defined $_[2]; $existing_file_no_ref_build; };
my $non_existing_file_sub = sub{ "/dev/null"; };

my $annotation_build = Test::MockObject->new;
$annotation_build->mock('rRNA_MT_file', $non_existing_file_sub);
$annotation_build->mock('annotation_file', $non_existing_file_sub);

my $reference_build =  Test::MockObject->new;
$reference_build->set_always('id', '1');

subtest 'file_from_annotation_build' => sub{
    plan tests => 3;

    throws_ok(
        sub{ $class->file_from_annotation_build(); },
        qr/but 4 were expected/,
        'does_annotation_build_have_required_files fails w/o params',
    );

    $annotation_build->mock('rRNA_MT_file', $existing_file_sub);
    is(
        $class->file_from_annotation_build($annotation_build, $reference_build, 'rRNA_MT_file'),
        $existing_file,
        'file_from_annotation_build rRNA_MT_file',
    );

    $annotation_build->mock('rRNA_MT_file', $existing_file_with_no_ref_build_sub);
    is(
        $class->file_from_annotation_build($annotation_build, $reference_build, 'rRNA_MT_file'),
        $existing_file_no_ref_build,
        'file_from_annotation_build rRNA_MT_file',
    );

};

subtest 'does_annotation_build_have_required_files' => sub{
    plan tests => 5;

    throws_ok(
        sub{ $class->does_annotation_build_have_required_files(); },
        qr/but 3 were expected/,
        'does_annotation_build_have_required_files fails w/o params',
    );

    $annotation_build->mock('rRNA_MT_file', $non_existing_file_sub);
    $annotation_build->mock('annotation_file', $non_existing_file_sub);
    throws_ok(
        sub{ $class->does_annotation_build_have_required_files($annotation_build, $reference_build); },
        qr/Cannot proceed\! Missing required files from annotation build: rRNA_MT_file annotation_file/,
        'does_annotation_build_have_required_files fails when annotation build does not have required files',
    );

    for my $file_method ( $class->required_files_from_annotation_build ) {
        $annotation_build->mock($file_method, $existing_file_sub);
        throws_ok(
            sub{ $class->does_annotation_build_have_required_files($annotation_build, $reference_build); },
            qr/Cannot proceed\! Missing required files from annotation build: /,
            "does_annotation_build_have_required_files fails when annotation build does not have $file_method",
        );
        $annotation_build->mock($file_method, $non_existing_file_sub);
    }

    $annotation_build->mock('rRNA_MT_file', $existing_file_sub);
    $annotation_build->mock('annotation_file', $existing_file_sub);
    lives_ok(
        sub{$class->does_annotation_build_have_required_files($annotation_build, $reference_build);},
        'does_annotation_build_have_required_files succeeds with annotation files',
    );

};

done_testing();
