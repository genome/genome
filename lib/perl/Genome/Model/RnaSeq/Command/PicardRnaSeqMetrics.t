#!/usr/bin/env genome-perl

use strict;
use warnings;

use Sub::Install;
use Test::Exception;
use Test::MockObject;
use Test::More tests => 4;

use above 'Genome';
use Genome::Test::Factory::Model::ImportedAnnotation;

my $class = 'Genome::Model::RnaSeq::Command::PicardRnaSeqMetrics';
use_ok($class) or die;

my ($build, $reference_build, $annotation_build);
subtest "setup" => sub{
    plan tests => 4;
    
    $build = Test::MockObject->new;
    $build->set_always('id', 1);
    $build->isa('Genome::Model::Build');
    $build->set_always('build_id', 1);
    $build->set_always('annotation_build', undef);
    ok($build, 'create mock build');

    $reference_build =  Test::MockObject->new;
    ok($reference_build, 'create mock reference build');
    $reference_build->set_always('id', '1');
    $build->set_always('reference_build', $reference_build);

    $annotation_build = Genome::Test::Factory::Model::ImportedAnnotation->create_mock_build;
    ok($annotation_build, 'create mock annotation build');
    $build->set_always('annotation_build', $annotation_build);

    my $alignment_result = Test::MockObject->new;
    ok($build, 'create mock alignment result');
    $build->set_always('alignment_result', $alignment_result);

    Sub::Install::reinstall_sub({
            into => $class,
            as => 'build',
            code => sub{ $build },
        });
};

subtest 'should_skip' => sub{
    plan tests => 5;

    my $cmd = $class->create(
        build_id => $build->id,
        build => $build,
    );

    # should skip w/o annotation build
    $build->set_always('annotation_build', undef);
    ok($cmd->should_skip, 'should_skip w/o annotation build');
    like($cmd->debug_message, qr/Skipping PicardRnaSeqMetrics since annotation build is not defined/, 'correct message');

    # should skip w/o required annotation files
    $annotation_build->rRNA_MT_file_does_not_exist;
    $annotation_build->annotation_file_does_not_exist;
    $build->set_always('annotation_build', $annotation_build);
    ok($cmd->should_skip, 'should_skip w/o required annotation files');
    like($cmd->debug_message, qr/Skipping PicardRnaSeqMetrics since annotation build is missing required files/, 'correct message');

    # should not skip
    $annotation_build->rRNA_MT_file_exists;
    $annotation_build->annotation_file_exists;
    ok(!$cmd->should_skip, 'should not skip w/ annotation build and required files');

};

subtest 'shortcut and execute should call should_skip' => sub{
    plan tests => 1;

    my $should_skip_sub = $class->can('shortcut');
    my $should_skip_cnt = 0;
    Sub::Install::reinstall_sub({
            into => $class,
            as => 'should_skip',
            code => sub{ $should_skip_cnt++; return 1; },
        });

    my $cmd = $class->create(
        build_id => $build->id,
        build => $build,
    );
    $cmd->shortcut;
    $cmd->execute;
    is($should_skip_cnt, 2, 'should_skip called twice');

    Sub::Install::reinstall_sub({
            into => $class,
            as => 'should_skip',
            code => $should_skip_sub,
        });
};

done_testing();
