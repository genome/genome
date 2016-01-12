#!/usr/bin/env genome-perl

use strict;
use warnings;

use Sub::Install;
use Test::Exception;
use Test::MockObject;
use Test::More tests => 3;

use above 'Genome';

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

    $annotation_build = Test::MockObject->new;
    ok($annotation_build, 'create mock annotation build');
    $annotation_build->set_always('rRNA_MT_file', '/dev/null');
    $annotation_build->set_always('annotation_file', '/dev/null');
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

subtest 'shortcut' => sub{
    plan tests => 4;

    my $cmd = $class->create(
        build_id => $build->id,
        build => $build,
    );

    # shortcut w/o annotation build
    $build->set_always('annotation_build', undef);
    lives_ok(sub{ $cmd->shortcut; }, 'shortcut w/o annotation build');
    like($cmd->debug_message, qr/Skipping PicardRnaSeqMetrics since annotation build is not defined/, 'correct message');

    # shortcut w/o required annotation files
    $build->set_always('annotation_build', $annotation_build);
    lives_ok(sub{ $cmd->shortcut; }, 'shortcut w/o required annotation files');
    like($cmd->debug_message, qr/Skipping PicardRnaSeqMetrics since annotation build is missing required files/, 'correct message');

};

done_testing();
