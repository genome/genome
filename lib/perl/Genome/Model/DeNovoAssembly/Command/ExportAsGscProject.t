#!/usr/bin/env genome-perl

BEGIN {
    $ENV{URI_DBI_NO_COMMIT} = 1;
}

use strict;
use warnings;

use above 'Genome';

use File::Temp;
use Genome::Test::Factory::Model::ReferenceAlignment;
use Test::Exception;
use Test::More tests => 7;

my %setup;
subtest 'setup' => sub{
    plan tests => 1;

    $setup{pkg} = 'Genome::Model::DeNovoAssembly::Command::ExportAsGscProject';
    use_ok($setup{pkg}) or die;

    $setup{project} = Genome::Project->__define__(name => 'DE_NOVO_WO-1999');
    $setup{tempdir} = File::Temp::tempdir(CLEANUP => 1);

};

subtest 'supported asemblers' => sub{
    plan tests => 5;

    ok($setup{pkg}->supported_assemblers, 'supported_assemblers');
    ok($setup{pkg}->subdirs_for_assembler('newbler de-novo-assemble'), 'subdirs_for_assembler');
    throws_ok(sub{ $setup{pkg}->subdirs_for_assembler}, qr/2 were expected/, 'subdirs_for_assembler w/o assembler');
    ok($setup{pkg}->assemblers_edit_dir('newbler de-novo-assemble'), 'subdirs_for_assembler');
    throws_ok(sub{ $setup{pkg}->subdirs_for_assembler}, qr/2 were expected/, 'subdirs_for_assembler w/o assembler');

};

subtest 'fails without models associated with project' => sub{
    plan tests => 1;

    throws_ok(
        sub{
            $setup{pkg}->execute(
                project => $setup{project},
                directory => $setup{tempdir},
            );
        },
        qr/No de novo models associated/,
        'Fails w/o project parts'
    );

};

subtest 'fails when de novo models do not exist' => sub{
    plan tests => 1;

    my $pp = Genome::ProjectPart->__define__(
        project => $setup{project},
        entity_id => 'BLAH',
        entity_class_name => 'Genome::Model::DeNovoAssembly',
    );

    throws_ok(
        sub{
            $setup{pkg}->execute(
                project => $setup{project},
                directory => $setup{tempdir},
            );
        },
        qr/Models associated with .+ do not exist/,
        'Fails w/o project parts'
    );

    $pp->delete;

};

subtest 'add models to project' => sub{
    plan tests => 3;

    my $denovo_newbler_pp = Genome::ProcessingProfile::DeNovoAssembly->__define__(name => 'DE_NOVO_PP', assembler_name => 'newbler de-novo-assemble');
    $setup{denovo_model} = Genome::Model::DeNovoAssembly->__define__(
        name => 'DENOVO_NEWBLER_MODEL',
        processing_profile => $denovo_newbler_pp,
        subject => Genome::Sample->create(name => 'HMPB-AAD13A05'),
    );
    ok($setup{denovo_model}, 'define newbler model');

    my $refalign_model = Genome::Model::ReferenceAlignment->__define__(name => 'REFALIGN_MODEL');
    ok($refalign_model, 'define refalign model');

    for my $m ( $setup{denovo_model}, $refalign_model ) {
        $setup{project}->add_part(
            entity_id => $m->id,
            entity_class_name => $m->__meta__->class_name,
        );
    }

    my @parts = $setup{project}->parts;
    is(@parts, 2, 'add models to project');

};

subtest 'fails when no succeeded builds' => sub{
    plan tests => 1;

    throws_ok(
        sub{
            $setup{pkg}->execute(
                project => $setup{project},
                directory => $setup{tempdir},
            );
        },
        qr/No succeeded builds found/,
        'Fails w/o project parts'
    );

};

subtest 'execute' => sub{
    plan tests => 8;

    my $data_directory = File::Temp::tempdir(CLEANUP => 1);
    my $build = Genome::Model::Build::DeNovoAssembly->__define__(
        model => $setup{denovo_model},
        data_directory => $data_directory,
        status => 'Succeeded',
    );
    my $phdball_dir = File::Spec->join($data_directory, 'consed', 'phdball_dir');
    Genome::Sys->create_directory($phdball_dir);
    my $phdball_ball_file = File::Spec->join($phdball_dir, 'phd.ball.1');
    Genome::Sys->write_file($phdball_ball_file, 'PHD');

    my $cmd = $setup{pkg}->execute(
        project => $setup{project},
        directory => $setup{tempdir},
    );
    ok($cmd->result, 'execute command');

    my $expected_project_directory = File::Spec->join($setup{tempdir}, $setup{denovo_model}->subject->name);
    ok(-d $expected_project_directory, 'created project dir');
    my $build_id_file = File::Spec->join($expected_project_directory, 'build_id='.$setup{denovo_model}->last_succeeded_build->id);
    ok(-e $build_id_file, 'created build id file');
    for my $subdir_name ( $cmd->subdirs_for_assembler( $setup{denovo_model}->processing_profile->assembler_name ) ) {
        my $subdir = File::Spec->join($expected_project_directory, $subdir_name);
        ok(-d $subdir, "created $subdir_name");
    }
    my $expected_phdball_file = File::Spec->join($expected_project_directory, 'phdball_dir', 'phd.ball.1');
    ok(-l $expected_phdball_file, 'linked phdball file');
    #diag("$data_directory\n$setup{tempdir}"); <STDIN>;

};

done_testing();
