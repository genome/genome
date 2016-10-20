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
use Test::MockObject;
use Test::More tests => 7;

my %setup;
subtest 'setup' => sub{
    plan tests => 3;

    $setup{pkg} = 'Genome::Model::DeNovoAssembly::Command::ExportAsGscProject';
    use_ok($setup{pkg}) or die;

    $setup{project} = Genome::Project->__define__(name => 'DE_NOVO_WO-1999');
    $setup{tempdir} = File::Temp::tempdir(CLEANUP => 1);

    my $denovo_newbler_pp = Genome::ProcessingProfile::DeNovoAssembly->__define__(name => 'DE_NOVO_PP', assembler_name => 'newbler de-novo-assemble');
    $setup{denovo_model} = Genome::Model::DeNovoAssembly->__define__(
        name => 'DENOVO_NEWBLER_MODEL',
        processing_profile => $denovo_newbler_pp,
        subject => Genome::Sample->create(name => 'HMPB-AAD13A05'),
    );
    ok($setup{denovo_model}, 'define newbler model');

    $setup{refalign_model} = Genome::Model::ReferenceAlignment->__define__(name => 'REFALIGN_MODEL');
    ok($setup{refalign_model}, 'define refalign model');

};

subtest 'supported assemblers' => sub{
    plan tests => 9;

    ok($setup{pkg}->supported_assemblers, 'supported_assemblers');

    my $assembler = 'newbler de-novo-assemble';
    ok($setup{pkg}->subdirs_for_assembler($assembler), 'subdirs_for_assembler');
    throws_ok(sub{ $setup{pkg}->subdirs_for_assembler; }, qr/2 were expected/, 'subdirs_for_assembler w/o assembler');

    ok($setup{pkg}->assemblers_edit_dir($assembler), 'assemblers_edit_dir');
    throws_ok(sub{ $setup{pkg}->assemblers_edit_dir; }, qr/2 were expected/, 'assemblers_edit_dir w/o assembler');

    my $model = Test::MockObject->new();
    $model->set_always('__display_name__', 'MOCK');
    $model->set_isa('Genome::Model::DeNovoAssembly');
    $model->set_always('processing_profile', $setup{denovo_model}->processing_profile);
    $model->set_list('instrument_data', 'INSTDATA1');
    ok($setup{pkg}->is_model_supported($model), "model ".$setup{denovo_model}->__display_name__." is supported");

    my $unsupported_denovo_pp = Test::MockObject->new();
    $unsupported_denovo_pp->set_always('assembler_name', 'unsupported');
    $model->set_always('processing_profile', $unsupported_denovo_pp);
    ok(!$setup{pkg}->is_model_supported($model), "model with unsupported assembler is NOT supported");

    $model->set_always('processing_profile', $setup{denovo_model}->processing_profile);
    $model->set_list('instrument_data', 'INSTDATA1', 'INSTDATA2');
    ok(!$setup{pkg}->is_model_supported($model), "model with more than one inst data is NOT supported");

    throws_ok(sub{ $setup{pkg}->is_model_supported}, qr/2 were expected/, 'is_model_supported w/o assembler');

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
    plan tests => 1;

    for my $m ( $setup{denovo_model}, $setup{refalign_model} ) {
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
    plan tests => 9;

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

    my $fastq_file = File::Spec->join($data_directory, 'input.fastq');
    Genome::Sys->write_file($fastq_file, 'FASTQ');

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
    my $expected_fastq_file = File::Spec->join($expected_project_directory, 'edit_dir', 'input.fastq');
    ok(-l $expected_phdball_file, 'linked fastq file');
    #diag("$data_directory\n$setup{tempdir}"); <STDIN>;

};

done_testing();
