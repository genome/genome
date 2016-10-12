#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use File::Temp;
use Genome::Test::Factory::Model::ReferenceAlignment;
use Test::Exception;
use Test::More tests => 5;

my %setup;
subtest 'setup' => sub{
    plan tests => 1;

    $setup{pkg} = 'Genome::Model::DeNovoAssembly::Command::ExportAsGscProject';
    use_ok($setup{pkg}) or die;

    $setup{tempdir} = File::Temp::tempdir(CLEANUP => 1);
    $setup{project} = Genome::Project->__define__(name => 'DE_NOVO_WO-1999');

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

    my $denovo_newbler_pp = Genome::ProcessingProfile::DeNovoAssembly->__define__(name => 'DE_NOVO_PP', assembler_name => 'newbler');
    my $denovo_newbler_model = Genome::Model::DeNovoAssembly->__define__(
        name => 'DENOVO_NEWBLER_MODEL',
        processing_profile => $denovo_newbler_pp,
    );
    ok($denovo_newbler_model, 'define newbler model');

    #my $refalign_model = Genome::Test::Factory::Model::ReferenceAlignment->generate_obj;
    my $refalign_model = Genome::Model::ReferenceAlignment->__define__(name => 'REFALIGN_MODEL');
    ok($refalign_model, 'define refalign model');

    for my $m ( $denovo_newbler_model, $refalign_model ) {
        $setup{project}->add_part(
            entity_id => $m->id,
            entity_class_name => $m->__meta__->class_name,
        );
    }

    my @parts = $setup{project}->parts;
    is(@parts, 2, 'add models to proejct');

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
    plan tests => 1;

    my $cmd = $setup{pkg}->execute(
        project => $setup{project},
        directory => $setup{tempdir},
    );
    ok($cmd->result, 'execute command');

};

done_testing();
