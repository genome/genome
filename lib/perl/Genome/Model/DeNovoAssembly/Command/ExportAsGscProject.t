#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use File::Temp;
use Genome::Test::Factory::Model::ReferenceAlignment;
use Test::Exception;
use Test::More tests => 4;

my %setup;
subtest 'setup' => sub{
    plan tests => 1;

    $setup{pkg} = 'Genome::Model::DeNovoAssembly::Command::ExportAsGscProject';
    use_ok($setup{pkg}) or die;

    $setup{tempdir} = File::Temp::tempdir(CLEANUP => 1);
    $setup{project} = Genome::Project->__define__(name => 'DE_NOVO_WO-1999');
    $setup{denovo_pp} = Genome::ProcessingProfile::DeNovoAssembly->__define__(name => 'DE_NOVO_PP');
    $setup{denovo_model} = Genome::Model::DeNovoAssembly->__define__(
        name => 'DE_NOVO_MODEL',
        processing_profile => $setup{denovo_pp},
    );
    $setup{refalign_model} = Genome::Test::Factory::Model::ReferenceAlignment->generate_obj;

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
        entity_class_name => $setup{denovo_model}->__meta__->class_name,
    );

    throws_ok(
        sub{
            $setup{pkg}->execute(
                project => $setup{project},
                directory => $setup{tempdir},
            );
        },
        qr/Models associated with .+ do not exist!/,
        'Fails w/o project parts'
    );

    $pp->delete;

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
