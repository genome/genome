#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use File::Temp;
use Test::More tests => 2;

my %setup;
subtest 'setup' => sub{
    plan tests => 1;

    $setup{pkg} = 'Genome::Model::DeNovoAssembly::Command::ExportAsGscProject';
    use_ok($setup{pkg}) or die;

    $setup{tempdir} = File::Temp::tempdir(CLEANUP => 1);
    $setup{project} = Genome::Project->__define__(name => 'DE_NOVO_WO-1999');

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
