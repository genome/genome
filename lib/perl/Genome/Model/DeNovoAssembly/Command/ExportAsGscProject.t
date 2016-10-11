#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use File::Temp;
use Test::More tests => 2;

my $pkg = 'Genome::Model::DeNovoAssembly::Command::ExportAsGscProject';
use_ok($pkg) or die;

subtest 'execute' => sub{
    plan tests => 1;

    my $tempdir = File::Temp::tempdir(CLEANUP => 1);
    my $cmd = $pkg->execute(
        directory => $tempdir,
    );
    ok($cmd->result, 'execute command');

};

done_testing();
