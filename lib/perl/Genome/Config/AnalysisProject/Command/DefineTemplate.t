#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

use above 'Genome';

use Test::More tests => 4;

my $pkg = 'Genome::Config::AnalysisProject::Command::DefineTemplate';
use_ok($pkg);

my $cmd = $pkg->create(
    name => 'test of ' . $pkg,
    environment => 'ad-hoc',
    no_config => 1,
);
isa_ok($cmd, $pkg, 'created command');
my $template = $cmd->execute;
isa_ok($template, 'Genome::Config::AnalysisProject');
is($template->status, 'Template', 'template is a template');
