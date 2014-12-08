#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More tests => 7;
use Test::Exception;

use_ok('Genome::Model::Tools::Joinx::VcfAnnotate');

my $annotate = Genome::Model::Tools::Joinx::VcfAnnotate->create(
    input_file => 'some_file',
    annotation_file => 'another_file',
    info_fields => 'FIELD',
    use_version => 1.5,
);
ok($annotate,"Able to create object with version 1.5");
ok($annotate->check_minimum_version($Genome::Model::Tools::Joinx::VcfAnnotate::MINIMUM_JOINX_VERSION), "Version 1.5 ok");

ok($annotate->use_version("1.10"), "Changed version to 1.10");
ok($annotate->check_minimum_version($Genome::Model::Tools::Joinx::VcfAnnotate::MINIMUM_JOINX_VERSION), "Version 1.10 ok");

ok($annotate->use_version("1.2"), "Changed version to 1.2");
dies_ok( sub { $annotate->check_minimum_version($Genome::Model::Tools::Joinx::VcfAnnotate::MINIMUM_JOINX_VERSION) }, "Version 1.2 fails minimum version check");

