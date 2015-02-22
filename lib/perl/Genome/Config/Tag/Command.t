#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Test::More tests => 5;

my $cmd_base = 'Genome::Config::Tag::Command';
use_ok($cmd_base);

for my $type ('List', 'Update::Description') {
    my $class = join('::', $cmd_base, $type);

    my $cmd = $class->create();
    isa_ok($cmd, 'Command', 'command can be created');

    ok($class->help_usage_complete_text, 'command produces help text');
}
