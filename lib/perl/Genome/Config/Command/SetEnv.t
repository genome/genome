#!/usr/bin/env genome-perl

use strict;
use warnings;

use Genome qw();
use Genome::Test::Config qw(setup_config);

use Test::More tests => 2;

use_ok('Genome::Config::Command::SetEnv');

my ($temp_dirs, $new_temp_dir) = Genome::Test::Config::temp_dir_helper();
local $ENV{XGENOME_CONFIG_SNAP} = $new_temp_dir->();
local $ENV{XGENOME_CONFIG_HOME} = $new_temp_dir->();
local $ENV{XGENOME_CONFIG_DIRS} = $new_temp_dir->();
setup_config(
    spec => {
        home_key => {
            env => 'XGENOME_HOME_KEY',
        },
    },
    global => {
        home_key => 'orig_value',
    },
);

my @commands = (
    Genome::Config::Command::SetEnv->create(key => 'home_key', value => 'foo'),
);
for my $command (@commands) {
    ok($command->execute);
}

