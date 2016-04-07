#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Genome::Test::Config qw(setup_config);

use Test::More tests => 4;

use_ok('Genome::Config::Command::List');

my ($temp_dirs, $new_temp_dir) = Genome::Test::Config::temp_dir_helper();
local $ENV{XGENOME_CONFIG_SNAP} = $new_temp_dir->();
local $ENV{XGENOME_CONFIG_HOME} = $new_temp_dir->();
local $ENV{XGENOME_CONFIG_DIRS} = $new_temp_dir->();
setup_config(
    spec => {
        home_key => {
            type => 'Str',
        },
        conf_key => {
            type => 'Str',
        },
    },
    global => {
        conf_key => 'conf_dir_value',
        home_key => 'conf_dir_value',
    },
);

my @commands = (
    Genome::Config::Command::List->create(),
    Genome::Config::Command::List->create(keys => [qw(home_key)]),
    Genome::Config::Command::List->create(keys => [qw(home_key conf_key)]),
);
for my $command (@commands) {
    ok($command->execute);
}

