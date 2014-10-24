#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

use above 'Genome';
use Genome::Test::Factory::Config::Profile::Item;

use Test::More tests => 6;

use_ok('Genome::Config::AnalysisProject::Command::TagConfigFile') or die 'Command cannot be used.';

my $profile_item = Genome::Test::Factory::Config::Profile::Item->setup_object();
isa_ok($profile_item, 'Genome::Config::Profile::Item', 'created profile item');

my $tag = Genome::Config::Tag->create(
    name => 'Genome::Config::AnalysisProject::Command::TagConfigFile test 1',
);
isa_ok($tag, 'Genome::Config::Tag', 'created tag');

my $cmd = Genome::Config::AnalysisProject::Command::TagConfigFile->create(
    tag => $tag,
    profile_items => [$profile_item],
);
isa_ok($cmd, 'Genome::Config::AnalysisProject::Command::TagConfigFile', 'created command');
ok($cmd->execute, 'executed command');

my $assigned_tag = $profile_item->tags;
is($assigned_tag,$tag, 'tag was added to profile item');
