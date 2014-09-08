#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

use above 'Genome';
use Genome::Test::Factory::Config::Profile::Item;

use Test::More tests => 10;

use_ok('Genome::Config::AnalysisProject::Command::UntagConfigFile') or die 'Command cannot be used.';

my $profile_item = Genome::Test::Factory::Config::Profile::Item->setup_object();
isa_ok($profile_item, 'Genome::Config::Profile::Item', 'created profile item');

my @tags = map {
    Genome::Config::Tag->create(
        name => 'Genome::Config::AnalysisProject::Command::UntagConfigFile test ' . $_,
    )
} (1..3);

for my $tag (@tags) {
    my $bridge = Genome::Config::Tag::Profile::Item->create(
        tag => $tag,
        profile_item => $profile_item
    );
    ok($bridge, 'tagged profile item for test');
}
is(scalar(@{[$profile_item->tags]}), 3, 'all tags assigned to profile item');

my $cmd = Genome::Config::AnalysisProject::Command::UntagConfigFile->create(
    tag => $tags[1],
    profile_items => [$profile_item],
);
isa_ok($cmd, 'Genome::Config::AnalysisProject::Command::UntagConfigFile', 'created command');
ok($cmd->execute, 'executed command');

my (@assigned_tags) = $profile_item->tags;
is(scalar(@assigned_tags), 2, 'one tag removed from assignment');
is_deeply(
    [sort map $_->id, @assigned_tags],
    [sort map $_->id, @tags[0,2]],
    'remaining tags are the expected tags'
);
