#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
    $ENV{UR_COMMAND_DUMP_STATUS_MESSAGES} = 1;
};

use Test::More tests => 5;
use above "Genome";

use Genome::Test::Factory::Sample;

my $class = 'Genome::Site::TGI::Synchronize::CheckForNameUpdates';

use_ok($class);

my @samples = map { Genome::Test::Factory::Sample->setup_object } 1..7;
my @organism_samples;
my @misc_updates;

for my $i (0..$#samples) {
    my $organism_sample = Genome::Site::TGI::Synchronize::Classes::OrganismSample->create(
        id => $samples[$i]->id,
        name => $samples[$i]->name . ($i % 2? 'xx' : ''),
    );

    my $misc_update = Genome::Site::TGI::Synchronize::Classes::MiscUpdate->create(
        description => 'UPDATE',
        subject_class_name => 'gsc.organism_sample',
        subject_property_name => 'full_name',
        old_value => $samples[$i]->name . ($i % 2? '' : 'xx'),
        new_value => $organism_sample->name,
        subject_id => $organism_sample->id,
        edit_date => '2010-01-01 01:01:01',
        is_reconciled => 0,
    );

    push @organism_samples, $organism_sample;
    push @misc_updates, $misc_update;
}

my $cmd = $class->create(
    updates => \@misc_updates,
);
isa_ok($cmd, $class, 'created command');

ok($cmd->execute(), 'executed command');

my $sample_count = scalar(@samples);
my $reconciled_count = grep { $_->is_reconciled } @misc_updates;
my $failed_count = grep { $_->has_failed } @misc_updates;

is($failed_count, int($sample_count / 2), 'found expected number of updates needing attention');
is($reconciled_count, $sample_count - $failed_count, 'found expected number of processed updates');

