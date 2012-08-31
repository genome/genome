#!/usr/bin/env genome-perl
use strict;
use warnings;
use above "Genome";
use Test::MockObject;
use Test::More;

my $i = Test::MockObject->new();
$i->set_isa('Genome::InstrumentData::AlignmentResult', 'UR::Object');
$i->set_always('__display_name__', "Mock Alignment Result 1");
$i->set_always('delete', 1);
ok($i, "created mock alignment result 1");

my $i2 = Test::MockObject->new();
$i2->set_isa('Genome::InstrumentData::AlignmentResult', 'UR::Object');
$i2->set_always('__display_name__', "Mock Alignment Result 2");
$i2->set_always('delete', 1);
ok($i2, "created mock alignment result 2");

my $remove_command = Genome::InstrumentData::Command::AlignmentResult::Remove->create(alignment_results => [$i, $i2]);
ok($remove_command->execute(), "command succeeded");

done_testing();
