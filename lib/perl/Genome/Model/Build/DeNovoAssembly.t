#!/usr/bin/env genome-perl
#
# TESTS ARE IN SUBCLASSES
#

use strict;
use warnings;

use above 'Genome';

use Test::More;

use_ok('Genome::Model::Build::DeNovoAssembly') or die;
use_ok('Genome::InstrumentData::InstrumentDataTestObjGenerator') or die;

my $id = Genome::InstrumentData::InstrumentDataTestObjGenerator::create_solexa_instrument_data("fake.bam");
my $id2 = Genome::InstrumentData::InstrumentDataTestObjGenerator::create_solexa_instrument_data("fake2.bam");
$id->read_orientation("reverse_forward");
$id2->read_orientation("reverse_forward");

my $read_orientation = Genome::Model::Build::DeNovoAssembly->resolve_attribute_for_instrument_data("read_orientation", 1, $id, $id2);
is($read_orientation, "reverse_forward", "Found read orientation when it was the same");

$read_orientation = Genome::Model::Build::DeNovoAssembly->resolve_attribute_for_instrument_data("read_orientation", 1, $id);
is($read_orientation, "reverse_forward", "Works for a single instrument data");

$id2->read_orientation("forward_reverse");
$read_orientation = Genome::Model::Build::DeNovoAssembly->resolve_attribute_for_instrument_data("read_orientation", 1, $id, $id2);
ok(!$read_orientation, "Didn't get read orientation when it was different");

my @read_orientations = Genome::Model::Build::DeNovoAssembly->resolve_attribute_for_instrument_data("read_orientation", 0, $id, $id2);
is_deeply(\@read_orientations, ["reverse_forward","forward_reverse"],"Got combined read orientation when it was different and not required to be the same");

@read_orientations = Genome::Model::Build::DeNovoAssembly->resolve_attribute_for_instrument_data("read_orientation", 1, $id);
is_deeply(\@read_orientations, ["reverse_forward"],"Got a list in a list context");

$id->original_est_fragment_size(200);
$id2->original_est_fragment_size(300);

my $average = Genome::Model::Build::DeNovoAssembly->resolve_average_for_attribute(attribute => "original_est_fragment_size", objects => [$id, $id2]);
is($average, 250, "Average calculated correctly");

my $disk_usage_result_subclass_names = Genome::Model::Build::DeNovoAssembly->_disk_usage_result_subclass_names;
is(@$disk_usage_result_subclass_names, 1, 'disk_usage_result_subclass_names');

done_testing();
