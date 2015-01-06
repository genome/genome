#!/usr/bin/env genome-perl

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
    $ENV{UR_COMMAND_DUMP_STATUS_MESSAGES} = 1;
}

use strict;
use warnings;

use above 'Genome';

use Test::More;

use_ok('Genome::Model::Tools::Picard::WithRequiredDownsampleRatio') or die;

class WithRequiredDownsampleRatio {
    is => 'Genome::Model::Tools::Picard::WithRequiredDownsampleRatio', 
};
my $downsample_ratio_property = WithRequiredDownsampleRatio->__meta__->property_meta_for_name('downsample_ratio');
ok($downsample_ratio_property, 'downsample_ratio property');
ok(!$downsample_ratio_property->is_optional, 'downsample_ratio is not optional')
    or diag explain $downsample_ratio_property;

my $tester = WithRequiredDownsampleRatio->create();
ok($tester, 'create');
my @errors = $tester->__errors__;
is(@errors, 1, 'errors for downsample_ratio of NA');
is($errors[0]->__display_name__, "INVALID: property 'downsample_ratio': No value specified for required property", 'correct error desc for downsample_ratio of NA');

done_testing();
