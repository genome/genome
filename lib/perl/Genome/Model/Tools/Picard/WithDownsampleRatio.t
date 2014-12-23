#!/usr/bin/env genome-perl

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
    $ENV{UR_COMMAND_DUMP_STATUS_MESSAGES} = 1;
}

use strict;
use warnings;

use above 'Genome';

use Test::More;

use_ok('Genome::Model::Tools::Picard::WithDownsampleRatio') or die;

class WithDownsampleRatio {
    is => 'Genome::Model::Tools::Picard::WithDownsampleRatio', 
};
my $downsample_ratio_property = WithDownsampleRatio->__meta__->property_meta_for_name('downsample_ratio');
ok($downsample_ratio_property, 'downsample_ratio property');
ok($downsample_ratio_property->is_optional, 'downsample_ratio is not optional')
    or diag explain $downsample_ratio_property;

my $tester = WithDownsampleRatio->create(
    downsample_ratio => 'NA',
);
ok($tester, 'create');
my @errors = $tester->__errors__;
is(@errors, 1, 'errors for downsample_ratio of NA');
is($errors[0]->__display_name__, "INVALID: property 'downsample_ratio': Invalid Float value.", 'correct error __display_name__ for downsample_ratio of NA');

$tester = WithDownsampleRatio->create(
    downsample_ratio => 0,
);
ok($tester, 'create');
@errors = $tester->__errors__;
is(@errors, 1, 'errors for downsample_ratio of 0');
is($errors[0]->__display_name__, "INVALID: property 'downsample_ratio': Must be greater than 0 and less than 1! 0", 'correct error __display_name__ for downsample_ratio of 0');

$tester = WithDownsampleRatio->create(
    downsample_ratio => 1,
);
ok($tester, 'create');
@errors = $tester->__errors__;
is(@errors, 1, 'errors for downsample_ratio of 1');
is($errors[0]->__display_name__, "INVALID: property 'downsample_ratio': Must be greater than 0 and less than 1! 1", 'correct error __display_name__ for downsample_ratio of 1');

$tester = WithDownsampleRatio->create(
    downsample_ratio => 0.333333,
);
ok($tester, 'create');
@errors = $tester->__errors__;
ok(!@errors, 'no errors for valid downsample_ratio');
is($tester->downsample_ratio, 0.333333, 'correctly set the down sample ratio from 0.333333 to 0.33');

done_testing();
