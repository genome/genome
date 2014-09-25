#! /gsc/bin/perl

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
    $ENV{UR_COMMAND_DUMP_STATUS_MESSAGES} = 1;
}

use strict;
use warnings;

use above 'Genome';

use Test::More;

use_ok('Genome::Model::Tools::Picard::WithDownsampleRatio') or die;

class DownsampleRatioTester {
    is => 'Genome::Model::Tools::Picard::WithDownsampleRatio', 
};
ok(DownsampleRatioTester->__meta__->property_meta_for_name('downsample_ratio')->is_optional, 'downsample_ratio is optional');

class DownsampleRatioIsRequiredTester {
    is => 'Genome::Model::Tools::Picard::WithDownsampleRatio', 
    has_optional_transient => { _make_downsample_ratio_required => {}, },
};
ok(!DownsampleRatioIsRequiredTester->__meta__->property_meta_for_name('downsample_ratio')->is_optional, 'downsample_ratio is not optional');

my $tester = DownsampleRatioTester->create(
    downsample_ratio => 'NA',
);
ok($tester, 'create');
my @errors = $tester->__errors__;
ok(@errors, 'errors for downsample_ratio of NA');
is($errors[0]->desc, 'Invalid Float value.', 'correct error desc for downsample_ratio of NA');

$tester = DownsampleRatioTester->create(
    downsample_ratio => 0,
);
ok($tester, 'create');
@errors = $tester->__errors__;
ok(@errors, 'errors for downsample_ratio of 0');
is($errors[0]->desc, 'Must be greater than 0 and less than 1! 0', 'correct error desc for downsample_ratio of 0');

$tester = DownsampleRatioTester->create(
    downsample_ratio => 1,
);
ok($tester, 'create');
@errors = $tester->__errors__;
ok(@errors, 'errors for downsample_ratio of 1');
is($errors[0]->desc, 'Must be greater than 0 and less than 1! 1', 'correct error desc for downsample_ratio of 1');

$tester = DownsampleRatioTester->create(
    downsample_ratio => 0.333333,
);
ok($tester, 'create');
@errors = $tester->__errors__;
ok(!@errors, 'no errors for valid downsample_ratio');
is($tester->downsample_ratio, 0.333333, 'correctly set the down sample ratio from 0.333333 to 0.33');

done_testing();
