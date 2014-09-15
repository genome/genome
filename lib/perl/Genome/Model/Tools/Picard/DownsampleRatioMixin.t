#! /gsc/bin/perl

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
    $ENV{UR_COMMAND_DUMP_STATUS_MESSAGES} = 1;
}

use strict;
use warnings;

use above 'Genome';

use Test::More;

use_ok('Genome::Model::Tools::Picard::DownsampleRatioMixin') or die;

class DownsampleRatioMixinTester {
    has => {
        downsample_ratio => Genome::Model::Tools::Picard::DownsampleRatioMixin->downsample_ratio_property,
    },
};
sub DownsampleRatioMixinTester::__errors__ {
    my $self = shift;
    return Genome::Model::Tools::Picard::DownsampleRatioMixin::__errors__($self);
}

my $tester = DownsampleRatioMixinTester->create(
    downsample_ratio => 'NA',
);
ok($tester, 'create');
my @errors = $tester->__errors__;
ok(@errors, 'errors for downsample_ratio of NA');
is($errors[0]->desc, 'Invalid number! NA', 'correct error desc for downsample_ratio of NA');

$tester = DownsampleRatioMixinTester->create(
    downsample_ratio => 0,
);
ok($tester, 'create');
@errors = $tester->__errors__;
ok(@errors, 'errors for downsample_ratio of 0');
is($errors[0]->desc, 'Must be greater than 0 and less than 1! 0', 'correct error desc for downsample_ratio of 0');

$tester = DownsampleRatioMixinTester->create(
    downsample_ratio => 1,
);
ok($tester, 'create');
@errors = $tester->__errors__;
ok(@errors, 'errors for downsample_ratio of 1');
is($errors[0]->desc, 'Must be greater than 0 and less than 1! 1', 'correct error desc for downsample_ratio of 1');

$tester = DownsampleRatioMixinTester->create(
    downsample_ratio => 0.333333,
);
ok($tester, 'create');
@errors = $tester->__errors__;
ok(!@errors, 'no errors for valid downsample_ratio');
is($tester->downsample_ratio, 0.33, 'correctly set the down sample ratio from 0.333333 to 0.33');

done_testing();
