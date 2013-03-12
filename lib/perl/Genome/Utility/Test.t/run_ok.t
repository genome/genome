use strict;
use warnings;

use Test::Builder::Tester;
use above "Genome";
use Genome::Utility::Test qw(run_ok);
use Test::More;

test_out('not ok 1 - false');
test_fail(+1);
run_ok('false');
test_test('run_ok fails on non-zero exit code');

test_out('ok 1 - true');
run_ok('true');
test_test('run_ok passes on zero exit code');

test_out('ok 1 - maybe');
run_ok('true', 'maybe');
test_test('run_ok sets test name correctly');

done_testing();
