use strict;
use warnings;

use Test::Builder::Tester;
use above "Genome";
use Genome::Utility::Test qw(capture_ok);
use Test::More;

test_out('not ok 1 - false');
test_err(q(/#\s+Failed test 'false'\n#\s+at .+ line \d+\.\n#\s+"false" unexpectedly returned exit value 1 at .* line \d+./));
capture_ok('false');
test_test('capture_ok fails on non-zero exit code');

test_out('ok 1 - true');
capture_ok('true');
test_test('capture_ok passes on zero exit code');

test_out('ok 1 - maybe');
capture_ok('true', 'maybe');
test_test('capture_ok sets test name correctly');

my $msg = 'hello';
test_out("ok 1 - echo $msg");
my ($ok, $output) = capture_ok(['echo', $msg]);
chomp $output;
test_test('ran capture_ok with expected output');
is($output, $msg, qq(capture_ok returned output, '$msg', as second argument));

done_testing();
