use strict;
use warnings;

use above "Genome";
use Test::More tests => 2;

use Term::ANSIColor;

BEGIN {
    use_ok 'Genome::Utility::Test', qw(strip_ansi);
}


my $test_string = color 'bold blue';
$test_string .= 'Hi there in blue';
$test_string .= color 'reset';
$test_string .= ' uncolored ';
$test_string .= colored 'more ugly colors', 'underline yellow on_magenta';

#note($test_string);
is(strip_ansi($test_string),
    'Hi there in blue uncolored more ugly colors',
    'Removed color markup');
