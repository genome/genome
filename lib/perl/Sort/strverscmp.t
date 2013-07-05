use strict;
use warnings;

use Test::More tests => 4;

use Sort::strverscmp;

subtest 'decompose_version' => sub {
    my @cases = (
        ['abc123def', 'abc'   , '123', 'def'],
        ['abc123'   , 'abc'   , '123', ''   ],
        ['abcdef'   , 'abcdef',    '', ''   ],
        ['123def'   , ''      , '123', 'def'],
    );
    plan tests => scalar(@cases);

    for my $case (@cases) {
        my ($orig, $estr, $enum, $erem) = @$case;

        subtest qq(Case: $orig) => sub {
            plan tests => 3;

            my ($str, $num, $rem) = Sort::strverscmp::decompose_version($orig);
            is($str, $estr, 'string matched');
            is($num, $enum, 'number matched');
            is($rem, $erem, 'remainder matched');
        };
    }
};

subtest 'decompose_fractional' => sub {
    my @cases = (
        ['009', '00',  '9'],
        [ '90', ''  , '90'],
        [  '9', ''  ,  '9'],
    );
    plan tests => scalar(@cases);

    for my $case (@cases) {
        my ($orig, $ezeroes, $enum, $erem) = @$case;

        subtest qq(Case: $orig) => sub {
            plan tests => 2;

            my ($zeroes, $num, $rem) = Sort::strverscmp::decompose_fractional($orig);
            is($zeroes, $ezeroes, 'zeroes matched');
            is($num, $enum, 'number matched');
        };
    }
};

subtest 'GNU strverscmp examples' => sub {
    plan tests => 5;

    is(strverscmp('no digit', 'no digit'), 0, q('no digit' == 'no digit'));
    is(strverscmp('item#99', 'item#100'), -1, q('item#99' < 'item#100'));
    is(strverscmp('alpha1', 'alpha001'), 1, q('alpha1' > 'alpha001'));
    is(strverscmp('part1_f012', 'part1_f01'), 1, q('part1_f012' > 'part1_f01'));
    is(strverscmp('foo.009', 'foo.0'), -1, q('foo.009' < 'foo.0'));
};

subtest 'custom examples' => sub {
    plan tests => 4;

    is(strverscmp('alpha1', 'beta1'), -1, q('alpha1' < 'beta1'));
    is(strverscmp('g', 'G'), 1, q('g' > 'G'));
    is(strverscmp('1.0.5', '1.0.50'), -1, q('1.0.5' < '1.0.50'));
    is(strverscmp('1.0.5', '1.05'), 1, q('1.0.5' > '1.05'));
};
