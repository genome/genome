#!/usr/bin/env genome-perl

use strict;
use warnings;
use Test::More tests => 44;
use above "Genome";
use Genome::Info::IUB;

#Test variant_alleles_for_iub
is(Genome::Info::IUB::variant_alleles_for_iub(undef, undef),undef, "variant_alleles_for_iub: undefined inputs return undef");  #test to make sure it handles undefined values
is(Genome::Info::IUB::variant_alleles_for_iub(undef, 'A'), undef, "variant_alleles_for_iub: undefined ref returns undef");
is(Genome::Info::IUB::variant_alleles_for_iub('A', undef),undef, "variant_alleles_for_iub: undefined iub returns undef");
is(Genome::Info::IUB::variant_alleles_for_iub('N', 'A'), undef, "variant_alleles_for_iub: ambiguous ref returns undef");
cmp_warn(sub { is(Genome::Info::IUB::variant_alleles_for_iub('Q', 'A'),
                undef, "variant_alleles_for_iub: invalid ref returns undef")
        },
        qr{Ambiguous reference bases \(Q\) not currently supported},
        'expected warning');
is(Genome::Info::IUB::variant_alleles_for_iub('G', 'Q'), undef, "variant_alleles_for_iub: invalid iub returns undef");
my @alleles = Genome::Info::IUB::variant_alleles_for_iub('G','W');
is_deeply(\@alleles, ['A','T'] , "variant_alleles_for_iub: returns values on proper input");

#test to make sure it handles undefined values
cmp_warn(sub {
            is(Genome::Info::IUB::iub_for_alleles(undef),
                undef, "iub_for_alleles: undefined inputs return undef")
        },
        qr{Conversion of more than 2 alleles to IUB code is currently unsupported \(1 passed\)},
        'expected warning');

cmp_warn(sub {
            is(Genome::Info::IUB::iub_for_alleles('A'),
                undef, "iub_for_alleles: invalid input returns undef");
        },
        qr{Conversion of more than 2 alleles to IUB code is currently unsupported \(1 passed\)},
        'expected warning');

cmp_warn(sub {
            is(Genome::Info::IUB::iub_for_alleles('A','T','G'),
                undef, "iub_for_alleles: unsupported input return undef");
        },
        qr{Conversion of more than 2 alleles to IUB code is currently unsupported \(3 passed\)},
        'expected warning');

is(Genome::Info::IUB::iub_for_alleles('A','N'), 'N', "iub_for_alleles: allele with 'N' returns 'N'");
is(Genome::Info::IUB::iub_for_alleles('a','T'), 'W', "iub_for_alleles: case insensitivity and valid result");


is(Genome::Info::IUB::iub_to_alleles(undef),undef, "iub_to_alleles: undefined inputs return undef");  #test to make sure it handles undefined values
is(Genome::Info::IUB::iub_to_alleles('Q'), undef, "iub_to_alleles: invalid input returns undef");
is_deeply([Genome::Info::IUB::iub_to_alleles('A')],['A','A'], "iub_to_alleles: returns two alleles for homozygote");
@alleles = sort(Genome::Info::IUB::iub_to_alleles('w'));
is_deeply(\@alleles, ['A','T'], "iub_to_alleles: case insensitivity and valid result");


is(Genome::Info::IUB::iub_to_bases(undef),undef, "iub_to_bases: undefined inputs return undef");  #test to make sure it handles undefined values
is(Genome::Info::IUB::iub_to_bases('Q'), undef, "iub_to_bases: invalid input returns undef");
is_deeply([Genome::Info::IUB::iub_to_bases('A')],['A'], "iub_to_bases: returns one base for homozygote");
my @bases = sort(Genome::Info::IUB::iub_to_bases('w'));
is_deeply(\@bases, ['A','T'], "iub_to_bases: case insensitivity and valid result");

is(Genome::Info::IUB::iub_to_string(undef),undef, "iub_to_string: undefined inputs return undef");  #test to make sure it handles undefined values
is(Genome::Info::IUB::iub_to_bases('Q'), undef, "iub_to_string: invalid input returns undef");
is(Genome::Info::IUB::iub_to_string('d'), 'AGT', 'iub_to_string: case insensitivity and valid result');

#Test reference_iub_to_base
is(Genome::Info::IUB::reference_iub_to_base(undef), undef, "reference_iub_to_base: undefined inputs return undef");
is(Genome::Info::IUB::reference_iub_to_base('A'), 'A', "reference_iub_to_base: A returns A");
is(Genome::Info::IUB::reference_iub_to_base('G'), 'G', "reference_iub_to_base: G returns G");
is(Genome::Info::IUB::reference_iub_to_base('C'), 'C', "reference_iub_to_base: C returns C");
is(Genome::Info::IUB::reference_iub_to_base('T'), 'T', "reference_iub_to_base: T returns T");
is(Genome::Info::IUB::reference_iub_to_base('R'), 'A', "reference_iub_to_base: R returns A");
is(Genome::Info::IUB::reference_iub_to_base('Y'), 'C', "reference_iub_to_base: Y returns C");
is(Genome::Info::IUB::reference_iub_to_base('M'), 'A', "reference_iub_to_base: M returns A");
is(Genome::Info::IUB::reference_iub_to_base('K'), 'G', "reference_iub_to_base: K returns G");
is(Genome::Info::IUB::reference_iub_to_base('S'), 'C', "reference_iub_to_base: S returns C");
is(Genome::Info::IUB::reference_iub_to_base('W'), 'A', "reference_iub_to_base: W returns A");
is(Genome::Info::IUB::reference_iub_to_base('B'), 'C', "reference_iub_to_base: B returns C");
is(Genome::Info::IUB::reference_iub_to_base('D'), 'A', "reference_iub_to_base: D returns A");
is(Genome::Info::IUB::reference_iub_to_base('H'), 'A', "reference_iub_to_base: H returns A");
is(Genome::Info::IUB::reference_iub_to_base('V'), 'A', "reference_iub_to_base: V returns A");
is(Genome::Info::IUB::reference_iub_to_base('N'), 'A', "reference_iub_to_base: N returns A");
is(Genome::Info::IUB::reference_iub_to_base('ACTGRYMKSWBDHVN'), 'ACTGACAGCACAAAA', "reference_iub_to_base: ACTGRYMKSWBDHVN returns ACTGACAGCACAAAA");

#$HeadURL$
#$Id$

sub cmp_warn {
    my $do = shift;
    my $warn = shift;
    my $msg = shift;

    my $got_warn;
    local $SIG{__WARN__} = sub { $got_warn = shift };
    $do->();

    like($got_warn, $warn, $msg);
}

