#!/usr/bin/env perl

use Test::More;
use above 'Genome';
use Genome::File::BedPe::Header;

use strict;
use warnings;

my $pkg = 'Genome::File::BedPe::Entry';

use_ok($pkg);

my $hdr = new Genome::File::BedPe::Header([]);

my $entry = $pkg->new($hdr, join("\t", qw(chr1 10 20 chr2 30 40 hello 45 + - a b c)));
ok($entry, "Created dummy entry");
is($entry->{chrom1}, "chr1");
is($entry->{start1}, 10);
is($entry->{end1}, 20);
is($entry->{chrom2}, "chr2");
is($entry->{start2}, 30);
is($entry->{end2}, 40);
is($entry->{name}, "hello");
is($entry->{score}, 45);
is($entry->{strand1}, '+');
is($entry->{strand2}, '-');
is($entry->{custom}[0], 'a');
is($entry->{custom}[1], 'b');
is($entry->{custom}[2], 'c');

done_testing();
