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
is($entry->{fields}{chrom1}, "chr1");
is($entry->{fields}{start1}, 10);
is($entry->{fields}{end1}, 20);
is($entry->{fields}{chrom2}, "chr2");
is($entry->{fields}{start2}, 30);
is($entry->{fields}{end2}, 40);
is($entry->{fields}{name}, "hello");
is($entry->{fields}{score}, 45);
is($entry->{fields}{strand1}, '+');
is($entry->{fields}{strand2}, '-');
is($entry->{fields}{custom}[0], 'a');
is($entry->{fields}{custom}[1], 'b');
is($entry->{fields}{custom}[2], 'c');

done_testing();
