#!/usr/bin/env perl

use Test::More;
use above 'Genome';
use Genome::File::BedPe::Header;

use strict;
use warnings;

my $pkg = 'Genome::File::BedPe::Entry';

use_ok($pkg);

my $header = new Genome::File::BedPe::Header([]);
my $entry_txt = join("\t", qw(chr1 10 20 chr2 30 40 hello 45 + - a b c));

subtest "create" => sub {
    my $entry = $pkg->new($header, $entry_txt);
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
};

subtest "to_string" => sub {
    my $entry = $pkg->new($header, $entry_txt);
    is($entry->to_string, $entry_txt, "to_string works");

    my $txt = join("\t", qw(chr1 10 20 chr2 30 40));
    $entry = $pkg->new($header, $txt);
    is($entry->to_string, $txt, "to_string works without optional fields");

    $txt = join("\t", qw(chr1 10 20 chr2 30 40 50));
    $entry = $pkg->new($header, $txt);
    is($entry->to_string, $txt, "to_string works without optional fields");

};

done_testing();
