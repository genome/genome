#!/usr/bin/env perl

use strict;
use warnings;

use above 'Genome';
use Test::More;
use Genome::File::BedPe::Header;
use Genome::File::BedPe::Entry;

my $pkg = 'Genome::File::BedPe::Entry';

use_ok($pkg);

my $entry_txt = join("\t", qw(chr1 10 20 chr2 30 40 hello 45 + - a b c));

subtest "create" => sub {
    my $header = new Genome::File::BedPe::Header([]);
    $header->set_custom_fields("c1", "c2", "c3");
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
    is($entry->custom_by_name("c1"), 'a', 'custom_by_name c1');
    is($entry->custom_by_name("c2"), 'b', 'custom_by_name c2');
    is($entry->custom_by_name("c3"), 'c', 'custom_by_name c3');

    eval {
        $entry->custom_by_name("banana");
    };
    ok($@, "asking for unknown custom field is an error");
};

subtest "to_string" => sub {
    my $header = new Genome::File::BedPe::Header([]);
    my $entry = $pkg->new($header, $entry_txt);
    is($entry->to_string, $entry_txt, "to_string works");

    my $txt = join("\t", qw(chr1 10 20 chr2 30 40));
    $entry = $pkg->new($header, $txt);
    is($entry->to_string, $txt, "to_string works without optional fields");

    $txt = join("\t", qw(chr1 10 20 chr2 30 40 50));
    $entry = $pkg->new($header, $txt);
    is($entry->to_string, $txt, "to_string works without optional fields");

};

subtest "new_from_fields" => sub {
    my $header = new Genome::File::BedPe::Header([]);
    $header->set_custom_fields("foo", "bar");
    my %fields = (
        chrom1 => 1, start1 => 2, end1 => 3,
        chrom2 => 4, start2 => 5, end2 => 6,
        name => 7,
        score => 8,
        strand1 => 9,
        strand2 => 10,
    );
    my %custom = (
        foo => "11",
        bar => "12",
    );


    my $entry = Genome::File::BedPe::Entry->new_from_fields($header, %fields, %custom);

    is(ref $entry, "Genome::File::BedPe::Entry", "new_from_fields yields entry");
    for my $k (keys %fields) {
        is($entry->{$k}, $fields{$k}, "$k set correctly");
    }

    is($entry->{custom}->[0], "11", "foo custom field");
    is($entry->{custom}->[1], "12", "bar custom field");

    is($entry->to_string, join("\t", 1..12), "to_string");
};

subtest "new_from_fields custom gap" => sub {
    my $header = new Genome::File::BedPe::Header([]);
    $header->set_custom_fields("foo", "bar", "baz");
    my %fields = (
        chrom1 => 1, start1 => 2, end1 => 3,
        chrom2 => 4, start2 => 5, end2 => 6,
        name => 7,
        score => 8,
        strand1 => 9,
        strand2 => 10,
    );
    my %custom = (
        baz => "13",
    );


    my $entry = Genome::File::BedPe::Entry->new_from_fields($header, %fields, %custom);

    is(ref $entry, "Genome::File::BedPe::Entry", "new_from_fields yields entry");
    for my $k (keys %fields) {
        is($entry->{$k}, $fields{$k}, "$k set correctly");
    }

    is($entry->{custom}->[0], ".", "foo custom field");
    is($entry->{custom}->[1], ".", "bar custom field");
    is($entry->{custom}->[2], "13", "baz custom field");

    is($entry->to_string, join("\t", 1..10, '.', '.', 13), "to_string");
};

subtest "new_from_fields missing required" => sub {
    my $header = new Genome::File::BedPe::Header([]);
    $header->set_custom_fields("foo");

    my %required = (
        chrom1 => 1, start1 => 2, end1 => 3,
        chrom2 => 4, start2 => 5, end2 => 6,
    );

    my %optional = (
        name => 7,
        score => 8,
        strand1 => 9,
        strand2 => 10,
    );

    for my $k (keys %required) {
        my %h = %required;
        delete $h{$k};
        eval {
            my $entry = Genome::File::BedPe::Entry->new_from_fields($header, %h);
        };
        ok($@, "omitting $k is an error");
    }

    for my $k (keys %optional) {
        my %h = %optional;
        my @expected = 1..11;
        $expected[$h{$k} - 1] = '.';
        delete $h{$k};

        # We use foo to make sure that field 10 gets set to . when missing
        my $entry = Genome::File::BedPe::Entry->new_from_fields($header, %required, %h, foo => "11");
        is(ref $entry, "Genome::File::BedPe::Entry", "omitting $k is not an error");
        is($entry->to_string, join("\t", @expected), "entry to_string");
    }
};


done_testing();
