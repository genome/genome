#!/usr/bin/env perl

use above 'Genome';
use Genome::File::BedPe::Reader;
use Genome::File::BedPe::Entry;

use Test::More;
use IO::String;

use strict;
use warnings;

my $pkg = 'Genome::File::BedPe::Writer';

use_ok($pkg);

my $bedpe_txt_with_header = <<EOF
#omg hi
#how are you?
1	2	3	4	5	6	hi	2	+	-	a	b	c	d	e
7	8	9	19	11	12	bye	3	-	+	f	g	h	i	j
EOF
;

my $bedpe_txt = <<EOF
1	2	3	4	5	6	hi	2	+	-	a	b	c	d	e
7	8	9	19	11	12	bye	3	-	+	f	g	h	i	j
EOF
;

subtest "Write entry" => sub {
    my $in_fh = new IO::String($bedpe_txt);
    my $out_str;
    my $out_fh = new IO::String($out_str);

    my $reader = Genome::File::BedPe::Reader->fhopen($in_fh, "test");
    my $writer = $pkg->fhopen($out_fh, "test", $reader->header);
    my @entries;
    while (my $entry = $reader->next) {
        push(@entries, $entry);
    }

    is(scalar @entries, 2, 'Got 2 entries');
    $writer->write($_) for (@entries);
    is($out_str, $bedpe_txt);
};


subtest "Write with header" => sub {
    my $in_fh = new IO::String($bedpe_txt_with_header);
    my $out_str;
    my $out_fh = new IO::String($out_str);

    my $reader = Genome::File::BedPe::Reader->fhopen($in_fh, "test");
    my $writer = $pkg->fhopen($out_fh, "test", $reader->header);
    my @entries;
    while (my $entry = $reader->next) {
        push(@entries, $entry);
    }

    is(scalar @entries, 2, 'Got 2 entries');
    $writer->write($_) for (@entries);
    is($out_str, $bedpe_txt_with_header);
};



done_testing();
