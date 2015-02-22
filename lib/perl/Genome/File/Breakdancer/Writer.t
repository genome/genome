#!/usr/bin/env perl

use above 'Genome';
use Genome::File::Breakdancer::Reader;
use Genome::File::Breakdancer::Entry;

use Test::More;
use IO::String;

use strict;
use warnings;

my $pkg = 'Genome::File::Breakdancer::Writer';

use_ok($pkg);

my $txt_with_header = <<EOF
#omg hi
#how are you?
1	10	2+0-	1	20	0+2-	DEL	10	35	2	lib1|2,2.89
1	11	16+5-	1	21	16+5-	INS	-10	32	3	lib1|1,0.01:lib2|2,NA
1	12	3+8-	1	22	3+8-	ITX	-10	43	3	lib1|1,NA:lib2|2,NA
EOF
;

my $txt = <<EOF
1	10	2+0-	1	20	0+2-	DEL	10	35	2	lib1|2,2.89
1	11	16+5-	1	21	16+5-	INS	-10	32	3	lib1|1,0.01:lib2|2,NA
1	12	3+8-	1	22	3+8-	ITX	-10	43	3	lib1|1,NA:lib2|2,NA
EOF
;

subtest "Write entry" => sub {
    my $in_fh = new IO::String($txt);
    my $out_str;
    my $out_fh = new IO::String($out_str);

    my $reader = Genome::File::Breakdancer::Reader->fhopen($in_fh, "test");
    my $writer = $pkg->fhopen($out_fh, "test", $reader->header);
    my @entries;
    while (my $entry = $reader->next) {
        push(@entries, $entry);
    }

    is(scalar @entries, 3, 'Got 2 entries');
    $writer->write($_) for (@entries);
    is($out_str, $txt);
};


subtest "Write with header" => sub {
    my $in_fh = new IO::String($txt_with_header);
    my $out_str;
    my $out_fh = new IO::String($out_str);

    my $reader = Genome::File::Breakdancer::Reader->fhopen($in_fh, "test");
    my $writer = $pkg->fhopen($out_fh, "test", $reader->header);
    my @entries;
    while (my $entry = $reader->next) {
        push(@entries, $entry);
    }

    is(scalar @entries, 3, 'Got 2 entries');
    $writer->write($_) for (@entries);
    is($out_str, $txt_with_header);
};



done_testing();
