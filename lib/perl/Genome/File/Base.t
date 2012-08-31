#!/usr/bin/env genome-perl
use strict;
use warnings;
use above "Genome";
use Test::More tests => 4;

my $path = Genome::Sys->create_temp_file_path;

class Genome::File::Foo { is => 'Genome::File::Base' };

my $f = Genome::File::Foo->get($path);
ok($f, "got a file object for a path");

my $fh = $f->open("w");
ok($fh, "open a write handle");

$fh->print(">abc\nAGCT\n-\n1234\n");
$fh->close;

my $fh2 = $f->open();
ok($fh2, "open a read handle");

my @lines = $fh2->getlines();
my $lines = join("",@lines);
is($lines, ">abc\nAGCT\n-\n1234\n", "content is correct");
