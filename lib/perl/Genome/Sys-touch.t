#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More tests => 7;

use File::Spec qw();
use File::Temp qw();
use File::stat qw(stat);

use Genome;

my $tmp_dir = File::Temp->newdir();
my $filename = File::Spec->join($tmp_dir->dirname, 'somefile');

ok(! -e $filename, 'somefile does not exist yet');
Genome::Sys->touch($filename);
ok(-f $filename, 'somefile was created');
is(stat($filename)->size, 0, 'newly created file is empty');

my $file = IO::File->new($filename, 'w') or die $!;
$file->print('foo');
$file->close();
ok(stat($filename)->size > 0, 'newly created file now has size');

my $atime = stat($filename)->atime;
my $mtime = stat($filename)->mtime;
my $size  = stat($filename)->size;

sleep(1);

Genome::Sys->touch($filename);

ok(stat($filename)->atime > $atime, 'touch advanced atime');
ok(stat($filename)->mtime > $mtime, 'touch advanced mtime');
is(stat($filename)->size,    $size, 'touch did not overwrite existing file');
