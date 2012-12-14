#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Data::Dumper 'Dumper';
use IO::String;
use Test::More;
use Storable 'retrieve';

use_ok('Genome::Utility::IO::Reader');

#< MAKE SOME READERS FOR TESTING >#
class Album::Reader {
    is => 'Genome::Utility::IO::Reader',
    has => [
    'next' => {
        calculate_from => [qw/ getline /],
        calculate => q| 
        return unless defined $getline;
        chomp $getline;
        my %album;
        @album{@{$self->headers}} = split(',', $getline);
        return \%album;
        |,
    },
    ],
    has_optional => [
    headers => {
        is => 'List',
    },
    ],
};

class NoNext::Reader {
    is => 'Genome::Utility::IO::Reader',
};

my $dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Utility-IO';
ok(-d $dir, "Test dir ($dir) exists");
my $albums_file = $dir.'/albums.csv';
ok(-f $albums_file, "Albums csv file ($albums_file) exists");
my $albums_stor = $dir.'/albums.stor';
ok(-f $albums_stor, "Albums stor file ($albums_stor) exists");
my $albums = retrieve($albums_stor);
ok($albums, "Got albums from stor file");

#< CLASS FOR INPUT THAT CAN'T GETLINE >#
class IO::NoGetline {
    is => 'UR::Object',
};

#< VALID READER W/ IO::FILE - TEST 'NEXT' >#
my $reader = Album::Reader->create(
    input => IO::File->new($albums_file, 'r'),
);
ok($reader, 'Created reader w/ IO::File');
$reader->headers($albums->{headers});
# test 'getline' and 'reset'
my $first_line = $reader->getline;
ok($first_line, "getline");
$reader->reset;
my $first_line_again = $reader->getline;
is($first_line, $first_line_again, "reset");
chomp $first_line;
is($first_line, join(',', @{$albums->{headers}}), 'First line matches headers');
# next
my @albums_from_next;
while ( my $album = $reader->next ) { push @albums_from_next, $album }
ok(@albums_from_next, 'next');
is_deeply(\@albums_from_next, $albums->{albums}, 'Albums from "next" match expected albums');

$reader->delete;

#< VALID READER W/ FILE - TEST 'ALL' >#
$reader = Album::Reader->create(
    input => $albums_file,
);
ok($reader, "Created reader w/ file ($albums_file)");
$reader->headers($albums->{headers});
$reader->getline; # header line
is($reader->get_original_input, $albums_file, 'get_original_input');
# all
my @albums_from_all = $reader->all;
ok(@albums_from_all, 'all');
is_deeply(\@albums_from_all, $albums->{albums}, 'Albums from "all" match expected albums');
$reader->delete;

#< INVALID PARAMS >#
$reader = Album::Reader->create();
ok(!$reader, 'Failed as expected - create reader w/o input');
$reader = Album::Reader->create(input => 'no_way_this_exists');
ok(!$reader, 'Failed as expected - create reader w/ invalid file');
$reader = Album::Reader->create(input => IO::NoGetline->create());
ok(!$reader, 'Failed as expected - create reader w/ object cant "getline"');

done_testing();
