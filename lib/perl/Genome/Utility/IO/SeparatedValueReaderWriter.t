#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

require File::Compare;
use Test::More;
use Test::Deep qw(cmp_bag);
use Test::Exception;

use_ok('Genome::Utility::IO::SeparatedValueReader') or die;
use_ok('Genome::Utility::IO::SeparatedValueWriter') or die;

my $test_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Utility-IO';
my $albums_csv = $test_dir.'/albums.csv';
my $albums_no_headers = $test_dir.'/albums.no_headers.csv';
my $albums_regexp = $test_dir.'/albums.test_regexp.csv';

my $tmpdir = File::Temp::tempdir(CLEANUP => 1);
my $output = $tmpdir.'/albums';

# NORMAL
my $reader = Genome::Utility::IO::SeparatedValueReader->create(
    input => $albums_csv,
);
ok($reader, 'create csv reader');
my $headers = $reader->headers;
my $writer = Genome::Utility::IO::SeparatedValueWriter->create(
    output => $output,
    headers => $reader->headers,
);
ok($writer, 'create csv writer');
while ( my $album = $reader->next ) {
    $writer->write_one($album);
}
$writer->output->flush;
is(File::Compare::compare($output, $albums_csv), 0, 'read/write file matches');

is($reader->line_number, 5, 'Line number incremented successfully');
ok($reader->reset, 'reset');
is($reader->line_number, 1, 'Reset to line 1');

# NO HEADERS
$reader = Genome::Utility::IO::SeparatedValueReader->create(
    input => $albums_no_headers,
    headers => $headers,
);
ok($reader, 'create reader w/ file w/o headers');
$headers = $reader->headers;
unlink $output;
$writer = Genome::Utility::IO::SeparatedValueWriter->create(
    output => $output,
    headers => $reader->headers,
);
ok($writer, 'create writer');
while ( my $album = $reader->next ) {
    $writer->write_one($album);
}
$writer->output->flush;
is(File::Compare::compare($output, $albums_csv), 0, 'read/write no headers file matches');

# REGEXP
$reader = Genome::Utility::IO::SeparatedValueReader->create(
    input => $albums_regexp,
    separator => ',+',
    is_regex => 1,
);
ok($reader, 'create reader w/ regexp file');
$headers = $reader->headers;
unlink $output;
$writer = Genome::Utility::IO::SeparatedValueWriter->create(
    output => $output,
    headers => $headers,
);
ok($writer, 'create writer');
while ( my $album = $reader->next ) {
    $writer->write_one($album);
}
is(File::Compare::compare($output, $albums_csv), 0, 'read/write regexp file matches');

# WRITER FAILS
my %fails = (
    'undef' => undef,
    'empty hash ref' => {},
    'not a hash ref' => 'a string',
    'different headers' => { different => 'headers' },
);
for my $desc ( keys %fails ) {
    dies_ok(sub {$writer->write_one($fails{$desc})}, "Failed as expected, tried to 'write one' w/ $desc");
}

# READER FAILS
$reader = Genome::Utility::IO::SeparatedValueReader->create(
    input => $albums_no_headers,
    headers => [qw/ not the right number of headers /],
);
ok($reader, 'create reader');
dies_ok(sub {$reader->next}, "line with too many columns causes die");

# This should fail because we ignore extra columns but we dont have the minimum
$reader = Genome::Utility::IO::SeparatedValueReader->create(
    input => $albums_no_headers,
    headers => [qw/ dont have data to fill all these columns /],
    allow_extra_columns => 1,
);
ok($reader, 'Created SVR to test too few columns while ignoring extra columns');
dies_ok(sub {$reader->next}, 'Failed as expected - next');

# This should succeed because we ignore extra columns
$reader = Genome::Utility::IO::SeparatedValueReader->create(
    input => $albums_no_headers,
    headers => [qw/ have enough data /],
    allow_extra_columns => 1,
);
ok($reader, 'Created SVR to test too many columns while ignoring extra columns');
ok($reader->next, 'Succeeded as expected');
my $extra = $reader->current_extra_columns;
ok(scalar(@$extra), 'extra columns set as expected');
cmp_bag($extra, ['legend', 'reggae'], 'extra columns are the trailing ones in the file');

#print "$tmpdir\n"; <STDIN>;
done_testing();
