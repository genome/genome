#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Cwd 'abs_path';
use Data::Dumper 'Dumper';
use File::Compare 'compare';
use File::Temp 'tempdir';
use Test::More tests => 15;
use Storable 'retrieve';

use_ok('Genome::Utility::IO::Writer');

#< MAKE SOME WRITERS FOR TESTING >#
class Album::Writer {
    is => 'Genome::Utility::IO::Writer',
    has_optional => [
    headers => {
        is => 'List',
    },
    ],
};
no warnings;
*Album::Writer::write_one = sub{ 
    my ($self, $album) = @_;
    $self->output->print(join(',', map { $album->{$_} } @{$self->headers}), "\n");
    return 1;
};
use warnings;

class NoWriteOne {
    is => 'Genome::Utility::IO::Writer',
};

my $dir = Genome::Config::get('test_inputs') . '/Genome-Utility-IO';
ok(-d $dir, "Test dir ($dir) exists");
my $albums_file = $dir.'/albums.csv';
ok(-f $albums_file, "Albums csv file ($albums_file) exists");
my $albums_stor = $dir.'/albums.stor';
ok(-f $albums_stor, "Albums stor file ($albums_stor) exists");
my $albums = retrieve($albums_stor);
ok($albums, "Got albums from stor file");
my $tmp_dir = tempdir(CLEANUP => 1);
ok(-d $tmp_dir, 'Got temp dir');

#< CLASS FOR OUPUT THAT CAN'T GETLINE >#
class IO::NoPrint {
    is => 'UR::Object',
};

#< VALID WRITER W/ OBJECT >#
my $writer = Album::Writer->create(
    output => IO::String->new,
);
ok($writer, 'Created writer w/ temp fh');
$writer->delete;

#< VALID WRITER W/ STDOUT >#
$writer = Album::Writer->create();
ok($writer, 'Created writer w/ STDOUT');
$writer->delete;

#< VALID WRITER W/ FILE >#
my $tmp_file = $tmp_dir.'/albums.csv';
$writer = Album::Writer->create(
    output => $tmp_file,
);
ok($writer, "Created writer w/ tmp file ($tmp_file)");
is($writer->get_original_output, abs_path($tmp_file), 'get_original_output');
$writer->headers($albums->{headers});
ok($writer->output->print(join(',', @{$writer->headers}),"\n"), "Wrote headers w/ 'print'");
for my $album ( @{$albums->{albums}} ) {
    $writer->write_one($album);
}
is(compare($tmp_file, $albums_file), 0, 'Generated and expected album files match');
$writer->delete;

#< INVALID PARAMS >#
$writer = Album::Writer->create(output => $albums_file);
ok(!$writer, 'Failed as expected - create writer w/ existing file');
$writer = Album::Writer->create(output => IO::NoPrint->create());
ok(!$writer, 'Failed as expected - create writer w/ object cant "print"');

#< INVALID writer >#
$writer = NoWriteOne->create();
ok(!$writer, 'Failed as expected - create writer w/o "write_one" method');



#########

=pod

=head1 Tests

=head1 Disclaimer

 Copyright (C) 2006 Washington University Genome Sequencing Center

 This script is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY
 or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
 License for more details.

=head1 Author(s)

 Eddie Belter <ebelter@watson.wustl.edu>

=cut
