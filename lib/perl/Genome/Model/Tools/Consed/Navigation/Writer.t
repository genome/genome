#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use File::Compare 'compare';
use File::Temp 'tempdir';
use Storable 'retrieve';
use Test::More tests => 8;

use_ok('Genome::Model::Tools::Consed::Navigation::Writer')
    or die;

my $dir = Genome::Config::get('test_inputs') . '/Genome-Consed';
ok(-d $dir, "Test dir ($dir) exists");
my $expected_nav = $dir.'/repeats.nav';
ok(-f $expected_nav, "Nav file ($expected_nav) exists");
my $tmp_dir = tempdir(CLEANUP => 1);
ok(-d $tmp_dir, "Temp dir ($tmp_dir) exists");
my $navs = retrieve($dir.'/navs.stor');
ok($navs, 'Got navs from stor file');

my $nav = $tmp_dir.'/repeats.nav';
my $writer = Genome::Model::Tools::Consed::Navigation::Writer->create(
    output => $nav,
    title => $navs->{title},
);
ok($writer, 'Created writer');
my $count = 0;
for my $nav ( @{$navs->{navs}} ) {
    $writer->write_one($nav) and $count++;
}
is($count, @{$navs->{navs}}, 'Wrote all navs');
is(compare($nav, $expected_nav), 0, 'Expected and generated nav file matches');



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
