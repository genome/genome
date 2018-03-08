#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use File::Compare 'compare';
use Storable 'retrieve';
use Test::More tests => 8;

use_ok('Genome::Model::Tools::Consed::Navigation::Reader')
    or die;

my $dir = Genome::Config::get('test_inputs') . '/Genome-Consed';
ok(-d $dir, "Test dir ($dir) exists");
my $nav = $dir.'/repeats.nav';
ok(-f $nav, "Nav file ($nav) exists");
my $navs = retrieve($dir.'/navs.stor');
ok($navs, 'Got navs from stor file');

my $reader = Genome::Model::Tools::Consed::Navigation::Reader->create(
    input => $nav,
);
ok($reader, 'Created nav reader');
is($reader->title, $navs->{title}, 'Title matches');
my @navs = $reader->all;
ok(@navs, 'Got navs');
is_deeply(\@navs, $navs->{navs}, 'Navs match');



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
