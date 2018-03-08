#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use File::Compare 'compare';
use File::Temp 'tempdir';
use Test::More tests => 7;

use_ok('Genome::Model::Tools::Consed::NavToAcenav')
    or die;

my $dir = Genome::Config::get('test_inputs') . '/Genome-Consed';
ok(-d $dir, "Test dir ($dir) exists");
my $nav = $dir.'/repeats.nav';
ok(-f $nav, "List file ($nav) exists");
my $expected_acenav = $dir.'/repeats.acenav';
my $tmp_dir = tempdir(CLEANUP => 1);
ok(-d $tmp_dir, "Temp dir ($tmp_dir) exists");
my $acenav = $tmp_dir.'/repeats.acenav';
my $nav2acenav = Genome::Model::Tools::Consed::NavToAcenav->create(
    nav => $nav,
    acenav => $acenav,
    acefile => '/home1/watson/seqmgr/M_BB0392D19/edit_dir/M_BB0392D19.fasta.screen.ace.old',
);
ok($nav2acenav, 'Created list to nav');
ok($nav2acenav->execute, 'execute');
is(compare($acenav, $expected_acenav), 0, 'Expected and generated acenav files match');



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
