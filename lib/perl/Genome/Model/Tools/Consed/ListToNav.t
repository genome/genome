#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use File::Compare 'compare';
use File::Temp 'tempdir';
use Test::More tests => 7;

use_ok('Genome::Model::Tools::Consed::ListToNav')
    or die;

my $dir = Genome::Config::get('test_inputs') . '/Genome-Consed';
ok(-d $dir, "Test dir ($dir) exists");
my $list = $dir.'/repeats.list';
ok(-f $list, "List file ($list) exists");
my $expected_nav = $dir.'/repeats.nav';
#ok(-f $expected_nav, "Nav file ($nav) exists");
my $tmp_dir = tempdir(CLEANUP => 1);
ok(-d $tmp_dir, "Temp dir ($tmp_dir) exists");

my $nav = $tmp_dir.'/repeats.nav';
my $list2nav = Genome::Model::Tools::Consed::ListToNav->create(
    list => $list,
    nav => $nav,
    title => 'Converted from list - repeats from /home1/watson/seqmgr/M_BB0392D19/edit_dir/M_BB0392D19.fasta.screen.ace.old',
);
ok($list2nav, 'Created list to nav');
ok($list2nav->execute, 'execute');
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
