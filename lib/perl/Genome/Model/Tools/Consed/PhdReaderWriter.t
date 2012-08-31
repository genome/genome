#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Data::Dumper;
use File::Compare;
use File::Temp;
use Test::More;

use_ok('Genome::Model::Tools::Consed::PhdReader') or die;
use_ok('Genome::Model::Tools::Consed::PhdWriter') or die;

my $phd_file = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Consed/L25990P6001H1.b1.phd.8';
my $tmpdir = File::Temp::tempdir(CLEANUP => 1);
my $phd_out = $tmpdir.'/phd';

my $phd = Genome::Model::Tools::Consed::PhdReader->read($phd_file);
ok($phd, 'read phd');
ok(Genome::Model::Tools::Consed::PhdWriter->write($phd_out, $phd), 'write phd');
is(File::Compare::compare($phd_out, $phd_file), 0, 'phd file matches');

#print "gvimdiff $phd_out $phd_file\n"; <STDIN>;
done_testing();
exit;

