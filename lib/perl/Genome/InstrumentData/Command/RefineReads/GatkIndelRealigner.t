#! /gsc/bin/perl

use strict;
use warnings;

use above 'Genome';

use Test::More;

# The running logic is tested inside GATK BP
my $class = 'Genome::InstrumentData::Command::RefineReads::GatkIndelRealigner';
use_ok($class) or die;
is_deeply([$class->result_names], ['indel realigner'], 'result_names');

done_testing();
