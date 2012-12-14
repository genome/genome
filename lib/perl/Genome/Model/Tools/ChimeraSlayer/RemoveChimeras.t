#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use Data::Dumper 'Dumper';
use File::Compare;
use File::Temp;
use Test::More;

use_ok('Genome::Model::Tools::ChimeraSlayer::RemoveChimeras') or die;

my $version = 1;
my $test_data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-ChimeraSlayer/RemoveChimeras/v'.$version;
my $sequences = $test_data_dir.'/sequences.fastq';
my $cpc_file = $test_data_dir.'/chimeras.CPS.CPC';
my $expected_output = $test_data_dir.'/sequences.chimera-free.fasta';

my $tmpdir = File::Temp::tempdir(CLEANUP => 1);
my $output = $tmpdir.'/out.fasta';

my $cmd = Genome::Model::Tools::ChimeraSlayer::RemoveChimeras->create(
    sequences => $sequences,
    chimeras => $cpc_file,
    output => $output,
);
ok($cmd, 'create');
$cmd->dump_status_messages(1);
    ok($cmd->execute, 'execute');
is(File::Compare::compare($output, $expected_output), 0, 'output mathces');
is_deeply($cmd->_metrics, { sequences => 7, chimeras => 3, output => 4, }, 'metrics match');

#print "$tmpdir\n"; <STDIN>;
done_testing();
