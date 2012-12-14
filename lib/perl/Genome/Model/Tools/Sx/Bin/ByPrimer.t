#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

require File::Compare;
use Test::More;

# Use
use_ok('Genome::Model::Tools::Sx::Bin::ByPrimer') or die;

# Create fails
ok(!Genome::Model::Tools::Sx::Bin::ByPrimer->execute(), 'create w/o primers');
ok(!Genome::Model::Tools::Sx::Bin::ByPrimer->execute(primers => [':seq']), 'create w/o name');
ok(!Genome::Model::Tools::Sx::Bin::ByPrimer->execute(primers => ['name:']), 'create w/o sequence');

# Files
my $dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Sx';
my $example_acgt = $dir.'/bin_by_primer.acgt.fastq';
ok(-s $example_acgt, 'example acgt');
my $example_tgca = $dir.'/bin_by_primer.tgca.fastq';
ok(-s $example_tgca, 'example tgca');
my $tmp_dir = File::Temp::tempdir(CLEANUP => 1);

# Ok
my $acgt_fastq = $tmp_dir.'/acgt.fastq';
my $tgca_fastq = $tmp_dir.'/tgca.fastq';
my $binner = Genome::Model::Tools::Sx::Bin::ByPrimer->create(
    input  => [ $dir.'/bin_by_primer.fastq' ],
    output => [ 'name=ACGT:'.$acgt_fastq, 'name=TGCA:'.$tgca_fastq ],
    primers => [ 'ACGT=ACGT', 'TGCA=TGCA' ],
);
ok($binner, 'create bin by primer');
$binner->dump_status_messages(1);
ok($binner->execute, 'execute binner');
is(File::Compare::compare($acgt_fastq, $dir.'/bin_by_primer.acgt.fastq'), 0, 'ACGT file matches');
is(File::Compare::compare($tgca_fastq, $dir.'/bin_by_primer.tgca.fastq'), 0, 'TGCA file matches');

#$acgt_fastq = $tmp_dir.'/acgt.rm.fastq';
#$tgca_fastq = $tmp_dir.'/tgca.rm.fastq';
#$binner = Genome::Model::Tools::Sx::Bin::ByPrimer->create(
#    input  => [ $dir.'/bin_by_primer.fwd.fastq', $dir.'/bin_by_primer.rev.fastq' ],
#    output => [ 'name=ACGT:'.$acgt_fastq, 'name=TGCA:'.$tgca_fastq ],
#    primers => [ 'ACGT=ACGT', 'TGCA=TGCA' ],
#    remove => 1,
#);
#ok($binner, 'create bin by primer');
#$binner->dump_status_messages(1);
#ok($binner->execute, 'execute binner');
#is(File::Compare::compare($acgt_fastq, $dir.'/bin_by_primer.acgt.rm.fastq'), 0, 'ACGT file mathces');
#is(File::Compare::compare($tgca_fastq, $dir.'/bin_by_primer.tgca.rm.fastq'), 0, 'TGCA file mathces');

my $pair1_fastq = $tmp_dir.'/PAIR1.rm.fastq';
my $pair2_fastq = $tmp_dir.'/PAIR2.rm.fastq';
$binner = Genome::Model::Tools::Sx::Bin::ByPrimer->create(
    input  => [ $dir.'/bin_by_primer.paired.fwd.fastq', $dir.'/bin_by_primer.paired.rev.fastq' ],
    output => [ 'name=PAIR1:'.$pair1_fastq, 'name=PAIR2:'.$pair2_fastq ],
    primers => [ 'PAIR1.F=ACGTA','PAIR1.R=ACGTT','PAIR2.F=TGCAC','PAIR2.R=TGCAG' ],
    remove => 1,
);
ok($binner, 'create bin by primer');
$binner->dump_status_messages(1);
ok($binner->execute, 'execute binner');
is(File::Compare::compare($pair1_fastq, $dir.'/bin_by_primer.PAIR1.rm.fastq'), 0, 'ACGT file matches');
is(File::Compare::compare($pair2_fastq, $dir.'/bin_by_primer.PAIR2.rm.fastq'), 0, 'TGCA file matches');

#print "$tmp_dir\n"; <STDIN>;

done_testing();
