#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

require File::Compare;
use Test::More;

# Use
use_ok('Genome::Model::Tools::Sx::RmDesc') or die;

# Files
my $dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Sx';
my $example_fastq = $dir.'/rm_desc.example.fastq';
ok(-s $example_fastq, 'example fastq');
my $in_fastq = $dir.'/rm_desc.in.fastq';
ok(-s $example_fastq, 'in fastq');
my $tmp_dir = File::Temp::tempdir(CLEANUP => 1);
my $out_fastq = $tmp_dir.'/out.fastq';

# Ok
my $binner = Genome::Model::Tools::Sx::RmDesc->create(
    input  => [ $in_fastq ],
    output => [ $out_fastq ],
);
ok($binner, 'create rm desc');
$binner->dump_status_messages(1);
ok($binner->execute, 'execute rm desc');
is(File::Compare::compare($out_fastq, $example_fastq), 0, 'ACGT file mathces');

#print "$tmp_dir\n"; <STDIN>;
done_testing();
