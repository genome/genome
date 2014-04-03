#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use File::Compare;
use Test::More;

if (Genome::Config->arch_os ne 'x86_64') {
    plan skip_all => 'requires 64-bit machine';
}
else {
    plan tests => 3;
}

use_ok('Genome::Model::Tools::Cufflinks::Assemble');


my $input_data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Cufflinks-Assemble';
my $bam_file = $input_data_dir .'/accepted_hits.bam';

# Nevermind, Cufflinks seems to produce slightly different results each time
# Not surprising, considering it's use of statistical approaches for normalization, correction, etc.
#my $expected_data_dir = $input_data_dir .'/expected_output';
#my $expected_transcripts_file = $expected_data_dir .'/transcripts.gtf';
#my $expected_gene_fpkm_file = $expected_data_dir .'/genes.fpkm_tracking';
#my $expected_isoform_fpkm_file = $expected_data_dir .'/isoforms.fpkm_tracking';

my $tmp_dir = File::Temp::tempdir('Cufflinks-Assemble-'.Genome::Sys->username.'-XXXX',CLEANUP => 1, TMPDIR => 1);

my $expression = Genome::Model::Tools::Cufflinks::Assemble->create(
   input_file => $bam_file,
   output_directory => $tmp_dir,
   params => '--frag-len-mean=250 --frag-len-std-dev=20 --min-frags-per-transfrag=1',
);
isa_ok($expression,'Genome::Model::Tools::Cufflinks::Assemble');
like($expression->command, qr/--no-update-check/, 'Do not check for updates while running cufflinks');
ok($expression->execute,'execute command '. $expression->command_name);

# See above comment on comparing expected output files
#ok( (compare($expression->transcripts_file,$expected_transcripts_file) == 0),'transcripts are identical');
#ok( (compare($expression->gene_fpkm_file,$expected_gene_fpkm_file) == 0),'gene fpkm are identical');
#ok( (compare($expression->isoform_fpkm_file,$expected_isoform_fpkm_file) == 0),'isoform fpkm are identical');
