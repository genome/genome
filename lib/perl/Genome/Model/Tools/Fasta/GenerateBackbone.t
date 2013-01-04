#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More tests => 4;

use File::Compare;

use above 'Genome';

BEGIN {
        use_ok('Genome::Model::Tools::Fasta::GenerateBackbone');
}
my $file_name = 'GENES.tsv';

my $tmp_dir = File::Temp::tempdir('Genome-Model-Tools-Fasta-GenerateBackbone-XXXXX', DIR => "$ENV{GENOME_TEST_TEMP}", CLEANUP => 1);
my $backbone_file = $tmp_dir .'/'. $file_name;

my $existing_data_dir = Genome::Config::reference_sequence_directory() . '/new_masked_ccds_ensembl_genbank_utr_nosv_all_transcriptome_quickfix';
my $fasta_file = $existing_data_dir .'/new_masked_ccds_ensembl_genbank_utr_nosv_all_transcriptome_quickfix.fa';
my $expected_backbone = $existing_data_dir .'/'. $file_name;

my $generate_backbone = Genome::Model::Tools::Fasta::GenerateBackbone->create(
                                                                              fasta_file => $fasta_file,
                                                                              output_file => $backbone_file,
                                                                          );
isa_ok($generate_backbone,'Genome::Model::Tools::Fasta::GenerateBackbone');

ok($generate_backbone->execute,'execute command '. $generate_backbone->command_name);

ok(!compare($backbone_file,$expected_backbone),'Backbone files are identical');
