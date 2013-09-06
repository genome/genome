use strict;
use warnings;

use above 'Genome';

use File::Temp 'tempdir';
use File::Basename;
use Test::More tests => 9;

BEGIN {
    use_ok('Genome::Model::Tools::GenePredictor');
    use_ok('Genome::Model::Tools::GenePredictor::Trnascan');
}

my $test_output_dir = tempdir('Genome-Model-Tools-GenePredictor-Trnascan-XXXXXX',
    TMPDIR => 1,
    CLEANUP => 1,
);
chmod(0755, $test_output_dir);
ok(-d $test_output_dir, "test output dir exists");

my $test_data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-GenePredictor';
ok(-d $test_data_dir, "test data directory exists at $test_data_dir");

my $fasta = $test_data_dir . '/Contig0a.masked.fasta';
ok(-e $fasta, "fasta file exists at $fasta");

my $command = Genome::Model::Tools::GenePredictor::Trnascan->create(
    fasta_file => $fasta, 
    raw_output_directory => $test_output_dir,
    prediction_directory => $test_output_dir,
);

isa_ok($command, 'Genome::Model::Tools::GenePredictor');
isa_ok($command, 'Genome::Model::Tools::GenePredictor::Trnascan');
ok($command->execute(), "executed trnascan command");

my @rna = Genome::Prediction::RNAGene->get(
    directory => $test_output_dir,
);
my $num_rna = scalar @rna;
ok($num_rna > 0, "able to retrieve $num_rna RNAGene objects");
