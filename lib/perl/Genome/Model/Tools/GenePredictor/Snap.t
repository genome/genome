use strict;
use warnings;

use above 'Genome';

use Bio::Seq;
use Bio::SeqIO;

use File::Temp qw(tempdir);
use File::Basename;
use Test::More tests => 12;

BEGIN {
    use_ok('Genome::Model::Tools::GenePredictor');
    use_ok('Genome::Model::Tools::GenePredictor::Snap');
}

my $test_output_dir = tempdir('Genome-Model-Tools-GenePredictor-Snap-XXXXXX',
    TMPDIR => 1,
    CLEANUP => 1,
);
chmod(0755, $test_output_dir);

my $test_data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-GenePredictor/';
ok(-d $test_data_dir, "test data directory exists at $test_data_dir");

my $fasta = $test_data_dir . 'Contig0a.masked.fasta';
ok(-e $fasta, "fasta file exists at $fasta");

my $model = $ENV{GENOME_SW} . '/snap/installed/HMM/C.elegans.hmm';
ok(-e $model, "model file exists at $model");

my $command = Genome::Model::Tools::GenePredictor::Snap->create(
    fasta_file => $fasta,
    model_files => $model, 
    version => '2010-07-28',
    raw_output_directory => $test_output_dir,
    prediction_directory => $test_output_dir,
);

isa_ok($command, 'Genome::Model::Tools::GenePredictor');
isa_ok($command, 'Genome::Model::Tools::GenePredictor::Snap');

ok($command->execute(), "executed snap command");

my @genes = Genome::Prediction::CodingGene->get(
    directory => $test_output_dir,
);
my $num_genes = scalar @genes;
ok($num_genes > 0, "able to retrieve $num_genes coding gene objects");

my @proteins = Genome::Prediction::Protein->get(
    directory => $test_output_dir,
);
my $num_proteins = scalar @proteins;
ok($num_proteins > 0, "able to retrieve $num_proteins protein objects");

my @transcripts = Genome::Prediction::Transcript->get(
    directory => $test_output_dir,
);
my $num_transcripts = scalar @transcripts;
ok($num_transcripts > 0, "able to retrieve $num_transcripts transcript objects");

my @exons = Genome::Prediction::Exon->get(
    directory => $test_output_dir,
);
my $num_exons = scalar @exons;
ok($num_exons > 0, "able to retrieve $num_exons exon objects");
