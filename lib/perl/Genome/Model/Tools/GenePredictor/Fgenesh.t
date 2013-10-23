use strict;
use warnings;

use above 'Genome';

use Bio::Seq;
use Bio::SeqIO;

use File::Temp qw(tempdir);
use File::Basename;
use Test::More tests => 18;

BEGIN {
    use_ok('Genome::Model::Tools::GenePredictor');
    use_ok('Genome::Model::Tools::GenePredictor::Fgenesh');
}

my $test_output_dir = tempdir('Genome-Model-Tools-GenePredictor-Fgenesh-XXXXXX',
    TMPDIR => 1,
    CLEANUP => 1,
);
chmod(0755, $test_output_dir);

my $test_data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-GenePredictor/';
ok(-d $test_data_dir, "test data directory exists at $test_data_dir");

my $fasta = $test_data_dir . 'Contig0a.masked.fasta';
ok(-e $fasta, "fasta file exists at $fasta");

my $model = $ENV{GENOME_SW} . '/softberry/fgenesh_installed/C_elegans';
ok(-e $model, "model file exists at $model");

my $command = Genome::Model::Tools::GenePredictor::Fgenesh->create(
    fasta_file => $fasta,
    model_file => $model,
    raw_output_directory => $test_output_dir,
    prediction_directory => $test_output_dir,
);

isa_ok($command, 'Genome::Model::Tools::GenePredictor');
isa_ok($command, 'Genome::Model::Tools::GenePredictor::Fgenesh');

ok($command->execute(), "executed fgenesh command");

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

my $gene = shift @genes;
my $protein = $gene->protein;
ok(defined $protein, "able to grab protein from coding gene via indirect property");

my $transcript = $gene->transcript;
ok(defined $transcript, "able to grab transcript from coding gene via indirect property");

@exons = $transcript->exons;
ok(@exons, "able to grab exons from transcript via indirect property");

$protein = $transcript->protein;
ok(defined $protein, "able to grab protein from transcript via indirect property");

my $exon = shift @exons;
$transcript = $exon->transcript;
ok(defined $transcript, "able to grab transcript from exon via indirect property");

$gene = $exon->coding_gene;
ok(defined $gene, "able to grab coding gene from exon via indirect property");
