use strict;
use warnings;

use above 'Genome';

use Bio::Seq;
use Bio::SeqIO;
use File::Temp;
use File::Basename;
use Test::More 'skip_all'; 

BEGIN {
    use_ok('Genome::Model::Tools::GenePredictor');
    use_ok('Genome::Model::Tools::GenePredictor::Glimmer3');
}

my $test_data_dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-GenePredictor';
ok(-d $test_data_dir, "test data directory exists at $test_data_dir");

my $fasta_file = $test_data_dir . '/Contig0a.20kb.masked.fasta';
ok(-e $fasta_file, "test fasta file exists at $fasta_file");
ok(-s $fasta_file, "test fasta file has size at $fasta_file");

my $model_file = $test_data_dir . "/HPAG1.glimmer3.icm";
ok(-e $model_file, "model file exists at $model_file");

my $pwm_file = $test_data_dir . "/HPAG1.glimmer3.pwm";
ok(-e $pwm_file, "glimmer pwm file exists at $pwm_file");

my $command = Genome::Model::Tools::GenePredictor::Glimmer3->create(
    fasta_file => $fasta_file,
    model_file => $model_file,
    pwm_file => $pwm_file,
    raw_output_directory => '/tmp',
    prediction_directory => '/tmp',
);

isa_ok($command, 'Genome::Model::Tools::GenePredictor::Glimmer3');
ok($command->execute(), 'executed command');

my @features = @{$command->{bio_seq_feature}};
ok(@features > 0, 'produced bio seq features');

foreach my $feature (@features) {
    isa_ok($feature, 'Bio::SeqFeature::Generic');
}

done_testing();
