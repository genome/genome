use strict;
use warnings;

use above "Genome";
use Test::More;

use Bio::Seq;

use Bio::SeqIO;

use Cwd;
use File::Temp;
use File::Basename;

unless ($ENV{UR_RUN_LONG_TESTS}) {
    plan skip_all => 'test takes twenty minutes to complete, run with UR_RUN_LONG_TESTS';
}

$ENV{UR_COMMAND_DUMP_STATUS_MESSAGES} = 1;

use_ok('Genome::Model::Tools::Predictor::Keggscan') or die;

my $test_data_base_dir = $ENV{GENOME_TEST_INPUTS} . '/';
my $fasta = join('/', $test_data_base_dir, 'Genome-Model-Tools-Predictor/medium.fasta');
ok(-e $fasta, "test fasta file found at $fasta") or die;

my $expected_dump_file = join('/', $test_data_base_dir, 'Genome-Model-Tools-Predictor/kegg.medium_fasta.dump.expected');
ok(-e $expected_dump_file, "expected dump file exists at $expected_dump_file") or die;

my $test_output_dir = File::Temp::tempdir(
    'Genome-Model-Tools-Predictor-XXXXX',
    TMPDIR => 1,
    CLEANUP => 1,
);
ok(-d $test_output_dir, "created output directory at $test_output_dir");

my $command = Genome::Model::Tools::Predictor::Keggscan->create(
    output_directory => $test_output_dir,
    input_fasta_file => $fasta,
    version => '56',
    parameters => '',
    dump_predictions_to_file => 1,
);
ok($command, 'created kegg command object');
ok($command->execute, 'executed kegg command');

ok(-e $command->raw_output_path, 'raw output file exists at ' . $command->raw_output_path);
ok(-e $command->dump_output_path, 'dump file exists at ' . $command->dump_output_path);

my $dump_diff = Genome::Sys->diff_file_vs_file($expected_dump_file, $command->dump_output_path);
ok(!$dump_diff, "no differences found between generated dump file " . $command->dump_output_path .
    " and expected dump file $expected_dump_file");

ok(-e $command->ace_file_path, 'ace file exists at ' . $command->ace_file_path);

done_testing();
