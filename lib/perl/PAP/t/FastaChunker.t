use strict;
use warnings;

use above "PAP";
use Workflow;

use Bio::Seq;
use Bio::SeqIO;

use File::Temp;
use File::Basename;
use Test::More;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
    $ENV{UR_COMMAND_DUMP_STATUS_MESSAGES} = 1;
}

use_ok('PAP::Command') or die;
use_ok('PAP::Command::FastaChunker') or die;

my $command = PAP::Command::FastaChunker->create(
    'fasta_file' => File::Basename::dirname(__FILE__) . '/data/large_test.fasta',
    'chunk_size' => 1000,
);
isa_ok($command, 'PAP::Command::FastaChunker') or die;

ok($command->execute(), 'successfully executed fasta chunker');

my @files = @{$command->fasta_files()};

ok(scalar(@files) == 7, 'found expected number of fasta chunks');

my $seq_count;
for my $file (@files) {
    $seq_count = 0;

    my $seqin = Bio::SeqIO->new(-file => $file, -format => 'Fasta');
    while (my $seq = $seqin->next_seq()) {
        $seq_count++;
    }

    ok($seq_count <= 1000, "sequence count for file $file is less than 1000");
}

ok(unlink(@files) == 7, 'removed fasta chunks successfully');
done_testing();

