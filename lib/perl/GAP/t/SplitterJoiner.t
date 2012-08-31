use strict;
use warnings;

#use Workflow;
use above 'GAP';

use Bio::Seq;
use Bio::SeqIO;

use File::Temp;
use File::Basename;
use Test::More tests => 15;

BEGIN {
    use_ok('GAP::Command');
    use_ok('GAP::Command::FastaSplitter');
    use_ok('GAP::Command::FastaJoiner');
}

my $splitter = GAP::Command::FastaSplitter->create(
                                                   'fasta_file' => File::Basename::dirname(__FILE__).'/data/C_elegans.ws184.fasta.bz2',
                                                   'chunk_size' => 1,
                                                  );
isa_ok($splitter, 'GAP::Command::FastaSplitter');

ok($splitter->execute(), 'execute splitter');

my $expected_file_count = 6;

my @files = @{$splitter->fasta_files()};

ok(scalar(@files) == $expected_file_count, 'file count');

foreach my $file (@files) {

    is(count_seq($file), 1, 'split file sequence count');
    
}

my $temp_fh = File::Temp->new();
my $temp_fn = $temp_fh->filename();

$temp_fh->close();

my $joiner = GAP::Command::FastaJoiner->create(
                                               'fasta_files' => \@files,
                                               'fasta_file'  => $temp_fn,
                                              );

ok($joiner->execute(), 'execute joiner');

is(count_seq($temp_fn), 6, 'join file sequence count');

ok(unlink(@files) == $expected_file_count, 'all files deleted');


sub count_seq {

    my ($seq_file) = (@_);


    my $seq_count = 0;

    my $seqin = Bio::SeqIO->new(-file => $seq_file, -format => 'Fasta');

    while (my $seq = $seqin->next_seq()) {
        $seq_count++;
    }

    return $seq_count;

}
