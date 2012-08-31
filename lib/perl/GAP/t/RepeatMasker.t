use strict;
use warnings;

use above 'GAP';

use Bio::Seq;
use Bio::SeqIO;
use File::Temp;
use File::Basename;
use Test::More;

BEGIN {
    use_ok('GAP::Command');
    use_ok('GAP::Command::RepeatMasker');
}

my $command = GAP::Command::RepeatMasker->create(
#    fasta_file => File::Basename::dirname(__FILE__).'/data/1183.dna',
    fasta_file => File::Basename::dirname(__FILE__).'/data/5_ctg.dna',
    repeat_library => '/gsc/var/lib/repeat/Trichinella_spiralis-3.7.1_070320.rep', 
    masked_fasta => '/gscuser/ssurulir/1183.dna.repeat_masker',
	xsmall => 1,
	ace => 1,
);
isa_ok($command, 'GAP::Command::RepeatMasker');

ok($command->execute(), 'execute');

my $seqio = Bio::SeqIO->new(
    -file => $command->masked_fasta(), 
    -format => 'Fasta'
);

while (my $seq = $seqio->next_seq()) {
    isa_ok($seq, 'Bio::SeqI');
}

done_testing();
