package GAP::Command::FastaJoiner;

use strict;
use warnings;

use Workflow;

use Bio::Seq;
use Bio::SeqIO;

use English;
use File::Temp;


class GAP::Command::FastaJoiner {
    is  => ['GAP::Command'],
    has => [
        fasta_files => { is => 'ARRAY',  doc => 'array of fasta file names' },
        fasta_file  => { is => 'SCALAR', doc => 'fasta file name'           },
    ],
};

operation GAP::Command::FastaJoiner {
    input  => [ 'fasta_file', 'fasta_files' ],
    output => [ ],
};

sub sub_command_sort_position { 10 }

sub help_brief {
    "Joint fasta files into one fasta file";
}

sub help_synopsis {
    return <<"EOS"
EOS
}

sub help_detail {
    return <<"EOS"
Need documenation here.
EOS
}

sub execute {

    my $self = shift;


    my $fasta_file = $self->fasta_file();

    my $seqout = Bio::SeqIO->new(-file => ">$fasta_file", -format => 'Fasta');
    
    foreach my $file (@{$self->fasta_files()}) {

        my $seqin = Bio::SeqIO->new(-file => $file, -format => 'Fasta');

        while (my $seq = $seqin->next_seq()) {
            $seqout->write_seq($seq);
        }

    }

    return 1;

}
 
1;
