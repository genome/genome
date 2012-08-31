package Genome::Model::Tools::Velvet::ToReadFasta;

use strict;
use warnings;

use Genome;
use IO::File;
use Bio::Seq;
use Bio::SeqIO;
use Genome::Model::Tools::Pcap::Ace::Reader;

class Genome::Model::Tools::Velvet::ToReadFasta {
    is           => 'Command::V2',
    has          => [
        ace_file    => {
            is      => 'String', 
            doc     => 'ace file name with path',
        }
    ],
    has_optional => [
        out_file    => {
            is      => 'String', 
            doc     => 'output fasta file name with path, default: ./reads.fasta',
            default => './reads.fasta',
        },
    ],
};
        

sub help_brief {
    'This tool grabs read_ids and their sequences from velvet_converted acefile',
}


sub help_detail {
    return <<EOS
This tool is needed to make read fasta for making fake Phds/Scfs for consed.
EOS
}


sub execute {
    my $self    = shift;
    my $acefile = $self->ace_file;
    
    unless (-s $acefile) {
        $self->error_message("Acefile $acefile not existing");
        return;
    }
    
    my $io = Bio::SeqIO->new(-format => 'fasta', -file => '>'.$self->out_file);
    my $fh = IO::File->new($acefile) or die "can't open $acefile\n";

    my $reader = Genome::Model::Tools::Pcap::Ace::Reader->new($fh);

    my %orient;

    while (my $obj = $reader->next_object) {
        if ($obj->{type} eq 'read_position') {
            $orient{$obj->{read_name}} = $obj->{u_or_c};
        }
        if ($obj->{type} eq 'read') {
            my $seq = $obj->{sequence};
            $seq =~ s/\*//g;
            my $sq = Bio::Seq->new(-seq => $seq, -id => $obj->{name});
            my $uc = delete $orient{$obj->{name}};
            $sq = $sq->revcom if $uc eq 'C';
            $io->write_seq($sq);
        }
    }
    $fh->close;
    
    return 1;
}

1;

