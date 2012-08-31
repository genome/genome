package Genome::Model::Tools::Hgmi::CalculateGcPercent;

use strict;
use warnings;

use Genome;

use Bio::Seq;
use Bio::SeqIO;

class Genome::Model::Tools::Hgmi::CalculateGcPercent {
    is => 'Command',
    has => [
        fasta_files => { is => 'ARRAY', doc => 'array of fasta file names',
                         is_input => 1, },
        gc_percent => { is => 'Float', is_optional => 1, doc => 'GC content',
                        is_output => 1, }
    ],
};

sub help_brief {
    "Write a set of fasta files for an assembly";
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
    
    my @files = @{$self->fasta_files()};

    my $gc_count;
    my $seq_length;
    
    foreach my $file (@files) {

        my $seqio = Bio::SeqIO->new(-file => $file, -format => 'Fasta');
        
        while (my $seq = $seqio->next_seq()) {
            
            my $seq_string = $seq->seq();
            
            ## Ns are usually gaps...
            my $n_count = $seq_string =~ tr/nN/nN/;
            
            $gc_count   += $seq_string =~ tr/gcGC/gcGC/;
            
            ## ...so don't count them when determining the sequence length
            $seq_length += ($seq->length() - $n_count); 
            
        }

    }

    $self->gc_percent(sprintf("%.1f", (($gc_count / $seq_length) * 100)));

    if($self->gc_percent < 30)
    {
        $self->gc_percent(30.0);
    }
    elsif($self->gc_percent > 70)
    {
        $self->gc_percent(70.0);
    }

    return 1;
    
}
