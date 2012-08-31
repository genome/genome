package GAP::Command::GenePredictor::tRNAscan;

use strict;
use warnings;

use Bio::SeqIO;
use GAP::Job::tRNAscan;

class GAP::Command::GenePredictor::tRNAscan {
    is  => ['GAP::Command::GenePredictor'],
    has => [
            domain => {
                        is  => 'TEXT',
                        doc => 'archaea/bacteria/eukaryota',
                        is_input => 1
                      },
           ],
};

sub sub_command_sort_position { 10 }

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

    my $seqio = Bio::SeqIO->new(-file => $self->fasta_file(), -format => 'Fasta');

    my $seq = $seqio->next_seq();
    
    ##FIXME: The last arg is the job_id, which is hardcoded here in 
    ##       a rather lame fashion.
    my $legacy_job = GAP::Job::tRNAscan->new(
                                             $seq,
                                             $self->domain(),
                                             2112,
                                         );
    
    $legacy_job->execute();
    
    my @features = $legacy_job->seq()->get_SeqFeatures();
    
    $self->bio_seq_feature(\@features);
           
    return 1;
    
}

1;
