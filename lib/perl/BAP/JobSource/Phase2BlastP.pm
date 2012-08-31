package BAP::JobSource::Phase2BlastP;

use strict;
use warnings;

use BAP::Job::Phase2BlastP;

use Bio::SeqIO;

use Carp;
use English;

use base qw(GAP::JobSource);


sub new {

    my ($class, @args) = @_; 


    my $self = $class->SUPER::new(@args);
    
    my $db = shift @args;
    unless (defined($db)) {
        croak 'missing db!';
    }

    my $fasta_file = shift @args;
    unless (defined($fasta_file)) {
        croak 'missing fasta file!';
    }
   
       
    my $core_num = shift @args;
    unless (defined($core_num)) {
        croak 'missing number of cores to run blast in JobSource!';
    }

    $self->{_evidence} = { };
    
    my $seqio = Bio::SeqIO->new(-file => $fasta_file, -format => 'Fasta');
    
    while (my $seq = $seqio->next_seq()) {
        push @{$self->{_jobs}}, BAP::Job::Phase2BlastP->new(
                                                            $seq, 
                                                            $db,
                                                            $self->next_job_id(),
                                                            $core_num,
                                                        );
    }
    
    return $self;
    
}

sub finish_job {

    my ($self, $job, $results) = @_;


    my $ref = $job->evidence();
    
    foreach my $gene (keys %{$ref}) {
        $self->{_evidence}->{$gene} = 1;
    }
    
}

sub evidence {

    my ($self) = @_;

    
    return $self->{_evidence};
    
}

1;
