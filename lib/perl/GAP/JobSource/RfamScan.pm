package GAP::JobSource::RfamScan;

use strict;
use warnings;

use Carp;
use English;

use GAP::Job::RfamScan;

use base qw(GAP::JobSource);


sub new {

    my ($class, @args) = @_;


    my $self = $class->SUPER::new(@args);
    
    my (
        $seq_source,
    ) = @args;
    
    unless (defined($seq_source)) {
        croak 'missing seq source!';
    }
    
    unless ($seq_source->can('next_seq')) {
        croak 'seq source does not implement next_seq()!';
    }
    
    $self->{_feature_ref} = [ ]; 
    
    while (my $seq = $seq_source->next_seq()) {
        
        push @{$self->{_jobs}}, GAP::Job::RfamScan->new(
                                                        $seq,
                                                        $self->next_job_id(),
                                                    );
        
    }
    
    return $self;
    
}

sub finish_job {

    my ($self, $job, $results) = @_;


    push @{$self->{_feature_ref}}, $job->seq()->get_SeqFeatures();
    
}


sub feature_ref {
    
    my ($self) = @_;
    
    return $self->{_feature_ref};
    
}

1;
