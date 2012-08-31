package GAP::JobSource::RNAmmer;

use strict;
use warnings;

use Carp;
use English;

use GAP::Job::RNAmmer;

use base qw(GAP::JobSource);


sub new {

    my ($class, @args) = @_;


    my $self = $class->SUPER::new(@args);
    
    my (
        $seq_source,
        $domain,
    ) = @args;
    
    unless (defined($seq_source)) {
        croak 'missing seq source!';
    }
    
    unless ($seq_source->can('next_seq')) {
        croak 'seq source does not implement next_seq()!';
    }
    
    unless (defined($domain)) {
        croak 'missing domain';
    }
    
    $self->{_feature_ref} = [ ]; 
    
    while (my $seq = $seq_source->next_seq()) {
        
        push @{$self->{_jobs}}, GAP::Job::RNAmmer->new(
                                                       $seq,
                                                       $domain,
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
