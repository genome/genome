package BAP::JobSource::Genemark;

use strict;
use warnings;

use Carp;
use English;

use BAP::Job::Genemark;

use base qw(GAP::JobSource);


sub new {

    my ($class, @args) = @_;


    my $self = $class->SUPER::new(@args);
    
    my (
        $seq_source,
        $gm_model,
        $circular_dna_flag
    ) = @args;
    
    unless (defined($seq_source)) {
        croak 'missing seq source!';
    }
    
    unless ($seq_source->can('next_seq')) {
        croak 'seq source does not implement next_seq()!';
    }
    
    unless (defined($gm_model)) {
        croak 'missing genemark model';
    }
    
    unless (defined($circular_dna_flag)) {
        $circular_dna_flag = 0;
    }
    
    $self->{_feature_ref} = [ ]; 
    
    while (my $seq = $seq_source->next_seq()) {
        
        push @{$self->{_jobs}}, BAP::Job::Genemark->new(
                                                         $seq,
                                                         $gm_model,
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

sub fail_job {
    
    my ($self, $job, $results, $error) = @_;

    
    my $job_id = $job->job_id();
    
    if ($self->{_fail_count}->{$job_id} >= 10) {
        warn "job $job_id continues to fail after multiple retries: $error";
        push @{$self->{_failed_jobs}}, $job;
    }
    else {
        $self->{_fail_count}->{$job_id}++;
        my $job_host = $job->execution_host();
        warn "job $job_id failed on $job_host - retrying";
        push @{$self->{_jobs}}, $job;
    }
    
}

1;
