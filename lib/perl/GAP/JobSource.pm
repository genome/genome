package GAP::JobSource;

use strict;
use warnings;


sub new {

    my ($class) = @_;

    
    my $self = { };
    bless $self, $class;
    
    $self->{_job_id}      = 0;
    $self->{_jobs}        = [ ];
    $self->{_fail_count}  = { };
    $self->{_failed_jobs} = [ ];
    
    return $self;
    
}

sub next_job_id {

    my ($self) = @_;


    my $next_job_id                      = $self->{_job_id}++;
    $self->{_fail_count}->{$next_job_id} = 0;
    
    return $next_job_id;
    
}

sub get_job {
    
    my ($self) = @_;
    
    
    return shift @{$self->{_jobs}};
    
}

sub fail_job {
    
    my ($self, $job, $results, $error) = @_;

    
    my $job_id = $job->job_id();
    
    if ($self->{_fail_count}->{$job_id} >= 3) {
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

sub failed_jobs {
    
    my ($self) = @_;

    
    return $self->{_failed_jobs};
    
}

1;
