package GAP::Job;

use strict;
use warnings;


sub job_id {

    my $self   = shift @_;
    my $job_id = shift @_;

    if (defined($job_id)) {
        $self->{_job_id} = $job_id; 
    }
    
    return $self->{_job_id};

}

sub seq {

    my $self = shift @_;

    return $self->{_seq};

}

sub execute {

    my $self = shift @_;

    $self->_determine_execution_host();

}

sub execution_host {

    my $self = shift @_;

    return $self->{_execution_host};

}

sub _determine_execution_host {

    my $self = shift @_;

    my $execution_host = 'unknown';

    if (exists($ENV{LSB_HOSTS})) { 
        $execution_host = $ENV{LSB_HOSTS};
    }
    elsif (exists($ENV{HOST})) {
        $execution_host = $ENV{HOST};
    }

    $self->{_execution_host} = $execution_host;

}

1;
