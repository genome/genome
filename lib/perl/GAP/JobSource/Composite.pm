package GAP::JobSource::Composite;

use strict;
use warnings;

use Carp;
use English;

use base qw(GAP::JobSource);


sub new {

    my ($class, @args) = @_;


    my $self = $class->SUPER::new(@args);
    
    my @job_sources = @args;
    
   
    $self->{_feature_ref} = [ ]; 
    $self->{_jobs}        = [ ];
    
    foreach my $job_source (@job_sources) {

        unless ($job_source->isa('GAP::JobSource')) {
            croak "all job sources must inherit from GAP::JobSource";
        }

        while (my $job = $job_source->get_job()) {

            $job->job_id($self->next_job_id());

            push @{$self->{_jobs}}, $job;
            
        }
        
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
