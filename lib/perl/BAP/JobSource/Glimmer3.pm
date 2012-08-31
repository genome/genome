package BAP::JobSource::Glimmer3;

use strict;
use warnings;

use Carp;
use English;

use BAP::Job::Glimmer;

use base qw(GAP::JobSource);


sub new {

    my ($class, @args) = @_;


    my $self = $class->SUPER::new(@args);
    
    my (
        $seq_source,
        $g3_model,
        $g3_pwm,
        $circular_dna_flag
    ) = @args;
    
    unless (defined($seq_source)) {
        croak 'missing seq source!';
    }
    
    unless ($seq_source->can('next_seq')) {
        croak 'seq source does not implement next_seq()!';
    }
    
    unless (defined($g3_model)) {
        croak 'missing glimmer3 model';
    }

     unless (defined($g3_pwm)) {
         croak 'missing glimmer3 pwm';
    }
    
    unless (defined($circular_dna_flag)) {
        $circular_dna_flag = 0;
    }
    
    $self->{_feature_ref} = [ ]; 
    
    while (my $seq = $seq_source->next_seq()) {
        
        push @{$self->{_jobs}}, BAP::Job::Glimmer->new(
                                                       'glimmer3',
                                                       $seq,
                                                       $g3_model,
                                                       $g3_pwm,
                                                       $circular_dna_flag,
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
