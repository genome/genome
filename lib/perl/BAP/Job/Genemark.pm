package BAP::Job::Genemark;

use strict;
use warnings;

use GAP::Job;
use Bio::Tools::Run::Genemark;
use Carp;

use base qw(GAP::Job);


sub new {

    my ($class, $seq, $model, $job_id) = @_;

    
    my $self = { };
    bless $self, $class;

    unless (defined($job_id)) {
        croak 'missing job id';
    }
    
    $self->job_id($job_id);
    
    unless (defined($seq)) {
        croak 'missing seq object!';
    }
    
    unless ($seq->isa('Bio::PrimarySeqI')) {
        croak 'seq object is not a Bio::PrimaySeqI!';
    }

    $self->{_seq} = $seq;
    
    unless (defined($model)) {
        croak 'missing model!';
    }

    $self->{_model} = $model;

    return $self;
    
}

sub execute {
    
    my ($self) = @_;

    $self->SUPER::execute(@_);


    my $seq     = $self->{_seq}; 
    my $model   = $self->{_model};
    
    my $factory = Bio::Tools::Run::Genemark->new(
                                                 '-program' => 'gmhmmp', 
                                                 '-m'       => $model,
                                             );
        
    my $parser = $factory->run($seq);
    
    while (my $gene = $parser->next_prediction()) {
        $seq->add_SeqFeature($gene);
    }            
    
}

1;
