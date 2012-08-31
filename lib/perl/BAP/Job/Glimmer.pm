package BAP::Job::Glimmer;

use strict;
use warnings;

use GAP::Job;
use Bio::Tools::Run::Glimmer;
use Carp;

use base qw(GAP::Job);


sub new {

    my ($class, $program, $seq, $icm, $pwm, $circular_flag, $job_id) = @_;

    
    my $self = { };
    bless $self, $class;

    unless (defined($job_id)) {
        croak 'missing job id';
    }
    
    $self->job_id($job_id);

    unless (defined($program)) {
        croak 'missing program!';
    }

    unless (($program eq 'glimmer2') || ($program eq 'glimmer3')) {
        croak "bad program '$program'!";
    }

    $self->{_program} = $program;
            
    unless (defined($seq)) {
        croak 'missing seq object!';
    }

    unless ($seq->isa('Bio::PrimarySeqI')) {
        croak 'seq object is not a Bio::PrimaySeqI!';
    }

    $self->{_seq} = $seq;
    
    unless (defined($icm)) {
        croak 'missing model/icm!';
    }

    $self->{_icm} = $icm;

    if ($program eq 'glimmer3') {

        unless (defined($pwm)) {
            croak 'missing pwm!';

        }

        $self->{_pwm} = $pwm;
        
    }

    if ($circular_flag) {
        $self->{_circular_dna} = 1;
    }
    else {
        $self->{_circular_dna} = 0;
    }
    
    return $self;
    
}

sub execute {
    
    my ($self) = @_;

    $self->SUPER::execute(@_);


    my $program = $self->{_program};
    my $seq     = $self->{_seq}; 
    my $model   = $self->{_icm};

    my %args = (
                -program => $program,
                -model   => $model,
            );

    if (exists($self->{_pwm})) {
        $args{'-b'} = $self->{_pwm};
    }

    unless ($self->{_circular_dna}) {
        $args{'-l'} = 1;
        $args{'-X'} = 1;
    }
   
    my $factory = Bio::Tools::Run::Glimmer->new(%args);
      
    my $parser = $factory->run($seq);
    
    while (my $gene = $parser->next_prediction()) {

        my $location = $gene->location();
        
        if ($location->isa('Bio::Location::SplitLocationI')) {
            
            foreach my $chunk ($location->sub_Location()) {
                
                my $new_location = Bio::Location::Fuzzy->new(
                                                             -start  => $chunk->start(),
                                                             -end    => $chunk->end(),
                                                             -strand => $chunk->strand(),
                                                            );
                
                my $new_gene = Bio::SeqFeature::Generic->new(
                                                             -seq_id      => $gene->seq_id(),
                                                             -location    => $new_location,
                                                             -strand      => $gene->strand(),
                                                             -source_tag  => $gene->source_tag(),
                                                             -primary_tag => $gene->primary_tag(),
                                                             -tag =>         { 'wraparound' => 1 },
                                                         );

                $seq->add_SeqFeature($gene);
                
            }
            
        }
        
        else {
            $seq->add_SeqFeature($gene);
        }
        
    }
    
}

1;
