package GAP::Job::tRNAscan;

use strict;
use warnings;

use GAP::Job;
use Bio::Tools::Run::tRNAscanSE;
use Carp;

use base qw(GAP::Job);


sub new {

    my ($class, $seq, $domain, $job_id) = @_;

    
    my $self = { };
    bless $self, $class;
    
    unless (defined($seq)) {
        croak 'missing seq object';
    }
    
    unless ($seq->isa('Bio::PrimarySeqI')) {
        croak 'seq object is not a Bio::PrimaySeqI!';
    }
    
    $self->{_seq} = $seq;
    
    unless (defined($domain)) {
        croak 'missing domain';
    }

    unless (
            ($domain eq 'archaea'  )
            ||
            ($domain eq 'bacteria' )
            ||
            ($domain eq 'eukaryota')
        ) {
        croak "invalid domain '$domain'"; 
    }
    
    $self->{_domain} = $domain;
    
    unless (defined($job_id)) {
        croak 'missing job id';
    }
    
    $self->job_id($job_id);
    
    
    return $self;
    
}

sub execute {
    
    my ($self) = @_;

    $self->SUPER::execute(@_);


    my $seq    = $self->{_seq};
    my $domain = $self->{_domain};
    
    my %switches = ( );

    if    ($domain eq 'archaea' ) { $switches{-A} = 1; }
    elsif ($domain eq 'bacteria') { $switches{-B} = 1; }
    
    my $factory = Bio::Tools::Run::tRNAscanSE->new(
                                                   '-program' => 'tRNAscan-SE',
                                                   %switches,
                                               );
    
    
    my $parser = $factory->run($seq);
    
    while (my $gene = $parser->next_prediction()) {
        $seq->add_SeqFeature($gene);
    }            
    
}

1;
