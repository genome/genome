package BAP::JobSource::InterGenicBlastX;

use strict;
use warnings;

use BAP::Job::InterGenicBlastX;

use Bio::SeqFeature::Collection;
use Carp;
use English;

use base qw(GAP::JobSource);


sub new {

    my ($class, @args) = @_;


    my $self = $class->SUPER::new(@args);
    
    my (
        $sequence_source,
        $feature_source,
        $blastx_db,
	    $core_num,
    ) = @args;
    
    unless (defined($sequence_source)) {
        croak 'missing sequence source!';
    }
    
    unless ($sequence_source->can('next_seq')) {
        croak 'sequence source does not implement next_seq()!';
    }
    
    unless (defined($feature_source)) {
        croak 'missing feature source!';
    }
    
    unless ($feature_source->can('next_feature')) {
        croak 'seq source does not implement next_seq()!';
    }

    unless (defined($blastx_db)) {
        croak 'missing blastx db!';
    }
    
    unless (defined($core_num)) {
	croak 'missing number of cores to run blast in JobSource!';
    }

    $self->{_feature_collection} = Bio::SeqFeature::Collection->new();
    
    my %features = ( );
    
    while (my $feat = $feature_source->next_feature) {

        my $seq_id = $feat->seq_id();

        unless (defined($seq_id)) {
            carp 'skipping feature with undef seq_id';
            next;
        }
        
        push @{$features{$seq_id}}, $feat;
        
    }
    
    while (my $seq = $sequence_source->next_seq()) {
        
        my $mask_ref = [ ];
        
        if (exists($features{$seq->display_id()})) {
            push @{$mask_ref}, @{$features{$seq->display_id}};
        }
        
        push @{$self->{_jobs}}, BAP::Job::InterGenicBlastX->new(
            $self->next_job_id(),
            $seq,
            $mask_ref,
            $blastx_db,
			$core_num,
        );
        
    }
    
    return $self;
    
}

sub finish_job {

    my ($self, $job, $results) = @_;

    
    $self->{_feature_collection}->add_features($job->genes());
    
}

sub feature_collection {

    my ($self) = @_;

    
    return $self->{_feature_collection};
    
}

1;
