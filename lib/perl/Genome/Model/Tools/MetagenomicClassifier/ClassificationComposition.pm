package Genome::Model::Tools::MetagenomicClassifier::ClassificationComposition;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::MetagenomicClassifier::ClassificationComposition {
    has => [
        confidence_threshold => {
            is => 'Number',
            doc => 'The confidence threshold that the classifications must be greter than or equal to to be considered confident.',
        },
    ],
};

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);
    return if not $self;

    my $threshold = $self->confidence_threshold;
    if ( not $threshold or $threshold > 1 or $threshold <= 0 ) {
        $self->error_message('Invalid confidence threshold: '.$threshold);
        return;
    }

    $self->{_classifications} = [];
    
    return $self;
}

sub add_classification {
    my ($self, $classification) = @_;

    Carp::confess('No root taxon found in classification: '.Data::Dumper::Dumper($classification)) if not $classification->{root}->{id};

    my $i = ( $classification->{root}->{confidence} >= $self->confidence_threshold ) ? 1 : 0;
    
    push @{$self->{_classifications}}, $classification;
    
    return 1;
}

sub confident_classifications {
    return $_[0]->{_classifications};
}

sub get_counts_for_domain_down_to_rank {
    my ($self, $domain, $to_rank) = @_;

    my $valid_domain = Genome::Model::Tools::MetagenomicClassifier->is_domain_valid($domain);
    return if not $valid_domain;

    my $valid_rank = Genome::Model::Tools::MetagenomicClassifier->is_rank_valid($to_rank);
    return if not $valid_rank;

    my @ranks;
    for my $rank ( Genome::Model::Tools::MetagenomicClassifier->taxonomic_ranks ) {
        push @ranks, $rank;
        last if $rank eq $to_rank;
    }

    my %counts;
    my $threshold = $self->confidence_threshold;
    for my $classification ( @{$self->confident_classifications} ) {
        next if lc($classification->{domain}->{id}) ne $domain;
        my $taxonomy = join(
            ':', 
            grep { defined } map { $classification->{$_}->{id} } @ranks
        );
        # Increment total
        $counts{$taxonomy}->{total}++;
        # Go thru the ranks
        for my $rank ( @ranks ) {
            my $confidence = $classification->{$rank}->{confidence}
                or next;
            if ( $confidence >= $threshold ) {
                $counts{$taxonomy}->{$rank}++;
            }
            elsif ( not defined $counts{$taxonomy}->{$rank} ) {
                $counts{$taxonomy}->{$rank} = 0;
            }
        }
    }

    return %counts;
}

1;

