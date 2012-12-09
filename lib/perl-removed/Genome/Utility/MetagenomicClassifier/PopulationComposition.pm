package Genome::Utility::MetagenomicClassifier::PopulationComposition;

use strict;
use warnings;

use Carp 'confess';
use Data::Dumper 'Dumper';
use Genome::Utility::MetagenomicClassifier;
use Regexp::Common;

sub new {
    my ($class, %params) = @_;

    my $self = bless {}, $class;

    if ( my $threshold = delete $params{confidence_threshold} ) {
        $self->_fatal_message(
            "Invalid confidence_threshold ($threshold) sent to 'new'"
        ) unless $threshold =~ /^$RE{num}{real}$/;
        $self->{confidence_threshold} = $threshold;
    }
    else {
        $self->{confidence_threshold} = 0.8;
    }

    $self->_fatal_message(
        "Unknown params sent to 'new': ".join(',', map { $_.' => '.$params{$_} } keys %params)
    ) if %params;
        
    $self->{_classifications} = [ [], [] ];
    
    # track stats?? would need to track when a classification is added, too
    
    return $self;
}

sub _fatal_message {
    my ($self, $msg) = @_;

    confess ref($self)." ERROR: $msg\n";
}
 
sub get_confidence_threshold {
    return $_[0]->{confidence_threshold};
}

sub add_classification {
    my ($self, $classification) = @_;

    my $i = ( ($classification->get_root_taxon->get_tag_values('confidence'))[0] >= $self->get_confidence_threshold ) 
    ? 1
    : 0;
    
    push @{$self->{_classifications}->[$i]}, $classification;
    
    return 1;
}

sub get_classifications {
    return map { @$_ } @{$_[0]->{_classifications}};
}

sub get_confident_classifications {
    return @{$_[0]->{_classifications}->[1]};
}

sub get_unconfident_classifications {
    return @{$_[0]->{_classifications}->[0]};
}

sub get_counts_for_domain_down_to_rank {
    my ($self, $domain, $to_rank) = @_;

    Genome::Utility::MetagenomicClassifier->validate_domain($domain)
        or return;

    Genome::Utility::MetagenomicClassifier->validate_taxonomic_rank($to_rank)
        or return;

    my $threshold = $self->get_confidence_threshold;

    my @ranks;
    for my $rank ( Genome::Utility::MetagenomicClassifier->taxonomic_ranks ) {
        push @ranks, $rank;
        last if $rank eq $to_rank;
    }

    my %counts;
    for my $classification ( $self->get_confident_classifications ) {
        next unless $classification->get_domain =~ /^$domain$/i;
        my $taxonomy = join(
            ':', 
            grep { defined } map { $self->_get_name_from_classification_for_rank($classification, $_) } @ranks
        );
        # Increment total
        $counts{$taxonomy}->{total}++;
        # Go thru the ranks
        for my $rank ( @ranks ) {
            my $confidence = $self->_get_confidence_from_classification_for_rank($classification, $rank)
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

sub _get_name_from_classification_for_rank {
    my ($self, $classification, $rank) = @_;

    my $method = 'get_'.$rank;

    return $classification->$method;
}

sub _get_confidence_from_classification_for_rank {
    my ($self, $classification, $rank) = @_;

    my $method = 'get_'.$rank.'_confidence';

    return $classification->$method;
}

sub _get_name_and_confidence_from_classification_for_rank {
    my ($self, $classification, $rank) = @_;

    my $method = 'get_'.$rank.'_name_and_confidence';

    return $classification->$method;
}

1;

=pod

=head1 Name

ModuleTemplate

=head1 Synopsis

=head1 Usage

=head1 Methods

=head2 

=over

=item I<Synopsis>

=item I<Arguments>

=item I<Returns>

=back

=head1 See Also

=head1 Disclaimer

Copyright (C) 2005 - 2008 Washington University Genome Sequencing Center

This module is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

=head1 Author(s)

B<Eddie Belter> I<ebelter@watson.wustl.edu>

=cut

#$HeadURL$
#$Id$

