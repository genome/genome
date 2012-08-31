package Genome::Utility::MetagenomicClassifier::ChimeraClassification;

use strict;
use warnings;

use Carp 'confess';
use Data::Dumper 'Dumper';

sub new {
    my ($class, %params) = @_;


    for my $req (qw/ name classification probe_classifications /) {
        _fatal_message("Required parameter ($req) not found.") unless $params{$req};
    }

    my $self =  bless \%params, $class;

    $self->{complemented} = 0 unless $self->{complemented};

    return $self;
}

sub name {
    my $self = shift;
    return $self->{name};
}
sub _fatal_message {
    my ($msg) = @_;

    confess __PACKAGE__." ERROR: $msg\n";
}

sub get_profile {
    my $self = shift;
    my @profile;
    push @profile, $self->name;
    push @profile, $self->maximum_common_depth;
    push @profile, $self->divergent_genera_count;
    push @profile, $self->divergent_probe_percent;
    push @profile, $self->classification->get_genus_confidence;
    push @profile, $self->maximum_divergent_confidence_difference;
    push @profile, $self->minimum_convergent_confidence_difference;
    return @profile;
}

sub probe_count {
    my $self = shift;
    my @probe_classifications = @{$self->probe_classifications};
    return scalar @probe_classifications;
}

sub probe_classifications {
    my $self = shift;
    return $self->{probe_classifications};
}

sub classification {
    my $self = shift;
    return $self->{classification};
}

sub maximum_divergent_confidence_difference {
    my $self = shift;
    my $classification = $self->classification;

    my $classification_confidence = $classification->get_genus_confidence;
    my $classification_genus = $classification->get_genus;

    my $max_divergent_confidence_diff = -1;
    
    foreach my $probe_class (@{$self->probe_classifications}) {
        my $probe_genus = $probe_class->get_genus;
        if ($probe_genus ne $classification_genus) {
            my $probe_confidence = $probe_class->get_genus_confidence;
            my $confidence_diff = $probe_confidence - $classification_confidence;
            $max_divergent_confidence_diff = $confidence_diff if $confidence_diff > $max_divergent_confidence_diff;
        }
    }
    return $max_divergent_confidence_diff;
}

sub minimum_convergent_confidence_difference {
    my $self = shift;
    my $classification = $self->classification;

    my $classification_confidence = $classification->get_genus_confidence;
    my $classification_genus = $classification->get_genus;

    my $min_convergent_confidence_diff = 1;
    
    foreach my $probe_class (@{$self->probe_classifications}) {
        my $probe_genus = $probe_class->get_genus;
        if ($probe_genus eq $classification_genus) {
            my $probe_confidence = $probe_class->get_genus_confidence;
            my $confidence_diff = $probe_confidence - $classification_confidence;
            $min_convergent_confidence_diff = $confidence_diff if $confidence_diff < $min_convergent_confidence_diff;
        }
    }
    return $min_convergent_confidence_diff;
}

sub convergent_probe_percent {
    my $self = shift;
    return 1 - $self->divergent_probe_percent;
}
sub divergent_probe_percent {
    my $self = shift;
    return $self->divergent_probe_count / $self->probe_count;
}

sub divergent_probe_count {
    my $self = shift;
    my $classification = $self->classification;
    my $classification_genus = $classification->get_genus;

    my $divergent_count = 0;
    
    foreach my $probe_class (@{$self->probe_classifications}) {
        my $probe_genus = $probe_class->get_genus;
        if ($probe_genus ne $classification_genus) {
            $divergent_count++;
        }
    }
    return $divergent_count;
}

sub divergent_classifications {
    my $self = shift;

    my $classification = $self->classification;
    my $classification_genus = $classification->get_genus;

    my %divergent_taxons;
    $divergent_taxons{$classification_genus} = $classification;

    foreach my $probe_class (@{$self->probe_classifications}) {
        my $probe_genus = $probe_class->get_genus;
        $divergent_taxons{$probe_genus} = $probe_class;
    }
    return values %divergent_taxons; 
}

sub divergent_genera_count {
    my $self = shift;
    my $classification = $self->classification;
    my $classification_genus = $classification->get_genus;

    my %seen_genera;
    $seen_genera{$classification_genus} = 1;

    foreach my $probe_class (@{$self->probe_classifications}) {
        my $probe_genus = $probe_class->get_genus;
        $seen_genera{$probe_genus}++;
    }

    return (scalar keys %seen_genera) - 1;
}

sub maximum_common_depth {
    my $self = shift;
    my $classification = $self->classification;
    my $max_depth = 10;
    my @probe_classifications = @{$self->probe_classifications};

    foreach my $probe_class (@probe_classifications) {
        my $current_depth = 0;
        my $classification_taxon = $classification->get_taxon;
        my $probe_taxon = $probe_class->get_taxon;
        while (defined $probe_taxon && $probe_taxon->id eq $classification_taxon->id) {
            $current_depth++;
            ($probe_taxon) = $probe_taxon->get_Descendents;
            ($classification_taxon) = $classification_taxon->get_Descendents;
        }

        if ($current_depth - 1 < $max_depth) {
            $max_depth = $current_depth - 1;
        }
    }
    return $max_depth;
}
 
#< NAME >#
sub get_name {
    return $_[0]->{name};
}

#< COMPLEMENTED >#
sub get_complemented { 
    return $_[0]->{complemented};
}

sub is_complemented { 
    return $_[0]->{complemented};
}

#< CLASSIFIER TYPE >#
sub get_classifier {
    return $_[0]->{classifier};
}

#< TAXON >#
sub get_taxon {
    return $_[0]->{taxon};
}

sub get_probe_taxons {
    return $_[0]->{taxon};
}

sub get_taxa {
    my $self = shift;

    my @taxa;
    unless ( $self->{taxa} ) {
        my $taxon = $self->get_taxon;
        do { 
            push @taxa, $taxon;
            ($taxon) = $taxon->get_Descendents;
        } until not defined $taxon;
        $self->{taxa} = \@taxa;
    }

    return @{$self->{taxa}};
}

sub taxa_count { 
    return scalar($_[0]->get_taxa);
}

sub _get_taxon_for_rank {
    return (grep { $_->rank eq $_[1] } $_[0]->get_taxa)[0];
}

sub _get_taxon_name_for_rank {
    my $taxon = $_[0]->_get_taxon_for_rank($_[1])
        or return 'none';
    return $taxon->id;
}

sub get_root_taxon {
    return $_[0]->_get_taxon_for_rank('root');
}

sub get_root {
    return $_[0]->_get_taxon_name_for_rank('root');
}

sub get_domain_taxon {
    return $_[0]->_get_taxon_for_rank('domain');
}

sub get_domain {
    return $_[0]->_get_taxon_name_for_rank('domain');
}

sub get_domain_confidence {
    my $self = shift;
    my $domain_taxon = $self->get_domain_taxon;
    unless ($domain_taxon) {
        return;
    }
    my ($confidence) = $domain_taxon->get_tag_values('confidence');
    return $confidence;
}

sub get_kingdom_taxon {
    return $_[0]->_get_taxon_for_rank('kingdom');
}

sub get_kingdom {
    return $_[0]->_get_taxon_name_for_rank('kingdom');
}

sub get_kingdom_confidence {
    my $self = shift;
    my $kingdom_taxon = $self->get_kingdom_taxon;
    unless ($kingdom_taxon) {
        return;
    }
    my ($confidence) = $kingdom_taxon->get_tag_values('confidence');
    return $confidence;
}

sub get_phylum_taxon {
    return $_[0]->_get_taxon_for_rank('phylum');
}

sub get_phylum {
    return $_[0]->_get_taxon_name_for_rank('phylum');
}

sub get_phylum_confidence {
    my $self = shift;
    my $phylum_taxon = $self->get_phylum_taxon;
    unless ($phylum_taxon) {
        return;
    }
    my ($confidence) = $phylum_taxon->get_tag_values('confidence');
    return $confidence;
}

sub get_class_taxon {
    return $_[0]->_get_taxon_for_rank('class');
}

sub get_class {
    return $_[0]->_get_taxon_name_for_rank('class');
}

sub get_class_confidence {
    my $self = shift;
    my $class_taxon = $self->get_class_taxon;
    unless ($class_taxon) {
        return;
    }
    my ($confidence) = $class_taxon->get_tag_values('confidence');
    return $confidence;
}

sub get_order_taxon {
    return $_[0]->_get_taxon_for_rank('order');
}

sub get_order {
    return $_[0]->_get_taxon_name_for_rank('order');
}

sub get_order_confidence {
    my $self = shift;
    my $order_taxon = $self->get_order_taxon;
    unless ($order_taxon) {
        return;
    }
    my ($confidence) = $order_taxon->get_tag_values('confidence');
    return $confidence;
}

sub get_family_taxon {
    return $_[0]->_get_taxon_for_rank('family');
}

sub get_family {
    return $_[0]->_get_taxon_name_for_rank('family');
}

sub get_family_confidence {
    my $self = shift;
    my $family_taxon = $self->get_family_taxon;
    unless ($family_taxon) {
        return;
    }
    my ($confidence) = $family_taxon->get_tag_values('confidence');
    return $confidence;
}

sub get_genus_taxon {
    return $_[0]->_get_taxon_for_rank('genus');
}

sub get_genus {
    return $_[0]->_get_taxon_name_for_rank('genus');
}

sub get_genus_confidence {
    my $self = shift;
    my $genus_taxon = $self->get_genus_taxon;
    unless ($genus_taxon) {
        return;
    }
    my ($confidence) = $genus_taxon->get_tag_values('confidence');
    return $confidence;
}

sub get_species_taxon {
    return $_[0]->_get_taxon_for_rank('species');
}

sub get_species {
    return $_[0]->_get_taxon_name_for_rank('species');
}

sub get_species_confidence {
    my $self = shift;
    my $species_taxon = $self->get_species_taxon;
    unless ($species_taxon) {
        return;
    }
    my ($confidence) = $species_taxon->get_tag_values('confidence');
    return $confidence;
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

#$HeadURL: svn+ssh://svn/srv/svn/gscpan/perl_modules/trunk/Genome/Utility/MetagenomicClassifier/SequenceClassification.pm $
#$Id: SequenceClassification.pm 43284 2009-02-04 22:15:30Z ebelter $

