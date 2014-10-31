package Genome::VariantReporting::Suite::Vep::AnnotationCategory;

use strict;
use warnings;
use Genome;


my @splice_sites = qw(
    splice_acceptor_variant
    splice_donor_variant
);


my @non_synonymous = qw(
    transcript_ablation
    stop_gained
    frameshift_variant
    stop_lost
    initiator_codon_variant
    transcript_amplification
    inframe_insertion
    inframe_deletion
    missense_variant
    incomplete_terminal_codon_variant
);


sub splice_site {
    return @splice_sites;
}


sub non_synonymous {
    return @non_synonymous;
}


sub is_category {
    my ($class, $category_name, @terms) = @_;

    unless (__PACKAGE__->can($category_name)) {
        Carp::confess("category $category_name is not implemented");
    }

    my $category = Set::Scalar->new(__PACKAGE__->$category_name);
    return !$category->intersection(Set::Scalar->new(@terms))->is_null;
}

1;
