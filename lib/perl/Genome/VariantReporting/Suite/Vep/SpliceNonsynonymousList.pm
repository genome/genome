package Genome::VariantReporting::Suite::Vep::SpliceNonsynonymousList;

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


sub splice_sites {
    return @splice_sites;
}


sub nonsynonymous_list {
    return @non_synonymous;
}


sub is_splice_site {
    my $splice_sites = Set::Scalar->new(@splice_sites);
    return !$splice_sites->intersection(Set::Scalar->new(@_))->is_null;
}


sub is_non_synonymous {
    my $non_synonymous = Set::Scalar->new(@non_synonymous);
    return !$non_synonymous->intersection(Set::Scalar->new(@_))->is_null;
}

1;
