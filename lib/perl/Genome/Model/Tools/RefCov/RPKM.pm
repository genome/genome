package Genome::Model::Tools::RefCov::RPKM;

use strict;
use warnings;

use Genome;

# RPKM (Reads Per Kilobase exon Model per million mapped reads)

class Genome::Model::Tools::RefCov::RPKM {
    has => [
        target_alignments_total => { },
        total_alignments => { },
        target_bp_total => { },
    ],
    has_optional => {
        rpkm => {},
    },
};

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);
    $self->_calculate_RPKM();
    return $self;
}

sub _calculate_RPKM {
    my $self = shift;

    #        10^9 x (the number of mappable reads aligned within a gene's exons)
    # RPKM = ---------------------------------------------------------------------------------------------------
    #        (the total number of mappable reads in the experiment) x (the sum of the exons in base pairs in Kb)

    my $rpkm = (10_000_000_000 * $self->target_alignments_total()) / ($self->total_alignments * $self->target_bp_total);
    $self->rpkm($rpkm);

    return $self;
}


1;  # End of package.
