package Genome::VariantReporting::Suite::Vep::DamagingFilter;

use strict;
use warnings;
use Genome;
use Genome::File::Vcf::VepConsequenceParser;
use Memoize qw();

class Genome::VariantReporting::Suite::Vep::DamagingFilter {
    is => 'Genome::VariantReporting::Framework::Component::Filter',
    doc => q{Filter out variants that aren't considered damaging},
};

sub name {
    return 'damaging';
}

sub requires_annotations {
    return ('vep');
}

sub filter_entry {
    my $self = shift;
    my $entry = shift;

    my %return_values;

    my $vep_parser = vep_parser($entry->{header});

    for my $alt_allele (@{$entry->{alternate_alleles}}) {
        my ($transcript) = $vep_parser->transcripts($entry, $alt_allele);
        $return_values{$alt_allele} = is_damaging($transcript);
    }

    return %return_values;
}

sub vep_parser {
    my $header = shift;
    return new Genome::File::Vcf::VepConsequenceParser($header);
}
Memoize::memoize('vep_parser', LIST_CACHE => 'MERGE');

sub is_damaging {
    my $transcript = shift;

    if (defined($transcript->{'consequence'})) {
        my $expression = always_damaging_expression();
        if ($transcript->{consequence} =~ qr($expression)) {
            return 1;
        }
        $expression = non_synonymous_expression();
        if ($transcript->{consequence} =~ qr($expression)) {
            $expression = non_synonymous_damaging_expression();
            if ($transcript->{polyphen} and ($transcript->{polyphen} =~ qr($expression))) {
                return 1;
            }
            if ($transcript->{sift} and ($transcript->{sift} =~ qr($expression))) {
                return 1;
            }
            if (!$transcript->{sift} and !$transcript->{polyphen}) {
                return 1;
            }
        }
    }

    return 0;
}

sub always_damaging_expression {
    my @always_damaging = qw(
        transcript_ablation
        transcript_amplification
        splice_donor_variant
        splice_acceptor_variant
        stop_gained
        frameshift_variant
        stop_lost
        mature_miRNA_variant);
    return join("|", @always_damaging);
}

sub non_synonymous_expression {
    my @non_synonymous = qw(
        initiator_codon_variant
        inframe_insertion
        inframe_deletion
        missense_variant
    );
    return join("|", @non_synonymous);
}

sub non_synonymous_damaging_expression {
    my @damaging = qw(
        damaging
        deleterious
    );
    return join("|", @damaging);
}

sub vcf_id {
    return 'DAMAGING';
}

sub vcf_description {
    return 'Filter out variants that aren\'t considered to be damaging';
}

1;

