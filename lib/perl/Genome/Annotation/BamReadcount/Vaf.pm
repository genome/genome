package Genome::Annotation::BamReadcount::Vaf;

use strict;
use warnings;
use Genome;

class Genome::Annotation::BamReadcount::Vaf {
    is => 'Genome::Annotation::BamReadcount::FilterBase',
    has => [
        min_vaf => {
            is => 'Number',
        },
    ],
};

sub process_entry {
    my $self = shift;
    my $entry = shift;

    my $allele = $self->get_allele_from_entry($entry);
    my $vaf = $self->calculate_vaf($entry, $allele);
    return $vaf >= $self->min_vaf;
}

sub get_allele_from_entry {
    my $self = shift;
    my $entry = shift;
    my $genotype = $entry->genotype_for_sample($self->sample_index);
    my @entry_allele_amino_acids = $entry->alleles;
    my @genotype_allele_pointers = $genotype->get_alleles;
    my @genotype_allele_amino_acids = map {$entry_allele_amino_acids[$_]} @genotype_allele_pointers;
    return $genotype_allele_amino_acids[1];
}

sub calculate_vaf {
    my $self = shift;
    my ($entry, $allele) = @_;
    my $bam_readcount_entry = $self->get_readcount_entry($entry);

    my $var_count = 0;
    for my $lib ($bam_readcount_entry->libraries) {
        $var_count += $lib->metrics_for($allele)->count;
    }

    return $var_count / $bam_readcount_entry->depth * 100;
}

1;

