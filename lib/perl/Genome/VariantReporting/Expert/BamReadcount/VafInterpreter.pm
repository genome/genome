package Genome::VariantReporting::Expert::BamReadcount::VafInterpreter;

use strict;
use warnings;
use Genome;
use Genome::VariantReporting::Expert::BamReadcount::VafCalculator;

class Genome::VariantReporting::Expert::BamReadcount::VafInterpreter {
    is => ['Genome::VariantReporting::Component::Interpreter', 'Genome::VariantReporting::Expert::BamReadcount::ComponentBase'],
};

sub name {
    return 'vaf';
}

sub requires_experts {
    return ('bam-readcount');
}

sub available_fields {
    return qw/
        vaf
        ref_count
        var_count
    /;
}

sub interpret_entry {
    my $self = shift;
    my $entry = shift;
    my $passed_alt_alleles = shift;

    my %return_values;
    my @sample_alt_alleles = sort $entry->alt_bases_for_sample($self->sample_index($entry->{header}));
    my %vafs = Genome::VariantReporting::Expert::BamReadcount::VafCalculator::calculate_vaf_for_all_alts(
        $entry, $self->get_readcount_entry($entry));


    for my $allele (@$passed_alt_alleles) {
        my $ref_allele;
        if (Genome::VariantReporting::Expert::BamReadcount::VafCalculator::is_deletion($entry->{reference_allele}, $allele)) {
            $ref_allele = substr($entry->{reference_allele}, 1, 1);
        }
        else {
            $ref_allele = substr($entry->{reference_allele}, 0, 1);
        }
        my $ref_count = Genome::VariantReporting::Expert::BamReadcount::VafCalculator::calculate_coverage_for_allele($self->get_readcount_entry($entry), $ref_allele, 'A');
        my $vaf;
        if (defined $vafs{$allele}) {
            $vaf = $vafs{$allele};
        }
        else {
            $vaf = undef;
        }
        $return_values{$allele} = {
            vaf => $vaf,
            var_count => Genome::VariantReporting::Expert::BamReadcount::VafCalculator::calculate_coverage_for_allele($self->get_readcount_entry($entry), $allele, $entry->{reference_allele}),
            ref_count => $ref_count,
        }
    }

    return %return_values;
}

1;

