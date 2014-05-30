package Genome::VariantReporting::Expert::BamReadcount::VafInterpreter;

use strict;
use warnings;
use Genome;
use Genome::VariantReporting::Expert::BamReadcount::VafCalculator;

class Genome::VariantReporting::Expert::BamReadcount::VafInterpreter {
    is => ['Genome::VariantReporting::Interpreter::Base', 'Genome::VariantReporting::Expert::BamReadcount::ComponentBase'],
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
    my %vafs = Genome::VariantReporting::Expert::BamReadcount::VafCalculator::calculate_vaf_for_multiple_alleles(
        $self->get_readcount_entry($entry), \@sample_alt_alleles);

    my $ref_count = Genome::Annotation::Expert::BamReadcount::VafCalculator::calculate_coverage_for_allele($self->get_readcount_entry($entry), $entry->{reference_allele});

    for my $allele (@$passed_alt_alleles) {
        my $vaf;
        if (defined $vafs{$allele}) {
            $vaf = $vafs{$allele};
        }
        else {
            $vaf = undef;
        }
        $return_values{$allele} = {
            vaf => $vaf,
            var_count => Genome::Annotation::Expert::BamReadcount::VafCalculator::calculate_coverage_for_allele($self->get_readcount_entry($entry), $allele),
            ref_count => $ref_count,
        }
    }

    return %return_values;
}

1;

