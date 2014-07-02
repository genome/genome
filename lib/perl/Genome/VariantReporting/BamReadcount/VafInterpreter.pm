package Genome::VariantReporting::BamReadcount::VafInterpreter;

use strict;
use warnings;
use Genome;
use Genome::VariantReporting::BamReadcount::VafCalculator;

class Genome::VariantReporting::BamReadcount::VafInterpreter {
    is => ['Genome::VariantReporting::Framework::Component::Interpreter', 'Genome::VariantReporting::BamReadcount::ComponentBase'],
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
        per_library_var_count
    /;
}

sub interpret_entry {
    my $self = shift;
    my $entry = shift;
    my $passed_alt_alleles = shift;

    my %return_values;
    my @sample_alt_alleles = sort $entry->alt_bases_for_sample($self->sample_index($entry->{header}));
    my %vafs = Genome::VariantReporting::BamReadcount::VafCalculator::calculate_vaf_for_all_alts(
        $entry, $self->get_readcount_entry($entry));

    for my $allele (@$passed_alt_alleles) {
        my $translated_reference_allele = $self->translate_ref_allele($entry->{reference_allele}, $allele);
        my $ref_count = Genome::VariantReporting::BamReadcount::VafCalculator::calculate_coverage_for_allele($self->get_readcount_entry($entry), $translated_reference_allele, 'A');
        my $vaf;
        if (defined $vafs{$allele}) {
            $vaf = $vafs{$allele};
        }
        else {
            $vaf = undef;
        }
        my $readcount_entry = $self->get_readcount_entry($entry);

        $return_values{$allele} = {
            vaf => $vaf,
            var_count => Genome::VariantReporting::BamReadcount::VafCalculator::calculate_coverage_for_allele($readcount_entry, $allele, $entry->{reference_allele}),
            ref_count => $ref_count,
            per_library_var_count => $self->per_library_coverage($readcount_entry, $allele, $entry->{reference_allele}),
        }
    }

    return %return_values;
}

# The $reference_allele must be untranslated. When we are interested in counts for reference, $allele should be the translated ref allele.
sub per_library_coverage {
    my ($self, $readcount_entry, $allele, $reference_allele) = @_;
    my $counts = Genome::VariantReporting::BamReadcount::VafCalculator::calculate_per_library_coverage_for_allele($readcount_entry, $allele, $reference_allele);
    return join(',', map {"$_:" . $counts->{$_} } keys %$counts);
}

sub translate_ref_allele {
    my ($self, $ref, $alt) = @_;
    if (Genome::VariantReporting::BamReadcount::VafCalculator::is_deletion($ref, $alt)) {
        return substr($ref, 1, 1);
    }
    else {
        return substr($ref, 0, 1);
    }
}

1;

