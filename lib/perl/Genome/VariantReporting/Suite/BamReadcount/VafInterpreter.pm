package Genome::VariantReporting::Suite::BamReadcount::VafInterpreter;

use strict;
use warnings;
use Genome;
use Genome::VariantReporting::Suite::BamReadcount::VafCalculator;
use Genome::VariantReporting::Suite::BamReadcount::VafInterpreterHelpers qw(
    many_samples_field_descriptions
    translate_ref_allele
);

class Genome::VariantReporting::Suite::BamReadcount::VafInterpreter {
    is => [
        'Genome::VariantReporting::Framework::Component::Interpreter',
        'Genome::VariantReporting::Suite::BamReadcount::ComponentBase',
    ],
    doc => 'Calculate the variant allele frequency, number of reads supporting the reference, and number of reads supporting variant for a sample',
};

sub name {
    return 'vaf';
}

sub requires_annotations {
    return ('bam-readcount');
}

sub field_descriptions {
    my $self = shift;
    return many_samples_field_descriptions($self);
}

sub _interpret_entry {
    my $self = shift;
    my $entry = shift;
    my $passed_alt_alleles = shift;

    my %return_values;

    my $readcount_entries = $self->get_readcount_entries($entry);
    return $self->null_interpretation($passed_alt_alleles) unless defined($readcount_entries);

    my %vafs = Genome::VariantReporting::Suite::BamReadcount::VafCalculator::calculate_vaf_for_all_alts(
        $entry, $readcount_entries);

    for my $allele (@$passed_alt_alleles) {
        my $readcount_entry = $readcount_entries->{$allele};
        if (!defined $readcount_entry) {
            $return_values{$allele} = {map {$_ => $self->interpretation_null_character} $self->available_fields};
        }
        else {
            my $translated_reference_allele = translate_ref_allele($entry->{reference_allele}, $allele);
            my $vaf = $vafs{$allele};

            $return_values{$allele} = {
                $self->create_sample_specific_field_name("vaf") => $vaf,
                $self->create_sample_specific_field_name("var_count") =>
                Genome::VariantReporting::Suite::BamReadcount::VafCalculator::calculate_coverage_for_allele(
                    $readcount_entry, $allele, $entry->{reference_allele}),
                $self->create_sample_specific_field_name("ref_count") =>
                Genome::VariantReporting::Suite::BamReadcount::VafCalculator::calculate_coverage_for_allele(
                    $readcount_entry, $translated_reference_allele, 'A'),
            }
        }
    }

    return %return_values;
}

1;

