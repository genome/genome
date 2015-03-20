package Genome::VariantReporting::Suite::BamReadcount::VafInterpreter;

use strict;
use warnings;
use Genome;
use Genome::VariantReporting::Suite::BamReadcount::VafInterpreterHelpers qw(
    many_samples_field_descriptions
    per_sample_field_descriptions
    translate_ref_allele
);

class Genome::VariantReporting::Suite::BamReadcount::VafInterpreter {
    is => [
        'Genome::VariantReporting::Framework::Component::Interpreter',
        'Genome::VariantReporting::Framework::Component::WithManySampleNames',
    ],
    has => [],
    doc => 'Calculate the variant allele frequency, number of reads supporting the reference, and number of reads supporting variant for multiple samples',
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
    for my $sample_name ($self->sample_names) {
        my $readcount_entries = $self->get_readcount_entries($entry, $sample_name);
        unless (defined($readcount_entries)) {
            for my $alt_allele (@$passed_alt_alleles) {
                $return_values{$alt_allele} = {map {$_ => $self->interpretation_null_character} per_sample_field_descriptions($sample_name)};
            }
        }

        my %vafs = Genome::VariantReporting::Suite::BamReadcount::VafCalculator::calculate_vaf_for_all_alts(
            $entry, $readcount_entries
        );

        for my $alt_allele (@$passed_alt_alleles) {
            my $readcount_entry = $readcount_entries->{$alt_allele};
            if (!defined $readcount_entry) {
                $return_values{$alt_allele} = {map {$_ => $self->interpretation_null_character} per_sample_field_descriptions($sample_name)};
            }
            else {
                my $translated_reference_allele = translate_ref_allele($entry->{reference_allele}, $alt_allele);
                my %results = (
                    $self->create_sample_specific_field_name('vaf', $sample_name) => $vafs{$alt_allele},
                    $self->create_sample_specific_field_name('var_count', $sample_name) =>
                        Genome::VariantReporting::Suite::BamReadcount::VafCalculator::calculate_coverage_for_allele(
                            $readcount_entry, $alt_allele, $entry->{reference_allele}
                        ),
                    $self->create_sample_specific_field_name('ref_count', $sample_name) =>
                        Genome::VariantReporting::Suite::BamReadcount::VafCalculator::calculate_coverage_for_allele(
                            $readcount_entry, $translated_reference_allele, 'A'
                        ),
                );
                for my $field_name (keys %results) {
                    $return_values{$alt_allele}->{$field_name} = $results{$field_name};
                }
            }
        }
    }
    return %return_values;
}

sub get_readcount_entries {
    my ($self, $entry, $sample_name) = @_;

    return Genome::File::Vcf::BamReadcountParser::get_bam_readcount_entries(
        $entry,
        $entry->{header}->index_for_sample_name($sample_name),
    );
}

1;
