package Genome::VariantReporting::Suite::BamReadcount::VafInterpreter;

use strict;
use warnings;
use Genome;
use Genome::VariantReporting::Suite::BamReadcount::VafCalculator;
use Genome::VariantReporting::Suite::BamReadcount::VafInterpreterHelpers qw(
    many_libraries_field_descriptions
    many_samples_field_descriptions
);

class Genome::VariantReporting::Suite::BamReadcount::VafInterpreter {
    is => ['Genome::VariantReporting::Framework::Component::Interpreter',
        'Genome::VariantReporting::Suite::BamReadcount::ComponentBase',
        'Genome::VariantReporting::Framework::Component::WithManyLibraryNames',],
    has_optional => [
        library_names => {
            is => 'Text',
            is_many => 1,
            is_translated => 1,
            doc => 'List of library names to be used in the report',
        },
    ],
    doc => 'Calculate the variant allele frequency, number of reads supporting the reference, and number of reads supporting variant for a sample and its libraries',
};

sub name {
    return 'vaf';
}

sub requires_annotations {
    return ('bam-readcount');
}

sub field_descriptions {
    my $self = shift;
    return (many_samples_field_descriptions($self),
    many_libraries_field_descriptions($self));
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
            my $translated_reference_allele = $self->translate_ref_allele($entry->{reference_allele}, $allele);
            my $vaf = $vafs{$allele};

            $return_values{$allele} = {
                $self->create_sample_specific_field_name("vaf") => $vaf,
                $self->create_sample_specific_field_name("var_count") =>
                Genome::VariantReporting::Suite::BamReadcount::VafCalculator::calculate_coverage_for_allele(
                    $readcount_entry, $allele, $entry->{reference_allele}),
                $self->create_sample_specific_field_name("ref_count") =>
                Genome::VariantReporting::Suite::BamReadcount::VafCalculator::calculate_coverage_for_allele(
                    $readcount_entry, $translated_reference_allele, 'A'),
                $self->flatten_hash($self->per_library_vaf($entry, $readcount_entry, $allele), "vaf"),
                $self->flatten_hash($self->per_library_coverage($readcount_entry, $allele, $entry->{reference_allele}), "var_count"),
                $self->flatten_hash($self->per_library_coverage($readcount_entry, $translated_reference_allele, 'A'), "ref_count"),
            }
        }
    }

    return %return_values;
}

sub per_library_vaf {
    my ($self, $entry, $readcount_entry, $allele) = @_;

    return Genome::VariantReporting::Suite::BamReadcount::VafCalculator::calculate_per_library_vaf_for_all_alts($entry, $readcount_entry)->{$allele};
}


# When checking for variant coverage: The $reference_allele must be untranslated
# When checking for reference coverage: The $reference_allele and $allele must both be the TRANSLATED reference
## This is because otherwise we will misinterpret the query as asking for insertion or deletion support inside the VafCalculator
sub per_library_coverage {
    my ($self, $readcount_entry, $allele, $reference_allele) = @_;
    return Genome::VariantReporting::Suite::BamReadcount::VafCalculator::calculate_per_library_coverage_for_allele($readcount_entry, $allele, $reference_allele);
}

sub flatten_hash {
    my ($self, $per_library_hash, $field_name) = @_;
    my %flattened_hash;
    for my $library_name (keys %$per_library_hash) {
        $flattened_hash{$self->create_library_specific_field_name($field_name, $library_name)} = $per_library_hash->{$library_name};
    }
    return %flattened_hash;
}

sub translate_ref_allele {
    my ($self, $ref, $alt) = @_;
    if (Genome::VariantReporting::Suite::BamReadcount::VafCalculator::is_deletion($ref, $alt)) {
        return substr($ref, 1, 1);
    }
    else {
        return substr($ref, 0, 1);
    }
}

1;

