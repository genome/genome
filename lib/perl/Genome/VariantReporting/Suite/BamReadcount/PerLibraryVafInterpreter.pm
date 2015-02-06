package Genome::VariantReporting::Suite::BamReadcount::PerLibraryVafInterpreter;

use strict;
use warnings;
use Genome;
use Genome::VariantReporting::Suite::BamReadcount::VafInterpreterHelpers qw(
    many_libraries_field_descriptions
    per_sample_field_descriptions
    translate_ref_allele
);
use Set::Scalar;
use List::MoreUtils qw(uniq);

class Genome::VariantReporting::Suite::BamReadcount::PerLibraryVafInterpreter {
    is => [
        'Genome::VariantReporting::Framework::Component::Interpreter',
        'Genome::VariantReporting::Framework::Component::WithManySampleNames',
        'Genome::VariantReporting::Framework::Component::WithManyLibraryNames',
    ],
    has => [],
    doc => 'Calculate the variant allele frequency, number of reads supporting the reference, and number of reads supporting variant for the libraries of multiple samples',
};

sub name {
    return 'per-library-vaf';
}

sub requires_annotations {
    return ('bam-readcount');
}

sub field_descriptions {
    my $self = shift;
    return many_libraries_field_descriptions($self);
}

sub _interpret_entry {
    my $self = shift;
    my $entry = shift;
    my $passed_alt_alleles = shift;

    if ($self->library_names) {
        my $available_libraries = Set::Scalar->new($self->available_libraries($entry));
        my $expected_libraries  = Set::Scalar->new($self->library_names);
        unless ($available_libraries->is_equal($expected_libraries)) {
            die $self->error_message(
                "Available libraries (%s) are not identical to planned libraries (%s)",
                join(', ', $available_libraries->members), join(', ', $expected_libraries->members)
            );
        }
    }
    my %return_values;

    for my $sample_name ($self->sample_names) {
        my $readcount_entries = $self->get_readcount_entries($entry, $sample_name);
        unless (defined($readcount_entries)) {
            for my $alt_allele (@$passed_alt_alleles) {
                $return_values{$alt_allele} = {map {$_ => $self->interpretation_null_character} per_sample_field_descriptions($sample_name)};
            }
        }

        for my $alt_allele (@$passed_alt_alleles) {
            my $readcount_entry = $readcount_entries->{$alt_allele};
            if (!defined $readcount_entry) {
                $return_values{$alt_allele} = {map {$_ => $self->interpretation_null_character} per_sample_field_descriptions($sample_name)};
            }
            else {
                my $translated_reference_allele = translate_ref_allele($entry->{reference_allele}, $alt_allele);
                my %results = (
                    $self->flatten_hash($self->per_library_vaf($entry, $readcount_entry, $alt_allele), "vaf"),
                    $self->flatten_hash($self->per_library_coverage($readcount_entry, $alt_allele, $entry->{reference_allele}), "var_count"),
                    $self->flatten_hash($self->per_library_coverage($readcount_entry, $translated_reference_allele, 'A'), "ref_count"),
                );
                for my $field_name (keys %results) {
                    $return_values{$alt_allele}->{$field_name} = $results{$field_name};
                }
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

sub available_libraries {
    my ($self, $entry) = @_;

    my %readcount_entries;
    for my $sample_name ($self->sample_names) {
        %readcount_entries = (%readcount_entries, %{$self->get_readcount_entries($entry, $sample_name)});
    }
    return uniq map {$_->name} map {$_->libraries} grep {defined($_)} values %readcount_entries;
}

sub get_readcount_entries {
    my ($self, $entry, $sample_name) = @_;

    return Genome::File::Vcf::BamReadcountParser::get_bam_readcount_entries(
        $entry,
        $entry->{header}->index_for_sample_name($sample_name),
    );
}

1;
