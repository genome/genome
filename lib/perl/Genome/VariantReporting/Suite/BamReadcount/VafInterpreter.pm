package Genome::VariantReporting::Suite::BamReadcount::VafInterpreter;

use strict;
use warnings;
use Genome;
use Genome::VariantReporting::Suite::BamReadcount::VafCalculator;
use Genome::VariantReporting::Suite::BamReadcount::VafInterpreterHelpers qw(basic_field_descriptions);

class Genome::VariantReporting::Suite::BamReadcount::VafInterpreter {
    is => ['Genome::VariantReporting::Framework::Component::Interpreter', 'Genome::VariantReporting::Suite::BamReadcount::ComponentBase'],
};

sub name {
    return 'vaf';
}

sub requires_annotations {
    return ('bam-readcount');
}

sub field_descriptions {
    my $self = shift;
    return basic_field_descriptions($self->sample_name);
}

sub _interpret_entry {
    my $self = shift;
    my $entry = shift;
    my $passed_alt_alleles = shift;

    my %return_values;

    my $readcount_entry = $self->get_readcount_entry($entry);
    return $self->null_interpretation($passed_alt_alleles) unless defined($readcount_entry);

    my %vafs = Genome::VariantReporting::Suite::BamReadcount::VafCalculator::calculate_vaf_for_all_alts(
        $entry, $readcount_entry);

    for my $allele (@$passed_alt_alleles) {
        my $translated_reference_allele = $self->translate_ref_allele($entry->{reference_allele}, $allele);
        my $vaf;
        if (defined $vafs{$allele}) {
            $vaf = $vafs{$allele};
        }
        else {
            $vaf = undef;
        }

        $return_values{$allele} = {
            vaf => $vaf,
            var_count => Genome::VariantReporting::Suite::BamReadcount::VafCalculator::calculate_coverage_for_allele($readcount_entry, $allele, $entry->{reference_allele}),
            per_library_var_count => $self->per_library_coverage($readcount_entry, $allele, $entry->{reference_allele}),
            ref_count => Genome::VariantReporting::Suite::BamReadcount::VafCalculator::calculate_coverage_for_allele($readcount_entry, $translated_reference_allele, 'A'),
            per_library_ref_count => $self->per_library_coverage($readcount_entry, $translated_reference_allele, 'A'),
            per_library_vaf => $self->per_library_vaf($entry, $readcount_entry, $allele),
        }
    }

    return %return_values;
}

sub per_library_vaf {
    my ($self, $entry, $readcount_entry, $allele) = @_;

    my $per_library_vafs = Genome::VariantReporting::Suite::BamReadcount::VafCalculator::calculate_per_library_vaf_for_all_alts($entry, $readcount_entry);

    return join(',', map {"$_:" . $per_library_vafs->{$allele}{$_} } keys %{$per_library_vafs->{$allele}});
}


# When checking for variant coverage: The $reference_allele must be untranslated
# When checking for reference coverage: The $reference_allele and $allele must both be the TRANSLATED reference
## This is because otherwise we will misinterpret the query as asking for insertion or deletion support inside the VafCalculator
sub per_library_coverage {
    my ($self, $readcount_entry, $allele, $reference_allele) = @_;
    my $counts = Genome::VariantReporting::Suite::BamReadcount::VafCalculator::calculate_per_library_coverage_for_allele($readcount_entry, $allele, $reference_allele);
    return join(',', map {"$_:" . $counts->{$_} } keys %$counts);
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

