package Genome::VariantReporting::BamReadcount::ManySamplesVafInterpreter;

use strict;
use warnings;
use Genome;

class Genome::VariantReporting::BamReadcount::ManySamplesVafInterpreter {
    is => ['Genome::VariantReporting::Framework::Component::Interpreter', 'Genome::VariantReporting::Framework::Component::WithManySampleNames'],
    has => [],
};

sub name {
    return 'many-samples-vaf';
}

sub requires_annotations {
    return ('bam-readcount');
}

sub available_fields {
    return qw/
        vaf
        ref_count
        var_count
        per_library_var_count
        per_library_ref_count
        per_library_vaf
    /;
}

sub interpret_entry {
    my $self = shift;
    my $entry = shift;
    my $passed_alt_alleles = shift;

    my %return_values;
    for my $sample_name ($self->sample_names) {
        my $interpreter = Genome::VariantReporting::BamReadcount::VafInterpreter->create(sample_name => $sample_name);
        my %result = $interpreter->interpret_entry($entry, $passed_alt_alleles);
        for my $alt_allele (@$passed_alt_alleles) {
            $return_values{$alt_allele}->{$sample_name} = $result{$alt_allele};
        }
    }
    return %return_values;
}

1;
