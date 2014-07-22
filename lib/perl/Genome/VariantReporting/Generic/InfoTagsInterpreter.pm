package Genome::VariantReporting::Generic::InfoTagsInterpreter;

use strict;
use warnings;
use Genome;

class Genome::VariantReporting::Generic::InfoTagsInterpreter {
    is => 'Genome::VariantReporting::Framework::Component::Interpreter',
};

sub name {
    return 'info-tags';
}

sub requires_annotations {
    return ();
}

sub available_fields {
    return qw/info_tags/;
}

sub process_interpret_entry {
    my $self = shift;
    my $entry = shift;
    my $passed_alt_alleles = shift;

    my %return_values;
    for my $alt_allele (@$passed_alt_alleles) {
        $return_values{$alt_allele} = {
            info_tags => join(",", keys %{$entry->info_for_allele($alt_allele)}),
        }
    }

    return %return_values;
}

1;
