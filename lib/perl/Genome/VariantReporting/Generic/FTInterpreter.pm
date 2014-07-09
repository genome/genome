package Genome::VariantReporting::Generic::FTInterpreter;

use strict;
use warnings FATAL => 'all';
use Genome;
use List::AllUtils qw/any/;

class Genome::VariantReporting::Generic::FTInterpreter {
    is => ['Genome::VariantReporting::Framework::Component::Interpreter', 'Genome::VariantReporting::Framework::Component::WithSampleName'],
};

sub name {
    return 'ft';
}

sub requires_experts {
    ();
}

sub available_fields {
    return qw/
        ft_string
    /;
}

sub interpret_entry {
    my ($self, $entry, $passed_alt_alleles) = @_;

    my $ft_string = $entry->sample_field($self->sample_index($entry->{header}), 'FT');

    my %return_values;
    for my $alt_allele (@$passed_alt_alleles) {
        $return_values{$alt_allele} = { ft_string => "" };
    }

    my @sample_alt_alleles = $entry->alt_bases_for_sample($self->sample_index($entry->{header}));
    for my $alt_allele (keys %return_values) {
        if (any { $_ eq $alt_allele } @sample_alt_alleles) {
            $return_values{$alt_allele} =  {
                ft_string => $ft_string
            };
        }
    }
    return %return_values;
}


1;
