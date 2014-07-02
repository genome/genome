package Genome::VariantReporting::Generic::FTInterpreter;

use strict;
use warnings FATAL => 'all';
use Genome;
use List::Util qw/first/;

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
    my ($self, $entry) = @_;

    my $ft_string = $entry->sample_field($self->sample_index($entry->{header}), 'FT');

    my %return_values;
    for my $alt_allele (@{$entry->{alternate_alleles}}) {
        $return_values{$alt_allele} =  {
            ft_string => $ft_string
        };
    }
    return %return_values;
}


1;
