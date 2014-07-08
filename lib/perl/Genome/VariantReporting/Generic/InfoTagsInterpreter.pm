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

sub requires_experts {
    return ();
}

sub available_fields {
    return qw/info_tags/;
}

sub interpret_entry {
    my $self = shift;
    my $entry = shift;

    my $info_types = $entry->{header}->info_types;

    my @matching_info_types;
    for my $info_type (keys %$info_types) {
        if (defined($entry->info($info_type))) {
            push @matching_info_types, $info_type;
        }
    }

    my %return_values;
    for my $alt_allele (@{$entry->{alternate_alleles}}) {
        $return_values{$alt_allele} = {
            info_tags => join(",", @matching_info_types),
        }
    }

    return %return_values;
}

1;
