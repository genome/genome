package Genome::VariantReporting::Interpreter::GmafInterpreter;

use strict;
use warnings;
use Genome;

class Genome::VariantReporting::Interpreter::GmafInterpreter {
    is => 'Genome::VariantReporting::Interpreter::Base',
    has => [
        dummy => {
            is_optional => 1,
        },
    ],
};

sub name {
    return 'gmaf'
}

sub requires_experts {
    return qw/
        dbsnp
    /;
}

sub available_fields {
    return qw/
        gmaf
    /;
}

sub interpret_entry {
    my $self = shift;
    my $entry = shift;
    my $passed_alt_alleles = shift;

    my %return_values;

    for my $variant_allele (@$passed_alt_alleles) {
        if (!defined $entry->info("GMAF")) {
            $return_values{$variant_allele}->{gmaf} = undef;
        }
        else {
            $return_values{$variant_allele}->{gmaf} = $entry->info("GMAF");
        }
    }

    return %return_values;
}

1;

