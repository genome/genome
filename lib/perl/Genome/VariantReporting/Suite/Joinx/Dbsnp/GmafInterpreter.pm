package Genome::VariantReporting::Suite::Joinx::Dbsnp::GmafInterpreter;

use strict;
use warnings;
use Genome;

class Genome::VariantReporting::Suite::Joinx::Dbsnp::GmafInterpreter {
    is => 'Genome::VariantReporting::Framework::Component::Interpreter',
    has => [
        dummy => {
            is_optional => 1,
        },
    ],
};

sub name {
    return 'gmaf'
}

sub requires_annotations {
    return qw/
        dbsnp
    /;
}

sub field_descriptions {
    return (
        gmaf => 'Global minor allele frequency from dbsnp',
    );
}

sub _interpret_entry {
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

