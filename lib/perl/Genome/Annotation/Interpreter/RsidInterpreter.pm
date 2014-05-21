package Genome::Annotation::Interpreter::RsidInterpreter;

use strict;
use warnings;
use Genome;

class Genome::Annotation::Interpreter::RsidInterpreter {
    is => 'Genome::Annotation::Interpreter::Base',
};

sub name {
    return 'rsid';
}

sub requires_experts {
    ();
}

sub available_fields {
    return qw/
        rsid
    /;
}

sub interpret_entry {
    my $self = shift;
    my $entry = shift;

    my $passed_alt_alleles = shift;

    my %return_values;
    for my $variant_allele (@$passed_alt_alleles) {
        $return_values{$variant_allele} = {
            rsid => join(",", @{$entry->{identifiers}}),
        };
    }
    return %return_values;
}

1;

