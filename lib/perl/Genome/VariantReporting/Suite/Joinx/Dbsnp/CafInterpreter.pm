package Genome::VariantReporting::Suite::Joinx::Dbsnp::CafInterpreter;

use strict;
use warnings;
use Genome;

use List::Util qw(max);

class Genome::VariantReporting::Suite::Joinx::Dbsnp::CafInterpreter {
    is => ['Genome::VariantReporting::Suite::Joinx::Dbsnp::ComponentBase', 'Genome::VariantReporting::Framework::Component::Interpreter'],
};

sub name {
    return 'caf'
}

sub requires_annotations {
    return qw/
        dbsnp
    /;
}

sub field_descriptions {
    return (
        caf => 'Allele frequency at this position based on 1000Genomes',
        max_alt_af => 'The highest allele frequency at this position'
    );
}

sub _interpret_entry {
    my $self = shift;
    my $entry = shift;
    my $passed_alt_alleles = shift;

    my %return_values;
    my $parser = $self->_caf_parser($entry->{header});
    my $caf = $parser->process_entry($entry);
    delete $caf->{$entry->{reference_allele}};
    my $highest_alt_af = _highest_af($caf);
    for my $variant_allele (@$passed_alt_alleles) {
        if (!defined $entry->info("CAF")) {
            $return_values{$variant_allele}->{caf} = undef;
            $return_values{$variant_allele}->{max_alt_af} = undef;
        }
        else {
            $return_values{$variant_allele}->{max_alt_af} = $highest_alt_af;
            $return_values{$variant_allele}->{caf} = $caf->{$variant_allele};
        }
    }

    return %return_values;
}

sub _highest_af {
    my $caf = shift;

    return max(grep {defined $_ and $_ ne "."} values %$caf);
}
1;

