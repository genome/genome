package Genome::VariantReporting::Suite::Joinx::Dbsnp::DbsnpInterpreter;

use strict;
use warnings;
use Genome;
use Genome::File::Vcf::DbsnpAFParser;

class Genome::VariantReporting::Suite::Joinx::Dbsnp::DbsnpInterpreter {
    is => 'Genome::VariantReporting::Framework::Component::Interpreter',
};

sub name {
    return 'dbsnp'
}

sub requires_annotations {
    return qw/
        dbsnp
    /;
}

sub field_descriptions {
    return (
        allele_frequency => 'Allele frequency from dbsnp CAF field',
    );
}

sub _caf_parser {
    my $header = shift;
    return Genome::File::Vcf::DbsnpAFParser->new($header);
}

Memoize::memoize('_caf_parser');

sub _interpret_entry {
    my $self = shift;
    my $entry = shift;
    my $passed_alt_alleles = shift;

    my $parser = _caf_parser($entry->{header});
    my $caf = $parser->process_entry($entry);
    my %return_values;

    for my $variant_allele (@$passed_alt_alleles) {
        if (!defined $caf) {
            $return_values{$variant_allele}->{allele_frequency} = undef;
        }
        elsif (defined $caf->{$variant_allele}) {
            $return_values{$variant_allele}->{allele_frequency} = $caf->{$variant_allele};
        }
        else {
            die sprintf("No allele frequency for allele ($variant_allele) given in caf: %s",
                Data::Dumper::Dumper($caf));
        }
    }

    return %return_values;
}

1;

