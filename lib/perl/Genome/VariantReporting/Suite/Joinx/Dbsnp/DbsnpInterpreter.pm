package Genome::VariantReporting::Suite::Joinx::Dbsnp::DbsnpInterpreter;

use strict;
use warnings;
use Genome;
use Genome::File::Vcf::DbsnpAFParser;
use Memoize qw();

class Genome::VariantReporting::Suite::Joinx::Dbsnp::DbsnpInterpreter {
    is => 'Genome::VariantReporting::Framework::Component::Interpreter',
    doc => 'Output the allele frequency from dbsnp CAF field',
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

Memoize::memoize('_caf_parser', LIST_CACHE => 'MERGE');

sub _interpret_entry {
    my $self = shift;
    my $entry = shift;
    my $passed_alt_alleles = shift;

    my $parser = _caf_parser($entry->{header});
    my $caf = $parser->process_entry($entry);
    my %return_values;

    for my $variant_allele (@$passed_alt_alleles) {
        if (defined($entry->info('CAF')) && !defined($caf->{$variant_allele})) {
            # This may be fine, but I want to hear about it.  I think it can
            # happen in two ways:
            #
            # 1) the dbsnp vcf has a '.' as the vaf for this allele
            # 2) the dbsnp vcf doesn't have this allele at all.  This entry
            #    can still have the CAF field, but only if there is another
            #    allele that IS in the dbsnp vcf.
            $self->debug_message("No allele frequency for allele (%s) in entry: %s",
                $variant_allele, $entry->to_string);
        }
        $return_values{$variant_allele}->{allele_frequency} = $caf->{$variant_allele};
    }

    return %return_values;
}

1;

