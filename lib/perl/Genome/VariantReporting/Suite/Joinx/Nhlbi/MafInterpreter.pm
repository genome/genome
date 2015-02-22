package Genome::VariantReporting::Suite::Joinx::Nhlbi::MafInterpreter;

use strict;
use warnings;
use Genome;

class Genome::VariantReporting::Suite::Joinx::Nhlbi::MafInterpreter {
    is => ['Genome::VariantReporting::Suite::Joinx::Nhlbi::ComponentBase', 'Genome::VariantReporting::Framework::Component::Interpreter'],
    doc => 'Output the European, African American, and overall minor allele frequencies from NHLBI',
};

sub name {
    return 'nhlbi'
}

sub requires_annotations {
    return qw/nhlbi/;
}

sub field_descriptions {
    return (
        NHLBI_EU_MAF => 'European minor allele frequency from NHLBI',
        NHLBI_AA_MAF => 'African American minor allele frequency from NHLBI',
        NHLBI_All_MAF => 'Overall minor allele frequency from NHLBI',
    );
}

sub _interpret_entry {
    my $self = shift;
    my $entry = shift;
    my $passed_alt_alleles = shift;

    my %return_values;
    my $eu_maf = $self->get_maf_for_entry($entry, "EU");
    my $aa_maf = $self->get_maf_for_entry($entry, "AA");
    my $all_maf = $self->get_maf_for_entry($entry, "All");
    for my $variant_allele (@$passed_alt_alleles) {
        $return_values{$variant_allele}->{NHLBI_EU_MAF} = $eu_maf;
        $return_values{$variant_allele}->{NHLBI_AA_MAF} = $aa_maf;
        $return_values{$variant_allele}->{NHLBI_All_MAF} = $all_maf;
    }
    return %return_values;
}
1;

