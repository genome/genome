package Genome::VariantReporting::Generic::GenotypeInterpreter;

use strict;
use warnings;

use Genome;

class Genome::VariantReporting::Generic::GenotypeInterpreter {
      is => ['Genome::VariantReporting::Framework::Component::Interpreter', 'Genome::VariantReporting::Framework::Component::WithSampleName'],
};

sub name {
    return 'genotype';
}

sub requires_experts {
    return qw/ /;
}

sub available_fields {
    return qw/
        genotype
    /;
}

sub interpret_entry {
    my $self = shift;
    my $entry = shift;

    my @sample_alt_alleles = sort $entry->alt_bases_for_sample($self->sample_index($entry->{header}));
    my $gt = $entry->genotype_for_sample($self->sample_index($entry->{header}));
    my %return_values = map { $_ => { genotype => 'not called' } } @{$entry->{alternate_alleles}};

    for my $alt_allele ( @sample_alt_alleles ) {
        if ( $gt->is_missing ) {
            # fail
            $return_values{$alt_allele} =  { genotype => 'missing' };
        }
        elsif ( $gt->is_homozygous ) {
            $return_values{$alt_allele} = { genotype => 'homozygous' };
        }
        else { # $gt->is_heterozygous
            $return_values{$alt_allele} = { genotype => 'heterozygous' };
        }
    }

    return %return_values;
}

1;

