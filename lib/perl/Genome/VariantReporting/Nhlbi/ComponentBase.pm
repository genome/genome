package Genome::VariantReporting::Nhlbi::ComponentBase;

use strict;
use warnings;
use Genome;

my %TYPES = (
    EU => 0,
    AA => 1,
    All => 2,
);

class Genome::VariantReporting::Nhlbi::ComponentBase {
};

sub get_maf_for_entry {
    my $self = shift;
    my $entry = shift;
    my $population_code = shift;

    unless (defined $TYPES{$population_code}) {
        die $self->error_message("Bad population code: ($population_code)");
    }

    my $maf = $entry->info("MAF");
    if (!defined $maf or $maf eq '.') {
        return;
    }
    my @population_mafs = split(",", $maf);
    unless (@population_mafs == 3) {
        die $self->error_message("MAF in unexpected format: $maf");
    }
    return $population_mafs[$TYPES{$population_code}];
}

1;

