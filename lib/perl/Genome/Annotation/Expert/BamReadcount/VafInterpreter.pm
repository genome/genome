package Genome::Annotation::Expert::BamReadcount::VafInterpreter;

use strict;
use warnings;
use Genome;
use Genome::Annotation::Expert::BamReadcount::VafCalculator;

class Genome::Annotation::Expert::BamReadcount::VafInterpreter {
    is => ['Genome::Annotation::Interpreter::Base', 'Genome::Annotation::Expert::BamReadcount::ComponentBase'],
};

sub name {
    return 'vaf';
}

sub requires_experts {
    return ('bam-readcount');
}

sub available_fields {
    return qw/
        vaf
    /;
}

sub interpret_entry {
    my $self = shift;
    my $entry = shift;
    my $passed_alt_alleles = shift;

    my %return_values;
    my @sample_alt_alleles = sort $entry->alt_bases_for_sample($self->sample_index($entry->{header}));
    my %vafs = Genome::Annotation::Expert::BamReadcount::VafCalculator::calculate_vaf_for_multiple_alleles(
        $self->get_readcount_entry($entry), \@sample_alt_alleles);
    for my $allele (@$passed_alt_alleles) {
        my $vaf;
        if (defined $vafs{$allele}) {
            $vaf = $vafs{$allele};
        }
        else {
            $vaf = undef;
        }
        $return_values{$allele} = {
            vaf => $vaf,
        }
    }

    return %return_values;
}

1;

