package Genome::VariantReporting::BamReadcount::CoverageInterpreter;

use strict;
use warnings FATAL => 'all';
use Genome;
use List::Util qw/first/;

class Genome::VariantReporting::BamReadcount::CoverageInterpreter {
    is => ['Genome::VariantReporting::Framework::Component::Interpreter', 'Genome::VariantReporting::BamReadcount::ComponentBase'],
};

sub name {
    return 'coverage';
}

sub requires_annotations {
    return ('bam-readcount');
}

sub available_fields {
    return qw/coverage/;
}

sub _interpret_entry {
    my ($self, $entry, $passed_alt_alleles) = @_;

    my $readcount_entry = $self->get_readcount_entry($entry);
    if ($readcount_entry) {
        my %return_hash;
        for my $alt (@$passed_alt_alleles) {
            $return_hash{$alt} = {
                coverage => $readcount_entry->depth,
            };
        }
        return %return_hash;
    } else {
        return $self->null_interpretation($passed_alt_alleles);
    }
}

1;
