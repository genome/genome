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

sub requires_experts {
    return ('bam-readcount');
}

sub available_fields {
    return qw/coverage/;
}

sub interpret_entry {
    my ($self, $entry) = @_;

    my $readcount_entry = $self->get_readcount_entry($entry);
    return return_hash($entry, ".") unless $readcount_entry;

    return return_hash($entry, $readcount_entry->depth);
}

sub return_hash {
    my $entry = shift;
    my $coverage = shift;

    my %return_hash;
    for my $alt (@{$entry->{alternate_alleles}}) {
        $return_hash{$alt} = {
            coverage => $coverage
        };
    }
    return %return_hash;
}

1;
