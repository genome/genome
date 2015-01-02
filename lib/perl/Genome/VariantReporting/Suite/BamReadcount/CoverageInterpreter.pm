package Genome::VariantReporting::Suite::BamReadcount::CoverageInterpreter;

use strict;
use warnings FATAL => 'all';
use Genome;
use List::Util qw/first/;

class Genome::VariantReporting::Suite::BamReadcount::CoverageInterpreter {
    is => ['Genome::VariantReporting::Framework::Component::Interpreter', 'Genome::VariantReporting::Suite::BamReadcount::ComponentBase'],
    doc => 'Determine the read depth at the variant position',
};

sub name {
    return 'coverage';
}

sub requires_annotations {
    return ('bam-readcount');
}

sub field_descriptions {
    return (
        coverage => 'Read depth at position'
    );
}

sub _interpret_entry {
    my ($self, $entry, $passed_alt_alleles) = @_;

    my $readcount_entries = $self->get_readcount_entries($entry);
    if ($readcount_entries) {
        my %return_hash;
        for my $alt (@$passed_alt_alleles) {
            my $readcount_entry = $readcount_entries->{$alt};
            if (defined $readcount_entry) {
                $return_hash{$alt} = {
                    coverage => $readcount_entry->depth,
                };
            }
            else {
                $return_hash{$alt} = { map { $_ => $self->interpretation_null_character } $self->available_fields };
            }
        }
        return %return_hash;
    } else {
        return $self->null_interpretation($passed_alt_alleles);
    }
}

1;
