package Genome::Model::Tools::Analysis::ResolveIgvReferenceName;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Analysis::ResolveIgvReferenceName {
    is => 'Command::V2',
    has_input => [
        reference_name => {
            is => 'Text',
            doc => "name of the reference to translate to the IGV reference name",
            example_values => ["GRCh37-lite-build37"],
        },
    ],
    has_optional_output => [
        igv_reference_name => {
            is => 'Text',
            is_transient => 1,
            default => 'reference',
            doc => "name of the IGV reference derived from the input reference name",
            example_values => ["reference_build36","b37","mm9"],
        },
    ],
};

sub help_detail {
  return <<HELP;
This tool will resolve a given GMS reference sequence name to an IGV reference name.
HELP
}

sub _doc_authors {
  return <<AUTHS;
 Susanna Siebert
AUTHS
}

sub execute {
    my $self = shift;

    my $igv_reference_name = resolve_igv_reference_name($self->reference_name);

    if (defined($igv_reference_name)) {
        $self->status_message("Found IGV reference name for (%s): $igv_reference_name", $self->reference_name);
        $self->igv_reference_name($igv_reference_name);
    }
    else {
        $self->status_message("Couldn't find IGV reference name for (%s). Defaulting to (%s)", $self->reference_name, $self->igv_reference_name);
    }

    return 1;
}

sub resolve_igv_reference_name {
    my $reference_name = shift;

    my %igv_reference_name_of = (
        'GRCh37-lite-build37'   => 'b37',
        'NCBI-human-build36'    => 'hg18',
        'g1k-human-build37'     => 'b37',
        'UCSC-mouse-buildmm9'   => 'mm9',
    );

    return $igv_reference_name_of{$reference_name} || undef;
}

1;
