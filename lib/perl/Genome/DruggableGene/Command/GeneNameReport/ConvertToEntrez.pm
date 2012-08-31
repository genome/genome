package Genome::DruggableGene::Command::GeneNameReport::ConvertToEntrez;

use strict;
use warnings;
use Genome;

class Genome::DruggableGene::Command::GeneNameReport::ConvertToEntrez {
    is => 'Genome::Command::Base',
    has => [
        gene_identifiers => {
            is => 'Text',
            shell_args_position => 1,
            doc => 'Gene identifiers to convert to entrez',
            is_many => 1,
        },
    ],
};

sub help_brief {
    'Translate a gene identifier to one or more Genome::DruggableGene::GeneNameReports';
}

sub help_synopsis {
    'genome druggable-gene gene-name-report convert-to-entrez ARK1D1 ARK1E1';
}

sub help_detail { help_brief() }

sub execute {
    my $self = shift;
    for my $gene_identifier ($self->gene_identifiers) {
        my ($entrez_genes, $intermediate_genes) = Genome::DruggableGene::GeneNameReport->convert_to_entrez($gene_identifier);
        my (@entrez, @intermediate);
        @entrez = @{$entrez_genes->{$gene_identifier}} if $entrez_genes->{$gene_identifier};
        @intermediate  = @{$intermediate_genes->{$gene_identifier}} if $intermediate_genes->{$gene_identifier};
        if(@entrez){
            $self->status_message($gene_identifier . " as entrez is:\n");
            $self->status_message($_->name . "\n") for (@entrez);
        }
        if(@intermediate) {
            $self->status_message($gene_identifier . " as intermediate is:\n");
            $self->status_message($_->name . "\n") for (@intermediate);
        }
    }
    return 1;
}

1;
