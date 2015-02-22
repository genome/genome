package Genome::Model::Tools::Gtf::LimitByBiotype;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Gtf::LimitByBiotype {
    is => 'Genome::Model::Tools::Gtf::Base',
    has => [
        transcript_info_tsv_file => {
            is => 'Text',
            doc => 'A tsv file generated through the Ensembl API that contains at least: gene_biotype, ensembl_transcript_id',
        },
        output_gtf_file => {
            is => 'Text',
            doc => 'The output GTF file.',
        },
        gene_biotypes => {
            doc => 'An array reference or comma-delimited string of biotypes.',
        },
    ],
};

sub execute {
    my $self = shift;

    my $biotypes = $self->gene_biotypes;
    my @biotypes;
    if (ref($biotypes) eq 'ARRAY') {
        @biotypes = @$biotypes;
    } else {
        @biotypes = split(',', $biotypes);
    }
    my %biotypes = map { $_ => 1 } @biotypes;

    my $transcript_reader = Genome::Utility::IO::SeparatedValueReader->create(
        input => $self->transcript_info_tsv_file,
        separator => "\t",
    );
    unless ($transcript_reader) {
        die('Failed to load transcript info: '. $self->transcript_info_tsv_file);
    }
    my %keeper_transcript_ids;
    while (my $data = $transcript_reader->next) {
        if ($biotypes{$data->{gene_biotype}}) {
            $keeper_transcript_ids{$data->{ensembl_transcript_id}} = 1;
        }
    }
    my @keeper_transcript_ids = keys %keeper_transcript_ids;
    unless (Genome::Model::Tools::Gtf::Limit->execute(
        input_gtf_file => $self->input_gtf_file,
        output_gtf_file => $self->output_gtf_file,
        id_type => 'transcript_id',
        ids => \@keeper_transcript_ids,
    )->result) {
        die('Failed to limit genes to only defined biotypes!');
    }
    return 1;
}
