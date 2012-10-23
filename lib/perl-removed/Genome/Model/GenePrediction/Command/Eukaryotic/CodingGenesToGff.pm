package Genome::Model::GenePrediction::Command::Eukaryotic::CodingGenesToGff;

use strict;
use warnings;
use Genome;
use Carp 'confess';

class Genome::Model::GenePrediction::Command::Eukaryotic::CodingGenesToGff {
    is => 'Genome::Command::Base',
    has => [
        prediction_directory => {
            is => 'DirectoryPath',
            doc => 'Directory containing gene predictions',
        },
        output_file => {
            is => 'FilePath',
            doc => 'File gff output is written to',
        },
    ],
};

sub help_detail {
    return <<EOS
Generates a gff file using gene predictions produced by the eukaryotic gene prediction pipeline
EOS
}

sub help_synopsis { help_detail() };
sub help_brief { help_detail() };

sub execute {
    my $self = shift;

    my $prediction_dir = $self->prediction_directory;
    confess "No prediction directory found at $prediction_dir!" unless -d $prediction_dir;

    my @coding_genes = Genome::Prediction::CodingGene->get(directory => $prediction_dir);
    $self->status_message("Found " . scalar @coding_genes . " coding gene predictions in $prediction_dir!");

    my $output_fh = IO::File->new($self->output_file, 'w');
    confess 'Could not get file handle for output file ' . $self->output_file unless $output_fh;

    for my $gene (@coding_genes) {
        my $transcript = $gene->transcript;
        my @exons = $transcript->exons;

        my ($start, $end) = (1, $transcript->spliced_length);
        my $desc = join(';', "ID=" . $gene->id, "Name=" . $gene->gene_name, "Target=" . 
            join(' ', $gene->gene_name, $start, $end, substr($gene->strand, 0, 1)));
        $output_fh->print(join("\t", $gene->sequence_name, $gene->source, 'match', $gene->start, 
                $gene->end, '.', substr($gene->strand, 0, 1), '.', $desc) . "\n");

        my $exon_count = 0;
        my $exon_length = 0;
        for my $exon (@exons) {
            my ($exon_start, $exon_end) = ($exon_length + 1, $exon_length + $exon->length);
            $exon_length += $exon->length;

            my $exon_desc = join(';', "ID=" . $gene->id . ":cds:$exon_count", "Parent=" . $gene->id, "Name=" . $gene->gene_name, 
                "Target=" . join(' ', $gene->id, $exon_start, $exon_end, substr($exon->strand, 0, 1)));
            $output_fh->print(join("\t", $exon->sequence_name, $gene->source, 'match_part', $exon->start,
                    $exon->end, $exon->score, substr($exon->strand, 0, 1), '.', $exon_desc) . "\n");
            $exon_count++;
        }
    }

    $output_fh->close;
    $self->status_message("Successfully completed conversion to gff, output file at " . $self->output_file);
    return 1;
}

1;

