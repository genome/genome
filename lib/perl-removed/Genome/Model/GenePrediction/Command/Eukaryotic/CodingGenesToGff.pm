package Genome::Model::GenePrediction::Command::Eukaryotic::CodingGenesToGff;

use strict;
use warnings;
use Genome;
use Carp 'confess';

class Genome::Model::GenePrediction::Command::Eukaryotic::CodingGenesToGff {
    is => 'Genome::Command::Base',
    has_input => [
        prediction_directory => {
            is => 'DirectoryPath',
            doc => 'Directory containing gene predictions',
            is_output => 1,
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
    $self->debug_message("Found " . scalar @coding_genes . " coding gene predictions in $prediction_dir!");

    my $coding_gene_reader = Genome::Utility::IO::SeparatedValueReader->create(
        input => $prediction_dir.'/coding_genes.csv',
        headers => Genome::DataSource::Predictions::CodingGenes->column_order,
    );

    my $transcript_reader = Genome::Utility::IO::SeparatedValueReader->create(
        input => $prediction_dir .'/transcripts.csv',
        headers => Genome::DataSource::Predictions::Transcripts->column_order,
    );

    my $exon_reader = Genome::Utility::IO::SeparatedValueReader->create(
        input => $prediction_dir .'/exons.csv',
        headers => Genome::DataSource::Predictions::Exons->column_order,
    );
    
    my $output_fh = IO::File->new($self->output_file, 'w');
    confess 'Could not get file handle for output file ' . $self->output_file unless $output_fh;
    my $exon = $exon_reader->next;
    while( my $gene = $coding_gene_reader->next) {
        my $transcript = $transcript_reader->next;
        unless ($transcript->{coding_gene_name} eq $gene->{gene_name}) {
            die('Coding Gene and Transcript files are out of order!');
        }

        my @exons;
        while ($exon && ($exon->{transcript_name} eq $transcript->{transcript_name}) ) {
            push @exons, $exon;
            $exon = $exon_reader->next;
        }
        unless (@exons) {
            die ('Coding Gene, Transcript and Exon files are out of order!');
        }
        my ($start, $end) = (1, $self->spliced_length(\@exons));
        my $gene_strand = substr($gene->{strand}, 0, 1);
        unless ($gene_strand =~ /^[\+\-]$/) {
            if ($gene_strand == 1) {
                $gene_strand = '+';
            } else {
                die('Unrecognized strand '. $gene_strand .' for gene '. $gene->{gene_name});
            }
        }
        my $desc = join(';', "ID=" . $gene->{gene_name}, "Name=" . $gene->{gene_name}, "Target=" . 
            join(' ', $gene->{gene_name}, $start, $end, $gene_strand));
        $output_fh->print(join("\t", $gene->{sequence_name}, $gene->{source}, 'match', $gene->{start}, 
                $gene->{end}, '.', $gene_strand, '.', $desc) . "\n");

        my $exon_count = 0;
        my $exon_length = 0;
        my $sort_order = 'sort_ascending';
        if ($gene_strand eq '-') {
            $sort_order = 'sort_descending';
        }
        for my $exon (sort $sort_order @exons) {
            my ($exon_start, $exon_end) = ($exon_length + 1, $exon_length + $self->exon_length($exon));
            $exon_length += $self->exon_length($exon);
            my $exon_strand = substr($exon->{strand}, 0, 1);
            unless ($exon_strand =~ /^[\+\-]$/) {
                if ($exon_strand == 1) {
                    $exon_strand = '+';
                } else {
                    die('Unrecognized strand '. $exon_strand .' for exon '. $exon->{exon_name});
                }
            }
            my $exon_desc = join(';', "ID=" . $gene->{gene_name} . ":cds:$exon_count", "Parent=" . $gene->{gene_name}, "Name=" . $gene->{gene_name}, 
                "Target=" . join(' ', $gene->{gene_name}, $exon_start, $exon_end, $exon_strand));
            $output_fh->print(join("\t", $exon->{sequence_name}, $gene->{source}, 'match_part', $exon->{start},
                    $exon->{end}, $exon->{score}, $exon_strand, '.', $exon_desc) . "\n");
            $exon_count++;
        }
    }

    $output_fh->close;
    $self->debug_message("Successfully completed conversion to gff, output file at " . $self->output_file);
    return 1;
}

sub exon_length {
    my ($self,$exon) = @_;
    return abs($exon->{start} - $exon->{end}) + 1;
}

sub spliced_length {
    my ($self,$exons) = @_;
    my $sum = 0;
    for my $exon (@$exons) {
        $sum += $self->exon_length($exon);
    }
    return $sum;
}

sub sort_ascending {
    $a->{start} <=> $b->{start};
}

sub sort_descending {
    $b->{start} <=> $a->{start};
}

1;

