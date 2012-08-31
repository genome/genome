package Genome::Model::Build::GenePrediction::Eukaryotic;

use strict;
use warnings;
use Genome;
use Carp 'confess';

use Bio::SeqIO;

class Genome::Model::Build::GenePrediction::Eukaryotic {
    is => 'Genome::Model::Build::GenePrediction',
};

# >75% of this comes from split fasta directories, which get cleaned
# up when the build successfully completes
sub calculate_estimated_kb_usage {
    return 2_048_000;
}

sub sorted_fasta_file {
    my $self = shift;
    return $self->data_directory . '/fasta.sorted';
}

sub workflow_name {
    my $self = shift;
    return 'eukaryotic gene prediction ' . $self->build_id;
}

sub repeat_masker_ace_file {
    my $self = shift;
    return $self->data_directory . "/repeat_masker.ace";
}

sub overly_masked_sequence_fasta_file {
    my $self = shift;
    return $self->data_directory . "/repeat_masker.overly_masked.fasta";
}

sub repeat_masker_gff_file {
    my $self = shift;
    return $self->data_directory . "/repeat_masker.gff";
}

sub predictions_ace_file {
    my $self = shift;
    return $self->data_directory . "/predictions.ace";
}

sub rna_predictions_ace_file {
    my $self = shift;
    return $self->data_directory . "/rna_predictions.ace";
}

sub log_directory {
    my $self = shift;
    return $self->data_directory . '/logs/';
}

sub split_fastas_output_directory {
    my $self = shift;
    return $self->data_directory . '/split_fastas/';
}

sub raw_output_directory {
    my $self = shift;
    return $self->data_directory . '/raw_predictor_output/';
}

sub prediction_directory {
    my $self = shift;
    return $self->data_directory;
}

sub dirs_ignored_by_diff {
    return qw(
        logs/
        raw_predictor_output/
    );
}

sub files_ignored_by_diff {
    return qw(
        split_fastas/(.*)rna_masked(.*)
        build.xml
        reports/Build_Initialized/report.xml
        reports/Build_Succeeded/report.xml
        repeat_masker.gff
    );
}

# Returns a list of sequence names in the assembly contigs file.
# TODO If this is a common request, may want to consider storing the
# sequence names somewhere
sub sequences { 
    my $self = shift;
    my $model = $self->model;
    my $fasta = $model->assembly_contigs_file;
    confess "No fasta file found at $fasta" unless -e $fasta;

    my $seq_obj = Bio::SeqIO->new(
        -file => $fasta,
        -format => 'Fasta',
    );
    confess "Could not create Bio::SeqIO object for $fasta" unless $seq_obj;

    my @sequences;
    while (my $seq = $seq_obj->next_seq) {
        push @sequences, $seq->display_id;
    }

    return \@sequences;
}
    
1;

