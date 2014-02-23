package Genome::Model::GenePrediction::Command::Eukaryotic::MaskRnaSequence;

use strict;
use warnings;

use Genome;
use Carp 'confess';
use File::Temp;
use File::Basename;
use File::Spec;
use Genome::Sys;

use Bio::SeqIO;

class Genome::Model::GenePrediction::Command::Eukaryotic::MaskRnaSequence {
    is => 'Command',
    has => [
        prediction_directory => {
            is => 'Path',
            is_input => 1,
            doc => 'Directory containing predictions, used to grab RNA predictions',
        },
        fasta_file => {
            is => 'Path',
            is_input => 1,
            doc => 'Fasta file with sequence that needs masking',
        },
    ],
    has_optional => [
        skip_if_no_rna_file => {
            is => 'Boolean',
            is_input => 1,
            doc => 'If this is set, not finding an rna file is not a fatal error and the masked fasta is the same as the input fasta',
            default => 1,
        },
        masked_fasta_file => {
            is => 'Path',
            is_input => 1,
            is_output => 1,
            doc => 'Fasta file that will contain masked sequence, defaults to fasta_file with a random suffix',
        },
    ],
};

sub help_brief { 
    return "Masks out the sequence associated with RNA predictions from a fasta file";
}

sub help_detail { 
    return "Masks out the sequence associated with RNA predictions from a fasta file";
}

sub execute {
    my $self = shift;

    $self->debug_message("Starting rna masking command, masking sequences in " . $self->fasta_file);

    unless (-e $self->fasta_file) {
        confess 'No fasta file found at ' . $self->fasta_file . '!';
    }

    unless (defined $self->masked_fasta_file) {
        my $full_fasta_path = File::Spec->rel2abs($self->fasta_file);
        my ($fasta_name, $fasta_dir) = fileparse($full_fasta_path);
        my $masked_fh = File::Temp->new(
            TEMPLATE => $fasta_name . ".rna_masked_XXXXXX",
            DIR => $fasta_dir,
            CLEANUP => 0,
            UNLINK => 0,
        );
        chmod(0666, $masked_fh->filename);
        $self->masked_fasta_file($masked_fh->filename);
        $masked_fh->close;
    }

    $self->debug_message("Masked sequences are being written to fasta file at " . $self->masked_fasta_file);

    # If no RNA genes are found, there won't be an rna file. This is a valid situation, so if no RNA gene file
    # and the skip flag is set, just copy the input fasta to the masked fasta location. To figure out the path
    # to the RNA file, need to use the data source and call its file resolver method.
    my $rna_data_source = Genome::Prediction::RNAGene->__meta__->data_source;
    my $file_resolver = $rna_data_source->can('file_resolver');
    my $rna_file = $file_resolver->($self->prediction_directory);
    unless (-e $rna_file) {
        if ($self->skip_if_no_rna_file) {
            $self->debug_message("No rna predictions file found at $rna_file" .
                " and skip_if_no_rna_file flag is set to true. Assuming that no rna predictions were" .
                " made and setting input fasta " . $self->fasta_file . " as masked output fasta!");

            my $cp_rv = Genome::Sys->copy_file($self->fasta_file, $self->masked_fasta_file);
            unless (defined $cp_rv and $cp_rv) {
                confess 'Could not copy input fasta ' . $self->fasta_file . ' to ' . $self->masked_fasta_file . "!";
            }

            $self->debug_message("Copy successful! Masked fasta file now at " . $self->masked_fasta_file . ", exiting!");
            return 1;
        }
        else {
            confess "No rna prediction file found at $rna_file";
        }
    }
    
    # This pre-loads all the predictions, which makes grabbing predictions per sequence faster below
    $self->debug_message("Found RNA predictions in $rna_file, now loading them into memory!");
    my @rna_predictions = Genome::Prediction::RNAGene->get(
        directory => $self->prediction_directory
    );

    my $masked_fasta = Bio::SeqIO->new(
        -file => '>' . $self->masked_fasta_file,
        -format => 'Fasta',
    );

    my $fasta = Bio::SeqIO->new(
        -file => $self->fasta_file,
        -format => 'Fasta',
    );

    $self->debug_message("Masking sequences in " . $self->fasta_file . " and writing to " . $self->masked_fasta_file);

    # Iterate through every sequence in the fasta and find the predictions associated with each one,
    # then mask out sequence associated with an rna prediction
    while (my $seq = $fasta->next_seq()) {
        $self->debug_message("Working on sequence " . $seq->display_id());
        my $seq_id = $seq->display_id();
        my $seq_length = $seq->length();
        my $seq_string = $seq->seq();
        
        my @predictions = Genome::Prediction::RNAGene->get(
            directory => $self->prediction_directory,
            sequence_name => $seq_id,
        );

        for my $prediction (@predictions) {
            my $start = $prediction->start;
            my $end = $prediction->end;

            # Make sure that start is less than end and within the bounds of the sequence
            ($start, $end) = ($end, $start) if $start > $end;
            $start = 1 if $start < 1;
            $end = $seq_length if $end > $seq_length;
            my $prediction_length = ($end - $start) + 1;

            # Any sequence within the rna prediction is replaced with an N
            substr($seq_string, $start - 1, $prediction_length, 'N' x $prediction_length);
        }

        my $masked_seq = Bio::Seq->new(
            -display_id => $seq_id,
            -seq => $seq_string
        );
        $masked_fasta->write_seq($masked_seq);
    }

    $self->debug_message("All sequences have been masked!");
    return 1;
}

1;

