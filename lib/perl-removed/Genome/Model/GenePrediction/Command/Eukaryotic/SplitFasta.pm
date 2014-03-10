package Genome::Model::GenePrediction::Command::Eukaryotic::SplitFasta;

use strict;
use warnings;

use Genome;
use Carp 'confess';
use File::Path 'make_path';
use POSIX 'ceil';
use Bio::SeqIO;

class Genome::Model::GenePrediction::Command::Eukaryotic::SplitFasta {
    is => 'Command',
    has => [
        fasta_file => {
            is => 'Path',
            is_input => 1,
            doc => 'Fasta file to be split up',
        },
        output_directory => {
            is => 'Path',
            is_input => 1,
            is_output => 1,
            doc => 'Directory in which split fastas are placed',
        },
    ],
    has_optional => [
        max_chunks => {
            is => 'Number',
            is_input => 1,
            default => 100,
            doc => 'Maximum number of chunks allowed, if this limit is reached then max bases per files is increased as needed',
        },
        max_bases_per_file => {
            is => 'Number',
            is_input => 1,
            default => 5000000,
            doc => 'Maximum number of bases allowed in each split fasta file',
        },
        fasta_files => {
            is => 'ARRAY',
            is_output => 1,
            doc => 'An array of split fasta files',
        },
        genome_size => {
            is => 'Number',
            is_output => 1,
            doc => 'Total number of bases in given fasta file, might be a useful metric',
        },
    ],
};

sub help_brief {
    return "Splits up a fasta into several chunks";
}

sub help_detail {
    return <<EOS
Given a fasta file, creates several smaller chunks in the given output_directory. Each fasta
chunk is no larger than the given max_bases_per_file parameter.
EOS
}

sub execute {
    my $self = shift;
    my $fasta_file_path = $self->fasta_file;
    unless (-e $fasta_file_path) {
        confess "$fasta_file_path does not exist!";
    }
    unless (-s $fasta_file_path) {
        confess "$fasta_file_path does not have size!";
    }

    my $fasta_file = Bio::SeqIO->new(
        -file => $fasta_file_path,
        -format => 'Fasta',
    );

    # The only way (I know of) to determine total sequence in the fasta file is to iterate through it. I hate
    # to iterate through the same file twice...
    my $total_bases = 0;
    while (my $sequence = $fasta_file->next_seq()) {
        $total_bases += $sequence->length;
    }

    # Now determine the chunk size... chunk size begins at the value provided by the user. If the total number of 
    # chunks (using this chunk size) exceeds the total allowable chunks, then chunk size is increased.
    my $chunk_size = $self->max_bases_per_file;
    my $total_chunks = ceil($total_bases / $chunk_size);
    unless ($total_chunks <= $self->max_chunks) {
        $chunk_size = ceil(($total_bases / $self->max_chunks) * 1.1);  # This is to allow some wiggle room, since each fasta won't
                                                                       # have exactly the maximum number the bases
        $self->debug_message("Increasing fasta chunk size from " . $self->max_bases_per_file . " to $chunk_size. This prevents the " .
            "total number of split fastas from exceeding " . $self->max_chunks . " (would have been $total_chunks).");
    }

    # Now split up the fasta file!
    my @filenames;
    my $current_fasta;
    my $counter = 0;
    my $current_chunk_size = 0;
    my $output_directory = $self->output_directory;

    $self->debug_message("Creating smaller fasta files in $output_directory containing sequence " .
        "from $fasta_file_path and each having no more than $chunk_size bases");

    unless (-d $output_directory) {
        my $rv = make_path($output_directory);
        confess "Could not make directory $output_directory!" unless defined $rv and $rv == 1;
    }

    # Need to reset the iterator before going through the fasta file again...
    $fasta_file = Bio::SeqIO->new(
        -file => $fasta_file_path,
        -format => 'Fasta',
    );

    # Now write sequence to each fasta chunk, making sure that the chunk size is never exceeded
    while (my $sequence = $fasta_file->next_seq()) {
        my $length = $sequence->length;

        if (not defined $current_fasta or ($current_chunk_size + $length) > $chunk_size) {
            my $filename = $output_directory . "/fasta_$counter";
            $current_fasta = Bio::SeqIO->new(
                -file => ">$filename",
                -format => 'Fasta',
            );

            $counter++;
            $current_chunk_size = 0;
            push @filenames, $filename;
        }

        $current_fasta->write_seq($sequence);
        $current_chunk_size += $length;
    }

    $self->fasta_files(\@filenames);
    $self->genome_size($total_bases);
    $self->debug_message("Created $counter fasta files in $output_directory.");
    $self->debug_message("Altogether, there are $total_bases bases of sequence!");
    return 1;
}
1;

