package Genome::Model::Tools::GenePredictor::Trnascan;

use strict;
use warnings;

use Genome;
use File::Temp;
use IO::File;
use Carp 'confess';
use File::Path 'make_path';

class Genome::Model::Tools::GenePredictor::Trnascan {
    is => 'Genome::Model::Tools::GenePredictor',
    has_optional => [
        domain => {
            is => 'Text',
            is_input => 1,
            valid_values => ['archaeal', 'bacterial', 'eukaryotic'],
            default => 'eukaryotic',
        },
        trnascan_install_path => {
            is => 'Path',
            default => '/gsc/bin/tRNAscan-SE',
            doc => 'Path to program executable',
        },
    ],
};

sub help_brief {
    return "Runs Trnascan on the provided fasta file";
}

sub help_synopsis {
    return "Runs Trnascan on the provided fasta file";
}

sub help_detail {
    return <<EOS
Runs tRNAsca on the provided fasta file, places raw output and prediction
output into provided directories
EOS
}

sub execute {
    my $self = shift;

    if ($self->skip_execution) {
        $self->status_message("Skip execution flag is set, exiting.");
        return 1;
    }

    $self->debug_message("Running trnascan on sequence in " . $self->fasta_file);

    unless (-d $self->raw_output_directory) {
        my $mk_rv = make_path($self->raw_output_directory);
        confess "Could not make raw ouput directory at " . $self->raw_output_directory unless defined $mk_rv and $mk_rv;
    }

    # Need a unique file name for raw output
    my $raw_output_fh = File::Temp->new(
        DIR => $self->raw_output_directory,
        TEMPLATE => "trnascan_raw_output_XXXXX",
        CLEANUP => 0,
        UNLINK => 0,
    );
    my $raw_output_file = $raw_output_fh->filename;
    $raw_output_fh->close;
    chmod(0666, $raw_output_file);
    $self->debug_message("Raw output being written to $raw_output_file");

    # Construct command and parameters/switches
    my @params;
    push @params, $self->fasta_file;
    push @params, "-B " if $self->domain eq 'bacterial';
    push @params, "-A " if $self->domain eq 'archaeal';
    push @params, "> $raw_output_file ";
    push @params, "2> $raw_output_file.error ";
   
    my $cmd = join(" ", $self->trnascan_install_path, @params);
    $self->debug_message("Preparing to run Trnascan-SE: $cmd");
    
    # FIXME Replace with Genome::Sys->shellcmd
    my $rv = system($cmd);
    confess 'Trouble executing tRNAscan!' unless defined $rv and $rv == 0;
    $self->debug_message("Done executing tRNAscan, now parsing output");

    # Parse output and create UR objects
    $raw_output_fh = IO::File->new($raw_output_file, 'r');
    for (1..3) { $raw_output_fh->getline };  # First three lines are headers
    while (my $line = $raw_output_fh->getline) {
        chomp $line;
        my ($seq_name, $trna_num, $begin, $end, $type, $codon, $intron_begin, $intron_end, $score) = split(/\s+/, $line);
        $self->debug_message("Parsing $seq_name");

        my $strand = 1;
        $strand = -1 if $begin > $end;
        ($begin, $end) = ($end, $begin) if $begin > $end;

        my $sequence = $self->get_sequence_by_name($seq_name);
        confess "Couldn't get sequence $seq_name!" unless $sequence;
        my $seq_string = $sequence->subseq($begin, $end);

        my $rna_gene = Genome::Prediction::RNAGene->create(
            directory => $self->prediction_directory,
            gene_name => $seq_name . '.t' . $trna_num,
            description => $type,
            start => $begin,
            end => $end,
            strand => $strand,
            source => 'trnascan',
            score => $score,
            sequence_name => $seq_name,
            sequence_string => $seq_string,
            codon => $codon,
            amino_acid => $type,
        );
    }

    $self->debug_message("Parsing done, getting locks and committing!");
    my @locks = $self->lock_files_for_predictions(qw/ Genome::Prediction::RNAGene /);
    UR::Context->commit;
    $self->release_prediction_locks(@locks);

    $self->status_message("Commit done, locks released, trnascan successfully completed!");
    return 1;
}

1;
