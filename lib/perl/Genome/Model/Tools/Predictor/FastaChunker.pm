package Genome::Model::Tools::Predictor::FastaChunker;

use strict;
use warnings;

use Genome;
use Bio::Seq;
use Bio::SeqIO;
use English;
use File::Temp;

# TODO There's another fasta chunk command used by EGAP... perhaps use that and remove this one?
class Genome::Model::Tools::Predictor::FastaChunker {
    is  => 'Command::V2',
    has => [
        input_fasta_file => { 
            is => 'FilePath',
            is_input => 1,
            doc => 'fasta file to be split',
        },
        chunk_size  => { 
            is => 'Number', 
            is_input => 1,
            doc => 'max number of sequences per fasta chunk',
        },
        output_directory => {
            is => 'DirectoryPath',
            is_optional => 1,
            is_input => 1,
            default => '/gscmnt/temp212/info/annotation/PAP_tmp',
            doc => 'Directory into which fasta chunks are placed',
        },
        lsf_queue => { 
            is_param => 1, 
            default_value => $ENV{GENOME_LSF_QUEUE_SHORT},
        },
        lsf_resource => { 
            is_param => 1, 
            default_value => 'rusage[tmp=100]',
        },
        fasta_files => { 
            is => 'ARRAY',  
            is_optional => 1,
            is_output => 1,
            doc => 'array of fasta files produced via chunking',
        },
    ],
};

sub help_brief {
    return 'split a fasta file into multiple smaller files';
}

sub help_detail {
    return help_brief();
}

sub execute {
    my $self = shift;

    my $input_file = $self->input_fasta_file;
    unless (-e $input_file) {
        die "No file found at $input_file!";
    }
    unless (-s $input_file) {
        die "File $input_file has no size!";
    }

    my $chunk_size = $self->chunk_size;
    unless ($chunk_size > 0) {
        die "Chunk size must be larger than 0, not $chunk_size!";
    }

    my $seq_in = Bio::SeqIO->new(-file => $input_file, -format => 'Fasta');
    unless ($seq_in) {
        die "Could not create reader for fasta file $input_file";
    }

    $self->status_message("Splitting up fasta file " . $self->input_fasta_file . " into chunks containing no more than $chunk_size sequences!");
    my $chunk_fh;
    my $seq_out;
    my $seq_count = 0;
    my @output_files;
    while (my $seq = $seq_in->next_seq()) {
        $seq_count++;
        if (($seq_count > $chunk_size) || (!defined($chunk_fh)))  {
            $seq_count = 1;
            $chunk_fh = File::Temp->new(
                'DIR'      => $self->output_directory,
                'SUFFIX'   => '.tmp',
                'TEMPLATE' => 'PAP_XXXXXXXX',
                'UNLINK'   => 0,
            );
            chmod 0666, $chunk_fh->filename();
            $seq_out = Bio::SeqIO->new(-fh => $chunk_fh, -format => 'Fasta');
            push @output_files, $chunk_fh->filename();
        }
        $seq_out->write_seq($seq);
    }

    $self->fasta_files(\@output_files);
    $self->status_messages("Done splitting up input fasta, created " . scalar @output_files . " chunks");
    return 1;
}
 
1;
