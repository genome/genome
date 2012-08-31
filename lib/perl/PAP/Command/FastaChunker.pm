package PAP::Command::FastaChunker;

use strict;
use warnings;

use Bio::Seq;
use Bio::SeqIO;
use English;
use File::Temp;

class PAP::Command::FastaChunker {
    is  => 'PAP::Command',
    has => [
        fasta_file => { 
            is => 'FilePath', 
            doc => 'fasta file name',
            is_input => 1,
        },
        chunk_size => { 
            is => 'Number', 
            doc => 'number of sequences per output file',
            is_input => 1,
        },
        fasta_files => { 
            is => 'ARRAY',  
            doc => 'array of fasta file names',
            is_optional => 1,
            is_output => 1,
        },
        lsf_queue => { 
            is_param => 1, 
            default_value => 'short',
        },
        lsf_resource => { 
            is_param => 1, 
            default_value => 'rusage[tmp=100]',
        },
    ],
};

sub sub_command_sort_position { 10 }

sub help_brief {
    "Chunk (split) a multi-fasta file into multiple smaller multi-fasta files";
}

sub help_synopsis {
    return help_brief();
}

sub help_detail {
    return help_brief();
}

sub execute {
    my $self = shift;

    my $input_file = $self->fasta_file;
    unless (-e $input_file) {
        die "No fasta file found at $input_file";
    }
    unless (-s $input_file) {
        die "Fasta file at $input_file has no size!";
    }

    my $chunk_size = $self->chunk_size;
    unless ($chunk_size > 0) {
        die "Chunk size must be greater than 0, $chunk_size is unacceptable!";
    }

    my $seq_in = Bio::SeqIO->new(
        -file => $input_file, 
        -format => 'Fasta'
    );
    unless ($seq_in) {
        die "Could not create Bio::SeqIO object for fasta file $input_file!";
    }

    my @output_files;
    my $chunk_fh;
    my $seq_out;
    my $seq_count = 0;
    
    $self->status_message("Starting fasta chunking, using chunk size $chunk_size and input fasta file $input_file");

    while (my $seq = $seq_in->next_seq()) {
        $seq_count++;
        if (($seq_count > $chunk_size) || (!defined($chunk_fh)))  {
            $seq_count = 1;

            ##FIXME: The temp dir location should not be hardcoded.  At least not here.
            $chunk_fh = File::Temp->new(
                'DIR'      => '/gscmnt/temp212/info/annotation/PAP_tmp',
                'SUFFIX'   => '.tmp',
                'TEMPLATE' => 'PAP_XXXXXXXX',
                'UNLINK'   => 0,
            );
            unless ($chunk_fh) {
                die "Could not create file handle for temp file!";
            }

            $seq_out = Bio::SeqIO->new(
                -fh => $chunk_fh, 
                -format => 'Fasta'
            );
            unless ($seq_out) {
                die "Could not create Bio::SeqIO object for output fasta chunk!";
            }

            push @output_files, $chunk_fh->filename();
            $self->status_message("Created fasta chunk " . $chunk_fh->filename);
        }

        $seq_out->write_seq($seq);
    }

    chmod 0666, @output_files;
    $self->fasta_files(\@output_files);

    # Okay, so this is a potential fix to a workflow problem. My hunch is that LSF jobs that finish too quickly
    # can't be tracked by workflow and leads to failures. This sleep should help prevent that. Maybe. We'll see.
    sleep(60) unless $ENV{UR_DBI_NO_COMMIT};
    $self->status_message("Done chunking!");
    return 1;
}
 
1;
