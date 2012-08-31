#$Id: FastaChunker.pm 39214 2008-09-30 19:29:40Z mjohnson $

package Genome::Model::Tools::Fasta::Chunker;

use strict;
use warnings;

use Genome;
use Workflow;

use Bio::Seq;
use Bio::SeqIO;

use English;
use File::Temp;


class Genome::Model::Tools::Fasta::Chunker {
    is  => ['Genome::Model::Tools::Fasta'],
    has_input => [                                  
        chunk_size  => { is => 'SCALAR', doc => 'number of sequences per output file', default => 10},  
        tmp_dir     => { is => 'SCALAR', doc => 'directory for saving temporary file chunks', default => Genome::Sys->create_temp_directory},
    ],
    has_output => [
        file_chunks => { is => 'ARRAY',  doc => 'array of fasta file names', is_optional => 1 },
    ],
};

sub sub_command_sort_position { 10 }

sub help_brief {
    "Chunk (split) a multi-fasta file into multiple smaller multi-fasta files";
}

sub help_synopsis {
    return <<"EOS"
EOS
}

sub help_detail {
    return <<"EOS"
Need documenation here.
EOS
}

sub execute {
    my $self = shift;
    my $input_file      = $self->fasta_file();
    unless (-s $input_file)
    {
        return $self->error_message("$input_file cannot be chunked (file is missing or empty)");
    }
    my $chunk_size      = int($self->chunk_size());     
    unless ($chunk_size > 0)
    {
       return $self->error_message("Integer value required for chunk size ($chunk_size)");
    }
    my $tmp_dir         = $self->tmp_dir();
    my $unique_template = $self->get_unique_id();
    my @output_files = ( );

    my $seq_in = Bio::SeqIO->new(-file => $input_file, -format => 'Fasta');

    my $chunk_fh;
    my $seq_out;
    
    my $seq_count = 0;
    my $tmp = 0;
    while (my $seq = $seq_in->next_seq()) {

        $seq_count++;
        
        if (($seq_count > $chunk_size) || (!defined($chunk_fh)))  {
        
            $seq_count = 1;

            $chunk_fh = File::Temp->new(
                                        'DIR'      => $tmp_dir,
                                        'SUFFIX'   => '.tmp',
                                        'TEMPLATE' => $tmp++ . "_" . $unique_template,
                                        'UNLINK'   => 0,
                                       );

            $seq_out = Bio::SeqIO->new(-fh => $chunk_fh, -format => 'Fasta');
        
            push @output_files, $chunk_fh->filename();

        }

        $seq_out->write_seq($seq);
    }

    $self->file_chunks(\@output_files);

    return 1;

}
 
sub get_unique_id()
{
    my $i = rand;
    return substr ($i,2,5) . 'XXXX';
}

1;
