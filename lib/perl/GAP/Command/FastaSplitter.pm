package GAP::Command::FastaSplitter;

use strict;
use warnings;

use Workflow;

use Bio::Seq;
use Bio::SeqIO;

use English;
use File::Temp;


class GAP::Command::FastaSplitter {
    is  => ['GAP::Command'],
    has => [
        fasta_file  => { is => 'SCALAR', doc => 'fasta file name'                             },
        chunk_size  => { is => 'SCALAR', doc => 'number of sequences per output file'         },
        fasta_files => { is => 'ARRAY',  doc => 'array of fasta file names', is_optional => 1 },
    ],
};

operation GAP::Command::FastaSplitter {
    input  => [ 'fasta_file', 'chunk_size' ],
    output => [ 'fasta_files'              ],
};

sub sub_command_sort_position { 10 }

sub help_brief {
    "Split a multi-fasta file into multiple smaller fasta files";
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


    my $input_file   = $self->fasta_file();     ##FIXME:  Should verify input_file is not empty
    my $chunk_size   = $self->chunk_size();     ##FIXME:  Should verify chunk_size > 0
    my @output_files = ( );

    if ($input_file =~ /\.bz2$/) {
        $input_file = "bzcat $input_file |";
    }

    my $seq_in = Bio::SeqIO->new(-file => $input_file, -format => 'Fasta');

    my $chunk_fh;
    my $seq_out;
    
    my $seq_count = 0;
    
    while (my $seq = $seq_in->next_seq()) {

        $seq_count++;

        if (($seq_count > $chunk_size) || (!defined($chunk_fh)))  {
        
            $seq_count = 1;

            ##FIXME: The temp dir location should not be hardcoded.  At least not here.
            $chunk_fh = File::Temp->new(
                                        'DIR'      => '/gscmnt/temp212/info/annotation/GAP_tmp',
                                        'SUFFIX'   => '.tmp',
                                        'TEMPLATE' => 'GAP_XXXXXXXX',
                                        'UNLINK'   => 0,
                                       );

            $seq_out = Bio::SeqIO->new(-fh => $chunk_fh, -format => 'Fasta');
        
            push @output_files, $chunk_fh->filename();

        }

        $seq_out->write_seq($seq);

    }

    $self->fasta_files(\@output_files);

    return 1;

}
 
1;
