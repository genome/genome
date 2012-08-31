package Genome::Model::Tools::Assembly::CreateOutputFiles::Gap;

use strict;
use warnings;

use Genome;
use Bio::SeqIO;
use IO::File;

class Genome::Model::Tools::Assembly::CreateOutputFiles::Gap {
    is => 'Genome::Model::Tools::Assembly::CreateOutputFiles',
    has => [
        directory => {
            is => 'Text',
            doc => 'Assembly build directory',
        },
    ],
};

sub help_brief {
    'Tool to create gap.txt file from velvet created contigs.fa file';
}

sub help_detail {
    "Tool to create a file containing estimated gap sizes between neighboring contigs";
}

sub execute {
    my $self = shift;

    unless (-s $self->directory.'/contigs.fa') {
	$self->error_message("Failed to find velvet generated contigs.fa file");
	return;
    }

    unlink $self->gap_sizes_file;
    my $gap_fh = Genome::Sys->open_file_for_writing( $self->gap_sizes_file );

    my $io = Bio::SeqIO->new(-format => 'fasta', -file => $self->directory.'/contigs.fa');

    while (my $seq = $io->next_seq) {
	my @bases = split (/N+/i, $seq->seq);
	my @n_s = split (/[ACGT]+/i, $seq->seq);
        #SINGLE CONTIG SCAFFOLD .. OR ALL NS .. NO GAP INFO NEEDED
	next if scalar @bases le 1;
	#RESOLVE BLANK ELEMENT IN SHIFTED ARRAY
	if ($seq->seq =~ /^N/i) {
	    shift @bases; #BLANK
	    shift @n_s;   #LEADING NS .. NOT NEEDED
	}
	else {
	    shift @n_s;   #BLANK
	}
	#GET RID OF TAILING NS .. NOT NEEDED
	pop @n_s if $seq->seq =~ /N$/i;
	#DOUBLE CHECK BASES AND GAP COUNTS
	unless (scalar @n_s == scalar @bases - 1) {
	    $self->error_message("Error: Unable to match sequences and gaps properly \n");
	    return;
	}
	#DECIPER SUPERCONTIG/CONTIG NAMES
	my ($node_num) = $seq->primary_id =~ /NODE_(\d+)_/;
	unless ($node_num) {
	    $self->error_message("Can not determine node number\n");
	    return;
	}
	my $supercontig_number = $node_num - 1;
	my $contig_number = 1;
	#ITERATE THROUGH BASES ARRAY AND ASSIGN GAP SIZE
	while (my $base_string = shift @bases) {
	    next if scalar @bases == 1; #LAST CONTIG OF SCAF .. NO GAP INFO
	    my $gap_size = length (shift @n_s);
	    #PRINT CONTIG NAME AND GAP SIZE: eg Contig855.26 91
	    $gap_fh->print ('Contig'.$supercontig_number.'.'.$contig_number.' '.$gap_size."\n");
	    $contig_number++;
	}
    }
    $gap_fh->close;
    return 1;
}

1;
