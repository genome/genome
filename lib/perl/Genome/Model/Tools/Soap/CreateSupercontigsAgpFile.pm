package Genome::Model::Tools::Soap::CreateSupercontigsAgpFile;

use strict;
use warnings;

use Genome;
use Bio::SeqIO;
use Data::Dumper;

class Genome::Model::Tools::Soap::CreateSupercontigsAgpFile {
    is => 'Genome::Model::Tools::Soap::Base',
    has => [
        assembly_directory => {
            is => 'Text',
            doc => 'Soap assembly directory',
        },
        min_contig_length => {
            is => 'Number',
            doc => 'Minimum contig length to process',
        },	
    ],
};

sub help_brief {
    'Tool to create supercontigs.agp file from soap created scaffold fasta file';
}

sub help_detail {
    return <<"EOS"
gmt soap create-supercontigs-agp-file --assembly-directory /gscmnt/111/soap_assembly --min-contig-length 200
EOS
}

sub execute {
    my $self = shift;

    unless ( $self->create_edit_dir ) {
	$self->error_message("Failed to create edit_dir");
	return;
    }

    unless (-d $self->assembly_directory) {
        $self->error_message("Failed to find assembly directory: ".$self->assembly_directory);
        return;
    }
    
    unless( -s $self->assembly_scaffold_sequence_file ) {
        $self->error_message("Failed to find soap scaffold sequence file: ".$self->assembly_scaffold_sequence_file);
        return;
    }

    unlink $self->supercontigs_agp_file;
    my $fh = Genome::Sys->open_file_for_writing( $self->supercontigs_agp_file );

    my $io = Bio::SeqIO->new( -format => 'fasta', -file => $self->assembly_scaffold_sequence_file );

    my $scaffold_number = 0;

    while (my $seq = $io->next_seq) {
        #remove lead/trail-ing Ns
        my $supercontig = $seq->seq;
        $supercontig =~ s/^N+//;
        $supercontig =~ s/N+$//;

        #skip if less than min contig length .. need to move this if all contigs in scaffold is < min contig length .. will skip an iteration
        next unless length $supercontig >= $self->min_contig_length;

	my $scaffold_name = 'Contig'.$scaffold_number++;
	my @bases = split (/N+/i, $supercontig);
	my @gaps = split (/[ACTG]+/i, $supercontig);

	shift @gaps; #empty string from split .. unless seq starts with Ns in which case it's just thrown out

	my $start_pos = 0;
	my $stop_pos = 0;
	my $fragment_order = 0;
	my $contig_order = 0;
        my $gap_length = 0;

        my $is_leading_contig = 1;

        for my $i ( 0 .. $#bases ) {
            #add contig length to gap when contig < min length
            $gap_length += length $gaps[$i - 1] unless $is_leading_contig;
            if ( length $bases[$i] < $self->min_contig_length ) {
                $gap_length += length $bases[$i];
                next;
            }
            #fragment info
            {
                #no frag info needed for first contig in scaffold
                next if $is_leading_contig;
                $start_pos = $stop_pos + 1;
                $stop_pos = $start_pos + $gap_length - 1;
                $fh->print( $scaffold_name."\t".$start_pos."\t".$stop_pos."\t".++$fragment_order."\tN\t".$gap_length."\tscaffold\tyes\tpaired-ends\n" );
                #reset gap length after printing a gap info
                $gap_length = 0;
            }
            #contig info
            {
                $start_pos = $stop_pos + 1;
                $stop_pos = $start_pos + ( length $bases[$i] ) - 1;
                my $contig_name = $scaffold_name.'.'.++$contig_order;
                my $contig_length = length $bases[$i];
                $fh->print( $scaffold_name."\t".$start_pos."\t".$stop_pos."\t".++$fragment_order."\tW\t".$contig_name."\t1\t".$contig_length."\t+\n" );
                #scaffold started .. no longer leading contig from this point
                $is_leading_contig = 0;
            }
	}
    }

    $fh->close;

    return 1;
}

1;
