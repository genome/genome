package Genome::Model::Tools::Soap::CreateSupercontigsFastaFile;

use strict;
use warnings;

use Genome;
use Bio::SeqIO;

class Genome::Model::Tools::Soap::CreateSupercontigsFastaFile {
    is => 'Genome::Model::Tools::Soap::Base',
    has => [
        assembly_directory => {
            is => 'Text',
            doc => 'Soap assembly directory',
        },
        min_contig_length => {
            is => 'Number',
            doc => 'Minimum contig length to filter',
        },
    ],
};

sub help_brief {
    'Tool to create supercontigs.fasta file from soap created scaffold fasta file';
}

sub help_detail {
    return <<"EOS"
gmt soap create-supercontigs-fasta-file --scaffold-sequence-file /gscmnt/111/soap_assembly/61EFS.cafSeq --assembly-directory /gscmnt/111/soap_assembly
EOS
}

sub execute {
    my $self = shift;

    unless ( $self->create_edit_dir ) {
	$self->error_message("Failed to creat edit_dir");
	return;
    }

    unless (-d $self->assembly_directory) {
        $self->error_message("Failed to find assembly directory: ".$self->assembly_directory);
        return;
    }

    unless( -s $self->assembly_scaffold_sequence_file ) {
        $self->error_message("Failed to find soap output file: ".$self->assembly_scaffold_sequence_file );
        return;
    }

    my $in = Bio::SeqIO->new(-format => 'fasta', -file => $self->assembly_scaffold_sequence_file);
    my $out = Bio::SeqIO->new(-format => 'fasta', -file => '>'.$self->supercontigs_fasta_file);

    my $supercontig_number = 0;
    while (my $seq = $in->next_seq) {
        my $old_fasta = $seq->seq;
        #shouldn't be heading or tailing Ns but could mess things up
        $old_fasta =~ s/^N+//;
        $old_fasta =~ s/N+$//;
        next unless length $old_fasta >= $self->min_contig_length;
        my @bases = split( /N+/i, $old_fasta );
        my @gaps = split ( /[ACGT]+/i, $old_fasta );
        shift @gaps; #split creates undef first array element
        my $pos = 0;
        my $new_fasta;
        for my $string ( @bases ) {
            if ( length $string >= $self->min_contig_length ) {
                $new_fasta .= $string;
                $new_fasta .= $gaps[$pos] if $gaps[$pos];#gap won't exist for last the string
            }
            else {
                #make strings < min_contig length gaps
                $new_fasta .= 'N' x length $string;
                $new_fasta .= $gaps[$pos] if $gaps[$pos];
            }
            $pos++;
        }
        unless( length $old_fasta eq length $new_fasta ) {
            $self->error_message( "Length of new supercontig fasta did not match the length of old: new: " . length ($new_fasta).", old: ".length $old_fasta );
            return;
        }
        $seq->id('Contig'.$supercontig_number++);
        $seq->seq( $new_fasta );

	$out->write_seq($seq);
    }

    return 1;
}

1;
