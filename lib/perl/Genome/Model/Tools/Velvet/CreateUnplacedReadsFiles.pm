package Genome::Model::Tools::Velvet::CreateUnplacedReadsFiles;

use strict;
use warnings;

use Genome;
use Bio::Seq;
use Bio::SeqIO;
use AMOS::AmosLib;

class Genome::Model::Tools::Velvet::CreateUnplacedReadsFiles {
    is => 'Genome::Model::Tools::Velvet::Base',
    has => [
	assembly_directory => {
	    is => 'Text',
	    doc => 'Assembly directory',
	},
        min_contig_length => {
            is => 'Number',
            doc => 'Minimum contig length to process',
        },
    ],
};

sub help_brief {
    'Tool to create pcap style reads.unplaced and reads.unplaced.fasta files'
}

sub help_detail {
    return <<EOS
gmt velvet create-unplaced-reads-files --assembly-directory /gscmnt/111/e_coli_velvet_asm --min-contig-length 200
EOS
}

sub execute {
    my $self = shift;

    #make edit_dir
    unless ( $self->create_edit_dir ) {
	$self->error_message("Assembly edit_dir does not exist and could not create one");
	return;
    }

    my $seek_positions = $self->load_sequence_seek_positions( $self->velvet_sequences_file );
    unless ($seek_positions) { #arrayref
	$self->error_message("Failed to get read names and seek_pos from Sequences file");
	return;
    }

    my $unplaced_reads = $self->_remove_placed_reads( $seek_positions );
    unless ($unplaced_reads) {
	#TODO - make sure this works for empty array ref too
	$self->error_message("Failed to remove placed reads from input reads");
	return;
    }

    unless ($self->_print_unplaced_reads($unplaced_reads, $self->velvet_sequences_file)) {
	$self->error_message("Failed to print unplaced reads");
	return;
    }

    return 1;
}

sub _print_unplaced_reads {
    my ($self, $unplaced_reads, $sequences_file) = @_;

    unlink $self->reads_unplaced_file;
    my $unplaced_fh = Genome::Sys->open_file_for_writing($self->reads_unplaced_file);
    my $fasta_out = Bio::SeqIO->new(-format => 'fasta', -file => '>'.$self->reads_unplaced_fasta_file);
    my $seq_fh = Genome::Sys->open_file_for_reading( $sequences_file );
    my $bio_seqio = Bio::SeqIO->new(-fh => $seq_fh, -format => 'fasta', -noclose => 1);
    for (0 .. $#$unplaced_reads) {
        next unless defined @$unplaced_reads[$_];
        my $seek_pos = ${$unplaced_reads}[$_];
        unless (defined $seek_pos) {
            $self->error_message("Failed to get read seek position for afg read index $_");
            return;
        }
        $seq_fh->seek($seek_pos, 0);
        my $seq = $bio_seqio->next_seq;
        $unplaced_fh->print( $seq->primary_id." unused\n" );
        my $seq_obj = Bio::Seq->new(-display_id => $seq->primary_id, -seq => $seq->seq);
        $fasta_out->write_seq($seq_obj);
    }
    $seq_fh->close;
    $unplaced_fh->close;

    return 1;
}

sub _remove_placed_reads {
    my ($self, $seek_pos) = @_;

    my $scaffold_info;
    unless( $scaffold_info = $self->get_scaffold_info_from_afg_file ) {
        $self->error_message( "Failed to get scaffolding info from afg file" );
        return;
    }

    my $afg_fh = Genome::Sys->open_file_for_reading($self->velvet_afg_file);

    while (my $record = getRecord($afg_fh)) {
	my ($rec, $fields, $recs) = parseRecord($record);
	if ($rec eq 'CTG') {

            my $contig_seq = $fields->{seq};
            $contig_seq =~ s/\n//g;

            my $contig_name = $fields->{eid};
            $contig_name =~ s/\-/\./;

            next if $scaffold_info->{$contig_name}->{is_removed};
            
	    for my $r (0 .. $#$recs) {
		my ($srec, $sfields, $srecs) = parseRecord($recs->[$r]);
		if ($srec eq 'TLE') {
		    #sfields:
		    #'src' => '19534',  #read id number
		    #'clr' => '0,90',   #read start, stop 0,90 = uncomp 90,0 = comp
		    #'off' => '75'      #read off set .. contig start position
		    @$seek_pos[$sfields->{src}] = undef;
		}
	    }
	}
    }

    $afg_fh->close;
    return $seek_pos;
}

1;
