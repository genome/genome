package Genome::Model::Tools::Assembly::CreateOutputFiles::ReadsUnplaced;

use strict;
use warnings;

use Genome;
use IO::File;
use Bio::SeqIO;
use Data::Dumper;

class Genome::Model::Tools::Assembly::CreateOutputFiles::ReadsUnplaced {
    is => 'Genome::Model::Tools::Assembly::CreateOutputFiles',
    has => [
	directory => {
	    is => 'Text',
	    doc => 'Main assembly directory',
	},
	reads_placed_file => {
	    is => 'Text',
	    doc => 'Assembly reads.placed file',
	    is_optional => 1,
	},
	reads_unplaced_out => {
	    is => 'Text',
	    doc => 'Reads unplaced output file',
	    is_optional => 1,
	},
	reads_unplaced_fasta_out => {
	    is => 'Text',
	    doc => 'Reads unplaced fasta out',
	    is_optional => 1,
	},
    ],
};

sub help_brief {
    'Tool to create assembly reads unplaced files'
}

sub help_synopsis {
}

sub help_detail {
}

sub execute {
    my $self = shift;

    my $placed_reads = $self->_get_placed_reads; #hashref of placed reads

    my $unplaced_out = ($self->reads_unplaced_out) ? $self->reads_unplaced : $self->reads_unplaced_file;
    my $unplaced_fasta_out = ($self->reads_unplaced_fasta_out) ? $self->reads_unplaced_fasta_out : $self->reads_unplaced_fasta_file;

    my $out = Bio::SeqIO->new(-format => 'fasta', -file => '>'.$unplaced_fasta_out);
    unlink $unplaced_out;
    my $fh_out = Genome::Sys->open_file_for_writing($unplaced_out);
    foreach ( $self->input_fasta_files ) {
	my $fh = IO::File->new("zcat $_ |") || die;
	my $io = Bio::SeqIO->new(-format => 'fasta', -fh => $fh);
	while (my $seq = $io->next_seq) {
	    next if exists $placed_reads->{$seq->id};
	    $out->write_seq($seq);
	    $fh_out->print($seq->id." unused\n");
	}
	$fh->close;
    }
    $fh_out->close;

    return 1;
}

sub _get_placed_reads {
    my $self = shift;
    my $reads_placed_file = ($self->reads_placed_file) ? $self->reads_placed_file : $self->directory.'/edit_dir/reads.placed';
    unless (-s $reads_placed_file) {
	$self->error_message("Failed to find reads.placed file in assembly edit_dir");
	return;
    }
    my $fh = Genome::Sys->open_file_for_reading($reads_placed_file);
    my %placed_reads;
    while (my $line = $fh->getline) {
	my ($read_name) = $line =~ /^\*\s+(\S+)\s+/;
	#newbler 454
	$read_name =~ s/(_left|_right)$//;
	#velvet solexa
	$read_name =~ s/\-\d+$//;
	$placed_reads{$read_name} = 1;
    }
    $fh->close;
    return \%placed_reads;
}

1;
