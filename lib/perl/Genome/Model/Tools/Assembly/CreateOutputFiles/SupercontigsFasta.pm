package Genome::Model::Tools::Assembly::CreateOutputFiles::SupercontigsFasta;

use strict;
use warnings;

use Genome;
use IO::File;
use Data::Dumper;

use Bio::SeqIO;

class Genome::Model::Tools::Assembly::CreateOutputFiles::SupercontigsFasta {
    is => 'Genome::Model::Tools::Assembly::CreateOutputFiles',
    has => [
	directory => {
	    is => 'Text',
	    doc => 'Assembly directory',
	},
    ],
    has_optional => [
	contigs_bases_file => {
	    is => 'Text',
	    doc => 'Assembly contigs.bases file',
	    is_mutable => 1,
	},
	gap_file => {
	    is => 'Text',
	    doc => 'Assembly gap.txt file',
	    is_mutable => 1,
	},
	output_file => {
	    is => 'Text',
	    doc => 'Assembly supercontigs.agp file',
	    is_mutable => 1,
	},
    ],
};

sub help_brief {
    'Tool to create supercontigs.fasta file for assemblies';
}

sub help_detail {
    "Tool to create supercontigs.fasta file for assemblies";
}

sub execute {
    my $self = shift;

    #validate input values
    unless ($self->_validate_files()) {
	$self->error_message("Failed to validate input values");
	return;
    }

    my $gap_sizes;
    unless ($gap_sizes = $self->_load_gap_sizes()) {
	$self->error_message("Failed to get gap sizes");
	return;
    }
    #TODO - names it get_ids_from_fasta and move to base class
    my $contig_names;
    unless ($contig_names = $self->_get_contig_names()) {
	$self->errror_message("Failed to get contig names from contigs.bases file");
    }

    my $sctg_out = Bio::SeqIO->new(-format => 'fasta', -file => '>'.$self->output_file);
    my $seek_positions =  $self->seek_pos_from_contigs_file($self->contigs_bases_file, 'fasta');
    #returns contig number, 1.3 (from name like Contig1.3) and fh seek position .. hash of array

    my $fasta;
    my $fasta_length = 0;
    my $fasta_and_gap_length = 0;

    foreach my $contig_number ( sort {$a<=>$b } keys %$seek_positions ) {
	my ($sctg_num, $ctg_num) = $contig_number =~ /(\d+)\.(\d+)/;
	
	my $seq_obj = $self->_get_seq_obj(@{$seek_positions->{$contig_number}}[0]);	
	my $sctg_name = $self->_derive_supercontig_name($seq_obj->primary_id);
	my $next_ctg_in_scaf = $self->_next_scaffold_contig_name($seq_obj->primary_id);
	
	if (! exists $contig_names->{$next_ctg_in_scaf}) {
	    #single contig scaffold or last contig in scaffold
	    $fasta .= $seq_obj->seq;
	    $fasta_length += length $seq_obj->seq;
	    $fasta_and_gap_length += length $seq_obj->seq;

	    $seq_obj->id($sctg_name);
	    $seq_obj->seq($fasta);
	    $seq_obj->desc("$fasta_length $fasta_and_gap_length");
	    $sctg_out->write_seq($seq_obj);

	    $fasta = '';
	    $fasta_length = 0;
	    $fasta_and_gap_length = 0;
	}
	else { #append with next contig/gap
	    unless (exists $gap_sizes->{$seq_obj->primary_id}) {
		$self->warning_message("Didn't find gap size for contig ".$seq_obj->primary_id."\n\t".
				       "Expected one becase $next_ctg_in_scaf exists .. setting it to default: 100 bp");
	    }
	    my $gap_size = (exists $gap_sizes->{$seq_obj->primary_id}) ? $gap_sizes->{$seq_obj->primary_id} : 100;

	    $fasta .= $seq_obj->seq;
	    $fasta .= 'N' x $gap_size;#$gap_sizes->{$seq_obj->primary_id};
	    $fasta_length += length $seq_obj->seq;
	    $fasta_and_gap_length += length $seq_obj->seq;
	    $fasta_and_gap_length += $gap_size;#$gap_sizes->{$seq_obj->primary_id};
	}
    }
    
    return 1;
}

sub _derive_supercontig_name {
    my ($self, $contig_name) = @_;

    my ($sctg_name) = $contig_name =~ /(Contig\d+)\.\d+/;
    unless ($sctg_name) {
	$self->error_message("Failed to derive supercontig name from contig name: $contig_name"..
			     "\n\tExpected name like Contig2.5");
	return;
    }
    return $sctg_name;
}

sub _validate_files {
    my $self = shift;

    $self->error_message("Failed to find directory: ".$self->directory) unless -d $self->directory;

    $self->contigs_bases_file($self->directory.'/edit_dir/contigs.bases') unless $self->contigs_bases_file;
    $self->error_message("Failed to find needed input file: ".$self->contigs_bases_file) and return
	unless -s $self->contigs_bases_file;;

    $self->gap_file($self->gap_sizes_file) unless $self->gap_file;
    $self->error_message("Failed to find needed input file: ".$self->gap_file) and return
	unless -e $self->gap_file; #can be zero size

    $self->output_file($self->supercontigs_fasta_file) unless $self->output_file;
    
    return 1;
}

sub _load_gap_sizes {
    my $self = shift;

    my %gap_sizes;
    my $fh = Genome::Sys->open_file_for_reading($self->gap_file);
    while (my $line = $fh->getline) {
	next if $line =~ /^\s+$/;
	my ($name, $size) = $line =~ /(Contig\S+)\s+(\d+)/;
	unless ($name and $size) {
	    $self->error_message("Failed to get contig name and gap size from line: $line".
				 "\n\tExpected lines like: Contig2.5 100");
	    return;
	}
	$gap_sizes{$name} = $size;
    }
    $fh->close;
    return \%gap_sizes;
}

sub _get_contig_names {
    my $self = shift;
    my %names;
    my $io = Bio::SeqIO->new(-format => 'fasta', -file => $self->contigs_bases_file);
    while (my $seq = $io->next_seq) {
	$names{$seq->primary_id} = 1;
    }
    return \%names;
}

sub _get_seq_obj {
    my ($self, $seek_pos) = @_;
    my $fh = Genome::Sys->open_file_for_reading($self->contigs_bases_file);
    $fh->seek($seek_pos, 0);
    my $io = Bio::SeqIO->new(-format => 'fasta', -fh => $fh);
    my $seq = $io->next_seq;
    unless ($seq) {
	$self->error_message("Failed to get seq object for seek position $seek_pos");
	return;
    }
    $fh->close;
    return $seq;
}

sub _next_scaffold_contig_name {
    my ($self, $contig_name) = @_;
    my ($sctg, $ctg) = $contig_name =~ /^Contig(\d+)\.(\d+)$/;
    unless (defined $sctg and defined $ctg) {
	$self->error_message("Can't determine supercontig and contig number from contig name: $contig_name".
			     "\n\tExpected name like Contig1.2");
	return;
    }
    return 'Contig'.$sctg.'.'.++$ctg;
}

1;
