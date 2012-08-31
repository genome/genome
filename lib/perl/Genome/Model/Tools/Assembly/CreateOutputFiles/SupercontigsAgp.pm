package Genome::Model::Tools::Assembly::CreateOutputFiles::SupercontigsAgp;

use strict;
use warnings;

use Genome;
use Bio::SeqIO;
use IO::File;
use Data::Dumper;

class Genome::Model::Tools::Assembly::CreateOutputFiles::SupercontigsAgp {
    is => 'Genome::Model::Tools::Assembly::CreateOutputFiles',
    has => [
	directory => {
	    is => 'Text',
	    doc => 'Assembly main build directory',
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
    'Tool to create supercontig.agp file for assemblies';
}

sub help_detail {
    "Tool to create supercontigs.agp file for assemblies";
}

sub execute {
    my $self = shift;

    unless ($self->_validate_files()) {
	$self->error_message("Failed to validate files needed to run this tool");
	return;
    }

    my $contig_names;
    unless ($contig_names = $self->_get_contig_names()) {
	$self->error_message("Failed to get contig names");
	return;
    }

    my $gap_sizes;
    unless ($gap_sizes = $self->_load_gap_sizes()) {
	$self->error_message("Failed to get gap sizes");
	return;
    }

    unlink $self->output_file;
    my $out_fh = Genome::Sys->open_file_for_writing($self->output_file);

    my $seek_positions = $self->seek_pos_from_contigs_file($self->contigs_bases_file, 'fasta');
    #returns contig number, 1.3 (from name like Contig1.3) and fh seek position .. hash of array

    #array to keep track of number of contigs in scaffold
    my @scaffold_data; #[0] = number of contig and gaps - ++ for each gap and contig
                       #[1] = next contig start position - incremented with prev contig length and gap length

    foreach my $contig_number (sort {$a<=>$b} keys %$seek_positions) {
	my ($sctg_num, $ctg_num) = $contig_number =~ /(\d+)\.(\d+)/;
	my $seq = $self->_get_seq_obj(@{$seek_positions->{$contig_number}}[0]);	

	my ($sctg_name) = $seq->primary_id =~ /(\S+)\.\d+/;
	my $start = (defined $scaffold_data[$sctg_num][1]) ? $scaffold_data[$sctg_num][1] : 0;
	$start++;
	my $stop = $start + length $seq->seq;
	$stop--;

	$scaffold_data[$sctg_num][0]++ unless defined $scaffold_data[$sctg_num][0];
	my $order = $scaffold_data[$sctg_num][0];
	
	#for contigs print:
	#sctg    start   stop    order   W       contig name     1       length  +
        #Contig1 1       380     1       W       Contig1.1       1       380     +
       
	$out_fh->print ("$sctg_name\t$start\t$stop\t$order\tW\t". $seq->primary_id ."\t1\t". length ($seq->seq) ."\t+\n");
	#Note - if contig is complemented, it would be notied with - at the end (instead of +)
	#but since we're getting sequences from fasta file, assume all are uncomplemented
	
	#Increment for next contig/gap position/order
	$scaffold_data[$sctg_num][0]++;
	$scaffold_data[$sctg_num][1] += length $seq->seq;

	#check .. if not last contig in scaffold, gap must exist
	my $next_contig = $sctg_num.'.'.++$ctg_num;
	if (exists $seek_positions->{$next_contig} and  ! exists $gap_sizes->{$seq->primary_id}) {
	    $self->warning_message("Next contig in scaffold Contig$next_contig exists but ". $seq->primary_id ." gap size does not".
				   "\n\tSetting gap size to default value of 100 bp");
	}

	#check .. if next contig does not exist gap size must not exist
	if (! exists $seek_positions->{$next_contig} and  exists $gap_sizes->{$seq->primary_id}) {
	    $self->warning_message($seq->primary_id ." is last contig in scaffold so gap size should not exist but does");
	}

	next unless exists $gap_sizes->{$seq->primary_id};

	#for gaps print:
	#sctg    start   stop    order   N       length  fragment        yes
	#Contig1 381     453     2       N       73      fragment        yes

	my $gap_size = (exists $gap_sizes->{$seq->primary_id}) ? $gap_sizes->{$seq->primary_id} : 100;
        my $status = ( exists $gap_sizes->{$seq->primary_id} ) ? 'N' : 'U';
	$start = $scaffold_data[$sctg_num][1];
	$start++;
	$stop = $start + $gap_size; #$gap_sizes->{$seq->primary_id};

	$stop--;
	$order = $scaffold_data[$sctg_num][0];
	$out_fh->print ("$sctg_name\t$start\t$stop\t$order\t$status\t". $gap_sizes->{$seq->primary_id}. "\tscaffold\tyes\tpaired-ends\n");

	#Increment fo next contig/gap position/order
	$scaffold_data[$sctg_num][0]++;
	$scaffold_data[$sctg_num][1] += $gap_size; #$gap_sizes->{$seq->primary_id};
    }

    $out_fh->close;

    return 1;
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

sub _validate_files {
    my $self = shift;

    $self->error_message("Failed to find or invalid directory: ".$self->directory) and return
	unless -d $self->directory;

    #contigs.bases file
    $self->contigs_bases_file($self->directory.'/edit_dir/contigs.bases') unless
	$self->contigs_bases_file;
    $self->error_message("Failed to find file or file is zero size: ".$self->contigs_bases_file) and return
	unless -s $self->contigs_bases_file;

    #gap.txt file
    $self->gap_file($self->gap_sizes_file) unless
	$self->gap_file;
    $self->error_message("Failed to find needed input file: ".$self->gap_file) and return
	unless -e $self->gap_file; #file can be zero size

    $self->output_file( $self->supercontigs_agp_file ) unless $self->output_file;

    return 1;
}


1;
