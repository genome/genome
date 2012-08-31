package Genome::Model::Tools::Assembly::CreateOutputFiles::ContigsFromAce;

use strict;
use warnings;

use Genome;
use IO::File;
use Data::Dumper;

use Bio::SeqIO;

class Genome::Model::Tools::Assembly::CreateOutputFiles::ContigsFromAce {
    is => 'Genome::Model::Tools::Assembly::CreateOutputFiles',
    has => [
	acefile => {
	    is => 'Text',
	    doc => 'Ace file to get fasta and qual from',
	},
	directory => {
	    is => 'Text',
	    doc => 'Assembly build directory, not edit_dir',
	},
    ],
    has_optional => [
	fasta_out => {
	    is => 'Text',
	    doc => 'Output fasta file name',
	    is_mutable => 1,
	},
	qual_out => {
	    is => 'Text',
	    doc => 'Output qual file name',
	    is_mutable => 1,
	},
    ],
    has_optional_transient => [
	_int_fasta_out => { is => 'Text', doc => 'Intermediate fasta out file' },
	_int_qual_out => { is => 'Text', doc => 'Intermediate qual out file' },
    ],
};

sub help_brief {
    'Tool to create contigs.bases and contigs.qual files from ace file';
}

sub help_detail {
    "Tool to to create a file containing fasta of contigs in ace files";
}

sub execute {
    my $self = shift;

    unless (-s $self->acefile) {
	$self->error_message("Failed to find ace file: ".$self->acefile);
	return;
    }

    #make intermediated unsorted bases and qual files
    unless ( $self->_get_fasta_qual_from_ace ) {
    	$self->error_message("Failed to get extrace fasta and qual from ace");
    	return;
    }

    #sort and print contigs bases and qual files
    unless ( $self->_sort_contigs_files ) {
    	$self->error_message("Failed to sort contigs fasta and qual files");
    	return;
    }

    #remove intermediate files
    unlink $self->_int_fasta_out, $self->_int_qual_out;

    return 1;
}

sub _sort_contigs_files {
    my $self = shift;

    #get fasta and qual file seek positioon
    my $fpos = $self->seek_pos_from_contigs_file($self->_int_fasta_out, 'fasta');
    my $qpos = $self->seek_pos_from_contigs_file($self->_int_qual_out, 'qual');

    unless ($self->fasta_out) {
	$self->fasta_out( $self->contigs_bases_file );
    }
    my $f_io = Bio::SeqIO->new(-format => 'fasta', -file => '>'.$self->fasta_out);

    unless ($self->qual_out) {
	$self->qual_out( $self->contigs_quals_file );
    }
    my $q_io = Bio::SeqIO->new(-format => 'qual', -file => '>'.$self->qual_out);

    #iterate through fasta seek pos and get fasta and qual
    foreach my $contig_number (sort {$a<=>$b} keys %$fpos) {
	my $f_seek_pos = @{$fpos->{$contig_number}}[0];
	my $q_seek_pos = @{$qpos->{$contig_number}}[0];
	unless (defined $f_seek_pos) {
	    $self->error_message("Failed to get fasta seek position");
	    return;
	}
	unless (defined $q_seek_pos) {
	    $self->error_message("Failed to get qual seek position");
	    return;
	}

	my $fbo = $self->_get_bio_obj($self->_int_fasta_out, 'fasta', $f_seek_pos);
	my $qbo = $self->_get_bio_obj($self->_int_qual_out, 'qual', $q_seek_pos);

	unless ($fbo->primary_id eq $qbo->primary_id) {
	    $self->error_message("Fasta and quality objects are not from the same contigs:\n\t".
				 "got from fasta: ".$fbo->primary_id." from quality ".$qbo->primary_id);
	    return;
	}

	$f_io->write_seq($fbo);
	$q_io->write_seq($qbo);
    }

    return 1;
}

sub _get_bio_obj {
    my ($self, $file, $format, $seek_pos) = @_;
    #this doesn't seek to work if fh is held open, ie, passed in
    my $fh = Genome::Sys->open_file_for_reading($file);
    $fh->seek($seek_pos, 0);
    my $io = Bio::SeqIO->new(-fh => $fh, -format => $format);
    my $seq = $io->next_seq;
    $fh->close;
    return $seq;
}

sub _get_fasta_qual_from_ace {
    my $self = shift;

    #existing ace parser are not used because it can't load really big ace files
    #need to clean this up a bit ..

    my $ace_fh = Genome::Sys->open_file_for_reading($self->acefile);

    #handle for intermediate fasta
    $self->_int_fasta_out($self->directory."/edit_dir/int.contigs.bases");
    unlink $self->_int_fasta_out;
    my $fasta_fh = Genome::Sys->open_file_for_writing($self->_int_fasta_out);

    #handle for intermediate qual
    $self->_int_qual_out($self->directory."/edit_dir/int.contigs.quals");
    unlink $self->_int_qual_out;
    my $qual_fh = Genome::Sys->open_file_for_writing($self->_int_qual_out);

    my $is_fasta = 0;
    my $is_qual = 0;
    my $contig_name;
    my $base_line_count = 0;
    my $comp; #Complemented or uncomplemented
    #this is a bit ugly .. just parsing through ace file line by line
    #and getting consensus fastas and quals
    
    while (my $line = $ace_fh->getline) {
	next if $line =~ /^\s+$/;
	chomp $line;
	if ($line =~ /^CO\s+/) {
	    ($contig_name, $comp) = $line =~ /^CO\s+(\S+)\s+\d+\s+\d+\s+\d+\s+(\w)/;
	    $fasta_fh->print (">$contig_name\n");
	    $is_fasta = 1;
	    next;
	}
	if ($line =~ /^BQ/) {
	    $is_fasta = 0;
	    $fasta_fh->print ("\n") if $base_line_count != 0;
	    $base_line_count = 0;
	    $is_qual = 1;
	    $qual_fh->print (">$contig_name\n");
	    next;
	}
	if ($is_fasta == 1) {
	    my @bases = split (//, $line);
	    foreach my $base (@bases) {
		if ($base =~ /^[acgtxn]$/i) {
		    $base_line_count++;
		    $fasta_fh->print (uc $base);
		    if ($base_line_count == 60) {
			$fasta_fh->print ("\n");
			$base_line_count = 0;
		    }
		}
		else {
		    next;
		}
	    }
	}
	if ($line =~ /^AF\s+/) {
	    $qual_fh->print ("\n") if $base_line_count != 0;
	    $base_line_count = 0;
	    $is_qual = 0;
	    next;
	}
	if ($is_qual == 1) {
	    my @quals = split (/\s+/, $line);
	    foreach my $qual (@quals) {
		next unless $qual =~ /^\d+$/;
		$base_line_count++;
		if ($base_line_count == 60) {
		    $qual_fh->print ("$qual\n");
		    $base_line_count = 0;
		}
		else {
		    $qual_fh->print ("$qual ");
		}
	    }
	}
    }

    $ace_fh->close;
    $fasta_fh->close;
    $qual_fh->close;
    
    return 1;
}

1;
