package Genome::Model::Tools::Assembly::CreateOutputFiles;

use strict;
use warnings;

use Genome;
use AMOS::AmosLib;

use Bio::SeqIO;

class Genome::Model::Tools::Assembly::CreateOutputFiles {
    is => 'Command',
    has => [
    ],
};

sub help_brief {
    'Tools for creating assembly output files'
}

sub help_detail {
    "Tools for creating assembly output files";
}

sub seek_pos_from_contigs_file {
    my ($self, $file, $type) = @_;

    unless (-s $file) {
	$self->error_message("Can't find file: $file");
	return;
    }
    unless ($type eq 'fasta' or $type eq 'qual') {
	$self->error_message("File type must be fasta or qual and not $type");
	return;
    }
    my %contig_infos;
    my $fh = Genome::Sys->open_file_for_reading($file) ||
        return;
    my $seek_pos = $fh->tell;
    my $io = Bio::SeqIO->new(-format => $type, -fh => $fh);
    while (my $seq = $io->next_seq) {
	my ($contig_number) = $seq->primary_id =~ /contig(\S+)/i;
	unless (defined $contig_number) {
	    $self->error_message("Expected header or id to like Contig1.1 or contig9 but got ".$seq->primary_id);
	    return;
	}
	#purpose of indexing by contig number is to be able
	#to sort by contig number and allow quicker access to seq or qual

	#temp fix for something that happened that's causing seek pos to be off by one
	#not sure what the problem is yet .. files did not change so likely due to Bio SeqIO changes
	$seek_pos = ( $seek_pos == 0 ) ? $seek_pos : $seek_pos - 1;

	push @{$contig_infos{$contig_number}}, $seek_pos;
        $seek_pos = $fh->tell;
    }
    $fh->close;
    return \%contig_infos;
}

sub get_gap_sizes {
    my $self = shift;
    
    my %gap_sizes;

    unless (-e $self->gap_sizes_file) {
	#file should exist 0 size even if assembly has no scaffolds
	#return blank hash if no gap sizes
	$self->error_message("Can't find gap.txt file: ".$self->gap_sizes_file);
	return;
    }

    my $fh = IO::File->new("<".$self->gap_sizes_file) ||
	die "Can not create file handle to read gap.txt file\n";

    while (my $line = $fh->getline) {
	chomp $line;
	my ($contig_name, $gap_size) = split (/\s+/, $line);
	unless ($contig_name =~ /Contig\d+\.\d+/ and $gap_size =~ /\d+/) {
	    $self->error_message("Gap.txt file lines should look like this: Contig4.1 125".
				 "\n\tbut it looks like this: ".$line);
	    return;
	}
	$gap_sizes{$contig_name} = $gap_size;
    }
    $fh->close;

    return \%gap_sizes;
}

sub get_contig_lengths {
    my ($self, $afg_file) = @_;
    my %contig_lengths;
    my $fh = Genome::Sys->open_file_for_reading($afg_file)
	or return;
    while (my $record = getRecord($fh)) {
	my ($rec, $fields, $recs) = parseRecord($record);
	if ($rec eq 'CTG') {
	    my $seq = $fields->{seq};
	    $seq =~ s/\n//g;

	    my ($sctg_num, $ctg_num) = split('-', $fields->{eid});
	    my $contig_name = 'Contig'.--$sctg_num.'.'.++$ctg_num;

	    $contig_lengths{$contig_name} = length $seq;
	}
    }
    $fh->close;
    return \%contig_lengths;
}

sub input_fasta_files {
    my $self = shift;
    my @input_fastas;
    foreach ( glob($self->directory."/edit_dir/*") ) {
	if ($_ =~ /(fasta|clip)\.gz$/) {
	    #check to make sure qual fine exists for fasta
	    #so it doesn't just grab any fastas
	    my $qual = $_;
	    $qual =~ s/\.gz$/\.qual\.gz/;
	    push @input_fastas, $_ if -s $qual;
	}
    }
    unless (@input_fastas) {
	$self->error_message("Failed to find any valid input fastas in dir: ".$self->directory.'/edit_dir');
	return;
    }
    return @input_fastas;
}


sub contigs_bases_file {
    return $_[0]->directory.'/edit_dir/contigs.bases';
}

sub contigs_quals_file {
    return $_[0]->directory.'/edit_dir/contigs.quals';
}

sub gap_sizes_file {
    return $_[0]->directory.'/edit_dir/gap.txt';
}

sub read_info_file {
    return $_[0]->directory.'/edit_dir/readinfo.txt';
}

sub reads_placed_file {
    return $_[0]->directory.'/edit_dir/reads.placed';
}

sub reads_unplaced_file {
    return $_[0]->directory.'/edit_dir/reads.unplaced';
}

sub reads_unplaced_fasta_file {
    return $_[0]->directory.'/edit_dir/reads.unplaced.fasta';
}

sub stats_file {
    return $_[0]->directory.'/edit_dir/stats.txt';
}

sub supercontigs_agp_file {
    return $_[0]->directory.'/edit_dir/supercontigs.agp';
}

sub supercontigs_fasta_file {
    return $_[0]->directory.'/edit_dir/supercontigs.fasta';
}

1;
