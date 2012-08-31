package Genome::Model::Tools::ViromeEvent::Assignment::Summary;

use strict;
use warnings;

use Genome;
use Workflow;
use Bio::SeqIO;
use IO::File;
use File::Basename;
use Data::Dumper;

class Genome::Model::Tools::ViromeEvent::Assignment::Summary{
    is => 'Genome::Model::Tools::ViromeEvent',
};

sub help_brief {
    "Tool to generate summary report for virome screening";
}

sub help_detail {
    return <<"EOS"
This script will read the assignment report files in the given 
directory and generate a summary report for a given library. It will report 
in each library, for each category, how many total sequence were 
assigned to this category, how many were assigned by BLASTN, how many
were assigned by TBLASTX.

It will also filter the virus lineage, leave out virus that are phage.
It will rank the virus lineage by range of percent ID from low to high. 

It will also generate a .InterestingReads report about the details of each lineage.
EOS
}

sub execute {
    my $self = shift;
    my $lib_name = basename($self->dir);
    $self->log_event("Assignment Summary starting on $lib_name");
    #ASSIGNMENT SUMMARY FILE
    my $assignment_sum_file = $self->dir.'/'.$lib_name.'.AssignmentSummary';
    my $asum_fh = IO::File->new("> $assignment_sum_file") || die
	"Can not create file handle for $assignment_sum_file";
    #GET READ COUNTS FROM BLAST IN/OUT FA FILES
    unless ($self->_print_seq_count_summary($asum_fh)) {
	$self->log_event("Failed to print seq count summary");
	return;
    }
    #GATHER LINEAGE AND RESULT METRICS FOR VIRAL HIT SEQS
    my $read_infos;
    unless ($read_infos = $self->_parse_assignment_report_file($asum_fh)) {
	$self->log_event("Failed to parse assignment report file");
	return;
    }
    $asum_fh->close;
    #GET DETAILED INFO ABOUT INTERESTING READS (VIRAL HIT READS)
    my $interesting_reads_file = $self->dir.'/'.$lib_name.'.InterestingReads';
    my $orig_fasta_file = $self->dir.'/'.$lib_name.'.fa';
    my $ir_fh = IO::File->new("> $interesting_reads_file") ||
	die "Can not create file handle for $interesting_reads_file";
    unless ($self->_parse_interesting_reads_file($ir_fh, $read_infos)) {
	$self->log_event("Failed to parse interesting reads files");
	return;
    }
    $ir_fh->close;
    $self->log_event("Assignment Summary completed for $lib_name");
    return 1;
}

sub _parse_interesting_reads_file {
    my $self = shift;
    my $fh = shift;
    my $read_infos = shift;
    my $lib_name = basename($self->dir);
    my $percent_id_cutoff = 100; #SO REALLY PRINTING EVERYTHING
    my $orig_fasta_file = $self->dir.'/'.$lib_name.'.fa';
    foreach my $lineage (sort keys %$read_infos) {
	my $reads = {};
	$fh->print("$lineage\ttotal number of reads ".
		      $read_infos->{$lineage}->{'read_count'}."\t".
		      '['.$read_infos->{$lineage}->{'low_id'}.', '.
		      $read_infos->{$lineage}->{'high_id'}.']%'."\n\n");
	delete $read_infos->{$lineage}->{'read_count'};
	delete $read_infos->{$lineage}->{'low_id'};
	delete $read_infos->{$lineage}->{'high_id'};
	foreach my $read_name (keys %{$read_infos->{$lineage}}) {
	    $reads->{$read_name} = 1;
	    foreach (@{$read_infos->{$lineage}->{$read_name}}) {
		$fh->print($_."\n");
	    }
	}
	$fh->print("\n");
	my $io = Bio::SeqIO->new(-format => 'fasta', -file => $orig_fasta_file);
	while (my $seq = $io->next_seq) {
	    if (exists $reads->{$seq->primary_id}) {
		$fh->print('>'.$seq->primary_id."\n".$seq->seq."\n\n");
	    }
	}
    }

    return 1;
}

sub _print_seq_count_summary {
    my $self = shift;
    my $fh = shift;
    #PRINT HEADER FOR THE TABLE
    $fh->print($self->dir."\n\n");
    $fh->printf("%-9s", '#total');
    foreach (qw/ uniq % Filtered total% good total% BNHG total% BNNT total% TBXNT total% TBXVG total% /) {
	$fh->printf("%-9s", $_);
    }
    $fh->print("\n");
    my @counts = $self->_get_read_counts_from_fastas();
    
    $fh->printf("%-9.0f", $counts[0]); #$counts[0] = total seq count
    $fh->printf("%-9.0f", $counts[1]); #$counts[1] = cd-hit filtered seqs
    $fh->printf("%-9.1f", $counts[1] * 100 / $counts[0]);
    $fh->printf("%-9.0f", $counts[2]); #$counts[2] = cd-hit filtered bad seqs
    $fh->printf("%-9.1f", $counts[2] * 100 / $counts[0]);
    $fh->printf("%-9.0f", $counts[3]); #$counts[3] = cd-hit filtered good seqs
    $fh->printf("%-9.1f", $counts[3] * 100 / $counts[0]);
    $fh->printf("%-9.0f", $counts[3] - $counts[4]); #$counts[4] = HG filtered seqs
    $fh->printf("%-9.1f", ($counts[3] - $counts[4]) * 100 / $counts[0]);
    $fh->printf("%-9.0f", $counts[4] - $counts[5]); #$counts[5] = NT blastN filtered seqs
    $fh->printf("%-9.1f", ($counts[4] - $counts[5]) * 100 / $counts[0]);
    $fh->printf("%-9.0f", $counts[5] - $counts[6]); #$counts[6] = NT blastX filtered seqs
    $fh->printf("%-9.1f", ($counts[5] - $counts[6]) * 100 / $counts[0]);
    $fh->printf("%-9.0f", $counts[6] - $counts[7]); #$counts[7] = Viral blastX filtered seqs
    $fh->printf("%-9.1f", ($counts[6] - $counts[7]) * 100 / $counts[0]);
    $fh->print("\n\n\n");

    return 1;
}

sub _parse_assignment_report_file {
    my $self = shift;
    my $fh = shift;
    my $lib_name = basename($self->dir);
    my $assignment_report_file = $self->dir.'/'.$lib_name.'.AssignmentReport';
    unless (-s $assignment_report_file) {
	$self->log_event("Failed to find ".basename($assignment_report_file));
	return;
    }
    my $report_fh = IO::File->new("< $assignment_report_file") ||
	die "Can not create file handle for $assignment_report_file";
    my @table_row_names = qw/ category Bacteria Fungi Homo Mus Phage Viruses other unassigned /;
    my $read_infos = {};  my $lineage;
    while (my $line = $report_fh->getline) {
	#PRINT THE TABLE
	if ($line =~ /^\s+(\S+)\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+$/) {
	    if (grep (/^$1$/, @table_row_names)) {
		$fh->print($line);
	    }
	}
	#CREATE A HASH OF VIRUS LINEAGE INFO => [READ INFO LINES]
	if ($line =~ /\s+total\s+number\s+of\s+reads\s+\d+$/) {
	    chomp $line;
	    $lineage = $line;
	    $lineage =~ s/\s+total.*$//;
	}
	if ($line =~ /^\S+\s+\d+\s+(\S+)\s+\d+\s+/) {
	    chomp $line;
	    if ($1 =~ /^gi/) {
		my ($read_name) = $line =~ /^(\S+)\s+/;
		push @{$read_infos->{$lineage}->{$read_name}}, $line;
	    }
	}
    }
    $report_fh->close;
    $fh->print("\n\n");

    #FIGURING OUT THE LONGEST WORDED LINEAGE FOR PRINT ALIGNMENT PURPOSES
    my $string_length = 0;
    foreach my $lineage (keys %$read_infos) {
	$string_length = (length $lineage > $string_length) ? length $lineage : $string_length;
    }
    #PRINT LINEAGE, NUMBER OF READS [LOW_ID, HIGH_ID]%
    $fh->printf("%-".$string_length."s", '#LINEAGE_NAME'); #PRINT HEADER
    $fh->print("\tREAD_COUNT\t\t\tID Low, High\n");
    foreach my $lineage (sort keys %$read_infos) {
	my $high_id = 0;  my $low_id = 101;  my $read_count = 0;
	foreach my $read_name (keys %{$read_infos->{$lineage}}) {
	    $read_count++;
	    foreach my $read_info_line (@{$read_infos->{$lineage}->{$read_name}}) {
		my @tmp = split (/\t+/, $read_info_line);
		my $read_id = $tmp[0];    #READ NAME
		my $percent_id = $tmp[6]; #PERCENT MATCH
		$high_id = ($percent_id > $high_id) ? $percent_id : $high_id;
		$low_id = ($percent_id < $low_id) ? $percent_id : $low_id;
	    }
	}
	#print $read_count.' '.$low_id.' '.$high_id."\n";
	$fh->printf("%-".$string_length."s", $lineage);
	$fh->print("\ttotal number of reads: $read_count\t".'['.$low_id.', '.$high_id.']%'."\n");
	$read_infos->{$lineage}->{'low_id'} = $low_id;
	$read_infos->{$lineage}->{'high_id'} = $high_id;
	$read_infos->{$lineage}->{'read_count'} = $read_count;
    }

    return $read_infos;
}

sub _get_read_counts_from_fastas {
    my $self = shift;
    my $lib_name = basename($self->dir);
    my @read_counts;
    #LIST OF FASTA FILE EXTENSIONS
    my @file_extensions = qw/ fa fa.cdhit_out fa.cdhit_out.masked.badSeq fa.cdhit_out.masked.goodSeq
                              HGfiltered.fa BNfiltered.fa TBXNTfiltered.fa unassigned.fa /;
    foreach my $ext (@file_extensions) {
	my $fasta_file = $self->dir.'/'.$lib_name.'.'.$ext;
	my $read_count = 0;
	if (-s $fasta_file) {
	    my $io = Bio::SeqIO->new(-format => 'fasta', -file => $fasta_file);
	    while (my $seq = $io->next_seq) {
		$read_count++;
	    }
	}
	elsif (-e $fasta_file) { #NO DATA TO PROCESS
	    $self->log_event("No reads to process in ".basename($fasta_file));
	}
	else { #FILE IS MISSING .. THIS SHOULDN'T HAPPEN
	    $self->log_event("Can not open file: ".basename($fasta_file));
	    return;
	}
	push @read_counts, $read_count;
    }
    return @read_counts;
}

1;
