package Genome::Model::Tools::ViromeEvent::Assignment::Generate;

use strict;
use warnings;

use Genome;
use Workflow;
use Switch;
use IO::File;
use File::Basename;
use Bio::SeqIO;
use Bio::SearchIO;

class Genome::Model::Tools::ViromeEvent::Assignment::Generate{
    is => 'Genome::Model::Tools::ViromeEvent',
};

sub help_brief {
    "gzhao's reporting for virome"
}
sub help_detail {
    "Tool to assign virome screen hits to organism catagories";
}

sub execute {
    my $self = shift;
    $self->log_event("Starting generate assignments");
    my $report_file = $self->dir.'/Analysis_Report_'.basename($self->dir);
    my $report_fh = IO::File->new("> $report_file") ||
	die "Can not create file handle for $report_file";
    my $sep = "***************************************************************************************\n\n";
    unless ($self->_print_intro($report_fh)) {
	$self->log_event("Failed to print introduction");
	return;
    }
    $report_fh->print($sep);
    unless ($self->_generate_sample_description($report_fh)) {
	$self->log_event("Failed to generate simple description");
	return;
    }
    $report_fh->print($sep);
    unless ($self->_generate_sequence_report($report_fh)) {
	$self->log_event("Failed generate sequence report");
	return;
    }
    $report_fh->print($sep);
    unless ($self->_generate_assignment_summary($report_fh)){
	$self->log_event("Failed to generate assignment summary");
	return;
    }
    $report_fh->print($sep);
    unless ($self->_report_interesting_reads($report_fh)) {
	$self->log_event("Failed to generate interesting reads report");
	return;
    }
    $report_fh->close;
    #GENERATE FASTA OF ALL UNASSIGNED READS
    unless ($self->_create_main_unassigned_reads_file()) {
	$self->log_event("Failed to create main unassigned reads file");
	return;
    }
    $self->log_event("Completed virome screening for ".basename($self->dir));
    return 1;
}

sub _create_main_unassigned_reads_file {
    my $self = shift;
    my $unassigned_reads = $self->dir.'/All_unassigned.fa';
    my $out = Bio::SeqIO->new(-format => 'fasta', -file => ">$unassigned_reads");
    foreach my $sample_name ($self->_get_sample_dir_names()) {
	my $file = $self->dir.'/'.$sample_name.'/'.$sample_name.'.unassigned.fa';
	if (-s $file) {
	    my $io = Bio::SeqIO->new(-format => 'fasta', -file => $file);
	    while (my $seq = $io->next_seq) {
		$out->write_seq($seq);
	    }
	}
	elsif (-e $file) { #NO DATA TO PROCESS
	    $self->log_event("No unassigned data to process for $sample_name");
	    next;
	}
	else { #FILE MISSING
	    $self->log_event("Failed to find unassigned.fa file for $sample_name");
	    return;
	}
    }
    return 1;
}

sub _report_interesting_reads {
    my $self = shift;
    my $fh = shift; #MAIN ANALYSIS FILE FH
    $fh->print("Interesting Reads\n\n");
    foreach my $dir_name($self->_get_sample_dir_names()) {
	$fh->print("Sample: $dir_name\n\n");
	my $file = $self->dir.'/'.$dir_name.'/'.$dir_name.'.InterestingReads';
	if (-s $file) {
	    my $reads_fh = IO::File->new("< $file") ||
		die "Can not create file handle for $file";
	    while (my $line = $reads_fh->getline) {
		$fh->print($line);
	    }
	    $reads_fh->close;
	}
	elsif (-e $file) {
	    $fh->print("No Interesting reads found\n\n");
	    $self->log_event("No interesting reads for $dir_name");
	}
	else {
	    $self->log_event("Missing interesting reads file ".basename($file));
	    return;
	}

	$fh->print("###############################################################\n\n");
    }
    $fh->print("End of Interesting Reads");
    return 1;
}

sub _generate_assignment_summary {
    my $self = shift;
    my $fh = shift;
    $fh->print("Assignments in each sample:\n\n");
    foreach my $dir_name ($self->_get_sample_dir_names()) {
	my $summary_file = $self->dir.'/'.$dir_name.'/'.$dir_name.'.AssignmentSummary';
	unless (-s $summary_file) {
	    $self->log_event("Missing or blank summary file ".basename($summary_file));
	    return;
	}
	my $sum_fh = IO::File->new("< $summary_file") ||
	    die "Can not open file: $summary_file";
	while (my $line = $sum_fh->getline) {
	    $fh->print($line);
	}
	$fh->print("###############################################################\n\n");
	$sum_fh->close;
    }
    $fh->print("End of Assignments\n\n");
    return 1;
}

sub _get_sample_dir_names {
    my $self = shift;
    my @sample_dirs;
    foreach (glob($self->dir."/*")) {
	next unless -d $_;
	my $lib_name = basename($_);
	#DETERMINE WHETHER DIR IS SAMPLE DIR BY CHECKING FOR CERTAIN FILES
	unless (-s $self->dir.'/'.$lib_name.'/'.$lib_name.'.fa') {
	    next;
	}
	#TODO - MAY HAVE TO CHECK FOR ADDITIONAL FILES/DIRS
	push @sample_dirs, $lib_name;
    }
    return @sample_dirs;
}

sub _generate_sequence_report {
    my $self = shift;
    my $fh = shift;
    $fh->print("Sequence Report\n\n".$self->dir."\n\n");
    $fh->printf("%40s%10s%10s%10s%10s%10s%10s%10s%10s%10s%10s%10s\n",
                 'sampleName','total','uniq','%','Filtered','%','good','%','BNassign','%','TBLASTX','%');
    my $counts;
    my $convert_table = $self->_convert_to_table_name();
    my @sample_dirs = $self->_get_sample_dir_names();
    foreach my $file_name ($self->_get_sample_dir_names()) {
	foreach my $ext ($self->_fasta_file_exts()) {
	    my $fa_file = $self->dir.'/'.$file_name.'/'.$file_name.'.'.$ext;
	    my $convert_name = $convert_table->{$ext};
	    my $seq_count = $self->_get_fasta_seq_count($fa_file);
	    $counts->{$file_name}->{$convert_name} = $seq_count;
	    $counts->{'total'}->{$convert_name}+= $seq_count;
	}
    }
    foreach my $lib_name (sort keys %$counts) {
	my $total_count = $counts->{$lib_name}->{'input_fa'};
	my $uniq_count = $counts->{$lib_name}->{'cd_hit'};
	my $uniq_ratio = sprintf("%.1f", $uniq_count * 100 / $total_count);
	my $filtered = $counts->{$lib_name}->{'cd_hit_bad'};
	my $filtered_ratio = sprintf("%.1f", $filtered * 100 / $total_count);
	my $good = $counts->{$lib_name}->{'cd_hit_good'};
	my $good_ratio = sprintf("%.1f", $good * 100 / $total_count);
	my $bn_assigned = $good - $counts->{$lib_name}->{'blastn'};
	my $bn_assigned_ratio = sprintf("%.1f", $bn_assigned * 100 / $total_count);
	#BLASTN FILTERED AND PASSED ONTO BLASTX
	my $bx_assigned = $counts->{$lib_name}->{'blastn'};
	my $bx_assigned_ratio = sprintf ("%.1f", $bx_assigned * 100 / $total_count);
	$fh->printf("%40s%10s%10s%10s%10s%10s%10s%10s%10s%10s%10s%10s\n",
	    $lib_name, $total_count, $uniq_count, $uniq_ratio, $filtered, $filtered_ratio, $good, $good_ratio, $bn_assigned, $bn_assigned_ratio, $bx_assigned, $bx_assigned_ratio);

    }
    $fh->printf("\nEnd of Sequence Report\n\n");

    return 1;
}

sub _fasta_file_exts {
    return qw/ fa fa.cdhit_out fa.cdhit_out.masked.badSeq fa.cdhit_out.masked.goodSeq BNfiltered.fa /;
}

sub _convert_to_table_name {
    return {
	'fa'                          => 'input_fa',
	'fa.cdhit_out'                => 'cd_hit',
	'fa.cdhit_out.masked.badSeq'  => 'cd_hit_bad',
	'fa.cdhit_out.masked.goodSeq' => 'cd_hit_good',
	'BNfiltered.fa'               => 'blastn',
    };
}

sub _generate_sample_description {
    my $self = shift;
    my $fh = shift;
    $fh->print("Summary:\n\n".$self->dir."\n");
    $fh->printf("%40s%25s%40s\n", 'ViralRead', 'PercentIDrange', 'IdentifiedVirus');
    foreach (glob($self->dir."/*")) {
	next unless -d $_; #ONLY SAMPLE DIRS ARE DIRS
	my $lib_name = basename($_);
	my $input_fa = $self->dir.'/'.$lib_name."/$lib_name.fa";
	unless (-s $input_fa) {
	    $self->log_event("No or blank input fasta found for $lib_name");
	    return;
	}
	my $seq_count = $self->_get_fasta_seq_count($input_fa);
	$fh->printf("\t%-40s%-25s\n", $lib_name, $seq_count);
	#PARSE ASSIGNMENT SUMMARY FILE
	my $asum_file = $self->dir.'/'.$lib_name."/$lib_name.AssignmentSummary";
	unless (-s $asum_file) {
	    $self->log_event("No or blank file ".basename($asum_file));
	    return;
	}
	my $asum_fh = IO::File->new("< $asum_file") ||
	    die "Can not create file handle for $asum_file";
	while (my $line = $asum_fh->getline) {
	    next unless $line =~ /\s+number\s+of\s+reads/;
	    chomp $line;
	    my ($range) = $line =~ /number\s+of\s+reads:\s+\d+\s+(.*)$/; #GET ID RANGE
	    my ($read_count) = $line =~ /number\s+of\s+reads:\s+(\d+)/;  #GET NUMBER OF READS
	    my ($lineage) = $line =~ /^(.*)\s+total\s+number\s+of/;      #GET LINEAGE INFO
	    $lineage =~ s/\s+$//;
	    my @tmp = split (';', $lineage);
	    #VIRUS NAME = $tmp[-1];
	    $fh->printf("%40s%25s%40s\n", $read_count, $range, $tmp[-1]);

	}
	$asum_fh->close;
    }
    $fh->print("End of summary\n\n");
    return 1;
}

sub _get_fasta_seq_count {
    my $self = shift;
    my $fa = shift;
    my $seq_count = 0;
    my $io = Bio::SeqIO->new(-format => 'fasta', -file => $fa);
    return unless $io;
    while (my $seq = $io->next_seq) {
	$seq_count++;
    }
    return $seq_count;
}

sub _print_intro {
    my $self = shift;
    my $fh = shift;
    $fh->print ("How to read this file:\n".
		"For the summary section:\n".
		"column 1: sample number and name\n".
		"column 2: sample description\n".
		"column 3: total number of sequences obtained for this sample\n\n".
		"If there is any viral sequence identified in this sample, it will show up under the information \n".
		"of this sample. There are 3 columns to describe the viral reads identified in this sample:\n".
		"column 1: number of viral reads\n".
		"column 2: range of percentage identity to blast hits. Some times one sequence can hit multiple \n".
		"sequence in the nt database and gives a range of percent identity.\n".
		"column 3: name of the virus\n\n".
		"For the sequence report section:\n".
		"column 1: sample number and name\n".
		"column 2: total number of reads obtained for this sample\n".
		"column 3: number of unique reads after removing redundancy\n".
		"column 4: percentage of unique reads (unique reads divided by total number of reads)\n".
		"column 5: number of Filtered reads (after repeat masking, reads does not have at least 50bp \n".
		"consecutive sequence without N)\n".
		"column 6: percentage of Filtered reads (Filtered reads divided by total number of reads)\n".
		"column 7: number of good reads\n".
		"column 8: percentage of good reads\n".
		"column 9: number of reads assigned by blastn (BNassign)\n".
		"column 10: percentage of reads assigned by blastn\n".
		"column 11: number of reads goes to TBLASTX\n".
		"column 12: percentage of sequences goes to TBLASTX\n\n".
		"For the Assignment in each sample section:\n".
		"Describes the number of sequences assigned to each category.\n\n".
		"If there is a Interesting read section, following are the description of the fields:\n".
		"Sample number and name\n".
		"Viral lineage\n".
		"Description of the reads and the top hit (the fields are: QueryName, Querylength, HitName, HitLen, \n".
		"HitDesc, Alignment Length, percent Identity, HitStart, HitEnd, e value)\n".
		"Sequence of the reads\n\n".
		"*The criteria for a viral lineage to appear in the \"Interesting Read\" section is the lowest \n".
		"percent identity is below 90% (Anellovirus is excluded)\n\n");
    return 1;
}

1;
