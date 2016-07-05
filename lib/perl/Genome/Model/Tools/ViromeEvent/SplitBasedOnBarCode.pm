package Genome::Model::Tools::ViromeEvent::SplitBasedOnBarCode;

use strict;
use warnings;

use Genome;
use Data::Dumper;
use IO::File;
use File::Basename;
use Bio::SeqIO;

class Genome::Model::Tools::ViromeEvent::SplitBasedOnBarCode{
    is => 'Genome::Model::Tools::ViromeEvent',
    has => [
        barcode_file => {
            doc => 'file of reads to be checked for contamination',
            is => 'String',
            is_input => 1,
        },
        fasta_file => {
            doc => 'file of reads to be checked for contamination',
            is => 'String',
            is_input => 1,
        },
    ],
};

sub help_brief {
    'Creates a list of directories and fasta files based on querying barcode file against fasta file';
}

sub help_detail {
    return <<"EOS"
This script accepts a 454 sequencing run output fasta file and sort the 
sequences based on their barcode. Each group of barcoded sequence 
represents a library. A "undecodable" file/dir will be created for 
sequences that do not have exact match to any of the used barcode.
One directory will be created for each library, which will hold a .fa file
holding all the sequences from this library.

In this version, primer B sequence at both ends were stripped off, check for
number of samples, duplication of primer B sequence.

In this version, calculate percentage of sequences in each length category
0-50, 51 - 100 etc.

In this version, both barcodes from 5 prime and 3 prime are used for 
decoding.
EOS
}

sub execute {
    my $self = shift;
    $self->log_event("Spliting based on barcode ".basename($self->fasta_file));

    if( not -s $self->barcode_file ) {
        $self->log_event("Failed to find barcode file or file is zero size: ".$self->barcode_file);
        return;
    }
    if( not -s $self->fasta_file ) {
        $self->log_event("Failed to find fasta file or file is zero size: ".$self->fasta_file);
        return;
    }
    if( not -d $self->dir ) {
        mkdir $self->dir;
        $self->log_event("Failed to find or create dir: ".$self->dir) and return if not
            -d $self->dir;
    }

    #CREATE A HASH LINKING SAMPLE NAME TO BARCODE SEQUENCE
    my $barcodes;
    unless ($barcodes = $self->_parse_barcode_file()) {
	$self->log_event("Failed to parse barcode file");
	return;
    }
    #CREATE GENERAL METRICS, CLIP PB SEQ FROM FASTA
    #FILTER FASTA BY BARCODE SEQUENCES
    my ($metrics, $fasta);
    unless (($metrics, $fasta) = $self->_filter_fasta_file($barcodes)) {
	$self->log_event("Failed to filter fasta file");
	return;
    }
    #SPARATE READS BASED ON BARCODE AND CREATE SAMPLES DIRS
    #FOR REST OF PIPELINE
    unless ($self->_create_run_samples($fasta, $barcodes)) {
	$self->log_event("Failed to create run samples");
	return;
    }
    #CREATE/PRINT GENERAL METRICS
    unless ($self->_print_metrics($metrics)) {
	$self->log_event("Failed to print metrics");
	return;
    }
    unless ($self->_print_distribution_metrics($fasta, $barcodes)) {
	$self->log_event("Failed to print dist metrics");
	return;
    }

    $self->log_event("Split based on barcode completed");
    return 1;
}

sub _print_metrics {
    my $self = shift;
    my $metrics = shift;
    my $analysis_out = $self->dir.'/Analysis_ReadStatistics_'.basename($self->dir);
    my $out_fh = IO::File->new("> $analysis_out") ||
	die "Can not create file handle for $analysis_out";
    $out_fh->print($self->fasta_file."\n".
		   "total number of samples: ".$metrics->{sample_count}."\n".
		   "total number of sequences: ".$metrics->{read_count}."\n" );

    #PRIMER B FOUND AT BOTH ENDS
    my $pb_at_both_ends_count = (exists $metrics->{reads_with_3and5prime_pb}) ? $metrics->{reads_with_3and5prime_pb} : '0';
    my $pb_at_both_ends_ratio = sprintf("%.1f", $pb_at_both_ends_count * 100 / $metrics->{read_count});
    $out_fh->print("\nNumber of sequences with 5' and 3' Primer B: $pb_at_both_ends_count ($pb_at_both_ends_ratio%)\n");
    #DECODED BY BOTH ENDS
    my $bc_at_both_ends_count = (exists $metrics->{pb_at_both_ends_decoded_by_5_and_3prime}) ?
	$metrics->{pb_at_both_ends_decoded_by_5_and_3prime} : '0';
    my $bc_at_both_ends_ratio = sprintf("%.1f", $bc_at_both_ends_count * 100 / $metrics->{read_count});
    $out_fh->print("\tDecoded with 5' and 3' ends: $bc_at_both_ends_count\n");
    #DECODED BY 5' END
    my $bc_at_5_end_count = (exists $metrics->{pb_at_both_ends_decoded_by_5prime}) ?
	$metrics->{pb_at_both_ends_decoded_by_5prime}: '0';
    $out_fh->print("\tDecoded with 5' end: $bc_at_5_end_count\n");
    #DECODED BY 3' END
    my $bc_at_3_end_count = (exists $metrics->{pb_at_both_ends_decoded_by_3prime}) ?
	$metrics->{pb_at_both_ends_decoded_by_3prime}: '0';
    $out_fh->print("\tDecoded with 3' end: $bc_at_3_end_count\n");
    #UNDECODED WITH NO VALID BARCODES:
    my $bc_at_both_ends_undecoded_count = (exists $metrics->{pb_at_both_ends_with_different_barcodes}) ?
	$metrics->{pb_at_both_ends_with_different_barcodes} : '0';
    $out_fh->print("\tUndecoded with no valid barcodes: $bc_at_both_ends_undecoded_count\n");
    #UNDECODED WITH DIFFERENT VALID BARCODES:
    my $bc_at_both_ends_dif_bc_count = (exists $metrics->{pb_at_both_ends_no_valid_bc}) ?
	$metrics->{pb_at_both_ends_no_valid_bc} : '0';
    $out_fh->print("\tUndecoded with different valid barcodes: $bc_at_both_ends_dif_bc_count\n");


    #PRIMER B FOUND AT 5' END ONLY
    my $pb_at_5_end_count = (exists $metrics->{reads_with_5prime_pb}) ? $metrics->{reads_with_5prime_pb} : '0';
    my $pb_at_5_end_ratio = sprintf("%.1f", $pb_at_5_end_count * 100 / $metrics->{read_count});
    $out_fh->print("Number of sequences with 5' Primer B: $pb_at_5_end_count ($pb_at_5_end_ratio%)\n");
    #DECODED BY 5' END
    my $pb_at_5_end_decoded_by_5_end_count = (exists $metrics->{pb_at_5end_decoded_by_5prime}) ?
	$metrics->{pb_at_5end_decoded_by_5prime} : '0';
    $out_fh->print("\tDecoded with 5' end: $pb_at_5_end_decoded_by_5_end_count\n");
    #UNDECODED
    my $pb_at_5_end_undecoded_count = (exists $metrics->{pb_at_5end_undecoded}) ? $metrics->{pb_at_5end_undecoded} : '0';
    $out_fh->print("\tUndecoded: $pb_at_5_end_undecoded_count\n");

    #PRIMER B FOUND AT 3' END ONLY
    my $pb_at_3_end_count = (exists $metrics->{reads_with_3prime_pb}) ? $metrics->{reads_with_3prime_pb} : '0';
    my $pb_at_3_end_ratio = sprintf("%.1f", $pb_at_3_end_count * 100 / $metrics->{read_count});
    $out_fh->print("Number of sequences with 3' Primer B: $pb_at_3_end_count ($pb_at_3_end_ratio%)\n");
    #DECODED BY 3' END
    my $pb_at_3_end_decoded_by_3_end_count = (exists $metrics->{pb_at_3end_decoded_by_3prime}) ?
	$metrics->{pb_at_3end_decoded_by_3prime} : '0';
    $out_fh->print("\tDecoded with 3' end: $pb_at_3_end_decoded_by_3_end_count\n");
    #UNDECODED
    my $pb_at_3_end_undecoded_count = (exists $metrics->{pb_at_3end_undecoded}) ? $metrics->{pb_at_3end_undecoded} : '0';
    $out_fh->print("\tUndecoded: $pb_at_3_end_undecoded_count\n");


    #NO PRIMER B AT EITHER ENDS
    my $pb_at_no_ends_count = (exists $metrics->{reads_with_no_pb}) ? $metrics->{reads_with_no_pb} : '0';
    my $pb_at_no_ends_ratio = sprintf("%.1f", $pb_at_no_ends_count * 100 / $metrics->{read_count});
    $out_fh->print("Number of sequences with no Primer B: $pb_at_no_ends_count ($pb_at_no_ends_ratio%)\n");
    #UNDECODED
    $out_fh->print("\tUndecoded: $pb_at_no_ends_count\n\n");

    #READS DECODED
    my $reads_decoded_count = 0;
    my $reads_decoded_ratio = 0;
    if (exists $metrics->{reads_decoded}) {
	$reads_decoded_count = $metrics->{reads_decoded};
	$reads_decoded_ratio = sprintf("%.1f", $reads_decoded_count * 100 / $metrics->{read_count});
    }
    $out_fh->print("number of sequences decoded: $reads_decoded_count ( $reads_decoded_ratio% )\n");
    #AVERAGE READ LENGTH
    #my $average_length = sprintf ("%.1f", $metrics->{total_sequence_length} / $metrics->{read_count});
    my $average_length = int ($metrics->{total_sequence_length} / $metrics->{read_count});
    $out_fh->print("average length of sequences: $average_length\n\n");

    $out_fh->close;
    return 1;
}

sub _print_distribution_metrics {
    my $self = shift;
    my $fasta = shift;
    my $barcodes = shift;
    my $analysis_out = $self->dir.'/Analysis_ReadStatistics_'.basename($self->dir);
    my $out_fh = (-s $analysis_out) ? IO::File->new(">> $analysis_out") : IO::File->new("> $analysis_out");
    unless ($out_fh) {
	die "Can not create file handle to append analysis file"
    }
    my $dist_increment = 50;
    my $dist_number = 8;

    #CALCUALTE DISTRIBUTION OF SEQUENCE LENGTHS IN INCREMENTS OF 50 BPS
    $out_fh->print("#######################################################\n".
                    "distribution of sequence lengths in the run:\n");
    my $str = "barcode\ttotal#\tAvgLen\t";
    #PRINT DISTRIBUTION SCALE
    #BARCODE TOTAL#  AvgLen          0       50      100     150     200     250     300     350     400     sample
    for (my $i = 0; $i <= $dist_number; $i++) {
	$str .= sprintf("%7.0f", $i * $dist_increment);
    }
    $out_fh->print($str."\t\tsample\n");
    my $dists = {};
    foreach my $barcode (sort keys %$fasta) {
	foreach my $read_name (keys %{$fasta->{$barcode}}) {
	    my $seq_length = length $fasta->{$barcode}->{$read_name};
	    my $dist_hit = int ($seq_length / $dist_increment);
	    #DISTRIBUTION FOR INDIVIDUAL SAMPLES
	    $dist_hit = ($dist_hit > $dist_number) ? $dist_number : $dist_hit;
	    $dists->{$barcode}->{$dist_hit}++;
	    $dists->{$barcode}->{'total_seq_lengths'} += $seq_length;
	    $dists->{$barcode}->{'total_read_count'}++;
	    $dists->{$barcode}->{'read_less_than_100bp'}++ if $seq_length <= 100;
	    #DISTRIBUTIONS FOR ALL DATA
	    $dists->{'total'}->{$dist_hit}++;
	    $dists->{'total'}->{'total_seq_lengths'} += $seq_length;
	    $dists->{'total'}->{'total_read_count'}++;
	}
    }
#   print Dumper $dists;
    #PRINT TOTAL DIST STATS FIRST
    my $string = "total\t".$dists->{'total'}->{'total_read_count'}."\t".
	         sprintf ("%.1f", $dists->{'total'}->{'total_seq_lengths'} / $dists->{'total'}->{'total_read_count'})."\t";
    for (my $i = 0; $i <= $dist_number; $i++) {
	my $val = (exists $dists->{'total'}->{$i}) ?
	    sprintf ("%.1f", $dists->{'total'}->{$i} * 100 / $dists->{'total'}->{'total_read_count'}) : '0.0';
	$string .= sprintf ("%7.1f", $val);
    }
    $out_fh->print($string."\n");
    #PRINT REST OF DIST STATS
    foreach my $bc (sort keys %$dists) {
	next if $bc eq 'total' || $bc eq 'undecodable';
	$string = "$bc\t".$dists->{$bc}->{'total_read_count'}."\t".
	          sprintf ("%.1f", $dists->{$bc}->{'total_seq_lengths'} / $dists->{$bc}->{'total_read_count'})."\t";
	for (my $i = 0; $i <= $dist_number; $i++) {
	    my $val = (exists $dists->{$bc}->{$i}) ?
		sprintf ("%.1f", $dists->{$bc}->{$i} * 100 / $dists->{$bc}->{'total_read_count'}) : '0.0';
	    $string .= sprintf ("%7.1f", $val);
	}
	$string .= "\t\t".$barcodes->{$bc} if exists $barcodes->{$bc};
	$out_fh->print($string."\n");
    }
    $out_fh->print("End of Distribution\n\n".
	           "#######################################################\n".
	           "Sample needs to be reamplified:\n");
    #SAMPLES FOR RE-SEQUENCING
    foreach my $bc (sort keys %$dists) {
	next if $bc eq 'total' || $bc eq 'undecodable';
	if (exists $dists->{$bc}->{'read_less_than_100bp'}) {
	    my $ratio = int ($dists->{$bc}->{'read_less_than_100bp'} * 100 / $dists->{$bc}->{'total_read_count'});
	    #IF >= 15% OF READS ARE LESS THAN 100 BP PUT ON RESEQUENCE LIST
	    $out_fh->print($barcodes->{$bc}."\n") if $ratio >= 15;
	}
    }
    $out_fh->print("End of sample list\n");
    $out_fh->close;
    return 1;
}

sub _create_run_samples {
    my $self = shift;
    my $fasta = shift;
    my $sample_names = shift;
    my $undecod_sample_name = basename $self->dir;
    $undecod_sample_name .= '_undecodable';
    foreach my $barcode (keys %$fasta) {
	my $sample_name = ($barcode eq 'undecodable') ? $undecod_sample_name : $sample_names->{$barcode};
	unless ($sample_name) {
	    $self->log_event("Can not get sample name for $barcode");
	    return;
	}
	my $sample_dir = $self->dir.'/'.$sample_name;
	unless (-d $sample_dir) {
	    system("mkdir $sample_dir");
	}
	my $sample_fa = $sample_dir.'/'.$sample_name.'.fa';
	if (-s $sample_fa) {
	    #VERIFY ALL READS ARE THERE IF ALL READS ARE THERE
	    if ($self->_verify_all_reads_present($sample_fa, $fasta->{$barcode})) {
		next;
	    }
	    unlink $sample_fa;
	    #ALSO REMOVE EVERYTHING ELSE THAT GETS CREATED AFTER THIS STEP
	}
	my $fa_out = Bio::SeqIO->new(-format => 'fasta', -file => ">$sample_fa");
	foreach my $read_name (keys %{$fasta->{$barcode}}) {
	    my $seq_obj = Bio::Seq->new(-display_id => $read_name, -seq => $fasta->{$barcode}->{$read_name});
	    $fa_out->write_seq($seq_obj);
	}
    }
    return 1;
}

sub _verify_all_reads_present {
    my $self = shift;
    my $fasta_file = shift;
    my $fasta = shift;
    my @existing_reads;
    my $io = Bio::SeqIO->new(-format => 'fasta', -file => $fasta_file);
    while (my $seq = $io->next_seq) {
	if (exists $fasta->{$seq->primary_id}) {
	    push @existing_reads, $seq->primary_id;
	}
	else {
	    $self->log_event("Failed to find ".$seq->primary_id." in existing fasta file");
	    return;
	}
    }
    unless (scalar (keys %$fasta) == scalar @existing_reads) {
	$self->log_event("Number of reads did not agree");
	return;
    }
    return 1;
}

sub _pb_seq {
    my $self = shift;
    return 'GTTTCCCAGTCACGATA';
}

sub _comp_pb_seq {
    my $self = shift;
    my $rev_pb_seq = reverse $self->_pb_seq();
    $rev_pb_seq =~ tr/ACGTacgt/TGCAtgca/;
    return $rev_pb_seq;
}

sub _get_barcode_lengths {
    my $self = shift;
    my $barcode = shift;
    my $barcode_length = 0;
    #BARCODE->{BARCODE} = SAMPLE_NAME
    foreach my $barcode (keys %$barcode) {
	#ALL BARCODES SHOULD BE THE SAME LENGTH
	if ($barcode_length > 0) {
	    unless ($barcode_length == length $barcode) {
		$self->log_event("$barcode is not sample length as previous one");
		return;
	    }
	}
	$barcode_length = length $barcode;
    }
    return $barcode_length;
}

sub _filter_fasta_file {
    my $self = shift;
    my $barcodes = shift;
    my $metrics = {}; my $fasta = {};
    $metrics->{sample_count} = 0;
    $metrics->{read_count} = 0;
    my $barcode_length = $self->_get_barcode_lengths($barcodes);
    my $pb_seq = $self->_pb_seq;
    my $comp_pb_seq = $self->_comp_pb_seq();
    my $io = Bio::SeqIO->new(-format => 'fasta', -file => $self->fasta_file);
    while (my $seq = $io->next_seq) {
	my $read_seq = $seq->seq; #WILL GET PROCESSED
	my $found_5prime_pb; my $found_3prime_pb;
	my $code_at_5end;    my $code_at_3end;
	my $sequence_decoded = 0;
	my $barcode; #BARCODE SEQ DERIVED FROM FASTA SEQ
	if ($read_seq =~ /$pb_seq/) {
	    $code_at_5end = substr("$`", -1*$barcode_length);
	    $read_seq = "$'";
	    $found_5prime_pb = 1;
	}
	if ($read_seq =~ /$comp_pb_seq/) {
	    $code_at_3end = substr("$'", 0, $barcode_length);
	    $read_seq = "$`";
	    $found_3prime_pb = 1;
	    #COMPLEMENT THE 3 END SEQ
	    $code_at_3end = reverse $code_at_3end;
	    $code_at_3end =~ tr/ACGTacgt/TGCAtgca/;
	}
	#IF THERE IS NO SEQUENCE TO BE DERIVED THEN NEXT
	unless ($read_seq) {
	    $self->log_event("No valid sequence left after clipping barcode sequence for ".$seq->primary_id);
	    next;
	}
        #NOT ENOUGH SEQUENCE LEFT FOR FURTHER PROCESSING
        if ( not length $read_seq >= 10 ) {
            $self->log_event('Skipping '.$seq->primary_id.' only '.(length $read_seq).' bases remain after clipping pb sequence');
            next;
        }
	#IF 5' AND 3' PRIMER B MATCHES
	if ($found_5prime_pb && $found_3prime_pb) {
	    $metrics->{reads_with_3and5prime_pb}++;
	     #IF VALID 3' AND 5' CODES EXIST
	    if (exists $barcodes->{$code_at_5end} && exists $barcodes->{$code_at_3end}) {
		#IF THEY'RE THE SAME .. SEQ IS DECODED
		if ($barcodes->{$code_at_5end} eq $barcodes->{$code_at_3end}) {
		    $sequence_decoded = 1;
		    $barcode = $code_at_5end;
		    $metrics->{pb_at_both_ends_decoded_by_5_and_3prime}++;
		}
		#IF THEY'RE DIFFERENT .. SEQ IS NOT DECODED
		else {
		    $metrics->{pb_at_both_ends_with_different_barcodes}++;
		}
	    }
	    #IF ONLY VALID 5' EXISTS
	    elsif (exists $barcodes->{$code_at_5end}) {
		$sequence_decoded = 1;
		$barcode = $code_at_5end;
		$metrics->{pb_at_both_ends_decoded_by_5prime}++;
	    }
	    #IF ONLY VALID 3' EXISTS
	    elsif (exists $barcodes->{$code_at_3end}) {
		$sequence_decoded = 1;
		$barcode = $code_at_3end;
		$metrics->{pb_at_both_ends_decoded_by_3prime}++;
	    }
	    #IF NEITHER ARE VALID BARCODES
	    else {
		#UNDECODED
		$metrics->{pb_at_both_ends_no_valid_bc}++;
	    }
	}
	#IF 5' PRIMER B MATCH ONLY
	elsif ($found_5prime_pb) {
	    $metrics->{reads_with_5prime_pb}++;
	    if (exists $barcodes->{$code_at_5end}) {
		$sequence_decoded = 1;
		$barcode = $code_at_5end;
		$metrics->{pb_at_5end_decoded_by_5prime}++;
	    }
	    else {
		$metrics->{pb_at_5end_undecoded}++;
	    }
	}
	#IF 3' PRIMER B MATCH ONLY
	elsif ($found_3prime_pb) {
	    $metrics->{reads_with_3prime_pb}++;
	    if (exists $barcodes->{$code_at_3end}) {
		$sequence_decoded = 1;
		$barcode = $code_at_3end;
		$metrics->{pb_at_3end_decoded_by_3prime}++;
	    }
	    else {
		$metrics->{pb_at_3end_undecoded}++
	    }
	}
	#NO MATCH
	else {
	    $metrics->{reads_with_no_pb}++;
	}
	#METRICS FOR LATER CALCULATIONS
	$metrics->{read_count}++;
	$metrics->{total_sequence_length} += length $seq->seq;
#	push @{$metrics->{lengths_of_seqs}}, length $seq->seq;
	#CREATE A FASTA HASH THAT CAN BE LOOKED UP BY BARCODE
	if ($sequence_decoded == 1) {
	    $fasta->{$barcode}->{$seq->primary_id} = $read_seq;
	    $metrics->{reads_decoded}++;
	}
	else {
	    $fasta->{undecodable}->{$seq->primary_id} = $read_seq;
	    $metrics->{undecoded}++;
	}
    }
    #GET A COUNT OF NUMBER OF SAMPLES
    foreach my $barcode (keys %$fasta) {
	next if $barcode eq 'undecodable';
	$metrics->{sample_count}++;
    }

    print Dumper $metrics;
    return $metrics, $fasta;
}

sub _parse_barcode_file {
    my $self = shift;
    my %barcode;
    unless (-s $self->barcode_file) {
	$self->log_event("Can not fine barcode file ".basename($self->barcode_file));
	return; #SHOULD EXIT .. BUT CAN'T IN WORK FLOW
    }
    my $pb_seq = $self->_pb_seq();
    my $fh = IO::File->new("< ".$self->barcode_file) ||
	die "Can not create file handle for barcode file";
    while (my $line = $fh->getline) {
	next if $line =~ /^\s+$/ || $line =~ /#/;
	chomp $line;
	my @tmp = split(/\s+/, $line);
	my $sample_number = $tmp[0];
	#$temp[1] ?? NOT SURE WHAT THIS IS 
	my $tag = $tmp[2]; #BARCODE SEQUENCE
        my $lib = $tmp[3]; #LIBRARY NAME
	$tag =~ s/$pb_seq//;
	next unless $tag; #SKIP IF TAG SEQ IS COMPOSED ENTIRELY OF PB SEQ
	$lib =~ s/\W+/_/g;
	my $sample_name = "S".$sample_number."_".$lib;
	$self->log_event("tag = $tag, sample_name = $sample_name");
	if (exists $barcode{$tag}) {
	    $self->log_event("$tag is duplicated");
	    return; #SHOULD EXIT
	}
	else {
	    $barcode{$tag} = $sample_name;
	}
    }
    $fh->close;
    return \%barcode;
}

1;
