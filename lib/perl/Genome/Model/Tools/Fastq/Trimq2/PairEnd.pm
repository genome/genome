package Genome::Model::Tools::Fastq::Trimq2::PairEnd;

use strict;
use warnings;

use Genome;
use File::Temp;
use File::Basename;

class Genome::Model::Tools::Fastq::Trimq2::PairEnd{
    is => 'Genome::Model::Tools::Fastq::Trimq2',
    has => [
        pair1_fastq_file  => {
            is  => 'Text',
            doc => 'the pair end fastq file 1',
        },
        pair2_fastq_file  => {
            is  => 'Text',
            doc => 'the pair end fastq file 2',
        },
    ],
    has_optional => [
        pair1_out_file    => {
            is  => 'Text',
            doc => 'the file path to use for the pair1 output file, default is xxx1.trimq2.fastq under same dir as pair1 fastq file',
        },
        pair2_out_file    => {
            is  => 'Text',
            doc => 'the file path to use for the pair2 output file, default is xxx2.trimq2.fastq under same dir as pair1 fastq file',
        },
        pair_as_frag_file => {
            is  => 'Text',
            doc => 'the file path to use for the pair as fragment fastq file, default is trimq2.pair_as_fragment.fastq under same dir as pair1 fastq file',
        },
    ],
};

sub help_synopsis {
    return <<EOS
gmt fastq trimq2 pair-end --pair1-fastq-file=lane1.fastq --pair1-out-file=lane1.trimmed.fastq
EOS
}

sub help_detail {
    return <<EOS
Trims fastq reads with phred quality 2 as either Illumina quality string as B (quality value 66) or sanger quality string as # (quality value 35)
EOS
}

sub execute {
    my $self = shift;

    unless (-s $self->pair1_fastq_file and -s $self->pair2_fastq_file) {
        $self->error_message('Need both pair1 and pair2 fastq file');
        return;
    }
    my $out_dir = $self->output_dir;
    
    my ($p1_base_name, $p1_base_dir) = fileparse($self->pair1_fastq_file);
    my ($p2_base_name, $p2_base_dir) = fileparse($self->pair2_fastq_file);
    
    $p1_base_name =~ s/\.fastq$// if $p1_base_name =~ /\.fastq$/;
    $p2_base_name =~ s/\.fastq$// if $p2_base_name =~ /\.fastq$/;

    $self->error_message("pair1 and 2 fastq file are not in the same directory") and return
        unless $p1_base_dir eq $p2_base_dir;
    $out_dir ||= $p1_base_dir;

    my $p1_in_fh = Genome::Sys->open_file_for_reading($self->pair1_fastq_file);
    unless ($p1_in_fh) {
        $self->error_message('Failed to open fastq file ' . $self->pair1_fastq_file . ": $!");
        return;
    }
    binmode $p1_in_fh, ":utf8";

    my $p2_in_fh = Genome::Sys->open_file_for_reading($self->pair2_fastq_file);
    unless ($p2_in_fh) {
        $self->error_message('Failed to open fastq file ' . $self->pair2_fastq_file . ": $!");
        return;
    }
    binmode $p2_in_fh, ":utf8";
    
    my $pair_filter_file = $out_dir . '/trimq2.pair_end.filtered.fastq'; 
    my $frag_filter_file = $out_dir . '/trimq2.pair_as_fragment.filtered.fastq';

    unless ($self->pair1_out_file){
        $self->pair1_out_file($out_dir ."/$p1_base_name.trimq2.fastq");
    }
    unless ($self->pair2_out_file){
        $self->pair2_out_file($out_dir ."/$p2_base_name.trimq2.fastq");
    }

    #create filehandles for trimmed output
    my $p1_out_fh = Genome::Sys->open_file_for_writing($self->pair1_out_file);
    unless ($p1_out_fh) {
        $self->error_message('Failed to open output file ' . $self->pair1_out_file . ": $!");
        return;
    }
    binmode $p1_out_fh, ":utf8";

    my $p2_out_fh = Genome::Sys->open_file_for_writing($self->pair2_out_file);
    unless ($p2_out_fh) {
        $self->error_message('Failed to open output file ' . $self->pair2_out_file . ": $!");
        return;
    }
    binmode $p2_out_fh, ":utf8";
    
    unless ($self->pair_as_frag_file){
        $self->pair_as_frag_file($out_dir . '/trimq2.pair_as_fragment.fastq');
    }

    my $frag_fh = Genome::Sys->open_file_for_writing($self->pair_as_frag_file);
    unless ($frag_fh) {
        $self->error_message('Failed to open ' . $self->pair_as_frag_file . ": $!");
        return;
    }
    binmode $frag_fh, ":utf8";

    #create filehandles for dumping filtered reads
    my $pair_filter_fh = Genome::Sys->open_file_for_writing($pair_filter_file);
    unless ($pair_filter_fh) {
        $self->error_message('Failed to open filtered file '. $pair_filter_file . ": $!");
        return;
    }
    binmode $pair_filter_fh, ":utf8";

    my $frag_filter_fh = Genome::Sys->open_file_for_writing($frag_filter_file);
    unless ($frag_filter_fh) {
        $self->error_message('Failed to open filtered file '. $frag_filter_file . ": $!");
        return;
    }
    binmode $frag_filter_fh, ":utf8";

    #create report fh
    unless ($self->report_file){
        $self->report_file( $out_dir . '/trimq2.report' );
    }
    if (-e $self->report_file) {
        $self->warning_message($self->report_file." already exist. Removing it");
        unlink $self->report_file;
    }

    my $report_fh = Genome::Sys->open_file_for_writing($self->report_file);
    unless ($report_fh) {
        $self->error_message("Failed to open report file " . $self->report_file . ": $!");
        return;
    }
    binmode $report_fh, ":utf8";

    #base counts
    my $ori_ct    = 0;
    my $trim_ct   = 0;
    my $filter_ct = 0;

    #read counts
    my $rd_ori_ct       = 0;
    my $rd_trim_ct      = 0;
    my $rd_filter_ct    = 0; 
    my $rd_pair_frag_ct = 0;
    
    my $qual_str = $self->trim_string;
    
    while (my $p1_header = $p1_in_fh->getline) {
        my $seq1  = $p1_in_fh->getline;
        my $sep1  = $p1_in_fh->getline;
        my $qual1 = $p1_in_fh->getline;

        my ($clean_header1) = $p1_header =~ /^@(\S+)\s+/;
        my $pair_basename1;
        
        if ($clean_header1 =~ /^(\S+)[12]$/) {
            $pair_basename1 = $1;
        }
        else {
            $self->error_message("Header $p1_header in ".$self->pair1_fastq_file." does not end with 1 or 2, not pair end data");
            return;
        }
        
        my $p2_header = $p2_in_fh->getline;
        my $seq2  = $p2_in_fh->getline;
        my $sep2  = $p2_in_fh->getline;
        my $qual2 = $p2_in_fh->getline;
    
        my ($clean_header2) = $p2_header =~ /^@(\S+)\s+/;
        my $pair_basename2;
        
        if ($clean_header2 =~ /^(\S+)[12]$/) {
            $pair_basename2 = $1;
        }
        else {
            $self->error_message("Header $p2_header in ".$self->pair2_fastq_file." does not end with 1 or 2, not pair end data");
            return;
        }

        unless ($pair_basename1 eq $pair_basename2) {
            $self->error_message("pair basename are different between pair 1 and 2: $pair_basename1 and $pair_basename2");
            return;
        }
        
        my $seq_length1 = (length $seq1) - 1; #account for new line
        my $seq_length2 = (length $seq2) - 1; #account for new line

        $ori_ct += $seq_length1 + $seq_length2;
        $rd_ori_ct += 2;

        my $trim_qual1;
        my $trim_length1;
        my $trim_seq1;
        my $trimmed_length1;
        my $was_trimmed1;
        if ($qual1 =~ /$qual_str/){
            $was_trimmed1 = 1;
            ($trim_qual1) = $qual1 =~ /^(\S*?)$qual_str/;
            $trim_length1 = length $trim_qual1; 
            $trim_seq1 = substr($seq1, 0, $trim_length1);
            $trimmed_length1 = $seq_length1 - $trim_length1
        }
        my $trim_qual2;
        my $trim_length2;
        my $trim_seq2;
        my $trimmed_length2;
        my $was_trimmed2;
        if ($qual2 =~ /$qual_str/){
            $was_trimmed2 = 1;
            ($trim_qual2) = $qual2 =~ /^(\S*?)$qual_str/;
            $trim_length2 = length $trim_qual2; 
            $trim_seq2 = substr($seq2, 0, $trim_length2);
            $trimmed_length2 = $seq_length2 - $trim_length2
        }

        my ($length1, $length2);
        if ($was_trimmed1){
            $length1 = $trim_length1 
        }else{
            $length1 =  $seq_length1;
        }
        if ($was_trimmed2){
            $length2 = $trim_length2 
        }else{
            $length2 =  $seq_length2;
        }

        if ($length1 >= $self->length_limit){
            if ($length2 >= $self->length_limit){

                if ($was_trimmed1){ #read was trimmed
                    $p1_out_fh->print($p1_header, $trim_seq1."\n", $sep1, $trim_qual1."\n"); 
                    $report_fh->print($clean_header1."\tT\t".$trimmed_length1."\n");  #in report t for trimmed
                    $rd_trim_ct++;
                    $trim_ct += $trimmed_length1;
                }else{ #read wasn't trimmed
                    $p1_out_fh->print($p1_header, $seq1, $sep1, $qual1);
                }

                if ($was_trimmed2){ #read was trimmed
                    $p2_out_fh->print($p2_header, $trim_seq2."\n", $sep2, $trim_qual2."\n"); 
                    $report_fh->print($clean_header2."\tT\t".$trimmed_length2."\n");  #in report t for trimmed
                    $rd_trim_ct++;
                    $trim_ct += $trimmed_length2;
                }else{ #read wasn't trimmed
                    $p2_out_fh->print($p2_header, $seq2, $sep2, $qual2);
                }
            }else{
                #print 1_as_frag, filter 2
                $rd_pair_frag_ct++;
                if ($was_trimmed1){ #read was trimmed
                    $frag_fh->print($p1_header, $trim_seq1."\n", $sep1, $trim_qual1."\n"); 
                    $report_fh->print($clean_header1."\tPF\t".$trimmed_length1."\n");  #In report PF for pair-end-as-fragment
                    $trim_ct += $trimmed_length1;
                    $rd_trim_ct++;
                }else{
                    $frag_fh->print($p1_header, $seq1, $sep1, $qual1);
                    $report_fh->print($clean_header1."\tPF\t50\n");  #In report PF for pair-end-as-fragment
                }

                $frag_filter_fh->print($p2_header, $seq2, $sep2, $qual2);
                $report_fh->print($clean_header2."\tF\t".$seq_length2."\n");  #In report F for filtered
                $filter_ct += $seq_length2;
                $rd_filter_ct++;  
            }
        }else{
            if ($length2 >= $self->length_limit){
                #filter 1, print 2_as_frag

                $frag_filter_fh->print($p1_header, $seq1, $sep1, $qual1);
                $report_fh->print($clean_header1."\tF\t".$seq_length1."\n");  #In report F for filtered
                $filter_ct += $seq_length1;
                $rd_filter_ct++;  
                
                $rd_pair_frag_ct++;
                if ($was_trimmed2){ #read was trimmed
                    $frag_fh->print($p2_header, $trim_seq2."\n", $sep2, $trim_qual2."\n"); 
                    $report_fh->print($clean_header2."\tPF\t".$trimmed_length2."\n");  #in report PF for pair-end-as-fragment
                    $trim_ct += $trimmed_length2;
                    $rd_trim_ct++;
                }else{
                    $frag_fh->print($p2_header, $seq2, $sep2, $qual2);
                    $report_fh->print($clean_header2."\tPF\t50\n");  #in report PF for pair-end-as-fragment
                }
            }else{
                #filterboth
                $pair_filter_fh->print($p1_header, $seq1, $sep1, $qual1);
                $report_fh->print($clean_header1."\tF\t".$seq_length1."\n");  #In report F for filtered
                $filter_ct += $seq_length1;
                $rd_filter_ct++;  

                $pair_filter_fh->print($p2_header, $seq2, $sep2, $qual2);
                $report_fh->print($clean_header2."\tF\t".$seq_length2."\n");  #In report F for filtered
                $filter_ct += $seq_length2;
                $rd_filter_ct++;  
            }
        }
    }

    $frag_fh->close;
    $p1_in_fh->close;
    $p2_in_fh->close;
    $p1_out_fh->close;
    $p2_out_fh->close;
    $pair_filter_fh->close;
    $frag_filter_fh->close;

    $report_fh->print("\nSummary:\n");
    my $rd_new_ct  = $rd_ori_ct - $rd_filter_ct;
    #my $rd_percent = 100*$rd_new_ct/$rd_ori_ct;
    my $new_ct  = $ori_ct - $trim_ct - $filter_ct;
    my $percent = 100*$new_ct/$ori_ct;

    $report_fh->print("\nNumberOfOriginalReads  NumberOfTrimmedReads  NumberOfFilteredReads  NumberOfRemainingReads  NumberOfPairAsFragmentReads\n");
    $report_fh->printf("%21s%22s%23s%24s%12s\n", $rd_ori_ct, $rd_trim_ct, $rd_filter_ct, $rd_new_ct, $rd_pair_frag_ct);
    $report_fh->print("\nNumberOfOriginalBases  NumberOfTrimmedBases  NumberOfFilteredBases  NumberOfResultingBases  Percentage\n");
    $report_fh->printf("%21s%22s%23s%24s%11.1f%%\n", $ori_ct, $trim_ct, $filter_ct, $new_ct, $percent);
    $report_fh->close;

    return 1;
}

1;

