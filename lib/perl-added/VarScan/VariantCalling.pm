package VarScan::VariantCalling;

use warnings;
use strict;
use Getopt::Long;

=head1 NAME

VarScan::VariantCalling - Functions for combining and calling variants

=head1 VERSION

    Version 1.03

=cut

our $VERSION = '1.03';
our $p_value_threshold = 1e-06;
our $seq_error_rate = 0.01;

=head1 SYNOPSIS

    This module contains numerous subroutines for variant calling.

=head1 FUNCTIONS

=cut

## Define shared variables ##

my $verbose, my $outfile, my $positions_file, my $min_coverage, my $min_avg_qual, my $min_reads1, my $min_reads2, my $min_var_freq, my $min_strands2;


################################################################################

=head2	combine_variants - combines overlapping variants

=cut
################################################################################

sub combine_variants
{                               
    (my $infile, my $outfile) = @_;

    if(!(-e $infile))
    {
	die "***ERROR*** $infile not found\n";
    }
    ## Determine type of file ##
    my $input = new FileHandle ($infile) or die "Can't open infile $infile\n!";
    my $first_line = <$input>;
    close($input);

    if($first_line)
    {
	my @lineContents = split(/\t/, $first_line);
	
	if($lineContents[3] && ($lineContents[3] eq "indel_type" || $lineContents[3] eq "INSERTION" || $lineContents[3] eq "DELETION"))
	{
	    combine_indels($infile, $outfile);
	}
	else
	{
	    combine_snps($infile, $outfile);
	}
    }
    else
    {
	system("touch $outfile");
    }

}


################################################################################

=head2	combine_snps - combines overlapping variants

=cut
################################################################################

sub combine_snps
{
    (my $infile, my $outfile) = @_;

    my %SNPcontexts = my %NumReads = my %SupportingStrands = my %VariantAlleles = ();
    
    my %SNPqualitySum = ();
    
    my %CombineStats = ();
    $CombineStats{'snp_events'} = 0;
    $CombineStats{'unique_snps'} = 0;
    
    
    ## Open the input file ##
    
    my $input = new FileHandle ($infile);
    my $lineCounter = 0;

    while (<$input>)
    {
	    chomp;
	    my $line = $_;
	    $lineCounter++;

	    my @lineContents = split(/\t/, $line);

	    if($lineContents[3] && $lineContents[0] ne "chrom" && $lineContents[0] ne "ref_name")
	    {
		    (my $chrom, my $position, my $allele1, my $allele2, my $read_name, my $read_pos, my $align_strand, my $base_qual_score, my $snp_in_context) = split(/\t/, $line);

		    $base_qual_score = 15 if(!$base_qual_score);
	    
		    if($chrom && $position && $align_strand && $allele1 && $allele2 && $allele2 =~ /[ACGT]/)
		    {
			my $snp_key = "$chrom\t$position\t$allele1\t$allele2";		
			    $CombineStats{'snp_events'}++;
			    $NumReads{$snp_key}++;
		    
			    $SNPqualitySum{$snp_key} = 0 if(!$SNPqualitySum{$snp_key});
			    $SNPqualitySum{$snp_key} += $base_qual_score;

			    $VariantAlleles{$snp_key} = "" if(!$VariantAlleles{$snp_key});
			    $SupportingStrands{$snp_key} = "" if(!$SupportingStrands{$snp_key});
			    
			    ## Save the strand information ##
			    
			    if(!$SupportingStrands{$snp_key} || (length($SupportingStrands{$snp_key}) < 2 && substr($SupportingStrands{$snp_key}, 0, 1) ne $align_strand))
			    {
				    $SupportingStrands{$snp_key} .= $align_strand;
			    }
			   
			    ## Save the SNP in context ##
			    
			    if($snp_in_context && (!$SNPcontexts{$snp_key} || length($snp_in_context) >= length($SNPcontexts{$snp_key})))
			    {
				$SNPcontexts{$snp_key} = $snp_in_context;
			    }
			    
			    ## Save the variant allele info ##
			    
			    if(!($VariantAlleles{$snp_key} =~ $allele2))
			    {
				    $VariantAlleles{$snp_key} .= "/" if($VariantAlleles{$snp_key});					
				    $VariantAlleles{$snp_key} .= $allele2 
			    }
		    }

	    }
    }
    
    open(OUTFILE, ">$outfile") or die "Can't open outfile: $!\n";
    print OUTFILE "chrom\tposition\tref\tvar\treads2\tavgQual\tstrands\tsnp_in_context\n";	
    foreach my $snp_key (sort keys %VariantAlleles)
    {
	    $CombineStats{'unique_snps'}++;
	    ## Determine average base quality ##
	    my $avg_base_qual = $SNPqualitySum{$snp_key} / $NumReads{$snp_key};
	    $avg_base_qual = sprintf("%d", $avg_base_qual);
	    my $num_reads = $NumReads{$snp_key};

	    my $num_strands = length($SupportingStrands{$snp_key});            

	    my $snp_in_context = "";
	    $snp_in_context = $SNPcontexts{$snp_key} if($SNPcontexts{$snp_key});

	    print OUTFILE "$snp_key\t$num_reads\t$avg_base_qual\t$num_strands\t$snp_in_context\n";
    }
    close(OUTFILE);


    print "$CombineStats{'snp_events'} substitution events\n";
    print "$CombineStats{'unique_snps'} unique SNPs\n";
    
    return 0;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}




################################################################################################
=head2 # CombineIndels - a subroutine that combines read indel events into indels
################################################################################################
=cut

sub combine_indels
{
        (my $infile, my $outfile) = @_;
#        my $outfile = $infile . ".combined";
	
	if(!($infile && -e $infile))
	{
		return("Indels file $infile not found");
	}

	my $default_qual_score = $VarScan::default_qual_score;

        ## Declare or reset variables ##
	my %IndelContext = my %IndelPrecisionScore = my %NumReads = my %SupportingReads = my %SupportingStrands = ();
	my %IndelQualitySum = my %IndelConfidenceSum = ();
	
	my %CombineStats = ();
	$CombineStats{'indel_events'} = 0;
	$CombineStats{'unique_indels'} = 0;
        
	## Open the input file ##
	
	my $input = new FileHandle ($infile);
	my $lineCounter = 0;

	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		my @lineContents = split(/\t/, $line);
		
		if($lineContents[4] && $lineContents[0] ne "chrom" && $lineContents[0] ne "ref_name")
		{
			my $chrom = $lineContents[0];
			my $chr_start = $lineContents[1];
			my $chr_stop = $lineContents[2];
			my $indel_type = $lineContents[3];
			my $indel_size = $lineContents[4];
			my $read_name = $lineContents[5];
			my $strand = $lineContents[8];
			my $conf_score = $lineContents[9];
			my $indel_context = $lineContents[10];
			my $qual_score = $lineContents[11];

                        ## Fix empty quality scores ##	
                        $qual_score = $default_qual_score if(!$qual_score);
			$conf_score = 0 if(!$conf_score);

			my $indel_key = "$chrom\t$chr_start\t$chr_stop\t$indel_type\t$indel_size";		
		
			$CombineStats{'indel_events'}++;

			$IndelQualitySum{$indel_key} += $qual_score;
			$IndelConfidenceSum{$indel_key} += $conf_score;
			$NumReads{$indel_key}++;

			## Build this supporting read ##
			
			my $supporting_read = "$read_name($conf_score:$strand:$qual_score)";
			$SupportingReads{$indel_key} .= "," if($SupportingReads{$indel_key});
			$SupportingReads{$indel_key} .= $supporting_read;
			
			## Use this event's indel context if it's more precise or more complete
			$IndelContext{$indel_key} = $indel_context if(!$IndelContext{$indel_key} || ($indel_context && $conf_score >= $IndelPrecisionScore{$indel_key} && length($indel_context) > length($IndelContext{$indel_key})));		
		
			## Take the highest precision score ##
			$IndelPrecisionScore{$indel_key} = $conf_score if(!$IndelPrecisionScore{$indel_key} || $conf_score > $IndelPrecisionScore{$indel_key});		
		}
	}
	
	open(OUTFILE, ">$outfile") or die "Can't open outfile: $!\n";
	print OUTFILE "chrom\tchr_start\tchr_stop\tindel_type\tindel_size\tindel_in_context\tnum_reads\tnum_strands\tavg_conf_score\tavg_base_qual\n";	
	foreach my $indel_key (sort keys %NumReads)
	{
		$CombineStats{'unique_indels'}++;
		my $num_reads = $NumReads{$indel_key};
		my $context = $IndelContext{$indel_key};
		my $num_strands = 1;
		if($SupportingReads{$indel_key} =~ '\:\+\:' && $SupportingReads{$indel_key} =~ '\:\-\:')
		{
			$num_strands = 2;
		}

		## Determine average base quality ##
		my $avg_base_qual = $IndelQualitySum{$indel_key} / $num_reads;
		$avg_base_qual = sprintf("%d", $avg_base_qual);

		my $avg_conf_score = $IndelConfidenceSum{$indel_key} / $num_reads;
		$avg_conf_score = sprintf("%d", $avg_conf_score);

#		print OUTFILE "$indel_key\t$context\t$num_reads\t$SupportingReads{$indel_key}\n"; #\t$avg_base_qual\t$num_strands\n";
		print OUTFILE "$indel_key\t$context\t$num_reads\t$num_strands\t$avg_conf_score\t$avg_base_qual\n"; #\t$SupportingReads{$indel_key}\n";	
	}
	close(OUTFILE);

	print "$CombineStats{'indel_events'} indel events\n";
	print "$CombineStats{'unique_indels'} unique indels\n";    

    return(0);
}



################################################################################################
=head2 limit_snps - Limit SNPs to chromosome positions from a file
#
################################################################################################
=cut

sub limit_snps
{
    (my $variants_file, my $positions_file, my $outfile) = @_;

    if(!$positions_file || !(-e $positions_file))
    {
            die "\n***Positions file not found!\n";
    }

    if(!$variants_file || !(-e $variants_file))
    {
            die "\n***Variants file not found!\n";
    }

	my %stats = ();
	$stats{'total_snps'} = $stats{'limit_snps'} = 0;
	my $input, my $lineCounter;

	print "Parsing ROI positions...\n";

	my %roi_positions = ();

	if($positions_file)
	{
		## Parse the alignment blocks ##
	
		$input = new FileHandle ($positions_file);
		$lineCounter = 0;
		
		while (<$input>)
		{
			chomp;
			my $line = $_;
			$lineCounter++;
			
			if($lineCounter >= 1)# && $lineCounter < 50000)
			{
				(my $chrom, my $position) = split(/\t/, $line);
				$roi_positions{$chrom . "\t" . $position} = 1;
			}
		}
		
		close($input);
		
		print "$lineCounter positions loaded\n";
	}

	print "Parsing SNPs...\n";
	## Parse the combined SNPs ##

	## Open the outfile ##
	if($outfile)
	{
		open(OUTFILE, ">$outfile") or die "Can't open outfile: $!\n";
	}
	
#	if($notfile)
#	{
#		open(NOTFILE, ">$notfile") or die "Can't open outfile: $!\n";
#	}

	$input = new FileHandle ($variants_file);
	$lineCounter = 0;
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		if($lineCounter == 1)
		{
			print OUTFILE "$line\n" if($outfile);
#			print NOTFILE "$line\n" if($notfile);
		}
		if($lineCounter > 1)# && $lineCounter < 1000)
		{
			$stats{'total_snps'}++;
#			(my $chrom, my $position, my $allele1, my $allele2, my $context, my $num_reads, my $reads) = split(/\t/, $line);			
			(my $chrom, my $position) = split(/\t/, $line);

			my $included_flag = 0;
			
			if(!$positions_file || $roi_positions{$chrom . "\t" . $position})
			{
				$stats{'limit_snps'}++;
				print OUTFILE "$line\n" if($outfile);
				$included_flag = 1;
			}
			
			if(!$included_flag)
			{
#				print NOTFILE "$line\n" if($notfile);
			}
		}
	}

	close($input);
	
	close(OUTFILE) if($outfile);
#	close(NOTFILE) if($notfile);

	print "$stats{'total_snps'} total SNPs\n";

	## Get percentage ##
	
        if($stats{'total_snps'})
        {
        	$stats{'limit_pct'} = sprintf("%.2f", ($stats{'limit_snps'} / $stats{'total_snps'} * 100)) . '%';
        	print "$stats{'limit_snps'} SNPs ($stats{'limit_pct'}) remained after filter\n";
        }
    
    return(0);
}


################################################################################################
=head2 remove_snps - Limit SNPs to chromosome positions from a file
#
################################################################################################
=cut

sub remove_snps
{
    (my $variants_file, my $positions_file, my $outfile) = @_;

    if(!$positions_file || !(-e $positions_file))
    {
            die "\n***Positions file not found!\n";
    }

    if(!$variants_file || !(-e $variants_file))
    {
            die "\n***Variants file not found!\n";
    }

	my %stats = ();
	$stats{'total_snps'} = $stats{'limit_snps'} = 0;
	my $input, my $lineCounter;

	print "Parsing ROI positions...\n";

	my %roi_positions = ();

	if($positions_file)
	{
		## Parse the alignment blocks ##
	
		$input = new FileHandle ($positions_file);
		$lineCounter = 0;
		
		while (<$input>)
		{
			chomp;
			my $line = $_;
			$lineCounter++;
			
			if($lineCounter >= 1)# && $lineCounter < 50000)
			{
				(my $chrom, my $position) = split(/\t/, $line);
				$roi_positions{$chrom . "\t" . $position} = 1;
			}
		}
		
		close($input);
		
		print "$lineCounter positions loaded\n";
	}

	print "Parsing SNPs...\n";
	## Parse the combined SNPs ##

	## Open the outfile ##
	if($outfile)
	{
		open(OUTFILE, ">$outfile") or die "Can't open outfile: $!\n";
	}
	
#	if($notfile)
#	{
#		open(NOTFILE, ">$notfile") or die "Can't open outfile: $!\n";
#	}

	$input = new FileHandle ($variants_file);
	$lineCounter = 0;
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		if($lineCounter == 1)
		{
			print OUTFILE "$line\n" if($outfile);
#			print NOTFILE "$line\n" if($notfile);
		}
		if($lineCounter > 1)# && $lineCounter < 1000)
		{
			$stats{'total_snps'}++;
#			(my $chrom, my $position, my $allele1, my $allele2, my $context, my $num_reads, my $reads) = split(/\t/, $line);			
			(my $chrom, my $position) = split(/\t/, $line);

			my $included_flag = 0;
			
			if(!$positions_file || !$roi_positions{$chrom . "\t" . $position})
			{
				$stats{'limit_snps'}++;
				print OUTFILE "$line\n" if($outfile);
				$included_flag = 1;
			}
			
			if(!$included_flag)
			{
#				print NOTFILE "$line\n" if($notfile);
			}
		}
	}

	close($input);
	
	close(OUTFILE) if($outfile);

	print "$stats{'total_snps'} total SNPs\n";

	## Get percentage ##
	
        if($stats{'total_snps'})
        {
        	$stats{'limit_pct'} = sprintf("%.2f", ($stats{'limit_snps'} / $stats{'total_snps'} * 100)) . '%';
        	print "$stats{'limit_snps'} SNPs ($stats{'limit_pct'}) remained after removal filter\n";
        }
    
    return(0);
}




################################################################################################
=head2 filter_variants - Filter SNPs based on numerous options
#
################################################################################################
=cut

sub filter_variants
{
    (my $variants_file) = @_;

    if(!$variants_file || !(-e $variants_file))
    {
            die "\n***Variants file not found!\n";
    }

    ## Get variable values ##
    
    $outfile = $VarScan::output_file;
    $min_coverage = $VarScan::min_coverage;
    $min_avg_qual = $VarScan::min_avg_qual;
    $min_reads1 = $VarScan::min_reads1;
    $min_reads2 = $VarScan::min_reads2;
    $min_var_freq = $VarScan::min_var_freq;
    $min_strands2 = $VarScan::min_strands2;
    $verbose = $VarScan::verbose;

    $min_var_freq = $min_var_freq / 100 if($min_var_freq > 1);

    my %stats = ();
    $stats{'total_snps'} = $stats{'limit_snps'} = 0;
    my $input, my $lineCounter;

    my %roi_positions = ();

    if($positions_file && -e $positions_file)
    {
	print "Parsing ROI positions...\n";

	## Parse the alignment blocks ##

	$input = new FileHandle ($positions_file);
	$lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		if($lineCounter >= 1)# && $lineCounter < 50000)
		{
			(my $chrom, my $position) = split(/\t/, $line);
			$roi_positions{$chrom . "\t" . $position} = 1;
		}
	}
	
	close($input);
	
	print "$lineCounter positions loaded\n";
    }

    print "Parsing SNPs...\n";
    ## Parse the combined SNPs ##

    ## Open the outfile ##
    if($outfile)
    {
	    open(OUTFILE, ">$outfile") or die "Can't open outfile: $!\n";
    }
    
    $input = new FileHandle ($variants_file);
    $lineCounter = 0;
    while (<$input>)
    {
	    chomp;
	    my $line = $_;
	    $lineCounter++;

	    my @lineContents = split(/\t/, $line);

	    my $chrom, my $position, my $coverage, my $reads1, my $reads2, my $avg_qual, my $strands2, my $var_freq;

	    if($lineContents[0] eq "chrom" || $lineContents[0] eq "ref_name")
	    {
		print OUTFILE "$line\n" if($outfile);
	    }
	    else
	    {
		if($lineContents[3] && $lineContents[3] && ($lineContents[3] eq "INSERTION" || $lineContents[3] eq "DELETION")) ## INDELS ##
		{
		    $position = $lineContents[1];
		    $coverage = $lineContents[5];
		    $reads1 = $lineContents[6];
		    $reads2 = $lineContents[7];
		    $avg_qual = $lineContents[10];
		    $strands2 = $lineContents[12];
		    $stats{'total_indels'}++;
		}
		elsif($lineContents[0] && $lineContents[1] && $lineContents[2] && $lineContents[3]) ## SNPS ##
		{
#		    print "$line\n";
		    $position = $lineContents[1];
		    $coverage = $lineContents[4];
		    $reads1 = $lineContents[5];
		    $reads2 = $lineContents[6];
		    $avg_qual = $lineContents[9];
		    $strands2 = $lineContents[11];
		    $strands2 = 1 if(!$strands2);
    
		    $stats{'total_snps'}++;
		}
    
		## Calculate var freq ##
		
		if(!$positions_file || $roi_positions{$chrom . "\t" . $position})
		{
			if(!$min_coverage || $coverage >= $min_coverage)
			{
			    if(!$min_reads2 || $reads2 >= $min_reads2)
			    {
				if(!$min_strands2 || $strands2 >= $min_strands2)
				{
				    if(!$min_avg_qual || $avg_qual >= $min_avg_qual)
				    {
					## Calculate variant frequency ##
					
					$var_freq = $reads2 / ($reads1 + $reads2) if($reads1 || $reads2);
					$var_freq = sprintf("%.4f", $var_freq);
					
					if(!$min_var_freq || $var_freq >= $min_var_freq)
					{
					    print OUTFILE "$line\n" if($outfile);
					    $stats{'limit_snps'}++;
					} # variant freq
				    }
				} # strands2
			    } # reads2
			} # coverage
		} # positions
	    }
		    

		    


    }

    close($input);
    
    close(OUTFILE) if($outfile);

    if($stats{'total_snps'})
    {
	print "$stats{'total_snps'} total SNPs\n";
	$stats{'limit_pct'} = sprintf("%.2f", ($stats{'limit_snps'} / $stats{'total_snps'} * 100)) . '%';
	print "$stats{'limit_snps'} SNPs ($stats{'limit_pct'}) remained after filter\n";
    }
    elsif($stats{'total_indels'})
    {
	print "$stats{'total_indels'} total indels\n";
	$stats{'limit_pct'} = sprintf("%.2f", ($stats{'limit_snps'} / $stats{'total_indels'} * 100)) . '%';
	print "$stats{'limit_snps'} indels ($stats{'limit_pct'}) remained after filter\n";

    }
    
    return(0);
}




################################################################################################
=head2 heuristic_indel_filter - Filter indels for homopolymer-associated artifacts

=cut
################################################################################################

sub heuristic_indel_filter
{
	my $indel = shift(@_);
	
	my $hf_min_quality_score = 15;
	my $hf_min_block_size = 11;
	my $hf_min_hp = 4;
	## Parse out the indel information ##
#	"$ref_name\t$ref_indel_start\t$ref_indel_stop\t$indel_type\t$indel_size\t$read_name\t$read_indel_start\t$read_indel_stop\t$align_strand\t$conf_score\t$indel_in_context\t$indel_qual";
	my @indelContents = split(/\t/, $indel);
	my $indel_type = $indelContents[3];
	my $indel_size = $indelContents[4];
	my $precision_score = $indelContents[9];
	my $indel_in_context = $indelContents[10];
	my $quality_score = $indelContents[11];

	
	## Split the indel-in-context ##
	
	(my $flank5, my $allele1, my $allele2, my $flank3) = split(/[\[\/\]]/, $indel_in_context);
	
	
	## HEURISTIC FILTER 1: SINGLE-BASE INDELS WITH LOW PRECISION SCORES ##	
	

	## HEURISTIC FILTER 2: INDELS WITH LOW BASE QUALITY SCORES ##		
	
	if($quality_score && $quality_score < $hf_min_quality_score)
	{
		return('low_quality_score');
	}
	
	
	## HEURISTIC FILTER 3: SINGLE-BASE INDELS WITH LOW PRECISION SCORES ##	
	
	if($indel_size == 1 && $precision_score < 70)
	{
		return('single_base_low_precision');
	}
	
	
	## HEURISTIC FILTER 4: SINGLE-BASE HOMOPOLYMERS ##	Example: TTAATGGAAA[-/A]CTCTTCCCTG
	
	if($indel_size == 1 && length($allele1) == 1 && length($allele2) == 1)
	{
		my $check_size = $hf_min_hp - 1;
		$check_size = 1 if($check_size < 1);
		my $check5 = substr($flank5, length($flank5) - $check_size, $check_size);
		my $check3 = substr($flank3, 0, $check_size);
		
		if($check5 eq ($allele1 x $check_size) || $check3 eq ($allele1 x $check_size) || $check5 eq ($allele2 x $check_size) || $check3 eq ($allele2 x $check_size))
		{
			return('single_base_homopolymer');		
		}
	}


	## HEURISTIC FILTER 4A: A/T HOMOPOLYMERS ##	Example: AA[-/A]

	if($indel_size <= 5)
	{
		my $check_allele;
		if($allele1 =~ '-')
		{
			$check_allele = $allele2;
		}
		elsif($allele2 =~ '-')
		{
			$check_allele = $allele1;
		}

		if($check_allele eq ('A' x length($check_allele)))
		{
		    if(substr($flank5, length($flank5) - 1, 1) eq 'A' || substr($flank3, 0, 1) eq 'A')
		    {
			return('A_homopolymer');			
		    }
		}
		if($check_allele eq ('T' x length($check_allele)))
		{
		    if(substr($flank5, length($flank5) - 1, 1) eq 'T' || substr($flank3, 0, 1) eq 'T')
		    {
			return('T_homopolymer');			
		    }
		}
	}


	## HEURISTIC FILTER 5: TWO-BASE HOMOPOLYMERS ##		Example: GTGAAATTTT[TT/--]GCTGTCTCAA
	
	if($indel_size == 2 && $precision_score == 100)
	{
		my $check_size = $hf_min_hp - 2;
		$check_size = 1 if($check_size < 1);
		my $check5 = substr($flank5, length($flank5) - $check_size, $check_size);
		my $check3 = substr($flank3, 0, $check_size);
		my $check_base;
		if($allele1 =~ '-' && substr($allele2, 0, 1) eq substr($allele2, 1, 1))
		{
			$check_base = substr($allele2, 0, 1);
		}
		elsif($allele2 =~ '-' && substr($allele1, 0, 1) eq substr($allele1, 1, 1))
		{
			$check_base = substr($allele1, 0, 1);
		}
		
		if($check_base)
		{
			if($check5 eq ($check_base x $check_size) || $check3 eq ($check_base x $check_size))
			{
				return('two_base_homopolymer');
			}
		}
	}


	## HEURISTIC FILTER 6: TWO-BASE SPLIT HOMOPOLYMERS ##	Example: TCC[--/CA]AAC
	
	if($indel_size == 2 && $precision_score == 100)
	{
		my $check5 = substr($flank5, length($flank5) - 2, 2);
		my $check3 = substr($flank3, 0, 2);				
		my $check_base1, my $check_base2;
		if($allele1 =~ '-' && substr($allele2, 0, 1) ne substr($allele2, 1, 1))
		{
			$check_base1 = substr($allele2, 0, 1);
			$check_base2 = substr($allele2, 1, 1);		
		}
		elsif($allele2 =~ '-' && substr($allele1, 0, 1) ne substr($allele1, 1, 1))
		{
			$check_base1 = substr($allele1, 0, 1);
			$check_base2 = substr($allele1, 1, 1);		
		}
		
		if($check_base1 && $check_base2 && $check5 eq ($check_base1 x 2) && $check3 eq ($check_base2 x 2))
		{
			return('two_base_split_homopolymer');
		}
		if($check_base1 && $check_base2 && $check5 eq ($check_base1 x 2) && $check3 eq ($check_base2 x 1))
		{
			return('two_base_split_homopolymer');
		}
		if($check_base1 && $check_base2 && $check5 eq ($check_base1 x 1) && $check3 eq ($check_base2 x 2))
		{
			return('two_base_split_homopolymer');
		}
	}


	return(0);
}






################################################################################################
=head2 combine_readcounts - Combine the results from multiple read count files
#
################################################################################################
=cut

sub combine_readcounts
{
    (my $readcounts_files, $outfile) = @_;

    ## Get variable values ##
    
    $min_coverage = $VarScan::min_coverage;
    $min_avg_qual = $VarScan::min_avg_qual;
    $min_reads1 = $VarScan::min_reads1;
    $min_reads2 = $VarScan::min_reads2;
    $min_var_freq = $VarScan::min_var_freq;
    $min_strands2 = $VarScan::min_strands2;
    $verbose = $VarScan::verbose;

    $min_var_freq = $min_var_freq / 100 if($min_var_freq > 1);

    my %stats = ();
    $stats{'total_snps'} = $stats{'limit_snps'} = 0;
    my $input, my $lineCounter;  
    my $header_line;

    ## Set up hashes to hold the relevant info ##
    
    my %read_coverage = my %reads1 = my %reads2 = my %readsN = my %strands1 = my %strands2 = my %quals1 = my %quals2 = ();
    my %num_libraries = my %variant_context = ();

    my @readcounts_files = split(/\,/, $readcounts_files);
    
    foreach my $variants_file (@readcounts_files)
    {
	if(!$variants_file || !(-e $variants_file))
	{
		die "\n***Variants file not found!\n";
	}

	print "Parsing file $variants_file...\n";

	## Parse the variant file ##

	$input = new FileHandle ($variants_file);
	$lineCounter = 0;
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
    
		my @lineContents = split(/\t/, $line);
    
		if($lineContents[0] eq "chrom" || $lineContents[0] eq "ref_name")
		{
		    $header_line = $line;
		}
		else
		{
		    if($lineContents[6] && $lineContents[3] && ($lineContents[3] eq "INSERTION" || $lineContents[3] eq "DELETION")) ## INDELS ##
		    {
			(my $chrom, my $chr_start, my $chr_stop, my $indel_type, my $indel_size, my $read_coverage, my $reads1, my $reads2, my $readsN, my $avg_qual1, my $avg_qual2, my $strands1, my $strands2, my $conf_score, my $context) = split(/\t/, $line);
			$strands1 = 1 if(!$strands1);
			$strands2 = 1 if(!$strands2);
			my $indel_key = "$chrom\t$chr_start\t$chr_stop\t$indel_type\t$indel_size";
			$read_coverage{$indel_key} += $read_coverage;
			$reads1{$indel_key} += $reads1;
			$reads2{$indel_key} += $reads2;
			$readsN{$indel_key} += $readsN;
			$quals1{$indel_key} .= " " if($quals1{$indel_key});
			$quals2{$indel_key} .= " " if($quals2{$indel_key});
			$quals1{$indel_key} .= $avg_qual1;
			$quals2{$indel_key} .= $avg_qual2;
			$strands1{$indel_key} = $strands1 if(!$strands1{$indel_key} || $strands1 > $strands1{$indel_key});
			$strands2{$indel_key} = $strands2 if(!$strands2{$indel_key} || $strands2 > $strands2{$indel_key});
			$num_libraries{$indel_key}++;
			$variant_context{$indel_key} = $context;
		    }
		    elsif($lineContents[0] && $lineContents[1] && $lineContents[2] && $lineContents[3]) ## SNPS ##
		    {
			(my $chrom, my $position, my $allele1, my $allele2, my $read_coverage, my $reads1, my $reads2, my $readsN, my $avg_qual1, my $avg_qual2, my $strands1, my $strands2) = split(/\t/, $line);
			$strands1 = 1 if(!$strands1);
			$strands2 = 1 if(!$strands2);
			my $snp_key = "$chrom\t$position\t$allele1\t$allele2";	
			$read_coverage{$snp_key} += $read_coverage;
			$reads1{$snp_key} += $reads1;
			$reads2{$snp_key} += $reads2;
			$readsN{$snp_key} += $readsN;
			$strands1{$snp_key} = $strands1 if(!$strands1{$snp_key} || $strands1 > $strands1{$snp_key});
			$strands2{$snp_key} = $strands2 if(!$strands2{$snp_key} || $strands2 > $strands2{$snp_key});
			$quals1{$snp_key} .= " " if($quals1{$snp_key});
			$quals2{$snp_key} .= " " if($quals2{$snp_key});
			$quals1{$snp_key} .= $avg_qual1;
			$quals2{$snp_key} .= $avg_qual2;
		    }
		}
	}
    
	close($input);	
    } # Go to next file 


    ## Open the outfile ##
    if($outfile)
    {
	    open(OUTFILE, ">$outfile") or die "Can't open outfile: $!\n";
	    print OUTFILE "$header_line";
	    print OUTFILE "\tnum_libraries" if(%num_libraries);
	    print OUTFILE "\tvariant_context" if(%variant_context);
	    print OUTFILE "\n";
    }

    foreach my $variant_key (sort byChrPos keys %read_coverage)
    {

	print OUTFILE "$variant_key\t";
	print OUTFILE $read_coverage{$variant_key} . "\t";
	print OUTFILE $reads1{$variant_key} . "\t";
	print OUTFILE $reads2{$variant_key} . "\t";
	print OUTFILE $readsN{$variant_key} . "\t";

	## Calculate the quals ##
	
	my $avg_qual1 = my $avg_qual2 = 0;

	my @quals = split(/\s+/, $quals1{$variant_key});
	my $sumq = my $numq = 0;
	foreach my $qual (@quals)
	{
	    $sumq += $qual;
	    $numq++;
	}
	$avg_qual1 = sprintf("%d", ($sumq / $numq)) if($numq);

	@quals = split(/\s+/, $quals2{$variant_key});
	$sumq = $numq = 0;
	foreach my $qual (@quals)
	{
	    $sumq += $qual;
	    $numq++;
	}
	$avg_qual2 = sprintf("%d", ($sumq / $numq)) if($numq);

	print OUTFILE "$avg_qual1\t";
	print OUTFILE "$avg_qual2\t";

	$strands1{$variant_key} = 0 if(!$strands1{$variant_key});
	$strands2{$variant_key} = 0 if(!$strands2{$variant_key});
	print OUTFILE $strands1{$variant_key} . "\t";
	print OUTFILE $strands2{$variant_key} . "\t";

	$num_libraries{$variant_key} = 0 if(!$num_libraries{$variant_key});
	print OUTFILE $num_libraries{$variant_key} . "\t" if(%num_libraries);

	print OUTFILE $variant_context{$variant_key} . "\t" if($variant_context{$variant_key});

	print OUTFILE "\n";
    }
    
    close(OUTFILE) if($outfile);

    
    return(0);
}





################################################################################################
=head2 compare_variants - compare variants and read counts between normal and tumor
#
################################################################################################
=cut

sub compare_variants
{
    (my $normal_readcounts, my $tumor_readcounts, $outfile) = @_;

    ## Get variable values ##
    
    $min_coverage = $VarScan::min_coverage;
    $min_avg_qual = $VarScan::min_avg_qual;
    $min_reads1 = $VarScan::min_reads1;
    $min_reads2 = $VarScan::min_reads2;
    $min_var_freq = $VarScan::min_var_freq;
    $min_strands2 = $VarScan::min_strands2;
    $verbose = $VarScan::verbose;

    $min_var_freq = $min_var_freq / 100 if($min_var_freq > 1);

    my %stats = ();
    $stats{'total_snps'} = $stats{'shared_snps'} = $stats{'normal_only_snps'} = $stats{'tumor_only_snps'} = 0;

    ## Parse the variants ##
    
    print "Loading Normal read counts...\n";
    my %normal_read_counts = load_read_counts($normal_readcounts);
    
    print "Loading Tumor read counts...\n";
    my %tumor_read_counts = load_read_counts($tumor_readcounts);


    ## Combine and sort variant keys ##
    
    print "Combining and sorting variant keys...\n";
    
    my %variant_keys = ();
    
    foreach my $normal_key (keys %normal_read_counts)
    {
	$variant_keys{$normal_key} = 1;
	if($tumor_read_counts{$normal_key})
	{
	    $stats{'shared_snps'}++;
	}
	else
	{
	    $stats{'normal_only_snps'}++;
	}
    }

    foreach my $tumor_key (keys %tumor_read_counts)
    {
	$variant_keys{$tumor_key} = 1;
	$stats{'tumor_only_snps'}++ if(!$normal_read_counts{$tumor_key});
    }

    $stats{'total_snps'} = $stats{'shared_snps'} + $stats{'normal_only_snps'} + $stats{'tumor_only_snps'};


    ## Open the outfile ##
    if($outfile)
    {
	    open(OUTSHARED, ">$outfile.Shared") or die "Can't open outfile: $!\n";
	    open(OUTNORMAL, ">$outfile.Normal-only") or die "Can't open outfile: $!\n";
	    open(OUTTUMOR, ">$outfile.Tumor-only") or die "Can't open outfile: $!\n";	    

	    print OUTSHARED "chrom\tposition\tref\tvar\treads2\tavgQual\tstrands\n";
	    print OUTNORMAL "chrom\tposition\tref\tvar\treads2\tavgQual\tstrands\n";
	    print OUTTUMOR "chrom\tposition\tref\tvar\treads2\tavgQual\tstrands\n";
    }
    
    foreach my $variant_key (sort byChrPos keys %variant_keys)
    {
	if($normal_read_counts{$variant_key} && $tumor_read_counts{$variant_key})
	{
	    ## Probably a germline SNP, so set aside ##
	    
	    print OUTSHARED "$variant_key\t";
	    print OUTSHARED "Shared\t";
	    print OUTSHARED $normal_read_counts{$variant_key} . "\t";
	    print OUTSHARED $tumor_read_counts{$variant_key} . "\n";
	}
	elsif($normal_read_counts{$variant_key} && !$tumor_read_counts{$variant_key})
	{
	    ## Probably a germline SNP, so set aside ##
	    
	    print OUTNORMAL "$variant_key\t0\t0\t0\n";    
	}
	if(!$normal_read_counts{$variant_key} && $tumor_read_counts{$variant_key})
	{
	    ## Probably a germline SNP, so set aside ##
	    
	    print OUTTUMOR "$variant_key\t0\t0\t0\n";
	}
    }
    
    close(OUTSHARED) if($outfile);
    close(OUTNORMAL) if($outfile);
    close(OUTTUMOR) if($outfile);

    
    print $stats{'total_snps'} . " total variants\n";
    print $stats{'normal_only_snps'} . " detected in Normal-only\n";
    print $stats{'tumor_only_snps'} . " detected in Tumor-only\n";
    print $stats{'shared_snps'} . " Shared by both Normal and Tumor\n";

    
    
    return(0);
}



################################################################################################
=head2 load_read_counts - Load a read counts file
#
################################################################################################
=cut

sub load_read_counts
{
	my $variants_file = shift(@_);
	my %read_counts = ();
        
	## Parse the variant file ##

	my $input = new FileHandle ($variants_file);
	my $lineCounter = 0;
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
    
		my @lineContents = split(/\t/, $line);
    
		if($lineContents[0] eq "chrom" || $lineContents[0] eq "ref_name")
		{

		}
		else
		{
		    if($lineContents[6] && $lineContents[3] && ($lineContents[3] eq "INSERTION" || $lineContents[3] eq "DELETION")) ## INDELS ##
		    {
			(my $chrom, my $chr_start, my $chr_stop, my $indel_type, my $indel_size, my $read_coverage, my $reads1, my $reads2, my $readsN, my $avg_qual1, my $avg_qual2, my $strands1, my $strands2, my $conf_score, my $context) = split(/\t/, $line);
			$strands1 = 1 if(!$strands1);
			$strands2 = 1 if(!$strands2);
			my $indel_key = "$chrom\t$chr_start\t$chr_stop\t$indel_type\t$indel_size";
			$read_counts{$indel_key} = "$read_coverage\t$reads1\t$reads2";
		    }
		    elsif($lineContents[0] && $lineContents[1] && $lineContents[2] && $lineContents[3]) ## SNPS ##
		    {
			(my $chrom, my $position, my $allele1, my $allele2, my $read_coverage, my $reads1, my $reads2, my $readsN, my $avg_qual1, my $avg_qual2, my $strands1, my $strands2) = split(/\t/, $line);
			$strands1 = 1 if(!$strands1);
			$strands2 = 1 if(!$strands2);
			my $snp_key = "$chrom\t$position\t$allele1\t$allele2";
			$read_counts{$snp_key} = "$read_coverage\t$reads1\t$reads2";
		    }
		}
	}
    
	close($input);
	
	return(%read_counts);
}




################################################################################################
=head2 call_somatic - Call SNPs as germline or somatic based on normal/tumor values
#
################################################################################################
=cut

sub call_somatic
{
    (my $normal_readcounts, my $tumor_readcounts, my $output_file) = @_;

    ## Set up hashes to hold the relevant info ##
    my %stats = ();
    my %normal_position_coverage = ();
    my %normal_read_coverage = my %normal_reads1 = my %normal_reads2 = my %normal_readsN = my %normal_strands1 = my %normal_strands2 = my %normal_quals1 = my %normal_quals2 = ();
    my %tumor_read_coverage = my %tumor_reads1 = my %tumor_reads2 = my %tumor_readsN = my %tumor_strands1 = my %tumor_strands2 = my %tumor_quals1 = my %tumor_quals2 = ();

    ## Get the NORMAL read counts ##

    print "Parsing Normal SNPs from $normal_readcounts...\n";

    my $input = new FileHandle ($normal_readcounts);
    my $lineCounter = 0;
    while (<$input>)
    {
	    chomp;
	    my $line = $_;
	    $lineCounter++;

	    my @lineContents = split(/\t/, $line);

	    if($lineContents[0] eq "chrom" || $lineContents[0] eq "ref_name")
	    {
#		$header_line = $line;
	    }
	    else
	    {
		if($lineContents[6] && $lineContents[3] && ($lineContents[3] eq "INSERTION" || $lineContents[3] eq "DELETION")) ## INDELS ##
		{
		    (my $chrom, my $chr_start, my $chr_stop, my $indel_type, my $indel_size, my $read_coverage, my $reads1, my $reads2, my $readsN, my $avg_qual1, my $avg_qual2, my $strands1, my $strands2, my $conf_score, my $context) = split(/\t/, $line);
		    $strands1 = 1 if(!$strands1);
		    $strands2 = 1 if(!$strands2);
		    my $indel_key = "$chrom\t$chr_start\t$chr_stop\t$indel_type\t$indel_size";
		    $normal_read_coverage{$indel_key} += $read_coverage;
		    $normal_reads1{$indel_key} += $reads1;
		    $normal_reads2{$indel_key} += $reads2;
		    $normal_readsN{$indel_key} += $readsN;
		    $normal_quals1{$indel_key} .= " " if($normal_quals1{$indel_key});
		    $normal_quals2{$indel_key} .= " " if($normal_quals2{$indel_key});
		    $normal_quals1{$indel_key} .= $avg_qual1;
		    $normal_quals2{$indel_key} .= $avg_qual2;
		    $normal_strands1{$indel_key} = $strands1 if(!$normal_strands1{$indel_key} || $strands1 > $normal_strands1{$indel_key});
		    $normal_strands2{$indel_key} = $strands2 if(!$normal_strands2{$indel_key} || $strands2 > $normal_strands2{$indel_key});
		}
		elsif($lineContents[0] && $lineContents[1] && $lineContents[2] && $lineContents[3]) ## SNPS ##
		{
		    (my $chrom, my $position, my $allele1, my $allele2, my $read_coverage, my $reads1, my $reads2, my $readsN, my $avg_qual1, my $avg_qual2, my $strands1, my $strands2) = split(/\t/, $line);
		    $strands1 = 1 if(!$strands1);
		    $strands2 = 1 if(!$strands2);
		    my $snp_key = "$chrom\t$position\t$allele1\t$allele2";	
		    $normal_read_coverage{$snp_key} += $read_coverage;
#		    $normal_position_coverage{$chrom . "\t" . $position} = $read_coverage if(!$normal_position_coverage{$chrom . "\t" . $position} || $read_coverage > $normal_position_coverage{$chrom . "\t" . $position});
		    $normal_position_coverage{$chrom . "\t" . $position} = $reads1 if(!$normal_position_coverage{$chrom . "\t" . $position} || $reads1 > $normal_position_coverage{$chrom . "\t" . $position});
		    $normal_reads1{$snp_key} += $reads1;
		    $normal_reads2{$snp_key} += $reads2;
		    $normal_readsN{$snp_key} += $readsN;
		    $normal_strands1{$snp_key} = $strands1 if(!$normal_strands1{$snp_key} || $strands1 > $normal_strands1{$snp_key});
		    $normal_strands2{$snp_key} = $strands2 if(!$normal_strands2{$snp_key} || $strands2 > $normal_strands2{$snp_key});
		    $normal_quals1{$snp_key} .= " " if($normal_quals1{$snp_key});
		    $normal_quals2{$snp_key} .= " " if($normal_quals2{$snp_key});
		    $normal_quals1{$snp_key} .= $avg_qual1;
		    $normal_quals2{$snp_key} .= $avg_qual2;
		    $stats{'snps_normal'}++;
		}
		else
		{
		    warn "Warning: Variant type not recognized: $line\n";
		}
	    }
    }

    close($input);	


    print $stats{'snps_normal'} . " SNPs called in Normal ($lineCounter lines)\n";

    print "Parsing Tumor SNPs from $tumor_readcounts...\n";
    
    ## Get the TUMOR read counts ##

    $input = new FileHandle ($tumor_readcounts);
    $lineCounter = 0;
    while (<$input>)
    {
	    chomp;
	    my $line = $_;
	    $lineCounter++;

	    my @lineContents = split(/\t/, $line);

	    if($lineContents[0] eq "chrom" || $lineContents[0] eq "ref_name")
	    {
#		$header_line = $line;
	    }
	    else
	    {
		if($lineContents[6] && $lineContents[3] && ($lineContents[3] eq "INSERTION" || $lineContents[3] eq "DELETION")) ## INDELS ##
		{
		    (my $chrom, my $chr_start, my $chr_stop, my $indel_type, my $indel_size, my $read_coverage, my $reads1, my $reads2, my $readsN, my $avg_qual1, my $avg_qual2, my $strands1, my $strands2, my $conf_score, my $context) = split(/\t/, $line);
		    $strands1 = 1 if(!$strands1);
		    $strands2 = 1 if(!$strands2);
		    my $indel_key = "$chrom\t$chr_start\t$chr_stop\t$indel_type\t$indel_size";
		    $tumor_read_coverage{$indel_key} += $read_coverage;
		    $tumor_reads1{$indel_key} += $reads1;
		    $tumor_reads2{$indel_key} += $reads2;
		    $tumor_readsN{$indel_key} += $readsN;
		    $tumor_quals1{$indel_key} .= " " if($tumor_quals1{$indel_key});
		    $tumor_quals2{$indel_key} .= " " if($tumor_quals2{$indel_key});
		    $tumor_quals1{$indel_key} .= $avg_qual1;
		    $tumor_quals2{$indel_key} .= $avg_qual2;
		    $tumor_strands1{$indel_key} = $strands1 if(!$tumor_strands1{$indel_key} || $strands1 > $tumor_strands1{$indel_key});
		    $tumor_strands2{$indel_key} = $strands2 if(!$tumor_strands2{$indel_key} || $strands2 > $tumor_strands2{$indel_key});
		}
		elsif($lineContents[0] && $lineContents[1] && $lineContents[2] && $lineContents[3]) ## SNPS ##
		{
		    (my $chrom, my $position, my $allele1, my $allele2, my $read_coverage, my $reads1, my $reads2, my $readsN, my $avg_qual1, my $avg_qual2, my $strands1, my $strands2) = split(/\t/, $line);
		    $strands1 = 1 if(!$strands1);
		    $strands2 = 1 if(!$strands2);
		    my $snp_key = "$chrom\t$position\t$allele1\t$allele2";	
		    $tumor_read_coverage{$snp_key} += $read_coverage;
		    $tumor_reads1{$snp_key} += $reads1;
		    $tumor_reads2{$snp_key} += $reads2;
		    $tumor_readsN{$snp_key} += $readsN;
		    $tumor_strands1{$snp_key} = $strands1 if(!$tumor_strands1{$snp_key} || $strands1 > $tumor_strands1{$snp_key});
		    $tumor_strands2{$snp_key} = $strands2 if(!$tumor_strands2{$snp_key} || $strands2 > $tumor_strands2{$snp_key});
		    $tumor_quals1{$snp_key} .= " " if($tumor_quals1{$snp_key});
		    $tumor_quals2{$snp_key} .= " " if($tumor_quals2{$snp_key});
		    $tumor_quals1{$snp_key} .= $avg_qual1;
		    $tumor_quals2{$snp_key} .= $avg_qual2;
		    $stats{'snps_tumor'}++;
		}
		else
		{
		    warn "Warning: Variant type not recognized: $line\n";
		}
	    }
    }

    close($input);

    print $stats{'snps_tumor'} . " SNPs called in Tumor ($lineCounter lines)\n";

    print "Comparing SNPs between Normal and Tumor...\n";

    if($output_file)
    {
	open(OUTFILE, ">$output_file") or die "can't open outfile: $!\n";
	print OUTFILE "chrom\tposition\tref\tvar\tnormal_reads1\tnormal_reads2\tnormal_freq\ttumor_reads1\ttumor_reads2\ttumor_freq\tsomatic_status\tp_value\n";
    }

    ## Go through the SNPs with Normal and Tumor coverage ##
    
    foreach my $snp_key (sort keys %tumor_read_coverage)
    {
	(my $chrom, my $position, my $allele1, my $allele2) = split(/\t/, $snp_key);
	
	## If we have coverage, but no SNP call, at this position, make the correction ##
	
#	if(!$normal_read_coverage{$snp_key} && $normal_position_coverage{$chrom . "\t" . $position})
#	{
#	    $normal_reads1{$snp_key} = $normal_read_coverage{$snp_key} = $normal_position_coverage{$chrom . "\t" . $position};
#	    $normal_reads2{$snp_key} = 0;
#	}
	
	if($normal_read_coverage{$snp_key})
	{
	    $stats{'snps_covered_both'}++;


	    ## Determine if coverage is sufficient for variant calling ##

	    my $normal_coverage = my $tumor_coverage = 0;
	    $normal_coverage = $normal_reads1{$snp_key} + $normal_reads2{$snp_key};
	    $normal_coverage = 0 if(!$normal_coverage);
	    $tumor_coverage = $tumor_reads1{$snp_key} + $tumor_reads2{$snp_key};
	    $tumor_coverage = 0 if(!$tumor_coverage);

	    if($normal_coverage >= 10 && $tumor_coverage >= 10)
	    {
		$stats{'snps_covered_both_10x'}++;
		
		## Determine variant frequency ##
		my $normal_freq = my $tumor_freq = 0;
		$normal_freq = $normal_reads2{$snp_key} / $normal_coverage if($normal_coverage);
		$tumor_freq = $tumor_reads2{$snp_key} / $tumor_coverage if($tumor_coverage);
		
		## Set default genotype to reference until proven otherwise ##	    
		
		my $reference_gt = my $normal_gt = my $tumor_gt = $allele1 . $allele1;
	
    
		## Determine normal genotype ##
    
		if($normal_reads2{$snp_key} > 0)
		{
		    ## Calculate a p-value ##
			
		    my $normal_p_value = VarScan::FisherTest::calculate_p_value($normal_coverage, 0, $normal_reads1{$snp_key}, $normal_reads2{$snp_key}, 0);
		    $normal_p_value = VarScan::FisherTest::format_p_value($normal_p_value);
	
		    ## Determine the normal genotype ##
		    
		    if($normal_freq >= 0.25 || $normal_p_value < 1e-06)
		    {
			$normal_gt = $allele1 . $allele2;
			$normal_gt = $allele2 . $allele2 if($normal_freq > 0.90);
		    }
		}
    
		if($tumor_reads2{$snp_key} > 0)
		{
		    ## Calculate a p-value ##
			
		    my $tumor_p_value = VarScan::FisherTest::calculate_p_value($tumor_coverage, 0, $tumor_reads1{$snp_key}, $tumor_reads2{$snp_key}, 0);
		    $tumor_p_value = VarScan::FisherTest::format_p_value($tumor_p_value);
	
		    ## Determine the tumor genotype ##
		    
		    if($tumor_freq >= 0.25 || $tumor_p_value < 1e-06)
		    {
			$tumor_gt = $allele1 . $allele2;
			$tumor_gt = $allele2 . $allele2 if($tumor_freq > 0.90);
		    }
		}
    
		## Determine the somatic status ##
		
		my $somatic_status = "Unknown";
		my $p_value = 1;
    
		## SCENARIO 1: Tumor = Normal = Reference ##
		
		if($tumor_gt eq $normal_gt && $normal_gt eq $reference_gt)
		{
		    $somatic_status = "Reference";
    
		    ## Compile the read counts ##
		    
		    my $reads1 = $normal_reads1{$snp_key} + $tumor_reads1{$snp_key};
		    my $reads2 = $normal_reads2{$snp_key} + $tumor_reads2{$snp_key};
		    my $coverage = $reads1 + $reads2;
		    
		    ## Calculate a p-value for the combined read counts supporting Germline ##
			
		    $p_value = VarScan::FisherTest::calculate_p_value($coverage, 0, $reads1, $reads2, 0);
		    $p_value = VarScan::FisherTest::format_p_value($p_value);
		}
		
		## SCENARIO 2: Tumor = Normal != Reference ##	    
    
		elsif($tumor_gt eq $normal_gt && $normal_gt ne $reference_gt)
		{
		    ## Somatic Status: Germline ##
		    
		    $somatic_status = "Germline";
		    
		    ## Compile the read counts ##
		    
		    my $reads1 = $normal_reads1{$snp_key} + $tumor_reads1{$snp_key};
		    my $reads2 = $normal_reads2{$snp_key} + $tumor_reads2{$snp_key};
		    my $coverage = $reads1 + $reads2;
		    
		    ## Calculate a p-value for the combined read counts supporting Germline ##
			
		    $p_value = VarScan::FisherTest::calculate_p_value($coverage, 0, $reads1, $reads2, 0);
		    $p_value = VarScan::FisherTest::format_p_value($p_value);		
		}
    
		## SCENARIO 3: Tumor != Normal ##	    
    
		elsif($tumor_gt ne $normal_gt)
		{
		    ## Determine type of change ##
		    
		    if($normal_gt eq ($allele1 . $allele2))
		    {
			if($tumor_gt eq ($allele2 . $allele2))
			{
			    $somatic_status = "LOH2";
			}
			else
			{
			    $somatic_status = "LOH1";
			}
		    }
		    else
		    {
			$somatic_status = "Somatic";
		    }
		    
		    ## Adjust read counts to the coverage level of the lowest sample ##
	
		    my $adj_normal_reads1 = my $adj_normal_reads2 = my $adj_tumor_reads1 = my $adj_tumor_reads2 = 0;		
    
		    ## Tumor has less coverage: adjust to Tumor ##
		    
		    if($tumor_coverage < $normal_coverage)
		    {
			$adj_tumor_reads1 = $tumor_reads1{$snp_key};
			$adj_tumor_reads2 = $tumor_reads2{$snp_key};
			$adj_normal_reads2 = sprintf("%d", $normal_freq * $tumor_coverage);
			$adj_normal_reads1 = $tumor_coverage - $adj_normal_reads2;
		    }
		    
		    ## Normal has less coverage: adjust to Normal ##
    
		    else
		    {
			$adj_normal_reads1 = $normal_reads1{$snp_key};
			$adj_normal_reads2 = $normal_reads2{$snp_key};
			$adj_tumor_reads2 = sprintf("%d", $tumor_freq * $normal_coverage);
			$adj_tumor_reads1 = $normal_coverage - $adj_tumor_reads2;
		    }
    
		    ## Calculate adjusted p-value for normal-to-tumor change ##
    
		    $p_value = VarScan::FisherTest::calculate_p_value($adj_normal_reads1, $adj_normal_reads2, $adj_tumor_reads1, $adj_tumor_reads2, 0);
		    $p_value = VarScan::FisherTest::format_p_value($p_value);
		    
		}
    
		## Reformat variant frequency ##
		
		$normal_freq = sprintf("%.2f", $normal_freq * 100) . '%';
		$tumor_freq = sprintf("%.2f", $tumor_freq * 100) . '%';
    
		## Print to the outfile ##
    
		print OUTFILE "$snp_key\t";
		print OUTFILE "$normal_reads1{$snp_key}\t$normal_reads2{$snp_key}\t$normal_freq\t$normal_gt\t";
		print OUTFILE "$tumor_reads1{$snp_key}\t$tumor_reads2{$snp_key}\t$tumor_freq\t$tumor_gt\t";
		print OUTFILE "$somatic_status\t";
		print OUTFILE "$p_value\t";
		print OUTFILE "\n";
		
		$stats{$somatic_status}++;
	    }
	    elsif($normal_coverage >= 10)
	    {
		$stats{'snps_missing_10x_tumor'}++;
	    }
	    elsif($tumor_coverage >= 10)
	    {
		$stats{'snps_missing_10x_normal'}++;
	    }
	}
	else
	{
	    print "No normal coverage for $snp_key\n";
	}
	
    }

    close(OUTFILE) if($output_file);

    print "$stats{'snps_covered_both'} had coverage data for both Normal and Tumor\n";
    print "$stats{'snps_missing_10x_normal'} had <10x coverage in Normal\n";
    print "$stats{'snps_missing_10x_tumor'} had <10x coverage in Tumor\n";
    print "$stats{'snps_covered_both_10x'} had 10x coverage in both samples\n";
    print "$stats{'Germline'} were called as Germline\n";
    print "$stats{'LOH1'} were called as LOH (loss of reference allele)\n";
    print "$stats{'LOH2'} were called as LOH (loss of variant allele)\n";
    print "$stats{'Somatic'} were called as Somatic\n";
}



sub byChrPos
{
    (my $chrom_a, my $pos_a) = split(/\t/, $a);
    (my $chrom_b, my $pos_b) = split(/\t/, $b);

    $chrom_a cmp $chrom_b
    or
    $pos_a <=> $pos_b;
    
#    $chrom_a = 23 if($chrom_a =~ 'X');
#    $chrom_a = 24 if($chrom_a =~ 'Y');
    
}

=head1 AUTHOR

    Daniel C. Koboldt, << <dkoboldt at genome.wustl.edu> >>
    The Genome Center at Washington University School of Medicine
    St. Louis, Missouri, USA

=head1 COPYRIGHT

    Copyright 2009 Daniel C. Koboldt and Washington University
    All rights reserved.

=head1 LICENSE

    This program is free for non-commercial use.

=cut

1; # End of VarScan::VariantCalling
