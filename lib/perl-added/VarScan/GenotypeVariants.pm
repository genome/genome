package VarScan::GenotypeVariants;

use warnings;
use strict;

=head1 NAME

VarScan::GenotypeVariants - The great new VarScan::GenotypeVariants!

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';

=head1 SYNOPSIS

Quick summary of what the module does.

Perhaps a little code snippet.

    use VarScan::GenotypeVariants;

    my $foo = VarScan::GenotypeVariants->new();
    ...

=head1 EXPORT

A list of functions that can be exported.  You can delete this section
if you don't export anything, such as for a purely object-oriented module.

=head1 FUNCTIONS

=cut

################################################################################################
=head2 GenotypeSNPs - a subroutine that determines read counts for SNPs
#
################################################################################################
=cut

sub genotype_snps
{
        (my $FileName, my $BlocksFileName) = @_;
	my $OutFileName = $FileName . ".readcounts";	
	
	my %GenotypeStats = ();
	my %VariantRegions = my %SNPsByPosition = my %CoverageByPosition = ();
	
#	print "Parsing SNPs...\n";
	## Parse the combined SNPs ##

	my $input = new FileHandle ($FileName);
	my $lineCounter = 0;
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		if($lineCounter > 1)
		{
#			(my $chrom, my $position, my $allele1, my $allele2, my $context, my $num_reads, my $reads) = split(/\t/, $line);			
			(my $chrom, my $position, my $allele1, my $allele2, my $num_reads, my $avg_qual, my $num_strands) = split(/\t/, $line);
			my $position_key = $chrom . ":" . $position; #substr($position, 0, 2);
			$SNPsByPosition{$position_key} += $num_reads;
			
			## Build a key using chromosome and first 2 bases of position for storing this SNP ##
			my $chrom_key = $chrom . ":" . substr($position, 0, 2);
			$VariantRegions{$chrom_key}++;			
			
			$GenotypeStats{'num_snps'}++;
		}
	}

	close($input);

#	print $GenotypeStats{'num_snps'} . " SNPs loaded\n";


	## Parse the alignment blocks ##

	$input = new FileHandle ($BlocksFileName);
	$lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		if($lineCounter > 1)  # && $lineCounter < 50000)
		{
			$GenotypeStats{'total_alignments'}++;
#			print "$GenotypeStats{'total_alignments'} alignments\n" if(!($GenotypeStats{'total_alignments'} % 1000));
		
			my @lineContents = split(/\s+/, $line);
			my $chrom = $lineContents[0];
			my $chr_start = $lineContents[1];
			my $chr_stop = $lineContents[2];
			my $strand = $lineContents[3];
			my $read_name = $lineContents[4];
	
			
			## Get possible chrom position keys ##
			
#			my $key_start = substr($chr_start, 0, 2);
#			my $key_stop = substr($chr_stop, 0, 2);
			
#			for(my $key_number = $key_start; $key_number <= $key_stop; $key_number++)
#			{
#				my $chrom_key = $chrom . ":" . $key_number;
#				if($VariantRegions{$chrom_key})
#				{		
					my $position_key;
					
					for(my $position = $chr_start; $position <= $chr_stop; $position++)
					{
						$position_key = $chrom . ":" . $position; #substr($position, 0, 2);
						if($SNPsByPosition{$position_key})
						{
							$CoverageByPosition{$position_key}++;
						}
					}
#				}			
#			}
		

		}
	}
	
	close($input);


	## Open the outfile ##
	
	open(OUTFILE, ">$OutFileName") or die "Can't open outfile: $!\n";
	print OUTFILE "chrom\tposition\tref\tvar\tcov\treads1\treads2\tavgQual\tstrands\n";

	$input = new FileHandle ($FileName);
	$lineCounter = 0;
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		if($lineCounter == 1)
		{
		}
		if($lineCounter > 1)
		{
			(my $chrom, my $position, my $allele1, my $allele2, my $num_reads, my $avg_qual, my $num_strands, my $context, my $supporting_reads) = split(/\t/, $line);
			my $position_key = $chrom . ":" . $position; #substr($position, 0, 2);
			
			my $read_coverage = $CoverageByPosition{$position_key};
			## Calculate $reads1 ##
			my $num_wt_reads = $CoverageByPosition{$position_key} - $SNPsByPosition{$position_key};
			
			if($num_wt_reads < 0)
			{
				print "Warning: Reads1 calculated to be less than zero: $line\n";
				exit(1);
			}

			print OUTFILE "$chrom\t$position\t$allele1\t$allele2\t$read_coverage\t$num_wt_reads\t$num_reads\t$avg_qual\t$num_strands\t$context\t$supporting_reads\n";			
		}
	}

	close($input);

	close(OUTFILE);

    return(0);
}




################################################################################################
=head2 EvaluateSNPs - evaluate the support for SNPs
#
################################################################################################
=cut

sub evaluate_snps
{
    my $min_p_threshold = 0;
    my $expected_error_rate = 0.1;
    
        (my $FileName) = @_;
	my $OutFileName = $FileName . ".p_value";	
	
	my %GenotypeStats = ();
	my %FisherSave = ();
        
        ## Open the output file ##
        
        open(OUTFILE, ">$OutFileName") or die "Can't open outfile: $!\n";
        print OUTFILE "chrom\tposition\tref\tvar\tcov\treads1\treads2\tavgQual\tstrands\tcontext\tp_value\n";
        
	my $input = new FileHandle ($FileName);
	my $lineCounter = 0;
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		if($lineCounter > 1)
		{
                    (my $chrom, my $position, my $allele1, my $allele2, my $coverage, my $reads1, my $reads2, my $avgQual, my $strands, my $context, my $read_support) = split(/\t/, $_);
                    
                    ## Calculate expected reads2/reads1 ##
                    
#                    my $exp_reads2 = 0;
                    my $exp_reads2 = sprintf("%d", ($coverage * $expected_error_rate));
                    my $exp_reads1 = ($reads1 + $reads2) - $exp_reads2;
                    
                   
                    ## Calculate p-value based on read counts ##
                    
                    my @FETtable;
                    push @FETtable, $reads2;
                    push @FETtable, $reads1;
                    push @FETtable, $exp_reads2;
                    push @FETtable, $exp_reads1;
                    
                    my $p_value;
                    
                    if(defined($FisherSave{$reads1 . '-' . $reads2}))
                    {
                        $p_value = $FisherSave{$reads1 . '-' . $reads2};
                    }
                    else
                    {
                        $p_value = VarScan::FisherTest::Test(@FETtable, $min_p_threshold);                    
                        $FisherSave{$reads1 . '-' . $reads2} = $p_value;
                    }

                    ## Calculate allele frequency ##
                    
                    my $variant_freq = $reads2 / ($reads1 + $reads2) * 100;
                    $variant_freq = sprintf("%.1f", $variant_freq) . '%';

                    ## Make a call of the variant genotype ##
                    
#                    my $variant_call = "het";
#my $chrom, my $position, my $allele1, my $allele2, my $coverage, my $reads1, my $reads2, my $avgQual, my $strands, my $context, my $read_support
                    print OUTFILE "$chrom\t$position\t$allele1\t$allele2\t$coverage\t$reads1\t$reads2\t$avgQual\t$strands\t$context\t";
                    print OUTFILE "$p_value\n";
                    
                }
        }
        close($input);

        close(OUTFILE);
        
        return(0);
}


################################################################################################
=head2 genotype_indels - get read counts for indels
#
################################################################################################
=cut

sub genotype_indels
{
    (my $indels_file, my $blocks_file) = @_;
    
    my $outfile = $indels_file . ".readcounts";

	## Parse the indels ##
	
	my %GenotypeStats = ();
	my %IndelsByPosition = my %VariantRegions = my %CoverageByPosition = ();
	$GenotypeStats{'num_indels'} = $GenotypeStats{'total_alignments'} = $GenotypeStats{'num_genotyped'} = 0;

	my $input = new FileHandle ($indels_file);
	my $lineCounter = 0;
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		my @lineContents = split(/\t/, $line);		
		
		if($lineContents[0] && $lineContents[0] ne "chrom" && $lineContents[0] ne "ref_name")
		{
			my $chrom = $lineContents[0];
			my $chr_start = $lineContents[1];
			my $chr_stop = $lineContents[2];
			my $indel_type = $lineContents[3];
			my $indel_size = $lineContents[4];
			my $indel_context = $lineContents[5];
			my $num_reads = $lineContents[6];
			my $chrom_key = $chrom . ":" . substr($chr_start, 0, 4);
			$VariantRegions{$chrom_key} .= "$chrom\t$chr_start\t$chr_stop\t$indel_type\t$indel_size\t$num_reads\n";
			
			$GenotypeStats{$indel_type}++;
			$GenotypeStats{'num_indels'}++;
		}
	}

	close($input);

	print "$GenotypeStats{'num_indels'} indels loaded\n";


	## Parse the alignment blocks file ##

	$input = new FileHandle ($blocks_file);
	$lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		if($lineCounter > 1)	# && $lineCounter < 10000)
		{
			$GenotypeStats{'total_alignments'}++;
			print "$GenotypeStats{'total_alignments'} alignments\n" if(!($GenotypeStats{'total_alignments'} % 10000));
		
			my @lineContents = split(/\s+/, $line);
			my $chrom = $lineContents[0];
			my $chr_start = $lineContents[1];
			my $chr_stop = $lineContents[2];
			my $strand = $lineContents[3];
			my $read_name = $lineContents[4];
	
			my $chrom_key = $chrom . ":" . substr($chr_start, 0, 4);
			if($VariantRegions{$chrom_key})
			{		
				my @regionIndels = split(/\n/, $VariantRegions{$chrom_key});
				my $numRegionIndels = @regionIndels;
				
				for(my $indelCounter = 0; $indelCounter < $numRegionIndels; $indelCounter++)
				{
					(my $indel_chrom, my $indel_chr_start, my $indel_chr_stop, my $indel_type, my $indel_size, my $indel_num_reads) = split(/\t/, $regionIndels[$indelCounter]);
									
					## Build a key for the indel ##
					
					my $indel_key = $indel_chrom . ":" . $indel_chr_start . "-" . $indel_chr_stop;
									
					## Check to see if our block totally spans the indel coordinates ##
					
					if($chr_start < $indel_chr_start && $chr_stop > $indel_chr_stop)
					{
						$CoverageByPosition{$indel_key}++;			
					}
					## Or if block is completely within a deleted region ##
					elsif($indel_chr_start < $chr_start && $indel_chr_stop > $chr_stop)
					{
						$CoverageByPosition{$indel_key}++;
					}
					## Or if they overlap but not completely ##
					elsif($indel_chr_stop > $chr_start && $indel_chr_start < $chr_stop)
					{
						## Calculate overlap bases ##
						
						my $overlap_bases = 0;
						
						## Indel Starts within block but goes beyond ##
						if($indel_chr_start > $chr_start && $indel_chr_start < $chr_stop) 
						{
							$overlap_bases = $chr_stop - $indel_chr_start;
						}
						## Starts before block but goes into it ##
						elsif($indel_chr_stop > $chr_start && $indel_chr_stop < $chr_stop)
						{
							$overlap_bases = $indel_chr_stop - $chr_start;
						}
												
						## Count read as supporting wild-type if it has at least 7 wt bases ##
						if($overlap_bases >= 7)
						{
							$CoverageByPosition{$indel_key}++;
							my $read_genotype = "WT";
#							print READGTS "$indel_chrom\t$indel_chr_start\t$indel_chr_stop\t$indel_type\t$indel_size\t$read_name\t$read_genotype\t$strand\t$overlap_bases\t$chrom:$chr_start-$chr_stop\n";						
						}
					}
				}
				
			}			

		}
	}
	
	close($input);
	
	print "$GenotypeStats{'total_alignments'} alignment blocks parsed\n";


	open(READCOUNTS, ">$outfile") or die "Can't open outfile: $!\n";
	print READCOUNTS "chrom\tchr_start\tchr_stop\tindel_type\tindel_size\tindel_context\treads1\treads2\tnum_strands\tconf_score\tavg_qual_score\n";

	## Now, re-parse the indel file and build genotype submission file ##

	$input = new FileHandle ($indels_file);
	$lineCounter = 0;
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		my @lineContents = split(/\t/, $line);
		
		if($lineContents[0] && $lineContents[0] ne "chrom")
		{
			my $chrom = $lineContents[0];
			my $chr_start = $lineContents[1];
			my $chr_stop = $lineContents[2];
			my $indel_type = $lineContents[3];
			my $indel_size = $lineContents[4];
			my $indel_context = $lineContents[5];
			my $num_reads = $lineContents[6];
			my $num_strands = $lineContents[7];
			my $conf_score = $lineContents[8];
			my @tempArray = split(/[\[\/\]]/, $indel_context);
			my $allele1 = $tempArray[1];
			my $allele2 = $tempArray[2];
			my $avg_qual_score = $lineContents[9];
			
			my $position_key = $chrom . ":" . $chr_start . "-" . $chr_stop; #substr($position, 0, 2);
			my $chr_name = $chrom;
#			$chr_name =~ s/[^0-9XYM]//g;
			
			## Get wt readcount ##
			my $reads_wt = $CoverageByPosition{$position_key};
			$reads_wt = 0 if(!$reads_wt);
			
			## Print to simple variants file with readcounts ##

			print READCOUNTS "$chrom\t$chr_start\t$chr_stop\t$indel_type\t$indel_size\t$indel_context\t";
			print READCOUNTS "$reads_wt\t$num_reads\t$num_strands\t$conf_score\t$avg_qual_score\n";

			my $discovery_string = "breakPointRead($allele1:$allele2:$avg_qual_score:reads1=$reads_wt:reads2=$num_reads)";

			$GenotypeStats{'num_genotyped'}++;		
		
		}
	}

	close(READCOUNTS);
	print "$GenotypeStats{'num_genotyped'} indels genotyped\n";

    return(0);
}




################################################################################################
=head2 EvaluateIndels - evaluate the support for indels
#
################################################################################################
=cut

sub evaluate_indels
{
    my $min_p_threshold = 0;
    my $expected_error_rate = 0.1;
    
        (my $FileName) = @_;
	my $OutFileName = $FileName . ".p_value";	
	
	my %GenotypeStats = ();
	my %FisherSave = ();
        
        ## Open the output file ##
        
        open(OUTFILE, ">$OutFileName") or die "Can't open outfile: $!\n";
        print OUTFILE "chrom\tchr_start\tchr_stop\tindel_type\tindel_size\tindel_context\treads1\treads2\tnum_strands\tconf_score\tavg_qual_score\tp_value\n";
        
	my $input = new FileHandle ($FileName);
	my $lineCounter = 0;
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		if($lineCounter > 1)
		{
                    (my $chrom, my $chr_start, my $chr_stop, my $indel_type, my $indel_size, my $indel_context, my $reads1, my $reads2, my $num_strands, my $conf_score, my $avg_qual_score) = split(/\t/, $_);
                    
                    ## Calculate expected reads2/reads1 ##
                    my $coverage = $reads1 + $reads2;                    
#                    my $exp_reads2 = 0;
                    my $exp_reads2 = sprintf("%d", ($coverage * $expected_error_rate));
                    my $exp_reads1 = ($reads1 + $reads2) - $exp_reads2;
                    
                   
                    ## Calculate p-value based on read counts ##
                    
                    my @FETtable;
                    push @FETtable, $reads2;
                    push @FETtable, $reads1;
                    push @FETtable, $exp_reads2;
                    push @FETtable, $exp_reads1;
                    
                    my $p_value;
                    
                    if(defined($FisherSave{$reads1 . '-' . $reads2}))
                    {
                        $p_value = $FisherSave{$reads1 . '-' . $reads2};
                    }
                    else
                    {
                        $p_value = VarScan::FisherTest::Test(@FETtable, $min_p_threshold);                    
                        $FisherSave{$reads1 . '-' . $reads2} = $p_value;
                    }

                    ## Calculate allele frequency ##
                    
                    my $variant_freq = $reads2 / ($reads1 + $reads2) * 100;
                    $variant_freq = sprintf("%.1f", $variant_freq) . '%';

                    ## Make a call of the variant genotype ##
                    
#                    my $variant_call = "het";
#my $chrom, my $position, my $allele1, my $allele2, my $coverage, my $reads1, my $reads2, my $avgQual, my $strands, my $context, my $read_support
                    print OUTFILE "$chrom\t$chr_start\t$chr_stop\t$indel_type\t$indel_size\t$indel_context\t$reads1\t$reads2\t$num_strands\t$conf_score\t$avg_qual_score\t";
                    print OUTFILE "$p_value\n";
                    
                }
        }
        close($input);

        close(OUTFILE);
        
        return(0);
}






=head1 AUTHOR

Dan Koboldt, C<< <dankoboldt at gmail.com> >>

=head1 BUGS

Please report any bugs or feature requests to
C<bug-varscan-genotypevariants at rt.cpan.org>, or through the web interface at
L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=VarScan>.
I will be notified, and then you'll automatically be notified of progress on
your bug as I make changes.

=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc VarScan

You can also look for information at:

=over 4

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/VarScan>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/VarScan>

=item * RT: CPAN's request tracker

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=VarScan>

=item * Search CPAN

L<http://search.cpan.org/dist/VarScan>

=back

=head1 ACKNOWLEDGEMENTS

=head1 COPYRIGHT & LICENSE

Copyright 2008 Dan Koboldt, all rights reserved.

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.

=cut

1; # End of VarScan::GenotypeVariants
