package VarScan::CallVariants;

use warnings;
use strict;

=head1 NAME

VarScan::CallVariants - The great new VarScan::CallVariants!

=head1 VERSION

Version 0.01

=cut

our $VERSION = '1.01';

our $min_average_quality =  15;
our $min_indel_size =       1;
our $max_indel_size =       100;
our $min_strand_support =   1;
our $min_read_coverage =    10;
our $min_reads2 =           2;
our $min_variant_freq =     0.01;
our $p_value_threshold =    1;#1e-04;

=head1 SYNOPSIS

Quick summary of what the module does.

Perhaps a little code snippet.

    use VarScan::CallVariants;

    my $foo = VarScan::CallVariants->new();
    ...

=head1 EXPORT

A list of functions that can be exported.  You can delete this section
if you don't export anything, such as for a purely object-oriented module.

=head1 FUNCTIONS

=cut

################################################################################################
=head2 GetSNPqualities - retrieve base qualities at SNP loci
#
################################################################################################
=cut

sub get_snp_qualities
{
        (my $infile, my $qual_db) = @_;
#	my $FileName = $output_dir . "/" . $sample_name . ".substitutions";
	
	
	## Open the output file ##
	
	open(OUTFILE, ">$infile.qual") or die "Can't open outfile: $!\n";
	
	## Open the input file ##
	
	my $input = new FileHandle ($infile);
	my $lineCounter = 0;

	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;

		if($lineCounter == 1)
		{
			print OUTFILE "$line\tqual_score\n";
		}
		if($lineCounter > 1)
		{
			my @lineContents = split(/\t/, $line);
			my $numContents = @lineContents;
			my $read_name = $lineContents[4];
			my $read_position = $lineContents[5];
			$read_position++;
			
			## Get the quality ##
			my $qpiece = $qual_db->seq($read_name, ($read_position - 0) => ($read_position + 0));			
			$qpiece = $qual_db->seq($read_name, ($read_position - 1) => ($read_position - 1)) if(!$qpiece);	
			## Convert to numeric ##
			
			my $qpiece_num = "";
			
			for(my $baseCounter = 0; $baseCounter < length($qpiece); $baseCounter++)
			{
				my $Qchar = substr($qpiece, $baseCounter, 1);
				my $Qnum = ord($Qchar) - 33;
				$qpiece_num .= " " if($qpiece_num);
				$qpiece_num .= $Qnum;
			}			
			
			
			## Determine qual score ##
			
			my $qual_score = $qpiece_num;
			
			## Append qual score to line ##
			
			my $newline = $line . "\t" . $qual_score;
			
			print OUTFILE "$newline\n";
#			print "$lineCounter\t$newline\n";		
#			return(0) if($lineCounter > 100);
		}
	}

	close($input);
	close(OUTFILE);
        
        return(0);
}



################################################################################################
=head2 Filter-and-call-SNPs - Filter SNP calls based on various criteria, and then output
#
################################################################################################
=cut

sub filter_and_call_snps
{
        (my $FileName, my $OutFileName) = @_;
	
	my %GenotypeStats = ();
        
        ## Open the output file ##
        
        open(OUTFILE, ">$OutFileName") or die "Can't open outfile: $!\n";
        print OUTFILE "chrom\tposition\tref\tvar\treads1\treads2\tp_value\tavgQual\tcontext\n";
        
	my $input = new FileHandle ($FileName);
	my $lineCounter = 0;
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		if($lineCounter > 1)
		{
                    $GenotypeStats{'unfiltered_snps'}++;
                    (my $chrom, my $position, my $allele1, my $allele2, my $coverage, my $reads1, my $reads2, my $avgQual, my $strands, my $context, my $p_value) = split(/\t/, $_);

                    ## Calculate allele frequency ##
                    
                    my $variant_freq = $reads2 / ($reads1 + $reads2);# * 100;
#                    $variant_freq = sprintf("%.1f", $variant_freq) . '%';

                    ## Make a call of the variant genotype ##

                    if($coverage >= $min_read_coverage)
                    {
                        if($reads2 >= $min_reads2)
                        {
                            if($variant_freq >= $min_variant_freq)
                            {
                                if($avgQual >= $min_average_quality)
                                {
                                    if($strands >= $min_strand_support)
                                    {
                                        if($p_value <= $p_value_threshold)
                                        {
                                            $GenotypeStats{'filtered_snps'}++;
                                            print OUTFILE "$chrom\t$position\t$allele1\t$allele2\t$reads1\t$reads2\t$p_value\t$avgQual\t$context\n";
                                        }
                                    }
                                }
                            }
                        }
                    }

                    
                }
        }
        close($input);

        close(OUTFILE);
        
        print "$GenotypeStats{'filtered_snps'} SNPs passed filters and were called\n";
        
        return(0);
}








################################################################################################
=head2 Filter-and-call-indels - Filter indel calls based on various criteria, and then output
#
################################################################################################
=cut

sub filter_and_call_indels
{
        (my $FileName, my $OutFileName) = @_;
	
	my %GenotypeStats = ();
        $GenotypeStats{'filtered_indels'} = 0;
        ## Open the output file ##
        
        open(OUTFILE, ">$OutFileName") or die "Can't open outfile: $!\n";
        print OUTFILE "chrom\tposition\tref\tvar\treads1\treads2\tp_value\tavgQual\tcontext\n";
        
        
	my $input = new FileHandle ($FileName);
	my $lineCounter = 0;
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		if($lineCounter > 1)
		{
                    $GenotypeStats{'unfiltered_indels'}++;
                    (my $chrom, my $chr_start, my $chr_stop, my $indel_type, my $indel_size, my $context, my $reads1, my $reads2, my $strands, my $conf_score, my $avgQual, my $p_value) = split(/\t/, $_);

                    ## Calculate coverage ##

                    my $coverage = $reads1 + $reads2;

                    ## Calculate allele frequency ##
                    
                    my $variant_freq = $reads2 / ($reads1 + $reads2); # * 100;
#                    $variant_freq = sprintf("%.1f", $variant_freq) . '%';

                    ## Make a call of the variant genotype ##
                    if($indel_size >= $min_indel_size && $indel_size <= $max_indel_size)
                    {
                        if($coverage >= $min_read_coverage)
                        {
                            if($reads2 >= $min_reads2 && $variant_freq >= $min_variant_freq)
                            {
                                if($avgQual >= $min_average_quality)
                                {
                                    if($strands >= $min_strand_support)
                                    {
                                        if($p_value <= $p_value_threshold)
                                        {
                                            $GenotypeStats{'filtered_indels'}++;
                                            print OUTFILE "$chrom\t$chr_start\t$chr_stop\t$indel_type\t$indel_size\t$reads1\t$reads2\t$p_value\t$avgQual\t$context\n";
                                        }
                                    }
                                }
                            }
                        }
                    }
                    
                }
        }
        close($input);

        close(OUTFILE);
        
        print "$GenotypeStats{'filtered_indels'} indels passed filters (size=$min_indel_size-$max_indel_size cov=$min_read_coverage reads2=$min_reads2 strands=$min_strand_support freq=$min_variant_freq min_pval=$p_value_threshold avgqual=$min_average_quality) and were called\n";
        
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




=head1 AUTHOR

Dan Koboldt, C<< <dankoboldt at gmail.com> >>

=head1 BUGS

Please report any bugs or feature requests to
C<bug-varscan-callvariants at rt.cpan.org>, or through the web interface at
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

Copyright 2009 Dan Koboldt, all rights reserved.

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.

=cut

1; # End of VarScan::CallVariants
