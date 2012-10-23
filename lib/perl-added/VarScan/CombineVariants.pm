package VarScan::CombineVariants;

use warnings;
use strict;

=head1 NAME

VarScan::CombineVariants - A VarScan module for combining variants detected in multiple reads

=head1 VERSION

Version 1.01

=cut

our $VERSION = '1.01';
our $min_average_quality =  20;
our $max_events_per_read = 5;   # The maximum # of translocation calls per read allowable
our $min_reads_per_translocation = 2;

=head1 SYNOPSIS

Combines variant predictions for SNPs, DNPs, insertions, deletions

    use VarScan::CombineVariants;

=head1 FUNCTIONS

=cut

################################################################################################
=head2 CombineSNPs - a subroutine that combines read substitution events into SNPs
################################################################################################
=cut

sub combine_substitutions
{
        (my $infile) = @_;
        my $outfile = $infile . ".combined";
	
	if(!($infile && -e $infile))
	{
		return("Substitutions file $infile not found");
	}	
	
	my %SNPcontexts = my %NumReads = my %SupportingReads = my %VariantAlleles = ();
	
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
		

		if($lineCounter > 1)
		{
			(my $chrom, my $position, my $allele1, my $allele2, my $read_name, my $read_pos, my $align_strand, my $snp_in_context, my $base_qual_score) = split(/\t/, $line);

                        $base_qual_score = $min_average_quality if(!$base_qual_score);
		
			my $snp_key = "$chrom\t$position\t$allele1\t$allele2";		
		
			if($allele2 =~ /[ACGT]/)
			{
				$CombineStats{'snp_events'}++;
				
				$NumReads{$snp_key}++;
				$SupportingReads{$snp_key} .= "," if($SupportingReads{$snp_key});			
				if($base_qual_score)
				{
					$SupportingReads{$snp_key} .= "$read_name($allele2:$align_strand:$base_qual_score)";				
				}
				else
				{
					$SupportingReads{$snp_key} .= "$read_name($allele2:$align_strand)";
				}

				$SNPcontexts{$snp_key} = $snp_in_context if(!$SNPcontexts{$snp_key} || length($snp_in_context) > length($SNPcontexts{$snp_key}));
			
                                $SNPqualitySum{$snp_key} = 0 if(!$SNPqualitySum{$snp_key});
				$SNPqualitySum{$snp_key} += $base_qual_score;
				if(!$base_qual_score)
                                {
                                    print "Somehow got $line\n";
                                }
				$VariantAlleles{$snp_key} = "" if(!$VariantAlleles{$snp_key});
	
				if(!($VariantAlleles{$snp_key} =~ $allele2))
				{
					$VariantAlleles{$snp_key} .= "/" if($VariantAlleles{$snp_key});					
					$VariantAlleles{$snp_key} .= $allele2 
				}
			}

		}
	}
	
	open(OUTFILE, ">$outfile") or die "Can't open outfile: $!\n";
#	print OUTFILE "chrom\tposition\tref_allele\tvar_allele\tcontext\tnum_reads\treads\n";
	print OUTFILE "chrom\tposition\tref\tvar\treads2\tavgQual\tstrands\n";	
	foreach my $snp_key (sort keys %SupportingReads)
	{
		$CombineStats{'unique_snps'}++;
		## Determine average base quality ##
		my $avg_base_qual = $SNPqualitySum{$snp_key} / $NumReads{$snp_key};
		$avg_base_qual = sprintf("%d", $avg_base_qual);
		my $num_reads = $NumReads{$snp_key};

		my $num_strands = 1;
		
		if($SupportingReads{$snp_key} =~ '\+' && $SupportingReads{$snp_key} =~ '\-')
		{
			$num_strands = 2;
		}

#		print OUTFILE "$snp_key\t$SNPcontexts{$snp_key}\t$NumReads{$snp_key}\t$SupportingReads{$snp_key}\n";
#		print OUTFILE "$snp_key\t$num_reads\t$avg_base_qual\t$num_strands\n";
		print OUTFILE "$snp_key\t$num_reads\t$avg_base_qual\t$num_strands\t$SNPcontexts{$snp_key}\t$SupportingReads{$snp_key}\n";
	}
	close(OUTFILE);


	print "$CombineStats{'snp_events'} substitution events\n";
	print "$CombineStats{'unique_snps'} unique SNPs\n";

    return(0);
}



################################################################################################
=head2 # CombineIndels - a subroutine that combines read indel events into indels
################################################################################################
=cut

sub combine_indels
{
        (my $infile) = @_;
        my $outfile = $infile . ".combined";
	
	if(!($infile && -e $infile))
	{
		return("Indels file $infile not found");
	}

        ## Declare or reset variables ##
	my %IndelContext = my %IndelPrecisionScore = my %NumReads = my %SupportingReads = ();
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
		

		if($lineCounter > 1)
		{
			my @lineContents = split(/\t/, $line);
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
                        $qual_score = $min_average_quality if(!$qual_score);
                
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
			$IndelContext{$indel_key} = $indel_context if(!$IndelContext{$indel_key} || ($conf_score >= $IndelPrecisionScore{$indel_key} && length($indel_context) > length($IndelContext{$indel_key})));		
		
			## Take the highest precision score ##
			$IndelPrecisionScore{$indel_key} = $conf_score if(!$IndelPrecisionScore{$indel_key} || $conf_score > $IndelPrecisionScore{$indel_key});		
		}
	}
	
	open(OUTFILE, ">$outfile") or die "Can't open outfile: $!\n";
	print OUTFILE "chrom\tchr_start\tchr_stop\tindel_type\tindel_size\tindel_in_context\tnum_reads\tnum_reads\tnum_strands\tavg_conf_score\tavg_base_qual\n";	
	foreach my $indel_key (sort keys %NumReads)
	{
		$CombineStats{'unique_indels'}++;
		my $num_reads = $NumReads{$indel_key};
		my $context = $IndelContext{$indel_key};
		my $num_strands = 1;
		if($SupportingReads{$indel_key} =~ '\+' && $SupportingReads{$indel_key} =~ '\-')
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
=head2 # combine_translocations - a subroutine that combines read translocations into translocs
################################################################################################
=cut

sub combine_translocations
{
        (my $infile, my $band_file) = @_;
        my $outfile = $infile . ".combined";
	
	if(!($infile && -e $infile))
	{
		return("Translocations file $infile not found");
	}	

        ## If a band file (BED formats for cyto bands, etc., was provided, use it ##
        my %genome_bands = ();        

        if($band_file && -e $band_file)
        {
            my $input = new FileHandle ($band_file);
            my $bandCounter = 0;
            while (<$input>)
            {
                    chomp;
                    my $line = $_;
                    $bandCounter++;
                    (my $chrom, my $start, my $stop, my $band) = split(/\t/, $line);
                    $genome_bands{$chrom} .= "$start\t$stop\t$band\n";
            }
            
            close($input);
            print "$bandCounter cyto bands loaded\n";
        }

	
	my %CombineStats = ();
	$CombineStats{'translocation_events'} = 0;
	$CombineStats{'unique_translocations'} = 0;
	
	my %NumReads = my %SupportingReads = my %EventsPerRead = ();
        
	## Open the input file ##
	
	my $input = new FileHandle ($infile);
	my $lineCounter = 0;

	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		

		if($lineCounter > 1)
		{
                        my @lineContents = split(/\t/, $line);
                        my $translocation = $lineContents[0];
                        my $read_name = $lineContents[1];
                        my $chrom1 = $lineContents[3];
                        my $chr_start1 = $lineContents[4];
#                        my $chr_stop1 = $lineContents[5];
                        my $chrom2 = $lineContents[10];
                        my $chr_start2 = $lineContents[11];

			my $translocation_key = "$translocation";		

                        ## If bands were provided, use them to generate a new key ##
                        
                        if($genome_bands{$chrom1} && $genome_bands{$chrom2})
                        {
                            my $band1 = coordinates_to_band($genome_bands{$chrom1}, $chr_start1);
                            my $band2 = coordinates_to_band($genome_bands{$chrom2}, $chr_start2);
                            $translocation_key = $chrom1 . $band1 . "-" . $chrom2 . $band2;
                        }

			$CombineStats{'translocation_events'}++;
			$NumReads{$translocation_key}++;
                        $SupportingReads{$translocation_key} .= "," if($SupportingReads{$translocation_key});			
                        $SupportingReads{$translocation_key} .= "$read_name";
                        $EventsPerRead{$read_name}++;
		}
	}
	
	open(OUTFILE, ">$outfile") or die "Can't open outfile: $!\n";
#	print OUTFILE "chrom\tposition\tref_allele\tvar_allele\tcontext\tnum_reads\treads\n";
	print OUTFILE "translocation\treads2\treads\n";	
	foreach my $translocation_key (sort keys %SupportingReads)
	{
		$CombineStats{'unique_translocations'}++;
		my $num_reads = $NumReads{$translocation_key};

                ## Get the individual reads ##
                
                my @supporting_reads = split(/\,/, $SupportingReads{$translocation_key});
                my $events_sum = 0;
                foreach my $read_name (@supporting_reads)
                {
                    if($read_name && $EventsPerRead{$read_name})
                    {
                        $events_sum += $EventsPerRead{$read_name};
                    }
                }

                my $avg_events_per_read = $events_sum / $num_reads;
                $avg_events_per_read = sprintf("%.2f", $avg_events_per_read);

                if($num_reads >= $min_reads_per_translocation && $avg_events_per_read < $max_events_per_read)
                {
                    print OUTFILE "$translocation_key\t$num_reads\t$avg_events_per_read\t$SupportingReads{$translocation_key}\n";
                    $CombineStats{'translocations_passed_filters'}++;
                }
	}
	close(OUTFILE);

	print "$CombineStats{'translocation_events'} translocation events\n";
	print "$CombineStats{'unique_translocations'} unique translocations\n";
        print "$CombineStats{'translocations_passed_filters'} passed events-per-read filters\n";
    return(0);
}



################################################################################################
=head2 # coordinates_to_band - match a chromosomal coordinate to a cytogenetic band
################################################################################################
=cut

sub coordinates_to_band
{
    (my $bands, my $position) = @_;
    
    my @bands = split(/\n/, $bands);
    
    foreach my $band (@bands)
    {
        (my $band_start, my $band_stop, my $band_name) = split(/\t/, $band);
        
        if($band_name && $position >= $band_start && $position <= $band_stop)
        {
            return($band_name);
        }
    }
    
    return(0);
}

=head1 AUTHOR

Dan Koboldt, C<< <dankoboldt at gmail.com> >>

=head1 BUGS

Please report any bugs or feature requests to
C<bug-varscan-combinevariants at rt.cpan.org>, or through the web interface at
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













1; # End of VarScan::CombineVariants
