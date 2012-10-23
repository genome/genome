package VarScan::ParseBowtie;

use warnings;
use strict;
use Getopt::Long;

=head1 NAME

VarScan::ParseBowtie - Parses Bowtie alignments

=head1 VERSION

    Version 1.03

=cut

our $VERSION = '1.03';

our $verbose = 1;

=head1 SYNOPSIS

    This module parses alignments in Bowtie format [http://bowtie-bio.sourceforge.net]

=head1 FUNCTIONS

=cut

## Define shared variables ##

my $fasta_file, my $quality_file, my $min_identity, my $min_align_score;
my $primer_trim, my $default_qual_score, my $min_qual_score;
my $output_dir, my $sample, my $output_alignments, my $output_snps, my $output_indels;

################################################################################

=head2	parse_alignments - parses the alignments file

=cut
################################################################################

sub parse_alignments
{
    (my $alignments_file) = @_;	#, my $output_file

    my %stats = ();
    $stats{'num_alignments'} = 0;

    ## Set default parameters ##
    my $min_align_score = 	25;
    my $min_identity = 		80;
    my $min_qual_score = 	15;
    my $default_qual_score = 	15;

    $fasta_file = $VarScan::fasta_file;
    $quality_file = $VarScan::quality_file;
    $min_identity = $VarScan::min_identity;
    $min_align_score = $VarScan::min_align_score;
    $primer_trim = $VarScan::primer_trim;
    $default_qual_score = $VarScan::default_qual_score;
    $min_qual_score = $VarScan::min_qual_score;
    $output_dir = $VarScan::output_dir;
    $sample = $VarScan::sample;
    
    $output_alignments = $VarScan::output_alignments;
    $output_snps = $VarScan::output_snps;
    $output_indels = $VarScan::output_indels;

    print "Minimum alignment score: $min_align_score\n" if($min_align_score);
    print "Minimum identity: $min_identity\n" if($min_identity);

    if($output_alignments)
    {
	    open(ALIGNMENTS, ">$output_alignments") or die "Can't open outfile: $!\n";
#	    print ALIGNMENTS "chrom\tposition\tref_allele\tvar_allele\tread_name\tread_pos\talign_strand\tqual_score\n";
    }

    if($output_snps)
    {
	    open(OUTFILE, ">$output_snps") or die "Can't open outfile: $!\n";
	    print OUTFILE "chrom\tposition\tref_allele\tvar_allele\tread_name\tread_pos\talign_strand\tqual_score\n";
    }
        
	my $input = new FileHandle ($alignments_file);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		
	
		if($lineCounter >= 1)
		{
			$stats{'num_alignments'}++;
                        
                        if($verbose)
                        {
                                print $stats{'num_alignments'} . " alignments parsed...\n" if(!($stats{'num_alignments'} % 10000))
                        }
			
			print ALIGNMENTS "$line\n" if($output_alignments);
			
			my @lineContents = split(/\t/, $line);
			my $read_name = $lineContents[0];
			my $align_strand = $lineContents[1];
			my $chromosome = $lineContents[2];
			my $chr_start = $lineContents[3];
			my $read_seq = $lineContents[4];
			my $read_qual = $lineContents[5];
			my $mismatch_list = $lineContents[7];
			
                        ## Adjust chromosome position to 1-based coordinates ##
                        
                        $chr_start++;
                        
			my $read_length = length($read_seq);
                        my $chr_stop = $chr_start + $read_length - 1;
                       
			if($mismatch_list)
			{
				my @mismatches = split(/\,/, $mismatch_list);
				foreach my $mismatch (@mismatches)
				{
					my $snp_chr_pos = my $qual_code = my $qual_score = 0;
					(my $read_pos, my $allele1, my $allele2) = split(/[\>\:]/, $mismatch);
				
                                	if($allele1 ne "N" && $allele2 ne "N")
					{
                                                if($align_strand eq "-")
                                                {
                                                        $snp_chr_pos = $chr_start + ($read_length - $read_pos - 1);
                                                        $qual_code = substr($read_qual, ($read_length - $read_pos - 1), 1);
                                                        ## Do NOT Correct alleles - bowtie gives read sequence in + orientation ##
                                                }
                                                else
                                                {
                                                        $snp_chr_pos = $chr_start + $read_pos;# - 1;
                                                        $qual_code = substr($read_qual, $read_pos, 1); #substr($read_qual, ($read_pos - 1), 1);
                                                }
                                                
                                                $qual_score = ord($qual_code) - 33;

                                                if($qual_score >= $min_qual_score)
                                                {
                                                        $stats{'num_snps'}++;
                                                        if($output_snps)
                                                        {
                                                                print OUTFILE "$chromosome\t$snp_chr_pos\t$allele1\t$allele2\t$read_name\t$read_pos\t$align_strand\t$qual_score\n";
                                                        }
                                                }
					}					
				}
			}
			
		}
	}
	
	close($input);

    ## Close the output file ##

    close(OUTFILE) if($output_snps);
    close(ALIGNMENTS) if($output_alignments);



    print "$stats{'num_alignments'} alignments were parsed\n";
    print "$stats{'num_snps'} substitution (SNP) events detected\n";

    return(0);
}




################################################################################

=head2	get_readcounts - gets read counts for a set of variants

=cut
################################################################################

sub get_readcounts
{
    (my $variants_file, my $alignments_file, my $output_file) = @_;

    my %stats = ();

    my %variant_regions = my %snps_by_position = my %indels_by_position = my %coverage_by_position = ();
    
    $min_qual_score = $VarScan::min_qual_score;
    
    print "Parsing SNPs...\n";
    ## Parse the combined SNPs ##
    
    my $input = new FileHandle ($variants_file);
    my $lineCounter = 0;
    while (<$input>)
    {
	    chomp;
	    my $line = $_;
	    $lineCounter++;

	    my @lineContents = split(/\t/, $line);

	    ## Determine type of variant entry ##

	    ## Header Line ##
	    if($lineContents[0] eq "chrom" || $lineContents[0] eq "ref_name")
	    {
		
	    }
	    ## INDEL ##
	    elsif($lineContents[6] && $lineContents[3] && ($lineContents[3] eq "INSERTION" || $lineContents[3] eq "DELETION")) ## INDELS ##
	    {
		    (my $chrom, my $chr_start, my $chr_stop, my $indel_type, my $indel_size, my $indel_context, my $num_reads, my $num_strands, my $avg_conf_score, my $avg_qual) = split(/\t/, $line);

		    my $indel_key = $chrom . ":" . $chr_start;
		    $indels_by_position{$indel_key} .= "$line\n";

		    ## Build a key using chromosome and first 2 bases of position for storing this indel ##
		    for(my $key_num = substr($chr_start, 0, 2); $key_num <= substr($chr_stop, 0, 2); $key_num++)
		    {
        		    my $chrom_key = $chrom . ":" . $key_num;
			    $variant_regions{$chrom_key}++;
		    }

		    $stats{$indel_type}++;
	    }
	    ## SNP ##
	    elsif($lineContents[0] && $lineContents[1] && $lineContents[2] && $lineContents[3]) ## SNPS ##
	    {
		    (my $chrom, my $position, my $allele1, my $allele2, my $num_reads, my $avg_qual, my $num_strands) = split(/\t/, $line);
    
		    my $position_key = $chrom . ":" . $position; #substr($position, 0, 2);
		    #$snps_by_position{$position_key} += $num_reads;
		    $snps_by_position{$position_key} .= $allele2;
		    
		    ## Build a key using chromosome and first 2 bases of position for storing this SNP ##
		    my $chrom_key = $chrom . ":" . substr($position, 0, 2);
		    $variant_regions{$chrom_key}++;			
		    
		    $stats{'SNP'}++;
	    }
	    else
	    {
		warn "Warning: Variant type not recognized: $line\n";
	    }


    }
    
    close($input);
    
    print $stats{'SNP'} . " SNPs loaded\n" if($stats{'SNP'});
    print $stats{'INSERTION'} . " insertions loaded\n" if($stats{'INSERTION'});
    print $stats{'DELETION'} . " deletions loaded\n" if($stats{'DELETION'});

    
    my %num_ref_reads = my %sum_ref_qual = my %ref_strands_seen = my %num_var_reads = my %num_N_reads = ();
    
    $input = new FileHandle ($alignments_file);
    $lineCounter = 0;
    
    while (<$input>)
    {
	    chomp;
	    my $line = $_;
	    $lineCounter++;		
    
	    if($lineCounter >= 1)
	    {
		    $stats{'num_alignments'}++;
		    
		    if($verbose)
		    {
			    print $stats{'num_alignments'} . " alignments parsed...\n" if(!($stats{'num_alignments'} % 100000))
		    }
		    
		    my @lineContents = split(/\t/, $line);
		    my $read_name = $lineContents[0];
		    my $align_strand = $lineContents[1];
		    my $chromosome = $lineContents[2];
		    my $chr_start = $lineContents[3];
		    my $read_seq = $lineContents[4];
		    my $read_qual = $lineContents[5];
		    my $mismatch_list = $lineContents[7];
		    
		    ## Adjust chromosome position to 1-based coordinates ##
		    
		    $chr_start++;
		    
		    my $read_length = length($read_seq);
		    my $chr_stop = $chr_start + $read_length - 1;
		    
		    ## Get any mismatch positions ##
		    my %mismatch_positions = ();
		    my %n_positions = ();
		    
		    if($mismatch_list)
		    {
			    my @mismatches = split(/\,/, $mismatch_list);
			    foreach my $mismatch (@mismatches)
			    {
				    (my $read_pos, my $allele1, my $allele2) = split(/[\>\:]/, $mismatch);
				    $read_pos++ if($align_strand eq "+");

				    if($allele1 eq "N" || $allele2 eq "N")
				    {
					    $n_positions{$read_pos} = 1;
				    }
				    else
				    {
					    $mismatch_positions{$read_pos} = 1;
				    }
			    }
		    }
		    
		    ## Go through every read position ##
		    
		    for(my $read_pos = 0; $read_pos < $read_length; $read_pos++)
		    {
			    my $ref_pos; my $qual_code; my $qual_score;
    
			    ## Determine base position quality and reference position ##
			    if($align_strand eq "-")
			    {
				    $ref_pos = $chr_start + ($read_length - $read_pos - 1);
			    }
			    else
			    {
				    $ref_pos = $chr_start + $read_pos - 1;
			    }
			    
			    my $position_key = $chromosome . ":" . $ref_pos;

			    my $read_base = "";
    
			    if($align_strand eq "-")
			    {
				    $qual_code = substr($read_qual, ($read_length - $read_pos - 1), 1);
				    $read_base = substr($read_seq, ($read_length - $read_pos - 1), 1);
			    }
			    else
			    {
				    $qual_code = substr($read_qual, ($read_pos - 1), 1);
				    $read_base = substr($read_seq, ($read_pos - 1), 1);
			    }

			    $qual_score = ord($qual_code) - 33;
			    
			    if($qual_score >= $min_qual_score)
			    {
				if($snps_by_position{$position_key} || $indels_by_position{$position_key})
				{
					$coverage_by_position{$position_key}++;
					
					## If base was N ##
					if($n_positions{$read_pos} || $qual_score < $min_qual_score)
					{
						$num_N_reads{$position_key}++;
					}
					
					## If base was mismatch ##
	
					elsif($mismatch_positions{$read_pos})
					{
						$num_var_reads{$position_key}++;
					}
					
					## Otherwise, base supports ref allele ##
					
					else
					{
						$ref_strands_seen{$position_key} .= $align_strand if(!$ref_strands_seen{$position_key} ||(length($ref_strands_seen{$position_key}) < 2 && $align_strand ne $ref_strands_seen{$position_key}));
						$num_ref_reads{$position_key}++;
						$sum_ref_qual{$position_key} += $qual_score;
					}
				}
			    }
		    }
		    
	    }
    }
    
    close($input);
	
	
    ## Open the outfile ##
    
    open(OUTFILE, ">$output_file") or die "Can't open outfile: $!\n";
    if($stats{'INSERTION'} || $stats{'DELETION'})
    {
	print OUTFILE "chrom\tchr_start\tchr_stop\tindel_type\tindel_size\t";
	print OUTFILE "read_coverage\treads1\treads2\treadsN\tavg_qual1\tavg_qual2\tstrands1\tstrands2\tindel_conf_score\tindel_in_context\n";
        print "Parsing indels again...\n";

    }
    else
    {
	print OUTFILE "chrom\tposition\tallele1\tallele2\t";
	print OUTFILE "read_coverage\treads1\treads2\treadsN\t";
	print OUTFILE "avg_ref_qual\tavg_var_qual\tnum_ref_strands\tnum_var_strands\tsnp_in_context\n";
        print "Parsing SNPs again...\n";
	$stats{'SNP'} = 0;

    }	    
    
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
		
	    }
	    elsif($lineContents[6] && $lineContents[3] && ($lineContents[3] eq "INSERTION" || $lineContents[3] eq "DELETION")) ## INDELS ##
	    {
		    (my $chrom, my $chr_start, my $chr_stop, my $indel_type, my $indel_size, , my $indel_context, my $reads2, my $num_var_strands, my $avg_conf_score, my $avg_var_qual) = split(/\t/, $line);
#		    my $indel_key = "$chrom\t$chr_start\t$chr_stop\t$indel_type\t$indel_size";
#		    my $reads1 = $num_ref_reads{$indel_key};
#		    my $readsN = $num_N_reads{$indel_key};
		    my $indel_key = "$chrom:$chr_start";
		    my $reads1 = $num_ref_reads{$indel_key};
		    my $readsN = $num_N_reads{$indel_key};
		    $reads1 = 0 if(!$reads1);
		    $readsN = 0 if(!$readsN);
		    my $read_coverage = $reads1 + $reads2 + $readsN;

		    my $num_ref_strands = 0;
		    $num_ref_strands = length($ref_strands_seen{$indel_key}) if($ref_strands_seen{$indel_key});

		    my $avg_ref_qual = $sum_ref_qual{$indel_key} / $reads1 if($reads1);
		    $avg_ref_qual = sprintf("%d", $avg_ref_qual) if($avg_ref_qual);
		    $avg_ref_qual = 0 if(!$avg_ref_qual);

		    print OUTFILE "$chrom\t$chr_start\t$chr_stop\t$indel_type\t$indel_size\t";
		    print OUTFILE "$read_coverage\t$reads1\t$reads2\t$readsN\t";
		    print OUTFILE "$avg_ref_qual\t$avg_var_qual\t$num_ref_strands\t$num_var_strands\t$avg_conf_score\t$indel_context\n";

	    }
	    elsif($lineContents[0] && $lineContents[1] && $lineContents[2] && $lineContents[3]) ## SNPS ##
	    {
		    (my $chrom, my $position, my $allele1, my $allele2, my $reads2, my $avg_var_qual, my $num_var_strands) = split(/\t/, $line);
		    my $position_key = $chrom . ":" . $position; #substr($position, 0, 2);
		    $num_ref_reads{$position_key} = 0 if(!$num_ref_reads{$position_key});
		    $num_var_reads{$position_key} = 0 if(!$num_var_reads{$position_key});
		    $num_N_reads{$position_key} = 0 if(!$num_N_reads{$position_key});
    
		    # Calculate read numbers ##
		    my $reads1 = $num_ref_reads{$position_key};
		    my $readsN = $num_N_reads{$position_key};				
		    my $read_coverage = $reads1 + $reads2 + $readsN;
    
		    my $avg_ref_qual = $sum_ref_qual{$position_key} / $reads1 if($reads1);
		    $avg_ref_qual = sprintf("%d", $avg_ref_qual) if($avg_ref_qual);
		    $avg_ref_qual = 0 if(!$avg_ref_qual);
    
		    my $num_ref_strands = 0;
		    $num_ref_strands = length($ref_strands_seen{$position_key}) if($ref_strands_seen{$position_key});
    
		    print OUTFILE "$chrom\t$position\t$allele1\t$allele2\t";
		    print OUTFILE "$read_coverage\t$reads1\t$reads2\t$readsN\t";
		    print OUTFILE "$avg_ref_qual\t$avg_var_qual\t$num_ref_strands\t$num_var_strands\n";
		    
		    $stats{'SNP'}++;
		    
	    }

    }
    
    close($input);
    
    close(OUTFILE);
	
    print "$stats{'SNP'} SNPs got read counts\n";	
    
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

1; # End of VarScan::ParseBowtie
