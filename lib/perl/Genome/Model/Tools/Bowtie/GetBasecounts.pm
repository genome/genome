
package Genome::Model::Tools::Bowtie::GetBasecounts;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# GetBasecounts.pm - 	Get unmapped/poorly-mapped reads by model id
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@watson.wustl.edu)
#
#	CREATED:	02/25/2009 by D.K.
#	MODIFIED:	02/25/2009 by D.K.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;

use FileHandle;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

class Genome::Model::Tools::Bowtie::GetBasecounts {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		variants_file	=> { is => 'Text', doc => "File containing combined Bowtie SNPs" },
		output_file	=> { is => 'Text', doc => "File to receive combined SNPs" },
		blocks_file	=> { is => 'Text', doc => "File containing Bowtie alignment blocks", is_optional => 1 },
		alignments_file	=> { is => 'Text', doc => "File containing Bowtie alignment blocks", is_optional => 1 },
                min_qual_score	=> { is => 'Text', doc => "Minimum quality score to include SNP", is_optional => 1 },
		verbose	=> { is => 'Text', doc => "Print interim progress updates", is_optional => 1 },
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Count the number of reads supporting every possible base at every position"                 
}

sub help_synopsis {
    return <<EOS
This command retrieves the locations of unplaced reads for a given genome model
EXAMPLE:	gmt bowtie parse-alignments --alignments-file bowtie.txt
EOS
}

sub help_detail {                           # this is what the user will see with the longer version of help. <---
    return <<EOS 

EOS
}


################################################################################################
# Execute - the main program logic
#
################################################################################################

sub execute {                               # replace with real execution logic.
	my $self = shift;

	## Get required parameters ##
	my $variants_file = $self->variants_file;
        my $blocks_file = $self->blocks_file;
        my $alignments_file = $self->alignments_file if(defined($self->alignments_file));
	my $outfile = $self->output_file;
	my $verbose = 1 if(defined($self->verbose));
        my $min_qual_score = 0;
        $min_qual_score = $self->min_qual_score if(defined($self->min_qual_score));

	my %GenotypeStats = ();
	my %VariantRegions = my %SNPsByPosition = my %CoverageByPosition = ();
	
	print "Parsing SNPs...\n";
	## Parse the combined SNPs ##

	my $input = new FileHandle ($variants_file);
	my $lineCounter = 0;
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		if($lineCounter > 1)# && $lineCounter < 1000)
		{
#			(my $chrom, my $position, my $allele1, my $allele2, my $context, my $num_reads, my $reads) = split(/\t/, $line);			
			(my $chrom, my $position, my $allele1, my $allele2, my $num_reads, my $avg_qual, my $num_strands) = split(/\t/, $line);

#                        if($num_reads >= 10 && $avg_qual >= 20)
 #                       {
                                my $position_key = $chrom . ":" . $position; #substr($position, 0, 2);
                                $SNPsByPosition{$position_key} += $num_reads;
                                
                                ## Build a key using chromosome and first 2 bases of position for storing this SNP ##
                                my $chrom_key = $chrom . ":" . substr($position, 0, 2);
                                $VariantRegions{$chrom_key}++;			
                                
                               $GenotypeStats{'num_snps'}++;
 #                       }
		}
	}

	close($input);

	print $GenotypeStats{'num_snps'} . " SNPs loaded\n";

	my %num_ref_reads = my %sum_ref_qual = my %ref_strands_seen = my %num_var_reads = my %num_N_reads = ();

	if($alignments_file)
	{
		my $input = new FileHandle ($alignments_file);
		my $lineCounter = 0;
		
		while (<$input>)
		{
			chomp;
			my $line = $_;
			$lineCounter++;		
		
			if($lineCounter >= 1)
			{
				$GenotypeStats{'num_alignments'}++;
				
				if($verbose)
				{
					print $GenotypeStats{'num_alignments'} . " alignments parsed...\n" if(!($GenotypeStats{'num_alignments'} % 10000))
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
						
						if($allele1 eq "N" || $allele2 eq "N")
						{
							$n_positions{$read_pos} = 1;
						}
						else
						{
							$mismatch_positions{$read_pos} = $allele2;
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
					my $base_code = "";

					if($SNPsByPosition{$position_key})
					{
#						$CoverageByPosition{$position_key}++;

						if($align_strand eq "-")
						{
							$qual_code = substr($read_qual, ($read_length - $read_pos - 1), 1);
							$base_code = substr($read_seq, ($read_length - $read_pos - 1), 1);
						}
						else
						{
							$qual_code = substr($read_qual, ($read_pos - 1), 1);
							$base_code = substr($read_seq, ($read_pos - 1), 1);
						}
	
						$qual_score = ord($qual_code) - 33;

						my $hash_name = $position_key . ":" . $base_code;
						$CoverageByPosition{$hash_name}++;
						
					}
				}
				
#				return(0) if($lineCounter > 1000);
			}
		}
		
		close($input);


		## Open the outfile ##
		
		open(OUTFILE, ">$outfile") or die "Can't open outfile: $!\n";
		print OUTFILE "chrom\tposition\tallele1\tallele2\t";
		print OUTFILE "read_coverage\treads1\treads2\treadsN\t";
		print OUTFILE "avg_ref_qual\tavg_var_qual\tnum_ref_strands\tnum_var_strands\n";		
	
		print "Parsing SNPs again...\n";
	
		$input = new FileHandle ($variants_file);
		$lineCounter = 0;
		while (<$input>)
		{
			chomp;
			my $line = $_;
			$lineCounter++;
			
			if($lineCounter == 1)
			{
			}
			if($lineCounter > 1)# && $lineCounter < 1000)
			{
				(my $chrom, my $position, my $allele1, my $allele2, my $reads2, my $avg_var_qual, my $num_var_strands) = split(/\t/, $line);
				my $position_key = $chrom . ":" . $position; #substr($position, 0, 2);
				my $reads_A = $CoverageByPosition{$position_key . ":A"};
				my $reads_C = $CoverageByPosition{$position_key . ":C"};
				my $reads_G = $CoverageByPosition{$position_key . ":G"};
				my $reads_T = $CoverageByPosition{$position_key . ":T"};
				my $reads_N = $CoverageByPosition{$position_key . ":N"};
				
				$reads_A = 0 if(!$reads_A);
				$reads_C = 0 if(!$reads_C);
				$reads_G = 0 if(!$reads_G);
				$reads_T = 0 if(!$reads_T);
				$reads_N = 0 if(!$reads_N);
				
				print "$chrom\t$position\t$allele1\t$reads_A\t$reads_C\t$reads_G\t$reads_T\t$reads_N\n";
				print OUTFILE "$chrom\t$position\t$allele1\t$reads_A\t$reads_C\t$reads_G\t$reads_T\t$reads_N\n";
#				print OUTFILE "$chrom\t$position\t$allele1\t$allele2\t";
#				print OUTFILE "$read_coverage\t$reads1\t$reads2\t$readsN\t";
#				print OUTFILE "$avg_ref_qual\t$avg_var_qual\t$num_ref_strands\t$num_var_strands\n";
				
			}
		}
	
		close($input);
	
		close(OUTFILE);


	}
	else
	{
		
		print "Parsing alignment blocks...\n";
	
		## Parse the alignment blocks ##
	
		$input = new FileHandle ($blocks_file);
		$lineCounter = 0;
		
		while (<$input>)
		{
			chomp;
			my $line = $_;
			$lineCounter++;
			
			if($lineCounter > 1)# && $lineCounter < 50000)
			{
				$GenotypeStats{'total_alignments'}++;
				print "$GenotypeStats{'total_alignments'} blocks parsed...\n" if($verbose && !($GenotypeStats{'total_alignments'} % 10000));
			
				my @lineContents = split(/\s+/, $line);
				my $chrom = $lineContents[0];
				my $chr_start = $lineContents[1];
				my $chr_stop = $lineContents[2];
				my $strand = $lineContents[3];
				my $read_name = $lineContents[4];
		
				
				## Get possible chrom position keys ##
				my $position_key;
				
				for(my $position = $chr_start; $position <= $chr_stop; $position++)
				{
					$position_key = $chrom . ":" . $position; #substr($position, 0, 2);
					if($SNPsByPosition{$position_key})
					{
						$CoverageByPosition{$position_key}++;
					}
				}
			
	
			}
		}
		
		close($input);
	
		my $keyCount = 0;
		foreach my $key (keys %CoverageByPosition)
		{
			$keyCount++;
		}
		
		print "$keyCount keys had coverage\n";
	
		## Open the outfile ##
		
		open(OUTFILE, ">$outfile") or die "Can't open outfile: $!\n";
		print OUTFILE "chrom\tposition\tref\tvar\tcov\treads1\treads2\tavgQual\tstrands\n";
	
		print "Parsing SNPs again...\n";
	
		$input = new FileHandle ($variants_file);
		$lineCounter = 0;
		while (<$input>)
		{
			chomp;
			my $line = $_;
			$lineCounter++;
			
			if($lineCounter == 1)
			{
			}
			if($lineCounter > 1)# && $lineCounter < 1000)
			{
				(my $chrom, my $position, my $allele1, my $allele2, my $num_reads, my $avg_qual, my $num_strands) = split(/\t/, $line);
				my $position_key = $chrom . ":" . $position; #substr($position, 0, 2);
	
				if($num_reads >= 10 && $avg_qual >= 20 && $CoverageByPosition{$position_key})
				{
					my $read_coverage = $CoverageByPosition{$position_key};
					## Calculate $reads1 ##
					my $num_wt_reads = $CoverageByPosition{$position_key} - $SNPsByPosition{$position_key};
					
					if($num_wt_reads < 0)
					{
	#                                        print "Warning: Reads1 calculated to be less than zero: $line\n";
	 #                                       exit(1);
					}
					else
					{
						print OUTFILE "$chrom\t$position\t$allele1\t$allele2\t$read_coverage\t$num_wt_reads\t$num_reads\t$avg_qual\t$num_strands\n";			
					}
				}
			}
		}
	
		close($input);
	
		close(OUTFILE);
	}

	
	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}



1;

