
package Genome::Model::Tools::Analysis::454::ParseNewblerVariants;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# ParseNewblerVariants - Load 454 reads from a sample-SFF tab-delimited file
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

class Genome::Model::Tools::Analysis::454::ParseNewblerVariants {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		diffs_file	=> { is => 'Text', doc => "The 454HCDiffs.txt or 454AllDiffs.txt file" },
                output_snps	=> { is => 'Text', doc => "Optional output file for SNPs", is_optional => 1 },
		output_indels	=> { is => 'Text', doc => "Optional output file for indels", is_optional => 1 },
                min_coverage	=> { is => 'Text', doc => "Minimum read depth to call variants", is_optional => 1 },
                min_var_freq	=> { is => 'Text', doc => "Minimum variant frequency to call variants", is_optional => 1 },
	],
};

#, is_optional => 1

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Run Newbler assembly of indel-supporting 454 reads"                 
}

sub help_synopsis {
    return <<EOS
This command retrieves the locations of unplaced reads for a given genome model
EXAMPLE:	gmt bowtie --query-file s_1_sequence.fastq --output-file s_1_sequence.Hs36.bowtie
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
	my $diffs_file = $self->diffs_file;
	my $output_snps = $self->output_snps if($self->output_snps);
	my $output_indels = $self->output_indels if($self->output_indels);
	my $min_coverage = 10;
	my $min_var_freq = 0.25;
	
	$min_coverage = $self->min_coverage if($self->min_coverage);
	$min_var_freq = $self->min_var_freq if($self->min_var_freq);

#	my $output_file = $self->output_file;

	if(!(-e $diffs_file))
	{
		die "Error: Diffs file not found!\n";
	}

	my %stats = ();
 
	## Open output files ##
	
	open(SNPS, ">$output_snps") or die "Can't open SNPs file: $!\n" if($output_snps);
	print SNPS "chrom\tposition\tref\tvar\tcoverage\treads1\treads2\n" if($output_snps);

	open(INDELS, ">$output_indels") or die "Can't open INDELs file: $!\n" if($output_indels);
	print INDELS "chrom\tchr_start\tchr_stop\tindel_type\tindel_size\tcoverage\treads1\treads2\tindel_in_context\n" if($output_indels);
 
 	my $input = new FileHandle ($diffs_file);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		
	
		if($line)
		{
			if(substr($line, 0, 1) eq ">" && !($line =~ 'Reference' || $line =~ 'Accno'))
			{
				(my $chrom, my $chr_start, my $chr_stop, my $allele1, my $allele2, my $read_depth, my $variant_freq) = split(/\t/, $line);
				$chrom = substr($chrom, 1, 99);
				$variant_freq =~ s/[^0-9\.]//g;
				$variant_freq = $variant_freq / 100;
				
				if($read_depth >= $min_coverage)
				{
					if($variant_freq >= $min_var_freq)
					{
						## Determine read counts ##
						
						my $reads2 = sprintf("%d", ($read_depth * $variant_freq));
						my $reads1 = $read_depth - $reads2;
						$reads1 = 0 if($reads1 < 0);
						
						## Determine variant type ##
						
						if($allele1 =~ '-' || $allele2 =~ '-')
						{
							my $indel_type = my $indel_size = "";
							## Determine indel type and size ##
							
							if($allele1 =~ '-')
							{
								$indel_type = 'INSERTION';
								$indel_size = length($allele2);
							}
							else
							{
								$indel_type = 'DELETION';
								$indel_size = length($allele1);								
							}

							print INDELS "$chrom\t$chr_start\t$chr_stop\t$indel_type\t$indel_size\t$read_depth\t$reads1\t$reads2\t[$allele1/$allele2]\n" if($output_indels);

							$stats{'num_indels'}++;
						}
						elsif(length($allele1) == 1 && length($allele2) == 1 && $allele1 ne 'N' && $allele2 ne 'N')
						{
							print SNPS "$chrom\t$chr_start\t$allele1\t$allele2\t$read_depth\t$reads1\t$reads2\n" if($output_snps);
							$stats{'num_snps'}++;							
						}
					}
				}
			}
		}
	}
	
	close($input);

	$stats{'num_snps'} = 0 if(!$stats{'num_snps'});
	$stats{'num_indels'} = 0 if(!$stats{'num_indels'});

	print "$stats{'num_snps'} SNPs parsed from file\n";
	print "$stats{'num_indels'} indels parsed from file\n"; 

 
	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}



1;

