
package Genome::Model::Tools::Bowtie::CallGenotypes;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# CallGenotypes.pm - 	Limit SNPs by various factors, like ROI positions.
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@watson.wustl.edu)
#
#	CREATED:	03/20/2009 by D.K.
#	MODIFIED:	03/20/2009 by D.K.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;

use FileHandle;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

class Genome::Model::Tools::Bowtie::CallGenotypes {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		variants_file	=> { is => 'Text', doc => "File containing SNPs" },
		output_file	=> { is => 'Text', doc => "Output file to receive limited SNPs" },
		verbose	=> { is => 'Text', doc => "Print interim progress updates", is_optional => 1 },
		min_coverage	=> { is => 'Text', doc => "Minimum read coverage to include SNP", is_optional => 1 },
		min_reads1	=> { is => 'Text', doc => "Minimum ref reads to include SNP", is_optional => 1 },
		min_reads2	=> { is => 'Text', doc => "Minimum var reads to include SNP", is_optional => 1 },
		min_strands2	=> { is => 'Text', doc => "Minimum var strands to include SNP", is_optional => 1 },
		min_avg_qual	=> { is => 'Text', doc => "Minimum average base quality to include SNP", is_optional => 1 },
		min_var_freq	=> { is => 'Text', doc => "Minimum variant frequency to include SNP", is_optional => 1 },
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Choose the most likely SNP call from readcounts"                 
}

sub help_synopsis {
    return <<EOS
This command parses alignments from Bowtie.
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
	my $outfile = $self->output_file;
	my $verbose = 1 if(defined($self->verbose));
	my $min_coverage = $self->min_coverage if(defined($self->min_coverage));
	my $min_reads1 = $self->min_reads1 if(defined($self->min_reads1));
	my $min_reads2 = $self->min_reads2 if(defined($self->min_reads2));
	my $min_strands2 = $self->min_strands2 if(defined($self->min_strands2));
	my $min_avg_qual = $self->min_avg_qual if(defined($self->min_avg_qual));
	my $min_var_freq = $self->min_var_freq if(defined($self->min_var_freq));

	my %stats = ();
	$stats{'total_snps'} = $stats{'limit_snps'} = 0;
	my $input, my $lineCounter;


	print "Parsing SNPs...\n";
	## Parse the combined SNPs ##

	my %variants_by_position = ();

	## Open the outfile ##
	open(OUTFILE, ">$outfile") or die "Can't open outfile: $!\n";
	
	$input = new FileHandle ($variants_file);
	$lineCounter = 0;
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		if($lineCounter == 1)
		{
			print OUTFILE "$line\n";
		}
		if($lineCounter > 1)# && $lineCounter < 1000)
		{
			$stats{'total_snps'}++;
#			(my $chrom, my $position, my $allele1, my $allele2, my $context, my $num_reads, my $reads) = split(/\t/, $line);			
			(my $chrom, my $position, my $allele1, my $allele2, my $coverage, my $reads1, my $reads2, my $readsN, my $avg_ref_qual, my $avg_var_qual, my $num_ref_strands, my $num_var_strands) = split(/\t/, $line);

			my $variant_freq = $reads2 / $coverage if($coverage);

			if(!$min_coverage || $coverage >= $min_coverage)
			{
				if(!$min_reads2 || $reads2 >= $min_reads2)
				{
					if(!$min_var_freq || $variant_freq >= $min_var_freq)
					{
						if(!$min_avg_qual || $avg_var_qual >= $min_avg_qual)
						{
							if(!$min_strands2 || $num_var_strands >= $min_strands2)
							{
								$stats{'limit_snps'}++;
								$variants_by_position{$chrom . "\t" . $position} .= "$line\n";
							}
						}
					}
				}
			}
		
			
		}
	}

	close($input);

	## Get percentage ##
	
	$stats{'limit_pct'} = sprintf("%.2f", ($stats{'limit_snps'} / $stats{'total_snps'} * 100)) . '%' if($stats{'total_snps'});
	print "$stats{'total_snps'} total SNPs\n";
	print "$stats{'limit_snps'} SNPs ($stats{'limit_pct'}) remained after filter\n";

	
	## Go through and sort positions ##
	print "Calling bases...\n";

	foreach my $chrom_key (keys %variants_by_position)
	{
		my @position_readcounts = split(/\n/, $variants_by_position{$chrom_key});
		@position_readcounts = sort byReadCount @position_readcounts;
		print OUTFILE $position_readcounts[0] . "\n";
	}

	sub byReadCount
	{
		(my $chrom_a, my $position_a, my $allele1_a, my $allele2_a, my $coverage_a, my $reads1_a, my $reads2_a, my $readsN_a, my $avg_ref_qual_a, my $avg_var_qual_a, my $num_ref_strands_a, my $num_var_strands_a) = split(/\t/, $a);
		(my $chrom_b, my $position_b, my $allele1_b, my $allele2_b, my $coverage_b, my $reads1_b, my $reads2_b, my $readsN_b, my $avg_ref_qual_b, my $avg_var_qual_b, my $num_ref_strands_b, my $num_var_strands_b) = split(/\t/, $b);

		## SORT ORDER: reads2 DESC, strands2 DESC, avg_var_qual DESC ##
		$reads2_b <=> $reads2_a
		or
		$num_var_strands_b <=> $num_var_strands_a
		or
		$avg_var_qual_b <=> $avg_var_qual_a;
	}
	
	close(OUTFILE);

	
	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}



1;

