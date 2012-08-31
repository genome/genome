
package Genome::Model::Tools::Analysis::SomaticPipeline::FilterGlfIndels;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# FilterGlfIndels - Merge glfSomatic/Varscan somatic calls in a file that can be converted to MAF format
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@watson.wustl.edu)
#
#	CREATED:	10/23/2009 by D.K.
#	MODIFIED:	10/23/2009 by D.K.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;
use FileHandle;
use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

class Genome::Model::Tools::Analysis::SomaticPipeline::FilterGlfIndels {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		variants_file	=> { is => 'Text', doc => "File of variants in indel format", is_optional => 0 },
		output_file     => { is => 'Text', doc => "Output file to receive filtered indels", is_optional => 0 },
		min_coverage     => { is => 'Text', doc => "Minimum coverage to allow indel [8]", is_optional => 1 },
		min_reads2     => { is => 'Text', doc => "Minimum read support to allow indel [2]", is_optional => 1 },
		min_var_freq     => { is => 'Text', doc => "Minimum allele frequency to allow indel [0.10]", is_optional => 1 },
		min_somatic_score     => { is => 'Text', doc => "Minimum somatic score to allow indel [1]", is_optional => 1 },
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Filters glfSomatic indel calls"                 
}

sub help_synopsis {
    return <<EOS
This command merges variant calls from the pipeline with their annotation information
EXAMPLE:	gmt analysis somatic-pipeline merge-snvs-with-annotation --variants-file [file] --annotation-file [file] --output-file [file]
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
	my $output_file = $self->output_file;
	my $min_coverage = 8;
	$min_coverage = $self->min_coverage if($self->min_coverage);

	my $min_reads2 = 2;
	$min_reads2 = $self->min_reads2 if($self->min_reads2);

	my $min_var_freq = 0.10;
	$min_var_freq = $self->min_var_freq if($self->min_var_freq);

	my $min_somatic_score = 1;
	$min_somatic_score = $self->min_somatic_score if($self->min_somatic_score);

	my %stats = ();
	$stats{'num_indels'} = $stats{'num_pass_filter'} = 0;
	
	## Open outfile ##
	
	open(OUTFILE, ">$output_file") or die "Can't open outfile: $!\n";
	
	## Parse the variants file ##
	
	my $input = new FileHandle ($variants_file);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;

		my @fields = split(/\t/, $line);

		my $chrom    = $fields[0];
		my $chr_start    = $fields[1];
		my $chr_stop = $fields[2];
		my $ref_allele = $fields[3];
		my $var_allele = $fields[4];
		my $indel_type = $fields[5];
		my $somatic_score = $fields[6];

		my $genotype_allele = "";
		
		if($indel_type =~ 'DEL')
		{
			$genotype_allele = '-' . $ref_allele;
		}
		elsif($indel_type =~ 'INS')
		{
			$genotype_allele = '+' . $var_allele;
		}

		#$my $tumor_coverage = $fields[12];
		my $tumor_reads1 = $fields[12];
		my $tumor_reads2 = $fields[13];
		my $normal_reads1 = $fields[25]; # number of reads that support indel 1 in normal
		my $normal_reads2 = $fields[26]; # number of reads that support indel 2 in normal 

		my $tumor_coverage = $tumor_reads1 + $tumor_reads2;
		my $normal_coverage = $normal_reads1 + $normal_reads2;		

		$stats{'num_indels'}++;

		if($tumor_coverage >= $min_coverage && $normal_coverage >= $min_coverage && $tumor_reads2 >= $min_reads2 && $somatic_score >= $min_somatic_score)
		{
			my $tumor_freq = ($tumor_reads2) / ($tumor_reads1 + $tumor_reads2);

			if($tumor_freq >= $min_var_freq)
			{
				## call genotype ##
				my $genotype = "";
				
				if($tumor_freq > 0.80)
				{
					$genotype = "$genotype_allele/$genotype_allele";
				}
				else
				{
					$genotype = "*/$genotype_allele";
				}
				
				$stats{'num_pass_filter'}++;

				$tumor_freq = sprintf("%.2f", ($tumor_freq * 100));

				print OUTFILE "$chrom\t$chr_start\t$chr_stop\t*\t$genotype\t$somatic_score\t$tumor_coverage\t$tumor_reads1\t$tumor_reads2\t$tumor_freq\%\t$normal_coverage\t$normal_reads1\t$normal_reads2\n";	
# if the next line ever become uncommented, please change ARGV usage to a property with shell_args_position
#				print OUTFILE "$line\n" if($ARGV[1]);
#				print "$chr\t$pos\t$indel1\t$indel2\t$tumor_reads1\t$tumor_reads2\t$tumor_freq\n";
			}
			
#			print "$tumor_reads1\t$tumor_reads2\t$tumor_freq\n";
		}
		elsif($normal_coverage < $min_coverage)
		{
#			warn "Normal coverage only $normal_coverage ($normal_reads1 + $normal_reads2)\n";
		}
	}

	close($input);
		
	close(OUTFILE);
	
	print "$stats{'num_indels'} indels\n";
	print "$stats{'num_pass_filter'} passed filter\n";
}





1;

