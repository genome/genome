
package Genome::Model::Tools::Analysis::SomaticPipeline::FormatIndelsForAnnotation;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# FormatIndelsForAnnotation - Merge glfSomatic/Varscan somatic calls in a file that can be converted to MAF format
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

class Genome::Model::Tools::Analysis::SomaticPipeline::FormatIndelsForAnnotation {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		variants_file	=> { is => 'Text', doc => "File of indel predictions", is_optional => 0 },
		output_file     => { is => 'Text', doc => "Output file to receive formatted lines", is_optional => 0 },
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Formats indels for the annotation pipeline"                 
}

sub help_synopsis {
    return <<EOS
This command formats indels for the annotation pipeline
EXAMPLE:	gmt analysis somatic-pipeline format-indels-for-annotation --variants-file [file] --output-file [file]
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
	
	## Open outfile ##
	
	open(OUTFILE, ">$output_file") or die "Can't open outfile: $!\n";
	
	## Parse the indels ##

	my $input = new FileHandle ($variants_file);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;

		my @lineContents = split(/\t/, $line);
		my $chrom = $lineContents[0];
	
		if(substr(uc($chrom), 0, 5) eq "CHROM" || substr(uc($chrom), 0, 3) eq "REF")
		{
			## Ignore header lines ##	
		}
		else
		{
			my $chr_start = $lineContents[1];
			my $chr_stop = my $allele1 = my $allele2 = "";
			my $ref = my $var = "";
			my $indel_type = my $indel_size = "";

			my $restColumn = 0;

			if($lineContents[2] =~ /[0-9]/)
			{
				$chr_stop = $lineContents[2];
				$ref = $lineContents[3];
				$var = $lineContents[4];
				$restColumn = 5;
			}
			else
			{
				$ref = $lineContents[2];
				$var = $lineContents[3];
				$restColumn = 4;
			}

			## Correct alleles ##

			if($ref eq '-' || $var eq '-')
			{
				$allele1 = $ref;
				$allele2 = $var;
				
				if($ref eq '-')
				{
					$indel_type = "INSERTION";
					$indel_size = length($var);
				}
				else
				{
					$indel_type = "DELETION";
					$indel_size = length($ref);
				}
			}
			elsif(substr($var, 0, 1) eq '+')
			{
				$allele1 = "-";
				$allele2 = uc($var);
				$allele2 =~ s/[^ACGTN]//g;
				$indel_type = "INSERTION";
				$indel_size = length($allele2);
			}
			elsif(substr($var, 0, 1) eq '-')
			{
				$allele2 = "-";
				$allele1 = uc($var);
				$allele1 =~ s/[^ACGTN]//g;
				$indel_type = "DELETION";
				$indel_size = length($allele1);
			}
			else
			{
				warn "Unable to format $line\n";
				$chrom = $chr_start = $chr_stop = $allele1 = $allele2 = "";
			}

			## If no chr stop, calculate it ##
			if(!$chr_stop)
			{
				if($indel_type eq "INSERTION")# || $indel_size == 1
				{
					$chr_stop = $chr_start + 1;
				}
				else
				{
					$chr_stop = $chr_start + $indel_size - 1;
				}
			}

			## If we have other information on line, output it ##
			my $numContents = @lineContents;
			my $rest_of_line = "";
			if($restColumn && $restColumn > 0 && $restColumn < $numContents)
			{
				for(my $colCounter = $restColumn; $colCounter < $numContents; $colCounter++)
				{
					$rest_of_line .= "\t" if($rest_of_line);
					$rest_of_line .= $lineContents[$colCounter];
				}

			}

			## If we have the necessary information, output line ##
			$chrom =~ s/chr// if(substr($chrom, 0, 3) eq "chr");
			$chrom = "" if($chrom =~ 'random');
			
			if($chrom && $chr_start && $chr_stop && $allele1 && $allele2)
			{
				print OUTFILE "$chrom\t$chr_start\t$chr_stop\t$allele1\t$allele2";
				print OUTFILE "\t$rest_of_line" if($rest_of_line);
				print OUTFILE "\n";
			}

		}
	}

	close($input);	
	
	close(OUTFILE);
}





1;

