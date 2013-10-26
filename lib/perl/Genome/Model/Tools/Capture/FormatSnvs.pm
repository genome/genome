
package Genome::Model::Tools::Capture::FormatSnvs;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# FormatSnvsForAnnotation - Merge glfSomatic/Varscan somatic calls in a file that can be converted to MAF format
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
use Genome::Model::Tools::Capture::Helpers qw/iupac_to_base fix_chrom sortByChrPos/;

class Genome::Model::Tools::Capture::FormatSnvs {
	is => 'Genome::Model::Tools::Capture',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		variants_file	=> { is => 'Text', doc => "File of SNV predictions", is_optional => 0, is_input => 1 },
		output_file     => { is => 'Text', doc => "Output file to receive formatted lines", is_optional => 0, is_input => 1, is_output => 1 },
		preserve_call	=> { is => 'Text', doc => "If set to 1, preserves the consensus call", is_optional => 1, is_input => 1 },
		append_line	=> { is => 'Text', doc => "If set to 1, appends extra columns in input lines to output lines", is_optional => 1, is_input => 1, default => 0 },
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Formats SNVs for the annotation pipeline"                 
}

sub help_synopsis {
    return <<EOS
This command formats indels for the annotation pipeline
EXAMPLE:	gmt analysis somatic-pipeline format-snvs-for-annotation --variants-file [file] --output-file [file]
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

	my @formatted = ();
	my $formatCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		

		my @lineContents = split(/\t/, $line);

		if(!(lc($lineContents[0]) =~ "chrom" || lc($lineContents[0]) =~ "ref_name"))
		{
			my $chrom = $lineContents[0];
			$chrom = fix_chrom($chrom);
			my $chr_start = $lineContents[1];
			my $chr_stop = my $allele1  = my $allele2 = "";

			my $restColumn = 0;

			if($lineContents[2] && $lineContents[2] =~ /[0-9]/)
			{
				$chr_stop = $lineContents[2];
				$allele1 = $lineContents[3];
				$allele2 = $lineContents[4];
				$restColumn = 5;
			}
			else
			{
				$chr_stop = $chr_start;
				$allele1 = $lineContents[2];
				$allele2 = $lineContents[3];
				$restColumn = 4;
			}

			if($chrom && $chr_start && $chr_stop)
			{
				$allele2 = iupac_to_base($allele1, $allele2) if(!$self->preserve_call);
	
				## If we have other information on line, output it ##
				my $numContents = @lineContents;
				my $rest_of_line = "";
				if($self->append_line && $restColumn && $restColumn > 0 && $restColumn < $numContents)
				{
					for(my $colCounter = $restColumn; $colCounter < $numContents; $colCounter++)
					{
						$rest_of_line .= "\t" if($colCounter > $restColumn);
						$rest_of_line .= $lineContents[$colCounter];
					}
	
				}
	
				$formatted[$formatCounter] = "$chrom\t$chr_start\t$chr_stop\t$allele1\t$allele2\t$rest_of_line";
				$formatCounter++;
			}
		}
	}

	close($input);
	
	## Sort the formatted indels by chr pos ##

	@formatted = sortByChrPos(@formatted);
	
	foreach my $snv (@formatted)
	{
		print OUTFILE "$snv\n";
	}
	
	
	close(OUTFILE);
}

1;
