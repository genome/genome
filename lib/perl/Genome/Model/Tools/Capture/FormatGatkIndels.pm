
package Genome::Model::Tools::Capture::FormatGatkIndels;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# FormatGatkIndelsForAnnotation - Merge glfSomatic/Varscan somatic calls in a file that can be converted to MAF format
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
use Genome::Model::Tools::Capture::Helpers qw(
    byChrPos
    fix_chrom
    iupac_to_base
);

class Genome::Model::Tools::Capture::FormatGatkIndels {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		variants_file	=> { is => 'Text', doc => "File of indel predictions", is_optional => 0, is_input => 1 },
		output_file     => { is => 'Text', doc => "Output file to receive formatted lines", is_optional => 0, is_input => 1, is_output => 1 },
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

	my @formatted = ();
	my $formatCounter = 0;
	
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
		my $chr_start = $lineContents[1];
		my $chr_stop = $lineContents[2];
		my $indel = $lineContents[3];

		if(substr(uc($chrom), 0, 5) eq "CHROM" || substr(uc($chrom), 0, 3) eq "REF")
		{
			## Ignore header lines ##	
			
		}
		else
		{
			$chrom = fix_chrom($chrom);
			
			my $ref = my $var = "";
			my $indel_type = my $indel_size = "";

			my $restColumn = 4;

			if(substr($indel, 0, 1) eq '+')
			{
				$indel_type = "INSERTION";
				$indel_size = length($indel) - 1;
				$chr_stop = $chr_start + 1;
				$ref = "0";
				$var = substr($indel, 1, 99);
			}
			else
			{
				$indel_type = "DELETION";
				$indel_size = length($indel) - 1;
				$chr_start++;
				$ref = substr($indel, 1, 99);
				$var = "0";
			}

			## If we have other information on line, output it ##
			my $numContents = @lineContents;
			my $rest_of_line = "";
			if($restColumn && $restColumn > 0 && $restColumn < $numContents)
			{
				for(my $colCounter = $restColumn; $colCounter < $numContents; $colCounter++)
				{
					$rest_of_line .= "\t" if($colCounter > $restColumn);
					$rest_of_line .= $lineContents[$colCounter];
				}
			}

			## Fix chromosome ##
			
			$chrom =~ s/chr//;
			$chrom = "MT" if($chrom eq "M");
			$chrom = "" if($chrom =~ 'random');

			## If we have the necessary information, output line ##
			
			if($chrom && $chr_start && $chr_stop && ($ref || $var))
			{
				$formatted[$formatCounter] = "$chrom\t$chr_start\t$chr_stop\t$ref\t$var\t$rest_of_line";
				$formatCounter++;

#				print OUTFILE "$chrom\t$chr_start\t$chr_stop\t$allele1\t$allele2";
#				print OUTFILE "\t$rest_of_line" if($rest_of_line);
#				print OUTFILE "\n";
			}

		}
	}

	close($input);	
	
	
	## Sort the formatted indels by chr pos ##

	@formatted = sort byChrPos @formatted;
	
	foreach my $snv (@formatted)
	{
		print OUTFILE "$snv\n";
	}
		
	
	close(OUTFILE);
}

1;
