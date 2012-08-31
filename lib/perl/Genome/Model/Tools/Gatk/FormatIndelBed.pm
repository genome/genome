
package Genome::Model::Tools::Gatk::FormatIndelBed;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# FormatIndelBed - Merge glfSomatic/Varscan somatic calls in a file that can be converted to MAF format
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

class Genome::Model::Tools::Gatk::FormatIndelBed {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		variants_file	=> { is => 'Text', doc => "BED file of indel predictions", is_optional => 0, is_input => 1 },
		output_file     => { is => 'Text', doc => "Output file to receive formatted lines", is_optional => 0, is_input => 1, is_output => 1 },
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Formats indels in the GATK Bed file for the annotation pipeline"                 
}

sub help_synopsis {
    return <<EOS
This command formats indels for the annotation pipeline
EXAMPLE:	gmt gatk format-indel-bed --variants-file GATK.bed --output-file GATK.bed.formatted
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
	
	if(-e $variants_file)
	{
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
			my $numContents = @lineContents;
			my $chrom = $lineContents[0];
			
			if($numContents > 5)
			{
				die "Warning: Input file $variants_file has $numContents columns, which is more than expected for BED format!\n";				
			}
			elsif($numContents <= 5) # 10 for germline, 17 for somatic 
			{
				my $chr_start = $lineContents[1];
				my $chr_stop = $lineContents[2];
				(my $indel, my $readcounts) = split(/\:/, $lineContents[3]);

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


				## Append readcounts to rest of line ##	
				my $rest_of_line = $readcounts;

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

}


sub byChrPos
{
    (my $chrom_a, my $pos_a) = split(/\t/, $a);
    (my $chrom_b, my $pos_b) = split(/\t/, $b);

	$chrom_a =~ s/X/23/;
	$chrom_a =~ s/Y/24/;
	$chrom_a =~ s/MT/25/;
	$chrom_a =~ s/M/25/;
	$chrom_a =~ s/[^0-9]//g;

	$chrom_b =~ s/X/23/;
	$chrom_b =~ s/Y/24/;
	$chrom_b =~ s/MT/25/;
	$chrom_b =~ s/M/25/;
	$chrom_b =~ s/[^0-9]//g;

    $chrom_a <=> $chrom_b
    or
    $pos_a <=> $pos_b;
    
#    $chrom_a = 23 if($chrom_a =~ 'X');
#    $chrom_a = 24 if($chrom_a =~ 'Y');
    
}

#############################################################
# ParseBlocks - takes input file and parses it
#
#############################################################

sub fix_chrom
{
	my $chrom = shift(@_);
	$chrom =~ s/chr// if(substr($chrom, 0, 3) eq "chr");
	$chrom =~ s/[^0-9XYMNT\_random]//g;	

	return($chrom);
}



1;

