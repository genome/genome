
package Genome::Model::Tools::Gatk::FormatVcf;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# FormatVcf - "Formats GATK UnifiedGenotyper indel detection vcf files"
#					
#	AUTHOR:		Will Schierding
#
#	CREATED:	03-Mar-2011 by W.S.
#	MODIFIED:	03-Mar-2011 by W.S.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;
use FileHandle;
use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

class Genome::Model::Tools::Gatk::FormatVcf {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		variants_file	=> { is => 'Text', doc => "File of indel predictions", is_optional => 0, is_input => 1 },
		output_file     => { is => 'Text', doc => "Output file to receive formatted lines", is_optional => 0, is_input => 1, is_output => 1 },
	],
};

sub sub_command_sort_position { 1 }

sub help_brief {                            # keep this to just a few words <---
    "Formats GATK UnifiedGenotyper indel detection vcf files"                 
}

sub help_synopsis {
    return <<EOS
This command runs the GATK indel detection pipeline
EXAMPLE:	gmt gatk format-vcf --variants-file input.vcf --output-file variant_format_for_pipeline.txt
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
		
		while (my $line = <$input>)
		{
			chomp($line);
			if ($line =~ m/^#/) {
				## Ignore header lines ##	
				next;
			}
			my @lineContents = split(/\t/, $line);
			my $numContents = @lineContents;
			my $chrom = $lineContents[0];
			my $pos = $lineContents[1] + 1;
			my $ref1 = $lineContents[3];
			my $alt = $lineContents[4];
			$chrom = fix_chrom($chrom);
				
			my ($ref, $var) = "";
			my ($indel_type, $indel_size, $chr_start, $chr_stop) = "";
#1	735233	.	TA	T	190.34	.	AC=1;AF=0.50;AN=2;DP=38;Dels=0.24;HRun=2;HaplotypeScore=2.4689;MQ=38.75;MQ0=0;QD=5.01;SB=-0.00;sumGLbyD=6.04	GT:AD:DP:GQ:PL	0/1:29,9:29:99:229,0,771
############SUBTRACT SO THAT ONE SIDE IS 0 AND THE OTHER NO LONGER HAS REF BASE AS START##############################	
			if(length($ref1) == 1) {
				$indel_type = "INSERTION";
				$ref = "0";
				$var = substr($alt, 1);
				$indel_size = length($var);
				$chr_start = $pos;
				$chr_stop = $chr_start + 1;


			}
			elsif(length($alt) == 1) {
				$indel_type = "DELETION";
				$ref = substr($ref1, 1);
				$var = "0";
				$indel_size = length($ref);
				$chr_start = $pos;
				$chr_stop = $chr_start + $indel_size - 1;
			}
			else {
				die "Cannot resolve insertion or deletion from ref:$ref alt:$alt\n";
			}
			## If we have other information on line, output it ##
			my $rest_of_line = "";
			my $restColumn = 5;
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
	else {
		die "Variants File Doesn't Exist";
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

