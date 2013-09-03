
package Genome::Model::Tools::Capture::ExportVcf;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# ExportVcf - Export mutation calls in VCF format
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@watson.wustl.edu)
#
#	CREATED:	11/17/2009 by D.K.
#	MODIFIED:	11/17/2009 by D.K.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;
use FileHandle;
use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

class Genome::Model::Tools::Capture::ExportVcf {
	is => 'Genome::Model::Tools::Capture',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		variants_file	=> { is => 'Text', doc => "File of predicted mutations", is_optional => 0, is_input => 1 },
		annotation_file	=> { is => 'Text', doc => "Annotation file(s) for the variants", is_optional => 0, is_input => 1 },
		output_file     => { is => 'Text', doc => "Outfile for annotation plus four columns", is_optional => 0, is_input => 1, is_output => 1 },
		vcf_file     => { is => 'Text', doc => "Outfile for final VCF records", is_optional => 0, is_input => 1, is_output => 1 },
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
	my $annotation_file = $self->annotation_file;
	my $output_file = $self->output_file;
	my $vcf_file = $self->vcf_file;
	
	my %stats = ();

	## Open outfile ##
	
	open(OUTFILE, ">$output_file") or die "Can't open outfile: $!\n";

	## Load the annotation ##
	
	my %annotation = load_annotation($annotation_file);
	
	## Parse the variants ##

	my $input = new FileHandle ($variants_file);
	my $lineCounter = 0;

	my @formatted = ();
	my $formatCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		

		(my $chromosome, my $chr_start, my $chr_stop, my $ref, my $cns) = split(/\t/, $line);
		my @lineContents = split(/\t/, $line);
		my $numContents = @lineContents;
		
		my $var = $self->iupac_to_base($ref, $cns);
		
		$stats{'num_variants'}++;
		$ref = "-" if($ref eq "0");
		$var = "-" if($var eq "0");
		my $key = join("\t", $chromosome, $chr_start, $chr_stop, $ref, $var);
		
		## PRoceed if we have annotation ##

		if($annotation{$key})
		{
			$stats{'num_with_annotation'}++;
			## The four columns we need to append to annotation are:
			# Score
			# Genotype
			# Tier 1, 2, 3, 4
			# Confidence
			
			my $score = my $genotype = my $tier = my $confidence = "";
			$confidence = "High";
			$tier = 1;
			
			if($numContents == 12 || $lineContents[5] =~ 'INS' || $lineContents[5] =~ 'DEL')
			{
				## Somatic Sniper only call ##
				$confidence = "Low";
				
				$score = $lineContents[6];
				
				if($cns eq "A" || $cns eq "C" || $cns eq "G" || $cns eq "T")
				{
					$genotype = $var . $var;				
				}
				elsif($cns =~ '\/')
				{
					$genotype = $cns;
				}
				else
				{
					$var = $self->iupac_to_base($ref, $cns);
					$genotype = $ref . $var;
				}
			}
			else
			{
				## Assume Varscan File ##
	
				my $tumor_cns = $lineContents[12];
				if($tumor_cns eq "A" || $tumor_cns eq "C" || $tumor_cns eq "G" || $tumor_cns eq "T")
				{
					$genotype = $var . $var;
				}
				elsif($tumor_cns =~ '\/')
				{
					$genotype = $tumor_cns;
				}
				else
				{
					$genotype = $ref . $var;
				}
				
				my $p_value = $lineContents[15];
				
				## Convert P-value to confidence score ##
				my $conf_score = 0;
				
				if($p_value == 0)
				{
					$conf_score = 255;
				}
				else
				{
					$conf_score = -10 * (log($p_value) / log(10));
					$conf_score = sprintf("%.2f", $conf_score);
					$conf_score = 255 if($conf_score > 255);					
				}
	
				$score = $conf_score;

			}
			
			$formatted[$formatCounter] = join("\t", $annotation{$key}, $score, $genotype, $tier, $confidence);
			$formatCounter++;			
		}
		else
		{
			warn "No annotation for $key\n";
		}
		

	}

	close($input);
	
	## Sort the formatted indels by chr pos ##

	@formatted = $self->sortByChrPos(@formatted);
	
	foreach my $snv (@formatted)
	{
		print OUTFILE "$snv\n";
	}

	
	
	close(OUTFILE);

	print $stats{'num_variants'} . " variants in file\n";
	print $stats{'num_with_annotation'} . " had annotation\n";

	## Run the VCF Conversion ##
	
	system("gmt tcga convert-mf-to-vcf --mf-file $output_file --output-file $vcf_file");
}





################################################################################################
# Execute - the main program logic
#
################################################################################################

sub load_annotation
{
	my $variant_file = shift(@_);
	my $input = new FileHandle ($variant_file);
	my $lineCounter = 0;
	
	my %annotation = ();
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		(my $chromosome, my $chr_start, my $chr_stop, my $ref, my $var) = split(/\t/, $line);
		$ref = "-" if($ref eq "0");
		$var = "-" if($var eq "0");
		
		my $key = join("\t", $chromosome, $chr_start, $chr_stop, $ref, $var);
		
		$annotation{$key} = $line;
	}
	
	close($input);	
	
	return(%annotation);
}

1;

