
package Genome::Model::Tools::Analysis::454::VerifyMaf;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# VerifyMaf - Align reads with SSAHA2 or other aligner
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

my %file1_results = my %file2_results = ();

class Genome::Model::Tools::Analysis::454::VerifyMaf {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		maf_file	=> { is => 'Text', doc => "Original MAF file" },
		annotation_file		=> { is => 'Text', doc => "Annotations that might be missing", is_optional => 1 },
		output_file		=> { is => 'Text', doc => "Output file for updated MAF" },
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Fixes inconsistent and/or missing information in maf file "                 
}

sub help_synopsis {
    return <<EOS
This command updates variants in a MAF file with 454 validation status
EXAMPLE:	gt analysis 454 verify-maf --maf-file original.maf --output-file fixed.maf
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
	my $maf_file = $self->maf_file;
	my $output_file = $self->output_file;

	if(!(-e $maf_file))
	{
		die "Error: MAF file 1 not found!\n";
	}

	my %stats = ();

	my %annotations = ();
	
	if($self->annotation_file)
	{
		%annotations = load_annotations($self->annotation_file);
	}


	open(OUTFILE, ">$output_file") or die "Can't open outfile: $!\n";

	my @output_lines = ();
	my $num_lines = 0;

	## Parse the MAF file ##
	
	my $input = new FileHandle ($maf_file);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		
	
		if($lineCounter == 1 && $line =~ "Chrom")
		{
			## Print the MAF header ##
			print OUTFILE "$line\n";
		}
		elsif($lineCounter == 1)
		{
			warn "Warning: No Header in MAF file!\n";
		}
		else
		{
			my @lineContents = split(/\t/, $line);			
			my $chromosome = $lineContents[4];
			my $chr_start = $lineContents[5];
			my $chr_stop = $lineContents[6];
			my $trv_type = $lineContents[8];
			my $var_type = $lineContents[9];
			my $ref_base = $lineContents[10];
			my $var_base = $lineContents[11];
			my $tumor_sample = $lineContents[15];
			my $normal_sample = $lineContents[16];
			my $current_val_status = $lineContents[24];
			my $val_method = $lineContents[28];
			my $tumor_val_allele1 = $lineContents[19];
			my $tumor_val_allele2 = $lineContents[20];
			my $normal_val_allele1 = $lineContents[21];
			my $normal_val_allele2 = $lineContents[22];
			my $source = $lineContents[27];
			my $sequencer = $lineContents[31];

			$var_base = $lineContents[12] if($lineContents[12] ne $ref_base);

			## Fix ambiguity codes in alleles ##
			
			if($var_type eq "SNP")
			{
				$tumor_val_allele1 = check_allele($ref_base, $var_base, $tumor_val_allele1, $tumor_val_allele2);
				$tumor_val_allele2 = check_allele($ref_base, $var_base, $tumor_val_allele2, $tumor_val_allele1);

				$normal_val_allele1 = check_allele($ref_base, $ref_base, $normal_val_allele1, $normal_val_allele2);
				$normal_val_allele2 = check_allele($ref_base, $ref_base, $normal_val_allele2, $normal_val_allele1);

				$lineContents[19] = $tumor_val_allele1;
				$lineContents[20] = $tumor_val_allele2;
				$lineContents[21] = $normal_val_allele1;
				$lineContents[22] = $normal_val_allele2;				
			}

			## Fix unknown validation status type ##
			
			if($current_val_status ne "unknown")
			{
				$val_method = check_platform($val_method);
				$lineContents[28] = $val_method;
				
			}

			## Check source ##
			
			if($source =~ "3730" || $source =~ "PCR")
			{
				$source = "PCR";
			}
			else
			{
				$source = "Capture";
			}
			$lineContents[27] = $source;


			## Check sequencer ##
			
			$sequencer = check_platform($sequencer);
			$lineContents[31] = $sequencer;

			
			## Fix unknown TRV type ##
			
			$trv_type = check_trv_type($trv_type);
			if($trv_type eq "unknown")
			{
				my $key = "$chromosome\t$chr_start\t$chr_stop\t$var_type";
				if($annotations{$key})
				{
					$trv_type = check_trv_type($annotations{$key});

					print "Updated to $trv_type\n";
				}
				
			}
			
			$lineContents[8] = $trv_type;
			print "$chromosome\t$chr_start\t$chr_stop\t$ref_base\t$var_base\t$var_type\n" if($trv_type =~ "unknown");

			## Compile a new line ##

			my $numContents = @lineContents;
			my $newline = "";
			for(my $colCounter = 0; $colCounter < $numContents; $colCounter++)
			{
				$newline .= "\t" if($newline);
				$newline .= "$lineContents[$colCounter]";
			}
			
			print OUTFILE "$newline\n";
			
		}

	}

	close($input);	
		
	close(OUTFILE);



	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}





#############################################################
# parse_file - takes input file and parses it
#
#############################################################

sub load_annotations
{
	my $annotation_file = shift(@_);
	
	my %annotations = ();
	
	my $input = new FileHandle ($annotation_file);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;	

		(my $chrom, my $chr_start, my $chr_stop, my $ref, my $var, my $var_type) = split(/\t/, $line);
		my @lineContents = split(/\t/, $line);
		my $annotation = $lineContents[13];
		my $key = "$chrom\t$chr_start\t$chr_stop\t$var_type";
		$annotations{$key} = $annotation;
	}
	
	close($input);
	
	return(%annotations);
}


#############################################################
# parse_file - takes input file and parses it
#
#############################################################

sub check_trv_type
{
	(my $trv_type) = @_;

	my $orig = $trv_type;

	# Valid values;
	#Frame_Shift_Del, Frame_Shift_Ins, In_Frame_Del, In_Frame_Ins, Missense_Mutation, Nonsense_Mutation, Silent, Splice_Site_Indel, Splice_Site_SNP, Targeted_Region 

	$trv_type = lc($trv_type);
	
	if($trv_type =~ 'missense')
	{
		return("Missense_Mutation");
	}
	elsif($trv_type =~ 'nonsense' || $trv_type =~ 'nonstop')
	{
		return("Nonsense_Mutation");
	}
	elsif($trv_type =~ 'silent')
	{
		return("Silent");
	}	
	elsif($trv_type =~ 'frame' && $trv_type =~ 'shift' && $trv_type =~ 'del')
	{
		return("Frame_Shift_Del");
	}
	elsif($trv_type =~ 'frame' && $trv_type =~ 'shift' && $trv_type =~ 'ins')
	{
		return("Frame_Shift_Ins");
	}
	elsif($trv_type =~ 'in' && $trv_type =~ 'frame' && $trv_type =~ 'del')
	{
		return("In_Frame_Del");
	}
	elsif($trv_type =~ 'in' && $trv_type =~ 'frame' && $trv_type =~ 'ins')
	{
		return("In_Frame_Ins");
	}
	elsif($trv_type =~ 'splice' && ($trv_type =~ "ins" || $trv_type =~ "del"))
	{
		return("Splice_Site_Indel");
	}
	elsif($trv_type =~ 'splice')
	{
		return("Splice_Site_SNP");
	}	
	elsif($trv_type =~ 'region')
	{
		return("Target_Region");
	}

	return($orig);
}




#############################################################
# parse_file - takes input file and parses it
#
#############################################################

sub check_allele
{
	(my $ref_base, my $var_base, my $this_base, my $other_base) = @_;

	if($this_base && ($this_base eq "A" || $this_base eq "C" || $this_base eq "G" || $this_base eq "T"))
	{
		return($this_base);
	}
	elsif($this_base)
	{
		if($other_base eq $ref_base)
		{
			return($var_base);
		}
		else
		{
			return($ref_base);
		}
	}
	else
	{
		return($var_base);
	}
		
}





#############################################################
# parse_file - takes input file and parses it
#
#############################################################

sub check_platform
{
	(my $platform) = @_;


	if($platform =~ "454" && $platform =~ "3730")
	{
		$platform = "Roche 454 / ABI 3730xl";
	}
	elsif($platform =~ "454")
	{
		$platform = "Roche 454";
	}
	elsif($platform =~ "Illumina")
	{
		$platform = "Illumina GAIIx";
	}
	elsif($platform =~ "3730")
	{
		$platform = "ABI 3730xl";
	}
	else
	{
		$platform = "Roche 454";
	}
	
	return($platform);
}







1;

