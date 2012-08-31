
package Genome::Model::Tools::Analysis::Mendelian::GenotypeTable;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# SearchRuns - Search the database for runs
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@watson.wustl.edu)
#
#	CREATED:	04/01/2009 by D.K.
#	MODIFIED:	04/01/2009 by D.K.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;

use FileHandle;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

my $num_affected = my $affecteds_missing = my $unaffecteds_variant = my $affecteds_variant = my $affecteds_ambiguous = 0;
my %genotypes = ();

class Genome::Model::Tools::Analysis::Mendelian::GenotypeTable {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		variant_file	=> { is => 'Text', doc => "List of variants to consider (annotation format)", is_optional => 0, is_input => 1},
		consensus_files	=> { is => 'Text', doc => "One or more consensus files for samples", is_optional => 0, is_input => 1},
		output_file	=> { is => 'Text', doc => "Output file for QC result", is_optional => 1, is_input => 1},
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Compiles a genotype table for all variants in a file"                 
}

sub help_synopsis {
    return <<EOS
This command compiles a genotype table for all variants in a file
EXAMPLE:	gmt analysis mendelian genotype-table --variant-file all.snvs --consensus-files myCons1,myCons2
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

	my $variant_file = $self->variant_file;

	my $consensus_files = $self->consensus_files;
	
	if($self->output_file)
	{
		open(OUTFILE, ">" . $self->output_file) or die "Can't open outfile: $!\n";		
	}

	my %stats = ();
	
	## Build an array of affected individuals' genotypes ##
		
	my @consensus_files = split(/\,/, $consensus_files);
		
	## Count the files of each type and print the header ##
	my $header = "variant";
	foreach my $consensus_file (@consensus_files)
	{
		$num_affected++;
		$header .= "\t" if($header);
		$header .= $consensus_file . "\treads1\treads2\tfreq";
		load_consensus($consensus_file);
	}

	if($self->output_file)
	{
		print OUTFILE "$header\n";
	}

	## Print the variants ##

	my $input = new FileHandle ($variant_file);
	my $lineCounter = 0;

	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		(my $chromosome, my $chr_start, my $chr_stop, my $ref, my $var) = split(/\t/, $line);

		if(length($ref) > 1 || $var eq "-")
		{
			## Undo the adjustment made when formatting deletions for annotation.
			$chr_start--;
		}

		my $variant_type = "SNP";
		
		if(length($ref) > 1 || $var eq "-")
		{
			$variant_type = "DEL";
		}
		elsif(length($var) > 1 || $ref eq "-")
		{
			$variant_type = "INS";
		}

		if($lineCounter >= 0)
		{
			$stats{'num_variants'}++;
						
			my $sample_genotype_string = "";
			
			## See how many affecteds carry it ##
			
			foreach my $consensus_file (@consensus_files)
			{
#				my %genotypes = load_consensus($consensus_file);
				my $sample_genotype = "-\t-\t-\t-";

				my $key = "$consensus_file\t$chromosome\t$chr_start";
				
				if($genotypes{$key})
				{
					(my $sample_call, my $sample_reads1, my $sample_reads2, my $sample_freq) = split(/\t/, $genotypes{$key});
					$sample_call = code_to_genotype($sample_call);					
					$sample_genotype = "$sample_call\t$sample_reads1\t$sample_reads2\t$sample_freq";
				}
				else
				{
					$affecteds_missing++;
				}

				$sample_genotype_string .= $sample_genotype . "\t";
			}
		
			if($self->output_file)
			{
				print OUTFILE "$line\t";
#						print OUTFILE "$affecteds_variant\t$unaffecteds_variant\t";
				print OUTFILE "$sample_genotype_string";
				
				
				print OUTFILE "\n";
			}					
		
			

		}		
		
	}
	
	close($input);
	
	if($self->output_file)
	{
		close(OUTFILE);
	}
	
	print $stats{'num_variants'} . " variants\n";

	
	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}




################################################################################################
# Load Genotypes
#
################################################################################################

sub load_consensus
{                               # replace with real execution logic.
	my $genotype_file = shift(@_);
#	my %genotypes = ();
	
	my $input = new FileHandle ($genotype_file);
	my $lineCounter = 0;
	my $gtCounter = 0;

	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;

		(my $chrom, my $position, my $ref, my $cns, my $reads1, my $reads2, my $var_freq, my $strands1, my $strands2, my $qual1, my $qual2, my $p_value) = split(/\t/, $line);

		if($ref =~ /[0-9]/)
		{
			($chrom, $position, my $stop, $ref, $cns, $reads1, $reads2, $var_freq, $strands1, $strands2, $qual1, $qual2, $p_value) = split(/\t/, $line);			
		}

		if(length($ref) > 1 || length($cns) > 1 || $ref eq "-" || $cns eq "-")
		{
			## If CNS is not formatted, do so ##

			if(!($cns =~ '/'))
			{
				$cns = "$ref/$cns";				
			}

		}
#		my $key = "$chrom\t$position";
		my $key = "$genotype_file\t$chrom\t$position";
		if($cns ne "N")
		{
			$genotypes{$key} = "$cns\t$reads1\t$reads2\t$var_freq";								
		}

	}
	close($input);
                            # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)

#	return(%genotypes);
}






sub code_to_genotype
{
	my $code = shift(@_);
	
	return("AA") if($code eq "A");
	return("CC") if($code eq "C");
	return("GG") if($code eq "G");
	return("TT") if($code eq "T");

	return("AC") if($code eq "M");
	return("AG") if($code eq "R");
	return("AT") if($code eq "W");
	return("CG") if($code eq "S");
	return("CT") if($code eq "Y");
	return("GT") if($code eq "K");

#	warn "Unrecognized ambiguity code $code!\n";
#	return("NN");
	return($code);
}


sub commify
{
	local($_) = shift;
	1 while s/^(-?\d+)(\d{3})/$1,$2/;
	return $_;
}


1;


