
package Genome::Model::Tools::Analysis::Hapmap::ExtractGenotypes;     # rename this when you give the module file a different name <--

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

my %stats = ();


class Genome::Model::Tools::Analysis::Hapmap::ExtractGenotypes {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		genotypes	=> { is => 'Text', doc => "Location of downloaded HapMap genotype file, comma-separated if multiple", is_optional => 0, is_input => 1},
		sample_name	=> { is => 'Text', doc => "Name of sample for which to extract genotypes", is_optional => 0, is_input => 1},
		output_file	=> { is => 'Text', doc => "Output file for extracted genotypes", is_optional => 1, is_input => 1},
		full_output	=> { is => 'Text', doc => "If set to 1, print extended information (RS#, alleles, strand)", is_optional => 1, is_input => 1},
		include_blanks	=> { is => 'Text', doc => "If set to 1, export blank (NN) genotypes", is_optional => 1, is_input => 1},
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Extracts genotypes for a single HapMap sample"                 
}

sub help_synopsis {
    return <<EOS
This command extracts genotypes for a single HapMap sample
EXAMPLE:	gmt analysis hapmap extract-genotypes --genotype-file genotypes_chr1_ASW_r27_nr.b36_fwd.txt --sample-name NA19654 --output-file NA19654.genotypes
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

	my $genotypes = $self->genotypes;
	my $sample_name = $self->sample_name;
	my $output_file = $self->output_file;

	if($self->output_file)
	{
		open(OUTFILE, ">" . $self->output_file) or die "Can't open outfile: $!\n";		
	}

	
	## Get all genotype files ##
	
	my @genotype_files = split(/\,/, $genotypes);
	
	foreach my $genotype_file (@genotype_files)
	{
		parse_genotypes($genotype_file, $sample_name, $self);
		$stats{'files_parsed'}++;
	}

	if($self->output_file)
	{
		close(OUTFILE);
	}
	
#	print $stats{'not_in_unaffected'} . " were NOT present in unaffected individuals\n";

	print $stats{'files_parsed'} . " files parsed\n";
	print commify($stats{'num_genotypes'}) . " genotypes parsed for $sample_name\n";
	print commify($stats{'genotypes_were_blank'}) . " were blank (NN)\n";
	print commify($stats{'genotypes_not_blank'}) . " were non-blank\n";
	print commify($stats{'genotypes_extracted'}) . " genotypes extracted\n";
	
	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}




################################################################################################
# Load Genotypes
#
################################################################################################

sub parse_genotypes
{
	my ($genotype_file, $sample_name, $self) = @_;

	warn "Parsing $genotype_file...\n";

	my $input = new FileHandle ($genotype_file);
	my $lineCounter = 0;
	
	my $sampleColumn = 0;

	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		my @lineContents = split(/\s+/, $line);
		my $numContents = @lineContents;
		
		if($lineCounter == 1)
		{
			for(my $colCounter = 1; $colCounter < $numContents; $colCounter++)
			{
				if($lineContents[$colCounter] eq $sample_name)
				{
					$sampleColumn = $colCounter;
					warn "Parsing genotypes for $sample_name from column $colCounter\n";
				}
			}
		}
		elsif($sampleColumn)
		{
			my $rs_number = $lineContents[0];
			my $alleles = $lineContents[1];
			my $chromosome = $lineContents[2];
			my $position = $lineContents[3];
			my $strand = $lineContents[4];
			
			my $genotype = $lineContents[$sampleColumn];
			$stats{'num_genotypes'}++;
			
			$chromosome =~ s/chr//;
		
			my $print_flag = 0;
			
			if($genotype ne "NN")
			{
				$print_flag = 1;
				$stats{'genotypes_not_blank'}++;				
			}
			else
			{
				$stats{'genotypes_were_blank'}++;
				
				if($self->include_blanks)
				{
					$print_flag = 1;
				}
			}
			
			if($print_flag)
			{
				if($self->full_output)
				{
					print OUTFILE join("\t", $chromosome, $position, $rs_number, $alleles, $strand, $genotype) . "\n";	
				}
				else
				{
					print OUTFILE join("\t", $chromosome, $position, $genotype) . "\n";					
				}

				$stats{'genotypes_extracted'}++;				
			}
		}
		else
		{
			close($input);
			return;
		}
	}
	
	close($input);
	
}



sub commify
{
	local($_) = shift;
	1 while s/^(-?\d+)(\d{3})/$1,$2/;
	return $_;
}

1;

