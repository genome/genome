package Genome::Model::Tools::Plink::BuildFiles;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# MutationRate - Calculate the mutation rate (per megabase) given a list of mutations (e.g. tier1 SNVs) and a set of regions (e.g. coding space)
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@genome.wustl.edu)
#
#	CREATED:	04/22/2011 by D.K.
#	MODIFIED:	04/22/2011 by D.K.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;

use FileHandle;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

## Declare global statistics hash ##

my %stats = ();

class Genome::Model::Tools::Plink::BuildFiles {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		input	=> { is => 'Text', doc => "Input basename for PLINK. There should be a basename.ped and a basename.map" , is_optional => 0, is_input => '1'},
		output	=> { is => 'Text', doc => "Output basename for PLINK files. By default, matches input basename" , is_optional => 1, is_input => '1'},
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Builds binary and sample-distance (genome) files for Plink analysis"                 
}

sub help_synopsis {
    return <<EOS
This command builds binary and sample-distance (genome) files for Plink analysis
EXAMPLE:	gmt plink build-files --input all-iscan-snps --output all-iscan-snps ...
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

	my $input = $self->input;
	my $output = $input;
	$output = $self->output if($self->output);


	if(!(-e "$input.ped"))
	{
		die "Error: Required file $input.ped not found\n";
	}
	if(!(-e "$input.map"))
	{
		die "Error: Required file $input.map not found\n";
	}

	## Build binary plink file ##
	
	if(!(-e "$input.bed"))
	{
		print "Building binary PED (BED) file...\n";
		my $cmd = Genome::Model::Tools::Plink->path_to_binary() . " --file $input --make-bed --out $output";
		print "$cmd\n";
		system($cmd);
		
	}

	## Build genome file if possible ##

	if(-e "$output.bed")
	{
		if(!(-e "$output.genome"))
		{
			print "Building cluster distances file to speed clustering...\n";
			my $cmd = Genome::Model::Tools::Plink->path_to_binary() . " --bfile $output --genome --out $output";
			print "$cmd\n";
			system($cmd);			
		}
	}

	## IF genome file exists, run MDS ##

	if((-e "$output.genome"))
	{
		if(!(-e "$output.mds"))
		{
			print "Running MDS clustering...\n";
			my $cmd = Genome::Model::Tools::Plink->path_to_binary() . "--bfile $output --read-genome $output.genome --cluster --out $output --mds-plot 5";
			print "$cmd\n";
			system($cmd);						
		}

		if(!(-e "$output.frq"))
		{
			print "Running MAF Analysis...\n";
			my $cmd = Genome::Model::Tools::Plink->path_to_binary() . "--bfile $output --freq --out $output ";
			print "$cmd\n";
			system($cmd);						
		}

		if(!(-e "$output.lmiss"))
		{
			print "Running Missingness Analysis...\n";
			my $cmd = Genome::Model::Tools::Plink->path_to_binary() . "--bfile $output --missing --out $output ";
			print "$cmd\n";
			system($cmd);						
		}


	}

	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}



1;
