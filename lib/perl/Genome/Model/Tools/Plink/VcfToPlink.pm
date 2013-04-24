package Genome::Model::Tools::Plink::VcfToPlink;     # rename this when you give the module file a different name <--

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

class Genome::Model::Tools::Plink::VcfToPlink {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		vcf_file	=> { is => 'Text', doc => "Input VCF file to be converted" , is_optional => 0, is_input => '1'},
		pedigree_file	=> { is => 'Text', doc => "A standard pedigree file with phenotype, gender, & relationships" , is_optional => 0, is_input => '1'},
		phenotype_file	=> { is => 'Text', doc => "Name for output phenotype file for PLINK" , is_optional => 0, is_input => '1'},
		plink_file	=> { is => 'Text', doc => "Name for output PLINK file" , is_optional => 0, is_input => '1'},
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Converts a VCF to PLINK format"                 
}

sub help_synopsis {
    return <<EOS
This command converts a VCF to PLINK format
EXAMPLE:	gmt plink vcf-to-plink --vcf-file [file] --pedigree-file [file] --phenotype-file [name] --plink-file [name]...
EOS
}

sub help_detail {                           # this is what the user will see with the longer version of help. <---
    return <<EOS 
The steps to convert a VCF to PLINK are as follows:
1.) Run VCFtools to convert to PLINK BED format
2.) Update individual & family IDs
3.) Update genders
4.) Update parent relationships
5.) Create a correct phenotype file to be used with --pheno command in PLINK
EOS
}


################################################################################################
# Execute - the main program logic
#
################################################################################################

sub execute {                               # replace with real execution logic.
	my $self = shift;

	my $input = $self->input;

#	my $cmd = Genome::Model::Tools::Plink->path_to_binary() . " --file $input --make-bed --out $output";



	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}



1;
