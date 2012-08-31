
package Genome::Model::Tools::Blat::AlignToGenome;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# AlignToGenome.pm - 	Align a query file (.fasta) to a reference genome on a per-chromosome basis
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@watson.wustl.edu)
#
#	CREATED:	10/20/2008 by D.K.
#	MODIFIED:	03/19/2009 by D.K.
#
#	NOTES:	
#			Default BLAT parameters are:
#				-mask=lower 	(takes advantage of soft-masking of repeats, HIGHLY recommended)
#				-out=pslx	(provides extended-PSL output with read/ref sequences)
#				-noHead		(specifies output without header)
#
#			Default BSUB parameters are: -q long -R"select[mem>3000] rusage[mem=3000]" 
#
#####################################################################################################################################

use strict;
use warnings;

use FileHandle;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

class Genome::Model::Tools::Blat::AlignToGenome {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		reference_dir	=> { is => 'Text', doc => "Directory containing chromosome files in FASTA format /gscmnt/839/info/medseq/reference_sequences/NCBI-human-build36/", is_optional => 1 },
		search_string	=> { is => 'Text', doc => "Seach string to identify chromosome files, [*.fasta]", is_optional => 1 },		
		query_file	=> { is => 'Text', doc => "Query file in FASTA format" },
		output_dir	=> { is => 'Text', doc => "Directory to store output files" },
		params		=> { is => 'Text', doc => "BLAT parameters [-mask=lower -out=pslx -noHead]", is_optional => 1 },
		lsf_queue	=> { is => 'Text', doc => "LSF queue if other than long [long]", is_optional => 1 },	
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Align query to reference sequence by chromosome"                 
}

sub help_synopsis {
    return <<EOS
This command runs BLAT alignments between a query and each chromosome of a reference genome
EOS
}

sub help_detail {                           # this is what the user will see with the longer version of help. <---
    return <<EOS 
This command would spawn BLAT alignments between myDeletions.fasta and each of the 24 human chrom refseqs:
 gmt blat align-to-genome --reference-dir /gscmnt/839/info/medseq/reference_sequences/NCBI-human-build36/ --search-string *.fasta --query-file myDeletions.fasta --output-dir ./ 

By default, each BLAT alignment will be launched in the LONG queue, and outputs will be in the ./ directory.

Output filenames are in the format [output_dir]/[query_filename].[reference_filename].psl.  For example, the above query
would have output files ./myDeletions.fasta.chr1.fa.psl ./myDeletions.fasta.chr2.fa.psl etc.
EOS
}


################################################################################################
# Execute - the main program logic
#
################################################################################################

sub execute {                               # replace with real execution logic.
	my $self = shift;

	## Get required parameters ##
	my $query_file = $self->query_file;
	my $output_dir = $self->output_dir;
	
	## Set defaults for optional parameters ##
	my $dir = "/gscmnt/839/info/medseq/reference_sequences/NCBI-human-build36/";
	my $search_string = "*.fasta";
	my $blat_params = "-mask=lower -out=pslx -noHead";
	my $lsf_queue = "long";

	## Get optional parameters if provided ##
	$dir = $self->reference_dir if($self->reference_dir);
	$search_string = $self->search_string if($self->search_string);
	$blat_params = $self->params if($self->params);
	$lsf_queue = $self->lsf_queue if($self->lsf_queue);

	## Verify that query file exists ##
	
	if(!(-e $query_file))
	{
		print "Input file does not exist. Exiting...\n";
		return(0);
	}	

	## Verify that output dir exists ##
	
	if(!(-d $output_dir))
	{
		print "Output directory does not exist. Exiting...\n";
		return(0);
	}	

	print "Running BLAT with parameters $blat_params\n";

	## Verify that directory exists ##

	if(-d $dir)
	{
		## Retrieve list of ref files ## 
		
		my $ref_files = `ls $dir/$search_string | grep -v all_sequences`;
	
		my $numRefFiles = 0;
	
		foreach my $ref_file (split(/\n/, $ref_files))
		{
			my @refFileContents = split(/\//, $ref_file);
			my $numFileContents = @refFileContents;
			my $chrom_name = $refFileContents[$numFileContents - 1];
			$numRefFiles++;
		
			my @queryFileContents = split(/\//, $query_file);
			$numFileContents = @queryFileContents;
			my $query_name = $queryFileContents[$numFileContents - 1];		
		
			## Build the output file name ##
			
			my $output_file = $output_dir . "/" . $query_name . "." . $chrom_name . ".psl";

			## Launch the BLAT alignment ##
			
			system("bsub -q $lsf_queue -R\"select[mem>3000] rusage[mem=3000]\" -oo $output_file.out blat $ref_file $query_file $blat_params $output_file");
		}
		
		print "BLAT alignments launched against $numRefFiles ref files\n";
	}
	else
	{
		print "Input directory does not exist. Exiting...\n";
		return(0);			
	}


	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}


1;

