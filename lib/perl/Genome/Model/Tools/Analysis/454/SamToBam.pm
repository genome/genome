
package Genome::Model::Tools::Analysis::454::SamToBam;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# SamToBam - Align reads with SSAHA2 or other aligner
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@watson.wustl.edu)
#
#	CREATED:	02/25/2009 by D.K.
#	MODIFIED:	02/25/2009 by D.K.
#
#	NOTES:	
#			
#####################################################################################################################################

my $ref_seq = Genome::Config::reference_sequence_directory() . '/NCBI-human-build36/all_sequences.fa';
my $ref_index = Genome::Config::reference_sequence_directory() . '/NCBI-human-build36/all_sequences.fa.fai';

use strict;
use warnings;

use FileHandle;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

class Genome::Model::Tools::Analysis::454::SamToBam {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		samples_file	=> { is => 'Text', doc => "Tab-delimited file of sample and SFF file(s)" },
		output_dir	=> { is => 'Text', doc => "Output directory" },
		aligner		=> { is => 'Text', doc => "Aligner that was used." },
		reference		=> { is => 'Text', doc => "Reference sequence [default=Hs36 ssaha2]", is_optional => 1 },
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Align 454 reads using FASTA locations from samples.tsv"                 
}

sub help_synopsis {
    return <<EOS
This command converts SAM files to BAM files
EXAMPLE:	gmt analysis 454 align-reads --samples-file data/samples.tsv --output-dir data --aligner ssaha2
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
	my $samples_file = $self->samples_file;
	my $output_dir = $self->output_dir;
	my $aligner = $self->aligner;

	if($self->reference)
	{
		$ref_seq = $self->reference;
		$ref_index = $ref_seq . ".fai";

		die "Reference files not found!\n$ref_seq\n$ref_index\n" if(!(-e $ref_seq && -e $ref_index));
	}

	if(!(-e $samples_file))
	{
		die "Error: Samples file not found!\n";
	}

	my $input = new FileHandle ($samples_file);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		
	
		my @lineContents = split(/\t/, $line);			
		my $sample_name = $lineContents[0];
		my $sff_files = $lineContents[1];				
		
		## Create sample output directories ##
		my $sample_output_dir = $output_dir . "/" . $sample_name;
		my $aligner_output_dir = $sample_output_dir . "/" . $aligner . "_out";
		my $aligner_scripts_dir = $sample_output_dir . "/" . $aligner . "_out/scripts";		
		my $aligner_output_file = "$aligner_output_dir/$sample_name.$aligner.sam";
		my $bam_file = "$aligner_output_dir/$sample_name.$aligner.bam";
		
		if(-e $aligner_output_file)
		{
			mkdir($aligner_scripts_dir) if(!(-d $aligner_scripts_dir));
			
			## open outfile ##
			my $script_filename = "$aligner_scripts_dir/sam_to_bam.sh";
			open(SCRIPT, ">$script_filename") or die "Can't open script file: $!\n";
			
			## Step 1: Convert SAM to BAM ##
			print SCRIPT "echo Converting SAM to BAM...\n";
			my $cmd = "samtools view -b -t $ref_index -o $bam_file $aligner_output_file";
			print SCRIPT "$cmd\n";

			## Step 2: Sort the BAM ##
			print SCRIPT "echo Sorting the BAM file...\n";			
			$cmd = "samtools sort $bam_file $bam_file.sorted";
			print SCRIPT "$cmd\n";
			
			## Step 3: Rename BAM file ##
			$cmd = "mv -f $bam_file.sorted.bam $bam_file";		
			print SCRIPT "$cmd\n";

			## Step 4: Index BAM file ##
			$cmd = "samtools index $bam_file";		
			print SCRIPT "$cmd\n";

			## Step 5: Generate Pileup ##
			print SCRIPT "echo Building the pileup file...\n";	
			$cmd = "samtools mpileup -f $ref_seq $bam_file >$bam_file.pileup";
			print SCRIPT "$cmd\n";

			## Update July 2011: Skip variant calling because pileup no longer supported ##
			## Step 6: Call SAMtools variants ##
#			print SCRIPT "echo Building the variants file...\n";	
#			$cmd = "samtools mpileup -c -v -f $ref_seq $bam_file >$bam_file.variants";
#			print SCRIPT "$cmd\n";

			## Step 7: Filter SAMtools variants ##
#			print SCRIPT "echo Building the pileup file...\n";	
#			$cmd = "samtools.pl varFilter -D 10000 $bam_file.variants >$bam_file.variants.filter";
#			print SCRIPT "$cmd\n";
			
			close(SCRIPT);
			
			## Change permissions ##
			system("chmod 755 $script_filename");
			
			## Run bsub ##
#			system("$script_filename");			
			system("bsub -q $ENV{GENOME_LSF_QUEUE_BUILD_WORKER_ALT} -R\"select[mem>2000] rusage[mem=2000]\" -oo $script_filename.out $script_filename");
		}
		else
		{
			die "Alignment output file $aligner_output_file not found!\n";
		}
		
	}
	
	close($input);

	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}



1;

