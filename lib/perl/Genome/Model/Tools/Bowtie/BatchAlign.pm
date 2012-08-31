
package Genome::Model::Tools::Bowtie::BatchAlign;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# BatchAlign.pm - 	Align reads to a reference genome using Bowtie
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@watson.wustl.edu)
#
#	CREATED:	04/22/2009 by D.K.
#	MODIFIED:	04/22/2009 by D.K.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;

use FileHandle;

## Bowtie Parameters ##
my $batch_size = 1000000;
my $num_cores = 1;
my $lsf_queue = "long";

my $novoalign_params = "-c $num_cores -a -l 36 -t 240 -k";	# -o SAM

my $path_to_novoalign = "/gscuser/dkoboldt/Software/NovoCraft/novocraftV2.05.13/novocraft/novoalign";
my $novoalign_reference = Genome::Config::reference_sequence_directory() . '/NCBI-human-build36/all_sequences.novoindex-k14-s3-v2.05.13';
#my $novoalign_reference = Genome::Config::reference_sequence_directory() . '/NCBI-human-build36/all_sequences.novoindex-k15-s2-v2.05.13';


#my $batch_dir = "/tmp/novoalign";

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

class Genome::Model::Tools::Bowtie::BatchAlign {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		query_file	=> { is => 'Text', doc => "Illumina/Solexa reads in FASTQ format" },
		batch_dir	=> { is => 'Text', doc => "Batch directory for intermediate files"},
		batch_name	=> { is => 'Text', doc => "Batch name for intermediate files", is_optional => 1 },
                reference	=> { is => 'Text', doc => "Path to bowtie-indexed reference [Defaults to Hs36]", is_optional => 1 },
		batch_size	=> { is => 'Text', doc => "Number of reads per batch [$batch_size]", is_optional => 1 },
		novo_params	=> { is => 'Text', doc => "Parameters for novoalign [$novoalign_params]", is_optional => 1 },
		num_cores	=> { is => 'Text', doc => "Number of cores for novoalign [$num_cores]", is_optional => 1 },
		lsf_queue	=> { is => 'Text', doc => "LSF queue to use for bsub [$lsf_queue]", is_optional => 1 },
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Batch-align reads to a reference genome using Bowtie"                 
}

sub help_synopsis {
    return <<EOS
This command breaks Illumina SE FASTQ files into batches and launches novoalign jobs
EXAMPLE:	gmt novoalign batch-align --query-file s_1_sequence.fastq --batch-dir novoalign_out
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
	my $query_file = $self->query_file;
	my $batch_dir = $self->batch_dir;
	$batch_size = $self->batch_size if($self->batch_size);
	$novoalign_params = $self->novo_params if($self->novo_params);
	$num_cores = $self->num_cores if($self->num_cores);
	$lsf_queue = $self->lsf_queue if($self->lsf_queue);

	if(!(-e $query_file))
	{
		die "Error: Query file not found!\n";
	}

	## Define Bowtie Reference (default to Hs36)

	my $reference = $novoalign_reference;

        if(defined($self->reference))
	{
		if(-e $self->reference)
		{
			$reference = $self->reference;
		}
		else
		{
			die "Error: Reference file not found!\n";
		}
	}


	## Get query file ##
	
	if(-e $query_file)
	{
		mkdir($batch_dir) if(!(-d $batch_dir));
		my $temp_fileroot;
		$temp_fileroot = $self->batch_name if($self->batch_name);
		$temp_fileroot = time() if(!$temp_fileroot);
		
		my $num_lines = $batch_size * 4;
		
		print "Split file basename: $temp_fileroot\n";
		print "Splitting the FASTQ file into batches of $batch_size reads...\n";
		system("split -l $num_lines $query_file $batch_dir/$temp_fileroot.");
#		system("split -l 4000000 $query_file $query_file.batch.");
		
		print "Obtaining split files...\n";
		
		my $file_list = `ls $batch_dir/$temp_fileroot.*`;
#		my $file_list = `ls $query_file.batch.*`;
		chomp($file_list);
		
		my @fastq_files = split(/\n/, $file_list);
		
		foreach my $fastq_file (@fastq_files)
		{
			## Get batch name @@
			
			my @temp = split(/\./, $fastq_file);
			my $batch_name = $temp[1];
			
			## Correct for files ending in .fa ##
			
			if($batch_name eq "fa")
			{
				system("mv $fastq_file " . $fastq_file . "q");
				$fastq_file .= "q";
				$batch_name .= "q";
			}
			
			my $batch_output_file = $fastq_file . ".novoalign";
			print "$fastq_file\t$batch_output_file\n";
			
			## Run the alignment ##
			system("bsub -q $lsf_queue -R\"select[type==LINUX64 && model != Opteron250 && mem>12000] rusage[mem=12000] span[hosts=1]\" -n $num_cores -M 12000000 -oo $batch_output_file.log \"$path_to_novoalign $novoalign_params -d $reference -f $fastq_file >$batch_output_file 2>$batch_output_file.err\"");
		}

		## Soft link to original fastq ##
		
		system("ln -s $query_file $batch_dir/$temp_fileroot.fastq");


	}

	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}



1;

