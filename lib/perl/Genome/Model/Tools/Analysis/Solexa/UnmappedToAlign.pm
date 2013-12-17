
package Genome::Model::Tools::Analysis::Solexa::UnmappedToAlign;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# UnmappedToAlign - 	Extract unmapped reads from a BAM file using SAMtools view with filters
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@watson.wustl.edu)
#
#	CREATED:	11/23/2009 by D.K.
#	MODIFIED:	11/23/2009 by D.K.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;
use FileHandle;
use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

my $aligner = "novoalign";
my $num_cores = 1;
my $lsf_queue = $ENV{GENOME_LSF_QUEUE_BUILD_WORKER};
my $bsub_cmd = "bsub -q $lsf_queue -R\"select[type==LINUX64 && model != Opteron250 && mem>8000] rusage[mem=8000] span[hosts=1]\" -n $num_cores -M 8000000";
my $novoalign_reference = Genome::Config::reference_sequence_directory() . '/NCBI-human-build36/all_sequences.novoindex-k14-s3-v2.05.13';
my $path_to_novoalign = "/gscuser/dkoboldt/Software/NovoCraft/novocraftV2.05.13/novocraft/novoalign";
#my $novoalign_params = "-c $num_cores -a -l 36 -t 240 -k -s 5 -d $novoalign_reference"; #-o SAM 
my $novoalign_params = "-c $num_cores -a -l 50 -t 240 -k -d $novoalign_reference"; #-o SAM 



class Genome::Model::Tools::Analysis::Solexa::UnmappedToAlign {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		output_dir	=> { is => 'Text', doc => "Directory to contain output files", is_optional => 0 },
		output_name	=> { is => 'Text', doc => "Name to use for output files, i.e. lane", is_optional => 0 },
		aligner	=> { is => 'Text', doc => "Aligner to use [novoalign]", is_optional => 1 },
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Maps unmapped reads from BAM files"                 
}

sub help_synopsis {
    return <<EOS
This command maps unmapped Illumina/Solexa reads
EXAMPLE:	gmt analysis solexa unmapped-to-align --output-dir unmapped --output-name myBam --aligner novoalign
			==> unmapped/myBam_1_sequence.pe.unmapped
			==> unmapped/myBam_2_sequence.pe.unmapped
			==> unmapped/myBam_1_sequence.se.unmapped
			==> unmapped/myBam_2_sequence.se.unmapped
			==> unmapped/myBam_sequence.se.unmapped			
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
	my $output_dir = $self->output_dir;
	my $name = $self->output_name;

	my $unmapped_pair_read1 = $output_dir . "/" . $name . "_1_sequence.pe.unmapped.fastq";
	my $unmapped_pair_read2 = $output_dir . "/" . $name . "_2_sequence.pe.unmapped.fastq";
	my $unmapped_frag_read1 = $output_dir . "/" . $name . "_1_sequence.unmapped.fastq";
	my $unmapped_frag_read2 = $output_dir . "/" . $name . "_2_sequence.unmapped.fastq";
	my $unmapped_single_read = $output_dir . "/" . $name . "_sequence.unmapped.fastq";
	
	## Specify alignment output dir ##
	
	my $alignment_output_dir = $output_dir . "/" . $aligner . "_out";
	mkdir($alignment_output_dir) if(!(-d $alignment_output_dir));
	
	my $unmapped_pair_out = $alignment_output_dir . "/" . $name . "noTrim_sequence.pe.unmapped.$aligner";
	my $unmapped_pair1_out = $alignment_output_dir . "/" . $name . "noTrim_sequence.read1.unmapped.$aligner";
	my $unmapped_pair2_out = $alignment_output_dir . "/" . $name . "noTrim_sequence.read2.unmapped.$aligner";
	my $unmapped_frag1_out = $alignment_output_dir . "/" . $name . "noTrim_sequence.frag1.unmapped.$aligner";
	my $unmapped_frag2_out = $alignment_output_dir . "/" . $name . "noTrim_sequence.frag2.unmapped.$aligner";
	my $unmapped_single_out = $alignment_output_dir . "/" . $name . "noTrim_sequence.single.unmapped.$aligner";

	## check for existence and size of files ##
	
	if(check_fastq($unmapped_pair_read1) && check_fastq($unmapped_pair_read2))
	{
		## Run Paired-end Alignment ##
		my $cmd = "$path_to_novoalign $novoalign_params -i 250 100 -f $unmapped_pair_read1 $unmapped_pair_read2 >$unmapped_pair_out";
		print "$unmapped_pair_out\n";
		system("$bsub_cmd -oo $unmapped_pair_out.log \"$cmd\"");		

		## Run Single-end Alignment for Read 1 ##
		$cmd = "$path_to_novoalign $novoalign_params -f $unmapped_pair_read1 >$unmapped_pair1_out";
		print "$unmapped_pair1_out\n";
		system("$bsub_cmd -oo $unmapped_pair1_out.log \"$cmd\"");		

		## Run Single-end Alignment for Read 2 ##
		$cmd = "$path_to_novoalign $novoalign_params -f $unmapped_pair_read2 >$unmapped_pair2_out";
		print "$unmapped_pair2_out\n";
		system("$bsub_cmd -oo $unmapped_pair2_out.log \"$cmd\"");
	}
	
	if(check_fastq($unmapped_frag_read1))
	{
		## Run Single-end Alignment ##
		my $cmd = "$path_to_novoalign $novoalign_params -f $unmapped_frag_read1 >$unmapped_frag1_out";
		print "$unmapped_frag1_out\n";
		system("$bsub_cmd -oo $unmapped_frag1_out.log \"$cmd\"");			
	}

	if(check_fastq($unmapped_frag_read2))
	{
		## Run Single-end Alignment ##
		my $cmd = "$path_to_novoalign $novoalign_params -f $unmapped_frag_read2 >$unmapped_frag2_out";
		print "$unmapped_frag2_out\n";
		system("$bsub_cmd -oo $unmapped_frag2_out.log \"$cmd\"");		
	}

	if(check_fastq($unmapped_single_read))
	{
		## Run Single-end Alignment ##
		my $cmd = "$path_to_novoalign $novoalign_params -f $unmapped_single_read >$unmapped_single_out";
		print "$unmapped_single_out\n";
		system("$bsub_cmd -oo $unmapped_single_out.log \"$cmd\"");
	}

	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}


################################################################################################
# Execute - the main program logic
#
################################################################################################

sub check_fastq
{                               
	(my $fastq_file) = @_;

	if(-e $fastq_file)
	{
		## Get file stats ##
		
		my @filestat = stat $fastq_file;
		
		my $file_size = $filestat[7];
		
		print "$fastq_file\t$file_size\n" if($file_size);
		return($file_size);
	}
	
	return(0);
}


1;

