
package Genome::Model::Tools::Analysis::454::LoadReads;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# LoadReads - Load 454 reads from a sample-SFF tab-delimited file
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

class Genome::Model::Tools::Analysis::454::LoadReads {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		samples_file	=> { is => 'Text', doc => "Tab-delimited file of sample and SFF file(s)" },
		output_dir	=> { is => 'Text', doc => "Output directory" },
		skip_if_present => { is => 'Text', doc => "Skip if SFF/Fasta/Qual files are present", is_optional => 1},
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Load 454 reads from samples.tsv and run BLAT alignment"                 
}

sub help_synopsis {
    return <<EOS
This command loads 454 data from a samples.tsv file and launches BLAT alignments
EXAMPLE:	gt analysis 454 load-reads --samples-file data/samples.tsv --output-dir data
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

	my $script_dir = $output_dir . "/scripts";

	if(!(-e $samples_file))
	{
		die "Error: Samples file not found!\n";
	}

	## Create output directories ##
	mkdir($output_dir) if(!(-d $output_dir));
	mkdir($script_dir) if(!(-d $script_dir));			



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
		my $sff_dir = $sample_output_dir . "/sff";
		my $fasta_dir = $sample_output_dir . "/fasta_dir";
		my $blat_dir = $sample_output_dir . "/blat_out";

		mkdir($sample_output_dir) if(!(-d $sample_output_dir));
		mkdir($sff_dir) if(!(-d $sff_dir));
		mkdir($fasta_dir) if(!(-d $fasta_dir));
		mkdir($blat_dir) if(!(-d $blat_dir));
		mkdir("$blat_dir/pslx") if(!(-d "$blat_dir/pslx"));

		## Open the sample setup script file ##
		
		my $ScriptFileName = "$script_dir/setup_" . $sample_name . ".sh";
		my $ScriptFileOut = "$script_dir/setup_" . $sample_name . ".sh.out";		
		open(SAMPLESCRIPT, ">$ScriptFileName") or die "Can't open outfile: $!\n";
		print SAMPLESCRIPT "#!/gsc/bin/sh\n";
		print SAMPLESCRIPT qq{date\n};	
		my $cmd;
		
		## Compile the SFF files into a single one ##
		
		print SAMPLESCRIPT qq{echo "Compiling SFF files..."\n};
		$cmd = "sfffile -o $sff_dir/$sample_name.sff $sff_files";
		print SAMPLESCRIPT "$cmd\n";
		
		
		## Extract sequence and quality ##
		print SAMPLESCRIPT qq{echo "Extracting FASTA sequence..."\n};
		$cmd = "sffinfo -s $sff_dir/$sample_name.sff >$fasta_dir/$sample_name.fasta";
		print SAMPLESCRIPT "$cmd\n";		
		
		print SAMPLESCRIPT qq{echo "Extracting QUALITY scores..."\n};
		$cmd = "sffinfo -q $sff_dir/$sample_name.sff >$fasta_dir/$sample_name.fasta.qual";
		print SAMPLESCRIPT "$cmd\n";
		
		## Run BLAT alignments ##	
		
		print SAMPLESCRIPT qq{echo "Running initial BLAT alignments..."\n};
		$cmd = "gmt blat align-to-genome --query-file $fasta_dir/$sample_name.fasta --output-dir $blat_dir/pslx";
#		print SAMPLESCRIPT "$cmd\n";			
		
		## Finish up and close the file ##
		
		print SAMPLESCRIPT qq{echo "Completed successfully!"\n};
		print SAMPLESCRIPT qq{date\n};	
		close(SAMPLESCRIPT);
		
		## Make the script executable ##
		system("chmod a+x $ScriptFileName");
		
		## Submit to the blades ##

		print "$sample_name\t$sample_output_dir\n";
		if($self->skip_if_present && -e "$sff_dir/$sample_name.sff" && -e "$fasta_dir/$sample_name.fasta" && -e "$fasta_dir/$sample_name.fasta.qual")
		{
			print "Skipping because output present...\n";
		}
		else
		{
			system(qq{bsub -q $ENV{GENOME_LSF_QUEUE_SHORT} -oo $ScriptFileOut -R "select[mem>2000] rusage[mem=2000]" $ScriptFileName});					
		}

	}
	
	close($input);

	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}



1;

