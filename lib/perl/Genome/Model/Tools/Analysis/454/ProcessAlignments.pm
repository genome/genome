
package Genome::Model::Tools::Analysis::454::ProcessAlignments;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# ProcessAlignments - Load 454 reads from a sample-SFF tab-delimited file
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

class Genome::Model::Tools::Analysis::454::ProcessAlignments {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		samples_file	=> { is => 'Text', doc => "Tab-delimited file of sample and SFF file(s)" },
		output_dir	=> { is => 'Text', doc => "Output directory" },
		regions_file	=> { is => 'Text', doc => "If provided, matches read alignments to regions", is_optional => 1 },
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Processes alignments of 454 data"                 
}

sub help_synopsis {
    return <<EOS
This command processes alignments for 454 datasets
	1.) Compiles individual per-chromosome BLAT alignments into single alignment file
	2.) Parses out the best alignments and their alignment blocks
	3.) Runs Varscan to detect SNPs and indels
	
EXAMPLE: gmt analysis 454 process-alignments --samples-file samples.tsv --output-dir data
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
	my $regions_file = $self->regions_file if($self->regions_file);

	my $script_dir = $output_dir . "/scripts";

	if(!(-e $samples_file))
	{
		die "Error: Samples file not found!\n";
	}
	
	## Create output directories ##
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
		my $varscan_dir = $sample_output_dir . "/varscan_out";

		if(-d $sample_output_dir)
		{
			## Open the sample setup script file ##
			
			my $ScriptFileName = "$script_dir/process2_" . $sample_name . ".sh";
			my $ScriptFileOut = "$script_dir/process2_" . $sample_name . ".sh.out";		
			open(SAMPLESCRIPT, ">$ScriptFileName") or die "Can't open outfile: $!\n";
			print SAMPLESCRIPT "#!/gsc/bin/sh\n";
			print SAMPLESCRIPT qq{date\n};	
			my $cmd;

			## Set up initial bsub command ##
			
			my $bsub_cmd = qq{bsub -q long -R \"select[type==LINUX64 && mem>4000] rusage[mem=4000]\" -oo $ScriptFileOut sh $ScriptFileName};


			if($regions_file)
			{
				## If blocks file exists, convert to layers ##
				if(-e "$blat_dir/$sample_name.psl.best-blocks.txt")
				{
					## Match reads to regions ##
					$cmd = "perl ~dkoboldt/Scripts/match_blat_to_regions.pl amplicons/regions.txt $blat_dir/$sample_name.psl.best-blocks.txt $blat_dir/$sample_name.psl.best.regions.layers";
					print SAMPLESCRIPT "$cmd\n";
				}
			}
			else
			{
				## STANDARD PROCESSING - COMPILE ALIGNMENTS AND VARSCAN ##
				
				## Create necessary output dirs ##
				
				mkdir($varscan_dir) if(!(-d $varscan_dir));
				
	
	
				## Count the number of existing BLAT files ##
	
				my $blatFiles = `find $blat_dir/pslx -name \"$sample_name*.psl\" -print | wc -l`;
				chomp($blatFiles);
				$blatFiles =~ s/[^0-9]//g;
				
				print "\t$blatFiles BLAT files\n";
				
				if($blatFiles >= 24 && $blatFiles <= 30)
				{
					## Compile the blat output for this sample ##
					print SAMPLESCRIPT qq{echo "Compiling BLAT PSLX files..."\n};
					$cmd = "cat $blat_dir/pslx/$sample_name*.psl >$blat_dir/$sample_name.psl";	
					print SAMPLESCRIPT "$cmd\n";
					## Add unmapped files ##
					#$cmd = "cat data/$patient_id/blat_out/pslx/unmapped*.psl >>data/$patient_id/blat_out/$sample_name.psl";	
					#print SAMPLESCRIPT "$cmd\n";				
				}
				else
				{
					warn "Incorrect # of BLAT files ($blatFiles) for this sample; skipping this sample!\n";
					$bsub_cmd = "";
				}
				
				## Parse the alignments ##
	
				print SAMPLESCRIPT qq{echo "Parsing BLAT alignments..."\n};
				$cmd = "gmt blat parse-alignments --output-blocks 1 --alignments-file $blat_dir/$sample_name.psl";
				print SAMPLESCRIPT "$cmd\n";
	
				print SAMPLESCRIPT qq{echo "Running Varscan..."\n};
				$cmd = "varscan easyrun $blat_dir/$sample_name.psl --fasta-file $fasta_dir/$sample_name.fasta --quality-file $fasta_dir/$sample_name.fasta.qual --output-dir $varscan_dir --sample $sample_name --ref-dir /gscmnt/sata180/info/medseq/biodb/shared/Hs_build36_mask1c/ --min-align-score 50 --min-identity 90 --primer-trim 10 --min-coverage 10 --min-reads2 2 --min-var-freq 0.25";
				print SAMPLESCRIPT "$cmd\n";					
			}

			print SAMPLESCRIPT qq{echo "Completed successfully!"\n};
			print SAMPLESCRIPT qq{date\n};	
			close(SAMPLESCRIPT);

			## Make the script executable ##
			system("chmod a+x $ScriptFileName");

			## Submit to the blades ##
	
			print "$sample_name\t$sample_output_dir\n";
			system($bsub_cmd);					

		}


	}
	
	close($input);

	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}



1;

