
package Genome::Model::Tools::Analysis::Solexa::AlignReads;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# LoadReads - Run maq sol2sanger on Illumina/Solexa files in a gerald_directory
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@watson.wustl.edu)
#
#	CREATED:	04/01/2009 by D.K.
#	MODIFIED:	04/01/2009 by D.K.
#
#	NOTES:	
#			
#####################################################################################################################################

my $bowtie_reference = Genome::Config::reference_sequence_directory() . '/NCBI-human-build36/all_sequences.bowtie';
my $bowtie_params = "-m 1 --best --strata -p 4";	# --trim3 25";

my $path_to_novoalign = "/gscuser/dkoboldt/Software/NovoCraft/novocraftV2.05.07/novocraft/novoalign";
my $novoalign_reference = Genome::Config::reference_sequence_directory() . '/NCBI-human-build36/all_sequences.novoindex-k14-s3-v2.05.13';
my $novoalign_params = "-a -l 36 -t 240";	# -o SAM

use strict;
use warnings;
use Cwd;

use FileHandle;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

class Genome::Model::Tools::Analysis::Solexa::AlignReads {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		flowcell_dir	=> { is => 'Text', doc => "Dir containing a fastq_dir", is_optional => 0 },
		include_lanes	=> { is => 'Text', doc => "Specify which lanes of a flowcell to include [e.g. 1,2,3]" , is_optional => 1},
		output_dir	=> { is => 'Text', doc => "Output dir for Bowtie; defaults to [dir]/bowtie_out" , is_optional => 1},
		aligner	=> { is => 'Text', doc => "Alignment tool to use (bowtie|maq|novoalign) [bowtie]" , is_optional => 1},
		trim3	=> { is => 'Text', doc => "Bases to trim from 3' end for Bowtie alignment [0]" , is_optional => 1},
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Obtains reads from a flowcell_id in FastQ format"                 
}

sub help_synopsis {
    return <<EOS
This command aligns reads to Hs36 (by default) after you've run load-reads
EXAMPLE 1:	gmt analysis solexa align-reads --flowcell_id 302RT --include-lanes 1,2,3,4 --output-dir output_dir --aligner bowtie
EXAMPLE 2:	gmt analysis solexa align-reads --sample-name H_GP-0365n --output-dir H_GP-0365n
EXAMPLE 3:	gmt analysis solexa align-reads --library-name H_GP-0365n-lib2 --output-dir H_GP-0365n
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
	my $flowcell_dir = $self->flowcell_dir;
	my $aligner = "bowtie";
	$aligner = $self->aligner if($self->aligner);
	my $output_dir = "./";
	$output_dir = $self->output_dir if($self->output_dir);

	if($self->trim3)
	{
		$bowtie_params .= " --trim3 " . $self->trim3;
	}

	## Handle include-lanes when specified ##

	my $include_lanes;
	$include_lanes = $self->include_lanes if($self->include_lanes);
	my %lanes_to_include = ();
	
	if($include_lanes)
	{
		my @lanes = split(/\,/, $include_lanes);
		foreach my $desired_lane (@lanes)
		{
			$lanes_to_include{$desired_lane} = 1;
		}
	}
	
	## Get all fastq files in flowcell dir ##
	
	if(-d $flowcell_dir)
	{
		my $fastq_dir = $flowcell_dir . "/fastq_dir";

		## Create alignment output directory ##
		my $alignment_dir = $flowcell_dir . "/" . $aligner . "_out";	
		mkdir($alignment_dir) if(!(-d $alignment_dir));

		## Approach 1: look for all possible lanes ##
		
		for(my $lane = 1; $lane <= 8; $lane++)
		{
			if(!$include_lanes || $lanes_to_include{$lane})
			{
				## Fragment-end read ##
				
				if(-e "$fastq_dir/s_" . $lane . "_sequence.fastq")
				{
					my $fastq_file1 = "$fastq_dir/s_" . $lane . "_sequence.fastq";

					## Get the read length ##
					
					my $seq = `head -2 $fastq_file1 | tail -1`;
					chomp($seq);
					my $read_len = length($seq);

					## Run the alignment ##
					
					if($aligner eq "bowtie")
					{
						my $reference = $bowtie_reference;

						## Launch SE ##
						print "$fastq_file1\tbowtie SE\n";
						my $alignment_outfile1 = $alignment_dir . "/s_" . $lane . "_sequence.$aligner";			# span[hosts=1]  -n 4
						system("bsub -q alignment -R\"select[type==LINUX64 && model != Opteron250 && mem>4000] rusage[mem=4000]\" -M 6000000 -oo $alignment_outfile1.log bowtie $bowtie_params --unfq $alignment_outfile1.unmapped.fastq --maxfq $alignment_outfile1.multiple.fastq $reference $fastq_file1 $alignment_outfile1");
					}
					elsif($aligner eq "novoalign")
					{
						## Adjust for shorter reads ##
						$novoalign_params = "-a" if($read_len <= 36);
						## Think about setting -l 50 for 75 bp reads
						## Launch SE ##
						my $reference = $novoalign_reference;
						#$novoalign_params = "-a -l 50 -t 240" if($read_length >= 70);	

						## Novoalign SE ##

						print "$fastq_file1 \t$read_len bp\tnovoalign SE\n";
						my $alignment_outfile = $alignment_dir . "/s_" . $lane . "_sequence.$aligner";

						system("bsub -q long -R\"select[type==LINUX64 && model != Opteron250 && mem>12000] rusage[mem=12000]\" -M 20000000 -oo $alignment_outfile.log \"$path_to_novoalign $novoalign_params -d $reference -f $fastq_file1 >$alignment_outfile\"");
					}					

				}
				
				## Paired-end reads ##
				
				elsif(-e "$fastq_dir/s_" . $lane . "_1_sequence.fastq")
				{
					my $fastq_file1 = "$fastq_dir/s_" . $lane . "_1_sequence.fastq";
					my $fastq_file2 = "$fastq_dir/s_" . $lane . "_2_sequence.fastq";

					## Get the read length ##
					
					my $seq1 = `head -2 $fastq_file1 | tail -1`;
					chomp($seq1);
					my $read_len1 = length($seq1);
                         
                         my $seq2 = `head -2 $fastq_file2 | tail -1`;
                         chomp($seq2);
                         my $read_len2 = length($seq2);
					
					## Run the alignment ##
					
					if($aligner eq "bowtie")
					{
						my $reference = $bowtie_reference;

						## Launch SE ##
						print "$fastq_file1\tbowtie SE\n";
						my $alignment_outfile1 = $alignment_dir . "/s_" . $lane . "_1_sequence.$aligner";
						system("bsub -q long -R\"select[type==LINUX64 && model != Opteron250 && mem>4000] rusage[mem=4000] span[hosts=1]\" -n 4 -M 6000000 -oo $alignment_outfile1.log bowtie $bowtie_params --unfq $alignment_outfile1.unmapped.fastq --maxfq $alignment_outfile1.multiple.fastq $reference $fastq_file1 $alignment_outfile1");

						print "$fastq_file2\tbowtie SE\n";
						my $alignment_outfile2 = $alignment_dir . "/s_" . $lane . "_2_sequence.$aligner";
						system("bsub -q long -R\"select[type==LINUX64 && model != Opteron250 && mem>4000] rusage[mem=4000] span[hosts=1]\" -n 4 -M 6000000 -oo $alignment_outfile2.log bowtie $bowtie_params --unfq $alignment_outfile2.unmapped.fastq --maxfq $alignment_outfile2.multiple.fastq $reference $fastq_file2 $alignment_outfile2");
					}
					elsif($aligner eq "novoalign")
					{
						## Adjust for shorter reads ##
						$novoalign_params = "-a" if($read_len1 <= 36 && $read_len2 <= 36);
						## Think about setting -l 50 for 75 bp reads
						## Launch SE ##
						my $reference = $novoalign_reference;
						#$novoalign_params = "-a -l 50 -t 240" if($read_length >= 70);	

#						print "$fastq_file1\tnovoalign SE\n";
#						my $alignment_outfile1 = $alignment_dir . "/s_" . $lane . "_1_sequence.$aligner";
#						system("bsub -q long -R\"select[type==LINUX64 && model != Opteron250 && mem>12000] rusage[mem=12000]\" -M 20000000 -oo $alignment_outfile1.log \"$path_to_novoalign $novoalign_params -d $reference -f $fastq_file1 >$alignment_outfile1\"");

#						print "$fastq_file2\tnovoalign SE\n";
#						my $alignment_outfile2 = $alignment_dir . "/s_" . $lane . "_2_sequence.$aligner";
#						system("bsub -q long -R\"select[type==LINUX64 && model != Opteron250 && mem>12000] rusage[mem=12000]\" -M 20000000 -oo $alignment_outfile2.log \"$path_to_novoalign $novoalign_params -d $reference -f $fastq_file2 >$alignment_outfile2\"");

						## Novoalign PE ##

						print "$fastq_file1 $fastq_file2\t$read_len1 bp\tnovoalign PE\n";
						my $alignment_outfile = $alignment_dir . "/s_" . $lane . "_sequence.$aligner";

						system("bsub -q long -R\"select[type==LINUX64 && model != Opteron250 && mem>12000] rusage[mem=12000]\" -M 20000000 -oo $alignment_outfile.log \"$path_to_novoalign $novoalign_params -d $reference -f $fastq_file1 $fastq_file2 >$alignment_outfile\"");
					}
				}
			}
		}

		## Approach 2: List all files in fastq_dir ##

		my $fastq_list = `find $fastq_dir -name "*.fastq" -print`;
		chomp($fastq_list);

		my @fastq_files = split(/\n/, $fastq_list);
		
		@fastq_files = sort @fastq_files;
		
		foreach my $fastq_file (@fastq_files)
		{
#			print "$fastq_file\n";
		}
	}

	my $sqlrun = my $rows_returned = "";

	if($sqlrun)
	{
#		print "$sqlrun\n"; exit(0);
		
		print "fcell\tlane\tlibrary_type\tfilt_reads\taln%\tsample_name\tlibrary_name\tstatus\n";
		
		my @lines = split(/\n/, $sqlrun);
		my %lane_pairs = ();
		
		foreach my $line (@lines)
		{
			if($line && (substr($line, 0, 4) eq "FLOW" || substr($line, 0, 1) eq "-"))
			{
				
			}
			elsif($line && $line =~ "Execution")
			{
				($rows_returned) = split(/\s+/, $line);
				print "$rows_returned rows returned\n";
			}
			elsif($line)
			{
				(my $flowcell, my $lane, my $sample, my $library, my $read_length, my $filt_clusters, my $seq_id, my $gerald_dir, my $insert_size, my $align_pct) = split(/\t/, $line);
				
				## Proceed if lane to be included ##
				if(!$include_lanes || $lanes_to_include{$lane})
				{
					## Get num reads ##
					
					my $num_reads = commify($filt_clusters);
					$align_pct = 0 if(!$align_pct);
					$align_pct = sprintf("%.2f", $align_pct) . '%';
					
					## Get SE or PE ##
					
					my $end_type = "SE";
					my $lane_name = $lane;
	
					if($insert_size)
					{
						$end_type = "PE";
						$lane_pairs{"$flowcell.$lane"} = 1 if(!$lane_pairs{"$flowcell.$lane"});
						$lane_name .= "_" . $lane_pairs{"$flowcell.$lane"};
						$lane_pairs{"$flowcell.$lane"}++;
					}
					
					## Create flowcell output dir and fastq output dir if necessary ##
					
					my $flowcell_dir = $output_dir . "/" . $flowcell;
					my $fastq_dir = $output_dir . "/" . $flowcell . "/fastq_dir";
					my $output_fastq = $fastq_dir . "/" . "s_" . $lane_name . "_sequence.fastq";
					
					## Create the output_dir ##
					
					my $alignment_dir = $flowcell_dir . "/" . $aligner . "_out";
					mkdir($alignment_dir) if(!(-d $alignment_dir));

					my $alignment_outfile = $flowcell_dir . "/" . $aligner . "_out/" . "s_" . $lane_name . "_sequence.$aligner";
					
					## Print result ##
					if(-e $output_fastq)
					{
						print "$flowcell \t$lane_name \t$read_length bp $end_type\t$num_reads \t$sample \t$output_fastq\t$alignment_dir\n";

						## Run the alignment ##
						
						if($aligner eq "bowtie")
						{
							my $reference = $bowtie_reference;

							## Launch SE ##
							system("bsub -q long -R\"select[type==LINUX64 && mem>4000] rusage[mem=4000]\" -M 6000000 -oo $alignment_outfile.log bowtie $bowtie_params --unfq $alignment_outfile.unmapped.fastq --maxfq $alignment_outfile.multiple.fastq $reference $output_fastq $alignment_outfile");

							## Launch PE ##
							if($end_type eq "PE" && $lane_pairs{"$flowcell.$lane"} eq "2")
							{
								my $file1 = $output_fastq;
								my $file2 = $output_fastq;
								$file2 =~ s/1\_sequence\.fastq/2\_sequence\.fastq/;
								my $paired_outfile = $flowcell_dir . "/" . $aligner . "_out/" . "s_" . $lane . "_paired.bowtie";
							#	system("bsub -q long -R\"select[type==LINUX64 && mem>4000] rusage[mem=4000]\" -oo $paired_outfile.log bowtie -m 1 -p 4 -I 50 -X 450 --unfq $paired_outfile.unmapped.fastq --maxfq $paired_outfile.multiple.fastq $reference -1 $file1 -2 $file2 $paired_outfile");
							}
						}
						elsif($aligner eq "novoalign")
						{
							## Think about setting -l 50 for 75 bp reads
							## Launch SE ##
							my $reference = $novoalign_reference;
							$novoalign_params = "-a -l 50 -t 240" if($read_length >= 70);	

							system("bsub -q alignment -R\"select[type==LINUX64 && mem>12000] rusage[mem=12000]\" -M 20000000 -oo $alignment_outfile.log \"$path_to_novoalign $novoalign_params -d $reference -f $output_fastq >$alignment_outfile\"");

							## Launch PE ##
							if($end_type eq "PE" && $lane_pairs{"$flowcell.$lane"} eq "2")
							{
								$novoalign_params .= " -i 250 50";
								my $file1 = $output_fastq;
								my $file2 = $output_fastq;
								$file2 =~ s/1\_sequence\.fastq/2\_sequence\.fastq/;
								my $paired_outfile = $flowcell_dir . "/" . $aligner . "_out/" . "s_" . $lane . "_paired.novoalign";
								#system("bsub -q alignment -R\"select[type==LINUX64 && mem>12000] rusage[mem=12000]\" -M 20000000 -oo $paired_outfile.log \"$path_to_novoalign $novoalign_params -d $reference -f $file1 $file2 >$paired_outfile\"");
							}
						}
					}
					else
					{
						print "*** File Missing: $output_fastq\n";
					}

				}
			}
		}
	}

	
	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}

sub commify
{
	local($_) = shift;
	1 while s/^(-?\d+)(\d{3})/$1,$2/;
	return $_;
}


1;

