
package Genome::Model::Tools::Analysis::454::ShowStats;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# ShowStats - Load 454 reads from a sample-SFF tab-delimited file
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

class Genome::Model::Tools::Analysis::454::ShowStats {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		samples_file	=> { is => 'Text', doc => "Tab-delimited file of sample and SFF file(s)" },
		output_dir	=> { is => 'Text', doc => "Output directory" },
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Align reads to a reference genome using Bowtie"                 
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

	my $script_dir = $output_dir . "/scripts";

	if(!(-e $samples_file))
	{
		die "Error: Samples file not found!\n";
	}
	
	print "SAMPLE  \tNUM_READS  \tTOTAL_BP  \tAVG_READLEN  \tREADS_ALIGNED \tALIGN_PCT \n";
	
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
			print $sample_name . "  \t";
			
			my $num_reads = my $average_rl = my $total_read_bp = my $total_read_kbp = my $num_aligned = my $align_pct = 0;
			
			$num_reads = `grep -c \">\" $fasta_dir/$sample_name.fasta`;
			chomp($num_reads);


			## Get the average read length ##
			
			my $read_lengths = `fastalength $fasta_dir/$sample_name.fasta`;
			chomp($read_lengths);
			my @readLens = split(/\n/, $read_lengths);
			$total_read_bp = 0;
			
			foreach my $read_length_line (@readLens)
			{
				(my $read_length) = split(/\s+/, $read_length_line);
				$total_read_bp += $read_length;
			}
			
			$average_rl = $total_read_bp / $num_reads;
			$average_rl = sprintf("%.1f", $average_rl);

			$total_read_kbp = sprintf("%d", ($total_read_bp / 1000));
			
			## Count the alignments ##
			
			my $blat_file = "$blat_dir/$sample_name.psl.best-alignments.txt";		
			
			if(-e $blat_file)
			{
				$num_aligned = `wc -l $blat_file`;
				chomp($num_aligned);
				$num_aligned--;
			}
			
			## Calculate align percent ##
						
			$align_pct = $num_aligned / $num_reads * 100 if($num_reads);
			$align_pct = sprintf("%.2f", $align_pct) . '%';



		
			print commify($num_reads) . "  \t";
			print commify($total_read_bp) . "  \t";
			print $average_rl . "    \t";
			print commify($num_aligned) . "  \t";
			print $align_pct . "\n";
		}


	}
	
	close($input);

	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}



sub commify
{
	local($_) = shift;
	1 while s/^(-?\d+)(\d{3})/$1,$2/;
	return $_;
}


1;

