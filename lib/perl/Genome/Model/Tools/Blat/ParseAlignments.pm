
package Genome::Model::Tools::Blat::ParseAlignments;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# ParseAlignments.pm - 	Parse the BLAT PSL output file; score and sort alignments
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@watson.wustl.edu)
#
#	CREATED:	10/20/2008 by D.K.
#	MODIFIED:	11/20/2008 by D.K.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;

use FileHandle;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

class Genome::Model::Tools::Blat::ParseAlignments {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		alignments_file	=> { is => 'Text', doc => "File containing all BLAT alignments in PSL/PSLX format" },
		min_identity	=> { is => 'Text', doc => "Minimum % identity to keep an alignment [0]", is_optional => 1},
		output_psl	=> { is => 'Text', doc => "Output scored PSL for best alignments [0]", is_optional => 1},	
		output_blocks	=> { is => 'Text', doc => "Output alignment blocks for best alignments [0]", is_optional => 1},		
		find_secondary	=> { is => 'Text', doc => "Find secondary alignments for multi-hit reads [0]", is_optional => 1},	
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Parse, score, and sort BLAT alignments"                 
}

sub help_synopsis {
    return <<EOS
This command parses a PSL-format BLAT output file, scores each alignment, and reports the best alignment for uniquely placed reads.
EXAMPLE:	gmt blat parse-alignments myBlatOutput.psl
Scored alignments for uniquely placed reads would be output to myBlatOutput.psl.best-alignments.txt
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
	my $alignments_file = $self->alignments_file;

	my $output_basename = $alignments_file;

	## Set defaults for optional parameters ##
	
	my $min_pct_identity = 0;
	my $output_psl = my $output_blocks = my $find_secondary = 0;
	my $scoreM = 1;	# Match points for scoring alignments
	my $scoreN = 2; # Mismatch penalty for scoring alignments
	my $scoreQ = 3; # Gap penalty for scoring alignments

	## Get optional parameters ##
	
	$min_pct_identity = $self->min_identity if($self->min_identity);
	$output_psl = $self->output_psl if($self->output_psl);
	$output_blocks = $self->output_blocks if($self->output_blocks);
	$find_secondary = $self->find_secondary if($self->find_secondary);

	## Verify that alignments file exists ##
	
	if(!(-e $alignments_file))
	{
		print "Input file does not exist. Exiting...\n";
		return(0);
	}	

	my %ReadAlignments = ();
	
	
	## Reset statistics ##
	
	my %AlignmentStats = ();
	
	## Open the infile ##

	print "Parsing PSL file...\n";
	my $input = new FileHandle ($alignments_file);
	my $lineCounter = my $pslFormatCounter = 0;

	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
	
		## Parse out only lines matching a BLAT-like result pattern ##
		
		if($line && $line =~ '\,') #$line=~/\d+\t\d+\d+\t\d+\t\d+\d+\t/)	
		{
			$pslFormatCounter++;
			
			my @lineContents = split(/\t/, $line);		
			my $numContents = @lineContents;
		
			$AlignmentStats{'num_hsps'}++;
			
			my $match_bases = $lineContents[0];
			my $mismatch_bases = $lineContents[1];
			my $rep_match_bases = $lineContents[2];
			my $query_gaps = $lineContents[4];
			my $subject_gaps = $lineContents[6];
			
			my $read_name = $lineContents[9];
			my $read_length = $lineContents[10];
			my $read_start = $lineContents[11];
			my $read_stop = $lineContents[12];				
			
			my $ref_name = $lineContents[13];
			my $num_blocks = $lineContents[17];
			my $alignment_gaps = $num_blocks - 1;
		
			## Calculate read-proportional length of the alignment ##
			
			my $pct_read_aligned = ($read_stop - $read_start) / $read_length * 100;
			
			## Calculate the sequence identity of the alignment ##
			my $alignment_identity = 0;
			$alignment_identity = ($match_bases + $rep_match_bases) / ($match_bases + $rep_match_bases + $mismatch_bases) * 100 if($match_bases || $rep_match_bases);			
		
			## Calculate a BLAST-like alignment score ##
			
			my $blast_score = (($match_bases + $rep_match_bases) * $scoreM) - ($mismatch_bases * $scoreN) - ($query_gaps * $scoreQ) - ($subject_gaps * $scoreQ);
			
			if($blast_score <= 0)
			{
				$AlignmentStats{'zero_score_hsp'}++;
			}
			elsif($alignment_identity < $min_pct_identity)
			{
				$AlignmentStats{'low_id_hsp'}++;
			}
			else
			{
				## Save the score along with the read alignment ##
				
				$ReadAlignments{$read_name} .= "\n" if($ReadAlignments{$read_name});
				$ReadAlignments{$read_name} .= $blast_score . "\t" . $line;
				
				$AlignmentStats{'scored_hsp'}++;
			}
		
		} # else line is not PSL format
	}
	
	close($input);


	## Warn if no lines parsed ##

	if($lineCounter > 0 && $pslFormatCounter == 0)
	{
		print "Input file contained 0 lines in PSL format...\n";
		return(0);	
	}

	print "Sorting and printing alignments...\n";

	## Define the alignment-sorting subroutine ##
	
	sub byAlignScore
	{
		my @temp = split(/\t/, $a);
		my $score_a = $temp[0];
		@temp = split(/\t/, $b);
		my $score_b = $temp[0];		
		$score_b <=> $score_a;
	}
	
	
	## Open output files ##

	open(BESTALIGNS, ">$output_basename.best-alignments.txt") or die "Can't open outfile $output_basename.best-alignments.txt: $!\n";
	print BESTALIGNS "score\tread_name\tread_start\tread_stop\tref_name\tref_start\tref_stop\talign_strand\tidentity\tpct_read_aligned\n";

	if($output_psl)
	{
		open(BESTPSL, ">$output_basename.best-alignments.psl") or die "Can't open outfile $output_basename.best-alignments.psl: $!\n";
		print BESTPSL "score\tmatch\tmismatch\trep\tN\tQgap_count\tQgap_bases\tTgap_count\tTgap_bases\tstrand\tname\tqSize\tqstart\tqend\ttname\ttsize\ttstart\ttend\tblock\tblockSizes\tblockStartQ\tblockStartT\tseq1\tseq2\n";
	}


	if($output_blocks)
	{
		open(BESTBLOCKS, ">$output_basename.best-blocks.txt") or die "Can't open outfile $output_basename.best-blocks.txt: $!\n";
		print BESTBLOCKS "ref_name\tref_start\tref_stop\talign_strand\tread_name\tread_start\tread_stop\tblock_no\n";
	}

	if($find_secondary)
	{
		open(SECONDARYALIGNS, ">$output_basename.secondary-alignments.txt") or die "Can't open outfile $output_basename.secondary-alignments.txt: $!\n";
		print SECONDARYALIGNS "score\tread_name\tread_start\tread_stop\tref_name\tref_start\tref_stop\talign_strand\tidentity\tpct_read_aligned\n";		
	}

	## Process alignments on read basis ##
	
	foreach my $ReadName (keys %ReadAlignments)
	{
		$AlignmentStats{'num_reads'}++;
	
		## Reset the read category ##
		
		my $read_category = "unknown";
	
		## Sort the HSPs for this read by their BLAST-like score ##
		
		my @ReadHSPs = split(/\n/, $ReadAlignments{$ReadName});
		@ReadHSPs = sort byAlignScore @ReadHSPs;
		my $numReadHSPs = @ReadHSPs;
	
		if($numReadHSPs == 1)
		{
			$read_category = "single_alignment";
		}
		elsif($numReadHSPs > 1)
		{
			(my $best_score) = split(/\t/, $ReadHSPs[0]);
			(my $second_best_score) = split(/\t/, $ReadHSPs[1]);
		
			## Determine how close the primary and secondary hit are ##
		
			my $score_diff = $best_score - $second_best_score;
			my $score_diff_pct = $score_diff / $best_score * 100;
		
		
			## Categorize based on score closeness ##
			
			if($score_diff_pct <= 1.00)
			{
				## Discard reads with competing best alignments ##
				$read_category = "competing_alignments";
#				print COMPETING "$ReadHSPs[0]\t$ReadHSPs[1]\n";	
			}
			else
			{
				$read_category = "multiple_alignments";
			}
		}	
	
		$AlignmentStats{$read_category}++;
	
		if($read_category eq "single_alignment" || $read_category eq "multiple_alignments")
		{
			$AlignmentStats{'best_alignments'}++;
			## Parse out the best alignment ##
			
			my @bestContents = split(/\t/, $ReadHSPs[0]);
			my $best_score = $bestContents[0];
			
			my $match_bases = $bestContents[1];
			my $mismatch_bases = $bestContents[2];
			my $rep_match_bases = $bestContents[3];
			my $query_gaps = $bestContents[5];
			my $subject_gaps = $bestContents[7];
						
			my $best_strand = $bestContents[9];
			my $best_read_name = $bestContents[10];
			my $best_read_length = $bestContents[11];			
			my $best_read_start = $bestContents[12];
			my $best_read_stop = $bestContents[13];
			my $best_ref_name = $bestContents[14];
			my $best_ref_start = $bestContents[16];
			my $best_ref_stop = $bestContents[17];
			
			## Adjust to 1-based positions ##
			$best_read_start++;
			$best_read_stop++;
			$best_ref_start++;
			$best_ref_stop++;
			
			## Calculate read-proportional length of the alignment ##
			
			my $pct_read_aligned = ($best_read_stop - $best_read_start) / $best_read_length * 100;
			$pct_read_aligned = sprintf("%.2f", $pct_read_aligned);
			
			## Calculate the sequence identity of the alignment ##
			my $alignment_identity = 0;
			$alignment_identity = ($match_bases + $rep_match_bases) / ($match_bases + $rep_match_bases + $mismatch_bases) * 100 if($match_bases || $rep_match_bases);
			$alignment_identity = sprintf("%.2f", $alignment_identity);
			
			## Print the best alignment ##
			print BESTALIGNS "$best_score\t$best_read_name\t$best_read_start\t$best_read_stop\t$best_ref_name\t$best_ref_start\t$best_ref_stop\t$best_strand\t$alignment_identity\t$pct_read_aligned\n";
			print BESTPSL "$ReadHSPs[0]\n" if($output_psl);
			
			## Get the best blocks ##
			if($output_blocks)
			{
				my $best_align_blocks = parseAlignmentBlocks($ReadHSPs[0]);
				print BESTBLOCKS $best_align_blocks if($best_align_blocks);			
			}
			
			
			
			for(my $secondCounter = 1; $secondCounter < $numReadHSPs; $secondCounter++)
			{
				my @secondContents = split(/\t/, $ReadHSPs[$secondCounter]);
				my $second_score = $secondContents[0];
				
				my $match_bases = $secondContents[1];
				my $mismatch_bases = $secondContents[2];
				my $rep_match_bases = $secondContents[3];
				my $query_gaps = $secondContents[5];
				my $subject_gaps = $secondContents[7];
							
				my $second_strand = $secondContents[9];
				my $second_read_name = $secondContents[10];
				my $second_read_length = $secondContents[11];			
				my $second_read_start = $secondContents[12];
				my $second_read_stop = $secondContents[13];
				my $second_ref_name = $secondContents[14];
				my $second_ref_start = $secondContents[16];
				my $second_ref_stop = $secondContents[17];
				
				## Adjust to 1-based positions ##
				$second_read_start++;
				$second_read_stop++;
				$second_ref_start++;
				$second_ref_stop++;
	
				## Calculate read-proportional length of the alignment ##
				
				my $pct_read_aligned = ($second_read_stop - $second_read_start) / $second_read_length * 100;
				$pct_read_aligned = sprintf("%.2f", $pct_read_aligned);
				
				## Calculate the sequence identity of the alignment ##
				my $alignment_identity = 0;
				$alignment_identity = ($match_bases + $rep_match_bases) / ($match_bases + $rep_match_bases + $mismatch_bases) * 100 if($match_bases || $rep_match_bases);
				$alignment_identity = sprintf("%.2f", $alignment_identity);

				print SECONDARYALIGNS "$second_score\t$second_read_name\t$second_read_start\t$second_read_stop\t$second_ref_name\t$second_ref_start\t$second_ref_stop\t$second_strand\t$alignment_identity\t$pct_read_aligned\n" if($find_secondary);

#				print SECONDARY "$ReadHSPs[$secondCounter]\n";
			}
		}
	}

	close(BESTALIGNS);
	close(BESTPSL) if($output_psl);
	close(BESTBLOCKS) if($output_blocks);

	## Print some statistics ##
		
	print "\t$AlignmentStats{'num_hsps'} alignments in the PSL file\n" if($AlignmentStats{'num_hsps'});	
	print "\t$AlignmentStats{'zero_score_hsp'} alignments had scores of zero\n" if($AlignmentStats{'zero_score_hsp'});		
	print "\t$AlignmentStats{'low_id_hsp'} low-identity alignments discarded\n" if($AlignmentStats{'low_id_hsp'});	
	print "\t$AlignmentStats{'scored_hsp'} scored alignments remained\n" if($AlignmentStats{'scored_hsp'});	
		
	print "\t$AlignmentStats{'num_reads'} reads had scored alignments\n" if($AlignmentStats{'num_reads'});
	print "\t$AlignmentStats{'competing_alignments'} ambigouously-mapped reads discarded\n" if($AlignmentStats{'competing_alignments'});		
	print "\t$AlignmentStats{'best_alignments'} uniquely mapped reads isolated\n" if($AlignmentStats{'best_alignments'});
	
	
	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}


################################################################################################
# parseAlignmentBlocks - get the blocks of gapped PSL alignments
#
################################################################################################

sub parseAlignmentBlocks
{
	my $alignmentLine = shift(@_);
	my @lineContents = split(/\t/, $alignmentLine);
	
	my $blocks = "";
	
	my $align_strand = $lineContents[9];
	my $read_name = $lineContents[10];
	my $ref_name = $lineContents[14];	
	
	my $num_blocks = $lineContents[18];
	my $block_sizes = $lineContents[19];
	my $read_block_starts = $lineContents[20];
	my $ref_block_starts = $lineContents[21];

	## Build arrays of the block info ##
	
	my @blockSizes = split(/,/, $block_sizes);
	my @readBlockStarts = split(/,/, $read_block_starts);
	my @refBlockStarts = split(/,/, $ref_block_starts);					

	## Loop through each block ##

	for(my $bCounter = 0; $bCounter < $num_blocks; $bCounter++)
	{
		my $block_number = $bCounter + 1;
		my $block_size = $blockSizes[$bCounter];
		
		my $read_block_start = $readBlockStarts[$bCounter];
		my $ref_block_start = $refBlockStarts[$bCounter];
	
		my $read_block_stop = $read_block_start + $block_size - 1;
		my $ref_block_stop = $ref_block_start + $block_size - 1;	
	
		## Fix zero-based positions ##
		$ref_block_start++;
		$ref_block_stop++;
		$read_block_start++;
		$read_block_stop++;

		## add to the list ##
		$blocks .= "$ref_name\t$ref_block_start\t$ref_block_stop\t$align_strand\t$read_name\t$read_block_start\t$read_block_stop\t$block_number\n";
	}

	return($blocks);
}



1;

