
package Genome::Model::Tools::Bowtie::DeDuplicate;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# GetUnmappedReads.pm - 	Get unmapped/poorly-mapped reads by model id
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

class Genome::Model::Tools::Bowtie::DeDuplicate {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		alignments_file	=> { is => 'Text', doc => "File containing Bowtie output" },
		output_file	=> { is => 'Text', doc => "Output substitution events to file" },
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "De-duplicate Bowtie alignments"                 
}

sub help_synopsis {
    return <<EOS
This command retrieves the locations of unplaced reads for a given genome model
EXAMPLE:	gmt bowtie parse-alignments --alignments-file bowtie.txt
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
	my $output_file = $self->output_file;

	my %stats = my %included = ();

	if($output_file)
	{
		open(OUTFILE, ">$output_file") or die "Can't open outfile: $!\n";
#		print OUTFILE "chrom\tposition\tref_allele\tvar_allele\tread_name\tread_pos\talign_strand\tqual_score\n";
	}
        
	my %best_read = my %best_score = ();
	
	my $input = new FileHandle ($alignments_file);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		
	
		if($lineCounter >= 1)
		{
			$stats{'num_alignments'}++;
                        
			print $stats{'num_alignments'} . " alignments parsed...\n" if(!($stats{'num_alignments'} % 10000));
 			
			my @lineContents = split(/\t/, $line);
			my $read_name = $lineContents[0];
			my $align_strand = $lineContents[1];
			my $chromosome = $lineContents[2];
			my $chr_start = $lineContents[3];
			my $read_seq = $lineContents[4];
			my $read_qual = $lineContents[5];
			
			## Calculate average base quality ##
			my $qual_sum = my $qual_num = 0;
			
			my @quals = split(//, $read_qual);
			foreach my $qual_code (@quals)
			{
				my $qual_score = ord($qual_code) - 33;
				$qual_sum += $qual_score;
				$qual_num++;
			}
			
			my $avg_qual = $qual_sum / $qual_num if($qual_num);

			
#			my $alignment_key = "$chromosome\t$chr_start\t$align_strand\t$read_seq";
			my $alignment_key = "$chromosome\t$chr_start\t$align_strand";

			if(!$best_read{$alignment_key} || $avg_qual > $best_score{$alignment_key})
			{
				$best_read{$alignment_key} = $read_name;
				$best_score{$alignment_key} = $avg_qual;
			}



			## Old method - doesn't try to find best match ##

#			if(!$included{$alignment_key})
#			{
#				$stats{'num_unique_alignments'}++;
#				print OUTFILE "$line\n";
#			}
#			else
#			{
#				$stats{'num_duplicate_alignments'}++;
#			}
#
#			$included{$alignment_key}++;
			

		}
	}
	
	close($input);

	print "Re-parsing alignments...\n";
	$stats{'num_alignments'} = 0;
	
	## Reparse the file, printing only the best reads ##

	$input = new FileHandle ($alignments_file);
	$lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		
	
		if($lineCounter >= 1)
		{
			$stats{'num_alignments'}++;
                        
			print $stats{'num_alignments'} . " alignments parsed...\n" if(!($stats{'num_alignments'} % 10000));
 			
			my @lineContents = split(/\t/, $line);
			my $read_name = $lineContents[0];
			my $align_strand = $lineContents[1];
			my $chromosome = $lineContents[2];
			my $chr_start = $lineContents[3];
			my $read_seq = $lineContents[4];
						
#			my $alignment_key = "$chromosome\t$chr_start\t$align_strand\t$read_seq";
			my $alignment_key = "$chromosome\t$chr_start\t$align_strand";

			if($best_read{$alignment_key} && $best_read{$alignment_key} eq $read_name)
			{
				$stats{'num_unique_alignments'}++;
				print OUTFILE "$line\n";
			}
			else
			{
				$stats{'num_duplicate_alignments'}++;
			}
			

		}
	}
	
	close($input);


	## Calculate duplication rate ##
	
	$stats{'duplication_rate'} = $stats{'num_duplicate_alignments'} / $stats{'num_alignments'} * 100;
	my $duplication_rate = sprintf("%.2f", $stats{'duplication_rate'}) . '%';

	print "$stats{'num_alignments'} alignments\n";
	print "$stats{'num_unique_alignments'} unique alignments\n";
	print "$stats{'num_duplicate_alignments'} duplicate alignments\n";	
	print "$duplication_rate duplication rate\n";
	
	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}



################################################################################################
# Execute - the main program logic
#
################################################################################################

sub flip_allele {
        my $allele = shift(@_);
        my $flipped = "";
        
        if(uc($allele) eq "A")
        {
                 $flipped = "T";    
        }
        elsif(uc($allele) eq "C")
        {
                 $flipped = "G";    
        }
        elsif(uc($allele) eq "G")
        {
                 $flipped = "C";    
        }
        elsif(uc($allele) eq "T")
        {
                 $flipped = "A";    
        }
        else
        {
            $flipped = "N";    
        }
        
        return($flipped);
}


1;

