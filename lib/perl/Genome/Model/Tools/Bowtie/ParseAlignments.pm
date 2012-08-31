
package Genome::Model::Tools::Bowtie::ParseAlignments;     # rename this when you give the module file a different name <--

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

class Genome::Model::Tools::Bowtie::ParseAlignments {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		alignments_file	=> { is => 'Text', doc => "File containing Bowtie output" },
		output_file	=> { is => 'Text', doc => "Output substitution events to file" },
                output_blocks	=> { is => 'Text', doc => "Output alignment blocks to file", is_optional => 1 },
                min_qual_score	=> { is => 'Text', doc => "Minimum quality score to include SNP", is_optional => 1 },
                verbose	=> { is => 'Text', doc => "Print interim progress updates", is_optional => 1 },
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Parse alignment output from Bowtie"                 
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
        my $output_blocks = $self->output_blocks if(defined($self->output_blocks));
        my $min_qual_score = 0;
        $min_qual_score = $self->min_qual_score if(defined($self->min_qual_score));
        my $verbose = 1 if(defined($self->verbose));


	my %stats = ();

	if($output_file)
	{
		open(OUTFILE, ">$output_file") or die "Can't open outfile: $!\n";
		print OUTFILE "chrom\tposition\tref_allele\tvar_allele\tread_name\tread_pos\talign_strand\tqual_score\n";
	}
        
        if($output_blocks)
        {
                open(BLOCKS, ">$output_blocks") or die "Can't open outfile: $!\n";
                print BLOCKS "ref_name\tref_start\tref_stop\talign_strand\tread_name\tread_start\tread_stop\tblock_number\n";
        }

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
                        
                        if($verbose)
                        {
                                print $stats{'num_alignments'} . " alignments parsed...\n" if(!($stats{'num_alignments'} % 10000))
                        }
			
			my @lineContents = split(/\t/, $line);
			my $read_name = $lineContents[0];
			my $align_strand = $lineContents[1];
			my $chromosome = $lineContents[2];
			my $chr_start = $lineContents[3];
			my $read_seq = $lineContents[4];
			my $read_qual = $lineContents[5];
			my $mismatch_list = $lineContents[7];
			
                        ## Adjust chromosome position to 1-based coordinates ##
                        
                        $chr_start++;
                        
			my $read_length = length($read_seq);
                        my $chr_stop = $chr_start + $read_length - 1;
			
                        if($output_blocks)
                        {
                                print BLOCKS "$chromosome\t$chr_start\t$chr_stop\t$align_strand\t$read_name\t1\t$read_length\t1\n";        
                        }
                        
                        
			if($mismatch_list)
			{
#				print "$line\n";
				my @mismatches = split(/\,/, $mismatch_list);
				foreach my $mismatch (@mismatches)
				{
					my $snp_chr_pos = my $qual_code = my $qual_score = 0;
					(my $read_pos, my $allele1, my $allele2) = split(/[\>\:]/, $mismatch);
				
                                	if($allele1 ne "N" && $allele2 ne "N")
					{
                                                if($align_strand eq "-")
                                                {
                                                        $snp_chr_pos = $chr_start + ($read_length - $read_pos - 1);
                                                        $qual_code = substr($read_qual, ($read_length - $read_pos - 1), 1);
                                                        ## Do NOT Correct alleles - bowtie gives read sequence in + orientation ##
#                                                        $allele1 = flip_allele($allele1);
 #                                                       $allele2 = flip_allele($allele2);
                                                }
                                                else
                                                {
                                                        $snp_chr_pos = $chr_start + $read_pos;# - 1;
                                                        $qual_code = substr($read_qual, $read_pos, 1); #substr($read_qual, ($read_pos - 1), 1);
                                                }
                                                
                                                $qual_score = ord($qual_code) - 33;

                                                if($qual_score >= $min_qual_score)
                                                {
                                                        $stats{'num_snps'}++;
                                                        if($output_file)
                                                        {
                                                                print OUTFILE "$chromosome\t$snp_chr_pos\t$allele1\t$allele2\t$read_name\t$read_pos\t$align_strand\t$qual_score\n";
                                                        }
                                                }
					}					
				}
			}
			
#			return(0) if($lineCounter > 10);
		}
	}
	
	close($input);

	print "$stats{'num_alignments'} reads aligned\n";
	print "$stats{'num_snps'} substitution events\n";
	
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

