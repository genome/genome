package Genome::Model::Tools::Analysis::MutationSpectrum;
package Genome::Model::Tools::Analysis::MutationRate;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# MutationRate - Calculate the mutation rate (per megabase) given a list of mutations (e.g. tier1 SNVs) and a set of regions (e.g. coding space)
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@genome.wustl.edu)
#
#	CREATED:	04/22/2011 by D.K.
#	MODIFIED:	04/22/2011 by D.K.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;

use FileHandle;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

## Declare global statistics hash ##

my %stats = ();

class Genome::Model::Tools::Analysis::MutationRate {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		sample_name	=> { is => 'Text', doc => "Descriptive name for sample" , is_optional => 1, default => 'Sample'},
		tier1_file	=> { is => 'Text', doc => "List of high-confidence tier 1 mutations" , is_optional => 0},
		tier2_file	=> { is => 'Text', doc => "List of high-confidence tier 2 mutations" , is_optional => 1},
		tier3_file	=> { is => 'Text', doc => "List of high-confidence tier 3 mutations" , is_optional => 1},
		tier1_space	=> { is => 'Text', doc => "BED file of tier 1 space" , is_optional => 0, default => '/gscmnt/sata921/info/medseq/make_tier_bed_files/NCBI-human-build36/tier1.bed'},
		tier2_space	=> { is => 'Text', doc => "BED file of tier 2 space" , is_optional => 0, default => '/gscmnt/sata921/info/medseq/make_tier_bed_files/NCBI-human-build36/tier2.bed'},
		tier3_space	=> { is => 'Text', doc => "BED file of tier 3 space" , is_optional => 0, default => '/gscmnt/sata921/info/medseq/make_tier_bed_files/NCBI-human-build36/tier3.bed'},
		coverage_factor	=> { is => 'Text', doc => "Fraction of space covered for mutation detection" , is_optional => 0, default => 1},
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Calculates the average mutation rate per megabase"                 
}

sub help_synopsis {
    return <<EOS
This command calculates the average mutation rate per megabase, given a list of mutations and a mutation-space.
The easiest input would be a tier 1 SNV file, which will calculate the coding mutation rate. By default, the "tier-spaces"
are taken from gmt fast-tier fast-tier BED files. If you want to extrapolate genome-wide mutation rate from exome, provide
a tier 1 SNV file, and use tier1+tier2+tier3 space with coverage_factor of something like 0.02
EXAMPLE:	gmt capture report-flagstat ...
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

        my $coverage_factor = $self->coverage_factor;
	## Get required parameters ##
        
        my $total_bases = my $total_mutations = 0;

        ## Use known values whenever possible ##
        my $tier1_bases = 43884962 if($self->tier1_space eq '/gscmnt/sata921/info/medseq/make_tier_bed_files/NCBI-human-build36/tier1.bed');
        my $tier2_bases = 248205479 if($self->tier2_space eq '/gscmnt/sata921/info/medseq/make_tier_bed_files/NCBI-human-build36/tier2.bed');
        my $tier3_bases = 1198993777 if($self->tier3_space eq '/gscmnt/sata921/info/medseq/make_tier_bed_files/NCBI-human-build36/tier3.bed');

        my $mutation_string = my $rate_string = "";

        ## If other files were provided, calculate bases within them ##        
        
        $tier1_bases = parse_regions_file($self->tier1_space) if(!$tier1_bases);
        $tier2_bases = parse_regions_file($self->tier2_space) if(!$tier2_bases);
        $tier3_bases = parse_regions_file($self->tier3_space) if(!$tier3_bases);
        my $non_tier1_bases = $tier2_bases + $tier3_bases;

        ## Load the tier 1 mutations and calculate its rate ##
            
        my $tier1_mutations = parse_mutations_file($self->tier1_file);      
        my $tier1_rate = ($tier1_mutations / $coverage_factor) / ($tier1_bases / 1000000);
        
#        print "TIER\tMUTS\tTOTAL_BASES\tMUTS_PER_MB\n";
#       print join("\t", "1", $tier1_mutations, commify($tier1_bases), $tier1_rate)  . "\n";
        $mutation_string .= "\t$tier1_mutations";
        $rate_string .= "\t$tier1_rate";
        
        $total_mutations = $tier1_mutations;        
        $total_bases = $tier1_bases;

        my $tier2_mutations = my $tier3_mutations = 0;

        if($self->tier2_file)
        {
            $tier2_mutations = parse_mutations_file($self->tier2_file);
            $total_bases += $tier2_bases;
            $total_mutations += $tier2_mutations;
            my $tier2_rate = ($tier2_mutations / $coverage_factor) / ($tier2_bases / 1000000);
#            print join("\t", "2", $tier2_mutations, commify($tier2_bases), $tier2_rate)  . "\n";
            $mutation_string .= "\t$tier2_mutations";
            $rate_string .= "\t$tier2_rate";

        }

        if($self->tier3_file)
        {
            $tier3_mutations = parse_mutations_file($self->tier3_file);
            $total_mutations += $tier3_mutations;
            $total_bases += $tier3_bases;
            my $tier3_rate = ($tier3_mutations / $coverage_factor) / ($tier3_bases / 1000000);
#            print join("\t", "3", $tier3_mutations, commify($tier3_bases), $tier3_rate)  . "\n";
            $mutation_string .= "\t$tier3_mutations";
            $rate_string .= "\t$tier3_rate";

        }

        ## Calculate non-tier1 rate ##
        
        my $non_tier1_mutations = $tier2_mutations + $tier3_mutations;
        my $non_tier1_rate = ($non_tier1_mutations / $coverage_factor) / ($non_tier1_bases / 1000000);
        $mutation_string .= "\t$non_tier1_mutations";
        $rate_string .= "\t$non_tier1_rate";

        ## Calculate overall mutation rate ##

        my $overall_rate = ($total_mutations / $coverage_factor) / ($total_bases / 1000000);
        $mutation_string .= "\t$total_mutations";
        $rate_string .= "\t$overall_rate";


#        print join("\t", "ALL", $total_mutations, commify($total_bases), $overall_rate)  . "\n";
        print $self->sample_name . $mutation_string . $rate_string . "\n";

	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}




#############################################################
# ParseFile - takes input file and parses it
#
#############################################################

sub parse_mutations_file
{
	(my $mutation_file) = @_;
	
        if(!(-e $mutation_file))
        {
            die "Mutation file not found: $mutation_file\n";
        }
        
        my $num_mutations = 0;
        
	## Parse the file ##

	my $input = new FileHandle ($mutation_file);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;

                my ($chrom, $chr_start, $chr_stop, $ref, $var) = split(/\t/, $line);
                
                if($chrom && lc($chrom) ne "chr" && lc($chrom) ne "chrom")
                {
                    $num_mutations++;   
                }

	}
	
	close($input);

	return($num_mutations);
	
}



#############################################################
# ParseFile - takes input file and parses it
#
#############################################################

sub parse_regions_file
{
	(my $regions_file) = @_;
	
        if(!(-e $regions_file))
        {
            die "Regions file not found: $regions_file\n";
        }
        
        my $num_bases = 0;
        
	## Parse the file ##

	my $input = new FileHandle ($regions_file);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;

                my ($chrom, $chr_start, $chr_stop) = split(/\t/, $line);
                
                if($chrom && lc($chrom) ne "chr" && lc($chrom) ne "chrom")
                {
                    for(my $position = $chr_start; $position <= $chr_stop; $position++)
                    {
                        $num_bases++;
                    }
                }

	}
	
	close($input);

	return($num_bases);
	
}


sub commify
{
	local($_) = shift;
	1 while s/^(-?\d+)(\d{3})/$1,$2/;
	return $_;
}

1;
