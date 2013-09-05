package Genome::Model::Tools::Analysis::Maf::Proximity;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# Proximity - Perform a proximity analysis on mutations in the MAF file.
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@watson.wustl.edu)
#
#	CREATED:	08/24/2010 by D.K.
#	MODIFIED:	08/24/2010 by D.K.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;

use FileHandle;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it
use Genome::Model::Tools::Analysis::Maf::Helpers qw/getMutationSamples byGeneTranscript load_aa_changes load_annotation/;

my %stats = ();
my $max_proximity = 0;

class Genome::Model::Tools::Analysis::Maf::Proximity {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		maf_file	=> { is => 'Text', doc => "Original MAF file" },
		annotation_file	=> { is => 'Text', doc => "Full annotation for variants in MAF file" },
		output_file	=> { is => 'Text', doc => "Output file for recurrence report", is_optional => 1 },
		output_maf	=> { is => 'Text', doc => "Output file of complete MAF with recurrent sites only", is_optional => 1 },		
		max_proximity	=> { is => 'Text', doc => "Maximum aa distance between mutations [10]", is_optional => 1, default => 10 },
		verbose		=> { is => 'Text', doc => "Print verbose output", is_optional => 1 },
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Performs a proximity analysis on mutations in a MAF file"                 
}

sub help_synopsis {
    return <<EOS
This command performs a proximity analysis on mutations in a MAF file
EXAMPLE:	gt analysis maf proximity --maf-file original.maf --output-file proximity-genes.tsv
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
	my $maf_file = $self->maf_file;
	$max_proximity = $self->max_proximity;# if(defined($self->max_proximity));

	if(!(-e $maf_file))
	{
		die "Error: MAF file not found!\n";
	}

	if($self->output_file)
	{
		open(OUTFILE, ">" . $self->output_file) or die "Can't open outfile: $!\n";
	}

	if($self->output_maf)
	{
		open(OUTMAF, ">" . $self->output_maf) or die "Can't open outfile: $!\n";
	}

	
	my %mutated_aa_positions = load_aa_changes($maf_file, $self->annotation_file, \%stats);

	$stats{'num_aa_positions'} = $stats{'num_recurrent_genes'} = $stats{'num_recurrent_positions'} = 0;
	

	## Declare hash to store min proximity for each mutation ##
	
	my %min_proximity_by_mutation = ();
	my %min_proximity_aa_key = ();

	foreach my $aa_key (sort byGeneTranscript keys %mutated_aa_positions)
	{
		$stats{'num_aa_positions'}++;

		(my $gene, my $transcript_name, my $aa_position) = split(/\t/, $aa_key);
		
		## Set the default min proximity to the max proximity plus 1 ##
		my $min_proximity = $max_proximity + 1;

		## Get Sample(s) with mutations at this position ##

		my $samples = getMutationSamples($mutated_aa_positions{$aa_key});
		my @samples = split(/\n/, $samples);
		my $num_samples = @samples;

		my @test = split(/\n/, $mutated_aa_positions{$aa_key});
		my $test_num = @test;

		## First, check for distance of zero ##
		
		if($num_samples > 1)
		{
			$min_proximity = 0.00;	
		}
		else
		{
			my $sample_name = $samples;
			## It's not recurrent at this particular residue, so look nearby for mutations not matching sample name##
			
			## Determine the minimum aa distance, from zero to maximum ##
	
			for(my $aa_distance = 1; $aa_distance <= $max_proximity; $aa_distance++)
			{
				## Check upstream ##
				my $distance_key = join("\t", $gene, $transcript_name, ($aa_position - $aa_distance));
				if($mutated_aa_positions{$distance_key})
				{
					my $nearby_samples = getMutationSamples($mutated_aa_positions{$distance_key});
					## If nearby mutation is from a different sample (or multiple ones), save it 
					if($nearby_samples ne $sample_name)
					{
						$min_proximity = $aa_distance;
					}
				}
				else
				{
					## Check downstream ##
					$distance_key = join("\t", $gene, $transcript_name, ($aa_position + $aa_distance));
					if($mutated_aa_positions{$distance_key})
					{
						my $nearby_samples = getMutationSamples($mutated_aa_positions{$distance_key});
						## If nearby mutation is from a different sample (or multiple ones), save it 
						if($nearby_samples ne $sample_name)
						{
							$min_proximity = $aa_distance;
						}
					}					
				}
			}
		}
		
		## Get each mutation line for this aa position ##
		
		my @aa_lines = split(/\n/, $mutated_aa_positions{$aa_key});
		my $num_lines = @aa_lines;

		## Save this proximity for each mutation ##
		
		foreach my $maf_line (@aa_lines)
		{
			if(!defined($min_proximity_by_mutation{$maf_line}) || $min_proximity < $min_proximity_by_mutation{$maf_line})
			{
				$min_proximity_by_mutation{$maf_line} = $min_proximity;
				$min_proximity_aa_key{$maf_line} = $aa_key;
			}
		}

	}

	## Go through each mutation in MAF to determine its distance ##
	my %proximity_counts = ();
	
	foreach my $maf_line (keys %min_proximity_by_mutation)
	{
		$stats{'num_mutations_examined'}++;
		my $proximity = $min_proximity_by_mutation{$maf_line};
		my $aa_key = $min_proximity_aa_key{$maf_line};
		if($proximity <= $max_proximity)
		{
			my @lineContents = split(/\t/, $maf_line);
			my $ref = $lineContents[10];
			my $var = $lineContents[11];
			$var = $lineContents[12] if($var eq $ref);
			
			
			print OUTFILE join("\t", $proximity, $aa_key, $lineContents[2], $lineContents[4], $lineContents[5], $lineContents[6], $ref, $var, $lineContents[15]) . "\n" if($self->output_file);
			print OUTMAF "$proximity\t$aa_key\t$maf_line\n" if($self->output_maf);
		}
		$proximity_counts{$proximity}++;
	}
	
	
	print $stats{'num_mutations_examined'} . " mutations evaluated\n";
	
	foreach my $proximity (sort keys %proximity_counts)
	{
		print $proximity_counts{$proximity} . " had a distance of $proximity\n";
	}

	close(OUTFILE) if($self->output_file);
	close(OUTMAF) if($self->output_maf);


}

1;
