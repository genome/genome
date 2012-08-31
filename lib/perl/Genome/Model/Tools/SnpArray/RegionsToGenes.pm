
package Genome::Model::Tools::SnpArray::RegionsToGenes;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# RegionsToGenes - retrieves the symbol(s) of Known Genes for a region based on coordinates using UCSC KnownGenes
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@watson.wustl.edu)
#
#	CREATED:	04/01/2009 by D.K.
#	MODIFIED:	04/01/2009 by D.K.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;

use FileHandle;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

my %stats = ();

class Genome::Model::Tools::SnpArray::RegionsToGenes {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		regions_file	=> { is => 'Text', doc => "Regions file in chromosome start stop format", is_optional => 0, is_input => 1 },
		ucsc_refgene	=> { is => 'Text', doc => "Path to the UCSC refGene.txt file", is_optional => 0, is_input => 1, default => '/gscmnt/sata809/info/medseq/dkoboldt/SNPseek2/ucsc/refGene.txt'},
		output_file	=> { is => 'Text', doc => "Output file", is_optional => 1, is_input => 1},
		output_gene_counts	=> { is => 'Text', doc => "Output file for the counts of affected regions by gene", is_optional => 1, is_input => 1},
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Appends the symbol(s) of known genes to a list of chromosomal regions"                 
}

sub help_synopsis {
    return <<EOS
This command appends the symbol(s) of known genes to a list of chromosomal regions
EXAMPLE:	gmt snp-array regions-to-genes --regions-file recurrent_cnas/recurrent.focal.amplification.1.tsv
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

sub execute
{                               # replace with real execution logic.
	my $self = shift;

	## Get required parameters ##
	my $regions_file = $self->regions_file;
	my $ucsc_refgene = $self->ucsc_refgene;
	my $output_file = $self->output_file;
	my $output_gene_counts = $self->output_gene_counts;

	my %genes_by_chrom = load_genes($ucsc_refgene);


	## Open output file ##
	
	if($output_file)
	{
		open(OUTFILE, ">$output_file") or die "can't open outfile: $!\n";
	}


	my %max_affected = ();

	## Parse the regions file ##

	my $input = new FileHandle ($regions_file);
	my $lineCounter = 0;

	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;

		if($lineCounter >= 1)
		{
			my ($chrom, $chr_start, $chr_stop, $num_hits) = split(/\t/, $line);
			
			if($chrom && lc(substr($chrom, 0, 5)) ne "chrom")
			{
				my $gene_list = "";
				## See if we have genes for this chromosome ##
				if($genes_by_chrom{$chrom})
				{
#					if($lineCounter < 100)
#					{						
						$gene_list = get_matching_genes($chrom, $chr_start, $chr_stop, $genes_by_chrom{$chrom});
						
						if($gene_list)
						{
							## Parse every gene symbol from the list ##
							my @genes = split(/\,/, $gene_list);
							
							foreach my $gene (@genes)
							{
								$max_affected{$gene} = $num_hits if(!$max_affected{$gene} || $num_hits > $max_affected{$gene});
							}
						}


#						print join("\t", $num_hits, $chrom, $chr_start, $chr_stop, $gene_list) . "\n";
#					}					
				}

				print OUTFILE join("\t", $line, $gene_list) . "\n" if($output_file);
			}
			elsif($chrom)
			{
				print OUTFILE "$line\tknown_gene(s)\n"; 
			}
		}


	}
	
	close($input);	
	close(OUTFILE);
	
	
	## Print out all genes ##
	
	if($output_gene_counts)
	{
		open (OUTFILE, ">$output_gene_counts") or die "Can't open outfile: $!\n";

		foreach my $gene (sort keys %max_affected)
		{
			print OUTFILE join("\t", $max_affected{$gene}, $gene) . "\n";
		}
	
		close(OUTFILE);
	}
	



}




################################################################################################
# Execute - the main program logic
#
################################################################################################

sub get_matching_genes
{
	my ($chrom, $chr_start, $chr_stop, $gene_list) = @_;
	
	my @genes = split(/\n/, $gene_list);
	
	my %gene_included = ();
	my $matching_genes = "";
	
	foreach my $gene_entry (@genes)
	{
		my ($gene_start, $gene_stop, $gene_name) = split(/\t/, $gene_entry);
		
		if(!$gene_included{$gene_name} && $gene_stop >= $chr_start && $gene_start <= $chr_stop)
		{
			$matching_genes .= "," if($matching_genes);
			$matching_genes .= $gene_name;
			$gene_included{$gene_name} = 1;
		}
	}
	
	return($matching_genes);
}



################################################################################################
# Execute - the main program logic
#
################################################################################################

sub load_genes
{
	my $FileName = shift(@_);
	my %genes_by_chrom = ();

	my $input = new FileHandle ($FileName);
	my $lineCounter = 0;

	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;

		if($lineCounter >= 1)
		{
			my @lineContents = split(/\t/, $line);
			my $chrom = $lineContents[2];
			my $tx_start = $lineContents[4];
			my $tx_stop = $lineContents[5];
			my $gene_symbol = $lineContents[12];
			
			$chrom =~ s/chr//g;
			$tx_start++;
			$tx_stop++;
			
			if($genes_by_chrom{$chrom})
			{
				$genes_by_chrom{$chrom} .= "\n";
			}
			
			$genes_by_chrom{$chrom} .= join("\t", $tx_start, $tx_stop, $gene_symbol);
#			print join("\t", $chrom, $tx_start, $tx_stop, $gene_symbol) . "\n";
		}


	}
	
	close($input);	
	
	return(%genes_by_chrom);
}




################################################################################################
# Execute - the main program logic
#
################################################################################################

sub run_recurrence
{
	my $self = shift(@_);
	my $chrom = shift(@_);
	my $event_size = shift(@_);
	my $event_type = shift(@_);

	## Get required parameters ##
	my $sample_events_file = $self->sample_events_file;
	my $output_basename = $self->output_basename;


	my %events_by_position = ();
	
	my $input = new FileHandle ($sample_events_file);
	my $lineCounter = 0;

	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;

		if($lineCounter <= 20)#($lineCounter > 0)#
		{
			my ($sample_name, $file_name) = split(/\t/, $line);


			my $events = load_events($file_name, $chrom, $event_size, $event_type);

			## Save the positions of this sample's events ##

			my @events = split(/\n/, $events);
			my $num_events = @events;
			
			foreach my $event (@events)
			{
				my ($chrom, $chr_start, $chr_stop) = split(/\t/, $event);
				
				## Count every base ##
				
				for(my $position = $chr_start; $position <= $chr_stop; $position++)
				{
					$events_by_position{$position}++;
				}
			}
			
			print "$chrom\t$sample_name\t$num_events\n";				
		}


	}
	
	close($input);
	
	
	## Look for regions of high recurrence ##
	
	my $outfile_name = $output_basename . "recurrent." . $event_size . "-" . $event_type . ".$chrom.tsv";
	open(OUTFILE, ">$outfile_name") or die "can't open outfile: $!\n";
	
	my $region_start = my $region_stop = my $region_events = 0;
	
	foreach my $position (sort numerically keys %events_by_position)
	{
		if(!$region_start)
		{
			## Start A new Region ##
			$region_start = $region_stop = $position;
			$region_events = $events_by_position{$position};
		}
		elsif($region_stop && ($position - $region_stop) == 1 && $events_by_position{$position} == $region_events)
		{
			## Extend the current region ##
			$region_stop = $position;
		}
		else
		{
			## End the current region because it's >1bp apart or has a different depth ##
			
			print OUTFILE join("\t", $chrom, $region_start, $region_stop, $region_events) . "\n";
			## Start A new Region ##
			$region_start = $region_stop = $position;
			$region_events = $events_by_position{$position};
		}
	}
	
	
	close(OUTFILE);	
}




sub numerically
{
	$a <=> $b;
}


################################################################################################
# Execute - the main program logic
#
################################################################################################

sub load_events
{
	my $FileName = shift(@_);
	my $desired_chrom = shift(@_);
	my $desired_class = shift(@_);
	my $desired_type = shift(@_);
	my $events = "";

	my $input = new FileHandle ($FileName);
	my $lineCounter = 0;

	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;

		if($lineCounter > 1)
		{
			my ($chrom, $chr_start, $chr_stop, $seg_mean, $num_segments, $num_markers, $p_value, $event_type, $event_size, $size_class) = split(/\t/, $_);

			if($chrom eq $desired_chrom)
			{
				if($size_class eq $desired_class)
				{
					if($event_type eq $desired_type)
					{
						$events .= $line . "\n";
					}
				}
			}
		}


	}
	
	close($input);	
	
	return($events);
}



1;


