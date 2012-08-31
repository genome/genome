
package Genome::Model::Tools::Bowtie::MatchToRegions;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# MatchToRegions.pm - 	Get unmapped/poorly-mapped reads by model id
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

class Genome::Model::Tools::Bowtie::MatchToRegions {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		alignments_file	=> { is => 'Text', doc => "File containing Bowtie output" },
		regions_file	=> { is => 'Text', doc => "BED-format file of target regions" },
		output_file	=> { is => 'Text', doc => "Output matching alignments to file", is_optional => 1 },
                output_layers	=> { is => 'Text', doc => "Output RefCov layers to file", is_optional => 1 },
                verbose	=> { is => 'Text', doc => "Print interim progress updates", is_optional => 1 },
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Match Bowtie alignments to a regions.tsv file, build RefCov layers"                 
}

sub help_synopsis {
    return <<EOS
This command matches genomic alignments to a list of target regions
EXAMPLE:	gmt bowtie match-to-regions --alignments myfile.bowtie --regions myregions.tsv --output-file myfile.regions.bowtie
EOS
}

sub help_detail {                           # this is what the user will see with the longer version of help. <---
    return <<EOS 

EOS
}


my %regions = ();

################################################################################################
# Execute - the main program logic
#
################################################################################################

sub execute {                               # replace with real execution logic.
	my $self = shift;

	## Get required parameters ##
	my $alignments_file = $self->alignments_file;
	my $regions_file = $self->regions_file;
	my $output_file = $self->output_file;
	my $output_layers = $self->output_layers;
        my $verbose = 1 if(defined($self->verbose));

	my %stats = ();

	if($output_file)
	{
		open(OUTFILE, ">$output_file") or die "Can't open outfile: $!\n";
	}
        
        if($output_layers)
        {
                open(LAYERS, ">$output_layers") or die "Can't open outfile: $!\n";
        }


	## Parse the Regions file ##

#	my %regions = ();
	my $input = new FileHandle ($regions_file);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		(my $chromosome, my $chr_start, my $chr_stop, my $region_name) = split(/\t/, $line);
	#	$regions{$region_name} = "$chromosome\t$chr_start\t$chr_stop";
#		$regions{$chromosome} = "" if(!$regions{$chromosome});
#		$regions{$chromosome} .= "$chr_start\t$chr_stop\t$region_name\n";
		for(my $position_key = substr($chr_start, 0, 3); $position_key <= substr($chr_stop, 0, 3); $position_key++)
		{
			my $region_key = $chromosome . ":" . $position_key;
			$regions{$region_key} = "" if(!$regions{$region_key});
			$regions{$region_key} .= "$chr_start\t$chr_stop\t$region_name\n";
		}
	}
	
	close($input);

	print "$lineCounter regions loaded\n";

	## Parse the Alignments ##

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
                        
                        if($verbose)
                        {
                                print $stats{'within_regions'} . " of " . $stats{'num_alignments'} . " alignments within targets...\n" if(!($stats{'num_alignments'} % 10000))
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

			my $region_layers = within_regions($read_name, $chromosome, $chr_start, $chr_stop);
			
			if($region_layers)
			{
				$stats{'within_regions'}++;
				print LAYERS "$region_layers" if($output_layers);
				print OUTFILE "$line\n" if($output_file);
			}

#			return(0) if($lineCounter > 10);
		}
	}
	
	close($input);

	print "$stats{'num_alignments'} reads aligned\n";
	
	if($stats{'num_alignments'})
	{
		$stats{'within_regions'} = 0 if(!$stats{'within_regions'});
		$stats{'pct_within'} = $stats{'within_regions'} / $stats{'num_alignments'} * 100;
		$stats{'pct_within'} = sprintf("%.2f", $stats{'pct_within'}) . '%';
		print $stats{'within_regions'} . " " . $stats{'pct_within'} . " were within regions\n";
	}
#	close(OUTFILE) if(OUTFILE);
#	close(LAYERS) if(LAYERS);
	
	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}



################################################################################################
# Within target regions - check if alignment is within one of our regions
#
################################################################################################

sub within_regions {
	
	(my $read, my $chrom, my $start, my $stop) = @_;
	
#	foreach my $region (keys %regions)
#	{
#		(my $region_chrom, my $region_start, my $region_stop) = split(/\t/, $regions{$region});
#		if($region_chrom eq $chrom)
#		{
#			if($position >= $region_start && $position <= $region_stop)
#			{
#				return(1);
#			}
#		}
#	}

#	if($regions{$chrom})
#	{
#		my @region_list = split(/\n/, $regions{$chrom});
#		
#		foreach my $region (@region_list)
#		{
#			(my $region_start, my $region_stop, my $region_name) = split(/\t/, $region);
#			if($position >= $region_start && $position <= $region_stop)
#			{
#				return($region);
#			}
#		}
#	}

	my $num_regions_overlapped = 0;
	my $region_layers = "";

	my $chrom_key = $chrom . ":" . substr($start, 0, 3);

	if($regions{$chrom_key})
	{
		my @region_list = split(/\n/, $regions{$chrom_key});
		
		foreach my $region (@region_list)
		{
			(my $region_start, my $region_stop, my $region_name) = split(/\t/, $region);
			if(($start >= $region_start && $start <= $region_stop) || ($stop >= $region_start && $stop <= $region_stop))
			{
				my $new_start = $start - $region_start + 1;
				my $new_stop = $stop - $region_start + 1;
				$new_start = 1 if($new_start < 1);
				$region_layers .= "$read\t$new_start\t$new_stop\t$region_name\n";
			}
		}
	}

	return($region_layers);
}


1;

