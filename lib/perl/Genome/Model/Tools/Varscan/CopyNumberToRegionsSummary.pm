
package Genome::Model::Tools::Varscan::CopyNumberToRegionsSummary;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# RunVarscan - Run Varscan somatic on two BAM files.
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@genome.wustl.edu)
#
#	CREATED:	12/09/2009 by D.K.
#	MODIFIED:	12/09/2009 by D.K.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;

use FileHandle;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

## SET DEFAULT PARAMS ##
my $undo_sd = 2;

class Genome::Model::Tools::Varscan::CopyNumberToRegionsSummary {
    is => 'Command',

    has => [                                # specify the command's single-value properties (parameters) <--- 
        input_file => {
            is => 'Text',
            doc => "Input file of regions followed by copy number",
            is_optional => 0,
        },
        output_file => {
            is => 'Text',
            doc => "Output file with summarized results",
            is_optional => 0,
        },
    ],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Output VarScan copy number mean for regions (e.g. exons) in a BED file"                 
}

sub help_synopsis {
    return <<EOS
This command outputs VarScan copy number mean at regions (e.g. exons) in a BED file
EXAMPLE:	gmt varscan copy-number-to-regions --bed-file myExons.bed --sample-dir myDir/ --output-file Sample.exons.copynumber ...
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
	my $input_file = $self->input_file;
	my $output_file = $self->output_file;

	my $amp_threshold = 1;
	my $del_threshold = -1;
	my $low_amp_threshold = 0.25;
	my $low_del_threshold = -0.25;


	my %gene_chrom = my %gene_start = my %gene_stop = ();
	my %gene_seg_mean_sum = my %gene_seg_mean_num = my %gene_amp = my %gene_del = my %gene_neutral = my %gene_na = ();
	my %gene_low_amp = my %gene_low_del = ();

	open(OUTFILE, ">$output_file") or die "Can't open outfile: $!\n";
	
	my $input = new FileHandle ($input_file);
	my $lineCounter = 0;

	my %lines_by_gene = ();
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		
		
		if($lineCounter > 1)
		{
			my ($chrom, $chr_start, $chr_stop, $name, $seg_mean) = split(/\s+/, $line);
			
#			if(!$ARGV[0]$name eq $ARGV[0] || $ARGV[0] eq 'all')
#			{
				$lines_by_gene{$name} .= "\n" if($lines_by_gene{$name});
				$lines_by_gene{$name} .= join("\t", $chrom, $chr_start, $chr_stop, $seg_mean);				
#			}

		}

	}
	
	close($input);	

	## Run through all of the exons seen for this gene ##

	foreach my $gene (keys %lines_by_gene)
	{
		my $seg_mean_sum = my $seg_mean_num = my $size_sum = 0;
		my @lines = split(/\n/, $lines_by_gene{$gene});
		foreach my $line (@lines)
		{
			my ($chrom, $chr_start, $chr_stop, $seg_mean) = split(/\t/, $line);
			my $size = $chr_stop - $chr_start + 1;
			$gene_chrom{$gene} = $chrom;
			$gene_start{$gene} = $chr_start if(!$gene_start{$gene} || $chr_start < $gene_start{$gene});
			$gene_stop{$gene} = $chr_stop if(!$gene_stop{$gene} || $chr_stop > $gene_stop{$gene});

			if($seg_mean ne "NA")
			{
				$seg_mean_sum += ($seg_mean * $size);
				$size_sum += $size;
				$seg_mean_num++;
			}
		}

		my $avg_seg_mean = "NA";

#		if($seg_mean_num)
		if($size_sum)
		{
#			$avg_seg_mean = $seg_mean_sum / $seg_mean_num;
			$avg_seg_mean = $seg_mean_sum / $size_sum;
			
			## Save the sample seg mean for this gene to the gene's total ##
			$gene_seg_mean_sum{$gene} += $avg_seg_mean;
			$gene_seg_mean_num{$gene}++;
			
			if($avg_seg_mean >= $amp_threshold)
			{
				$gene_amp{$gene}++;
			}
			elsif($avg_seg_mean >= $low_amp_threshold)
			{
				$gene_low_amp{$gene}++;	
			}
			elsif($avg_seg_mean <= $del_threshold)
			{
				$gene_del{$gene}++;
			}
			elsif($avg_seg_mean <= $low_del_threshold)
			{
				$gene_low_del{$gene}++;
			}
			else
			{
				$gene_neutral{$gene}++;
			}
		}
		else
		{
			$gene_na{$gene}++;
		}
		
		print OUTFILE join("\t", $gene_chrom{$gene}, $gene_start{$gene}, $gene_stop{$gene}, $gene, $seg_mean_num, $avg_seg_mean) . "\n";

	}

	close(OUTFILE);
	

}





#############################################################
# parse_array - parses the file
#
#############################################################

sub get_copy_result
{
	my ($chrom, $chr_start, $chr_stop, $events) = @_;
	
	my $size = $chr_stop - $chr_start + 1;
	my $best_overlap = my $best_seg_mean = "";
	my $best_event = "";
	if($events)
	{
		my @events = split(/\n/, $events);		
		my $overlap = 0;

		foreach my $check_event (@events)
		{
			my ($event_chrom, $event_chr_start, $event_chr_stop, $event_seg_mean) = split(/\t/, $check_event);
			my $event_size = $event_chr_stop - $event_chr_start + 1;
			if($event_chrom eq $chrom)
			{
				if($event_chr_stop >= $chr_start && $event_chr_start <= $chr_stop)
				{
					## Determine the overlap ##

					if($event_chr_start >= $chr_start && $event_chr_stop <= $chr_stop)
					{
						## Check event completely within sought event #
						$overlap = $event_chr_stop - $event_chr_start + 1;
					}
					elsif($chr_start >= $event_chr_start && $chr_stop <= $event_chr_stop)
					{
						## Sought event completely within check event #
						$overlap = $chr_stop - $chr_start + 1;
					}
					elsif($event_chr_start < $chr_start && $event_chr_stop >= $chr_start)
					{
						## Check event overlaps at front edge ##
						$overlap = $event_chr_stop - $chr_start;
					}
					elsif($event_chr_stop > $chr_stop && $event_chr_start <= $chr_stop)
					{
						## Check event overlaps at back edge ##
						$overlap = $chr_stop - $event_chr_start;
					}

					if(!$best_overlap || $overlap >= $best_overlap)
					{
						$best_seg_mean = $event_seg_mean;
						$best_overlap = $overlap;
						$best_event = $check_event;
					}
				}
			}
		}
		

		
	}
	
	if($best_overlap)
	{
		return($best_seg_mean . "\t" . $best_event);		
	}
	else
	{
		return();
	}

	
}



#############################################################
# load_cds_exons - load the coordiantes of the BED file
#
#############################################################

sub load_cds_exons
{
	my $FileName = shift(@_);

	my %exons = ();

	my $input = new FileHandle ($FileName);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		my ($chrom, $chr_start, $chr_stop, $gene) = split(/\t/, $line);
		
		$exons{$chrom} .= "\n" if($exons{$chrom});
		$exons{$chrom} .= "$line";
	}
	
	close($input);

	return(%exons);
}




#############################################################
# load_regions - load a sample's segmented CBS calls
#
#############################################################

sub load_regions
{
	my $FileName = shift(@_);

	my $regions = ();

	my $input = new FileHandle ($FileName);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		
		
		if($lineCounter > 1)
		{
			my ($num, $id, $chrom, $chr_start, $chr_stop, $num_mark, $seg_mean, $bstat, $pval, $lcl, $ucl) = split(/\s+/, $line);
			$chrom =~ s/[^0-9XYMT]//g;
			my $key = join("\t", $chrom, $chr_start, $chr_stop);
			$regions .= "\n" if($regions);
			$regions .= join("\t", $chrom, $chr_start, $chr_stop, $seg_mean);
		}
	}
	
	close($input);

	return($regions);	
}



#############################################################
# get_copy_result - retrieve a sample's CBS result
#
#############################################################

sub old_get_copy_result
{
	my ($chrom, $chr_start, $regions) = @_;
	
	my @regions = split(/\n/, $regions);
	
	foreach my $region (@regions)
	{
		my ($region_chrom, $region_chr_start, $region_chr_stop, $region_seg_mean) = split(/\t/, $region);
		
		if($region_chrom eq $chrom && $region_chr_stop >= $chr_start && $region_chr_start <= $chr_start)
		{
			return($region_seg_mean);
		}
	}
	
	return("");
}



1;

