
package Genome::Model::Tools::Varscan::CopyNumberSegments;     # rename this when you give the module file a different name <--

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

class Genome::Model::Tools::Varscan::CopyNumberSegments {
    is => 'Command',

    has => [                                # specify the command's single-value properties (parameters) <--- 
        regions_file => {
            is => 'Text',
            doc => "Path to copy number regions from Varscan copyCaller",
            is_optional => 0,
        },
        output_basename  => {
            is => 'Text',
            doc => "Output file basename for cnv plots",
            is_optional => 0,
        },
        min_depth  => {
            is => 'Text',
            doc => "Minimum depth for a region (in one sample) to include it",
            is_optional => 0,
            default => 8,
        },
        min_points_to_plot  => {
            is => 'Text',
            doc => "Minimum number of points for a chromosome to plot it",
            is_optional => 0,
            default => 100,
        },
        undo_sd  => {
            is => 'Text',
            doc => "Remove change-points of less than this number of standard deviations",
            is_optional => 0,
            default => 4,
        },
        lsf_command  => {
            is => 'Text',
            doc => "If set to something like bsub -q $ENV{GENOME_LSF_QUEUE_BUILD_WORKER}, will run bsub",
            is_optional => 1,
        },
        array_data  => {
            is => 'Text',
            doc => "If set to 1, expect array data in id, chrom, pos, value format",
            is_optional => 1,
        },
        plot_y_min  => {
            is => 'Text',
            doc => "The minimum value on y-axis in CN plots",
            is_optional => 0,
            default => -5,
        },
        plot_y_max  => {
            is => 'Text',
            doc => "The minimum value on y-axis in CN plots",
            is_optional => 0,
            default => 5,
        },
    ],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Generate plots of exome copy number from Varscan copyCaller or SNP array calls"                 
}

sub help_synopsis {
    return <<EOS
Generate plots of exome copy number from Varscan copyCaller calls
EXAMPLE:	gmt capture copy-number-plots ...
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
	my $regions_file = $self->regions_file;
	my $output_basename = $self->output_basename;
	my $min_depth = $self->min_depth;
	$undo_sd = $self->undo_sd;

	## Open the index HTML file ##
	my @tempArray = split(/\//, $output_basename);
	my $numElements = @tempArray;
	my $image_basename = $tempArray[$numElements - 1];
	

	my $num_printed_in_column = 0;

	## Parse the variant file ##

	my $input = new FileHandle ($regions_file);
	my $lineCounter = 0;
	
	my $current_chrom = "";
	my $current_chrom_results = "";
	my $metMinDepth = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		my @lineContents = split(/\t/, $line);
		my $numContents = @lineContents;
		
		my ($chrom, $chr_start, $chr_stop, $normal, $tumor, $log_value) = split(/\t/, $line);
		## Parse newer output ##
		($chrom, $chr_start, $chr_stop, my $num_positions, $normal, $tumor, $log_value) = split(/\t/, $line) if($numContents > 6);
		
		
		if($lineCounter <= 2 && $line =~ 'REF')
		{
			# Skip SNP array headers #
		}
		elsif($lineCounter > 1 || $chrom ne "chrom")
		{
			if($self->array_data)
			{
				(my $id, $chrom, $chr_start, $log_value) = split(/\t/, $line);
				$chr_stop = $chr_start + 1;
				$num_positions = 1;
				$normal = $min_depth + 1;
				$tumor = $min_depth + 1;
			}


			## Process the previous chromosome ##
			if($current_chrom && $chrom ne $current_chrom)
			{
				print "Chromosome $chrom...\n";
				process_results($self, $current_chrom, $current_chrom_results);
				$current_chrom_results = "";	
			}
	
			if($normal >= $min_depth || $tumor >= $min_depth)
			{
				$metMinDepth++;
				$current_chrom = $chrom;		
				$current_chrom_results .= "\n" if($current_chrom_results);

# Current pipeline prior to 7/13/11: uses only start of region ##
#$current_chrom_results .= $line;

				## Determine region size. ##
				
				my $region_size = $chr_stop - $chr_start + 1;

				## If region size is less than 200 bp, report just the midpoint ##
				
				if($region_size <= 1000)
				{
					my $midpoint = sprintf("%d", ($chr_start + $chr_stop) / 2);

					if($midpoint >= $chr_start && $midpoint <= $chr_stop)
					{
						$current_chrom_results .= join("\t", $chrom, $chr_stop, $num_positions, $normal, $tumor, $log_value);
					}
					else
					{
						warn "no Midpoint $midpoint $chr_start $chr_stop\n";
					}
				}
				
				## Otherwise, report both start and stop ##
				
				else
				{
					$current_chrom_results .= join("\t", $chrom, $chr_start, $num_positions, $normal, $tumor, $log_value) . "\n";
					## Add the stop position ##
					$current_chrom_results .= join("\t", $chrom, $chr_stop, $num_positions, $normal, $tumor, $log_value);					
				}




			}

		}
		

	}
	
	close($input);
	
	process_results($self, $current_chrom, $current_chrom_results);
	
	print "$lineCounter lines parsed from Varscan output\n";
	print "$metMinDepth met minimum depth of $min_depth\n";


	open(SEGMENTS, ">$output_basename.segments.tsv") or die "Can't open outfile: $!\n";
	print SEGMENTS join("\t", "ID", "sample", "chrom", "loc.start", "loc.end", "num.mark", "seg.mean", "bstat", "pval", "lcl", "ucl") . "\n";

	open(INDEX, ">$output_basename.index.html") or die "Can't open outfile: $!\n";
	print INDEX "<HTML><BODY><TABLE CELLSPACING=0 CELLPADDING=5 BORDER=0 WIDTH=\"100%\">\n";
	print INDEX "<TR>\n";

	for(my $chrom = 1; $chrom <= 24; $chrom++)
	{
		my $chrom_name = $chrom;
		$chrom_name = "X" if($chrom == 23);
		$chrom_name = "Y" if($chrom == 24);

		my $chrom_filename = $output_basename . ".$chrom_name.infile";
		my $segments_filename = "$chrom_filename.segments.p_value";
		my $image_filename = $image_basename . "." . $chrom_name . ".jpg";
		print INDEX "<TD><A HREF=\"$image_filename\"><IMG SRC=\"$image_filename\" HEIGHT=240 WIDTH=320 BORDER=0></A></TD>\n";

		print SEGMENTS parse_segments($segments_filename) if(-e $segments_filename);

		$num_printed_in_column++;

		if($num_printed_in_column >= 4)
		{
			print INDEX "</TR><TR>\n";
			$num_printed_in_column = 0;
		}


	}

	print INDEX "</TR></TABLE></BODY></HTML>\n";
	close(INDEX);

	close(SEGMENTS);
	
	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}




################################################################################################
# Execute - the main program logic
#
################################################################################################

sub parse_segments
{                               # replace with real execution logic.
	my $segments_file = shift(@_);

	my $result = "";

	my $input = new FileHandle ($segments_file);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		if($lineCounter > 1)
		{
			my @lineContents = split(/\s+/, $line);
			$line =~ s/\"X\"/X/;
			$line =~ s/\"Y\"/Y/;
			$result .= join("\t", @lineContents) . "\n";
		}
	}
	
	close($input);
	
	return($result);
}


################################################################################################
# Process results - filter variants by type and into high/low confidence
#
################################################################################################

sub process_results
{
	my ($self, $chrom, $lines) = @_;

	my $output_basename = $self->output_basename;

	my $chrom_filename = $output_basename . ".$chrom.infile";
	my $script_filename = $output_basename . ".$chrom.R";
	my $image_filename = $output_basename . "." . $chrom . ".jpg";

	my @lines = split(/\n/, $lines);
	my $num_lines = @lines;

	if($num_lines >= $self->min_points_to_plot)
	{
		print "CHROMOSOME $chrom: $num_lines lines\n";

		## Print all lines to file ##
		open(OUTFILE, ">$chrom_filename") or die "Can't open outfile $chrom_filename: $!\n";
	
		my $num_columns = 0;
		foreach my $line (@lines)
		{		
			if($line =~ 'Infinity' || $line =~ '∞')
			{
				## Don't try to plot sites with log2 infinity 	
			}
			else
			{
				my @lineContents = split(/\t/, $line);
				$num_columns = @lineContents;
				print OUTFILE "$line\n";			
			}
	
		}
	
		close(OUTFILE);
		

		## Begin R Script ##
	
		open(SCRIPT, ">$script_filename") or die "Can't open script $script_filename: $!\n";	
		print SCRIPT "library(DNAcopy)\n";
		print SCRIPT "regions <- read.table(\"$chrom_filename\")\n";
		print SCRIPT "png(\"$image_filename\", height=600, width=800)\n";
		print SCRIPT qq{
CNA.object <- CNA(regions\$V$num_columns, regions\$V1, regions\$V2, data.type="logratio", sampleid=c("Chromosome $chrom"))\n
smoothed.CNA.object <- smooth.CNA(CNA.object)\n
};

		print SCRIPT qq{
segment.smoothed.CNA.object <- segment(smoothed.CNA.object, undo.splits="sdundo", undo.SD=$undo_sd, verbose=1)
p.segment.smoothed.CNA.object <- segments.p(segment.smoothed.CNA.object)
x <- c($undo_sd, length(p.segment.smoothed.CNA.object\$pval[p.segment.smoothed.CNA.object\$num.mark>=20]))
write.table(x, file="$chrom_filename.segments.sd", append=TRUE)
};

## Allow multiple runs for multiple SDs ##
my $this_undo_sd = $undo_sd;
while($this_undo_sd > 0.5)
{
	$this_undo_sd -= 0.5;
	print SCRIPT qq|
if(length(p.segment.smoothed.CNA.object\$pval[p.segment.smoothed.CNA.object\$num.mark>=20]) < 50)
{
x <- c($this_undo_sd, length(p.segment.smoothed.CNA.object\$pval[p.segment.smoothed.CNA.object\$num.mark>=20]))
write.table(x, file="$chrom_filename.segments.sd", append=TRUE)
segment.smoothed.CNA.object <- segment(smoothed.CNA.object, undo.splits="sdundo", undo.SD=$this_undo_sd, verbose=1)
p.segment.smoothed.CNA.object <- segments.p(segment.smoothed.CNA.object)	
}
|;	
}




## Original using their plotting mechanism ##

#		print SCRIPT qq{
#plot(segment.smoothed.CNA.object, type="w", cex=0.5, cex.axis=1.5, cex.lab=1.5, ylim=c(-4,4))
#write.table(p.segment.smoothed.CNA.object, file="$chrom_filename.segments.p_value")
#};
		my $ymin = $self->plot_y_min;
		my $ymax = $self->plot_y_max;
		print SCRIPT qq{
detach(package:DNAcopy)
par(mar=c(4,4,2,2))
plot(segment.smoothed.CNA.object\$data\$maploc, segment.smoothed.CNA.object\$data\$Chromosome.$chrom, pch=19, cex=0.25, cex.axis=1.25, cex.lab=1.5, col="cornflowerblue", ylim=c($ymin,$ymax), main="Chromosome $chrom", xlab="Position", ylab="Copy Number Change (log2)")
segments(segment.smoothed.CNA.object\$output\$loc.start, segment.smoothed.CNA.object\$output\$seg.mean, segment.smoothed.CNA.object\$output\$loc.end, segment.smoothed.CNA.object\$output\$seg.mean, col="red", lwd=2)
write.table(p.segment.smoothed.CNA.object, file="$chrom_filename.segments.p_value")
};
		print SCRIPT "dev.off()\n";
		close(SCRIPT);
		
		if($self->lsf_command)
		{
#			print "Running $script_filename\n";
			system($self->lsf_command . " \"R --no-save < $script_filename\"");			
		}
		else
		{
			print "Running $script_filename\n";
			system("R --no-save < $script_filename");			
		}

		

	}


}


################################################################################################
# Process results - filter variants by type and into high/low confidence
#
################################################################################################

sub old_process_results
{
	my ($self, $chrom, $lines) = @_;

	my $output_basename = $self->output_basename;

	my $chrom_filename = $output_basename . ".$chrom.infile";
	my $script_filename = $output_basename . ".R";
	my $image_filename = $output_basename . "." . $chrom . ".jpg";

	my @lines = split(/\n/, $lines);
	my $num_lines = @lines;

	print "CHROMOSOME $chrom: $num_lines lines\n";

	## Print all lines to file ##
	open(OUTFILE, ">$chrom_filename") or die "Can't open outfile $chrom_filename: $!\n";

	my $num_columns = 0;
	foreach my $line (@lines)
	{		
		if($line =~ 'Infinity')
		{
			## Don't try to plot sites with log2 infinity 	
		}
		else
		{
			my @lineContents = split(/\t/, $line);
			$num_columns = @lineContents;
			print OUTFILE "$line\n";			
		}

	}

	close(OUTFILE);
	
	print "Chrom: $chrom lines $num_lines\n";

	if($num_lines >= $self->min_points_to_plot)
	{
		## Begin R Script ##
	
		open(SCRIPT, ">$script_filename") or die "Can't open script $script_filename: $!\n";
	
		print SCRIPT "regions <- read.table(\"$chrom_filename\")\n";
		print SCRIPT "png(\"$image_filename\", height=600, width=800)\n";
		print SCRIPT "plot(regions\$V2, regions\$V$num_columns, col=\"blue\", cex=0.5, cex.axis=1.5, cex.lab=1.5, pch=19, ylim=c(-4,4), main=\"Chromosome $chrom\", xlab=\"Position on chromosome $chrom\", ylab=\"Log2 Ratio (Tumor/Normal)\")\n";
		print SCRIPT "points(regions\$V3, regions\$V$num_columns, col=\"blue\", cex=0.5, cex.axis=1.5, cex.lab=1.5, pch=19, ylim=c(-4,4))\n";
		print SCRIPT "dev.off()\n";
		close(SCRIPT);
		
		print "Running $script_filename\n";
		system("R --no-save < $script_filename");		
	}

}


1;

