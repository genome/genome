
package Genome::Model::Tools::Varscan::CopyNumberPlots;     # rename this when you give the module file a different name <--

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
my $min_loh_size = 10;
my $min_loh_snps = 3;
my $num_loh_regions = 0;

class Genome::Model::Tools::Varscan::CopyNumberPlots {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		regions_file		=> { is => 'Text', doc => "Path to copy number regions from Varscan copyCaller", is_optional => 0 },
		output_basename 	=> { is => 'Text', doc => "Output file basename for cnv plots", is_optional => 0 },
		min_points_to_plot 	=> { is => 'Text', doc => "Minimum number of points for a chromosome to plot it", is_optional => 0, default => 100 },
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Generate plots of exome copy number from Varscan copyCaller calls"                 
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
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		(my $chrom) = split(/\t/, $line);
		
		if($current_chrom && $chrom ne $current_chrom)
		{
			process_results($self, $current_chrom, $current_chrom_results);
			$current_chrom_results = "";	
		}

		$current_chrom = $chrom;		
		$current_chrom_results .= "\n" if($current_chrom_results);
		$current_chrom_results .= $line;
	}
	
	close($input);
	
	process_results($self, $current_chrom, $current_chrom_results);


	open(INDEX, ">$output_basename.index.html") or die "Can't open outfile: $!\n";
	print INDEX "<HTML><BODY><TABLE CELLSPACING=0 CELLPADDING=5 BORDER=0 WIDTH=\"100%\">\n";
	print INDEX "<TR>\n";

	for(my $chrom = 1; $chrom <= 24; $chrom++)
	{
		my $chrom_name = $chrom;
		$chrom_name = "X" if($chrom == 23);
		$chrom_name = "Y" if($chrom == 24);

		my $image_filename = $image_basename . "." . $chrom_name . ".jpg";
		print INDEX "<TD><A HREF=\"$image_filename\"><IMG SRC=\"$image_filename\" HEIGHT=240 WIDTH=320 BORDER=0></A></TD>\n";

		$num_printed_in_column++;

		if($num_printed_in_column >= 4)
		{
			print INDEX "</TR><TR>\n";
			$num_printed_in_column = 0;
		}


	}

	print INDEX "</TR></TABLE></BODY></HTML>\n";
	close(INDEX);
	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
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

