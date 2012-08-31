
package Genome::Model::Tools::SnpArray::ProcessCmdsCalls;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# SearchRuns - Search the database for runs
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

class Genome::Model::Tools::SnpArray::ProcessCmdsCalls {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		map_file	=> { is => 'Text', doc => "Three-column file of genotype calls chrom, pos, genotype", is_optional => 0, is_input => 1 },
		cmds_file	=> { is => 'Text', doc => "Three-column file of genotype calls chrom, pos, genotype", is_optional => 0, is_input => 1 },
		num_sd	=> { is => 'Text', doc => "Minimum number of SD to call a change-point", is_optional => 1, is_input => 1, default => 1},
		output_basename	=> { is => 'Text', doc => "Output file for QC result", is_optional => 1, is_input => 1},
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Merges CMDS calls with map file and generates per-chromosome data"                 
}

sub help_synopsis {
    return <<EOS
This command merges CMDS calls with map file and generates per-chromosome data
EXAMPLE:	gmt snp-array process-cmds-calls --map-file /gscmnt/sata181/info/medseq/llin/Ovarian/SNP/CMDS_444_samples/map.csv --cmds-file /gscmnt/sata181/info/medseq/llin/Ovarian/SNP/CMDS_444_samples/NORMALIZATION/TCGA-10-0933.csv
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
	my $map_file = $self->map_file;
	my $cmds_file = $self->cmds_file;
	my $output_basename = $self->output_basename;
	my $num_sd = $self->num_sd;

	## Load the map file ##
	print "Loading map file and opening outfiles...\n";

	my $input = new FileHandle ($map_file);
	my @map_lines = ();
	my $lineCounter = 0;
	my %file_handles = ();

	while (<$input>)
	{
		chomp;
		my $line = $_;
		my ($snp, $chrom, $position) = split(/\t/, $line);
		
		$chrom = "X" if($chrom eq "23");
		$chrom = "Y" if($chrom eq "24");
		
		if($chrom ne "CHR" && !$file_handles{$chrom})
		{
			open($file_handles{$chrom}, ">$output_basename.$chrom.infile") or die "Can't open outfile: $!\n";
		}
		
		$lineCounter++;
		$map_lines[$lineCounter] = $line;
	}
	close($input);

	print "$lineCounter lines loaded\n";


	## Parse the CMDS values file ##

	my %position_printed = ();

	## Load the map file ##
	print "Parsing CMDS file...\n";

	$input = new FileHandle ($cmds_file);
	$lineCounter = 0;

	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;

		## Skip header ##

		if($lineCounter > 1)
		{
			if($map_lines[$lineCounter])
			{
				my $log2 = $line;
				my ($snp, $chrom, $position) = split(/\t/, $map_lines[$lineCounter]);

				$chrom = "X" if($chrom eq "23");
				$chrom = "Y" if($chrom eq "24");

				if(!$position_printed{"$chrom\t$position"})
				{
					if($file_handles{$chrom})
					{
						my $outfile = $file_handles{$chrom};
						print $outfile "$chrom\t$position\t$log2\n";					
					}
					else
					{
						warn "NO Filehandles for $chrom\n";
					}
					$position_printed{"$chrom\t$position"} = 1;
				}


			}
			else
			{
				die "No map result for line $lineCounter; make sure $map_file is same length as $cmds_file\n";
			}			
		}


	}
	close($input);


	## Close output file ##
	
	foreach my $chrom (sort keys %file_handles)
	{
		close($file_handles{$chrom}) if($file_handles{$chrom});
		my $chrom_filename = "$output_basename.$chrom.infile";
		my $script_filename = $chrom_filename . ".R";
		my $image_filename = "$output_basename.$chrom.png";
		open(SCRIPT, ">$script_filename") or die "Can't open script $script_filename: $!\n";
	
		print SCRIPT "library(DNAcopy)\n";

		print SCRIPT "regions <- read.table(\"$chrom_filename\")\n";
		print SCRIPT "png(\"$image_filename\", height=600, width=800)\n";

		print SCRIPT qq{
CNA.object <- CNA(regions\$V3, regions\$V1, regions\$V2, data.type="logratio", sampleid=c("Chromosome $chrom"))\n
smoothed.CNA.object <- smooth.CNA(CNA.object)\n
segment.smoothed.CNA.object <- segment(smoothed.CNA.object, undo.splits="sdundo", undo.SD=$num_sd, verbose=1)
p.segment.smoothed.CNA.object <- segments.p(segment.smoothed.CNA.object)
plot(segment.smoothed.CNA.object, type="w", cex=0.5, cex.axis=1.5, cex.lab=1.5)
write.table(p.segment.smoothed.CNA.object, file="$chrom_filename.segments.p_value")
};
		print SCRIPT "dev.off()\n";
		close(SCRIPT);
		
		print "Running $script_filename\n";
		system("R --no-save < $script_filename");# if($chrom eq "X" || $chrom eq "Y");	
	}

	## Parse out segments and build an index HTML file ##

	open(SEGMENTS, ">$output_basename.segments.tsv") or die "Can't open outfile: $!\n";
	print SEGMENTS join("\t", "ID", "sample", "chrom", "loc.start", "loc.end", "num.mark", "seg.mean", "bstat", "pval", "lcl", "ucl") . "\n";

	open(INDEX, ">$output_basename.index.html") or die "Can't open outfile: $!\n";
	print INDEX "<HTML><BODY><TABLE CELLSPACING=0 CELLPADDING=5 BORDER=0 WIDTH=\"100%\">\n";
	print INDEX "<TR>\n";

	## Figure out the relative path to image files ##
	my @tempArray = split(/\//, $output_basename);
	my $numElements = @tempArray;
	my $image_basename = $tempArray[$numElements - 1];
	my $num_printed_in_column = 0;
	
	for(my $chrom = 1; $chrom <= 24; $chrom++)
	{
		my $chrom_name = $chrom;
		$chrom_name = "X" if($chrom == 23);
		$chrom_name = "Y" if($chrom == 24);

		my $chrom_filename = $output_basename . ".$chrom_name.infile";
		my $segments_filename = "$chrom_filename.segments.p_value";
		my $image_filename = $image_basename . "." . $chrom_name . ".png";
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
			$result .= join("\t", @lineContents) . "\n";
		}
	}
	
	close($input);
	
	return($result);
}


###############################################################################
# commify - add appropriate commas to long integers
###############################################################################

sub commify
{
	local($_) = shift;
	1 while s/^(-?\d+)(\d{3})/$1,$2/;
	return $_;
}


1;

