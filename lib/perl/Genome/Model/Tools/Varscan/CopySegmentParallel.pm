
package Genome::Model::Tools::Varscan::CopySegmentParallel;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# Varscan::Somatic	Runs Varscan somatic pipeline on Normal/Tumor BAM files
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@genome.wustl.edu)
#
#	CREATED:	12/09/2009 by D.K.
#	MODIFIED:	12/29/2009 by D.K.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;

use FileHandle;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

class Genome::Model::Tools::Varscan::CopySegmentParallel {
    is => 'Genome::Model::Tools::Varscan',

    has => [                                # specify the command's single-value properties (parameters) <--- 
        basename => {
            is => 'Text',
            doc => "Path and basename to varScan copynumber files e.g. dir/varScan.output",
            is_optional => 0,
            is_input => 1,
            is_output => 0,
        },
        output => {
            is => 'Text',
            doc => "Output file for copy number results",
            is_optional => 0,
            is_input => 1,
            is_output => 1,
        },
        chromosome => {
            is => 'Text',
            doc => "Optional chromosome to do one only",
            is_optional => 1,
            is_input => 1,
        },
        reference => {
            is => 'Text',
            doc => "Reference FASTA file for BAMs",
            is_optional => 0,
            example_values => [(Genome::Config::reference_sequence_directory() . '/NCBI-human-build36/all_sequences.fa')]
        },
        skip_if_output_present => {
            is => 'Text',
            doc => "If set to 1, skip execution if output files exist",
            is_optional => 1,
            is_input => 1,
        },
        min_depth => {
            is => 'Text',
            doc => "Minimum depth in both samples to include for segmentation",
            is_optional => 0,
            default => 20,
        },
        undo_sd => {
            is => 'Text',
            doc => "Remove change-points of less than this number of standard deviations",
            is_optional => 0,
            default => 4,
        },
    ],
    has_param => [
        lsf_resource => { default_value => 'select[model!=Opteron250 && type==LINUX64 && tmp>1000] rusage[mem=4000]'},
    ],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Segments copy number files from VarScan copy-number-parallel"                 
}

sub help_synopsis {
    return <<EOS
This command segments copy number files from VarScan copy-number-parallel and builds a summary HTML file
EXAMPLE:	gmt varscan copy-segment-parallel --basename varscan_out/varScan.output.copynumber --output varscan_out/varScan.output.copynumber.cbs ...
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
	my $basename = $self->basename;
	my $output = $self->output;
	my $reference = $self->reference;
	my $index_file = "$reference.fai";
	
	if(!(-e $index_file))
	{
		die "Index file for reference ($index_file) not found!\n";
	}

	## Check skip if output present ##
	
	if($self->skip_if_output_present)
	{

	}

	if(-e $index_file)
	{	
		open(HTML, ">$output.index.html") or die "Can't open outfile: $!\n";
		print HTML "<HTML><BODY><TABLE CELLSPACING=0 CELLPADDING=5 BORDER=0 WIDTH=\"100%\">\n";
		print HTML "<TR>\n";
		my $num_printed_in_column = 0;

		## Get the image basename ##
		my @tempArray = split(/\//, $output);
		my $numElements = @tempArray;
		my $image_basename = $tempArray[$numElements - 1];


		my $input = new FileHandle ($index_file);
		my $lineCounter = 0;
	
		while (<$input>)
		{
			chomp;
			my $line = $_;
			$lineCounter++;
			
			my ($chrom) = split(/\t/, $line);

			if($chrom =~ 'NT_')
			{
#				print "Skipping $chrom\n";								
			}
			else
			{
				if(!$self->chromosome  || $chrom eq $self->chromosome)
				{
					if(substr($chrom, 0, 2) ne "GL")
					{
						print "$chrom\t";
						my $input_file = $basename . ".$chrom.copynumber";
						
						if(-e $input_file)
						{
							my $output_basename = $output . ".$chrom";

							parse_regions($input_file, $output_basename, $chrom, $self);
							my $image_filename = $image_basename . ".$chrom.jpg";				
							print HTML "<TD><A HREF=\"$image_filename\"><IMG SRC=\"$image_filename\" HEIGHT=240 WIDTH=320 BORDER=0></A></TD>\n";					
							$num_printed_in_column++;

							if($num_printed_in_column >= 4)
							{
								print HTML "</TR><TR>\n";
								$num_printed_in_column = 0;
							}

						}
						
					}
				}


			}

		}
	
		close($input);
		
		close(HTML);

	}
	else
	{
		die "Error: One of your BAM files doesn't exist!\n";
	}

	
	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}




################################################################################################
# Execute - the main program logic
#
################################################################################################

sub parse_regions
{                               # replace with real execution logic.
	my $segments_file = shift(@_);
	my $output_basename = shift(@_);
	my $chrom = shift(@_);
	my $self = shift(@_);
	
	my $result = "";

	my $chrom_filename = $output_basename . ".infile";
	open(OUTFILE, ">$chrom_filename") or die "Can't open outfile: $!\n";

	my $input = new FileHandle ($segments_file);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		if($lineCounter > 1)
		{
			my ($chrom, $chr_start, $chr_stop, $num_positions, $normal_depth, $tumor_depth, $log2_ratio) = split(/\t/, $line);
			if($normal_depth >= $self->min_depth && $tumor_depth >= $self->min_depth)
			{
				my $midpoint = sprintf("%d", ($chr_start + $chr_stop) / 2);
				print OUTFILE join("\t", $chrom, $midpoint, $log2_ratio) . "\n";				
			}

		}
	}
	
	close($input);

	close(OUTFILE);
	

	if($lineCounter >= 100)
	{
		my $script_filename = $output_basename . ".R";
		my $undo_sd = $self->undo_sd;

		my $image_filename = $output_basename . ".jpg";				
		
		## Begin R Script ##
	
		open(SCRIPT, ">$script_filename") or die "Can't open script $script_filename: $!\n";	
		print SCRIPT "library(DNAcopy)\n";
		print SCRIPT "regions <- read.table(\"$chrom_filename\")\n";
		print SCRIPT "png(\"$image_filename\", height=600, width=800)\n";
		print SCRIPT qq{
CNA.object <- CNA(regions\$V3, regions\$V1, regions\$V2, data.type="logratio", sampleid=c("Chromosome $chrom"))\n
smoothed.CNA.object <- smooth.CNA(CNA.object)\n
};

		print SCRIPT qq{
segment.smoothed.CNA.object <- segment(smoothed.CNA.object, undo.splits="sdundo", undo.SD=$undo_sd, alpha=0.001, verbose=1)
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
if(length(p.segment.smoothed.CNA.object\$pval[p.segment.smoothed.CNA.object\$num.mark>=40]) < 50)
{
x <- c($this_undo_sd, length(p.segment.smoothed.CNA.object\$pval[p.segment.smoothed.CNA.object\$num.mark>=20]))
write.table(x, file="$chrom_filename.segments.sd", append=TRUE)
segment.smoothed.CNA.object <- segment(smoothed.CNA.object, undo.splits="sdundo", undo.SD=$this_undo_sd, alpha=0.001, verbose=1)
p.segment.smoothed.CNA.object <- segments.p(segment.smoothed.CNA.object)	
}
|;	
}



## Original using their plotting mechanism ##

#		print SCRIPT qq{
#plot(segment.smoothed.CNA.object, type="w", cex=0.5, cex.axis=1.5, cex.lab=1.5, ylim=c(-4,4))
#write.table(p.segment.smoothed.CNA.object, file="$chrom_filename.segments.p_value")
#};

		print SCRIPT qq{
detach(package:DNAcopy)
par(mar=c(4,4,2,2))
plot(segment.smoothed.CNA.object\$data\$maploc, segment.smoothed.CNA.object\$data\$Chromosome.$chrom, pch=19, cex=0.25, cex.axis=1.25, cex.lab=1.5, col="cornflowerblue", ylim=c(-5,5), main="Chromosome $chrom", xlab="Position", ylab="Copy Number Change (log2)")
segments(segment.smoothed.CNA.object\$output\$loc.start, segment.smoothed.CNA.object\$output\$seg.mean, segment.smoothed.CNA.object\$output\$loc.end, segment.smoothed.CNA.object\$output\$seg.mean, col="red", lwd=2)
write.table(p.segment.smoothed.CNA.object, file="$chrom_filename.segments.p_value")
};

		print SCRIPT "dev.off()\n";
		close(SCRIPT);
		print "Running $script_filename\n";
		if($chrom eq "20" || $chrom eq "21" || $chrom eq "22" || $chrom eq "Y" || $chrom eq "MT")
		{
			system("bsub -q $ENV{GENOME_LSF_QUEUE_SHORT} -R\"select[mem>2000] rusage[mem=2000]\" \"R --no-save < $script_filename\"");			
		}
		else
		{
			system("bsub -q $ENV{GENOME_LSF_QUEUE_BUILD_WORKER} -R\"select[mem>2000] rusage[mem=2000]\" \"R --no-save < $script_filename\"");						
		}

		
	}
}


1;

