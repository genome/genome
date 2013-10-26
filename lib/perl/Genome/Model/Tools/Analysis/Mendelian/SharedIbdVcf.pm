
package Genome::Model::Tools::Analysis::Mendelian::SharedIbdVcf;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# SearchRuns - Search the database for runs
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@genome.wustl.edu)
#
#	CREATED:	07/12/2013 by D.K.
#	MODIFIED:	07/12/2013 by D.K.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;


use FileHandle;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it


class Genome::Model::Tools::Analysis::Mendelian::SharedIbdVcf {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		vcf_file	=> { is => 'Text', doc => "Input in VCF format", is_optional => 0, is_input => 1},
		output_basename	=> { is => 'Text', doc => "Output basename for files", is_optional => 0, is_input => 1},
		control_samples	=> { is => 'Text', doc => "Comma-separated list of control sample names", is_optional => 1, is_input => 1},
		male_samples	=> { is => 'Text', doc => "Comma-separated list of male sample names", is_optional => 1, is_input => 1},		
		inheritance_model	=> { is => 'Text', doc => "Presumed model of mendelian inheritance [autosomal-dominant]", is_optional => 1, is_input => 1, default => 'autosomal-dominant'},
		min_coverage	=> { is => 'Text', doc => "Minimum coverage to refute a possible variant in an affected", is_optional => 1, is_input => 1, default => 20},
		min_call_rate	=> { is => 'Text', doc => "Minimum callrate for affecteds to include a variant", is_optional => 1, is_input => 1, default => 0.50},
		plot_resolution	=> { is => 'Text', doc => "Bin size for plotting shared IBD segments", is_optional => 1, is_input => 1, default => 100000},
		path_to_beagle	=> { is => 'Text', doc => "Path to the BEAGLE executable required for fastIBD, i.e. /gsc/pkg/bio/beagle/beagle-3.3.12/beagle.jar", is_optional => 0, is_input => 1},
		centromere_file	=> { is => 'Text', doc => "A UCSC 0-based BED file of centromere locations per chromosome", is_optional => 0, is_input => 1, default => '/gscmnt/sata809/info/medseq/dkoboldt/SNPseek2/ucsc/hg19/gapTable.centromere.nochr.txt'},
		hapmap_file_dir	=> { is => 'Text', doc => "Directory containing genetic distances from HapMap, e.g. /gscmnt/sata809/info/medseq/dkoboldt/SNPseek2/ucsc/geneticMap/", is_optional => 0},
		hapmap_file_name => { is => 'Text', doc => "Search string to use for chromosome file with CHROM to replace chromosome name, e.g. genetic_map_GRCh37_chrCHROM.txt", is_optional => 0},
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Attempts to identify regions of identity-by-descent shared among affecteds"                 
}

sub help_synopsis {
    return <<EOS
This tool attempts to identify regions of identity-by-descent (IBD) shared among affecteds
EXAMPLE:	gmt analysis mendelian shared-ibd-vcf --vcf-file myVCF.vcf --output-basename myVCF.shared-ibd
EOS
}

sub help_detail {                           # this is what the user will see with the longer version of help. <---
    return <<EOS 
This tool attempts to identify regions of identity-by-descent (IBD) shared among affecteds
EOS
}


################################################################################################
# Execute - the main program logic
#
################################################################################################

sub execute {                               # replace with real execution logic.
	my $self = shift;
	my $vcf_file = $self->vcf_file;

	my %stats = ();

	my %control_sample = ();
	if($self->control_samples)
	{
		my @samples = split(/\,/, $self->control_samples);
		foreach my $sample (@samples)
		{
			$control_sample{$sample} = 1;
			$stats{'num_samples_control'}++;
		}
	}

	my %male_sample = ();
	if($self->male_samples)
	{
		my @samples = split(/\,/, $self->male_samples);
		foreach my $sample (@samples)
		{
			$male_sample{$sample} = 1;
			$stats{'num_samples_male'}++;
		}
	}


	my %centromeres = load_centromeres($self);

	my $num_aff_pairs = parse_affected_samples($self);

	if($num_aff_pairs)
	{
		print "$num_aff_pairs affected sample pairs parsed from " . $self->vcf_file . "\n";
	}
	else
	{
		die "Unable to parse affected samples from VCF\n";
	}
	
	my $min_num_aff_pairs = $num_aff_pairs;
	for(my $aff_pairs = $num_aff_pairs; $aff_pairs > 0; $aff_pairs--)
	{
		if(($aff_pairs / $num_aff_pairs) > 0.80)
		{
			$min_num_aff_pairs = $aff_pairs			
		}

	}

	my $affected_pairs_file = $self->output_basename . ".pairs.affected";
	my $output_basename = $self->output_basename;
	my @temp = split(/\//, $output_basename);
	my $numContents = @temp;
	my $true_basename = $temp[$numContents - 1];
	
	open(INDEX, ">$output_basename.plot.index.html") or die "Can't open outfile: $!\n";
	print INDEX "<HTML><BODY><TABLE CELLSPACING=0 CELLPADDING=5 BORDER=0 WIDTH=\"100%\">\n";
	print INDEX "<TR>\n";
	my $num_printed_in_column = 0;


#	for(my $chrCounter = 2; $chrCounter <= 2; $chrCounter++)
	for(my $chrCounter = 1; $chrCounter <= 23; $chrCounter++)	
	{
		my $chrom = $chrCounter;
		$chrom = "X" if($chrCounter == 23);
		$chrom = "Y" if($chrCounter == 24);

		my $cen_start = my $cen_stop = -1;
		($cen_start, $cen_stop) = split(/\t/, $centromeres{$chrom}) if($centromeres{$chrom});
		
		warn "CHROMOSOME $chrom\n";
		
		## Find the file ##
		my $query = $self->hapmap_file_name;
		$query =~ s/CHROM/$chrom/;
		my $hapmap_file = $self->hapmap_file_dir . "/" . $query;
		
		
		my $beagle_input_file = $output_basename . ".input." . $chrom;
		my $beagle_markers_file = $output_basename . ".markers." . $chrom;
		my $beagle_output_base = $output_basename . ".output." . $chrom;
		
		if(-e $hapmap_file)
		{
			warn "Converting to BEAGLE format...\n";
			## Convert to beagle format ##
			my $cmd = "gmt vcf convert-to-beagle --vcf-file " . $self->vcf_file . " --output-file $beagle_input_file --markers-file $beagle_markers_file --chromosome $chrom --hapmap-file $hapmap_file";
			system($cmd);

			if(-e $beagle_input_file && -e $beagle_markers_file)
			{
				$cmd = "java -Xmx14000m -jar " . $self->path_to_beagle . " fastibd=true ibdpairs=" . $affected_pairs_file . " unphased=$beagle_input_file out=$beagle_output_base omitprefix=true markers=$beagle_markers_file missing=? ";
				system($cmd);
				
				## Convert FIBD file ##
				
				my $beagle_fibd_file = $beagle_input_file . ".fibd";
				system("gunzip $beagle_fibd_file.gz") if -e "$beagle_fibd_file.gz";
				## REmove unnecessary beagle files ##
				system("rm -rf $beagle_input_file.gprobs.gz");
				system("rm -rf $beagle_input_file.phased.gz");
				system("rm -rf $beagle_input_file.dose.gz");
				system("rm -rf $beagle_input_file.r2");
				system("rm -rf $beagle_input_file.ibd");
				
				if(-e $beagle_fibd_file)
				{
					warn "Converting FIBD file to coordinates...\n";
					convert_fibd_file($beagle_fibd_file, $beagle_markers_file);
					$beagle_fibd_file .= ".coords";
					
					if(-e $beagle_fibd_file)
					{
						## Build the IBD plot ##
						build_chrom_plot($beagle_fibd_file, $output_basename . ".$chrom", $chrom, $cen_start, $cen_stop, $self, $min_num_aff_pairs, $num_aff_pairs);# if($chrom eq "8");
						my $image_filename = "$true_basename.$chrom.plot.png";					
						print INDEX "<TD><A HREF=\"$image_filename\"><IMG SRC=\"$image_filename\" HEIGHT=100 WIDTH=400 BORDER=0></A></TD>\n";
					
						$num_printed_in_column++;
					
						if($num_printed_in_column >= 4)
						{
							print INDEX "</TR><TR>\n";
							$num_printed_in_column = 0;
						}
					}
				}
				else
				{
					warn "Warning: Beagle FastIBD output file $beagle_fibd_file not found!\n";
				}
				
			}
			else
			{
				warn "***ERROR*** Missing converted BEAGLE file(s) $beagle_input_file $beagle_markers_file\n";
			}
		}
		else
		{
			warn "$hapmap_file not found\n";
		}
	}
	
	close(INDEX);

	
	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}


################################################################################################
# LoadAnnotation - load the VEP annotation 
#
################################################################################################

sub parse_affected_samples
{
	my $self = shift(@_);
	## Parse the VCF file ##
	my %aff_samples = ();

	my %control_sample = ();
	if($self->control_samples)
	{
		my @samples = split(/\,/, $self->control_samples);
		foreach my $sample (@samples)
		{
			$control_sample{$sample} = 1;
		}
	}
	
	my $input = new FileHandle ($self->vcf_file);
	my $lineCounter = 0;

	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
			
		my @lineContents = split(/\t/, $line);
		my $numContents = @lineContents;
		
		if(substr($line, 0, 1) eq '#')
		{
			if(substr($line, 0, 6) eq '#CHROM')
			{
				for(my $colCounter = 9; $colCounter < $numContents; $colCounter++)
				{
					my $sample = $lineContents[$colCounter];
					if($control_sample{$sample})
					{
						## Skip ##
					}
					else
					{
						## include ##
						$aff_samples{$sample} = 1;
					}
				}
				
				## Output the sample pairs file ##
				my $num_aff_pairs = 0;
				
				my $outfile_name = $self->output_basename . ".pairs.affected";
				open(OUTFILE, ">$outfile_name") or die "Can't open outfile: $!\n";
				print OUTFILE "sample1\tsample2\n";
				my %included = ();
				foreach my $sample1 (sort keys %aff_samples)
				{
					foreach my $sample2 (sort keys %aff_samples)
					{
						if($sample1 ne $sample2)
						{
							my $key = join("\t", $sample1, $sample2);
							my $key2 = join("\t", $sample2, $sample1);
							
							if(!$included{$key} && !$included{$key2})
							{
								$num_aff_pairs++;
								print OUTFILE "$key\n";
								$included{$key} = 1;
								$included{$key2} = 1;
							}
						}
					}
				}
				
				## Close input and return ##
				close($input);
				return($num_aff_pairs);
				
			}
		}
	}
	
	close($input);
	return(0);
}


################################################################################################
# LoadAnnotation - load the VEP annotation 
#
################################################################################################

sub load_centromeres
{
	my $self = shift(@_);
	my %centromeres = ();
	## Parse the VCF file ##
	
	my $input = new FileHandle ($self->centromere_file);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		my ($chrom, $chr_start, $chr_stop) = split(/\t/, $line);
		$centromeres{$chrom} = join("\t", $chr_start, $chr_stop);
	}
	
	close($input);

	return(%centromeres);
}



#############################################################
# parse_file - parses the file
#
#############################################################

sub convert_fibd_file
{
	my ($fibd_file, $markers_file) = @_;
	
	my @markers = load_markers($markers_file);

	## Create new fibd output_file ##
	
	my $output_file = $fibd_file . ".coords";
	open(OUTFILE, ">$output_file") or die "Can't open output file: $!\n";

	my $input = new FileHandle ($fibd_file);
	my $lineCounter = 0;
	
	my %sample_combo = ();
	my $combo_counter = 1;
	
	my $size_sum = my $score_sum = 0;
	my $num_segments = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		
	
		my ($sample1, $sample2, $index1, $index2, $score) = split(/\t/, $line);
		
		my $marker_start = $markers[$index1];
		my $marker_stop = $markers[$index2];
		$marker_stop = $markers[$index2 - 1] if(!$markers[$index2]);

		my $num_mark = $index2 - $index1;

		my ($chrom) = split(/\:/, $marker_start);		
#		my ($start_name, $start_pos) = split(/\t/, $marker_start);
#		my ($stop_name, $stop_pos) = split(/\t/, $marker_stop);
		my ($start_name) = split(/\t/, $marker_start);
		my ($stop_name) = split(/\t/, $marker_stop);
		my ($start_chrom, $start_pos) = split(/\:/, $start_name);
		my ($stop_chrom, $stop_pos) = split(/\:/, $stop_name);

		$size_sum += $num_mark;
		$score_sum += $score;
		$num_segments++;

		my $combo_id = 0;
		if($sample_combo{$sample1 . "\t" . $sample2})
		{
			$combo_id = $sample_combo{$sample1 . "\t" . $sample2};
		}
		else
		{
			$combo_id = $combo_counter;
			$sample_combo{$sample1 . "\t" . $sample2} = $combo_id;
			$combo_counter++;
		}

		print OUTFILE join("\t", $sample1, $sample2, $index1, $index2, $score, $num_mark, $chrom, $start_pos, $stop_pos, $combo_counter) . "\n";
	}
	
	close($input);
	
	close(OUTFILE);
}



#############################################################
# parse_file - parses the file
#
#############################################################

sub build_chrom_plot
{
	my ($infile, $outbase, $chrom_name, $cen_start, $cen_stop, $self, $min_num_aff_pairs, $num_aff_pairs) = @_;
	
	my %segments = parse_ibd_segments($infile, $self);

	my $resolution = $self->plot_resolution;

	my @regions = ();
	my $numRegions = 0;
	my $num_pairs = 0;
	my %ibd_positions = ();
		
	foreach my $pair (sort keys %segments)
	{
		$num_pairs++;
		my @segments = split(/\n/, $segments{$pair});
		my $numSegments = @segments;
	
		my %position_counted = ();
	
		print "$pair\t$numSegments\n";
	
		foreach my $segment (@segments)
		{
			my ($chrom, $start_pos, $stop_pos, $num_mark, $score) = split(/\t/, $segment);
	
			$start_pos = sprintf("%d", $start_pos / $resolution);
			$stop_pos = sprintf("%d", $stop_pos / $resolution);
	
			for(my $position = $start_pos; $position <= $stop_pos; $position++)
			{
				my $key = join("\t", $chrom, $position);
				
				if(!$position_counted{$key})
				{
					if($ibd_positions{$key})
					{
						$ibd_positions{$key}++;
					}
					else
					{
						$ibd_positions{$key} = 1;
					}
					$position_counted{$key} = 1;
				}
				
	
			}
		}
	
	
	
	}
	
	warn "$num_pairs unique affected pairs\n";

	my $shared_ibd_positions_outfile = "$outbase.shared.positions";
	my $shared_ibd_regions_outfile = "$outbase.shared.regions";
	open(OUTFILE, ">$shared_ibd_positions_outfile") or die "Can't open outfile: $!\n";
	print OUTFILE "chrom\tposition\tibd_samples\n";

	open(COMBINEDREGIONS, ">$shared_ibd_regions_outfile") or die "Can't open outfile: $!\n";
	print COMBINEDREGIONS "chrom\tchr_start\tchr_stop\tsize_bp\tibd_samples\n";
	
#	print "See $outbase.ibd\n";
#	print "See $outbase.ibdregions\n";
	
	my $region_chrom = my $region_start = my $region_stop = my $region_ibd = 0;
	my $region_positions = 0;
	
	foreach my $key (sort byChrPos keys %ibd_positions)
	{
		my ($chrom, $position) = split(/\t/, $key);
		$position = $position * $resolution;
		my $ibd_count = $ibd_positions{$key};
	
		if(!$region_ibd)
		{
			## Start new region ##
			$region_chrom = $chrom;
			$region_start = $region_stop = $position;
			$region_ibd = $ibd_count;
			$region_positions = 1;
		}
		elsif($chrom ne $region_chrom || $ibd_count != $region_ibd)
		{
			## End current region and start a new one ##
	
			print COMBINEDREGIONS join("\t", $region_chrom, $region_start, $region_stop, $region_positions, $region_ibd) . "\n";
			$region_chrom = $chrom;
			$region_start = $region_stop = $position;
			$region_ibd = $ibd_count;
			$region_positions = 1;
		}
		else
		{
			## Add to current region	
			$region_stop = $position;
			$region_positions++;
		}
	
	
		print OUTFILE join("\t", $chrom, $position, $ibd_count) . "\n";
	}

	close(OUTFILE);
	close(COMBINEDREGIONS);

	
	open(OUTSCRIPT, ">$outbase.R") or die "Can't open outfile: $!\n";
	
#	print OUTSCRIPT qq{ibd <- read.table("$outbase.ibd", header=T)\n
#	png("$outbase.plot.png", height=600, width=800)\n
#	plot(ibd\$position, ibd\$ibd_samples, pch=19, cex=0.25, col="blue", type="h", ylim=c(0,$num_pairs), xlab="Position", ylab="Samples IBD", main="$chrom_name")\n
#	};

	print OUTSCRIPT "library(ggplot2)\n";

	## REmoved this x-axis label: + xlab("Position on chr$chrom_name")  ##

	if($cen_start && $cen_stop)
	{
		print OUTSCRIPT qq{ibd <- read.table("$shared_ibd_positions_outfile", header=T)\n
		png("$outbase.plot.png", height=200, width=800)\n
		ggplot(ibd, aes(x = position, y = ibd_samples)) + geom_area(fill="blue") + xlab("") + ylab("Shared IBD") + opts(title="chr$chrom_name") + scale_y_continuous(limits = c(0,$num_aff_pairs)) + geom_vline(xintercept = $cen_start, linetype = \"longdash\") + geom_vline(xintercept = $cen_stop, linetype = \"longdash\") + geom_hline(yintercept = $min_num_aff_pairs, color = \"red\")\n
		};
		
	}
	else
	{
		print OUTSCRIPT qq{ibd <- read.table("$shared_ibd_positions_outfile", header=T)\n
		png("$outbase.plot.png", height=200, width=800)\n
		ggplot(ibd, aes(x = position, y = ibd_samples)) + geom_area(fill="blue") + xlab("") + ylab("Shared IBD") + opts(title="chr$chrom_name") + scale_y_continuous(limits = c(0,$num_aff_pairs)) + geom_hline(yintercept = $min_num_aff_pairs, color = \"red\")\n
		};
		
	}


	
	print OUTSCRIPT "dev.off()\n";

	## Build 150x150 plot ##

	if($cen_start && $cen_stop)
	{
		print OUTSCRIPT qq{
		png("$outbase.smallplot.png", height=150, width=150)\n
		ggplot(ibd, aes(x = position, y = ibd_samples)) + geom_area(fill="blue") + xlab("") + ylab("") + scale_y_continuous(limits = c(0,$num_aff_pairs)) + opts(title="chr$chrom_name", plot.margin = unit(c(0.01,0.01,0.01,0.01), \"cm\"), axis.ticks=theme_blank(), axis.text.x=theme_blank(),axis.title.x=theme_blank(), axis.text.y=theme_blank(),axis.title.y=theme_blank()) + geom_vline(xintercept = $cen_start, linetype = \"longdash\") + geom_vline(xintercept = $cen_stop, linetype = \"longdash\") + geom_hline(yintercept = $min_num_aff_pairs, color = \"red\")\n
		};		
	}
	else
	{
		print OUTSCRIPT qq{
		png("$outbase.smallplot.png", height=150, width=150)\n
		ggplot(ibd, aes(x = position, y = ibd_samples)) + geom_area(fill="blue") + xlab("") + ylab("") + scale_y_continuous(limits = c(0,$num_aff_pairs)) + opts(title="chr$chrom_name", plot.margin = unit(c(0.01,0.01,0.01,0.01), \"cm\"), axis.ticks=theme_blank(), axis.text.x=theme_blank(),axis.title.x=theme_blank(), axis.text.y=theme_blank(),axis.title.y=theme_blank()) + geom_hline(yintercept = $min_num_aff_pairs, color = \"red\")\n
		};		
	}



	
	print OUTSCRIPT "dev.off()\n";



	close(OUTSCRIPT);
	
	system("R --no-save <$outbase.R");



}



#############################################################
# parse_file - parses the file
#
#############################################################

sub parse_ibd_segments
{
	my %segments = ();
	my $FileName = shift(@_);
	my $self = shift(@_);

	my %control_sample = ();
	if($self->control_samples)
	{
		my @samples = split(/\,/, $self->control_samples);
		foreach my $sample (@samples)
		{
			$control_sample{$sample} = 1;
		}
	}


	my $input = new FileHandle ($FileName);
	my $lineCounter = 0;
	

	my $combo_counter = 1;
	
	my $size_sum = my $score_sum = 0;
	my $num_segments = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		
	
		my ($sample1, $sample2, $index1, $index2, $score, $num_mark, $chrom, $start_pos, $stop_pos) = split(/\t/, $line);
		
		if($control_sample{$sample1} || $control_sample{$sample2})
		{
			## Ignore segments shared with control samples ##	
		}
		elsif($score >= 10e-10)
		{
			## These are less likely to be in IBD
		}
		else
		{
			my $key = join("\t", $sample1, $sample2);
			$segments{$key} .= "\n" if($segments{$key});
			$segments{$key} .= join("\t", $chrom, $start_pos, $stop_pos, $num_mark, $score);			
		}

	}
	
	close($input);
	
	return(%segments);
}



#############################################################
# parse_file - parses the file
#
#############################################################

sub load_markers
{
	my $FileName = shift(@_);

	my @markers = ();
	my $input = new FileHandle ($FileName);
	my $lineCounter = 0;
	
	
	while (<$input>)
	{
		chomp;
		my $line = $_;

		$markers[$lineCounter] = $line;
		$lineCounter++;		
	
	}
	
	close($input);	

	return(@markers);
}


sub byChrPos
{
	my ($chrom_a, $pos_a) = split(/\t/, $a);
	if($chrom_a eq "X")
	{
		$chrom_a = 23;
	}
	elsif($chrom_a eq "Y")
	{
		$chrom_a = 24;		
	}
	
	my ($chrom_b, $pos_b) = split(/\t/, $b);

	if($chrom_b eq "X")
	{
		$chrom_b = 23;
	}
	elsif($chrom_b eq "Y")
	{
		$chrom_b = 24;		
	}
	
	$chrom_a <=> $chrom_b
	or
	$pos_a <=> $pos_b;

}


sub numericallyDesc
{
	$b = 0 if($b eq '.');
	$a = 0 if($a eq '.');
	$b <=> $a;
}

sub commify
{
	local($_) = shift;
	1 while s/^(-?\d+)(\d{3})/$1,$2/;
	return $_;
}


1;


