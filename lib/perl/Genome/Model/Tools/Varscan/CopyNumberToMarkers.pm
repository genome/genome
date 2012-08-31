
package Genome::Model::Tools::Varscan::CopyNumberToMarkers;     # rename this when you give the module file a different name <--

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

class Genome::Model::Tools::Varscan::CopyNumberToMarkers {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		sample_dir	=> { is => 'Text', doc => "Sample directory with VarScan output", is_optional => 0 },
		basename	=> { is => 'Text', doc => "Base name for copy number files", is_optional => 0, default => 'varScan.*.copynumber.cbs' },
		bed_file 	=> { is => 'Text', doc => "BED file of exon (or tier 1) definitions", is_optional => 0 },
		output_file 	=> { is => 'Text', doc => "Output file for marker-based copynumber", is_optional => 0 },
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Output VarScan copy number mean at positions (e.g. exons) in a BED file"                 
}

sub help_synopsis {
    return <<EOS
This command outputs VarScan copy number mean at positions (e.g. exons) in a BED file
EXAMPLE:	gmt varscan copy-number-to-markers --bed-file myExons.bed --sample-dirs-file Sample-VarScan-CopyNumber-Dirs.tsv --output-dir exon_segments/ ...
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
	my $sample_dir = $self->sample_dir;
	my $bed_file = $self->bed_file;	
	my $output_file = $self->output_file;
	my $basename = $self->basename;
	
	print "Loading BED Targets...\n";
	my %cds_exons = load_cds_exons($bed_file);

	open(OUTFILE, ">$output_file") or die "Can't open outfile: $!\n";
	
	my $num_exons = my $num_exons_matched = 0;
	
	for(my $chrCounter = 1; $chrCounter <= 23; $chrCounter++)
	{			
		my $chrom = $chrCounter;
		$chrom = "X" if($chrCounter == 23);

		print "$chrom...";
		
		my $exome_file = `ls $sample_dir/$basename.$chrom.infile.segments.p_value`;
		chomp($exome_file);

		if(-e $exome_file)
		{
			## Load the exon targets on this chrom ##
			my $exons = $cds_exons{$chrom};
			
			## Load the sample's CBS segments for this chrom ##
			my $copy_regions = load_regions($exome_file);

			my @exons = split(/\n/, $exons);
			
			foreach my $exon (@exons)
			{
				$num_exons++;
				my $matched_flag = 0;
				
				## Parse out exon coordinates ##
				
				my ($chrom, $chr_start, $chr_stop, $gene) = split(/\t/, $exon);
				my $exon_result = "";
				
				## Get Result for Start of Exon ##

				$exon_result = get_copy_result($chrom, $chr_start, $copy_regions);
				
				if($exon_result)
				{
					$matched_flag = 1;
					print OUTFILE join("\t", $chrom, $chr_start, $exon_result) . "\n";
				}

				## Get Result for Stop of Exon ##

				$exon_result = get_copy_result($chrom, $chr_stop, $copy_regions);
				
				if($exon_result)
				{
					$matched_flag = 1;
					print OUTFILE join("\t", $chrom, $chr_stop, $exon_result) . "\n";
				}
				
				$num_exons_matched++ if($matched_flag);
			}
		}
		else
		{
			warn "File not found: $sample_dir/varScan.*.copynumber.cbs.$chrom.infile.segments.p_value\n";
		}
	}
		

	close(OUTFILE);
	print "\n";
	
	print "$num_exons_matched of $num_exons exons matched\n";			

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

sub get_copy_result
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

