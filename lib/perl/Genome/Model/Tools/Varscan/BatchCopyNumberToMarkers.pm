
package Genome::Model::Tools::Varscan::BatchCopyNumberToMarkers;     # rename this when you give the module file a different name <--

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

class Genome::Model::Tools::Varscan::BatchCopyNumberToMarkers {
    is => 'Command',

    has => [                                # specify the command's single-value properties (parameters) <--- 
        sample_dirs_file => {
            is => 'Text',
            doc => "Tab-delimited file of sample and VarScan-CopyNumber dir",
            is_optional => 0,
        },
        bed_file => {
            is => 'Text',
            doc => "BED file of exon (or tier 1) definitions",
            is_optional => 0,
        },
        output_dir => {
            is => 'Text',
            doc => "Output directory where per-sample files will be created",
            is_optional => 0,
        },
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
	my $sample_dirs_file = $self->sample_dirs_file;
	my $bed_file = $self->bed_file;	
	my $output_dir = $self->output_dir;

	print "Loading BED Targets...\n";
	my %cds_exons = load_cds_exons($bed_file);

	## Parse the sample Dirs file ##

	my $input = new FileHandle ($sample_dirs_file);
	my $lineCounter = 0;

	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		my ($sample, $dir) = split(/\t/, $line);

		my $num_exons = my $num_exons_matched = 0;

		## parse out patient ID ##
		
		my @sampleContents = split(/\-/, $sample);
#		my $patient = join("-", "TCGA", $sampleContents[1], $sampleContents[2]);
		my $patient = $sample;
#		my $patient = join("-", "TCGA", $sampleContents[1], $sampleContents[2], $sampleContents[3], $sampleContents[4], $sampleContents[5]);
		print "$patient\n";
		my $output_file = $output_dir . "/$patient.varScan.cbs.exon.tsv";
		
		if(!(-e $output_file))
		{
			my $cmd = "gmt varscan copy-number-to-markers --sample-dir $dir --bed-file $bed_file --output-file $output_file";
			system("bsub -q $ENV{GENOME_LSF_QUEUE_BUILD_WORKER} -R\"select[mem>4000] rusage[mem=4000]\" $cmd");
			sleep(1);
		}
		else
		{
			print "Skipping $patient...\n";
		}
	}
	
	close($input);

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

