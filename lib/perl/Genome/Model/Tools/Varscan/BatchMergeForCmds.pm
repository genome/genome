
package Genome::Model::Tools::Varscan::BatchMergeForCmds;     # rename this when you give the module file a different name <--

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

class Genome::Model::Tools::Varscan::BatchMergeForCmds {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		sample_dirs_file		=> { is => 'Text', doc => "Tab-delimited file of sample and VarScan-CopyNumber dir", is_optional => 0 },
		bed_file 	=> { is => 'Text', doc => "BED file of exon (or tier 1) definitions", is_optional => 0 },
		markers_dir 	=> { is => 'Text', doc => "Directory where CopyNumberToMarkers results are stored", is_optional => 0 },
		output_dir 	=> { is => 'Text', doc => "Output directory where per-sample files will be created", is_optional => 0 },
		chromosome 	=> { is => 'Text', doc => "Name of chromosome to include", is_optional => 0 },
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Merges sample copy-number-markers files into per-chromosome files for CMDS"                 
}

sub help_synopsis {
    return <<EOS
This command merges sample copy-number-markers files into per-chromosome files for CMDS
EXAMPLE:	gmt varscan batch-merge-for-cmds --bed-file myExons.bed --sample-dirs-file Sample-VarScan-CopyNumber-Dirs.tsv --markers-dir exon_segments/ --output-dir cmds_input/ ...
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
	my $markers_dir = $self->markers_dir;
	my $output_dir = $self->output_dir;

	print "Loading BED Targets...\n";
	my %cds_exons = load_cds_exons($bed_file);
	my %patient_seg_means = ();
	my %patients = ();

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
		$patients{$patient} = 1;
		print "$patient\n";
		my $markers_file = $markers_dir . "/$patient.varScan.cbs.exon.tsv";
		
		if(-e $markers_file)
		{
			my %seg_means = load_regions($markers_file, $self);

			## Save each seg mean for this patient ##
			foreach my $key (keys %seg_means)
			{
				my $new_key = $patient . "\t" . $key;
				$patient_seg_means{$new_key} = $seg_means{$key};
			}
		}
		else
		{
			print "Skipping $patient because no markers file $markers_file...\n";
		}
	}
	
	close($input);


	## Parse the CDS Exons or Targets File to Build the Output File ##
	
	$input = new FileHandle ($bed_file);
	$lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		my ($chrom, $chr_start, $chr_stop, $gene) = split(/\t/, $line);

		## Only output regions on this chromosome ##
		if($self->chromosome eq $chrom)
		{
			## Save results for chr_start and/or chr_stop ##
			
			my $key = join("\t", $chrom, $chr_start);
			
			my $chr_start_result = $key;
			my $chr_stop_result = $key;
			my $samples_called = my $samples_total = 0;
			
			foreach my $patient (sort keys %patients)
			{
				my $patient_key = join("\t", $patient, $key);
				$samples_total++;
				
				if($patient_seg_means{$patient_key})
				{
					$chr_start_result .= "\t$patient_seg_means{$patient_key}";
					$samples_called++;
				}
				else
				{
					$chr_start_result .= "\tNA";
				}
			}

			my $call_rate = $samples_called / $samples_total;
#			if($call_rate >= 0.50)
#			{
				print OUTFILE $chr_start_result . "\n";				
#			}
			$chr_start_result = "";


			## Chr Stop ##
			
			$key = join("\t", $chrom, $chr_stop);
					
			$samples_called = $samples_total = 0;
			
			foreach my $patient (sort keys %patients)
			{
				my $patient_key = join("\t", $patient, $key);
				$samples_total++;
				
				if($patient_seg_means{$patient_key})
				{
					$chr_stop_result .= "\t$patient_seg_means{$patient_key}";
					$samples_called++;
				}
				else
				{
					$chr_stop_result .= "\tNA";
				}
			}
			

			$call_rate = $samples_called / $samples_total;

#			if($call_rate >= 0.50)
#			{
				print OUTFILE $chr_stop_result  . "\n";				
#			}
			


		}



#		return(0) if($lineCounter > 50);
	}
	
	close($input);	

}





#############################################################
# parse_file - parses the file
#
#############################################################

sub load_regions
{
	my $FileName = shift(@_);
	my $self = shift(@_);
	my $NORMALIZE_MEANS = my $AMP_ONLY = my $DEL_ONLY = 0;
	
	my %seg_means = ();
	
	my $input = new FileHandle ($FileName);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		
		
		my ($chrom, $position, $seg_mean) = split(/\t/, $line);
		
		if($self->chromosome eq $chrom)
		{
			my $key = join("\t", $chrom, $position);

			if($NORMALIZE_MEANS)
			{
				if($seg_mean <= -0.25)
				{
					## Deletion, set to -1 ##
					$seg_mean = -1;
				}
				elsif($seg_mean > -0.25 && $seg_mean < 0.25)
				{
					$seg_mean = 0;
				}
				elsif($seg_mean >= 0.25)
				{
					## Gain, so set to 1##
					$seg_mean = 1 if($seg_mean < 1);
				}				
			}
			elsif($AMP_ONLY)
			{
				if($seg_mean >= 0.25)
				{
					## Gain, so leave it ##
				}
				else
				{
					$seg_mean = 0;
				}
			}
			elsif($DEL_ONLY)
			{
				if($seg_mean <= -0.25)
				{
					## Gain, so leave it ##
				}
				else
				{
					$seg_mean = 0;
				}
			}

			$seg_means{$key} = $seg_mean;			
		}

	}
	
	close($input);

	return(%seg_means);
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






1;

