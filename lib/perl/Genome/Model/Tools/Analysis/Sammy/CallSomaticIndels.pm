
package Genome::Model::Tools::Analysis::Sammy::CallSomaticIndels;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# CallSomatic - Call somatic variants from normal/tumor BAM files
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@watson.wustl.edu)
#
#	CREATED:	07/28/2009 by D.K.
#	MODIFIED:	07/28/2009 by D.K.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;

use FileHandle;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

class Genome::Model::Tools::Analysis::Sammy::CallSomaticIndels {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		output_dir	=> { is => 'Text', doc => "Output directory for somatic calls", is_optional => 0 },
		sample_name	=> { is => 'Text', doc => "Sample name for file naming purposes", is_optional => 0 },
		regions_file	=> { is => 'Text', doc => "Tab-delimited file of target regions", is_optional => 0 },
		normal_bam	=> { is => 'Text', doc => "BAM file for normal sample", is_optional => 1 },
		tumor_bam	=> { is => 'Text', doc => "BAM file for tumor sample" , is_optional => 1},
		normal_pileup	=> { is => 'Text', doc => "Pileup file for normal sample", is_optional => 1 },
		tumor_pileup	=> { is => 'Text', doc => "Pileup file for tumor sample" , is_optional => 1},
		reference	=> { is => 'Text', doc => "Reference file for alignments" , is_optional => 1},
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Calls somatic variants from normal and tumor BAM files"                 
}

sub help_synopsis {
    return <<EOS
This command calls somatic variants from Normal and Tumor alignments files using Dan Koboldt's Sammy package
EXAMPLE:	gmt analysis sammy
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
#	my $sample_name = $self->sample_name;
	my $reference_file = $self->reference;
	$reference_file = Genome::Config::reference_sequence_directory() . '/NCBI-human-build36/all_sequences.fa' if(!$self->reference);

	## Create directory ##
	mkdir($self->output_dir) if(!(-d $self->output_dir));
	my $cmd;

	## Ready pileup files ##
	
	my $normal_pileup = $self->normal_pileup;
	my $tumor_pileup = $self->tumor_pileup;
	
	$normal_pileup = $self->output_dir . "/" . $self->sample_name . ".normal.pileup" if(!$normal_pileup);
	$tumor_pileup = $self->output_dir . "/" . $self->sample_name . ".tumor.pileup" if(!$tumor_pileup);
	
	my $normal_indel = $self->output_dir . "/" . $self->sample_name . ".normal.indel";
	my $tumor_indel = $self->output_dir . "/" . $self->sample_name . ".tumor.indel";

	$normal_indel = "$normal_pileup.indel" if(-e "$normal_pileup.indel");
	$tumor_indel = "$tumor_pileup.indel" if(-e "$tumor_pileup.indel");

	## Convert BAM to ROI Pileup ##
	if(!(-e $normal_pileup))
	{
		if($self->normal_bam)
		{
			die "Normal BAM file $self->normal_bam not found!\n" if(!(-e $self->normal_bam));
	
			$cmd = "samtools pileup -f " . $reference_file . " " . $self->normal_bam;
			$cmd .= " | " . call_sammy() . "limit-indels --regions-file " . $self->regions_file . " --output-file $normal_pileup";
			
			print "Generating Normal Pileup file...\n";
			system($cmd);
		}
		else
		{
			die "Error: Normal BAM file " . $self->normal_bam . " not found!\n";
		}		
	}

	
	## Call indels on ROI pileup ##
	if(-e $normal_pileup && !(-e $normal_indel))
	{
		$cmd = call_sammy() . "pileup2indel " . $normal_pileup . " --min-coverage 10 --min-reads2 2 --min-var-freq 0.25 --p-value 1.0E-06 >$normal_indel";
		print "Calling indels in Normal...\n";
		system($cmd);
	}


	## Convert BAM to ROI Pileup ##

	if(!(-e $tumor_pileup))
	{
		if($self->tumor_bam)
		{
			die "Normal BAM file $self->tumor_bam not found!\n" if(!(-e $self->tumor_bam));

			## If no pileup file exists, extract it from the BAM ##
			$cmd = "samtools pileup -f " . $reference_file . " " . $self->tumor_bam;
			$cmd .= " | " . call_sammy() . "limit-indels --regions-file " . $self->regions_file . " --output-file $tumor_pileup";
			
			print "Generating Tumor Pileup file...\n";
			system($cmd);
		}
		else
		{
			die "Error: Tumor BAM file " . $self->tumor_bam . " not found!\n";	
		}		
	}

	
	## Call indels on ROI pileup ##
	if(-e $tumor_pileup && !(-e $tumor_indel))
	{
		$cmd = call_sammy() . "pileup2indel " . $tumor_pileup . " --min-coverage 10 --min-reads2 2 --min-var-freq 0.25 --p-value 1.0E-06 >$tumor_indel";
		print "Calling indels in Tumor...\n";
		system($cmd);
	}


	## Compare indels between normal and tumor ##
	
	my $compared_indels = $self->output_dir . "/" . $self->sample_name . ".indels.compared";

	if(-e $normal_indel && -e $tumor_indel)
	{	
		## Compare the indels ##
		
		$cmd = call_sammy() . "compare " . $normal_indel . " " . $tumor_indel . " " . $compared_indels;

		print "Comparing indels between Normal and Tumor...\n";
		system($cmd);
	}
	else
	{
		die "Error: Normal or Tumor indels file missing...\n";
	}


	## Make the somatic call ##

	my $compared_somatic = $self->output_dir . "/" . $self->sample_name . ".indels.compared.status";
	
	if(-e $compared_indels && -e $normal_pileup && -e $tumor_pileup)
	{
		## Run the somatic ##

		$cmd = call_sammy() . "somatic-indel " . $normal_pileup . " " . $tumor_pileup . " " . $compared_indels . " " . $compared_somatic;		

		print "Calling variants as Germline or Somatic...\n";
		system($cmd);
	}

	
	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}

sub call_sammy
{
	my $classpath = "/gscuser/dkoboldt/Software/Sammy3";
	my $cmd = "java -Xms2000m -Xmx2000m -classpath $classpath Sammy ";
	return($cmd);
}



sub commify
{
	local($_) = shift;
	1 while s/^(-?\d+)(\d{3})/$1,$2/;
	return $_;
}

1;

