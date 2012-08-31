
package Genome::Model::Tools::Capture::SomaticToRefAlign;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# SomaticToRefAlign - Run exome-based copy number analysis in parallel on somatic-variation group models
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@genome.wustl.edu)
#
#	CREATED:	02/09/2012 by D.K.
#	MODIFIED:	03/15/2012 by D.K.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;

use FileHandle;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

## Declare global statistics hash ##

my %stats = ();


class Genome::Model::Tools::Capture::SomaticToRefAlign {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		group_id		=> { is => 'Text', doc => "ID of somatic-variation model group" , is_optional => 0},
		output_bam_files	=> { is => 'Text', doc => "Output sample-bam files (one per sample) to this file" , is_optional => 1},
		output_bam_list		=> { is => 'Text', doc => "Output tumor, normal BAM, tumor BAM, normal list for MuSiC to this file" , is_optional => 1},
		output_somatic_dirs	=> { is => 'Text', doc => "Output last succeeded build directory per model to this file" , is_optional => 1},
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Fetches refalign models for a somatic-variation group"                 
}

sub help_synopsis {
    return <<EOS
This command fetches refalign models for a somatic-variation group
EXAMPLE:	gmt capture somatic-to-ref-align --group-id 3328 
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

	my $group_id = $self->group_id;

	## Get the models in each model group ##

	my $model_group = Genome::ModelGroup->get($group_id);
	my @models = $model_group->models;
	my $debug_counter = 0;

	if($self->output_bam_files)
	{
		open(OUTBAMFILES, ">" . $self->output_bam_files) or die "Can't open outfile: $!\n";
	}

	if($self->output_bam_list)
	{
		open(OUTBAMLIST, ">" . $self->output_bam_list) or die "Can't open outfile: $!\n";
	}	

	if($self->output_somatic_dirs)
	{
		open(OUTDIRS, ">" . $self->output_somatic_dirs) or die "Can't open outfile: $!\n";
	}

	foreach my $model (@models)
	{
		$stats{'models_in_group'}++;
		
		my $model_id = $model->genome_model_id;
		my $model_name = $model->name;
		my $subject_name = $model->subject_name;
		$subject_name = "Model" . $model_id if(!$subject_name);

		## Get Builds ##
		my @builds = $model->builds;
		my $last_build_dir = "";
		foreach my $build ( @builds )
		{
			my $build_id = $build->id;
			my $build_status = $build->status;
			if($build_status eq "Succeeded")
			{
				$last_build_dir = $build->data_directory;
			}
		}
		

		## Get normal and tumor model ##
		
		my $normal_model = $model->normal_model;
		my $tumor_model = $model->tumor_model;

		## Get Model IDs ##
		
		my $normal_model_id = $normal_model->id;
		my $tumor_model_id = $tumor_model->id;

		## Get TCGA-Compliant Subject Names ##

		my $normal_subject_name = $normal_model->subject_name;
		my $tumor_subject_name = $tumor_model->subject_name;

		my $tumor_sample = $tumor_model->subject_name;
		my $normal_sample = $normal_model->subject_name;
		my $tumor_model_dir = $tumor_model->last_succeeded_build_directory;
		my $normal_model_dir = $normal_model->last_succeeded_build_directory;
		my $tumor_bam = `ls $tumor_model_dir/alignments/*.bam | head -1`; chomp($tumor_bam);
		my $normal_bam = `ls $normal_model_dir/alignments/*.bam | head -1`; chomp($normal_bam);

		print join("\t", $model_id, $model_name, $last_build_dir) . "\n";

		## Output somatic dirs ##
		
		if($self->output_somatic_dirs && $last_build_dir)
		{
			print OUTDIRS join("\t", $tumor_sample, $normal_sample, $last_build_dir) . "\n";
		}

		## Output bam files (one sample per line)

		if($self->output_bam_files)
		{
			print OUTBAMFILES join("\t", $tumor_sample, $tumor_bam) . "\n";
			print OUTBAMFILES join("\t", $normal_sample, $normal_bam) . "\n";
		}

		## Output bam list for music (tumor name, one pair per line ##
		
		if($self->output_bam_list)
		{
			print OUTBAMLIST join("\t", $tumor_sample, $normal_bam, $tumor_bam, $normal_sample) . "\n";
		}
	}
	
	close(OUTDIRS) if($self->output_somatic_dirs);
	close(OUTBAMFILES) if($self->output_bam_files);
	close(OUTBAMLIST) if($self->output_bam_list);
	
	foreach my $key (sort keys %stats)
	{
		print $stats{$key} . " " . $key . "\n";
	}
	
	return 1;
}




###################################################
# get_flagstat - 
#
###################################################

sub get_flagstat
{
	my $bam_file = shift(@_);
	my $flagstat = "";
	
	if(-e "$bam_file.flagstat")
	{
		$flagstat = `cat $bam_file.flagstat`;	
	}
	else
	{
		$flagstat = `samtools flagstat $bam_file`;
	}


	## IF we got it, parse it ##

	if($flagstat)
	{
		my %cov_stats = ();	
		my @lines = split(/\n/, $flagstat);
		
		foreach my $line (@lines)
		{
			## Erase the plus zero ##
			
			$line =~ s/\s\+\s0//;
			
			(my $num_reads) = split(/\s+/, $line);
			my $category = $line;
			$category =~ s/$num_reads\s//;
			
			## Remove stuff with parentheses ##
			my $split_char = " \\(";
			($category) = split(/$split_char/, $category);
			
			$cov_stats{$category} = $num_reads if($category);			
		}
		
		return(%cov_stats);
	}
	
	return();
}






#############################################################
# ParseFile - takes input file and parses it
#
#############################################################

sub avg_read_len
{
	my $FileName = shift(@_);
	return(100);
	my $len_sum = my $len_num = 0;
	my $read_seqs = `samtools view $FileName 2>/dev/null | head -10000 | cut -f 10`;

	my @lines = split(/\n/, $read_seqs);
	foreach my $line (@lines)
	{
		my $read_len = length($line);
		$len_sum += $read_len;
		$len_num++;
	}

	my $avg_readlen = $len_sum / $len_num if($len_num);
	$avg_readlen = sprintf("%.2f", $avg_readlen);
	return($avg_readlen);
}


1;

