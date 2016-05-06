
package Genome::Model::Tools::Capture::ExomeCopyNumber;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# ExomeCopyNumber - Run exome-based copy number analysis in parallel on somatic-variation group models
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


class Genome::Model::Tools::Capture::ExomeCopyNumber {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		group_id		=> { is => 'Text', doc => "ID of somatic-variation model group" , is_optional => 0},
		output_dir	=> { is => 'Text', doc => "An output directory to hold copy number results" , is_optional => 1},
		varscan_params	=> { is => 'Text', doc => "Parameters to pass to VarScan copynumber" , is_optional => 1, default => '--min-coverage 20 --min-segment-size 25 --max-segment-size 100'},
		reference	=> 	{ is => 'Text', doc => "Reference to use for bam-readcounts-based filters" , is_optional => 1, example_values => ['/gscmnt/gc4096/info/model_data/2857786885/build102671028/all_sequences.fa']},
		test_only	=> { is => 'Text', doc => "Set to 1 to see how the command would run without doing anything" , is_optional => 1},
		run_limit	=> { is => 'Text', doc => "If specified, will only launch this many copy number in parallel jobs" , is_optional => 1},		
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Runs exome-based copy number on a model group of somatic-variation models"                 
}

sub help_synopsis {
    return <<EOS
This command runs exome-based copy number on a model group of somatic-variation models
EXAMPLE:	gmt capture exome-copy-number --group-id 3328 --output-dir varscan_copynumber
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
	my $output_dir = $self->output_dir;

	## Get the models in each model group ##

	my $model_group = Genome::ModelGroup->get($group_id);
	my @models = $model_group->models;
	my $debug_counter = 0;

	foreach my $model (@models)
	{
		$stats{'models_in_group'}++;
		
		my $model_id = $model->genome_model_id;
		my $model_name = $model->name;
		my $subject_name = $model->subject_name;
		$subject_name = "Model" . $model_id if(!$subject_name);

		## Get normal and tumor model ##
		
		my $normal_model = $model->normal_model;
		my $tumor_model = $model->tumor_model;

		## Get Model IDs ##
		
		my $normal_model_id = $normal_model->id;
		my $tumor_model_id = $tumor_model->id;

		## Get TCGA-Compliant Subject Names ##

		my $normal_subject_name = $normal_model->subject_name;
		my $tumor_subject_name = $tumor_model->subject_name;

		my ($sample_prefix) = split(/\-/, $normal_subject_name);
		$normal_subject_name =~ s/$sample_prefix/TCGA/;
		
		($sample_prefix) = split(/\-/, $tumor_subject_name);
		$tumor_subject_name =~ s/$sample_prefix/TCGA/;

		my @temp = split(/\-/, $tumor_subject_name);
		my $patient_id = join("-", $temp[0], $temp[1], $temp[2]);

		my $tumor_sample = $tumor_model->subject_name;
		my $normal_sample = $normal_model->subject_name;
		my $tumor_model_dir = $tumor_model->last_succeeded_build_directory;
		my $normal_model_dir = $normal_model->last_succeeded_build_directory;
		my $tumor_bam = `ls $tumor_model_dir/alignments/*.bam | head -1`; chomp($tumor_bam);
		my $normal_bam = `ls $normal_model_dir/alignments/*.bam | head -1`; chomp($normal_bam);

		my $output_dir = $self->output_dir . "/" . $tumor_sample . "-" . $normal_sample;
	
		print join("\t", $tumor_sample, $normal_sample) . "\n";
		
		if($self->test_only)
		{
			## Get the flagstat ##
#			print "Getting Flagstat...\n";	
			my %normal_flagstat = get_flagstat($normal_bam);
			my %tumor_flagstat = get_flagstat($tumor_bam);
	
#			print "Computing Average Read Length...\n";
			my $normal_readlen = avg_read_len($normal_bam);
			my $tumor_readlen = avg_read_len($tumor_bam);
			
			## Determine the total unique GBP ##
			
			my $normal_unique_bp = ($normal_flagstat{'mapped'} - $normal_flagstat{'duplicates'}) * $normal_readlen;
			my $tumor_unique_bp = ($tumor_flagstat{'mapped'} - $tumor_flagstat{'duplicates'}) * $tumor_readlen;
	
			my $normal_tumor_ratio = $normal_unique_bp / $tumor_unique_bp;
			$normal_tumor_ratio = sprintf("%.4f", $normal_tumor_ratio);
	
#			print "Normal: $normal_readlen ==> $normal_unique_bp\n";
#			print "Tumor: $tumor_readlen ==> $tumor_unique_bp\n";
#			print "Ratio: $normal_tumor_ratio\n";
			print join("\t", $tumor_sample, $normal_sample, $tumor_unique_bp, $normal_unique_bp, $normal_tumor_ratio) . "\n";


		}
		else
		{
			my $file_count = 0;

			## Check output dir ##
			if(!(-d $output_dir))
			{
				mkdir($output_dir) if(!(-d $output_dir));				
			}
			else
			{
				## If output dir already existed, count files ##
				$file_count = `ls $output_dir/varScan.output*copynumber 2>/dev/null | wc -l`;
				chomp($file_count);
			}

			## Proceed if not all files are present ##
	
			if($file_count < 20)
			{
				my $cmd = "gmt varscan copy-number-parallel --output $output_dir/varScan.output --normal-bam $normal_bam --tumor-bam $tumor_bam --varscan-params=\"" . $self->varscan_params . "\"";
				$cmd .= " --reference " . $self->reference if($self->reference);

				$stats{'models_to_run'}++;

				if(!$self->run_limit || $stats{'models_launched'} < $self->run_limit)
				{
					system($cmd);
					$stats{'models_launched'}++;					
				}

			}
			else
			{
				warn "Skipping $tumor_sample because it has $file_count files...\n";
				$stats{'models_skipped'}++;
			}
			
		}

#		exit(0);

	}	
	
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

