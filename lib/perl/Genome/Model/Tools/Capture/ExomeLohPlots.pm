
package Genome::Model::Tools::Capture::ExomeLohPlots;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# ExomeLohPlots - Run exome-based copy number analysis in parallel on somatic-variation group models
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


class Genome::Model::Tools::Capture::ExomeLohPlots {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		group_id		=> { is => 'Text', doc => "ID of somatic-variation model group" , is_optional => 0},
		output_dir	=> { is => 'Text', doc => "An output directory to hold LOH results (can be same as copy number dir)" , is_optional => 0},
		min_coverage	=> { is => 'Text', doc => "Minimum coverage to include variant" , is_optional => 0, default => 20},
		min_var_freq	=> { is => 'Text', doc => "Minimum variant allele frequency in normal for het" , is_optional => 0, default => 40},
		max_var_freq	=> { is => 'Text', doc => "Maximum variant allele frequency in normal for het" , is_optional => 0, default => 60},
		reference	=> 	{ is => 'Text', doc => "Reference to use for bam-readcounts-based filters; defaults to build 37" , is_optional => 1, example_values => ['/gscmnt/sata420/info/model_data/2857786885/build102671028/all_sequences.fa']},
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Runs LOH analysis using VarScan calls on a model group of somatic-variation models"                 
}

sub help_synopsis {
    return <<EOS
This command runs LOH analysis using VarScan calls on a model group of somatic-variation models
EXAMPLE:	gmt capture exome-loh-plots --group-id 3328 --output-dir varscan_copynumber
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
	
		if(!(-d $output_dir))
		{
			mkdir($output_dir) if(!(-d $output_dir));			
		}
		
		print "Processing $output_dir\n";
		print "Tumor: $tumor_model_dir\n";
		print "Normal: $normal_model_dir\n";
		## Get succeeded build dir ##
		my $last_build_dir = "";
		my @builds = $model->builds;
		
		foreach my $build (@builds)
		{
			if($build->status eq 'Succeeded')
#			if($build->status ne 'Running')
			{
				$last_build_dir = $build->data_directory;
			}
		}

		if($last_build_dir)
		{
			print "Somatic: $last_build_dir\n";
		
			my $germline_hc_file = `ls $last_build_dir/variants/snv/varscan-somatic*/varscan-high-confidence*/snvs.Germline.hc`;
			chomp($germline_hc_file);

			print "Parsing Normal Hets...\n";

			## Process germline file ##			
			my $germline_outfile = $output_dir . "/varScan.snp.Germline.het";
			if($germline_hc_file)
			{
				my $num_germline_hets = parse_hets($self, $germline_hc_file, $germline_outfile);
				print "$num_germline_hets Germline hets included\n";
			}
			
			my $loh_hc_file = `ls $last_build_dir/variants/snv/varscan-somatic*/varscan-high-confidence*/snvs.LOH.hc`;
			chomp($loh_hc_file);

			my $loh_outfile = $output_dir . "/varScan.snp.LOH.het";			

			if($loh_hc_file)
			{
				my $num_loh_hets = parse_hets($self, $loh_hc_file, $loh_outfile);
				print "$num_loh_hets LOH calls included\n";
			}

			## Combine into single file ##
			my $combined_outfile = $output_dir . "/varscan.snp.Normal.het";

			if(!(-e $combined_outfile))
			{
                            #my $current_len = `cat $combined_outfile | wc -l`; chomp($current_len);
                            print "Combining and sorting...\n";
                            system("cat $germline_outfile $loh_outfile >$combined_outfile");
                            system("gmt capture sort-by-chr-pos --input $combined_outfile --output $combined_outfile");
                            #my $new_len = `cat $combined_outfile | wc -l`; chomp($new_len);
			}


			## Run false-positive filter ##
			if(!(-e "$combined_outfile.fpfilter"))
			{
				print "Running false-positive filter...\n";
				my $cmd = "gmt somatic filter-false-positives --variant-file $combined_outfile --output-file $combined_outfile.fpfilter --bam-file $normal_bam --reference " . $self->reference;
				system("bsub -q $ENV{GENOME_LSF_QUEUE_BUILD_WORKER} -R\"select[mem>8000] rusage[mem=8000]\" -M 8000000 -oo $output_dir/fperr.log $cmd");
			}
			else
			{
				my $cmd = "gmt varscan loh-segments --variant-file $combined_outfile.fpfilter --output-basename $combined_outfile.fpfilter.plot --varscan-cn-basename varScan.output.copynumber.called.cbs";
				system("bsub -q $ENV{GENOME_LSF_QUEUE_BUILD_WORKER} $cmd");				
			}



#			exit(0);
		}
		else
		{
			warn "No Succeeded build dir for $model_id ($model_name)\n";
		}

#		exit(0);

	}	
	
	return 1;
}


#############################################################
# parse_germline_hets
#
#############################################################

sub parse_hets
{
	my ($self, $FileName, $OutFile) = @_;

	open(OUTFILE, ">$OutFile") or die "Can't open outfile: $!\n";
	
	my $input = new FileHandle ($FileName);
	my $lineCounter = 0;
	
	my $num_hets_included = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		
	
		my @lineContents = split(/\t/, $line);
		my $numContents = @lineContents;
		
		my $normal_cov = my $normal_freq = my $tumor_cov = 0;
		
		if($lineContents[2] =~ /[0-9]/)
		{
			## Already in chrom-start-stop-format (this is unexpected) ##
			$normal_cov = $lineContents[5] + $lineContents[6];
			$normal_freq = $lineContents[7];
			$normal_freq =~ s/\%//;
			$tumor_cov = $lineContents[9] + $lineContents[10];
			
			if($normal_cov >= $self->min_coverage)
			{
				if($normal_freq >= $self->min_var_freq && $normal_freq <= $self->max_var_freq)
				{
					print OUTFILE "$line\n";
					$num_hets_included++;
				}
			}
		}
		else
		{
			$normal_cov = $lineContents[4] + $lineContents[5];
			$normal_freq = $lineContents[6];
			$normal_freq =~ s/\%//;
			$tumor_cov = $lineContents[8] + $lineContents[9];
			
			if($normal_cov >= $self->min_coverage)
			{
				if($normal_freq >= $self->min_var_freq && $normal_freq <= $self->max_var_freq)
				{
					## Reformat the line ##
					my $new_line = join("\t", $lineContents[0], $lineContents[1], $lineContents[1]);
					for(my $colCounter = 2; $colCounter < $numContents; $colCounter++)
					{
						$new_line .= "\t" . $lineContents[$colCounter];
					}

					print OUTFILE "$new_line\n";
					$num_hets_included++;
				}
			}
		}
		
#		return($num_hets_included) if($lineCounter > 1000);

	}
	
	close($input);	

	close(OUTFILE);
	
	return($num_hets_included);
}


1;

