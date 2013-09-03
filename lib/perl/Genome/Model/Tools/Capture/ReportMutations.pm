
package Genome::Model::Tools::Capture::ReportMutations;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# ReportMutations - Build Genome Models for Capture Datasets
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
use Genome::Model::Tools::Capture::Helpers 'iupac_to_base';

## Declare global statistics hash ##

my %stats = ();

class Genome::Model::Tools::Capture::ReportMutations {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		group_id		=> { is => 'Text', doc => "ID of model group" , is_optional => 0},
		show_builds		=> { is => 'Text', doc => "If set to 1, show all builds of somatic models" , is_optional => 1},
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Report mutations for a capture somatic model group"                 
}

sub help_synopsis {
    return <<EOS
This command reports mutations for a capture somatic model group
EXAMPLE:	gmt capture report-mutations --group-id 2788
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

	
	## Reset counters ##
	my %num_mutated_samples = my %num_mutations = ();	
	
	## Keep stats in a single hash ##
	
	my %stats = ();
	
	## Get the models in each model group ##

	my $model_group = Genome::ModelGroup->get($group_id);
	my @models = $model_group->models; 

	foreach my $model (@models)
	{
		$stats{'models_in_group'}++;
		
		my $model_id = $model->genome_model_id;
		my $subject_name = $model->subject_name;
		
		my $build_dir = "";
		my $model_status = "New";
		my $lsb_stats = "";

		my $num_builds = 0;		

		my $build_ids = my $build_statuses = "";
		my @builds = $model->builds;

		## Get last succeeded build ##
		my $last_succeeded_build;
		if($self->show_builds)
		{
			foreach my $build (@builds)
			{
				$num_builds++;
				$build_ids .= "," if($build_ids);
				$build_ids .= $build->id;
				$build_statuses .= "," if($build_statuses);
				$build_statuses .= $build->status;
				if($build->status eq "Succeeded")
				{
					$last_succeeded_build = $build;
				}
			}
		}
		
		if(@builds)
		{
			$model_status = "Run";

			if($model->last_succeeded_build_directory)
			{
				my %sample_stats = ();
				
				$model_status = "Done";

				$stats{'models_finished'}++;
				
				$build_dir = $model->last_succeeded_build_directory;


				## Process the results ##
				
				$sample_stats{'tier1_snvs'} = $sample_stats{'tier1_indels'} = $sample_stats{'tier1_gatk'} = 0;
				my $tier1_snv_file = "$build_dir/merged.somatic.snp.filter.novel.tier1";
				my $tier1_indel_file = "$build_dir/merged.somatic.indel.filter.tier1";
				my $tier1_gatk_file = "$build_dir/gatk.output.indel.formatted.Somatic.tier1";				
				my $annotation_snv_file = "$build_dir/annotation.somatic.snp.transcript";
				my $annotation_indel_file = "$build_dir/annotation.somatic.indel.transcript";
				my $annotation_gatk_file = "$build_dir/annotation.somatic.gatk-indel.transcript";

				my $mutation_file_list = my $annotation_file_list = "";

				if(-e $tier1_snv_file)
				{
					$sample_stats{'tier1_snvs'} = `cat $tier1_snv_file | wc -l`;
					chomp($sample_stats{'tier1_snvs'});
					## Append SNV files to lists ##
					$mutation_file_list .= "," if($mutation_file_list);
					$mutation_file_list .= $tier1_snv_file;
					$annotation_file_list .= "," if($annotation_file_list);
					$annotation_file_list .= $annotation_snv_file;
				}

				if(-e $tier1_indel_file)
				{
					$sample_stats{'tier1_indels'} = `cat $tier1_indel_file | wc -l`;
					chomp($sample_stats{'tier1_indels'});
					## Append indel files to lists ##
					$mutation_file_list .= "," if($mutation_file_list);
					$mutation_file_list .= $tier1_indel_file;
					$annotation_file_list .= "," if($annotation_file_list);
					$annotation_file_list .= $annotation_indel_file;

				}

				if(-e $tier1_gatk_file)
				{
					$sample_stats{'tier1_gatk'} = `cat $tier1_gatk_file | wc -l`;
					chomp($sample_stats{'tier1_gatk'});
					## Append SNV files to lists ##
					$mutation_file_list .= "," if($mutation_file_list);
					$mutation_file_list .= $tier1_gatk_file;
					$annotation_file_list .= "," if($annotation_file_list);
					$annotation_file_list .= $annotation_gatk_file;

				}

				my %gene_mutation_counts = match_mutations_to_genes($mutation_file_list, $annotation_file_list);

				foreach my $gene (keys %gene_mutation_counts)
				{
					$num_mutated_samples{$gene}++;
					$num_mutations{$gene} += $gene_mutation_counts{$gene};
				}

				$lsb_stats = join("\t", $sample_stats{'tier1_snvs'}, $sample_stats{'tier1_indels'}, $sample_stats{'tier1_gatk'});
			}
			else
			{
				$stats{'models_running'}++;
			}
		}
		
		if($self->show_builds)
		{
			print join("\t", $subject_name, $model_id, $model_status, $num_builds, $build_ids, $build_statuses, $build_dir) . "\n";
		}
		else
		{
			print join("\t", $subject_name, $model_id, $model_status, $lsb_stats, $build_dir) . "\n";
		}


	}	
	
#	print $stats{'models_in_group'} . " models in group\n" if($stats{'models_in_group'});
#	print $stats{'models_running'} . " models running\n" if($stats{'models_running'});
#	print $stats{'models_finished'} . " models finished\n" if($stats{'models_finished'});


	## Build an array of gene mutation counts ##
	my @gene_sample_mutations = ();
	my $numEntries = 0;

	my $brc50_gene_results = "";

	foreach my $gene (keys %num_mutated_samples)
	{
		my $num_samples = $num_mutated_samples{$gene};
		my $num_mutations = $num_mutations{$gene};
		
		$gene_sample_mutations[$numEntries] = join("\t", $gene, $num_samples, $num_mutations);

		if($gene eq "PIK3CA" || $gene eq "MT-ND5" || $gene eq "TP53" || $gene eq "TTN" || $gene eq "SYNE1" || $gene eq "MAP3K1" || $gene eq "ATR" || $gene eq "MYST3" || $gene eq "BIRC6" || $gene eq "RB1")
		{
			$brc50_gene_results .= join("\t", $gene, $num_samples, $num_mutations) . "\n";			
		}


		$numEntries++;
	}

	print "RESULTS FOR BRC50 TOP MUTATED GENES\n";
	print "$brc50_gene_results\n";

	## Sort by number of samples ##

	@gene_sample_mutations = sort byNumSamples (@gene_sample_mutations);

	sub byNumSamples
	{
		my ($gene_a, $num_a) = split(/\t/, $a);
		my ($gene_b, $num_b) = split(/\t/, $b);
		$num_b <=> $num_a;
	}

	## Print the top 20 ##
	
	print "TOP 20 MUTATED GENES\n";
	for(my $printCounter = 0; $printCounter < 20; $printCounter++)
	{
		print "$gene_sample_mutations[$printCounter]\n";
	}

}




################################################################################################
# Execute - the main program logic
#
################################################################################################

sub match_mutations_to_genes
{                               # replace with real execution logic.
	my ($mutation_files, $annotation_files) = @_;

	my %gene_mutation_counts = ();
	
	my @mutation_files = split(/\,/, $mutation_files);
	my @annotation_files = split(/\,/, $annotation_files);
	
	my $numFiles = @mutation_files;

	my %mutation_counted = ();
	
	for (my $fileCounter = 0; $fileCounter < $numFiles; $fileCounter++)
	{
		my $mutation_file = $mutation_files[$fileCounter];
		my $annotation_file = $annotation_files[$fileCounter];
		my %mutation_list = parse_both_files($mutation_file, $annotation_file);

		## Parse the mutations that were returned ##
		
		foreach my $gene_mutation (keys %mutation_list)
		{
			my ($chrom, $chr_start, $chr_stop, $ref, $var, $gene) = split(/\t/, $gene_mutation);
			my $key = join("\t", $chrom, $chr_start, $chr_stop, $ref, $var);
			
			if(!$mutation_counted{$key})
			{
				$gene_mutation_counts{$gene}++;
				$mutation_counted{$key} = 1;
			}
		}
	}


	return(%gene_mutation_counts);
}

################################################################################################
# Execute - the main program logic
#
################################################################################################

sub parse_both_files
{
	my ($mutation_file, $annotation_file) = @_;
	my %gene_mutations = ();

	## Load annotation ##
	my %annotation = ();
	my $input = new FileHandle ($annotation_file);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		my ($chrom, $chr_start, $chr_stop, $ref, $var) = split(/\t/, $line);
		$var = iupac_to_base($ref, $var);

		my @lineContents = split(/\t/, $line);
		my $var_type = $lineContents[5];
		my $gene_name = $lineContents[6];
		my $trv_type = $lineContents[13];
		my $tx_pos = $lineContents[14];
		my $aa_change = $lineContents[15];
		my $key = join("\t", $chrom, $chr_start, $chr_stop, $ref, $var);
		$annotation{$key} = join("\t", $var_type, $gene_name, $trv_type, $tx_pos, $aa_change);
	}

	close($input);


	## Parse mutations ##
	$input = new FileHandle ($mutation_file);
	$lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;
		
		my ($chrom, $chr_start, $chr_stop, $ref, $var) = split(/\t/, $line);
		$var = iupac_to_base($ref, $var);

		my $key = join("\t", $chrom, $chr_start, $chr_stop, $ref, $var);

		if($annotation{$key})
		{
			my ($var_type, $gene_name, $trv_type, $tx_pos, $aa_change) = split(/\t/, $annotation{$key});
			$gene_mutations{$key . "\t" . $gene_name}++;
		}
	}

	close($input);

	return(%gene_mutations);
}

1;

