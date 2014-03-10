
package Genome::Model::Tools::Analysis::454::SomaticPipeline;     # rename this when you give the module file a different name <--

#####################################################################################################################################
# SomaticPipeline - Runs the Varscan somatic pipeline on matched tumor-normal data
#					
#	AUTHOR:		Dan Koboldt (dkoboldt@watson.wustl.edu)
#
#	CREATED:	02/25/2009 by D.K.
#	MODIFIED:	02/25/2009 by D.K.
#
#	NOTES:	
#			
#####################################################################################################################################

use strict;
use warnings;

use FileHandle;

use Genome;                                 # using the namespace authorizes Class::Autouse to lazy-load modules under it

class Genome::Model::Tools::Analysis::454::SomaticPipeline {
	is => 'Command',                       
	
	has => [                                # specify the command's single-value properties (parameters) <--- 
		patients_file	=> { is => 'Text', doc => "Tab-delimited file of normal and tumor" },
		data_dir	=> { is => 'Text', doc => "Output directory for 454 data. " },
		output_dir	=> { is => 'Text', doc => "Output directory for VarScan-somatic results e.g. normal_tumor. " },
		varscan_dir	=> { is => 'Text', doc => "Output directory for Varscan in the sample dir", default => "varscan_somatic"},
		varscan_params	=> { is => 'Text', doc => "Default params for varscan somatic", default => "--min-coverage 8 --min-var-freq 0.08 --p-value 0.10 --somatic-p-value 0.05 --strand-filter 1"},
		aligner		=> { is => 'Text', doc => "Aligner to use" },
		reference		=> { is => 'Text', doc => "Reference sequence [default=Hs36 ssaha2]", is_optional => 1 },
		skip_if_output_present	=> { is => 'Text', doc => "Skip if output present", is_optional => 1 },		
		# Make workflow choose 64 bit blades
		lsf_resource => {
		    is_param => 1,
		    default_value => 'rusage[mem=6000,tmp=10000] select[type==LINUX64 && model != Opteron250 && mem>6000 && tmp>10000] span[hosts=1]',
		},
		lsf_queue => {
		    is_param => 1,
		    default_value => $ENV{GENOME_LSF_QUEUE_BUILD_WORKER},
		},
	],
};

sub sub_command_sort_position { 12 }

sub help_brief {                            # keep this to just a few words <---
    "Runs the Varscan somatic pipeline on matched normal-tumor samples"                 
}

sub help_synopsis {
    return <<EOS
This command runs the Varscan somatic pipeline on matched normal-tumor samples
EXAMPLE:	gmt analysis 454 somatic-pipeline --patients-file data/paired-normal-tumor.tsv --output-dir data --aligner ssaha2
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
	my $patients_file = $self->patients_file;
	my $data_dir = $self->data_dir;
	my $varscan_dir = $self->varscan_dir;
	my $aligner = $self->aligner;

	if(!(-e $patients_file))
	{
		die "Error: Samples file not found!\n";
	}

	my $input = new FileHandle ($patients_file);
	my $lineCounter = 0;
	
	while (<$input>)
	{
		chomp;
		my $line = $_;
		$lineCounter++;		
	
		my ($normal_sample_name, $tumor_sample_name) = split(/\t/, $line);			
		
		## Identify BAM files ##
		my $normal_sample_data_dir = $data_dir . "/" . $normal_sample_name;
		my $normal_bam_file = $normal_sample_data_dir . "/" . $aligner . "_out/$normal_sample_name.$aligner.bam";

		my $tumor_sample_data_dir = $data_dir . "/" . $tumor_sample_name;
		my $tumor_bam_file = $tumor_sample_data_dir . "/" . $aligner . "_out/$tumor_sample_name.$aligner.bam";

		if(-e $normal_bam_file && -e $tumor_bam_file)
		{
			print "$tumor_sample_name\t$normal_bam_file\t$tumor_bam_file\n";
			
#			my $varscan_data_dir = $tumor_sample_data_dir . "/$varscan_dir";
			my $varscan_data_dir = $self->output_dir . "/$tumor_sample_name-$normal_sample_name";
			mkdir($varscan_data_dir) if(!(-d $varscan_data_dir));

			open(SCRIPT, ">$varscan_data_dir/script_pipeline.sh") or die "Can't open outfile: $!\n";
			print SCRIPT "#!/gsc/bin/sh\n";
			my $cmd = "";
			
			if(!($self->skip_if_output_present && -e "$varscan_data_dir/varScan.output.snp"))
			{
				print SCRIPT "echo \"\nRunning Varscan somatic...\"\n";
				$cmd = "gmt varscan somatic --normal-bam $normal_bam_file --tumor-bam $tumor_bam_file --output $varscan_data_dir/varScan.output --varscan-params=\"" . $self->varscan_params . "\"";
				print SCRIPT "$cmd\n";				
			}

			if(!($self->skip_if_output_present && -e "$varscan_data_dir/varScan.output.snp.formatted"))
			{
				print SCRIPT "echo \"\nFormatting SNVs...\"\n";
				$cmd = "gmt capture format-snvs --variant $varscan_data_dir/varScan.output.snp --output $varscan_data_dir/varScan.output.snp.formatted";
				print SCRIPT "$cmd\n";
				
			}

			# Format indels
			if(!($self->skip_if_output_present && -e "$varscan_data_dir/varScan.output.indel.formatted"))
			{
				print SCRIPT "echo \"\nFormatting Indels...\"\n";
				$cmd = "gmt capture format-indels --variant $varscan_data_dir/varScan.output.indel --output $varscan_data_dir/varScan.output.indel.formatted";
				print SCRIPT "$cmd\n";				
			}
			
			## Process SNVs ##
			if(!($self->skip_if_output_present && -e "$varscan_data_dir/varScan.output.snp.formatted.Somatic"))
			{
				print SCRIPT "echo \"\nProcessing SNVs...\"\n";
				$cmd = "gmt varscan process-somatic --status-file $varscan_data_dir/varScan.output.snp.formatted";
				print SCRIPT "$cmd\n";	
			}
			
			## Process Indels ##
			if(!($self->skip_if_output_present && -e "$varscan_data_dir/varScan.output.indel.formatted.Somatic"))
			{
				print SCRIPT "echo \"\nProcessing Indels...\"\n";
				$cmd = "gmt varscan process-somatic --status-file $varscan_data_dir/varScan.output.indel.formatted";
				print SCRIPT "$cmd\n";	
			}


			## Filter homopolymer-associated indels ##
			if(!($self->skip_if_output_present && -e "$varscan_data_dir/varScan.output.indel.formatted.Somatic.filter"))
			{
				print SCRIPT "echo \"\nFiltering indels for homopolymers...\"\n";
				$cmd = "gmt somatic monorun-filter --tumor-bam $tumor_bam_file --variant-file $varscan_data_dir/varScan.output.indel.formatted.Somatic --output-file $varscan_data_dir/varScan.output.indel.formatted.Somatic.filter";
				print SCRIPT "$cmd\n";				
			}

			## Filter false positive SNVs ##
			if(!($self->skip_if_output_present && -e "$varscan_data_dir/varScan.output.snp.formatted.Somatic.hc.filter"))
			{
				print SCRIPT "echo \"\nFiltering SNVs for false positives...\"\n";
				$cmd = "gmt somatic filter-false-positives --bam-file $tumor_bam_file --variant-file $varscan_data_dir/varScan.output.snp.formatted.Somatic.hc --output-file $varscan_data_dir/varScan.output.snp.formatted.Somatic.hc.filter --filtered-file $varscan_data_dir/varScan.output.snp.formatted.Somatic.hc.filter.removed";
				print SCRIPT "$cmd\n";	
			}

			## Identify novel SNVs ##			
			if(!($self->skip_if_output_present && -e "$varscan_data_dir/varScan.output.snp.formatted.Somatic.hc.filter.novel"))
			{
				print SCRIPT "echo \"\nIsolating novel SNVs...\"\n";
				$cmd = "gmt annotate lookup-variants --variant $varscan_data_dir/varScan.output.snp.formatted.Somatic.hc.filter --output $varscan_data_dir/varScan.output.snp.formatted.Somatic.hc.filter.novel --filter-out-submitters=\"SNP500CANCER,OMIMSNP,CANCER-GENOME,CGAP-GAI,LCEISEN,ICRCG\"";
				print SCRIPT "$cmd\n";				
			}			

			## Identify novel Indels ##
			if(!($self->skip_if_output_present && -e "$varscan_data_dir/varScan.output.indel.formatted.Somatic.filter.novel"))
			{
				print SCRIPT "echo \"\nIsolating novel indels...\"\n";
				$cmd = "gmt annotate lookup-variants --variant $varscan_data_dir/varScan.output.indel.formatted.Somatic.filter --output $varscan_data_dir/varScan.output.indel.formatted.Somatic.filter.novel --filter-out-submitters=\"SNP500CANCER,OMIMSNP,CANCER-GENOME,CGAP-GAI,LCEISEN,ICRCG\"";
				print SCRIPT "$cmd\n";				
			}

			## Annotate novel SNVs ##			
			if(!($self->skip_if_output_present && -e "$varscan_data_dir/varScan.output.snp.formatted.Somatic.hc.filter.novel.transcript-annotation"))
			{
				print SCRIPT "echo \"\nAnnotating novel SNVs...\"\n";
				$cmd = "gmt annotate transcript-variants --annotation-filter top --variant-file $varscan_data_dir/varScan.output.snp.formatted.Somatic.hc.filter.novel --output-file $varscan_data_dir/varScan.output.snp.formatted.Somatic.hc.filter.novel.transcript-annotation";
				print SCRIPT "$cmd\n";				

				print SCRIPT "echo \"\nUCSC-annotating novel SNVs...\"\n";
				$cmd = "gmt somatic ucsc-annotator --input-file $varscan_data_dir/varScan.output.snp.formatted.Somatic.hc.filter.novel --output-file $varscan_data_dir/varScan.output.snp.formatted.Somatic.hc.filter.novel.ucsc-annotation";
				print SCRIPT "$cmd\n";				
			}			

			## Annotate novel Indels ##
			if(!($self->skip_if_output_present && -e "$varscan_data_dir/varScan.output.indel.formatted.Somatic.filter.novel.transcript-annotation"))
			{
				print SCRIPT "echo \"\nAnnotating novel indels...\"\n";
				$cmd = "gmt annotate transcript-variants --annotation-filter top --variant-file $varscan_data_dir/varScan.output.indel.formatted.Somatic.filter.novel --output-file $varscan_data_dir/varScan.output.indel.formatted.Somatic.filter.novel.transcript-annotation";
				print SCRIPT "$cmd\n";				

				print SCRIPT "echo \"\nUCSC-annotating novel Indels...\"\n";
				$cmd = "gmt somatic ucsc-annotator --input-file $varscan_data_dir/varScan.output.indel.formatted.Somatic.filter.novel --output-file $varscan_data_dir/varScan.output.indel.formatted.Somatic.filter.novel.ucsc-annotation";
				print SCRIPT "$cmd\n";				
			}

			## Tier the SNVs ##
			if(!($self->skip_if_output_present && -e "$varscan_data_dir/varScan.output.snp.formatted.Somatic.hc.filter.novel.tier1"))
			{
				print SCRIPT "echo \"\nTiering novel SNVs...\"\n";
				$cmd = "gmt somatic tier-variants --variant-file $varscan_data_dir/varScan.output.snp.formatted.Somatic.hc.filter.novel --tier1-file $varscan_data_dir/varScan.output.snp.formatted.Somatic.hc.filter.novel.tier1 --tier2-file $varscan_data_dir/varScan.output.snp.formatted.Somatic.hc.filter.novel.tier2 --tier3-file $varscan_data_dir/varScan.output.snp.formatted.Somatic.hc.filter.novel.tier3 --tier4-file $varscan_data_dir/varScan.output.snp.formatted.Somatic.hc.filter.novel.tier4 --transcript-annotation-file $varscan_data_dir/varScan.output.snp.formatted.Somatic.hc.filter.novel.transcript-annotation --ucsc-file $varscan_data_dir/varScan.output.snp.formatted.Somatic.hc.filter.novel.ucsc-annotation";
				print SCRIPT "$cmd\n";				
			}			

			## Tier the Indels ##
			if(!($self->skip_if_output_present && -e "$varscan_data_dir/varScan.output.indel.formatted.Somatic.filter.novel.tier1"))
			{
				print SCRIPT "echo \"\nTiering novel indels...\"\n";
				$cmd = "gmt somatic tier-variants --variant-file $varscan_data_dir/varScan.output.indel.formatted.Somatic.filter.novel --tier1-file $varscan_data_dir/varScan.output.indel.formatted.Somatic.filter.novel.tier1 --tier2-file $varscan_data_dir/varScan.output.indel.formatted.Somatic.filter.novel.tier2 --tier3-file $varscan_data_dir/varScan.output.indel.formatted.Somatic.filter.novel.tier3 --tier4-file $varscan_data_dir/varScan.output.indel.formatted.Somatic.filter.novel.tier4 --transcript-annotation-file $varscan_data_dir/varScan.output.indel.formatted.Somatic.filter.novel.transcript-annotation --ucsc-file $varscan_data_dir/varScan.output.indel.formatted.Somatic.filter.novel.ucsc-annotation";
				print SCRIPT "$cmd\n";				
			}	


			close(SCRIPT);

			system("bsub -q $ENV{GENOME_LSF_QUEUE_BUILD_WORKER} -R\"select[type==LINUX64 && model != Opteron250 && mem>6000 && tmp>20000] rusage[mem=6000]\" -M 6000000 -u \'" . Genome::Sys->username . "\' -oo $varscan_data_dir/script_pipeline.sh.out sh $varscan_data_dir/script_pipeline.sh");		
		}
	}
	
	close($input);

	return 1;                               # exits 0 for true, exits 1 for false (retval/exit code mapping is overridable)
}



1;

