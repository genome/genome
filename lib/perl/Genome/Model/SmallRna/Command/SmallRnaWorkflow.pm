package Genome::Model::SmallRna::Command::SmallRnaWorkflow;

use strict;
use warnings;
use Genome;
use Workflow;
use Workflow::Simple;

my $DEFAULT_CLUSTERS = '5000';
my $DEFAULT_CUTOFF = '2';
my $DEFAULT_ZENITH = '5';
my $DEFAULT_MIN_DEPTH = '1';

class Genome::Model::SmallRna::Command::SmallRnaWorkflow {
	is        => ['Genome::Model::SmallRna::Command::Base'],
	has_input => [
		bam_file => {
			is  => 'Text',
			doc => 'Input BAM File',
		},
		output_base_dir => {
			is  => 'Text',
			doc => 'Path of the directory to write output file',
		},

		annotation_files => {
			is => 'Text',
			doc =>'Comma separated list of input BED files',
		},
		
		annotation_name => {
        		    is => 'String',
            		    doc => 'Comma delimited list of the Annotation Tracks. Should be in the same order as the list of annotation bed files.',
        	},
		
		normalized_bam_file => {
			    is =>'Text',
			    doc => 'Bam file for head bin normalization',
		},

	   	minimum_zenith => {
            		is => 'String',
            		doc => 'Minimum zenith depth for generating clusters',
            		default_value => $DEFAULT_ZENITH,
        	},
		minimum_depth => {
            		is => 'String',
            		doc => 'Minimum depth to filter coverage',
            		default_value => $DEFAULT_MIN_DEPTH,
        },

        
		read_size_bin => {
			is => 'Text',
			doc =>'Min_max read length bin: eg 17_75',
		},
		
		subcluster_min_mapzero => {
			is        => 'Text',
			is_optional => 1,
			doc       =>'Minimum %MapZero Alignments to call subclusters',
			default_value => $DEFAULT_CUTOFF,

		},
		
		input_cluster_number => {
            is => 'Text',
            is_optional => 1,
            doc => 'Number of TOP Clusters to calculate statistcs',
            default_value => $DEFAULT_CLUSTERS,
	   },
	   
	   

	],
};

sub execute {
	my $self            		= shift;
	my $base_output_dir 		= $self->output_base_dir;
	my $clusters_number		= $self->input_cluster_number;
	my $bin 			= $self->read_size_bin;
	my $bam_file        		= $self->bam_file;
	my $output_dir 			= $base_output_dir .'/'.$bin;

	#print $output_dir."\n";
	Genome::Sys->create_directory($output_dir);

#### FILTER-BAM PARAMETERS####
	my $filtered_bam    		= $output_dir . '/'.$bin.'.bam';
	my $xa_tag_value		= '1';
	
#### CLUSTER-COVERAGE PARAMETERS## 	

	my $zenith 			= $self->minimum_zenith;
	my $coverage_stats  		= $output_dir . '/coverage_stats.tsv';
	my $regions_bed			= $output_dir . '/all_regions.bed'; 

#### STATS-GENERATOR PARAMETERS###

	my $sorted_clusters_bed		= $output_dir . '/top_sorted_clusters.bed';
	my $alignment_stats		= $output_dir . '/alignment_stats.tsv';
	my $subclusters			= $output_dir . '/subclusters.bed'; 
	my $subclusters_intersect	= $output_dir . '/subclusters_intersect.tsv'; 
	my $min_mapscore		= $self->subcluster_min_mapzero;
        my $flagstat_17_70_file		= $self->normalized_bam_file.'.flagstat';
 
### ANNOTATE-CLUSTER PARAMETERS###

	my $annotations			= $self->annotation_files; 
	my $anno_name			= $self->annotation_name;
	my $annotation_output		= $output_dir . '/annotation_intersect.tsv'; 
	
### SPREADSHEET PARAMETERS###
	
	my $spreadsheet			= $output_dir . '/Final_spreadsheet.tsv';
	my $number 			= $self->input_cluster_number;
	
####

	my %params = (
		input_bam_file    					=> $bam_file,
		filtered_bam_file					=> $filtered_bam,
#		xa_tag			  				=> $xa_tag_value,
		read_size_bin						=> $bin,
		zenith_depth						=> $zenith,
		minimum_depth						=> $self->minimum_depth,
		stats_file 						=> $coverage_stats,
		bed_file						=> $regions_bed,
		flagstat_17_70_file					=> $flagstat_17_70_file,
		output_stats_file					=> $alignment_stats,
		output_clusters_file					=> $sorted_clusters_bed,
		output_subclusters_file					=> $subclusters,
		output_subcluster_intersect_file			=> $subclusters_intersect,
		subcluster_min_mapzero					=> $min_mapscore,
		annotation_bed_file					=> $annotations,
		annotation_name						=> $anno_name,
		output_tsv_file						=> $annotation_output,
		output_spreadsheet					=> $spreadsheet,
		input_cluster_number					=> $number
	);

	my $module_path = $self->get_class_object->module_path;
	my $xml_path    = $module_path;
	$xml_path =~ s/\.pm/\.xml/;
	my $workflow = Workflow::Operation->create_from_xml($xml_path);
	
	my @errors   = $workflow->validate;
	unless ( $workflow->is_valid ) {
		die(    'Errors encountered while validating workflow '
			  . $xml_path . "\n"
			  . join( "\n", @errors ) );
	}
	### THIS IS NOT REQUIRED OTHERWISE, BUT ONLY WHEN THIS WORKFLOW IS RUN PARALLEL'
	
	unless (Workflow::Simple->can('run_workflow_lsf'))
	{
		die('Unable to load method');
	}
	my $output = Workflow::Simple::run_workflow_lsf( $xml_path, %params );
	unless ( defined $output ) {
		my @errors = @Workflow::Simple::ERROR;
		for (@errors) {
			print STDERR $_->error . "\n";
		}
		return;
	}
	return 1;
}

1;

__END__

