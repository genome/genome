package Genome::Model::SmallRna::Command::ParallelWorkflow;

use strict;
use warnings;
use Genome;
use Workflow::Simple;

my $DEFAULT_CLUSTERS = '5000';
my $DEFAULT_CUTOFF = '2';
my $DEFAULT_ZENITH = '5';

class Genome::Model::SmallRna::Command::ParallelWorkflow {
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
		
		minimum_zenith => {
            is => 'String',
            doc => 'Minimum zenith depth for generating clusters',
            default_value => $DEFAULT_ZENITH,
        },
        
		size_bins => {
			is => 'Text',
			doc =>'comma separated list of Min_max read length bins: eg 17_75,17_25',
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
	my $self            = shift;

	my $base_output_dir = $self->output_base_dir;
	my $clusters_number	= $self->input_cluster_number;
	my @size_array		= split (',',$self->size_bins);
	
	
	my $bin 			= \@size_array;

	my $bam_file        = $self->bam_file;

	my $zenith 			= $self->minimum_zenith;
	my $min_mapscore	= $self->subcluster_min_mapzero;
	my $annotations		= $self->annotation_files; 
	my $anno_name		= $self->annotation_name;
	my $number 			= $self->input_cluster_number;
	
####

	my %params = (
		bam_file    						=> $bam_file,
		output_base_dir						=> $base_output_dir,
		annotation_files					=> $annotations,
		annotation_name						=> $anno_name,
		minimum_zenith						=> $zenith,
		read_size_bins						=> $bin,		
		subcluster_min_mapzero				=> $min_mapscore,	
		input_cluster_number				=> $number,
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
	my $output = run_workflow_lsf( $xml_path, %params );
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

