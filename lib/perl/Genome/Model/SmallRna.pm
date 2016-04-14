package Genome::Model::SmallRna;
use strict;
use warnings;
BEGIN { $INC{"Genome/Model/Build/SmallRna.pm"} = 1; $INC{"Genome/ProcessingProfile/SmallRna.pm"} = 1; $INC{"Genome/Model/Command/Define/SmallRna.pm"} = 1; };
use Genome;

# DEFAULTS
my $DEFAULT_CLUSTERS = '5000';
my $DEFAULT_CUTOFF = '2';
my $DEFAULT_ZENITH = '5';
my $DEFAULT_MIN_DEPTH = '1';
my $DEFAULT_BIN 	= '17_70';


class Genome::Model::SmallRna {
    is  => 'Genome::ModelDeprecated',
    has => [
        ref_model_id => {
            via => 'ref_model',
            to => 'id',
        },
    ],
    has_input => [
        ref_model => {
            is    => 'Genome::Model::ReferenceAlignment',
            doc   => 'Ref model with/without Coverage ROI for downstream smallrna analysis ',
        },
        ref_model_coverage => {
        	is_optional => 1,
            	is    => 'Genome::Model::ReferenceAlignment',
            	doc   => 'Ref model with coverage against second ROI list   ',
        },        
    ],
    has_param => [
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
            is_optional => 1,
            doc => 'Minimum zenith depth for generating clusters',
            default_value => $DEFAULT_ZENITH,
        },
        minimum_depth => {
            is => 'String',
            is_optional => 1,
            doc => 'Minimum depth to filter coverage',
            default_value => $DEFAULT_MIN_DEPTH,
        },
        normalization_bin => {
			is => 'Text',
			doc =>'Head bin to normalize by: eg 17_70 ',
            default_value => $DEFAULT_BIN,
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

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
TO DO
EOS
}



sub _resolve_subject {
    my $self = shift;
   # $DB::single = 1;
    my $subject = $self->_infer_subjects_from_ref_model();
    return $subject;
}


sub _infer_subjects_from_ref_model {
    my $self = shift;
    my $subject;
  #  my $input_model = $self->ref_model_id,
	$subject = $self->ref_model->subject;
	return $subject;
}

sub help_detail_for_create_profile {
    return <<EOS
  TO DO
EOS
}

sub help_manual_for_define_model {
    return <<EOS
TO DO
EOS
}

sub _resolve_workflow_for_build {
    my $self = shift;
    my $build = shift;

    my $workflow = Genome::WorkflowBuilder::DAG->from_xml_filename(__FILE__ . '.xml');

    my $log_directory = $build->log_directory;
    $workflow->recursively_set_log_dir($log_directory);
    
    $workflow->name($build->workflow_name);
  	
    my $data_directory = $build->data_directory;    
    my $ref_build      = $build->ref_build;
  
    if ($self->ref_model->region_of_interest_set_name)
     {
		my $coverage_dir_name=$ref_build->reference_coverage_directory;
		unless (-e $coverage_dir_name) {
        $self->error_message("Coverage dir $coverage_dir_name does not exist!");
        die $self->error_message;
    	}			
		my $symlink = "ln -s ".$coverage_dir_name." ".$data_directory."/Coverage_".$self->ref_model->id;
		print "Creating Symbolic Link to the Coverage dir for Model 1 "."\n";
		print $symlink."\n";
		system  ($symlink);
     	}
 		else{
                print "No Coverage info found for reference model 1"."\n";
        	}
        
    if ($self->ref_model_coverage && $self->ref_model_coverage->region_of_interest_set_name)
     {			
     			my $second_ref_build =$self->ref_model_coverage->last_succeeded_build;
                my $second_cov_dir = $second_ref_build->reference_coverage_directory;
                
                unless (-e $second_cov_dir) {
       		    $self->error_message("Coverage dir $second_cov_dir does not exist!");
        		die $self->error_message;
    			}	
                
                my $symlink_2 = "ln -s ".$second_cov_dir." ".$data_directory."/Coverage_".$self->ref_model_coverage->id ;
				print "Creating Symbolic Link to the Coverage dir for Model 2 "."\n";
				print $symlink_2."\n";
				system  ($symlink_2);
     }
 	else
 	{
               	print "No Coverage info found for reference model 2"."\n";
    }
    return $workflow;
}

sub map_workflow_inputs {
    my $self = shift;
    my $build = shift;

    my @inputs = ();
    unless ($self->subject)
    {
    my $subject = $self->_resolve_subject();
    }
	
    my $data_directory = $build->data_directory;
    my $ref_build      = $build->ref_build;
    my $bam_file       = $ref_build->whole_rmdup_bam_file;

    unless (-e $bam_file) {
        $self->error_message("Bam file $bam_file does not exist!");
        die $self->error_message;
    }
	
    my @size_array = split (',', $self->size_bins);
    my $bin = \@size_array;
    my $normalized_output_dir = $data_directory .'/'.$self->normalization_bin;
    Genome::Sys->create_directory($normalized_output_dir);
    my $normalized_filtered_bam = $normalized_output_dir . '/'.$self->normalization_bin . '.bam';

    push @inputs, normalization_bin       => $self->normalization_bin;
    push @inputs, normalized_filtered_bam => $normalized_filtered_bam;
    push @inputs, bam_file                => $bam_file;
    push @inputs, output_base_dir         => $data_directory;
    push @inputs, annotation_files        => $self->annotation_files;
    push @inputs, annotation_name         => $self->annotation_name;
    push @inputs, minimum_zenith          => $self->minimum_zenith;
    push @inputs, minimum_depth           => $self->minimum_depth;
    push @inputs, size_bins               => $bin;
    push @inputs, subcluster_min_mapzero  => $self->subcluster_min_mapzero;
    push @inputs, input_cluster_number    => $self->input_cluster_number;

    return @inputs;
}

1;
