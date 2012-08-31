package Genome::Model::Tools::HmpShotgun::MergeAlignmentsMulti;

use strict;
use warnings;

use Genome;
use Workflow;
use IO::File;

class Genome::Model::Tools::HmpShotgun::MergeAlignmentsMulti {
    is  => ['Command'],
    has => [
        working_directory => {
            is  => 'ARRAY',
            is_input => '1',
            doc => 'The working directory.',
        },
        concise_files => {
        	is  => 'ARRAY',
            is_input => '1',
            doc => 'The reads to align.',
        },
        paired_end1_concise_file => {
        	is  => 'String',
            is_output => '1',
            is_optional => '1',
            doc => 'The resulting alignment summary of multiple hits.',
        },
        paired_end2_concise_file => {
        	is  => 'String',
            is_output => '1',
            is_optional => '1',
            doc => 'The resulting alignment summary of multiple hits.',
        },
        
	],
    has_param => [
           lsf_resource => {
           default_value => 'select[model!=Opteron250 && type==LINUX64] rusage[mem=4000]',
           },
    ],
};

sub help_brief {
    'Align reads against a given metagenomic reference.';
}

sub help_detail {
    return <<EOS
    Align reads against a given metagenomic reference.
EOS
}

sub execute {
    my $self = shift;

    $self->dump_status_messages(1);
    $self->status_message(">>>Running MergeAlignmentsMulti at ".UR::Time->now);
    #my $model_id = $self->model_id;
    my $concise_files_ref = $self->concise_files;
    my @concise_files = @$concise_files_ref;
    
    #get working directory from  parallelized inputs
    my $working_directory_ref = $self->working_directory;
    my @working_directory_list = @$working_directory_ref;
    my $working_directory = $working_directory_list[0];
    
    #my $working_directory = $self->working_directory;
    $self->status_message("Working directory: ".$working_directory);
    
    #grep the items based on the sequence name
    my @pe1 = grep(/^.*s_\d_1.*/,@concise_files);
    my @pe2 = grep(/^.*s_\d_2.*/,@concise_files);
   
	my $pe1_output_file = $working_directory."/pe1_combined_concise_file.txt";
	my $pe2_output_file = $working_directory."/pe2_combined_concise_file.txt";
    
    $self->status_message("Merging ".join("\n",@pe1)." into $pe1_output_file");
    $self->merge($pe1_output_file, \@pe1);
    $self->status_message("Merging ".join("\n",@pe2)." into $pe2_output_file");
    $self->merge($pe2_output_file, \@pe2);
    
    $self->paired_end1_concise_file($pe1_output_file);
    $self->paired_end2_concise_file($pe2_output_file);
    $self->status_message("<<<Completed MergeAlignmentsMulti for testing at at ".UR::Time->now);
    return 1;
}

sub merge {
    
    my $self = shift;
    my $combined_file = shift;
    my $files_to_merge_ref = shift;
    my @concise_files = @$files_to_merge_ref;  
    
    my @expected_output_files = ($combined_file);
    
    my $rv_check = Genome::Sys->are_files_ok(input_files=>\@expected_output_files);
    
    if (defined($rv_check)) {
	    if ($rv_check == 1) {
	    	#shortcut this step, all the required files exist.  Quit.
	    	$self->status_message("Skipping this step.  File exists: $combined_file.  If you would like to regenerate these files, remove them and rerun.");
	   	    return 1;
	    }
    }
	
    if (scalar(@concise_files) == 0) {
	 $self->error_message("*** No files to merge.  Quitting.");
	 return;
    } elsif ( scalar(@concise_files) == 1) {
		$self->status_message("Only one alignment file is present.  Not merging, only copying.");
		my $cp_cmd = "cp ".$concise_files[0]." ".$combined_file;
		my $rv_cp = Genome::Sys->shellcmd(cmd=>$cp_cmd);
		if ($rv_cp != 1) {
			$self->error_message("<<<Failed MergeAligments. Copy failed.  Return value: $rv_cp");
			return;
		} 
    } else {
		
	    $self->status_message("Merging files: ".join("\n",@concise_files) );
	    $self->status_message("Destination file: ".$combined_file);    
        
            my $rv_cat = Genome::Sys->cat(input_files=>\@concise_files,output_file=>$combined_file);        
		if ($rv_cat != 1) {
	    	$self->error_message("<<<Failed MergeAlignments on cat.  Return value: $rv_cat");
	    	return;
	    }	    
    
	}
    
    Genome::Sys->mark_files_ok(input_files=>\@expected_output_files);
    return 1;   
}
 
1;
