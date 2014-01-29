package Genome::Model::Tools::HmpShotgun::RefCov;

use strict;
use warnings;

use Genome;
use Workflow;
use IO::File;

class Genome::Model::Tools::HmpShotgun::RefCov {
    is  => ['Command'],
    has => [
        working_directory => {
            is  => 'String',
            is_input => '1',
            doc => 'The working directory.',
        },
        aligned_bam_file => {
        	is  => 'String',
                is_input => '1',
                doc => 'The reference sequence.',
        },
        regions_file => {
        	is  => 'String',
                is_input => '1',
                doc => 'The reads to align.',
        },
        read_count_file => {
        	is  => 'String',
                is_input => '1',
                doc => 'The reads/contig summary.',
        },
        other_hits_file => {
        	is  => 'String',
                is_input => '1',
                doc => 'The other_hits summary.',
        },
        combined_file => {
        	    is  => 'String',
                is_output => '1',
                is_optional => '1',
                doc => 'The resulting alignment.',
        },
        
	],
    has_param => [
           lsf_resource => {
           default_value => 'select[model!=Opteron250 && type==LINUX64] rusage[mem=4000]',
           },
    ],
};

sub help_brief {
    'Run the reference coverage report.';
}

sub help_detail {
    return <<EOS
    Runs the reference coverage report.
EOS
}

sub execute {
    my $self = shift;

    $self->dump_status_messages(1);
    
    $self->debug_message(">>>Running HMP RefCov at ".UR::Context->current->now);
    #my $model_id = $self->model_id;
    $self->debug_message("Aligned Bam File: ".$self->aligned_bam_file);
    $self->debug_message("Regions file: ".$self->regions_file);
    $self->debug_message("Read count file: ".$self->read_count_file);
    
    #$self->debug_message("<<<Completed HMP RefCov for testing at ".UR::Context->current->now);
    #return 1;
    
    #expected output files
    my $stats_file = $self->working_directory."/reports/refcov_stats.txt";
    
    my $readcount_file = $self->read_count_file;
    my $combined_file = $self->working_directory."/reports/combined_refcov.txt";
    
    $self->combined_file($combined_file);
    
    
    my @expected_refcov_output_files = ($stats_file);
    
    $self->debug_message("Output stats file: ".$stats_file);
    
    my $rv_check = Genome::Sys->are_files_ok(input_files=>\@expected_refcov_output_files);
    if ($rv_check) {
    	$self->debug_message("Expected output files exist.  Skipping generation of ref cov stats file.");
    } else {
  
    	my $cmd = "genome-perl5.10 -S gmt ref-cov standard ".$self->aligned_bam_file." ".$self->regions_file." ".$stats_file;    
     														
    	$self->debug_message("Running ref cov report at ".UR::Context->current->now);
    	my $rv = Genome::Sys->shellcmd(cmd=>$cmd);
    	if ($rv == 1) {
    		Genome::Sys->mark_files_ok(input_files=>\@expected_refcov_output_files);
    	}
    	$self->debug_message("RefCov file generated at ".UR::Context->current->now);
    }
    	
    	
    my @expected_output_files = ($combined_file);
    my $rv_output_check = Genome::Sys->are_files_ok(input_files=>\@expected_output_files);
    #The re-check of the ref cov check value (rv_check) is necessary because if the refcov stats file has been regenerated,
    #then the subsequent report files should be regenerated.
    if ($rv_output_check && $rv_check) {
   		$self->debug_message("Expected output files exist.  Skipping generation of the combined file.");
   		$self->debug_message("<<<Completed RefCov at ".UR::Context->current->now);
   		return 1;
    } else {
    	$self->debug_message("The previous ref cov stats file may have been regenerated.  Attempting to explicitly delete: $combined_file" );
    	unlink($combined_file);
    }
    
    $self->debug_message("Now combining ref cov stats at ".UR::Context->current->now);
    
    my $refcov_headers_file = "/gscmnt/sata409/research/mmitreva/databases/Bacterial_assemblies.Dec2009.headers_for_refcov.txt";
    
    my $cmd_combine = "perl /gscmnt/sata409/research/mmitreva/sabubuck/HMP_CLINICAL_SAMPLES_JAN_2010/SCRIPTS/combine_refcov_results.pl ".
    					"-refcov $stats_file -db_headers  $refcov_headers_file " .
    					"-ref_counts $readcount_file  -output $combined_file ";
    
    my $rv_combine = Genome::Sys->shellcmd(cmd=>$cmd_combine);
    
    $self->debug_message("Done combining ref cov stats with read counts at ".UR::Context->current->now);
    
    if ($rv_combine) {
    	Genome::Sys->mark_files_ok(input_files=>\@expected_output_files);
    }
    
    $self->debug_message("<<<Completed RefCov at ".UR::Context->current->now);
    
    return 1;
}
1;
