package Genome::Model::Tools::HmpShotgun::MergeAlignments;

use strict;
use warnings;

use Genome;
use Workflow;
use IO::File;

class Genome::Model::Tools::HmpShotgun::MergeAlignments {
    is  => ['Command'],
    has => [
        working_directory => {
            is  => 'ARRAY',
            is_input => '1',
            doc => 'The working directory.',
        },
        alignment_files => {
        	is  => 'ARRAY',
            is_input => '1',
            doc => 'The reads to align.',
        },
        unaligned_files => {
        	is  => 'ARRAY',
            is_input => '1',
            doc => 'The unaligned files.',
        },
        reference1_aligned_file => {
        	is  => 'String',
            is_output => '1',
            is_optional => '1',
            doc => 'The resulting alignment.',
        },
        reference2_aligned_file => {
        	is  => 'String',
            is_output => '1',
            is_optional => '1',
            doc => 'The resulting alignment.',
        },
         _working_directory => {
        	is  => 'String',
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
    $self->status_message(">>>Running MergeAlignments at ".UR::Time->now);
    #my $model_id = $self->model_id;
    my $alignment_files_ref = $self->alignment_files;
    my @alignment_files = @$alignment_files_ref;
    
    my $unaligned_files_ref = $self->unaligned_files;
    my @unaligned_files = @$unaligned_files_ref;
    
    #get parallelized inputs 
    #all of the alignment jobs are sending in the same working directory.  
    #pick the first one
    my $working_directory_ref = $self->working_directory;
    my @working_directory_list = @$working_directory_ref;
    my $working_directory = $working_directory_list[0];
    $self->_working_directory($working_directory);
    
    #my $working_directory = $self->working_directory."/alignments/";
    $self->status_message("Working directory: ".$working_directory);
    
    #first cat the unaligned reads
    my $unaligned_combined = $working_directory."/unaligned_combined.sam";
    my $rv_unaligned_check = Genome::Sys->is_file_ok($unaligned_combined);
    if ($rv_unaligned_check) {
    	$self->status_message("Output file exists.  Skipping the generation of the unaligned reads files.  If you would like to regenerate these files, remove them and rerun.");  	
    } else {
    	my $rv_cat = Genome::Sys->cat(input_files=>\@unaligned_files,output_file=>$unaligned_combined);
    	if ($rv_cat) {
    		Genome::Sys->mark_file_ok($unaligned_combined);
    	} else {
    		$self->status_message("There was a problem generating the combined unaligned file: $unaligned_combined");
    		#may want to return here.
    	}
    }
    
    
    #This is super hacky.  Need the name of where the bam came from 
    #probably pass this in from the workflow.
    #Right now, it is grep-ing the path of the input files to determine 
    #which reference sequence is which.
    
    my $refseq1 = "metagenome1";
    my $refseq2 = "metagenome2";
    
    my @refseq1 = grep(/^.*$refseq1.*/, @alignment_files); 
    my @refseq2 = grep(/^.*$refseq2.*/, @alignment_files); 
    
    my @refseq_list = ( \@refseq1, \@refseq2 );
    my @refseq_name_list = ( $refseq1, $refseq2 );
    
    my $refseq1_name_sorted_bam = $self->resolve_name_sorted_file_name($refseq1,"bam");  
    my $refseq1_name_sorted_sam = $self->resolve_name_sorted_file_name($refseq1,"sam");
    my $refseq2_name_sorted_bam = $self->resolve_name_sorted_file_name($refseq2,"bam");
    my $refseq2_name_sorted_sam = $self->resolve_name_sorted_file_name($refseq2,"sam");
    #my $refseq2_name_sorted_sam = $working_directory."/".$refseq2."_name_sorted.sam"; 
    
    my @expected_output_files = ( $refseq1_name_sorted_bam, $refseq1_name_sorted_sam, $refseq2_name_sorted_bam, $refseq2_name_sorted_sam);
    my $rv_check = Genome::Sys->are_files_ok(input_files=>\@expected_output_files);
    
    if (defined($rv_check)) {
	    if ($rv_check == 1) {
	    	#shortcut this step, all the required files exist.  Quit.
	    	$self->reference1_aligned_file($refseq1_name_sorted_sam);
            $self->reference2_aligned_file($refseq2_name_sorted_sam);
	    	$self->status_message("Skipping this step.  If you would like to regenerate these files, remove them and rerun.");
	   	    $self->status_message("<<<Completed MergeAlignments at ".UR::Time->now);
	   	    return 1;
	    }
    }
	
    if (scalar(@alignment_files) < 2) {
             $self->error_message("*** Invalid number of files to merge: ".scalar(@alignment_files).". Must have 2 or more.  Quitting.");
             return;
    } else {
            for my $refseq_list_item (@refseq_list) {
                my @alignment_files_per_refseq = @$refseq_list_item;
                my $refseq_name = shift(@refseq_name_list);
                my $merged_alignment_files_per_refseq = $working_directory."/".$refseq_name."_coord_sorted.bam";
                my $merged_alignment_files_per_refseq_name_sorted =   $self->resolve_name_sorted_file_name($refseq_name,"bam");
                my $merged_alignment_files_per_refseq_sam =  $self->resolve_name_sorted_file_name($refseq_name,"sam");
                
                $self->status_message("Merging files: ".join("\n",@alignment_files_per_refseq) );
                $self->status_message("Destination file: ".$merged_alignment_files_per_refseq);
        
                #get from pp eventually
                my $picard_path = "/gsc/scripts/lib/java/samtools/picard-tools-1.07/";
               
                my $merge_tool = "java -Xmx2g -cp $picard_path/MergeSamFiles.jar net.sf.picard.sam.MergeSamFiles MSD=true SO=coordinate AS=true tmp_dir=$working_directory VALIDATION_STRINGENCY=SILENT O=$merged_alignment_files_per_refseq ";
                my $list_of_files = join(' I=',@alignment_files_per_refseq);                              
                my $cmd_merge = $merge_tool." I=".$list_of_files;		                   
                my $rv_merge = Genome::Sys->shellcmd(cmd=>$cmd_merge);											 
            
                if ($rv_merge != 1) {
                        $self->error_message("<<<Failed MergeAlignments on picard merge.  Return value: $rv_merge");
                        return;
                }
                $self->status_message("Merge complete.");
                
                #name sorting
			    my $sorter = Genome::Model::Tools::Sam::SortBam->create(file_name=>$merged_alignment_files_per_refseq,
										    							name_sort=>1,
			    														output_file=>$merged_alignment_files_per_refseq_name_sorted);
			   	my $rv_sort = $sorter->execute;
			   	if ($rv_sort != 1) {
			   		$self->error_message("Sort failed.  Return value: $rv_sort");
			    	return;
			   	}
                
                ####end sort and move                
                
                $self->status_message("Converting from bam to sam file: $merged_alignment_files_per_refseq_name_sorted to $merged_alignment_files_per_refseq_sam");
                #my $cmd_convert = "java -Xmx2g -cp $picard_path/SamFormatConverter.jar net.sf.picard.sam.SamFormatConverter VALIDATION_STRINGENCY=SILENT I=$merged_alignment_files_per_refseq O=$merged_alignment_files_per_refseq_sam";  
                my $cmd_convert = "samtools view $merged_alignment_files_per_refseq_name_sorted > $merged_alignment_files_per_refseq_sam";
                my $rv_convert = Genome::Sys->shellcmd(cmd=>$cmd_convert);											 
            
                if ($rv_convert != 1) {
                        $self->error_message("<<<Failed MergeAlignments on bam to sam conversion.  Return value: $rv_merge");
                        return;
                }
                    
            }#end for loop
            
        
    
	}
    
    Genome::Sys->mark_files_ok(input_files=>\@expected_output_files);
    $self->reference1_aligned_file($refseq1_name_sorted_sam);
    $self->reference2_aligned_file($refseq2_name_sorted_sam);
    $self->status_message("<<<Completed MergeAlignments for testing at at ".UR::Time->now);
    return 1;
 
}

sub resolve_name_sorted_file_name {
	my $self = shift;
	my $refseq_name = shift;
	my $extension = shift;
	return $self->_working_directory."/".$refseq_name."_name_sorted.".$extension;
}

1;
