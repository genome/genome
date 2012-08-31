# boberkfe: seems like the bam merge logic here could be folded in w/ DedupSam
# boberkfe: to eliminate some duplicated code
#

package Genome::Model::Event::Build::ReferenceAlignment::DeduplicateLibraries::Maq;

use strict;
use warnings;

use Genome;
use File::Basename;
use File::Copy;
use IO::File;
use Data::Dumper;
use File::stat;

class Genome::Model::Event::Build::ReferenceAlignment::DeduplicateLibraries::Maq {
    is => ['Genome::Model::Event::Build::ReferenceAlignment::DeduplicateLibraries'],
};

sub execute {
    my $self = shift;
    my $now = UR::Time->now;
  
    $self->status_message("Starting DeduplicateLibraries::Maq");
    my $alignments_dir = $self->resolve_accumulated_alignments_path;
    $self->status_message("Accumulated alignments directory: ".$alignments_dir);
   
    unless (-e $alignments_dir) { 
       $self->error_message("Alignments dir didn't get allocated/created, can't continue '$alignments_dir':  $!");
       return;
    }

    #get the instrument data
    my $build = $self->build;
    my @instrument_data = $build->instrument_data;
    my %library_alignments;
    my @all_alignments;

    #accumulate the maps per library
    for my $instrument_data (@instrument_data) {
        my $library = $instrument_data->library_name;
        #for RT#51519 library name contains white space, those white-space lib names should be fixed upstream when putting into DB.
        if ($library =~ /\s+/) {
            $self->warning_message("Library name: $library contains space. will replace with _. Note the lib_tag in bam file will be changed too");
            $library =~ s/\s+/\_/g;
        }

        my @alignments = $build->alignment_results_for_instrument_data($instrument_data);
        for my $alignment (@alignments) {
            my @maps = $alignment->alignment_file_paths;
            push @{$library_alignments{$library}}, @maps;  #for the dedup step
            push @all_alignments, @maps;                   #for the whole genome map file
        }
    }
    
    $self->status_message("Starting dedup workflow with params:");
    #prepare the input for parallelization
    my @list_of_library_alignments;
    for my $library_key ( keys %library_alignments ) {
	    my @read_set_list = @{$library_alignments{$library_key}};	
        $self->status_message("Library: ".$library_key." Read sets count: ". scalar(@read_set_list) ."\n");
        if (scalar(@read_set_list)>0) {
            my %library_alignments_item = ( $library_key => \@read_set_list );  
            push @list_of_library_alignments, \%library_alignments_item;
        } 
        else {
            $self->status_message("Not including library: $library_key because it is empty.");
        } 
    }
    $self->status_message("Size of library alignments: ".@list_of_library_alignments ); 

    if (scalar(@list_of_library_alignments)==0) {
        $self->status_message("None of the libraries contain data.  Quitting.");
        return; 
    }

    #parallelization starts here
    require Workflow::Simple;
        
    my $op = Workflow::Operation->create(
        name => 'Deduplicate libraries.',
        operation_type => Workflow::OperationType::Command->get('Genome::Model::Event::Build::ReferenceAlignment::DeduplicateLibraries::Dedup')
    );
    $op->parallel_by('library_alignments');

    # db disconnect prior to long operation
    if (Genome::DataSource::GMSchema->has_default_handle) {
        $self->status_message("Disconnecting GMSchema default handle.");
        Genome::DataSource::GMSchema->disconnect_default_dbh();
    }
    
    my %extra_params = ();
    $extra_params{dedup_version} = $self->model->duplication_handler_version if $self->model->duplication_handler_version;
    $extra_params{dedup_params}  = $self->model->duplication_handler_params  if $self->model->duplicaiton_handler_params;

    my $output = Workflow::Simple::run_workflow_lsf(
        $op,
        'ref_list'  => $self->model->reference_sequence_build->full_consensus_sam_index_path($self->model->samtools_version),
        'accumulated_alignments_dir' => $alignments_dir, 
        'library_alignments' => \@list_of_library_alignments,
        'aligner_version' => $self->model->read_aligner_version,
        %extra_params,
    );

    #check workflow for errors 
    if (!defined $output) {
       foreach my $error (@Workflow::Simple::ERROR) {
           $self->error_message($error->error);
       }
       return;
    }
    else {
       my $results = $output->{result};
       my $result_libraries = $output->{library_name};
       for (my $i = 0; $i < scalar(@$results); $i++) {
           my $rv = $results->[$i];
            if ($rv != 1) {
                $self->error_message("Workflow had an error while rmdup'ing library: ". $result_libraries->[$i]); 
                die "Workflow had an error while rmdup'ing library: ". $result_libraries->[$i];
            }
        }
    }
 
   #merge those Bam files...BAM!!!
   $now = UR::Time->now;
   $self->status_message(">>> Beginning Bam merge at $now.");
   my $sam_path = Genome::Model::Tools::Sam->path_for_samtools_version($self->model->merger_version);
   my $bam_merge_tool = $sam_path.' merge';
   my $bam_index_tool = $sam_path.' index';
   my $bam_merged_output_file = $self->build->whole_rmdup_bam_file; 
   my $bam_final;
 
   if (-s $bam_merged_output_file )  {
        $self->error_message("The bam file: $bam_merged_output_file already exists.  Skipping bam processing.  Please remove this file and rerun to generate new bam files.");
   }  
   else {
       #get the bam files from the alignments directory
       my @bam_files = <$alignments_dir/*.bam>;

       #remove previously merged/rmdup bam files from the list of files to merge... 
       my $i=0;
       for my $each_bam (@bam_files) {
            #if the bam file name contains the string '_rmdup.bam', remove it from the list of files to merge
            my $substring_index = index($each_bam, "_rmdup.bam");
            unless ($substring_index == -1) {
                    $self->status_message($bam_files[$i]. " will not be merged.");
                    delete $bam_files[$i];
            }		
            $i++;
       }

       if (scalar(@bam_files) == 0 ) {
            $self->error_message("No bam files have been found at: $alignments_dir");
       } 
       elsif (scalar(@bam_files) == 1) {
            my $single_file = shift(@bam_files);
            $self->status_message("Only one bam file has been found at: $alignments_dir. Not merging, only renaming.");
            #my $rename_cmd = "mv ".$single_file." ".$bam_non_merged_output_file;
            $self->status_message("Renaming Bam file from $single_file to $bam_merged_output_file");
            my $bam_rename_rv = move($single_file,$bam_merged_output_file); 
            unless ($bam_rename_rv==1) {
                    $self->error_message("Bam file rename error!  Return value: $bam_rename_rv");
            } 
            else {
                    #renaming success
                    $bam_final = $bam_merged_output_file; 
            } 
       } 
       else {
            $self->status_message("Multiple Bam files found.  Bam files to merge: ".join(",",@bam_files) );
            my $merger_params = $self->model->merger_params;
            $bam_merge_tool .= " $merger_params" if $merger_params;
            my $bam_merge_cmd = "$bam_merge_tool $bam_merged_output_file ".join(" ",@bam_files); 
            $self->status_message("Bam merge command: $bam_merge_cmd");

            #my $bam_merge_rv = system($bam_merge_cmd);
            my $bam_merge_rv = Genome::Sys->shellcmd(cmd=>$bam_merge_cmd,
                                                                     input_files=>\@bam_files,
                                                                     output_files=>[$bam_merged_output_file],
                                                                    );
            $self->status_message("Bam merge return value: $bam_merge_rv");
            unless ($bam_merge_rv == 1) {
                    $self->error_message("Bam merge error!  Return value: $bam_merge_rv");
            } 
            else {
                    #merging success
                    $bam_final = $bam_merged_output_file;
            }
       }

       my $bam_index_rv;
       if (defined $bam_final) {
            $self->status_message("Indexing bam file: $bam_final");
            my $bam_index_cmd = $bam_index_tool ." ". $bam_final;
            #$bam_index_rv = system($bam_index_cmd);
            $bam_index_rv = Genome::Sys->shellcmd(cmd=>$bam_index_cmd,
                                                                  input_files=>[$bam_final],
                                                                  output_files=>[$bam_final.".bai"],
                                                                 );
            unless ($bam_index_rv == 1) {
                    $self->error_message("Bam index error!  Return value: $bam_index_rv");
            } 
            else {
                    #indexing success
                    $self->status_message("Bam indexed successfully.");
            }
       }  
       else {
            #no final file defined, something went wrong	
            $self->error_message("Bam index error!  Return value: $bam_index_rv");
       }

       $now = UR::Time->now;
       $self->status_message("<<< Completing Bam merge at $now.");

       #remove intermediate files
       $now = UR::Time->now;
       $self->status_message(">>> Removing intermediate files at $now");
       
       #remove bam files 
       for my $each_bam_file (@bam_files) {
            $self->status_message("Executing unlink command on $each_bam_file and $each_bam_file.bai");
            my $rm_rv1 = unlink($each_bam_file);
            my $rm_rv2 = unlink($each_bam_file.".bai"); #remove each index as well
            unless ($rm_rv1 == 1) {
                    $self->error_message("There was a problem with the bam remove command: $rm_rv1");
            }  
            unless ($rm_rv2 == 1) {
                    $self->error_message("There was a problem with the bam index remove command: $rm_rv2");
            }
       } 

   } #end else for Bam merge process

   $now = UR::Time->now;
   $self->status_message("<<< Completed removing intermediate files at $now");

   #starting map merge of all library maps 
   $now = UR::Time->now;
   $self->status_message(">>> Beginning mapmerge at $now .");
   my $out_filepath = $alignments_dir;

   my @libraries =  keys %library_alignments; 
   $self->status_message("Libraries: ".join(",",@libraries));
   my @maps_to_merge;
   my $cmd;
   for my $library (@libraries) {
       my $library_file = $alignments_dir .'/'. $library.'.map';
       if (-e $library_file) {
           push @maps_to_merge, $library_file;
       }
   }

   if (@maps_to_merge) {
       $now = UR::Time->now;
       my $maq_pathname = Genome::Model::Tools::Maq->path_for_maq_version($self->model->read_aligner_version);                                
       $cmd ="$maq_pathname mapmerge ". $self->build->whole_rmdup_map_file ." ".join(" ",@maps_to_merge);
   }

   $self->status_message("Executing $cmd at $now.");
   #my $rv = system($cmd);
   my $rv = Genome::Sys->shellcmd(  cmd=>$cmd,
                                                    output_files=>[$self->build->whole_rmdup_map_file],
                                                    input_files=>\@maps_to_merge,
                                                  );
   unless ($rv == 1) {
       $self->error_message("Unexpected return value($rv) from command: $cmd");
       die($self->error_message);
   }

   $self->create_bam_md5;

    for my $file (grep {-f $_} glob($build->accumulated_alignments_directory . "/*")) {
        $self->status_message("Setting $file to read-only");
        chmod 0444, $file;
    }

    $now = UR::Time->now;
    $self->status_message("<<< Completed mapmerge at $now .");
    $self->status_message("*** All processes completed. ***");

    return $self->verify_successful_completion();
}


sub verify_successful_completion {
    my $self = shift;

    my $return_value = 1;
    my $build = $self->build;
            
    unless (-e $build->whole_rmdup_map_file) {
	    $self->error_message("Can't verify successful completeion of Deduplication step. ".$build->whole_rmdup_map_file." does not exist!");	  	
	    return 0;
    } 
    return $return_value;
}

sub calculate_required_disk_allocation_kb {
    my $self = shift;

    $self->status_message("calculating how many map files will get incorporated");
    my $build = $self->build;
    my @instrument_data = $build->instrument_data;
    my @build_maps;
    for my $instrument_data (@instrument_data) {
        my @alignments = $build->alignment_results_for_instrument_data($instrument_data);
        for my $alignment (@alignments) {
            my @aln_maps = $alignment->alignment_file_paths;
            push @build_maps, @aln_maps;
        }
    }
    my $total_size;
    
    print Dumper(\@build_maps);

    for (@build_maps) {
        $total_size += stat($_)->size;
    }

    #take the total size plus a 5 gb safety margin
    # the total size accounts for 4x the map size, including all intermediary files
    ###stage 1
    # - deduped per library map files 1x
    # - map -> sam conversion (allow 5x map size for uncompressed sam file) 5x
    # - sam -> bam conversion (allow 1.25 map size for bam) 1.25x
    ## total burst: 7.25x initial map
    ## remaining:
    # - deduped per library map files 1x
    # - per library map -> bam conversion (1.25 map size) 1.25x
    ## remaining: 2.25x
    
    # stage 2
    # merge all map files:
    # adds 1x total bam size
    ## remaining: 3.25x
    
    $total_size = sprintf("%.0f", ($total_size/1024) * 7.25); 
    return $total_size;
}


1;
