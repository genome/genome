package Genome::Model::Event::Build::ReferenceAlignment::DeduplicateLibraries::Samtools;

#REVIEW fdu 11/18/2009
#This module reflects old way to do dedup, which is running "samtools
#rmdup" on each merged per library bam by workflow, and then merging them 
#together to a final merged_rmdup.bam. The new way (see Picard.pm) is
#to run picard Markduplicate on the whole bam mixed with all libraries. 
#Once the new way is proved ok by quite a few of model builds, this
#module along with DedupSam.pm will soon become obsolete and be removable.

use strict;
use warnings;

use Genome;
use File::Basename;
use File::Copy;
use IO::File;
use File::stat;
use Data::Dumper;

class Genome::Model::Event::Build::ReferenceAlignment::DeduplicateLibraries::Samtools {
    is => ['Genome::Model::Event::Build::ReferenceAlignment::DeduplicateLibraries'],
};

sub execute {
    my $self = shift;
    my $now = UR::Time->now;
  
    $self->status_message("Starting DeduplicateLibraries::Samtools");
    my $alignments_dir = $self->resolve_accumulated_alignments_path;
    $self->status_message("Accumulated alignments directory: ".$alignments_dir);
   
    unless (-e $alignments_dir) { 
       $self->error_message("Alignments dir didn't get allocated/created, can't continue '$alignments_dir':  $!");
       return;
    }

    #my $bam_merged_output_file = $alignments_dir."/".$self->model->subject_name."_merged_rmdup.bam";
    my $bam_merged_output_file = $self->build->whole_rmdup_bam_file; 
    if (-e $bam_merged_output_file) {
        $self->status_message("A merged and rmdup'd bam file has been found at: $bam_merged_output_file");
        $self->status_message("If you would like to regenerate this file, please delete it and rerun.");
        $now = UR::Time->now;
        $self->status_message("Skipping the rest of DeduplicateLibraries::Samtools at $now");
        $self->status_message("*** All processes skipped. ***");
        return 1;
    } 

    #get the instrument data
    my $build = $self->build;
    my @instrument_data = $build->instrument_data;
    my %library_alignments;
    my @all_alignments;

    #accumulate the maps per library
    for my $instrument_data (@instrument_data) {
        my $library = $instrument_data->library_name;
        #for RT#51519 library name contains white space, those white-space lib names need be fixed upstream when putting into DB.
        if ($library =~ /\s+/) {
            $self->warning_message("Library name: $library contains space. will replace with _.");
            $library =~ s/\s+/\_/g;
        }
        my @alignments = $build->alignment_results_for_instrument_data($instrument_data);
        for my $alignment (@alignments) {
            my @bams = $alignment->alignment_bam_file_paths;
            $self->status_message("bam file paths: ". @bams);

            push @{$library_alignments{$library}}, @bams;  #for the dedup step
            push @all_alignments, @bams;                   #for the whole genome map file
        }
    }
    
    $self->status_message("Starting SAM dedup workflow with params:");
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

    my $merger_name    = $self->model->merger_name || 'picard';
    my $merger_version = $self->model->merger_version;
    my $merger_params  = $self->model->merger_params;

    #parallelization starts here
    require Workflow::Simple;
        
    my $op = Workflow::Operation->create(
            name => 'Deduplicate libraries.',
            operation_type => Workflow::OperationType::Command->get('Genome::Model::Event::Build::ReferenceAlignment::DeduplicateLibraries::DedupSam')
    );

    $op->parallel_by('library_alignments');

    # db disconnect prior to long operation
    if (Genome::DataSource::GMSchema->has_default_handle) {
        $self->status_message("Disconnecting GMSchema default handle.");
        Genome::DataSource::GMSchema->disconnect_default_dbh();
    }

    my $output = Workflow::Simple::run_workflow_lsf(
            $op,
            'accumulated_alignments_dir' => $alignments_dir, 
            'library_alignments' => \@list_of_library_alignments,
            'dedup_version' => $self->model->duplication_handler_version,
            'dedup_params'  => $self->model->duplication_handler_params,
            'merger_name' => $merger_name,
            'merger_version' => $merger_version,
            'merger_params'  => $merger_params,
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
   
   #remove original library input files
   my @original_to_remove_files = grep {$_ !~ m/rmdup/ }<$alignments_dir/*.bam>;
   for (@original_to_remove_files) {
        $self->status_message("Removing intermediate library file $_");
        unlink($_);
   }
 
   #merge those Bam files...BAM!!!
   $now = UR::Time->now;
   $self->status_message(">>> Beginning Bam merge at $now.");
   #my $bam_merged_output_file = $alignments_dir."/".$self->model->subject_name."_merged_rmdup.bam";

   my @bam_files = <$alignments_dir/*_merged_rmdup.bam>;
   
    #remove previously merged/rmdup bam files from the list of files to merge... 
    #    my $i=0;
    #   for my $each_bam (@bam_files) {
    #        #if the bam file name contains the string '_rmdup.bam', remove it from the list of files to merge
    #        my $substring_index = index($each_bam, "_rmdup.bam");
    #        unless ($substring_index == -1) {
    #                $self->status_message($bam_files[$i]. " will not be merged.");
   #                delete $bam_files[$i];
   #         }		
   #         $i++;
   #    }

   # these are already sorted coming out of the initial merge, so don't bother re-sorting
   my $merge_rv = Genome::Model::Tools::Sam::Merge->execute(
       files_to_merge => \@bam_files,
       merged_file => $bam_merged_output_file,
       is_sorted => 1,
       #software => $merge_software,
       merger_name => $merger_name,
       merger_version => $merger_version,
       merger_params => $merger_params
   ); 

   $now = UR::Time->now;
   $self->status_message("<<< Completing Bam merge at $now.");

   #remove intermediate files
   $now = UR::Time->now;
   $self->status_message(">>> Removing intermediate files at $now");

    #delete everything except big dedup bam file and index 
    my @all_files = <$alignments_dir/*>;
    for my $each_bam_file (@all_files) {
        if ( ($each_bam_file eq $bam_merged_output_file) || ($each_bam_file eq $bam_merged_output_file.".bai" ) ) {   
            $self->status_message("Keeping $each_bam_file");
        } 
        else {
            $self->status_message("Executing unlink command on $each_bam_file");
            my $rm_rv1 = unlink($each_bam_file);
            unless ($rm_rv1 == 1) {
                $self->error_message("There was a problem with the bam remove command: $rm_rv1");
            }
        }
    }

    $self->create_bam_md5;

    for my $file (grep {-f $_} glob($build->accumulated_alignments_directory . "/*")) {
        $self->status_message("Setting $file to read-only");
        chmod 0444, $file;
    }

    $now = UR::Time->now;
    $self->status_message("<<< Completed removing intermediate files at $now");
    $self->status_message("*** All processes completed. ***");

    return $self->verify_successful_completion();
}


sub verify_successful_completion {
    my $self = shift;
    my $return_value = 1;
    my $build = $self->build;
            
    unless (-e $build->whole_rmdup_bam_file) {
	    $self->error_message("Can't verify successful completeion of Deduplication step. ".$build->whole_rmdup_bam_file." does not exist!");	  	
	    return 0;
    } 

    return $return_value;

}

sub calculate_required_disk_allocation_kb {
    my $self = shift;

    $self->status_message("calculating how many bam files will get incorporated...");
    
    my $build = $self->build;
    my @instrument_data = $build->instrument_data;
    my @build_bams;
    for my $instrument_data (@instrument_data) {
        my @alignments = $build->alignment_results_for_instrument_data($instrument_data);
        for my $alignment (@alignments) {
            my @aln_bams = $alignment->alignment_bam_file_paths;
            push @build_bams, @aln_bams;
        }
    }
    my $total_size;
    
    for (@build_bams) {
        $total_size += stat($_)->size;
    }

    #take the total size plus a 10% safety margin
    # 3x total size; individual/deduped per-lib bams, full build deduped bam
    $total_size = sprintf("%.0f", ($total_size/1024)*1.1); 
    $total_size = ($total_size * 2);

    return $total_size;
}


1;
