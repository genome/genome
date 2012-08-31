#:boberkfe maybe it's not the totally best place that this event holds the
#:boberkfe disk allocation for the accumulated alignments path.
#:boberkfe but until this becomes a software result it seems like it's the
#:boberkfe best place i could come up with

package Genome::Model::Event::Build::ReferenceAlignment::DeduplicateLibraries;

use strict;
use warnings;

use Genome;
use Command; 

class Genome::Model::Event::Build::ReferenceAlignment::DeduplicateLibraries {
    is => ['Genome::Model::Event'],
    has => [
         disk_allocation     => {
                                calculate_from => [ 'class', 'id' ],
                                calculate => q|
                                    my $disk_allocation = Genome::Disk::Allocation->get(
                                                          owner_class_name => $class,
                                                          owner_id => $id,
                                                      );
                                    return $disk_allocation;
                                |,
        },
    ]
};

sub command_subclassing_model_property {
    return 'duplication_handler_name';
}

sub is_not_to_be_run_by_add_reads {
    return 1;
}

sub resolve_accumulated_alignments_path {
    my $self = shift;

    my $build_accumulated_alignments_path = $self->build->accumulated_alignments_directory;

    if (-d $build_accumulated_alignments_path || -l $build_accumulated_alignments_path ) {
        print "$build_accumulated_alignments_path already exists as a directory or symlink, not creating an allocation";
        return $build_accumulated_alignments_path;
    }

    my $kb_needed = $self->calculate_required_disk_allocation_kb; 

    my $allocation_path = sprintf("build_merged_alignments/build%s",$self->build->id);

    my $allocation = $self->disk_allocation;

    unless ($allocation) {

        $allocation = Genome::Disk::Allocation->allocate(
                                                                  disk_group_name => 'info_genome_models',
                                                                  allocation_path => $allocation_path,
                                                                  kilobytes_requested => $kb_needed,
                                                                  owner_class_name => $self->class,
                                                                  owner_id => $self->id,
                                                                  );
        unless ($allocation) {
             $self->error_message('Failed to get disk allocation for accumulated alignments.  This dedup event needed $kb_needed kb in order to run.');
             die $self->error_message;
        }
    }

    unless (symlink($allocation->absolute_path, $build_accumulated_alignments_path)) {
            $self->error_message("Failed to symlink " . $allocation->absolute_path . " to the accumulated alignments path in the build dir."); 
            die;
    }

    return $build_accumulated_alignments_path;

}

sub calculate_required_alignment_disk_allocation_kb {
    die "abstract, must be defined in your dedup subclass!";
}
  

sub create_bam_md5 {
    my $self = shift;

    my $bam_merged_output_file = $self->build->whole_rmdup_bam_file;
    my $md5_file = $bam_merged_output_file.'.md5';
    my $cmd = "md5sum $bam_merged_output_file > $md5_file";

    $self->status_message("Creating md5 file for the whole rmdup BAM file...");

    my $md5_rv  = Genome::Sys->shellcmd(
        cmd                        => $cmd, 
        input_files                => [$bam_merged_output_file],
        #output_files               => [$md5_file],
        skip_if_output_is_present  => 0,
    ); 

    $self->warning_message("Failed to create bam md5 for $bam_merged_output_file") unless $md5_rv == 1;

    return 1;
}


1;

