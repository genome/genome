package Genome::Model::Tools::SimpleAlignment::MergeAlignments;

use strict;
use warnings;

use Genome;
use Workflow;
use IO::File;

class Genome::Model::Tools::SimpleAlignment::MergeAlignments {
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
        merged_file => {
            is  => 'String',
            is_input => '1',
            is_output => '1',
            doc => 'The resulting alignment.',
        },
        lsf_resource => {
            is_param  => 1,
            value => "-R 'select[mem>4000 model!=Opteron250 && type==LINUX64] rusage[mem=4000]' -M 4000000",
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
    
    #get parallelized inputs 
    my $working_directory = $self->working_directory;
    my $merged_file = $self->merged_file;
    
    if (scalar(@alignment_files) == 0) {
             $self->error_message("*** Invalid number of files to merge: ".scalar(@alignment_files).". Must have 1 or more.  Quitting.");
             return;
    } else {
                
                $self->status_message("Merging files: ".join("\n",@alignment_files) );
                $self->status_message("Destination file: ".$merged_file);
         
                #get from pp eventually
                #my $picard_path = "/gsc/scripts/lib/java/samtools/picard-tools-1.07/";
                #my $merge_tool = "java -Xmx3g -cp $picard_path/MergeSamFiles.jar net.sf.picard.sam.MergeSamFiles MSD=true SO=coordinate AS=true tmp_dir=$working_directory VALIDATION_STRINGENCY=SILENT O=$merged_file ";
                #my $list_of_files = join(' I=',@alignment_files);
                #my $cmd_merge = $merge_tool." I=".$list_of_files;		
                #my $rv_merge = Genome::Sys->shellcmd(cmd=>$cmd_merge);											 
           
                my $cmd_merge = Genome::Model::Tools::Sam::Merge->create(
                    files_to_merge => \@alignment_files,
                    merged_file    => $merged_file,
                    bam_index      => 1,
                    is_sorted      => 1,
                );

                my $rv_merge = $cmd_merge->execute;
 
                if ($rv_merge != 1) {
                        $self->error_message("<<<Failed MergeAlignments on picard merge.  Return value: $rv_merge");
                        return;
                }
                $self->status_message("Merge complete.");
                
    }
    
    #Genome::Sys->mark_files_ok(input_files=>\@expected_output_files);
    $self->status_message("<<<Completed MergeAlignments at ".UR::Time->now);
    return 1;
 
}
1;
