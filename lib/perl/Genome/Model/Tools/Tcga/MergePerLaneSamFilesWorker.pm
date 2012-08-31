package Genome::Model::Tools::Tcga::MergePerLaneSamFilesWorker;

use strict;
use warnings;

use Genome;
use Workflow;
use IO::File;

class Genome::Model::Tools::Tcga::MergePerLaneSamFilesWorker {
    is  => ['Command'],
    has => [
        seq_ids => {
            is  => 'String',
            is_input => 1,
            doc => 'The directory containing the Maq map files.',
        },
        working_directory => {
            is => 'String',
            is_input => 1,
            doc => 'The working directory of the tool.',
        },
        seq_dict_sam_file => {
            is => 'String',
            is_input => 1,
            doc => 'The sequence dictionary for the reference.',
        },
        ref_list => {
            is => 'String',
            is_input => 1,
            doc => 'The ref_list used for SamToBam.',
            default => Genome::Config::reference_sequence_directory() . '/NCBI-human-build36/all_sequences.fa.fai',
            is_optional => 1,
        },
    ],
    has_param => [
           lsf_resource => {
           default_value => 'select[model!=Opteron250 && type==LINUX64] rusage[mem=2000]',
           },
    ],

};

sub help_brief {
    'Convert Maq map files into the TCGA format.';
}

sub help_detail {
    return <<EOS
    Convert Maq map files into the TCGA format.
EOS
}

sub execute {
    my $self = shift;
    my $working_dir = $self->working_directory;
    
    my $seq_id = $self->seq_ids; 
    
    my $pid = getppid();
    my $log_path = "$working_dir/logs/merge_per_lane_$seq_id"."_$pid.txt";
    my $log_fh = Genome::Sys->open_file_for_writing($log_path);
    unless($log_fh) {
       $self->error_message("Failed to open output filehandle for: " .  $log_path );
       die "Could not open file ".$log_path." for writing.";
    } 

    print $log_fh "Merging sams for seq id: $seq_id\n";

    my $sam_directory = $self->working_directory."/sams/";
    my $bam_directory = $self->working_directory."/bams/";

    my $seqdict = $self->seq_dict_sam_file; 
    my $rg_file = $self->working_directory."/header/$seq_id.rg"; 
    my $pg_file = $self->working_directory."/header/$seq_id.pg"; 
    my $aligned_file = $self->working_directory."/aligned/$seq_id.sam"; 
    my $unaligned_file = $self->working_directory."/unaligned/$seq_id.sam"; 

    if (!-s $unaligned_file) {
        print $log_fh "$unaligned_file not found for $seq_id.  Returning.\n";
        return 1;
    }

    if (!-s $aligned_file) {
        print $log_fh "$aligned_file not found for $seq_id.  Returning.\n";
        return 1;
    }

    my @files_to_merge;

    my $per_lane_sam_file = "$sam_directory/$seq_id.sam";
   
    if (!-s $per_lane_sam_file) {
        push(@files_to_merge, $seqdict, $rg_file, $pg_file, $aligned_file, $unaligned_file); 
        print $log_fh "\nCat-ing together: ".join("\n",@files_to_merge);
        my $cat_rv = Genome::Sys->cat(input_files=>\@files_to_merge,output_file=>$per_lane_sam_file);
        if ($cat_rv ne 1) {
            print $log_fh "\nFor $seq_id, error during cat! Return value $cat_rv";
            die "For $seq_id, error cat-ing all files together.  Return value: $cat_rv";
        } else {
            print $log_fh "\nCat successful.";
        } 
    } else {
        print log_fh "\nThe per lane sam file: $per_lane_sam_file already exists.  Skipping the generation of this file.";
    }

    #do bams
    my $per_lane_bam_file = $bam_directory."/$seq_id.bam";
    
    print $log_fh "\nConverting sam $per_lane_sam_file to bam $per_lane_bam_file";

    if (!-s $per_lane_bam_file) { 
        my $to_bam = Genome::Model::Tools::Sam::SamToBam->create(
            bam_file => $per_lane_bam_file, 
            sam_file => $per_lane_sam_file,                                                      
            ref_list => $self->ref_list,
            keep_sam => 1,
            fix_mate => 1,
            index_bam => 0,
        );
        my $rv_to_bam = $to_bam->execute();
        if ($rv_to_bam ne 1) { 
            $self->error_message("There was an error converting the Sam file $per_lane_sam_file to $per_lane_bam_file.  Return value was: $rv_to_bam");
            die "There was an error converting the Sam file $per_lane_sam_file to $per_lane_bam_file.  Return value was: $rv_to_bam";
        } else {
            print $log_fh "\nConversion successful.";
        }
    } else {
        print log_fh "\nThe fixed mate bam file: $per_lane_bam_file already exists.  Skipping the generation of this file.";
    }

    #cleanup 
    #my @tmp_files = <$per_lane_bam_file_directory/*>;
    #print log_fh "Directory contents: ".join("\n",@tmp_files);
    #for my $tmp_file (@tmp_files) {
    #    if ($tmp_file != $per_lane_bam_file) {
    #        unlink($tmp_file);
    #        print log_fh "Unlinked $tmp_file";
    #    }
    #}
  
    $log_fh->close;
    return 1; 

}
1;
