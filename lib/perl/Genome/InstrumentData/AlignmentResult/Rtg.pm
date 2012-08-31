package Genome::InstrumentData::AlignmentResult::Rtg;

use strict;
use warnings;
use File::Basename;

use Genome;

#  So, you want to build an aligner?  Follow these steps.
#
#  1) set aligner name in the UR class def
#  2) set required_rusage
#  3) Implement run_aligner
#  4) Implement aligner_params_for_sam_header
#
#  You also will want to create a Genome::InstrumentData::Command::Align::YOURCLASSHERE
#  so you can align from the command line too!

class Genome::InstrumentData::AlignmentResult::Rtg {
    is => 'Genome::InstrumentData::AlignmentResult',
    
    # TODO: Put your aligner name here
    has_constant => [
        aligner_name => { value => 'rtg', is_param=>1 },
    ],
};

sub required_arch_os { 'x86_64' }

#TODO: Put the LSF resources required to run the alignment here.
sub required_rusage { 
    "-R 'select[type==LINUX64 && tmp>90000 && mem>30000] span[hosts=1] rusage[tmp=90000, mem=30000]' -M 30000000 -n 4";
}

#
#  Implement this method with the actual logic to run your aligner.
#  The pathnames for input files are passed in.
#
#  The expected output is an "all_sequences.sam" file in the scratch directory.
#  YOU MUST APPEND TO THIS FILE, NOT OVERWRITE.  This is because the calling code
#  may run multiple passes of run_aligner.  For example, when running trimming,
#  there may be two invocations - one with paired data, and one with a set of remaining
#  singleton reads whose mates were clipped too far in to be usable.
#
#  The sam file also needs to NOT have any SAM headers (i.e. lines starting with "@").
#  The pipeline adds its own detailed headers that are appropriate to the data.
#  Having headers here already will cause issues in the pipeline
#  downstream with picard merge.
#

sub _run_aligner {
    my $self = shift;
    my @input_pathnames = @_;


    # get refseq info
    my $reference_build = $self->reference_build;
    
    #TODO: This'll give you the path to a reference index in the reference directory ending in .fa
    #If your flavor of aligner uses a different file extension for its index, put it here.
    my $reference_fasta_path = $reference_build->full_consensus_path('fa');
    
    # This is your scratch directory.  Whatever you put here will be wiped when the alignment
    # job exits.
    my $scratch_directory = $self->temp_scratch_directory;
    
    # This is the alignment output directory.  Whatever you put here will be synced up to the
    # final alignment directory that gets a disk allocation.
    my $staging_directory = $self->temp_staging_directory;
 
    # This is the SAM file you should be appending to.  Dont forget, no headers!
    my $sam_file = $scratch_directory . "/all_sequences.sam";

    
    # TODO: implement your aligner logic here.  If you need to condition on whether the
    # data is single-ended or paired-ended, key off the number of files passed in (1=SE, 2=PE)
    # Under no circumstances should you ever get more than 2 files, if you do then that's bad and
    # you should die.
     
    my $output_dir = $self->temp_scratch_directory . "/rtg_default";
    my $sdf_dir = $self->temp_scratch_directory . "/rtg_sdf";
    my $format_cmd = "/gscmnt/gc2146/info/medseq/rtg_software/rtg-BETA-2011-04-29-34609/rtg format -o $sdf_dir";
     $format_cmd .=" --format=fastq";
     $format_cmd .=" --quality-format=sanger";  

    my $finished_sam; 
    if (@input_pathnames == 1) {
        $self->status_message("_run_aligner called in single-ended mode.");
        my $fastq_name = $input_pathnames[0];
        $format_cmd .= " $fastq_name";
        $finished_sam = $output_dir . "/alignments.sam.gz";
    } elsif (@input_pathnames == 2) {
        $self->status_message("_run_aligner called in paired-end mode.");
        my $left_fastq = $input_pathnames[0];
        my $right_fastq = $input_pathnames[1];
         $format_cmd .= " -l $left_fastq";
         $format_cmd .= " -r $right_fastq";
         $finished_sam = $output_dir . "/mated.sam.gz";
    } else {
        $self->error_message("_run_aligner called with " . scalar @input_pathnames . " files.  It should only get 1 or 2!");
        die $self->error_message;
    }
     my $cmd = "/gscmnt/gc2146/info/medseq/rtg_software/rtg-BETA-2011-04-29-34609/rtg map";
       $cmd .= " -o $output_dir";
       $cmd .=" -t /gscmnt/gc2146/info/medseq/rtg_software/rtg-BETA-2011-04-29-34609/build_36.sdf";
       $cmd .=" --legacy-cigars";
       #   $cmd .=" --read-names";
       $cmd .=" --report-unmapped";
       #   $cmd .=" --format=fastq";
       #  $cmd .=" --quality-format=sanger";
       $cmd .=" -n 1";
       $cmd .=" -T 4";
       $cmd .=" -i $sdf_dir";
  
    my $finished_unmated_sam = $output_dir . "/unmated.sam.gz";
    my $unmated_cmd = "zcat $finished_unmated_sam | grep -v '^\@' >> $sam_file";
    my $finished_unmapped_sam = $output_dir . "/unmapped.sam.gz";
    my $unmapped_cmd = "zcat $finished_unmapped_sam | grep -v '^\@' >> $sam_file"; 
    my $bam_cmd="zcat $finished_sam | grep -v '^\@' >> $sam_file";
    $self->status_message("Trying $format_cmd");
    Genome::Sys->shellcmd(  cmd => $format_cmd, input_files =>[@input_pathnames], output_files => [] );
    $self->status_message("Trying: $cmd");
    Genome::Sys->shellcmd(  cmd => $cmd, input_files =>[], output_files => [$finished_sam] );
    $self->status_message("Trying: $bam_cmd");
    Genome::Sys->shellcmd(  cmd => $bam_cmd, input_files => [$finished_sam], output_files => [$sam_file], skip_if_output_is_present=>0);
     $self->status_message("Trying: $unmapped_cmd");
    Genome::Sys->shellcmd(  cmd => $unmapped_cmd, input_files => [$finished_unmapped_sam], output_files => [$sam_file], skip_if_output_is_present=>0);
    if(@input_pathnames == 2) {
        $self->status_message("Trying: $unmated_cmd");
        Genome::Sys->shellcmd(  cmd => $unmated_cmd, input_files => [$finished_unmated_sam], output_files => [$sam_file], skip_if_output_is_present=>0);
    }
    
    #  `cp $sam_file /gscmnt/sata921/info/medseq/rtg_evaluation/`;
    # `cp $output_dir /gscmnt/sata921/info/medseq/rtg_evaluation/`;
   

    # confirm that at the end we have a nonzero sam file, this is what'll get turned into a bam and
    # copied out.
    unless (-s $sam_file) {
        die "The sam output file $sam_file is zero length; something went wrong.";
    }
    
    # TODO:
    # If you have any log files from the aligner you're wrapping that you'd like to keep,
    # copy them out into the staging directory here.
    
    
    # TODO:
    # Do any last minute checks on the staged data and die if they fail.
    
    
    # If we got to here, everything must be A-OK.  AlignmentResult will take over from here
    # to convert the sam file to a BAM and copy everything out.
    
    return 1;
}

# TODO: This should be the command line used to run your aligner, not just the params.
# sorry for the mis-named method name, it'll get fixed soon.
#
# This will end up in our re-composed sam/bam header.

sub aligner_params_for_sam_header {
    my $self = shift;
    
    return "rtg map" . $self->aligner_params;
}

sub fillmd_for_sam {
    return 0;
}

# If your aligner adds read group tags, or you are handling it in the wrapper, this needs to be 0
# otherwise 1.  If you are streaming output straight to BAM then you need to take care of adding RG
# tags with either the wrapper or the aligner itself, and this needs to be 0.
sub requires_read_group_addition {
    return 1;
}

# if you are streaming to bam, set this to 1.  Beware of read groups.
sub supports_streaming_to_bam {
    return 0;
}

# If your aligner accepts BAM files as inputs, return 1.  You'll get a set of BAM files as input, with
# suffixes to define whether it's paired or unparied.
# [input_file.bam:0] implies SE, [input_file.bam:1, input_file.bam:2] implies paired.
sub accepts_bam_input {
    return 0;
}

# Use this to prep the index.  Indexes are saved for each combo of aligner params & version, as runtime
# params for some aligners also call for specific params when building the index.
sub prepare_reference_sequence_index {
    my $class = shift;
    my $refindex = shift;

    # If you need the parameters the aligner is being run with, in order to customize the index, here they are.
    my $aligner_params = $refindex->aligner_params;

    my $staging_dir = $refindex->temp_staging_directory;
    my $staged_fasta_file = sprintf("%s/all_sequences.fa", $staging_dir);

    my $actual_fasta_file = $staged_fasta_file;

    if (-l $staged_fasta_file) {
        $class->status_message(sprintf("Following symlink for fasta file %s", $staged_fasta_file));
        $actual_fasta_file = readlink($staged_fasta_file);
        unless($actual_fasta_file) {
            $class->error_message("Can't read target of symlink $staged_fasta_file");
            return;
        } 
    }

    # Do whatever it takes here to put the index into $staging_dir.

    # If you return 1 then this reference will be copied up from $staging_dir into a software result for future use.
    return 1;
}
