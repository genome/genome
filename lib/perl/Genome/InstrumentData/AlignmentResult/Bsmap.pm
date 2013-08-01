package Genome::InstrumentData::AlignmentResult::Bsmap;

use strict;
use warnings;
use IO::File;
use File::Basename;
use File::Copy;
use File::Temp;
use Sort::Naturally;
use Genome;

#  So, you want to build an aligner?  Follow these steps.
#
#  1) set aligner name in the UR class def
#  2) set required_rusage
#  3) Implement run_aligner
#  4) Implement aligner_params_for_sam_header

class Genome::InstrumentData::AlignmentResult::Bsmap {
    is => 'Genome::InstrumentData::AlignmentResult',
    has_constant => [
        aligner_name => { value => 'bsmap', is_param=>1 }
    ]
};

sub required_arch_os { 'x86_64' }

# LSF resources required to run the alignment.
#"-R 'select[model!=Opteron250 && type==LINUX64 && mem>16000 ** tmp > 150000] span[hosts=1] rusage[tmp=150000, mem=16000]' -M 16000000 -n 1";
sub required_rusage {
    my $self = shift;
    my $cores = 4; # default to 4, although decomposed_aligner_params should always include "-p N" even if the processing-profile does not

    my $params = $self->decomposed_aligner_params();
    if ($params =~ /-p\s+(\d+)/) {
        $cores = $1;
    }

    return "-R 'select[type==LINUX64 && mem>16000 && tmp > 300000] span[hosts=1] rusage[tmp=300000, mem=16000]' -M 16000000 -n $cores";
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

    # get refseq info and fasta files
    my $reference_build = $self->reference_build;
    my $reference_fasta_path = $reference_build->full_consensus_path('fa');

    # get the index directory
    #my $reference_index = $self->get_reference_sequence_index();
    #my $reference_index_directory = dirname($reference_index->full_consensus_path());

    #my $reference_index_directory = $reference_index->data_directory(); # better way to do this?
    #print "Ref index dir: $reference_index_directory\n";
    # example dir /gscmnt/sata921/info/medseq/cmiller/methylSeq/bratIndex

    # This is your scratch directory.  Whatever you put here will be wiped when the alignment
    # job exits.
    my $scratch_directory = $self->temp_scratch_directory;

    # This is (another) temporary directory. The difference between this and the scratch directory is that
    # this is blown away between calls of _run_aligner when running in force_fragment mode, while
    # the scratch directory will stick around.
    my $temporary_directory = File::Temp->tempdir("_run_aligner_XXXXX", DIR => $scratch_directory);

    # This is the alignment output directory.  Whatever you put here will be synced up to the
    # final alignment directory that gets a disk allocation.
    my $staging_directory = $self->temp_staging_directory;

    # This is the SAM file you should be appending to.  Dont forget, no headers!
    my $sam_file = $scratch_directory . "/all_sequences.sam";
    # import format
    #my $import_format = $self->instrument_data->import_format; # TODO what is this supposed to be?

    # decompose aligner params for each stage of alignment
    my $aligner_params = $self->decomposed_aligner_params;

    # get the command path
    my $bsmap_cmd_path = Genome::Model::Tools::Bsmap->path_for_bsmap_version($self->aligner_version);

    # Data is single-ended or paired-ended: key off the number of files passed in (1=SE, 2=PE)
    # Under no circumstances should you ever get more than 2 files, if you do then that's bad and
    # you should die.
    my $paired_end = 0;

    my $min_insert_size = -1;
    my $max_insert_size = -1;


    # seems dubious
    for (@input_pathnames) {
        $_ =~ s/^(.+\.(?:bam|sam))(?::[12])?$/$1/g;
    }

    if (@input_pathnames == 1) {
        $self->status_message("_run_aligner called in single-ended mode.");
    } elsif (@input_pathnames == 2) {
        $self->status_message("_run_aligner called in paired-end mode.");
        $paired_end = 1;

        # get the instrument-data so we can calculate the minimum and maximum insert size
        my $instrument_data = $self->instrument_data();

        if (defined($instrument_data->median_insert_size) && defined($instrument_data->sd_above_insert_size)){
            my $median_insert_size = $instrument_data->resolve_median_insert_size();
            my $sd_above_insert_size = $instrument_data->resolve_sd_insert_size();

            # TODO this may be an area for improvement
            if(defined($instrument_data->sd_below_insert_size)){
                my $sd_below_insert_size = $instrument_data->sd_below_insert_size();
                $min_insert_size = $median_insert_size - (3*$sd_below_insert_size);
            } else {
                $min_insert_size = $median_insert_size - (3*$sd_above_insert_size);
            }
            $max_insert_size = $median_insert_size + (3*$sd_above_insert_size);
        } else { #use defaults from bsmap
            $self->error_message("Unable to get insert size info from instrument data, using bsmap defaults: 28/500");
            $max_insert_size = 500;
            $min_insert_size = 28;
        }

    } else {
        die $self->error_message("_run_aligner called with " . scalar @input_pathnames . " files.  It should only get 1 or 2!");
    }

    ###################################################
    # run the bsmap aligner, output into temp sam file
    ###################################################
    $self->status_message("Running bsmap aligner.");
    my $temp_sam_output = $temporary_directory . "/mapped_reads.sam";

    #add no-header option for later versions
    if($self->aligner_version =~ /(\d+.\d+)/){
        if($1 >= 2.7 ){
            $aligner_params .= " -H";
        }
    }

    my $align_cmd = sprintf("%s %s -d %s %s -o %s",
        $bsmap_cmd_path,
        $paired_end
            ? "-a $input_pathnames[0] -b $input_pathnames[1] -m $min_insert_size -x $max_insert_size"
            : "-a $input_pathnames[0]",
        $reference_fasta_path,
        $aligner_params, # example: -p 4 (4 cores) -z 33 (initial qual char) -v 4 (max mismatches) -q 20 (qual trimming)
        $temp_sam_output
    );


    my $rv = Genome::Sys->shellcmd(
        cmd => $align_cmd,
        input_files => [@input_pathnames, $reference_fasta_path],
        output_files => [$temp_sam_output]
    );
    unless($rv) { die $self->error_message("Alignment failed."); }

    ###################################################
    # append temp sam file to all_sequences.sam
    ###################################################
    $self->status_message("Appending mapped_reads.sam to all_sequences.sam.");

    # This will also remove @header lines from the raw sam file. The duplicate
    # headers cause Picard crashes.
    my $append_cmd = sprintf(
        "grep -v -E '%s' %s >> %s", '^@', $temp_sam_output, $sam_file);

    $rv = Genome::Sys->shellcmd(
        cmd => $append_cmd,
        input_files  => [$temp_sam_output],
        output_files => [$sam_file],
        skip_if_output_is_present => 0 # because there will already be an all_sequences.sam we're appending to
        );

    unless($rv) { die $self->error_message("Appending failed."); }

    ###################################################
    # Sort
    ###################################################
    $self->status_message("Resorting all_sequences.sam by coordinate.");
    $self->_sort_sam($sam_file);

    ###################################################
    # clean up
    ###################################################

    # confirm that at the end we have a nonzero sam file, this is what'll get turned into a bam and copied out.
    unless (-s $sam_file) { die $self->error_message("The sam output file $sam_file is zero length; something went wrong."); }

    # TODO Any log files for staging directory?
    # TODO Any last minute checks?

    # If we got to here, everything must be A-OK.  AlignmentResult will take over from here
    # to convert the sam file to a BAM and copy everything out.
    return 1;
}

sub _disconnect_from_db {
    my ($self) = @_;

    $self->status_message("Closing data source db handle...");
    if ($self->__meta__->data_source->has_default_handle) {
        if ($self->__meta__->data_source->disconnect_default_handle) {
            $self->status_message("Disconnected data source db handle (as expected).");
        } else {
            $self->status_message("Unable to disconnect data source db handle.");
        }
    } else {
        $self->status_message("Data source db handle already closed.");
    }
}

sub _check_db_connection {
    my ($self) = @_;

    if ($self->__meta__->data_source->has_default_handle) {
        $self->status_message("Data source db handle unexpectedly reconnected itself.");
    } else {
        $self->status_message("Data source db handle still closed (as expected).");
    }
}

# Sort a sam file.
sub _sort_sam {
    my ($self, $given_sam) = @_;

    my $unsorted_sam = "$given_sam.unsorted";

    # Prepare sort command
    unless (move($given_sam, $unsorted_sam)) {
        die $self->error_message(
            "Unable to move $given_sam to $unsorted_sam. " .
            "Cannot proceed with sorting.");
    }

    my $picard_sort_cmd = Genome::Model::Tools::Picard::SortSam->create(
        sort_order             => 'coordinate',
        input_file             => $unsorted_sam,
        output_file            => $given_sam,
        max_records_in_ram     => 2000000,
        maximum_memory         => 8,
        maximum_permgen_memory => 256,
        temp_directory         => $self->temp_scratch_directory,
        use_version            => $self->picard_version,
    );

    # Disconnect from db
    $self->_disconnect_from_db();

    # Run sort command
    unless ($picard_sort_cmd and $picard_sort_cmd->execute) {
        die $self->error_message(
            "Failed to create or execute Picard sort command.");
    }

    # Hopefully we're still disconnected
    $self->_check_db_connection();

    # Clean up
    unless (unlink($unsorted_sam)) {
        $self->status_message("Could not unlink $unsorted_sam.");
    }

    return $given_sam;
}

# TODO: This should be the command line used to run your aligner, not just the params.
# sorry for the mis-named method name, it'll get fixed soon.
#
# This will end up in our re-composed sam/bam header.
sub aligner_params_for_sam_header {
    my $self = shift;

    my $aligner_params = $self->decomposed_aligner_params();

    # the number of processors does not affect the sam file
    $aligner_params =~ s/-p\s+\d+//;

    # recompact white space
    $aligner_params =~ s/(^)?(?(1)\s+|\s+(?=\s|$))//g;

    return "bsmap $aligner_params";
}

sub decomposed_aligner_params {
    # TODO if the user specifies params in the processing profile but neglects to require a necessary param
    #       this will not automatically fill it in
    my $self = shift;

    my $full_params;

    if (ref($self)) { # if this is an instance of AlignmentResult
        $full_params = $self->aligner_params();
    } else {
        $full_params = shift;
    }

    # split a colon-delimited list of arguments
    #my @params = split(":", $full_params || "::");

    my $defaults = ("-p 4 -q 20 -v 4 -z ! -R");
    # -p is the number of processors to give it
    # -z ! is required to get the correct quality score output
    # -v is the number of mismatches to allow
    # -q is the minimum quality score used in trimming


    # create our params hash, using default arguments if none were supplied
    my $aligner_params = $full_params || $defaults;

    # if # of processors is not specified, override it to 4 (bsmap defaults to 1)
    $aligner_params .= " -p 4" unless $aligner_params =~ /-p\s+\d/;

    # attemp to compact and sort command-line arguments for consistency:
    # compacts strings of whitespace down to a single character; strips all white space from beginning and end of string
    $aligner_params =~ s/(^)?(?(1)\s+|\s+(?=\s|$))//g;
    # split by each argument, sort, rejoin
    $aligner_params = join(" ",sort(split(/\s(?=-)/, $aligner_params)));

    return $aligner_params;
}

sub prepare_reference_sequence_index {
    my $class = shift;

    $class->status_message("BSMAP doesn't require index preparation, doing nothing.");

    return 0;
}

sub fillmd_for_sam {
    return 1;
}

sub requires_read_group_addition {
    return 1;
}

# only accept bam input if we're not running in force fragment mode
sub accepts_bam_input {
    my $self = shift;
    return $self->force_fragment ? 0 : 1;
}

sub supports_streaming_to_bam {
    return 0;
}
