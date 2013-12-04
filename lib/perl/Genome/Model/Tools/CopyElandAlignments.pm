package Genome::Model::Tools::CopyElandAlignments;

use strict;
use warnings;
use Genome;
use Carp;

# AUTHOR : Todd Wylie <twylie@wustl.edu>
# DATE   : Mon Mar 29 12:54:28 CDT 2010

class Genome::Model::Tools::CopyElandAlignments {
    is  => 'Command',
    has => [
            flowcell_id => {
                            is  => 'Text',
                            doc => "Provide a flow cell id for archiving Eland s_N_sorted.txt files.",
                           },
            outdir      => {
                            is  => 'Text',
                            doc => "Output directory where files will be copied.",
                           },
           ],
  has_optional => [
                   fastq => {
                             is  => 'Boolean',
                             doc => "[OPTIONAL] Copy in corresponding FASTQ files as well.",
                            },
                  ],
      has_many => [
                   lanes => {
                             is_optional => 1,
                             is          => 'Text',
                             doc         => "[OPTIONAL] Choose specific lanes to copy... default is all.",
                            },
                  ],
              };


sub help_brief { "This script will copy (non-multiplexed) Eland s_N_sorted.txt files to an archive directory if they still exist in a Gerald run directory structure." }


sub help_synopsis {
    return <<EOS
gmt copy-eland-alignments --flowcell-id=<string> --output=<string> --fastq --lanes=1,2,3
EOS
}


sub help_detail {
    return <<EOS
Eland alignment information is transient within Gerald run directories; the
information usually lasts about 1 week post-processing before being removed. To
archive (i.e., copy) pertinent Eland alignment information (s_N_sorted.txt)
files, this script will copy files to a given out directory given a single FLOW
CELL id. CURRENTLY, ONLY WORKS ON NON-MULTIPLEXED RUNS.
** WARNING ** You will probably need ~100Gb of disk space for the copy.
EOS
}


sub execute {
    my $self = shift;

    # Navigate to pipeline directory given the passed FLOW CELL id. Does the
    # Eland alignment still exist? Nothing fancy here.

    my @possible_files = qw(
                               s_1_1_sorted.txt
                               s_1_2_sorted.txt
                               s_2_1_sorted.txt
                               s_2_2_sorted.txt
                               s_3_1_sorted.txt
                               s_3_2_sorted.txt
                               s_4_1_sorted.txt
                               s_4_2_sorted.txt
                               s_5_1_sorted.txt
                               s_5_2_sorted.txt
                               s_6_1_sorted.txt
                               s_6_2_sorted.txt
                               s_7_1_sorted.txt
                               s_7_2_sorted.txt
                               s_8_1_sorted.txt
                               s_8_2_sorted.txt
                          );

    # File EVAL.
    my (@eland_files, @fastq_files);
    my @lane_paths = Genome::InstrumentData::Solexa->get( flow_cell_id => $self->flowcell_id() );
    foreach my $lane_path (@lane_paths) {
        foreach my $possible_file (@possible_files) {
            unless ($self->lanes()) {
                my $eval_file = $lane_path->gerald_directory() . '/' . $possible_file;
                if (-e $eval_file) {
                    push (@eland_files, $eval_file);
                    print 'EVAL: ' . $eval_file . "\n";
                }
            }
            else {
                foreach my $lane ($self->lanes()) {
                    my $lane_eval = 's_' . $lane;
                    if ($possible_file =~ /$lane_eval/) {
                        my $eval_file = $lane_path->gerald_directory() . '/' . $possible_file;
                        if (-e $eval_file) {
                            push (@eland_files, $eval_file);
                            print 'EVAL: ' . $eval_file . "\n";
                        }
                    }
                }
            }
        }
    }
    my $file_count = scalar( @eland_files );
    if ($file_count == 0) {
        croak "Only $file_count number of files found. Too few, please investigate.";
    }
    unless ($self->lanes()) {
        if ($file_count < 8) {
            croak "Only $file_count number of files found. Too few, please investigate.";
        }
    }
    print "(FILE COUNT = $file_count)\n";


    # Output directory.
    if (-d $self->outdir()) { croak "Output directory already exists." }
    my $outdir_path = $self->outdir() . '/' . $self->flowcell_id();
    mkdir( $self->outdir(), 0776) or croak "Could not make output directory: " . $self->outdir();
    mkdir( $outdir_path,    0776) or croak "Could not make output directory: $outdir_path";

    # README.txt file logging.
    my $readme_file = $outdir_path . '/' . 'README.txt';
    open (README, ">$readme_file") or croak "Could not write README.txt file to output directory.";
    print README join ("\n", @eland_files);
    print README join ("\n", @fastq_files) if ($self->fastq());
    close (README);

    # LSF bsub script to copy files into directory.
    foreach my $eland_file (@eland_files) {
        my @fields      = split (/\//, $eland_file);
        my $sorted_file = pop( @fields );
        my $bsub_script = $outdir_path . '/' . $sorted_file . '.sh';
        my $error_log   = $outdir_path . '/' . $sorted_file . '.error';
        open (SCRIPT, ">$bsub_script") or croak "Could not write the bsub script file: $bsub_script";
        if ($self->fastq()) {
            my $fastq_file = $eland_file;
            $fastq_file    =~ s/sorted/sequence/;
            print SCRIPT join (
                               q/ /,
                               'cp ',
                               $fastq_file,
                               $outdir_path . '/',
                              ) . "\n";
            print SCRIPT join (
                               q/ /,
                               'cp ',
                               $eland_file,
                               $outdir_path . '/',
                               "\ntouch " . $outdir_path . '/' . $sorted_file . '.FINISHED',
                              ) . "\n";
            close (SCRIPT);
        }
        else {
            print SCRIPT join (
                               q/ /,
                               'cp ',
                               $eland_file,
                               $outdir_path . '/',
                               "\ntouch " . $outdir_path . '/' . $sorted_file . '.FINISHED',
                              ) . "\n";
            close (SCRIPT);
        }
        system ("bsub -q $ENV{GENOME_LSF_QUEUE_BUILD_WORKER} -e $error_log \'sh $bsub_script\'\n");  # launch the bsub job
    }

    return 1;
}


1;  # end of package
