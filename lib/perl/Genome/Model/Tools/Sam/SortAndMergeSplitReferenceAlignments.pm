package Genome::Model::Tools::Sam::SortAndMergeSplitReferenceAlignments;

use strict;
use warnings;

use Genome;
use Command;
use IO::File;
use File::Basename;
use File::Spec;
use Sys::Hostname;
use Genome::Utility::AsyncFileSystem qw(on_each_line);


class Genome::Model::Tools::Sam::SortAndMergeSplitReferenceAlignments {
    is  => 'Genome::Model::Tools::Sam',
    has => [
        input_files => {
            is => 'Text',
            is_many => 1,
            doc => 'BAM input files to merge'
        },
        output_file => {
            is => 'Text',
            doc => 'merged BAM file to write'
        },
        keep_sorted_files => {
            is => 'Boolean',
            default => 1,
        },
        type => {
            is => 'Text',
            is_optional => 1,
            valid_values => ['unsorted', 'fragment', 'paired_end', 'paired_end_and_fragment'],
            default => 'paired_end_and_fragment',
            doc => 'note that fragment, paired_end, and paired_end_and_fragment require that input be sorted by read name'
        },
        softclip_mismatch => {
            is => 'Boolean',
            is_optional => 1,
            default=>0,
            doc => 'If this parameter is specified, softclipped bases will be added to the mismatch count from the NM: tag',
        },
    ],
};

sub help_brief {
    'Tool to sort and merge BAM or SAM files aligned against a split reference (eg one that was larger than 4GiB before being broken up)';
}

sub help_detail {
    return 'Tool to sort and merge BAM or SAM files aligned against a split reference (eg one that was larger than 4GiB before being broken up).  ' . 
           'Writes a .merge_ok file to the output directory.  If this file exists, this command is a no-op.  Note that temporary/intermediate files ' .
           'made by samtools sort are kept in the same directory as the output file.  If the input file is compressed, budget approximately 12x ' . 
           'the input filesize for the duration of the merging process.';
}

sub execute {
    my $self = shift;

    $self->dump_status_messages(1);
    $self->dump_error_messages(1);

    my @in_files = $self->input_files;
    my $out_file = $self->output_file;

    # sort bam files by name for merge
    my $counter;
    my @bams_sorted;
    for my $bam (@in_files) {
        $counter++;
        my $sorted_bam = $self->output_file.".sorted_part_".$counter;

        #TODO this shortcutting step should be removed, as it's not guaranteed correct, but it's in place for hmp pipeline shortcutting for now

        my $rv;
        eval{
            my $skip;
            #if (-e $sorted_bam){
            #    my $bam_lc = `samtools view $bam | wc -l`;
            #    my $sort_lc = `samtools view $sorted_bam | wc -l`;
            #    $skip = 1 if $bam_lc == $sort_lc;
            #}
            unless ($skip){
                my $command = Genome::Model::Tools::Sam::SortBam->create(
                    file_name => $bam,
                    name_sort => 1,
                    output_file => $sorted_bam,
                    maximum_memory => 3000000000,
                );
                $rv=$command->execute;
            }
            push @bams_sorted, $sorted_bam . '.bam';
        };
        if ($@ or !$rv){
            $self->error_message("Fail to sorted bam $bam: $@");
            for (@bams_sorted){
                unlink $_;
            }
            die "Fail to sorted bam $bam: $@";
        }
    }

    print "DEBUG: Starting merge...\n";
    # merge bam files
    my $rv;
    eval{ 
        my $command = Genome::Model::Tools::Sam::MergeSplitReferenceAlignments->create(
            input_files => \@bams_sorted,
            output_file => $out_file,
            softclip_mismatch => $self->softclip_mismatch,
        );
        $rv = $command->execute;
    };
    if ($@ or !$rv){
        $self->error_message("Fail to merge sorted bams: $@");
        for (@bams_sorted){
            unlink $_ unless $self->keep_sorted_files;
        }
        die "Fail to merge sorted bams: $@";
    }
    for (@bams_sorted){
        unlink $_ unless $self->keep_sorted_files;
    }
    return 1;
}

1;
