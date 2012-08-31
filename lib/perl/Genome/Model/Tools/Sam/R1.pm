package Genome::Model::Tools::Sam::R1;

use strict;
use warnings;

use Genome;
use Command;
use IO::File;
use File::Basename;
use Sys::Hostname;
use Genome::Utility::AsyncFileSystem qw(on_each_line);

class Genome::Model::Tools::Sam::R1 {
    is  => 'Genome::Model::Tools::Sam',
    has => [
        input_files => {
            is => 'Text',
            is_many => 1,
            doc => 'SAM or BAM input files to merge'
        },
        input_format => {
            is => 'Text',
            valid_values => ['SAM', 'BAM']
        },
        output_file => {
            is => 'Text',
            doc => 'merged SAM or BAM file to write'
        },
        output_format => {
            is => 'Text',
            valid_values => ['SAM', 'BAM']
        },
        type => {
            is => 'Text',
            is_optional => 1,
            valid_values => ['unsorted', 'fragment', 'paired_end', 'paired_end_and_fragment'],
            default => 'paired_end_and_fragment',
            doc => 'note that fragment, paired_end, and paired_end_and_fragment require that input be sorted by read name'
        }
    ],
};

sub help_brief {
    'Tool to merge BAM or SAM files aligned against a split reference (eg one that was larger than 4GiB before being broken up)';
}

sub help_detail {
    return 'Tool to merge BAM or SAM files aligned against a split reference (eg one that was larger than 4GiB before being broken up)';
}

sub _make_bamsam_arg {
    my ($arg, $file, $format);
    ($arg, $file, $format) = @_;
    return '--' . $arg . '-' . lc($format) . '=\"' . $file . '\"';
}

sub _merge_command {
    my $type = shift;
    my $in_files = shift;
    my $in_format = shift;
    my $out_file = shift;
    my $out_format = shift;

    my $cmd = 'bash -c "LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/gsc/var/tmp/genome/lib:/gsc/pkg/boost/boost_1_42_0/lib /gsc/var/tmp/genome/bin/merge-split-reference-alignments-v1.1.2';
    foreach my $in_file (@$in_files) {
        $cmd .= ' ' . _make_bamsam_arg('split', $in_file, $in_format);
    }
    $cmd .= ' ' . _make_bamsam_arg('merged', $out_file, $out_format);
    if($type eq 'unsorted') {
    }
    elsif($type eq 'fragment') {
        $cmd .= ' --select-best';
    }
    elsif($type eq 'paired_end') {
        $cmd .= ' --select-best-pair';
    }
    elsif($type eq 'paired_end_and_fragment') {
        $cmd .= ' --select-best-pair-and-fragment';
    }
    else {
        die;
    }
    $cmd .= '"';
    return $cmd;
}

sub execute {
    my $self = shift;

    $self->dump_status_messages(1);
    $self->dump_error_messages(1);

    my $type = $self->type;
    my @in_files = $self->input_files;
    my $in_format = $self->input_format;
    my $out_file = $self->output_file;
    my $out_format = $self->output_format;

    my $cmd = _merge_command($type, \@in_files, $in_format, $out_file, $out_format);

    return Genome::Sys->shellcmd(cmd => $cmd, input_files => \@in_files, output_files => [$out_file], skip_if_output_is_present => 0);
}

1;
