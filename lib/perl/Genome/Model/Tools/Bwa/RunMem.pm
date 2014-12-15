package Genome::Model::Tools::Bwa::RunMem;

use strict;
use warnings;
use Genome;
use File::Slurp qw(read_file write_file);
use List::AllUtils qw(any);
use Text::ParseWords qw(shellwords);

my $REQUIRED_SAMTOOLS_VERSION = '0.1.19';
my %VALID_SAMTOOLS_VERSIONS = ($REQUIRED_SAMTOOLS_VERSION => 1);

class Genome::Model::Tools::Bwa::RunMem {
    is => "Command::V2",
    has_input => [
        samtools_version => {
            is => "String",
            doc => "The version of samtools to use for sorting and calmd",
            valid_values => [sort keys %VALID_SAMTOOLS_VERSIONS],
            default_value => $REQUIRED_SAMTOOLS_VERSION,
        },

        temp_dir => {
            is => "String",
            doc => "Temp directory to use (one is generated if omitted)",
            is_optional => 1,
        },

        bwa_version => {
            is => "String",
            doc => "Bwa version to use",
        },

        input_fastqs => {
            is => "String",
            is_many => 1,
            doc => "Unaligned read input fastq files",
        },

        output_file => {
            is_output => 1,
            is => "String",
            doc => "Output file path",
        },

        aligner_params => {
            is => "String",
            doc => "Additional params for bwa mem",
            default_value => "",
        },

        aligner_log_path => {
            is_output => 1,
            is => "String",
            doc => "Path where aligner log (bwa mem stderr) will be written. " .
                "All other commands log to stderr.",
        },

        # Note: seems like our bwa mem aligner index files are not
        # faidx indexed.
        indexed_fasta => {
            is => "String",
            doc => "Path to faidx indexed reference sequence",
        },

        aligner_index_fasta => {
            is => "String",
            doc => "Path to aligner index fasta to align against",
        },

        sam_header_path => {
            is => "String",
            doc => "Path to sam header to use",
        },

        max_sort_memory => {
            is => "String",
            doc => "Maximum amount of memory to use while sorting. ".
                "Suffixes K/M/G are recognized",
            default_value => "4G",
        },

        num_threads => {
            is => "Integer",
            doc => "Number of alignment/sorting threads",
            default_value => 1,
        },
    ]
};

sub help_brief {
    "Tool to run the bwa mem aligner";
}

sub help_synopsis {
<<EOS
PAIRED END ALIGNMENT
    gmt bwa run-mem \\
        --bwa-version 0.7.10 \\
        --aligner-log-path bwa.log \\
        --indexed-fasta /path/to/ref.fa \\
        --aligner-index-fasta /path/to/aligner_index_ref.fa \\
        --input-files r1.fq,r2.fq \\
        --sam-header-path header.sam \\
        --aligner-params '-R "\@RG\\tID:rg1\\tLB:lb1\\tSM:sm1"'

SINGLE END ALIGNMENT
    gmt bwa run-mem \\
        --bwa-version 0.7.10 \\
        --aligner-log-path bwa.log \\
        --indexed-fasta /path/to/ref.fa \\
        --aligner-index-fasta /path/to/aligner_index_ref.fa \\
        --input-files reads.fq \\
        --sam-header-path header.sam \\
        --aligner-params '-R "\@RG\\tID:rg1\\tLB:lb1\\tSM:sm1"'
EOS
}

sub help_detail {
<<EOS
 Aligns FASTQ format reads in --input-files with bwa mem producing a
 coordinate-sorted output BAM file with MD tags corrected by samtools calmd.

 There are two fasta parameters, --indexed-fasta and --aligner-index-fasta.
 This is because our bwa mem aligner indexes are not faidx indexed, and
 samtools calmd requires such.

 Versions of bwa prior to 0.7.7 produced very incorrect MD tags. All current
 versions still produce tags that are not exactly like what samtools calmd
 generates. This is due to the way BWA encodes the reference sequence; its
 2bit encoding will replace Ns in the reference with random characters in
 "ACGT". In the future, we may decide we don't care and stop running calmd
 as part of this tool for bwa 0.7.7+.
EOS
}

sub _samtools_path {
    my $self = shift;

    # Note that the interface for sort is different in each of 0.1.18, 0.1.19, and 1.0
    # 1.0 performs a little better, but it is not installed.
    return Genome::Model::Tools::Sam->path_for_samtools_version("0.1.19");
}

sub _param_join {
    return join(" \\\n        ", @_);
}

sub _validate_params {
    my $self = shift;

    my @params = shellwords($self->aligner_params);
    # Check for -t or -t#
    my @thread_params = grep {
        my $x = $_;
        $x =~ s/['"]//g;
        $x =~ /^-t[0-9]*/
    } @params;

    if (@thread_params) {
        die $self->error_message(
            "This module does not support '-t <n>' as part of the " .
            "--aligner-params. Use --num-threads instead."
            );
    }

    # samtools 0.1.19 is required for this module to work properly
    if (!exists $VALID_SAMTOOLS_VERSIONS{$self->samtools_version}) {
        die $self->error_message(
            sprintf "Samtools version %s is invalid, supported versions: ",
            $self->samtools_version,
            join(", ", sort keys %VALID_SAMTOOLS_VERSIONS)
            );
    }
}

sub _aligner_command {
    my $self = shift;
    my $ver = $self->bwa_version;
    if (!Genome::Model::Tools::Bwa->supports_mem($ver)) {
        die $self->error_message("Bwa version '$ver' does not support bwa mem");
    }
    my $bwa = Genome::Model::Tools::Bwa->path_for_bwa_version($ver);
    my $log = $self->aligner_log_path;
    my @extra_params = ($self->aligner_params, "-t", $self->num_threads);

    my @in = $self->input_fastqs;
    my $fa = $self->aligner_index_fasta;
    return _param_join("$bwa mem", @extra_params, $fa, @in, "2> $log");
}

sub _sam_replace_header_cmdline {
    my $self = shift;
    my $hdr = $self->sam_header_path;
    my $perl_interp = $^X;

    # grind off any SAM header lines (that start with @)
    my $body = qq!while (<>) { /^[^@]/ && last } print; print while(<>)!;
    return "(cat $hdr; $perl_interp -e '$body')";
}

sub _sam_to_uncompressed_bam_cmdline {
    my $self = shift;
    my $samtools_exe_path = $self->_samtools_path;
    return "$samtools_exe_path view -Su -";
}

sub _sort_cmdline {
    my $self = shift;
    my $samtools_exe_path = $self->_samtools_path;
    my $comp_level = 0;
    my $threads = $self->num_threads;
    my $max_mem = $self->max_sort_memory;
    my $tmp_path = sprintf("%s/samtools-sort", $self->temp_dir);

    return "$samtools_exe_path sort -l $comp_level -\@ $threads -m $max_mem -o -f - $tmp_path";
}

sub _calmd_cmdline {
    my $self = shift;
    my $fasta = $self->indexed_fasta;
    my $samtools_exe_path = $self->_samtools_path;
    return "$samtools_exe_path calmd -b - $fasta";
}

sub _make_pipeline {
    my @cmds = @_;
    return join(" \\\n    | ", @cmds);
}

sub _pipeline_commands {
    my $self = shift;
    return (
        $self->_aligner_command,
        $self->_sam_replace_header_cmdline,
        $self->_sam_to_uncompressed_bam_cmdline,
        $self->_sort_cmdline,
        $self->_calmd_cmdline
        );
}

sub _script_text {
    my ($self, $pipestatus_path) = @_;

    my @pipeline = $self->_pipeline_commands;

    my $cmd = sprintf("%s > %s", _make_pipeline(@pipeline), $self->output_file);
    return <<EOS
#!/bin/bash

set -o pipefail
RV=0
$cmd
PSTAT=\${PIPESTATUS[@]}
echo "Pipe status: \$PSTAT"

echo \$PSTAT > $pipestatus_path

for x in \$PSTAT
do
    if [ "\$x" != "0" ]
    then
        exit 1
    fi
done

EOS
;
}

# This sub tries to parse a file containing a single line containing the
# result of ${PIPESTATUS[@]} from bash after executing the relevant pipeline.
# It should not be called unless an error is actually detected.
sub _parse_pipestatus {
    my ($self, $path) = @_;

    return "PIPESTATUS error file does not exist" unless -s $path;

    my @lines = read_file($path);
    chomp @lines;

    # there should only be one line in the file
    return "Too many lines in PIPESTATUS error file" unless scalar @lines == 1;

    my @status = split(/\s/, $lines[0]);

    my @commands = $self->_pipeline_commands;
    # the number of commands we intended to execute should be equal
    # to the number of exit codes we received.
    return "Wrong number of status codes in PIPESTATUS error file"
        unless scalar @status == scalar @commands;

    return "no failures" unless any { $_ != 0 } @status;

    # who dun goofed?
    my @failures = grep {$status[$_] != 0} 0..$#status;
    return sprintf "The following commands crashed:\n  %s",
        join("\n  ",
            map {sprintf "%s: %s", $_, $commands[$_]} @failures
            );
}

sub execute {
    my $self = shift;
    $self->_validate_params;

    my @input_fastqs = $self->input_fastqs;
    my $n_inputs = scalar @input_fastqs;
    if ($n_inputs < 1 || $n_inputs > 2) {
        die $self->error_message(
            "Don't know what to do with $n_inputs input files! " .
            "Give one for single or two for paired end alignment."
            );
    }

    if (!defined $self->temp_dir) {
        $self->temp_dir(Genome::Sys->create_temp_directory);
    }

    $self->status_message(sprintf("Temp directory is: %s", $self->temp_dir));

    my $script_path = sprintf("%s/run.sh", $self->temp_dir);
    my $pipestatus_path = sprintf("%s/pipestatus.txt", $self->temp_dir);
    write_file($script_path, $self->_script_text($pipestatus_path));

    $self->debug_message(
        sprintf
            "Executing script:\n\n" .
            "-------- BEGIN SCRIPT --------\n%s" .
            "--------- END SCRIPT ---------\n",
            join("", read_file($script_path))
        );

    $self->status_message(
        "Warnings below from samtools sort/calmd about truncated bam files " .
        "may be safely ignored. They are caused by piping data to samtools " .
        "(it can't seek to check for the bam EOF marker)."
        );

    eval {
        Genome::Sys->shellcmd(cmd => "bash $script_path");
    };

    if ($@) {
        my $msg = sprintf "Pipeline failed: %s",
            $self->_parse_pipestatus($pipestatus_path);

        die $self->error_message($msg);
    }

    return 1;
}

1;
