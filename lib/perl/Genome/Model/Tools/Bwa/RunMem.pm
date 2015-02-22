package Genome::Model::Tools::Bwa::RunMem;

use strict;
use warnings;
use Genome;
use File::Slurp qw(read_file write_file);
use List::AllUtils qw(any);
use Text::ParseWords qw(shellwords);

my $REQUIRED_SAMTOOLS_VERSION = '0.1.19';
my %VALID_SAMTOOLS_VERSIONS = ($REQUIRED_SAMTOOLS_VERSION => 1);

my $REQUIRED_BEDTOOLS_VERSION = '2.17.0';
my %VALID_BEDTOOLS_VERSIONS = ($REQUIRED_BEDTOOLS_VERSION => 1);

my $MIN_PER_THREAD_SORT_MEMORY_MB = 50;

# Samtools likes to use quite a bit more than memory you tell it to.
# I evaluated this on very short reads (36bp) which means a high bam record/file size ratio.
# I think part of the problem in samtools is not accounting for memory used to track the
# records themselves, so this may be overly conservative for longer read lengths.
my $SAMTOOLS_MEMORY_CORRECTION = 0.75;

class Genome::Model::Tools::Bwa::RunMem {
    is => "Command::V2",
    has_input => [
        bwa_version => {
            is => "String",
            doc => "Bwa version to use",
        },

        samtools_version => {
            is => "String",
            doc => "The version of samtools to use for sorting and calmd",
            valid_values => [sort keys %VALID_SAMTOOLS_VERSIONS],
            default_value => $REQUIRED_SAMTOOLS_VERSION,
        },

        bedtools_version => {
            is => "String",
            doc => "Bedtools version to use (for bam => fastq conversion)",
            valid_values => [sort keys %VALID_BEDTOOLS_VERSIONS],
            default_value => $REQUIRED_BEDTOOLS_VERSION,
        },

        input_files => {
            is => "String",
            is_many => 1,
            doc => "Unaligned read input fastq files (1 or 2), or a single bam file",
        },

        is_bam => {
            is => "Boolean",
            doc => "The input is a single bam file (add -p to --aligner-params for PE bams)",
            default_value => 0,
        },

        bam_exclude_flags => {
            is => "Boolean",
            doc => "Flag mask to exclude when streaming from bam files",
            default_value => 0x100 | 0x200 | 0x800,
        },

        limit_to_read_group => {
            is => "String",
            doc => "Only align reads in the given read group (requires --is-bam)",
            is_optional => 1,
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

        max_sort_memory_mb => {
            is => "String",
            doc => "Maximum amount of memory to use while sorting (MB). ",
            default_value => "4096",
        },

        num_threads => {
            is => "Integer",
            doc => "Number of alignment/sorting threads",
            default_value => 1,
        },

        temp_dir => {
            is => "String",
            doc => "Temp directory to use (one is generated if omitted)",
            is_optional => 1,
        },
    ]
};

sub help_brief {
    "Tool to run the bwa mem aligner";
}

sub help_synopsis {
<<EOS
PAIRED END ALIGNMENT
    FASTQ:
        gmt bwa run-mem \\
            --bwa-version 0.7.10 \\
            --aligner-log-path bwa.log \\
            --indexed-fasta /path/to/ref.fa \\
            --aligner-index-fasta /path/to/aligner_index_ref.fa \\
            --input-files r1.fq,r2.fq \\
            --sam-header-path header.sam \\
            --aligner-params '-R "\@RG\\tID:rg1\\tLB:lb1\\tSM:sm1"'

    BAM (note the -p in aligner params):
        gmt bwa run-mem \\
            --bwa-version 0.7.10 \\
            --aligner-log-path bwa.log \\
            --indexed-fasta /path/to/ref.fa \\
            --aligner-index-fasta /path/to/aligner_index_ref.fa \\
            --input-files reads.bam \\
            --sam-header-path header.sam \\
            --aligner-params '-p -R "\@RG\\tID:rg1\\tLB:lb1\\tSM:sm1"'

SINGLE END ALIGNMENT
    FASTQ:
        gmt bwa run-mem \\
            --bwa-version 0.7.10 \\
            --aligner-log-path bwa.log \\
            --indexed-fasta /path/to/ref.fa \\
            --aligner-index-fasta /path/to/aligner_index_ref.fa \\
            --input-files reads.fq \\
            --sam-header-path header.sam \\
            --aligner-params '-R "\@RG\\tID:rg1\\tLB:lb1\\tSM:sm1"'

    BAM (note the lack of -p in aligner params):
        gmt bwa run-mem \\
            --bwa-version 0.7.10 \\
            --aligner-log-path bwa.log \\
            --indexed-fasta /path/to/ref.fa \\
            --aligner-index-fasta /path/to/aligner_index_ref.fa \\
            --input-files reads.bam \\
            --sam-header-path header.sam \\
            --aligner-params '-R "\@RG\\tID:rg1\\tLB:lb1\\tSM:sm1"'

EOS
}

sub help_detail {
<<EOS
 Aligns reads (FASTQ or BAM) in --input-files with bwa mem producing a
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
    return Genome::Model::Tools::Sam->path_for_samtools_version($self->samtools_version);
}

sub _bedtools_path {
    my $self = shift;

    # The bedtools module returns the path to the bedtools dir, not the executable
    return File::Spec->catfile(
        Genome::Model::Tools::BedTools->path_for_bedtools_version($self->bedtools_version),
        "bin",
        "bedtools")
        ;
}

sub _param_join {
    return join(" \\\n        ", @_);
}

sub _validate_params {
    my $self = shift;

    if ($self->num_threads < 1) {
        die $self->error_message("--num-threads must be >= 1!");
    }

    if ($self->limit_to_read_group && !$self->is_bam) {
        die $self->error_message("--limit-to-read-group requires --is-bam!");
    }

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

sub _script_pre_lines {
    my $self = shift;

    return <<EOS;
exec 3>&1 # Save stdout to fd3 (for bedtools)
EOS
}

sub _script_post_lines {
    my $self = shift;

    return <<EOS;
exec 3>&- # Close fd3
EOS
}

sub _bam_input_file {
    my $self = shift;

    die $self->error_message("Can't stream something that isn't a bam!") unless $self->is_bam;
    my @input_files = $self->input_files;

    die $self->error_message("When aligning bam files, only one input is allowed (use -p in --aligner params for PE)")
        unless @input_files == 1;

    return $input_files[0];
}

sub _limit_by_read_group_and_flags_commands {
    my $self = shift;

    my $bam_path = $self->_bam_input_file;

    my $rg_limiting_opts = '';

    if ($self->limit_to_read_group) {
        $rg_limiting_opts = sprintf "-r%s", $self->limit_to_read_group;
    }

    return join " ",
        $self->_samtools_path,
        "view",
        "-u",
        "-F", $self->bam_exclude_flags,
        $rg_limiting_opts,
        $bam_path,
        ;
}

sub _bam_to_fastq_commands {
    my $self = shift;

    my $bam_path = $self->_bam_input_file;

    my $bad_sort_msg = sprintf "ERROR: Input bam %s is not in name sorted order, abort!",
        $bam_path;

    my $path = $self->_bedtools_path;

    # Bedtools does not /enforce/ that the input bam is in name sorted order.
    # It does, however, warn about such. We'll juggle file descriptors around a
    # bit to enable detecting these warnings on stderr and failing early.

    return <<EOS;
{
    # swap stdin and stderr
    $path bamtofastq -i /dev/stdin -fq /dev/stdout -fq2 /dev/stdout 3>&1 1>&2 2>&3 \\
        | {
            tee >(
                if grep -q 'WARNING: Query'
                then
                    echo '$bad_sort_msg'
                    exit 1
                fi
                ) \\
            ;
        } 2>&1 \\
    ;
# swap stdin and stderr again (back to their original configuration)
} 3>&1 1>&2 2>&3 \\
EOS
}

sub _aligner_commands {
    my $self = shift;
    my $ver = $self->bwa_version;
    if (!Genome::Model::Tools::Bwa->supports_mem($ver)) {
        die $self->error_message("Bwa version '$ver' does not support bwa mem");
    }
    my $bwa = Genome::Model::Tools::Bwa->path_for_bwa_version($ver);
    my $log = $self->aligner_log_path;
    my @extra_params = ($self->aligner_params, sprintf("-t %d", $self->num_threads));

    my @cmds;

    my @in = $self->input_files;
    if ($self->is_bam) {
        push @cmds,
            $self->_limit_by_read_group_and_flags_commands,
            $self->_bam_to_fastq_commands
            ;

        @in = ("-");
    }

    my $fa = $self->aligner_index_fasta;
    push @cmds, _param_join("$bwa mem", @extra_params, $fa, @in, "2> $log");

    return @cmds;
}

sub _sam_replace_header_commands {
    my $self = shift;
    my $hdr = $self->sam_header_path;
    my $perl_interp = $^X;

    # grind off any SAM header lines (that start with @)
    my $body = qq!while (<>) { /^[^@]/ && last } print; print while(<>)!;
    return "{ cat $hdr && $perl_interp -e '$body' ; }";
}

sub _sam_to_uncompressed_bam_commands {
    my $self = shift;
    my $samtools_path = $self->_samtools_path;
    return "$samtools_path view -Su -";
}

sub _sort_memory_per_thread_mb {
    my $self = shift;
    my $mem = $self->max_sort_memory_mb;
    my $threads = $self->num_threads;

    my $total = $SAMTOOLS_MEMORY_CORRECTION * $mem;
    my $missing = $mem - $total;
    my $per_thread = int($total / $threads + 0.5);

    if ($per_thread < $MIN_PER_THREAD_SORT_MEMORY_MB) {
        die $self->error_message("Too many threads ($threads) to divide ${total}MB " .
            "sort memory amongst (${missing}MB reserved for samtools overhead)");
    }

    return $per_thread;
}

sub _sort_commands {
    my $self = shift;
    my $samtools_path = $self->_samtools_path;
    my $comp_level = 0;
    my $threads = $self->num_threads;
    my $per_thread_mem = $self->_sort_memory_per_thread_mb;
    my $tmp_path = $self->temp_dir;

    $tmp_path = Genome::Sys->create_temp_directory unless $tmp_path;
    $tmp_path = File::Spec->catfile($tmp_path, "samtools-sort");

    $self->debug_message("max per-thread sorting memory is ${per_thread_mem}MB");

    return "$samtools_path sort -l $comp_level -\@ $threads -m ${per_thread_mem}M -o - $tmp_path";
}

sub _calmd_commands {
    my $self = shift;
    my $fasta = $self->indexed_fasta;
    my $samtools_path = $self->_samtools_path;
    return "$samtools_path calmd -u - $fasta";
}

sub _final_bam_conversion_commands {
    my $self = shift;
    my $samtools_path = $self->_samtools_path;
    my $num_threads = $self->num_threads;
    return "$samtools_path view -b -\@ $num_threads -";
}

sub _pipeline_commands {
    my $self = shift;
    return (
        $self->_aligner_commands,
        $self->_sam_replace_header_commands,
        $self->_sam_to_uncompressed_bam_commands,
        $self->_sort_commands,
        $self->_calmd_commands,
        $self->_final_bam_conversion_commands,
        );
}

sub execute {
    my $self = shift;
    $self->_validate_params;

    $self->debug_message(sprintf("bwa version: %s", $self->bwa_version));
    $self->debug_message(sprintf("samtools version: %s", $self->samtools_version));
    $self->debug_message(sprintf("bedtools version: %s", $self->bedtools_version));

    my @input_files = $self->input_files;
    my $n_inputs = scalar @input_files;
    if ($n_inputs < 1 || $n_inputs > 2) {
        die $self->error_message(
            "Don't know what to do with $n_inputs input files! " .
            "Give one for single or two for paired end alignment."
            );
    }

    $self->status_message(
        "Warnings below from samtools sort/calmd about truncated bam files " .
        "may be safely ignored. They are caused by piping data to samtools " .
        "(it can't seek to check for the bam EOF marker)."
        );

    my %params = (
        interpreter => "/bin/bash",
        pre_commands => [$self->_script_pre_lines],
        post_commands => [$self->_script_post_lines],
        pipe_commands => [$self->_pipeline_commands],
        redirects => sprintf("> %s", $self->output_file),
        );

    my $pipe = Genome::Sys::ShellPipeline->create(%params);
    return $pipe->execute;
}

1;
