package Genome::Model::RnaSeq::DetectFusionsResult::ChimerascanResult;

use strict;
use warnings;

use above 'Genome';

class Genome::Model::RnaSeq::DetectFusionsResult::ChimerascanResult {
    is => "Genome::Model::RnaSeq::DetectFusionsResult::ChimerascanBase",
    has_input => [
        original_bam_paths => {
            is => "Text",
            is_many => 1,
            doc => "The path(s) to the original instrument_data BAM files."
        }
    ],
};


sub _resolve_original_files {
    my ($self, $reuse_bam) = @_;

    my @fastq_files;
    if ($reuse_bam) {
        return $self->_resolve_original_files_reusing_bam();
    } else {
        unless ($self->original_bam_paths) {
            die("Couldn't find 'original_bam_paths' to make fastq files!");
        }
        # get fastq1/2 from the BAMs
        my (@fastq1_files, @fastq2_files);
        for my $bam_path ($self->original_bam_paths) {
            my $tmp_dir = File::Temp::tempdir('tempXXXXX',
                DIR => $self->temp_staging_directory,
                CLEANUP => 1
            );
            my $queryname_sorted_bam = File::Spec->join($tmp_dir,
                    'original_queryname_sorted.bam');
            $self->_qname_sort_bam($bam_path, $queryname_sorted_bam);

            # make fastqs from the qname sorted bam
            my $fastq1 = File::Spec->join($tmp_dir, "original_fastq1");
            my $fastq2 = File::Spec->join($tmp_dir, "original_fastq2");

            $self->_convert_bam_to_fastqs($queryname_sorted_bam, $fastq1, $fastq2);

            push @fastq1_files, $fastq1;
            push @fastq2_files, $fastq2;
        }

        # concatinate forward/reverse fastqs together
        my $fastq1 = File::Spec->join($self->temp_staging_directory, 'fastq1');
        my $cmd = sprintf('cat %s > %s', join(" ", @fastq1_files), $fastq1);
        Genome::Sys->shellcmd(
            cmd => $cmd,
            input_files => [@fastq1_files],
            output_files => [$fastq1],
        );

        my $fastq2 = File::Spec->join($self->temp_staging_directory, 'fastq2');
        $cmd = sprintf('cat %s > %s', join(" ", @fastq2_files), $fastq2);
        Genome::Sys->shellcmd(
            cmd => $cmd,
            input_files => [@fastq2_files],
            output_files => [$fastq2],
        );
        return ($fastq1, $fastq2, undef);
    }
}

sub _run_chimerascan {
    my ($self, $bowtie_version, $c_pargs, $c_opts) = @_;

    my $python_path = $self->_python_path_for_version($self->version);
    local $ENV{PYTHONPATH} = ($ENV{PYTHONPATH} ? $ENV{PYTHONPATH} . ":" : "") .
            $python_path;
    $self.status_message("Added '$python_path' to PYTHONPATH");

    my $cmd_path = $self->_path_for_version($self->version);
    unless ($cmd_path) {
        die $self->error_message("Failed to find a path for chimerascan for " .
                "version " . $self->version . "!");
    }

    my $bowtie_path = Genome::Model::Tools::Bowtie->base_path($bowtie_version);

    my $arguments_str = join(" ", @{$c_pargs});
    my $output_directory = $c_pargs->[-1];
    my $cmd = "python $cmd_path/chimerascan_run.py -v $c_opts " .
              "--bowtie-path=" . $bowtie_path .
              " $arguments_str > $output_directory/chimera_result.out";

    Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files => $c_pargs,
    );
}

sub _path_for_version {
    my ($class, $version) = @_;

    my ($base_path, $bin_sub, $python_sub) =
            _get_chimerascan_path_for_version($version);
    my $result = join('/', $base_path, $bin_sub);

    if (-e $result) {
        return $result;
    } else {
        die("Binary path ($bin_sub) does not exist under chimerascan path " .
                "($base_path)!");
    }
}

sub _python_path_for_version {
    my ($class, $version) = @_;

    my ($base_path, $bin_sub, $python_sub) =
            _get_chimerascan_path_for_version($version);
    my $result = join('/', $base_path, $python_sub);

    if (-e $result) {
        return $result;
    } else {
        die("Python path ($python_sub) does not exist under chimerascan " .
            "path ($base_path)!");
    }
}

sub _get_chimerascan_path_for_version {
    my ($version) = @_;

    my $path = sprintf("/usr/lib/chimerascan%s", $version);
    if (-e $path) {
        # installed via deb-package method
        my $bin_sub = 'bin';
        my $python_sub = join('/', 'lib', 'python2.6', 'site-packages');
        return $path, $bin_sub, $python_sub;
    } else {
        # fall back to old installation method
        $path = $ENV{GENOME_SW} . "/chimerascan/chimerascan-$version";
        if (-e $path) {
            my $bin_sub = 'chimerascan';
            my $python_sub = join('/', 'build', 'lib.linux-x86_64-2.6');
            return $path, $bin_sub, $python_sub;
        } else {
            die("You requested an unavailable version of Chimerascan. " .
                "Requested: $version");
        }
    }
}

sub _chimerascan_index_cmd {
    return 'chimerascan-index';
}

sub _chimerascan_result_class {
    return 'Genome::Model::RnaSeq::DetectFusionsResult::ChimerascanResult';
}

1;
