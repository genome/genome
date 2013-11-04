package Genome::Model::RnaSeq::DetectFusionsResult::Chimerascan::FixedReadLength::Result;

use strict;
use warnings;

use above 'Genome';

class Genome::Model::RnaSeq::DetectFusionsResult::Chimerascan::FixedReadLength::Result {
    is => "Genome::Model::RnaSeq::DetectFusionsResult::Chimerascan::ResultBase",
};

sub _run_chimerascan {
    my ($self, $bowtie_version, $c_pargs, $c_opts) = @_;

    my $python_path = $self->_python_path_for_version($self->version);
    local $ENV{PYTHONPATH} = ($ENV{PYTHONPATH} ? $ENV{PYTHONPATH} . ":" : "") .
            $python_path;
    $self->status_message("Added '$python_path' to PYTHONPATH");

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
    return 'Genome::Model::RnaSeq::DetectFusionsResult::Chimerascan::FixedReadLength::Result';
}

1;
