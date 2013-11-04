package Genome::Model::RnaSeq::DetectFusionsResult::Chimerascan::VariableReadLength::Result;

use strict;
use warnings;

use above 'Genome';

class Genome::Model::RnaSeq::DetectFusionsResult::Chimerascan::VariableReadLength::Result {
    is => "Genome::Model::RnaSeq::DetectFusionsResult::Chimerascan::ResultBase",
};

sub get_executable_path {
    my ($self, $version) = @_;

    my $executable = File::Which::which("chimerascan-vrl$version");

    unless (-e $executable) {
        die $self->error_message("Failed to find a path to installed " .
                "chimerascan-vrl for version '" . $self->version . "'!");
    }
    return $executable;
}

sub _run_chimerascan {
    my ($self, $bowtie_version, $c_pargs, $c_opts) = @_;

    my $executable = $self->get_executable_path($self->version);
    my $bowtie_path = Genome::Model::Tools::Bowtie->base_path($bowtie_version);

    my $arguments_str = join(" ", @{$c_pargs});
    my $output_directory = $c_pargs->[-1];
    my $cmd = "$executable chimerascan_run.py -v $c_opts " .
              "--bowtie-path=" . $bowtie_path .
              " $arguments_str > $output_directory/chimera_result.out";

    Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files => $c_pargs,
    );
}

sub _chimerascan_index_cmd {
    return 'chimerascan-vrl-index';
}

sub _chimerascan_result_class {
    return 'Genome::Model::RnaSeq::DetectFusionsResult::Chimerascan::VariableReadLength::Result';
}

1;
