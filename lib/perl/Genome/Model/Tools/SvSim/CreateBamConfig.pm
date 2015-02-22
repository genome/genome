package Genome::Model::Tools::SvSim::CreateBamConfig;

use Genome;
use File::Basename qw(dirname);
use File::chdir;

use strict;
use warnings;

my $DEFAULT_VERSION = "1.4.2";

class Genome::Model::Tools::SvSim::CreateBamConfig {
    is => "Command::V2",
    has_input => [
        breakdancer_version => {
            is => "Text",
            doc => "The version of breakdancer to use",
            default_value => $DEFAULT_VERSION,
            is_optional => 1,
        },
        bam_files => {
            is => "Text",
            doc => "The list of bam files to examine",
            is_many => 1,
        },
        output_file => {
            is => "Text",
            is_output => 1,
            doc => "The output file. Note that histograms will also be placed in the containing directory",
        },
    ],
};

sub execute {
    my $self = shift;
    my $bd = Genome::Model::Tools::Breakdancer->create(
        use_version => $self->breakdancer_version
        );

    my @bams = $self->bam_files;
    my $exe = $bd->breakdancer_config_command;
    my $dir = dirname($self->output_file);
    my $cmd = join(" ", $exe, "-g", "-h", @bams, ">", $self->output_file);

    {
        local $CWD = $dir;
        return Genome::Sys->shellcmd(
            cmd => $cmd,
            input_files => \@bams,
            output_files => [$self->output_file],
            allow_zero_size_output_files => 1,
            );
    }
}

1;
