package Genome::Model::Tools::CompleteGenomics;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::CompleteGenomics {
    is => 'Command::Tree',
};

sub path_to_cgatools {
    my $class = shift;

    #if/when versioning support added, just change this here
    return '/gsc/bin/cgatools';
}

sub run_command {
    my $class = shift;
    my $subcommand = shift;
    my %shellcmd_opts = @_;

    Genome::Sys->shellcmd(
        cmd => join(' ', $class->path_to_cgatools, $subcommand),
        %shellcmd_opts,
    );

    return 1;
}

1;
