package Genome::Model::Tools::Manta::Base;

use strict;
use warnings;

use Genome;
use Carp qw(confess);

class Genome::Model::Tools::Manta::Base {
    is => 'Genome::Model::Tools::Base',
    is_abstract => 1,
};

sub help_detail {
    "Manta calls structural variants (SVs) and indels from mapped paired-end sequencing reads."
}

sub versions {
    my $class = shift;
    my %VERSIONS = (
        '0.29.6'     => '/gscmnt/gc13001/info/model_data/jwalker_scratch/src/manta-0.29.6.centos5_x86_64/bin',
    );
    return %VERSIONS;
}

sub tool_path {
    my $self = shift;
    return File::Spec->join($self->path_for_version($self->version),$self->_tool_subcommand_name);
}

# tool subcommand name as a string (e.g., 'configManta.py')
sub _tool_subcommand_name {
    confess "sub _tool_subcommand_name must be overridden in child class";
}


1;
