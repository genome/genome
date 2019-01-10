package Genome::Model::Tools::Sambamba::Base;

use strict;
use warnings;

use Genome;
use Carp qw(confess);
use File::Spec;

class Genome::Model::Tools::Sambamba::Base {
    is => 'Genome::Model::Tools::Base',
    is_abstract => 1,
};

sub versions {
    my $class = shift;

    my @versions = ('0.6.5');

    my %versions;
    for my $version (@versions) {
        my $alloc = Genome::Disk::Allocation->get(allocation_path => 'tools/sambamba_v'.$version);
        unless ($alloc) {
            $class->fatal_message('failed to locate expected version %s', $version);
        }

        $versions{$version} = File::Spec->join($alloc->absolute_path, 'sambamba_v'.$version);
    }

    return %versions;
}

sub tool_path {
    my $self = shift;
    return ($self->path_for_version($self->version),$self->_tool_subcommand_name);
}

# tool subcommand name as a string (e.g., 'view')
sub _tool_subcommand_name {
    confess ('sub _tool_subcommand_name must be overridden in child class');
}


1;
