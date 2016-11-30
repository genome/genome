package Genome::Model::Tools::Sambamba::Base;

use strict;
use warnings;

use Genome;
use Carp qw(confess);

class Genome::Model::Tools::Sambamba::Base {
    is => 'Genome::Model::Tools::Base',
    is_abstract => 1,
};

sub versions {
    my $class = shift;
    my %versions = (
        '0.6.5' => '/gscmnt/gc13001/info/tools/sambamba_v0.6.5/sambamba_v0.6.5'
    );
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
