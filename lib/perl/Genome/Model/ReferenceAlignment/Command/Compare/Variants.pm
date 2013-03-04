package Genome::Model::ReferenceAlignment::Command::Compare::Variants;

use strict;
use warnings;

use Genome;

class Genome::Model::ReferenceAlignment::Command::Compare::Variants {
    #is => 'Genome::Model::Compare::Base',
    is => 'Command::V2',
    has_input => [
        from_build => {
            is => 'Genome::Model::Build',
            shell_args_position => 0,
        },
        to_build => {
            is => 'Genome::Model::Build',
            shell_args_position => 1,
        },
        output_dir => {
            is => 'FilesystemDirectory',
        },
    ],
};

sub execute {
    my $self = shift;
    $self->status_message('Please compare variants from build '. $self->from_build->__display_name__ .' to build '. $self->to_build->__display_name__ .'\!');
    Genome::Sys->shellcmd(cmd => 'touch '. $self->output_dir .'/variants.diff');
    return 1;
}

1;
