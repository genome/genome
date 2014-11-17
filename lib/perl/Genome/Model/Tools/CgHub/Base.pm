package Genome::Model::Tools::CgHub::Base;

use strict;
use warnings;

use Genome;

require File::Temp;
require File::Spec;

class Genome::Model::Tools::CgHub::Base {
    is => "Command::V2",
    is_abstract => 1,
    has_optional_output => {
        output_file => {
            is => 'Text', 
            doc => 'Save standard output/error to this file.',
        },
    },
};

sub execute {
    my $self = shift;

    my $cmd = $self->_build_command;
    die $self->error_message('Failed to build CG Hub command!') if not $cmd;

    my $run_ok = $self->_run_command($cmd);
    die $self->error_message('Failed to run CG Hub command!') if not $run_ok;

    my $success_ok = $self->_verify_success;
    die $self->error_message('Ran CG Hub command, but it was determined that it was not successful!') if not $success_ok;

    return 1;
}

sub _run_command {
    my ($self, $cmd) = @_;

    my $output_file = $self->output_file;
    if ( not $output_file ) {
        my $tmp_dir = File::Temp::tempdir(CLEANUP => 1);
        $output_file = $self->output_file(
            File::Spec->catfile($tmp_dir, 'cghub.out')
        ); 
    }
    $cmd .= " | tee $output_file";

    local $ENV{'PATH'} = $ENV{'PATH'} . ':/cghub/bin';
    my $rv = eval{ Genome::Sys->shellcmd(cmd => $cmd); };
    if ( not $rv ) {
        $self->error_message($@) if $@;
        $self->error_message('Failed to execute CG Hub command!');
        return;
    }

    return 1;
}

1;

