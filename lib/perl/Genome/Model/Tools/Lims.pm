package Genome::Model::Tools::Lims;
use strict;
use warnings;
use Genome;

class Genome::Model::Tools::Lims {
    is => 'Command::Tree',
    has_many_optional_input => [
        args => { is => 'Text', is_many => 1, shell_args_position => 1, },
    ],
    doc => 'all tools which interface directly with TGI LIMS should go here'
};

sub sub_command_sort_position { -1 }

sub resolve_class_and_params_for_argv {
    my $class = shift;
    if ($^X eq '/usr/bin/perl') {
        #warn "delegate to LIMS interpreter: running LIMS mod on $class with @_\n";
        return ($class, { args => [@_] });
    }
    else {
        #warn "already on the LIMS interpreter: running normally\n";
        return $class->SUPER::resolve_class_and_params_for_argv(@_);
    }
}

sub execute {
    my $self = shift;
    my @args = $self->args;
    my $cmd = "/gsc/bin/perl `which gmt` lims @args";
    #warn "on non-LIMS interpreter, shelling out to run this: $cmd\n";
    my $exit_code = system $cmd;
    $exit_code /= 256;
    return !$exit_code;
}

1;

