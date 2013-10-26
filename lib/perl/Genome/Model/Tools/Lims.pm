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
        $ENV{'LIMS_PERL'} = 1;
        return $class->SUPER::resolve_class_and_params_for_argv(@_);
    }
}

sub execute {
    my $self = shift;
    my $exit_code = system('/gsc/bin/perl', '-S', 'gmt', 'lims', $self->args);;
    return ($exit_code == 0);
}

1;
