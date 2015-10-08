package Genome::Model::ReferenceVariation::Command::Base;

use strict;
use warnings;

use Genome;

class Genome::Model::ReferenceVariation::Command::Base {
    is => 'Command::V2',
    has_input => [
        build => {
            is => 'Genome::Model::Build::ReferenceVariation',
            doc => 'build for which to run the command',
            is_output => 1,
            shell_args_position => 1,
        },
    ],
    doc => 'Base class for reference variation pipeline commands.',
    is_abstract => 1,
};

sub shortcut {
    my $self = shift;

    my $cmd = $self->_command;
    my $retval = $cmd->shortcut;

    my $accessor = $self->_result_accessor;

    if ($retval and my $result = $cmd->$accessor) {
        $self->status_message('Found existing result: %s', $result->__display_name__);
        $self->_postprocess_result($result);
    }

    return $retval;
}

sub execute {
    my $self = shift;

    my $cmd = $self->_command;
    my $retval = $cmd->execute;

    my $accessor = $self->_result_accessor;

    if ($retval and my $result = $cmd->$accessor) {
        $self->status_message('Generated result: %s', $result->__display_name__);
        $self->_postprocess_result($result);
    } else {
        die $self->error_message('Failed to produce result.');
    }

    return $retval;
}

sub _command {
    my $self = shift;
    my $command_class = $self->_command_class;
    my @params = $self->_params_for_command;

    my $cmd = $command_class->create(@params);
    unless ($cmd) {
        die $self->error_message('Failed to create command.');
    }

    return $cmd;
}

sub _result_accessor {
    Carp::confess('subclass must implement _result_accessor');
}

sub _command_class {
    Carp::confess('subclass must implement _command_class');
}

sub _params_for_command {
    Carp::confess('subclass must implement _params_for_command');
}

sub _label_for_result {
    return;
}

sub _postprocess_result {
    my $self = shift;
    my $result = shift;

    return unless $self->_label_for_result;

    $result->add_user(user => $self->build, label => $self->_label_for_result);
}

sub sub_command_category {
    return 'pipeline steps';
}

1;
