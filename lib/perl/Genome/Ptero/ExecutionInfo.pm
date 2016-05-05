package Genome::Ptero::ExecutionInfo;

use strict;
use warnings;

use Ptero::HTTP qw();
use Ptero::Proxy::Workflow::Execution;

use Scalar::Util;

use Genome;

use Genome::Utility::Text;
use Genome::Utility::Inputs;

class Genome::Ptero::ExecutionInfo {
    is => 'Command::V2',
    roles => ['Genome::Role::CommandWithColor'],
    has_input => [
        execution_id => {
            is => 'Text',
            shell_args_position => 1,
            doc => 'ID of the execution for which to retrieve information',
        },
        display_inputs => {
            is => 'Boolean',
            default_value => 0,
            doc => 'Display information about the inputs to this step',
        },
        display_outputs => {
            is => 'Boolean',
            default_value => 1,
            doc => 'Display the standard output and error from the execution',
        }
    ],
    has_transient_optional => {
        execution => {
            is => 'Ptero::Concrete::Workflow::Execution',
            doc => 'The PTero Perl SDK entity representing an execution',
        },
        job_info => {
            is => 'HASH',
            doc => 'Information about the job for the execution',
        },
    },
    doc => 'retrieve information about a PTero execution',
};

sub execute {
    my $self = shift;

    $self->_display_basic_info;
    $self->_display_inputs if $self->display_inputs;
    $self->_display_outputs if $self->display_outputs;

    return 1;
}

sub job_info {
    my $self = shift;

    unless($self->__job_info) {
        my $execution = $self->execution;

        my $url = $execution->{data}{jobUrl};
        unless ($url) {
            $self->fatal_message('No job information available for execution %s.', $execution->{id});
        }

        $self->debug_message('Retrieving data from <%s>.', $url);

        my $data = Ptero::HTTP::make_request_and_decode_response(method => 'GET', url => $url);
        unless ($data and $data->{jobId}) {
            $self->fatal_message('No job information found at <%s>.', $url);
        }

        $self->__job_info($data);
    }

    return $self->__job_info;
}

sub execution {
    my $self = shift;

    unless($self->__execution) {
        my $url = Genome::Config::get('ptero_workflow_submit_url');
        $url =~ s/workflows$/executions/;
        $url .= '/' . $self->execution_id;

        my $proxy = Ptero::Proxy::Workflow::Execution->new(url => $url);
        unless($proxy) {
            $self->fatal_message('No execution found for id %s.', $self->execution_id);
        }

        $self->__execution($proxy->concrete_execution);
    }

    return $self->__execution;
}

sub _display_basic_info {
    my $self = shift;
    my $job_info = $self->job_info;

    for my $key ('jobId', 'status', 'lsfJobId', 'commandLine', 'command') {
        my $value = $job_info->{$key};
        next unless $value;

        if (ref($value) eq 'ARRAY') { $value = join(' ', @$value); }

        my $padded_key = Genome::Utility::Text::justify(uc $key, 'right', 15);
        print $self->_color_pair($padded_key, $value), "\n";
    }
}

sub _display_inputs {
    my $self = shift;
    my $execution = $self->execution;

    print $self->_color_heading('INPUTS'), "\n";

    my $inputs = Genome::Utility::Inputs::decode($execution->{inputs});

    for my $key (sort keys %$inputs) {
        my $padded_key = '  ' . $key;
        my $value = $self->_display_value_for_input($inputs->{$key});
        print $self->_color_pair( $padded_key, $value ), "\n";
    }
}

sub _display_value_for_input {
    my $self = shift;
    my $input = shift;

    if (Scalar::Util::blessed($input)) {
        if ($input->can('__display_name__')) {
            return $input->__display_name__;
        } else {
            return ref($input);
        }
    } elsif(ref $input eq 'ARRAY') {
        return join('; ', map $self->_display_value_for_input($_), @$input);
    } elsif(ref $input eq 'HASH') {
        return join('; ',
            map
                join(': ', $_, $self->_display_value_for_input($input->{$_})),
                sort keys %$input
        );
    } else {
        return $input;
    }
}

sub _display_outputs {
    my $self = shift;
    my $job_info = $self->job_info;

    for my $output ('stdout', 'stderr') {
        if (defined $job_info->{$output}) {
            print
                $self->_color_heading(uc $output),"\n",
                $job_info->{$output},
                "\n";
        }
    }
}

1;
