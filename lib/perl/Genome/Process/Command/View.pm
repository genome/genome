package Genome::Process::Command::View;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::Process::Command::View {
    is => [
        'Genome::Command::Viewer',
        'Genome::Command::WorkflowMixin',
    ],
    has => [
        process => {
            is => 'Genome::Process',
            shell_args_position => 1,
            doc => 'Genome::Process',
        },
        status_events => {
            is => 'Boolean',
            is_optional => 1,
            default_value => 0,
            doc => 'Display process status-events.'
        },
        inputs => {
            is => 'Boolean',
            is_optional => 1,
            default_value => 0,
            doc => 'Display process inputs.'
        },
        full => {
            is => 'Boolean',
            is_optional => 1,
            default_value => 0,
            doc => 'Display all process information. Equivalent to "--events --inputs --workflow".'
        },
    ],
    doc => 'Displays basic information about a process.',
};

sub help_detail {
    my $self = shift;
    return <<EOP;
    Displays basic information about a process.
EOP
}

sub help_synopsis {
    my $self = shift;
    return <<EOP;
Displays information about the events, inputs, and workflow for a process.

To display a full report about a process try:

    genome process view --full <process_id>
EOP
}

sub write_report {
    my ($self, $width, $handle) = @_;

    $self->_display_process($handle);

    if ($self->full) {
        $self->inputs(1);
        $self->status_events(1);
        $self->workflow(1);
    }

    for my $thing (["inputs", "Inputs", "_display_input"],
                   ["status_events", "Events", "_display_event"]) {
        my ($item, $section_name, $method_name) = @{$thing};
        if($self->$item) {
            my @items = $self->process->$item;
            $self->_display_many($handle, $section_name,
                    $method_name, @items);
        }
    }

    my $workflow = $self->process->newest_workflow_instance;
    $self->_display_workflow($handle, $workflow);
    $self->_display_logs($handle, $workflow);

    1;
}

sub _display_process {
    my ($self, $handle) = @_;
    my $process = $self->process;

    my $format_str = <<EOS;
%s
%s %s
%s %s
%s

%s
%s
%s
%s

EOS

    print $handle sprintf($format_str,
        $self->_color_heading('Process'),
        map {$self->_pad_right($_, $self->COLUMN_WIDTH)} (
            $self->_color_pair('Process ID', $process->id),
            $self->_color_pair('Process Status',
                $self->_status_color($process->status)),
            $self->_color_pair('Process Started',
                    $self->_clean_up_timestamp($process->started_at)),
            $self->_color_pair('Process Ended',
                $self->_clean_up_timestamp($process->ended_at)),
            $self->_color_pair('Run by', $process->created_by),
        ),
        $self->_color_pair('Process Class', $process->class),
        $self->_color_pair('Software Revision', $process->software_revision),
        $self->_color_pair('SoftwareResult Test Name(s)',
            $self->_software_result_test_names($process->unique_results)),
        $self->_color_pair('MetaData Directory', $process->metadata_directory));
}

sub _display_input {
    my ($self, $handle, $input) = @_;

    print $handle $input->name . "\n";
    print $handle $self->_color_pair('    ID', $input->value_id) . "\n";
    print $handle $self->_color_pair('    Type', $input->value_class_name) . "\n";
    if ($input->array_index != -1) {
        print $handle $self->_color_pair('    Index', $input->array_index) . "\n";
    }
}

sub _display_event {
    my ($self, $handle, $event) = @_;

    my $format_str = <<EOS;
%s
    %s
    %s
EOS
    print $handle sprintf($format_str,
        $event->timestamp,
        $self->_color_pair('Old Status', $event->old_status),
        $self->_color_pair('New Status',
            $self->_clean_up_timestamp($event->new_status)),
    );
}


1;
