package Genome::Model::Build::Command::View;

use strict;
use warnings;

use Genome;

class Genome::Model::Build::Command::View {
    is => [
        'Genome::Command::Viewer',
        'Genome::Command::WorkflowMixin',
    ],
    has => [
        build => {
            is => 'Genome::Model::Build',
            shell_args_position => 1,
            doc => 'Genome::Model::Build',
        },
        events => {
            is => 'Boolean',
            is_optional => 1,
            default_value => 0,
            doc => 'Display build events.'
        },
        full => {
            is => 'Boolean',
            is_optional => 1,
            default_value => 0,
            doc => 'Display all build information. Equivalent to "--events --inputs --workflow --notes".'
        },
        inputs => {
            is => 'Boolean',
            is_optional => 1,
            default_value => 0,
            doc => 'Display build inputs.'
        },
        notes => {
            is => 'Boolean',
            is_optional => 1,
            default_value => 0,
            doc => 'Display build notes.'
        },
        workflow => {
            is => 'Boolean',
            is_optional => 1,
            default_value => 1,
            doc => 'Display workflow.',
        },
    ],
};

sub help_synopsis {
    my $self = shift;
    my $result .= <<EOP;
    Displays basic information about a build.
EOP
    return $result;
}

sub help_detail {
    my $self = shift;
    my $result .= <<EOP;
Displays information about the events, inputs, and workflow for a build.

To display a full report about a build try:

    genome model build view --full <build_id>

This will display almost all the information about the build, excluding input and output connectors.
To also include input and ouput connectors:

    genome model build view --full --connectors <build_id>
EOP
    return $result;
}

sub write_report {
    my ($self, $width, $handle) = @_;

    my $build = $self->build;
    $self->_display_build($handle, $build);

    for my $thing (["inputs", "Inputs", "_display_input"],
                   ["events", "Events", "_display_event"],
                   ["notes", "Notes", "_display_note"]) {
        my ($item, $section_name, $method_name) = @{$thing};
        if($self->full || $self->$item) {
            my @items = $build->$item;
            $self->_display_many($handle, $section_name,
                    $method_name, @items);
        }
    }

    if($self->full || $self->workflow) {
        my $workflow = $self->build->newest_workflow_instance;
        $self->_display_workflow($handle, $workflow);
    }

    1;
}

sub _display_build {
    my ($self, $handle, $build) = @_;

    my $processing_profile = $build->model->processing_profile;

    my $format_str = <<EOS;
%s
%s %s
%s %s
%s %s
%s %s

%s
%s%s
%s

EOS

    print $handle sprintf($format_str,
        map {$self->_pad_right($_, $self->COLUMN_WIDTH)} (
            $self->_color_heading('Build'),
            $self->_color_pair('Build ID', $build->id),
            $self->_color_pair('Build Status',
                $self->_status_color($build->status)),
            $self->_color_pair('Model ID', $build->model->id),
            $self->_color_pair('Model Name', $build->model->name),
            $self->_color_pair('Run by', $build->run_by),
            $self->_color_pair('Processing Profile ID', $processing_profile->id),
            $self->_color_pair('Build Scheduled',
                    $self->_clean_up_timestamp($build->date_scheduled)),
            $self->_color_pair('Build Completed',
                $self->_clean_up_timestamp($build->date_completed)),
        ),
        $self->_color_pair('Build Class', $build->class),
        $self->_color_pair('Software Revision', $build->software_revision),
        $self->_software_result_test_name,
        $self->_color_pair('Data Directory', $build->data_directory));
}

sub _software_result_test_name {
    my ($self) = @_;
    my $build = $self->build;
    my @results = $build->all_results;
    my $label = 'SoftwareResult Test Name(s)';

    unless (scalar(@results)) {
        return "\n" . $self->_color_pair($label, "No results found for this build");
    }
    return '' unless scalar(@results);

    my $UNDEF = '***undef***';
    my @test_names = map {$_->test_name ? $_->test_name : $UNDEF} @results;

    my %counts;
    $counts{$_}++ for @test_names;
    my @sorted_counts = (
        sort { $b->[1] <=> $a->[1] }
        map { [ $_, $counts{$_} ] } keys %counts
    );

    my @entries;
    while (my $entry = shift(@sorted_counts)) {
        my ($name, $count) = @{$entry};
        if ($name eq $UNDEF) {
            $name = "undef";
        } else {
            $name = '"' . $name . '"';
        }

        push @entries, sprintf("%s (%d)", $name, $count);
    }
    my $value = join(", ", @entries);

    return "\n" . $self->_color_pair($label, $value);
}

sub _display_many {
    my ($self, $handle, $section_name, $method_name, @items) = @_;

    print $handle $self->_color_heading($section_name) . "\n";
    if(@items) {
        for my $item (@items) {
            $self->$method_name($handle, $item);
        }
    } else {
        print $handle "None\n";
    }
    print $handle "\n";
}

sub _display_note {
    my ($self, $handle, $note) = @_;

    print $handle $self->_color_pair("Editor", $note->editor_id);
    print $handle $self->_color_pair("    Date/Time", $note->entry_date) . "\n";
    print $handle $self->_color_pair("    Header", $note->header_text) . "\n";
    print $handle $self->_color_pair("    Body", $note->body_text) . "\n";
}

sub _display_input {
    my ($self, $handle, $input) = @_;

    print $handle $input->name . "\n";
    print $handle $self->_color_pair('    ID', $input->value_id) . "\n";
    print $handle $self->_color_pair('    Type', $input->value_class_name) . "\n";
}

sub _display_event {
    my ($self, $handle, $event) = @_;

    my $format_str = <<EOS;
%s %s
    %s
    %s   %s
    %s
EOS
    print $handle sprintf($format_str,
        $self->_color_pair('ID', $self->_pad_right($event->id, $self->COLUMN_WIDTH)),
        $self->_color_pair('Status',
            $self->_status_color($event->event_status)),
        $self->_color_pair('Name', $event->event_type),
        $self->_color_pair('Scheduled',
            $self->_clean_up_timestamp($event->date_scheduled)),
        $self->_color_pair('Completed',
            $self->_clean_up_timestamp($event->date_completed)),
        $self->_color_pair('LSF ID', $event->lsf_job_id));
}

sub _resolve_running_child_end_time {
    my ($self) = @_;
    unless ('Abandoned' eq $self->build->status) {
        return DateTime->now(time_zone=>'local');
    }
    return;
}

1;
