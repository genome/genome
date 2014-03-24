package Genome::Model::Build::Command::DetermineError;

use warnings;
use strict;

use Genome;
use File::Slurp 'read_file';
use IPC::System::Simple qw(capture);
use Try::Tiny;

class Genome::Model::Build::Command::DetermineError {
    is => 'Genome::Command::Base',
    has => [
        build => {
            is => 'Genome::Model::Build',
            shell_args_position => 1,
            doc => 'The build',
        },
        display_results => {
            is => 'Boolean',
            default => 1,
        },
        assume_build_status => {
            is => 'Text',
            is_optional => 1,
            doc => 'Treat the build as if it has this status when determining the error. (Useful for finding errors in already-abandoned builds.)',
            valid_values => ['Failed', 'Unstartable'],
        },
    ],
    has_optional_output => [
        error_type => {
            is => "Text",
        },
        error_source_file => {
            is => "Text",
            default => "Unknown",
        },
        error_source_line => {
            is => "Text",
            default => "Unknown",
        },
        error_log => {
            is => "Text",
            default => "Unknown",
        },
        error_text => {
            default => "Unknown",
            is => "Text",
        },
        error_date => {
            default => "Unknown",
            is => "Text",
        },
        error_host => {
            default => "Unknown",
            is => "Text",
        },
    ],
};

sub help_detail {
    return <<HELP;
This command determines information about a Failed or Unstartable build.
HELP
}

sub execute {
    my $self = shift;

    my $status = $self->assume_build_status || $self->build->status;

    if ($status eq 'Failed') {
        $self->handle_failed();
    } elsif ($status eq 'Unstartable') {
        $self->handle_unstartable();
    } else {
        $self->error_message(sprintf("Build (%s) has an unexpected status: %s",
                $self->build->id, $self->build->status));
        die $self->error_message;
    }

    if ($self->display_results) {
        $self->print_results;
    }

    return 1;
}

sub print_results {
    my $self = shift;

    print $self->formatted_results;
}

sub formatted_results {
    my $self = shift;

    my $results = join('',
        "\n",
        sprintf("Build: %s\t\t", $self->_color($self->build->id, 'bold')),
        sprintf("Error Type: %s\t\t", $self->_color($self->error_type, 'bold')),
        sprintf("Error Date: %s\n", $self->_color($self->error_date, 'bold')),
        sprintf("Error Source File: %s\n", $self->_color($self->error_source_file, 'bold')),
        sprintf("Error Source Line: %s\n", $self->_color($self->error_source_line, 'bold')),
        sprintf("Error Log: %s\n", $self->_color($self->error_log, 'bold')),
        sprintf("Error Host: %s\n", $self->_color($self->error_host, 'bold')),
        sprintf("Error Message: %s\n", $self->_color($self->error_text, 'bold')),
    );

    return $results;
}

sub handle_failed {
    my $self = shift;

    $self->error_type("Failed");

    try {
        my $workflow = $self->build->newest_workflow_instance;
        $self->handle_failed_from_logs($workflow);
    } catch {
        my $error_log = find_error_log($self->build->log_directory);
        if (defined $error_log)  {
            $self->set_status_from_log_file($error_log);
        }
    };
}

sub find_error_log {
    my ($log_dir) = @_;
    try {
        my @logs = capture(qw(grep -l ERROR), $log_dir);
        chomp @logs;
        return $logs[0];
    } catch {
        return;
    };
}

sub handle_failed_from_logs {
    my ($self, $workflow) = @_;

    my @failed_steps = failed_workflow_steps($workflow);
    for my $failed_step (@failed_steps) {
        if ($failed_step->current->can('stderr') and my $error_log = $failed_step->current->stderr) {
            if (-e $error_log and -s $error_log) {

                $self->set_status_from_log_file($error_log);

                # if we can't determine the date from the log file, try the workflow step.
                unless ($self->error_date) {
                    if ($failed_step->end_time) {
                        $self->error_date($failed_step->end_time);
                    } elsif ($failed_step->start_time) {
                        $self->error_date("Sometime After " . $failed_step->start_time);
                    }
                }
                return;
            }
        }
    }
}

sub set_status_from_log_file {
    my $self = shift;
    my ($error_log) = @_;

    my ($error_source_file, $error_source_line, $error_host, $error_date, $error_text) = parse_error_log($error_log);
    unless($error_source_file || $error_source_line || $error_host || $error_date || $error_text) {
        my $output_log = $error_log;
        $output_log =~ s/\.err$/.out/;

        if(-e $output_log and -s _) {
            #Look for LSF errors
            ($error_source_file, $error_source_line, $error_host, $error_date, $error_text) = parse_output_log($output_log);
            if($error_text) {
                $error_log = $output_log; #only prefer .out files if we got a result
            }
        }
    }


    $self->error_log($error_log);

    $self->error_source_file($error_source_file) if $error_source_file;
    $self->error_source_line($error_source_line) if $error_source_line;
    $self->error_text($error_text) if $error_text;
    $self->error_host($error_host) if $error_host;
    $self->error_date($error_date) if $error_date;
}

sub parse_error_log {
    my $filename = shift;

    my ($error_source_file, $error_source_line, $error_host, $error_date, $error_text);

    my $file_text = read_file($filename);

    ($error_date, $error_host, my $date_removed_text) = get_error_date($file_text);
    if ($error_date) {
        my $text = '(ERROR[\s\:]+.{1,400}?)';
        my $file = 'at\s([^\s]*?\.pm)';
        my $line = 'line\s(\d+)';
        my $query = join('\s+', $text, $file, $line);
        if ($date_removed_text =~ m/$query/s) {
            $error_text = $1;
            $error_source_file = $2;
            $error_source_line = $3;
        } elsif ($date_removed_text =~ m/(ERROR.*)/) {
            $error_text = $1;
        }
    }

    return ($error_source_file, $error_source_line, $error_host, $error_date, $error_text);
}

sub parse_output_log {
    my $filename = shift;

    my ($error_source_file, $error_source_line, $error_host, $error_date, $error_text);
    my $file_text = Genome::Sys->read_file($filename);

    if($file_text =~ /(TERM_\w+:.+$)/m) {
        $error_text = $1;

        if($file_text =~ /Results reported at (.*)$/m) {
            $error_date = $1;
        }

        if($file_text =~ /Job was executed on host\(s\) <([^>]+)>/) {
            $error_host = $1;
        }

        $error_source_file = 'n/a';
        $error_source_line = 'n/a';
    }

    return ($error_source_file, $error_source_line, $error_host, $error_date, $error_text);
}

sub get_error_date {
    my $file_text = shift;

    my $date = '(\d{4}-\d{2}-\d{2}\s\d{2}:\d{2}:\d{2})-\d{4}';
    my $host = '([^\:]+)\:';

    my ($error_date, $error_host, $formatted_text);
    if ($file_text =~ m/$date\s$host\s+ERROR/) {
        $error_date = $1;
        $error_host = $2;
        ($formatted_text = $file_text) =~ s/$date\s+$host\s//g;
    } else {
        ($error_date, $error_host, $formatted_text) = get_error_date_from_ptero_log($file_text);
    }

    return $error_date, $error_host, $formatted_text;
}

sub get_error_date_from_ptero_log {
    my $file_text = shift;

    my ($error_date, $error_host, $formatted_text);

    my $date = '\[(\d{4}/\d{2}/\d{2}\s\d{2}\:\d{2}\:\d{2}).\d+\]';
    if ($file_text =~ m/$date\s+ERROR/) {
        $error_date = $1;

        # remove everything after the error and search backwards through the file
        # for the host
        (my $error_truncated_text = $file_text) =~ s/ERROR.*//s;

        my @lines = split(/\n/, $error_truncated_text);
        for my $line (reverse @lines) {
            if ($line =~ m/Starting log annotation on host:\s(.*)/) {
                $error_host = $1;
                last;
            }
        }

        ($formatted_text = $file_text) =~ s/$date\s//g;
    }

    return $error_date, $error_host, $formatted_text;
}

sub failed_workflow_steps {
    my $workflow = shift;

    my $failed_steps = [];
    workflow_visitor($workflow, $failed_steps);

    return sort {step_time($a) cmp step_time($b)} @$failed_steps;
}

sub step_time {
    my $step = shift;
    return $step->end_time if $step->end_time;
    return $step->start_time if $step->start_time;
    return 'Unknown';
}

sub workflow_visitor {
    my ($step, $failed_steps) = @_;

    my $status = $step->status;

    if($status eq 'failed' or $status eq 'crashed' or
       $status eq 'running' or $status eq 'running*') {
        push(@{$failed_steps}, $step);
    }

    if ($step->can('related_instances')) {
        for my $sub_step ($step->related_instances) {
            workflow_visitor($sub_step, $failed_steps);
        }
    }
}

sub set_if_defined {
    my ($self, $field, $error, $accessor_base) = @_;

    # prefer $accessor_base then $inferred_accessor
    my $inferred_accessor = "inferred_$accessor_base";
    my $value = $error->$accessor_base;
    unless ($value) {
        $value = $error->$inferred_accessor;
    }

    if ($value) {
        $self->$field($value);
    }
}

sub handle_unstartable {
    my $self = shift;

    $self->error_type("Unstartable");

    my $note = $self->error_note;
    if ($note) {
        $self->error_date($note->entry_date);
        $self->error_text($note->body_text);
        $self->set_file_and_line($note->body_text);
    } else {
        $self->warning_message(sprintf('Could not find any notes on build (%s) which seem related to errors',
                $self->build->id));
    }
}

sub error_note {
    my $self = shift;

    my @notes = grep {$_->header_text !~ m/Build Created/} $self->build->notes;
    return unless (@notes);

    my @sorted_notes = sort {$a->entry_date cmp $b->entry_date} @notes;

    # prefer 'Unstartable" notes
    my @unstartable_notes = grep {$_->header_text eq 'Unstartable'} @sorted_notes;
    if (@unstartable_notes) {
        return shift @unstartable_notes;
    } else {
        return shift @sorted_notes;
    }
}

sub set_file_and_line {
    my ($self, $text) = @_;
    if ($text =~ m/\sat\s(.*\.pm)\sline\s(\d+)/) {
        $self->error_source_file($1);
        $self->error_source_line($2);
    }
}
