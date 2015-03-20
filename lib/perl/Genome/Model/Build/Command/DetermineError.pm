package Genome::Model::Build::Command::DetermineError;

use warnings;
use strict;

use Genome;
use IPC::System::Simple qw(capture);
use File::ReadBackwards;
use Try::Tiny;

# This is normally loaded automagically when Perl sees you're using %+.
# but due to some interaction with Class::Autouse (when it installs its
# UNIVERSAL::AUTOLOAD handler), this magic stops working.  Using the
# module explicitly fixes the problem in perl 5.10.1.  The bug seems to
# be fixed by 5.18
use Tie::Hash::NamedCapture;

our $WORKFLOW_DATE_AND_HOST = qr{
        (?<date>\d{4}-\d{2}-\d{2}\s\d{2}:\d{2}:\d{2})-\d{4}  # workflow date
        \s                                                   # and
        (?<host>[^\:]+)\:                                    # host
    }x;
our $PTERO_DATE = qr{
        \[(?<date>\d{4}/\d{2}/\d{2}\s\d{2}\:\d{2}\:\d{2}).\d+\]
    }x;
our $FLOW_DATE = qr{
        (?<date>\d{4}-\d{2}-\d{2}\s\d{2}\:\d{2}\:\d{2}).\d+
    }x;

our $ERROR_LOCATION = qr{
        at \s (?<error_source_file>\S+\.pm) \s line \s (?<error_source_line>\d+)
    }x;

our $ERROR_FINDING_REGEX = qr{
                (?:
                    $WORKFLOW_DATE_AND_HOST
                    | $PTERO_DATE
                    | $FLOW_DATE
                )
                \s
                (?:
                    (?:(?<error_text>ERROR:? .*?) \s $ERROR_LOCATION)
                    |
                    (?<error_text>ERROR:? .*)
                )
            }x;

our $EXCEPTION_FINDING_REGEX = qr{
                (?:
                    $WORKFLOW_DATE_AND_HOST
                    | $PTERO_DATE
                    | $FLOW_DATE
                )
                \s
                (?<error_text>.*) \s (?<!called\s) $ERROR_LOCATION
        }x;

our $PTERO_HOST_FINDING_REGEX = qr{Starting log annotation on host:\s(.*)};

class Genome::Model::Build::Command::DetermineError {
    is => 'Genome::Command::WithColor',
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
            } else {
                #if we still don't know, try to find the last message printed by die() or warn()
                ($error_source_file, $error_source_line, $error_host, $error_date, $error_text) = find_die_or_warn_in_log($error_log);
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

    my($get_one_more_line, $found_host);
    no warnings 'exiting';
    SCAN_FILE:
    for(1) {
        Genome::Sys->iterate_file_lines(
            $filename,
            $PTERO_HOST_FINDING_REGEX,
                sub { $found_host = $1 },

            $ERROR_FINDING_REGEX,
                sub {
                    ($error_date, $error_text, $error_source_file, $error_source_line)
                        = @+{'date','error_text','error_source_file','error_source_line'};
                    $error_host = $found_host || $+{'host'};

                    last SCAN_FILE if ($error_source_file);
                    $get_one_more_line = 1;
                },

            sub {
                if ($get_one_more_line) {
                    # sometimes we die immediately after the error
                    ($error_source_file, $error_source_line) = shift =~ m/at (\S+\.pm) line (\d+)/;
                    last SCAN_FILE;
                }
            }
        );
    }

    return unless $error_date;
    return ($error_source_file, $error_source_line, $error_host, $error_date, $error_text);
}

sub parse_output_log {
    my $filename = shift;

    my ($error_source_file, $error_source_line, $error_host, $error_date, $error_text);

    $error_source_file = 'n/a';
    $error_source_line = 'n/a';

    no warnings 'exiting';
    SCAN_FILE:
    for (1) {
        Genome::Sys->iterate_file_lines(
                $filename,
                qr{Results reported at (.*)$/},
                    sub { $error_date = $1 },

                qr{Job was executed on host\(s\) <([^>]+)>},
                    sub { $error_host = $1 },

                qr{(TERM_\w+:.+$)},
                    sub { $error_text = $1; last SCAN_FILE }, # This appears third in the LSF output, so we've located all three fields.
        );
    }

    return unless $error_text;
    return ($error_source_file, $error_source_line, $error_host, $error_date, $error_text);
}

sub find_die_or_warn_in_log {
    my $filename = shift;

    my ($error_source_file, $error_source_line, $error_host, $error_date, $error_text);

    my $backwards_fh = File::ReadBackwards->new($filename);

    no warnings 'exiting';
    SCAN_FILE:
    for(1) {
        Genome::Sys->iterate_file_lines(
            $backwards_fh,
            $PTERO_HOST_FINDING_REGEX,
                sub {
                    $error_host = $1;
                    last SCAN_FILE if $error_text;
                },

            $EXCEPTION_FINDING_REGEX,
                sub {
                    unless($error_text) {
                        ($error_date, $error_text, $error_source_file, $error_source_line, $error_host)
                            = @+{'date','error_text','error_source_file','error_source_line','host'};
                        last SCAN_FILE if ($error_host);
                    }
                },
        );
    }

    $backwards_fh->close;

    return ($error_source_file, $error_source_line, $error_host, $error_date, $error_text);
}

sub failed_workflow_steps {
    my $workflow = shift;

    my $failed_steps = [];
    workflow_visitor($workflow, $failed_steps);

    return sort {
        normalize_parent_id($a->parent_instance_id) <=> normalize_parent_id($b->parent_instance_id) ||
        step_time($a) cmp step_time($b)
    } @$failed_steps;
}

sub step_time {
    my $step = shift;
    return $step->end_time if $step->end_time;
    return 's' . $step->start_time if $step->start_time;
    return 'Unknown';
}

sub normalize_parent_id {
    my $parent_id = shift;
    return $parent_id ? -1 : 1;
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
