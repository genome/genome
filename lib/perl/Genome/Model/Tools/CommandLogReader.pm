package Genome::Model::Tools::CommandLogReader;

use strict;
use warnings;
use Genome;
use IO::File;
use DateTime;

class Genome::Model::Tools::CommandLogReader{
    is => 'Command',
    has =>[
        start_date => {
            is => 'Text',
            doc => "Logs written after this date (YYYY-MM-DD) will be included in the reader's output.",
        },
    ],
    has_optional =>[
        report_type => {
            is => 'Text',
            valid_values => ["full","abbreviated","summary"],
            default => "summary",
            doc => "A full report lists every user/command/parameter combination and counts number of times each combination occurs. An abbreviated report excludes the parameter column. A summary report includes how many commands each user ran (if group by is user, does not include commands) or how many times a command was executed (if group by is command, does not include user). Defaults to summary.",
        },
        output_file => {
            is => 'Text',
            doc => "Store log information in the specified file, defaults to STDOUT.",
            default => "STDOUT",
        },
        end_date => {
            is => 'Text',
            doc => "Logs written before this date (YYYY-MM-DD) will be included in the reader's output, defaults to current date.",
        },
        group_by => {
            is => 'Text',
            valid_values => ["user","command_class"],
            default => "user",
            doc => "Groups output by user name or command class, defaults to user.",
        },
        print_headers => {
            is => 'Boolean',
            default => 1,
            doc => "Print column headers on first line of output.",
        },
    ],
};

sub help_brief {
    "Generates a report on command line usage of genome and gmt scripts";
}

sub help_synopsis {
    return <<"EOS"
EOS
}

sub help_detail {
    return <<"EOS"
This module generates a report on command line usage of genome scripts, grouping by either user name or command
EOS
}

sub output_attributes {
    my $self = shift;
    if ($self->group_by eq "user") {
        if ($self->report_type eq "full") {
            return (qw/ user command_class command_with_args times_called/);
        }
        elsif ($self->report_type eq "abbreviated") {
            return (qw/ user command_class times_called/);
        }
        elsif ($self->report_type eq "summary") {
            return (qw/ user num_commands/);
        }
    }
    elsif ($self->group_by eq "command_class") {
        if ($self->report_type eq "full") {
            return (qw/ command_class user times_called command_with_args/);
        }
        elsif ($self->report_type eq "abbreviated") {
            return (qw/ command_class user times_called/);
        }
        elsif ($self->report_type eq "summary") {
            return (qw/ command_class times_called/);
        }
    }
    else {
        $self->error_message("Unknown group by attribute or report type, unable to provide column headers");
        return;
    }
}

sub _create_file {
    my ($self, $output_file) = @_;
    my $output_fh;

    if (-e $output_file) {
        $self->warning_message("found previous output file, removing $output_file");
        unlink $output_file;
        if (-e $output_file) {
            die "failed to remove previous file: $! ($output_file)";
        }
    }
    $output_fh = IO::File->new("> $output_file");
    unless ($output_fh) {
        die "Can't open file ($output_file) for writing: $!";
    }

    return $output_fh;
}

sub execute {
    my $self = shift;

    my $base_log_dir = Genome::Site::TGI::Security->base_log_dir();
    unless (-e $base_log_dir and -d $base_log_dir) {
        $self->error_message("Log directory $base_log_dir does not exist, exiting");
        return;
    }

    if ($self->start_date !~ /\d{4}-\d{1,2}-\d{1,2}/) {
        $self->error_message("Invalid start date format, use YYYY-MM-DD");
        return;
    }

    if (defined $self->end_date and $self->end_date !~ /\d{4}-\d{1,2}-\d{1,2}/) {
        $self->error_message("Invalid end date format, use YYYY-MM-DD");
        return;
    }

    my ($start_year, $start_month, $start_day) = split(/-/, $self->start_date);
    my $start_date = DateTime->new(year => $start_year, month => $start_month, day => $start_day);

    my $end_date = DateTime->now(time_zone => 'America/Chicago');
    if (defined $self->end_date) {
        my ($end_year, $end_month, $end_day) = split(/-/, $self->end_date);
        $end_date = DateTime->new(year => $end_year, month => $end_month, day => $end_day);
    }

    my %info;
    my @log_columns = Genome::Site::TGI::Security->log_columns();;
    until (DateTime->compare($start_date, $end_date) >= 0) {
        my $log_dir = $base_log_dir . "/" . $start_date->year . "/";
        my $log_file = $log_dir . $start_date->month . "-" . $start_date->day . ".log";

        if (-e $log_file) {
            my $log_svr = Genome::Utility::IO::SeparatedValueReader->create(
                input => $log_file,
                headers => \@log_columns,
                separator => "\t",
                is_regex => 1,
                ignore_extra_columns => 1,
            );

            while (my $line = $log_svr->next) {
                my $user = $line->{user} || '-';
                my $command_class = $line->{command_class} || '-';
                my $command_with_args = $line->{command_with_args} || '-';
                if ($self->group_by eq "command_class") {
                    $info{$command_class}->{$user}->{$command_with_args}++;
                }
                elsif ($self->group_by eq "user") {
                    $info{$user}->{$command_class}->{$command_with_args}++;
                }
            }
        } else{
            $self->warning_message("Could not find $log_file, continuing.");
        }

        $start_date->add(days => 1);
    }

    # Print to output file, tab delimited
    # Format (grouped by user): user command params(space delimited) times_called
    # Format (grouped by command): command user params(space delimited) times_called
    my $output_fh;
    my $output_file = $self->output_file;
    if ($self->output_file =~ /STDOUT/i) {
        $output_fh = 'STDOUT';
    }
    else {
        $output_fh = $self->_create_file($output_file);
    }

    $output_fh->print(join("\t", $self->output_attributes) . "\n") if $self->print_headers;

    for my $group (sort keys %info) {
        if ($self->report_type eq "summary") {
            my $sum = 0;
            for my $subgroup (sort keys %{$info{$group}}) {
                for my $params (keys %{$info{$group}->{$subgroup}}) {
                    $sum += $info{$group}->{$subgroup}->{$params};
                }
            }
            $output_fh->print($group . "\t");
            $output_fh->print($sum . "\n");
        }
        else {
            for my $subgroup (sort keys %{$info{$group}}) {
                if ($self->report_type eq "full") {
                    for my $params (keys %{$info{$group}->{$subgroup}}) {
                        $output_fh->print($group . "\t");
                        $output_fh->print($subgroup . "\t");
                        $output_fh->print($params . "\t");
                        $output_fh->print($info{$group}->{$subgroup}->{$params} . "\n");
                    }
                }

                elsif ($self->report_type eq "abbreviated") {
                    my $sum = 0;
                    for my $params (keys %{$info{$group}->{$subgroup}}) {
                        $sum += $info{$group}->{$subgroup}->{$params};
                    }
                    $output_fh->print($group . "\t");
                    $output_fh->print($subgroup . "\t");
                    $output_fh->print($sum . "\n");
                }
            }
        }
    }
    close $output_fh unless $output_fh eq 'STDOUT';
}
1;

