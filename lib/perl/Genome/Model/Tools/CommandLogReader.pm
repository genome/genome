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
            valid_values => ["user","command"],
            default => "user",
            doc => "Groups output by user name or command name, defaults to user.",
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

sub log_attributes {
    my $self = shift;
    return (qw/ date time user command params/);
}

sub log_directory_base {
    my $self = shift;
    return "/gsc/var/log/genome/command_line_usage";
}

sub output_attributes {
    my $self = shift;
    if ($self->group_by eq "user") {
        if ($self->report_type eq "full") {
            return (qw/ user command params times_called/);
        }
        elsif ($self->report_type eq "abbreviated") {
            return (qw/ user command times_called/);
        }
        elsif ($self->report_type eq "summary") {
            return (qw/ user num_commands/);
        }
    }
    elsif ($self->group_by eq "command") {
        if ($self->report_type eq "full") {
            return (qw/ command user params times_called/);
        }
        elsif ($self->report_type eq "abbreviated") {
            return (qw/ command user times_called/);
        }
        elsif ($self->report_type eq "summary") {
            return (qw/ command times_called/);
        }
    }
    else {
        $self->error_message("Unknown group by attribute or report type, unable to provide column headers");
        return;
    }
}

sub _sort_params {
    my $self = shift;
    my $params = shift;

    my @params = sort(split('--', $params));
    @params = sort(split('/', $params)) unless @params;
    return unless @params;

    my $str;
    for my $p (@params) {
        next if length $p == 0;
        $p =~ s/^\s+|\s+$//g; # Remove leading and trailing spaces
        unless ($p =~ /=/) {  # Replace first space with = unless there's already an =
            $p =~ s/\s+/=/;
        }
        $str .= $p . " ";
    }
    chop $str;
    return $str;
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

    unless (-e $self->log_directory_base and -d $self->log_directory_base) {
        $self->error_message("Log directory $self->log_directory_base does not exist, exiting");
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

    my $end_date = DateTime->now;
    if (defined $self->end_date) {
        my ($end_year, $end_month, $end_day) = split(/-/, $self->end_date);
        $end_date = DateTime->new(year => $end_year, month => $end_month, day => $end_day);
    }

    my %info;
    my @log_columns = $self->log_attributes();
    until (DateTime->compare($start_date, $end_date) == 1) {
        my $log_dir = $self->log_directory_base . "/" . $start_date->year . "/";
        my $log_file = $log_dir . $start_date->month . "-" . $start_date->day . ".log";

        $start_date->add(days => 1);

        unless (-e $log_file) {
            $self->warning_message("Could not find $log_file, continuing.");
            next;
        }

        my $log_svr = Genome::Utility::IO::SeparatedValueReader->create(
            input => $log_file,
            headers => \@log_columns,
            separator => "\t",
            is_regex => 1,
            ignore_extra_columns => 1,
        );

        while (my $line = $log_svr->next) {
            my $user = $line->{user};
            my $command = $line->{command};
            my $params = $self->_sort_params($line->{params});
            if ($params) {
                if ($self->group_by eq "command"){
                    $info{$command}->{$user}->{$params}++;
                }
                elsif ($self->group_by eq "user") {
                    $info{$user}->{$command}->{$params}++;
                }
            }
            else {
                if ($self->group_by eq "command") {
                    $info{$command}->{$user}->{'none'}++;
                }
                elsif ($self->group_by eq "user") {
                    $info{$user}->{$command}->{'none'}++;
                }
            }
        }
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

