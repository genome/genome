package Genome::Sys::ShellPipeline;

use strict;
use warnings;

use Genome;
use Try::Tiny qw(try catch finally);
use Carp qw(croak);
use File::Spec qw();
use List::AllUtils qw(any);
use File::Slurp qw(read_file write_file);

class Genome::Sys::ShellPipeline {
    is => "UR::Object",
    attributes_have => [qw(is_input is_output is_optional)],
    has_input => {
        pre_commands => {
            is => "ARRAY",
            doc => "Arrayref of commands to execute before the main pipe commands",
            default_value => [],
        },

        post_commands => {
            is => "ARRAY",
            doc => "Arrayref of commands to execute after the main pipe commands",
            default_value => [],
        },

        pipe_commands => {
            is => "ARRAY",
            doc => "Arrayref of commands to be piped together",
        },

        redirects => {
            is => "Text",
            doc => "If specified, this text is appended to the concatenated pipe " .
                "commands without a pipe inbetween (e.g., > /tmp/stdout.txt 2>> /tmp/log)",
            is_optional => 1,
        },

        interpreter => {
            is => "Text",
            doc => "The interpreter to use to run the pipe script",
            default_value => "/bin/bash",
        },

        script_path => {
            is => "Text",
            doc => "If specified, write the pipe script to this path. " .
                "By default, it is written to a temporary location.",
            is_optional => 1,
        },
    },

    has_output => [
        return_codes => {
            is => "ARRAY",
            doc => "The return code of each command in pipe_commands",
        },
    ]
};

sub _validate_params {
    my $self = shift;

    my $commands = $self->pipe_commands;
    croak $self->error_message("pipe_commands argument must be a nonempty arrayref")
        unless (ref $commands eq "ARRAY" && scalar @$commands > 0);

    for my $x (qw(pre post)) {
        my $var = sprintf "%s_commands", $x;
        my $value = $self->$var;
        croak $self->error_message("$var argument must be undef or an arrayref")
            unless !defined $value || ref $value eq "ARRAY";
    }

}

sub execute {
    my $self = shift;

    $self->_validate_params;

    my $script_path = $self->script_path;
    my $tmpdir = Genome::Sys->create_temp_directory;
    if (!defined $script_path) {
        $script_path = File::Spec->catfile($tmpdir, "pipe.sh");
    }

    my $pipestatus_path = File::Spec->catfile($tmpdir, "pipe.status");

    write_file($script_path, $self->script_text($pipestatus_path));

    $self->debug_message(
        sprintf
            "Executing script $script_path:\n\n" .
            "-------- BEGIN SCRIPT --------\n%s" .
            "--------- END SCRIPT ---------\n",
            join("", read_file($script_path))
        );


    chmod oct(755), $script_path;

    my @errors;

    try {
        Genome::Sys->shellcmd(cmd => $script_path);
    }
    catch {
        push @errors, $_;
    };

    push @errors, $self->_check_pipestatus($pipestatus_path);
    $self->error_message($_) for @errors;
    croak $self->error_message() if @errors;

    return 1;
}

sub _make_pipe {
    my @cmds = @_;
    return join " \\\n    | ", @cmds;
}

sub script_text {
    my ($self, $pipestatus_path) = @_;

    my $commands = $self->pipe_commands;
    my $pre = join "\n", @{$self->pre_commands};
    my $post = join "\n", @{$self->post_commands};

    my $cmd = _make_pipe @$commands;
    $cmd .= " " . $self->redirects if $self->redirects;

    my $interp = $self->interpreter;
    croak $self->error_message("Script interpreter $interp is not executable!") unless -x $interp;
    return <<EOS;
#!$interp

$pre

set -o pipefail
$cmd
PSTAT=\${PIPESTATUS[@]}

echo \$PSTAT > $pipestatus_path

for x in \$PSTAT
do
    if [ "\$x" != "0" ]
    then
        exit 1
    fi
done

$post
EOS
}

# This sub tries to parse a file containing a single line containing the
# result of ${PIPESTATUS[@]} (a space delimited list of exit codes from
# each command in the pipe expression.
sub _check_pipestatus {
    my ($self, $path) = @_;

    return "PIPESTATUS error file does not exist" unless -s $path;

    my @lines = read_file $path;
    chomp @lines;

    # there should only be one line in the file
    return "Too many lines in PIPESTATUS error file" unless scalar @lines == 1;

    my @status = split(/\s/, $lines[0]);
    $self->return_codes(\@status);

    my @commands = @{$self->pipe_commands};
    # the number of commands we intended to execute should be equal
    # to the number of exit codes we received.
    return "Wrong number of status codes in PIPESTATUS error file"
        unless scalar @status == scalar @commands;

    return unless any { $_ != 0 } @status;

    # which command failed?
    my @failures = grep {$status[$_] != 0} 0..$#status;
    return sprintf "The following commands crashed:\n  %s",
        join("\n  ",
            map {sprintf "%s: %s", $_, $commands[$_]} @failures
            );
}
