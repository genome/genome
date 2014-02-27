package Genome::Model::Tools::Park::Base;

use strict;
use warnings;

use Genome;
use File::Basename qw(dirname);
use Carp qw(confess);
use Params::Validate qw(:types validate_pos);
use IPC::Run qw(start pump finish);

class Genome::Model::Tools::Park::Base {
    is => 'UR::Object',
    is_abstract => 1,
};

sub _template_path {
    die 'abstract',
}

sub _link_process {
    my ($self, $process_uri, $target) = @_;

    my $link_script = $self->_rex_script_path('process link');
    Genome::Sys->shellcmd(cmd => "$link_script --process $process_uri --target $target");
    return;
}

sub _run_rex_process {
    my ($self, $source_path, $inputs_filename) = validate_pos(@_, 1,
        {type => SCALAR}, {type => SCALAR});

    my $script_path = $self->_rex_script_path('process start');
    my $cmd = [$script_path, '--source-path', $source_path,
        '--inputs-file', $inputs_filename];

    return $self->_run_generic_process($cmd);
}

sub _run_generic_process {
    my ($self, $cmd) = validate_pos(@_, 1, {type => ARRAYREF});

    my $in;
    my $out = "x" x 100_000; # preallocate $out to avoid performance penalty
    $out = "";

    my $process_uri;
    eval {
        my $h = start $cmd, \$in, \$out;

        pump $h until $out =~ /Launching Process (\S*) /;
        $process_uri = $1;
        $self->status_message("$out\n");
        $out = "";

        $self->status_message("\nView the progress of this process with:\n");
        $self->status_message("    %s %s\n", $self->_rex_script_path('process view'),
            $process_uri);

        finish $h;
        print "Process completed successfully.\n";
    };
    if ($@) {
        $self->status_message("$out\n");
        die "Process FAILED: $@\n";
    }
    return $process_uri;
}

my %REX_COMMAND_SCRIPTS = (
    'process start' => 'rex_process_start.sh',
    'process link' => 'rex_process_link.sh',
    'process view' => 'rex_process_view.sh',
);

sub _rex_script_path {
    my $self = shift;
    my $command = shift;

    my $script_filename = $REX_COMMAND_SCRIPTS{$command};
    unless (defined($script_filename)) {
        confess "Couldn't find script for command ($command)";
    }

    my $this_dir = dirname(__FILE__);
    my $full_path = File::Spec->join($this_dir, $script_filename);
    unless (-x $full_path) {
        confess "Couldn't execute script ($full_path) for command ($command)";
    }
    return $full_path;
}


sub _generate_inputs_file {
    my $self = shift;

    my $template_fh = Genome::Sys->open_file_for_reading($self->_template_path);
    my ($inputs_fh, $inputs_filename) = Genome::Sys->create_temp_file();

    my $line_number = 0;
    while (my $line = $template_fh->getline) {
        $line_number++;
        my ($workflow_input_name, $accessor) = _parse_line($line);

        unless ($self->can($accessor)) {
            confess(sprintf("Class (%s) has no accessor (%s) that the template (%s) expected on line (%d)",
                    ref $self, $accessor, $self->_template_path, $line_number));
        }
        if (!defined($self->$accessor)) {
            confess(sprintf("Attribute (%s) is undefined, which is disallowed in an inputs file.", $accessor));
        }
        my $value = [$self->$accessor] if ref($self->$accessor) ne 'ARRAY';

        my $string = join("\t", @$value);
        printf $inputs_fh "%s\t%s\n", $workflow_input_name, $string;
    }
    return $inputs_filename;
}

sub _parse_line {
    my $line = shift;

    chomp($line);
    my @parts = split(/\t/, $line);
    my $workflow_input_name = shift(@parts);
    my $accessor = shift(@parts);

    unless (defined($accessor)) {
        $accessor = (split(/\./, $workflow_input_name))[-1];
    }
    return ($workflow_input_name, $accessor);
}

1;
