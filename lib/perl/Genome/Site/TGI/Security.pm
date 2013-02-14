use strict;
use warnings;

package Genome::Site::TGI::Security;

use File::Spec qw();
use File::Basename;
use DateTime;
use IO::File;
use Sys::Hostname;

sub command_class {
    my @argv = @ARGV;

    my $command = basename($0);
    my %command_class = (
        genome => 'Genome::Command',
        gmt    => 'Genome::Model::Tools',
    );
    my $command_class = $command_class{$command};

    while ($command_class && @argv) {
        last unless ($command_class->is_sub_command_delegator);
        last if ($argv[0] && index($argv[0], '-') >= 0);
        $command_class = $command_class->class_for_sub_command(shift @argv) || $command_class;
    }

    return $command_class;
}

sub base_log_dir {
    return File::Spec->join($ENV{GENOME_LOG_DIR}, 'command_line_usage');
}

sub log_columns {
    return qw(date time user command_class command command_with_args);
}

# Get info about command and write to log file.
# Send error information if log file isn't found or can't be accessed, which is usually
# a good indication that there are disk issues afoot.
sub log_command {
    my @argv = @ARGV;

    # Fork twice... grandchild process will be an orphan, so parent won't wait for it to complete
    my $pid = (defined(&DB::DB)? -1 : UR::Context::Process->fork()); #don't fork in the debugger
    if ($pid) {
        return 1;
    }

    my $command_class = command_class() || '-';
    my $command = $0;
    my $command_with_args = join(' ', basename($command), @argv);

    my $params = join(" ", @argv);
    my $dt = DateTime->now;
    $dt->set_time_zone('America/Chicago');
    my $date = $dt->ymd;
    my $time = $dt->hms;
    my $host = hostname;

    my $log_dir = join('/', base_log_dir(), $dt->year);
    my $log_file = $dt->month . "-" . $dt->day . ".log";
    my %log_data = (
        date => $date,
        time => $time,
        user => Genome::Sys->username,
        command_class => $command_class,
        command => $command,
        command_with_args => $command_with_args,
    );
    my $log_msg;
    for my $column (log_columns()) {
        unless ($log_data{$column}) {
            warn "missing $column!";
        }
        $log_msg = (defined($log_msg) ? join("\t", $log_msg, $log_data{$column}) : $log_data{$column});
    }

    unless (-e $log_dir and -d $log_dir) {
        mkdir $log_dir;
        chmod(0777, $log_dir);
    }

    my $log_path = join('/', $log_dir, $log_file);
    my $log_fh = IO::File->new($log_path, 'a');
    unless ($log_fh) {
        print STDERR "Could not get file handle for log file at $log_path, command execution will continue!\n";

        exit 0;
    }

    flock($log_fh, 2);
    chmod(0666, $log_path) unless -s $log_path;
    print $log_fh "$log_msg\n";
    close $log_fh;
    exit 0;
}

1;

