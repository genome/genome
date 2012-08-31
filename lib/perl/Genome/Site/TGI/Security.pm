
package Genome::Site::TGI::Security;

use File::Basename;
use DateTime;
use IO::File;
use Sys::Hostname;

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
    
    my $command = basename($0) . " ";
    while (@argv) {
        last unless defined $argv[0] and $argv[0] !~ /^-/;
        $command .= (shift @argv) . " ";
    }

    my $params = join(" ", @argv);
    my $dt = DateTime->now;
    $dt->set_time_zone('America/Chicago');
    my $date = $dt->ymd;
    my $time = $dt->hms;
    my $host = hostname;

    my $base_log_dir = $ENV{GENOME_LOG_DIR};
    my $log_dir = join('/', $base_log_dir, 'command_line_usage', $dt->year);
    my $log_file = $dt->month . "-" . $dt->day . ".log";
    my $log_msg = join("\t", $date, $time, Genome::Sys->username, $command, $params);

    unless (-e $log_dir and -d $log_dir) {
        mkdir $log_dir;
        chmod(0777, $log_dir);
    }

    my $log_path = join('/', $log_dir, $log_file);
    my $log_fh = IO::File->new($log_path, 'a');
    unless ($log_fh) {
        print STDERR "Could not get file handle for log file at $log_path, command execution will continue!\n";

        my $email_msg = "User: " . Genome::Sys->username . "\n" .
                        "Date: $date\n" .
                        "Time: $time\n" .
                        "Host: $host\n" .
                        "Command: $command\n" .
                        "Could not write to log file at $log_path : $!\n";

        App::Mail->mail(
            To      => 'jeldred@genome.wustl.edu',
            From    => Genome::Config->user_email,
            Subject => "Error writing to log file",
            Message => $email_msg
        );
        exit 0;
    }

    flock($log_fh, 2);
    chmod(0666, $log_path) unless -s $log_path;
    print $log_fh "$log_msg\n";
    close $log_fh;
    exit 0;
}

1;

