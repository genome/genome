package Genome::DataSource::CommonRDBMS;

use Exporter qw(import);
our @EXPORT_OK = qw(log_error log_commit_time);

use Genome;
use UR::DataSource::Pg;


class Genome::DataSource::CommonRDBMS {
    doc => 'Mixin class to implement pausing access to the database',
};

my $query_pause = _make_db_pause_function('query_pause_sentry_file_path');
foreach my $signal ( qw( query precreate_handle sequence_nextval ) ){
    __PACKAGE__->create_subscription(
        method => $signal,
        callback => $query_pause,
    );
}

my $commit_pause = _make_db_pause_function('commit_pause_sentry_file_path');
UR::Context->current->create_subscription(
    method => 'precommit',
    callback => sub {
        __PACKAGE__->$commit_pause;
    }
);

sub query_pause_sentry_file_path {
    $ENV{GENOME_DB_QUERY_PAUSE};
}

sub commit_pause_sentry_file_path {
    $ENV{GENOME_DB_PAUSE};
}

sub pause_sleep_length_seconds { 30 }

sub _make_db_pause_function {
    my $pause_check_method = shift;

    return sub {
        my $self = shift;

        my $sentry = $self->$pause_check_method;
        return 1 unless $sentry and -e $sentry;

        print STDERR "Database querying has been paused; disconnecting db handles.  "
                    . "Please wait until the query pause is released...\n";

        my @data_sources = UR::Context->all_objects_loaded(__PACKAGE__);
        for my $ds (@data_sources) {
            $ds->disconnect_default_handle if $ds->has_default_handle;
        }

        while (1) {
            sleep pause_sleep_length_seconds;
            last unless $sentry and -e $sentry;
        }

        print STDERR "Database updating has been resumed, continuing query!\n";
        return 1;
    };
}


sub _check_pg_version {
    my $class = shift;

    my $required_pg_version = '2.19.3';
    require DBD::Pg;
    my $pg_version = $DBD::Pg::VERSION;
    if (($pg_version ne $required_pg_version) && !defined $ENV{'LIMS_PERL'}) {
        my $error_message = "**** INCORRECT POSTGRES DRIVER VERSION ****\n" .
            "You are using a Perl version that includes an incorrect DBD::Pg driver.\n" .
            "You are running $pg_version and need to be running $required_pg_version.\n" .
            "Your sync has been aborted to protect data integrity in the Postgres database.\n" .
            "Please be sure you are using 'genome-perl' and not /gsc/bin/perl.\n\n\n" .
            "This event has been logged with apipe; if you are unsure of why you received this message\n" .
            "open an apipe-support ticket with the date/time of occurrence and we will assist you.\n";
        UR::Object->error_message($error_message);
        log_error($error_message);
        die $error_message;
    }
}


sub log_error {
    my $error = shift;
    my $log_string = create_log_message($error);
    my $log_fh = open_error_log();
    # If we can't get a file handle to the log file, no worries. Just continue without making a peep.
    if ($log_fh) {
        my $lock_status = File::lockf::lock($log_fh);
        # this returns 0 on success
        unless ($lock_status != 0) {
            $log_fh->print("$log_string\n");
            File::lockf::ulock($log_fh);
        }
        $log_fh->close;
    }
}

sub log_commit_time {
    my($db_name, $time) = @_;

    my $original_cmdline = get_command_line();
    my $execution_id = $ENV{'GENOME_EXECUTION_ID'} || '';

    my $path = _determine_base_log_pathname();
    $path .= "-${db_name}.timing";
    my $fh = IO::File->new($path,'>>');
    unless ($fh) {
        print STDERR "Couldn't open $path for appending: $!\n";
        return;
    }
    chmod(0666,$path);
    my $lock_status = File::lockf::lock($fh);
    # lock gives us a zero return value on success
    unless ($lock_status != 0) {
        $fh->print(join("\t",$THIS_COMMIT_ID, $execution_id, $time, $original_cmdline) . "\n");
        File::lockf::ulock($fh);
    }
    $fh->close();
}

sub get_command_line {
    my @commands;
    if ($INC{'Command/V2.pm'}) {
        push @commands, Command::V2->get('original_command_line true' => 1);
    }
    if ($INC{'Command/V1.pm'}) {
        push @commands, Command::V1->get('original_command_line true' => 1);
    }
    @commands = sort { $a->id cmp $b->id } @commands;
    my $original_cmdline = $commands[0] ? $commands[0]->original_command_line : $0;

    return $original_cmdline;
}


sub create_log_message {
    my $error = shift;

    require DateTime;
    my $dt = DateTime->now(time_zone => 'America/Chicago');
    my $date = $dt->ymd;
    my $time = $dt->hms;
    my $user = Genome::Sys->username;
    my $acting_user = getlogin();
    my $perl_path = $^X;

    my $original_cmdline = get_command_line();

    my $path = _determine_base_log_pathname();

    require Sys::Hostname;
    my $host = Sys::Hostname::hostname();

    my $string = join("\n", join(',', $date, $time, $host, $user, $acting_user, $perl_path, $original_cmdline, $THIS_COMMIT_ID), $error);
    return $string;
}

sub _determine_base_log_pathname {
    require DateTime;
    my $dt = DateTime->now(time_zone => 'America/Chicago');
    my $date = $dt->ymd;

    my $base_dir = '/gsc/var/log/genome/postgres/';
    my $dir = join('/', $base_dir, $dt->year);
    unless (-d $dir) {
        mkdir $dir;
        chmod(0777, $dir);
    }

    my $file = join('-', $dt->month, $dt->day);
    return join('/',$dir, $file);
}

sub open_error_log {
    my $path = _determine_base_log_pathname();
    $path .= '.log';
    my $fh = IO::File->new($path, '>>');
    unless ($fh) {
        print STDERR "Couldn't open $path for appending: $!\n";
        return;
    }
    chmod(0666,$path);
    return $fh;
}


BEGIN {
    __PACKAGE__->_check_pg_version();
}


1;
