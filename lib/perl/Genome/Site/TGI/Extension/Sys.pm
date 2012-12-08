package Genome::Sys;

# Short Term: shellcmd() should probably be rewritten, it does not correctly use $! after the system call.  Would also be nice if
# it could support IO wrapping of command being executed.  shellcmd() might be better off in its own module, since its not strictly
# a filesystem function.

use strict;
use warnings;

use Genome;
use Genome::Sys;  # ensure our overrides take precedence

use Time::HiRes;
use Genome::Utility::Instrumentation;

use Data::Dumper;
require Carp;
require IO::Dir;
require IO::File;
require File::Basename;
require File::Path;
require File::Copy;
require Genome::Utility::Text;
use Sys::Hostname;
use File::Find;
#use Archive::Extract;

require MIME::Lite;

#####
# Methods useful for bsubbing jobs and checking their status
# FIXME There's a lower level API to LSF that doesn't rely on a the command line 
# interface. That should be used here, but this works well enough for now.
#####

sub bsub_and_wait {
    my $class = shift;
    my %params = @_;
    my $job_id = Genome::Sys->bsub(%params);
    my $status = Genome::Sys->wait_for_lsf_job($job_id);
    return ($job_id, $status);
}

# FIXME This is incomplete and should be expanded to accept more bsub parameters
sub bsub {
    my $class = shift;
    my %params = @_;
    my $queue = $params{queue} || 'long';
    my $job_group = $params{job_group};
    my $cmd = $params{cmd};
    my $log_file = $params{log_file};
    unless ($cmd) {
        die "Must be given a command to bsub!";
    }

    my $bsub_cmd = "bsub";
    $bsub_cmd .= " -q $queue";
    $bsub_cmd .= " -g $job_group" if $job_group;
    $bsub_cmd .= " -o $log_file" if $log_file;
    $bsub_cmd .= " $cmd";

    my $bsub_output = `$bsub_cmd`;
    my ($job_id) = $bsub_output =~ /Job <(\d+)> is submitted to/;
    unless ($job_id) {
        die "Could not get job id from bsub output!";
    }
    return $job_id;
}

sub get_lsf_job_status {
    my ($class, $job_id) = @_;
    unless ($job_id) {
        die "Must be given a job id!";
    }
    my $cmd = "bjobs -aw $job_id 2>&1";

    my $status;
    for (1..3) {
        my $output = `$cmd`;
        my ($headers, $line) = split("\n", $output);
        if ($headers and $headers =~ /^JOBID/ and $line) {
            (undef, undef, $status) = split(/\s+/, $line); 
            last if $status;
        }
        sleep 10;
    }
    unless ($status) {
        die "Could not get status for job $job_id";
    }
    return $status;
}

sub wait_for_lsf_job {
    my ($class, $job_id) = @_;
    unless ($job_id) {
        die "Must be given job id!";
    }

    my $status;
    while (1) {
        $status = Genome::Sys->get_lsf_job_status($job_id);
        last if $status eq 'DONE' or $status eq 'EXIT';
        sleep 10;
    }
    return $status;
}

sub wait_for_lsf_jobs {
    my ($class, @job_ids) = @_;
    unless (@job_ids) {
        die "Must be given job ids!";
    }

    my %job_statuses;
    my $all_jobs_complete = 1;
    $DB::single = 1;
    do {
        $all_jobs_complete = 1;
        JOB_ID: for my $job_id (@job_ids) {
            if (exists $job_statuses{$job_id} and ($job_statuses{$job_id} eq 'DONE'
                    or $job_statuses{$job_id} eq 'EXIT')) {
                next JOB_ID;
            }
            my $status = $class->get_lsf_job_status($job_id);
            $job_statuses{$job_id} = $status;
            $all_jobs_complete = 0 unless $status eq 'DONE' or $status eq 'EXIT';
        }
    } while (!$all_jobs_complete);

    return %job_statuses;
}

sub kill_lsf_job {
    my ($class, $job_id) = @_;
    unless ($job_id) {
        die "Must be given job id!";
    }

    my $cmd = "bkill $job_id";
    my $rv = system($cmd);
    $rv = $rv >> 8;
    unless ($rv == 0) {
        die "Could not kill LSF job $job_id, received return status $rv";
    }
    return 1 if $rv == 0;
}


# used by Genome::Sys->snapshot_revision when present
sub _simplify_inc {
    my $self = shift;
    my @inc = @_;
    # if the only path is like /gsc/scripts/opt/genome/snapshots/genome-1213/lib/perl then just call it genome-1213
    # /gsc/scripts/opt/genome/snapshots/genome-1213/lib/perl -> genome-1213
    # /gsc/scripts/opt/genome/snapshots/custom/genome-foo/lib/perl -> custom/genome-foo
    if (@inc == 1 and $inc[0] =~ /^\/gsc\/scripts\/opt\/genome\/snapshots\//) {
        $inc[0] =~ s/^\/gsc\/scripts\/opt\/genome\/snapshots\///;
        $inc[0] =~ s/\/lib\/perl$//;
    }
    return @inc;
}

# this helps us clean-up locks

my %SYMLINKS_TO_REMOVE;

# disk usage management

sub disk_usage_for_path {
    my $self = shift;
    my $path = shift;

    unless (-d $path) {
        $self->error_message("Path $path does not exist!");
        return;
    }

    return unless -d $path;
    my $cmd = "du -sk $path 2>&1";
    my $du_output = qx{$cmd};
    my $kb_used = ( split( ' ', $du_output, 2 ) )[0];
    unless (Scalar::Util::looks_like_number($kb_used)) {
        $self->error_message("du output is not a number: $kb_used");
        return;
    }

    return $kb_used;
}


#< Files >#

sub diff_text_vs_text {
    my ($self,$t1,$t2) = @_;
    my $p1 = $self->create_temp_file_path();
    $self->write_file($p1, $t1);
    my $p2 = $self->create_temp_file_path();
    $self->write_file($p2, $t2);

    return $self->diff_file_vs_file($p1, $p2);
}

sub diff_file_vs_text {
    my ($self,$f1,$t2) = @_;
    my $p2 = $self->create_temp_file_path();
    $self->write_file($p2, $t2);

    return $self->diff_file_vs_file($f1, $p2);
}

sub diff_file_vs_file {
    my ($self,$f1,$f2) = @_;

    my $diff_fh = IO::File->new("sdiff -s $f1 $f2 2>&1 |");
    unless ($diff_fh) {
        Carp::croak("Can't run 'sdiff -s $f1 $f2' for diff_file_vs_file(): $!");
    }
    my $diff_output = do { local( $/ ) ; <$diff_fh> };
    return $diff_output;
}

sub open_file_for_appending {
    my ($self, $file) = @_;

    unless ( defined $file ) {
        Carp::croak("No append file given");
    }

    if ( -d $file ) {
        Carp::croak("Append file ($file) is a directory and cannot be opend as a file");
    }

    my ($name, $dir) = File::Basename::fileparse($file);
    unless ( $dir ) {
        Carp::croak("Can't determine directory from append file ($file)");
    }

    unless ( -w $dir ) {
        Carp::croak("Do not have WRITE access to directory ($dir) for append file ($name)");
    }

    return $self->_open_file($file, 'a');
}

sub validate_file_for_execution {
    my ($self, $file) = @_;

    my $rv = eval { $self->validate_file_for_reading($file) };
    if ($@ or not $rv) {
        my $msg = "Cannot read file $file and therefore cannot execute it";
        $msg .= ", reason: $@" if $@;
        Carp::croak $msg;
    }

    unless (-x $file) {
        Carp::croak "File $file cannot be executed!";
    }

    return 1;
}


sub bzip {
    my $self = shift;
    my $file = shift;

    # validate_file_for_reading throws its own exceptions when there are problems
    $self->validate_file_for_reading($file)
        or return;

    my $bzip_cmd = "bzip2 -z $file";
    my $result_file = $file.".bz2";
    # shellcmd throws its own exceptions when there are problems, including checking existence of output file
    $self->shellcmd(cmd=>$bzip_cmd,
                    output_files=>[$result_file]
                    );

    return $result_file;

}

sub bunzip {
    my $self = shift;
    my $file = shift;

    $self->validate_file_for_reading($file)
        or return;

    if ($file=~m/\.bz2$/) {

        #the -k option will keep the bzip file around
        my $bzip_cmd = "bzip2 -dk $file";

        #get the unzipped file name by removing the .bz2 extension.
        $file=~m/(\S+).bz2/;
        my $result_file = $1;

        $self->shellcmd(cmd=>$bzip_cmd,
                        output_files=>[$result_file],
                        );

        return $result_file;

    } else {
        Carp::croak("Input file ($file) does not have .bz2 extension.  Not unzipping.");
    }

}

sub extract_archive {
    my $self = shift;
    my %params = @_;
    my $from = delete $params{from};
    my $to = delete $params{to};

    if(%params) {
        my @crap = %params;
        Carp::confess("Unknown params passed to extract_archive: @crap");
    }

    unless($from) {
        Carp::croak("No 'from' passed to extract_archive");
    }

    unless($to) {
        Carp::croak("No 'to' passed to extract_archive");
    }

    $self->validate_file_for_reading($from);

    #use binaries rather than perl modules to extract archives.
    #the perl modules will often load the entire archive into RAM
    #we operate on files way too big to be doing that
    $Archive::Extract::PREFER_BIN = 1;

    require Archive::Extract;

    my $archive = Archive::Extract->new(archive => $from);
    unless($archive) {
        Carp::croak("Can't create extract object for archive $from");
    }

    my $rv = $archive->extract( to => $to);
    unless($rv) {
        Carp::croak("Unable to extract $from into $to");
    }

    return $archive;
}



sub open_gzip_file_for_writing {
    my ($self, $file) = @_;

    $self->validate_file_for_writing($file)
        or return;

    if (-e $file) {
        unless (unlink $file) {
            Carp::croak("Can't unlink $file: $!");
        }
    }
    my $pipe = "| bgzip -c > $file";

    return $self->_open_file($pipe);
}

sub open_file_for_overwriting {
    my ($self, $file) = @_;

    if ( not defined $file ) {
        Carp::croak('Cannot open file for over writing. No file given.');
    }

    if ($file eq '-') {
        Carp::croak("Cannot open STDOUT (-) for over writing.");
    }

    if ( -d $file ) {
        Carp::croak("Cannot open file ($file) for over writing. It is a directory.");
    }

    if ( -e $file ) {
        unlink $file;
    }

    my ($name, $dir) = File::Basename::fileparse($file);
    unless ( $dir ) {
        Carp::croak("Cannot open file ($file) for over writing. Failed to get directory from file ($file).");
    }

    unless ( -w $dir ) {
        Carp::croak("Cannot open file ($file) for over writing. Do not have write access to directory ($dir).");
    }

    my $fh = IO::File->new($file, 'w');
    return $fh if $fh;

    Carp::croak("Failed to open file ($file) for over write: $!");
}

sub copy_directory {
    my ($self, $source, $dest) = @_;

    $self->status_message("copying directory: $source to $dest...\n");

    $self->shellcmd(
        cmd => "cp -r '$source' '$dest'",
        input_directories => [$source],
        output_directories => [$dest],
    );

    return 1;
}


#< Dirs >#
sub validate_existing_directory {
    my ($self, $directory) = @_;

    unless ( defined $directory ) {
        Carp::croak("Can't validate_existing_directory: No directory given");
    }

    unless ( -e $directory ) {
        Carp::croak("Can't validate_existing_director: $directory: Path does not exist");
    }


    unless ( -d $directory ) {
        Carp::croak("Can't validate_existing_director: $directory: path exists but is not a directory");
    }

    return 1;
}

sub validate_directory_for_read_access {
    my ($self, $directory) = @_;

    # Both underlying methods throw their own exceptions
    $self->validate_existing_directory($directory)
        or return;

    return $self->_can_read_from_directory($directory);
}

sub validate_directory_for_write_access {
    my ($self, $directory) = @_;

    # Both underlying methods throw their own exceptions
    $self->validate_existing_directory($directory)
        or return;

    return $self->_can_write_to_directory($directory);
}

sub validate_directory_for_read_write_access {
    my ($self, $directory) = @_;

    # All three underlying methods throw their own exceptions
    $self->validate_existing_directory($directory)
        or return;

    $self->_can_read_from_directory($directory)
        or return;

    return $self->_can_write_to_directory($directory);
}

sub recursively_validate_directory_for_read_write_access {
    my ($self, $directory) = @_;

    my $wanted = sub {
        my $full_path = $File::Find::name;
        if (-f $full_path) {
            eval {
                Genome::Sys->validate_file_for_reading($full_path);
                Genome::Sys->validate_file_for_writing_overwrite($full_path);
            };

            if ($@) {
                Carp::croak "Cannot read or write $full_path in $directory!";
            }
        }
    };

    find($wanted, $directory);
    return 1;
}

sub _can_read_from_directory {
    my ($self, $directory) = @_;

    unless ( -r $directory ) {
        Carp::croak("Directory ($directory) is not readable");
    }

    return 1;
}

sub _can_write_to_directory {
    my ($self, $directory) = @_;

    unless ( -w $directory ) {
        Carp::croak("Directory ($directory) is not writable");
    }

    return 1;
}

sub open_directory {
    my ($self, $directory) = @_;

    my $dh = IO::Dir->new($directory);

    unless ($dh) {
        $directory ||= '';
        Carp::croak("Can't open_directory $directory: $!");
    }
    return $dh;
}


# FIXME there are several places where it excludes named pipes explicitly...
# These may not always be appropriate in the general sense, but may be
# for things under Genome::*

sub cat {
    my ($self,%params) = @_;
    my $input_files = delete $params{input_files};
    my $output_file = delete $params{output_file};
    my $mode = ($params{append_mode} ? ">>" : ">");

    my $skip_if_output_is_present = ($params{append_mode}? 0 : 1);

    my @input_files = @$input_files;
    while(my @next = splice(@input_files,0,10)) {
        my $rv = $self->shellcmd(
            cmd => "cat @next $mode $output_file",
            input_files => \@next,
            output_files => [$output_file],
            skip_if_output_is_present => $skip_if_output_is_present,
        );

        unless($rv) {
            die('Failed to cat ' . @next);
        }

        $mode = '>>';
        $skip_if_output_is_present = 0;
    }

    return 1;
}

sub lock_resource {
    my ($self,%args) = @_;

    my $total_lock_start_time = Time::HiRes::time();

    my $resource_lock = delete $args{resource_lock};
    my ($lock_directory,$resource_id,$parent_dir);
    if ($resource_lock) {
        $parent_dir = File::Basename::dirname($resource_lock);
        $self->create_directory($parent_dir);
        unless (-d $parent_dir) {
            Carp::croak("failed to make parent directory $parent_dir for lock $resource_lock!: $!");
        }
    }
    else {
        $lock_directory =  delete $args{lock_directory} || Carp::croak('Must supply lock_directory to lock resource');
        $self->create_directory($lock_directory);
        $resource_id = $args{'resource_id'} || Carp::croak('Must supply resource_id to lock resource');
        $resource_lock = $lock_directory . '/' . $resource_id . ".lock";
        $parent_dir = $lock_directory
    }
    my $wait_on_self = delete $args{wait_on_self} || 0;
    my $basename = File::Basename::basename($resource_lock);

    my $block_sleep = delete $args{block_sleep};
    $block_sleep = 60 unless defined $block_sleep;
    my $max_try = delete $args{max_try};
    $max_try = 7200 unless defined $max_try;
    my $wait_announce_interval = delete $args{wait_announce_interval};
    $wait_announce_interval = 0 unless defined $wait_announce_interval;

    my ($my_host, $my_pid, $my_lsf_id, $my_user) = (hostname, $$, ($ENV{'LSB_JOBID'} || 'NONE'), Genome::Sys->username);
    my $job_id = (defined $ENV{'LSB_JOBID'} ? $ENV{'LSB_JOBID'} : "NONE");
    my $lock_dir_template = sprintf("lock-%s--%s_%s_%s_%s_XXXX",$basename,$my_host,$ENV{'USER'},$$,$job_id);
    my $tempdir =  File::Temp::tempdir($lock_dir_template, DIR=>$parent_dir, CLEANUP=>1);

    unless (-d $tempdir) {
        Carp::croak("Failed to create temp lock directory ($tempdir)");
    }

    # make this readable for everyone
    chmod(0777, $tempdir) or Carp::croak("Can't chmod 0777 path ($tempdir): $!");

    # drop an info file into here for compatibility's sake with old stuff.
    # put a "NOKILL" here on LSF_JOB_ID so an old process doesn't try to snap off the job ID and kill me.
    my $lock_info = IO::File->new("$tempdir/info", ">");
    unless ($lock_info) {
        Carp::croak("Can't create info file $tempdir/info: $!");
    }
    $lock_info->printf("HOST %s\nPID $$\nLSF_JOB_ID_NOKILL %s\nUSER %s\n",
                       $my_host,
                       $ENV{'LSB_JOBID'},
                       $ENV{'USER'},
                     );
    $lock_info->close();

    my $initial_time = time;
    my $last_wait_announce_time = $initial_time;

    my $lock_attempts = 1;
    my $ret;
    while(!($ret = symlink($tempdir,$resource_lock))) {
        # TONY: The only allowable failure is EEXIST, right?
        # If any other error comes through, we end up in bigger trouble.
        use Errno qw(EEXIST ENOENT :POSIX);
        if ($! != EEXIST) {
            $self->error_message("Can't create symlink from $tempdir to lock resource $resource_lock because: $!");
            Carp::croak($self->error_message());
        }
        my $symlink_error = $!;
        chomp $symlink_error;
        return undef unless $max_try--;

        my $target = readlink($resource_lock);

        #If symlink no longer exists
        if ($! == ENOENT || !$target) {
            sleep $block_sleep;
            next;
        }

        #symlink could have disappeared between the top of the while loop and now, so we check to see if it's changed
        my $target_exists = (-e $target);
        unless($target eq readlink($resource_lock)){
            next;
        }

        if (!$target_exists) {
            # TONY: This means the lock symlink points to something that's been deleted
            # That's _really_ bad news and should probably get an email like below.
            $self->error_message("Lock ($resource_lock) exists but target ($target) does not exist.");
            Carp::croak($self->error_message);
        }
        my $target_basename = File::Basename::basename($target);


        $target_basename =~ s/lock-.*?--//;;
        my ($host, $user, $pid, $lsf_id) = split /_/, $target_basename;

        if (not $wait_on_self and $host eq $my_host and $user eq $my_user and $pid eq $my_pid and $lsf_id eq $my_lsf_id) {
            $self->warning_message("Looks like I'm waiting on my own lock, forcing unlock...");
            $self->unlock_resource(resource_lock => $resource_lock, force => 1);
            next;
        }

        my $info_content=sprintf("HOST %s\nPID %s\nLSF_JOB_ID %s\nUSER %s",$host,$pid,$lsf_id,$user);

        my $time = time;
        my $elapsed_time = $time - $last_wait_announce_time;
        if ($elapsed_time >= $wait_announce_interval) {
            $last_wait_announce_time = $time;
            my $total_elapsed_time = $time - $initial_time;
            $self->status_message("waiting (total_elapsed_time = $total_elapsed_time seconds) on lock for resource '$resource_lock': $symlink_error. lock_info is:\n$info_content");
        }

        if ($lsf_id ne "NONE") {
            my ($job_info,$events) = Genome::Model::Event->lsf_state($lsf_id);
                 unless ($job_info) {
                     $self->warning_message("Invalid lock for resource $resource_lock\n"
                                            ." lock info was:\n". $info_content ."\n"
                                            ."Removing old resource lock $resource_lock\n");
                     unless ($Genome::Sys::IS_TESTING) {
                        my $message_content = <<END_CONTENT;
Hey Apipe,

This is a lock attempt on %s running as PID $$ LSF job %s and user %s.

I'm about to remove a lock file held by a LSF job that I think is dead.

The resource is:

%s

Here's info about the job that I think is gone.

%s

I'll remove the lock in an hour.  If you want to save the lock, kill me
before I unlock the process!

Your pal,
Genome::Utility::Filesystem

END_CONTENT

                        my $msg = MIME::Lite->new(From    => sprintf('"Genome::Utility::Filesystem" <%s@genome.wustl.edu>', $ENV{'USER'}),
                                              To      => 'apipe-run@genome.wustl.edu',
                                              Subject => 'Attempt to release a lock held by a dead process',
                                              Data    => sprintf($message_content, $my_host, $ENV{'LSB_JOBID'}, $ENV{'USER'}, $resource_lock, $info_content),
                                            );
                        $msg->send();
                        $self->status_message('Sleeping for one hour...');
                        sleep 60 * 60;
                 }
                     $self->unlock_resource(resource_lock => $resource_lock, force => 1);
                     #maybe warn here before stealing the lock...
               }
           }
        sleep $block_sleep;
        $lock_attempts += 1;
       }
    $SYMLINKS_TO_REMOVE{$resource_lock} = 1;

    # do we need to activate a cleanup handler?
    $self->cleanup_handler_check();

    my $total_lock_stop_time = Time::HiRes::time();
    my $lock_time_miliseconds = 1000 * ($total_lock_stop_time - $total_lock_start_time);

    my $caller_name = _resolve_caller_name(caller());

    Genome::Utility::Instrumentation::timing("lock_resource.$caller_name", $lock_time_miliseconds);
    Genome::Utility::Instrumentation::timing("lock_resource_attempts.$caller_name",
        $lock_attempts);

    Genome::Utility::Instrumentation::timing('lock_resource.total', $lock_time_miliseconds);
    Genome::Utility::Instrumentation::timing('lock_resource_attempts.total',
        $lock_attempts);

    return $resource_lock;
}

sub _resolve_caller_name {
    my ($package, $filename, $line) = @_;
    $package =~ s/::/./g;
    return $package;
}

sub unlock_resource {
    my ($self,%args) = @_;
    my $resource_lock = delete $args{resource_lock};
    my $force = delete $args{force};

    my ($lock_directory,$resource_id);
    unless ($resource_lock) {
        $lock_directory =  delete $args{lock_directory} || Carp::croak('Must supply lock_directory to lock resource');
        $resource_id = $args{'resource_id'} || Carp::croak('Must supply resource_id to lock resource');
        $resource_lock = $lock_directory . '/' . $resource_id . ".lock";
    }

    my $target = readlink($resource_lock);
    if (!$target) {
        if ($! == ENOENT) {
            $self->error_message("Tried to unlock something that's not locked -- $resource_lock.");
            Carp::croak($self->error_message);
        } else {
            $self->error_message("Couldn't readlink $resource_lock: $!");
        }
    }
    unless (-d $target) {
        $self->error_message("Lock symlink '$resource_lock' points to something that's not a directory - $target. ");
        Carp::croak($self->error_message);
    }
    my $basename = File::Basename::basename($target);
    $basename =~ s/lock-.*?--//;;
    my ($thost, $tuser, $tpid, $tlsf_id) = split /_/, $basename;
    my $my_host = hostname;
    my $my_job_id = (defined $ENV{'LSB_JOBID'} ? $ENV{'LSB_JOBID'} : "NONE");

    unless ($force) {
        unless ($thost eq $my_host
             && $tuser eq $ENV{'USER'}
             && $tpid eq $$
             && $tlsf_id eq $my_job_id) {

             $self->error_message("This lock does not look like it belongs to me.  $basename does not match $my_host $ENV{'USER'} $$ $my_job_id.");
             delete $SYMLINKS_TO_REMOVE{$resource_lock}; # otherwise the lock would be forcefully cleaned up when process exits
             Carp::croak($self->error_message);
        }
    }

    my $unlink_rv = unlink($resource_lock);
    if (!$unlink_rv) {
        $self->error_message("Failed to remove lock symlink '$resource_lock':  $!");
        Carp::croak($self->error_message);
    }

    my $rmdir_rv = File::Path::rmtree($target);
    if (!$rmdir_rv) {
        $self->error_message("Failed to remove lock symlink target '$target', but we successfully unlocked.");
        Carp::croak($self->error_message);
    }

    delete $SYMLINKS_TO_REMOVE{$resource_lock};
    $self->cleanup_handler_check();
    return 1;
}

# This method does _not_ throw exceptions since it seems like a non-critical method
sub check_for_path_existence {
    my ($self,$path,$attempts) = @_;

    unless (defined $attempts) {
        $attempts = 5;
    }

    while ($attempts-- > 0) {
        return 1 if -e $path;
        sleep(1);
    }
    return;
}


# FIXME - I think this is a private function to Filesystem.pm
sub cleanup_handler_check {
    my $self = shift;

    my $symlink_count = scalar keys %SYMLINKS_TO_REMOVE;

    if ($symlink_count > 0) {
        $SIG{'INT'} = \&INT_cleanup;
        $SIG{'TERM'} = \&INT_cleanup;
        $SIG{'HUP'} = \&INT_cleanup;
        $SIG{'ABRT'} = \&INT_cleanup;
        $SIG{'QUIT'} = \&INT_cleanup;
        $SIG{'SEGV'} = \&INT_cleanup;
    } else {
        delete $SIG{'INT'};
        delete $SIG{'TERM'};
        delete $SIG{'HUP'};
        delete $SIG{'ABRT'};
        delete $SIG{'QUIT'};
        delete $SIG{'SEGV'};
    }
}

END {
    exit_cleanup();
};

sub INT_cleanup {
    exit_cleanup();
    print STDERR "INT/TERM cleanup activated in Genome::Utility::Filesystem\n";
    Carp::confess;
}

sub exit_cleanup {
    for my $sym_to_remove (keys %SYMLINKS_TO_REMOVE) {
        if (-l $sym_to_remove) {
            warn("Removing remaining resource lock: '$sym_to_remove'");
            unlink($sym_to_remove) or warn "Can't unlink $sym_to_remove: $!";
        }
    }
}


sub get_classes_in_subdirectory {
    my ($subdirectory) = @_;

    unless ( $subdirectory ) {
        Carp::croak("No subdirectory given to get_classes_in_subdirectory");
    }

    my $genome_dir = Genome->get_base_directory_name();
    my $inc_directory = substr($genome_dir, 0, -7);
    unless ( $inc_directory ) {
        Carp::croak("Could not get inc directory for Genome.\n");
    }

    my $directory = $inc_directory.'/'.$subdirectory;
    return unless -d $directory;

    my @classes;
    for my $module ( glob("$directory/*pm") ) {
        $module =~ s#$inc_directory/##;
        $module =~ s#\.pm##;
        $module =~ s#/#::#g;
        push @classes, $module;
    }

    return @classes;
}

sub get_classes_in_subdirectory_that_isa {
    my ($subdirectory, $isa) = @_;

    unless ( $isa ) {
        Carp::confess("No isa given to get classes in directory that isa\n");
    }

    my @classes;
    for my $class ( get_classes_in_subdirectory($subdirectory) ) {
        next unless $class->isa($isa);
        push @classes, $class;
    }

    return @classes;
}


sub directory_size_recursive {
    my ($self,$directory) = @_;
    my $size;
    unless (-e $directory) {
        Carp::croak("directory $directory does not exist");
    }
    find(sub { $size += -s if -f $_ }, $directory);
    return $size;
}

sub is_file_ok {
    my ($self, $file) = @_;

    my $ok_file = $file.".ok";

    #if the file exists and is ok, return 1
    my $rv;
    eval{$rv = $self->validate_file_for_reading($file)};
    if ($rv) {
        if (-e $ok_file) {
            return 1;
        } else {
            #if the file exists, but is not ok, erase the file, return
            my $unlink_rv = unlink($file);
            $self->status_message("File $file not ok.  Deleting.");
            if ($unlink_rv ne 1) {
               Carp::croak($self->error_message("Can't unlink $file.  No ok file found."));
            }
            return;
        }
    } else {
        #if the file doesn't exist, but the ok file does, unlink the ok file.
        if (-e $ok_file) {
        	$self->status_message("File $ok_file exists but does not have an original file.  Deleting.");
            my $unlink_rv = unlink($ok_file);
            if ($unlink_rv ne 1) {
               Carp::croak($self->error_message("Can't unlink $ok_file.  No original file found."));
            }
            return;
        }
    }

    return;

}

sub mark_file_ok {
    my ($self, $file) = @_;

    my $ok_file = $file.".ok";

    if (-f $file ) {
        my $touch_rv = $self->shellcmd(cmd=>"touch $ok_file");
        if ($touch_rv ne 1) {
            Carp::croak($self->error_message("Can't touch ok file $ok_file."));
        } else {
            return 1;
        }
    } else {
    	$self->status_message("Not touching.  Cannot validate file for reading: ".$file);
    }
    return;
}

sub mark_files_ok {
	my ($self,%params) = @_;
	my $input_files = delete $params{input_files};
	for my $input_file (@$input_files) {
		$self->status_message("Marking file: ".$input_file);
		$self->mark_file_ok($input_file);
	}
	return 1;
}


sub are_files_ok {

	my ($self,%params) = @_;
	my $input_files = delete $params{input_files};
	my $all_ok = 1;
	for my $input_file (@$input_files) {
		if (!$self->is_file_ok($input_file) ) {
			$all_ok = 0;
		}
	}

	if ($all_ok != 1) {
    	#delete all the files and start over
    	$self->status_message("Files are NOT OK.  Deleting files: ");
    	$self->status_message(join("\n",@$input_files));
    	for my $file (@$input_files) {
            if (-e $file){
                unlink($file) or Carp::croak("Can't unlink $file: $!");
            }
            if (-e "$file.ok"){
                unlink("$file.ok") or Carp::croak("Can't unlink ${file}.ok: $!");
            }
    	}
    	return;
    } else {
    	#shortcut this step, all the required files exist.
    	$self->status_message("Expected output files already exist.");
   	    return 1;
    }

	return;
}

sub remove_directory_tree {
    my ($self, $directory) = @_;
    unless (-d $directory) {
        $self->warning_message("No directory found at $directory, cannot remove");
        return;
    }

    File::Path::remove_tree($directory, { error => \my $remove_errors });
    if (@$remove_errors) {
        my $error_summary;
        for my $error (@$remove_errors) {
            my ($file, $message) = %$error;
            if ($file eq '') {
                $error_summary .= "General error encountered, message: $message\n";
            }
            else {
                $error_summary .= "File $file, message: $message\n";
            }
        }

        if ($error_summary) {
            $self->error_message("Problems encountered removing $directory\n$error_summary");
            return;
        }
    }
    return 1;
}

sub get_mem_total_from_proc {
    my $mem_total;
    my $meminfo_fh = IO::File->new('/proc/meminfo', 'r');
    if ($meminfo_fh) {
        while (my $meminfo = $meminfo_fh->getline) {
            ($mem_total) = $meminfo =~ /MemTotal:\s+(\d+)/;
            last if ($mem_total);
        }
    }
    return $mem_total;
}

sub get_mem_limit_from_bjobs {
    my $mem_limit;
    my $LSB_JOBID = $ENV{LSB_JOBID};
    my $bjobs_cmd = qx(which bjobs);
    if ($bjobs_cmd && $LSB_JOBID) {
        chomp $bjobs_cmd;
        my $bjobs = qx($bjobs_cmd -l $LSB_JOBID);
        my ($bjobs_mem_limit_kb) = $bjobs =~ /MEMLIMIT\s+(\d+)/;
        $mem_limit = $bjobs_mem_limit_kb if ($bjobs_mem_limit_kb);
    }
    return $mem_limit;
}


# detect maximum available memory
# would probably be better to not do this this way but it's a start
sub mem_limit_kb {
    my $class = shift;

    my $mem_limit_kb;

    # get physical total memory
    if (-e '/proc/meminfo') {
        my $mem_total = $class->get_mem_total_from_proc;
        $mem_limit_kb = $mem_total if $mem_total;
    }

    # get LSF memory limit
    if ($ENV{LSB_JOBID}) {
        my $mem_limit = $class->get_mem_limit_from_bjobs;
        $mem_limit_kb = $mem_limit if $mem_limit;
    }

    return $mem_limit_kb;
}

1;

=pod

=head1 Name

Genome::Sys;

=head1 Synopsis

Houses some generic file and directory methods

=head1 Usage

 require Genome::Sys;

 # Call methods directly:
 Genome::Sys->create_directory($new_directory);

=head1 Methods for Files

=head2 validate_file_for_reading

 Genome::Sys->validate_file_for_reading('/tmp/users.txt')
    or ...;

=over

=item I<Synopsis>   Checks whether the given file is defined, exists, and is readable

=item I<Arguments>  file (string)

=item I<Returns>    true on success, false on failure

=back

=head2 open_file_for_reading

 Genome::Sys->open_file_for_reading('/tmp/users.txt')
    or die;

=over

=item I<Synopsis>   First validates the file for reading, then creates a IO::File for it.

=item I<Arguments>  file (string)

=item I<Returns>    IO::File object

=back

=head2 validate_file_for_writing

 Genome::Sys->validate_file_for_writing('/tmp/users.txt')
    or die;

=over

=item I<Synopsis>   Checks whether the given file is defined, does not exist, and that the directory it is in is writable

=item I<Arguments>  file (string)

=item I<Returns>    true on success, false on failure

=back

=head2 open_file_for_writing

 Genome::Sys->open_file_for_writing('/tmp/users.txt')
    or die;

=over

=item I<Synopsis>   First validates the file for writing, then creates a IO::File for it.

=item I<Arguments>  file (string)

=item I<Returns>    IO::File object

=back

=head2 copy_file

 Genome::Sys->copy_file($FROM, $TO)
    or ...;

=over

=item I<Synopsis>   Validates the $FROM as a file for reading, and validates the $TO as a file for writing.

=item I<Arguments>  from file (string), to file (string)

=item I<Returns>    true on success, false on failure

=back

=head1 Methods for Directories

=head2 validate_existing_directory

 Genome::Sys->validate_existing_directory('/tmp/users')
    or die;

=over

=item I<Synopsis>   Checks whether the given directory is defined, and is a directory (does not check permissions)

=item I<Arguments>  directory (string)

=item I<Returns>    true on success, false on failure

=back

=head2 open_directory

 Genome::Sys->open_directory('/tmp/users')
    or die;

=over

=item I<Synopsis>   First validates the directory, the creates a IO::Dir handle for it

=item I<Arguments>  IO::Dir (object)

=item I<Returns>    true on success, false on failure

=back

=head2 create_directory

 Genome::Sys->create_directory('/tmp/users')
    or die;

=over

=item I<Synopsis>   Creates the directory with the default permissions 02775

=item I<Arguments>  directory (string)

=item I<Returns>    true on success, false on failure

=back

=head1 Methods for Locking

=head2 lock_resource

Document me!

=head2 unlock_resource

Document me!

=head1 See Also

=head1 Disclaimer

Copyright (C) 2005 - 2008 Washington University Genome Sequencing Center

This module is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

=head1 Author(s)

B<Eddie Belter> I<ebelter@watson.wustl.edu>

=cut

#$HeadURL$
#$Id$
