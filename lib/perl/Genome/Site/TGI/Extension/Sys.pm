package Genome::Sys;

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
require Genome::Utility::Text; #TODO remove, unused
use Sys::Hostname;
use File::Find;
use Params::Validate qw(:types);
use IO::Socket;
use Sys::Hostname qw(hostname);
#use Archive::Extract;

require MIME::Lite;

use Genome::Sys::LSF::bsub qw();

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
    return Genome::Sys::LSF::bsub::bsub(@_);
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

sub bsub_and_wait_for_completion {
    my($class, %bsub_args) = @_;

    my $special_args;
    ($special_args, %bsub_args) = _bsub_and_wait_for_completion__validate_args(%bsub_args);
    my($commands, $on_submit_cb, $on_complete_cb) = @$special_args{'commands','on_submit','on_complete'};

    my %seq_to_job_id;
    for(my $seq = 0; $seq < @$commands; $seq++) {
        my $cmd = $commands->[$seq];
        my %this_cmd_bsub_args;
        if (ref($cmd) and ref($cmd) eq 'HASH') {
            %this_cmd_bsub_args = %$cmd;
        } else {
            %this_cmd_bsub_args = ( cmd => $cmd );
        }

        my $job_id = $class->bsub(
                                  %bsub_args,
                                  %this_cmd_bsub_args,
                                );
        $on_submit_cb->($seq, $job_id) if $on_submit_cb;
        $seq_to_job_id{$seq} = $job_id;
    }

    $class->_bsub_and_wait_for_completion__wait_on_jobs_with_lsf_dependency([values %seq_to_job_id], { queue => $bsub_args{queue} || 'short' } );

    # immediately after being reaped in _bsub_and_wait_for_completion_wait_on_jobs(), bjobs
    # reports their status as still "RUN".  Instead, we go into the traditional polling until
    # LSF updates their status to EXIT or DONE
    my %job_statuses = $class->wait_for_lsf_jobs(values %seq_to_job_id);
    my @statuses = map { $job_statuses{$_} }
                   map { $seq_to_job_id{$_} }
                   ( 0 .. $#$commands );
    return @statuses;
}

sub _bsub_and_wait_for_completion__validate_args {
    my %bsub_args = @_;

    my %special_args;
    $special_args{commands} = my $commands = delete $bsub_args{cmds};
    unless ($commands and ref($commands) and ref($commands) eq 'ARRAY') {
        Carp::croak(q(arg 'commands' is required and must be an arrayref of commands'));
    }
    if (exists $bsub_args{cmd}) {
        Carp::croak(q(arg 'cmd' is not allowed, use 'cmds' instead, an arrayref of commands));
    }
    $special_args{on_submit} = delete $bsub_args{on_submit};
    $special_args{on_complete} = delete $bsub_args{on_complete};
    return(\%special_args, %bsub_args);
}

sub _bsub_and_wait_for_completion__wait_on_jobs_with_lsf_dependency {
    my ($class, $job_ids, $bsub_args) = @_;
    local $ENV{LSF_NIOS_JOBSTATUS_INTERVAL} = 1;
    $class->bsub(%$bsub_args,
                 wait_for_completion => 1,
                 depend_on => (join "&&", map { "ended($_)" } @$job_ids) ,
                 cmd => "echo DONE");
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
        sleep 10 unless $all_jobs_complete;
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
            allow_zero_size_input_files => 1,
        );

        unless($rv) {
            die('Failed to cat ' . @next);
        }

        $mode = '>>';
        $skip_if_output_is_present = 0;
    }

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
    my $LSB_JOBID = $ENV{LSB_JOBID};
    my $bjobs_cmd = qx(which bjobs);
    if ($bjobs_cmd && $LSB_JOBID) {
        chomp $bjobs_cmd;
        my $bjobs = qx($bjobs_cmd -l $LSB_JOBID);
        my ($bjobs_mem_limit, $bjobs_mem_limit_units) = $bjobs =~ /MEMLIMIT\s+(\d+\.?\d*)\s+([TGMK])\b/;
        if ($bjobs_mem_limit_units eq 'T') {
            return int($bjobs_mem_limit * 1024 * 1024 * 1024);
        }
        elsif ($bjobs_mem_limit_units eq 'G') {
            return int($bjobs_mem_limit * 1024 * 1024);
        }
        elsif ($bjobs_mem_limit_units eq 'M') {
            return int($bjobs_mem_limit * 1024);
        }
        elsif ($bjobs_mem_limit_units eq 'K') {
            return int($bjobs_mem_limit);
        }
        else {
            # failed to parse
        }
    }
    return undef;
}


# detect maximum available memory
# would probably be better to not do this this way but it's a start
sub mem_limit_kb {
    my $class = shift;

    # get LSF memory limit
    if ($ENV{LSB_JOBID}) {
        my $mem_limit = $class->get_mem_limit_from_bjobs;
        return $mem_limit if $mem_limit;
    }

    # get physical total memory
    if (-e '/proc/meminfo') {
        my $mem_total = $class->get_mem_total_from_proc;
        return $mem_total if $mem_total;
    }

    return undef;
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

=head1 Disclaimer

Copyright (C) 2005 - 2008 Washington University Genome Sequencing Center

This module is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY or the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

=head1 Author(s)

B<Eddie Belter> I<ebelter@watson.wustl.edu>

=cut

#$HeadURL$
#$Id$
