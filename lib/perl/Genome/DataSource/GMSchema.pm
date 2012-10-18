package Genome::DataSource::GMSchema;

use strict;
use warnings;
use Genome;
use Carp;
use File::lockf;


class Genome::DataSource::GMSchema {
    is => ['UR::DataSource::Oracle'],
};

sub table_and_column_names_are_upper_case { 1; }

sub server {
    "dwrac";
}

sub login {
    "mguser";
}

sub auth {
    "mguser_prd";
}

sub owner {
    "MG";
}


sub clone_db_handles_for_child_process {
    my $self = shift;

    Genome::DataSource::GMSchemaOracle->get()->clone_db_handles_for_child_process;
    Genome::DataSource::PGTest->get()->clone_db_handles_for_child_process;

    return $self->SUPER::clone_db_handles_for_child_process;
}

our $THIS_COMMIT_ID = 'not within _sync_database';

# This datasource now commits to both Oracle AND postgres. The postgres commit is
# done within an eval so its result does not in any way affect the Oracle commit.
sub _sync_database {
    my $self = shift;
    my %params = @_;

    local $THIS_COMMIT_ID = UR::Object::Type->autogenerate_new_object_id_uuid();

    # Not disconnecting/forking with no commit on to prevent transactions from being
    # closed, which can cause failures in tests that have multiple commits.
    if ($ENV{UR_DBI_NO_COMMIT}) {
        my $oracle_sync_rv = Genome::DataSource::GMSchemaOracle->_sync_database(@_);
        unless ($oracle_sync_rv) {
            Carp::confess "Could not sync to oracle!";
        }
        return 1;
    }

    $self->pause_db_if_necessary;


    # fork if we don't skip.
    my $skip_postgres = (defined $ENV{GENOME_DB_SKIP_POSTGRES} && -e $ENV{GENOME_DB_SKIP_POSTGRES}); 
    my $use_postgres = !$skip_postgres;

    if ($ENV{GENOME_QUERY_POSTGRES}) {
        Genome::Site::TGI->undo_table_name_patch;
        my %classes = map { $_->class => 1 } @{$params{changed_objects}};
        for my $class (sort keys %classes) {
            my $meta = UR::Object::Type->get($class);
            next unless $meta;
            my @metas = ($meta, $meta->ancestry_class_metas);
            for my $meta (@metas) {
                if ($meta->class_name =~ /::Ghost$/) {
                    my $non_ghost_class = $meta->class_name;
                    $non_ghost_class =~ s/::Ghost$//;
                    my $non_ghost_meta = UR::Object::Type->get($non_ghost_class);
                    if ($non_ghost_meta) {
                        $meta->table_name($non_ghost_meta->table_name);
                    }
                    else {
                        Carp::confess "Could not find meta object for non-ghost class $non_ghost_class!";
                    }
                }
            }
        }
    }
	
	

	my ($pg_signal_reader, $pg_signal_writer);

	my $pid;
	
	if ($use_postgres) {
            pipe($pg_signal_reader, $pg_signal_writer);
            $pid = UR::Context::Process->fork();
	} else {
            $pid = $$;
	}
	
    if ($pid) {
        # we're in the parent, close the reader.
        if ($use_postgres) {
                close($pg_signal_reader);
                $pg_signal_writer->autoflush(1);
        }

		
        my $sync_time_start = Time::HiRes::time();
        my $oracle_sync_rv = Genome::DataSource::GMSchemaOracle->_sync_database(@_);

        my $sync_time_duration = Time::HiRes::time() - $sync_time_start;
        unless ($oracle_sync_rv) {
            Carp::confess "Could not sync to oracle!";
        }
        if ($use_postgres) {
            log_commit_time('oracle',$sync_time_duration);
            waitpid($pid, -1);
			
			print $pg_signal_writer("1\n");
			close($pg_signal_writer);
        }		
        if ($ENV{GENOME_QUERY_POSTGRES}) {
            Genome::Site::TGI->redo_table_name_patch;
        }
        return 1;
    } elsif (defined $pid) {
        # Fork twice so parent (process doing Oracle commit) doesn't wait for child
        # to finish.
        # Ignoring SIG_CHLD prevents "Child process #### reaped" from appearing in logs
        $SIG{CHLD} = 'IGNORE';
        my $second_pid = fork();
        Carp::confess "Can't fork" unless defined $second_pid;
        if ($second_pid) {
            # Using POSIX exit prevents END and DESTROY blocks from executing.
            POSIX::_exit(0);
        }

        # close this to stop us from blocking on the read even when our parent exits.
        close($pg_signal_writer);

        # builds will bomb out unless we tell POE that we forked.
        eval {
                POE::Kernel->has_forked();
        };
    
        # Turtles all the way down... the logging logic can potentially bomb and emit warnings that the user
        # shouldn't see, so eval everything!
        eval { 
            my $stderr = '';;
            local *STDERR;
            open STDERR, '>', \$stderr;
            my $sync_time_start = Time::HiRes::time();
			
            eval {
                $DB::single = 1;
                my $pg_commit_rv;
                my $pg_sync_rv = Genome::DataSource::PGTest->_sync_database(@_);
				
												
                my $pg_signal = <$pg_signal_reader>;
                if (defined $pg_signal) {
                    $pg_commit_rv = Genome::DataSource::PGTest->commit;
                }
            };
            my $sync_time_duration = Time::HiRes::time() - $sync_time_start;
            if ($stderr ne '' || $@) {
                my $error = '';
                $error .= "EXCEPTION:" . $@ if $@;
                $error .= "STDERR: " . $stderr if $stderr;
				print $error, "\n";
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
            log_commit_time('pg',$sync_time_duration);
        };
        POSIX::_exit(0);
    }
    else {
        Carp::confess "Problem forking for postgres commit!";
    }

    return 1;
}

sub log_commit_time {
    my($db_name, $time) = @_;

    # See if this process was started from the commandline
    my @commands;
    if ($INC{'Command/V2.pm'}) {
        push @commands, Command::V2->get('original_command_line true' => 1);
    }
    if ($INC{'Command/V1.pm'}) {
        push @commands, Command::V1->get('original_command_line true' => 1);
    }
    @commands = sort { $a->id cmp $b->id } @commands;
    my $original_cmdline = $commands[0] ? $commands[0]->original_command_line : $0;
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


sub create_log_message {
    my $error = shift;

    require DateTime;
    my $dt = DateTime->now;
    $dt->set_time_zone('America/Chicago');
    my $date = $dt->ymd;
    my $time = $dt->hms;
    my $user = Genome::Sys->username;

    require Sys::Hostname;
    my $host = Sys::Hostname::hostname();

    my $string = join("\n", join(',', $date, $time, $host, $user, $THIS_COMMIT_ID), $error);
    return $string;
}

sub _determine_base_log_pathname {
    require DateTime;
    my $dt = DateTime->now;
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


sub _init_created_dbh {
    my ($self, $dbh) = @_;
    return unless defined $dbh;

    $self->SUPER::_init_created_dbh($dbh);

    $dbh->do('alter session set "_hash_join_enabled"=TRUE');

    # stores program name as "MODULE" and user name as "ACTION"
    $self->set_userenv(
        'dbh' => $dbh,
        'module' => substr(Cwd::abs_path($0), -48, 48)
    ); # our oracle module variable is 48 characters

    return $dbh;
}

sub _get_sequence_name_for_table_and_column {
    my ($self, $table_name, $column_name) = @_;
    if ($table_name =~ /PROCESSING_PROFILE/) {
        return 'PROCESSING_PROFILE_SEQ';
    }
    elsif($table_name =~ /GENOME_MODEL_BUILD/) {
        return 'GENOME_MODEL_EVENT_SEQ';
    }
    elsif($table_name =~ /SOFTWARE_RESULT/) {
        return 'GENOME_MODEL_EVENT_SEQ';
    }
    elsif($table_name =~ /MISC_NOTE/) {
        return 'GENOME_MODEL_SEQ';
    }
    elsif ($table_name =~ /INSTRUMENT_DATA/) {
        return 'SEQ_SEQ';
    }
    elsif ($column_name eq 'ID') {
        return $table_name . '_SEQ';
    }
    elsif($table_name =~ /GSC./) {
        return 'IMPORTED_INSTRUMENT_DATA_SEQ';
    }
    else {
        $self->SUPER::_get_sequence_name_for_table_and_column($table_name, $column_name);
    }
}

sub pause_db_if_necessary {
    my $self = shift;
    return 1 unless -e $ENV{GENOME_DB_PAUSE};

    my @o = grep { ref($_) eq 'UR::DeletedRef' } UR::Context->all_objects_loaded('UR::Object');
    if (@o) {
        print Data::Dumper::Dumper(\@o);
        Carp::confess();
    }

    # Determine what has changed.
    my @changed_objects = (
        UR::Context->all_objects_loaded('UR::Object::Ghost'),
        grep { $_->__changes__ } UR::Context->all_objects_loaded('UR::Object')
        #UR::Util->mapreduce_grep(sub { $_[0]->__changes__ },$self->all_objects_loaded('UR::Object'))
    );

    my @real_changed_objects = grep {UR::Context->resolve_data_source_for_object($_)} @changed_objects;

    return 1 unless (@real_changed_objects);


    print "Database updating has been paused, please wait until updating has been resumed...\n";

    my @data_sources = UR::Context->all_objects_loaded('UR::DataSource::RDBMS');
    for my $ds (@data_sources) {
        $ds->disconnect_default_handle if $ds->has_default_handle;
    }

    while (1) {
        sleep sleep_length();
        last unless -e $ENV{GENOME_DB_PAUSE};
    }

    print "Database updating has been resumed, continuing commit!\n";
    return 1;
}

sub sleep_length {
    return 30;
}


1;

