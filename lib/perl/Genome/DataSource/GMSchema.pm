package Genome::DataSource::GMSchema;

use strict;
use warnings;
use Genome;
use Carp;
use File::lockf;
use DBD::Pg;
use IO::Socket;
use List::MoreUtils qw(any);


class Genome::DataSource::GMSchema {
    is => ['UR::DataSource::RDBMSRetriableOperations', 'UR::DataSource::Oracle'],
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

    # This @changed_objects "shortcut" is to avoid an issue with autocommit on
    # the Pg driver where it would cause an unwanted lock.  This was previously
    # in Command::Dispatch::Shell from UR commit 7fd973b.
    my @changed_objects = @{$params{changed_objects}};
    unless (@changed_objects) {
        return 1;
    }

    local $THIS_COMMIT_ID = UR::Object::Type->autogenerate_new_object_id_uuid();

    my $required_pg_version = '2.19.3';

    my $pg_version = $DBD::Pg::VERSION;
    if (($pg_version ne $required_pg_version) && !defined $ENV{'LIMS_PERL'}) {
        $self->error_message("**** INCORRECT POSTGRES DRIVER VERSION ****\n" .
            "You are using a Perl version that includes an incorrect DBD::Pg driver.\n" .
            "You are running $pg_version and need to be running $required_pg_version.\n" .
            "Your sync has been aborted to protect data integrity in the Postgres database.\n" .
            "Please be sure you are using 'genome-perl' and not /gsc/bin/perl.\n\n\n" .
            "This event has been logged with apipe; if you are unsure of why you received this message\n" .
            "open an apipe-support ticket with the date/time of occurrence and we will assist you.\n");
        log_error($self->error_message);
        die $self->error_message;
    }


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

    # Attempt to get a meta db handle first.  This way, if the meta db doesn't exist,
    # then we'll create it before forking.  Otherwise, it'll attempt to create in the fork,
    # fail, and then we'll lose the synced data.
    my $meta_dbh = Genome::DataSource::Meta->get_default_handle();

    if ($ENV{GENOME_QUERY_POSTGRES}) {
        Genome::Site::TGI->undo_table_name_patch;
        my %classes = map { $_->class => 1 } @changed_objects;
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


    my ($parent_oracle_control_sock, $child_pg_control_sock);

    my $pid;

    if ($use_postgres) {
        ($parent_oracle_control_sock, $child_pg_control_sock) = IO::Socket->socketpair(AF_UNIX, SOCK_STREAM, PF_UNSPEC);
        unless ($parent_oracle_control_sock && $child_pg_control_sock) {
            die $self->error_message("Uh-oh. Couldn't prepare oracle/postgres sync control socket pair.");
        }

        $pid = UR::Context::Process->fork();
    } else {
        $pid = $$;
    }

    if ($pid) {
        my $sync_time_start = Time::HiRes::time();
        my $oracle_sync_rv = Genome::DataSource::GMSchemaOracle->_sync_database(@_);
        my $sync_time_duration = Time::HiRes::time() - $sync_time_start;
        unless ($oracle_sync_rv) {
            Carp::confess "Could not sync to oracle!";
        }

        if ($use_postgres) {
            log_commit_time('oracle', $sync_time_duration);
            close $parent_oracle_control_sock;

            print $child_pg_control_sock "1\n";
            eval {
                local $SIG{'ALRM'} = sub {
                    log_error('Timed out waiting for Postgres to sync! Databases have possibly diverged.');
                    die; # to exit eval
                };
                alarm 30;
                my $pg_signal = <$child_pg_control_sock>;
                alarm 0;
            };
        }

        if ($ENV{GENOME_QUERY_POSTGRES}) {
            Genome::Site::TGI->redo_table_name_patch;
        }

        return 1;
    }
    elsif (defined $pid) {
        # close this to stop us from blocking on the read even when our parent exits.
        close $child_pg_control_sock;

        # Fork twice so parent (process doing Oracle commit) doesn't wait for
        # child to finish.  Ignoring SIG_CHLD prevents "Child process ####
        # reaped" from appearing in logs.
        $SIG{CHLD} = 'IGNORE';
        my $second_pid = fork();
        Carp::confess "Can't fork" unless defined $second_pid;
        if ($second_pid) {
            POSIX::_exit(0); # avoids END and DESTROY blocks
        }

        # builds will bomb out unless we tell POE that we forked.
        eval { POE::Kernel->has_forked() };

        # Turtles all the way down... the logging logic can potentially bomb
        # and emit warnings that the user shouldn't see, so eval everything!
        eval {
            my $stderr = '';;
            local *STDERR;
            open STDERR, '>', \$stderr;
            my $sync_time_start = Time::HiRes::time();

            eval {
                my $pg_commit_rv;
                my $pg_sync_rv = Genome::DataSource::PGTest->_sync_database(@_);

                my $pg_signal = <$parent_oracle_control_sock>;
                if (defined $pg_signal) {
                    $pg_commit_rv = Genome::DataSource::PGTest->SUPER::commit;
                }
            };
            my $sync_time_duration = Time::HiRes::time() - $sync_time_start;
            if ($stderr ne '' || $@) {
                my $error = '';
                $error .= "EXCEPTION:" . $@ if $@;
                $error .= "STDERR: " . $stderr if $stderr;
                print "***** POSTGRES SYNC ERROR *****\n";
                print "The postgres write failed, for the following reason.  However, the Oracle sync succeeded.\n\n";
                print "This is a non-fatal error and only affects Postgres testing. Your task completed successfully.\n";
                print $error, "\n";
                log_error($error);
            }
            log_commit_time('pg', $sync_time_duration);
        };
        print $parent_oracle_control_sock "1\n";
        POSIX::_exit(0); # avoids END and DESTROY blocks
    }
    else {
        Carp::confess "Problem forking for postgres commit!";
    }

    return 1;
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
    return 1 unless $ENV{GENOME_DB_PAUSE} and -e $ENV{GENOME_DB_PAUSE};

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
        last unless $ENV{GENOME_DB_PAUSE} and -e $ENV{GENOME_DB_PAUSE};
    }

    print "Database updating has been resumed, continuing commit!\n";
    return 1;
}

sub sleep_length {
    return 30;
}

my @retriable_operations = (
    qr(ORA-25408), # can not safely replay call
);
sub should_retry_operation_after_error {
    my($self, $sql, $dbi_errstr) = @_;
    return any { $dbi_errstr =~ /$_/ } @retriable_operations;
}

sub create_iterator_closure_for_rule {
    my $self = shift;
    my @create_iter_params = @_;

    $self->_retriable_operation(sub {
        $self->pause_db_queries_if_necessary();

        $self->SUPER::create_iterator_closure_for_rule(@create_iter_params);
    });
}

sub create_dbh {
    my $self = shift;
    my @create_dbh_params = @_;

    $self->_retriable_operation(sub {
        $self->pause_db_queries_if_necessary();

        $self->SUPER::create_dbh(@create_dbh_params);
    });
}

sub pause_db_queries_if_necessary {
    my $self = shift;
    return 1 unless $ENV{GENOME_DB_QUERY_PAUSE} and -e $ENV{GENOME_DB_QUERY_PAUSE};

    print "Database querying has been paused; disconnecting db handles.  Please wait until the query pause is released...\n";

    my @data_sources = UR::Context->all_objects_loaded('UR::DataSource::RDBMS');
    for my $ds (@data_sources) {
        $ds->disconnect_default_handle if $ds->has_default_handle;
    }

    while (1) {
        sleep 30;
        last unless $ENV{GENOME_DB_QUERY_PAUSE} and -e $ENV{GENOME_DB_QUERY_PAUSE};
    }

    print "Database updating has been resumed, continuing query!\n";
    return 1;
}

# A list of the old GM schema tables that Genome::Model should ignore
my @OLD_GM_TABLES = qw(
ALL_ALLELE_TYPE
ANALYSIS_METHOD
CHROMOSOME
COLLABORATOR
COLLABORATOR_SAMPLE
CONSERVATION_SCORE
EGI_TYPE
EXTERNAL_GENE_ID
GENE
GENE_EXPRESSION
GENE_GENE_EXPRESSION
GENE_GENOMIC_FEATURE
GENE_GENOTYPE
GENOME_MODEL_1
GENOME_UPDATE_HISTORY
GENOMIC_FEATURE
GENOTYPE
GENOTYPE_VARIATION
GE_DETECTION
GE_TECH_TYPE
GF_FEATURE_TYPE
GO_XREF
GROUP_INFO
GROUP_TYPES
HISTOLOGY
IPRO_GENE_TRANSCRIPT_XREF
IPRO_RESULTS
MAF
META_GROUP
META_GROUP_STATUS
PF_FEATURE_TYPE
PP_DETECTION_SOFT
PP_MAPPING_REFERENCE
PP_TECH_TYPE
PROCESS_PROFILE
PROCESS_PROFILE_SOURCE
PROTEIN
PROTEIN_FEATURE
PROTEIN_VARIATION
PROTEIN_VARIATION_SCORE
PVS_SCORE_TYPE
PVS_SOFTWARE
READ_GROUP
READ_GROUP_GENOTYPE
READ_GROUP_GENOTYPE_OLD
READ_GROUP_INFO
READ_INFO
REPEAT_INFO
RGGI_INFO_TYPE
RGG_INFO
RGG_INFO_OLD
RI_REPEAT_TYPE
SAMPLE
SAMPLE_GENE
SAMPLE_GENE_EXPRESSION
SAMPLE_GENOTYPE
SAMPLE_GROUP_INFO
SAMPLE_HISTOLOGY
SAMPLE_SUBTYPE
SAMPLE_TISSUE
SAMPLE_TYPE
SEQUENCE_DIFF
SEQUENCE_DIFF_EVAL
SEQUENCE_DIFF_PART
SGIAM_SCORE_TYPE
SGI_ANALYSIS_METHOD
SG_INFO_TYPE
STRAND
SUBMITTER
SUBMITTER_METHOD
TERM
TISSUE
TRANSCRIPT
TRANSCRIPT_SOURCE
TRANSCRIPT_STATUS
TRANSCRIPT_SUB_STRUCTURE
TRANSCRIPT_VARIATION
TSS_STRUCTURE_TYPE
TV_OLD
TV_TYPE
VARIANT_REVIEW_DETAIL_OLD
VARIANT_REVIEW_LIST_FILTER
VARIANT_REVIEW_LIST_MEMBER
VARIATION
VARIATION_FREQUENCY
VARIATION_GROUP
VARIATION_INSTANCE
VARIATION_ORIG
VARIATION_ORIG_VARIATION_TYPE
VARIATION_SCORE
VS_SCORE_TYPE
VS_SOFTWARE
);

sub _ignore_table {
    my($self,$table_name) = @_;

    return 1 if $self->SUPER::_ignore_table($table_name);

    return scalar(grep { $_ eq $table_name } @OLD_GM_TABLES);
}

1;

