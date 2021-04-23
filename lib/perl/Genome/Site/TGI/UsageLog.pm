package Genome::Site::TGI::UsageLog;

use strict;
use warnings;

use Cwd qw(abs_path);
use DBI qw();
use File::Basename qw(basename);
use Sys::Hostname qw(hostname);
use Carp;

require Genome::Config;

my $command = -e $0 ? join(' ', abs_path($0), @ARGV) : join(' ', $0, @ARGV);

sub import {
    if (should_record_usage()) {
        my $hostname  = hostname();
        my $is_blade  = $hostname =~ /blade[0-9-]+\.gsc\.wustl\.edu/;
        my $is_gmt    = index(basename($0), 'gmt')    == 0;
        my $is_genome = index(basename($0), 'genome') == 0;
        if (  !$is_blade
            || $is_blade && ($is_gmt || $is_genome)
        ) {
            record_usage(
                recorded_at  => 'now',
                hostname     => $hostname,
                username     => $<,
                genome_path  => genome_path(),
                perl_path    => abs_path($^X),
                command      => $command,
                perl5lib     => $ENV{PERL5LIB},
                git_revision => git_revision(),
            );
        }
    }
}

my $should_record_usage;
sub should_record_usage {
    return if $should_record_usage++;
    my $log_usage = Genome::Config::get('log_usage');
    return $log_usage eq '' || $log_usage;
}

sub record_usage {
    my %values = @_;

    my $dbh = log_dbi('DBI', 'connect', dsn(), db_username(), db_password(), db_opts());
    if (! $dbh) {
        my $connect_error = $DBI::errstr;
        eval qq( END { print STDERR "Warning: Failed to connect to logging DB when starting up: \Q$connect_error\E\n" } );

    } else {
        my @columns = qw(recorded_at hostname username perl_path genome_path git_revision command perl5lib);
        my $sth = log_dbi($dbh, 'prepare', q{INSERT INTO usage_log (} . join(', ', @columns) . q{) VALUES (} . join(', ', map { '?' } @columns) . q{)});
        log_dbi($sth, 'execute', map { $values{$_} } @columns);

        log_dbi($dbh, 'disconnect');
    }
}

sub dsn { 'dbi:Pg:dbname=genome_usage;host=gms-postgres.gsc.wustl.edu;sslcert=/gscmnt/gc2560/core/gms-database-certs/gms-postgres.crt;sslkey=/gscmnt/gc2560/core/gms-database-certs/gms-postgres.pem;sslrootcert=/gscmnt/gc2560/core/gms-database-certs/gms-postgres.server.crt' }
sub db_username { Genome::Config::get('usage_log_db_username') }
sub db_password { Genome::Config::get('usage_log_db_password') }
sub db_opts { { AutoCommit => 1, PrintError => 0, RaiseError => 0 } }

sub log_dbi {
    my $receiver = shift;
    my $method = shift;
    my @args = @_;

    return $receiver->$method(@args);
}

sub genome_path {
    unless ($INC{'Genome.pm'}) {
        return '';
    }
    return abs_path($INC{'Genome.pm'});
}

sub git_revision {
    my $genome_path = genome_path();
    return '' unless $genome_path;

    my $work_tree = ($genome_path =~ /^(.*)\/lib\/perl\/Genome.pm$/)[0];
    return '' unless $work_tree;

    my $stdout = qx(git --git-dir "$work_tree/.git" rev-parse HEAD 2> /dev/null);
    return '' unless $stdout;

    chomp $stdout;
    return $stdout;
}

1;
