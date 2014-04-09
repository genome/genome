package Genome::Site::TGI::UsageLog;

use strict;
use warnings;

use Cwd qw(abs_path);
use DBI qw();
use File::Basename qw(basename);
use Sys::Hostname qw(hostname);

use Genome::Utility::Instrumentation qw();

my $command = -e $0 ? join(' ', abs_path($0), @ARGV) : join(' ', $0, @ARGV);

if (!exists $ENV{GENOME_LOG_USAGE} || $ENV{GENOME_LOG_USAGE}) {
    record_usage(
        recorded_at  => 'now',
        hostname     => hostname(),
        username     => (getpwuid($<))[0],
        genome_path  => genome_path(),
        perl_path    => abs_path($^X),
        command      => $command,
        perl5lib     => $ENV{PERL5LIB},
        git_revision => git_revision(),
    );
}

sub record_usage {
    my %values = @_;

    my $dbh = log_dbi('DBI', 'connect', dsn(), db_username(), db_password(), db_opts());

    my @columns = qw(recorded_at hostname username perl_path genome_path git_revision command perl5lib);
    my $sth = log_dbi($dbh, 'prepare', q{INSERT INTO usage_log (} . join(', ', @columns) . q{) VALUES (} . join(', ', map { '?' } @columns) . q{)});
    log_dbi($sth, 'execute', map { $values{$_} } @columns);

    log_dbi($dbh, 'disconnect');
}

sub dsn { 'dbi:Pg:dbname=genome_usage;host=gms-postgres.gsc.wustl.edu' }
sub db_username { 'genome_usage' }
sub db_password { 'r7Z1QU4ZVez2CR' }
sub db_opts { { AutoCommit => 1, PrintError => 0, RaiseError => 0 } }

sub log_dbi {
    my $receiver = shift;
    my $method = shift;
    my @args = @_;
    my $rv = $receiver->$method(@args);
    if ($rv) {
        Genome::Utility::Instrumentation::inc(join('.', 'genome_usage', $method));
    } else {
        Genome::Utility::Instrumentation::inc(join('.', 'genome_usage', $method . '_failed'));
    }
    return $rv;
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
