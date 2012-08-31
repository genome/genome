
use strict;
use warnings;

use Test::More tests => 8;
use DBI;
use File::Temp qw/ tempfile /;

BEGIN 
{
    use_ok("BAP::Utils::ActivityLog");
}

# something more permanent and dynamic....
#my $test_db_file = '/gscuser/josborne/src/pm-bap/BAP/Utils/mgap_activity.db';
my ($tfh, $test_db_file) = tempfile();
# possibly set up a temporary sqlite db here?
if( -f $test_db_file )
{
    unlink( $test_db_file );
}
my $create_dbh = DBI->connect("dbi:SQLite:dbname=$test_db_file",'','',
                           { RaiseError => 1, AutoCommit => 1 });
$create_dbh->do(qq/ create table activity_log (
                       row_id        INTEGER PRIMARY KEY,
                       activity      varchar2(20),
                       sequence_id   number,
                       sequence_name varchar2(25),
                       organism_name varchar2(25),
                       host          varchar2(40),
                       user          varchar2(20),
                       started       TIMESTAMP,
                       finished      TIMESTAMP
                            )/, {});
$create_dbh->do(qq/CREATE INDEX activity_log_finished_idx ON activity_log(finished)/, {} );
$create_dbh->do(qq/CREATE INDEX activity_log_started_idx ON activity_log(started)/, {} );

# always always always set this before logging anything
$ENV{ACTIVITY_LOG} = $test_db_file;
# get test host and test user
my $testhost = undef;
if (exists($ENV{LSB_HOSTS}) )
{
    $testhost = $ENV{LSB_HOSTS};
}
elsif (exists($ENV{HOST}) )
{
    $testhost = $ENV{HOST};
}

my $testuser = undef;
if (exists($ENV{USER}) )
{
    $testuser = $ENV{USER};
}
elsif (exists($ENV{LOGIN}) )
{
    $testuser = $ENV{LOGIN};
}
elsif (exists($ENV{USERNAME}) )
{
    $testuser = $ENV{USERNAME};
}



ok(BAP::Utils::ActivityLog::log_action(10, 'test','test_sequence','test_org' ),
   'logging works');

my $testdbh = DBI->connect("dbi:SQLite:dbname=$test_db_file",'','',
                           { RaiseError => 1, AutoCommit => 1 });

my $testsql = qq(select * from activity_log);
my $teststh = $testdbh->prepare($testsql);
$teststh->execute();
my $testptr = $teststh->fetchrow_hashref();
ok($testptr->{activity} eq 'test', 'activity matches');
ok($testptr->{sequence_id} == 0, 'seq id matches');
ok($testptr->{sequence_name} eq 'test_sequence', 'seq name matches');
ok($testptr->{organism_name} eq 'test_org', 'organism name matches');
ok($testptr->{host} eq $testhost, 'host matches');
ok($testptr->{user} eq $testuser, 'user matches');

$testdbh->do(qq(delete from activity_log where activity = 'test'), {});
unlink($test_db_file);
