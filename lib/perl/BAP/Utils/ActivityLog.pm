package BAP::Utils::ActivityLog;

use strict;
use warnings;

use BAP::Config;
use DBI;
use DateTime;


# should have just a bare method for logging times for different things.
sub log_action
{
    my ($elapsed_seconds,$step,
        $sequence_name, $organism) = @_;
    my $sequence_id = 0;
    my $db_file = BAP::Config->new()->activity_db_file();
    if(exists($ENV{ACTIVITY_LOG}))
    {
        $db_file = $ENV{ACTIVITY_LOG};
        if( ! -f $db_file )
        {
            print STDERR "in test mode, not logging action\n";
            return 1;
        }
    }
    my $dbh = DBI->connect("dbi:SQLite:dbname=$db_file",'','',
                           {RaiseError => 1, AutoCommit => 1});
    unless(defined($dbh))
    {

        return 0;
    }
    my $sql = <<SQL;
    INSERT INTO activity_log (activity,
                              sequence_id,
                              sequence_name,
                              organism_name,
                              host,
                              user,
                              started,
                              finished)
        VALUES (?,?,?,?,?,?,
                strftime('%s', 'now') - $elapsed_seconds,
                strftime('%s', 'now')
        );
SQL

    my $host = undef;
    my $user = undef;

    if (exists($ENV{LSB_HOSTS}) )
    {
        $host = $ENV{LSB_HOSTS};
    }
    elsif (exists($ENV{HOST}) )
    {
        $host = $ENV{HOST};
    }

    if (exists($ENV{USER}) )
    {
        $user = $ENV{USER};
    }
    elsif (exists($ENV{LOGIN}) )
    {
        $user = $ENV{LOGIN};
    }
    elsif (exists($ENV{USERNAME}) )
    {
        $user = $ENV{USERNAME};
    }

    # a diagnostic
#    print $db_file,"\n";
    $dbh->do($sql, {}, 
             $step, $sequence_id, $sequence_name, 
             $organism,
             $host, $user);
   

    return 1;
}

sub mark_time 
{
    return DateTime->now(time_zone => 'America/Chicago');
}

1;
