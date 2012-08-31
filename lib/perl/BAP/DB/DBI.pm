package BAP::DB::DBI;

use strict;
use warnings;

use base 'Class::DBI::Oracle';

my $db_options = {
                   __PACKAGE__->_default_attributes(),
                   AutoCommit  => 0,
                   LongReadLen => 10000000, 
                 };
                 
our $db_env = 'prod';

my $db_auth = {
               'prod' => {
                          'sid'  => 'DWRAC',
                          'user' => 'mgapuser',
                          'pass' => 'mg_dw',
                      },
               'dev' => {
                         'sid'  => 'DWDEV',
                         'user' => 'mgapuser',
                         'pass' => 'mg_dev',
                     },
           };

__PACKAGE__->_remember_handle('Main');
__PACKAGE__->autoupdate(0);

my $current_dbh;

sub db_Main {

    if (
        defined($current_dbh) &&
        $current_dbh->FETCH('Active') &&
        $current_dbh->ping() 
       ) {
    
        return $current_dbh;
    
    }
    else {
    
        $current_dbh = DBI->connect_cached(
                                          "dbi:Oracle:$db_auth->{$db_env}{'sid'}", 
                                          $db_auth->{$db_env}{'user'},
                                          $db_auth->{$db_env}{'pass'},
                                          $db_options,
                                          );

        return $current_dbh;

    }

}

1;
