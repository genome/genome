
use strict;
use warnings;

use Time::HiRes qw(sleep);
use Test::WWW::Selenium;
use Test::More tests => 9;
use Test::Exception;

use Sys::Hostname;
use POSIX ":sys_wait_h";

my $PID_FILE = "/tmp/selenium.pid";
my $HOST = hostname();

END {
    print "Trying to shutdown selenium\n";
    exit if ! -e $PID_FILE;

    open(my $fh, $PID_FILE);
    while(<$fh>) {
        my $pid = $_;
        print "killing pid $pid\n";
        kill(9, $pid);
    }
    close($fh);
    unlink($PID_FILE);
}

main();
exit;


sub main {

    die "Error: Expect LDAP user/pass as first two arguments" if ! $ARGV[1];

    start_selenium();
    sleep 10; # child proc is up but might not be ready yet

    #parent proc
    test_imp(@ARGV);
}


sub start_selenium {

    print "Trying to start selenium\n";
    my $pid = fork();
    die 'Error: cant fork()' if ! defined $pid;

    if ($pid == 0) {
        # child proc
        my $server_cmd = 'java -jar /gsc/scripts/lib/java/selenium-server.jar'
                        . ' -firefoxProfileTemplate /gscuser/jlolofie/.mozilla/firefox/hjpikpi9.selenium/ &'
                        . ' echo $! >>'
                        . " $PID_FILE";
        exec($server_cmd);
    }
}


sub test_imp {

    my ($user, $pass) = @_;
    chomp($pass);

    #https://sso.gsc.wustl.edu/idp/Authn/UserPassword
    print "Running the tests...\n";
    my ($sel, $try, $err);
    while (++$try <= 3) {

        eval {
            $sel = Test::WWW::Selenium->new(
                host        => "$HOST",
                port        => 4444,
                browser     => "*chrome",
                browser_url => "https://imp.gsc.wustl.edu/"
            );
        };

        if ($@) {
            $err = $@;
            undef $@;
            sleep 10;
        } else {
            $try = 9;
            last;
        }
    }


    $sel->open_ok("/search");
    $sel->set_speed(2000);
    $sel->type_ok("name=j_username", $user);
    $sel->type_ok("name=j_password", $pass);
    $sel->click_ok("submit");
    $sel->wait_for_page_to_load_ok("30000");

    $sel->type_ok("name=query", "brc59");
    $sel->click_ok("id=searchButton");
    $sel->wait_for_page_to_load_ok("30000");
    $sel->is_element_present_ok("link=BRC59.tumor.somatic.pipeline (2869425631)");

}




