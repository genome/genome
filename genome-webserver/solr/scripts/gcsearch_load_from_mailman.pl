#!/usr/bin/env genome-perl

use above "Genome";

use strict;
use warnings;

use Data::Dumper;
use Email::Simple;
use Genome::Sys::LockProxy qw();
use LWP::UserAgent;

STDOUT->autoflush(1);
STDERR->autoflush(1);

my $DEBUG = $ENV{'DEBUG_MAIL_LOADER'};
my $total_emails;

#Just clear the cache for each entry instead of building all the search result views for now
Genome::Search->get()->refresh_cache_on_add(0);


#emails are in: http://gscsmtp.wustl.edu/pipermail/[list]/[year]-[monthname].txt
my $HTTP_PATH = 'http://gscsmtp.wustl.edu/pipermail'; 
my $MONTH_NAMES = month_names();


my $lock_resource = 'gcsearch/mailman_loader';
if ($ENV{UR_DBI_NO_COMMIT}) {
    $lock_resource .= '_dev';
}

my $lock = Genome::Sys::LockProxy->new(
    resource => $lock_resource,
    scope => 'site',
)->lock(max_try => 0);
unless ($lock) {
    die "could not lock, another instance must be running.";
}

main();

$lock->unlock();

exit;

sub main {

    my $startyear = shift(@ARGV);
    my $startmonth = shift(@ARGV);
    my @lists = @ARGV;

    undef $startyear unless $startyear and $startyear =~ /\d{4}/ and $startyear > 1990;
    undef $startmonth unless defined $startmonth and $startmonth = int($startmonth) and $startmonth > 0 and $startmonth < 13;
    $startmonth-- if $startmonth;  #Subtract one to map traditional counting to Perl gmtime/array indexing
    unless (@lists) {
        @lists = map (lc $_, @{list_names()} );
    }

    print "LISTS:\n";
    print Dumper \@lists;
    print "start year: $startyear\nstart month: $startmonth\n" if $DEBUG;

    collect_docs(\@lists, $startyear, $startmonth);

    print "total emails added: $total_emails\n";
}

sub collect_docs {
	
    my ($lists, $startyear, $startmonth) = @_;

    #By default just update the current month's data.
    my @today = gmtime();
    my $currentyear = 1900 + $today[5];
    my $currentmonth = $today[4];
    $startyear ||= $currentyear;
    $startmonth = $currentmonth unless defined $startmonth;

    my $ua = LWP::UserAgent->new();

    for my $list (@$lists) {

        print "list: $list\n" if $DEBUG;
        for my $year ($startyear..$currentyear) {

            print "year: $year\n" if $DEBUG;
            for my $month (($year == $startyear ? $startmonth : 0)..($year == $currentyear ? $currentmonth : 11)) {

                print "month: $month\n" if $DEBUG;
                #process the month index to find the ids for links to the individual messages
                my $summary_response = $ua->get($HTTP_PATH . '/' . $list . '/'  . $year . '-' . $MONTH_NAMES->[$month] . '/date.html');

            	my @message_ids;
            	if($summary_response->is_success) {
            		@message_ids = extract_message_ids( $summary_response->content );
            	} else {
            	    if($summary_response->code eq 404) {
            	        next; #Some mailing lists are low-traffic and don't get messages every month.
            	    } else {
            		   die "Couldn't get message ids for $list for $year " . $MONTH_NAMES->[$month] . ": " . $summary_response->status_line;
            	    }
            	}

                #process the text dump of the months' messages; add them to the index
            	my $response = $ua->get($HTTP_PATH . '/' . $list . '/'  . $year . '-' . $MONTH_NAMES->[$month] . '.txt');

                if ($response->is_success) {
                    my @emails = get_emails( $response->content );                   
                    eval {
                        process_emails( $list, $year, $MONTH_NAMES->[$month], \@emails, \@message_ids);
                    };

                    if ($@) {
                        print "Error! skipped $list $year-$month\n" . $@ . "\n";
                    }
                } else {
                    die "Couldn't process $list for $year " . $MONTH_NAMES->[$month] . ": " . $response->status_line;
                    next;
                }

            }
        }
    }
}

sub extract_message_ids {
	
	my ($month_list) = @_;

	my @message_ids = $month_list =~ /<A HREF="(\d+).html"/gi;

	return @message_ids;
}

sub process_emails {

    my ($list_name, $year, $month_name, $emails, $message_ids) = @_;

    my @sorted_message_ids = sort( { $a cmp $b } @$message_ids); #The monthly summary puts them in ID order, not date order

    my $num_emails = scalar(@$emails);
    my $num_message_ids = scalar(@sorted_message_ids);

    if ($num_emails != $num_message_ids) {
        die "parsed $num_emails emails from file but $num_message_ids html files exist ($list_name, $year-$month_name)";
    }

    my $i = 0;
    my @sys_emails;
    for my $email (@$emails) {
        my $msg_id = $sorted_message_ids[$i];
        my $search_id = join('/', $list_name, $year . '-' . $month_name, $msg_id);

        if ($DEBUG) {
            my $url = join('/', $HTTP_PATH, $search_id . '.html');
            print "$url\n";
        }
        
        $email->header_set('X-Genome-Search-ID', $search_id);
        
        my $sys_email = Genome::Sys::Email->get($search_id);
#  If the observer in G:S:Email works, initialize already happens after get()
#        $sys_email->initialize($email);

        push @sys_emails, $sys_email;
        
        $i++;

    }

    $total_emails += $i;

    return Genome::Search->add(@sys_emails);
}



sub get_emails {

    my ($text) = @_;

    my @email_texts = split(/(?m)(?>^From\s[^:\n]*(?::\d\d){2}.*\n)/,$text);
    shift @email_texts;

    my @emails = map( Email::Simple->new($_), @email_texts );

    #Remove duplicates
    my %seen = ();
    @emails = grep { ! $seen{ $_->header('Message-ID') }++ } @emails;

    return @emails;
}

sub month_names {
    return [
       'January', 'February', 'March',
       'April', 'May', 'June',
       'July', 'August', 'September',
       'October', 'November', 'December'
    ];
}

sub list_names {
    my $ua = LWP::UserAgent->new();
    
    my $response = $ua->get($HTTP_PATH . '/');
    
    if ($response->is_success) {
        my $content = $response->content;
        
        my @lists = $content =~ m!<A HREF="([\w\-]+)/">\1!gi; #Match links to list subdirectories
        
        return \@lists;
    } else {
        die "Couldn't load list of mailing lists from server: " . $response->status_line;
    }
}
