#!/usr/bin/perl -w

use strict;
use warnings;
use Data::Dumper;
use LWP;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
binmode(STDOUT, ":utf8");

my $outdir = '';
my $idfile = '';
my $debug = '0'; #For testing purposes download only a subset of data
my $refresh = '0';

GetOptions ('outdir=s'=>\$outdir, 'idfile=s'=>\$idfile, 'debug=i'=>\$debug, 'refresh=i'=>\$refresh);

my $usage=<<INFO;

  Example usage: 
  
  download_data.pl  --outdir=/gscmnt/sata132/techd/mgriffit/DruggableGenes/KnownDruggable/ClinicalTrialsGov/xml_records_singles/  --idfile=/gscuser/ogriffit/Projects/DruggableGenes/KnownDruggable/ClinicalTrialsGov/record_ids.txt  --debug=0 --refresh=0

  
  Details:
  --outdir                    PATH.  Location to dump individual xml files
  --idfile                    PATH.  File name. Will be used to store a list of clinical trial IDs found by scraping the website
  --debug=1                 Set to 1 to test on a small number of records
  --refresh=0               Set to 1 to force fresh download even when file already exists

INFO

unless ($outdir && $idfile){
  print GREEN, "\n\n$usage\n\n", RESET;
  exit();
}

unless ($outdir =~ /\/$/){
  $outdir .= "/";
}
unless (-e $outdir && -d $outdir){
  print RED, "\n\nOutdir is not valid: $outdir\n\n", RESET;
  exit();
}
 
#Create an LWP object
my $browser = LWP::UserAgent->new;

#Base URL for clinicaltrials.gov
my $base_url = "http://clinicaltrials.gov/";

#Page with all record links for crawlers
my $crawler_url = $base_url."ct2/crawl";

#See if this page exists
my $response = $browser->get($crawler_url);

#Die if it could not be reached
die "Can't get $crawler_url -- ", $response->status_line unless $response->is_success;

#Grab the actual html content of this page:
my $html = $response->content;

#Grab all the URLs on the page - use URI to infer the full path
my @urls;
while( $html =~ m/<A HREF=\"(.*?)\"/ig ) {    
  #print URI->new_abs( $1, $response->base ) ,"\n";
  my $url_object = URI->new_abs( $1, $response->base );
  my $url = $url_object->as_string;
  if ($url =~ /ct2\/crawl\/\d+/){
    push(@urls, $url);
  }
}

#Now go through all the actual pages and grab the record IDs
$| = 0;
my $url_count = scalar(@urls);
my %records;
print "\n\nFound $url_count pages of URLs to record ID pages";
my $c=0;
foreach my $url (@urls){
  $c++;
  if ($debug){
    if ($c > 1){
      last();
    }
  }

  my %recordsl;
  my $response = $browser->get($url);
  die "Can't get $url -- ", $response->status_line unless $response->is_success;
  my $html = $response->content;
  print "\n\tScraping: $url";

  my @record_ids;
  while( $html =~ m/(NCT\d+)/g ) {
    my $rid = $1;
    $recordsl{$rid}=1;
    $records{$rid}=1;
  }
  my $records_on_page = keys %recordsl;
  my $cum_record_count = keys %records;
  print "\n\t\tFound $records_on_page records on this page ($cum_record_count distinct records so far)";
}
$| = 1;

#Now go through all records and get an xml file for each
my $r = 0;
foreach my $rid (sort keys %records){
  $r++;
  my $outfile = "$outdir"."$rid".".xml";
  my $path = "http://clinicaltrials.gov/show/$rid"."?displayxml=true";
  #Check for existence of the xml record.  If it is missing download it.
  if (-e $outfile && $refresh==0){
    print "\n\t\t\t$r: file already exists";
  }else{
    my $wget_cmd = "wget $path -O $outfile 1>/dev/null 2>/dev/null";
    print "\n\t\t\t$r: $wget_cmd";
    system($wget_cmd);
  }
}

#Dump a complete list of the clinical trial records
open (OUT, ">$idfile") || die "\n\nCould not open idfile for output: $idfile\n\n";
foreach my $rid (sort keys %records){
  print OUT "$rid\n";
}
close(OUT);

print "\n\nPrinted complete list of clinical trial record IDs to: $idfile\n\n";

exit();
