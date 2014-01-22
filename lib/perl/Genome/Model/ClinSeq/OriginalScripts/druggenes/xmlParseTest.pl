#!/usr/bin/env genome-perl -w
#Written by Malachi Griffith

#Parse an XML database file from DrugBank
#

use strict;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use Data::Dumper;
use XML::Simple;

my $infile = '';
my $verbose = '';

GetOptions ('infile=s'=>\$infile, 'verbose=s'=>\$verbose);


my $usage=<<INFO;

  Example usage: 
  
  parseDrugBankXml.pl  --infile=drugbank.xml

  
  Details:
  --infile                    PATH.  XML data file dowloaded from www.drugbank.ca
  --verbose=1                 Print more output while running

INFO

unless ($infile){
  print GREEN, "\n\n$usage\n\n", RESET;
  exit();
}


my $xs1 = XML::Simple->new();

#my $xml = $xs1->XMLin($infile, ForceArray => 1);
#my $xml = $xs1->XMLin($infile);
my $xml = $xs1->XMLin($infile, KeyAttr => ['drugbank-id', 'id'] );

print Dumper $xml;

exit();


