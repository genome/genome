#!/usr/bin/env genome-perl
#Written by Jason Walker, modified by Malachi Griffith
#Get a list of patient common names from the user.  
#Use the Genome API to list information about each of these patients relating to exome or other capture data sets
  

use strict;
use warnings;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use Data::Dumper;
use above 'Genome';

#Required
my $common_names = '';

GetOptions ('common_names=s'=>\$common_names);

my $usage=<<INFO;
  Example usage: 
  
  listExomeDatasets.pl  --common_names='BRC18,BRC36,BRC38'

INFO

if ($common_names){
  print GREEN, "\n\nAttempting to find exome (and other capture) datasets for: $common_names\n\n", RESET;
}else{
  print RED, "\n\nRequired parameter missing", RESET;
  print GREEN, "\n\n$usage", RESET;
}


#my @common_names = qw/BRC18 BRC36 BRC38/;
my @common_names = split(",", $common_names);

for my $common_name (@common_names) {

  #Get an 'individual object using the patient common name
  print BLUE, "\n\n$common_name", RESET;
  my $individual = Genome::Individual->get(
    common_name => $common_name,
  );
  #Get sample objects associated with the individual object
  my @samples = $individual->samples;
  my $scount = scalar(@samples);
  print BLUE, "\n\tFound $scount samples", RESET;
  
  #Get additional info for each sample 
  for my $sample (@samples) {
    #Display basic sample info
    my $sample_name = $sample->name || "UNDEF";
    my $extraction_type = $sample->extraction_type || "UNDEF";
    my $sample_common_name = $sample->common_name || "UNDEF";
    my $tissue_desc = $sample->tissue_desc || "UNDEF";
    my $cell_type = $sample->cell_type || "UNDEF";
    #print BLUE, "\n\t\tSAMPLE\t". $common_name ."\t". $sample->name ."\t". $sample->extraction_type ."\t". $sample->common_name ."\t". $sample->tissue_desc ."\t". $sample->cell_type, RESET;
    print MAGENTA, "\n\t\tSAMPLE\tCN: $common_name\tSN: $sample_name\tET: $extraction_type\tSCN: $sample_common_name\tTD: $tissue_desc\tCT: $cell_type", RESET;

    #Get libraries associated with each sample
    my @libraries = $sample->libraries;

    for my $library (@libraries) {
      #Get instrument data for each sample
      # or my @instrument_data = $sample->solexa_lanes;
      my @instrument_data = Genome::InstrumentData::Solexa->get( library_name => $library->name, );
      my %by_trsn;

      #Get the target region set name for each instrument data record - and create a list of these for the library
      for my $id (@instrument_data) {
        push @{$by_trsn{$id->target_region_set_name || ''}}, $id;
      }

      #Now foreach target region set name, if it is defined, display the instrument IDs.  WGS data does not have a target region set name, so these will be skipped
      for my $trsn (keys %by_trsn) {
        # Comment out if you want WGS
        #Skip this 
        if ($trsn =~ /^\s*$/) { next; }
        my @trsn_ids = @{$by_trsn{$trsn}};
        print BLUE, "\n\t\t\tLIBRARY/TargetSet\t". $library->name ."\t". $trsn ."\t". scalar(@trsn_ids), RESET;
        for my $trsn_id (@trsn_ids) {
          print BLUE, "\n\t\t\t\tINSTRUMENT_DATA\t". $trsn_id->id, RESET;
        }
      }
    }
  }
}

print "\n\n";

exit();
