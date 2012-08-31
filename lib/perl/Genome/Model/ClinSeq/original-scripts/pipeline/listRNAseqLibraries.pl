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
  
  listRNAseqLibraries.pl  --common_names='BRC18,BRC36,BRC38'

INFO

if ($common_names){
  print GREEN, "\n\nAttempting to find RNA-seq libraries for: $common_names\n\n", RESET;
}else{
  print RED, "\n\nRequired parameter missing", RESET;
  print GREEN, "\n\n$usage", RESET;
}

print YELLOW, "LEGEND:\nCN = common name\nSN = sample name\nET = extraction_type\nSCN = sample common name\nTD = tissue description\nCT = cell type\n", RESET;

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
    my $sample_name = $sample->name;
    my $extraction_type = $sample->extraction_type || "UNDEF";
    my $sample_common_name = $sample->common_name || "UNDEF";
    my $tissue_desc = $sample->tissue_desc || "UNDEF";
    my $cell_type = $sample->cell_type || "UNDEF";
    if ($extraction_type eq "rna"){
      print MAGENTA, "\n\t\tSAMPLE\tCN: $common_name\tSN: $sample_name\tET: $extraction_type\tSCN: $sample_common_name\tTD: $tissue_desc\tCT: $cell_type", RESET;
      #Get libraries associated with each sample
      my @libraries = $sample->libraries;
      for my $library (@libraries) {
        print BLUE, "\n\t\t\tLIBRARY\t". $library->name, RESET;

        #Get instrument data ids for the RNA-seq libraries only
        my @instrument_data = Genome::InstrumentData::Solexa->get( library_name => $library->name );
        my $id_count = scalar (@instrument_data);
        my @iid_list;

        #In case the RNA-seq library was captured, get the target region set name(s)
        #Get the target region set name for each instrument data record - and create a list of these for the library
        for my $id (@instrument_data) {
          my $run_name = $id->run_name;
          my $iid = $id->id;
          push(@iid_list, $iid);
          my $trsn = $id->target_region_set_name || 'n/a';
        }
        if ($id_count){
          print BLUE, "\n\t\t\t\tFound $id_count datasets (@iid_list)", RESET;
        }else{
          print BLUE, "\n\t\t\t\tFound $id_count datasets", RESET;
        }

        for my $id (@instrument_data) {
          my $run_name = $id->run_name;
          my $iid = $id->id;
          push(@iid_list, $iid);
          my $trsn = $id->target_region_set_name || 'n/a';
          print BLUE, "\n\t\t\t\t$iid\t$trsn", RESET;

        }



      }
    }
  }
}

print "\n\n";

exit();
