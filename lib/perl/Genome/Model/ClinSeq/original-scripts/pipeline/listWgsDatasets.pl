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
my $sample_names = '';
my $report_only = '';

GetOptions ('common_names=s'=>\$common_names, 'sample_names=s'=>\$sample_names, 'report_only=s'=>\$report_only);

my $usage=<<INFO;
  Example usage: 
  
  listWgsDatasets.pl  --common_names='BRC18,BRC36,BRC38'

  OR

  listWgsDatasets.pl  --sample_names='H_MF-Vaco432-1106984,H_MF-Vaco432-50.1R-1106985'

  OR 

  listWgsDatasets.pl  --sample_names='H_MF-Vaco432-1106984,H_MF-Vaco432-50.1R-1106985'  --report_only=1


INFO


if ($common_names){
  unless ($report_only){
    print GREEN, "\n\nAttempting to find WGS datasets for: $common_names\n\n", RESET;
  }
}elsif($sample_names){
  unless ($report_only){
    print GREEN, "\n\nAttempt to find WGS datasets for: $sample_names", RESET;
  }
}else{
  print RED, "\n\nRequired parameter missing", RESET;
  print GREEN, "\n\n$usage", RESET;
}


if ($common_names){
  my @common_names = split(",", $common_names);

  for my $common_name (@common_names) {

    #Get an 'individual object using the patient common name
    unless ($report_only){
      print BLUE, "\n\n$common_name", RESET;
    }
    my $individual = Genome::Individual->get(
      common_name => $common_name,
    );
    #Get sample objects associated with the individual object
    my @samples = $individual->samples;
    my $scount = scalar(@samples);
    unless ($report_only){
      print BLUE, "\n\tFound $scount samples", RESET;
    }
    #Get additional info for each sample 
    for my $sample (@samples) {
      #Display basic sample info
      my $sample_name = $sample->name || "UNDEF";
      my $extraction_type = $sample->extraction_type || "UNDEF";
      my $sample_common_name = $sample->common_name || "UNDEF";
      my $tissue_desc = $sample->tissue_desc || "UNDEF";
      my $cell_type = $sample->cell_type || "UNDEF";
      unless ($report_only){
        print MAGENTA, "\n\t\tSAMPLE\tCN: $common_name\tSN: $sample_name\tET: $extraction_type\tSCN: $sample_common_name\tTD: $tissue_desc\tCT: $cell_type", RESET;
      }
      #Get libraries associated with each sample
      my @libraries = $sample->libraries;

      for my $library (@libraries) {
        #Get instrument data for each sample
        # or my @instrument_data = $sample->solexa_lanes;
        my @instrument_data = Genome::InstrumentData::Solexa->get( library_name => $library->name, );
        my @id_list;
        for my $id (@instrument_data) {
          #If a target region set name is defined, this is not WGS data, do not list
          unless ($id->target_region_set_name){
            if ($extraction_type eq "genomic dna"){
              push (@id_list, $id->id);
            }
          }
        }
        
        #If a target region set name is defined, this is not WGS data, do not list
        unless ($report_only){
          print BLUE, "\n\t\t\tLIBRARY\t". $library->name, RESET;
          print BLUE, "\n\t\t\t\tINSTRUMENT_DATA\t@id_list", RESET;
        }
      }
    }
  }
}


if ($sample_names){
  my @sample_names = split(",", $sample_names);

  for my $sample_name (@sample_names) {

    my $sample = Genome::Sample->get( name => $sample_name);
    my $patient = $sample->patient;
    my $patient_name = $patient->name;
    my $patient_common_name = $patient->common_name || "UNDEF";

    #Display basic sample info
    my $sample_name = $sample->name || "UNDEF";
    my $extraction_type = $sample->extraction_type || "UNDEF";
    my $sample_common_name = $sample->common_name || "UNDEF";
    my $tissue_desc = $sample->tissue_desc || "UNDEF";
    my $cell_type = $sample->cell_type || "UNDEF";
    unless ($report_only){
      print MAGENTA, "\nSAMPLE\tSN: $sample_name\tET: $extraction_type\tSCN: $sample_common_name\tTD: $tissue_desc\tCT: $cell_type", RESET;
    }
    #Get libraries associated with each sample
    my @libraries = $sample->libraries;
    my @grand_id_list;
    for my $library (@libraries) {
      #Get instrument data for each sample
      # or my @instrument_data = $sample->solexa_lanes;
      my @instrument_data = Genome::InstrumentData::Solexa->get( library_name => $library->name, );
      my @id_list;
      for my $id (@instrument_data) {
        #If a target region set name is defined, this is not WGS data, do not list
        unless ($id->target_region_set_name){
          if ($extraction_type eq "genomic dna"){
            push (@id_list, $id->id);
            push (@grand_id_list, $id->id);
          }
        }
      }
      my $ii_count = scalar(@id_list);
      if ($ii_count > 0){
        unless ($report_only){
          print BLUE, "\n\tLIBRARY\t". $library->name, RESET;
          print BLUE, "\n\t\tINSTRUMENT_DATA\t@id_list", RESET;
        }
      }
    }
    my $grand_ii_count = scalar(@grand_id_list);
    #If a target region set name is defined, this is not WGS data, do not list
    if ($grand_ii_count > 0){
      if ($report_only){
        print "$patient_common_name\t$patient_name\t$sample_common_name\t$sample_name\t$extraction_type\t$sample_common_name\t$tissue_desc\t$cell_type\n";
      }
    }
  }
}

print "\n\n";

exit();
