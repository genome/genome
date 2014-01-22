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
  
  listRNAseqBuildsByCommonName.pl  --common_names='BRC18,BRC36,BRC38'

INFO

if ($common_names){
  print GREEN, "\n\nAttempting to find exome (and other capture) datasets for: $common_names\n\n", RESET;
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

  #my @samples = $individual->samples;
  #my @models = Genome::Model::RnaSeq->get("subject"=>\@samples);
  #my @models = Genome::Model::RnaSeq->get("subject.source.common_name"=>$common_name);
  
  my @models = Genome::Model::RnaSeq->get("subject.patient.common_name"=>$common_name);
  
  foreach my $model (@models){
    my $model_id = $model->id;
    my $model_name = $model->name;
    my $model_subject = $model->subject;
    my $model_subject_name = $model_subject->name;
    my $processing_profile = $model->processing_profile;
    my $processing_profile_name = $processing_profile->name;

    print BLUE, "\n\tModel. ID: $model_id\tSubject: $model_subject_name\tModel Name: $model_name", RESET;
    print BLUE, "\n\t\tPP: $processing_profile_name", RESET;

  }


  #print Dumper @models;

}

print "\n\n";

exit();
