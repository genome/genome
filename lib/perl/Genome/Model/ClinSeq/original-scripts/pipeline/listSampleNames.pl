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
my $model_group = '';
my $extraction_type_filter = '';
my $outfile = '';

GetOptions ('common_names=s'=>\$common_names, 'model_group=i'=>\$model_group, 'extraction_type_filter=s'=>\$extraction_type_filter, 'outfile=s'=>\$outfile);

my $usage=<<INFO;
  Example usage: 
  
  listSampleNames.pl  --common_names='BRC18,BRC36,BRC38'

  listSampleNames.pl  --model_group='44014'

  listSampleNames.pl  --model_group='44014'  --extraction_type_filter='rna'  --outfile=ClinSeq_RNA_Samples.tsv

INFO

if ($common_names){
  print GREEN, "\n\nAttempting to find sample names from common names: $common_names\n\n", RESET;
}elsif ($model_group){
  print GREEN, "\n\nAttempting to find sample names from model group id: $model_group\n\n", RESET;
}else{
  print RED, "\n\nRequired parameter missing", RESET;
  print GREEN, "\n\n$usage", RESET;
}

my @outlines;

if ($common_names){
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
      my $extraction_type = $sample->extraction_type;
      my $sample_common_name = $sample->common_name || "UNDEF";
      my $tissue_desc = $sample->tissue_desc || "UNDEF";
      my $cell_type = $sample->cell_type || "UNDEF";

      #Apply filters
      if ($extraction_type_filter){
        next unless ($extraction_type =~ /$extraction_type_filter/);
      }

      #Print summary line
      print MAGENTA, "\nSAMPLE\tCN: $common_name\tSN: $sample_name\tET: $extraction_type\tSCN: $sample_common_name\tTD: $tissue_desc\tCT: $cell_type", RESET;
      push(@outlines, "$common_name\t$sample_name\t$extraction_type\t$sample_common_name\t$tissue_desc\t$cell_type\n");
    }
  }
}

if ($model_group){
  my $model_group = Genome::ModelGroup->get($model_group);
  my @models= $model_group->models;
  foreach my $model (@models){
    my $model_id = $model->id;
    my $subject = $model->subject;
    my $common_name = $subject->common_name || "NA";
    my @samples = $subject->samples;
    foreach my $sample (@samples){
      my $sample_name = $sample->name;
      my $extraction_type = $sample->extraction_type;
      my $sample_common_name = $sample->common_name || "UNDEF";
      my $tissue_desc = $sample->tissue_desc || "UNDEF";
      my $cell_type = $sample->cell_type || "UNDEF";

      #Apply filters
      if ($extraction_type_filter){
        next unless ($extraction_type =~ /$extraction_type_filter/);
      }

      print MAGENTA, "\nSAMPLE\tCN: $common_name\tSN: $sample_name\tET: $extraction_type\tSCN: $sample_common_name\tTD: $tissue_desc\tCT: $cell_type", RESET;
      push(@outlines, "$common_name\t$sample_name\t$extraction_type\t$sample_common_name\t$tissue_desc\t$cell_type");
    }
  }
}

if ($outfile){
  print BLUE, "\n\nDumping samples to file: $outfile", RESET;
  my @outlines_sort = sort @outlines;
  my $out_string = join("\n", @outlines_sort);
  open (OUT, ">$outfile") || die "\n\nCould not open output file: $outfile\n\n";
  print OUT "$out_string";
  close(OUT);
}

print "\n\n";

exit();
