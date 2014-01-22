#!/usr/bin/env genome-perl
#Written by Malachi Griffith and Gabe Sanderson

use Term::ANSIColor qw(:constants);
use strict;
use warnings;
use above 'Genome';
use Getopt::Long;

#This script attempts to find cases suitable for complete clinseq analysis (i.e. has WGS, Exome, and RNA-seq data)

#Want to be able to query for cases starting with:
#Single common name: e.g. AML103
#Common names query: e.g. AML%
#Individual id (or list of them)
#Individual name (or list of them)

#Get individuals by query
#my @ids = (2872001431);
#my @individuals = Genome::Individual->get(\@ids);

#my $search_string = 'AML%';
#my $search_string = 'AML103';

my $search_string = '';
my $outfile = '';

GetOptions ('search_string=s'=>\$search_string, 'outfile=s'=>\$outfile);

my $usage=<<INFO;
  Example usage: 
  
  identifyClinseqCandidatesByCommonName.pl  --search_string='AML%'  --outfile='AML_Individuals_Summary.tsv'

INFO

unless ($search_string && $outfile){
  print GREEN, "\n\n$usage\n\n", RESET;
  exit(1);
}


my @individuals = Genome::Individual->get(
    common_name => {operator => 'like', value => $search_string},
    #id => {operator => '>=', value => '100000'},
);

#Sort individuals objects by their common name
my %tmp;
my $c = 0;
for my $individual (@individuals) {
  $c++;
  my $common_name = $individual->common_name;
  $tmp{$c}{common_name} = $common_name; 
  $tmp{$c}{individual} = $individual;
}
my @individuals_sort;
foreach my $c (sort {$tmp{$a}->{common_name} cmp $tmp{$b}->{common_name}} keys %tmp){
  push(@individuals_sort, $tmp{$c}{individual});
}

my $individual_count = scalar @individuals;
print "\n\nFound $individual_count individuals\n\n";

my %cc_common_names;
my %cc_names;
my %cc_ids;

my %mc_common_names;
my %mc_names;
my %mc_ids;

#Keep track of target region set names that do not have a feature list type definition
my %undefined_fls;

#Go through each individual ID found
my %ind;
for my $individual (@individuals_sort) {
  my $individual_name = $individual->name || "Unknown";
  my $individual_id = $individual->id;
  my $common_name = $individual->common_name;
  my @samples = $individual->samples;
  my $sample_count = scalar @samples;
  
  $ind{$individual_id}{name} = $individual_name;
  $ind{$individual_id}{common_name} = $common_name;
  $ind{$individual_id}{sample_count} = $sample_count;

  print "\n\n$common_name - $individual_name ($individual_id) has $sample_count samples";

  my ($tumor_wgs, $normal_wgs, $tumor_exome, $normal_exome, $tumor_rna, $normal_rna) = (0,0,0,0,0,0);

  for my $sample (@samples) {
    my $sample_name = $sample->name;
    my $sample_id = $sample->id;
    my $sample_common_name = $sample->common_name || "Unknown";
    my $sample_extraction_type = $sample->extraction_type || "Unknown";

    my @instrument_data;
    my @new_instrument_data = $sample->instrument_data;
 
    my @solexa_instrument_data = grep {$_->class eq "Genome::InstrumentData::Solexa"} @new_instrument_data;
    my $sid_count = scalar @solexa_instrument_data;

    push @instrument_data, @solexa_instrument_data;
    unless (@instrument_data) {
      #print "No instrument data for sample . " . $sample->id;
      next;
    }

    my (@exome_data, @rna_data, @wgs_data, @unknown, @other);
    my %trsns;
    my %feature_list_types;
    for my $instrument_data (@instrument_data) {
      my $sample_type = $instrument_data->sample_type;
      if ($sample_type && $sample_type =~ /rna/) {
        push @rna_data, $instrument_data;
      } else {
        my $target_region = $instrument_data->target_region_set_name;
        if ($target_region){ 
          $trsns{$target_region} = 1;
          my $fl = Genome::FeatureList->get(name => $target_region);
          if (not $fl or not $fl->content_type) {
            push @unknown, $instrument_data;
          }elsif ($fl->content_type eq 'exome') {
            push @exome_data, $instrument_data;
          }else {
            push @other, $instrument_data;        
          }      
        }else{   
          push @wgs_data, $instrument_data;      
        }    
      }  
    } 
    my $trsn_list = '';
    my @trsns = sort keys %trsns;
    if (scalar @trsns){
      foreach my $trsn (@trsns){
        my $fl = Genome::FeatureList->get(name => $trsn);
        my $content_type_def = 0;
        if ($fl && $fl->content_type){
          $content_type_def = 1;
        }else{
          $undefined_fls{$trsn}=1;
        }
        $trsn_list .= "$trsn($content_type_def),";
      }
    }

    my $wgs_data_count = scalar(@wgs_data);
    my $exome_data_count = scalar(@exome_data);
    my $rna_data_count = scalar(@rna_data);
    print "\n\tSample: $sample_id $sample_name $sample_common_name $sample_extraction_type ($sid_count lanes of data)\twgs:$wgs_data_count\texome:$exome_data_count\trna:$rna_data_count\t[trsns: $trsn_list]";

    if ($sample_common_name =~ /tumor|met/i && $sample_extraction_type eq 'genomic dna' && $wgs_data_count){
      $tumor_wgs = 1;
    }
    if ($sample_common_name =~ /tumor|met/i && $sample_extraction_type eq 'genomic dna' && $exome_data_count){
      $tumor_exome = 1;
    }
    if ($sample_common_name =~ /tumor|met/i && $sample_extraction_type eq 'rna' && $rna_data_count){
      $tumor_rna = 1;
    }
    if ($sample_common_name =~ /normal|germline/i && $sample_extraction_type eq 'genomic dna' && $wgs_data_count){
      $normal_wgs = 1;
    }
    if ($sample_common_name =~ /normal|germline/i && $sample_extraction_type eq 'genomic dna' && $exome_data_count){
      $normal_exome = 1;
    }
    if ($sample_common_name =~ /normal|germline/i && $sample_extraction_type eq 'rna' && $rna_data_count){
      $normal_rna = 1;
    }
  }

  $ind{$individual_id}{tumor_wgs} = $tumor_wgs;
  $ind{$individual_id}{tumor_exome} = $tumor_exome;
  $ind{$individual_id}{tumor_rna} = $tumor_rna;
  $ind{$individual_id}{normal_wgs} = $normal_wgs;
  $ind{$individual_id}{normal_exome} = $normal_exome;
  $ind{$individual_id}{normal_rna} = $normal_rna;
  my $sum = $tumor_wgs + $tumor_exome + $tumor_rna + $normal_wgs + $normal_exome + $normal_rna;
  $ind{$individual_id}{sum} = $sum;

  #Is this sample a candidate for ClinSeq analysis?
  if ($tumor_wgs && $tumor_exome && $tumor_rna && $normal_wgs && $normal_exome){
    print BLUE, "\n\tClinSeq candidate: YES (tumor_wgs:$tumor_wgs, normal_wgs:$normal_wgs, tumor_exome:$tumor_exome, normal_exome:$normal_exome, tumor_rna:$tumor_rna, normal_rna:$normal_rna) = $sum", RESET;
    $cc_common_names{$common_name} = 1;
    $cc_names{$individual_name} = 1;
    $cc_ids{$individual_id} = 1;
  }elsif ($tumor_wgs && $tumor_exome && $tumor_rna && $normal_wgs){
    print YELLOW, "\n\tClinSeq candidate: MAYBE (tumor_wgs:$tumor_wgs, normal_wgs:$normal_wgs, tumor_exome:$tumor_exome, normal_exome:$normal_exome, tumor_rna:$tumor_rna, normal_rna:$normal_rna) = $sum", RESET;
    $mc_common_names{$common_name} = 1;
    $mc_names{$individual_name} = 1;
    $mc_ids{$individual_id} = 1;
  }else{
    print RED, "\n\tClinSeq candidate: NO (tumor_wgs:$tumor_wgs, normal_wgs:$normal_wgs, tumor_exome:$tumor_exome, normal_exome:$normal_exome, tumor_rna:$tumor_rna, normal_rna:$normal_rna) = $sum", RESET;
  }
}
my $cc_common_names_string = join(",", sort keys %cc_common_names);
my $cc_names_string = join(",", sort keys %cc_names);
my $cc_ids_string = join(",", sort keys %cc_ids);

my $mc_common_names_string = join(",", sort keys %mc_common_names);
my $mc_names_string = join(",", sort keys %mc_names);
my $mc_ids_string = join(",", sort keys %mc_ids);

#Summarize those target region sets that do not have a feature list or feature list type defined
if (scalar keys %undefined_fls){
  my $undefined_fls_string = join("\n", keys %undefined_fls);
  print YELLOW, "\n\nThe following target region set names do not have a feature list or feature list type defined:\n$undefined_fls_string", RESET;
}else{
  print BLUE, "\n\nAll target region set names do have a feature list and feature list type defined", RESET;
}

#Summarize candidate clinseq cases found
if (scalar keys %cc_ids){
  print BLUE, "\n\nCOMPLETE CANDIDATES", RESET;
  print BLUE, "\nComplete candidate common names: $cc_common_names_string", RESET;
  print BLUE, "\nComplete candidate names: $cc_names_string", RESET;
  print BLUE, "\nComplete candidate ids: $cc_ids_string", RESET;
}else{
  print YELLOW, "\n\nCOMPLETE CANDIDATES", RESET;
  print YELLOW, "\nFound 0 Complete ClinSeq candidates", RESET;
}

#Summarize candidate clinseq cases found
if (scalar keys %mc_ids){
  print BLUE, "\n\nPARTIAL CANDIDATES", RESET;
  print BLUE, "\nPartial candidate common names: $mc_common_names_string", RESET;
  print BLUE, "\nPartial candidate names: $mc_names_string", RESET;
  print BLUE, "\nPartial candidate ids: $mc_ids_string", RESET;
}else{
  print YELLOW, "\n\nPARTIAL CANDIDATES", RESET;
  print YELLOW, "\nFound 0 Partial ClinSeq candidates", RESET;
}

#Print a summary of all individuals found
open (OUT, ">$outfile") || die "\n\nCould not open output file: $outfile\n\n";
print OUT "individual_id\tname\tcommon_name\tsample_count\ttumor_wgs\tnormal_wgs\ttumor_exome\tnormal_exome\ttumor_rna\tnormal_rna\tsum\n";
foreach my $individual_id (sort {$a <=> $b} keys %ind){
  print OUT "$individual_id\t$ind{$individual_id}{name}\t$ind{$individual_id}{common_name}\t$ind{$individual_id}{sample_count}\t$ind{$individual_id}{tumor_wgs}\t$ind{$individual_id}{normal_wgs}\t$ind{$individual_id}{tumor_exome}\t$ind{$individual_id}{normal_exome}\t$ind{$individual_id}{tumor_rna}\t$ind{$individual_id}{normal_rna}\t$ind{$individual_id}{sum}\n";
}
close(OUT);

print "\n\n";


