#!/usr/bin/env genome-perl -w
#Written by Malachi Griffith

use strict;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use Data::Dumper;
use above 'Genome';
use Genome::Model::ClinSeq::Util qw(:all);

my $infile = '';
my $option = '';

GetOptions ('infile=s'=>\$infile, 'option=i'=>\$option);

unless ($infile && $option){
  print RED, "\nRequired input parameter(s) missing or incorrect format!\n\n", RESET;
  print GREEN, "\nExample usage: summarizeSnvIndelFile.pl  --infile=/gscmnt/ams1180/info/model_data/2875816457/build111571998/effects/snvs.hq.tier1.v1.annotated  --option=1\n\n", RESET;
  print GREEN, "\nExample usage: summarizeSnvIndelFile.pl  --infile=/gscmnt/ams1180/info/model_data/2875816457/build111571998/effects/indels.hq.tier1.v1.annotated  --option=2\n\n", RESET;
  print GREEN, "\n\nOptions:", RESET;
  print GREEN, "\n1 -> SNVs.  Missense/nonsense/splice.  Unique on coord with gene names and amino acid effects concatenated for multiple transcripts", RESET;
  print GREEN, "\n2 -> SNVs.  Missense/nonsense/splice.  Unique on coord with gene names and amino acid effects and sift/polyphen results concatenated for multiple transcripts", RESET;
  print GREEN, "\n3 -> SNVs.  All.  Unique on coord with gene names and amino acid effects concatenated for multiple transcripts", RESET;
  print GREEN, "\n4 -> Indels.  Gene affecting.  Unique on coord with gene names, amino acid effects and indel types concatenated for multiple transcripts", RESET;

  print GREEN, "\n\n", RESET;
  exit();
}

unless (-e $infile){
  print RED, "\n\nCould not find file: $infile\n\n", RESET;
  exit();
}


#Option 1.  Produce a simple TSV with the location of each SNV, and the gene symbols associated with that location.  Limit to missense/nonsense/splice mutations only
#Make unique on coordinate position and supply a comma separated list of genes/transcripts reported for that coord
#Similarly concatenate the amino acid effects reported
if ($option == 1){
  open (IN, "$infile") || die "\n\nCould not open input file: $infile\n\n";
  my %data;
  my $d = 0;
  while(<IN>){
    chomp($_);
    unless ($_ =~ /missense|nonsense|splice/){
      next();
    }
    $d++;
    my @line = split("\t", $_);
    my $coord = "$line[0]:$line[1]";
    my $symbol = $line[6];
    my $aa_change = $line[15];

    $data{$coord}{order} = $d;
    if ($data{$coord}{symbols}){
      my $symbol_ref = $data{$coord}{symbols};
      $symbol_ref->{$symbol}=1;
      my $aa_ref = $data{$coord}{aa_changes};
      $aa_ref->{$aa_change}=1;
    }else{
      my %sym;
      my %aa;
      $sym{$symbol}=1;
      $aa{$aa_change}=1;
      $data{$coord}{symbols} = \%sym;
      $data{$coord}{aa_changes} = \%aa;
    }
  }
  close(IN);
  foreach my $pos (sort {$data{$a}->{order} <=> $data{$b}->{order}} keys %data){
    my %symbols = %{$data{$pos}{symbols}};
    my @symbols = keys %symbols;
    my $symbol_string = join(", ", @symbols);
    my %aa_changes = %{$data{$pos}{aa_changes}};
    my @aa_changes = sort keys %aa_changes;
    my $aa_changes_string = join(", ", @aa_changes);
    print "$pos\t$symbol_string\t$aa_changes_string\n";
  }
}


#Option 2.  Produce a simple TSV with the location of each SNV, and the gene symbols associated with that location.  Limit to missense/nonsense/splice mutations only
#Make unique on coordinate position and supply a comma separated list of genes/transcripts reported for that coord
#Similarly concatenate the amino acid effects reported
#Similarly concatenate SIFT/polyphen results
if ($option == 2){
  open (IN, "$infile") || die "\n\nCould not open input file: $infile\n\n";
  my %data;
  my $d = 0;
  while(<IN>){
    chomp($_);
    unless ($_ =~ /missense|nonsense|splice/){
      next();
    }
    $d++;
    my @line = split("\t", $_);
    my $coord = "$line[0]:$line[1]";
    my $symbol = $line[6];
    my $aa_change = $line[15];
    my $sift = $line[21];
    my $polyphen = $line[22];

    $data{$coord}{order} = $d;
    if ($data{$coord}{symbols}){
      my $symbol_ref = $data{$coord}{symbols};
      $symbol_ref->{$symbol}=1;
      my $aa_ref = $data{$coord}{aa_changes};
      $aa_ref->{$aa_change}=1;
      my $sift_ref = $data{$coord}{sift};
      $sift_ref->{$sift}=1;
      my $polyphen_ref = $data{$coord}{polyphen};
      $polyphen_ref->{$polyphen}=1;
    }else{
      my %sym;
      my %aa;
      my %sift;
      my %polyphen;
      $sym{$symbol}=1;
      $aa{$aa_change}=1;
      $sift{$sift}=1;
      $polyphen{$polyphen}=1;
      $data{$coord}{symbols} = \%sym;
      $data{$coord}{aa_changes} = \%aa;
      $data{$coord}{sift} = \%sift;
      $data{$coord}{polyphen} = \%polyphen;
    }
  }
  close(IN);
  foreach my $pos (sort {$data{$a}->{order} <=> $data{$b}->{order}} keys %data){
    my %symbols = %{$data{$pos}{symbols}};
    my @symbols = keys %symbols;
    my $symbol_string = join(", ", @symbols);
    my %aa_changes = %{$data{$pos}{aa_changes}};
    my @aa_changes = keys %aa_changes;
    my $aa_changes_string = join(", ", @aa_changes);
    my %sift = %{$data{$pos}{sift}};
    my @sift = keys %sift;
    my $sift_string = join(", ", @sift);
    my %polyphen = %{$data{$pos}{polyphen}};
    my @polyphen = keys %polyphen;
    my $polyphen_string = join(", ", @polyphen);
    print "$pos\t$symbol_string\t$aa_changes_string\t$sift_string\t$polyphen_string\n";
  }
}


#Option 3.  Produce a simple TSV with the location of each SNV, and the gene symbols associated with that location.  No filtering of SNV mutations
#Make unique on coordinate position and supply a comma separated list of genes/transcripts reported for that coord
#Similarly concatenate the amino acid effects reported
if ($option == 3){
  open (IN, "$infile") || die "\n\nCould not open input file: $infile\n\n";
  my %data;
  my $d = 0;
  while(<IN>){
    chomp($_);
    $d++;
    my @line = split("\t", $_);
    my $coord = "$line[0]:$line[1]";
    my $symbol = $line[6];
    my $aa_change = $line[15];

    $data{$coord}{order} = $d;
    if ($data{$coord}{symbols}){
      my $symbol_ref = $data{$coord}{symbols};
      $symbol_ref->{$symbol}=1;
      my $aa_ref = $data{$coord}{aa_changes};
      $aa_ref->{$aa_change}=1;
    }else{
      my %sym;
      my %aa;
      $sym{$symbol}=1;
      $aa{$aa_change}=1;
      $data{$coord}{symbols} = \%sym;
      $data{$coord}{aa_changes} = \%aa;
    }
  }
  close(IN);
  foreach my $pos (sort {$data{$a}->{order} <=> $data{$b}->{order}} keys %data){
    my %symbols = %{$data{$pos}{symbols}};
    my @symbols = keys %symbols;
    my $symbol_string = join(", ", @symbols);
    my %aa_changes = %{$data{$pos}{aa_changes}};
    my @aa_changes = sort keys %aa_changes;
    my $aa_changes_string = join(", ", @aa_changes);
    print "$pos\t$symbol_string\t$aa_changes_string\n";
  }
}


#Option 4.  Produce a simple TSV with the location of each Indel, and the gene symbols associated with that location.  Ignore: 5prime, 3prime, intronic indels
#Make unique on coordinate position and supply a comma separated list of genes/transcripts reported for that coord
#Similarly concatenate the amino acid effects reported, and the type of mutation
if ($option == 4){
  open (IN, "$infile") || die "\n\nCould not open input file: $infile\n\n";
  my %data;
  my $d = 0;
  while(<IN>){
    chomp($_);
    if ($_ =~ /3_prime_flanking_region|3_prime_untranslated_region|5_prime_flanking_region|5_prime_untranslated_region|intronic/){
      next();
    }
    $d++;
    my @line = split("\t", $_);
    my $coord = "$line[0]:$line[1]-$line[2]";
    my $symbol = $line[6];
    my $aa_change = $line[15];
    my $indel_type = $line[13];

    $data{$coord}{order} = $d;
    if ($data{$coord}{symbols}){
      my $symbol_ref = $data{$coord}{symbols};
      $symbol_ref->{$symbol}=1;
      my $aa_ref = $data{$coord}{aa_changes};
      $aa_ref->{$aa_change}=1;
      my $it_ref = $data{$coord}{indel_types};
      $it_ref->{$indel_type}=1;
    }else{
      my %sym;
      my %aa;
      my %it;
      $sym{$symbol}=1;
      $aa{$aa_change}=1;
      $it{$indel_type}=1;
      $data{$coord}{symbols} = \%sym;
      $data{$coord}{aa_changes} = \%aa;
      $data{$coord}{indel_types} = \%it;
    }
  }
  close(IN);
  foreach my $pos (sort {$data{$a}->{order} <=> $data{$b}->{order}} keys %data){
    my %symbols = %{$data{$pos}{symbols}};
    my @symbols = keys %symbols;
    my $symbol_string = join(", ", @symbols);
    my %aa_changes = %{$data{$pos}{aa_changes}};
    my @aa_changes = keys %aa_changes;
    my $aa_changes_string = join(", ", @aa_changes);
    my %it_types = %{$data{$pos}{indel_types}};
    my @it_types = keys %it_types;
    my $it_types_string = join(", ", @it_types);
    print "$pos\t$symbol_string\t$aa_changes_string\t$it_types_string\n";
  }
}



