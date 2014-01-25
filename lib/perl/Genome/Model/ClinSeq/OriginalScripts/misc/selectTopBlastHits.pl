#!/usr/bin/env genome-perl -w
#Written by Malachi Griffith

#Example usage:
#BLAST COMMAND or cat of blast file | selectTopHits.pl

use strict;
use Getopt::Long;
use Term::ANSIColor qw(:constants);

my $top_hits = '';

GetOptions ('top_hits=i'=>\$top_hits);

unless ($top_hits =~ /^\d+/){
  print RED, "\nRequired input parameter(s) missing or incorrect format!\n\n", RESET;
  print RED, "\nExample usage: BLAST COMMAND or cat of blast file | selectTopHits.pl  --top_hits=1 [number to display]\n", RESET;
  exit();
}

my %hits;
my $c = 1;
my $hc = 0;
my $last_query;
while(<STDIN>){
  chomp($_);
  my $line = $_;
  my @line = split("\t", $_);
  my $query_id = $line[0];	
  my $bitscore = $line[11];
  $hc++;

  if ($c == 1){
    $last_query = $query_id;
  }


  if ($query_id eq $last_query){
    #Store the hit in the growing hit list
    $hits{$hc}{line} = $line;
    $hits{$hc}{bitscore} = $bitscore;

    if (eof){
      &printHits('-hits'=>\%hits);
    }

  }else{
    &printHits('-hits'=>\%hits);
    
    #Wipe the hit list and store the first hit of a new list
    %hits = ();
    $hc = 1;
    $hits{$hc}{line} = $line;
    $hits{$hc}{bitscore} = $bitscore;
 
    if (eof){
      &printHits('-hits'=>\%hits);
    }

  }
  $c++;
  $last_query = $query_id;
}

exit();


sub printHits{
  my %args = @_;
  my %hits = %{$args{'-hits'}};
  my $oc = 0;
  foreach my $h (sort {$hits{$b}->{bitscore} <=> $hits{$a}->{bitscore}} keys %hits){
    $oc++;
    if ($oc <= $top_hits){
      print "$hits{$h}{line}\n"
    }
  }
  return();
}




