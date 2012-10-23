#####################################################################################
# Copyright (C) 
#    2007 Ken Chen (kchen@watson.wustl.edu)
# All rights reserved!
#####################################################################################
#####################################################################################
# MG::IO::Polyscan::Contig
#####################################################################################

use MG::IO::Polyscan::SNP;
use MG::IO::Polyscan::INDEL;

package MG::IO::Polyscan::Contig;
#------------------------------------------------
our $VERSION = '1.0';
#------------------------------------------------
 use strict;
 use warnings;
 use Log::Log4perl;


=head1 NAME

 MG::IO::Polyscan -- read and write information from/to a Polyscan file
 MG::IO::Polyscan::SNP
 MG::IO::Polyscan::INDEL
 MG::IO::Polyscan::Contig

=head1 SYNOPSIS

 my $pp=new MG::IO::Polyscan();
 $pp->readIn($inPolyscanfile);  #read in original polyscan.out

or

 my $pp=new MG::IO::Polyscan(polyscan => $inPolyscanfile);

=head2 Synopsis subheading

 MG::IO::Polyscan::SNP
 MG::IO::Polyscan::INDEL
 MG::IO::Polyscan::Contig

=head1 DESCRIPTION

This module reads and writes from a Polyscan file.

=head1 Constructor and Initialization

=item (object) new (polyscan => Polypred_filename)

=head2 arguments (optional)

=item yamlfile

 (Default: $0.yaml)     The YAML file takes the same arguments as given here
                        (except the 'yamlfile' argument).  See the synopsis
                        for examples.

=item polyscan

                        Polyscan (3.0) filename

=head2 Methods

=item examethod(argument)

Description of example method.

=cut

sub new {
  my ($class)=@_;
  my $self={
	    _logger=>undef,
	    _NAME=>undef,
	    _SNP_READSITES=>undef,
	    _INDEL_READSITES=>undef
	   };
  bless $self, $class;
  return $self;
}

sub getNAME{
  my ($self)=@_;
  return $self->{_NAME};
}

sub getSNP_READSITES{
  my ($self)=@_;
  return $self->{_SNP_READSITES};
}

sub getINDEL_READSITES{
  my ($self)=@_;
  return $self->{_INDEL_READSITES};
}

sub ReadIn{
  #parse in PolyScan SNP output
  my ($self, $fin)=@_;
  $_=<$fin>; chomp;
  $self->{_NAME}=$_;

  until($_=~/END\_CONTIG/ || eof($fin)) {
    do {
      $_=<$fin>;
    } until ($_=~/BEGIN/ ||eof($fin));

    if($_=~/\_INDEL/){
      my %indels;
      do {
	$_=<$fin>; chomp;
	if(/^(\d+)\s+(\d+)\s+(\S+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\S+)/){
	  my $cpos=$1;
	  my $read=$3;
	  my $psindel=new MG::IO::Polyscan::INDEL();
	  $psindel->setValue($_);
	  $indels{$cpos}{$read}=$psindel;
	}
      } until ($_=~/END\_INDEL/||eof($fin));
      $self->{_INDEL_READSITES}=\%indels;
    }
    elsif($_=~/\_SNP/){
      my %snps;
      do {
	$_=<$fin>; chomp;
	if(/^(\d+)\s+(\d+)\s+(\S+)\s+(\D+)/){
	  my $cpos1=$1;
	  my $read=$3;
	  my $psgt=new MG::IO::Polyscan::SNP();
	  $psgt->setValue($_);
	  $snps{$cpos1}{$read}=$psgt;
	}
      } until ($_=~/END\_SNP/||eof($fin));
      $self->{_SNP_READSITES}=\%snps;
    }
    $_=<$fin>;
  }
}

sub setIndels{
  my ($self, $indels)=@_;
  $self->{_INDEL_READSITES}=$indels;
}

sub getIndels{
  my ($self)=@_;
  return $self->{_INDEL_READSITES};
}

sub printOut{
  #print out PolyScan SNP output
  my ($self, $fout)=@_;

  printf $fout "<BEGIN_CONTIG>\n%s\n\n",$self->{_NAME};

  print $fout "<BEGIN_INDEL>\n";
  my %indels=%{$self->{_INDEL_READSITES}};
  foreach my $cpos(sort {$a <=> $b} keys %indels){
    foreach my $rd(sort keys %{$indels{$cpos}}){
      $indels{$cpos}{$rd}->printOut($fout);
    }
  }
  print $fout "<END_INDEL>\n\n";

  print $fout "<BEGIN_SNP>\n";
  my %snps=%{$self->{_SNP_READSITES}};
  foreach my $cpos(sort {$a <=> $b} keys %snps){
    foreach my $rd(sort keys %{$snps{$cpos}}){
      $snps{$cpos}{$rd}->printOut($fout);
    }
  }
  print $fout "<END_SNP>\n";


  print $fout "<END_CONTIG>\n";
}



1;
