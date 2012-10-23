#####################################################################################
# Copyright (C) 
#    2007 Ken Chen (kchen@watson.wustl.edu)
# All rights reserved!
#####################################################################################

#####################################################################################
# MG::IO::Polyscan::SNP
#####################################################################################

package MG::IO::Polyscan::SNP;
#------------------------------------------------
our $VERSION = '1.0';
#------------------------------------------------
 use strict;
 use warnings;
 use Log::Log4perl;

#23      120     H_P013-PCR018_0013.g1   CT      45      heterozygous

sub new {
  my ($class)=@_;
  my $self={
	    _refpos=>undef,
	    _rdpos=>undef,
	    _read=>undef,
	    _allel1=>undef,
	    _allel2=>undef,
	    _score=>undef,
	    _var_tag=>undef,
	    _logger=>undef
	   };
  bless $self, $class;
  return $self;
}

sub setValue{
  my ($self, $line)=@_;
  my @u=split /\s+/,$line;
  if($#u<3 || $#u>5){
    my($log) = $self->{_logger};
    $log->error("SNP annotation in polyscan file needs to have 4 to 6 columns");
    return;
  }
  else{
    $self->{_refpos}=$u[0];
    $self->{_rdpos}=$u[1];
    $self->{_read}=$u[2];
    my $gt=$u[3];
    my ($allel1,$allel2)=split //, uc($gt);
    $self->{_allel1}=$allel1;
    $self->{_allel2}=$allel2;
    $self->{_score}=$u[4];
    if($#u>4){
      $self->{_var_tag}=$u[5];
    }
  }
}

sub get{
  my ($self)=@_;
  return ($self->{_refpos},$self->{_rdpos}, $self->{_read}, $self->{_allel1}, $self->{_allel2},$self->{_score}, $self->{_var_tag});
}

sub printOut{
  my ($self, $fout) = @_;
  printf $fout "%s\t%s\t%s\t%s\t%s",$self->{_refpos},$self->{_rdpos},$self->{_read},$self->{_allel1}.$self->{_allel2},$self->{_score};
  if(defined $self->{_score}){
    printf $fout "\t%s", $self->{_score};
    if(defined $self->{_var_tag}){
      printf $fout "\t%s\n", $self->{_var_tag};
    }
    else{
      printf $fout "\n";
    }
  }
  else{
    printf $fout "\n";
  }

}
