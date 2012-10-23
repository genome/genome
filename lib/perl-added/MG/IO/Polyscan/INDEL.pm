#####################################################################################
# Copyright (C) 
#    2007 Ken Chen (kchen@watson.wustl.edu)
# All rights reserved!
#####################################################################################

#####################################################################################
# MG::IO::Polyscan::INDEL
#####################################################################################

package MG::IO::Polyscan::INDEL;
#------------------------------------------------
our $VERSION = '1.0';
#------------------------------------------------
 use strict;
 use warnings;
 use Log::Log4perl;

#6874    519     H_P015-PCR008_0121.g1   20      deletion        TTCATTTGCTCATTCATTCA    36      0.89

sub new {
  my ($class)=@_;
  my $self={
	    _refpos=>undef,
	    _rdpos=>undef,
	    _read=>undef,
	    _size=>undef,
	    _type=>undef,
	    _seqindel=>undef,
	    _score=>undef,
	    _spratio=>undef,
	    _comment=>undef,
	    _refpos2=>undef, #this helps define a region where indel pos is ambiguous
	    _logger=>undef
	   };
  bless $self, $class;
  return $self;
}

sub setValue{
  my ($self, $line)=@_;
  my($log) = $self->{_logger};
  my @u=split /\s+/,$line;
  if($#u<7){
    $log->error("INDEL annotation in polyscan.out file should have at least 8 columns");
    return;
  }
  else{
    $self->{_refpos}=$u[0];
    $self->{_rdpos}=$u[1];
    $self->{_read}=$u[2];
    $self->{_size}=$u[3];
    $self->{_type}=$u[4];
    $self->{_seqindel}=$u[5];
    $self->{_score}=$u[6];
    $self->{_spratio}=$u[7];
    $self->{_comment}=$u[8] if (defined $u[8]);
    $self->{_refpos2}=$u[9] if (defined $u[9]);
  }
}

sub get{
  my ($self)=@_;
  return ($self->{_refpos},$self->{_rdpos},$self->{_read},$self->{_size},$self->{_type},$self->{_seqindel},$self->{_score},$self->{_spratio},$self->{_comment},$self->{_refpos2});

}

sub getIndelSize{
  my ($self)=@_;
  return $self->{_size};
}

sub getIndelPos{
  my ($self)=@_;
  return $self->{_refpos};
}

sub getIndelType{
  my ($self)=@_;
  return $self->{_type};
}

sub printOut{
  my ($self, $fout)=@_;
  printf $fout "%d\t%d\t%s\t%d\t%s\t%s\t%d\t%.2f",$self->{_refpos},$self->{_rdpos},$self->{_read},$self->{_size},$self->{_type},$self->{_seqindel},$self->{_score},$self->{_spratio};
  if(defined $self->{_comment}){
    print $fout "\t%s\n", $self->{_comment};
  }
  else{
    print $fout "\n";
  }
}
