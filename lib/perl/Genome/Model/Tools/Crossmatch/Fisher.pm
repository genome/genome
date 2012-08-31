# Copyright (C) 2008 Washington University in St. Louis
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

########################################################################################
# Author: Ken Chen (kchen@genome.wustl.edu)
# Date: July, 2008
# This module can used to detect SNP and indels based on the cross-match alignment
# of 454 reads to the reference sequences using Fisher Exact Test
########################################################################################

use strict;
use warnings;
package Genome::Model::Tools::Crossmatch::Fisher;

sub new{
  my ($class, %arg) = @_;
  my $self={
	    fin=>$arg{fin} || undef,
	    dcpos=>undef,
	    align=>undef,
	    readpos=>undef,
	    min_base_qual=>$arg{MinBaseQual} || 15,
	    homopolymer_size=>$arg{HomoPolymerIndelSize} || 2
	   };
  bless($self, $class || ref($class));
  ($self->{align},$self->{dcpos})=$self->Read($self->{fin});
  return $self;
}

sub Read{
  my ($self,$f_cm)=@_;
  my %DCReads;
  my %DCPoses;
  my @discrepancies;
  my ($paln_readname,$tmpstr,$score,$perc_sub,$perc_del,$perc_ins,$readname,$r1,$r2,$r3,$comp,$refseq,$g1,$g2,$g3);
  open(CM,"<$f_cm") || die "unable to open $f_cm\n";
  while(<CM>){
    chomp;
    if (/^ALIGNMENT/){
      if($#discrepancies>=0){
	if(defined $paln_readname){
	  foreach my $disc(@discrepancies){
	    push @{$DCReads{$paln_readname}->{discrepancy}},$disc;
	  }
	  undef @discrepancies;
	}
      }
      my @u=split /\s+/;
      pop @u if($u[$#u] eq '*');  #include both higher and lower scoring match

      if($#u==12){
	($tmpstr,$score,$perc_sub,$perc_del,$perc_ins,$readname,$r1,$r2,$r3,$refseq,$g1,$g2,$g3)=@u;
	undef $comp;
      }
      elsif($#u>=13){
	($tmpstr,$score,$perc_sub,$perc_del,$perc_ins,$readname,$r1,$r2,$r3,$comp,$refseq,$g1,$g2,$g3)=@u;
      }
      else{
      }
      $paln_readname=$readname;

      my ($r_start,$r_end,$r_rest);
      if($r1=~/\(/){
	$r_start=$r2;$r_end=$r3;
	($r_rest)=($r1=~/\((\d+)\)/);
      }
      else{
	$r_start=$r1;$r_end=$r2;
	($r_rest)=($r3=~/\((\d+)\)/);
      }

      my $rlen=$r_end+$r_rest;
      my ($ref_start,$ref_end);
      if($g1=~/\(/){
	$ref_start=$g2; $ref_end=$g3;
      }
      else{
	$ref_start=$g1; $ref_end=$g2;
      }

      my $aln;
      ($aln->{score},$aln->{perc_sub},$aln->{perc_del},$aln->{perc_ins},$aln->{r_start},$aln->{r_end},$aln->{refseq},$aln->{ref_start},$aln->{ref_end},$aln->{orientation})=($score,$perc_sub,$perc_del,$perc_ins,$r_start,$r_end,$refseq,$ref_start,$ref_end, $comp ||'U');
      push @{$DCReads{$readname}->{aln}},$aln;
    }
    elsif(/DISCREPANCY/){
      my @u=split /\s+/;
      if($u[1]=~/([SDI])\-*(\d*)/){  #substitutions and indels
	my ($dctype,$dcsize)=($1,$2);
	$dcsize=1 if(!defined $dcsize || length($dcsize)<=0);
	shift @u;
	push @discrepancies,join("\t",@u);
	
	my $dcpos=$u[3];
	my ($base,$qual)=($u[2]=~/(\S)\((\d+)\)/);
	if(defined $comp && $comp eq 'C'){  #reverse complement
	  $base=~tr/ACGT/TGCA/;
	  #complement indel position
	  if($dctype eq 'D'){
	    $dcpos=$dcpos-$dcsize+1;
	  }
	  elsif($dctype eq 'I'){
	    $dcpos=$dcpos-1;
	  }
	  else{}
	}

	my $dcinfo;
	($dcinfo->{type},
	 $dcinfo->{rpos},
	 $dcinfo->{base},
	 $dcinfo->{qual},
	 $dcinfo->{read},
	 $dcinfo->{aln_orient},
	 $dcinfo->{score}
	)=($u[0],$u[1],$base,$qual,$readname,$comp || 'U',$score);
	
	my $maxlen_homopolymer=&Indel_454_homopolymerErr($u[$#u],$comp);
	if( $qual>=$self->{min_base_qual} && #Require minimal quality score to be trusted
	    (
	     ($dctype eq 'S') ||
	     ($dctype=~/[DI]/ && $dcsize>1) ||
	     ($dctype=~/[DI]/ && $dcsize==1 && $maxlen_homopolymer<$self->{homopolymer_size})
	    )
	  )
	  {
	    my @DClist;
	    if(defined $DCPoses{$dcpos}){
	      @DClist=@{$DCPoses{$dcpos}};
	    }
	    push @DClist,$dcinfo;
	    $DCPoses{$dcpos}=\@DClist;
	  }
      }
    }
    elsif(/SCORE_HISTOGRAM/){
      my @u=split /\s+/;
      my ($bestscore,$nbest)=($u[$#u]=~/(\d+)\((\d+)\)/);
      my ($sec_bestscore,$nsec_best)=($u[$#u-1]=~/(\d+)\((\d+)\)/);
      ($DCReads{$u[1]}->{best},$DCReads{$u[1]}->{nbest},$DCReads{$u[1]}->{sec_best},$DCReads{$u[1]}->{nsec_best},$DCReads{$u[1]}->{associate})=($bestscore,$nbest,$sec_bestscore||0,$nsec_best||0,$readname);
      }
    else{}
  }
  if($#discrepancies>=0){
    $DCReads{$paln_readname}->{discrepancy}=\@discrepancies if(defined $paln_readname);
  }
  close(CM);
  return (\%DCReads,\%DCPoses);
}

sub Indel_454_homopolymerErr{
  #Test find if a 1bp indel is caused by upstream homopolymer errors
  my ($string,$comp)=@_;
  my @bases=split //,$string;
  my $midpos=int($#bases/2);
  my $maxlen_homopolymer=0;

  foreach my $start($midpos,$midpos+1){
    my $pbase=$bases[$start];
    my $end;
    for($end=$start+1; $end<=$#bases; $end++){
      my $base=$bases[$end];
      last if($base ne $pbase);
    }
    my $len=$end-$start;
    $maxlen_homopolymer=($maxlen_homopolymer<$len)?$len:$maxlen_homopolymer;
  }

  foreach my $start($midpos,$midpos-1){
    my $pbase=$bases[$start];
    my $end;
    for($end=$start-1; $end>=0; $end--){
      my $base=$bases[$end];
      last if($base ne $pbase);
    }
    my $len=$start-$end;
    $maxlen_homopolymer=($maxlen_homopolymer<$len)?$len:$maxlen_homopolymer;
  }

  return $maxlen_homopolymer;
}

sub GetAlleleInfo{
  my ($self,$refseq,$f_qual,$refpos,$Min_Reads)=@_;
  my $refbase=substr $refseq, $refpos-1, 1;
  my %Allele_info;

  my $DCreadsCount=0;
  if(defined $self->{dcpos}{$refpos}){
    my @DCreads=@{$self->{dcpos}{$refpos}};
    $DCreadsCount=$#DCreads+1;
    if($DCreadsCount>=$Min_Reads){
      foreach my $read(@DCreads){
	my $base=$read->{base} || $refbase;
	my $DCtype=$read->{type};
	$base='-' if($DCtype=~/[DI]/);
	$Allele_info{$DCtype}{$base}{$read->{aln_orient}}{SumQual}+=$read->{qual};
	if(!defined $Allele_info{$DCtype}{$base}{$read->{aln_orient}}{MaxQual} ||
	   $Allele_info{$DCtype}{$base}{$read->{aln_orient}}{MaxQual}<$read->{qual}){
	  $Allele_info{$DCtype}{$base}{$read->{aln_orient}}{MaxQual}=$read->{qual};
	}
	$Allele_info{$DCtype}{$base}{$read->{aln_orient}}{Rpos}{$read->{rpos}}++;
      }
    }
  }

  my %readpos;
  #Map reference positions to unpadded read positions
  foreach my $read(keys %{$self->{align}}){  #all the reads that are aligned
    my $align=${$self->{align}}{$read};
    if(($align->{ref_start}<=$refpos) && ($refpos<=$align->{ref_end})){
      my $rpos;
      if($align->{orientation} eq 'U'){  #forward alignment
	$rpos=$refpos-$align->{ref_start}+$align->{r_start};
      }
      else{  #reverse alignment
	$rpos=$align->{ref_end}-$refpos+$align->{r_start};
      }
      foreach my $dc(@{$align->{discrepancy}}){
	my ($dctype,$pos)=($dc=~/(\S+)\s+(\d+)/);
	next if($pos>$rpos);  #doesn't affect alignment coordinates
	if($dctype=~/D\-*(\d*)/){
	  $rpos-=$1 || 1;
	}
	elsif($dctype=~/I\-*(\d*)/){
	  $rpos+=$1 || 1;
	}
	else{}
      }
      $readpos{$read}=$rpos;
    }
  }

  #Get Quality from file
  open(FQIN,"<$f_qual") || die "$f_qual is not available\n";
  do{
    $_=<FQIN>;
  } until (/^>(\S+)\s/ || eof(FQIN));
  my $readname=$1;
  while(!eof(FQIN)){
    my $qualstr='';
    do{
      $_=<FQIN>; chomp;
      $qualstr=join(' ',$qualstr,$_) if($_ !~ /^>(\S+)\s/);
    } until (/^>(\S+)\s/ || eof(FQIN));
    if(!eof(FQIN)){
      my @quals=split /\s+/,$qualstr;
      if(defined $readpos{$readname}){
	my $align=${$self->{align}}{$readname};
	my $rpos=$readpos{$readname};
	warn "rpos=$rpos out of scope $#quals in getAlleleInfo in CrossMatch.pm" if($rpos<1 || $rpos>$#quals);
	my $qual=$quals[$rpos];
	if($qual>$self->{min_base_qual}){
	  my $base=$refbase;
	  my $DCtype='W';
	  if(!defined $base){
	    print "";
	  }

	  $Allele_info{$DCtype}{$base}{$align->{orientation}}{SumQual}+=$qual;
	  if(!defined $Allele_info{$DCtype}{$base}{$align->{orientation}}{MaxQual} ||
	     $Allele_info{$DCtype}{$base}{$align->{orientation}}{MaxQual}<$qual){
	    $Allele_info{$DCtype}{$base}{$align->{orientation}}{MaxQual}=$qual;
	  }
	  $Allele_info{$DCtype}{$base}{$align->{orientation}}{Rpos}{$rpos}++;
	}
      }
      ($readname)=($_=~/^\>(\S+)\s/);
      #print "$readname\n";
    }
  }
  close(FQIN);
  return (\%Allele_info);
}

1;
