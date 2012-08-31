#!/usr/bin/env genome-perl
# Merge call sets produced by AssemblyValidation.pl scripts
use strict;
use warnings;
use Getopt::Std;
use Bio::SeqIO;
use Bio::Seq;

my %opts = (d=>20,p=>20,s=>10,z=>0.1,w=>50, l=>200, r=>0.8, C=>0.5, n=>1000);
getopts('p:d:f:w:s:z:hL:ce:l:r:C:n:',\%opts);
die("
Usage:   MergeAssembledCallsets.pl <the result index file>\n
         the tab/space delimited index file must have 3 columns: set name, .csv file, .fasta file
Options:
         -d INT     do not merge SVs unless they differ less than [$opts{d}] bp in start and end positions, subject to 1 duplication
         -p INT     do not merge SVs unless they differ less than [$opts{p}] bp in size
         -z FLOAT   ignore SVs with more than $opts{z} size difference from the expected size
         -n INT     ignore SVs with assembled breakpoint [$opts{n}] bp off the predicted breakpoint
         -s INT     bp minimal size to trust [$opts{s}]
         -f FILE    dump breakpoint sequences to a fasta file
         -L STRING  ignore SVs that are detected in sets with STRING in the name
         -w INT     filter out prediction with weighted assembly score smaller than [$opts{w}]
         -h         apply a high quality alignment filter on the calls
         -c         sort result by assembly score
         -C         always include SVs with CNA > $opts{C}
         -e FILE    exclude assembled SV overlap with Breakdancer SVs in FILE
         -l INT     use in conjunction with -e options: ambiguity allowed at breakpoints [$opts{l}] bp
         -r FLOAT   use in conjunction with -e options: fraction of overlap bewteen 2 events [$opts{r}]
") unless (@ARGV);

my %f_csv;
my %f_fasta;
my %diff_size;
my %diff_loc;
open(INX,"<$ARGV[0]") || die "unable to open $ARGV[0]\n";
while(<INX>){
  chomp;
  next if(/^\#/);
  my ($set,$f_indel,$f_fasta,$diff_size,$diff_loc)=split;
  push @{$f_csv{$set}},$f_indel;
  #push @{$f_fasta{$set}},$f_fasta;
  $f_fasta{$set}=$f_fasta;
  $diff_size{$set}=$diff_size if(defined $diff_size);
  $diff_loc{$set}=$diff_loc if(defined $diff_loc);
}

my %SVs;
foreach my $set(keys %f_csv){
  my @files=@{$f_csv{$set}};
  foreach my $file(@files){
    &AddSVs($set,$file);
  }
}

my %EXSVs;
if($opts{e}){  #exclude SVs
  open(SVIN,"<$opts{e}") || die "unable to open $opts{e}\n";
  while(<SVIN>){
    next if(/^\#/);
    chomp;
    my $sv;
    my @extra;
    ($sv->{chr1},$sv->{start},$sv->{ori1},$sv->{chr2},$sv->{end},$sv->{ori2},$sv->{type},$sv->{size},$sv->{qual},$sv->{nrps},$sv->{rp_sample},$sv->{allele_frequency},@extra)=split /\s+/;
    next unless(defined $sv->{type} && defined $sv->{chr1} && defined $sv->{chr2});
    push @{$EXSVs{$sv->{type}}{$sv->{chr1}}},$sv;
  }
  close(SVIN);
  foreach my $chr(sort keys %SVs){
    my $idx=1;
    foreach my $chr2(sort keys %{$SVs{$chr}}){
      foreach my $pos(sort {$a<=>$b} keys %{$SVs{$chr}{$chr2}}){
	my $SV=$SVs{$chr}{$chr2}{$pos};
	my $hit=&Overlap($chr,$chr2,$pos,$SV);
	if($hit){
#	  printf "exclude: %s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\n",$chr,$SV->{OuterStart},$SV->{InnerStart},$chr2,$SV->{InnerEnd},$SV->{OuterEnd},$SV->{type},join(',',@{$SV->{ori}}),$SV->{minsize},$SV->{maxsize},join(',',@{$SV->{group}}),join(',',@{$SV->{wAsmscore}});
	  delete $SVs{$chr}{$chr2}{$pos};
	}
      }
    }
  }
}

#print table
#printf "#ID\tCHR\tPOS\tTYPE\tSIZE\tGROUPS\tHET\n";
#printf "#ID\tCHR1\tOUTER_START\tINNER_START\tCHR2\tINNER_END\tOUTER_END\tTYPE\tORIENTATION\tMINSIZE\tMAXSIZE\tSOURCE\tSCORES\tCopy_Number\tGene\tKnown\n";
printf "#ID\tCHR1\tOUTER_START\tINNER_START\tCHR2\tINNER_END\tOUTER_END\tTYPE\tORIENTATION\tMINSIZE\tMAXSIZE\tSOURCE\tSCORES\tCopy_Number\n";

my $fout;
if($opts{f}){
  if( -s $opts{f}){
    `rm -f $opts{f}`;
  }
  $fout = Bio::SeqIO->new(-file => ">>$opts{f}" , '-format' => 'Fasta');
}

my %SVscores;
foreach my $chr(sort keys %SVs){
  my $idx=1;
  foreach my $chr2(sort keys %{$SVs{$chr}}){
    foreach my $pos(sort {$a<=>$b} keys %{$SVs{$chr}{$chr2}}){
      my $SV=$SVs{$chr}{$chr2}{$pos};
      my $cid="$chr.$idx";
      #my $pos=($SV->{OuterStart}==$SV->{InnerStart})?$SV->{OuterStart}:join(',',$SV->{OuterStart},$SV->{InnerStart});
      my $size=($SV->{minsize}==$SV->{maxsize})?$SV->{minsize}:join('-',$SV->{minsize},$SV->{maxsize});
      #printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$cid,$chr,$pos,"Deletion",$size,join(',',@{$SV->{group}}),join(',',@{$SV->{het}});
      my $maxscore=-1;
      my $winningGroup;
      my $winningContig;
      my $ignore=0;
      for(my $i=0;$i<=$#{$SV->{wAsmscore}};$i++){
	if($maxscore<${$SV->{wAsmscore}}[$i]){
	  $maxscore=${$SV->{wAsmscore}}[$i];
	  $winningGroup=${$SV->{group}}[$i];
	  $winningContig=${$SV->{prestr}}[$i];
	}
	$ignore=1 if(defined $opts{L} && ${$SV->{group}}[$i]=~/$opts{L}/i);
      }
      next if($ignore);

      my $keys=join('|',$cid,$chr,$chr2,$pos,$winningGroup,$winningContig);
      $SVscores{$keys}=$maxscore;
      $idx++;
    }
  }
}

my @SVkeys;
if($opts{c}){
  @SVkeys=sort {$SVscores{$b}<=>$SVscores{$a}} keys %SVscores;
}
else{
  @SVkeys=sort keys %SVscores;
}

foreach my $key(@SVkeys){
  my ($cid,$chr,$chr2,$pos,$winningGroup,$winningContig)=split /\|/,$key;
  my $SV=$SVs{$chr}{$chr2}{$pos};
  printf "%s\t%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\t%s\n",$cid,$chr,$SV->{OuterStart},$SV->{InnerStart},$chr2,$SV->{InnerEnd},$SV->{OuterEnd},$SV->{type},join(',',@{$SV->{ori}}),$SV->{minsize},$SV->{maxsize},join(',',@{$SV->{group}}),join(',',@{$SV->{wAsmscore}}),$SV->{extra};

  if($opts{f}){
    my $seqobj;
    my $contigfile=$f_fasta{$winningGroup};
    my $in  = Bio::SeqIO->newFh(-file => "$contigfile" , '-format' => 'Fasta');
    while ( my $seq = <$in> ) {
      # do something with $seq
      my ($id,$var,$ins,$strand,$score)=split /\,/,$seq->id;
      $id=~s/ID\://;
      next unless($id eq $winningContig);
      my @u=split /\,/,$seq->id;
      shift @u;
      my $sequence=$seq->seq();
      my $dispid=join(',',$cid,$winningGroup,@u);
      $seqobj = Bio::Seq->new( -display_id => "$dispid",
			       -seq => $sequence);
      last;
    }
    if(defined $seqobj){
      $fout->write_seq($seqobj);
    }
    else{
      print STDERR "Contig $winningContig not found in $contigfile\n";
    }
  }

}

sub AddSVs{
  my ($set,$fin)=@_;
  open(FIN,"<$fin") || die "unable to open $fin\n";
  while(<FIN>){
    chomp;
    my ($chr,$start,$chr2,$end,$ori,$size,$type,$het,$wAsmscore,$read_len,$perc_aligned,$n_seg,$n_sub,$n_indel,$nbp_indel,$microhomology,$scarstr,$prestr,$asm_parm,@extra)=split /\s+/;
    next unless(defined $prestr && $size=~/^\d/);
    my ($pchr1,$pre_pos1,$pchr2,$pre_pos2,$pre_type,$pre_size,$pre_ori)=split /\./,$prestr;
    my $size_diff_cutoff=$diff_size{$set} || $opts{z};
    $size_diff_cutoff=1e10 if($type && $type=~/ctx/i);
    my $loc_diff_cutoff=$diff_loc{$set} || $opts{n};
    $start=~s/\(\d+\)//g;
    $end=~s/\(\d+\)//g;
    $size=~s/\(\d+\)//g;
    $type=~s/\(\w+\)//g;
    my $extra_info=join("\t",@extra);
    my $f_cna=0;

    if(@extra && $extra[0]=~/\.bam/){
      my @cnstr=split /\,/,$extra[0];
      my ($ncn)=($cnstr[0]=~/\:(\S+)$/);
      if($ncn!~/NA/i){
	for(my $i=1;$i<=$#cnstr;$i++){
	  my ($cn)=($cnstr[$i]=~/\:(\S+)$/);
	  $f_cna=$cn-$ncn if($cn!~/NA/i);
	}
      }
    }

    #filtering criteria

    $perc_aligned=~s/\%//;
    my $flanksize=$read_len*$perc_aligned/100;
    my $sub_rate=$n_sub/$flanksize;
    my $subindel_rate=$n_indel/$flanksize;

    if(defined $opts{h} &&  ($n_seg>2
			     || $n_seg==2 && (($sub_rate>0.006 && $n_sub>2) || ($subindel_rate> 0.002 && $n_indel>1)|| $nbp_indel>20)
			     || $n_seg==1 && ($sub_rate>0.006 || $subindel_rate>0.001 || $nbp_indel>5)
			     #|| $n_seg==2 && (($sub_rate>0.006 && $n_sub>2) || ($subindel_rate> 0.002 && $n_indel>1)|| $nbp_indel>20) && ($size<=99999) && (abs($f_cna)<$opts{C})
			     #|| $n_seg==1 && ($sub_rate>0.006 || $subindel_rate>0.001 || $nbp_indel>5)&& ($size<=99999) && (abs($f_cna)<$opts{C})
			     || ($type ne 'INS' ) && ($perc_aligned<80)
			     #|| ($type eq 'DEL' ) && ($size>99999) && ($f_cna>-$opts{C})
			     #|| ($type eq 'ITX' ) && ($size>99999) && ($f_cna<$opts{C})
			     #|| ($type eq 'INV' ) && ($size>99999)
			     #|| $microhomology>=200
			     || $size<$opts{s} || $type ne $pre_type || ($pre_size>0 && abs($size-$pre_size)/$pre_size>$opts{z}) || abs($start-$pre_pos1)>$loc_diff_cutoff
			     || $wAsmscore >0 && $wAsmscore<$opts{w}
			    )
      ){
      printf STDERR  "%s\t%s\t%d\t%s\t%d\t%s\t%d\t%s\t%s\t%d\t%d\t%d%%\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\n",$set,$chr,$start,$chr2,$end,$ori,$size,$type,$het,$wAsmscore,$read_len,$perc_aligned,$n_seg,$n_sub,$n_indel,$nbp_indel,$microhomology,$scarstr,$prestr,$asm_parm,$extra_info;
      next;
    }

    #printf "%s\t%d\t%s\t%d\t%s\t%d\t%s\t%s\t%d\t%d\t%d%%\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\n", $chr,$start,$chr2,$end,$ori,$size,$type,$het,$wAsmscore,$read_len,$perc_aligned,$n_seg,$n_sub,$n_indel,$nbp_indel,$microhomology,$scarstr,$prestr,$asm_parm;

    my $found;
    my @tds=($type && $type=~/ctx/i)?(0):(-$size,0,$size);
    foreach my $td(@tds){
      my $locrange=int($size*0.01);
      my $sizerange=$locrange;
      #$locrange=($locrange>$opts{d})?$locrange:$opts{d};
      #$sizerange=($sizerange>$opts{p})?$sizerange:$opts{p};
      $locrange=$opts{d};
      $sizerange=$opts{p};
      for (my $i=-$locrange+$td;$i<=$locrange+$td;$i++){
	my $SV=$SVs{$chr}{$chr2}{$start+$i};
	if(defined $SV
	   && ($SV->{maxsize}<=$size+$sizerange && $SV->{minsize}>=$size-$sizerange || $type=~/ctx/i)
	   && $end>=$SV->{OuterEnd}-$locrange && $end<=$SV->{InnerEnd}+$locrange
	  ){
	  $found=$i;
	  last;
	}
      }
    }

    my $SVregion;

    if(defined $found){
      $SVregion=$SVs{$chr}{$chr2}{$start+$found};
      $SVregion->{minsize}=($SVregion->{minsize}>$size)?$size:$SVregion->{minsize};
      $SVregion->{maxsize}=($SVregion->{maxsize}<$size)?$size:$SVregion->{maxsize};
      $SVregion->{OuterStart}=($SVregion->{OuterStart}>$start)?$start:$SVregion->{OuterStart};
      $SVregion->{InnerStart}=($SVregion->{InnerStart}<$start)?$start:$SVregion->{InnerStart};
      $SVregion->{OuterEnd}=($SVregion->{OuterEnd}<$end)?$end:$SVregion->{OuterEnd};
      $SVregion->{InnerEnd}=($SVregion->{InnerEnd}>$end)?$end:$SVregion->{InnerEnd};
    }
    else{
      $SVregion->{minsize}=$size;
      $SVregion->{maxsize}=$size;
      $SVregion->{OuterStart}=$start;
      $SVregion->{InnerStart}=$start;
      $SVregion->{OuterEnd}=$end;
      $SVregion->{InnerEnd}=$end;
      $SVregion->{type}=$type;
      $SVregion->{extra}=$extra_info;
    }

    my $exist=0;
    for(my $ki=0;$ki<=$#{$SVregion->{group}};$ki++){
      my $group=${$SVregion->{group}}[$ki];
      if($group eq $set){
	$exist=1;
	if($wAsmscore>${$SVregion->{wAsmscore}}[$ki]){  #select the one with higher assembly score
	  ${$SVregion->{wAsmscore}}[$ki]=$wAsmscore;
	  ${$SVregion->{prestr}}[$ki]=$prestr;
	  ${$SVregion->{het}}[$ki]=$het;
	  ${$SVregion->{ori}}[$ki]=$ori;
	}
	last;
      }
    }
    if(! $exist){
      push @{$SVregion->{group}},$set;
      push @{$SVregion->{het}},$het;
      push @{$SVregion->{wAsmscore}},$wAsmscore;
      push @{$SVregion->{prestr}},$prestr;
      push @{$SVregion->{ori}},$ori;
      $SVs{$chr}{$chr2}{$start+($found||0)}=$SVregion;
    }
  }
  close(FIN);
}


sub Overlap{
  my ($chr1,$chr2,$pos,$SV)=@_;
  my ($pt,$pc1,$ps,$pc2,$pe)=($SV->{type},$chr1,($SV->{OuterStart}+$SV->{InnerStart})/2,$chr2,($SV->{InnerEnd}+$SV->{OuterEnd})/2);
  my $hit=0;
  if(defined $EXSVs{$SV->{type}}{$chr1}){
    my @bsvs=@{$EXSVs{$SV->{type}}{$chr1}};
    foreach my $bsv(@bsvs){
      my ($t,$c1,$s,$c2,$e)=($bsv->{type},$bsv->{chr1},$bsv->{start},$bsv->{chr2},$bsv->{end});
      if(defined $pe && defined $e && $pt eq $t){ # same SV type
	my $breakpoints_overlap=($pc1 eq $c1 && $pc2 eq $c2 && abs($ps-$s)<=$opts{l} && abs($pe-$e)<=$opts{l})?1:0;
	my $interval_overlap=0;
	if($c1 eq $c2){  #same chromosome
	  ($ps,$pe)=sort {$a<=>$b} ($ps,$pe);
	  ($s,$e)=sort {$a<=>$b} ($s,$e);
	  my $length=$pe-$ps+1;
	  my $size=$e-$s+1;
	  my $overlap=(($e<$pe)?$e:$pe)-(($s>$ps)?$s:$ps)+1;
	  my $perc1=$overlap/$size;
	  my $perc2=$overlap/$length;
	  $interval_overlap=(defined $opts{r} && $perc1>=$opts{r} && $perc2>=$opts{r})?1:0;
	}
	if($breakpoints_overlap || $interval_overlap){
	  $hit=1;
	  last;
	}
      }
    }
  }

  return $hit;
}
