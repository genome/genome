#!/usr/bin/env genome-perl
# Annotate Assembled SVs based UCSC gene table
use strict;
use warnings;
use Getopt::Std;
use Bio::SeqIO;
use Bio::Seq;
use DBI;

my $dir = '/gsc/scripts/share/BreakAnnot_file';
my $build_dir = $dir . '/human_build36';
my $mousetable       = $dir . '/Mouse.July2007.RefSeqgene.tab';
my $flist_cancergene = $dir . '/Cancer_genes.csv';
my %opts = (
	    g=> $build_dir.'/Human.Mar2006.RefSeqGenes.tab',
	    s=> $build_dir.'/Human.Mar2006.SegDups.tab',
	    p=> $build_dir.'/dbsnp130.indel.named.csv',
	    v=> $build_dir.'/ncbi36_submitted.gff',
	    L=>50,
	    l=>200,
	    w=>200,
	    r=>0.5,
	    k=>5,
	    m=>100
	   );
getopts('g:l:Ma:b:c:d:o:r:AL:w:m:k:p:v:s:',\%opts);
die("
Add annotation (UCSC gene, dbSNP and others) to a text file
Usage:   BreakAnnot.pl <a BreakDancer output file>\n
Options:
         -g FILE  Use UCSC Human gene table at $opts{g}
         -M       Use UCSC Mouse gene table at $mousetable
         -p FILE  Use UCSC dbSNP file at $opts{p}
         -v FILE  Use dbVar file at $opts{v}
         -s FILE  Use UCSC segmental duplication file at $opts{s}
         -o STR   Only analysis the specified chromosome
         -l INT   Proximity of breakpoints to gene annotations [$opts{l}] bp
         -L INT   Proximity of breakpoints to segmental duplication [$opts{L}] bp
         -w INT   Look for repeat annotation within +- [$opts{w}] bp of breakpoint
         -m INT   annotate breakpoint +- [$opts{k}] bp overlapping more than [$opts{m}] bp of masked repeat
         -k INT   See -m [$opts{k}]
         -A       Annotate Merged SV assembly file
         -a INT   Which column contains the chromosome of the breakpoint 1 (1-based)
         -b INT   Which column contains the position of the  breakpoint1 (1-based)
         -c INT   Which column contains the chromosome of the breakpoint 2 (1-based)
         -d INT   Which column contains the position of the breakpoint 2 (1-based)
         -r FLOAT Fraction of overlap (reciprocal) required to hit a dbSNP/dbVar SVs [$opts{r}]
\n") unless (@ARGV);

$opts{g}=$mousetable if(defined $opts{M});  #switch to mouse
if($opts{A}){
  $opts{a}=2; $opts{b}=4; $opts{c}=5; $opts{d}=6;  #SV assembly format
}

my $cancergenelst=&ReadCancerGenelst($flist_cancergene);
my $db = "ucsc";
my $user = "mgg_admin";
my $password = "c\@nc3r";
my $dataBase = "DBI:mysql:$db:mysql2";
my $dbh = DBI->connect($dataBase, $user, $password) ||
  die "ERROR: Could not connect to database: $! \n";

my (%BK1s,%BK2s,%ROIs,%chrs);
my @SVs;
my $header;
foreach my $fin(@ARGV){
  my $FIN;
  if($fin eq '-'){
    $FIN=\*STDIN;
  }
  else{
    open(IN,"<$fin") || die "unable to open $fin\n";
    $FIN=\*IN;
  }
  while(<$FIN>){
    chomp;
    if(/^\#/){
      $header=$_ if(/^\#/ && /Chr/i);
      next;
    }
    chomp;
    my @u=split /\s+/;
    my ($chr,$start,$ori1,$chr2,$end,$ori2,$type,$size,$score,$nreads,$nread_lib,$AF,@extra)=@u;
    next if(defined $opts{o} && $chr ne $opts{o});

    $chr=$u[$opts{a}-1] if(defined $opts{a});
    $start=$u[$opts{b}-1] if(defined $opts{b});
    $start=~s/\D.*//gi;
    $chr2=$u[$opts{c}-1] if(defined $opts{c});
    $end=$u[$opts{d}-1] if(defined $opts{d});
    $end=~s/\D.*//gi;
    next unless($start=~/^\d+/ && $end=~/^\d+/ && defined $chr);
    push @SVs,$_;
    $BK1s{$chr}{$start}++;
    $BK2s{$chr2}{$end}++;
    $ROIs{$chr}{$start}=$end if($chr eq $chr2);
    $chrs{$chr}++; $chrs{$chr2}++;
  }
}

my %AROIs;
my %ABKs;
my %SBKs;
my %dbBK2s;
my %dbVarBK2s;
my %RPMK;
foreach my $chr(keys %chrs){
  next if(defined $opts{o} && $chr ne $opts{o});
  my $annot=&ReadUCSCGeneAnnotation($chr);

  #Select all transcripts overlapping with ROI
  my @starts=sort {$a<=>$b} (keys %{$ROIs{$chr}});
  my @txEnds=sort {$a<=>$b} keys %{$$annot{$chr}};
  #push @{$Annot{$e->{chrom}}{$e->{txEnd}}{$e->{txStart}}},$e;
  foreach my $start(@starts){
    my $end=$ROIs{$chr}{$start};
    while (@txEnds>0 && $start>$txEnds[0]){
      shift @txEnds;
    }
    next unless(@txEnds>0);
    foreach my $txEnd(@txEnds){
      my $hit=0;
      foreach my $txStart(keys %{$$annot{$chr}{$txEnd}}){
	if($end>=$txStart){
	  foreach my $transcript(@{$$annot{$chr}{$txEnd}{$txStart}}){
	    push @{$AROIs{$chr}{$start}{$end}},$transcript;
	  }
	  $hit=1;
	}
      }
      last if($hit == 0);
    }
  }

  #Select all transcripts that overlap with the breakpoint
  my @poses=sort {$a<=>$b} (keys %{$BK1s{$chr}}, keys %{$BK2s{$chr}});
  @txEnds=sort {$a<=>$b} keys %{$$annot{$chr}};
  foreach my $pos(@poses){
    while (@txEnds>0 && $pos>$txEnds[0]+$opts{l}){
      shift @txEnds;
    }
    next unless(@txEnds>0);
    foreach my $txEnd(@txEnds){
      foreach my $start(keys %{$$annot{$chr}{$txEnd}}){
	if($pos>$start-$opts{l}){
	  foreach my $transcript(@{$$annot{$chr}{$txEnd}{$start}}){
	    push @{$ABKs{$chr}{$pos}},$transcript;
	  }
	}
      }
    }
  }

  next if(defined $opts{M});  #For mouse, only apply gene annotation and ignore everything else

  #Select all SegDup that overlaps the breakpoint
  my $segdup=&ReadUCSCSegDupAnnotation($chr);
  my @segEnds=sort {$a<=>$b} keys %{$$segdup{$chr}};
  foreach my $pos(@poses){
    while (@segEnds>0 && $pos>$segEnds[0]+$opts{l}){
      shift @segEnds;
    }
    next unless(@segEnds>0);
    foreach my $start(keys %{$$segdup{$chr}{$segEnds[0]}}){
      if($pos>$start-$opts{l}){
	foreach my $dup(@{$$segdup{$chr}{$segEnds[0]}{$start}}){
	  push @{$SBKs{$chr}{$pos}},$dup;
	}
      }
    }
  }

  #Obtain dbSNP annotation
  @poses=sort {$a<=>$b} (keys %{$BK2s{$chr}});
  my $dbSNP=&Read_dbSNPAnnotation($chr);
  my @chromEnds=sort {$a<=>$b} keys %{$$dbSNP{$chr}};
  foreach my $pos(@poses){
    while (@chromEnds>0 && $pos>$chromEnds[0]+$opts{l}){
      shift @chromEnds;
    }
    next unless(@chromEnds>0);
    foreach my $start(keys %{$$dbSNP{$chr}{$chromEnds[0]}}){
      if($pos>=$start-$opts{l}){
	foreach my $var(@{$$dbSNP{$chr}{$chromEnds[0]}{$start}}){
	  push @{$dbBK2s{$chr}{$pos}},$var;
	}
      }
    }
  }

  my $dbVar=&Read_dbVarAnnotation($chr);
  @chromEnds=sort {$a<=>$b} keys %{$$dbVar{$chr}};
  foreach my $pos(@poses){
    while (@chromEnds>0 && $pos>$chromEnds[0]+$opts{l}){
      shift @chromEnds;
    }
    next unless(@chromEnds>0);
    foreach my $start(keys %{$$dbVar{$chr}{$chromEnds[0]}}){
      if($pos>=$start-$opts{l}){
	foreach my $var(@{$$dbVar{$chr}{$chromEnds[0]}{$start}}){
	  push @{$dbVarBK2s{$chr}{$pos}},$var;
	}
      }
    }
  }

  #Prepare repeatMasker table
  my $table = "chr$chr"."_rmsk";
  #  my $query = "SELECT genoStart, genoEnd, repClass
  #              	FROM $table
  #              	WHERE (repClass = 'Satellite' || repClass = 'Low_complexity'|| repClass = 'Simple_repeat' || repClass = 'SINE' || repClass = 'LINE' || repClass = 'LTR')
  #                    && genoEnd >= ? && genoStart <= ?
  #              	ORDER BY genoStart";
  my $query = "SELECT genoStart, genoEnd, repClass
              	FROM $table
              	WHERE genoEnd >= ? && genoStart <= ?
              	ORDER BY genoStart";

  $RPMK{$chr} = $dbh->prepare($query) || die "Could not prepare statement '$query': $DBI::errstr \n";

}


print "$header\tRefseqGene\tDataBases\tSegDup\tRepeat\tShortIndex\n" if(defined $header);
foreach my $sv(@SVs){
  my @u=split /\s+/,$sv;
  my ($chr,$start,$ori1,$chr2,$end,$ori2,$type,$size,$score,$nreads,$nread_lib,$AF,@extra)=@u;
  $chr=$u[$opts{a}-1] if(defined $opts{a});
  $start=$u[$opts{b}-1] if(defined $opts{b});
  $chr2=$u[$opts{c}-1] if(defined $opts{c});
  $end=$u[$opts{d}-1] if(defined $opts{d});

  my $geneAnnot=&GetGeneAnnotation($chr,$start,$chr2,$end,$type);
  print "$sv\t$geneAnnot";

  if(!defined $opts{M}){  #For human, apply the following crap
    my $dbSNPAnnot=&GetVarAnnotation($chr,$start,$chr2,$end,\%dbBK2s);
    my $dbVarAnnot=&GetVarAnnotation($chr,$start,$chr2,$end,\%dbVarBK2s);
    my $dbSegDupAnnot=&GetSegDupAnnotation($chr,$start,$chr2,$end);
    my $repeatAnnot=&GetRepeatMaskerAnnotation($chr,$start,$chr2,$end);

    printf "\t%s",join(',',$dbSNPAnnot,$dbVarAnnot);
    printf "\t%s",$dbSegDupAnnot;
    printf "\t%s",$repeatAnnot;
    printf "\tchr%s\:%d\-%d,chr%s\:%d\-%d",$chr,$start-500,$start+500,$chr2,$end-500,$end+500;
    if($chr eq $chr2){
      print ",chr$chr:$start\-$end";
    }
  }
  print "\n";
}

sub GetGeneAnnotation{
  my ($chr,$start,$chr2,$end,$type)=@_;

  #Overlapping genes
  my %Cancergenes;
  my %OverlapGenes;
  foreach my $e(@{$AROIs{$chr}{$start}{$end}}){
    $OverlapGenes{$e->{name2}}++;
    if(defined $$cancergenelst{uc($e->{name2})}){
      $Cancergenes{$e->{name2}}++;
    }
  }

  #Annotate Breakpoint
  my (@e1s,@e2s);
  if(defined $ABKs{$chr}{$start}){
    @e1s=@{$ABKs{$chr}{$start}};
  }
  if(defined $ABKs{$chr2}{$end}){
    @e2s=@{$ABKs{$chr2}{$end}};
  }

  #select the transcript (from multiple ones)
  my ($e1,$e2,$struct1,$struct2,$annot1,$annot2);
  for(my $i=0;$i<=$#e1s;$i++){
    for(my $j=0;$j<=$#e2s;$j++){
      if($e1s[$i]->{name} eq $e2s[$j]->{name}){  #same carrying transcript
	$e1=$e1s[$i];
	$e2=$e2s[$j];
      }
    }
  }
  if(!defined $e1){
    for(my $i=0;$i<=$#e1s;$i++){
      if(!defined $e1 || $e1->{txSize} < $e1s[$i]->{txSize}){  #longest isoform
	$e1=$e1s[$i];
      }
    }
  }
  if(!defined $e2){
    for(my $i=0;$i<=$#e2s;$i++){
      if(!defined $e2 || $e2->{txSize} < $e2s[$i]->{txSize}){  #longest isoform
	$e2=$e2s[$i];
      }
    }
  }
  $OverlapGenes{$e1->{name2}}++ if(defined $e1);
  $OverlapGenes{$e2->{name2}}++ if(defined $e2);

  $struct1=&AnnotExon($start,$e1) if(defined $e1);
  $annot1=(defined $struct1)?sprintf "%s%d:%d\/%d",$struct1->{unit},$struct1->{id},$struct1->{pos},$struct1->{size}:'NA';

  $struct2=&AnnotExon($end,$e2) if(defined $e2);
  $annot2=(defined $struct2)?sprintf "%s%d:%d\/%d",$struct2->{unit},$struct2->{id},$struct2->{pos},$struct2->{size}:'NA';

  my @overlapgenes=keys %OverlapGenes;
  my $gene=(@overlapgenes)?'Gene:'.join('|',@overlapgenes):'-';

  if(defined $e1 && defined $e2){
    if($e1->{name} eq $e2->{name}){  #same transcript
      $gene.=sprintf ",%s\|%s\:%s\-%s",$e1->{name}||'NA',$e1->{name2}||'NA',$annot1||'NA',$annot2||'NA';
      if(defined $struct1 && defined $struct2 && (abs($struct2->{id}-$struct1->{id})>0 ||
						  $struct1->{unit}=~/exon/i)
	){
	$gene.=',AffectCoding';
      }
      if(defined $struct1 && defined $struct2 && (abs($struct2->{id}-$struct1->{id})>1)
	){
	$gene.=',novelSplice';
      }
    }
    else{
      $gene.=sprintf ",%s\|%s\:%s\-%s\|%s\:%s",$e1->{name}||'NA',$e1->{name2}||'NA',$annot1||'NA',$e2->{name}||'NA',$e2->{name2}||'NA',$annot2||'NA';
      $gene.=',AffectCoding,Fusion';
    }
  }
  elsif(defined $e1 || defined $e2){
    $gene.=sprintf ",%s\|%s\:%s\-%s\|%s\:%s",$e1->{name}||'NA',$e1->{name2}||'NA',$annot1||'NA',$e2->{name}||'NA',$e2->{name2}||'NA',$annot2||'NA';
    $gene.=',AffectCoding';
  }
  else{}

  my @cgenes=keys %Cancergenes;
  $gene.=',Cancer:'.join('|',@cgenes) if($#cgenes>=0);
  return $gene;
}

sub GetRepeatMaskerAnnotation{
  my ($chr1,$pos1,$chr2,$pos2)=@_;
  my ($nbrpt1,$rep1)=&GetBKRepeatMaskerAnnotation($chr1,$pos1);
  my ($nbrpt2,$rep2)=&GetBKRepeatMaskerAnnotation($chr2,$pos2);
  my $repeatannot='-';
  if($nbrpt1>$opts{m} || $nbrpt2>$opts{m}){
    $repeatannot=sprintf "Repeat:%s-%s", $rep1 || 'NA',$rep2 ||'NA';
  }
  return $repeatannot;
}

sub GetBKRepeatMaskerAnnotation{
  my ($chr1,$pos)=@_;
  if($chr1=~/[MN]/){
    return (0,undef);
  }
  my $start=$pos-$opts{w};
  my $stop=$pos+$opts{w};
  $RPMK{$chr1}->execute($start, $stop) ||
    die "Could not execute statement for repeat masker table with (chr$chr1, $start, $stop): $DBI::errstr \n";
#  my %repeatCoords;
  my %repCount;
  while ( my ($chrStart, $chrStop, $repClass) =  $RPMK{$chr1}->fetchrow_array() ) {
    my $start_last = ($chrStart > $start) ? $chrStart : $start;
    my $stop_last = ($chrStop < $stop) ? $chrStop : $stop;
    #foreach ($start_last..$stop_last) { $repeatCoords{$_} = 1; $repCount{$repClass}++;}
    $repCount{$repClass}=$stop_last-$start_last+1 if($start_last-$opts{k}<=$pos && $pos<=$stop_last+$opts{k});
  }
#  my $repeatbase=0;
  my @sortedrepClass;
#  foreach ($start..$stop) {    if ( defined $repeatCoords{$_}  ) { $repeatbase++; }}
  my $maxClass;
  my $maxClassCount=0;
  foreach(keys %repCount){
    if(defined $repCount{$_} && $repCount{$_}>$maxClassCount){
      $maxClass=$_;
      $maxClassCount=$repCount{$_};
    }
  }
  return ($maxClassCount,$maxClass);
}

sub GetSegDupAnnotation{
  my ($chr,$start,$chr2,$end)=@_;

  #Annotate Breakpoint
  my (@e1s,@e2s);
  if(defined $SBKs{$chr}{$start}){
    @e1s=@{$SBKs{$chr}{$start}};
  }
  if(defined $SBKs{$chr2}{$end}){
    @e2s=@{$SBKs{$chr2}{$end}};
  }

  #select the transcript (from multiple ones)
  my ($e1,$e2,$struct1,$struct2,$annot1,$annot2);
  for(my $i=0;$i<=$#e1s;$i++){
    for(my $j=0;$j<=$#e2s;$j++){
      if($e1s[$i]->{name} eq $e2s[$j]->{name}){  #same carrying transcript
	$e1=$e1s[$i];
	$e2=$e2s[$j];
      }
    }
  }
  if(!defined $e1){
    for(my $i=0;$i<=$#e1s;$i++){
      if(!defined $e1 || $e1->{Size} < $e1s[$i]->{Size}){
	$e1=$e1s[$i];
      }
    }
  }
  if(!defined $e2){
    for(my $i=0;$i<=$#e2s;$i++){
      if(!defined $e2 || $e2->{Size} < $e2s[$i]->{Size}){
	$e2=$e2s[$i];
      }
    }
  }

  my $segdup='-';
  if(defined $e1 || defined $e2){
    $segdup=sprintf "SegDup:%s\-%s",$e1->{name}||'NA',$e2->{name}||'NA';
  }
  return $segdup;
}


sub GetVarAnnotation{
  my ($chr,$start,$chr2,$end,$db)=@_;
  return '-' if($chr ne $chr2);
  my @vars;
  if(defined $$db{$chr2}{$end}){
    @vars=@{$$db{$chr2}{$end}};
  }
  my $bestvar;
  my ($maxratio1,$maxratio2)=(0,0);
  foreach my $var(@vars){
    my $pos1=($end<$var->{chromEnd})?$end:$var->{chromEnd};
    my $pos2=($start>$var->{chromStart})?$start:$var->{chromStart};
    my $overlap=$pos1-$pos2+1;
    my $ratio1=$overlap/(abs($end-$start)+1);
    my $ratio2=$overlap/(abs($var->{chromEnd}-$var->{chromStart})+1);
    if($ratio1>=$opts{r} && $ratio2>=$opts{r} && ($ratio1>=$maxratio1 || $ratio2>=$maxratio2)){
      $bestvar=$var;
      $maxratio1=$ratio1;
      $maxratio2=$ratio2;
    }
  }

  my $varreport='-';
  if(defined $bestvar){
    $varreport=$bestvar->{name};
  }
  return $varreport;
}

sub AnnotExon{
  my ($pos,$e)=@_;
  my @exonStarts=split /\,/,$e->{exonStarts};
  my @exonEnds=split /\,/,$e->{exonEnds};

  my $report;
  for(my $i=0;$i<=$#exonStarts;$i++){
    if($exonStarts[$i]<=$pos && $pos<=$exonEnds[$i]){
      my $id=$i+1;
      my $rpos=$pos-$exonStarts[$i]+1;
      if($e->{strand} eq '-'){
	$id=$#exonStarts-$i+1;
	$rpos=$exonEnds[$i]-$pos+1;
      }
      ($report->{unit},$report->{id},$report->{pos},$report->{size})=('Exon',$id,$rpos,$exonEnds[$i]-$exonStarts[$i]+1);
    }
    elsif($i>0 && $exonEnds[$i-1]<$pos && $pos<$exonStarts[$i]){
      my $id=$i;
      my $rpos=$pos-$exonEnds[$i-1]+1;
      if($e->{strand} eq '-'){
	$id=$#exonStarts-$i+1;
	$rpos=$exonStarts[$i]-$pos+1;
      }
      ($report->{unit},$report->{id},$report->{pos},$report->{size})=('Intron',$id,$rpos,$exonStarts[$i]-$exonEnds[$i-1]+1);
    }
  }
  return $report;
}

sub ReadUCSCGeneAnnotation{
  my ($chr)=@_;
  my %Annot;
  open(AN,"<$opts{g}") || die "Unable to open $opts{g}\n";
  while(<AN>){
    chomp;
    next if(/^\#/);
    my $e;
    ($e->{bin},$e->{name},$e->{chrom},$e->{strand},$e->{txStart},$e->{txEnd},$e->{cdsStart},$e->{cdsEnd},$e->{exonCount},$e->{exonStarts},$e->{exonEnds},$e->{id},$e->{name2},$e->{cdsStartStat},$e->{cdsEndStat},$e->{exonFrames})=split;
    $e->{chrom}=~s/chr//;
    next unless($e->{chrom} eq $chr);
    $e->{txSize}=abs($e->{txStart}-$e->{txEnd}+1);
    push @{$Annot{$e->{chrom}}{$e->{txEnd}}{$e->{txStart}}},$e;
  }
  close(AN);
  return \%Annot;
}

sub ReadUCSCSegDupAnnotation{
  my ($chr)=@_;
  my %SegDup;
  open(SEGDUP,"<$opts{s}") || die "Unable to open $opts{s}\n";
  while(<SEGDUP>){
    chomp;
    next if(/^\#/);
    my $e;
    my @extra;
    ($e->{bin},$e->{chrom},$e->{Start},$e->{End},$e->{name},$e->{score},$e->{strand},@extra)=split;
    $e->{chrom}=~s/chr//;
    next unless($e->{chrom} eq $chr);
    $e->{Size}=abs($e->{Start}-$e->{End}+1);
    push @{$SegDup{$e->{chrom}}{$e->{End}}{$e->{Start}}},$e;
  }
  close(SEGDUP);
  return \%SegDup;
}

sub Read_dbSNPAnnotation{
  my ($chr)=@_;
  my %dbSNP;
  open(DBSNP,"<$opts{p}") || die "Unable to open $opts{p}\n";
  while(<DBSNP>){
    chomp;
    next if(/^\#/);
    my $p;
    ($p->{bin},$p->{chrom},$p->{chromStart},$p->{chromEnd},$p->{name},$p->{score},$p->{strand},$p->{refNCBI},$p->{refUCSC},$p->{observed},$p->{molType},$p->{class},$p->{valid},$p->{avHet},$p->{avHetSE},$p->{func},$p->{locType},$p->{weight})=split /\t+/;
    $p->{chrom}=~s/chr//;
    next unless($p->{chrom} eq $chr);
    push @{$dbSNP{$p->{chrom}}{$p->{chromEnd}}{$p->{chromStart}}},$p;
  }
  close(DBSNP);
  return \%dbSNP;
}

sub Read_dbVarAnnotation{
  my ($chr)=@_;
  my %dbVar;
  open(DBVAR,"<$opts{v}") || die "Unable to open $opts{v}\n";
  while(<DBVAR>){
    chomp;
    next if(/^\#/);
    my ($p,$db,$tmp1,$tmp2,$tmp3);
    ($p->{chrom},$db,$p->{var},$p->{chromStart},$p->{chromEnd},$tmp1,$tmp2,$tmp3,$p->{name})=split /\t+/;
    $p->{chrom}=~s/chr//;
    next unless($p->{chrom} eq $chr);
    push @{$dbVar{$p->{chrom}}{$p->{chromEnd}}{$p->{chromStart}}},$p;
  }
  close(DBVAR);
  return \%dbVar;
}

sub ReadCancerGenelst{
  my @fin=@_;
  my %genes;
  foreach my $f(@fin){
    open(FIN,"<$f") || die "unable to open $f\n";
    $_=<FIN>;
    while(<FIN>){
      next if(/^\#/);
      chomp;
      my @t=split;
      $genes{$t[0]}++;
    }
  }
  return \%genes;
}
