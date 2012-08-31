package Genome::Model::Tools::CopyNumber::PurityHmm;

use strict;
use warnings;

my $LZERO=-10000000000;
my $LSMALL=$LZERO/2;
my $minLogExp = -log(-$LZERO);
my $PI=3.1415927;
my $f_TransProb=1e-6;

class Genome::Model::Tools::CopyNumber::PurityHmm {
    is => 'Command',
    has_input => [
        purity => {
            is => 'String',
            doc => 'purity',
            default => '0.8',
        },
        ploidy => {
            is => 'String',
            default => '2',
            doc => 'ploidy',
        },
        mean2x => {
            is => 'String',
            doc => 'mean 2x',
            default => "1",
        },
        var => {
            is => 'String',
            doc => 'var',
            default => 1e4,
        },
        max_cn => {
            is => 'String',
            doc => 'mx cn',
            default => '4',
        },
        cnv_rate => {
            is => 'String',
            doc => 'cnv rate',
            default => 0.0001,
        },
        cnv_size => {
            is => 'String',
            doc => 'cnv size',
            default => 0.1,
        },
        copy_number_file => {
            is => 'String',
            doc => 'Copy number file to operate on',
        },
    ],
    has_transient_optional => [
        ploidy_p => {
            doc => 'class variable',
        },
        ploidy_q => {
            doc => 'class variable',
        },
        avgState => {
            doc => 'class variable',
        },
        transP => {
            doc => 'class variable',
        },
        ceiling => {
            doc => 'class variable',
        },
    ],
};

sub printParm {
  my ($self,$fout,$ploidy)=@_;
  if(defined $self->ploidy_p && defined $self->ploidy_q && (!defined $ploidy)){
    printf $fout "purity=%.6f\tavg_state=%.6f\tp_arm_ploidy=%.6f\tq_arm_ploidy=%.6f\tmean2x=%.6f\tvar=%.6f\n",$self->purity,$self->avgState||-1,$self->ploidy_p,$self->ploidy_q,$self->mean2x,$self->var;
  }
  else{
    printf $fout "purity=%.6f\tavg_state=%.6f\tploidy=%.6f\tmean2x=%.6f\tvar=%.6f\n",$self->purity,$self->avgState||-1,$self->ploidy,$self->mean2x,$self->var;
  }
}

sub LogLikelihoodRatio {
  my ($self,$data,$i_start,$i_end,$CN)=@_;
  my ($logProbNeutral,$logProbCNV)=(0,0);

  for(my $t=$i_start;$t<=$i_end;$t++){
    my $dy=$self->Ceiling(${$data->{y}}[$t]);
    my $logprob=$self->LogGaussianPDF($dy,2);
    $logProbNeutral+=$logprob;
    $logprob=$self->LogGaussianPDF($dy,$CN);
    $logProbCNV+=$logprob;
  }
  my $log10Ratio=($logProbCNV-$logProbNeutral)/log(10);
  return $log10Ratio;
}

sub LogLikelihood {
  my ($self,$data,$i_start,$i_end,$CN)=@_;
  my $logProbCNV=0;

  for(my $t=$i_start;$t<=$i_end;$t++){
    my $dy=$self->Ceiling(${$data->{y}}[$t]);
    my $logprob=$self->LogGaussianPDF($dy,$CN);
    $logProbCNV+=$logprob;
  }
  $logProbCNV=$logProbCNV/log(10);
  return $logProbCNV;
}

sub Fit {
  my ($self,$data,$alg,$TPupdate,$PurityUpdate, $ofh)=@_;
  my $PLScore;
  my $Diff;
  my $LScore=$LZERO;
  my $iter=1;
  my $state;
  do{
    $PLScore=$LScore;
    if($alg && $alg=~/vit/i){
      ($LScore,$state)=$self->Viterbi($data,$TPupdate);
      $self->ViterbiUpdate($data,$state,$TPupdate,$PurityUpdate);
    }
    else{
      ($LScore)=$self->ForwardBackword($data,$TPupdate);
    }

    $Diff=$LScore-$PLScore;
    printf $ofh "Iter=%d\tLOGPROB=%.6f\tDIFF=%.6f\n",$iter,$LScore,$Diff;
    $iter++;
  } until(abs($Diff)<0.0001 || $iter>20);

  #$state=$self->SegPloidy($data);
  return $state;
}

sub SegPloidy {
  #Compute ploidy
  my ($self,$ofh,$data,$minsize,$Chr_Median_Readcount,$cmStart,$cmEnd)=@_;
  my ($LScore,$state)=$self->Viterbi($data);
  my $T=$#{$data->{x}} || 1;
  my $adjusted_meanreadcount=$self->purity*$self->mean2x;
  my $contamination_rate=1-$self->purity;
  my @adjusted_cns;
  for(my $i=0; $i<=$T;$i++){
    #my $adjusted_cn=2*(${$data->{y}}[$i]-$contamination_rate*$self->{mean2x})/$adjusted_meanreadcount;  #Li's equation
    my $tn_readcount_ratio=1;
    if(defined $$Chr_Median_Readcount{normal} && $$Chr_Median_Readcount{tumor} && $$Chr_Median_Readcount{normal}>0){
      $tn_readcount_ratio=$$Chr_Median_Readcount{tumor}/$$Chr_Median_Readcount{normal};
    }
    my $adjusted_cn;
    if(defined $data->{z}){
      $adjusted_cn=2*(${$data->{y}}[$i]-$contamination_rate*${$data->{z}}[$i]*$tn_readcount_ratio)/$adjusted_meanreadcount;
    }
    else{
      $adjusted_cn=2*${$data->{y}}[$i]/$adjusted_meanreadcount;
    }
    push @adjusted_cns,$adjusted_cn;
  }
  if(defined $cmStart && defined $cmEnd){  #compute Ploidy by chromosome arm
    my ($ploidy_stat_p,$ploidy_stat_q)=(0,0);
    my ($count_p,$count_q)=(0,0);
    my ($interval_p,$interval_q)=(0,0);
    my ($p_pos,$q_pos);
    for(my $i=0; $i<=$T;$i++){
      if(${$data->{x}}[$i]<$cmStart){
	#$ploidy_stat_p+=$$state[$i];
	$ploidy_stat_p+=$adjusted_cns[$i];
	$count_p++;
        $interval_p+=${$data->{x}}[$i];
        $interval_p-=$p_pos if(defined $p_pos);
        $p_pos=${$data->{x}}[$i];
      }
      elsif(${$data->{x}}[$i]>$cmEnd){
	#$ploidy_stat_q+=$$state[$i];
	$ploidy_stat_q+=$adjusted_cns[$i];
	$count_q++;
        $interval_q+=${$data->{x}}[$i];
	$interval_q-=$q_pos if(defined $q_pos);
 	$q_pos=${$data->{x}}[$i];
      }
      else{}
    }
    $self->ploidy_p(($interval_p>$minsize)?$ploidy_stat_p/$count_p:-1);
    $self->ploidy_q(($interval_q>$minsize)?$ploidy_stat_q/$count_q:-1);
    printf $ofh "p_interval: %d\tq_interval: %d\n",$interval_p,$interval_q;
  }
  else{
    my $ploidy_stat=0;
    my $ploidy_interval=0;
    for(my $i=0; $i<=$T;$i++){
      #$ploidy_stat+=$$state[$i];
      $ploidy_stat+=$adjusted_cns[$i];
      $ploidy_interval+=${$data->{x}}[$i];
      $ploidy_interval-=${$data->{x}}[$i-1] if($i>0);
    }
    $self->ploidy(($ploidy_interval>$minsize)?$ploidy_stat/$T:-1);
  }

  my $avg_stat=0;
  for(my $i=0; $i<=$T;$i++){
    $avg_stat+=$$state[$i];
  }
  $self->avgState(($T>10)?$avg_stat/$T:-1);
  return ($state,\@adjusted_cns);
}

sub ViterbiUpdate {
  my ($self,$data,$state,$TPupdate,$PurityUpdate)=@_;
  my $T=$#{$data->{y}}+1;

  my @transCount;
  my @mean_gc_state;
  my @num_state;
  my $pstate;

  #update purity
  my $purity_stat=0;
  my $nmarkers=0;
  for(my $i=0;$i<$T;$i++){
    my $j=$$state[$i];
    if($TPupdate){
      if(defined $pstate){
	$transCount[$pstate][$j]++;
      }
      $pstate=$j;
    }
    next if($j>=2);  #skip the copy number neutral and gain states
    #estimate purity using the deletions (see documentation purity_ploidy.lyx)
    $purity_stat+=($self->Ceiling(${$data->{y}}[$i])-$self->mean2x)/($self->mean2x) * 2/($j-2);
    $nmarkers++;
  }

  if($PurityUpdate){
    if($nmarkers<10){
      $self->error_message("Insufficient somatic alteration for purity estimation!\n");
    }
    else{
      $self->purity(($nmarkers>10)?$purity_stat/$nmarkers:$self->purity);
    }
  }

  #update mean
  my $mean_stat=0;
  $nmarkers=0;
  for(my $i=0;$i<$T;$i++){
    my $j=$$state[$i];
    my $r=1+($j/2-1)*$self->purity;
    next if($r<=0);
    $mean_stat+=$self->Ceiling(${$data->{y}}[$i])/$r;
  }
  $self->mean2x(($nmarkers>10)?$mean_stat/$nmarkers:$self->mean2x);

  #update var
  my $var_stat=0;
  $nmarkers=0;
  for(my $i=0;$i<$T;$i++){
    my $j=$$state[$i];
    my $m=$self->mean2x*($j/2*$self->purity+1-$self->purity);
    $var_stat+=($self->Ceiling(${$data->{y}}[$i])-$m)**2;
    $nmarkers++;
  }
  $self->var(($nmarkers>10)?$var_stat/($nmarkers-1):1);

  #update transition
  if($TPupdate){
    my @transP;
    my @transP_old=@{$self->transP};
    for(my $j=0;$j<=$self->max_cn;$j++){
      my $count=0;
      for(my $k=0;$k<=$self->max_cn;$k++){
	$count+=$transCount[$j][$k] || 0;
      }
      if($count>0){
	$transP[$j]=($count-$transCount[$j][$j])/$count;
      }
      else{
	$transP[$j]=$transP_old[$j];
      }
      $transP[$j]=($transP[$j]<$f_TransProb)?$f_TransProb:$transP[$j];   #floor transition prob
    }
    $self->transP(\@transP);
  }
}

sub Viterbi {
  my ($self,$data,$TPupdate)=@_;
  my @phi;
  my @delta;

  my @delta_t;
  for(my $j=0;$j<=$self->max_cn;$j++){
    #my $logObsProb=$self->LogGaussianPDF($y[0],$j);
    my $logObsProb=$self->LogGaussianPDF(${$data->{y}}[0],$j);
    push @delta_t,$logObsProb;
  }
  push @delta,\@delta_t;

  my @pdelta_t=@delta_t;
  my @transP=@{$self->transP};
  for(my $t=1;$t<=$#{$data->{y}};$t++){
    @delta_t=();
    my @phi_t;

    my $dist=(${$data->{x}}[$t]-${$data->{x}}[$t-1])/1000000;
    $dist=($dist<=0)?10:$dist;  #Assuming 10Mb distance if abnormal
    my $theta=exp(-2*$dist);

    for(my $i=0;$i<=$self->max_cn;$i++){
      my $logObsProb=$self->LogGaussianPDF(${$data->{y}}[$t],$i);
      my $maxlogprob=$LZERO;
      my $lab=0;

      for(my $j=0;$j<=$self->max_cn;$j++){  #previous time point
	my $transprob;
	if($TPupdate){
	  $transprob=($j==$i)?(1-$transP[$j]):$transP[$j]/($self->max_cn);
	}
	else{
	  $transprob=($j==$i)?$theta:(1-$theta)/($self->max_cn);
	}

	my $logprob=$pdelta_t[$j]+log($transprob);
	if($logprob>$maxlogprob){
	  $maxlogprob=$logprob;
	  $lab=$j;
	}
      }
      push @delta_t,$maxlogprob+$logObsProb;
      push @phi_t,$lab;
    }
    my @delta_tin=@delta_t;
    push @delta,\@delta_tin;

    push @phi,\@phi_t;
    @pdelta_t=@delta_t;
  }

  #trace back
  my @labs;
  my $maxlogprob=$LZERO;
  my $lab=0;
  for(my $j=0;$j<=$self->max_cn;$j++){  #previous time point
    my $logprob=$delta_t[$j];
    if($logprob>$maxlogprob){
      $maxlogprob=$logprob;
      $lab=$j;
    }
  }
  push @labs,$lab;
  my $LOGScore=$maxlogprob;

  my $mean_stat=0;
  my $var_stat=0;
  for(my $t=$#phi;$t>=0;$t--){
    $lab=$phi[$t][$lab];
    push @labs,$lab;
  }
  @labs=reverse(@labs);
  return ($LOGScore,\@labs);
}

sub ForwardBackword {
  my ($self,$data,$TPupdate)=@_;
  my @alpha;
  my @alpha_t0;
  my @beta;
  my @beta_t0;
  my @transP=@{$self->transP};

  #Forward procedure
  for(my $i=0;$i<=$self->max_cn;$i++){
    push @alpha_t0,log(1/$self->max_cn);  #uniform initial prob
  }

  @beta_t0=@alpha_t0;
  push @alpha,\@alpha_t0;  #t=0;

  my $T=$#{$data->{y}}+1;
  my @palpha_t=@alpha_t0;
  for(my $t=1;$t<=$T;$t++){
    my @alpha_t;

    my $dist=($t<$T)?(${$data->{x}}[$t]-${$data->{x}}[$t-1])/1000000:0;
    $dist=($dist<=0)?10:$dist;  #Assuming 10Mb distant in between if abnormal
    my $theta=exp(-2*$dist);

    for(my $j=0;$j<=$self->max_cn;$j++){
      my $logObsProb=$self->LogGaussianPDF(${$data->{y}}[$t-1],$j);
      my $sumlogprob=$LZERO;
      for(my $i=0;$i<=$self->max_cn;$i++){  #previous time point
	my $transprob;
	if($TPupdate){
	  $transprob=($j==$i)?(1-$transP[$i]):$transP[$i]/($self->{max_cn});
	}
	else{
	  $transprob=($j==$i)?$theta:(1-$theta)/($self->max_cn);
	}
	$transprob=($transprob<$f_TransProb)?$f_TransProb:$transprob;   #floor transition prob

	my $logprob=$palpha_t[$i]+log($transprob);
	$sumlogprob=&LAdd($sumlogprob,$logprob);
      }
      push @alpha_t,$sumlogprob+$logObsProb;
    }
    push @alpha,\@alpha_t;
    @palpha_t=@alpha_t;
  }

  #Backward procedure
  push @beta,\@beta_t0;  #t=T

  my @pbeta_t=@beta_t0;
  for(my $t=$T-1;$t>=0;$t--){
    my @beta_t;
    my $dist=($t==$T-1)?0:(${$data->{x}}[$t+1]-${$data->{x}}[$t])/1000000;
    $dist=($dist<=0)?10:$dist;  #Assuming 10Mb distant in between if abnormal
    my $theta=exp(-2*$dist);

    for(my $i=0;$i<=$self->max_cn;$i++){
      my $sumlogprob=$LZERO;
      for(my $j=0;$j<=$self->max_cn;$j++){  #previous time point
	my $logObsProb=$self->LogGaussianPDF(${$data->{y}}[$t],$j);

	my $transprob;
	if($TPupdate){
	  $transprob=($j==$i)?(1-$transP[$i]):$transP[$i]/($self->max_cn);
	}
	else{
	  $transprob=($j==$i)?$theta:(1-$theta)/($self->max_cn);
	}
	$transprob=($transprob<$f_TransProb)?$f_TransProb:$transprob;   #floor transition prob

	my $logprob=$pbeta_t[$j]+log($transprob)+$logObsProb;
	$sumlogprob=&LAdd($sumlogprob,$logprob);
      }
      push @beta_t,$sumlogprob;
    }
    push @beta,\@beta_t;
    @pbeta_t=@beta_t;
  }
  @beta=reverse(@beta);

  my $LOGScore=$LZERO;
  for(my $i=0;$i<=$self->max_cn;$i++){
    $LOGScore=&LAdd($LOGScore,$palpha_t[$i]);
  }

  #update parm;
  my ($varstat,$nsample,$purity_stat)=(0,0,0);
  my @gamma;
  my @eta_ij;
  my @eta_i;

  if($TPupdate){
    for(my $i=0;$i<=$self->max_cn;$i++){
      for(my $j=0;$j<=$self->max_cn;$j++){
	$eta_ij[$i][$j]=$LZERO;
      }
      $eta_i[$i]=$LZERO;
    }
  }

  for(my $t=1;$t<=$T;$t++){
    my $LogProbSum=$LZERO;
    my @gamma_t;

    my ($dist,$theta);
    if($TPupdate){
      $dist=($t==$T)?0:(${$data->{x}}[$t]-${$data->{x}}[$t-1])/1000000;
      $theta=exp(-2*$dist);
    }

    for(my $i=0;$i<=$self->max_cn;$i++){
      push @gamma_t,$alpha[$t][$i]+$beta[$t][$i];
      $LogProbSum=&LAdd($LogProbSum,$gamma_t[$i]);
      if($TPupdate && ($t<$T)){
	for(my $j=0;$j<=$self->max_cn;$j++){
	  my $logObsProb=$self->LogGaussianPDF(${$data->{y}}[$t],$j);
	  my $transprob=($j==$i)?(1-$transP[$i]):$transP[$i]/($self->max_cn);
	  $transprob=($transprob<$f_TransProb)?$f_TransProb:$transprob;   #floor transition prob

	  $eta_ij[$i][$j]=&LAdd($eta_ij[$i][$j],$alpha[$t][$i]+log($transprob)+$logObsProb+$beta[$t+1][$j]);
	}
	$eta_i[$i]=&LAdd($eta_i[$i],$gamma_t[$i]);
      }
    }
    for(my $i=0;$i<=$self->max_cn;$i++){
      $gamma_t[$i]=exp($gamma_t[$i]-$LogProbSum);
      next if($i==2 || abs(${$data->{y}}[$t-1])>2.5);  #skip copy-number neutral and outliers

      $purity_stat+=$gamma_t[$i]*${$data->{y}}[$t-1]/($i-2);
      $nsample+=$gamma_t[$i];
    }
    push @gamma,\@gamma_t;
  }

  #update purity
  if($nsample<10){
    $self->error_message("Insufficient data for purity estimation!\n");
  }
  else{
    $self->purity($purity_stat/$nsample);
  }

  for(my $t=1;$t<=$T;$t++){
    for(my $i=0;$i<=$self->max_cn;$i++){
      my $mean=$self->purity*($i-2);
      $varstat+=$gamma[$t-1][$i]*(${$data->{y}}[$t-1]-$mean)**2;
    }
  }
  #update var
  $self->var($varstat/$T);
  $self->var(($self->var<0.01)?0.01:$self->var);

  if($TPupdate){
    #update transP
    undef @transP;
    for(my $i=0;$i<=$self->max_cn;$i++){
      my $Sum;
      my $transSum=0;
      for(my $j=0;$j<=$self->max_cn;$j++){
	my $prob=exp($eta_ij[$i][$j]-$eta_i[$i]);
	$Sum+=$prob;
      next if($i==$j);
	$transSum+=$prob;
      }
      my $transprob=$transSum/$Sum;
      $transprob=($transprob<$f_TransProb)?$f_TransProb:$transprob;   #floor transition prob
      push @transP,$transprob;
    }

    $self->transP(\@transP);
  }

  return ($LOGScore);
}

sub LogGaussianPDF {
  my ($self,$x,$j)=@_;
  $x=$self->Ceiling($x);
  my $m=$self->mean2x*($j/2*$self->purity+1-$self->purity);
  my $v=$self->var;
  my $logp=-1/2*log(2*$PI*$v)-($x-$m)**2/(2*$v);
  return $logp;
}

sub Ceiling {
  my ($self,$y)=@_;
  $y=($y>$self->ceiling)?$self->ceiling:$y;
  return $y;
}

sub LAdd {
  my ($x, $y)=@_;
  my ($temp,$diff,$z);
  if ($x<$y) {
    $temp = $x; $x = $y; $y = $temp;
  }
  $diff = $y-$x;
  if ($diff<$minLogExp){
    return  ($x<$LSMALL)?$LZERO:$x;
  }
  else {
    $z = exp($diff);
    return $x+log(1.0+$z);
  }
  return $z;
}

sub init {
    my $self = shift;
    $self->ceiling($self->mean2x*$self->max_cn*2);
    my @transP;
    for(my $i=0;$i<=$self->max_cn;$i++){
        if($i==2){
            push @transP,$self->cnv_rate;
        }
        else{
            push @transP,$self->cnv_size;
        }
    }
    $self->transP(\@transP);
    return 1;
}

1;
