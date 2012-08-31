#!/usr/bin/env genome-perl
#Compare the first file with the rest of the files in terms of overlap in coordinates
use strict;
use warnings;
use Getopt::Std;
use DBI;
use lib "/gscuser/jwallis/svn/perl_modules/test_project/jwallis";

my %opts = (s=>0,q=>0,l=>0,m=>0,x=>0);
getopts('s:z:a:b:c:r:q:f:l:e:m:dpnhx:y:o:g', \%opts);
die("
Usage:   CompareInterval.pl <cfg file contains a list of files and parsing instruction (name,filename,skip,delimiter,chr_col,start_col,end_col, [score, size])> | <or a breakdancer file when -d is used>\n
Options:
         -f STRING specify the first file to be compared
         -d        compare deletions in the BreakDancer file to a set of references (DGV etc.)
         -s INT    which set to compare against the rest [$opts{s}]
         -z INT    which set to be compared to
         -l INT    allow +-[$opts{l}] bp ambiguity in positioning
         -m INT    ignore SVs of size shorter than this bp
         -e INT    size of the SVs (when specified) can not differ by more than this bp
         -c STRING chromosome
         -a INT    start position
         -b INT    end position
         -q INT    confidence score cutoff [$opts{q}]
	 -p	   print only those overlapped, by default off
         -r FLOAT  overlap ratio
	 -n	   print only overlapped samples
	 -h	   to print in simple VCR version
	 -y FLOAT  to print only those with at least one |t_cn-n_cn| >= this num [$opts{y}]
	 -x INT	   to print ctx and the allowed ambiguity is [$opts{x}]
	 -o	   one sided for ctx
	 -g		if to eliminate centromere
") unless (@ARGV);

my @sets;
my @tags;
my $f_cfg=$ARGV[0];
my ($selset,$seltag);
my $idx=0;
if($opts{d}){
  $seltag="BreakDancer";
  $selset=&ReadRegions($ARGV[0],1,'\t',0,1,4,7,8);
  $f_cfg='/gscmnt/sata194/info/sralign/kchen/1000genomes/analysis/scripts/compareLargeIndel.cmg';
  $idx++;
}

my $centro_read = 0;
my %CenStart;
my %CenEnd;
my $assembly_hit_cen = 0;
my $breakdancer_hit_cen = 0;
my $BreakpointBuffer = 500;
my $maxFraction = 0.8;
my %satelliteToChrStatements;

my $span;

my $num_col_ = 0;
my $num_col_1 = 0;
my @num_col;

open(CONFIG,"<$f_cfg") || die "unable to open $f_cfg\n";
while(<CONFIG>){
  next unless (/\S+/);
#  print $_."\n";
  my ($tag,$file,$skip,$delimiter,$chr1,$start,$chr2,$stop,$type,$nreads1,$nreads2,$normal_cn,$tumor_cn,$size,$score,$support_reads,$span_size)=split;
  if($idx==0 && $opts{f}){   #replace the first file, format specification kept
    $file=$opts{f};
  }
  if($idx==$opts{s}){
    $seltag=$tag;
    $selset=&ReadRegions($file,$skip,$delimiter,$chr1,$start,$chr2,$stop,$type,$nreads1,$nreads2,$normal_cn,$tumor_cn,$size,$score,$support_reads,$span_size);
    $num_col_1 = $num_col_;
  }
  elsif($opts{z}){
    if($opts{z}=~/$idx/){
      push @tags,$tag;
      push @sets,&ReadRegions($file,$skip,$delimiter,$chr1,$start,$chr2,$stop,$type,$nreads1,$nreads2,$normal_cn,$tumor_cn,$size,$score,$support_reads,$span_size);
      push @num_col, $num_col_;
    }
  }
  else{
    push @tags,$tag;
    push @sets,&ReadRegions($file,$skip,$delimiter,$chr1,$start,$chr2,$stop,$type,$nreads1,$nreads2,$normal_cn,$tumor_cn,$size,$score,$support_reads,$span_size);
    push @num_col, $num_col_;
  }
  $idx++;
}
close(CONFIG);

print "Assembly hit centromere: $assembly_hit_cen\tBreakDancer hit centromere: $breakdancer_hit_cen\n" if($opts{g});
unshift @tags,$seltag;
unshift @sets,$selset;
unshift @num_col, $num_col_1;

my $nsets=$#sets+1;
my %concordance;
my @nSVs;

my @record_all;
my @record_all_to_print;
my @idx;
my @set_chr;
foreach my $chr(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y'){
  for(my $i=0;$i<$nsets;$i++){
    push @idx,0;
    $set_chr[$i]=${$sets[$i]}{$chr};
#    if($i == 1){
#    for(my $j = 0; $j < $#{$set_chr[$i]}; $j++){
#    print $j. "\t". ${$set_chr[$i]}[$j]->{chr1}."\t". ${$set_chr[$i]}[$j]->{start}. "\t".${$set_chr[$i]}[$j]->{end}."\n" if(${$set_chr[$i]}[$j]->{start} &&  ${$set_chr[$i]}[$j]->{start}=~ /536740/);
#    print "wrong!\n" if(${$set_chr[$i]}[$j]->{chr1} eq 2);
#    }
#    }
    $nSVs[$i]+=$#{$set_chr[$i]}+1;
  }
  #print $#{$set_chr[0]}."\n";
  my $i = 0;
  while($i < $nsets){
  while($idx[$i]<=$#{$set_chr[$i]}){
    my $Sreg=${$set_chr[$i]}[$idx[$i]];
#    print "Original:.\n";
#    &RegprintAll($tags[0],$Sreg);
    if(!$opts{x}){
	&iterate_func($i, $Sreg->{start}, $Sreg->{chr2});
    }
    else{
	&iterate_func($i, $Sreg->{end} - $opts{x}, $Sreg->{chr2});
    }
    
    #&print_iterate();
#    $idx[$i]++;
#    print "$idx[0]\n";
#    print $#{$set_chr[0]};
  }
  $i++;
  }
  

  &post_process_combine();
  &hierarchy_span_regions();
#  &print_iterate_no_check();
 # &print_iterate();
  @record_all = ();
  @record_all_to_print = ();
  @idx = ();
  @set_chr = ();
}

#for(my $i=0;$i<$nsets;$i++){
#  for(my $j = $i+1; $j < $nsets; $j++){
#        my $num1=(defined $concordance{$i}{$j} && defined $nSVs[$i])?$concordance{$i}{$j}*100/$nSVs[$i]:0;
#	my $num2=(defined $concordance{$i}{$j} && defined $nSVs[$i])?$concordance{$i}{$j}*100/$nSVs[$j]:0;
#	printf "# %s\:%s\t%d\/%d \(%.2f%%\)\t%d\/%d \(%.2f%%\)\n",$tags[$i],$tags[$j],$concordance{$i}{$j}||0,$nSVs[$i],$num1,$concordance{$i}{$j}||0,$nSVs[$j],$num2;
#  }
#}

# hierarchily span the regions by the 12th column and 13th column
sub hierarchy_span_regions{
  my $num = $#record_all_to_print + 1;
  for(my $i = 0; $i < $num; $i++){
  	my $index = $record_all_to_print[$i]->{ind};
  	chomp $index;
  	my @u = split(/\|/, $index);
  	
	my $to_print_tmp;
  	for(my $jj = 0; $jj < $#u+1; $jj++){
  	  my $reg = $record_all[$u[$jj]];
  	  
  	  for(my $j = 0; $j < $nsets; $j++){
  	    next if(!$reg->{$j});
  	    my $tmp = $reg->{$j};
		my $for_num = $to_print_tmp->{$j};
		my $num = scalar keys %$for_num;

		for my $key (sort {$a<=>$b} keys %$tmp){
			my $tmp_ = $tmp->{$key};
			$to_print_tmp->{$j}->{$num++} = $tmp_;
		}
	  }
	}
	
	# hierarchy determine outer start, inner start, inner stop, outer stop
	my ($outer_start, $inner_start, $inner_stop, $outer_stop);
	my $first_time = 1;
	my $level_processed;
	for(my $j = 0; $j < $nsets; $j++){
	  if($to_print_tmp->{$j}){
	    $level_processed = $j;
	    my $tmp = $to_print_tmp->{$j};
	    for my $key (sort {$a <=> $b} keys %$tmp){
	      my $tmp_ = $tmp->{$key};
	      my $start = $tmp_->{start};
	      my $end = $tmp_->{end};
	      if($first_time == 0){
	        $outer_start = $start if($outer_start > $start);
	        $inner_start = $start if($inner_start < $start);
	        $inner_stop = $end if($inner_stop > $end);
	        $outer_stop = $end if($outer_stop < $end);
	      }
	      if($first_time == 1){
	        $outer_start = $start;
	        $inner_start = $start;
	        $inner_stop = $end;
	        $outer_stop = $end;
	        $first_time = 0;
	      }
	    }
	    last;
	  }
	}
	
	# span the region
	$record_all_to_print[$i]->{outer_start} = $outer_start - $span->{$level_processed};
	$record_all_to_print[$i]->{inner_start} = $inner_start + $span->{$level_processed};
	$record_all_to_print[$i]->{inner_stop} = $inner_stop - $span->{$level_processed};
	$record_all_to_print[$i]->{outer_stop} = $outer_stop + $span->{$level_processed};			
	
	# print key (the first column with the range)
	&RegprintKey_span($record_all_to_print[$i]); 	      
	
	# print	the others
	for(my $j = 0; $j < $nsets; $j++){
	  if($to_print_tmp->{$j}){
		my $tmp = $to_print_tmp->{$j};
		for my $key (sort {$a<=>$b} keys %$tmp){
	 		my $tmp_ = $tmp->{$key};
			if($opts{h}){
	  			&RegprintVCR($tags[$j], $tmp_);
			}
			else{
	  			&RegprintAll($tags[$j], $tmp_);
			}
		}
	  }
	}	  	  
	print "\n";
	$to_print_tmp = ();
  }
}  	

# to print record_all from record_all_to_print
sub print_iterate_no_check{
  my $print_num = $#record_all_to_print + 1;
  for(my $i = 0; $i < $print_num; $i++){
	my $index = $record_all_to_print[$i]->{ind};
	chomp $index;
	my @u=split(/\|/, $index);
	&RegprintKey($record_all_to_print[$i]);
	# to organize before printing out (sort by lists)
	my $to_print_tmp;

	for(my $jj = 0; $jj < $#u+1; $jj++){
	  my $reg = $record_all[$u[$jj]];

	  for(my $j = 0; $j < $nsets; $j++){
		next if(!$reg->{$j});
		my $tmp = $reg->{$j};
		my $for_num = $to_print_tmp->{$j};
		my $num = scalar keys %$for_num;

		for my $key (sort {$a<=>$b} keys %$tmp){
			my $tmp_ = $tmp->{$key};
			$to_print_tmp->{$j}->{$num++} = $tmp_;
		}
	  }
	}

	for(my $j = 0; $j < $nsets; $j++){
	  if($to_print_tmp->{$j}){
		my $tmp = $to_print_tmp->{$j};
		for my $key (sort {$a<=>$b} keys %$tmp){
	 		my $tmp_ = $tmp->{$key};
			if($opts{h}){
	  			&RegprintVCR($tags[$j], $tmp_);
			}
			else{
	  			&RegprintAll($tags[$j], $tmp_);
			}
		}
	  }
	}

	print "\n";
	$to_print_tmp = ();
  }
}

sub post_process_combine{
  my $print_num = $#record_all + 1;
  my $push = 0;
  for(my $i = 0; $i < $print_num; $i++){
	my $reg = $record_all[$i];
	# not applying the copy number last filter now for BD
	#if($opts{y}){
 	#  if(&check_cn_alternative($reg) == 0){
	#	next;
	#  }
	#}
	if($reg->{num} == 1 && $opts{p}){
	  next;
	}
	my $tmp; # record everything to combine
	$tmp->{start} = $reg->{start};
	$tmp->{end} = $reg->{end};
	$tmp->{chr1} = $reg->{chr1};
	$tmp->{chr2} = $reg->{chr2};
	$tmp->{ind} = $i;
	if($push == 0){
	  push @record_all_to_print, $tmp;
	  $push = 1;
	}
	else{
	  my $record_num = $#record_all_to_print + 1;	
	  my $num = 0;
	  my $overlap = 0;
	  for(my $j = $record_num - 1; $j >= 0; $j--){
		$num ++;
		if($num > 100){
		  last;
		}
		if(&Overlap($tmp, $record_all_to_print[$j])){
		  $record_all_to_print[$j]->{ind} = $record_all_to_print[$j]->{ind}."|".$tmp->{ind};
		  $overlap = 1;
		  last;
		}
	  }
	  if($overlap == 0){
		push @record_all_to_print, $tmp;
	  }

	}
  }
}

sub check_cn_alternative{
  my ($reg) = @_;
  my $ret = 0;
  for(my $j = 0; $j < $nsets; $j++){
	if(defined $reg->{$j}){
	  my $Xreg = $reg->{$j};
	  for my $key (keys %$Xreg){
		my $tmp = $Xreg->{$key};
	  	my $normal_cn = $tmp->{normal_cn};
	  	my $tumor_cn = $tmp->{tumor_cn};
	 	#printf "%s\t%s\n", $normal_cn, $tumor_cn;
	  	if($normal_cn =~ /^\d+/ && $tumor_cn =~ /^\d+/){
			if(abs($tumor_cn - $normal_cn) >= $opts{y}){
		  		$ret = 1;
		  		last;
			}
	  	}
	  }
	  if($ret == 1){
		last;
	  }
	}
  }
  return $ret;
}

sub print_iterate{
  my $print_num = $#record_all+1;
  for(my $i = 0; $i < $print_num; $i++){
	my $reg = $record_all[$i];
	if($opts{y}){
	  if(&check_cn_alternative($reg) == 0){
		next;
	  }
	}
	
	# to see if to print key 
	if($reg->{num} == 1 && $opts{p}){
	  next;
	}	
	&RegprintKey($reg);

	for(my $j = 0; $j < $nsets; $j++){
	  if(! defined $reg->{$j}){
		if(! $opts{n}){
		  &RegprintAll_none($tags[$j], $j);
		}
	  }
	  else{
		my $tmp = $reg->{$j};
		for my $key (sort {$a<=>$b} keys %$tmp){
			my $tmp_ = $tmp->{$key};
			if($opts{h}){
				&RegprintVCR($tags[$j], $tmp_);
			}
			else{
				&RegprintAll($tags[$j], $tmp_);
			}
		}
	  }
	}
	print "\n";
  }
  @record_all = ();
}

# recursive function begin here
sub iterate_func{
  my ($i, $start, $chr) = @_;
  my $Sreg=${$set_chr[$i]}[$idx[$i]];
  my $Qreg;
  my $ever_come = 0;
  my $record;
  for(my $j=$i+1;$j<$nsets;$j++){ # for a particular 'parent' subnode, compare with every following subnode 
	my $Qreg=${$set_chr[$j]}[$idx[$j]];
	my $overlap = 0;
	if($idx[$j]<=$#{$set_chr[$j]}){
	  $overlap = &Overlap($Sreg,$Qreg);
	  my $start_cutoff;
	  if(!$opts{x}){
	  	$start_cutoff = $Sreg->{start}<$start ? $Sreg->{start}:$start;
	  }
	  else{
		$start_cutoff = $Sreg->{end} - $opts{x} < $start ? $Sreg->{end} - $opts{x}:$start;
	  }
	  $start_cutoff = $start_cutoff < 0 ? 0:$start_cutoff; # start_curoff is to prevent the recursive one going to down. it will stop when it touches the start_cutoff. start_cutoff is the smallest start of each one launching the following recursion for per chromosome. and it is smallest end - opts{x} for each one launching the following recursion for ctx. 
	  while(!$overlap && $idx[$j]<=$#{$set_chr[$j]} && $Qreg->{end}<$start_cutoff){
		my $start_;
		if(!$opts{x}){
		  $start_ = $start_cutoff < $Qreg->{start} ? $start_cutoff : $Qreg->{start};
		}
		else{
		  $start_ = $start_cutoff < $Qreg->{end}-$opts{x}?$start_cutoff:$Qreg->{end}-$opts{x};
		}
		$start_ = $start_<0?0:$start_; 
		&iterate_func($j, $start_, $Qreg->{chr2}); # launch another recursion
		#$idx[$j]++;
		$Qreg=${$set_chr[$j]}[$idx[$j]];
		$overlap = &Overlap($Sreg,$Qreg);
	  }
	  if($overlap){
		if($ever_come == 0){
		  $ever_come = 1;
		  $record->{chr1} = $Sreg->{chr1};
		  $record->{start} = $Sreg->{start};
		  $record->{chr2} = $Sreg->{chr2};
		  $record->{end} = $Sreg->{end};
	          $record->{$i}->{0} = $Sreg;
		  $record->{num} = 1;
		  $record->{parent} = $i;
#$count_all ++ if($record->{chr2} eq "1" && $i == 1);		  
		}
		$idx[$j]++;
		      #print $idx[$j]."\t".$Qreg->{start}."\t".$Qreg->{chr1}."\t". $Qreg->{chr2}."\n" if($Qreg->{chr1} && $Qreg->{chr2});# && $Qreg->{chr1} =~ /X/ && $Qreg->{chr2} =~ /10/);
		$record->{$j}->{0} = $Qreg;
		$record->{num} ++;	 
		$concordance{$i}{$j} ++; 
#$count_all ++ if($record->{chr2} eq "1" && $j == 1);		  		
	  }
	}
  }
  if($ever_come == 0){# && !$opts{p}){ not apply -p here but in print process
	my $exist = &combine_record_simple($Sreg, $i);
	if($exist){
	  $record->{$i}->{0} = $Sreg;
	  $record->{chr1} = $Sreg->{chr1};
	  $record->{start} = $Sreg->{start};	
	  $record->{chr2} = $Sreg->{chr2};
	  $record->{end} = $Sreg->{end};
	  #print $record->{start}."\t".$record->{end}."\t".$record->{chr1}."\n";
	  $record->{num} = 1;
	  $record->{parent} = $i;
#	  		  $count_all ++ if($record->{chr2} eq "1");
	  push @record_all, $record;
	}
  }
  else{
	push @record_all, $record;
  }
  $idx[$i] ++; # add up no matter combined or not  
}

sub combine_record_simple{
	my ($record, $j) = @_;
	my $exist = 1;
	# never overlapped with others, see if can be combined with any
	my $record_num = $#record_all + 1;
	my $i = $record_num - 1;
	my $iterate = 0;
	while($i >= 0 && $exist == 1){
	  $iterate ++;
	  #if(!$opts{x} && ($opts{r} && $record_all[$i]->{end} < $record->{start} + ($record->{end} - $record->{start})*$opts{r} && $record_all[$i]->{parent} == $j || !$opts{r} && $record_all[$i]->{end} < $record->{start} && $record_all[$i]->{parent} == $j || $record->{end} - $record->{start} > 1000000 && $iterate > 1000)){
	  #	last;
	  #}
	  #if($opts{x} && ($record->{end} - $opts{x} > $record_all[$i]->{end} || $record->{chr2} ne $record_all[$i]->{chr2} || $record->{chr1} ne $record_all[$i]->{chr1}) ){
	#	last;
	#  }
	  #if($record_all[$i]->{num}>1){
		if(&Overlap($record, $record_all[$i])){
		  my $record_j_hash_ref = $record_all[$i]->{$j};
		  my $record_j_hash_ref_size = scalar keys %$record_j_hash_ref;
		  $record_all[$i]->{$j}->{$record_j_hash_ref_size} = $record;
		  $exist = 0;
		}
	  #}
	  $i--;
	}
	return $exist;
}
	
	
sub RegprintKey{
	my ($reg) = @_;
	printf "%s\t%s\t%s\t%s\t%d\t",join('.',$reg->{chr1},$reg->{start},$reg->{chr2},$reg->{end}), $reg->{chr1}, $reg->{start}, $reg->{chr2}, $reg->{end};
}

sub RegprintKey_span{
	my ($reg) = @_;
	printf "%s\t",join('.',$reg->{chr1},$reg->{outer_start},$reg->{inner_start},$reg->{chr2},$reg->{inner_stop},$reg->{outer_stop});#, $reg->{chr1}, $reg->{start}, $reg->{chr2}, $reg->{end};
}

sub RegprintAll_none{
	my ($tag, $i) = @_;
	my $num = $num_col[$i];
	printf "%s\t", $tag;
	for(my $j = 0; $j < $num; $j++){
		printf "-\t";
	}
}  


sub Regprint{
  my ($tag,$reg)=@_;
  printf "%s\t%s\t%d\t%s\t%d\t",$tag,$reg->{chr1},$reg->{start},$reg->{chr2},$reg->{end};
}

sub RegprintAll{
  my ($tag,$reg)=@_;
  printf "%s\t%s\t",$tag,$reg->{all};
}

sub RegprintVCR{
  my ($tag,$reg)=@_;
  if(!$opts{x}){
	printf "%s\|%s\:%s\:%s\:%s",$tag,$reg->{chr1}, $reg->{start},$reg->{chr2},$reg->{end};
	printf "\:tp%s",$reg->{type} if(defined $reg->{type});
	printf "\:sz%s",$reg->{size} if(defined $reg->{size});
	printf "\:sc%s",$reg->{score} if(defined $reg->{score});
	printf "\:nrdF%s", $reg->{nreads1} if(defined $reg->{nreads1});# && !defined $reg->{normal});
	printf "\:nrdS%s", $reg->{nreads2} if(defined $reg->{nreads2});# && !defined $reg->{tumor});
	printf "\:Nsp%s", $reg->{normal} if(defined $reg->{normal});
	printf "\:Tsp%s", $reg->{tumor} if(defined $reg->{tumor});
	printf "\:Ncn%s", $reg->{normal_cn} if(defined $reg->{normal_cn});
	printf "\:Tcn%s", $reg->{tumor_cn} if(defined $reg->{tumor_cn});
	printf "\t";
  }
  else{
	printf "%s\|%s\:%s\:%s\:%s", $tag, $reg->{chr1}, $reg->{start}, $reg->{chr2}, $reg->{end};
	printf "\:nrdF%s", $reg->{nreads1} if(defined $reg->{nreads1});# && !defined $reg->{normal});
	printf "\:nrdS%s", $reg->{nreads2} if(defined $reg->{nreads2});# && !defined $reg->{tumor});
	printf "\:Nsp%s", $reg->{normal} if(defined $reg->{normal});
	printf "\:Tsp%s", $reg->{tumor} if(defined $reg->{tumor});
	printf "\t";
  }
}

sub Overlap{
  my ($reg1,$reg2)=@_;
  my $overlap=0;
  if(!$reg2->{end} || !$reg1->{start} || !$reg1->{end} || !$reg2->{start}){
    return $overlap;
  }
  if($opts{x} && ((!$reg1->{chr1} || !$reg1->{chr2} || !$reg2->{chr1} || !$reg2->{chr2}) ||  $reg1->{chr1} eq $reg1->{chr2} || $reg2->{chr1} eq $reg2->{chr2} ) ){ # to get rid of no chromosomes if ctx
    return $overlap;
  }
  if(!$opts{x} && ($reg1->{chr1} ne $reg1->{chr2} || $reg2->{chr1} ne $reg2->{chr2})){
      return $overlap;
  }

  if(!$opts{x}){ # check the overlap of per chromosome
  my $a=$reg2->{end}-$reg1->{start}+$opts{l};
  my $b=$reg1->{end}-$reg2->{start}+$opts{l};

  if($a>=0 && $b>=0){
    my $size1=abs($reg1->{end}-$reg1->{start})+1+$opts{l};
    my $size2=abs($reg2->{end}-$reg2->{start})+1+$opts{l};
    $overlap=($a<$b)?$a:$b;
    $overlap=($overlap<$size1)?$overlap:$size1;
    $overlap=($overlap<$size2)?$overlap:$size2;
    if(defined $opts{r}){
      my $r1=($size1>0)?$overlap/$size1:1;
      my $r2=($size2>0)?$overlap/$size2:1;
      $overlap=($r1>=$opts{r} && $r2>=$opts{r})?1:0;
      #$overlap=($r1>=$opts{r} && $r2>=$opts{r})?1:0;
    }
    if(defined $opts{e} && defined $reg1->{size} && defined $reg2->{size}){
      my $size_diff=abs($reg1->{size}-$reg2->{size});
      $overlap=0 if($size_diff > $opts{e});
    }
  }
  }
  else{ # check the overlap of the ctx
    if($reg1->{chr1} eq $reg2->{chr1} && $reg1->{chr2} eq $reg2->{chr2}){		
		if(!$opts{o} && abs($reg1->{start} - $reg2->{start}) < $opts{x} && abs($reg1->{end} - $reg2->{end}) < $opts{x} || $opts{o} && (abs($reg1->{start} - $reg2->{start}) < $opts{x} || abs($reg1->{end} - $reg2->{end}) < $opts{x})){
		  $overlap = 1;
		}
    }
  }
  return $overlap;
}

# from the header pass from library to bam file
sub parse_header_sp_reads{
  my ($header, $reg) = @_;
  # parse sp reads to libraries
  my $sp_reads = $reg->{sp_reads};
  my @u=split(/\:/, $sp_reads);
  for(my $i = 0; $i < $#u + 1; $i ++ ){
  	my @v = split(/\|/, $u[$i]);
  	if($#v + 1 != 2){
            #print "parsing sp reads failed: ". $sp_reads. "\n";
  		next;
  	}
  	#print "hello ".$v[0]. ": ". $v[1];
  	#print "\n";
  	$reg->{$header->{$v[0]}} += $v[1]   if(defined $header->{$v[0]} && defined $reg->{$header->{$v[0]}});
  	$reg->{$header->{$v[0]}} = $v[1]	if(defined $header->{$v[0]} && !defined $reg->{$header->{$v[0]}});
  }
  return $reg;
}

sub read_header{
  my ($string, $header) = @_;
#  print $string."\n";
  return $header if($string !~ /^#/ || $string =~ /^#Chr1/); # if not beginning with a header # or it's the header of the columns, return, nothing added to the header
  my @u = split(/\t/, $string);
  my $bam;
  my $library;
  for(my $i = 0; $i < $#u+1; $i++){
#  print "$u[$i]\n";
  	$bam = "tumor" if($u[$i] =~ /tumor.bam/);
  	$bam = "normal" if($u[$i] =~ /normal.bam/);
  	$library = $u[$i] if($u[$i] =~ /library:/);
#  	print "yes$library\n" if($u[$i] =~ /library:/);
  }
  $library =~ s/library:// if($library);
  return $header if(!$bam || !$library);
  if(!defined $header->{$library}){  	
	$header->{$library} = $bam;
  }
  else{
	if($header->{$library} ne $bam){
	  print "wrong! library ". $library . " belong to both bam file " . $bam . " and ". $header->{$library} . "\n";
	  return $header;
	}
  }
  return $header;
}

sub ReadRegions{
  my ($fin,$n_skip,$delimiter,@c)=@_;
#  print $fin. "\n";
  my %Regions;
  my %hRegions;
  open(SV,"<$fin") || die "unable to open $fin\n";
  my $BD = 0; # indicate if it's a breakdancer file
  my $record = 0; # indicate whether need to record the library and bam file
  my $header; # library to bam file
  for(my $i=0;$i<$n_skip;$i++){
  	$_=<SV>;    #skip $n_skip lines
  	$header = &read_header($_, $header) if($record == 1);
  	$record = 1 if($BD == 1 && $_=~ /#Library Statistics:/);
  	$BD = 1 if($_ =~ /^#Software: BreakDancerMax/);  		
  }
  if($record == 1 && (scalar keys %$header) == 0){
  	# BD but no header : novo ctx, read the configure file
  	my $source_config = $ARGV[1];
  	open(config_source, "<$source_config") || die "unable to open $source_config\n";
  	while(<config_source>){
	  	$header = &read_header($_, $header);
	}
  }
  my %chrcount;
  my $chr_hash;
  my $chr_num = 1;
  while(<SV>){
    chomp;
    my $reg;
    #my @u=split $delimiter;
    my @u = split(/\s+/,$_); 
    if($#u < 2 || $u[0] =~ /^#/){
	next;
    }
    
    $num_col_ = @u;
    $reg->{all}=$_;

    # basic forms
    ($reg->{chr1},$reg->{start},$reg->{chr2},$reg->{end})=($u[$c[0]],$u[$c[1]],$u[$c[2]],$u[$c[3]]);
    next if($reg->{start} !~ /^\d+$/ ||
	    $reg->{end} !~ /^\d+$/);
    if($opts{x} && $reg->{chr1} eq $reg->{chr2}){
        next;
    }
    if((!$opts{x}) && $reg->{chr1} ne $reg->{chr2}){
        next;
    }
    if($opts{h}){
	# extra ones
	for(my $p = 4; $p <= 12; $p++){
	#print $u[$c[$p]];
	#print "\n";
		if(defined $c[$p] && $c[$p]=~/^\d+/){
		    	$reg->{type} = $u[$c[$p]] if($p == 4);
			$reg->{nreads1} = $u[$c[$p]] if($p == 5);
			$reg->{nreads2} = $u[$c[$p]] if($p == 6);
			$reg->{normal_cn} = $u[$c[$p]] if($p == 7);
			$reg->{tumor_cn} = $u[$c[$p]] if($p == 8);
			$reg->{score} = $u[$c[$p]] if($p == 10);
			$reg->{size} = $u[$c[$p]] if($p == 9);
			$reg->{size} = $reg->{end} - $reg->{start} + 1 		if(!defined $reg->{size} && $reg->{start} =~ /^\d+/ && $reg->{end} =~ /^\d+/);
			if($reg->{score} && $reg->{score} =~ /,/){
				my @scores = split(/,/,$reg->{score});
				my $max_score = 0;
				for(my $i_sc = 0; $i_sc <= $#scores; $i_sc ++){
					$max_score = $scores[$i_sc] if($max_score < $scores[$i_sc]);
				}
				$reg->{score} = $max_score;
			}	
		}
		if(defined $c[$p] && $p == 11 && $c[$p] ne "NA" && defined $u[$c[$p]]) { # if NA, like CNA
			if($u[$c[$p]] !~ /^\d+/ && $u[$c[$p]] !~ /tumor/i && $u[$c[$p]] !~ /normal/i){ # like BD
				$reg->{sp_reads} = $u[$c[$p]];	
				$reg = &parse_header_sp_reads($header, $reg) if($reg->{sp_reads});
				$reg->{tumor} = 0 if(!$reg->{tumor});
				$reg->{normal} = 0 if(!$reg->{normal});			
			}
			elsif($u[$c[$p]] =~ /^\d+/) { # like Pindel
				$reg->{tumor} = $u[$c[$p]] if($fin =~ /tumor/);
				$reg->{normal} = $u[$c[$p]] if($fin =~ /normal/);
			}
			elsif($u[$c[$p]] =~ /tumor/i || $u[$c[$p]] =~ /normal/i){ # like SD and AS
				($reg->{tumor}) = ($u[$c[$p]] =~ /tumor.*?(\d+)/i);
				($reg->{normal}) = ($u[$c[$p]] =~ /normal.*?(\d+)/i);
			}
		}
		if(defined $c[$p] && $p == 12 && $c[$p] ne "NA" && ! defined $span->{$idx}){
			$span->{$idx} = $c[$p];
			#print "$idx\t" . $span->{$idx} . "\n";
		}
		
	}
    }
    
    #print $reg->{sp_reads}."\t".$reg->{tumor}."\t".$reg->{normal}."\n";
    $reg->{chr1}=~s/chr//;
    $reg->{chr1}=~s/\"//gi;
    $reg->{chr2}=~s/chr//;
    $reg->{chr2}=~s/\"//gi;

    next if(
    	! defined $reg->{start} ||
    	! defined $reg->{end} ||
	    $reg->{start} !~ /^\d+$/ ||
	    $reg->{end} !~ /^\d+$/ ||
	    defined $reg->{score} && $reg->{score} < $opts{q} ||
	    $opts{c} && ($reg->{chr2} ne $opts{c}) ||
	    $opts{a} && $reg->{start} < $opts{a} ||
	    $opts{b} && $reg->{end} >  $opts{b} ||
	    defined $reg->{size} && defined $opts{m} && $reg->{size} < $opts{m} && !$opts{x}
	   );
	   #print $reg->{score}. "\t" . $reg->{size} . "\n" if($reg->{score} && $reg->{size} && $fin =~ /Pindel/);
	   
	if($opts{g}){
		my $hit_cen = 0;
		if($centro_read == 0){
			$centro_read = 1;
			my $chr_cen;
			#read data table
			foreach $chr_cen (1..22,"X","Y") {
				my $db = "ucsc";
				my $user = "mgg_admin";
				my $password = "c\@nc3r";
				my $dataBase = "DBI:mysql:$db:mysql2";
				my $dbh = DBI->connect($dataBase, $user, $password) ||
				    die "ERROR: Could not connect to database: $! \n";
    
			    my $table = "chr$chr_cen"."_rmsk";
    			my $query = "SELECT genoStart, genoEnd
              	FROM $table
              	WHERE (repClass = 'Satellite' || repClass = 'Low_complexity'|| repClass = 'Simple_repeat')
                    && genoEnd >= ? && genoStart <= ?
              	ORDER BY genoStart";

			    $satelliteToChrStatements{$chr_cen} = $dbh->prepare($query) ||
				die "Could not prepare statement '$query': $DBI::errstr \n";
			}
		}
		# make a region around breakpoints
		my %regions = ();
		${$regions{$reg->{start}-$BreakpointBuffer}}{$reg->{start}+$BreakpointBuffer} = 1;
   	    ${$regions{$reg->{end}-$BreakpointBuffer}}{$reg->{end}+$BreakpointBuffer} = 1;
		my $satelliteRegion = 0;
	    my %outputRegion = ();
	    my $chrStart;
	    my $chrStop;
	    my $start;
	    my $stop;
	    my %satelliteCoords;
	    my $totalSatellite;
	    my $fraction;
	    foreach $start ( keys %regions ) {
			foreach $stop ( keys %{$regions{$start}} ) {
			    %satelliteCoords = ();
    			$totalSatellite = 0;
    	
    			$satelliteToChrStatements{$reg->{chr1}}->execute($start, $stop) ||
					die "Could not execute statement for repeat masker table with (chr$reg->{chr1}, $start, $stop): $DBI::errstr \n";
			    while ( ($chrStart, $chrStop) =  $satelliteToChrStatements{$reg->{chr1}}->fetchrow_array() ) {
   			    	my $start_last = ($chrStart > $start) ? $chrStart : $start;
			    	my $stop_last = ($chrStop < $stop) ? $chrStop : $stop;
					foreach ($start_last..$stop_last) { $satelliteCoords{$_} = 1; }
   				}
    			foreach ($start..$stop) {
					if ( defined $satelliteCoords{$_}  ) { $totalSatellite++; }
    			}
    			$fraction = $totalSatellite/($stop - $start + 1);
    			if ( $fraction > $maxFraction ) { $satelliteRegion = 1; ${$outputRegion{$start}}{$stop} = $fraction;  last; }#print "$reg->{chr1}\t$reg->{start}\t$reg->{end}\t$start\t$stop\t$fraction\n"; last;}
			}
			last if($satelliteRegion == 1);
    	}
		$assembly_hit_cen ++ if($satelliteRegion == 1 && $fin =~ /assembled/);
		$breakdancer_hit_cen ++ if($satelliteRegion == 1 && $fin =~ /.sv/);    	
    	next if($satelliteRegion == 1);
    }
	
    if(!$opts{x}){
	$reg->{size} = $reg->{end} - $reg->{start} + 1 if(! defined $reg->{size});
    	$hRegions{$reg->{chr2}}{join('.',$reg->{start},$reg->{end})}=$reg;
    }
    else{
	    $chr_hash->{$reg->{chr1}} = $chr_num ++ if(!defined $chr_hash->{chr1});
	    my $tmp_chr = $chr_hash->{$reg->{chr1}};
	    $tmp_chr = "0".$chr_hash->{$reg->{chr1}} if($tmp_chr =~ /\d+/ && $tmp_chr < 10);
		$hRegions{$reg->{chr2}}{join('.',$reg->{end},$tmp_chr.$reg->{start})}=$reg;
    }
    $chrcount{$reg->{chr2}}++;
#    print $reg->{chr2}."\n" if($reg->{chr2});
#      	print "hello\n" if($reg->{start} && $reg->{end} && $reg->{start} =~ /79789716/);
  }

#  print "done with $fin\n";
#  print $chrcount{1}."\n" if($chrcount{1});
	
  foreach my $chr(keys %hRegions){
    my @Regs;
    if(!$opts{x}){
    foreach my $coor(sort {$a <=> $b} keys %{$hRegions{$chr}}){
      push @Regs,$hRegions{$chr}{$coor};
    }
    }
    else{
	foreach my $coor(sort {$a <=> $b} keys %{$hRegions{$chr}}){
	  push @Regs, $hRegions{$chr}{$coor};
	}
    }
    $Regions{$chr}=\@Regs;
  }
  close(SV);
  return \%Regions;
}
