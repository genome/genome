#!/bin/tcsh

#set num = `echo $argv[1] | perl -ane '@a = split(":",$F[0]); print $#a+1'`;
#echo $num;
#cd /gscmnt/sata861/info/medseq/Xian/comparison/BRC50/perl_sample/assembly_added/final/no_filters;
#set j = 0;
#while($j < $num)
#	set i = `echo $argv[1] $j | perl -ane '@a = split(":", $F[0]); print $a[$F[1]]'`;

#set result_dir = $argv[1]
set no_filters = "/no_filters";
set result_dir = $argv[1]${no_filters};
echo ${result_dir}
set project_name = $argv[2];
set i = $argv[3];
set user = $argv[4];
	set f = ${result_dir}/${project_name}${i};
	if(-e $f.CR2.tmp) then
	    rm -f $f.CR2.tmp
	endif
	# CR2
	`grep BD $f | grep -v AS | perl -ane '$a=0; foreach ( @F[1..$#F] ) { ( $ncn,$tcn ) = ( $_=~/Ncn(\S+)\:Tcn(\S+)/ ) ; ( $tp ) = ($_=~/tp(\S+):/); $a = 1 if($ncn - $tcn > 0.5 && $tp =~ /DEL/ || $tcn-$ncn>0.5 && $tp =~ /ITX/);} print "$_" if($a);' > $f.CR2.tmp`
	# filter
	bsub -u xfan@genome.wustl.edu -J ${project_name}${i} -q 'short' "perl filters.pl $f.CR2.tmp ${f}.CR2 ${f}.CR2.filteredout";
#	@ j = ${j} + 1;		
#end
