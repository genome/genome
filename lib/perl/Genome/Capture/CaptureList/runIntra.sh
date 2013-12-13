#!/bin/tcsh

set dir = $argv[1];
set name = $argv[2];
set i = $argv[3];
set output_dir = $argv[4];
#set input = $argv[1];
#set num = `echo $argv[1] | perl -ane '@a = split(":",$F[0]); print $#a+1'`;
#echo $num;
#set j = 0;
#while($j < $num)
	#set i = `echo $argv[1] $j | perl -ane '@a = split(":", $F[0]); print $a[$F[1]]'`;
#	echo $i;
	#if($argv[2] == 1) then
        #generate the configure file
		./file2cfg_sample.sh ${dir} ${name} ${i} ${output_dir}
#	endif
if(!(-e "${output_dir}/no_filters")) then
    mkdir ${output_dir}/no_filters;
endif
	bsub -u $argv[5]@${GENOME_EMAIL_DOMAIN} -J ${name}${i} -q 'short' "perl bigComparison_final_simple_more.pl -h -n -r 0.75 ${output_dir}/${name}${i} > ${output_dir}/no_filters/${name}${i}"
#	@ j = ${j} ;	
#end	

#foreach i (18 20 22 35 38 48 49)
#	if($argv[2] == 1) then
#		~/myscript/SVsummary/file2cfg_sample.sh ${i} > BRC${i}
#	endif
#	bsub -u xfan@genome.wustl.edu -J BRC${i} -q 'short' "perl ~/myscript/SVsummary/bigComparison_final_simple_more.pl -h -n -r 0.75 BRC${i} > no_filters/BRC${i}.1"
#end
