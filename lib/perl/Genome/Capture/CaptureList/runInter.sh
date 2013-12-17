#!/bin/tcsh
#set num = `echo $argv[1] | perl -ane '@a = split(":",$F[0]); print $#a+1'`;

set dir = $argv[1];
set name = $argv[2];
set i = $argv[3];
set output_dir = $argv[4];
#set j = 0;
#while($j < $num)
	#set i = `echo $argv[1] $j | perl -ane '@a = split(":", $F[0]); print $a[$F[1]]'`;
#	if($argv[2] == 1) then
# generate the configure file
		./file2cfg_ctx_sample.sh ${dir} ${name} ${i} ${output_dir}
#	endif
if(!(-e "${output_dir}/no_filters")) then
    mkdir ${output_dir}/no_filters;
endif
	bsub -u $argv[5]@${GENOME_EMAIL_DOMAIN} -J ${name}${i}_ctx -q 'short' "perl bigComparison_final_simple_more.pl -h -n -x 1 ${output_dir}/${name}${i}_ctx ${dir}/${name}${i}/BreakDancer/${name}${i}.chr1.sv > ${output_dir}/no_filters/${name}${i}_ctx"
#	@ j = $j + 1;	
#end
