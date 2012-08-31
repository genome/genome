#!/bin/tcsh

#set num = `echo $argv[1] | perl -ane '@a = split(":",$F[0]); print $#a+1'`;
#echo $num;
#cd /gscmnt/sata861/info/medseq/Xian/comparison/BRC50/perl_sample/assembly_added/final/no_filters;
#set j = 0;
#while($j < $num)
#	set i = `echo $argv[1] $j | perl -ane '@a = split(":", $F[0]); print $a[$F[1]]'`;
#	set f = "BRC"$i;

#set result_dir = $argv[1];
set no_filters = "/no_filters";
set result_dir = $argv[1]${no_filters};
set project_name = $argv[2]
set number = $argv[3]
set f = ${project_name}${number};
cd ${result_dir}
set output_file_name = $f".capture";
set output_file = ${result_dir}/${output_file_name};
echo ${output_file}
if(-e ${output_file}) then
    rm ${output_file}
endif
echo "" > ${output_file};
	set f_ctx = $f"_ctx";
	echo "${f_ctx}";
	echo $f;
	set all = 0;
	# CR1
	set g1 = `grep AS $f | wc -l`;
	echo "\tCR1:" $g1;
        grep AS $f >> ${output_file};
	# CR2
	set g2 = `more ${f}.CR2 | wc -l`;
	echo "\tCR2:" $g2;
        cat ${f}.CR2 >> ${output_file};
	# CR3
	set g3 = `grep PD $f | grep CN | wc -l`
	echo "\tCR3:" $g3;
        grep PD $f | grep CN >> ${output_file};
	# CR4
	set g4 = `grep PD $f | grep -v CN | grep -v BD | grep -v AS | grep -v PDN | wc -l`;
	echo "\tCR4:" $g4;
        grep PD $f | grep -v CN | grep -v BD | grep -v AS | grep -v PDN >> ${output_file};
	# CR5
	set g5 = `grep CN $f | grep -v BD | grep -v PD | grep -v AS | perl -ane '($s)=($F[1]=~/sc(\S+)\:/); print "$_" if($s>=100);' | wc -l`
	echo "\tCR5:" $g5;
        grep CN $f | grep -v BD | grep -v PD | grep -v AS | perl -ane '($s)=($F[1]=~/sc(\S+)\:/); print "$_" if($s>=100);' >> ${output_file};
    # HD
    set hd_intra = `grep HD $f | grep -v BD | grep -v PD | grep -v AS | grep -v CN | wc -l`;
    echo "\tHD:" $hd_intra;
	    grep HD $f | grep -v BD | grep -v PD | grep -v AS | grep -v CN >> ${output_file};
	# CTX
	set ctx = `more ${f_ctx} | wc -l`;
	echo "\tCTX:" $ctx;
        cat ${f_ctx} >> ${output_file};
	@ all = $all + $g5 + $g4 + $g3 + $g2 + $g1 + $ctx;
	echo "\tIn all:" $all;
#	@ j = $j + 1;	
#end
