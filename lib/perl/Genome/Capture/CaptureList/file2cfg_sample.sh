#!/bin/tcsh

#set num = `echo $argv[1] | perl -ane '@a = split(":",$F[0]); print $#a+1'`;
	
set project_dir = $argv[1];
set project_name = $argv[2];
set project_chr = ".chr";
set number = $argv[3];

set assembly = "BreakDancer/";
set assembly_file_end = ".allchr.Q40.somatic.assembled.HQfiltered.csv";

set BreakDancer = "BreakDancer/";
set BreakDancer_file_end1 = ".chr";
set BreakDancer_file_end2 = ".sv"


set Pindel_end = "pindel/";
set Pindel_normal_end = "normal/";
set Pindel_tumor_end = "tumor/";
set Pindel_file = "adapted_deletions.csv.big_deletions";

set Hydra = "Hydra/";
set Hydra_file_end = "_Hydra.sv"

set CNA = "CNA/";
set CNA_file_end = "_CNA.seg";

set col_AS = "	\\t	1	2	4	6	7	NA	NA	NA	NA	10	12	NA  0";
set col_BD = "	\\t	0	1	3	4	6	2	5	11	12	7	8	10  100";
set col_PD = "	\\t	0	1	0	2	NA	NA	NA	NA	NA	3	NA	4   0";
set col_CNA = "	\\t 0	1	0	2	NA	NA	NA	5	7	3	9	NA  10000";
#									type n1	n2 n_cn	t_cn sz	sc	sp_r
set col_HD = "  \\t 0   1   2   3   NA  NA  NA  NA  NA  NA  NA  NA  175";

set AS = ".AS";
set BD = ".BD";
set PD = ".PD";
set CN = ".CN";
set HD = ".HD";

set output_dir = $argv[4];

set j = 0;
#while($j < $num)
#	set number = `echo $argv[1] $j | perl -ane '@a = split(":", $F[0]); print $a[$F[1]]'`;
#@ j = $j + 1;

#set number = $argv[1];
#echo ${number}
	set output_file = ${output_dir}/${project_name}${number};
	if (-e ${output_file}) then
		rm ${output_file}
	endif
	
	@ ever_come = 0
		
		# Assembly
		set file_name_AS = ${project_dir}/${BreakDancer}${project_name}${number}${assembly_file_end}
		echo ${file_name_AS}
		if(-e ${file_name_AS}) then
			set skip = `more ${file_name_AS} | perl -ane 'print "$_" if($F[0] =~ /^#/)' | wc -l`;
			if(${skip} == 0) then
				echo "${file_name_AS}";
			endif
			set tag = ${project_name}${number}${AS};
			if (${ever_come} == 0) then
				echo "${tag}\t${file_name_AS}\t${skip}${col_AS}" > ${output_file}
			else
				echo "${tag}\t${file_name_AS}\t${skip}${col_AS}" >> ${output_file}
			endif
			@ ever_come = 1
		else
			echo ${file_name_AS}
		endif
		

		# PD normal
		set file_name_PD_normal = ${project_dir}/${Pindel_end}${Pindel_normal_end}${Pindel_file};
            
		if(-e ${file_name_PD_normal}) then
			set tag = ${project_name}${number}${PD}N;
			if(${ever_come} == 0) then
				echo "${tag}\t${file_name_PD_normal}\t0${col_PD}" > ${output_file}
			else
				echo "${tag}\t${file_name_PD_normal}\t0${col_PD}" >> ${output_file}
			endif
			@ ever_come = 1
		else
		    echo $file_name_PD_normal
		endif
		
		# PD tumor
		set file_name_PD_tumor = ${project_dir}/${Pindel_end}${Pindel_tumor_end}${Pindel_file};
            
		if(-e ${file_name_PD_tumor}) then
			set tag = ${project_name}${number}${PD}T;		
			if(${ever_come} == 0) then
				echo "${tag}\t${file_name_PD_tumor}\t0${col_PD}" > ${output_file}
			else
				echo "${tag}\t${file_name_PD_tumor}\t0${col_PD}" >> ${output_file}
			endif
			@ ever_come = 1
		else
		    echo $file_name_PD_tumor
		endif

		# BD
	foreach i (1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y)
		set file_name_BD = ${project_dir}/${BreakDancer}${project_name}${number}${BreakDancer_file_end1}${i}${BreakDancer_file_end2}
		if (-e ${file_name_BD}) then
			set skip = `more ${file_name_BD} | perl -ane 'print "$_" if($F[0] =~ /^#/)' | wc -l`;
			if(${skip} == 0) then
				echo "${file_name_BD}";
			endif
			set tag = ${project_name}${number}${BD}${project_chr}${i};
			if (${ever_come} == 0) then
				echo "${tag}\t${file_name_BD}\t${skip}${col_BD}" > ${output_file}
			else
				echo "${tag}\t${file_name_BD}\t${skip}${col_BD}" >> ${output_file}
			endif
			@ ever_come = 1
		endif
	end				
			
		# CNA
		set file_name_CNA = ${project_dir}/${CNA}${project_name}${number}${CNA_file_end}
		if (-e ${file_name_CNA}) then
                    set skip = `grep "#CHR" -n -m 1 ${file_name_CNA} | perl -ane 'my @u = split(":", $F[0]); print($u[0])'`
			if(${skip} == 0) then
				echo "${file_name_CNA}";
			endif
			set tag = ${project_name}${number}${CN};
			if(${ever_come} == 0) then
				echo "${tag}\t${file_name_CNA}\t${skip}${col_CNA}" > ${output_file}
			else
				echo "${tag}\t${file_name_CNA}\t${skip}${col_CNA}" >> ${output_file}
			endif
			@ ever_come = 1	
		else
			echo ${file_name_CNA}		
		endif

                # Hydra
                set file_name_HD = ${project_dir}/${Hydra}${project_name}${number}${Hydra_file_end}
                if (-e ${file_name_HD}) then
                    echo ${file_name_HD};
                    set skip = `more ${file_name_HD} | perl -ane 'print "$_" if($F[0] =~ /^#/)' | wc -l`;
                    set tag = ${project_name}${number}${HD};
                    if(${ever_come} == 0) then
                        echo "${tag}\t${file_name_HD}\t${skip}${col_HD}" > ${output_file}
                    else
                        echo "${tag}\t${file_name_HD}\t${skip}${col_HD}" >> ${output_file}
                    endif
                    @ ever_come = 1
                else
                	echo ${file_name_HD}
                endif
                        
#end
	
