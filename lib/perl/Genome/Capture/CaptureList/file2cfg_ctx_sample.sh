#!/bin/tcsh

set ctx = "_ctx";

set project_dir = $argv[1];
set project_name = $argv[2];
set project_chr = ".chr";
set number = $argv[3];

set assembly = "BreakDancer/";
set assembly_file_end = ".ctx.Q40.somatic.assembled.HQfiltered.csv";

set BreakDancer = "BreakDancer/";
set BreakDancer_file_end1 = ".novo.ctx";

set Hydra = "Hydra/";
set Hydra_file_end = "_Hydra.ctx"

set col_AS = "	\\t	1	2	4	6	7	NA	NA	NA	NA	10	12	NA  0";
set col_BD = "	\\t	0	1	3	4	6	2	5	NA	NA	7	8	10  100";
#									type n1	n2 n_cn	t_cn sz	sc	
set col_HD = "  \\t 0   1   2   3   NA  NA  NA  NA  NA  NA  NA  NA  175";

set AS = ".AS";
set BD = ".BD";
set HD = ".HD";

set output_dir = $argv[4];

#@ number = 0

#while(${number} <= 52 )
#	@ number = ${number} + 1
#set number = $argv[1];
#echo ${number}
	set output_file = ${output_dir}/${project_name}${number}${ctx};
	if (-e ${output_file}) then
		rm ${output_file}
	endif
	
	@ ever_come = 0
		
		# Assembly
		set file_name_AS_all = ${project_dir}/${BreakDancer}${project_name}${number}${assembly_file_end}
		echo ${file_name_AS_all}
                set file_name_AS = ${project_dir}/${BreakDancer}${project_name}${number}${assembly_file_end}".somaticCap";
                if(-e ${file_name_AS}) then
                    rm ${file_name_AS}
                endif
                `more ${file_name_AS_all} | perl -ane 'print "$_" if($F[11]!~/normal/)' > ${file_name_AS}`;
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
		
		# BD
	#foreach i (1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y)
		set file_name_BD = ${project_dir}/${BreakDancer}${project_name}${number}${BreakDancer_file_end1}
		if (-e ${file_name_BD}) then
#			set skip = `more ${file_name_BD} | perl -ane 'print "$_" if($F[0] =~ /^#/)' | wc -l`;
			if(${skip} == 0) then
#				echo "${file_name_BD}";
			endif
			set tag = ${project_name}${number}${BD};
			if (${ever_come} == 0) then
#				echo "${tag}\t${file_name_BD}\t${skip}${col_BD}" > ${output_file}
			else
#				echo "${tag}\t${file_name_BD}\t${skip}${col_BD}" >> ${output_file}
			endif
			@ ever_come = 1
		endif
	#end				

        set file_name_HD = ${project_dir}/${Hydra}${project_name}${number}${Hydra_file_end}
        if(-e ${file_name_HD}) then
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
	
