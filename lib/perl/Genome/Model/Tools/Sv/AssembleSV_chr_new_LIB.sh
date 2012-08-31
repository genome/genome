#!/gsc/bin/sh
ID=$1;  #Sample ID, e.g., BRC6, AML2
LIB=$2;  #Normal Libraries
NOVO=$3;
DIR=$4;
SCRIPT_DIR="$( cd "$( dirname "$0" )" && pwd )"
if [ $# == 4 ]; then
    MOUSE=0;
else
    MOUSE=$5;
fi

cd ${DIR};
PATH_=`pwd`;
echo $PATH_;

if [ ! -e ${ID}$NOVO.ctx.annot ]; then
    ${SCRIPT_DIR}/BreakAnnot.pl ${ID}$NOVO.ctx > ${ID}$NOVO.ctx.annot
fi
verbalCheck="AssemblyValidation finished ok.";
for s in normal tumor; do 
    mkdir assembly_intra_$s;
    mkdir assembly_inter_$s;
    for i in `seq 1 22` X Y ; do 
	f="$ID.chr$i.Q40.somatic.assembled.$s.NEW"; 
        # put the check of jobs ahead of logs
        J=`bjobs -J cf${ID}c${i}_${s}_NEW -u all`;
        if [ "$J" == "" ]; then
    	    if [ -e $f.log ]; then
	        e=`tail -n 1  $f.log`; 
    	        if [[ "$e" != "$verbalCheck" ]]; then 
		    rm -f $f.log ; 
	        fi;
	    fi;
#	if [ -e $f.csv ]; then
#	    w=`wc -l $f.csv | awk '{print $1;}'`;
#	    if [[ "$w" -lt "2" ]]; then
#		rm -f $f.log;
#	    fi
#	fi
	    if [ ! -e $f.log ]; then
	        if [ -e $f.csv ]; then
		    rm -f $f.csv;
	        fi
	        if [ -e $f.fasta ]; then
		    rm -f $f.fasta;
	        fi
	        echo ${f}.csv;
	        if [ ! -e $ID.chr$i.sv.annot ]; then
		    ${SCRIPT_DIR}/BreakAnnot.pl $ID.chr$i.sv > $ID.chr$i.sv.annot
	        fi
                #J=`bjobs -J cf${ID}c${i}_${s}_NEW -u all`;
                #if [ "$J" == "" ]; then
                if [ $MOUSE == 1 ]; then
                    bsub -q apipe -N -u $USER@watson.wustl.edu -M 8000000 -R "select[type==LINUX64] select[mem>8000] rusage[mem=8000]" -e $f.log -J cf${ID}c${i}_${s}_NEW "gmt sv assembly-validation --bam-files $s.bam --sv-file $ID.chr$i.sv.annot --breakpoint-seq-file ${f}.fasta --cm-aln-file ${f}.cm --min-breakdancer-score 40 --intermediate-read-dir ${PATH_}/assembly_intra_$s --specify-chr $i --output-file ${f}.csv --min-size-of-confirm-asm-sv 10 --skip-libraries ${LIB} --assemble-mouse --reference-file /gscmnt/839/info/medseq/reference_sequences/NCBI-mouse-build37/all_sequences.fa";
                fi
                if [ $MOUSE != 1 ]; then
                    bsub -q apipe -N -u $USER@watson.wustl.edu -M 8000000 -R "select[type==LINUX64] select[mem>8000] rusage[mem=8000]" -e $f.log -J cf${ID}c${i}_${s}_NEW "gmt sv assembly-validation --bam-files $s.bam --sv-file $ID.chr$i.sv.annot --breakpoint-seq-file ${f}.fasta --cm-aln-file ${f}.cm --min-breakdancer-score 40 --intermediate-read-dir ${PATH_}/assembly_intra_$s --specify-chr $i --output-file ${f}.csv --min-size-of-confirm-asm-sv 10 --skip-libraries ${LIB}";
                fi
            fi
        fi

        f="$ID.ctx.chr$i.Q40.somatic.assembled.$s.NEW"; 
        # put the check of jobs ahead of logs
        K=`bjobs -J cf${ID}ctx${i}_${s}_NEW -u all`;
        if [ "$K" == "" ]; then
            if [ -e $f.log ]; then
                e=`tail -n 1  $f.log`; 
                if [[ "$e" != "$verbalCheck" ]]; then 
                    rm -f $f.log ; 
                fi;
            fi;
            #	if [ -e $f.csv ]; then
            #	    w=`wc -l $f.csv | awk '{print $1;}'`;
            #	    if [[ "$w" -lt "2" ]]; then
            #		rm -f $f.log;
            #	    fi
            #	fi
            if [ ! -e $f.log ]; then
                if [ -e $f.csv ]; then
                    rm $f.csv;
                fi
                if [ -e $f.fasta ]; then
                    rm $f.fasta;
                fi
                echo $f.csv
                #K=`bjobs -J cf${ID}ctx${i}_${s}_NEW -u all`;
                #if [ "$K" == "" ]; then
                if [ $MOUSE == 1 ]; then
                    bsub -q apipe -N -u $USER@watson.wustl.edu -M 8000000 -R "select[type==LINUX64] select[mem>8000] rusage[mem=8000]" -e $f.log -J cf${ID}ctx${i}_${s}_NEW "gmt sv assembly-validation --bam-files $s.bam --sv-file ${ID}$NOVO.ctx.annot --specify-chr $i --breakpoint-seq-file ${f}.fasta --cm-aln-file ${f}.cm --min-breakdancer-score 40 --intermediate-read-dir ${PATH_}/assembly_inter_$s --output-file ${f}.csv --min-size-of-confirm-asm-sv 10 --skip-libraries ${LIB} --assemble-mouse --reference-file /gscmnt/839/info/medseq/reference_sequences/NCBI-mouse-build37/all_sequences.fa";
                fi
                if [ $MOUSE != 1 ]; then
                    bsub -q apipe -N -u $USER@watson.wustl.edu -M 8000000 -R "select[type==LINUX64] select[mem>8000] rusage[mem=8000]" -e $f.log -J cf${ID}ctx${i}_${s}_NEW "gmt sv assembly-validation --bam-files $s.bam --sv-file ${ID}$NOVO.ctx.annot --specify-chr $i --breakpoint-seq-file ${f}.fasta --cm-aln-file ${f}.cm --min-breakdancer-score 40 --intermediate-read-dir ${PATH_}/assembly_inter_$s --output-file ${f}.csv --min-size-of-confirm-asm-sv 10 --skip-libraries ${LIB}";
                fi

            fi
        fi
    done
done
cd ..;
