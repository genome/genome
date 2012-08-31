#!/gsc/bin/sh

ID=$1;
#scripts/fixCTXswapXY_chr.pl $ID
DIR=$2;
SCRIPT_DIR="$( cd "$( dirname "$0" )" && pwd )"
cd ${DIR};
if [ -s merge.NEW.index ]; then
    rm -f merge.NEW.index;
fi
if [ -s merge_ctx.NEW.index ]; then
    rm -f merge_ctx.NEW.index;
fi

for s in tumor normal; do 
    for i in `seq 1 22` X Y; do 
	if [ -s $ID.chr$i.Q40.somatic.assembled.$s.NEW.csv ]; then
	    echo ${s}$i $ID.chr$i.Q40.somatic.assembled.$s.NEW.csv $ID.chr$i.Q40.somatic.assembled.$s.NEW.fasta >> merge.NEW.index; 
	fi
	if [ -s $ID.ctx.chr$i.Q40.somatic.assembled.$s.NEW.csv ]; then
	    echo ${s}$i $ID.ctx.chr$i.Q40.somatic.assembled.$s.NEW.csv $ID.ctx.chr$i.Q40.somatic.assembled.$s.NEW.fasta >> merge_ctx.NEW.index; 
	fi
    done
done

${SCRIPT_DIR}/MergeAssembledCallsets.pl -c -f $ID.allchr.Q40.somatic.assembled.HQfiltered.NEW.fasta -d 200 -h merge.NEW.index 1> $ID.allchr.Q40.somatic.assembled.HQfiltered.NEW.csv 2>$ID.allchr.Q40.somatic.assembled.HQfiltered_out.NEW.csv
${SCRIPT_DIR}/BreakAnnot.pl -A $ID.allchr.Q40.somatic.assembled.HQfiltered.NEW.csv > $ID.allchr.Q40.somatic.assembled.HQfiltered.NEW.csv.annot
${SCRIPT_DIR}/MergeAssembledCallsets.pl -c -f $ID.ctx.Q40.somatic.assembled.HQfiltered.NEW.fasta -d 500 -h merge_ctx.NEW.index 1> $ID.ctx.Q40.somatic.assembled.HQfiltered.NEW.csv 2> $ID.ctx.Q40.somatic.assembled.HQfiltered_out.NEW.csv
${SCRIPT_DIR}/BreakAnnot.pl -A $ID.ctx.Q40.somatic.assembled.HQfiltered.NEW.csv > $ID.ctx.Q40.somatic.assembled.HQfiltered.NEW.csv.annot
cd ..
