#!/gsc/bin/sh

ID=$1;
DIR=$2;
SCRIPT_DIR="$( cd "$( dirname "$0" )" && pwd )"
#scripts/fixCTXswapXY_chr.pl $ID
cd ${DIR};
if [ -s merge.index ]; then
    rm -f merge.index;
fi
#if [ -s merge_ctx.index ]; then
#    rm -f merge_ctx.index;
#fi

for file in tumor normal; do
   echo ${file} 
    if [ -s assembly_${file}/${file}.csv ]; then
        if [ -s assembly_${file}/${file}.fasta ]; then
            cp assembly_${file}/${file}.csv $ID.${file}.csv;
            cp assembly_${file}/${file}.fasta $ID.${file}.fasta;
        fi
    fi
    #    for i in `seq 1 22` X Y; do 
    #if [ -s $ID.chr$i.Q40.somatic.assembled.$s.csv ]; then
    if [ -s $ID.${file}.csv ]; then
        echo ${file} $ID.${file}.csv $ID.${file}.fasta >> merge.index; 
    fi
    #    done
done

MAC_SCRIPT="${SCRIPT_DIR}/MergeAssembledCallsets.pl"
if [ -x "$MAC_SCRIPT" ]; then
    $MAC_SCRIPT -c -f $ID.HQfiltered.fasta -d 200 -h merge.index 1> $ID.HQfiltered.csv 2>$ID.HQfiltered_out.csv
else
    echo "${MAC_SCRIPT}: No such file or directory" 1>&2
    exit 1
fi

${SCRIPT_DIR}/BreakAnnot.pl -A $ID.HQfiltered.csv > $ID.HQfiltered.csv.annot
cd ..
