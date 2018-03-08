TEST_CHR=20
TEST_START=44860117
TEST_END=44863487
NORMAL_SAMPLE=H_NJ-HCC1395-HCC1395_BL
NORMAL_BAM=/gscmnt/gc9022/info/build_merged_alignments/merged-alignment-blade13-4-5.gsc.wustl.edu-jwalker-25191-5e2355d1bb8c4ff78279ddb3947fc0d3/5e2355d1bb8c4ff78279ddb3947fc0d3.bam
TUMOR_SAMPLE=H_NJ-HCC1395-HCC1395
TUMOR_BAM=/gscmnt/gc8002/info/build_merged_alignments/merged-alignment-blade13-2-11.gsc.wustl.edu-jwalker-27691-f1e64d1802b0493095d70b837a467f47/f1e64d1802b0493095d70b837a467f47.bam
REF_FASTA=/gscmnt/ams1102/info/model_data/2869585698/build106942997/all_sequences.fa

samtools1.2 view ${NORMAL_BAM} ${TEST_CHR}:${TEST_START}-${TEST_END} | cut -f 1 | sort | uniq > ${NORMAL_SAMPLE}.${TEST_CHR}_${TEST_START}-${TEST_END}.txt

gmt picard filter-sam-reads --use-version=1.123 --input=${NORMAL_BAM}  --read-list-file=${NORMAL_SAMPLE}.${TEST_CHR}_${TEST_START}-${TEST_END}.txt --filter=includeReadList --output-file=${NORMAL_SAMPLE}.${TEST_CHR}_${TEST_START}-${TEST_END}.bam --nowrite-reads-files

samtools1.2 view ${TUMOR_BAM} ${TEST_CHR}:${TEST_START}-${TEST_END} | cut -f 1 | sort | uniq > ${TUMOR_SAMPLE}.${TEST_CHR}_${TEST_START}-${TEST_END}.txt

gmt picard filter-sam-reads --use-version=1.123 --input=${TUMOR_BAM} --read-list-file=${TUMOR_SAMPLE}.${TEST_CHR}_${TEST_START}-${TEST_END}.txt --filter=includeReadList --output-file=${TUMOR_SAMPLE}.${TEST_CHR}_${TEST_START}-${TEST_END}.bam --nowrite-reads-files

samtools1.2 faidx ${REF_FASTA} ${TEST_CHR}:${TEST_START}-${TEST_END} > ${TEST_CHR}_${TEST_START}-${TEST_END}.fa

samtools1.2 faidx ${TEST_CHR}_${TEST_START}-${TEST_END}.fa

/gscmnt/gc2560/core/speedseq_freeze/v3/speedseq/bin/bwa index ${TEST_CHR}_${TEST_START}-${TEST_END}.fa

/gscmnt/gc2560/core/speedseq_freeze/v3/speedseq/bin/speedseq realign -o ${TUMOR_SAMPLE}.${TEST_CHR}_${TEST_START}-${TEST_END}.realigned ${TEST_CHR}_${TEST_START}-${TEST_END}.fa ${TUMOR_SAMPLE}.${TEST_CHR}_${TEST_START}-${TEST_END}.bam

/gscmnt/gc2560/core/speedseq_freeze/v3/speedseq/bin/speedseq realign -o ${NORMAL_SAMPLE}.${TEST_CHR}_${TEST_START}-${TEST_END}.realigned ${TEST_CHR}_${TEST_START}-${TEST_END}.fa ${NORMAL_SAMPLE}.${TEST_CHR}_${TEST_START}-${TEST_END}.bam
