samtools1.2 view /gscmnt/gc9022/info/build_merged_alignments/merged-alignment-blade13-4-5.gsc.wustl.edu-jwalker-25191-5e2355d1bb8c4ff78279ddb3947fc0d3/5e2355d1bb8c4ff78279ddb3947fc0d3.bam 20:42220611-42542245 | cut -f 1 | sort | uniq > /tmp/H_NJ-HCC1395-HCC1395_BL.20_42220611-42542245.txt

gmt picard filter-sam-reads --use-version=1.123 --input=/gscmnt/gc9022/info/build_merged_alignments/merged-alignment-blade13-4-5.gsc.wustl.edu-jwalker-25191-5e2355d1bb8c4ff78279ddb3947fc0d3/5e2355d1bb8c4ff78279ddb3947fc0d3.bam  --read-list-file=/tmp/H_NJ-HCC1395-HCC1395_BL.20_42220611-42542245.txt --filter=includeReadList --output-file=H_NJ-HCC1395-HCC1395_BL.20_42220611-42542245.bam --nowrite-reads-files

samtools1.2 view /gscmnt/gc8002/info/build_merged_alignments/merged-alignment-blade13-2-11.gsc.wustl.edu-jwalker-27691-f1e64d1802b0493095d70b837a467f47/f1e64d1802b0493095d70b837a467f47.bam 20:42220611-42542245 | cut -f 1 | sort | uniq > /tmp/H_NJ-HCC1395-HCC1395.20_42220611-42542245.txt

gmt picard filter-sam-reads --use-version=1.123 --input=/gscmnt/gc8002/info/build_merged_alignments/merged-alignment-blade13-2-11.gsc.wustl.edu-jwalker-27691-f1e64d1802b0493095d70b837a467f47/f1e64d1802b0493095d70b837a467f47.bam --read-list-file=/tmp/H_NJ-HCC1395-HCC1395.20_42220611-42542245.txt --filter=includeReadList --output-file=H_NJ-HCC1395-HCC1395.20_42220611-42542245.bam --nowrite-reads-files

samtools1.2 faidx /gscmnt/ams1102/info/model_data/2869585698/build106942997/all_sequences.fa 20:42220611-42542245 > 20_42220611-42542245.fa

samtools1.2 faidx 20_42220611-42542245.fa

/gscmnt/sata849/info/speedseq_freeze/v3/speedseq/bin/bwa index 20_42220611-42542245.fa

/gscmnt/sata849/info/speedseq_freeze/v3/speedseq/bin/speedseq realign -o H_NJ-HCC1395-HCC1395.20_42220611-42542245.realigned 20_42220611-42542245.fa H_NJ-HCC1395-HCC1395.20_42220611-42542245.bam

/gscmnt/sata849/info/speedseq_freeze/v3/speedseq/bin/speedseq realign -o H_NJ-HCC1395-HCC1395_BL.20_42220611-42542245.realigned 20_42220611-42542245.fa H_NJ-HCC1395-HCC1395_BL.20_42220611-42542245.bam
