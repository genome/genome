#!/usr/bin/tcsh


set srr_ids = $argv[1]
set sra_samples = $argv[2]
set picard_dir = $argv[3]
set tmp_dir = $argv[4]
set ascp_user = $argv[5]
set ascp_pw = $argv[6]
set pwd = `pwd`

echo "Checking for appropriate scripts..."
set fastq_dump = `which fastq-dump | awk '{print NF}'`
set trimBWA = `which trimBWAstyle.usingBam.pl | awk '{print NF}'`
set samtools = `which samtools | awk '{print NF}'`
setenv ASPERA_SCP_PASS $ascp_pw

if ( $trimBWA > 1 ) then
	echo "Could not find trimBWAstyle.usingBAM.pl.  Please make sure this is in your path and try again."
	exit
endif

if ( $fastq_dump > 1 ) then 
        echo "Could not find fastq-dump.  Please make sure this is in your path and try again."
        exit
endif

if ( $samtools > 1 ) then 
        echo "Could not find samtools.  Please make sure this is in your path and try again."
        exit
endif
echo ""

foreach sample ( `grep -f $srr_ids $sra_samples | awk '{print $2}' | sort | uniq` )

	#Number of runs in this sample
	set runs = `grep $sample $sra_samples | awk '{print $1}' | sort | uniq | wc | awk '{print $1}'`
	
       	echo "`date` ${sample}: Processing $runs run(s)..."

	if ( -d $sample ) then
	    echo "`date` ${sample}: Processing already started..."
	else
	    mkdir $sample
	endif

	set sn = "sn"
	set ri = "ri"
	set sum = 0
	set failed = 0
	
	foreach srr_id ( `grep $sample $sra_samples | awk '{print $1}' | sort | uniq` )

		if ($failed) continue

		#Pull meta data if it does not exist
		if (! -e $sample/$srr_id.xml) then
		    echo "\t`date` ${srr_id}: Pulling meta data..."		
		    wget 'http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?run='$srr_id'&retmode=xml' -O $sample/$srr_id.xml >& /dev/null
		else
		    echo "\t`date` ${srr_id}: Meta data exists."
		endif

		#Check to see if SRA has been converted to fastq
		if (-e $sample/$srr_id.fastq-dump) then		
			set dumped = `grep Success $sample/$srr_id.fastq-dump | wc | awk '{print $1}'`
		else
			set dumped = 0
		endif

		#If fastqs have not been created, run fastq-dump
		if ( $dumped > 0 && -e ${sample}/${srr_id}_1.masked && -e ${sample}/${srr_id}_1.masked) then
			#Fastq dump was sucessful
			set reads = `grep Success $sample/${srr_id}.fastq-dump | awk '{print $7}'`
                        @ sum += `echo $reads`
                        set masked = `cat $sample/${srr_id}*masked | awk '{print $1}' | perl -ne '$sum += $_; if (eof) {print $sum;}'`
                        set mbases = `cat $sample/${srr_id}*masked | awk '{print $2}' | perl -ne '$sum += $_; if (eof) {print $sum;}'`
                        echo "\t`date` ${srr_id}: $masked total reads masked $mbases total bases masked."
			echo "\t`date` ${srr_id}: $reads pairs already dumped successfully."
		else
			#Check that data has been downloaded
			if (! -e ${srr_id}) then
				echo "\tERROR: No SRA data for ${srr_id}"
				set failed = 1
				continue
				#exit
			endif

			#Check SRA files for download integrity
			echo "\t`date` ${srr_id}: Checking download integrity... "
			set orig_md5 = `cat $srr_id/col/READ/md5 | awk '$2 ~ /data/ {print $1}'`
			set new_md5 = `md5sum $srr_id/col/READ/data | awk '{print $1}'`
			if ($orig_md5 != $new_md5) then
				echo "\tERROR ${srr_id}: md5 does not check out.  Please remove this directory and retry your download."
				set failed = 1
				continue
				#exit
			endif

			#Dump fastq files from SRA
			echo "\t`date` ${srr_id}: Pulling Fastqs from SRA data..."
			fastq-dump -E -A $srr_id -D $srr_id/ -DB '@$sn/$ri' -DQ '+$sn/$ri' -O $sample/ >& $sample/$srr_id.fastq-dump
     
			#Check that reads were dumped correctly
			set reads = `grep Success $sample/${srr_id}.fastq-dump | awk '{print $7}'`
			if ($reads < 1) then
				echo "\tERROR ${srr_id}: fastq-dump failed"
				set failed = 1
				continue
				#exit
			else
				echo "\t`date` ${srr_id}: $reads pairs dumped successfully."
				@ sum += `echo $reads`
				if (! -e $sample/${srr_id}_1.fastq || ! -e $sample/${srr_id}_2.fastq) then
					echo "\tERROR ${srr_id}: Read one and read two files were not downloaded properly.  Check that NCBI has both ends for this run."
					set failed = 1
					continue
					#exit
				endif
			endif
			
			#Count number of human masked reads
			echo "\t`date` ${srr_id}: Counting masked reads..."
 			perl -ne 'chomp; $count = tr/N/n/; if ($count == length $_) { $masked++; $b += $count;} if (eof){ print "$masked\t$b\n";}' $sample/${srr_id}_1.fastq > $sample/${srr_id}_1.masked
 			perl -ne 'chomp; $count = tr/N/n/; if ($count == length $_) { $masked++; $b += $count;} if (eof){ print "$masked\t$b\n";}' $sample/${srr_id}_2.fastq > $sample/${srr_id}_2.masked
                        set masked = `cat $sample/${srr_id}*masked | awk '{print $1}' | perl -ne '$sum += $_; if (eof) {print $sum;}'`
                        set mbases = `cat $sample/${srr_id}*masked | awk '{print $2}' | perl -ne '$sum += $_; if (eof) {print $sum;}'`
			echo "\t`date` ${srr_id}: $masked total reads masked $mbases total bases masked."

			echo "\t`date` ${srr_id}: Appending fastqs to sample file..."
			#Concatenate the fastqs into one sample
			#___Modified to touch files before appending into them ... jmartin 100907
			touch ${sample}/${sample}_1.fastq
			touch ${sample}/${sample}_2.fastq
			cat ${sample}/${srr_id}_1.fastq >> ${sample}/${sample}_1.fastq
			cat ${sample}/${srr_id}_2.fastq >> ${sample}/${sample}_2.fastq

			#Remove the fastqs
			echo "\t`date` ${srr_id}: Removing fastqs..."
			rm ${sample}/${srr_id}_1.fastq
			rm ${sample}/${srr_id}_2.fastq
		endif
	end

	if ($failed) then
		echo "`date` ${sample}: Some runs are not available.  SKIPPING"
		echo ""
		continue
	endif

	#Converting fastqs into a BAM
	echo "`date` ${sample}: Converting Fastq to BAM..."
	if (! -e $sample/$sample.bam) then
		java -jar $picard_dir/FastqToSam.jar F1=$sample/{$sample}_1.fastq F2=$sample/${sample}_2.fastq O=$sample/$sample.bam V=Standard SAMPLE_NAME=$sample TMP_DIR=$tmp_dir >& ${sample}/FastqToSam.out
	else
		echo "\tSKIPPING: $sample/$sample.bam exists."
	endif

	#Checking that BAM is complete
	if (-e $sample/FastqToSam.out) then
		set bamReads = `grep Processed $sample/FastqToSam.out | awk '{print $2}'`
	else 
		set bamReads = 0
	endif

	if ( $sum == $bamReads ) then
		if ( -e $sample/${sample}_1.fastq ) then		
			echo "\tRemoving sample fastqs: $sample/${sample}_1.fastq $sample/${sample}_2.fastq"
			rm $sample/*.fastq
		endif
	else
		echo "\tERROR: BAM file reads ($bamReads) do not match SRA dumped reads ($sum)."
		exit
	endif

	#Removing Duplicates
	echo "`date` ${sample}: Removing duplicates..."

	if (-e ${sample}/$sample.denovo_duplicates_marked.bam ) then
		if (-e  $sample/EstimateLibraryComplexity.out) then
		    if ( `grep "done" $sample/EstimateLibraryComplexity.out | wc | awk '{print $1}'` > 0 ) then
			echo "\tSKIPPING: $sample/$sample.denovo_duplicates_marked.bam exists"
		    else
			echo "\tERROR: EstimateLibraryComplexity did not complete.  Please delete ${sample}/$sample.denovo_duplicates_marked.bam and $sample/EstimateLibraryComplexity.out and try again"
			exit
		    endif
	        else 
		    echo "\tERROR: EstimateLibraryComplexity did not complete.  Please delete ${sample}/$sample.denovo_duplicates_marked.bam and $sample/EstimateLibraryComplexity.out and try again"
		    exit
		endif
	else
		#___Modified to force specific memory allocation, and to turn off java garbage collection timeout ... jmartin 100907
		java -Xmx43g -XX:-UseGCOverheadLimit -jar $picard_dir/EstimateLibraryComplexity.jar I=${sample}/$sample.bam O=${sample}/$sample.denovo_duplicates_marked.bam METRICS_FILE=${sample}/$sample.denovo_duplicates_marked.metrics REMOVE_DUPLICATES=true TMP_DIR=$tmp_dir >& ${sample}/EstimateLibraryComplexity.out
		samtools flagstat $sample/$sample.denovo_duplicates_marked.bam > $sample/$sample.denovo_duplicates_marked.counts
        endif
	
	set nondupReads = `awk '$2 ~ /paired/ {print $1}' $sample/$sample.denovo_duplicates_marked.counts`
	echo "`date` ${sample}: $nondupReads left after duplication removal"
	
	#Trimming low quality and masked reads
	echo "`date` ${sample}: Trimming Q2 bases..."
	if (! -e ${sample}/${sample}.denovo_duplicates_marked.trimmed.1.fastq.bz2 ) then
		trimBWAstyle.usingBam.pl -o 33 -q 3 -f ${sample}/$sample.denovo_duplicates_marked.bam > ${sample}/trimBWAstyle.out
		set qualReads = `grep reads ${sample}/trimBWAstyle.out | awk '{print $3}'`
		set qualBases = `grep bases ${sample}/trimBWAstyle.out | awk '{print $3}'`
		echo "`date` ${sample}: $qualReads reads ($qualBases bases) left after trimming" 
		echo "`date` ${sample}: Compressing trimmed files..."
		bzip2 ${sample}/${sample}.denovo_duplicates_marked.trimmed.*.fastq
	else 
		echo "\tSKIPPING: ${sample}/$sample.duplicates_removed.Q2trimmed.fastq.bz2 exists."
                set qualReads = `grep reads ${sample}/trimBWAstyle.out | awk '{print $3}'`
                set qualBases = `grep bases ${sample}/trimBWAstyle.out | awk '{print $3}'`
                echo "`date` ${sample}: $qualReads reads ($qualBases bases) left after trimming"
	endif

	echo "`date` ${sample}: Creating checksum files..."
	if (! -e ${sample}/${sample}.md5) then
		md5sum $sample/*.bz2 > $sample/$sample.md5
	else
		echo "\tSKIPPING: $sample/$sample.md5 exists"
	endif

	echo "`date` ${sample}: Uploading files to DACC..."
	if (! -e ${sample}/${sample}.ascp) then
#_______________TEMPORARY CHANGE JUST FOR TESTING ... jmartin 100902
		####ascp -QTd $sample/*.md5 $sample/*.xml $sample/*.bz2 $ascp_user@aspera.hmpdacc.org:/WholeMetagenomic/02-ScreenedReads/ProcessedForAssembly/$sample/ > $sample/$sample.ascp
		ascp -QTd -l100M $sample/*.md5 $sample/*.xml $sample/*.bz2 $ascp_user@aspera.hmpdacc.org:/WholeMetagenomic/02-ScreenedReads/ProcessedForAssembly/$sample/ > $sample/$sample.ascp
	else
		set files = `grep files ${sample}/${sample}.ascp | awk '{print $4}'`
		if ( $files > 4 ) then
			echo "\tSKIPPING: $sample transfered."
		else
#_______________________TEMPORARY CHANGE JUST FOR TESTING ... jmartin 100902
	                ####ascp -QTd $sample/*.md5 $sample/*.xml $sample/*.bz2 $ascp_user@aspera.hmpdacc.org:/WholeMetagenomic/02-ScreenedReads/ProcessedForAssembly/$sample/ > $sample/$sample.ascp
			ascp -QTd -l100M $sample/*.md5 $sample/*.xml $sample/*.bz2 $ascp_user@aspera.hmpdacc.org:/WholeMetagenomic/02-ScreenedReads/ProcessedForAssembly/$sample/ > $sample/$sample.ascp
		endif		
	endif
	
	echo "`date` ${sample}: Processing complete."
	echo ""
end
