#!/usr/bin/tcsh


set srr_ids = $argv[1]      #This is a FILE containing a list of 'subset_names' for vervet (subset_name == lane/index specification)
set sra_samples = $argv[2]  #This is a FILE containing a 2-column list showing '<subset_name>	 <full_sample_name>'
set picard_dir = $argv[3]   #This should point to the version of picard to use...due to some custom modifications made for us by the author, this needs to be: $ENV{GENOME_SW_LEGACY_JAVA}/samtools/picard-tools-1.27
set tmp_dir = $argv[4]      #This points to a temporary working directory created by: Genome::Sys->create_temp_directory
set pwd = `pwd`

echo "Checking for appropriate scripts..."
set trimBWA = `which trimBWAstyle.usingBam.pl | awk '{print NF}'`
set samtools = `which samtools | awk '{print NF}'`

if ( $trimBWA > 1 ) then
	echo "Could not find trimBWAstyle.usingBAM.pl.  Please make sure this is in your path and try again."
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






#### This part replaced other code ... this is ONLY for VERVET project ... jmartin 111030
	#Locate BAM files for current sample
	#_NOTE: In the 'sra_samples' file, I have the subset_name (Erica's ID) & full_sample_name (eg. WFAA-1080-0104013252) stored.......I should be able to get 'bam_path' using full_sample_name
        #       BUT, in this case, I will have TWO (or more?) bam files for each full_sample_name, since for normal HMP that merging happened during the code I removed for this version.   So
	#       I need to MERGE the multiple bam files per sample name here
	#
	# ----> This merges all the individual lane/index (i.e. subset_name) bam files into a SINGLE bam file at the location ${sample}/$sample.bam
	#
	# NOTE: Need to be sure to always ONLY use subset names that are requested (since some samples may have additional subset names)
	#
        set merged_bam_file = ${sample}/$sample.bam
        set cmd_merge = "samtools merge $merged_bam_file "
	set component_bam_count = 0
	foreach line ("`genome instrument-data list imported -filter sample_name="$sample" -show bam_path,description`")
	        if ("$line" =~ "BAM_PATH") then
		    set var = 1
		else if ("$line" =~ "---*") then
		    set var = 1
		else
		    set bam_file_path = `echo $line | awk '{print $1;}'`
		    set current_subset_name = `echo $line | awk 'BEGIN{FS="(";}{print $2;}' | awk '{print $1;}' | awk 'BEGIN{FS=":";}{print $2;}'`
		    set want_this_subset = 0
		    foreach subset (`cat $srr_ids`)
			if ("$subset" == "$current_subset_name") then
			    set want_this_subset = 1
			endif
		    end
		    if ($want_this_subset > 0) then
			set cmd_merge = "$cmd_merge $bam_file_path "
			@ component_bam_count = $component_bam_count + 1
			set last_component_bam = $bam_file_path
		    endif
		endif
	end

        #Run samtools merge for the current sample if there are 2+ component bam files, if just 1 component bam file just make a copy ... jmartin 111030
	if ($component_bam_count > 1) then
	    echo "`date` ${sample}: Merging bam files using this command =>$cmd_merge<="
	    $cmd_merge
	else if ($component_bam_count == 1) then
	    echo "`date` ${sample}: Only 1 component bam file found ... no need to merge"
	    cp $last_component_bam $merged_bam_file
	else
	    echo "`date` ${sample}: Did not find any component bam files for sample =>$sample<="
	    exit
	endif

	#Make sure the merged bam file is sitting in the right location
	if (! -e $merged_bam_file) then
	        echo "`date` ${sample}: Failed to build merged bam file..."
		exit
        else
		echo "`date` ${sample}: Merging bam files complete and file exists"
	endif

	#Now check to make sure the merged bam file contains the expected number of vervet-host-free reads
        set input_read_count = 0
	foreach line ("`genome instrument-data list imported -filter sample_name="$sample" -show description`")
	        if ("$line" =~ "DESCRIPTION") then
		    set var = 1
		else if ("$line" =~ "---*") then
		    set var = 1
		else
		    set current_subset_name = `echo $line | awk 'BEGIN{FS="(";}{print $2;}' | awk '{print $1;}' | awk 'BEGIN{FS=":";}{print $2;}'`
		    set want_this_subset = 0
		    foreach subset (`cat $srr_ids`)
			if ("$subset" == "$current_subset_name") then
			    set want_this_subset = 1
			endif
		    end
		    if ($want_this_subset > 0) then
			set reads_in_this_subset = `echo $line | awk '{print $NF;}' | awk 'BEGIN{FS=":";}{print $2;}' | sed 's/)//'`
			@ input_read_count = $input_read_count + $reads_in_this_subset
		    endif
		endif
	end
	echo "`date` ${sample}: counted merge input reads=>$input_read_count<="

	set tmp_sam = ${sample}/$sample.tmp_sam
	samtools view $merged_bam_file > $tmp_sam
	set reads_in_merged_bam = `wc $tmp_sam | awk '{print $1;}'`
	echo "`date` ${sample}: counted reads in merge output=>$reads_in_merged_bam<="
	rm -f $tmp_sam

	#Here I test to make sure that the merged bam is not truncated...because some bam files were built using older sam, its possible that there may be an error message stating "[bam_header_read] EOF marker is absent. The input is probably truncated." ... this step is meant to confirm that is just an artifact of using different versions of samtools, or an artifact of the fact that these bams are NOT from alignments, but are just a convenient way to store sequence (thus no header is present)
	if ($input_read_count != $reads_in_merged_bam) then
	        echo "`date` ${sample}: merged bam does not have same number reads as the input bam files merged when using this command=>$cmd_merge<=  ...  input bam read count=>$input_read_count<= merged bam read count=>$reads_in_merged_bam<="
		exit
	endif
	echo "`date` ${sample}: Merge of bam files successful"

	############### NEW ON 111118 ... jmartin
	#Here I do a READ-NAME BASED sort the merged bam file
	set cmd_nameBased_sort = "java -jar $picard_dir/SortSam.jar I=${sample}/$sample.bam O=${sample}/$sample.SSsort.bam TMP_DIR=$tmp_dir SORT_ORDER=queryname"
	$cmd_nameBased_sort
	if (-e ${sample}/$sample.SSsort.bam) then
	    echo "`date` ${sample}: Name-based sort of bam file (via SortSam.jar) successful"
	else
	    echo "`date` ${sample}: ERROR: failed to do name-based sort (via SortSam.jar) on merged bam file"
	    exit
	endif
	############### NEW ON 111118 ... jmartin

#### This part replaced other code ... this is ONLY for VERVET project







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
####		java -Xmx43g -XX:-UseGCOverheadLimit -jar $picard_dir/EstimateLibraryComplexity.jar I=${sample}/$sample.bam O=${sample}/$sample.denovo_duplicates_marked.bam METRICS_FILE=${sample}/$sample.denovo_duplicates_marked.metrics REMOVE_DUPLICATES=true TMP_DIR=$tmp_dir >& ${sample}/EstimateLibraryComplexity.out
		#___Modified to use 60Gb memory allocation to see if this helped eliminate the truncated bam files being produced ... jmartin 111108
		####java -Xmx60g -XX:-UseGCOverheadLimit -jar $picard_dir/EstimateLibraryComplexity.jar I=${sample}/$sample.bam O=${sample}/$sample.denovo_duplicates_marked.bam METRICS_FILE=${sample}/$sample.denovo_duplicates_marked.metrics REMOVE_DUPLICATES=true TMP_DIR=$tmp_dir >& ${sample}/EstimateLibraryComplexity.out
		#___Modified to use SORTED bam file ... jmartin 111114
		java -Xmx60g -XX:-UseGCOverheadLimit -jar $picard_dir/EstimateLibraryComplexity.jar I=${sample}/$sample.SSsort.bam O=${sample}/$sample.denovo_duplicates_marked.bam METRICS_FILE=${sample}/$sample.denovo_duplicates_marked.metrics REMOVE_DUPLICATES=true TMP_DIR=$tmp_dir >& ${sample}/EstimateLibraryComplexity.out
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
	
	echo "`date` ${sample}: Processing complete."
end
