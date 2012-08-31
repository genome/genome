#!/bin/tcsh

set sddir=$1;
set sdfile=$2;
set libfile=$3;
set bamfile=$4;
set scriptpath=`dirname $0`;

# do with bamfile: dir/tumor.bam,normal.bam..., separated by comma, will create dirctory $sddir/assembly_tumor; $sddir/assembly_tumor_old; $sddir/assembly_normal; $sddir/assembly_normal_old

echo $bamfile $sddir | perl -ane '@a = split(" "); ($b)=($a[0]=~/\/([^\/]+)$/); $b=~s/.bam//g; @c=split(",",$b); foreach $d (@c) {$command="mkdir ".$a[1]."/assembly_".$d."\n";system($command);}'

# get the library for this sample; it should have the normal anyway
set g="NA";
echo $g;
if(-e ${libfile}) then
	set g=`cut -f3,5 ${libfile} | grep normal | cut -f2 | perl -ane 'chomp; @x = split ":"; print "$x[1]\n"' | sort -u | awk '{if($NR==1) a=$0; if($NF>1) $a=$a "," $0;} END{print $a}'`;
	#join(",", @x);'`;
	#$g=`echo $g | perl -ane 'chomp; $_=~s/ /,/g; print $_'`;
	echo $g;
endif
#set g="HAHA";
	
set project=`echo $libfile | perl -ane '($x)=($F[0]=~/\/([^\/]+).cfg$/); print $x'`;
#set project="TEST";

set gmt='`which gmt`';

foreach directory (`ls -d $sddir/assembly_*`)

#	set directory=$sddir/assembly_tumor_old;
#echo $directory $bamfile;

# annotation
	if(! -e ${sdfile}.annot) then
	    ${scriptpath}/BreakAnnot.pl $sdfile > ${sdfile}.annot
	endif

	set bam_complete=`echo $directory $bamfile | perl -ane '($x)=($F[0]=~/assembly_([^\_]+)/); ($y)=($F[1]=~/^(.+)\/[^\/]+$/); print $y . "/". $x.".bam"'`;
#echo $bam_complete;
	set job=`echo $directory | perl -ane '($x)=($F[0]=~/assembly_(\S).+$/); $y=uc($x); print $y'`;
	
	set bam = `echo $directory | perl -ane '($x)=($F[0]=~/assembly_([^\_]+)/); print $x;'`;
#echo $bam;	
	# tumor
	
	set R1='"select[type==LINUX64 && mem>8000]';
	set R2='rusage[mem=8000]"';
	
	if(! `echo $directory | perl -ane 'print $F[0]=~/old/'`) then
		# new
		set new="bsub -q apipe -N -u someone@genome.wustl.edu -M 8000000 -R ${R1} ${R2} -e $directory/$bam.log -J ${project}_${job}_SD_AS 'perl -I ${scriptpath} $gmt sv assembly-validation --bam-files $bam_complete --sv-file ${sdfile} --breakpoint-seq-file $directory/$bam.fasta --cm-aln-file $directory/$bam.cm --min-breakdancer-score 40 --intermediate-read-dir $directory --output-file $directory/$bam.csv --skip-libraries $g --min-size-of-confirm-asm-sv 10'";
		if(-e ${directory}/${bam}.csv) then
			rm -f ${directory}/${bam}.csv;
		endif
		bsub -q apipe -N -u someone@genome.wustl.edu -M 8000000 -R "select[type==LINUX64 && mem>8000] rusage[mem=8000]" -e $directory/$bam.log -J ${project}_${job}_SD_AS "perl -I ${scriptpath} $gmt sv assembly-validation --bam-files $bam_complete --sv-file $sdfile --breakpoint-seq-file $directory/$bam.fasta --cm-aln-file $directory/$bam.cm --min-breakdancer-score 40 --intermediate-read-dir $directory --output-file $directory/$bam.csv --skip-libraries $g --min-size-of-confirm-asm-sv 10";
	endif
	
end	

