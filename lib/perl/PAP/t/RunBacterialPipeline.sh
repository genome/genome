#! /bin/bash
## This is specific to running bacterial pipeline on a specified test genome from Sasi's Repo
set -e

NOW="$(date +"%Y%m%d")$$"

INPUT_VAR=$1

if [ ${INPUT_VAR} == "klebsiella" ]
then
	# Generate LOCUS_NAME - this will be the argument for the perl script that generates the config file
	LOCUS_NAME="KLEBSIELLA"$NOW"MTST"
	TEMPLATE="/gscuser/ssurulir/test-genome/KLEBSIELLA/KLEBSIELLA_template_config_date"
	CONFIG_DIR="/gscuser/ssurulir/test-genome/KLEBSIELLA/"
	CORE_GENE_CHECK=" "
elif [ ${INPUT_VAR} == "bifidobacterium" ]
then
	# Generate LOCUS_NAME - this will be the argument for the perl script that generates the config file
	LOCUS_NAME="BIFBER"$NOW"MTST"
	TEMPLATE="/gscuser/ssurulir/test-genome/BIFBRE2011TST/BIFBRE2011TST_template_config_date"
	CONFIG_DIR="/gscuser/ssurulir/test-genome/BIFBRE2011TST/"
	CORE_GENE_CHECK="--skip-core-check"
else
	echo "Jibberish"
	exit
fi

# We will open file handle start write to fh
exec 6>&1
README_FNAME="/gscmnt/278/analysis/HGMI/config_files/READMEs/batch_runs/README_ss_"$LOCUS_NAME"_"$NOW
exec > $README_FNAME

echo $(date)
echo "** Generating config-file for this batch run"

# Generate config file - this config file gets written to /gscuser/ssurulir/test-genome/BIFBRE2011TST
perl /gscuser/bdericks/repos/genome/pap/lib/perl/PAP/t/generate_config.pl --locus-name $LOCUS_NAME --template $TEMPLATE

cd /gscuser/ssurulir/workspace/SandBox
# Generated file ought to be at..error checking is done before bsub'ing the command
CONFIG_FNAME=$CONFIG_DIR$LOCUS_NAME"_config_"$NOW

LOCUS_ID=`grep "locus_id" ${CONFIG_FNAME} |cut -d" " -f2 `
LOCUS_TAG=`grep "locus_tag" ${CONFIG_FNAME} |cut -d" " -f2 `
ASSEMBLY_NAME=`grep "assembly_name:" ${CONFIG_FNAME} |cut -d" " -f2 `
DIR_NAME=`grep "org_dirname:" ${CONFIG_FNAME} |cut -d" " -f2 `
ANNO_DIR="/gscmnt/278/analysis/HGMI/"$DIR_NAME"/"$ASSEMBLY_NAME"/Version_1.0/Genbank_submission/Version_1.0/Annotated_submission"

## Output and err file for this run
SCREENOUT="/gscmnt/278/analysis/HGMI/config_files/Screenoutput/batch_runs/"$LOCUS_NAME"_Screenoutput_"$NOW".bsub.out"
SCREENERR="/gscmnt/278/analysis/HGMI/config_files/Screenoutput/batch_runs/"$LOCUS_NAME"_Screenoutput_"$NOW".bsub.err"

## bsub command
BSUB="bsub -o $SCREENOUT -e $SCREENERR -q long -n 2 -R 'span[hosts=1] rusage[mem=4098]' -N -u ssurulir@genome.wustl.edu,vbhonagi@genome.wustl.edu gmt hgmi hap --config $CONFIG_FNAME $CORE_GENE_CHECK"

if [ -f $CONFIG_FNAME ];
then
	echo "** Config file: "$CONFIG_FNAME
	cat $CONFIG_FNAME
	echo -e
	echo "** Running bacterial pipeline from Sasi's repo"
	cd /gscuser/ssurulir/workspace/genome-stable/lib/perl

	echo $(pwd) "> "$BSUB
	eval ${BSUB}

	JOBID=`grep "is submitted to queue" ${README_FNAME} |cut -d'<' -f2|cut -d'>' -f1 `
	#echo ${JOBID}
	#echo $(date)
	echo "Sleeping for 1hr..z z z"
	sleep 3600
	#JOBSTATUS=`bjobs -a ${JOBID}|cut -d' ' -f3 `
	JOBSTATUS=`perl /gscuser/ssurulir/workspace/SandBox/lsfHistory.pl --job-id ${JOBID}`
	while [ $JOBSTATUS != "DONE" ]; do
		echo "Job not complete yet...sleeping for 1 more hour..z z z"
		sleep 3600
		#JOBSTATUS=`bjobs -a ${JOBID}|cut -d' ' -f3 `
		JOBSTATUS=`perl /gscuser/ssurulir/workspace/SandBox/lsfHistory.pl --job-id ${JOBID}`
	done

	cd $ANNO_DIR
	DATFILE=`ls | grep "sqlite-$LOCUS_NAME-*.dat" `
	if [ -f $DATFILE ]
	then 
		echo "Found $DATFILE"
		echo "Now running submission scrips"
		if [ -d /gscmnt/239/info/submissions/assembly_submissions/human_microbes/Klebsiella_spMS92-3/draft/ ]; then
			/gscmnt/277/analysis/personal_dirs/vbhonagi/HGMI_submission_scripts/run_parse_scripts_for_FNLsubmission_svn_SeqAnnot_4Testgenome.sh /gscmnt/278/analysis/HGMI/Acedb/Development ${LOCUS_ID}_  Klebsiella_sp $LOCUS_TAG 1 /gscmnt/239/info/submissions/assembly_submissions/human_microbes/Klebsiella_spMS92-3/draft/ -nornaScr HMP:1111  >&submn_stdout
		else
			echo "Assembly dir /gscmnt/239/info/submissions/assembly_submissions/human_microbes/Klebsiella_spMS92-3/draft/ does not exist"
		fi
		exec 1>&6 6>&-
	else
		echo "Missing dat file from BER naming. Something doesn't look right. Perhaps it's Sunday/Tuesday/Thursday/Saturday ;)"
		exec 1>&6 6>&-
		exit
	fi
else 
	echo "** Config file: "$(pwd)"/"$CONFIG_FNAME " not found."
	exec 1>&6 6>&-
	exit
fi
