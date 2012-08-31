#!/bin/bash

set -o errexit   # exit if any command fails
set -o pipefail  # fail if any command in a pipe fails
set -o nounset   # fail if an env var is used but unset


# updating solr-dev
# TODO: disk_group, disk_volume, flowcell



bsub -u jlolofie@genome.wustl.edu -q apipe GENOME_DEV_MODE=1 perl /gsc/scripts/bin/gcsearch_load_from_db --add user,individual,sample,population_group,taxon --lock 2
bsub -u jlolofie@genome.wustl.edu -q apipe GENOME_DEV_MODE=1 perl /gsc/scripts/bin/gcsearch_load_from_db --add processing_profile,model_group,library,work_order --lock 3



for i in `seq 0 20`; do
echo "loading model chunk number $i (21 total)"
GENOME_DEV_MODE=1 perl /gsc/scripts/bin/gcsearch_load_from_db --add model --lock 1 --chunk $i
done



for i in `seq 0 20`; do
echo "loading flowcell chunk number $i (21 total)"
GENOME_DEV_MODE=1 perl /gsc/scripts/bin/gcsearch_load_from_db --add flowcell --lock 4 --chunk $i
done





