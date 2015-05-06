#! /bin/bash

eval "$(genome config get --format=bash park_beta_root)"
source $GENOME_PARK_BETA_ROOT/workon_park_user.sh
rex process start $@
