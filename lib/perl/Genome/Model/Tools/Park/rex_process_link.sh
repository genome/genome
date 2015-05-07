#! /bin/bash

PARK_BETA_ROOT="$(genome config get park_beta_root)"
source $PARK_BETA_ROOT/workon_park_user.sh
rex process link $@
