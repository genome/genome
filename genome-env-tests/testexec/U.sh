#!/bin/bash

if test -z "$WORKSPACE"
then
    echo "WORKSPACE is not set" >&2
    exit 1
fi

source "$BATS_TEST_DIRNAME/test_helper.bash"

for M in sqitch/genome ur ; do
    submodule_is_clean $M
done
submodule_is_not_initialized ur
module_not_loaded_from_submodule UR
apipe_test_db_is_not_used
