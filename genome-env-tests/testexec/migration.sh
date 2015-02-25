#!/bin/bash

type -p genome-env-next

if test -z "$WORKSPACE"
then
    echo "WORKSPACE is not set" >&2
    exit 1
fi

for M in sqitch/genome ur workflow ; do
    submodule_is_clean $M
done

submodule_is_not_initialized sqitch/genome

if echo "$GENOME_DS_GMSCHEMA_SERVER" | grep -qv 'apipe-test-db'
then
    echo "GENOME_DS_GMSCHEMA_SERVER should refer to apipe-test-db" >&2
    exit 1
fi
