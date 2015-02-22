#!/bin/bash

type -p genome-env-next

if test -z "$WORKSPACE"
then
    echo "WORKSPACE is not set" >&2
    exit 1
fi

for M in sqitch/genome ur workflow ; do
    if ! git diff --exit-code $M
    then
        echo "submodule should be clean: $M" >&2
        exit 1
    fi
done

if git submodule status sqitch/genome | grep -qv '^-'
then
    echo "submodule should not be initialized: sqitch/genome" >&2
    exit 1
fi

if echo "$GENOME_DS_GMSCHEMA_SERVER" | grep -qv 'apipe-test-db'
then
    echo "GENOME_DS_GMSCHEMA_SERVER should refer to apipe-test-db" >&2
    exit 1
fi
