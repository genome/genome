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

submodule_is_not_initialized workflow

if test "$WORKSPACE/genome/ur/lib/UR.pm" != "$(perl -MUR -e 'print $INC{q(UR.pm)}, qq(\n)')"
then
    echo "UR should be loaded from submodule" >&2
    exit 1
fi

if test "$WORKSPACE/genome/workflow/lib/Workflow.pm" == "$(perl -MWorkflow -e 'print $INC{q(Workflow.pm)}, qq(\n)')"
then
    echo "Workflow should not be loaded from submodule" >&2
    exit 1
fi

if echo "$GENOME_DS_GMSCHEMA_SERVER" | grep -q 'apipe-test-db'
then
    echo "GENOME_DS_GMSCHEMA_SERVER should not refer to apipe-test-db" >&2
    exit 1
fi
