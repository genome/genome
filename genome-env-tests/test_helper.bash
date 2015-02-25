function abs_path {
    echo "$( cd -P "$1" && pwd )"
}

function tempdir {
    set -o errexit

    local TEMPDIR;

    if man mktemp | grep -q BSD
    then
        TEMPDIR="$(mktemp -d -t genome-env-testing)"
    fi

    if man mktemp | grep -q GNU
    then
        TEMPDIR="$(mktemp --tmpdir -d genome-env-testing.XXXXXXXX)"
    fi

    echo "$TEMPDIR"
}

function cache_repo {
    local NAME="$1"
    local URL="$2"
    local CLONE_BASE_DIR="$3"

    local CACHE_BASE_DIR="/tmp/genome-env-cache"
    local CACHE_DIR="$CACHE_BASE_DIR/$NAME.git"
    local CLONE_DIR="$CLONE_BASE_DIR/$NAME"

    if test -d $CACHE_DIR
    then
        git --git-dir $CACHE_DIR fetch > /dev/null 2>&1
    else
        git clone --bare $URL $CACHE_DIR > /dev/null 2>&1
    fi

    git clone --reference $CACHE_DIR $URL $CLONE_DIR > /dev/null 2>&1
}

function init_workspace {
    export WORKSPACE="$(tempdir)"
    cache_repo genome https://github.com/genome/genome.git "$WORKSPACE"
    cache_repo ur https://github.com/genome/UR.git "$WORKSPACE"
    cache_repo workflow https://github.com/genome/tgi-workflow.git "$WORKSPACE"
    cache_repo genome-sqitch https://github.com/genome/genome-sqitch.git "$WORKSPACE"
}

function rm_workspace {
    if test -n "$WORKSPACE" -a -d "$WORKSPACE"
    then
        rm -rf "$WORKSPACE"
    fi
}

function submodule_is_clean {
    if ! git diff --exit-code $1
    then
        echo "submodule should be clean: $1" >&2
        exit 1
    fi
}

function submodule_is_initialized {
    if git submodule status $1 | grep -q '^-'
    then
        echo "submodule should be initialized: $1" >&2
        exit 1
    fi
}

function submodule_is_not_initialized {
    if git submodule status $1 | grep -qv '^-'
    then
        echo "submodule should be not initialized: $1" >&2
        exit 1
    fi
}

function apipe_test_db_is_used {
    if echo "$GENOME_DS_GMSCHEMA_SERVER" | grep -qv 'apipe-test-db'
    then
        echo "GENOME_DS_GMSCHEMA_SERVER should refer to apipe-test-db" >&2
        exit 1
    fi
}

function apipe_test_db_is_not_used {
    if echo "$GENOME_DS_GMSCHEMA_SERVER" | grep -q 'apipe-test-db'
    then
        echo "GENOME_DS_GMSCHEMA_SERVER should not refer to apipe-test-db" >&2
        exit 1
    fi
}

function module_loaded_from_submodule {
    local SUBMODULE="${1,,}"
    if test "$WORKSPACE/genome/$SUBMODULE/lib/$1.pm" != "$(perl -M$1 -e "print \$INC{q($1.pm)}, qq(\\n)")"
    then
        echo "$1 should be loaded from submodule" >&2
        exit 1
    fi
}

function module_not_loaded_from_submodule {
    local SUBMODULE="${1,,}"
    if test "$WORKSPACE/genome/$SUBMODULE/lib/$1.pm" == "$(perl -M$1 -e "print \$INC{q($1.pm)}, qq(\\n)")"
    then
        echo "$1 should not be loaded from submodule" >&2
        exit 1
    fi
}
