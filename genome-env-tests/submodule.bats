#!/usr/bin/env bats

load test_helper

function setup {
    init_workspace
    export GE_NO_REDIRECT=1
    cd $WORKSPACE/genome
}

function teardown {
    rm_workspace
}

@test "unclean submodule fails" {
    git submodule update --init ur
    echo >> ur/lib/UR.pm
    run $BATS_TEST_DIRNAME/../bin/genome-env-next
    test "$status" -ne 0
}

@test "submodule with commits fails" {
    git submodule update --init ur
    (
        cd ur
        echo >> lib/UR.pm
        git commit -m whitespace lib/UR.pm
    )
    run $BATS_TEST_DIRNAME/../bin/genome-env-next
    test "$status" -ne 0
}
