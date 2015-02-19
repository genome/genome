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

@test "basic test: genome-env-next -D" {
    $BATS_TEST_DIRNAME/../bin/genome-env-next -D "$BATS_TEST_DIRNAME/testexec/D.sh"
}

@test "basic test: genome-env-next -U" {
    export PERL5LIB="$WORKSPACE/ur/lib:$PERL5LIB"
    $BATS_TEST_DIRNAME/../bin/genome-env-next -U "$BATS_TEST_DIRNAME/testexec/U.sh"
}

@test "basic test: genome-env-next -u" {
    $BATS_TEST_DIRNAME/../bin/genome-env-next -u $WORKSPACE/ur "$BATS_TEST_DIRNAME/testexec/U.sh"
}

@test "basic test: genome-env-next -W" {
    export PERL5LIB="$WORKSPACE/workflow/lib:$PERL5LIB"
    $BATS_TEST_DIRNAME/../bin/genome-env-next -W "$BATS_TEST_DIRNAME/testexec/W.sh"
}

@test "basic test: genome-env-next -w" {
    $BATS_TEST_DIRNAME/../bin/genome-env-next -w $WORKSPACE/workflow "$BATS_TEST_DIRNAME/testexec/W.sh"
}
