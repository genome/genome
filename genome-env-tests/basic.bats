#!/usr/bin/env bats

load test_helper
export BATS_TEST_DIRNAME

function setup {
    init_workspace
    export GE_NO_REDIRECT=1
    cd $WORKSPACE/genome
}

function teardown {
    rm_workspace
}

@test "basic test: genome-env" {
    $BATS_TEST_DIRNAME/../bin/genome-env "$BATS_TEST_DIRNAME/testexec/default.sh"
}

@test "basic test: genome-env -D" {
    $BATS_TEST_DIRNAME/../bin/genome-env -D "$BATS_TEST_DIRNAME/testexec/D.sh"
}

@test "basic test: genome-env -U" {
    export PERL5LIB="$WORKSPACE/ur/lib:$PERL5LIB"
    $BATS_TEST_DIRNAME/../bin/genome-env -D -U "$BATS_TEST_DIRNAME/testexec/U.sh"
}

@test "basic test: genome-env -u" {
    $BATS_TEST_DIRNAME/../bin/genome-env -D -u $WORKSPACE/ur "$BATS_TEST_DIRNAME/testexec/U.sh"
}

@test "basic test: genome-env -M" {
    $BATS_TEST_DIRNAME/../bin/genome-env -D -M "$BATS_TEST_DIRNAME/testexec/M.sh"
}

@test "basic test: genome-env -m" {
    $BATS_TEST_DIRNAME/../bin/genome-env -m $WORKSPACE/genome-sqitch "$BATS_TEST_DIRNAME/testexec/migration.sh"
}
