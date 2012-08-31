export LESS=-eiMXR
GIT_PS1_SHOWDIRTYSTATE=1

last_perl_tests_status() {
    local IS_GENOME_GIT=$(git describe --tags 2> /dev/null | grep -P '^genome-' | wc -l)
    if (( ${IS_GENOME_GIT} > 0 )); then
        if [ -n "$(which curl)" ]; then
            local RESULT=$(curl -sk https://apipe-ci.gsc.wustl.edu/job/1%20Genome%20Perl%20Tests/lastCompletedBuild/api/xml?xpath=/freeStyleBuild/result | grep -vP '^XPath' | sed -r 's/^<\w+>(.*)<\/\w+>/\1/')
            local DESC=$(curl -sk https://apipe-ci.gsc.wustl.edu/job/1%20Genome%20Perl%20Tests/lastCompletedBuild/api/xml?xpath=/freeStyleBuild/description | grep -vP '^XPath' | sed -r 's/^<\w+>(.*)<\/\w+>/\1/')
            local FAILCOUNT=$(curl -sk https://apipe-ci.gsc.wustl.edu/job/1%20Genome%20Perl%20Tests/lastCompletedBuild/api/xml?xpath=/freeStyleBuild/action/failCount | grep -vP '^XPath' | sed -r 's/^<\w+>(.*)<\/\w+>/\1/')
            local NO_COLOR="\033[0m"
            local BLUE="\033[0;34m"
            local CYAN="\033[0;36m"
            local RED="\033[0;31m"
            local YELLOW="\033[1;33m"
            if [ "${RESULT}" = "FAILURE" ]; then
                echo " (RED)"
            fi
            if [ "${RESULT}" = "UNSTABLE" ]; then
                echo " (YELLOW)"
            fi
            if [ "${RESULT}" = "SUCCESS" ]; then
                echo -e " (${CYAN}*${NO_COLOR})"
            fi
        fi
    fi
}

last_perl_tests_color() {
    local IS_GENOME_GIT=$(git describe --tags 2> /dev/null | grep -P '^genome-' | wc -l)
    local NO_COLOR="\033[0m"
    if (( ${IS_GENOME_GIT} > 0 )); then
        if [ -n "$(which curl)" ]; then
            local RESULT=$(curl -sk https://apipe-ci.gsc.wustl.edu/job/1%20Genome%20Perl%20Tests/lastCompletedBuild/api/xml?xpath=/freeStyleBuild/result | grep -vP '^XPath' | sed -r 's/^<\w+>(.*)<\/\w+>/\1/')
            local BLUE="\033[0;34m"
            local CYAN="\033[0;36m"
            local RED="\033[0;31m"
            local YELLOW="\033[1;33m"
            if [ "${RESULT}" = "FAILURE" ]; then
                echo $RED
            fi
            if [ "${RESULT}" = "UNSTABLE" ]; then
                echo $YELLOW
            fi
            if [ "${RESULT}" = "SUCCESS" ]; then
                echo $CYAN
            fi
        else
            # fallback to no color if not able to determine color
            echo $NO_COLOR
        fi
    fi
}

parse_git_branch() {
    local IS_GENOME_GIT=$(git remote -v 2> /dev/null | grep origin | grep genome | wc -l)
    local BRANCH=$(git branch 2> /dev/null | sed -e '/^[^*]/d' -e 's/^\*\ //')
    local NO_COLOR="\033[0m"
    if (( ${IS_GENOME_GIT} )); then
        local BRANCH_COLOR=$(last_perl_tests_color)
    fi
    if [ -n "${BRANCH}" ]; then
        echo -e " (${BRANCH_COLOR}${BRANCH}${NO_COLOR})"
    fi
}

export PS1="${PS1}\$(parse_git_branch) $ "
