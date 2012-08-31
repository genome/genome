#!/bin/bash

GIT_SRV=git
GIT_PATH=/srv/git
SKEL_URL=ssh://$GIT_SRV$GIT_PATH/example-c-project.git
DO_PUSH=1

[ -z "$TMPDIR" ] && TMPDIR=/tmp

function bail() {
    echo $1
    exit 1
}

function usage() {
    echo "Usage: `basename $0` [-n] <project name>"
    echo "  Initialize a new c/c++ project"
    echo ""
    echo "Options:"
    echo "  -h|--help    - this message"
    echo "  -n|--no-push - do not automatically commit and push the new subrepo"
}

function enforce_clean_master() {
    branch=`git branch | perl -ne 's/^\* // && print'`
    if [ "$branch" != "master" ]
    then
        bail "This program must be run on branch master, you are on branch $branch"
    fi
    git fetch || bail "git fetch failed"
    git status -uno | grep -A1 '^# Your branch ' &&
        bail "Your branch is not in sync with the remote repo, pull/push first"

    git push --dry-run > /dev/null 2>&1 ||
        bail "git push --dry-run failed."
}

function create_git_repo() {
    name=$1 
    repo_path=$GIT_PATH/$name.git
    ssh $GIT_SRV "if [ ! -e $repo_path ]; then mkdir $repo_path && cd $repo_path && git init --bare --shared; else echo $repo_path already exists on server $GIT_SRV; exit 1; fi"
    [ $? -ne 0 ] && {
        echo "Aborted."
        exit 1
    }

    (pushd $TMPDIR && git clone ssh://$GIT_SRV/$repo_path &&
    pushd $name && touch .gitignore &&
    git add .gitignore && git commit -m 'added .gitignore' && git push origin master &&
    popd && rm -rf $name) || {
        echo "Failed to set up repository $name"
        exit 1
    }

    git submodule add ssh://$GIT_SRV/$repo_path || {
        echo "Failed to set up submodule for $name"
        exit 1
    }
}

function create_project_skel() {
    git archive --remote $SKEL_URL master | tar xf -
}

# main -------------------------------------------------------------------

enforce_clean_master

project_name=
while [ $# -gt 0 ]; do
    case $1 in
        -h|--help) usage; exit 0 ;;
        -n|--no-push) DO_PUSH=0 ;;
        -*) bail "Unexpected argument: $1" ;;
        *) [ "$project_name" != "" ] && {
                bail "Multiple project names specified: '$project_name' and '$1'"
            }
            project_name=$1
            ;;
    esac
    shift
done

if [ -z "$project_name" ]; then
    echo "Error: no project name specified"
    usage
    exit 1
fi

[ -e $project_name ] && {
    echo "Error: something named $project_name already exists in `pwd`"
    ls -ld $project_name
    exit 1
}

echo ""
echo "Create new project named $project_name? (y/n)"
while read ans; do
    case $ans in
        y*) break;;
        n*) echo "Aborted"; exit 1;;
        *) echo "Please anwser y or n";;
    esac
done

create_git_repo $project_name
(pushd $project_name && create_project_skel $project_name && ls -l)
git add $project_name
[ $DO_PUSH -eq 1 ] && {
    git pull
    git commit -m "Added submodule $project_name"
    git push
    echo "** pushed changes"
}
git status

echo "Completed."
