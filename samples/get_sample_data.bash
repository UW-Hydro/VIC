#! /usr/bin/env bash

set -e

# save a few directory names
CWD=$PWD
SAMPLES="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

git submodule init

# Get the sample data
cd $SAMPLES/data
git fetch origin
# TODO - get data branch/tag that aligns with the state of the VIC repo.
git checkout master

# back to the original directory
cd $CWD
