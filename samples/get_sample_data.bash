#! /usr/bin/env bash

set -e

# save a few directory names
CWD=$PWD
SAMPLES="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Get the sample data
cd $SAMPLES/VIC_sample_data
git fetch origin
# TODO - get data branch/tag that aligns with the state of the VIC repo.
git checkout master

# back to the original directory
cd $CWD

# Check that the
if [ ! -d "~/workdir" ]; then
  echo "WARNING: `~/workdir` does not exist. Sample data will not run without modifying the global parameter file."
fi
