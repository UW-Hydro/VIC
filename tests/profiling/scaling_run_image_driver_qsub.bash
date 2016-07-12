#!/usr/bin/env bash
set -e
set -x

# A POSIX variable
OPTIND=1         # Reset in case getopts has been used previously in the shell.

# Initialize our own variables:
vic_exe=""
vic_global=""
max_cores="64"
qsub_template=""

dir="$(dirname $0)"   # Returns "/from/hear/to"
dt=$(date '+%Y%m%d');
timing_table_file="$dir/timing.$dt.txt"


function usage {
  echo "Usage: `basename $0` -e vic_executable -g vic_global -n $max_cores -t qsub_template template_filename --h for help";
}

while getopts "h?e::g::t::n:" opt; do
    case "$opt" in
    h|\?)
        usage
        exit 0
        ;;
    n)  max_cores=$OPTARG ;;
    e)  vic_exe=$OPTARG ;;
    g)  vic_global=$OPTARG ;;
    t)  qsub_template=$OPTARG ;;
    esac
done

shift $((OPTIND-1))

[ "$1" = "--" ] && shift

GIT_VERSION=$(git describe --abbrev=4 --dirty --always --tags)
VIC_VERSION=$($vic_exe -v)

echo "Scaling Output --> $timing_table_file"

cat >$timing_table_file <<EOL
----------------- START VIC DISTRIBUTED SCALING PROFILE -----------------

Date                      : $(date)
Machine                   : $HOSTNAME
User                      : $USER
VIC Test Git Version      : $GIT_VERSION
VIC Executable            : $vic_exe
VIC Global Parameter File : $vic_global
Test Maximum MPI Cores    : $max_cores

VIC Executable Version Info
---------------------------
$VIC_VERSION

Cores | Time (Seconds)
----------------------
EOL


np="1"
while [ $np -le $max_cores ]
do
  run_script="$dir/vic_scaling_test.$np.qsub"
  qsub -V $qsub_template
  np=$[$np*2]
done

exit 0

# End of file
