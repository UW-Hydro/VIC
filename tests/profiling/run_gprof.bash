#!/usr/bin/env bash
set -e
# set -x

# A POSIX variable
OPTIND=1         # Reset in case getopts has been used previously in the shell.

# Initialize our own variables:
vic_exe=""
vic_global=""

function usage {
  echo "Usage: `basename $0` -e vic_executable -g vic_global -h for help";
}


while getopts "h?e::g::" opt; do
    case "$opt" in
    h|\?)
        usage
        exit 0
        ;;
    e)  vic_exe=$OPTARG ;;
    g)  vic_global=$OPTARG ;;
    esac
done

shift $((OPTIND-1))

[ "$1" = "--" ] && shift

GIT_VERSION=$(git describe --abbrev=4 --dirty --always --tags)

echo "----------------------- START VIC GPROF PROFILING -----------------------"
echo ""
echo "Date                      : $(date)"
echo "Machine                   : $HOSTNAME"
echo "User                      : $USER"
echo "VIC Test Git Version      : $GIT_VERSION"
echo "VIC Executable            : $vic_exe"
echo "VIC Global Parameter File : $vic_global"
echo ""

echo "VIC Executable Version Info"
echo "---------------------------"
$vic_exe -v
echo ""

$vic_exe -g $vic_global 2>&1 | tail -n 55

now=`date +"%y%m%d"`
gprof $vic_exe | gprof2dot | dot -Tpng -o vic_call_graph_$now.png

echo ""
echo "------------------------ END VIC GPROF PROFILING ------------------------"
exit 0

# End of file
