#!/usr/bin/env bash
set -e
# set -x

# A POSIX variable
OPTIND=1         # Reset in case getopts has been used previously in the shell.

# Initialize our own variables:
vic_exe=""
vic_global=""
max_cores="8"

function usage {
  echo "Usage: `basename $0` -e vic_executable -g vic_global -n $max_cores -h for help";
}


while getopts "h?e::g::n:" opt; do
    case "$opt" in
    h|\?)
        usage
        exit 0
        ;;
    n)  max_cores=$OPTARG ;;
    e)  vic_exe=$OPTARG ;;
    g)  vic_global=$OPTARG ;;
    esac
done

shift $((OPTIND-1))

[ "$1" = "--" ] && shift

GIT_VERSION=$(git describe --abbrev=4 --dirty --always --tags)

echo "-------------------- START VIC LOCAL SCALING PROFILE --------------------"
echo ""
echo "Date                      : $(date)"
echo "Machine                   : $HOSTNAME"
echo "User                      : $USER"
echo "VIC Test Git Version      : $GIT_VERSION"
echo "VIC Executable            : $vic_exe"
echo "VIC Global Parameter File : $vic_global"
echo "Test Maximum MPI Cores    : $max_cores"
echo ""

echo "VIC Executable Version Info"
echo "---------------------------"
$vic_exe -v
echo ""

printf "%5s | %s\n" "Cores" "Time (Seconds)"
echo "----------------------"

START=$(date +%s)
$vic_exe -g $vic_global &> /dev/null
END=$(date +%s)
DIFF=$(echo "$END - $START" | bc)
printf "%5s | %f\n" "--" $DIFF

np="1"
while [ $np -le $max_cores ]
do
  START=$(date +%s)
  mpiexec -np $np $vic_exe -g $vic_global &> /dev/null
  END=$(date +%s)
  DIFF=$(echo "$END - $START" | bc)
  printf "%5s | %f\n" $np $DIFF
  np=$[$np*2]
done

echo ""
echo "--------------------- END VIC LOCAL SCALING PROFILE ---------------------"
exit 0

# End of file
