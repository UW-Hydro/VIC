#!/usr/bin/env bash
# Download forcings for 0.5 degree global simulation

set -e

machine='other'

if [ -n "${WORKDIR+1}" ]; then
    forcings_path=$WORKDIR/global_0.5_degree_forcings
    mkdir -p $forcings_path
    echo "forcings path:  $forcings_path"
else
    echo "WORKDIR is not defined"
    exit 1
fi

forcings_path=$WORKDIR/global_0.5_degree_forcings/
mkdir -p $forcings_path

if [[ "$machine" == 'hydra' ]]; then
    cp /home/ftp/pub/global/cells/*tgz $forcings_path
elif [[ "$machine" == 'hydro' ]]; then
    cp /nfs/dynamo/ftp/pub/global/cells/*tgz $forcings_path
else
    data_path="ftp://ftp.hydro.washington.edu/pub/global/cells/"
    wget --mirror -P $forcings_path -nH --cut-dirs=3 -A "*tgz" $data_path  
fi

# unar the forcing blocks
for tarfile in $forcings_path/*tgz;
do
    tar -C $forcings_path -z -x -v -f $tarfile
done

# remove tarfiles
rm $forcings_path/*tgz

for dir in $forcings_path/*;
do 
    mv $dir/* $forcings_path
    rm -r $dir
done

echo "Done getting the global forcings"
echo "Forcing files:  $forcings_path"
echo "Forcings readme file:  http://www.hydro.washington.edu/SurfaceWaterGroup/Data/met_global_0.5deg.html"
