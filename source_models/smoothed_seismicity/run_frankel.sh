#PBS -P w84
#PBS -q express
#PBS -l walltime=16:00:00
#PBS -l ncpus=9
#PBS -l mem=32GB
#PBS -l wd

module load intel-cc/12.1.9.293
module load intel-fc/12.1.9.293
module load geos/3.5.0
module load hdf5/1.8.10
#module load openmpi/1.6.3
module load gdal
module unload python
module load python/2.7.11
module load python/2.7.11-matplotlib

# To get rtree to run
export SPATIALINDEX_C_LIBRARY=/short/n74/jdg547/spatialindex-src-1.8.5/lib/libspatialindex_c.so.4
export LD_LIBRARY_PATH=/short/n74/jdg547/spatialindex-src-1.8.5/lib:$LD_LIBRARY_PATH
# Python paths for local openquake installs and dependencies
export PYTHONPATH=.:/home/547/jdg547/.local/lib/python2.7/site-packages:${PYTHONPATH}
export PYTHONPATH=.:/short/w84/NSHA18/sandpit/jdg547/oq-hazardlib:${PYTHONPATH}
export PYTHONPATH=.:/short/w84/NSHA18/sandpit/jdg547/oq-engine:${PYTHONPATH}
export PYTHONPATH=.:/short/w84/NSHA18/sandpit/jdg547/hmtk:${PYTHONPATH}
#export PYTHONPATH=.::/short/w84/NSHA18/sandpit/jdg547/NSHA2018/:${PYTHONPATH}
export PYTHONPATH=.::/short/w84/NSHA18/sandpit/jdg547/:${PYTHONPATH}

counter=0
one=1
bvals=(0.819 0.835 1.208 0.747 0.727 1.062 0.892 0.944 1.355)
for bval in ${bvals[*]}; do
    python fixed_smoothing.py $bval > frankel_b$bval.log &
    counter=$(($counter+$one));
    # Once we have submitted jobs for all 9 bvalues, break from this loop.
    if [ $counter = 9 ];
    then
        break
    fi
done      

wait