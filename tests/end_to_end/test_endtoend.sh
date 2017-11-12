#!/bin/bash

export DUST_DIR=$HOME/myrepo/dust
export DECALS_SIM_DIR=$HOME/mydata/obiwan/obiwan
export LEGACY_SURVEY_DIR=$HOME/mydata/obiwan/legacy_survey_dir
#export UNWISE_COADDS_DIR=/global/cscratch1/sd/desiproc/unwise-coadds/fulldepth:/global/cscratch1/sd/desiproc/unwise-coadds/w3w4
#export UNWISE_COADDS_TIMERESOLVED_DIR=/global/cscratch1/sd/desiproc/unwise-coadds/time_resolved_neo2
#export UNWISE_COADDS_TIMERESOLVED_INDEX /global/cscratch1/sd/desiproc/unwise-coadds/time_resolved_neo2/time_resolved_neo2-atlas.fits

#export PATH=$HOME/myinstall/astrometry/lib/python/astrometry/util:$PATH
#export PATH=$HOME/myinstall/astrometry/lib/python/astrometry/blind:$PATH
#export PYTHONPATH=$HOME/myinstall/astrometry/lib/python:$PYTHONPATH
#export PYTHONPATH=$HOME/myinstall/tractor/lib/python2.7/site-packages:$PYTHONPATH

#export PYTHONPATH=$HOME/myrepo/legacypipe/py:$PYTHONPATH
#export PYTHONPATH=$HOME/myrepo/theValidator:$PYTHONPATH
#export PYTHONPATH=.:${PYTHONPATH}

python obiwan/kenobi.py -b 1238p245 -n 2 --DR 5 -o elg --add_sim_noise --zoom 1550 1650 1550 1650

num_fils=`find ${DECALS_SIM_DIR}/elg -type f|wc -l`
# 0 is True, 1 if False
if [ "${num_fils}" == "34" ]; then
    exit 0
else
    exit 42
fi
