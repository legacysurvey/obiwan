#!/bin/bash

source /srv/py3_venv/bin/activate
export PYTHONPATH=$CSCRATCH/repos_for_docker/obiwan/py:$CSCRATCH/repos_for_docker/legacypipe/py:$PYTHONPATH

python $CSCRATCH/repos_for_docker/obiwan/tests/test_200x200_pixel_regions.py 


