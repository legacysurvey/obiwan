#!/bin/bash

export QDO_BATCH_PROFILE=cori
#qdo_dir=/global/cscratch1/sd/kaylanb/test/qdo
#qdo_dir=/global/cscratch1/sd/kaylanb/test/qdo_install
#export PATH=${qdo_dir}/bin:${PATH}
#export PYTHONPATH=${qdo_dir}:${PYTHONPATH}

#export QDORC_LOCATION=$qdo_dir
export QDO_BACKEND=postgres
export QDO_DB_NAME=desirun
export QDO_DB_HOST=scidb2.nersc.gov
export QDO_DB_USER=desirun_admin
#export QDO_DB_PASS  --> see bashrc.ext

