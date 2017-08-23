#!/bin/bash

# RUN: ./run_fits2db <db_table_name> <fits_data_table>
# COMPILE: fitsdb.c
#Cori: gcc -o fits2db -I/global/cscratch1/sd/kaylanb/software/desi/desiconda/2017-04-13_aux/include fits2db.c -L/global/cscratch1/sd/kaylanb/software/desi/desiconda/2017-04-13_aux/lib -lcfitsio -lm
#Edison: ?

db_table=${1}
fits_table=${2}
echo db_table=$db_table
echo fits_table=$fits_table
#Make table
/global/cscratch1/sd/kaylanb/test/fits2db/fits2db_cori --sql=postgres --create --noload -C -X -t ${db_table} ${fits_table}| /usr/bin/psql -U desi_admin -d desi -h scidb2.nersc.gov

#Load data
#for fn in table1.fits table2.fits;/global/cscratch1/sd/kaylanb/test/fits2db/fits2db_cori --sql=postgres -C -X -t heyo2 $fn| /usr/bin/psql -U desi_admin -d desi -h scidb2.nersc.gov;done
/global/cscratch1/sd/kaylanb/test/fits2db/fits2db_cori --sql=postgres -C -X -t ${db_table} ${fits_table}| /usr/bin/psql -U desi_admin -d desi -h scidb2.nersc.gov

echo Done


