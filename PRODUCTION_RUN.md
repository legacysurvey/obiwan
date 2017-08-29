# Instructions for doing a production run

### Randoms + PSQL DB

Given KDE trained models, saved as .pkl files
```sh
ls $obiwan_out/randoms
elg-kde.pickle lrg-kde.pickle
```

Generate "ndraws" random ra,dec points on the unit sphere each point also sampling the KDE model. Note the followinging can be run with mpi4py if millions of randoms. For ELGs over a 9deg2 region,
```sh
python $obiwan_code/obiwan/py/obiwan/draw_radec_color_z.py --dowhat sample --obj elg --ra1 122.3 --ra2 125.3 --dec1 23.0 --dec2 26.0 --ndraws 240000 --outdir $CSCRATCH/obiwan_out/randoms
```

Load the randoms into the desi database at scidb2.nersc.gov
```
cd $obiwan_code/bin
bash $obiwan_code/bin/run_fits2db.sh obiwan_elg_9deg $obiwan_out/randoms/elg_randoms/rank_1_seed_1.fits 
db_table=obiwan_elg_9deg
```

Index and cluster the db table for fast ra,dec querying
```sh
psql -U desi_admin -d desi -h scidb2.nersc.gov
desi=> \i /global/cscratch1/sd/kaylanb/obiwan_code/obiwan/etc/db_cluster 
```

### Prepare QDO runs

Get list of bricks touching the regions of randoms, also turn this into a task list for qdo 
```sh
python obiwan/py/obiwan/runmanager/qdo_tasks.py --ra1 123.3 --ra2 125.3 --dec1 23.0 --dec2 26.0 --nobj_total 240000 --nobj_per_run 300
```
Note, there are 240000 randoms in the db and we are configuring the qdo task list so that 300 objects are injected every time obiwan runs

Create the qdo queue
```sh
qdo create obiwan_9deg
qdo load obiwan_9deg tasks_inregion.txt
```

Make a copy of "qdo_job.sh" and "slurm_job.sh" and edit the copies for the specified run case
```sh
cd $obiwan_code/bin
cp qdo_job.sh qdo_job_9deg.sh
cp slurm_job.sh slurm_job_9deg.sh
```
WARNING make sure to tell obiwan which table in the desi db to use
```sh
python obiwan/kenobi.py --dataset dr5 ... --db_table obiwan_elg_9deg
```

### Run a 5 min debug job to test the setup
This is where the "slurm_job.sh" copy comes in. Lets say the desi db randoms table is called "obiwan_elg_9deg", then run the slurm_job.sh job as so
```sh
cd $obiwan_out/obiwan_elg_9deg
sbatch $obiwan_code/obiwan/bin/slurm_job_9deg.sh
```
when it finishes and it appears to have worked, remove its outputs
```sh
rm -r $obiwan_out/obiwan_elg_9deg/elg/123/1238p245/rs0
```

### Finally, launch with QDO
```sh
cd $obiwan_out/obiwan_elg_9deg
qdo launch obiwan 40 --cores_per_worker 4 --batchqueue debug --walltime 00:30:00 --script $obiwan_code/obiwan/bin/qdo_job_9deg.sh --keep_env
qdo launch obiwan 40 --cores_per_worker 4 --batchqueue regular --walltime 05:00:00 --script $obiwan_code/obiwan/bin/qdo_job_9deg.sh --keep_env
```

### Managing your qdo production run
The scripts are in `obiwan/py/obiwan/runmanager/*.py`. I ran the 9 deg2 test using two qdo ques, named `obiwan_9deg` for all `rs*/` jobs and (when those all completed) `obiwan_9deg_doskip` for all `skip_rs*/` jobs.

(1) To get a list of all log.<brickname> and slurm-<slurmid>.out files, sorted by status of "succeeded, failed, running" in the qdo db, and a tally of each error that occurred, do
```sh
cd $obiwan_code
python $obiwan_code/obiwan/py/obiwan/runmanager/status.py --qdo_quename obiwan_9deg --outdir /global/cscratch1/sd/kaylanb/obiwan_out/obiwan_elg_9deg --obj elg
```
replace "obiwan_9deg" by "obiwan_9deg_doskip" to get the same for the `skip_rs*` jobs.

(2) If you have a region on the sky where you want to run obiwan, you can get a list of all the tasks you need to do so with
```py
from obiwan.runmanager.qdo_tasks import write_qdo_tasks_normal
write_qdo_tasks_normal(ra1=123.3,ra2=124.3,dec1=24.0,dec2=25.0,
                       nobj_total=2400, nobj_per_run=500)
```
where ra1,ra2,dec1,dec2 are the corners of your region. nobj_total is the total number of randoms in the psql randoms db for your region, and nobj_per_run is how many randoms to inject per obiwan run.

Once you finish all the above runs, you need to inject the randoms that were skipped. You can get a list of all these tasks with a text file listing the bricknames to inject the skipped randoms into
```py
from obiwan.runmanager.qdo_tasks import write_qdo_tasks_skipids
write_qdo_tasks_skipids(brick_list_fn, nobj_per_run=500)
```


