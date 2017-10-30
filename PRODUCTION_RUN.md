# Instructions for doing a production run

### Completed production runs

(1) elg_9deg2_ra125:
9deg2 region centered at ra,dec= 124.8,24.5
* name_for_run: elg_9deg2_ra125
* region: ra1,ra2, dec1,dec2= 122.3,125.3, 23.0,26.0
* kdedir for randoms: $CSCRATCH/obiwan_out/randoms
* outdir for randoms: $CSCRATCH/obiwan_out/randoms/elg_9deg2_ra125
* psql db, desi, table for randoms: obiwan_elg_9deg
* psql db, qdo, table for tasks: obiwan_9deg, obiwan_9deg_doskip

(2) elg_9deg2_ra175:
9deg2 region centered at ra,dec= 175.0,24.5
* name_for_run: elg_9deg2_ra175
* region: ra1,ra2, dec1,dec2= 173.5,176.5, 23.0,26.0
* kdedir for randoms: $CSCRATCH/obiwan_out/randoms
* outdir for randoms: $CSCRATCH/obiwan_out/randoms/elg_9deg2_ra175
* psql db, desi, table for randoms: obiwan_elg_ra175
* psql db, qdo, table for tasks: obiwan_ra175, obiwan_ra175_doskip

The rest of this markdown is for `elg_9deg2_ra175`.

### Randoms + PSQL DB

Given KDE trained models, saved as .pkl files
```sh
ls $obiwan_out/randoms
elg-kde.pickle lrg-kde.pickle
```

Generate "ndraws" random ra,dec points on the unit sphere each point also sampling the KDE model. Note the followinging can be run with mpi4py if millions of randoms. 
```sh
export ra1=173.5
export ra2=176.5
export dec1=23.0
export dec2=26.0
export kdedir=$CSCRATCH/obiwan_out/randoms
export outdir=$CSCRATCH/obiwan_out/randoms/elg_9deg2_ra175
python $obiwan_code/obiwan/py/obiwan/draw_radec_color_z.py --dowhat sample --obj elg --ra1  --ra2 ${ra2} --dec1 ${dec1} --dec2 ${dec2} --ndraws 240000 --kdedir ${kdedir} --outdir ${outdir}
```

Load the randoms into the desi database at nerscdb03 (the scidb2.nersc.gov server is no more)
```
bash $obiwan_code/obiwan/bin/run_fits2db.sh obiwan_elg_ra175 /global/cscratch1/sd/kaylanb/obiwan_out/randoms/elg_9deg2_ra175/rank_1_seed_1.fits
```

Index and cluster the db table for fast ra,dec querying
```sh
cd $obiwan_code/obiwan/etc
cat cluster_randoms | sed s#name#obiwan_elg_ra175#g > cluster_temp
psql -U desi_admin -d desi -h scidb2.nersc.gov
desi=> \i /global/cscratch1/sd/kaylanb/obiwan_code/obiwan/etc/cluster_temp
```

Also make sure the bricks file is in the desi db. Do
```sh
psql_desi 
desi=> \d obiwan_bricks
```
If its not there do
```sh
cd $HOME/
rsync -av /global/project/projectdirs/desi/www/users/kburleigh/obiwan/legacysurveydir/survey-bricks.fits.gz .
gunzip survey-bricks.fits.gz
bash $obiwan_code/bin/run_fits2db.sh obiwan_bricks survey-bricks.fits 
```
Index and cluster it
```sh
psql -U desi_admin -d desi -h scidb2.nersc.gov
desi=> \i /global/cscratch1/sd/kaylanb/obiwan_code/obiwan/etc/cluster_bricks
```

### Prepare QDO runs

*1) Make the QDO task list*
Generate the qdo tasks to run, which includes the list of bricks taht touch your radec region and how many randoms to inject per task 
```sh
>>> from obiwan.runmanager.qdo_tasks import TaskList
>>> T= TaskList(ra1=173.5,ra2=176.5, dec1=23.0,dec2=26.0,
                nobj_total=240000, nobj_per_run=300)
>>> T.bricks()
>>> T.tasklist(do_skipid='no',do_more='no')
```
where ra1,ra2,dec1,dec2 are the corners of your region. nobj_total is the total number of randoms in the psql randoms db for your region, and nobj_per_run is how many randoms to inject per obiwan run.

Note, there are 240000 randoms in the db and we are configuring the qdo task list so that 300 objects are injected every time obiwan runs

Create the qdo queue
```sh
qdo create obiwan_ra175 
qdo load obiwan_ra175 tasks_inregion.txt
```

Then create that qdo queue
```sh
qdo create obiwan_ra175_domore
qdo load obiwan_ra175_dormore tasks_skipid_no_more_yes_minid_240001.txt
```

*2) configure the slurm job*

Make a copy of "qdo_job.sh" and "slurm_job.sh",
```sh
cd $obiwan_code/bin
cp qdo_job.sh qdo_job_9deg.sh
cp slurm_job.sh slurm_job_9deg.sh
```
and the following fields in both job files need to be updated (`do_skipids` only needs to be updated for slurm_job.sh, its auotmatic for qdo_job.sh)
```sh
export name_for_run=elg_9deg2_ra175
export randoms_db=elg_9deg2_ra175
export do_skipids=no
export do_more=no
```
Additonally, the following fields may need to be updated depending on your run
```sh
export dataset=dr5
export brick=1750p225
export object=elg
export nobj=300
```

### Run a 5 min debug job to test that everything works on the compute nodes
```sh
export outdir=$obiwan_out/${name_for_run}
cd $outdir
sbatch $obiwan_code/obiwan/bin/slurm_job_9deg.sh
```
when it finishes, grep the open the resulting `slurm*.out` file, find the file it says it is "logging to", and grep that file for the success string
```sh
grep "decals_sim:All done!" <logging to file>
```
If it finds it, then it work and you can remove the output directory
```
rm -r $obiwan_out/${name_for_run}/elg/${bri}/${brick}/rs0
```

### Finally, launch with QDO
```sh
cd $obiwan_out/${name_for_run}
export qdo_quename=obiwan_ra175
qdo launch ${qdo_quename} 40 --cores_per_worker 4 --batchqueue regular --walltime 05:00:00 --script $obiwan_code/obiwan/bin/qdo_job_9deg.sh --keep_env
```

### Add more randoms mid-run
Eventually you'll need to add more randoms. For instance if after finishing all QDO runs, the randoms you recover in the simulated tractor catalogues have less than 10x target density.

To add more randoms repeat previous instructions but add the "--startid" option. `240,000` randoms were added initially. Each gets a primary key from 1 to the number of randoms. So the randoms you add mid-run need to have primary keys that start at `240,001`. Lets make `720,000` more
```sh
export startid=240001
python $obiwan_code/obiwan/py/obiwan/draw_radec_color_z.py --dowhat sample --obj elg --ra1  --ra2 ${ra2} --dec1 ${dec1} --dec2 ${dec2} --ndraws 240000 --kdedir ${kdedir} --outdir ${outdir} --startid ${startid}
```

It is easiest to load the additional randoms to a new temporary table in the DB then insert that table's rows into the randoms DB. If your new randoms fits table is `/global/cscratch1/sd/kaylanb/obiwan_out/randoms/elg_9deg2_ra175/more.fits` then 
```sh
bash $obiwan_code/obiwan/bin/run_fits2db.sh obiwan_test /global/cscratch1/sd/kaylanb/obiwan_out/randoms/elg_9deg2_ra175/more.fits
```

Now add `obiwan_test` to the end of the randoms table and delete `obiwan_test`
```sh
desi=> insert into obiwan_elg_ra175 select * from obiwan_test;
desi=> drop table obiwan_test;
```

Make a QDO task list for your additional randoms. Specify `minid` to skip all primary keys below your 240,001
```sh
>>> from obiwan.runmanager.qdo_tasks import TaskList
>>> T= TaskList(ra1=173.5,ra2=176.5, dec1=23.0,dec2=26.0,
                nobj_per_run=300,
                nobj_total=240000 + 720000)
>>> T.bricks()
>>> T.tasklist(do_skipid='no',do_more='yes',minid=240001)
```
then run these randoms from a new QDO queue
```sh
qdo create obiwan_ra175_domore
qdo load obiwan_ra175_dormore tasks_skipid_no_more_yes_minid_240001.txt
```

### Managing your qdo production run
Manage your qdo production run with `obiwan/py/obiwan/runmanager/status.py`. To get a list of all log.<brickname> and slurm-<slurmid>.out files, sorted by status of "succeeded, failed, running" in the qdo db, and a tally of each error that occurred, do
```sh
cd $obiwan_code
python $obiwan_code/obiwan/py/obiwan/runmanager/status.py --qdo_quename ${qdo_quename} --outdir /global/cscratch1/sd/kaylanb/obiwan_out/${name_for_run} --obj elg
```

Once you finish all the above runs, you will make a second qdo que. We previously made `obiwan_ra175` to do the usual obiwan runs, and now we make `obiwan_ra175_doskip` to do the randoms that we skipped. You can get a list of these tasks with
```sh
cat tasks_inregion | awk '{print $1}'|sort|uniq > brick_list.txt
python
>>> from obiwan.runmanager.qdo_tasks import write_qdo_tasks_skipids
>>> write_qdo_tasks_skipids('brick_list.txt', nobj_per_run=300)
```
which outputs a file `tasks_skipids.txt`

Create a new qdo queue_name for the skipid runs and load the new tasks
```sh
export qdo_quename=obiwan_ra175_doskip
qdo create ${qdo_quename} 
qdo load ${qdo_quename} tasks_skipids.txt
```

The qdo tasks automatically set the `do_skipid` flag, so you dont need to edit the `qdo_job_9deg.sh` file. Just run it with your new qdo `que_name`
```sh
cd $obiwan_out/${name_for_run}
qdo launch ${qdo_quename} 40 --cores_per_worker 4 --batchqueue regular --walltime 05:00:00 --script $obiwan_code/obiwan/bin/qdo_job_9deg.sh --keep_env
```


