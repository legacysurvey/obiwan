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
This is where the "slurm_job.sh" copy comes in. Run the job
```sh
cd $obiwan_out/9deg
sbatch $obiwan_code/obiwan/bin/slurm_job_9deg.sh
```
when it finishes and it appears to have worked, remove its outputs
```sh
rm -r $obiwan_out/9deg/elg/123/1238p245/rs0
```

### Finally, launch with QDO
```sh
cd $obiwan_out/9deg
qdo launch obiwan 40 --cores_per_worker 4 --batchqueue debug --walltime 00:30:00 --script $obiwan_code/obiwan/bin/qdo_job_9deg.sh --keep_env
```


