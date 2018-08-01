Run on 1000 deg2 of sky
==========================================

On-going production runs
^^^^^^^^^^^^^^^^^^^^^^^^^

.. _elg-dr5:

#. ELGs for 1/2 of DR5

    - /global/cscratch1/sd/kaylanb/obiwan_out/elg_dr5
    - name "elg_dr5"
    - ra [109.0,278.5], dec [-11.1,35.4]
    - all dr5 bricks in Middle region (90 < ra < 280) with grz coverage
    - randoms dir: /global/cscratch1/sd/kaylanb/obiwan_out/randoms/elg_dr5

        * 5x ELG density: 2400 * 10 * area * 4
        * 170*47*2400*5*4/1e9 = 0.383e9
        * times 4 because 0.42 of 10k in ELG box and expect 1/2 not recovered

    * psql db "desi", table for randoms "obiwan_elg_dr5"
    * psql db "qdo", table for tasks "obiwan_elg_dr5"

#. Same as #1 but Stars

Completed production runs
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. _elg-9deg2-ra175:

#. ELGs for 9deg2 region centered at ra,dec= 175.0,24.5

    * /global/cscratch1/sd/kaylanb/obiwan_out/elg_9deg2_ra175
    * name "elg_9deg2_ra175"
    * region: ra1,ra2, dec1,dec2= 173.5,176.5, 23.0,26.0
    * kdedir for randoms: $CSCRATCH/obiwan_out/randoms
    * outdir for randoms: $CSCRATCH/obiwan_out/randoms/elg_9deg2_ra175
    * psql db "desi", table for randoms "obiwan_elg_ra175"
    * psql db "qdo", table for tasks "obiwan_ra175", "obiwan_ra175_doskip"

#. ELGs for 9deg2 region centered at ra,dec= 124.8,24.5

    * /global/cscratch1/sd/kaylanb/obiwan_out/obiwan_elg_9deg
    * name_for_run: elg_9deg2_ra125
    * region: ra1,ra2, dec1,dec2= 122.3,125.3, 23.0,26.0
    * kdedir for randoms: $CSCRATCH/obiwan_out/randoms
    * outdir for randoms: $CSCRATCH/obiwan_out/randoms/elg_9deg2_ra125
    * psql db, desi, table for randoms: obiwan_elg_9deg
    * psql db, qdo, table for tasks: obiwan_9deg, obiwan_9deg_doskip

Data Sets that Obiwan Supports
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

#. DR5

.. * DR5 zeropoints (from legacy_zeropoints, minimum set of columns needed for legacypip)e)
.. * DR5 CCDs

#. DR3_eBOSS

    * DR3 zeropoints (from IDL zeropoints)many more columns than legacypipe need, and comput
    * DR3 CCDs
    * additional eBOSS CCDs (e.g. beyond MJD cutoff of March 2016)

The "dataset" keyword tells obiwan which dataset to run. For example::

    python obiwan/kenobi.py --dataset {dr5,dr3_eboss}

The corresponding config in Legacypipe:

#. DR5

    * --run dr5
    * only use survey-ccds files ignore "survey-ccds-kd" files

#. DR3_eBOSS

    * --run dr3_eboss
    * --no-rex --use-simp: turn off REX model, use SIMP instead
    * --nsigma 6: default is 6 but enforce this
    * only use the survey-ccds files I made that exactly includes the `DR3_eBOSS` image list using idl zeropoints, also ignore "survey-ccds-kd" files

The NERSC PostgreSQL DB
^^^^^^^^^^^^^^^^^^^^^^^^^

The psql db at NERSC is called "desi". It's hosted at scidb2.nersc.gov and you sign in with user "desi_user". You'll need the ``.pgpass`` password file. Then put the ``.pgpass`` file in ``$HOME/`` on Cori and give it user rw permissions.::

  cp <path/to/kaylans/.pgpass $HOME/
  chmod u+rw $HOME/.pgpass


Make sure the bashrc_obiwan loaded the ``postresql`` module. Then do::

  psql -U desi_user -d desi -h scidb2.nersc.gov

It worked if that brings you to a ``desi=>`` prompt

Randoms
^^^^^^^^^^^^^^^^^^^^^^^^^

We fit a Mixture of Gaussians to the desired n(z) and seperately joint-sample the color,shape,redshift (10k row fits table) from Tractor catalogues matched to a Truth region. The 10k row fits table is here::

    ls $obiwan_out/randoms/elg_100deg2
    elg_sample_5dim_10k.fits

Generate "ndraws" random ra,dec points on the unit sphere each point also sampling the KDE model. Run with multiple core if > 1 Million randoms. Either way fill in and submit this slurm job:
https://github.com/legacysurvey/obiwan/blob/master/bin/slurm_job_randoms.sh

Load the randoms into the desi database at nerscdb03 (the scidb2.nersc.gov server is no more). First create the table::

    bash $obiwan_out/obiwan/bin/fits2db_create.sh obiwan_elg_dr5 /global/cscratch1/sd/kaylanb/obiwan_out/randoms/elg_dr5/randoms_rank_0.fits

then load all the randoms fits tables::

    name=elg_dr5
    num=`find $obiwan_out/randoms/$name -name "randoms*.fits"|wc -l`
    let num=$num-1
    for i in `seq 0 $num`;do bash $obiwan_code/obiwan/bin/fits2db_load.sh obiwan_$name $obiwan_out/randoms/$name/randoms_rank_${i}.fits;done


Index and cluster the db table for fast ra,dec querying::

    cd $obiwan_code/obiwan/etc
    cat cluster_randoms | sed s#name#obiwan_elg_100deg2#g > cluster_temp
    psql -U desi_admin -d desi -h nerscdb03.nersc.gov
    desi=> \i /global/cscratch1/sd/kaylanb/obiwan_code/obiwan/etc/cluster_temp


Also make sure the bricks file is in the desi db. Do::

    psql_desi
    desi=> \d obiwan_bricks

If its not there do::

    cd $HOME/
    rsync -av /global/project/projectdirs/desi/www/users/kburleigh/obiwan/legacysurveydir/survey-bricks.fits.gz .
    gunzip survey-bricks.fits.gz
    bash $obiwan_code/bin/run_fits2db.sh obiwan_bricks survey-bricks.fits

Index and cluster it::

    psql -U desi_admin -d desi -h scidb2.nersc.gov
    desi=> \i /global/cscratch1/sd/kaylanb/obiwan_code/obiwan/etc/cluster_bricks


Prepare QDO runs
^^^^^^^^^^^^^^^^^^^^^^^^^

1. Make the QDO task list

Generate the qdo tasks to run, which includes the list of bricks that are in your radec region. It is too slow to query the randoms db for each brick's number of randoms, so instead estimate as expectation number + 2 StdErros per brick. Run the script like this::

    python $obiwan_code/obiwan/py/obiwan/runmanager/qdo_tasks.py --obj elg --radec 109.0 278.5 -11.1 35.4 --nobj_total 383000000 --survey_bricks /home/kaylan/mydata/survey-bricks.fits.gz --bricks_fn elg_dr5/dr5_bricks_inMid_grz.txt

which writes out the task file. Now create the qdo queue::

    qdo create obiwan_elg_dr5
    qdo load obiwan_elg_dr5 <task-file.txt>


2. Run a single brick to test that everything works

Go to your output directory and copy over the template slurm job script::

    export outdir=$CSCRATCH/obiwan_out/elg_100deg2
    cd $outdir
    cp $obiwan_code/obiwan/bin/slurm_job.sh ./slurm_job_100deg2.sh

Modify the relevant fields for the run, namely::

    export name_for_run=elg_9deg2_ra175
    export randoms_db=elg_9deg2_ra175
    export do_skipids=no
    export do_more=no
    export dataset=dr5
    export brick=1750p225
    export object=elg
    export nobj=300

Run it as a 30 min debug job::

    sbatch slurm_job_100deg2.sh

When it finishes, grep the open the resulting `slurm*.out` file, find the file it says it is "logging to", and grep that file for the success string::

    grep "decals_sim:All done!" <logging to file>

If the success string is there, cleanup the testrun outputs, add the new slurm job script to the obiwan repo, and being the production run::

    rm -r $obiwan_out/${name_for_run}/elg/${bri}/${brick}/rs0
    cp slurm_job_100deg2.sh $obiwan_code/obiwan/bin/
    # cd to obiwan repo and git add, git commit

**3) Production runs with QDO**
Copy over the template qdo job script,::

    cd $outdir
    cp $obiwan_code/obiwan/bin/qdo_job.sh ./qdo_job_100deg2.sh

and edit the relevant fileds as before. Now launch the qdo jobs::

    export qdo_quename=obiwan_elg_100deg
    qdo launch ${qdo_quename} 40 --cores_per_worker 4 --batchqueue regular --walltime 05:00:00 --script $outdir/qdo_job_100deg2.sh --keep_env

Once you see successful runs,::

    cp qdo_job_100deg2.sh $obiwan_code/obiwan/bin/
    # cd to obiwan repo and git add, git commit


Add more randoms mid-run
^^^^^^^^^^^^^^^^^^^^^^^^^
Eventually you'll need to add more randoms. For instance if after finishing all QDO runs, the randoms you recover in the simulated tractor catalogues have less than 10x target density.

To add more randoms repeat previous instructions but add the "--startid" option. `240,000` randoms were added initially. Each gets a primary key from 1 to the number of randoms. So the randoms you add mid-run need to have primary keys that start at `240,001`. Lets make `720,000` more::

    export startid=240001
    python $obiwan_code/obiwan/py/obiwan/draw_radec_color_z.py --dowhat sample --obj elg --ra1  --ra2 ${ra2} --dec1 ${dec1} --dec2 ${dec2} --ndraws 240000 --kdedir ${kdedir} --outdir ${outdir} --startid ${startid}

It is easiest to load the additional randoms to a new temporary table in the DB then insert that table's rows into the randoms DB. If your new randoms fits table is `/global/cscratch1/sd/kaylanb/obiwan_out/randoms/elg_9deg2_ra175/more.fits` then::

    bash $obiwan_code/obiwan/bin/run_fits2db.sh obiwan_test /global/cscratch1/sd/kaylanb/obiwan_out/randoms/elg_9deg2_ra175/more.fits


Now add `obiwan_test` to the end of the randoms table and delete `obiwan_test`::

    desi=> insert into obiwan_elg_ra175 select * from obiwan_test;
    desi=> drop table obiwan_test;

Make a QDO task list for your additional randoms. Specify `minid` to skip all primary keys below your 240,001::

    from obiwan.runmanager.qdo_tasks import TaskList
    T= TaskList(ra1=173.5,ra2=176.5, dec1=23.0,dec2=26.0,
                nobj_per_run=300,
                nobj_total=240000 + 720000)
    T.bricks()
    T.tasklist(do_skipid='no',do_more='yes',minid=240001)

then run these randoms from a new QDO queue::

    qdo create obiwan_ra175_domore
    qdo load obiwan_ra175_dormore tasks_skipid_no_more_yes_minid_240001.txt

Managing your qdo production run
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Manage your qdo production run with `obiwan/py/obiwan/runmanager/status.py`. To get a list of all log.<brickname> and slurm-<slurmid>.out files, sorted by status of "succeeded, failed, running" in the qdo db, and a tally of each error that occurred, do::

    cd $obiwan_code
    python $obiwan_code/obiwan/py/obiwan/runmanager/status.py --qdo_quename ${qdo_quename} --outdir /global/cscratch1/sd/kaylanb/obiwan_out/${name_for_run} --obj elg

Once you finish all the above runs, you will make a second qdo que. We previously made `obiwan_ra175` to do the usual obiwan runs, and now we make `obiwan_ra175_doskip` to do the randoms that we skipped. You can get a list of these tasks with::

    cat tasks_inregion | awk '{print $1}'|sort|uniq > brick_list.txt

::

    from obiwan.runmanager.qdo_tasks import write_qdo_tasks_skipids
    write_qdo_tasks_skipids('brick_list.txt', nobj_per_run=300)

which outputs a file `tasks_skipids.txt`. Now create a new qdo queue_name for the skipid runs and load the new tasks::

    export qdo_quename=obiwan_ra175_doskip
    qdo create ${qdo_quename}
    qdo load ${qdo_quename} tasks_skipids.txt


The qdo tasks automatically set the `do_skipid` flag, so you dont need to edit the `qdo_job_9deg.sh` file. Just run it with your new qdo `que_name`::

    cd $obiwan_out/${name_for_run}
    qdo launch ${qdo_quename} 40 --cores_per_worker 4 --batchqueue regular --walltime 05:00:00 --script $obiwan_code/obiwan/bin/qdo_job_9deg.sh --keep_env
