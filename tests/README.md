# Tests

### End to End 
There are three __datasets__ 
 1) DR5: DR5 zeropoints (which are the minimum set of columns needed for legacypipe, and computed by legacy zeropoints) and DR5 imaging
 2) DR3: DR3 zeropoints (many more columns than legacypipe need, and computed by idl zeropoints) and DR3 imaging
 3) `DR3_eBOSS`: DR3 zeropoints (many more columns than legacypipe need, and computed by idl zeropoints) but DR3 + additional eBOSS imaging

The "dataset" keyword tells obiwan which dataset to run on, e.g.
```sh 
python obiwan/kenobi.py --dataset {dr5,dr3,`dr3_eboss`}
```

When obiwan calss "runbrick" it has to specify the following options (unless runbrick is run from an older git checkout) to actually use each "dataset" as described above 
 1) DR5
 * --run dr5: use only survey-ccds files as named in DR5 and ignore "survey-ccds-kd" files
 2) DR3
 * --run dr3: use only survey-ccds files as named in DR3 and ignore "survey-ccds-kd" files
 * --no-rex --use-simp: turn off REX model, use SIMP instead 
 * --nsigma 6: default is 6 but enforce this 
 * __others__?
 3) `DR3_eBOSS`
 * --run `dr3_eboss`: use the survey-ccds files I made that exactly includes the `DR3_eBOSS` image list using idl zeropoints, also ignore "survey-ccds-kd" files
 * --no-rex --use-simp: turn off REX model, use SIMP instead 
 * --nsigma 6: default is 6 but enforce this 

These legacypipedir test dirs were made as follows
 * `obiwan/test/end_to_end/legacypipedir_1238p245_dataset_DR5`
 * `obiwan/test/end_to_end/legacypipedir_1238p245_dataset_DR3`
  - export brick=1238p245;export bri=`echo $brick|head -c 3`
  - read /global/project/projectdirs/cosmo/data/legacysurvey/dr3/coadd/${bri}/${brick}/legacysurvey-${brick}-ccds.fits, cut to 1 ccd (I chose the one with best seeing) and do `a.set('ccd_x1',a.ccd_x0+200)` and `a.set('ccd_y1',a.ccd_y0+200)`, then save to a new fit table "ccds.fits" 
  - `python legacyanalysis/create_testcase.py ccds.fits ${obiwan_code}/obiwan/py/obiwan/test/end_to_end/legacypipedir_${brick}_dataset_DR3 ${brick} --fpack`

### Run the end to end test problem for __obiwan/kenobi.py__
__DR5__
```sh
export dataset=DR5
export brick=1238p245
```

__DR3__
```sh
export dataset=DR3
export brick=1238p245
```

`__DR3_eBOSS__`
```sh
export dataset=DR3_eBOSS
export brick=TODO
```

now run it
```sh
export LEGACY_SURVEY_DIR=/global/cscratch1/sd/kaylanb/obiwan_code/obiwan/py/obiwan/test/end_to_end/legacypipedir_${brick}_dataset_${dataset}
python obiwan/kenobi.py --dataset ${dataset} -b ${brick} -n 2 -o elg --outdir ${dataset}_outdir --add_sim_noise
```
### Run an end test problem for __legacypipe/runbrick.py_

__DR5__
* full setup
 - `export LEGACY_SURVEY_DIR=/global/cscratch1/sd/kaylanb/obiwan_data/legacysurveydir_dr5;python legacypipe/runbrick.py --skip --no-write --threads 16 --skip-calibs --brick 1238p245 --outdir dr5_test --nsigma 6 --zoom 1100 1200 1100 1200 --run dr5`
* test problem
 - `export LEGACY_SURVEY_DIR=/global/cscratch1/sd/kaylanb/obiwan_code/obiwan/py/obiwan/test/end_to_end/test_legacypipedir_dr5ccds_1238p245;python legacypipe/runbrick.py --skip --no-write --skip-calibs --brick 1238p245 --outdir dr5_test --nsigma 6 --zoom 1 200 2900 3100 --run dr5`
__DR3__
* test problem

__DR3 eBOSS__
* can't use brick 1238p245 b/c not in eBOSS dr3 footprint

##### obiwan/kenobi.py
__DR5__
 - `export LEGACY_SURVEY_DIR=/global/cscratch1/sd/kaylanb/obiwan_code/obiwan/py/obiwan/test/end_to_end/test_legacypipedir_dr5ccds_1238p245;python legacypipe/runbrick.py --skip --no-write --skip-calibs --brick 1238p245 --outdir dr5_obiwan_test --nsigma 6 --zoom 1 200 2900 3100`




