# Description of Obiwan Outputs
I completed a full obiwan production run on a 9 deg2 region, which is a box having ra = (122.3, 125.3) and dec = (23, 26). Below I'm using the outputs from this production run to describe what the outputs from Obiwan look like and how they are organized. Everything is on NERSC is in this directory
`/global/cscratch1/sd/kaylanb/obiwan_out`
and has "chmod a+r" permissions so you can open and inspect the files yourself.

For info on how I created the files or ran the production run, see [PRODUCTION_RUN.md](https://github.com/legacysurvey/obiwan/blob/master/PRODUCTION_RUN.md)

### Randoms
`randoms/`
* `*kde.pickle`: saved Kernel Density Estimates (KDE) for g,r,z,redshift,rhalf,sersicn,etc.
* `elg_randoms/*.fits`: fits tables of random ra,dec from unit sphere and random draw from above KDE 

### Obiwan Outputs
obiwan_elg_9deg/
* `slurm-*.out`: top level logs from slurm job scripts
* `elg/bri/brick/`: Tractor Catalogues and Obiwan meta data for each brick
where bri is the first three letters of each brick.

For example, brick=1234p230 has outputs in
`obiwan_elg_9deg/elg/123/1238p245/`
* `rs{0,300,600,...}/`: outputs from injecting first 300 randoms, next 300, etc.
* `skip_rs{0,300}/`: outputs from injecting first 300 randoms **that were skipped** because withthin 5 arcsec of another random in all `rs{0,300,600,...}/` directories above

The `rs*/` and `skip_rs*/` directories have the same files, for example
`obiwan_elg_9deg/elg/123/1238p245/rs0/`
* `log.1238p245`: log file for this run
* `obiwan/`: obiwan metadata
* `coadd/, metrics/, tractor/, tractor-i/`: usual imaging pipeline outputs

Note, `coadd/` is the largest output so only the rs0/coadd/ directory is kept. All other `rs*` and all `skip_rs*` only have the .jpg files in their coadd/ directory

The new directory to the usual imaging pipeline is `obiwan/`, for example
`obiwan_elg_9deg/elg/123/1238p245/rs0/obiwan/`
* `metacat-elg-1238p245.fits`: info about how obiwan was run (e.g. brickname)
* `simcat-elg-1238p245.fits`: binary table with 1 row for each random that was injected, colums are: id, seed, ra, dec, grz flux, and other source properties 
* `skippedids-elg-1238p245.fits`: binary table with 1 row for each random that was skipped because within 5 arcsec of another. Only 1 column: id, which is how the random is looked up in the randoms PSQl db 
