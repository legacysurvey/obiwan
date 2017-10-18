# Description of Obiwan Outputs
I completed a full obiwan production run on a 9 deg2 region, which is a box having ra = (173.5,176.5) and dec = (23.0,26.0). Below Im using the outputs from this production run to describe what the outputs from Obiwan look like and how they are organized. All outputs from the run are in this directory
`/global/cscratch1/sd/kaylanb/obiwan_out/elg_9deg2_ra175`
and has "chmod a+r" permissions so you can open and inspect the files yourself. The initial set of randoms are here `/global/cscratch1/sd/kaylanb/obiwan_out/randoms/elg_9deg2_ra175`, and also stored in a NERSC psql db.

For info on how I created the files or ran the production run, see [PRODUCTION_RUN.md](https://github.com/legacysurvey/obiwan/blob/master/PRODUCTION_RUN.md)

### Randoms
`randoms/elg_9deg2_ra175/`
* `elg_randoms/*.fits`: fits tables of random ra,dec from unit sphere and random draw from above KDE 

### Obiwan Outputs
`obiwan_out/elg_9deg2_ra175/`
* `slurm-*.out`: top level logs from slurm job scripts
* `elg/bri/brick/`: Tractor Catalogues and Obiwan meta data for each brick
where bri is the first three letters of each brick.

For example, brick=1757p240 has outputs in
`obiwan_out/elg_9deg2_ra175/elg/175/1757p240/`
* `rs{0,300,600,...}/`: outputs from injecting first 300 randoms, next 300, etc.
* `skip_rs{0,300}/`: outputs from injecting first 300 randoms **that were skipped** because withthin 5 arcsec of another random in all `rs{0,300,600,...}/` directories above

The `rs*/` and `skip_rs*/` directories have the same files, for example
`obiwan_out/elg_9deg2_ra175/elg/175/1757p240/rs0`
* `log.1757p240`: log file for this run
* `obiwan/`: obiwan metadata
* `coadd/, metrics/, tractor/, tractor-i/`: usual imaging pipeline outputs

Note, `coadd/` is the largest output so only the rs0/coadd/ directory is kept. All other `rs*` and all `skip_rs*` only have the .jpg files in their coadd/ directory

The new directory to the usual imaging pipeline is `obiwan/`, for example
`obiwan_out/elg_9deg2_ra175/elg/175/1757p240/rs0/obiwan/`
* `metacat-elg-1757p240.fits`: info about how obiwan was run (e.g. brickname)
* `simcat-elg-1757p240.fits`: binary table with 1 row for each random that was injected, colums are: id, seed, ra, dec, grz flux, and other source properties 
* `skippedids-elg-1757p240.fits`: binary table with 1 row for each random that was skipped because within 5 arcsec of another. Only 1 column: id, which is how the random is looked up in the randoms PSQl db 
