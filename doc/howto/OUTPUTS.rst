Data Model for Obiwan
======================

I completed a full obiwan production run on a 9 deg2 region, which is a box having ra = (173.5,176.5) and dec = (23.0,26.0). Below I describe what the outputs from Obiwan look like and how they are organized, for the :ref:`9 deg2 <elg-9deg2-ra175>` and :ref:`half of DR5 <elg-dr5>` production runs. 

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
* `skip_rs{0,300}/`: outputs from injecting first 300 randoms *that were skipped* because withthin 5 arcsec of another random in all `rs{0,300,600,...}/` directories above
* `more_rs{0,300,600,...}/`: same as rs*/ but for randoms that were added to the randoms DB *after* the rs*/ runs completed. These "more_" prefixes will appear whenever there are less randoms that 10x target density.
* `more_skip_rs{0,300,600,...}/`: same as skip_rs*/ but for the more_rs*/ directories.

All `*rs*/` directories have the same data model. For example,
`obiwan_out/elg_9deg2_ra175/elg/175/1757p240/rs0`
* `log.1757p240`: log file for this run
* `obiwan/`: obiwan metadata
* `coadd/, metrics/, tractor/, tractor-i/`: usual imaging pipeline outputs

Note, `coadd/` is the largest output so only the `rs0/coadd/` directory has model, chi2, etc. images. All other `*rs*/coadd/` directories do not have these; they only have the much smaller .jpg files.

Within the `*rs*/` directories, the data model is *identical* to the Data Releases except for one additional directory `obiwan/`. For example,
`obiwan_out/elg_9deg2_ra175/elg/175/1757p240/rs0/obiwan/`
* `metacat-elg-1757p240.fits`: info about how obiwan was run (e.g. brickname)
* `simcat-elg-1757p240.fits`: binary table with 1 row for each random that was injected, colums are: id, seed, ra, dec, grz flux, and other source properties 
* `skippedids-elg-1757p240.fits`: binary table with 1 row for each random that was skipped because within 5 arcsec of another. Only 1 column: id, which is how the random is looked up in the randoms PSQl db 
