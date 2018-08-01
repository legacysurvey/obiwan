**********************
Data Model
**********************

The data and directory structure for the :ref:`9 deg2 <elg-9deg2-ra175>` and :ref:`half of DR5 <elg-dr5>` production runs is described below.

Monte Carlo simulation files
------------------------------------------------
The root directory for all Monte Carlo simulation products is, on NERSC:

- **/global/cscratch1/sd/kaylanb/obiwan_out/**

Randoms
------------------------------------------------
The randoms are stored in FITS tables here:

- **.../randoms/elg_9deg2_ra175/elg_randoms/**

and in the PostreSQL DB

- **-h nerscdb03.nersc.gov -d desi -U desi_user**

Each random has a unique id, positon (Ra, DEC), and realistic properties (brightness, shape). Each position is a random draw from the unit sphere and each set of properties is a random draw from a Gaussian Mixture model.

Obiwan Outputs
------------------------------------------------
Results from running our :ref:`pipeline <https://github.com/legacysurvey/legacypipe>` on the modified images are here:

- **.../elg_9deg2_ra175/**

The top level directory has these files

* The usual six directories: *tractor,tractor-i,coadd,metrics,checkpoint,logs*

    * The coadd/ directory has a subset of the usual outputs because it is the directory with the largest file sizes

* Monte Carlo simulation metadata directory: *obiwan*
* slurm-[0-9]+.out: stdout, stderr from the SLURM job scripts

Subdirectories follow the usual Data Relase format of **.../bri/brick/**, where *bri* is the first three letters of each brick. The multiple iteractions per brick are identified by directories named **rs[0-9]+**, where the *rs* stands for the *Row* of the unique id-sorted table of randoms in the brick to *Start* from. The *[0-9]+* is the index of that row, which would be 0, 500, 1000, etc. when 500 fake galaxies are added to the images in each iteration.

For example, the tractor catalogues containing the first 1500 fake sources in brick *1757p240* are in

- .../tractor/175/1757p240/**rs0**/
- .../tractor/175/1757p240/**rs500**/
- .../tractor/175/1757p240/**rs1000**/

However, there will be additional directories to **rs[0-9]+**, such as **skip_rs[0-9]+**. These are the sources that were *skipped* because they were within 5 arcsec of another random. There are actually four sets of directories like this:

- **rs[0-9]+**: initial set of iterations per brick
- **skip_rs[0-9]+**: skipped sources from the initial set of iterations
- **more_rs[0-9]+**: new set of iterations per brick for randoms that were added to the database after the run completed (e.g. more randoms were needed)
- **more_skip_rs[0-9]+**: skipped sources from the new set of iterations

Simulation Metadata
------------------------------------------------
The metadata for each Monte Carlo simulation is stored in

- **.../obiwan/bri/brick/rs[0-9]+/**

It consists of four files:

- **metacat-elg-brick.fits**: how obiwan was run (e.g. brickname)
- **simcat-elg-brick.fits**: truth table of the final properties of the sources added to the images (e.g. guassian noise and galactic extinction make the fluxes actually added different from those in the DB)
- **skippedids-elg-brick.fits**: ids of the randoms that were skipped because the within 5 arcsec of another random
- **sim_ids_added.fits**: unique ids of the randoms that overlap with at least one CCD so are actually added to the images
