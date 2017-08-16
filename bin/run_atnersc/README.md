# Running Obiwan @ NERSC

The idea is for any NERSC user to easily do optimized production runs of Obiwan using a Docker Image. The steps are basically
 - unpack some tar.gz files to the user's scratch space
 - git clone https://github.com/kaylanb/obiwan.git
 - submit the included slurm job script that will load the Docker Image and run the obiwan repo

### Setup
For eBOSS DR3 do
```sh
$ export obiwan_dir=$CSCRATCH/obiwan
$ mkdir $obiwan_dir
$ cd $obiwan_dir
$ wget <path>/obiwan_eBOSS_DR3.tar.gz
$ tar -xzv obiwan_eBOSS_DR3.tar.gz
```
### Possible Production Runs
You can do the following runs with obiwan

| Run | Sources | Docs |
| ------ | ------ | ------ |
| eBOSS DR3 | ELG | fill in |
| DR5 | ELG | fill in |


