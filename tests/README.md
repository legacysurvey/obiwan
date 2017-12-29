# Tests

The "testcase*" directories were made as follows
 * `obiwan/test/end_to_end/legacypipedir_1238p245_dataset_DR5`
 * `obiwan/test/end_to_end/legacypipedir_1238p245_dataset_DR3`
  - export brick=1238p245;export bri=`echo $brick|head -c 3`
  - read /global/project/projectdirs/cosmo/data/legacysurvey/dr3/coadd/${bri}/${brick}/legacysurvey-${brick}-ccds.fits, cut to 1 ccd (I chose the one with best seeing) and do `a.set('ccd_x1',a.ccd_x0+200)` and `a.set('ccd_y1',a.ccd_y0+200)`, then save to a new fit table "ccds.fits" 
  - `python legacyanalysis/create_testcase.py ccds.fits ${obiwan_code}/obiwan/py/obiwan/test/end_to_end/legacypipedir_${brick}_dataset_DR3 ${brick} --fpack`



