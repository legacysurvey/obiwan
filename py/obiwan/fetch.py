"""Download and unpack tarballs"""

from six.moves import urllib
import tarfile
import os

def fetch_targz(targz_url, outdir):
    """downloads targz_url file to outdir, then untars it

    Args:
        targz_url: like 'http://www.google.com/path/to/data.tar.gz'
        outdir: where to download the tar.gz file to
    """
    name= os.path.basename(targz_url)
    outfn= os.path.join(outdir, name)
    assert('.tar.gz' in name)
    #local_dir= os.path.dirname(local_fn)
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    if not os.path.exists(outfn):
        print('Grabbing: %s\n Putting here: %s' % (targz_url,outfn))
        urllib.request.urlretrieve(targz_url, outfn)
        tgz = tarfile.open(outfn)
        tgz.extractall(path= outdir)
        tgz.close()
    else:
        print('Already exists: %s' % (outfn,))

    
