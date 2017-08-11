from six.moves import urllib
import tarfile
import os

def fetch_targz(remote_fn,local_fn):
    """grabs tar.gz file from web, puts on local dir and unpacks it

    Args:
      remote_fn: web link ending in .tar.gz
      local_fn: where to put it, ending in .tar.gz
    """
    local_dir= os.path.dirname(local_fn)
    if not os.path.isdir(local_dir):
        os.makedirs(local_dir)
    if not os.path.exists(local_fn):
        print('Grabbing: %s\n Putting here: %s' % (remote_fn,local_fn))
        urllib.request.urlretrieve(remote_fn, local_fn)
        tgz = tarfile.open(local_fn)
        tgz.extractall(path= local_dir)
        tgz.close()
    else:
        print('Already exists: %s' % (local_fn,))
