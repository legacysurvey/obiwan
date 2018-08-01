#!/bin/bash

#stop immediately if any of these commands have an error
source /srv/py3_venv/bin/activate
export PYTHONPATH=/srv/repos_for_docker/obiwan/py:/srv/repos_for_docker/legacypipe/py:$PYTHONPATH
cd /
exec "$@"
