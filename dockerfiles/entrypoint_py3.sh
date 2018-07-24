#!/bin/bash

#stop immediately if any of these commands have an error
source /srv/py3_venv/bin/activate
exec "$@"
