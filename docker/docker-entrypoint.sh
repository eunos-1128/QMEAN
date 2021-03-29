#!/usr/bin/bash

# Exit immediately on commands with a non-zero exit status.
set -eo pipefail


# In case there are no command line arguments, run Celery. (To be done)
# In case of command line arguments, check if it is QMEAN and run it
# In case of command line arguments but not QMEAN, just execute the command line
if [ $# -eq 0 ]; then
    : #set -- celery -A qm.clryschdlr worker --loglevel=info
else
    if [ $1 == "qmeandisco" ]; then
        set -- /qmean/run-qmeandisco.py "${@:2}"
    fi
fi

exec "$@"
