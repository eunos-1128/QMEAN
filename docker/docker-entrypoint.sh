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
        set -- /qmean/run_qmean.py --method QMEANDisCo "${@:2}"
    elif [ $1 == "qmean" ]; then
        set -- /qmean/run_qmean.py --method QMEAN "${@:2}"
    elif [ $1 == "qmeanbrane" ]; then
        set -- /qmean/run_qmean.py --method QMEANBrane "${@:2}"
    fi
fi

exec "$@"
