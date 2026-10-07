#!/usr/bin/env bash
# run a list of sos_probe_sdp.py argument lines (stdin, one job per line) with at most $P parallel jobs
P=${P:-6}
xargs -P "$P" -I{} bash -c "$(dirname "$0")/sos_probe_run.sh {}"
