#!/bin/bash
# Bytes of the given files/dirs currently in this node's page cache.
find "$@" -type f -print0 | xargs -0 fincore -b -n -o RES | awk '{s+=$1} END{print s+0}'
