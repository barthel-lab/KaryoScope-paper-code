#!/usr/bin/env python3
"""Drop every file under the given paths from this node's page cache with
posix_fadvise(DONTNEED). Works on the NFS client without root."""
import os, sys
for p in sys.argv[1:]:
    walk = os.walk(p) if os.path.isdir(p) else [(os.path.dirname(p), [], [os.path.basename(p)])]
    for d, _, fs in walk:
        for f in fs:
            fd = os.open(os.path.join(d, f), os.O_RDONLY)
            os.posix_fadvise(fd, 0, 0, os.POSIX_FADV_DONTNEED)
            os.close(fd)
