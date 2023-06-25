#!/bin/bash
# if possible remove .afg from $1?
# check for existence of -reads, -contig and -tle
export basename="$1"
cp ${basename}-reads.afg ${basename}.afg
sed \$d < ${basename}-contig.afg >> ${basename}.afg
tr ' ' '\n' < ${basename}-tle.afg >> ${basename}.afg
echo \} >> ${basename}.afg
echo "You can now view ${basename}.afg in tablet."
exit 0
