#!/bin/bash

perl $WORKDIR/software/parseBAM.pl --input APOYTHmutmerge.bam --output $WORKDIR/matrix/APOYTHmut.matrix --verbose --removeDuplicates --minCoverage 3