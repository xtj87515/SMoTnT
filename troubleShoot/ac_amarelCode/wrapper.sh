#!/bin/bash
module use /projects/community/modulefiles/
module load gsl
./../source/SMOTNT "$@"
