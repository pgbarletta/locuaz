#!/bin/csh
#
#Input pdbs are read from $filelist
#In $dispfile information on AIR energy and violation +is written .

setenv HADDOCK /home/abonvin/haddock2.0_devel
set cnsExec="/home/software/cns_solve_1.1.2/opteron-x86_64-linux/bin/cns"
set outfile="ana_airs.out"

$cnsExec < ana_airs.inp > $outfile
gzip $outfile
