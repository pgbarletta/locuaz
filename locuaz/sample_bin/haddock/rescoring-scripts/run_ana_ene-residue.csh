#!/bin/csh
#
#Input pdbs are read from $filelist
#In $dispfile information on per residue interaction energies

setenv HADDOCK /home/abonvin/haddock2.0_devel
set cnsExec="/home/software/cns_solve_1.1.2/opteron-x86_64-linux/bin/cns"
set outfile="ana_ene-residue.out"

$cnsExec < ana_ene-residue.inp > $outfile
gzip $outfile
