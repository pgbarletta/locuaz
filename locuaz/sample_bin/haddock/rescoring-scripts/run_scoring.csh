#!/bin/csh
#
## scoring.inp uses the following cns input files:
#-calc_free-ene.cns
#-print_coorheader.cns
#-separate.cns
#-def_solv_param.cns
#-scale_intra_only.cns
#-scale_inter_only.cns
#-symmultimer.cns
#Parameters are read from $HADDOCK/toppar/
#Input pdbs are read from $filelist
#Output pdbs (with HADDOCK-header) are written to $Filenames.fileroot _$count
#In $dispfile some information is written from generate.inp etc.

setenv HADDOCK /m100_work/AIRC_Fortun21/SF/haddock/haddock2.1 
set cnsExec="/m100_work/AIRC_Fortun21/SF/haddock/cns_solve_1.3/ibm-ppc64le-linux/bin/cns"
set outfile="scoring.out"

source /m100_work/AIRC_Fortun21/SF/haddock/cns_solve_1.3/cns_solve_env
source /m100_work/AIRC_Fortun21/SF/haddock/haddock2.1/haddock_configure.csh

$cnsExec < scoring.inp > $outfile
