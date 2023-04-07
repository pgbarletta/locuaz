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

setenv HADDOCK /home/abonvin/haddock2.1
set cnsExec="/home/abonvin/software/cns_solve_1.2/intel-x86_64bit-linux/bin/cns"
set outfile="scoring.out"

$cnsExec < scoring-air.inp > $outfile
gzip $outfile
