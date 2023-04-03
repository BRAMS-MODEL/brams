#!/bin/sh
cd .
##### NEW
cat << Eof > gf.inp
 &run
  runname   = 'gf-3d-clean',  
  runlabel  = 'gf-3d-clean',  
  version   = 1,
 &end
Eof

mk ; gf.x > new.txt
echo " compare "
cmp gf-3d-ref.gra gf-3d-clean.gra

exit
##### OLD 
cat << Eof > gf.inp
 &run
  runname   = 'gf-1d-ok',  
  runlabel  = 'gf-1d-ok',  
  version   = 2,
 &end
Eof
gf.x > old.txt
