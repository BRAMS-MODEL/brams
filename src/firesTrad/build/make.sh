#!/bin/bash
gfortran -ffree-form -ffree-line-length-none -cpp -O3 -frecord-marker=4 -o dsa23bem.exe dsa23bem.f90
gfortran -ffree-form -ffree-line-length-none -cpp -O3 -frecord-marker=4 -o goes23bem.exe goes23bem.f90
gfortran -ffree-form -ffree-line-length-none -cpp -O3 -frecord-marker=4 -o rviirs23bem.exe rviirs23bem.f90
gfortran -m64 -ffree-form -ffree-line-length-none -cpp -O3 -frecord-marker=4 -o rmeteosat23bem.exe hdf5utils.f90 rmeteosat23bem.f90 -I/home/luiz.flavio/apps/include -L/home/luiz.flavio/apps/lib/ -lhdf5_fortran -lhdf5hl_fortran 



#####/opt/hdf5/lib/libhdf5hl_fortran.a /opt/hdf5/lib/libhdf5_hl.a /opt/hdf5/lib/libhdf5_fortran.a /opt/hdf5/lib/libhdf5.a -lz -lrt -ldl -lm -Wl,-rpath -Wl,/opt/hdf5/lib/ #-L/mnt/d/Ariane/LIBS/jpeg9/lib -ljpeg -L/mnt/d/Ariane/LIBS/zlib-1.2.4/lib -lz


################# no used yet ##########
#gfortran -m64 -ffree-form -ffree-line-length-none -cpp -O3 -frecord-marker=4 -o rmodis23bem.exe hdf4utils.f90 rmodis23bem.f90 -L/lustre_xc50/luiz_flavio/hdf4/lib -lmfhdf -ldf #-L/mnt/d/Ariane/LIBS/jpeg6/lib -ljpeg -L/mnt/d/Ariane/LIBS/zlib-1.2.4/lib -lz
mv *.exe ../bin  		
