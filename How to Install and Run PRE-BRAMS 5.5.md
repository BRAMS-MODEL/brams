# How to Install and Run PRE-BRAMS 5.5



1. Download the BRAMS 5.5 source from ftp:
   
   ```bash
   wget http://ftp.cptec.inpe.br/pesquisa/bramsrd/BRAMS/releases/stable/brams_v5.5.tar.xz
   ```

2. Unzip and unpack it
   
   ```bash
   xz -d brams_v5.5.tar.xz
   tar -xvf brams_v5.5.tar
   ```

3. Goto PRE_BRAMS build directory:
   
   ```bash
   cd BRAMS-5.5/auxProgs/PRE-BRAMS/build/
   ```

4. Edit the paths.mk file to put your paths:
   You must setup the paths for GRIB2 and NETCDF installation and the paths for BRAMS instalation. Do not change the FIXED VALUES. 
   
   ```makefile
   # USER RELATED
   CHEM_MECH=RELACS_TUV
   GRIB2_PATH=/lustre_xc50/luiz_flavio/models/BRAMS-5.5/auxProgs/PRE-BRAMS/install/grib2/
   NETCDFC_PATH=/opt/cray/pe/netcdf/4.4.1.1.6/gnu/51/
   NETCDFF_PATH=/opt/cray/pe/netcdf/4.4.1.1.6/gnu/51/
   #
   #
   #
   RAMS_ROOT=/lustre_xc50/luiz_flavio/models/BRAMS-5.5
   # FIXED VALUES
   SRC_PATH=../src
   CHEM_MODEL_PATH=$(RAMS_ROOT)/src/brams/ccatt/$(CHEM_MECH)
   DUMP_PATH=$(RAMS_ROOT)/src/utils/dump/
   BRAMS_INCLUDE=$(RAMS_ROOT)/src/utils/include
   ```
   
   <mark>If You have no grib2 or NetCDF installed see the section 11 and following.</mark>

5. compile the code:
   
   ```bash
   make OPT=gnu clean
   
   make OPT=gnu
   
   ```
   
   > The compilation is prepared only to gfortran. If You use another compiler copy the include file and creates another one for the new compiler.
   > 
   > If everything is ok a executable (<mark>prep_1.0</mark>) will be created in bin diretory on PRE-BRAMS bin directory

6. Goto bin directory and create a data directory for test
   
   ```bash
   cd ../bin
   mkdir data
   cd data
   ```
   
   > You must create the data directory in other directory if You prefer.

7. Download a set of data from NCEP/GFS for test:
   
   ```bash
   wget https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/gfs.20210203/00/gfs.t00z.pgrb2.0p25.f000
   
   wget https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/gfs.20210203/00/gfs.t00z.pgrb2.0p25.f006
   
   wget https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/gfs.20210203/00/gfs.t00z.pgrb2.0p25.f012
   
   wget https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/gfs.20210203/00/gfs.t00z.pgrb2.0p25.f018
   
   wget https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/gfs.20210203/00/gfs.t00z.pgrb2.0p25.f024
   
   ```
   
   > The example above is for date 2021 02 03, just one day of data.

8. Download the CAMS aerosol/chem climatology:
   
   if You will run the model for South America locations: 
   
   ```bash
   wget http://ftp.cptec.inpe.br/pesquisa/bramsrd/BRAMS/files/CAMS-SouthAmerica.tar
   ```
   
   if You will run the model for others locations:
   
   ```bash
   wget http://ftp.cptec.inpe.br/pesquisa/bramsrd/BRAMS/files/CAMS-Globe.tar
   ```
   
   > The climatology files must be downloaded just once. Keep it in some directory for future uses.

9. Setup the namelist file:
   
   go to your bin directory (`cd ..`)
   
   edit the **pre.nml** file
   
   ```bash
   $ARGS_INPUT
   
   !!!!! DATE !!!!!
   
   init_year  = 2021,
   init_month = 02
   init_day   = 03
   init_hour  = 0,
   final_year = 2021
   final_month= 02
   final_day = 04
   final_hour = 00,
   
   !!!!! TIME STEP !!!!!!!
   
   step  = 6, !Timestep in hours
   
   !!!!! ATMOS !!!!!
   atmos_type   = 1, !0=DP, 1=GFS Grib2
   atmos_prefix ='gfs.t00z.pgrb2.0p25.f',
   atmos_sufix  ='.2021020300.grib2',
   atmos_idir   ='./data/GFS_0p25/',
   levels = 23,
   initial_latitude = -70., !initial latitude for domain of model (-90 to 90)  
   final_latitude  = 20., !final latitude for domain of model(-90 to 90)
   initial_longitude = 250., !initial longitude for domain of model (0 to 360)
   final_longitude = 358., !Final longitude for domain of model (0 to 360)
   
   !!!!! CHEM !!!!!! 
   chem_type     = 1, !0 = no Chem, 1 = CAMS 
   chem_idir  = "./data/",
   chem1_prefix ='',
   chem1_sufix  ='-CAMS-EC-2010-2019-AMS',
   
   !!!!! OUTPUT !!!!!
   out_type   = 2, !0=text, 1=VFM, 2=Grads
   out_prefix = 'IC',
   out_sufix  = '',
   out_dir    = './',
   	
   $END
   ```
   
   > See that initial and final latitude and longitude must be greater then your domain area of interest. In the case above the area is over South America. For GFS input keep 23 levels. The output is pointed for Grads file. Use it! The model now works with bin file in input.

10. Run the code:
    
    ```bash
    ./prep_1.0
    ```
    
    > A bunch of information will be displayed in screen. If is all ok the IC/CC will be write in IC file at out_dir directory pointed in namelist.

11. Install WGRIB2 API package from NCEP
    
    In order to read grib2 (pattern data file), the model uses wgrib2 lib & API. 
    
    Create a install directory, download and install it.
    
    ```bash
    mkdir ~/~install
    cd ~/install
    wget https://www.ftp.cpc.ncep.noaa.gov/wd51we/wgrib2/wgrib2.tgz
    tar -zxvf wgrib2.tgz
    cd grib2
    ```
    
    Edit (use Your editor) the makefile, search the directives bellow and change accordingly the values:
    
    ```bash
    USE_NETCDF3=0
    USE_NETCDF4=0
    USE_REGEX=1
    USE_TIGGE=1
    USE_MYSQL=0
    USE_IPOLATES=3
    USE_SPECTRAL=0
    USE_UDF=0
    USE_OPENMP=0
    USE_PROJ4=0
    USE_WMO_VALIDATION=0
    DISABLE_TIMEZONE=0
    MAKE_FTN_API=1
    DISABLE_ALARM=0
    USE_G2CLIB=0
    USE_PNG=0
    USE_JASPER=0
    USE_AEC=0
    ```
    
    Now continue the instalation:
    
    ```bash
    make CC=gcc FC=gfortran
    make CC=gcc FC=gfortran lib
    cd ..
    sudo mv grib2 /opt
    ```
    
    > If You prefer or don't have  sudo permission You can use another directory to put the code. 

12. Install the NetCDF C libraries and packages
    
    The model use NetCDF Unidata file format to read and write files. Get the C NetCDF C library and files at Unidata [NetCDF Download site]([Unidata | Downloads | NetCDF](https://www.unidata.ucar.edu/downloads/netcdf/))
    
    
    Get the Stable Release of NETCDF-C Releases (Gziped format). After download it:
    
    
    ```bash
    cd ~/install
    mv ~/Downloads/netcdf-c-4.7.3.tar.gz .
    tar -xzvf netcdf-c-4.7.3.tar.gz
    cd netcdf-c-4.7.3
    ./configure --prefix=/opt/netcdfc FC=/opt/mpich3/bin/mpif90 CC=/opt/mpich3/bin/mpicc LT_SYS_LIBRARY_PATH=/opt/hdf5/lib CFLAGS='-I/opt/hdf5/include' LIBS='-L/opt/hdf5/lib'
    make
    sudo make install
    
    ```
    
    > The version 4.7.3 above is only an example. You must try newest if available.
    > Don't use a version < 4.7.3
    > 
    > You will need a MPI instaled
    
    

13. Install the NetCDF Fortran libraries and packages
    
    Get the Fortran NetCDF Fortran library and files at Unidata NetCDF Download site (the same above). Roll the page some lines.
    Get the Stable Release of NetCDF-Fortran Releases (Gziped format). After download it:
    
    ```bash
    
    cd ~/install
    mv ~/Downloads/netcdf-fortran-4.5.2.tar.gz .
    tar -xzvf netcdf-fortran-4.5.2.tar.gz
    cd netcdf-fortran-4.5.2
    export LD_LIBRARY_PATH=/opt/netcdfc/lib:/opt/hdf5/lib:$LD_LIBRARY_PATH
    export CC=/opt/mpich3/bin/mpicc
    export FC=/opt/mpich3/bin/mpif90
    export CPPFLAGS="-I/opt/netcdfc/include -I/opt/hdf5/include"
    export LDFLAGS="-L/opt/netcdfc/lib -L/opt/hdf5/lib"
    ./configure --prefix=/opt/netcdff
    make
    sudo make install
    ```
    
    
    
    


