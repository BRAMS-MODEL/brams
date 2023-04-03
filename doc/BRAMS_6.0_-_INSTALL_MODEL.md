# BRAMS 6.0 - INSTALL GUIDE

Before You install the model BRAMS-6.0  You must have the prerequisites installed. See the document  [BRAMS_6.0_-_INSTALL_PREREQUISITES](http://ftp.cptec.inpe.br/pesquisa/bramsrd/BRAMS-6.0/docs/BRAMS_6.0_-_INSTALL_PREREQUISITES.html) (extension md or html) to see how to install the requisites. **Even You believe the system are ready, please, read the prerequisites document.**

1. Building PATHS and linking compilers.
   
   In order to compiler the BRAMS code using the compilers and libraries You install in prerequisites make the commands bellow.
   
   ```batch
   export PATH={YOUR_DIR}/bin:$PATH
   export LD_LIBRARY_PATH={YOUR_DIR}/lib:$LD_LIBRARY_PATH
   sudo ln -s /usr/bin/gfortran {YOUR_DIR}/bin/gfortran
   sudo ln -s /usr/bin/gcc {YOUR_DIR}/bin/gcc
   ```
   
   > Notice: See the prerequisites document to use the correct {YOUR_DIR} . <mark>Pay attention in correct compiler You are using!!!!</mark>
   
   Please, check if Your path have in first part **{YOUR_DIR}** and if the correct library is in first part of LD_LIBRARY_PATH.  
   
   If You have doubt please, read the prerequisites' doc and use the alias You build (see 9. in prerequisites)
   
   ```batch
   echo $PATH
   echo $LD_LIBRARY_PATH
   ```
   
   Please check if Fortran and C version is the correct. (The example is only for Gnu. You must change if use another compiler) 
   
   ```
   gfortran --version
   gcc --version
   ```
   
   The results must be something like. See that in this case we use 8.4.0.
   
   ```
   $ gfortran --version
   GNU Fortran (Ubuntu 8.4.0-4ubuntu1) 8.4.0
   Copyright (C) 2018 Free Software Foundation, Inc.
   This is free software; see the source for copying conditions.  There is NO
   warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
   $ gcc --version
   gcc (Ubuntu 8.4.0-4ubuntu1) 8.4.0
   Copyright (C) 2018 Free Software Foundation, Inc.
   This is free software; see the source for copying conditions.  There is NO
   warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
   ```

2. ## Configure the model
   
   Now You must configure the model passing all libreries You will be use. The **{YOUR_BIN_AREA}** is the folder you want to put the binary and namelists of model. 
   
   ```bash
   ./configure --program-prefix=BRAMS_6.0 --prefix={YOUR_BIN_AREA} --enable-jules    --with-chem=RELACS_TUV --with-aer=SIMPLE --with-fpcomp={YOUR_DIR}/bin/mpif90    --with-cpcomp={YOUR_DIR}/bin/mpicc --with-fcomp={your_ortran_compiler} --with-ccomp={your_C_compiler} --with-netcdff={YOUR_DIR} --with-netcdfc={YOUR_DIR} --with-wgrib2={YOUR_DIR}
   ```

Bellow an example of use using Gnu, where {YOUR_BIN_AREA} is /home/oscar.

```bash
./configure --program-prefix=BRAMS_6.0 --prefix=/home/oscar --enable-jules    --with-chem=RELACS_TUV --with-aer=SIMPLE --with-fpcomp=/opt/gnu8/bin/mpif90    --with-cpcomp=/opt/gnu8/bin/mpicc --with-fcomp=gfortran --with-ccomp=gcc --with-netcdff=/opt/gnu8 --with-netcdfc=/opt/gnu8 --with-wgrib2=/opt/gnu8
```

You can configure the model without GRIB2 library, in this case don't use the `--with-wgrib2` directive. Remember that without this option You can not read grib2 files, but You can use for test purposes.

3. ## Make and Make install
   
   The make command will create the brams-6.0 executable. After creation it is necessary to run the make install command so that the basic files for the run are copied to the area set in {YOUR_BIN_AREA}
   
   ```bash
   make
   make install
   ```
   
   The output of BRAMS's model is presented with colors in terminal. If You save the output log (by ">&" use) may be a problem if You want to edit the file. Is possibles to make a filter to extract all special symbols used to color the output. To make the filter just
   
   ```bash
   make filter
   make install-filter
   ```
   
   BRAM's models needs a lot of input files. One of them is the Initial and boundary conditions (IC). This files is created from a global model data. Two diferents data input must be created, from NCEP model GFS or from ERA5 data. To create the IC data You must create the PRE-BRAMS utility.
   
   ```bash
   make pre-brams
   make install-pre-brams
   ```

4. ## Download tables, fixed files an a small test case
   
   A series of fixed files and tables are required to run the model. You can get them by downloading the compressed file available in the ftp area. This file has 16GB and may take a while to download. It depends on your network speed. We advise you to check the checksum of the file to make sure that the download has not broken it.
   
   ```bash
   cd {YOUR_BIN_AREA}
   cd ..
   wget http://ftp.cptec.inpe.br/pesquisa/bramsrd/BRAMS-6.0/test_set_data/MD5SUM
   wget http://ftp.cptec.inpe.br/pesquisa/bramsrd/BRAMS-6.0/test_set_data/brams6.0_test_bin.tar.xz
   md5sum brams6.0_test_bin.tar.xz
   cat MD5SUM
   tar -xvf brams6.0_test_bin.tar.xz
   ```

5. ## Download tables, fixed files an a Big test case
   
   The big test case is complete. It run about all South America  with 8km of resolution, 1017x993 points, 45 levels, and simulate by 216 hours. This case do not run in small machines and need systems with big numbers of processors/cores.
   
   ```bash
   cd {your bin area}
   wget http://ftp.cptec.inpe.br/pesquisa/bramsrd/toEgeon/bm/bm.tar.bz2
   bunzip2 bm.tar.bz2
   tar -xvf bm.tar
   ```

6. ## Intall GRADS file
   
   The output of model as You read above is in GRADS file. May be used another formats but we will show how in expert users guide. If You are using Ubuntu Linux or some Debia derivated Linux You can install grads in a simple way, using apt-get command
   
   ```bash
   sudo apt-get install grads
   ```
   
   To install grads right from source, please, see the information on site [GrADS Downloads](http://cola.gmu.edu/grads/downloads.php). Pay attention on the necessary libraries explained on [GrADS Supplibs](http://cola.gmu.edu/grads/gadoc/supplibs2.html)
   
   To learn about grads we recommend to read the document [Grads Manual from NCEP](https://www.cpc.ncep.noaa.gov/products/international/grads/Advanced_GrADS_Manual.pdf).
   
    
