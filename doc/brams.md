project_dir: ./src
output_dir: ./html
project: BRAMS MODEL DOCUMENTATION
include: ../src/utils/include
project_github: https://projetos.cptec.inpe.br/projects/brams/repository
project_website:http://brams.cptec.inpe.br/
summary: Brazilian developments on the Regional Atmospheric Modeling System ![](http://brams.cptec.inpe.br/wp-content/uploads/2015/11/logo-brams.navigation.png "")
author: The Brams' Team
author_description: Brams is under Cathedral Model's Development. A small team based at CPTEC/INPE can do the commit of new versions in an SubVersion's (version control system). But all the code is free for modifications under CC-GPL License, except the 3rd parts (notice inside the code).
github: https://projetos.cptec.inpe.br/
email: brams-help@cptec.inpe.br
fpp_extensions: fpp
predocmark: >
media_dir: ./media
docmark_alt: #
predocmark_alt: <
display: public
         protected
         private
source: false
graph: true
search: true
macro: TEST
       LOGIC=.true.
extra_mods: json_module: http://jacobwilliams.github.io/json-fortran/
            futility: http://cmacmackin.github.io
license: CC-GPL
extra_filetypes: sh #

BRAMS is based on the Regional Atmospheric Modeling System (RAMS) originally developed at CSU/USA. BRAMS software is under a free license (CC-GPL). Currently, it is developed and maintained by CPTEC/INPE, USP and others institutions in Brazil and abroad. BRAMS/RAMS are multipurpose numerical weather prediction models that were designed to simulate atmospheric circulations spanning from hemispheric scales down to large eddy simulations (LES) of the planetary boundary layer. Several new functionalities and modifications were included to improve the numerical representation of key physical processes over tropical and subtropical regions (Freitas et al., 2005, 2009), a complete in-line module for atmospheric chemistry and aerosol processes (Longo et al., 2013), as well as a state-of-the-art model for simulation of atmosphere-surface exchange of water, energy, momentum, carbon and others biogeochemical tracers (Moreira et al., 2013).

The current version of BRAMS contains the following additional physical parameterizations:

    + Radiation:  CARMA (Toon et al., 1988) and RRTMG (Iacono et al., 2008) schemes for long- and short-wave, including aerosols effects and coupled with microphysics and convection schemes
    + Microphysics: a double moment from RAMS CSU version, Thompson double moment and aerosol aware (Thompson and Eidhammer, 2014).
    + Convection schemes:  Souza (1999) for shallow convection, Grell and Deveny (2002) for deep convection and Grell and Freitas (2014) scale and aerosol aware for deep and shallow convection including convective transport and wet removal of tracers.
    + Turbulence parameterizations: Nakanishi & Nino (2004) TKE based formulation, Taylor’s theory based formulation (Campos Velho, 1998).
    + Surface interaction and carbon cycle: TEB (Town Energy Budget) scheme to simulate urban areas (Freitas et al., 2007) and the Joint UK Land Environment Simulator (JULES) model (Moreira et al. 2013).
    + In-line emission, deposition, transport and chemistry model (CCATT) for atmospheric composition and air pollution studies (Longo et al., 2013).

Besides, BRAMS includes the following features:

    + Highly accurate and monotonic scheme (Walcek 2000, Freitas et al. 2012) for advection of scalars.
    + Complete, mass conservative formulation for the Exner function prognostic equation (Medvigy et al. 2005).
    + mComputational Scalability: high efficiency up to ∼10000 cores with the MPI approach.
    + Digital filter for model initialization.
    + Model output at GrADS format during the model run time.
    + Coupling with the STILT Lagrangian Particle Dispersion Modelling.
    + Output to drive air parcels trajectory model calculation.


@Note
** A BRAMS guide for Beginners (First Time users) - Version 5.2.5**<br/>
** How to Install ** <br/>
---
**1. Introduction**<br/>
This page will guide You to install software infrastructure and BRAMS model. BRAMS code is distributed under CC-GPL License. You can download, copy, modify and redistribute the code. Some rights are reserved. We recommend the user to read the license terms on site Creative Commons. BRAMS works with Linux or Unix operational’s systems. As a first approach, we recommend the Linux UBUNTU flavours distribution: UBUNTU <https://www.ubuntu.com/download/desktop?>, XUBUNTU <https://xubuntu.org/download/>, LUBUNTU <https://lubuntu.net/downloads/> or MATE <https://ubuntu-mate.org/download/>. <br/>
Please, read carefully each of the steps below in order to install and run the model.
---
**2. Install Fortran and C Compilers** <br/>
BRAMS has been tested with the compilers: INTEL compilers <https://software.intel.com/en-us/fortran-compilers/try-buy>, PGI compilers <https://www.pgroup.com/products/community.htm> and GNU Fortran compiler (GPL) <https://gcc.gnu.org/fortran/> or <https://askubuntu.com/questions/358907/how-do-i-install-gfortran>. Follow the instructions of each site to install the compilers.
---
**3. Install MPI Libraries and Software** <br/>
BRAMS works only in parallel mode. One can run the model using a single processor/core, but must to do it using MPI with the MPIRUN command. We recommend download and install the last version of MPICH stable release. The last working version is 3.2.1. You may have problems with new versions. Also take care to choose the correct version to your OS.

After download the MPICH in <https://www.mpich.org/downloads/>, please, proceed the installation:

    Uncompress: ~\> tar -zxvf mpich-3.1.4.tar.gz
    goto mpich directory: ~\> cd mpich-3.1.4
    Configure mpich’s makefile: ~\> ./configure -disable-fast CFLAGS=-O2 FFLAGS=-O2 CXXFLAGS=-O2 FCFLAGS=-O2 -prefix=/opt/mpich3 CC=icc FC=ifort F77=ifort
    make and install: ~\> make; ~\> sudo make install

Notes:

    -prefix= directory where mpich will be installed when You make the install command;
     CC= C compiler according You install on step (2) above;
    FC= Fortran compiler according You install on step (2) above;
    F77= Use the same as FC
---
**4. Download and compile the BRAMS Model** <br/>
To download BRAMS, you must fill a form with your identification. Go to BRAMS page in <http://brams.cptec.inpe.br/get-started/download-requisition/> fill the form and get the model.
After download the code, the instalation of BRAMS code is the following:

    extract files: ~\>tar -zxvf brams-5.2.5-src.tgz
    go to build directory: ~\>cd BRAMS/build/
    configure the model’s makefile: ~\> ./configure -program-prefix=BRAMS -prefix=/home/xuser -enable-jules -with-chem=RELACS_TUV -with-aer=SIMPLE -with-fpcomp=/opt/mpich3/bin/mpif90 -with-cpcomp=/opt/mpich3/bin/mpicc -with-fcomp=ifort -with-ccomp=icc
    make and make install: ~\>make; sudo make install

Notes:

    -prefix= directory where BRAMS will be installed when You make the install command;

    -with-chem= The chemical to be used (RELACS_TUV, RELACS_MX, CB07_LUT, CB07_TUV, RACM_TUV)

    -with-aer= aerosol mechanism (SIMPLE or MATRIX(under test))

    -with-fpcomp= and –with-cpcomp= mpich parallel compiler installed on step 3 above.

    -with-fcomp= and –with-ccomp= compiler installed on step 2 above.

@endnote

@Bug
If You find any bugs please, send an email to:

- <mailto:brams-help@cptec.inpe.br>

@endbug

**How to Run the BRAMS Model**:

Before run the model, You must get the input data. There are two packages for initial tests:

- Small meteorological case for laptops and desktops (722MB)Click to get it:
<ftp://ftp.cptec.inpe.br/brams/BRAMS/data/meteo-only.tgz>

- Small chemical case (using RELACS_TUV) for laptops and desktops (945 MB):
<ftp://ftp.cptec.inpe.br/brams/BRAMS/data/meteo-chem.tgz>

Both cases are just for testing and learning processes. To get data for any different runs, visit the Input Data page in <http://brams.cptec.inpe.br/input-data/>.

    After downloading a test case, uncompress it using the command “tar”: ~> tar -zxvf meteo-only.tgz
    Goto the test directory: ~>cd meteo-only
    Make a tmp directory (necessary to Jules): ~>mkdir ./tmp
    Export the tmp directory (especially if you are running NOT locally !): ~>export TMPDIR=./tmp
    Take care about stack size of your machine. Make it at least 65536 or unlimited : ~>ulimit -s 65536
    Run the model using mpirun installed on step 3: ~>/opt/mpich3/bin/mpirun -np 4 ./brams-5.2.5

Notes:

    -np= the numbers of cores You will use.
    If You do not know try the command: $ lscpu (see the information about CPU(s))

    ./brams-5.2.5 – The BRAM’s executable code. Please, point to the directory of code.

    When the model start a lot of logs will be printed on your screen, Pay attention for errors that will be printed during the run.

    If you got an error and needs for help, please, send a e-mail to brams_help@cptec.inpe.br and inform the error. Please, attach the log printed on the screen.

**Visualize the results with GrADS**:

BRAMS, by default, writes the output in subdirectory dataout/POSPROCESS. Later you can set others directories and features on namelist file (expert users). The format of the output in POS is in GrADS software (COLA/IGES). You must install GrADS in your computer.

Get the software from OpenGrads and install it on your computer, or, if You use UBUNTU You can get grads using: ~\> sudo apt-get install grads

After install, please:

    Goto pos directory: ~\> cd dataout/POSPROCESS
    Run grads software: ~\> grads -l
    When grads prompt appears on terminal (ga->), you can choose one of the output files, check they by listing: ga-> !ls -latr *.ctl
    Choose one and open it: ga->open METEO-ONLY-A-2015-08-27-030000-g1.ctl (the name is just an example)
    after the file is open You will see information about the file, LON, LAT, LEV, etc.
    list all the variables available in output: ga->q file
    Choose one of them and proceed the plot: ga->d tempc (in this example plotting tempc – temperature)

You can see more information about Grads on Cola: Grads User’s Guide
