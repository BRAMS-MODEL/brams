! %W% %P% %G% %U%

! Code converted using TO_F90 by Alan Miller
! Date: 2011-09-15  Time: 16:28:24

!***********************************************************************
!***********************************************************************
! Development of this code was sponsored by EPRI, 3412 Hillview Ave.,  *
! Palo Alto, CA 94304 under Contract WO8221-01                         *
!                                                                      *
! Developed by Yang Zhang, Betty Pun, and Christian Seigneur,          *
! Atmospheric and Environmental Research, Inc., 2682 Bishop Drive,     *
! Suite 120, San Ramon, CA 94583                                       *

! Development of previously available modules are listed at the        *
! begining of the code of the corresponding module.  Some of these     *
! modules may be copyrighted                                           *
!***********************************************************************
! RCS file, release, date & time of last delta, author, state, [and locker]
! $Header: /usr/gchm/dwright/MODEL/MADRID/models/CCTM/src/aero/aero_MADRID2_qssa/ISRPIA.EXT,v 1.1 2002/10/24 17:03:55 models3 Exp $

!=======================================================================

! *** ISORROPIA PLUS CODE
! *** INCLUDE FILE 'ISRPIA.EXT'
! *** THIS FILE CONTAINS THE DECLARATIONS OF THE GLOBAL CONSTANTS
!     AND VARIABLES.

! *** COPYRIGHT 1996-98, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES

!=======================================================================
!***********************************************************************
!  REVISION HISTORY:                                                   *
!     Modified by YZ of AER to be used in MADRID2, Oct. 10, 2002       *
!                                                                      *
!***********************************************************************

IMPLICIT DOUBLE PRECISION (a-h, o-z)
INTEGER, PARAMETER :: ncomp=5, nions=7, ngasaq=3, nslds=9, npair=13, nzsr=100, nerrmx=25

! *** INPUT VARIABLES **************************************************

INTEGER :: metstbl
COMMON /inpt/ w(ncomp), waer(ncomp), temp, rh, iprob, metstbl

! *** WATER ACTIVITIES OF PURE SALT SOLUTIONS **************************

COMMON /zsr / awas(nzsr), awss(nzsr), awac(nzsr), awsc(nzsr),  &
    awan(nzsr), awsn(nzsr), awsb(nzsr), awab(nzsr), awsa(nzsr), awlc(nzsr)

! *** DELIQUESCENCE RELATIVE HUMIDITIES ********************************

INTEGER :: wftyp
COMMON /drh / drh2so4,  drnh42s4, drnahso4, drnacl,   drnano3,  &
    drna2so4, drnh4hs4, drlc,     drnh4no3, drnh4cl
COMMON /mdrh/ drmlcab,  drmlcas,  drmasan,  drmg1,    drmg2,  &
    drmg3,    drmh1,    drmh2,    drmi1,    drmi2,  &
    drmi3,    drmq1,    drmr1,    drmr2,    drmr3,  &
    drmr4,    drmr5,    drmr6,    drmr7,    drmr8,  &
    drmr9,    drmr10,   drmr11,   drmr12,   drmr13, wftyp

! *** VARIABLES FOR LIQUID AEROSOL PHASE *******************************

DOUBLE PRECISION :: molal, molalr, m0
REAL :: ionic
LOGICAL :: calaou, calain, frst, dryf
COMMON /ions/ molal(nions), molalr(npair), gama(npair), zz(npair),  &
    z(nions),     gamou(npair),  gamin(npair),m0(npair), gasaq(ngasaq),  &
    epsact,       coh,           chno3,       chcl,  &
    water,        ionic,         iacalc,  &
    frst,         calain,        calaou,      dryf

! *** VARIABLES FOR SOLID AEROSOL PHASE ********************************

COMMON /salt/ ch2so4,  cnh42s4, cnh4hs4, cnacl,   cna2so4,  &
    cnano3,  cnh4no3, cnh4cl,  cnahso4, clc

! *** VARIABLES FOR GAS PHASE ******************************************

COMMON /gas / gnh3, ghno3, ghcl

! *** EQUILIBRIUM CONSTANTS ********************************************

COMMON /equk/ xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10,  &
    xk11, xk12, xk13, xk14, xkw, xk21, xk22, xk31, xk32, xk41, xk42

! *** MOLECULAR WEIGHTS ************************************************

DOUBLE PRECISION :: imw
COMMON /othr/ r, imw(nions), wmw(ncomp), smw(npair)

! *** SOLUTION/INFO VARIABLES ******************************************

CHARACTER (LEN=15) :: scase
COMMON /case/ sulratw, sulrat, sodrat, scase

COMMON /soln/ eps, maxit, nsweep, ndiv, iclact

! *** ERROR SYSTEM *****************************************************

CHARACTER (LEN=40) :: errmsg
INTEGER :: errstk, nofer
LOGICAL :: stkofl
COMMON /eror/ stkofl, nofer, errstk(nerrmx), errmsg(nerrmx)

! *** GENERIC VARIABLES ************************************************

CHARACTER (LEN=14) :: version
COMMON /cgen/ great, tiny, tiny2, zero, one, version

! *** FOR INTERACTION WITH HYDROPHILIC ORGANIC COMPOUNDS
DOUBLE PRECISION :: organion, watorg
COMMON /org/ organion, watorg


! *** END OF INCLUDE FILE **********************************************

