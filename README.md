# BRAMS
THE BRAZILIAN DEVELOPMENTS ON THE REGIONAL ATMOSPHERIC MODELING SYSTEM

BRAMS IS A NUMERICAL MODELING SYSTEM DESIGNED FOR REGIONAL SCALE ATMOSPHERIC FORECASTING AND RESEARCH, WITH FOCUS ON ATMOSPHERIC CHEMISTRY, AIR QUALITY AND BIOGEOCHEMICAL CYCLES.
BRAMS is based on the Regional Atmospheric Modeling System (RAMS) originally developed at CSU/USA. BRAMS software is under a free license (CC-GPL). Currently, it is developed and maintained by CPTEC/INPE, USP and others institutions in Brazil and abroad. BRAMS/RAMS are multipurpose numerical weather prediction models that were designed to simulate atmospheric circulations spanning from hemispheric scales down to large eddy simulations (LES) of the planetary boundary layer. Several new functionalities and modifications were included to improve the numerical representation of key physical processes over tropical and subtropical regions (Freitas et al., 2005, 2009), a complete in-line module for atmospheric chemistry and aerosol processes (Longo et al., 2013), as well as a state-of-the-art model for simulation of atmosphere-surface exchange of water, energy, momentum, carbon and others biogeochemical tracers (Moreira et al., 2013).
The current version of BRAMS contains the following additional physical parameterizations:

Radiation: CARMA (Toon et al., 1988) and RRTMG (Iacono et al., 2008) schemes for long- and short-wave, including aerosols effects and coupled with microphysics and convection schemes
Microphysics: a double moment from RAMS CSU version, Thompson double moment and aerosol aware (Thompson and Eidhammer, 2014).
Convection schemes: Souza (1999) for shallow convection, Grell and Deveny (2002) for deep convection and Grell and Freitas (2014) scale and aerosol aware for deep and shallow convection including convective transport and wet removal of tracers.
Turbulence parameterizations: Nakanishi & Nino (2004) TKE based formulation, Taylor’s theory based formulation (Campos Velho, 1998).
Surface interaction and carbon cycle: TEB (Town Energy Budget) scheme to simulate urban areas (Freitas et al., 2007) and the Joint UK Land Environment Simulator (JULES) model (Moreira et al. 2013).
In-line emission, deposition, transport and chemistry model (CCATT) for atmospheric composition and air pollution studies (Longo et al., 2013).
Besides, BRAMS includes the following features:

Highly accurate and monotonic scheme (Walcek 2000, Freitas et al. 2012) for advection of scalars.
Complete, mass conservative formulation for the Exner function prognostic equation (Medvigy et al. 2005).
Computational Scalability: high efficiency up to ∼10000 cores with the MPI approach.
Digital filter for model initialization.
Model output at GrADS format during the model run time.
Coupling with the STILT Lagrangian Particle Dispersion Modelling.
Output to drive air parcels trajectory model calculation.
And ongoing work includes:

Fire spread model.
MATRIX (Bauer et al., 2008) aerosol model.
Data assimilation with GSI / 3D-VAR technique.
Time integration schemes (Runge-Kutta 3rd order, Adams-Bashforth-Moulton 3rd order, leapfrog with Robert-Asselin-Williams filter 2nd order).
Advection schemes (from 2nd to 6th order) with positivity and monotonicity constraints and WENO (Weighted Essentially Non-Oscillatory) formulation.
Code reformulation for fully compliant with the current FORTRAN standards.
