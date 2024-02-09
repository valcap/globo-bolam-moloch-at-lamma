!    Pochva:
!    Numerical scheme for the parametrization of thermo-hydrous processes
!    in soil and at its surface (own soil, vegetation, snow cover)
!
! Scheme is elaborated by O.Drofa, CNR-ISAC, Italy (o.drofa-at-isac.cnr.it)
!
! ---------------------------------------------------------------------------------------------------------------------------------
MODULE POCHVA_SCHEME

USE MODULE_POCHVA
IMPLICIT NONE

CONTAINS

! ---------------------------------------------------------------------------------------------------------------------------------

SUBROUTINE POCHVA(LEV_SOIL_EXT, TS_EXT, T_SURF_EXT, QS_EXT, QV_SURF_EXT, &
 SNOW, FSNOW, SNOW_ALBEDO, &
 FLUX_HEAT_BOTTOM, FLUX_WATER_BOTTOM, RUN_OFF_EXT, RUN_OFF_TOT_EXT, &
 MASK_SOIL_EXT, H_ATM_BOTT_LEV, PSURF, PATM_BOTT_LEV, TATM_BOTT_LEV, QATM_BOTT_LEV, &
 FLUX_RAIN, FLUX_SNOW, FLUX_RAD_TOT, FLUX_RAD_VIS, FLUX_HEAT_SPEC, FLUX_HEAT_LAT, FLUX_Q, KTURBH, KTURBWV, &
 QSOIL_MAX, QSOIL_MIN, CSOIL, RHOSOIL, PSISOIL, WCONDSOIL, PARBSOIL, PARCSOIL, &
 QSOIL_REL_WILT, QSOIL_REL_REF, LAIVEG, LAIVEG_MAX, FRACVEG, ROOTVEG, IND_LEV_TBOTTOM, IND_LEV_QBOTTOM, &
 VALMISS_EXT, DTIME, STEP_NUM, FLAG_TURB_FLUX, FLAG_INI, FLAG_PAR_CHANGE, PROC_IND_EXT, &
 SNOW_DIRT_EXT, TS_SURF_EXT, QS_SURF_EXT, FICE_SOIL_EXT, FICE_SOIL_SURF_EXT, &
 FLAG_SNOW_INIT, LEV_SNOW_EXT, TSNOW_EXT, FICE_SNOW_EXT, SNOW_AGE_EXT, SNOW_MELT_AGE_EXT, RHO_SNOW_EXT)

! Main procedure of Pochva scheme

! --- Input/output description ----------------------------------------------------------------------------------------------------

! Input and output variables:
!      TS_EXT : temperature at the levels inside of the soil - main prognostic variable (K), 1-st level is the upper inside, 3D array
!      QS_EXT : specific volumetric water content at the levels inside of the soil - main prognostic variable (m^3/m^3), 1-st level is the upper inside, 3D array
!      T_SURF_EXT : temperature at the surface - main prognostic variable (K), 2D array
!      SNOW : water content (kg/m^2) of snow cover, 2D array

! Input variables:
!      NPOINT_X and NPOINT_Y : number of grid point along horizontal axises of model domain
!      LEV_SOIL_EXT : depth of each levels inside of the soil (m), at this levels main
!                     variables (temperature and water content) are definied, 1D array
!      MASK_SOIL_EXT : Flag for application of the Pochva scheme (1 or 0), 2D array (integer)
!      H_ATM_BOTT_LEV : altitude of the bottom atmospheric level (m), 2D array
!      PSURF : pressure at the surface (Pa), 2D array
!      PATM_BOTT_LEV : pressure at the bottom atmospheric level (Pa), 2D array
!      TATM_BOTT_LEV : temperature at the bottom atmospheric level (K), 2D array
!      QATM_BOTT_LEV : specific humidity at the bottom atmospheric level (kg/kg), 2D array
!      FLUX_RAIN     : flux of liquid precipitation at the surface (kg/m^2/s), 2D array
!      FLUX_SNOW     : flux of solid (ice) precipitation at the surface (kg/m^2/s), 2D array
!      FLUX_RAD_TOT  : flux of total radiation at the surface (Joule/m^2/s), 2D array
!      FLUX_RAD_VIS  : flux of visible radiation at the surface (Joule/m^2/s), 2D array
!      KTURBH and KTURBWV : coefficients of turbulent exchange of heat and water vapour
!                           in the atmosphere surface layer (m^2/s), 2D array
!      DTIME : time step (s), real variable
!      STEP_NUM : number of current time step, integer variable
!      FLAG_TURB_FLUX : flag of status of turbulent fluxes (3 fluxes: specific and latent heat, water vapour)
!                       between the surface anf the lowset atmospheric lower,
!                       if FLAG_TURB_FLUX=0, then this fluxes will be defined by Pochva scheme 
!                       as output fields using input value of coefficients of turbulent exchange,
!                       if FLAG_TURB_FLUX=1, then this fluxes must be defined as input fields
!                       and will be used by Pochva scheme without changing of them;  
!      FLAG_INI : flag for initialization procedure application, integer variable: 1 - true, not 1 - false
!      FLAG_PAR_CHANGE : flag for long (climate) run application, it permits to change LAI and cover faction of vegetation,
!                        and, to change the status of grid point for the case of thiсk sea ice,
!                        thick sea ice point is sumulated as a glacier point, and thus points may appear and disappear,
!                        integer variable: 0 - false, not 0 - true
! Physical pamameters of the soil:
!      QSOIL_MAX      : Maximum specific volumetric soil water content (m^3/m^3) at all soil levels, 1-st level is the upper inside, 3D array
!      QSOIL_MIN      : Minimum specific volumetric soil water content (m^3/m^3) at all soil levels, 1-st level is the upper inside, 3D array
!      CSOIL          : Dry soil thermal capacity (J/K/m^3) at all soil levels, 1-st level is the upper inside, 3D array
!      RHOSOIL        : Dry soil density (kg/m^3) at all soil levels, 1-st level is the upper inside, 3D array
!      PSISOIL        : Idraulic potential of saturated soil (m) at all soil levels, 1-st level is the upper inside, 3D array
!      WCONDSOIL      : Idraulic conductivity of saturated soil (m/s) at all soil levels, 1-st level is the upper inside, 3D array
!      PARBSOIL       : Parameter in soil water transport equation (dimensionless) at all soil levels, 1-st level is the upper inside, 3D array
!      PARCSOIL       : Parameter in soil water transport equation (dimensionless) at all soil levels, 1-st level is the upper inside, 3D array
!      QSOIL_REL_WILT :  Relative soil water content (proportion) at which the vegetation wilt at all soil levels, 1-st level is the upper inside, 3D array
!      QSOIL_REL_REF  : Relative soil water content (proportion) at which the vegetation stop evapotraspiration increase at all soil levels, 1-st level is the upper inside, 3D array
!      LAIVEG         : Leaf Area Index of vegetation (dimensionless), 2D array
!      LAIVEG_MAX     : LAI maximum value permited for used vegetation dataset (refer value), real variable
!      FRACVEG        : Fraction of vegetation (proportion), 2D array
!      ROOTVEG        : Root depth of vegetation (m), 2D array
!      VALMISS_EXT    : Flag-value (soil temperature or soil water content) to defined
!                       soil bottom level, real variable
!      IND_LEV_TBOTTOM : index of bottom (climate) soil level for heat trasport at each grid point, 2D array
!      IND_LEV_QBOTTOM : index of bottom soil level for water trasport at each grid point, 2D array

! Output variables:
!      QV_SURF_EXT : specific humidfity at the surface - main prognostic variable (kg/kg), 2D array
!      FSNOW : fraction of grib box covered by snow, 2D array
!      SNOW_ALBEDO : albedo of snow, 2D array
!      FLUX_HEAT_BOTTOM : Heat flux (J/s/m/m) at soil bottom level, 2D array
!      FLUX_WATER_BOTTOM : Water flux (kg/s/m/m) at soil bottom level, 2D array
!      RUN_OFF_EXT : Water flux (kg/s/m/m) of soil run-off in 1,5 m surface layer, 2D array
!      RUN_OFF_TOT_EXT : Water flux (kg/s/m/m) of soil run-off, 2D array

! Input or output variables (dipending on FLAG_TURB_FLUX, see upper):
!      FLUX_HEAT_SPEC : Surface specific heat flux (J/s/m/m), 2D array
!      FLUX_HEAT_LAT  : Surface latent heat flux (J/s/m/m), 2D array
!      FLUX_Q         : Water vapour surface flux (kg/s/m/m), 2D array

! Input optional variables:
!      SNOW_DIRT_EXT  : snow_dirt: snow "dirtibility", that is coefficient
! (proportion 0-1) of dirt growth of snow surface dirty due to vegetation
! waste, aerosol deposition, etc., it is used in snow albedo definition
!      FLAG_SNOW_INIT : flag of status of snow variables (levels depth, temperature, ice faction, age, period of melting)
!                       if FLAG_SNOW_INIT=0, then snow variables are defined in input and will be used in initializition of Pochva scheme,
!                       if FLAG_SNOW_INIT=1, then snow variables are not defined in input and will be defined by Pochva scheme as output only
!                       Attention: it is recommended to set FLAG_SNOW_INIT=0

! Input and output optional variables:
!      TS_SURF_EXT  : temperature (K) at the soil top, in the case of snow presence it is not surface temperature, but is temperature at snow-soil boundary, 2D array
!      QS_SURF_EXT  : specific volumetric water content at the soil top (m^3/m^3), 2D array
!      FICE_SOIL_EXT: fraction of ice phase in soil water content (proportion) at the levels inside of the soil, 1-st level is the upper inside, 3D array
!      FICE_SOIL_SURF_EXT: fraction of ice phase in soil water content (proportion) at the soil top, 2D array

! Input and output or output only optional variables dipending on FLAG_SNOW_INIT (see over):
!      LEV_SNOW_EXT : depth of snow levels (kg/m/m), must have third array dimension equal 11, 
!                     1-st level is the snow top, 11-th level is surface between snow and soil, 3D array 
!      TSNOW_EXT    : snow temperature (K), must have third array dimension equal 11, 
!                     1-st level is the snow top, 11-th level is surface between snow and soil, 3D array 
!      FICE_SNOW_EXT: fraction of ice phase in snow (proportion) at the snow levels, must have third array dimension equal 11, 
!                     1-st level is the snow top, 11-th level is surface between snow and soil, 3D array 
!      SNOW_AGE_EXT : snow age (total cumulated period of existence, day), must have third array dimension equal 11, 
!                     1-st level is the snow top, 11-th level is surface between snow and soil, 3D array 
!      SNOW_MELT_AGE_EXT : melting snow age (total cumulated period when snow was underwent melting, day), must have third array dimension equal 11, 
!                     1-st level is the snow top, 11-th level is surface between snow and soil, 3D array 

! Output optional variables:
!      RHO_SNOW_EXT : snow density (kg/m/m/m) - diagnostic variable, must have third array dimension equal 11, 
!                     1-st level is the snow top, 11-th level is surface between snow and soil, 3D array 

! ---------------------------------------------------------------------------------------------------------------------------------

REAL, DIMENSION(:), INTENT(in) :: LEV_SOIL_EXT
REAL, DIMENSION(:,:,:), INTENT(inout) :: TS_EXT, QS_EXT
REAL, DIMENSION(:,:,:), INTENT(in) :: QSOIL_MAX, QSOIL_MIN, CSOIL, RHOSOIL, PSISOIL, WCONDSOIL, PARBSOIL, PARCSOIL,&
 QSOIL_REL_WILT, QSOIL_REL_REF
INTEGER, DIMENSION(:,:), INTENT(in) :: MASK_SOIL_EXT
REAL, DIMENSION(:,:), INTENT(inout) :: T_SURF_EXT, QV_SURF_EXT, SNOW
REAL, DIMENSION(:,:), INTENT(out) :: FSNOW, SNOW_ALBEDO, FLUX_HEAT_BOTTOM, FLUX_WATER_BOTTOM, RUN_OFF_EXT, RUN_OFF_TOT_EXT
REAL, DIMENSION(:,:), INTENT(in) :: H_ATM_BOTT_LEV,&
 PSURF, PATM_BOTT_LEV, TATM_BOTT_LEV, QATM_BOTT_LEV,&
 FLUX_RAIN, FLUX_SNOW, FLUX_RAD_TOT, FLUX_RAD_VIS, KTURBH, KTURBWV,&
 LAIVEG, FRACVEG, ROOTVEG
INTEGER, DIMENSION(:,:), INTENT(in) :: IND_LEV_TBOTTOM, IND_LEV_QBOTTOM
REAL, DIMENSION(:,:) :: FLUX_HEAT_SPEC, FLUX_HEAT_LAT, FLUX_Q
INTEGER, INTENT(in) :: FLAG_TURB_FLUX, FLAG_INI, STEP_NUM, PROC_IND_EXT
INTEGER, INTENT(inout) :: FLAG_PAR_CHANGE
REAL, INTENT(in) :: VALMISS_EXT, DTIME, LAIVEG_MAX

REAL, DIMENSION(:,:), INTENT(in), OPTIONAL :: SNOW_DIRT_EXT
REAL, DIMENSION(:,:), INTENT(inout), OPTIONAL :: TS_SURF_EXT, QS_SURF_EXT, FICE_SOIL_SURF_EXT
REAL, DIMENSION(:,:,:), INTENT(inout), OPTIONAL :: FICE_SOIL_EXT
INTEGER, INTENT(in), OPTIONAL :: FLAG_SNOW_INIT
REAL, DIMENSION(:,:,:), INTENT(inout), OPTIONAL :: LEV_SNOW_EXT, TSNOW_EXT, &
 FICE_SNOW_EXT, SNOW_AGE_EXT, SNOW_MELT_AGE_EXT
REAL, DIMENSION(:,:,:), INTENT(out), OPTIONAL :: RHO_SNOW_EXT

INTEGER :: NPOINT_X_EXT, NPOINT_Y_EXT, NLEV_SOIL_EXT, NLEV_SNOW_EXT, &
 NLEV, I, J, K, IUNIT
INTEGER, DIMENSION(3) :: SHAPE_3D
REAL :: SNALBED

! For testing output "test" --->
!INTEGER :: N_OUTPOINT=0, N_OUTPOINT2=7
!INTEGER, DIMENSION(10) :: &
! I_OUTPOINT=    (/ 25,   24,    9,   34,    6,    8,   13,   -1,   -1,   -1/), &
! J_OUTPOINT=    (/ 51,   16,   16,   38,   27,   20,   48,   -1,   -1,   -1/), &
! IPROC_OUTPOINT=(/ 59,   52,   66,   61,  111,   71,   40,   -1,   -1,   -1/)
INTEGER :: N_OUTPOINT=0, N_OUTPOINT2=0
INTEGER, DIMENSION(10) :: &
 I_OUTPOINT=    (/  8,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1/), &
 J_OUTPOINT=    (/ 42,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1/), &
 IPROC_OUTPOINT=(/ 61,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1/)
CHARACTER (LEN=30) :: FILE_OUT
! <---

 SHAPE_3D=SHAPE(TS_EXT)

 NPOINT_X_EXT=SHAPE_3D(1)
 NPOINT_Y_EXT=SHAPE_3D(2)
 NLEV_SOIL_EXT=SHAPE_3D(3)

 PROC_IND=PROC_IND_EXT

! Initialization of all variables using in the Pochva scheme

 IF (FLAG_INI == 1) THEN
  
   IF (PRESENT(FICE_SOIL_EXT)) THEN
     SHAPE_3D=SHAPE(FICE_SOIL_EXT)
     NLEV=SHAPE_3D(3)
     IF (NLEV /= NLEV_SOIL_EXT) THEN
       PRINT *,'In POCHVA scheme input arrays fice_soil_ext do not have correct dimension for 3 index,'
       PRINT *,' that must be 1:',NLEV_SOIL_EXT,', stop'
       STOP
     ENDIF
   ENDIF
  
   IF (PRESENT(LEV_SNOW_EXT)) THEN
     SHAPE_3D=SHAPE(LEV_SNOW_EXT)
     NLEV_SNOW_EXT=SHAPE_3D(3)
     IF (NLEV_SNOW_EXT /= NLEV_SNOW+1) THEN
       PRINT *,'In POCHVA scheme input arrays lev_snow_ext do not have correct dimension for 3 index,'
       PRINT *,' that must be 1:',NLEV_SNOW+1,', stop'
       STOP
     ENDIF
   ENDIF
  
   IF (PRESENT(TSNOW_EXT)) THEN
     SHAPE_3D=SHAPE(TSNOW_EXT)
     NLEV=SHAPE_3D(3)
     IF (NLEV /= NLEV_SNOW+1) THEN
       PRINT *,'In POCHVA scheme input arrays tsnow_ext do not have correct dimension for 3 index,'
       PRINT *,' that must be 1:',NLEV_SNOW+1,', stop'
       STOP
     ENDIF
   ENDIF
  
   IF (PRESENT(FICE_SNOW_EXT)) THEN
     SHAPE_3D=SHAPE(FICE_SNOW_EXT)
     NLEV=SHAPE_3D(3)
     IF (NLEV /= NLEV_SNOW+1) THEN
       PRINT *,'In POCHVA scheme input arrays fice_snow_ext do not have correct dimension for 3 index,'
       PRINT *,' that must be 1:',NLEV_SNOW+1,', stop'
       STOP
     ENDIF
   ENDIF
  
   IF (PRESENT(SNOW_AGE_EXT)) THEN
     SHAPE_3D=SHAPE(SNOW_AGE_EXT)
     NLEV=SHAPE_3D(3)
     IF (NLEV /= NLEV_SNOW+1) THEN
       PRINT *,'In POCHVA scheme input arrays snow_age_ext do not have correct dimension for 3 index,'
       PRINT *,' that must be 1:',NLEV_SNOW+1,', stop'
       STOP
     ENDIF
   ENDIF
  
   IF (PRESENT(SNOW_MELT_AGE_EXT)) THEN
     SHAPE_3D=SHAPE(SNOW_MELT_AGE_EXT)
     NLEV=SHAPE_3D(3)
     IF (NLEV /= NLEV_SNOW+1) THEN
       PRINT *,'In POCHVA scheme input arrays snow_melt_age_ext do not have correct dimension for 3 index,'
       PRINT *,' that must be 1:',NLEV_SNOW+1,', stop'
       STOP
     ENDIF
   ENDIF
  
   IF (PRESENT(RHO_SNOW_EXT)) THEN
     SHAPE_3D=SHAPE(RHO_SNOW_EXT)
     NLEV=SHAPE_3D(3)
     IF (NLEV /= NLEV_SNOW+1) THEN
       PRINT *,'In POCHVA scheme input arrays rho_snow_ext do not have correct dimension for 3 index,'
       PRINT *,' that must be 1:',NLEV_SNOW+1,', stop'
       STOP
     ENDIF
   ENDIF

   CALL POCHVA_INIT(NPOINT_X_EXT, NPOINT_Y_EXT, NLEV_SOIL_EXT, &
 LEV_SOIL_EXT, MASK_SOIL_EXT, TS_EXT, T_SURF_EXT, QS_EXT, SNOW, &
 QSOIL_MAX, QSOIL_MIN, CSOIL, RHOSOIL, PSISOIL, WCONDSOIL, PARBSOIL, PARCSOIL,&
 QSOIL_REL_WILT, QSOIL_REL_REF, LAIVEG, LAIVEG_MAX, FRACVEG, ROOTVEG, &
 IND_LEV_TBOTTOM, IND_LEV_QBOTTOM, VALMISS_EXT, DTIME, &
 TS_SURF_EXT, QS_SURF_EXT, FICE_SOIL_EXT, FICE_SOIL_SURF_EXT, &
 FLAG_SNOW_INIT, LEV_SNOW_EXT, TSNOW_EXT, FICE_SNOW_EXT, SNOW_AGE_EXT, SNOW_MELT_AGE_EXT, SNOW_DIRT_EXT)

! test --->
  DO K=1,N_OUTPOINT2
    IF (PROC_IND==IPROC_OUTPOINT(K)) THEN
    IUNIT=100+K
    WRITE (FILE_OUT,'(A,I2.2,A)') "initial_",K,".dat2"
    OPEN (IUNIT,FILE=FILE_OUT,STATUS="UNKNOWN")
      I=I_OUTPOINT(K)
      J=J_OUTPOINT(K)
      WRITE (IUNIT,*) NLEV_SOIL_EXT,LEV_SOIL_EXT, TS_EXT(I,J,:), &
 T_SURF_EXT(I,J), QS_EXT(I,J,:), SNOW(I,J), &
 QSOIL_MAX(I,J,:), QSOIL_MIN(I,J,:), &
 CSOIL(I,J,:), RHOSOIL(I,J,:), PSISOIL(I,J,:), WCONDSOIL(I,J,:), PARBSOIL(I,J,:), PARCSOIL(I,J,:),&
 QSOIL_REL_WILT(I,J,:), QSOIL_REL_REF(I,J,:), LAIVEG(I,J), LAIVEG_MAX, FRACVEG(I,J), ROOTVEG(I,J), VALMISS_EXT, &
 IND_LEV_TBOTTOM(I,J), IND_LEV_QBOTTOM(I,J), &
 QS_SURF_EXT(I,J), LEV_SNOW_EXT(I,J,:), TSNOW_EXT(I,J,:), SNOW_AGE_EXT(I,J,:), SNOW_MELT_AGE_EXT(I,J,:)
    ENDIF
  ENDDO
  DO K=1,N_OUTPOINT2
    IF (PROC_IND==IPROC_OUTPOINT(K)) THEN
      IUNIT=100+K
      CLOSE (IUNIT)
    ENDIF
  ENDDO
  DO K=1,N_OUTPOINT2
    IF (PROC_IND==IPROC_OUTPOINT(K)) THEN
      IUNIT=200+K
      WRITE (FILE_OUT,'(A,I2.2,A)') "atm_forcing_",K,".dat2"
      OPEN (IUNIT,FILE=FILE_OUT,STATUS="UNKNOWN")
    ENDIF 
  ENDDO
  DO K=1,N_OUTPOINT
    IF (PROC_IND==IPROC_OUTPOINT(K)) THEN
      IUNIT=300+(K-1)*20
      IUNIT=IUNIT+1
      WRITE (FILE_OUT,'(A,I2.2,A)') "tsoil_",K,".dat"
      OPEN (IUNIT,FILE=FILE_OUT,STATUS="UNKNOWN")
      IUNIT=IUNIT+1
      WRITE (FILE_OUT,'(A,I2.2,A)') "frac_sice_",K,".dat"
      OPEN (IUNIT,FILE=FILE_OUT,STATUS="UNKNOWN")
      IUNIT=IUNIT+1
      WRITE (FILE_OUT,'(A,I2.2,A)') "flux_entropy_soil_",K,".dat"
      OPEN (IUNIT,FILE=FILE_OUT,STATUS="UNKNOWN")
      IUNIT=IUNIT+1
      WRITE (FILE_OUT,'(A,I2.2,A)') "qsoil_",K,".dat"
      OPEN (IUNIT,FILE=FILE_OUT,STATUS="UNKNOWN")
      IUNIT=IUNIT+1
      WRITE (FILE_OUT,'(A,I2.2,A)') "runoff_",K,".dat"
      OPEN (IUNIT,FILE=FILE_OUT,STATUS="UNKNOWN")
      IUNIT=IUNIT+1
      WRITE (FILE_OUT,'(A,I2.2,A)') "flux_water_soil_",K,".dat"
      OPEN (IUNIT,FILE=FILE_OUT,STATUS="UNKNOWN")
      IUNIT=IUNIT+1
      WRITE (FILE_OUT,'(A,I2.2,A)') "tsnow_",K,".dat"
      OPEN (IUNIT,FILE=FILE_OUT,STATUS="UNKNOWN")
      IUNIT=IUNIT+1
      WRITE (FILE_OUT,'(A,I2.2,A)') "lev_snow_",K,".dat"
      OPEN (IUNIT,FILE=FILE_OUT,STATUS="UNKNOWN")
      IUNIT=IUNIT+1
      WRITE (FILE_OUT,'(A,I2.2,A)') "fice_snow_",K,".dat"
      OPEN (IUNIT,FILE=FILE_OUT,STATUS="UNKNOWN")
      IUNIT=IUNIT+1
      WRITE (FILE_OUT,'(A,I2.2,A)') "snow_age_",K,".dat"
      OPEN (IUNIT,FILE=FILE_OUT,STATUS="UNKNOWN")
      IUNIT=IUNIT+1
      WRITE (FILE_OUT,'(A,I2.2,A)') "snow_melt_age_",K,".dat"
      OPEN (IUNIT,FILE=FILE_OUT,STATUS="UNKNOWN")
      IUNIT=IUNIT+1
      WRITE (FILE_OUT,'(A,I2.2,A)') "flux_entropy_snow_",K,".dat"
      OPEN (IUNIT,FILE=FILE_OUT,STATUS="UNKNOWN")
      IUNIT=IUNIT+1
      WRITE (FILE_OUT,'(A,I2.2,A)') "surf_",K,".dat"
      OPEN (IUNIT,FILE=FILE_OUT,STATUS="UNKNOWN")
      IUNIT=IUNIT+1
      WRITE (FILE_OUT,'(A,I2.2,A)') "atm_",K,".dat"
      OPEN (IUNIT,FILE=FILE_OUT,STATUS="UNKNOWN")
    ENDIF 
  ENDDO
! <---

 ENDIF

! Setting to zero of some additional variables

 FLUX_ENTROPY_FOREST_SOIL(:,:)=0.

! Definition of some internal variables

 DTSTEP=DTIME
 FLUX_PREC_LIQ(1:NPOINT_X,1:NPOINT_Y)=FLUX_RAIN(1:NPOINT_X,1:NPOINT_Y)
 FLUX_PREC_ICE(1:NPOINT_X,1:NPOINT_Y)=FLUX_SNOW(1:NPOINT_X,1:NPOINT_Y)

 IF (FLAG_PAR_CHANGE /= 0) THEN
   CALL POCHVA_PAR_CHANGE(MASK_SOIL_EXT, TS_EXT, T_SURF_EXT, SNOW, &
 LAIVEG, LAIVEG_MAX, FRACVEG, IND_LEV_TBOTTOM, IND_LEV_QBOTTOM, VALMISS_EXT, &
 FLAG_SNOW_INIT, LEV_SNOW_EXT, TSNOW_EXT, FICE_SNOW_EXT, SNOW_AGE_EXT, SNOW_MELT_AGE_EXT)
   FLAG_PAR_CHANGE=0
 ENDIF

!IF (PROC_IND==0) PRINT *,' Pochva step ',STEP_NUM

! test --->
DO K=1,N_OUTPOINT2
  IF (PROC_IND==IPROC_OUTPOINT(K)) THEN
    IUNIT=200+K
    I=I_OUTPOINT(K)
    J=J_OUTPOINT(K)
    WRITE (IUNIT,1002) STEP_NUM,H_ATM_BOTT_LEV(I,J),PSURF(I,J),PATM_BOTT_LEV(I,J),TATM_BOTT_LEV(I,J),QATM_BOTT_LEV(I,J),&
   FLUX_RAIN(I,J),FLUX_SNOW(I,J),FLUX_RAD_TOT(I,J),FLUX_RAD_VIS(I,J),KTURBH(I,J),KTURBWV(I,J),&
   LAI_VEGLOW(I,J),FVEGLOW(I,J),DTSTEP 
  ENDIF
ENDDO
! <---

 DO J=1,NPOINT_Y
 DO I=1,NPOINT_X
 IF (MASK_SOIL(I,J) == 1) THEN

   NLEV=NLEV_SOIL_WATER(I,J)

! For eventual evepotraspitation processes

   DO K=0,NLEV-1
     IF (QS_REL(I,J,K) > QS_REL_VEGWILT(I,J,K).AND.QS_REL(I,J,K) < QS_REL_VEGREF(I,J,K)) THEN
       QS_REL_VEG(I,J,K)=(QS_REL(I,J,K)-QS_REL_VEGWILT(I,J,K))/(QS_REL_VEGREF(I,J,K)-QS_REL_VEGWILT(I,J,K))
     ELSEIF (QS_REL(I,J,K) <= QS_REL_VEGWILT(I,J,K)) THEN
        QS_REL_VEG(I,J,K)=0.
     ELSE
        QS_REL_VEG(I,J,K)=1.
     ENDIF
   ENDDO

 ENDIF
 ENDDO
 ENDDO

! Determination of temperature and air specific humidity at the complex terrain surface
! Definition of turbolent fluxe of water vapour between the bottom atmospheric levels ans terrain surface
! Definition of turbolent and radiation fluxes of entropy between the bottom atmospheric levels ans terrain surface
! Parametrization of vegetation evapotraspiration - forecast of soil water content increment

 CALL POCHVA_ATM_FLUX (H_ATM_BOTT_LEV, PSURF, PATM_BOTT_LEV, TATM_BOTT_LEV, QATM_BOTT_LEV, KTURBH, KTURBWV, &
 FLUX_RAD_TOT, FLUX_RAD_VIS, FLUX_HEAT_SPEC, FLUX_HEAT_LAT, FLUX_Q, FLAG_TURB_FLUX, STEP_NUM)

! Parametrization of soil water transport processes - forecast of soil water content increment

 CALL POCHVA_SOIL_WATER(STEP_NUM)

! Definition of soil water decrease becuase of evapotraspiration of low anf forest (high) vegetation

 CALL POCHVA_VEGETATION(STEP_NUM)

! Parametrization of soil thermal exchange - forecast of soil temperature

 CALL POCHVA_SOIL_TEMPER(STEP_NUM)

! Parametrization of water transport and termodynamical processes in the snow -
! forecast of the temperature, water content and cover fraction of snow

 CALL POCHVA_SNOW(SNOW, STEP_NUM)

! Definition of output variables

 DO I=1,NPOINT_X
 DO J=1,NPOINT_Y
 IF (MASK_SOIL(I,J) == 1) THEN
   T_SURF_EXT(I,J)=TSURF(I,J)
   QV_SURF_EXT(I,J)=QSURF(I,J)
   IF (FLAG_TURB_FLUX == 0) THEN
     FLUX_HEAT_SPEC(I,J)=FLUX_ENTROPY_ATM_SPEC(I,J)*TSURF(I,J)
     FLUX_HEAT_LAT(I,J)=FLUX_ENTROPY_ATM_LAT(I,J)*TSURF(I,J)
     FLUX_Q(I,J)=FLUX_WV(I,J)
   ENDIF
   NLEV=MIN(NLEV_SOIL_EXT, IND_LEV_TBOTTOM(I,J))
   DO K=1,NLEV
     TS_EXT(I,J,K)=EXP(TS(I,J,K))*T0
   ENDDO
   IF (PRESENT(TS_SURF_EXT)) THEN
     TS_SURF_EXT(I,J)=EXP(TS(I,J,0))*T0
   ENDIF
   IF (PRESENT(FICE_SOIL_EXT)) THEN
     FICE_SOIL_EXT(I,J,1:NLEV)=FRAC_SICE(I,J,1:NLEV)
   ENDIF
   IF (PRESENT(FICE_SOIL_SURF_EXT)) THEN
     FICE_SOIL_SURF_EXT(I,J)=FRAC_SICE(I,J,0)
   ENDIF
   NLEV=MIN(NLEV_SOIL_EXT, IND_LEV_QBOTTOM(I,J))
   DO K=1,NLEV
     QS_EXT(I,J,K)=QS(I,J,K)
   ENDDO
   IF (PRESENT(QS_SURF_EXT)) THEN
     QS_SURF_EXT(I,J)=QS(I,J,0)
   ENDIF
   FSNOW(I,J)=FSNCOVER(I,J)
   IF (PRESENT(LEV_SNOW_EXT)) LEV_SNOW_EXT(I,J,:)=VALMISS_EXT
   IF (PRESENT(TSNOW_EXT)) TSNOW_EXT(I,J,:)=VALMISS_EXT
   IF (PRESENT(FICE_SNOW_EXT)) FICE_SNOW_EXT(I,J,:)=VALMISS_EXT
   IF (PRESENT(SNOW_AGE_EXT)) SNOW_AGE_EXT(I,J,:)=VALMISS_EXT
   IF (PRESENT(SNOW_MELT_AGE_EXT)) SNOW_MELT_AGE_EXT(I,J,:)=VALMISS_EXT
   IF (PRESENT(RHO_SNOW_EXT)) RHO_SNOW_EXT(I,J,:)=VALMISS_EXT
   DO K=0,NLEV_SNOW
     IF (TSNOW(I,J,K) == VAL_MISSING) EXIT
   ENDDO
   NLEV=K-2
   IF (NLEV >= 0) THEN
     IF (PRESENT(LEV_SNOW_EXT)) LEV_SNOW_EXT(I,J,1:NLEV+2)=LEV_SNOW(I,J,0:NLEV+1)
     IF (PRESENT(TSNOW_EXT)) TSNOW_EXT(I,J,1:NLEV+2)=EXP(TSNOW(I,J,0:NLEV+1))*T0
     IF (PRESENT(FICE_SNOW_EXT)) FICE_SNOW_EXT(I,J,1:NLEV+2)=FICE_SNOW(I,J,0:NLEV+1)
     IF (PRESENT(SNOW_AGE_EXT)) SNOW_AGE_EXT(I,J,1:NLEV+2)=SNOW_AGE(I,J,0:NLEV+1)
     IF (PRESENT(SNOW_MELT_AGE_EXT)) SNOW_MELT_AGE_EXT(I,J,1:NLEV+2)=SNOW_MELT_AGE(I,J,0:NLEV+1)
     IF (PRESENT(RHO_SNOW_EXT)) RHO_SNOW_EXT(I,J,1:NLEV+2)=RHOSNOW(I,J,0:NLEV+1)
   ENDIF
   CALL SNOW_ALBEDO_DEF(SNALBED, I, J) ! Snow albedo definition
   SNOW_ALBEDO(I,J)=SNALBED
   NLEV=MAX(IND_LEV_TBOTTOM(I,J), 1)
   FLUX_HEAT_BOTTOM(I,J)=FLUX_SOIL_ENTROPY(I,J,NLEV)*EXP(0.5*(TS(I,J,NLEV)+TS(I,J,NLEV-1)))*T0
   NLEV=MAX(IND_LEV_QBOTTOM(I,J),1)
   FLUX_WATER_BOTTOM(I,J)=FLUX_SOIL_WATER_BOTTOM_DIAG(I,J)*RHOW ! m/s ---> kg/m^2/s
   RUN_OFF_TOT_EXT(I,J)=SUM(RUN_OFF(I,J,0:NLEV))
   RUN_OFF_EXT(I,J)=0.
   DO K=0,NLEV
     RUN_OFF_EXT(I,J)=RUN_OFF_EXT(I,J)+RUN_OFF(I,J,K)
     IF (LEV_SOIL(K) > 1.5) EXIT
   ENDDO
   RUN_OFF_TOT_EXT(I,J)=RUN_OFF_TOT_EXT(I,J)/DTSTEP
   RUN_OFF_EXT(I,J)=RUN_OFF_EXT(I,J)/DTSTEP

 ENDIF
 ENDDO
 ENDDO

! test --->
1001 FORMAT (100E16.8)
1002 FORMAT (I8,100E16.8)
 DO K=1,N_OUTPOINT
!IF (PROC_IND==IPROC_OUTPOINT(K).AND.MOD((STEP_NUM-1)*INT(DTSTEP), 3600) == 0) THEN
   IF (PROC_IND==IPROC_OUTPOINT(K)) THEN
   I=I_OUTPOINT(K)
   J=J_OUTPOINT(K)
   NLEV=NLEV_SOIL_HEAT(I,J)
   IUNIT=300+(K-1)*20
   IUNIT=IUNIT+1 ! tsoil, 301, 321, 341, 361, 381, 401, 421, 441, 461, 481
   WRITE (IUNIT,1001) STEP_NUM*DTSTEP/3600.,&
 EXP(TS(I,J,0:NLEV))*T0-T0
   IUNIT=IUNIT+1 ! frac_sice, 302, 322, 342, 362, 382, 402, 422, 442, 462, 482
   WRITE (IUNIT,1001) STEP_NUM*DTSTEP/3600.,&
 FRAC_SICE(I,J,0:NLEV)
   IUNIT=IUNIT+1 ! flux_entropy_soil, 303, 323, 343, 363, 383, 403, 423, 443, 463, 483
   WRITE (IUNIT,1001) STEP_NUM*DTSTEP/3600.,&
 FLUX_SOIL_ENTROPY(I,J,0:NLEV)
NLEV=NLEV_SOIL_WATER(I,J)
   IUNIT=IUNIT+1 ! qsoil, 304, 324, 344, 364, 384, 404, 424, 444, 464, 484
   WRITE (IUNIT,1001) STEP_NUM*DTSTEP/3600.,&
 (QS(I,J,0:NLEV)-QSMIN(I,J,0:NLEV))/(QSMAX(I,J,0:NLEV)-QSMIN(I,J,0:NLEV))
   IUNIT=IUNIT+1 ! runoff, 305, 325, 345, 365, 385, 405, 425, 445, 465, 485
   WRITE (IUNIT,1001) STEP_NUM*DTSTEP/3600.,&
 RUN_OFF_EXT(I,J),RUN_OFF_TOT_EXT(I,J),RUN_OFF(I,J,0:NLEV)
   IUNIT=IUNIT+1 ! flux_water_soil, 306, 326, 346, 366, 386, 406, 426, 446, 466, 486
   WRITE (IUNIT,1001) STEP_NUM*DTSTEP/3600.,&
 FLUX_SOIL_WATER(I,J,0:NLEV)
   IUNIT=IUNIT+1 ! tsnow, 307, 327, 347, 367, 387, 407, 427, 447, 467, 487
   WRITE (IUNIT,1001) STEP_NUM*DTSTEP/3600.,&
 EXP(TSNOW(I,J,0:NLEV_SNOW))*T0-T0
   IUNIT=IUNIT+1 ! lev_snow, 308, 328, 348, 368, 388, 408, 428, 448, 468, 488
   WRITE (IUNIT,1001) STEP_NUM*DTSTEP/3600.,&
 LEV_SNOW(I,J,0:NLEV_SNOW)
   IUNIT=IUNIT+1 ! fice_snow, 309, 329, 349, 369, 389, 409, 429, 449, 469, 489
   WRITE (IUNIT,1001) STEP_NUM*DTSTEP/3600.,&
 FICE_SNOW(I,J,0:NLEV_SNOW)
   IUNIT=IUNIT+1 ! snow_age, 310, 330, 350, 370, 390, 410, 430, 450, 470, 490
   WRITE (IUNIT,1001) STEP_NUM*DTSTEP/3600.,&
 SNOW_AGE(I,J,0:NLEV_SNOW)
   IUNIT=IUNIT+1 ! snow_melt_age, 311, 331, 351, 371, 391, 411, 431, 451, 471, 491
   WRITE (IUNIT,1001) STEP_NUM*DTSTEP/3600.,&
 SNOW_MELT_AGE(I,J,0:NLEV_SNOW)
   IUNIT=IUNIT+1 ! flux_entropy_snow, 312, 332, 352, 372, 392, 412, 432, 452, 472, 492
   WRITE (IUNIT,1001) STEP_NUM*DTSTEP/3600.,&
 FLUX_SNOW_ENTROPY(I,J,0:NLEV_SNOW)
   IUNIT=IUNIT+1 ! surf, 313, 333, 353, 373, 393, 413, 433, 453, 473, 493
   WRITE (IUNIT,1001) STEP_NUM*DTSTEP/3600.,&
 FSNCOVER(I,J),&  !  2
 TSURF(I,J)-T0,&  !  3
 QSURF(I,J),&  !  4
 FLUX_ENTROPY_ATM(I,J),&  !  5
 FLUX_ENTROPY_ATM_SOIL(I,J),&  !  6
 FLUX_ENTROPY_ATM_SNOW(I,J),&  !  7
 FLUX_ENTROPY_SNOW_SOIL(I,J),&  !  8
 FLUX_WV(I,J),&  !  9
 FLUX_WV_SOIL(I,J),&  ! 10
 FLUX_WV_SNOW(I,J),&  ! 11
 FLUX_WV_VEGLOW_DRY(I,J),&  ! 12
 FLUX_WV_VEGLOW_WET(I,J),&  ! 13
 FLUX_PREC_LIQ(I,J),&  ! 14
 FLUX_PREC_ICE(I,J),&  ! 15
 FLUX_W_LIQ_ATM_SOIL(I,J),&  ! 16
 FLUX_W_LIQ_ATM_SNOW(I,J),&  ! 17
 FLUX_W_ICE_ATM_SNOW(I,J),&  ! 18
 FLUX_W_LIQ_ATM_VEGLOW(I,J),&  ! 19
 FLUX_WATER_SNOW_SOIL(I,J),&  ! 20
 FVEGLOW(I,J),&  ! 21
 FVEGLOWWET(I,J),&  ! 22
 LAI_VEGLOW(I,J),&  ! 23
 QW_VEGLOW(I,J),&  ! 24
 FLUX_RAD_TOT(I,J),& ! 25
 FLUX_RAD_VIS(I,J),& ! 26
 FLUX_ENTROPY_ATM_SPEC(I,J),& ! 27
 FLUX_ENTROPY_ATM_LAT(I,J),&  ! 28
 FLUX_ENTROPY_ATM_RAD(I,J),&  ! 29
 SNOW_ALBEDO(I,J),&           ! 30
 SNOW_DIRT(I,J)               ! 31
   IUNIT=IUNIT+1 ! atm, 314, 334, 354, 374, 394, 414, 434, 454, 474, 494
   WRITE (IUNIT,1001) STEP_NUM*DTSTEP/3600.,&
 TATM_BOTT_LEV(I,J)-T0,& ! 2
 QATM_BOTT_LEV(I,J),&    ! 3
 PATM_BOTT_LEV(I,J),&    ! 4
 H_ATM_BOTT_LEV(I,J),&   ! 5
 KTURBH(I,J),&           ! 6
 KTURBWV(I,J)            ! 7
ENDIF
 ENDDO
! <---

RETURN
END SUBROUTINE POCHVA
! ---------------------------------------------------------------------------------------------------------------------------------
!    Pochva:
!    Numerical scheme for the parametrization of thermo-hydrous processes
!    in soil and at its surface (own soil, vegetation, snow cover)
!
! Scheme is elaborated by O.Drofa, CNR-ISAC, Italy (o.drofa-at-isac.cnr.it)
!
! ---------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE POCHVA_ATM_FLUX (DZA, PSURF, PATM1, TATM1, QATM1, KTURBH, KTURBWV, &
 RADIATION_FLUX, RAD_VIS_FLUX, FLUX_HEAT_SPEC, FLUX_HEAT_LAT, FLUX_Q, FLAG_TURB_FLUX, ISTEP)

! Determination of temperature and air specific humidity at the complex terrain surface (bare soil, dry and wet "leaf", snow)
! Definition of turbolent flux of water vapour between the bottom atmospheric levels ans terrain surface
! Definition of turbolent and radiation fluxes of entropy between the bottom atmospheric levels ans terrain surface

!USE MODULE_POCHVA
!IMPLICIT NONE

! Input variables: NPOINT_X, NPOINY_Y, MASK_SOIL, MASK_SOIL_AXED, 
! DZA, PSURF, PATM1, TATM1, QATM1, KTURBH, KTURBWV, 
! FLUX_PREC_LIQ, FLUX_PREC_ICE, RADIATION_FLUX, RAD_VIS_FLUX, 
! FLUX_HEAT_SPEC, FLUX_HEAT_LAT, FLUX_Q, FLAG_TURB_FLUX,
! QS, QS_REL, QS_REL_VEG, TS, TSNOW, TFOREST, FSNCOVER, FVEGLOW, FFOREST,
! FVEGLOWWET, FFORESTWET, FLAG_SNOW_THICK 

! Output variables: QSURF, TSURF, FLUX_ENTROPY_ATM, FLUX_ENTROPY_ATM_SPEC,
! FLUX_ENTROPY_ATM_LAT, FLUX_ENTROPY_ATM_RAD, FLUX_ENTROPY_ATM_SOIL,
! FLUX_ENTROPY_ATM_FOREST,
! FLUX_WV, FLUX_WV_SOIL, FLUX_WV_SNOW,
! FLUX_WV_VEGLOW_DRY, FLUX_WV_VEGLOW_WET, FLUX_WV_FOREST_DRY, FLUX_WV_FOREST_WET,
! FLUX_SOIL_WATER_VEG, FLUX_W_LIQ_ATM_VEGLOW, FLUX_W_LIQ_ATM_FOREST,
! FLUX_W_LIQ_ATM_SNOW, FLUX_W_ICE_ATM_SNOW

REAL, DIMENSION(:,:) :: DZA, PSURF, PATM1, TATM1, QATM1, KTURBH, KTURBWV, &
 RADIATION_FLUX, RAD_VIS_FLUX, FLUX_HEAT_SPEC, FLUX_HEAT_LAT, FLUX_Q
REAL, DIMENSION(NPOINT_X,NPOINT_Y) :: KTURBWV_1M
INTEGER :: FLAG_TURB_FLUX, ISTEP, NLEV, I, J, K, FLAG_EVAP_VEGLOW, FLAG_EVAP_FOREST
REAL :: FRAC_SOIL, FRAC_LEAF_DRY, FRAC_LEAF_WET, FRAC_FOREST_DRY, FRAC_FOREST_WET,&
 TSURF_OLD, TSNOW_SURF, ESW, ESI, QSAT, QSNOW, QSAT_FOREST, QVSOIL, QVLEAF, ALFA_LEAF,&
 ROA0, ZTSOIL0, ZT0, ZT1, PD0, PD1, PV0, PV1, SDRY0, SDRY1, SWV0, SWV1,& 
 ZFR1, ZFLUX_SW_ATM, ZFLUX_SI_ATM, ZFLUX1, ZFLUX2
 
 FLUX_ENTROPY_ATM(:,:)=0.
 FLUX_ENTROPY_ATM_SPEC(:,:)=0.
 FLUX_ENTROPY_ATM_LAT(:,:)=0.
 FLUX_ENTROPY_ATM_RAD(:,:)=0.
 FLUX_ENTROPY_ATM_SOIL(:,:)=0.
 FLUX_ENTROPY_ATM_SNOW(:,:)=0.
 FLUX_ENTROPY_ATM_FOREST(:,:)=0.
 FLUX_WV_SOIL(:,:)=0.
 FLUX_WV_SNOW(:,:)=0.
 FLUX_WV_VEGLOW_DRY(:,:)=0.
 FLUX_WV_VEGLOW_WET(:,:)=0.
 FLUX_WV_FOREST_DRY(:,:)=0.
 FLUX_WV_FOREST_WET(:,:)=0.
 FLUX_W_LIQ_ATM_SOIL(:,:)=0.
 FLUX_W_LIQ_ATM_SNOW(:,:)=0.
 FLUX_W_ICE_ATM_SNOW(:,:)=0.
 FLUX_W_LIQ_ATM_VEGLOW(:,:)=0.
 FLUX_W_LIQ_ATM_FOREST(:,:)=0.

! Turbulent exchange coefficient at 1 m over surface:
 KTURBWV_1M(:,:) = KTURBWV(:,:)/DZA(:,:)

 DO J=1,NPOINT_Y
 DO I=1,NPOINT_X
 IF (MASK_SOIL(I,J) == 1) THEN

   TSURF_OLD=TSURF(I,J)

   ZTSOIL0=EXP(TS(I,J,0))*T0 ! Temperature at the bare soil surface

! Fractions of various surface type

! Grid box surface fraction of bare soil
   FRAC_SOIL=(1.-FVEGLOW(I,J))*(1.-FFOREST(I,J))

! Grid box surface fraction of dry leaf low vegetation
   ZFR1=(1.-FFOREST(I,J))*FVEGLOW(I,J)
   FRAC_LEAF_DRY=ZFR1*(1.-FVEGLOWWET(I,J))

! Grid box surface fraction of wet leaf low vegetation
   FRAC_LEAF_WET=ZFR1*FVEGLOWWET(I,J)

! Grid box surface fraction of dry leaf of the forest vegetation:
   FRAC_FOREST_DRY=FFOREST(I,J)*(1.-FFORESTWET(I,J))

! Grid box surface fraction of wet leaf of the forest vegetation:
   FRAC_FOREST_WET=FFOREST(I,J)*FFORESTWET(I,J)

! ---- Snow surface ----------------------------------------------------

   IF (TSNOW(I,J,0) > VAL_MISSING) THEN ! Snow exists

!!!     TSNOW_SURF=EXP(TSNOW(I,J,0))*T0
     IF (FLAG_SNOW_THICK(I,J) == 1) THEN ! Snow is thick
       TSNOW_SURF=EXP(TSNOW(I,J,0))*T0
     ELSE ! Snow is thin
       TSNOW_SURF=TSURF_OLD
     ENDIF

     CALL PRES_VAP_SAT(TSNOW_SURF, ESW, ESI)

     IF (TSNOW(I,J,0) >= 0.) THEN
       QSNOW=ESW*EPS/(PSURF(I,J)-ESW*(1.-EPS))
     ELSE
       QSNOW=ESI*EPS/(PSURF(I,J)-ESI*(1.-EPS))
     ENDIF

   ELSE

     TSNOW_SURF=EXP(TS(I,J,0))*T0
     QSNOW=0.

   ENDIF

   NLEV=NLEV_SOIL_WATER(I,J)

! ---- Bare soil and low vegetation surface ----------------------------
   
   CALL PRES_VAP_SAT(ZTSOIL0, ESW, ESI)

   IF (ZTSOIL0 >= T0) THEN

! Soil (with vegetation low vegetaion) surface temperature is upper then T0

     QSAT=ESW*EPS/(PSURF(I,J)-ESW*(1.-EPS))

! Evapotraspiration by low vegetation

! Control of possibilty for this process

     FLAG_EVAP_VEGLOW=1
     DO K=0,NLEV-1
       IF (DZ_DZSUM_VEGLOW(I,J,K) <= 0.) EXIT
       IF (TS(I,J,K) < 0.) FLAG_EVAP_VEGLOW=0
     ENDDO
     IF (QATM1(I,J) > QSAT) FLAG_EVAP_VEGLOW=0
     IF (LAI_VEGLOW(I,J) <= 0.) FLAG_EVAP_VEGLOW=0
     IF (RAD_VIS_FLUX(I,J) <= 0.) FLAG_EVAP_VEGLOW=0

     IF (FLAG_EVAP_VEGLOW == 1) THEN

! Low vegetation evapotraspiration is possible

       CALL QV_LEAF(QVLEAF, ALFA_LEAF, QSAT, QATM1(I,J), KTURBWV_1M(I,J), RAD_VIS_FLUX(I,J), LAI_VEGLOW(I,J), LAI_MAX,&
 QS_REL_VEG(I,J,0:NLEV-1), DZ_DZSUM_VEGLOW(I,J,0:NLEV-1), QV_LEAF_BETA, NLEV)

     ELSE

! Low vegetation evapotraspiration is impossible

       QVLEAF=QATM1(I,J)

     ENDIF

   ELSE

! Soil (with vegetation low vegetaion) surface temperature is below then T0

     FLAG_EVAP_VEGLOW=0
     QSAT=ESI*EPS/(PSURF(I,J)-ESI*(1.-EPS))
     QVLEAF=QATM1(I,J)

   ENDIF

   CALL QV_BARE_SOIL(QVSOIL, QSAT, QATM1(I,J), KTURBWV_1M(I,J), QS(I,J,0), QSMAX(I,J,0), QSMIN(I,J,0), &
 PAR_B(I,J,0), MASK_SOIL_AXED(I,J))

! ---- Forest (high) vegetation surface --------------------------------
   
   CALL PRES_VAP_SAT(TFOREST(I,J), ESW, ESI)

   IF (TFOREST(I,J) >= T0) THEN

! Forest (high) vegetation surface temperature is upper then T0

     QSAT_FOREST=ESW*EPS/(PSURF(I,J)-ESW*(1.-EPS))

! Evapotraspiration by forest (high) vegetation

! Control of possibilty for this process

     FLAG_EVAP_FOREST=1
     DO K=0,NLEV-1
       IF (DZ_DZSUM_FOREST(I,J,K) <= 0.) EXIT
       IF (TS(I,J,K) < TS0) FLAG_EVAP_FOREST=0
     ENDDO
     IF (QATM1(I,J) < QSAT_FOREST) FLAG_EVAP_FOREST=0
     IF (LAI_FOREST(I,J) <= 0.) FLAG_EVAP_FOREST=0
     IF (RAD_VIS_FLUX(I,J) <= 0.) FLAG_EVAP_FOREST=0

     IF (FLAG_EVAP_FOREST == 1) THEN

! Forest (high) vegetation evapotraspiration is possible

       CALL QV_LEAF(QFOREST(I,J), ALFA_LEAF, QSAT_FOREST, QATM1(I,J), KTURBWV_1M(I,J), RAD_VIS_FLUX(I,J), &
 LAI_FOREST(I,J), LAI_MAX,&
 QS_REL_VEG(I,J,0:NLEV-1), DZ_DZSUM_FOREST(I,J,0:NLEV-1), QV_FOREST_BETA, NLEV)

     ELSE

! Forest (high) vegetation evapotraspiration is impossible

       QFOREST(I,J)=QATM1(I,J)

     ENDIF

   ELSE

! Forest (high) vegetation surface temperature is below then T0

     QFOREST=QATM1(I,J)
     QSAT_FOREST=ESI*EPS/(PSURF(I,J)-ESI*(1.-EPS))

   ENDIF

! ----------------------------------------------------------------------

! Temperature and water vapour mixing ratio at the complex terrain surface (representative for a whole grib box)

   IF (TSNOW(I,J,0) <= VAL_MISSING) THEN ! Snow not exists

     TSURF(I,J) = ZTSOIL0
     QSURF(I,J) = QVSOIL*FRAC_SOIL + QVLEAF*FRAC_LEAF_DRY + QSAT*FRAC_LEAF_WET

   ELSE ! Snow exists

     IF (FLAG_SNOW_THICK(I,J) == 1) THEN ! Snow is thick
       TSURF(I,J) = TSNOW_SURF
     ELSE ! Snow is thin
       TSURF(I,J) = ZTSOIL0
     ENDIF
     QSURF(I,J) = QSNOW

   ENDIF

   ROA0=PSURF(I,J)/(RD*TSURF(I,J)*(1.+EP*QSURF(I,J)))

!------------------ FLUXES -------------------------------------------

! Turbulent flux of water vapour between the lowest atmespheric layer and the complex terrain surface (representative for a whole grib box)

! This total value must be devide between various "water storage":
! bare soil;
! water evaporated by the low vegetation dry leafs;
! water deposited on the low vegetation wet leafs;
! snow.

   IF (FLAG_TURB_FLUX == 0) THEN
! Turbulent fluxes are defined here
!     FLUX_WV(I,J)=ROA0*KTURBWV(I,J)/DZA(I,J)*(QATM1(I,J)-QSURF(I,J))
     ZFLUX1=ROA0*KTURBWV(I,J)/DZA(I,J)*(QATM1(I,J)-QSURF(I,J))
     ZFLUX2=MIN(ABS(ZFLUX1), 2.E-4)
     FLUX_WV(I,J)=SIGN(ZFLUX2,ZFLUX1)
   ELSE
! Turbulent fluxes are defined by input data
     FLUX_WV(I,J)=FLUX_Q(I,J)
   ENDIF

   IF (TSNOW(I,J,0) <= VAL_MISSING.OR.FLAG_SNOW_THICK(I,J) == 0) THEN ! Snow not exists or snow is not think
! If snow is not thick, then water vapour phase change is simulated at the soil surface,
! in common with simulation of other water phase change

     IF (FLUX_WV(I,J) < 0.) THEN

! Water vapour flux upward - evaporation

! Bare soil
     FLUX_WV_SOIL(I,J)=FLUX_WV(I,J)*FRAC_SOIL

! Dry leaf surface of the low vegetation
       IF (FLAG_EVAP_VEGLOW /= 0 ) THEN
         FLUX_WV_VEGLOW_DRY(I,J)=FLUX_WV(I,J)*FRAC_LEAF_DRY ! Evaportraspiration is possible
       ELSE
         FLUX_WV_VEGLOW_DRY(I,J)=0. ! Evaportraspiration is impossible
       ENDIF

! Wet leaf surface of the low vegetation
       FLUX_WV_VEGLOW_WET(I,J)=FLUX_WV(I,J)-FLUX_WV_VEGLOW_DRY(I,J)-FLUX_WV_SOIL(I,J)
       IF (QW_VEGLOW(I,J)+FLUX_WV_VEGLOW_WET(I,J)*DTSTEP < 0.) THEN
         FLUX_WV_VEGLOW_WET(I,J)=-QW_VEGLOW(I,J)/DTSTEP
         FLUX_WV_SOIL(I,J)=FLUX_WV(I,J)-FLUX_WV_VEGLOW_WET(I,J)-FLUX_WV_VEGLOW_DRY(I,J)
       ENDIF

     ELSE

! Water vapour flux downward - deposition

! Dry leaf surface of the low vegetation
       FLUX_WV_VEGLOW_DRY(I,J)=0.

! Wet leaf surface of the low vegetation
       FLUX_WV_VEGLOW_WET(I,J)=FLUX_WV(I,J)*(FRAC_LEAF_DRY+FRAC_LEAF_WET)
       IF (QW_VEGLOW(I,J)+FLUX_WV_VEGLOW_WET(I,J)*DTSTEP > QW_VEGLOW_MAX(I,J)) &
 FLUX_WV_VEGLOW_WET(I,J)=(QW_VEGLOW_MAX(I,J)-QW_VEGLOW(I,J))/DTSTEP

! Bare soil
       FLUX_WV_SOIL(I,J)=FLUX_WV(I,J)-FLUX_WV_VEGLOW_WET(I,J)

      ENDIF

! Snow
     FLUX_WV_SNOW(I,J)=0.

   ELSE ! Snow exists

! Bare soil
     FLUX_WV_SOIL(I,J)=0.
! Dry leaf surface of the low vegetation
     FLUX_WV_VEGLOW_DRY(I,J)=0.
! Wet leaf surface of the low vegetation
     FLUX_WV_VEGLOW_WET(I,J)=0.

! Snow
     FLUX_WV_SNOW(I,J)=FLUX_WV(I,J)

   ENDIF

! Flux of atmospheric precipitation liquid water and ice on the complex terrain surface (representative for a whole grib box)

! The total value of liquid water flux must be devide between various "water storage":
! water deposited on the low vegetation wet leafs;
! soil;
! snow.
! The whole value of ice water flux is directed to the unique "water storage" that is snow

   IF (TSNOW(I,J,0) <= VAL_MISSING) THEN ! Snow not exists

! Liquid precipitation on the leaf surface of the low vegetation
     FLUX_W_LIQ_ATM_VEGLOW(I,J)=FLUX_PREC_LIQ(I,J)*FVEGLOW(I,J)
     IF (QW_VEGLOW(I,J)+FLUX_W_LIQ_ATM_VEGLOW(I,J)*DTSTEP > QW_VEGLOW_MAX(I,J)) &
     FLUX_W_LIQ_ATM_VEGLOW(I,J)=(QW_VEGLOW_MAX(I,J)-QW_VEGLOW(I,J))/DTSTEP

! Liquid precipitation on the soil surface
     FLUX_W_LIQ_ATM_SOIL(I,J)=MAX( FLUX_PREC_LIQ(I,J)-FLUX_W_LIQ_ATM_VEGLOW(I,J), 0.)

! Liquid precipitation on the snow surface
     FLUX_W_LIQ_ATM_SNOW(I,J)=0.

   ELSE ! Snow exists

! Liquid precipitation on the leaf surface of the low vegetation
     FLUX_W_LIQ_ATM_VEGLOW(I,J)=0.

! Liquid precipitation on the soil surface
     FLUX_W_LIQ_ATM_SOIL(I,J)=0.

! Liquid precipitation on the snow surface
     FLUX_W_LIQ_ATM_SNOW(I,J)=FLUX_PREC_LIQ(I,J)

   ENDIF

! Solid (ice) precipitation on the snow surface
   FLUX_W_ICE_ATM_SNOW(I,J)=FLUX_PREC_ICE(I,J)

! Flux of total entropy on the complex terrain surface (representative for a whole grib box)

! Turbulent flux of water vapour between the lowest atmespheric layer and the complex terrain surface (representative for a whole grib box)
! Flux of total entropy consists from 5 components:
! 1) turbulent flux of dry air entropy between the lowest atmespheric layer and the surface;
! 2) turbulent flux of water vopour entropy between the lowest atmespheric layer and the surface;
! 3) flux of the entropy of liquid water precipitation;
! 4) flux of the entropy of solid water (ice) precipitation;
! 5) flux of the total radiation entropy

   IF (FLAG_TURB_FLUX == 0) THEN

! Turbulent fluxes are defined here

     ZT0=ALOG(TSURF(I,J)/T0)
     ZT1=ALOG(TATM1(I,J)/T0)
     PV0=PSURF(I,J)*QSURF(I,J)/(QSURF(I,J)*(1.-EPS)+EPS)
     PV1=PATM1(I,J)*QATM1(I,J)/(QATM1(I,J)*(1.-EPS)+EPS)
     PD0=PSURF(I,J)-PV0
     PD1=PATM1(I,J)-PV1

! 1) Turbulent flux of dry air entropy between the lowest atmespheric layer and the surface

     SDRY0=(1.-QSURF(I,J))*(CPD*ZT0-RD*ALOG(PD0/P0))
     SDRY1=(1.-QATM1(I,J))*(CPD*ZT1-RD*ALOG(PD1/P0))

!     FLUX_ENTROPY_ATM_SPEC(I,J)=ROA0*KTURBH(I,J)/DZA(I,J)*(SDRY1-SDRY0)
     ZFLUX1=ROA0*KTURBH(I,J)/DZA(I,J)*(SDRY1-SDRY0)
     ZFLUX2=MIN(ABS(ZFLUX1), 1.)
     FLUX_ENTROPY_ATM_SPEC(I,J)=SIGN(ZFLUX2,ZFLUX1)

! 2) Turbulent flux of water vapour entropy between the lowest atmespheric layer and the surface

     SWV0=QSURF(I,J)*(CPV*ZT0-RV*ALOG(PV0/E0)+LIVDT0)
     SWV1=QATM1(I,J)*(CPV*ZT1-RV*ALOG(PV1/E0)+LIVDT0)

!     FLUX_ENTROPY_ATM_LAT(I,J)=ROA0*KTURBWV(I,J)/DZA(I,J)*(SWV1-SWV0)
     ZFLUX1=ROA0*KTURBWV(I,J)/DZA(I,J)*(SWV1-SWV0)
     ZFLUX2=MIN(ABS(ZFLUX1), 2.)
     FLUX_ENTROPY_ATM_LAT(I,J)=SIGN(ZFLUX2,ZFLUX1)

   ELSE

! Turbulent fluxes are defined by input data

! 1) Turbulent flux of dry air entropy between the lowest atmespheric layer and the surface
     FLUX_ENTROPY_ATM_SPEC(I,J)=FLUX_HEAT_SPEC(I,J)/TSURF_OLD

! 2) Turbulent flux of water vapour entropy between the lowest atmespheric layer and the surface
     FLUX_ENTROPY_ATM_LAT(I,J)=FLUX_HEAT_LAT(I,J)/TSURF_OLD

   ENDIF ! FLAG_TURB_FLUX

!! 3) Flux of the entropy of liquid water precipitation;
!
!   ZFLUX_SW_ATM=FLUX_PREC_LIQ(I,J)*(CW*(ZT0*0.3+ZT1*0.7)+LIWDT0)
!
!! 4) Flux of the entropy of solid water (ice) precipitation;
!
!   Z1=MIN(ZT0*0.3+ZT1*0.7, 0.)
!   ZFLUX_SI_ATM=FLUX_PREC_ICE(I,J)*CI*Z1

! 5) Flux of the total radiation entropy

   FLUX_ENTROPY_ATM_RAD(I,J)=RADIATION_FLUX(I,J)/TSURF_OLD

! Total entropy flux

!   FLUX_ENTROPY_ATM(I,J) = FLUX_ENTROPY_ATM_SPEC(I,J) + FLUX_ENTROPY_ATM_LAT(I,J) + ZFLUX_SW_ATM + ZFLUX_SI_ATM + FLUX_ENTROPY_ATM_RAD(I,J) 
   FLUX_ENTROPY_ATM(I,J) = FLUX_ENTROPY_ATM_SPEC(I,J) + FLUX_ENTROPY_ATM_LAT(I,J) + FLUX_ENTROPY_ATM_RAD(I,J) 

! The total entropy flux flow on the surface of soil with low vegetation or on the snow surface

   IF (TSNOW(I,J,0) > VAL_MISSING.AND.FLAG_SNOW_THICK(I,J) == 1) THEN ! Snow exists and snow is think
     FLUX_ENTROPY_ATM_SOIL(I,J) = 0.
     FLUX_ENTROPY_ATM_SNOW(I,J) = FLUX_ENTROPY_ATM(I,J)
   ELSE
     FLUX_ENTROPY_ATM_SOIL(I,J) = FLUX_ENTROPY_ATM(I,J)
     FLUX_ENTROPY_ATM_SNOW(I,J) = 0.
   ENDIF

   FLUX_ENTROPY_ATM_FOREST(I,J) = 0.

 ENDIF
 ENDDO
 ENDDO

RETURN
END SUBROUTINE POCHVA_ATM_FLUX
! ---------------------------------------------------------------------------------------------------------------------------------
!    Pochva:
!    Numerical scheme for the parametrization of thermo-hydrous processes
!    in soil and at its surface (own soil, vegetation, snow cover)
!
! Scheme is elaborated by O.Drofa, CNR-ISAC, Italy (o.drofa-at-isac.cnr.it)
!
! ---------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE POCHVA_SOIL_WATER(ISTEP)
! Procedure of soil water transport because of сapillary and gravity forces action
! ---> forecast of water soil content increment

!USE MODULE_POCHVA
!IMPLICIT NONE

! Input variables: NPOINT_X, NPOINY_Y, MASK_SOIL, MASK_SOIL_AXED, NLEV_SOIL_MAX, NLEV_SOIL_WATER, D_LEV_SOIL, D_LEV_SOIL_H,
! FRAC_SICE, QSMAX, QSMIN, PSIS_SAT, PAR_B, PAR_C, KAPPA_H, FLUX_WV_SOIL, FLUX_W_LIQ_ATM_SOIL

! Output variables: PSIS, FLUX_SOIL_WATER, QS_REL

! Input and output variables: QS

INTEGER :: ISTEP, NLEV, I, J, K, KH
REAL :: Z01, Z02, Z03, Z04, Z05, Z06, Z11, Z12, Z13, Z14, Z15, Z16, ERROR, &
 SUM_DQS, WEIGHT, ZMIN
REAL, DIMENSION(0:NLEV_SOIL_MAX) :: KAPPA, DFLUX_SOIL_WATER_DZ, ASW
REAL, DIMENSION(3) :: ZLIMIT

 FLUX_SOIL_WATER(:,:,:)=0.
 DQS(:,:,:)=0.
 RUN_OFF(:,:,:)=0.

!open (31,file="oshibka_pochva.txt",status="unknown",position="append")
 DO J=1,NPOINT_Y
 DO I=1,NPOINT_X

 IF (MASK_SOIL(I,J) == 1.AND.MASK_SOIL_AXED(I,J) /= 1) THEN

  NLEV=NLEV_SOIL_WATER(I,J)

!do k=0,nlev
!if (qs(i,j,k)>qsmax(i,j,k).or.qs(i,j,k)<qsmin(i,j,k)) write(31,*) 'POCHVA_SOIL_WATER start oshibka ',&
! proc_ind,i,j,k,qs(i,j,k),qsmax(i,j,k),qsmin(i,j,k),istep
!enddo

  QS_OLD(I,J,:)=QS(I,J,:)
  QSI(I,J,:)=QS(I,J,:)*FRAC_SICE(I,J,:)

! Soil hydraulic potential, hydraulic conductivity at all soil full-levels
! and the denominator (ASW) in the right part of prognostic equation at all full soil full- levels
! (ASW is the function of QS, FRAC_SICE, PSIS, QSMAX)

   DO K=0,NLEV
     IF (FRAC_SICE(I,J,K) < 0.999) THEN
!!!       Z01=QSMAX(I,J,K)-QS(I,J,K)*FRAC_SICE(I,J,K)
       Z02=QS(I,J,K)*(1.-FRAC_SICE(I,J,K))
!!!       Z03=Z01/Z02
!!!       Z04=QSMAX(I,J,K)/Z01
!!!       Z05=Z03**PAR_B(I,J,K)
!!!       Z06=Z04**PAR_C(I,J,K)
!!!       PSIS(I,J,K)=PSIS_SAT(I,J,K)*Z05*Z06
       Z03=QSMAX(I,J,K)/Z02
       Z04=Z03**PAR_B(I,J,K)
       PSIS(I,J,K)=PSIS_SAT(I,J,K)*Z04
!!!       KAPPA(K)=KAPPA_SAT(I,J,K)*((Z02/QSMAX(I,J,K))**(2.*PAR_B(I,J,K)+3.))
       KAPPA(K)=KAPPA_SAT(I,J,K)*((Z02/QSMAX(I,J,K))**(PAR_B(I,J,K)+3.))
       ASW(K)=-QSMAX(I,J,K)*PSIS_SAT(I,J,K)*PAR_B(I,J,K)/QS(I,J,K)*Z04
!!!       IF (FRAC_SICE(I,J,K) > 0.) THEN
!!!         KAPPA(K)=KAPPA_SAT(I,J,K)*((Z02/Z01)**(2.*PAR_B(I,J,K)+3.))
!!!         ASW(K)=QSMAX(I,J,K)*PSIS_SAT(I,J,K)*&
!!! ( -PAR_B(I,J,K)*Z03**(PAR_B(I,J,K)-1)*(FRAC_SICE(I,J,K)/Z02+Z01/(Z02**2))*Z06 +&
!!! Z05*PAR_C(I,J,K)*Z04**(PAR_C(I,J,K)-1)*QSMAX(I,J,K)*FRAC_SICE(I,J,K)/(Z01**2) )
!!!       ELSE
!!!         KAPPA(K)=KAPPA_SAT(I,J,K)*((QS(I,J,K)/QSMAX(I,J,K))**(2.*PAR_B(I,J,K)+3.))
!!!         ASW(K)=QSMAX(I,J,K)*PSIS_SAT(I,J,K)*&
!!! ( -PAR_B(I,J,K)*Z03**(PAR_B(I,J,K)-1)*(Z01/(Z02**2)) )
!!!       ENDIF
     ELSE
       PSIS(I,J,K)=VAL_MISSING
!!!       KAPPA(K)=0.
       ASW(K)=1.E8
     ENDIF
   ENDDO 

! Soil water flux at half-levels 

   DO KH=1,NLEV
     K=KH
     IF (PSIS(I,J,K-1) /= VAL_MISSING.AND.PSIS(I,J,K) /= VAL_MISSING) THEN
       KAPPA_H(I,J,KH)=0.5*(KAPPA(K-1)+KAPPA(K))
       FLUX_SOIL_WATER(I,J,KH)=KAPPA_H(I,J,KH)*((-PSIS(I,J,K)+PSIS(I,J,K-1))/D_LEV_SOIL(K))
! Limitation of water fluxes
       IF (FLUX_SOIL_WATER(I,J,KH) /= 0.) THEN
!         ZLIMIT(1)=ABS(0.5*(PSIS(I,J,K)+PSIS(I,J,K-1))*D_LEV_SOIL(K)/DTSTEP*0.1) ! 10% is maximum change of average layer hydraulic potential
         ZLIMIT(1)=ABS(PSIS(I,J,K)*D_LEV_SOIL(K)/DTSTEP*0.1) ! 10% is maximum change of hydraulic potential at level k
         ZLIMIT(2)=ABS(PSIS(I,J,K-1)*D_LEV_SOIL(K)/DTSTEP*0.1) ! 10% is maximum change of hydraulic potential at level k-1
         ZLIMIT(3)=ABS(FLUX_SOIL_WATER(I,J,KH))
         ZMIN=MINVAL(ZLIMIT)
         FLUX_SOIL_WATER(I,J,KH)=SIGN(ZMIN,FLUX_SOIL_WATER(I,J,KH))
       ENDIF
     ELSE
       KAPPA_H(I,J,KH)=0.
       FLUX_SOIL_WATER(I,J,KH)=0.
     ENDIF
   ENDDO

  FLUX_SOIL_WATER(I,J,0)=(FLUX_WV_SOIL(I,J)+FLUX_W_LIQ_ATM_SOIL(I,J)+FLUX_WATER_SNOW_SOIL(I,J))/RHOW

! Flux of soil water at bottom layer, because of flux of soil hydraulic potential
! is not correct for water flux control in the case of freesing soil

   IF (ANY (FRAC_SICE(I,J,NLEV-1:NLEV) > 0.) ) THEN
     Z01=LEV_SOIL_H(NLEV)-LEV_SOIL(NLEV-1)
     Z02=LEV_SOIL(NLEV)-LEV_SOIL_H(NLEV)
     FLUX_SOIL_WATER_BOTTOM_DIAG(I,J)=KAPPA_H(I,J,NLEV)/D_LEV_SOIL(NLEV)* &
 ((QS(I,J,NLEV-1)-QSI(I,J,NLEV-1))*Z01-(QS(I,J,NLEV)-QSI(I,J,NLEV))*Z02)
   ELSE
     FLUX_SOIL_WATER_BOTTOM_DIAG(I,J)=FLUX_SOIL_WATER(I,J,NLEV)
   ENDIF

! Definition of the increment of soil water content at all full-levels

   DO K=0,NLEV-1
     DFLUX_SOIL_WATER_DZ(K)=(-FLUX_SOIL_WATER(I,J,K+1)+FLUX_SOIL_WATER(I,J,K))/D_LEV_SOIL_H(K)
     DQS(I,J,K)=DTSTEP/ASW(K)*DFLUX_SOIL_WATER_DZ(K)

! Forecast of soil water content

     QS(I,J,K)=QS(I,J,K)+DQS(I,J,K)
     QS_REL(I,J,K)=(QS(I,J,K)-QSMIN(I,J,K))/(QSMAX(I,J,K)-QSMIN(I,J,K)) 

   ENDDO

! Conservation of column water

   ERROR=SUM(QS(I,J,0:NLEV-1)*D_LEV_SOIL_H(0:NLEV-1))-SUM(QS_OLD(I,J,0:NLEV-1)*D_LEV_SOIL_H(0:NLEV-1))
   ERROR=ERROR+(-FLUX_SOIL_WATER(I,J,0)+FLUX_SOIL_WATER_BOTTOM_DIAG(I,J))*DTSTEP
!   PRINT *,'ERROR 1',ERROR

   IF (ABS(ERROR) > 1.E-6 ) THEN

     IF (ERROR > 0.) THEN

      Z01=SUM(QS_REL(I,J,0:NLEV-1))
      IF (Z01 > 1.E-6) THEN
        DO K=0,NLEV-1
          WEIGHT=(QS_REL(I,J,K))/Z01
          QS(I,J,K)=(QS(I,J,K)*D_LEV_SOIL_H(K)-ERROR*WEIGHT)/D_LEV_SOIL_H(K)
          DQS(I,J,K)=QS(I,J,K)-QS_OLD(I,J,K)
        ENDDO 
      ENDIF

     ELSE

      Z01=SUM(1.-QS_REL(I,J,0:NLEV-1))
      IF (Z01 > 1.E-6) THEN
        DO K=0,NLEV-1
          WEIGHT=(1.-QS_REL(I,J,K))/Z01
          QS(I,J,K)=(QS(I,J,K)*D_LEV_SOIL_H(K)-ERROR*WEIGHT)/D_LEV_SOIL_H(K)
          DQS(I,J,K)=QS(I,J,K)-QS_OLD(I,J,K)
        ENDDO 
      ENDIF

     ENDIF

   ENDIF

!   ERROR=SUM(QS(I,J,0:NLEV-1)*D_LEV_SOIL_H(0:NLEV-1))-SUM(QS_OLD(I,J,0:NLEV-1)*D_LEV_SOIL_H(0:NLEV-1))
!   ERROR=ERROR+(-FLUX_SOIL_WATER(I,J,0)+FLUX_SOIL_WATER_BOTTOM_DIAG(I,J))*DTSTEP
!   PRINT *,'ERROR 2',ERROR

! Definition of run-off

   DO K=0,NLEV-1
     IF (QS(I,J,K) > QSMAX(I,J,K)) THEN
       RUN_OFF(I,J,K)=RUN_OFF(I,J,K)+(QS(I,J,K)-QSMAX(I,J,K))*RHOW*D_LEV_SOIL_H(K)
       QS(I,J,K)=QSMAX(I,J,K)
     ENDIF
     QS(I,J,K)=MAX(QS(I,J,K), QSI(I,J,K))
     QS(I,J,K)=MAX(QS(I,J,K), QSMIN(I,J,K))
     QS_REL(I,J,K)=(QS(I,J,K)-QSMIN(I,J,K))/(QSMAX(I,J,K)-QSMIN(I,J,K)) 
     DQS(I,J,K)=QS(I,J,K)-QS_OLD(I,J,K)
   ENDDO

! Flux of soil water (kg/m^2/s) and it's derivative respect to z (kg/m^3/s)

   FLUX_SOIL_WATER(I,J,:)=FLUX_SOIL_WATER(I,J,:)*RHOW

!do k=0,nlev
!if (qs(i,j,k)>qsmax(i,j,k).or.qs(i,j,k)<qsmin(i,j,k)) write (31,*) 'POCHVA_SOIL_WATER finish oshibka ',&
! proc_ind,i,j,k,qs(i,j,k),qsmax(i,j,k),qsmin(i,j,k),istep
!enddo

 ENDIF

 IF (MASK_SOIL(I,J) == 1.AND.MASK_SOIL_AXED(I,J) == 1) THEN

   NLEV=NLEV_SOIL_WATER(I,J)
   FLUX_SOIL_WATER(I,J,0)=FLUX_WV_SOIL(I,J)+FLUX_W_LIQ_ATM_SOIL(I,J)+FLUX_WATER_SNOW_SOIL(I,J)
   RUN_OFF(I,J,0)=FLUX_SOIL_WATER(I,J,0)*DTSTEP
   FLUX_SOIL_WATER(I,J,1:NLEV)=0.
   RUN_OFF(I,J,1:NLEV)=0.

 ENDIF

 ENDDO
 ENDDO
!close (31)

RETURN
END SUBROUTINE POCHVA_SOIL_WATER
! ---------------------------------------------------------------------------------------------------------------------------------
!    Pochva:
!    Numerical scheme for the parametrization of thermo-hydrous processes
!    in soil and at its surface (own soil, vegetation, snow cover)
!
! Scheme is elaborated by O.Drofa, CNR-ISAC, Italy (o.drofa-at-isac.cnr.it)
!
! ---------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE POCHVA_VEGETATION(ISTEP)

! Dinamics of water beacuse of vegetation activity:
! soil water decrease becuase of evapotraspiration of low anf forest (high) vegetation;
! definition of water content on low vegetation leaf and on hight (forest) vegetatiob leaf,
! definition of wet leaf fraction for low ang hight (forest) vegetation.

!USE MODULE_POCHVA
!IMPLICIT NONE

! Input variables: NPOINT_X, NPOINY_Y, NLEV_SOIL_MAX, NLEV_SOIL_WATER,
! MASK_SOIL, MASK_SOIL_AXED, FLUX_WV_VEGLOW_DRY, FLUX_WV_VEGLOW_WET,
! FLUX_WV_FOREST_DRY, FLUX_WV_FOREST_WET, QS_REL_VEG,&
! DTSTEP, RHOW, QSMAX, QSMIN, D_LEV_SOIL_H, DZSUM_VEGLOW, DZSUM_FOREST 

! Output variables: FLUX_SOIL_WATER_VEG

! Input and output variables: DQS

INTEGER :: ISTEP, NLEV, I, J, K
REAL :: ZFLUX_WV, ZDQSREL, ZDQS, ZFRAC_K, ZDQW1, ZDQW2, Z1
REAL, DIMENSION(0:NLEV_SOIL_MAX) :: ZFLUX_SW_LOWVEG, ZFLUX_SW_FOREST

 FLUX_SOIL_WATER_VEG(:,:,:)=0.

!open (31,file="oshibka_pochva.txt",status="unknown",position="append")
 DO J=1,NPOINT_Y
 DO I=1,NPOINT_X
 IF (MASK_SOIL(I,J) == 1.AND.MASK_SOIL_AXED(I,J) /= 1) THEN

   NLEV=NLEV_SOIL_WATER(I,J)

!do k=0,nlev
!if (qs(i,j,k)>qsmax(i,j,k).or.qs(i,j,k)<qsmin(i,j,k)) write (31,*) 'POCHVA_VEGETATION start oshibka ',&
! proc_ind,i,j,k,qs(i,j,k),qsmax(i,j,k),qsmin(i,j,k),istep
!enddo

! Definition of soil water decrease becuase of evapotraspiration of low anf forest (high) vegetation

! Evapotraspiration by low vegetation

   ZFLUX_SW_LOWVEG(:)=0.
   IF (FLUX_WV_VEGLOW_DRY(I,J) /= 0..AND.QV_LEAF_BETA /= 0.) THEN

     ZFLUX_WV=FLUX_WV_VEGLOW_DRY(I,J)*DTSTEP/(RHOW*DZSUM_VEGLOW(I,J))
     ZFLUX_SW_LOWVEG(NLEV)=0.
     DO K=NLEV-1,0,-1
       ZFRAC_K=QS_REL_VEG(I,J,K)*DZ_DZSUM_VEGLOW(I,J,K)/QV_LEAF_BETA
       ZDQSREL=ZFLUX_WV*ZFRAC_K
       ZDQS=ZDQSREL*(QSMAX(I,J,K)-QSMIN(I,J,K))
       Z1=MIN( MAX( QS(I,J,K)+ZDQS, QSMIN(I,J,K)), QSMAX(I,J,K))
       DQS(I,J,K)=DQS(I,J,K)+Z1-QS(I,J,K)
       QS(I,J,K)=Z1
       ZFLUX_SW_LOWVEG(K)=ZFLUX_SW_LOWVEG(K+1)-FLUX_WV_VEGLOW_DRY(I,J)*ZFRAC_K
     ENDDO

   ENDIF

! Evapotraspiration by high vegetation (forest)

   ZFLUX_SW_FOREST(:)=0.
   IF (FLUX_WV_FOREST_DRY(I,J) /= 0..AND.QV_FOREST_BETA /= 0.) THEN

     ZFLUX_WV=FLUX_WV_FOREST_DRY(I,J)*DTSTEP/(RHOW*DZSUM_FOREST(I,J))
     ZFLUX_SW_FOREST(NLEV)=0.
     DO K=NLEV-1,0,-1
       ZFRAC_K=QS_REL_VEG(I,J,K)*DZ_DZSUM_FOREST(I,J,K)/QV_FOREST_BETA
       ZDQSREL=ZFLUX_WV*ZFRAC_K
       ZDQS=ZDQSREL*(QSMAX(I,J,K)-QSMIN(I,J,K))
       DQS(I,J,K)=DQS(I,J,K)+ZDQS
       QS(I,J,K)=QS(I,J,K)+ZDQS
       ZFLUX_SW_FOREST(K)=ZFLUX_SW_FOREST(K+1)-FLUX_WV_FOREST_DRY(I,J)*ZFRAC_K
     ENDDO

!do k=0,nlev
!if (qs(i,j,k)>qsmax(i,j,k).or.qs(i,j,k)<qsmin(i,j,k)) write (31,*) 'POCHVA_VEGETATION finish oshibka ',&
! proc_ind,i,j,k,qs(i,j,k),qsmax(i,j,k),qsmin(i,j,k),istep
!enddo

   ENDIF

   FLUX_SOIL_WATER_VEG(I,J,1:NLEV)=ZFLUX_SW_LOWVEG(1:NLEV)+ZFLUX_SW_FOREST(1:NLEV)

! Definition of water content on low vegetation leaf and wet fraction of this leafs

   ZDQW1=FLUX_W_LIQ_ATM_VEGLOW(I,J)*DTSTEP
   ZDQW2=FLUX_WV_VEGLOW_WET(I,J)*DTSTEP
   QW_VEGLOW(I,J)=MIN(MAX(QW_VEGLOW(I,J)+ZDQW1+ZDQW2, 0.), QW_VEGLOW_MAX(I,J))
   IF (QW_VEGLOW(I,J) > 0..AND.QW_VEGLOW_MAX(I,J) > 0.) THEN
     FVEGLOWWET(I,J)=(QW_VEGLOW(I,J)/QW_VEGLOW_MAX(I,J))**0.7
   ELSE
     FVEGLOWWET(I,J)=0.
   ENDIF

! Definition of water content on hight (forest) vegetation leaf and wet fraction of this leafs

   ZDQW1=FLUX_W_LIQ_ATM_FOREST(I,J)*DTSTEP
   ZDQW2=FLUX_WV_FOREST_WET(I,J)*DTSTEP
   QW_FOREST(I,J)=MIN(MAX(QW_FOREST(I,J)+ZDQW1+ZDQW2, 0.), QW_FOREST_MAX(I,J))
   IF (QW_FOREST(I,J) > 0..AND.QW_FOREST_MAX(I,J) > 0.) THEN
     FFORESTWET(I,J)=(QW_FOREST(I,J)/QW_FOREST_MAX(I,J))**0.7
   ELSE
     FFORESTWET(I,J)=0.
   ENDIF

 ENDIF
 ENDDO
 ENDDO
!close (31)

RETURN
END SUBROUTINE POCHVA_VEGETATION
! ---------------------------------------------------------------------------------------------------------------------------------
!    Pochva:
!    Numerical scheme for the parametrization of thermo-hydrous processes
!    in soil and at its surface (own soil, vegetation, snow cover)
!
! Scheme is elaborated by O.Drofa, CNR-ISAC, Italy (o.drofa-at-isac.cnr.it)
!
! ---------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE POCHVA_SOIL_TEMPER(ISTEP)

! Simulation of soil thermal processes ---> forecast of the soil temperature

!USE MODULE_POCHVA
!IMPLICIT NONE

! Input variables: NPOINT_X, NPOINT_Y, MASK_SOIL, MASK_SOIL_AXED, NLEV_SOIL_MAX, NLEV_SOIL_HEAT, D_LEV_SOIL, D_LEV_SOIL_H,
! QS, QSI, RHOG, CG, RHOW, CW, RHOI, CI, RHOWCW, RHOICI, LIWDT0,
! FLUX_ENTROPY_ATM_SOIL, FLUX_ENTROPY_SNOW_SOIL, FLUX_ENTROPY_FOREST_SOIL,
! FLAG_SNOW_THICK, DLEV_SNOW_BOTTOM, NLEV_SNOW

! Output variables: DFRAC_SICE_DT, FLUX_SOIL_ENTROPY

! Input and output variables: TS, FRAC_SICE, FICE_SNOW

INTEGER :: ISTEP, NLEV, I, J, K, KH, K_TOP, NLEVS, ITER_MAX=20, ITER
REAL, DIMENSION(0:NLEV_SOIL_MAX) :: CSOIL, RHOSOIL, CRHOSOIL, TSTAR, S, S_INI, S_FIN, DS
REAL, DIMENSION(NLEV_SOIL_MAX) :: LAMBDASDCS
REAL :: F_OLD, F_NEW, T_OLD, T_NEW, DFDT_OLD, FICE, DFICE, DELTAT=1.E-7, &
 FS, DFS, Z1, Z2, Z3, Z4, SUM_DS, ERROR, WEIGHT, ZMIN
REAL, DIMENSION(2) :: ZLIMIT

 FLUX_SOIL_ENTROPY(:,:,:)=0.

! Definition of soil thermal conductivity at all full-levels

 CALL THERM_CONDUCT

!open (31,file="oshibka_pochva.txt",status="unknown",position="append")
 DO J=1,NPOINT_Y
 DO I=1,NPOINT_X
 IF (MASK_SOIL(I,J) == 1) THEN

   NLEV=NLEV_SOIL_HEAT(I,J)

!do k=0,nlev
!if (ts(i,j,k)>0.4.or.ts(i,j,k)<-0.45) write (31,*) 'POCHVA_SOIL_TEMPER start oshibka ',&
! proc_ind,i,j,k,ts(i,j,k),exp(ts(i,j,k))*t0-t0,istep
!enddo

! ------- Simulation of conductive thermal transport --------

! Top boundary condition (FLUX_ENTROPY_FOREST_SOIL not considerated!)

   IF (FLAG_SNOW_THICK(I,J) == 0) THEN !
     FLUX_SOIL_ENTROPY(I,J,0)=FLUX_ENTROPY_ATM_SOIL(I,J)
   ELSE
     FLUX_SOIL_ENTROPY(I,J,0)=FLUX_ENTROPY_SNOW_SOIL(I,J)
   ENDIF

! Total specific heat and total density of the humid soil at full-levels
! and wet soil entropy at the start 

   IF (MASK_SOIL_AXED(I,J) /= 1) THEN
     DO K=0,NLEV
       RHOSOIL(K)=RHOG(I,J,K)+QS(I,J,K)*((1.-FRAC_SICE(I,J,K))*RHOW+FRAC_SICE(I,J,K)*RHOI)
       CSOIL(K)=CG(I,J,K)+QS(I,J,K)*((1.-FRAC_SICE(I,J,K))*CW+FRAC_SICE(I,J,K)*CI)
       CRHOSOIL(K)=CSOIL(K)*RHOSOIL(K) 
       S_INI(K)=CRHOSOIL(K)*TS(I,J,K)
     ENDDO
   ELSE
     RHOSOIL(0:NLEV)=RHOG(I,J,0:NLEV)
     CSOIL(0:NLEV)=CG(I,J,0:NLEV)
     CRHOSOIL(0:NLEV)=CSOIL(0:NLEV)*RHOSOIL(0:NLEV) 
     S_INI(0:NLEV)=CRHOSOIL(0:NLEV)*TS(I,J,0:NLEV)
   ENDIF

! Soil thermal conductivity devided by soil total specific heat Total at half-levels

   DO KH=1,NLEV
     K=KH
     LAMBDASDCS(KH)=(LAMBDAG(I,J,K-1)/CSOIL(K-1)+LAMBDAG(I,J,K)/CSOIL(K))*0.5
   ENDDO

! Conductive thermal flux at half-levels

   Z1=0.05
   DO KH=1,NLEV
     K=KH
     FLUX_SOIL_ENTROPY(I,J,KH)=LAMBDASDCS(KH)*(-CSOIL(K)*TS(I,J,K)+CSOIL(K-1)*TS(I,J,K-1))/D_LEV_SOIL(K)
     IF (FLUX_SOIL_ENTROPY(I,J,KH) /= 0.) THEN
! Limitation of entropy fluxes
       ZLIMIT(1)=CRHOSOIL(K)*TS_REF*D_LEV_SOIL(K)/DTSTEP*Z1 ! 1% (5% for uppest layer) is maximum change of layer entropy with reference temperature TS_REF
       ZLIMIT(2)=ABS(FLUX_SOIL_ENTROPY(I,J,KH))
       ZMIN=MINVAL(ZLIMIT)
       FLUX_SOIL_ENTROPY(I,J,KH)=SIGN(ZMIN,FLUX_SOIL_ENTROPY(I,J,KH))
     ENDIF
   Z1=0.01
   ENDDO

! Soil temperature (ln(T/T0)) after conductive thermal transport simulation at the full-levels

   DO KH=0,NLEV-1
     K=KH
     DS(K)=DTSTEP*(-FLUX_SOIL_ENTROPY(I,J,KH+1)+FLUX_SOIL_ENTROPY(I,J,KH))/D_LEV_SOIL_H(KH)
     S_FIN(K)=S_INI(K)+DS(K)
     TSTAR(K)=TS(I,J,K)+DS(K)/CRHOSOIL(K)
   ENDDO

! Conservation of column entropy

   ERROR=SUM((S_FIN(0:NLEV-1)-S_INI(0:NLEV-1))*D_LEV_SOIL_H(0:NLEV-1))
   ERROR=ERROR+(-FLUX_SOIL_ENTROPY(I,J,0)+FLUX_SOIL_ENTROPY(I,J,NLEV))*DTSTEP

   IF (ERROR /= 0. ) THEN

     WEIGHT=1./FLOAT(NLEV)
     DO K=0,NLEV-1
       Z1=D_LEV_SOIL_H(K)*CRHOSOIL(K)
       TSTAR(K)=(TSTAR(K)*Z1-ERROR*WEIGHT)/Z1
       S_FIN(K)=TSTAR(K)*CRHOSOIL(K)
     ENDDO 

   ENDIF

 IF (MASK_SOIL_AXED(I,J) /= 1) THEN

! ------- Simulation of soil water phase changing --------

   IF (DLEV_SNOW_BOTTOM(I,J) > 0..AND.LEV_SNOW(I,J,0) /= VAL_MISSING) THEN ! Snow exists

! Soil surface (soil level number 0) with snow

     K_TOP=1

! Actual snow level number

     DO K=0,NLEV_SNOW
       IF (LEV_SNOW(I,J,K) == VAL_MISSING) EXIT
     ENDDO
     NLEVS=K-1

     K=0

     FICE=QSI(I,J,K)/QS(I,J,K)
     FS=FICE_SNOW(I,J,NLEVS)

     IF (TSTAR(K) > TS_FICE_CRIT) THEN

! Phase chinging is possible
 
       Z1=RHOWCW-RHOICI
       Z2=RHOW*LIWDT0

! Entropy of soil water (must be conserved during phase changing)

       S(K)=QS(I,J,K)*(RHOWCW*TSTAR(K)+Z2*(1.-FICE)-FICE*TSTAR(K)*Z1)*D_LEV_SOIL_H(0)+ & ! soil part
 (CI*TSTAR(K)+(1./FS-1.)*(CW*TSTAR(K)+LIWDT0))*DLEV_SNOW_BOTTOM(I,J) ! snow part

! Iteraction procedure

       T_OLD=TSTAR(K)
       CALL FRAC_SOIL_ICE_POINT(T_OLD, FICE, DFICE, I, J, K)
       FICE=QSI(I,J,K)/QS(I,J,K)
       F_OLD=FS
       CALL FRAC_SNOW_ICE(T_OLD, FS, DFS)
       FS=F_OLD
       F_OLD=QS(I,J,K)*(RHOWCW*T_OLD+Z2*(1.-FICE)-FICE*T_OLD*Z1)*D_LEV_SOIL_H(0)+ & ! soil part
 (CI*T_OLD+(1./FS-1.)*(CW*T_OLD+LIWDT0))*DLEV_SNOW_BOTTOM(I,J) ! snow part

       DO ITER=1,ITER_MAX

         DFDT_OLD=QS(I,J,K)*(RHOWCW-Z2*DFICE-(DFICE*T_OLD+FICE)*Z1)*D_LEV_SOIL_H(0)+ & ! soil part
 (CI+(1./FS-1.)*CW+(-1./FS**2*DFS)*(CW*T_OLD)+LIWDT0)*DLEV_SNOW_BOTTOM(I,J) ! snow part
         T_NEW=T_OLD+(S(K)-F_OLD)/DFDT_OLD
         CALL FRAC_SOIL_ICE_POINT(T_NEW, FICE, DFICE, I, J, K)
         CALL FRAC_SNOW_ICE(T_NEW, FS, DFS)
         F_NEW=QS(I,J,K)*(RHOWCW*T_NEW+Z2*(1.-FICE)-FICE*T_NEW*Z1)*D_LEV_SOIL_H(0)+ & ! soil part
 (CI*T_NEW+(1./FS-1.)*(CW*T_NEW+LIWDT0))*DLEV_SNOW_BOTTOM(I,J) ! snow part

         IF (ABS(T_NEW-T_OLD) < DELTAT) EXIT ! Itaractions are converged
         T_OLD=T_NEW
         F_OLD=F_NEW

       ENDDO

       TS(I,J,K)=T_NEW
       FRAC_SICE(I,J,K)=FICE
       FICE_SNOW(I,J,NLEVS)=FS
       IF (FICE_SNOW(I,J,NLEVS) < 1.) SNOW_MELT_AGE(I,J,NLEVS)=SNOW_MELT_AGE(I,J,NLEVS)+DTSTEP/86400.

     ELSE

! Phase chinging is impossible, temperature is the too low

       TS(I,J,K)=TSTAR(K)
       FRAC_SICE(I,J,K)=1.
       FICE_SNOW(I,J,NLEVS)=1.

     ENDIF

   ELSE ! Snow do not exist

     K_TOP=0

   ENDIF

! Soil levels

   DO K=K_TOP,NLEV-1

     FICE=QSI(I,J,K)/QS(I,J,K)

     IF ((TSTAR(K) < 0..OR.(TSTAR(K) >= 0..AND.FICE > 0.)).AND.TSTAR(K) > TS_FICE_CRIT) THEN

! Phase changing is possible
 
       Z1=RHOWCW-RHOICI
       Z2=RHOW*LIWDT0

! Entropy of soil water (must be conserved during phase changing)

       S(K)=QS(I,J,K)*(RHOWCW*TSTAR(K)+Z2*(1.-FICE)-FICE*TSTAR(K)*Z1)

! Iteraction procedure

       T_OLD=TSTAR(K)
       CALL FRAC_SOIL_ICE_POINT(T_OLD, FICE, DFICE, I, J, K)
       FICE=QSI(I,J,K)/QS(I,J,K)
       F_OLD=QS(I,J,K)*(RHOWCW*T_OLD+Z2*(1.-FICE)-FICE*T_OLD*Z1)
       DO ITER=1,ITER_MAX

         DFDT_OLD=QS(I,J,K)*(RHOWCW-Z2*DFICE-(DFICE*T_OLD+FICE)*Z1)
         T_NEW=T_OLD+(S(K)-F_OLD)/DFDT_OLD
         CALL FRAC_SOIL_ICE_POINT(T_NEW, FICE, DFICE, I, J, K)
         F_NEW=QS(I,J,K)*(RHOWCW*T_NEW-FICE*T_NEW*Z1+Z2*(1.-FICE))

         IF (ABS(T_NEW-T_OLD) < DELTAT) EXIT ! Itaractions are converged
         T_OLD=T_NEW
         F_OLD=F_NEW

       ENDDO

       TS(I,J,K)=T_NEW
       FRAC_SICE(I,J,K)=FICE

     ELSE

! Phase chinging is impossible

       TS(I,J,K)=TSTAR(K)
       IF (TSTAR(K) >= 0.) THEN
         FRAC_SICE(I,J,K)=0.
       ELSE
         FRAC_SICE(I,J,K)=1.
       ENDIF

     ENDIF

   ENDDO

 ELSE ! MASK_SOIL_AXED(I,J) == 1 ! Glacier

   TS(I,J,0:NLEV-1)=TSTAR(0:NLEV-1)

   IF (DLEV_SNOW_BOTTOM(I,J) > 0..AND.LEV_SNOW(I,J,0) /= VAL_MISSING) THEN ! Snow exists
! Eventual water phase changing in thw bottom half-layer of snow

! Actual snow level number

     DO K=0,NLEV_SNOW
       IF (LEV_SNOW(I,J,K) == VAL_MISSING) EXIT
     ENDDO
     NLEVS=K-1

     IF (TS(I,J,0) > 0..OR.FICE_SNOW(I,J,NLEVS) < 1.) THEN

! Phase chinging is possible

       K=0

       FS=FICE_SNOW(I,J,NLEVS)

       S(K)=DLEV_SNOW_BOTTOM(I,J)*FS*CI*TS(I,J,K)+DLEV_SNOW_BOTTOM(I,J)*(1.-FS)*(CW*TS(I,J,K)+LIWDT0)

       FS=1.
       T_NEW=S(K)/(DLEV_SNOW_BOTTOM(I,J)*CI)

       IF (T_NEW < 0.) THEN
! Solution is found
         TS(I,J,0)=T_NEW
         FICE_SNOW(I,J,NLEVS)=FS
       ELSE
! Snow is melting
         TS(I,J,0)=0.
         FS=1.-S(K)/(DLEV_SNOW_BOTTOM(I,J)*LIWDT0)
         FICE_SNOW(I,J,NLEVS)=MAX( MIN( FS, 1.), 1.E-8)
         SNOW_MELT_AGE(I,J,NLEVS)=SNOW_MELT_AGE(I,J,NLEVS)+DTSTEP/86400.
       ENDIF

     ENDIF ! Phase chinging is possible

   ENDIF ! Snow exists

   DO K=0,NLEV-1
     TS(I,J,K)=MIN(TS(I,J,K),0.)
   ENDDO

 ENDIF ! MASK_SOIL_AXED(I,J)

!do k=0,nlev
!if (ts(i,j,k)>0.4.or.ts(i,j,k)<-0.45) write (31,*) 'POCHVA_SOIL_TEMPER finish oshibka ',&
! proc_ind,i,j,k,ts(i,j,k),exp(ts(i,j,k))*t0-t0,istep
!enddo

 ENDIF
 ENDDO
 ENDDO
!close (31)

RETURN
END SUBROUTINE POCHVA_SOIL_TEMPER
! ---------------------------------------------------------------------------------------------------------------------------------
!    Pochva:
!    Numerical scheme for the parametrization of thermo-hydrous processes
!    in soil and at its surface (own soil, vegetation, snow cover)
!
! Scheme is elaborated by O.Drofa, CNR-ISAC, Italy (o.drofa-at-isac.cnr.it)
!
! ---------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE POCHVA_SNOW(SNOW, ISTEP)
! Procedure of termodinamical process  of snow soil water transport
! because of transport and phase transition of water mass and entropy
! ---> forecast of snow temperature, total water content, density and cover fraction

!USE MODULE_POCHVA
!IMPLICIT NONE

! Input variables: NPOINT_X, NPOINY_Y, MASK_SOIL, NLEV_SNOW, D_LEV_SNOW, D_LEV_SNOW_MIN,
! FLUX_ENTROPY_ATM_SNOW, FLUX_WV_SNOW, FLUX_W_LIQ_ATM_SNOW, FLUX_W_ICE_ATM_SNOW

! Output variables: FLUX_ENTROPY_SNOW_SOIL, FLUX_WATER_SNOW_SOIL, 
! FLUX_SNOW_ENTROPY, FLUX_SNOW_WATER, DLEV_SNOW_BOTTOM, SNOW, FSNCOVER

! Input and output variables: TSNOW, LEV_SNOW, FICE_SNOW, RHOSNOW,
! SNOW_AGE, SNOW_MELT_AGE, FLAG_SNOW_THICK

REAL, DIMENSION(:,:) :: SNOW

INTEGER :: ISTEP, I, J, K, KH, NLEV, NLEV_NEW, N
REAL :: ZRHOH, LAMBDA_SNOW_H, C_TOT_H, ZDLEV, ZA, ZB, ZNUMER, ZDENOM, FICE_NEW, ZT_NEW, FLUX_CONDUCT, &
 ZMIN, Z01, Z02, Z03, ZS, &
 DSNOW_ATM, DSNOW_ATM_LIQ, DSNOW_ATM_ICE, DSNOW_MELT, SNOW_OLD, SNOW_NEW, &
 SNOW_WATER_SUM, ERROR, SUM_DS, SUM_D_LEV, WEIGHT
REAL, DIMENSION(0:NLEV_SNOW+1) :: RHO_TOT, C_TOT, SNOW_ENTROPY, S_INI, S_FIN, DS,  TSTAR, &
 ZD_LEV_SNOW, SNOW_ICE=1., SNOW_WATER=0., ZFICE_NEW, SNOW_ENTROPY_NEW, &
 ZLEV_SNOW_NEW, ZD_LEV_SNOW_NEW, RHO_OLD, RHO_NEW, F_OLD, F_NEW
REAL, DIMENSION(2) :: ZLIMIT

!open (31,file="oshibka_pochva.txt",status="unknown",position="append")
 DO J=1,NPOINT_Y
 DO I=1,NPOINT_X
 IF (MASK_SOIL(I,J) ==1) THEN

! Actual snow level number at time step start

   DO K=0,NLEV_SNOW
     IF (LEV_SNOW(I,J,K) == VAL_MISSING) EXIT
   ENDDO
   NLEV=K-2

! NLEV is level number in snow from its top (number 0) to the last inside level
! over (!) bottom surface on the snow-soil boundary

   FLAG_SNOW_THICK(I,J)=0
   DLEV_SNOW_BOTTOM(I,J)=0.

   IF (NLEV >= 0) THEN

     TSNOW(I,J,NLEV+1)=TS(I,J,0)

     ZD_LEV_SNOW(:)=D_LEV_SNOW(I,J)
     ZD_LEV_SNOW(0)=LEV_SNOW(I,J,1)*0.5
     ZD_LEV_SNOW(1)=ZD_LEV_SNOW(0)+D_LEV_SNOW(I,J)*0.5
     IF (NLEV > 0) THEN
       ZD_LEV_SNOW(NLEV+1)=D_LEV_SNOW(I,J)*0.5
     ELSE
       ZD_LEV_SNOW(NLEV+1)=LEV_SNOW(I,J,1)*0.5
     ENDIF

     IF (LEV_SNOW(I,J,NLEV+1) >= D_LEV_SNOW_MIN) THEN
       FLAG_SNOW_THICK(I,J)=1
       DLEV_SNOW_BOTTOM(I,J)=ZD_LEV_SNOW(NLEV+1)
     ELSE
       DLEV_SNOW_BOTTOM(I,J)=LEV_SNOW(I,J,1)
     ENDIF

   ENDIF

! Entropy balance with "old" mass configuration

! If snow is less than threshold value (D_LEV_SNOW_MIN), then
! snow is "transparent" for heat flux, so, entripy exchange is direct between
! the atmospere and the soil surface, in this case snow water phase change 
! is simulated in soil thermical processes scheme: equation of the temperature
! in soil leve number 0

   IF (NLEV >= 0) THEN

     SNOW_OLD=LEV_SNOW(I,J,NLEV+1)

     IF (FLAG_SNOW_THICK(I,J) == 1) THEN ! Entropy process scheme in snow beacuase snow is thick enough

! ------- Simulation of conductive thermal transport --------

! Total specific heat and total density and wet entropy at the time step start

       DO K=0,NLEV+1
         Z01=1.-FICE_SNOW(I,J,K)
         RHO_TOT(K)=FICE_SNOW(I,J,K)*RHOI+Z01*RHOW
         C_TOT(K)=FICE_SNOW(I,J,K)*CI+Z01*CW
         S_INI(K)=C_TOT(K)*TSNOW(I,J,K)
       ENDDO 

! Top boundary condition
! Entropy flux in snow with mass axis using has unit (J*m)/(kg*s*K),
! but entropy flux it the surface and in soil has unit J/(K*s*m**2)
  
       FLUX_SNOW_ENTROPY(I,J,0)=FLUX_ENTROPY_ATM_SNOW(I,J)/RHO_TOT(0)

! Entropy flux at half levels:
  
       DO KH=1,NLEV+1
         K=KH
         ZRHOH=(RHOSNOW(I,J,K-1)+RHOSNOW(I,J,K))*0.5
         C_TOT_H=(C_TOT(K-1)+C_TOT(K))*0.5
         ZDLEV=LEV_SNOW(I,J,K)-LEV_SNOW(I,J,K-1) 
!! Z02 is the coefficient to simulate of "convective mixing" in the fresh snow 
!         Z01=(-TSNOW(I,J,K)+TSNOW(I,J,K-1))/ZDLEV
!         IF (ZRHOH <= RHOSNOW_FRESH*2.) THEN
!           Z02=MAX( MIN(1.E4*(Z01)**2., 10.), 1.)
!         ELSE
           Z02=1.
!         ENDIF
         LAMBDA_SNOW_H=2.60E-6*(ZRHOH)**2*Z02 ! J/(s*m*K)
         FLUX_CONDUCT=LAMBDA_SNOW_H/C_TOT_H*(-S_INI(K)+S_INI(K-1))/ZDLEV
! Limitation of conductive entropy fluxes
         IF (FLUX_CONDUCT /= 0. ) THEN
           ZLIMIT(1)=CI*TS_REF*ZDLEV/(DTSTEP*RHOI)*0.05 ! 5% is maximum change of layer entropy with reference temperature TS_REF 
           ZLIMIT(2)=ABS(FLUX_CONDUCT)
           ZMIN=MINVAL(ZLIMIT)
           FLUX_CONDUCT=SIGN(ZMIN,FLUX_CONDUCT)
         ENDIF
         FLUX_SNOW_ENTROPY(I,J,KH)=FLUX_CONDUCT
       ENDDO

! Entropy flux at separating surface between the snow and the soil:
! Entropy flux in snow with mass axis using has unit (J*m)/(kg*s*K),
! but entropy flux at the surface and in soil has unit J/(K*s*m**2)
  
       FLUX_ENTROPY_SNOW_SOIL(I,J)=FLUX_SNOW_ENTROPY(I,J,NLEV+1)*RHO_TOT(NLEV+1)
  
! Snow temperature after conductive thermal transport simulation at the full-levels
  
       DO K=0,NLEV
         KH=K 
         DS(K)=DTSTEP*RHO_TOT(K)*(-FLUX_SNOW_ENTROPY(I,J,KH+1)+FLUX_SNOW_ENTROPY(I,J,KH))/ZD_LEV_SNOW(KH)
         S_FIN(K)=S_INI(K)+DS(K)
         TSTAR(K)=TSNOW(I,J,K)+DS(K)/C_TOT(K)
       ENDDO

! Conservation of column entropy

       ERROR=SUM(DS(0:NLEV)*ZD_LEV_SNOW(0:NLEV))
       ERROR=ERROR+(-FLUX_SNOW_ENTROPY(I,J,0)*RHO_TOT(0)+FLUX_SNOW_ENTROPY(I,J,NLEV+1)*RHO_TOT(NLEV+1))*DTSTEP

       IF (ERROR /= 0. ) THEN

         Z01=SUM(ZD_LEV_SNOW(0:NLEV))
         DO K=0,NLEV
           WEIGHT=(ZD_LEV_SNOW(K))/Z01
           Z02=ZD_LEV_SNOW(K)*C_TOT(K)
           TSTAR(K)=MIN((TSTAR(K)*Z02-ERROR*WEIGHT)/Z02, 0.02)
           DS(K)=TSTAR(K)*C_TOT(K)-S_INI(K)
         ENDDO 

       ENDIF

! ------- Simulation of snow water phase changing --------

       DO K=0,NLEV

         IF (TSTAR(K) <= 0..AND.FICE_SNOW(I,J,K) == 1.) THEN
           TSNOW(I,J,K)=TSTAR(K)
           CYCLE
         ENDIF

         SNOW_ENTROPY(K)=ZD_LEV_SNOW(K)*(FICE_SNOW(I,J,K)*CI*TSTAR(K)+(1.-FICE_SNOW(I,J,K))*(CW*TSTAR(K)+LIWDT0))

         FICE_NEW=1.
         ZT_NEW=SNOW_ENTROPY(K)/(ZD_LEV_SNOW(K)*CI)

         IF (ZT_NEW < 0.) THEN
! Solution is found

           TSNOW(I,J,K)=ZT_NEW
           FICE_SNOW(I,J,K)=FICE_NEW

         ELSE

! Snow is melting

           TSNOW(I,J,K)=0.
           FICE_NEW=1.-SNOW_ENTROPY(K)/(ZD_LEV_SNOW(K)*LIWDT0)
           FICE_NEW=MAX( MIN( FICE_NEW, 1.), 1.E-8)
           FICE_SNOW(I,J,K)=FICE_NEW
           SNOW_MELT_AGE(I,J,K)=SNOW_MELT_AGE(I,J,K)+DTSTEP/86400.

         ENDIF

       ENDDO

     ELSE ! No entropy processes scheme in snow
  
       FLUX_ENTROPY_SNOW_SOIL(I,J)=0.
       TSNOW(I,J,0)=MIN(TS(I,J,0), -1.E-3)
       TSNOW(I,J,1)=TS(I,J,0)
       FICE_SNOW(I,J,0)=FICE_SNOW(I,J,1) 
  
     ENDIF ! SNOW >= SNOW_MIN

     DO K=0,NLEV
       SNOW_AGE(I,J,K)=SNOW_AGE(I,J,K)+DTSTEP/86400.
     ENDDO

   ELSE ! NLEV < 0

     FLUX_ENTROPY_SNOW_SOIL(I,J)=0.

   ENDIF ! NLEV >= 0

! Run-off of melted water downward snow profile, and
! Renewing of snow mass configration due to the interaction with atmosphere
! and snow melting

   DSNOW_ATM_LIQ=FLUX_W_LIQ_ATM_SNOW(I,J)*DTSTEP
   DSNOW_ATM_ICE=(FLUX_WV_SNOW(I,J)+FLUX_W_ICE_ATM_SNOW(I,J))*DTSTEP
   DSNOW_ATM=DSNOW_ATM_LIQ+DSNOW_ATM_ICE

! Snow surface is refreshed by atmosferic precipitation:

   IF (SNOW_AGE(I,J,0) > 0..AND.FLUX_WV_SNOW(I,J)+FLUX_W_ICE_ATM_SNOW(I,J) > 1.E-4) THEN
     SNOW_AGE(I,J,0)=0.
     SNOW_MELT_AGE(I,J,0)=0.
   ENDIF

   FLUX_WATER_SNOW_SOIL(I,J)=0.

   IF (NLEV >= 0) THEN ! Snow exists at the time step start

     DO K=0,NLEV+1
       SNOW_WATER(K)=ZD_LEV_SNOW(K)*(1.-FICE_SNOW(I,J,K))
!if (PROC_IND==71.and.i==16.and.j==2) print *,'poisk 5',K,SNOW_WATER(K),ZD_LEV_SNOW(K),FICE_SNOW(I,J,K),ISTEP
     ENDDO
     SNOW_WATER_SUM=SUM( SNOW_WATER(0:NLEV+1) )
     SNOW_WATER(0)=SNOW_WATER(0)+DSNOW_ATM_LIQ
!if (PROC_IND==71.and.i==16.and.j==2) print *,'poisk 4',SNOW_WATER(8),"-",SNOW_WATER(:)

     DSNOW_MELT=MAX( SNOW_WATER(NLEV)+SNOW_WATER(NLEV+1), 0.)
     FLUX_WATER_SNOW_SOIL(I,J)=DSNOW_MELT/DTSTEP

     FLUX_SNOW_WATER(I,J,0)=FLUX_W_LIQ_ATM_SNOW(I,J)
     DO K=1,NLEV+1
       FLUX_SNOW_WATER(I,J,K)=SNOW_WATER(K-1)/DTSTEP
     ENDDO 

     SNOW_WATER_SUM=MAX( SNOW_WATER_SUM+DSNOW_ATM_LIQ-DSNOW_MELT, 0.)
     SNOW_NEW=MAX( SNOW_OLD+DSNOW_ATM-DSNOW_MELT, 0.)

     NLEV_NEW=NLEV

     IF (SNOW_NEW < DELTA_SNOW) SNOW_NEW=0.

     IF (SNOW_NEW > 0.) THEN ! Snow exists at the time step finish

       SNOW_WATER(NLEV:NLEV+1)=0. ! Water transits to the soil by melting snow water flux (see below)
!if (PROC_IND==71.and.i==16.and.j==2) print *,'poisk 3',SNOW_WATER(8),"-",SNOW_WATER(:)
! Melting water transits downward along the profile
       DO K=NLEV,1,-1
         SNOW_WATER(K)=SNOW_WATER(K-1)
       ENDDO
       SNOW_WATER(0)=0.
!if (PROC_IND==71.and.i==16.and.j==2) print *,'poisk 2',SNOW_WATER(8),"-",SNOW_WATER(:)

       LEV_SNOW(I,J,:)=VAL_MISSING
       LEV_SNOW(I,J,0)=0.
       IF (SNOW_NEW >= D_LEV_SNOW(I,J)+D_LEV_SNOW_MIN) THEN
         N=INT(SNOW_NEW/D_LEV_SNOW(I,J))
         IF (SNOW_NEW-D_LEV_SNOW(I,J)*FLOAT(N) >= D_LEV_SNOW_MIN) N=N+1
         DO K=1,N
           LEV_SNOW(I,J,K)=SNOW_NEW-FLOAT(N-K)*D_LEV_SNOW(I,J)
         ENDDO
         NLEV_NEW=N-1
       ELSE
         LEV_SNOW(I,J,1)=SNOW_NEW
         NLEV_NEW=0
       ENDIF

       ZD_LEV_SNOW(:)=D_LEV_SNOW(I,J)
       ZD_LEV_SNOW(0)=LEV_SNOW(I,J,1)*0.5
       ZD_LEV_SNOW(1)=ZD_LEV_SNOW(0)+D_LEV_SNOW(I,J)*0.5
       IF (NLEV_NEW > 0) THEN
         ZD_LEV_SNOW(NLEV_NEW+1)=D_LEV_SNOW(I,J)*0.5
       ELSE
         ZD_LEV_SNOW(NLEV_NEW+1)=LEV_SNOW(I,J,1)*0.5
       ENDIF

       IF (NLEV_NEW < NLEV) SNOW_WATER(NLEV_NEW+1)=SNOW_WATER(NLEV+1)+SNOW_WATER(NLEV)
!if (PROC_IND==71.and.i==16.and.j==2) print *,'poisk 1',SNOW_WATER(8),"-",SNOW_WATER(:)

       DO K=0,NLEV_NEW+1
! poisk
         FICE_SNOW(I,J,K)=1.
!         IF (ZD_LEV_SNOW(K) > 0.) THEN
if (isnan(SNOW_WATER(K)).or.isnan(ZD_LEV_SNOW(K))) then
write (703,*) SNOW_WATER(K),ZD_LEV_SNOW(K),K,I,J,PROC_IND
print *,'poisk 0',SNOW_WATER(K),ZD_LEV_SNOW(K),K,I,J,PROC_IND
endif
if (isnan(SNOW_WATER(K)).or.isnan(ZD_LEV_SNOW(K)).or.abs(SNOW_WATER(K))>1.e4.or.abs&
(ZD_LEV_SNOW(K))>1.e4.or.abs(ZD_LEV_SNOW(K))<1.e-8) then
write (702,*) SNOW_WATER(K),ZD_LEV_SNOW(K)
endif
           FICE_SNOW(I,J,K)=MIN( MAX( 1.-SNOW_WATER(K)/ZD_LEV_SNOW(K), 1.E-4), 1.)
!         ENDIF
       ENDDO

       IF (NLEV_NEW > NLEV) FICE_SNOW(I,J,0)=1.

       IF (NLEV_NEW /= NLEV) THEN ! Snow level number changes

         IF (NLEV_NEW > NLEV) THEN ! Snow levels number increased

! The upper inside level appears, temperature at air-snow surface changes to
! conserv the snow entropy at the upper 0.5 layer

           DO K=NLEV_NEW+1,1,-1
             TSNOW(I,J,K)=TSNOW(I,J,K-1)
             FICE_SNOW(I,J,K)=FICE_SNOW(I,J,K-1)
             SNOW_AGE(I,J,K)=SNOW_AGE(I,J,K-1)
             SNOW_MELT_AGE(I,J,K)=SNOW_MELT_AGE(I,J,K-1)
             RHOSNOW(I,J,K)=RHOSNOW(I,J,K-1)
           ENDDO

         ELSE ! Snow levels number decreased

! The lower inside level disappears, temperature at soil-snow surface changes to
! conserv the snow entropy at the lower 1.5 layer

           SNOW_ICE(NLEV:NLEV+1)=ZD_LEV_SNOW(NLEV:NLEV+1)*FICE_SNOW(I,J,NLEV:NLEV+1)
           SNOW_WATER(NLEV:NLEV+1)=ZD_LEV_SNOW(NLEV:NLEV+1)-SNOW_ICE(NLEV:NLEV+1)

           ZS=SNOW_ICE(NLEV+1)*CI*TSNOW(I,J,NLEV+1)+SNOW_ICE(NLEV)*CI*TSNOW(I,J,NLEV)+ &
 SNOW_WATER(NLEV+1)*(CW*TSNOW(I,J,NLEV+1)+LIWDT0)+SNOW_WATER(NLEV)*(CW*TSNOW(I,J,NLEV)+LIWDT0)
           ZNUMER=ZS-LIWDT0*(SNOW_WATER(NLEV+1)+SNOW_WATER(NLEV))
           ZDENOM=CI*(SNOW_ICE(NLEV+1)+SNOW_ICE(NLEV))+CW*(SNOW_WATER(NLEV+1)+SNOW_WATER(NLEV))
           TSNOW(I,J,NLEV_NEW+1)=MIN( ZNUMER/ZDENOM, 0.)
           TS(I,J,0)=TSNOW(I,J,NLEV_NEW+1)
           FICE_SNOW(I,J,NLEV_NEW+1)=(SNOW_ICE(NLEV+1)+SNOW_ICE(NLEV))/(ZD_LEV_SNOW(NLEV+1)+ZD_LEV_SNOW(NLEV))
           FICE_SNOW(I,J,NLEV_NEW+1)=MAX( MIN(FICE_SNOW(I,J,NLEV_NEW+1), 1.), 1.E-6)
           SNOW_AGE(I,J,NLEV_NEW+1)=SUM( SNOW_AGE(I,J,NLEV_NEW:NLEV_NEW+1) )*0.5
           SNOW_MELT_AGE(I,J,NLEV_NEW+1)=SUM( SNOW_MELT_AGE(I,J,NLEV_NEW:NLEV_NEW+1) )*0.5
           RHOSNOW(I,J,NLEV_NEW+1)=SUM( RHOSNOW(I,J,NLEV_NEW:NLEV_NEW+1) )*0.5

         ENDIF ! Snow level number increased or decreased

       ENDIF ! Snow level number changes

       IF (NLEV_NEW+1 > NLEV_SNOW-1.OR.(D_LEV_SNOW(I,J)>D_LEV_SNOW_BASE+0.1.AND.&
LEV_SNOW(I,J,NLEV_NEW+1) <= (D_LEV_SNOW(I,J)-D_LEV_SNOW_BASE)*FLOAT(NLEV_SNOW-1))) THEN
! It is necessary to change step along vertical axis and to redistribute snow leves
         NLEV=NLEV_NEW

         IF (NLEV_NEW+1 > NLEV_SNOW-1) THEN ! Snow cover if thick too
           D_LEV_SNOW(I,J)=D_LEV_SNOW(I,J)+D_LEV_SNOW_BASE
         ELSE ! Step along vertical axis may bee decreased
           D_LEV_SNOW(I,J)=D_LEV_SNOW(I,J)-D_LEV_SNOW_BASE
         ENDIF

         ZLEV_SNOW_NEW(:)=VAL_MISSING
         ZLEV_SNOW_NEW(0)=0.
         IF (LEV_SNOW(I,J,NLEV+1) >= D_LEV_SNOW(I,J)+D_LEV_SNOW_MIN) THEN
           N=INT(LEV_SNOW(I,J,NLEV+1)/D_LEV_SNOW(I,J))
           IF (LEV_SNOW(I,J,NLEV+1)-D_LEV_SNOW(I,J)*FLOAT(N) >= D_LEV_SNOW_MIN) N=N+1
           DO K=1,N
             ZLEV_SNOW_NEW(K)=LEV_SNOW(I,J,NLEV+1)-FLOAT(N-K)*D_LEV_SNOW(I,J)
           ENDDO
           NLEV_NEW=N-1
         ELSE
           ZLEV_SNOW_NEW(1)=LEV_SNOW(I,J,NLEV+1)
           NLEV_NEW=0
         ENDIF

         RHO_OLD(0:NLEV+1)=FICE_SNOW(I,J,0:NLEV+1)*RHOI+(1.-FICE_SNOW(I,J,0:NLEV+1))*RHOW
         SNOW_ENTROPY(0:NLEV+1)=RHO_OLD(0:NLEV+1)*(FICE_SNOW(I,J,0:NLEV+1)*CI*TSNOW(I,J,0:NLEV+1)+ &
 (1.-FICE_SNOW(I,J,0:NLEV+1))*(CW*TSNOW(I,J,0:NLEV+1)+LIWDT0))

! Interpolation along the vertical axis into new axis points with conservation of 
! integral of interpolating function

         CALL INTERP_CONSERV_1D(NLEV+2, LEV_SNOW(I,J,0:NLEV+1), SNOW_ENTROPY(0:NLEV+1), &
 NLEV_NEW+2, ZLEV_SNOW_NEW(0:NLEV_NEW+1), SNOW_ENTROPY_NEW(0:NLEV_NEW+1))

         IF (MINVAL(FICE_SNOW(I,J,0:NLEV+1)) < 0.9999) THEN
           CALL INTERP_CONSERV_1D(NLEV+2, LEV_SNOW(I,J,0:NLEV+1), FICE_SNOW(I,J,0:NLEV+1), &
 NLEV_NEW+2, ZLEV_SNOW_NEW(0:NLEV_NEW+1), ZFICE_NEW(0:NLEV_NEW+1))
           FICE_SNOW(I,J,0:NLEV_NEW+1)=MIN( MAX( ZFICE_NEW(0:NLEV_NEW+1), 0.01), 1.)
           CALL INTERP_CONSERV_1D(NLEV+2, LEV_SNOW(I,J,0:NLEV+1), RHO_OLD(0:NLEV+1), &
 NLEV_NEW+2, ZLEV_SNOW_NEW(0:NLEV_NEW+1), RHO_NEW(0:NLEV_NEW+1))
           RHO_NEW(0:NLEV_NEW+1)=MIN(MAX(RHO_NEW(0:NLEV_NEW+1), RHOSNOW_FRESH), RHOSNOW_FIRN)
         ELSE
           FICE_SNOW(I,J,0:NLEV_NEW+1)=1.
           RHO_NEW(0:NLEV_NEW+1)=RHOI
         ENDIF

         DO K=0,NLEV_NEW
           Z01=1.-FICE_SNOW(I,J,K)
           ZNUMER=SNOW_ENTROPY_NEW(K)/RHO_NEW(K)-LIWDT0*Z01
           ZDENOM=CI*FICE_SNOW(I,J,K)+CW*Z01
           TSNOW(I,J,K)=ZNUMER/ZDENOM
           IF (TSNOW(I,J,K) > 0.) THEN
             TSNOW(I,J,K)=0.
!!!             Z01=CW*TSNOW(I,J,K)+LIWDT0
!!!             ZNUMER=Z01-SNOW_ENTROPY_NEW(K)/RHO_NEW(K)
!!!             ZDENOM=Z01-CI*TSNOW(I,J,K)
!!!             FICE_SNOW(I,J,K)=ZNUMER/ZDENOM
             FICE_SNOW(I,J,K)=1.-SNOW_ENTROPY_NEW(K)/(RHO_NEW(K)*LIWDT0)
             FICE_SNOW(I,J,K)=MIN( MAX( FICE_SNOW(I,J,K), 0.01), 1.)
           ENDIF
         ENDDO

         F_OLD(0:NLEV+1)=SNOW_AGE(I,J,0:NLEV+1)
         CALL INTERP_1D(NLEV+2, LEV_SNOW(I,J,0:NLEV+1), F_OLD(0:NLEV+1), &
 NLEV_NEW+2, ZLEV_SNOW_NEW(0:NLEV_NEW+1), F_NEW(0:NLEV_NEW+1))
         SNOW_AGE(I,J,0:NLEV_NEW+1)=MAX(F_NEW(0:NLEV_NEW+1), 0.)

         F_OLD(0:NLEV+1)=SNOW_MELT_AGE(I,J,0:NLEV+1)
         CALL INTERP_1D(NLEV+2, LEV_SNOW(I,J,0:NLEV+1), F_OLD(0:NLEV+1), &
 NLEV_NEW+2, ZLEV_SNOW_NEW(0:NLEV_NEW+1), F_NEW(0:NLEV_NEW+1))
         SNOW_MELT_AGE(I,J,0:NLEV_NEW+1)=MAX(F_NEW(0:NLEV_NEW+1), 0.)

         F_OLD(0:NLEV+1)=RHOSNOW(I,J,0:NLEV+1)
         CALL INTERP_1D(NLEV+2, LEV_SNOW(I,J,0:NLEV+1), F_OLD(0:NLEV+1), &
 NLEV_NEW+2, ZLEV_SNOW_NEW(0:NLEV_NEW+1), F_NEW(0:NLEV_NEW+1))
         RHOSNOW(I,J,0:NLEV_NEW+1)=MIN(MAX(F_NEW(0:NLEV_NEW+1), RHOSNOW_FRESH), RHOSNOW_FIRN)

         LEV_SNOW(I,J,0:NLEV_NEW+1)=ZLEV_SNOW_NEW(0:NLEV_NEW+1)

       ENDIF ! Step along vertical axis changed

     ELSE ! Snow does not exist at the time step finish

       NLEV_NEW=-2

     ENDIF ! Snow at the time step finish

   ELSE ! Snow do not exist at the time step start

     SNOW_NEW=DSNOW_ATM

     IF (SNOW_NEW >= DELTA_SNOW) THEN ! Snow exists at the time step finish - snow appears

       D_LEV_SNOW(I,J)=D_LEV_SNOW_BASE
       NLEV_NEW=0
       LEV_SNOW(I,J,0)=0.
       LEV_SNOW(I,J,1)=SNOW_NEW
       TSNOW(I,J,0:1)=MIN (TS(I,J,0), 0.)
       RHOSNOW(I,J,0:1)=RHOSNOW_FRESH
       FICE_SNOW(I,J,0:1)=MAX( MIN( 1.-DSNOW_ATM_LIQ/DSNOW_ATM, 1.), 0.1)
       SNOW_AGE(I,J,0:1)=DTSTEP/86400.
       IF (FICE_SNOW(I,J,0) < 1.) THEN
         SNOW_MELT_AGE(I,J,0:1)=SNOW_AGE(I,J,0)
       ELSE
         SNOW_MELT_AGE(I,J,0:1)=0.
       ENDIF

     ELSE ! Snow do not exist at the time step finish

       NLEV_NEW=-2

     ENDIF ! Snow at the time step finish

   ENDIF ! Snow at the time step start

! Actual snow level number at time step finish

   DO K=NLEV_NEW+2,NLEV_SNOW
     LEV_SNOW(I,J,K)=VAL_MISSING
     TSNOW(I,J,K)=VAL_MISSING
     FICE_SNOW(I,J,K)=VAL_MISSING
     RHOSNOW(I,J,K)=VAL_MISSING
     SNOW_AGE(I,J,K)=0.
     SNOW_MELT_AGE(I,J,K)=0.
   ENDDO

! Snow density at the time step finish

   CALL SNOW_DENSITY(I, J, NLEV_NEW+1)

! Snow cover fraction and water content of the snow cover spreaded on the whole grid box 

   IF (NLEV_NEW >= 0) THEN
     SNOW_NEW=LEV_SNOW(I,J,NLEV_NEW+1)
     FSNCOVER(I,J)=MIN( ((SNOW_NEW/SNOW_FRMAX)**0.5), 1.0)
     SNOW(I,J)=SNOW_NEW
   ELSE
     SNOW(I,J)=0.
     FSNCOVER(I,J)=0.
   ENDIF

!do k=0,nlev_new
!if (tsnow(i,j,k)>1.E-6.or.tsnow(i,j,k)<-0.45) write (31,*) 'POCHVA_SOIL_SNOW finish oshibka tsnow',&
!!if (tsnow(i,j,k)>1.E-6.or.tsnow(i,j,k)<-1.20) write (31,*) 'POCHVA_SOIL_SNOW finish oshibka tsnow',&
! proc_ind,i,j,k,tsnow(i,j,k),exp(tsnow(i,j,k))*t0-t0,istep,lev_snow(i,j,k),nlev,nlev_new
!enddo
!do k=0,nlev_new+1
!if (fice_snow(i,j,k)>1.or.fice_snow(i,j,k)<0.1) write (31,*) 'POCHVA_SOIL_SNOW finish oshibka fice_snow',&
! proc_ind,i,j,k,fice_snow(i,j,k),istep
!enddo
!do k=0,nlev_soil_water(i,j)
!if (qs(i,j,k)>qsmax(i,j,k).or.qs(i,j,k)<qsmin(i,j,k)) write (31,*) 'POCHVA_SOIL_SNOW finish oshibka qs ',&
! proc_ind,i,j,k,qs(i,j,k),qsmax(i,j,k),qsmin(i,j,k),istep
!enddo

 ENDIF
 ENDDO
 ENDDO
!close(31)

RETURN
END SUBROUTINE POCHVA_SNOW
! ---------------------------------------------------------------------------------------------------------------------------------
!    Pochva:
!    Numerical scheme for the parametrization of thermo-hydrous processes
!    in soil and at its surface (own soil, vegetation, snow cover)
!
! Scheme is elaborated by O.Drofa, CNR-ISAC, Italy (o.drofa-at-isac.cnr.it)
!
! ---------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE POCHVA_INIT(NPOINT_X_EXT, NPOINT_Y_EXT, NLEV_SOIL_EXT, &
 LEV_SOIL_EXT, MASK_SOIL_EXT, TS_EXT, T_SURF_EXT, QS_EXT, SNOW, &
 QSOIL_MAX, QSOIL_MIN, CSOIL, RHOSOIL, PSISOIL, WCONDSOIL, PARBSOIL, PARCSOIL,&
 QSOIL_REL_WILT, QSOIL_REL_REF, LAIVEG, LAIVEG_MAX, FRACVEG, ROOTVEG, &
 IND_LEV_TBOTTOM, IND_LEV_QBOTTOM, VALMISS_EXT, DTIME, &
 TS_SURF_EXT, QS_SURF_EXT, FICE_SOIL_EXT, FICE_SOIL_SURF_EXT, &
 FLAG_SNOW_INIT, LEV_SNOW_EXT, TSNOW_EXT, FICE_SNOW_EXT, SNOW_AGE_EXT, SNOW_MELT_AGE_EXT, SNOW_DIRT_EXT)

! Procedure of initialition of all variables used in the scheme "Pochva"

!USE MODULE_POCHVA
IMPLICIT NONE

INTEGER :: NPOINT_X_EXT, NPOINT_Y_EXT, NLEV_SOIL_EXT
REAL, DIMENSION(:) :: LEV_SOIL_EXT
REAL, DIMENSION(:,:,:) :: TS_EXT, QS_EXT, &
 QSOIL_MAX, QSOIL_MIN, CSOIL, RHOSOIL, PSISOIL, WCONDSOIL, PARBSOIL, PARCSOIL,&
 QSOIL_REL_WILT, QSOIL_REL_REF
INTEGER, DIMENSION(:,:) :: MASK_SOIL_EXT
REAL, DIMENSION(:,:) :: T_SURF_EXT, SNOW, LAIVEG, FRACVEG, ROOTVEG
INTEGER, DIMENSION(:,:) :: IND_LEV_TBOTTOM, IND_LEV_QBOTTOM
REAL :: LAIVEG_MAX, VALMISS_EXT, DTIME

REAL, DIMENSION(:,:), OPTIONAL :: SNOW_DIRT_EXT, TS_SURF_EXT, QS_SURF_EXT, FICE_SOIL_SURF_EXT
REAL, DIMENSION(:,:,:), OPTIONAL :: FICE_SOIL_EXT
INTEGER, OPTIONAL :: FLAG_SNOW_INIT
REAL, DIMENSION(:,:,:), OPTIONAL :: LEV_SNOW_EXT, TSNOW_EXT, &
 FICE_SNOW_EXT, SNOW_AGE_EXT, SNOW_MELT_AGE_EXT

INTEGER :: NLEV, I, J, K, N, NN, FLAG_SNOW
REAL :: Z01, Z02, ZFICE, ZDFICE, ZQS

! Time step

 DTSTEP=DTIME

! Points number of the numerical domain

 NPOINT_X=NPOINT_X_EXT
 NPOINT_Y=NPOINT_Y_EXT

! Soil level number

 NLEV_SOIL_MAX=NLEV_SOIL_EXT

! Snow (maximum) level number

! NLEV_SNOW=10

! Allocations

 IF (ALLOCATED(NLEV_SOIL_HEAT)) DEALLOCATE(NLEV_SOIL_HEAT)
 ALLOCATE(NLEV_SOIL_HEAT(NPOINT_X,NPOINT_Y))
 NLEV_SOIL_HEAT(:,:)=0
 IF (ALLOCATED(NLEV_SOIL_WATER)) DEALLOCATE(NLEV_SOIL_WATER)
 ALLOCATE(NLEV_SOIL_WATER(NPOINT_X,NPOINT_Y))
 NLEV_SOIL_WATER(:,:)=0
 IF (ALLOCATED(LEV_SOIL)) DEALLOCATE(LEV_SOIL)
 ALLOCATE(LEV_SOIL(0:NLEV_SOIL_MAX))
 LEV_SOIL(:)=0.
 IF (ALLOCATED(LEV_SOIL_H)) DEALLOCATE(LEV_SOIL_H)
 ALLOCATE(LEV_SOIL_H(NLEV_SOIL_MAX))
 LEV_SOIL_H(:)=0.
 IF (ALLOCATED(D_LEV_SOIL)) DEALLOCATE(D_LEV_SOIL)
 ALLOCATE(D_LEV_SOIL(NLEV_SOIL_MAX))
 D_LEV_SOIL(:)=0.
 IF (ALLOCATED(D_LEV_SOIL_H)) DEALLOCATE(D_LEV_SOIL_H)
 ALLOCATE(D_LEV_SOIL_H(0:NLEV_SOIL_MAX-1))
 D_LEV_SOIL_H(:)=0.

 IF (ALLOCATED(TS)) DEALLOCATE(TS)
 ALLOCATE(TS(NPOINT_X,NPOINT_Y,0:NLEV_SOIL_MAX))
 TS(:,:,:)=0.
 IF (ALLOCATED(QS)) DEALLOCATE(QS)
 ALLOCATE(QS(NPOINT_X,NPOINT_Y,0:NLEV_SOIL_MAX))
 QS(:,:,:)=0.
 IF (ALLOCATED(QS_REL)) DEALLOCATE(QS_REL)
 ALLOCATE(QS_REL(NPOINT_X,NPOINT_Y,0:NLEV_SOIL_MAX))
 QS_REL(:,:,:)=0.
 IF (ALLOCATED(PSIS)) DEALLOCATE(PSIS)
 ALLOCATE(PSIS(NPOINT_X,NPOINT_Y,0:NLEV_SOIL_MAX))
 PSIS(:,:,:)=0.
 IF (ALLOCATED(KAPPA_H)) DEALLOCATE(KAPPA_H)
 ALLOCATE(KAPPA_H(NPOINT_X,NPOINT_Y,NLEV_SOIL_MAX))
 KAPPA_H(:,:,:)=0.
 IF (ALLOCATED(FRAC_SICE)) DEALLOCATE(FRAC_SICE)
 ALLOCATE(FRAC_SICE(NPOINT_X,NPOINT_Y,0:NLEV_SOIL_MAX))
 FRAC_SICE(:,:,:)=0.
 IF (ALLOCATED(DFRAC_SICE_DT)) DEALLOCATE(DFRAC_SICE_DT)
 ALLOCATE(DFRAC_SICE_DT(NPOINT_X,NPOINT_Y,0:NLEV_SOIL_MAX))
 DFRAC_SICE_DT(:,:,:)=0.
 IF (ALLOCATED(QSURF)) DEALLOCATE(QSURF)
 ALLOCATE(QSURF(NPOINT_X,NPOINT_Y))
 QSURF(:,:)=0.
 IF (ALLOCATED(TSURF)) DEALLOCATE(TSURF)
 ALLOCATE(TSURF(NPOINT_X,NPOINT_Y))
 TSURF(:,:)=0.
 IF (ALLOCATED(TFOREST)) DEALLOCATE(TFOREST)
 ALLOCATE(TFOREST(NPOINT_X,NPOINT_Y))
 TFOREST(:,:)=0.
 IF (ALLOCATED(QFOREST)) DEALLOCATE(QFOREST)
 ALLOCATE(QFOREST(NPOINT_X,NPOINT_Y))
 QFOREST(:,:)=0.
 IF (ALLOCATED(QS_REL_VEG)) DEALLOCATE(QS_REL_VEG)
 ALLOCATE(QS_REL_VEG(NPOINT_X,NPOINT_Y,0:NLEV_SOIL_MAX))
 QS_REL_VEG(:,:,:)=0.

 IF (ALLOCATED(D_LEV_SNOW)) DEALLOCATE(D_LEV_SNOW)
 ALLOCATE(D_LEV_SNOW(NPOINT_X,NPOINT_Y))
 D_LEV_SNOW(:,:)=D_LEV_SNOW_BASE ! kg/m^2
 IF (ALLOCATED(TSNOW)) DEALLOCATE(TSNOW)
 ALLOCATE(TSNOW(NPOINT_X,NPOINT_Y,0:NLEV_SNOW))
 TSNOW(:,:,:)=-1.E-2
 IF (ALLOCATED(LEV_SNOW)) DEALLOCATE(LEV_SNOW)
 ALLOCATE(LEV_SNOW(NPOINT_X,NPOINT_Y,0:NLEV_SNOW))
 LEV_SNOW(:,:,:)=VAL_MISSING
 IF (ALLOCATED(FICE_SNOW)) DEALLOCATE(FICE_SNOW)
 ALLOCATE(FICE_SNOW(NPOINT_X,NPOINT_Y,0:NLEV_SNOW))
 FICE_SNOW(:,:,:)=1.
 IF (ALLOCATED(RHOSNOW)) DEALLOCATE(RHOSNOW)
 ALLOCATE(RHOSNOW(NPOINT_X,NPOINT_Y,0:NLEV_SNOW))
 RHOSNOW(:,:,:)=MIN( RHOSNOW_FRESH*2., RHOSNOW_OLD )
 IF (ALLOCATED(SNOW_AGE)) DEALLOCATE(SNOW_AGE)
 ALLOCATE(SNOW_AGE(NPOINT_X,NPOINT_Y,0:NLEV_SNOW))
 SNOW_AGE(:,:,:)=0.
 IF (ALLOCATED(SNOW_MELT_AGE)) DEALLOCATE(SNOW_MELT_AGE)
 ALLOCATE(SNOW_MELT_AGE(NPOINT_X,NPOINT_Y,0:NLEV_SNOW))
 SNOW_MELT_AGE(:,:,:)=0.
 IF (ALLOCATED(DLEV_SNOW_BOTTOM)) DEALLOCATE(DLEV_SNOW_BOTTOM)
 ALLOCATE(DLEV_SNOW_BOTTOM(NPOINT_X,NPOINT_Y))
 DLEV_SNOW_BOTTOM(:,:)=0.
 IF (ALLOCATED(FLAG_SNOW_THICK)) DEALLOCATE(FLAG_SNOW_THICK)
 ALLOCATE(FLAG_SNOW_THICK(NPOINT_X,NPOINT_Y))
 FLAG_SNOW_THICK(:,:)=0
 IF (ALLOCATED(SNOW_DIRT)) DEALLOCATE(SNOW_DIRT)
 ALLOCATE(SNOW_DIRT(NPOINT_X,NPOINT_Y))
 SNOW_DIRT(:,:)=0.5

 IF (ALLOCATED(QS_OLD)) DEALLOCATE(QS_OLD)
 ALLOCATE(QS_OLD(NPOINT_X,NPOINT_Y,0:NLEV_SOIL_MAX))
 QS_OLD(:,:,:)=0.
 IF (ALLOCATED(QSI)) DEALLOCATE(QSI)
 ALLOCATE(QSI(NPOINT_X,NPOINT_Y,0:NLEV_SOIL_MAX))
 QSI(:,:,:)=0.
 IF (ALLOCATED(DQS)) DEALLOCATE(DQS)
 ALLOCATE(DQS(NPOINT_X,NPOINT_Y,0:NLEV_SOIL_MAX))
 DQS(:,:,:)=0.
 IF (ALLOCATED(RUN_OFF)) DEALLOCATE(RUN_OFF)
 ALLOCATE(RUN_OFF(NPOINT_X,NPOINT_Y,0:NLEV_SOIL_MAX))
 RUN_OFF(:,:,:)=0.
 IF (ALLOCATED(FLUX_ENTROPY_ATM)) DEALLOCATE(FLUX_ENTROPY_ATM)
 ALLOCATE(FLUX_ENTROPY_ATM(NPOINT_X,NPOINT_Y))
 FLUX_ENTROPY_ATM(:,:)=0.
 IF (ALLOCATED(FLUX_ENTROPY_ATM_SPEC)) DEALLOCATE(FLUX_ENTROPY_ATM_SPEC)
 ALLOCATE(FLUX_ENTROPY_ATM_SPEC(NPOINT_X,NPOINT_Y))
 FLUX_ENTROPY_ATM_SPEC(:,:)=0.
 IF (ALLOCATED(FLUX_ENTROPY_ATM_LAT)) DEALLOCATE(FLUX_ENTROPY_ATM_LAT)
 ALLOCATE(FLUX_ENTROPY_ATM_LAT(NPOINT_X,NPOINT_Y))
 FLUX_ENTROPY_ATM_LAT(:,:)=0.
 IF (ALLOCATED(FLUX_ENTROPY_ATM_RAD)) DEALLOCATE(FLUX_ENTROPY_ATM_RAD)
 ALLOCATE(FLUX_ENTROPY_ATM_RAD(NPOINT_X,NPOINT_Y))
 FLUX_ENTROPY_ATM_RAD(:,:)=0.
 IF (ALLOCATED(FLUX_ENTROPY_ATM_SOIL)) DEALLOCATE(FLUX_ENTROPY_ATM_SOIL)
 ALLOCATE(FLUX_ENTROPY_ATM_SOIL(NPOINT_X,NPOINT_Y))
 FLUX_ENTROPY_ATM_SOIL(:,:)=0.
 IF (ALLOCATED(FLUX_ENTROPY_ATM_SNOW)) DEALLOCATE(FLUX_ENTROPY_ATM_SNOW)
 ALLOCATE(FLUX_ENTROPY_ATM_SNOW(NPOINT_X,NPOINT_Y))
 FLUX_ENTROPY_ATM_SNOW(:,:)=0.
 IF (ALLOCATED(FLUX_ENTROPY_ATM_FOREST)) DEALLOCATE(FLUX_ENTROPY_ATM_FOREST)
 ALLOCATE(FLUX_ENTROPY_ATM_FOREST(NPOINT_X,NPOINT_Y))
 FLUX_ENTROPY_ATM_FOREST(:,:)=0.
 IF (ALLOCATED(FLUX_ENTROPY_SNOW_SOIL)) DEALLOCATE(FLUX_ENTROPY_SNOW_SOIL)
 ALLOCATE(FLUX_ENTROPY_SNOW_SOIL(NPOINT_X,NPOINT_Y))
 FLUX_ENTROPY_SNOW_SOIL(:,:)=0.
 IF (ALLOCATED(FLUX_ENTROPY_FOREST_SOIL)) DEALLOCATE(FLUX_ENTROPY_FOREST_SOIL)
 ALLOCATE(FLUX_ENTROPY_FOREST_SOIL(NPOINT_X,NPOINT_Y))
 FLUX_ENTROPY_FOREST_SOIL(:,:)=0.
 IF (ALLOCATED(FLUX_WV)) DEALLOCATE(FLUX_WV)
 ALLOCATE(FLUX_WV(NPOINT_X,NPOINT_Y))
 FLUX_WV(:,:)=0.
 IF (ALLOCATED(FLUX_WV_SOIL)) DEALLOCATE(FLUX_WV_SOIL)
 ALLOCATE(FLUX_WV_SOIL(NPOINT_X,NPOINT_Y))
 FLUX_WV_SOIL(:,:)=0.
 IF (ALLOCATED(FLUX_WV_SNOW)) DEALLOCATE(FLUX_WV_SNOW)
 ALLOCATE(FLUX_WV_SNOW(NPOINT_X,NPOINT_Y))
 FLUX_WV_SNOW(:,:)=0.
 IF (ALLOCATED(FLUX_WV_VEGLOW_DRY)) DEALLOCATE(FLUX_WV_VEGLOW_DRY)
 ALLOCATE(FLUX_WV_VEGLOW_DRY(NPOINT_X,NPOINT_Y))
 FLUX_WV_VEGLOW_DRY(:,:)=0.
 IF (ALLOCATED(FLUX_WV_VEGLOW_WET)) DEALLOCATE(FLUX_WV_VEGLOW_WET)
 ALLOCATE(FLUX_WV_VEGLOW_WET(NPOINT_X,NPOINT_Y))
 FLUX_WV_VEGLOW_WET(:,:)=0.
 IF (ALLOCATED(FLUX_WV_FOREST_DRY)) DEALLOCATE(FLUX_WV_FOREST_DRY)
 ALLOCATE(FLUX_WV_FOREST_DRY(NPOINT_X,NPOINT_Y))
 FLUX_WV_FOREST_DRY(:,:)=0.
 IF (ALLOCATED(FLUX_WV_FOREST_WET)) DEALLOCATE(FLUX_WV_FOREST_WET)
 ALLOCATE(FLUX_WV_FOREST_WET(NPOINT_X,NPOINT_Y))
 FLUX_WV_FOREST_WET(:,:)=0.
 IF (ALLOCATED(FLUX_PREC_LIQ)) DEALLOCATE(FLUX_PREC_LIQ)
 ALLOCATE(FLUX_PREC_LIQ(NPOINT_X,NPOINT_Y))
 FLUX_PREC_LIQ(:,:)=0.
 IF (ALLOCATED(FLUX_PREC_ICE)) DEALLOCATE(FLUX_PREC_ICE)
 ALLOCATE(FLUX_PREC_ICE(NPOINT_X,NPOINT_Y))
 FLUX_PREC_ICE(:,:)=0.
 IF (ALLOCATED(FLUX_W_LIQ_ATM_SOIL)) DEALLOCATE(FLUX_W_LIQ_ATM_SOIL)
 ALLOCATE(FLUX_W_LIQ_ATM_SOIL(NPOINT_X,NPOINT_Y))
 FLUX_W_LIQ_ATM_SOIL(:,:)=0.
 IF (ALLOCATED(FLUX_W_LIQ_ATM_SNOW)) DEALLOCATE(FLUX_W_LIQ_ATM_SNOW)
 ALLOCATE(FLUX_W_LIQ_ATM_SNOW(NPOINT_X,NPOINT_Y))
 FLUX_W_LIQ_ATM_SNOW(:,:)=0.
 IF (ALLOCATED(FLUX_W_ICE_ATM_SNOW)) DEALLOCATE(FLUX_W_ICE_ATM_SNOW)
 ALLOCATE(FLUX_W_ICE_ATM_SNOW(NPOINT_X,NPOINT_Y))
 FLUX_W_ICE_ATM_SNOW(:,:)=0.
 IF (ALLOCATED(FLUX_W_LIQ_ATM_VEGLOW)) DEALLOCATE(FLUX_W_LIQ_ATM_VEGLOW)
 ALLOCATE(FLUX_W_LIQ_ATM_VEGLOW(NPOINT_X,NPOINT_Y))
 FLUX_W_LIQ_ATM_VEGLOW(:,:)=0.
 IF (ALLOCATED(FLUX_W_LIQ_ATM_FOREST)) DEALLOCATE(FLUX_W_LIQ_ATM_FOREST)
 ALLOCATE(FLUX_W_LIQ_ATM_FOREST(NPOINT_X,NPOINT_Y))
 FLUX_W_LIQ_ATM_FOREST(:,:)=0.
 IF (ALLOCATED(FLUX_WATER_SNOW_SOIL)) DEALLOCATE(FLUX_WATER_SNOW_SOIL)
 ALLOCATE(FLUX_WATER_SNOW_SOIL(NPOINT_X,NPOINT_Y))
 FLUX_WATER_SNOW_SOIL(:,:)=0.
 IF (ALLOCATED(FLUX_SOIL_ENTROPY)) DEALLOCATE(FLUX_SOIL_ENTROPY)
 ALLOCATE(FLUX_SOIL_ENTROPY(NPOINT_X,NPOINT_Y,0:NLEV_SOIL_MAX))
 FLUX_SOIL_ENTROPY(:,:,:)=0.
 IF (ALLOCATED(FLUX_SOIL_WATER)) DEALLOCATE(FLUX_SOIL_WATER)
 ALLOCATE(FLUX_SOIL_WATER(NPOINT_X,NPOINT_Y,0:NLEV_SOIL_MAX))
 FLUX_SOIL_WATER(:,:,:)=0.
 IF (ALLOCATED(FLUX_SOIL_WATER_VEG)) DEALLOCATE(FLUX_SOIL_WATER_VEG)
 ALLOCATE(FLUX_SOIL_WATER_VEG(NPOINT_X,NPOINT_Y,NLEV_SOIL_MAX))
 FLUX_SOIL_WATER_VEG(:,:,:)=0.
 IF (ALLOCATED(FLUX_SOIL_WATER_BOTTOM_DIAG)) DEALLOCATE(FLUX_SOIL_WATER_BOTTOM_DIAG)
 ALLOCATE(FLUX_SOIL_WATER_BOTTOM_DIAG(NPOINT_X,NPOINT_Y))
 FLUX_SOIL_WATER_BOTTOM_DIAG(:,:)=0.
 IF (ALLOCATED(LAMBDAG)) DEALLOCATE(LAMBDAG)
 ALLOCATE(LAMBDAG(NPOINT_X,NPOINT_Y,0:NLEV_SOIL_MAX))
 LAMBDAG(:,:,:)=0.
 IF (ALLOCATED(FLUX_SNOW_ENTROPY)) DEALLOCATE(FLUX_SNOW_ENTROPY)
 ALLOCATE(FLUX_SNOW_ENTROPY(NPOINT_X,NPOINT_Y,0:NLEV_SNOW))
 FLUX_SNOW_ENTROPY(:,:,:)=0.
 IF (ALLOCATED(FLUX_SNOW_WATER)) DEALLOCATE(FLUX_SNOW_WATER)
 ALLOCATE(FLUX_SNOW_WATER(NPOINT_X,NPOINT_Y,0:NLEV_SNOW))
 FLUX_SNOW_WATER(:,:,:)=0.

 IF (ALLOCATED(QSMAX)) DEALLOCATE(QSMAX)
 ALLOCATE(QSMAX(NPOINT_X,NPOINT_Y,0:NLEV_SOIL_MAX))
 QSMAX(:,:,:)=0.
 IF (ALLOCATED(QSMIN)) DEALLOCATE(QSMIN)
 ALLOCATE(QSMIN(NPOINT_X,NPOINT_Y,0:NLEV_SOIL_MAX))
 QSMIN(:,:,:)=0.
 IF (ALLOCATED(CG)) DEALLOCATE(CG)
 ALLOCATE(CG(NPOINT_X,NPOINT_Y,0:NLEV_SOIL_MAX))
 CG(:,:,:)=0.
 IF (ALLOCATED(RHOG)) DEALLOCATE(RHOG)
 ALLOCATE(RHOG(NPOINT_X,NPOINT_Y,0:NLEV_SOIL_MAX))
 RHOG(:,:,:)=0.
 IF (ALLOCATED(PSIS_SAT)) DEALLOCATE(PSIS_SAT)
 ALLOCATE(PSIS_SAT(NPOINT_X,NPOINT_Y,0:NLEV_SOIL_MAX))
 PSIS_SAT(:,:,:)=0.
 IF (ALLOCATED(KAPPA_SAT)) DEALLOCATE(KAPPA_SAT)
 ALLOCATE(KAPPA_SAT(NPOINT_X,NPOINT_Y,0:NLEV_SOIL_MAX))
 KAPPA_SAT(:,:,:)=0.
 IF (ALLOCATED(PAR_B)) DEALLOCATE(PAR_B)
 ALLOCATE(PAR_B(NPOINT_X,NPOINT_Y,0:NLEV_SOIL_MAX))
 PAR_B(:,:,:)=0.
 IF (ALLOCATED(PAR_C)) DEALLOCATE(PAR_C)
 ALLOCATE(PAR_C(NPOINT_X,NPOINT_Y,0:NLEV_SOIL_MAX))
 PAR_C(:,:,:)=0.
 IF (ALLOCATED(PAR_FB)) DEALLOCATE(PAR_FB)
 ALLOCATE(PAR_FB(NPOINT_X,NPOINT_Y,0:NLEV_SOIL_MAX))
 PAR_FB(:,:,:)=0.
 IF (ALLOCATED(QS_REL_VEGWILT)) DEALLOCATE(QS_REL_VEGWILT)
 ALLOCATE(QS_REL_VEGWILT(NPOINT_X,NPOINT_Y,0:NLEV_SOIL_MAX))
 QS_REL_VEGWILT(:,:,:)=0.
 IF (ALLOCATED(QS_REL_VEGREF)) DEALLOCATE(QS_REL_VEGREF)
 ALLOCATE(QS_REL_VEGREF(NPOINT_X,NPOINT_Y,0:NLEV_SOIL_MAX))
 QS_REL_VEGREF(:,:,:)=0.

 IF (ALLOCATED(LAI_VEGLOW)) DEALLOCATE(LAI_VEGLOW)
 ALLOCATE(LAI_VEGLOW(NPOINT_X,NPOINT_Y))
 LAI_VEGLOW(:,:)=0.
 IF (ALLOCATED(LAI_FOREST)) DEALLOCATE(LAI_FOREST)
 ALLOCATE(LAI_FOREST(NPOINT_X,NPOINT_Y))
 LAI_FOREST(:,:)=0.
 IF (ALLOCATED(FVEGLOW)) DEALLOCATE(FVEGLOW)
 ALLOCATE(FVEGLOW(NPOINT_X,NPOINT_Y))
 FVEGLOW(:,:)=0.
 IF (ALLOCATED(FFOREST)) DEALLOCATE(FFOREST)
 ALLOCATE(FFOREST(NPOINT_X,NPOINT_Y))
 FFOREST(:,:)=0.
 IF (ALLOCATED(FVEGLOWWET)) DEALLOCATE(FVEGLOWWET)
 ALLOCATE(FVEGLOWWET(NPOINT_X,NPOINT_Y))
 FVEGLOWWET(:,:)=0.
 IF (ALLOCATED(FFORESTWET)) DEALLOCATE(FFORESTWET)
 ALLOCATE(FFORESTWET(NPOINT_X,NPOINT_Y))
 FFORESTWET(:,:)=0.
 IF (ALLOCATED(ROOT_DEPTH_VEGLOW)) DEALLOCATE(ROOT_DEPTH_VEGLOW)
 ALLOCATE(ROOT_DEPTH_VEGLOW(NPOINT_X,NPOINT_Y))
 ROOT_DEPTH_VEGLOW(:,:)=0.
 IF (ALLOCATED(QW_VEGLOW)) DEALLOCATE(QW_VEGLOW)
 ALLOCATE(QW_VEGLOW(NPOINT_X,NPOINT_Y))
 QW_VEGLOW(:,:)=0.
 IF (ALLOCATED(QW_FOREST)) DEALLOCATE(QW_FOREST)
 ALLOCATE(QW_FOREST(NPOINT_X,NPOINT_Y))
 QW_FOREST(:,:)=0.
 IF (ALLOCATED(ROOT_DEPTH_FOREST)) DEALLOCATE(ROOT_DEPTH_FOREST)
 ALLOCATE(ROOT_DEPTH_FOREST(NPOINT_X,NPOINT_Y))
 ROOT_DEPTH_FOREST(:,:)=0.
 IF (ALLOCATED(QW_VEGLOW_MAX)) DEALLOCATE(QW_VEGLOW_MAX)
 ALLOCATE(QW_VEGLOW_MAX(NPOINT_X,NPOINT_Y))
 QW_VEGLOW_MAX(:,:)=0.
 IF (ALLOCATED(QW_FOREST_MAX)) DEALLOCATE(QW_FOREST_MAX)
 ALLOCATE(QW_FOREST_MAX(NPOINT_X,NPOINT_Y))
 QW_FOREST_MAX(:,:)=0.
 IF (ALLOCATED(DZSUM_VEGLOW)) DEALLOCATE(DZSUM_VEGLOW)
 ALLOCATE(DZSUM_VEGLOW(NPOINT_X,NPOINT_Y))
 DZSUM_VEGLOW(:,:)=0.
 IF (ALLOCATED(DZSUM_FOREST)) DEALLOCATE(DZSUM_FOREST)
 ALLOCATE(DZSUM_FOREST(NPOINT_X,NPOINT_Y))
 DZSUM_FOREST(:,:)=0.

 IF (ALLOCATED(FSNCOVER)) DEALLOCATE(FSNCOVER)
 ALLOCATE(FSNCOVER(NPOINT_X,NPOINT_Y))
 FSNCOVER(:,:)=0.

 IF (ALLOCATED(DZ_DZSUM_VEGLOW)) DEALLOCATE(DZ_DZSUM_VEGLOW)
 ALLOCATE(DZ_DZSUM_VEGLOW(NPOINT_X,NPOINT_Y,0:NLEV_SOIL_MAX-1))
 DZ_DZSUM_VEGLOW(:,:,:)=0.
 IF (ALLOCATED(DZ_DZSUM_FOREST)) DEALLOCATE(DZ_DZSUM_FOREST)
 ALLOCATE(DZ_DZSUM_FOREST(NPOINT_X,NPOINT_Y,0:NLEV_SOIL_MAX-1))
 DZ_DZSUM_FOREST(:,:,:)=0.

 IF (ALLOCATED(MASK_SOIL)) DEALLOCATE(MASK_SOIL)
 ALLOCATE(MASK_SOIL(NPOINT_X,NPOINT_Y))
 MASK_SOIL(:,:)=0
 IF (ALLOCATED(MASK_SOIL_AXED)) DEALLOCATE(MASK_SOIL_AXED)
 ALLOCATE(MASK_SOIL_AXED(NPOINT_X,NPOINT_Y))
 MASK_SOIL_AXED(:,:)=0

! Definitions

! Full-levels depth (in m, reading from soil surface)

 LEV_SOIL(0)=0.
 LEV_SOIL(1:NLEV_SOIL_MAX)=LEV_SOIL_EXT(1:NLEV_SOIL_MAX)

! Half-levels depth (in m, reading from soil surface)

 DO K=1,NLEV_SOIL_MAX
  LEV_SOIL_H(K)=LEV_SOIL(K-1)+0.5*(LEV_SOIL(K)-LEV_SOIL(K-1))
 ENDDO

! Difference between full-levels depth (m)

 DO K=1,NLEV_SOIL_MAX
   D_LEV_SOIL(K)=LEV_SOIL(K)-LEV_SOIL(K-1)
 ENDDO

! Difference between half-levels depth (m)

 D_LEV_SOIL_H(0)=LEV_SOIL_H(1)
 DO K=1,NLEV_SOIL_MAX-1
   D_LEV_SOIL_H(K)=LEV_SOIL_H(K+1)-LEV_SOIL_H(K)
 ENDDO
 SUM_D_LEV_SOIL_H=0.
 DO K=0,NLEV_SOIL_MAX-1
   SUM_D_LEV_SOIL_H=SUM_D_LEV_SOIL_H+D_LEV_SOIL_H(K)
 ENDDO

! Constant (in time) soil and vegetation physical parameters

! Maximum specific volumetric soil water content (m^3/m^3) at whole soil levels (0:NLEV_SOIL_MAX)

 QSMAX(1:NPOINT_X,1:NPOINT_Y,1:NLEV_SOIL_MAX)=QSOIL_MAX(1:NPOINT_X,1:NPOINT_Y,1:NLEV_SOIL_MAX)
 QSMAX(1:NPOINT_X,1:NPOINT_Y,0)=QSMAX(1:NPOINT_X,1:NPOINT_Y,1)

! Minimum specific volumetric soil water content (m^3/m^3) at whole soil levels (0:NLEV_SOIL_MAX)

 QSMIN(1:NPOINT_X,1:NPOINT_Y,1:NLEV_SOIL_MAX)=QSOIL_MIN(1:NPOINT_X,1:NPOINT_Y,1:NLEV_SOIL_MAX)
 QSMIN(1:NPOINT_X,1:NPOINT_Y,0)=QSMIN(1:NPOINT_X,1:NPOINT_Y,1)

! Dry soil thermal capacity (J/K/m^3) at whole soil levels (0:NLEV_SOIL_MAX)
 CG(1:NPOINT_X,1:NPOINT_Y,1:NLEV_SOIL_MAX)=CSOIL(1:NPOINT_X,1:NPOINT_Y,1:NLEV_SOIL_MAX)
 CG(1:NPOINT_X,1:NPOINT_Y,0)=CG(1:NPOINT_X,1:NPOINT_Y,1)

! Dry soil density (kg/m^3) at whole soil levels (0:NLEV_SOIL_MAX)
 RHOG(1:NPOINT_X,1:NPOINT_Y,1:NLEV_SOIL_MAX)=RHOSOIL(1:NPOINT_X,1:NPOINT_Y,1:NLEV_SOIL_MAX)
 RHOG(1:NPOINT_X,1:NPOINT_Y,0)=RHOG(1:NPOINT_X,1:NPOINT_Y,1)

! Idraulic potential of saturated soil (m) at whole soil levels (0:NLEV_SOIL_MAX)
 PSIS_SAT(1:NPOINT_X,1:NPOINT_Y,1:NLEV_SOIL_MAX)=PSISOIL(1:NPOINT_X,1:NPOINT_Y,1:NLEV_SOIL_MAX)
 PSIS_SAT(1:NPOINT_X,1:NPOINT_Y,0)=PSIS_SAT(1:NPOINT_X,1:NPOINT_Y,1)

! Idraulic conductivity of saturated soil (m/s) at whole soil levels (0:NLEV_SOIL_MAX)
 KAPPA_SAT(1:NPOINT_X,1:NPOINT_Y,1:NLEV_SOIL_MAX)=WCONDSOIL(1:NPOINT_X,1:NPOINT_Y,1:NLEV_SOIL_MAX)
 KAPPA_SAT(1:NPOINT_X,1:NPOINT_Y,0)=WCONDSOIL(1:NPOINT_X,1:NPOINT_Y,1)

! Parameter "b" in soil water transport equation (dimensionless) at whole soil levels (0:NLEV_SOIL_MAX)
 PAR_B(1:NPOINT_X,1:NPOINT_Y,1:NLEV_SOIL_MAX)=PARBSOIL(1:NPOINT_X,1:NPOINT_Y,1:NLEV_SOIL_MAX)
 PAR_B(1:NPOINT_X,1:NPOINT_Y,0)=PAR_B(1:NPOINT_X,1:NPOINT_Y,1)

! Parameter "c" in soil water transport equation (dimensionless) at whole soil levels (0:NLEV_SOIL_MAX)
 PAR_C(1:NPOINT_X,1:NPOINT_Y,1:NLEV_SOIL_MAX)=PARCSOIL(1:NPOINT_X,1:NPOINT_Y,1:NLEV_SOIL_MAX)
 PAR_C(1:NPOINT_X,1:NPOINT_Y,0)=PAR_C(1:NPOINT_X,1:NPOINT_Y,1)

! Relative soil water content (proportion) at which the vegetation wilt at whole soil levels (0:NLEV_SOIL_MAX)
 QS_REL_VEGWILT(1:NPOINT_X,1:NPOINT_Y,1:NLEV_SOIL_MAX)=QSOIL_REL_WILT(1:NPOINT_X,1:NPOINT_Y,1:NLEV_SOIL_MAX)
 QS_REL_VEGWILT(1:NPOINT_X,1:NPOINT_Y,0)=QS_REL_VEGWILT(1:NPOINT_X,1:NPOINT_Y,1)

! Relative soil water content (proportion) at which the vegetation stop evapotraspiration increase at whole soil levels (0:NLEV_SOIL_MAX)
 QS_REL_VEGREF(1:NPOINT_X,1:NPOINT_Y,1:NLEV_SOIL_MAX)=QSOIL_REL_REF(1:NPOINT_X,1:NPOINT_Y,1:NLEV_SOIL_MAX)
 QS_REL_VEGREF(1:NPOINT_X,1:NPOINT_Y,0)=QS_REL_VEGREF(1:NPOINT_X,1:NPOINT_Y,1)

! Leaf Area Index of low vegetation (dimensionless)
 LAI_VEGLOW(1:NPOINT_X,1:NPOINT_Y)=LAIVEG(1:NPOINT_X,1:NPOINT_Y)
 LAI_MAX=LAIVEG_MAX

! Leaf Area Index of forest (high) vegetation (dimensionless)
 LAI_FOREST(1:NPOINT_X,1:NPOINT_Y)=0.

! Fraction of low vegetation (not forest) (proportion)
 FVEGLOW(1:NPOINT_X,1:NPOINT_Y)=FRACVEG(1:NPOINT_X,1:NPOINT_Y)

! Fraction of high vegetation (forest) (proportion)
 FFOREST(1:NPOINT_X,1:NPOINT_Y)=0.

! Fraction of wet leaf surface of low vegetation (not forest) (proportion)
 FVEGLOWWET(1:NPOINT_X,1:NPOINT_Y)=0.

! Fraction of wet leaf surface of high vegetation (forest) (proportion)
 FFORESTWET(1:NPOINT_X,1:NPOINT_Y)=0.

! Root depth of low vegetation (m)
 ROOT_DEPTH_VEGLOW(1:NPOINT_X,1:NPOINT_Y)=ROOTVEG(1:NPOINT_X,1:NPOINT_Y)

! Root depth of forest (high) vegetation (m)
 ROOT_DEPTH_FOREST(1:NPOINT_X,1:NPOINT_Y)=0.

! Water deposited on the leaf surface of the low vegetation (kg/m^2)
 QW_VEGLOW(1:NPOINT_X,1:NPOINT_Y)=0.

! Water deposited on the leaf surface of the hight (forest) vegetation (kg/m^2)
 QW_FOREST(1:NPOINT_X,1:NPOINT_Y)=0.

! Maximum value of water deposited on the leaf surface of the low vegetation (kg/m^2)
 IF (LAIVEG_MAX > 0.) THEN
!   QW_VEGLOW_MAX(1:NPOINT_X,1:NPOINT_Y)=LAI_VEGLOW(1:NPOINT_X,1:NPOINT_Y)/LAIVEG_MAX*10. ! 10 kg/m**2 at LAI of tropical forest Run4
   QW_VEGLOW_MAX(1:NPOINT_X,1:NPOINT_Y)=LAI_VEGLOW(1:NPOINT_X,1:NPOINT_Y)/LAIVEG_MAX*5. ! 5 kg/m**2 at LAI of tropical forest
 ELSE
   QW_VEGLOW_MAX(1:NPOINT_X,1:NPOINT_Y)=0.
 ENDIF 

! Maximum value of water deposited on the leaf surface of the hight (forest) vegetation (kg/m^2)
 QW_FOREST_MAX(1:NPOINT_X,1:NPOINT_Y)=0.

! Parameter for determination of soil ice fraction

 PAR_FA=-4./ALOG((T0+T_FICE_CRIT)/T0)
 DO J=1,NPOINT_Y
 DO I=1,NPOINT_X
 DO K=0,NLEV_SOIL_MAX
   Z01=MIN(MAX(PAR_B(I,J,K), 4.), 12.)
   Z02=(Z01-4.)/8.
   PAR_FB(I,J,K)=2.-Z02
 ENDDO
 ENDDO
 ENDDO

! Number of soil levels used in each domain point

 NLEV_SOIL_HEAT(:,:)=IND_LEV_TBOTTOM(:,:)
 NLEV_SOIL_WATER(:,:)=IND_LEV_QBOTTOM(:,:)
  
! For vegetation scheme: fraction of each soil layer with respect to summary soil layer affacted by vegetation evapotraspiration processe

 DZSUM_VEGLOW(:,:)=0.
 DZSUM_FOREST(:,:)=0.
 DO J=1,NPOINT_Y
 DO I=1,NPOINT_X
   DO K=0,NLEV_SOIL_MAX-1
     IF (LEV_SOIL(K) <= ROOT_DEPTH_VEGLOW(I,J)) DZSUM_VEGLOW(I,J)=DZSUM_VEGLOW(I,J)+D_LEV_SOIL_H(K)
     IF (LEV_SOIL(K) <= ROOT_DEPTH_FOREST(I,J)) DZSUM_FOREST(I,J)=DZSUM_FOREST(I,J)+D_LEV_SOIL_H(K)
   ENDDO
   IF (DZSUM_VEGLOW(I,J) > 0.) THEN
     DO K=0,NLEV_SOIL_MAX-1
       IF (LEV_SOIL(K) <= ROOT_DEPTH_VEGLOW(I,J)) THEN
         DZ_DZSUM_VEGLOW(I,J,K)=D_LEV_SOIL_H(K)/DZSUM_VEGLOW(I,J)
       ELSE
         DZ_DZSUM_VEGLOW(I,J,K)=0.
       ENDIF
     ENDDO
   ELSE
     DZ_DZSUM_VEGLOW(I,J,:)=0.
   ENDIF
   IF (DZSUM_FOREST(I,J) > 0.) THEN
     DO K=0,NLEV_SOIL_MAX-1
       IF (LEV_SOIL(K) <= ROOT_DEPTH_FOREST(I,J)) THEN
         DZ_DZSUM_FOREST(I,J,K)=D_LEV_SOIL_H(K)/DZSUM_FOREST(I,J)
       ELSE
         DZ_DZSUM_FOREST(I,J,K)=0.
       ENDIF
     ENDDO
   ELSE
     DZ_DZSUM_FOREST(I,J,:)=0.
   ENDIF
 ENDDO
 ENDDO

! Flag for application of the Pochva scheme (1 or 0)

 MASK_SOIL(:,:)=MASK_SOIL_EXT(:,:)

! Flag for application of the Pochva scheme in the axed version (1 or 0)

 MASK_SOIL_AXED(:,:)=0
 DO J=1,NPOINT_Y
 DO I=1,NPOINT_X
   IF (PSIS_SAT(I,J,0) > 0.) MASK_SOIL_AXED(I,J)=1
 ENDDO
 ENDDO

! Initialization of main prognostic variables

 DO J=1,NPOINT_Y
 DO I=1,NPOINT_X
 IF (MASK_SOIL(I,J)==1) THEN

! Definition of soil temperature at full-levels (K)

   NLEV=NLEV_SOIL_HEAT(I,J)
   TSURF(I,J)=T_SURF_EXT(I,J)
   IF (PRESENT(TS_SURF_EXT)) THEN
     TS(I,J,0)=ALOG(TS_SURF_EXT(I,J)/T0)
   ELSE
     TS(I,J,0)=ALOG(T_SURF_EXT(I,J)/T0)
   ENDIF
   TS(I,J,1:NLEV)=ALOG(TS_EXT(I,J,1:NLEV)/T0)
   TS(I,J,NLEV+1:NLEV_SOIL_MAX)=TS(I,J,NLEV)

! Definition of specific and relative soil water content at full-levels (m^3/m^3)

   NLEV=NLEV_SOIL_WATER(I,J)
   QS(I,J,1:NLEV)=MIN( MAX( QS_EXT(I,J,1:NLEV),QSMIN(I,J,1:NLEV)), QSMAX(I,J,1:NLEV))
   IF (PRESENT(QS_SURF_EXT)) THEN
     QS(I,J,0)=MIN( MAX( QS_SURF_EXT(I,J),QSMIN(I,J,0)), QSMAX(I,J,0))
   ELSE
     QS(I,J,0)=QS(I,J,1)
   ENDIF
   QS_REL(I,J,0:NLEV)=(QS(I,J,0:NLEV)-QSMIN(I,J,0:NLEV))/(QSMAX(I,J,0:NLEV)-QSMIN(I,J,0:NLEV))
   QS(I,J,NLEV+1:NLEV_SOIL_MAX)=QS(I,J,NLEV)
   QS_REL(I,J,NLEV+1:NLEV_SOIL_MAX)=QS_REL(I,J,NLEV)

! Definition of soil ice fraction and it's derivative respect to temperature at full-levels

   NLEV=NLEV_SOIL_HEAT(I,J)
   IF (PRESENT(FICE_SOIL_EXT)) THEN
     FRAC_SICE(I,J,1:NLEV)=FICE_SOIL_EXT(I,J,1:NLEV)
     DO K=1,NLEV
       QSI(I,J,K)=QS(I,J,K)*FRAC_SICE(I,J,K)
       CALL FRAC_SOIL_ICE_POINT(TS(I,J,K), ZFICE, ZDFICE, I, J, K)
       DFRAC_SICE_DT(I,J,K)=ZDFICE
     ENDDO
   ELSE
     DO K=1,NLEV
       CALL FRAC_SOIL_ICE_POINT(TS(I,J,K), ZFICE, ZDFICE, I, J, K)
       QSI(I,J,K)=QS(I,J,K)*ZFICE
       FRAC_SICE(I,J,K)=ZFICE
       DFRAC_SICE_DT(I,J,K)=ZDFICE
     ENDDO
   ENDIF

   IF (PRESENT(FICE_SOIL_SURF_EXT)) THEN
     FRAC_SICE(I,J,0)=FICE_SOIL_SURF_EXT(I,J)
     QSI(I,J,0)=QS(I,J,0)*FRAC_SICE(I,J,0)
     CALL FRAC_SOIL_ICE_POINT(TS(I,J,0), ZFICE, ZDFICE, I, J, 0)
     DFRAC_SICE_DT(I,J,0)=ZDFICE
   ELSE
     CALL FRAC_SOIL_ICE_POINT(TS(I,J,0), ZFICE, ZDFICE, I, J, 0)
     QSI(I,J,0)=QS(I,J,0)*ZFICE
     FRAC_SICE(I,J,0)=ZFICE
     DFRAC_SICE_DT(I,J,0)=ZDFICE
   ENDIF

   QSI(I,J,NLEV+1:NLEV_SOIL_MAX)=QSI(I,J,NLEV)
   FRAC_SICE(I,J,NLEV+1:NLEV_SOIL_MAX)=FRAC_SICE(I,J,NLEV)
   DFRAC_SICE_DT(I,J,NLEV+1:NLEV_SOIL_MAX)=DFRAC_SICE_DT(I,J,NLEV)

   IF (MASK_SOIL_AXED(I,J) == 1) THEN ! Glacier
     TS(I,J,0:NLEV_SOIL_MAX)=MIN(TS(I,J,0:NLEV_SOIL_MAX), -1.E-3)
     QSI(I,J,0:NLEV_SOIL_MAX)=0.
     FRAC_SICE(I,J,0:NLEV_SOIL_MAX)=0.
     DFRAC_SICE_DT(I,J,0:NLEV_SOIL_MAX)=0.
     QS(I,J,0:NLEV_SOIL_MAX)=0.
     QS_REL(I,J,0:NLEV_SOIL_MAX)=0.
   ENDIF

! Definition of snow cover variables

   FSNCOVER(I,J)=MIN( ((SNOW(I,J)/SNOW_FRMAX)**0.5), 1.0)

   TSNOW(I,J,:)=VAL_MISSING
   FICE_SNOW(I,J,:)=VAL_MISSING
   RHOSNOW(I,J,:)=VAL_MISSING
   SNOW_AGE(I,J,:)=0.
   SNOW_MELT_AGE(I,J,:)=0.

   IF (PRESENT(FLAG_SNOW_INIT)) THEN
     FLAG_SNOW=FLAG_SNOW_INIT
   ELSE
     FLAG_SNOW=0
   ENDIF

   IF (PRESENT(LEV_SNOW_EXT).AND.FLAG_SNOW == 1) THEN

     DO K=0,NLEV_SNOW
       IF (INT(LEV_SNOW_EXT(I,J,K+1)) /= INT(VALMISS_EXT)) THEN
         LEV_SNOW(I,J,K)=LEV_SNOW_EXT(I,J,K+1)
       ELSE
         EXIT
       ENDIF
     ENDDO
     N=K-2
     D_LEV_SNOW(I,J)=D_LEV_SNOW_BASE
     IF (N > 1) THEN
       NN=NINT((LEV_SNOW(I,J,2)-LEV_SNOW(I,J,1))/D_LEV_SNOW_BASE)
       D_LEV_SNOW(I,J)=D_LEV_SNOW_BASE*FLOAT(NN)
     ENDIF

   ELSE

     D_LEV_SNOW(I,J)=D_LEV_SNOW_BASE
     DO WHILE (.TRUE.)
       IF (D_LEV_SNOW(I,J)*FLOAT(NLEV_SNOW) >= SNOW(I,J)) EXIT
       D_LEV_SNOW(I,J)=D_LEV_SNOW(I,J)+D_LEV_SNOW_BASE
     ENDDO
  
     IF (SNOW(I,J) > 0.) THEN
       LEV_SNOW(I,J,0)=0.
       IF (SNOW(I,J) >= D_LEV_SNOW(I,J)+D_LEV_SNOW_MIN) THEN  
         N=INT(SNOW(I,J)/D_LEV_SNOW(I,J))
         IF (SNOW(I,J)-D_LEV_SNOW(I,J)*FLOAT(N) >= D_LEV_SNOW_MIN) N=N+1
         DO K=1,N
           LEV_SNOW(I,J,K)=SNOW(I,J)-FLOAT(N-K)*D_LEV_SNOW(I,J)
         ENDDO
         N=N-1
       ELSE
         N=0
         LEV_SNOW(I,J,1)=SNOW(I,J)
       ENDIF
     ELSE
       N=-1
     ENDIF

   ENDIF

   IF (N >= 0) THEN

     IF (PRESENT(TSNOW_EXT).AND.FLAG_SNOW == 1) THEN
       TSNOW(I,J,0:N)=ALOG(TSNOW_EXT(I,J,1:N+1)/T0)
       TSNOW(I,J,N+1)=MIN (TS(I,J,0), 0.)
     ELSE
       IF (N >= 1 ) THEN
         TSNOW(I,J,0)=MIN(ALOG(TSURF(I,J)/T0), -1.E-3)
         TSNOW(I,J,N+1)=MIN(ALOG(TSURF(I,J)/T0), -1.E-3)
         Z01=1./FLOAT(N+1)
         DO K=1,N
           TSNOW(I,J,K)=TSNOW(I,J,0)*(1.-Z01*FLOAT(K))+TSNOW(I,J,N+1)*Z01*FLOAT(K)
         ENDDO
       ELSEIF (N == 0 ) THEN
         TSNOW(I,J,0)=MIN(ALOG(TSURF(I,J)/T0), -1.E-3)
         TSNOW(I,J,1)=MIN(TS(I,J,0), -1.E-3)
       ENDIF
     ENDIF

     IF (PRESENT(FICE_SNOW_EXT).AND.FLAG_SNOW == 1) THEN
       FICE_SNOW(I,J,0:N+1)=FICE_SNOW_EXT(I,J,1:N+2)
     ELSE
       FICE_SNOW(I,J,0:N+1)=1.
     ENDIF
  
     IF (PRESENT(SNOW_AGE_EXT).AND.FLAG_SNOW == 1) THEN
       SNOW_AGE(I,J,0:N+1)=SNOW_AGE_EXT(I,J,1:N+2)
     ELSE
       SNOW_AGE(I,J,0:N+1)=10.
     ENDIF
  
     IF (PRESENT(SNOW_MELT_AGE_EXT).AND.FLAG_SNOW == 1) THEN
       SNOW_MELT_AGE(I,J,0:N+1)=SNOW_MELT_AGE_EXT(I,J,1:N+2)
     ELSE
       SNOW_MELT_AGE(I,J,0:N+1)=0.
     ENDIF
  
     RHOSNOW(I,J,0:N+1)=MIN( RHOSNOW_FRESH*2., RHOSNOW_OLD )

   ENDIF

   IF (N >= 0) THEN
     IF (LEV_SNOW(I,J,N+1) >= D_LEV_SNOW_MIN) THEN
       FLAG_SNOW_THICK(I,J)=1
       DLEV_SNOW_BOTTOM(I,J)=(LEV_SNOW(I,J,N+1)-LEV_SNOW(I,J,N))*0.5
     ELSE
       FLAG_SNOW_THICK(I,J)=0
       DLEV_SNOW_BOTTOM(I,J)=LEV_SNOW(I,J,1)
     ENDIF
   ELSE
     DLEV_SNOW_BOTTOM(I,J)=0.
   ENDIF

   IF (PRESENT(SNOW_DIRT_EXT)) THEN
     SNOW_DIRT(I,J)=SNOW_DIRT_EXT(I,J)
   ELSE
     SNOW_DIRT(I,J)=0.5
   ENDIF

! Definition of forest vegetation variables

  TFOREST(I,J)=EXP(TS(I,J,0))*T0

 ENDIF
 ENDDO
 ENDDO

RETURN
END SUBROUTINE POCHVA_INIT
! ---------------------------------------------------------------------------------------------------------------------------------
!    Pochva:
!    Numerical scheme for the parametrization of thermo-hydrous processes
!    in soil and at its surface (own soil, vegetation, snow cover)
!
! Scheme is elaborated by O.Drofa, CNR-ISAC, Italy (o.drofa-at-isac.cnr.it)
!
! ---------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE POCHVA_PAR_CHANGE(MASK_SOIL_EXT, TS_EXT, T_SURF_EXT, SNOW, &
 LAIVEG, LAIVEG_MAX, FRACVEG, IND_LEV_TBOTTOM, IND_LEV_QBOTTOM, VALMISS_EXT, &
 FLAG_SNOW_INIT, LEV_SNOW_EXT, TSNOW_EXT, FICE_SNOW_EXT, SNOW_AGE_EXT, SNOW_MELT_AGE_EXT)

! Procedure changes of internal array of fisical parameters in the case of
! changing (external) of LAI and cover faction of vegetation,
! and, in the case of changhing (external) of the status of grid point for the case of thiсk sea ice,
! thick sea ice point is sumulated as a glacier point, and thus points may appear and disappear.

!USE MODULE_POCHVA
!IMPLICIT NONE

INTEGER, DIMENSION(:,:) :: MASK_SOIL_EXT, IND_LEV_TBOTTOM, IND_LEV_QBOTTOM
REAL, DIMENSION(:,:,:) :: TS_EXT
REAL, DIMENSION(:,:) :: T_SURF_EXT, SNOW, LAIVEG, FRACVEG
REAL :: LAIVEG_MAX, VALMISS_EXT

INTEGER, OPTIONAL :: FLAG_SNOW_INIT
REAL, DIMENSION(:,:,:), OPTIONAL :: LEV_SNOW_EXT, TSNOW_EXT, &
 FICE_SNOW_EXT, SNOW_AGE_EXT, SNOW_MELT_AGE_EXT

INTEGER :: I, J, K, NLEV, N, NN, FLAG_SNOW
REAL :: QW_VEGLOW_MAX_OLD, Z01

! Eventual changing of low vegetation LAI and fraction

 DO J=1,NPOINT_Y
 DO I=1,NPOINT_X
 IF (MASK_SOIL(I,J) == 1) THEN

! Leaf Area Index of low vegetation (dimensionless)
   LAI_VEGLOW(I,J)=LAIVEG(I,J)

! Fraction of low vegetation (not forest) (proportion)
   FVEGLOW(I,J)=FRACVEG(I,J)

! Maximum value of water deposited on the leaf surface of the low vegetation (kg/m^2)
   QW_VEGLOW_MAX_OLD=QW_VEGLOW_MAX(J,J)
   IF (LAIVEG_MAX > 0.) THEN
!     QW_VEGLOW_MAX(I,1:J)=LAI_VEGLOW(I,J)/LAIVEG_MAX*10. ! 10 kg/m**2 at LAI of tropical forest Run4
     QW_VEGLOW_MAX(J,J)=LAI_VEGLOW(I,J)/LAIVEG_MAX*5. ! 5 kg/m**2 at LAI of tropical forest
   ELSE
     QW_VEGLOW_MAX(I,J)=0.
   ENDIF 

! Water deposited on the leaf surface of the low vegetation (kg/m^2)
   
   IF (QW_VEGLOW_MAX_OLD > 0..AND.LAI_VEGLOW(I,J) < QW_VEGLOW_MAX_OLD) THEN
     Z01=QW_VEGLOW(I,J)*(1.-QW_VEGLOW_MAX(I,J)/QW_VEGLOW_MAX_OLD)
     QW_VEGLOW(I,J)=MAX(QW_VEGLOW(I,J)-Z01, 0.)
     IF (T_SURF_EXT(I,J) >= T0) THEN
       FLUX_PREC_LIQ(I,J)=FLUX_PREC_LIQ(I,J)+Z01/DTSTEP
     ELSE
       FLUX_PREC_ICE(I,J)=FLUX_PREC_ICE(I,J)+Z01/DTSTEP
     ENDIF
   ENDIF

 ENDIF
 ENDDO
 ENDDO

! Eventual changing of grid point status: sea point with thick ice cover becames! "land" point with glacier,
! "land" point with glacier that is thick sea ice becames "normal" sea point 

 DO J=1,NPOINT_Y
 DO I=1,NPOINT_X

   IF (MASK_SOIL_EXT(I,J) == 1.AND.MASK_SOIL(I,J) == 0) THEN
! Sea pont becames land point with glacier

     MASK_SOIL(I,J)=1
     MASK_SOIL_AXED(I,J)=1 ! ! Flag for application of the Pochva scheme in the axed version
     NLEV_SOIL_HEAT(I,J)=IND_LEV_TBOTTOM(I,J)
     NLEV_SOIL_WATER(I,J)=IND_LEV_QBOTTOM(I,J)
 
! Constant (in time) soil and vegetation physical parameters

     QSMAX(I,J,:)=1.
     QSMIN(I,J,:)=0.
     CG(I,J,:)=CI
     RHOG(I,J,:)=915.
     PSIS_SAT(I,J,:)=2.2
     KAPPA_SAT(I,J,:)=0.
     PAR_B(I,J,:)=0.
     PAR_C(I,J,:)=0.
     QS_REL_VEGWILT(I,J,:)=0.
     QS_REL_VEGREF(I,J,:)=0.
     LAI_VEGLOW(I,J)=0.
     LAI_FOREST(I,J)=0.
     FVEGLOW(I,J)=0.
     FFOREST(I,J)=0.
     FVEGLOWWET(I,J)=0.
     FFORESTWET(I,J)=0.
     ROOT_DEPTH_VEGLOW(I,J)=0.
     ROOT_DEPTH_FOREST(I,J)=0.
     QW_VEGLOW(I,J)=0.
     QW_FOREST(I,J)=0.
     PAR_FB(I,J,:)=2.
     DZSUM_VEGLOW(I,J)=0.
     DZSUM_FOREST(I,J)=0.

! Initialization of main prognostic variables

     NLEV=NLEV_SOIL_HEAT(I,J)
     TSURF(I,J)=MIN(T_SURF_EXT(I,J), T0-0.3)
     TS(I,J,1:NLEV)=MIN(ALOG(TS_EXT(I,J,1:NLEV)/T0), -1.E-3)
     IF (SNOW(I,J) > 0.) THEN
       IF (PRESENT(TSNOW_EXT).AND.INT(TSNOW_EXT(I,J,1))==INT(VALMISS_EXT)) THEN
         DO K=0,NLEV_SNOW
           TS(I,J,K)=ALOG(TSNOW_EXT(I,J,K+1)/T0)
           IF (INT(TSNOW_EXT(I,J,K+2))==INT(VALMISS_EXT)) EXIT
         ENDDO
       ELSE
         TS(I,J,0)=TS(I,J,1)
       ENDIF
     ELSE
       TS(I,J,0)=MIN(T_SURF_EXT(I,J), -1.E-3)
     ENDIF
     QSI(I,J,0:NLEV)=0.
     FRAC_SICE(I,J,0:NLEV)=0.
     DFRAC_SICE_DT(I,J,0:NLEV)=0.
     NLEV=NLEV_SOIL_WATER(I,J)
     QS(I,J,0:NLEV)=0.
     QS_REL(I,J,0:NLEV)=0.

     FSNCOVER(I,J)=MIN( ((SNOW(I,J)/SNOW_FRMAX)**0.5), 1.0)

     LEV_SNOW(I,J,:)=VAL_MISSING
     TSNOW(I,J,:)=VAL_MISSING
     FICE_SNOW(I,J,:)=VAL_MISSING
     SNOW_AGE(I,J,:)=0.
     SNOW_MELT_AGE(I,J,:)=0.
     RHOSNOW(I,J,:)=VAL_MISSING

     IF (PRESENT(FLAG_SNOW_INIT)) THEN
       FLAG_SNOW=FLAG_SNOW_INIT
     ELSE
       FLAG_SNOW=0
     ENDIF
  
     IF (PRESENT(LEV_SNOW_EXT).AND.FLAG_SNOW == 1) THEN
  
       DO K=0,NLEV_SNOW
         IF (INT(LEV_SNOW_EXT(I,J,K+1)) /= INT(VALMISS_EXT)) THEN
           LEV_SNOW(I,J,K)=LEV_SNOW_EXT(I,J,K+1)
         ELSE
           EXIT
         ENDIF
       ENDDO
       N=K-2
       D_LEV_SNOW(I,J)=D_LEV_SNOW_BASE
       IF (N > 1) THEN
         NN=NINT((LEV_SNOW(I,J,2)-LEV_SNOW(I,J,1))/D_LEV_SNOW_BASE)
         D_LEV_SNOW(I,J)=D_LEV_SNOW_BASE*FLOAT(NN)
       ENDIF
  
     ELSE
  
       D_LEV_SNOW(I,J)=D_LEV_SNOW_BASE
       DO WHILE (.TRUE.)
         IF (D_LEV_SNOW(I,J)*FLOAT(NLEV_SNOW) >= SNOW(I,J)) EXIT
         D_LEV_SNOW(I,J)=D_LEV_SNOW(I,J)+D_LEV_SNOW_BASE
       ENDDO
    
       IF (SNOW(I,J) > 0.) THEN
         LEV_SNOW(I,J,0)=0.
         IF (SNOW(I,J) >= D_LEV_SNOW(I,J)+D_LEV_SNOW_MIN) THEN  
           N=INT(SNOW(I,J)/D_LEV_SNOW(I,J))
           IF (SNOW(I,J)-D_LEV_SNOW(I,J)*FLOAT(N) >= D_LEV_SNOW_MIN) N=N+1
           DO K=1,N
             LEV_SNOW(I,J,K)=SNOW(I,J)-FLOAT(N-K)*D_LEV_SNOW(I,J)
           ENDDO
           N=N-1
         ELSE
           N=0
           LEV_SNOW(I,J,1)=SNOW(I,J)
         ENDIF
       ELSE
         N=-1
       ENDIF
  
     ENDIF
  
     IF (N >= 0) THEN

       IF (PRESENT(TSNOW_EXT).AND.FLAG_SNOW == 1) THEN
         TSNOW(I,J,0:N+1)=ALOG(TSNOW_EXT(I,J,1:N+2)/T0)
       ELSE
         IF (N >= 1 ) THEN
           TSNOW(I,J,0)=MIN(ALOG(TSURF(I,J)/T0), -1.E-3)
           TSNOW(I,J,N+1)=MIN(ALOG(TSURF(I,J)/T0), -1.E-3)
           Z01=1./FLOAT(N+1)
           DO K=1,N
             TSNOW(I,J,K)=TSNOW(I,J,0)*(1.-Z01*FLOAT(K))+TSNOW(I,J,N+1)*Z01*FLOAT(K)
           ENDDO
         ELSEIF (N == 0 ) THEN
           TSNOW(I,J,0)=MIN(ALOG(TSURF(I,J)/T0), -1.E-3)
           TSNOW(I,J,1)=MIN(TS(I,J,0), -1.E-3)
         ENDIF
       ENDIF
    
       IF (PRESENT(FICE_SNOW_EXT).AND.FLAG_SNOW == 1) THEN
         FICE_SNOW(I,J,0:N+1)=FICE_SNOW_EXT(I,J,1:N+2)
       ELSE
         FICE_SNOW(I,J,0:N+1)=1.
       ENDIF
    
       IF (PRESENT(SNOW_AGE_EXT).AND.FLAG_SNOW == 1) THEN
         SNOW_AGE(I,J,0:N+1)=SNOW_AGE_EXT(I,J,1:N+2)
       ELSE
         SNOW_AGE(I,J,0:N+1)=10.
       ENDIF
    
       IF (PRESENT(SNOW_MELT_AGE_EXT).AND.FLAG_SNOW == 1) THEN
         SNOW_MELT_AGE(I,J,0:N+1)=SNOW_MELT_AGE_EXT(I,J,1:N+2)
       ELSE
         SNOW_MELT_AGE(I,J,0:N+1)=0.
       ENDIF
    
       RHOSNOW(I,J,0:N+1)=MIN( RHOSNOW_FRESH*2., RHOSNOW_OLD )

     ENDIF
  
     IF (N >= 0) THEN
       IF (LEV_SNOW(I,J,N+1) >= D_LEV_SNOW_MIN) THEN
         FLAG_SNOW_THICK(I,J)=1
         DLEV_SNOW_BOTTOM(I,J)=(LEV_SNOW(I,J,N+1)-LEV_SNOW(I,J,N))*0.5
       ELSE
         FLAG_SNOW_THICK(I,J)=0
         DLEV_SNOW_BOTTOM(I,J)=LEV_SNOW(I,J,1)
       ENDIF
     ELSE
       DLEV_SNOW_BOTTOM(I,J)=0.
     ENDIF

     TFOREST(I,J)=EXP(TS(I,J,0))*T0

   ELSE IF (MASK_SOIL_EXT(I,J) == 0.AND.MASK_SOIL(I,J) == 1) THEN
! Land point with glacier becames sea point

     MASK_SOIL(I,J)=0
     TSURF(I,J)=VAL_MISSING
     QSURF(I,J)=VAL_MISSING
     TS(I,J,:)=VAL_MISSING
     FRAC_SICE(I,J,:)=VAL_MISSING
     QS(I,J,:)=1.
     SNOW(I,J)=0.
     FSNCOVER(I,J)=0.
     LEV_SNOW(I,J,:)=VAL_MISSING
     TSNOW(I,J,:)=VAL_MISSING
     FICE_SNOW(I,J,:)=VAL_MISSING
     SNOW_AGE(I,J,:)=0.
     SNOW_MELT_AGE(I,J,:)=0.
     RHOSNOW(I,J,:)=VAL_MISSING
     FLUX_SOIL_ENTROPY(I,J,:)=0.
     RUN_OFF(I,J,:)=0.
     FLUX_SOIL_WATER_BOTTOM_DIAG(I,J)=0.
     FLUX_ENTROPY_ATM_SPEC(I,J)=0.
     FLUX_ENTROPY_ATM_LAT(I,J)=0.
     FLUX_WV(I,J)=0.

   ENDIF 

 ENDDO
 ENDDO

RETURN
END SUBROUTINE POCHVA_PAR_CHANGE
! ---------------------------------------------------------------------------------------------------------------------------------

END MODULE POCHVA_SCHEME
