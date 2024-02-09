!    Pochva:
!    Numerical scheme for the parametrization of thermo-hydrous processes
!    in soil and at its surface (own soil, vegetation, snow cover)
!
! Scheme is elaborated by O.Drofa, CNR-ISAC, Italy (o.drofa-at-isac.cnr.it)
!
! Versione 21.1.1: August 2021, application of surface layer turbulent koefficient
! that is defined for the air layer of 1 m thickness
! (SUBROUTINE QV_BARE_SOIL, SUBROUTINE QV_LEAF) 
!
! ---------------------------------------------------------------------------------------------------------------------------------
MODULE MODULE_POCHVA

 INTEGER, PARAMETER ::&
 NLEV_SNOW=10    ! Snow level number (including surface)

 INTEGER ::&
 NLEV_SOIL_MAX,& ! Maximum soil levels number that can be used in the scheme (including surface)
 NPOINT_X,&      ! Points number along X-asis of the model grid
 NPOINT_Y        ! Points number along Y-asis of the model grid

 INTEGER :: PROC_IND ! Index of precessor used in the parallel calculation only

 REAL, DIMENSION(:), ALLOCATABLE ::&
 LEV_SOIL,&    ! Soil levels depth (in m, reading from soil surface)
 LEV_SOIL_H,&  ! Soil half-levels depth (in m, reading from soil surface)
 D_LEV_SOIL,&  ! Difference between full-levels depth (m)
 D_LEV_SOIL_H  ! Difference between half-levels depth (m)

 REAL, DIMENSION(:,:,:), ALLOCATABLE :: &
 QS,&          ! Specific volumetric soil water content (m^3/m^3)
 QS_REL,&      ! Relative soil water content (proportion)
 PSIS,&        ! Soil hydraulic potential (m)
 KAPPA_H,&     ! Soil hydraulic conductivity at half-levels (m/s)
 TS,&          ! ln(T/T0), where T is Soil temperature (K), and T0=273.15(K) is triple point temperature
 FRAC_SICE,&   ! Fraction of ice (crystal phase) of soil water (proportion)
 DFRAC_SICE_DT,&! Derivative with respect to temperature of fraction of ice (freeze fase) of soil water
 DQS,&         ! Increment of total soil water content during time step (m^3/m^3)
 QS_OLD,&      ! Specific volumetric soil water content (m^3/m^3) at time step start
 QSI,&         ! Specific volumetric content of soil water in ice phase (m^3/m^3)
 RUN_OFF,&     ! Soil water run-off (kg/m^2)
 FLUX_SOIL_ENTROPY,&  ! Flux of soil entropy (J/K/m^2/s)
 FLUX_SOIL_WATER,&    ! Flux of soil water (kg/m^2/s)
 FLUX_SOIL_WATER_VEG ,&    ! Flux of soil water (kg/m^2/s) because of evapotraspiration
 LAMBDAG,&     ! Soil thermal conductivity (J/s/K/m)
 QSMAX,&       ! Maximum specific volumetric soil water content (m^3/m^3)
 QSMIN,&       ! Minimum specific volumetric soil water content (m^3/m^3)
 CG,&          ! Dry soil thermal capacity (J/K/m^3) 
 RHOG,&        ! Dry soil density (kg/m^3) 
 PSIS_SAT,&    ! Idraulic potential of saturated soil (m) 
 KAPPA_SAT,&   ! Idraulic conductivity of saturated soil (m/s)
 PAR_B,&       ! Parameter in soil water transport equation (dimensionless)
 PAR_C,&       ! Parameter in soil water transport equation (dimensionless)
 PAR_FB,&      ! Parameter for determination of soil ice fraction (dimensionless) >=1 (low sloped function) and <=2 (great sloped function function)
 QS_REL_VEGWILT,&  ! Relative soil water content (proportion) at which the vegetation wilt
 QS_REL_VEGREF,&   ! Relative soil water content (proportion) at which the vegetation stop evapotraspiration increase
 QS_REL_VEG,&      ! Relative, rispect to QS_REL_VEGWILT and QS_REL_VEGREF, value of soil water content (used in the evapotraspiration's simulation)
 DZ_DZSUM_VEGLOW,& ! Fraction of a soil layer with respect to summary soil layer affacted by low vegetation evapotraspiration processe 
 DZ_DZSUM_FOREST,& ! Fraction of a soil layer with respect to summary soil layer affacted by forest (high) vegetation evapotraspiration processe 
 TSNOW,&       ! ln(Tsnow/T0), where Tsnow is Snow temperature (K),  and T0=273.15(K) is triple point temperature
 FICE_SNOW,&   ! Fraction (in mass) of ice phase in the snow (proportion)
 LEV_SNOW,&    ! Snow levels depth (in kg/m^2, starting from snow surface)
 RHOSNOW,&     ! Snow density (kg/m^3)
 SNOW_AGE,&    ! Total age of snow (day)
 SNOW_MELT_AGE,&      ! Period when snow was be exposed to melting (day)
 FLUX_SNOW_ENTROPY,&  ! Flux of snow entropy (J/K/m^2/s)
 FLUX_SNOW_WATER      ! Flux of snow water (kg/m^2/s)

 REAL, DIMENSION(:,:), ALLOCATABLE :: &
 FLUX_SOIL_WATER_BOTTOM_DIAG, &    ! Flux of soil water (kg/m^2/s) at soil bottom layer,
! used for diagnostic variable and variable for water conservation control
 FSNCOVER,&           ! Fraction of snow cover (proportion)
 TSURF,&              ! Temperature at the complex terrain surface (K) (representative for a whole grib box)
 QSURF,&              ! Specific humidity at the complex terrain surface (K) (representative for a whole grib box)
 TFOREST,&            ! Temperature at the forest (high vegetation) top (K)
 QFOREST,&            ! Spesific humidity at the forest (high vegetation) top (kg/kg)
 FLUX_ENTROPY_ATM,&        ! Total flux of entropy between the lowset atmospheric layer and the complex terrain surface (J/K/m^2/s)
 FLUX_ENTROPY_ATM_SPEC,&   ! Flux of entropy defined by specipic heat flux between the lowest atmospheric layer and the complex terrain surface (J/K/m^2/s)
 FLUX_ENTROPY_ATM_LAT,&    ! Flux of entropy defined by latent heat flux between the lowest atmospheric layer and the complex terrain surface (J/K/m^2/s)
 FLUX_ENTROPY_ATM_RAD,&    ! Flux of entropy defined by radiaton net flux between the lowest atmospheric layer and the complex terrain surface (J/K/m^2/s)
 FLUX_ENTROPY_ATM_SOIL,&   ! Flux of entropy between the lowest atmospheric lvele and the surface of the soil covered by low vegetation (J/K/m^2/s)
 FLUX_ENTROPY_ATM_SNOW,&   ! Flux of entropy between the lowest atmospheric lvele and the surface of the snow cover (J/K/m^2/s)
 FLUX_ENTROPY_ATM_FOREST,& ! Flux of entropy between the lowest atmospheric lvele and the surface of the forest (high) vegetation (J/K/m^2/s)
 FLUX_ENTROPY_SNOW_SOIL,&  ! Flux of entropy between the lowest snow level and the surface of the soil covered by low vegetation (J/K/m^2/s)
 FLUX_ENTROPY_FOREST_SOIL,&! Flux of entropy between the forest head surface and the surface of the soil covered by low vegetation (J/K/m^2/s)
 FLUX_WV,&                 ! Total flux of water vapour between the lowset atmospheric layer and the complex terrain surface (kg/m^2/s)
 FLUX_WV_SOIL,&            ! Flux of water vapour between the lowset atmospheric layer and surface of bare soil (kg/m^2/s)
 FLUX_WV_SNOW,&            ! Flux of water vapour between the lowset atmospheric layer and surface of the snow cover (kg/m^2/s)
 FLUX_WV_VEGLOW_DRY,&      ! Flux of water vapour between the lowset atmospheric layer and dry leaf surface of low vegetation (kg/m^2/s)
 FLUX_WV_VEGLOW_WET,&      ! Flux of water vapour between the lowset atmospheric layer and wet leaf surface of low vegetation (kg/m^2/s)
 FLUX_WV_FOREST_DRY,&      ! Flux of water vapour between the lowset atmospheric layer and dry leaf surface of forest (high) vegetation (kg/m^2/s)
 FLUX_WV_FOREST_WET,&      ! Flux of water vapour between the lowset atmospheric layer and wet leaf surface of forest (high) vegetation (kg/m^2/s)
 FLUX_PREC_LIQ,&           ! Flux of liquid water precipitation at the terrain surface (kg/m^2/s)
 FLUX_PREC_ICE,&           ! Flux of soil water precipitation at the terrain surface (kg/m^2/s)
 FLUX_W_LIQ_ATM_SOIL,&     ! Flux of liquid water from atmpsphere to the soil surface (kg/m^2/s)
 FLUX_W_LIQ_ATM_SNOW,&     ! Flux of liquid water from atmpsphere to the snow surface (kg/m^2/s)
 FLUX_W_ICE_ATM_SNOW,&     ! Flux of ice water from atmpsphere to the snow surface (kg/m^2/s)
 FLUX_W_LIQ_ATM_VEGLOW,&   ! Flux of liquid water from atmpsphere to the leaf surface of the low vegetation (kg/m^2/s)
 FLUX_W_LIQ_ATM_FOREST,&   ! Flux of liquid water from atmpsphere to the leaf surface of the high (forest) vegetation (kg/m^2/s)
 FLUX_WATER_SNOW_SOIL,&!Water flux between the lowest snow level and the surface of the soil covered by low vegetation (kg/m^2/s)
 FVEGLOW,&            ! Fraction of low vegetation (not forest) (proportion)
 FFOREST,&            ! Fraction of high vegetation (forest) (proportion)
 FVEGLOWWET,&         ! Fraction of wet leaf surface of low vegetation (not forest) (proportion)
 FFORESTWET,&         ! Fraction of wet leaf surface of high vegetation (forest) (proportion)
 LAI_VEGLOW,&         ! Leaf Area Index of low vegetation (dimensionless)
 LAI_FOREST,&         ! Leaf Area Index of forest (high) vegetation (dimensionless)
 ROOT_DEPTH_VEGLOW,&  ! Root depth of low vegetation (m)
 ROOT_DEPTH_FOREST,&  ! Root depth of forest (high) vegetation (m)
 QW_VEGLOW,&          ! Water deposited on the leaf surface of the low vegetation (kg/m^2) 
 QW_FOREST,&          ! Water deposited on the leaf surface of the hight (forest) vegetation (kg/m^2) 
 QW_VEGLOW_MAX,&      ! Maximum value of water deposited on the leaf surface of the low vegetation (kg/m^2) 
 QW_FOREST_MAX,&      ! Maximum value of water deposited on the leaf surface of the hight (forest) vegetation (kg/m^2) 
 DZSUM_VEGLOW,&       ! Summary soil layer affacted by low vegetation evapotraspiration processe 
 DZSUM_FOREST,&       ! Summary soil layer affacted by forest (high) vegetation evapotraspiration processe 
 D_LEV_SNOW,&         ! Step along vertical axis in the snow (kg/m^2)
 DLEV_SNOW_BOTTOM,&   ! Mass (in kg/m^2) of bottom half-layer of snow
 SNOW_DIRT            ! "Dirtibility" of snow surface, used for snow albedo definition, coefficient (0-1)

 REAL, DIMENSION(:,:), ALLOCATABLE :: PREC_ACCUM, FLUX_Q_ACCUM, RUN_OFF_ACCUM, FLUX_WATER_BOTTOM_ACC

 REAL, PARAMETER :: & 
 RHOW=1.E3,&           ! Water density
 RHOI=0.9E3,&          ! Ice density
 P0=1.E5,&             ! Refer thermodynamical total pressure
 E0=611.,&             ! Refer thermodynamical water vapour partial pressure
 T0=273.15,&           ! Refer thermodynamical termerature (triple point)
 TS0=0.,&              ! ALOG(T/T0) T=T0
 T_FICE_CRIT=-30.,&    ! Critical temperature (C deg.) when all soil water is ice phase
 TS_FICE_CRIT=ALOG((T_FICE_CRIT+T0)/T0), &
 TS_REF=0.03,&         ! Soil reference temperature (ln(T/T0)), 0.03 ---> 8.3 deg.C
 RD=287.05,&           ! Dry air gas constant
 RV=461.51,&           ! Water vapour gas constant
 EPS=RD/RV, EP=RV/RD-1., &
 CPD=1004.6,&          ! Specific heat of air for constant pressure
 CPV=1869.46,&         ! Specific heat of water vapour for constant pressure
 CW=4186.8,&           ! Specific heat of water
 CI=CW/2.,&            ! Specific heat of ice
 RHOWCW=RHOW*CW, RHOICI=RHOI*CI,&
 LIW=333560.5,&        ! Specific heat of phase transition from ice to water
 LIV=2834170.5,&       ! Specific heat of phase transition from ice to vapour
 LWV=LIV-LIW,&         ! Specific heat of phase transition from water to vapour
 LIWDT0=LIW/T0, LIVDT0=LIV/T0, &
 CCW1=(CPV-CW)/RV, CCW2=LWV/(T0*RV)-CCW1, CCI1=(CPV-CI)/RV, CCI2=LIV/(T0*RV)-CCI1, &
 LAMBDA_WATER=0.6,&    ! Heat conductivity coefficient of water
 LAMBDA_ICE=2.,&       ! Heat conductivity coefficient of ice
! RHOSNOW_FRESH=100.,&  ! Density of fresh snow (kg/m^3)
! RHOSNOW_OLD=700.,&    ! Density of old show (kg/m^3)
! RHOSNOW_FIRN=800.,&   ! Density of firn (kg/m^3)
 RHOSNOW_FRESH=150.,&  ! Density of fresh snow (kg/m^3)
 RHOSNOW_OLD=800.,&    ! Density of old show (kg/m^3)
 RHOSNOW_FIRN=850.,&   ! Density of firn (kg/m^3)
 DELTA_SNOW=1.E-4,&    ! Snow "depth" that is eqivalent to zero
 D_LEV_SNOW_MIN=5.0,&  ! Minimum step along vertical axis in the snow (kg/m^2)
 D_LEV_SNOW_BASE=10.,& ! Generic value of step along vertical axis in the snow (kg/m^2)
 SNOW_FRMAX=5.,&       ! Water content of snow cover when snow cover fraction became equal 1.
 VAL_MISSING=-9999.

 REAL :: &
 DTSTEP, &    ! Time step
 PAR_FA, LAI_MAX, QV_LEAF_BETA, QV_FOREST_BETA, &
 SUM_D_LEV_SOIL_H     ! Sum of difference between half-levels depth (m)

 INTEGER, DIMENSION(:,:), ALLOCATABLE :: &
 NLEV_SOIL_HEAT,& ! Number of soil levels used in each domain point (including surface) for heat transport
 NLEV_SOIL_WATER,&! Number of soil levels used in each domain point (including surface) for water transport
 MASK_SOIL,&      ! Flag for application of the Pochva scheme
 MASK_SOIL_AXED,& ! Flag for application of the Pochva scheme in the axed version
                  ! with thermal conduction only (for glaciers etc.)
 FLAG_SNOW_THICK  ! Flag for application of heat process scheme in snow
                  ! (if snow is thick, then it applies =1, otherwise
                  ! does not apply =0 )

CONTAINS

! ##################################################################
SUBROUTINE PRES_VAP_SAT(T, ESW, ESI)

! Partial pressure of saturated water vapour over liquid water (ESW) and ice (ESI)
! as the function of the temperature (T) using an analytic formula (based entropy)

IMPLICIT NONE

REAL :: T, ESW, ESI, ZT0T

 ZT0T=T0/T
 ESW=E0*EXP(-CCW1*ALOG(ZT0T)+CCW2*(1.-ZT0T)) ! Partial pressure over water
 ESI=E0*EXP(-CCI1*ALOG(ZT0T)+CCI2*(1.-ZT0T)) ! Partial pressure over ice

END SUBROUTINE PRES_VAP_SAT
! ##################################################################
 SUBROUTINE FRAC_SOIL_ICE

! Fraction (FRAC_SICE) of freez water (ice) in total soil water content
! as the function of TS= ln(T/T0), where TS is soil temperature, and soil physical pameters (via PAR_FB),
! and its first derivate (DFRAC_SICE_DT)

 IMPLICIT NONE
 INTEGER :: I, J, K
 REAL :: Z, ZTANH, ZSECH

 DO J=1,NPOINT_Y
 DO I=1,NPOINT_X
 IF (MASK_SOIL(I,J) ==1 ) THEN
   DO K=0,NLEV_SOIL_HEAT(I,J)
     IF (TS(I,J,K) < 0..AND.TS(I,J,K) > TS_FICE_CRIT) THEN
       Z=TS(I,J,K)*PAR_FA*PAR_FB(I,J,K)
       ZTANH=TANH(Z)
       ZSECH=1./COSH(Z)
       FRAC_SICE(I,J,K)=-ZTANH
       DFRAC_SICE_DT(I,J,K)=-PAR_FA*PAR_FB(I,J,K)*(ZSECH)**2
     ELSEIF (TS(I,J,K) >= 0.) THEN
        FRAC_SICE(I,J,K)=0.
        DFRAC_SICE_DT(I,J,K)=0.
     ELSE
        FRAC_SICE(I,J,K)=1.
        DFRAC_SICE_DT(I,J,K)=0.
      ENDIF
   ENDDO
 ENDIF
 ENDDO
 ENDDO

 END SUBROUTINE FRAC_SOIL_ICE
! ##################################################################
 SUBROUTINE FRAC_SOIL_ICE_POINT(T, FICE, DFICE, I, J, K)

! Fraction (FICE) of freez water (ice) in total soil water content
! as the function of T=ln(T/T0), where T is soil temperature, and soil physical pameters (via PAR_FB),
! its first derivate (DFICE)

 IMPLICIT NONE
 REAL :: T, FICE, DFICE
 INTEGER :: I, J, K
 REAL :: Z, ZTANH, ZSECH

 IF (T < 0..AND.T > TS_FICE_CRIT) THEN
   Z=T*PAR_FA*PAR_FB(I,J,K)
   ZTANH=TANH(Z)
   ZSECH=1./COSH(Z)
   FICE=-ZTANH
   DFICE=-PAR_FA*PAR_FB(I,J,K)*ZSECH**2
 ELSEIF (T >= 0.) THEN
   FICE=0.
   DFICE=0.
 ELSE
   FICE=1.
   DFICE=0.
 ENDIF

 END SUBROUTINE FRAC_SOIL_ICE_POINT
! ##################################################################
 SUBROUTINE INVERT_FRAC_SOIL_ICE_POINT(FICE, T, I, J, K)

! Inverted function of fraction (FICE) of freez water (ice) in total soil water content
! dipending on T=ln(T/T0), where T is soil temperature, and soil physical pameters (via PAR_FB).
! ln(T/T0)=(-tanh^-1(Fice))/a/b

 IMPLICIT NONE
 REAL :: FICE, T
 INTEGER :: I, J, K

 IF (FICE < 1.AND.FICE >= 0.) THEN
   T=-0.5/(PAR_FA*PAR_FB(I,J,K))*ALOG((1.+FICE)/(1.-FICE))
   T=MIN(T, -1.E-3)
 ELSE
   T=TS_FICE_CRIT+1.E-2
 ENDIF

 END SUBROUTINE INVERT_FRAC_SOIL_ICE_POINT
! ##################################################################
 SUBROUTINE FRAC_SNOW_ICE(T, FICE, DFICE)

! Fraction (FICE) of freez water (ice) in total snow water content
! as the function of T=ln(T/T0), where T is snow temperature
! its first derivate (DFICE)

 IMPLICIT NONE
 REAL :: T, FICE, DFICE
 REAL :: Z, ZTANH, ZSECH
 REAL, PARAMETER :: A=10000., A_DEV_2=-5000., B=-0.000366, T2=0.000732 ! B=-ln((0.2+t0)/t0)

 DFICE=0.
 IF (T < T2.AND.T > 0.) THEN ! 0>t>0.2 degree
   Z=(T+B)*A
   ZTANH=TANH(Z)
   ZSECH=1./COSH(Z)
   FICE=(-ZTANH+1.)*0.5
   DFICE=A_DEV_2*ZSECH**2
 ELSEIF (T >= T2) THEN
   FICE=1.E-8 ! at +0.2 degree, it should not be equal to zero 
 ELSE ! t < 0
   FICE=1.
 ENDIF

 END SUBROUTINE FRAC_SNOW_ICE
! ##################################################################
 SUBROUTINE THERM_CONDUCT

! Thermal conductivity (LAMBDAG) as the function of the soil hydraulic potential(PSI)
! and content of freeze water fase in total soil water conent (FRAC_SICE and QS)

 IMPLICIT NONE
 INTEGER :: I, J, K
 REAL :: Z1

 DO J=1,NPOINT_Y
 DO I=1,NPOINT_X
 IF (MASK_SOIL(I,J) == 1) THEN
   IF (MASK_SOIL_AXED(I,J) == 0) THEN
     DO K=0,NLEV_SOIL_HEAT(I,J)
!       IF (PSIS(I,J,K) /= VAL_MISSING ) THEN
!!         LAMBDAG(I,J,K)=MIN( 0.2+300.*EXP(-6.-ALOG10(ABS(PSIS(I,J,K)))), 2.5)+QS(I,J,K)*FRAC_SICE(I,J,K)*LAMBDA_ICE
!         LAMBDAG(I,J,K)=MIN( 0.5+1400.*EXP(-6.-ALOG10(ABS(PSIS(I,J,K)))), 2.5)+QS(I,J,K)*FRAC_SICE(I,J,K)*LAMBDA_ICE
!       ELSE
!         LAMBDAG(I,J,K)=QS(I,J,K)*FRAC_SICE(I,J,K)*LAMBDA_ICE
!       ENDIF
!       LAMBDAG(I,J,K)=MIN(LAMBDAG(I,J,K), 2.5)
       Z1=RHOG(I,J,K)*1.E-3
!       LAMBDAG(I,J,K)=MIN( Z1*0.7*SQRT(QS_REL(I,J,K))+Z1*1.1, 3.2)+QS(I,J,K)*FRAC_SICE(I,J,K)*LAMBDA_ICE
!       LAMBDAG(I,J,K)=MIN(LAMBDAG(I,J,K), 3.2)
!       LAMBDAG(I,J,K)=MIN( Z1*1.0*SQRT(QS_REL(I,J,K))+Z1*0.5, 3.0)+QS(I,J,K)*FRAC_SICE(I,J,K)*LAMBDA_ICE
       LAMBDAG(I,J,K)=MIN( Z1*1.0*SQRT(QS_REL(I,J,K))+Z1*0.3, 3.0)+QS(I,J,K)*FRAC_SICE(I,J,K)*LAMBDA_ICE
       LAMBDAG(I,J,K)=MIN(LAMBDAG(I,J,K), 3.0)
!       LAMBDAG(I,J,K)=MIN( Z1*1.2*QS_REL(I,J,K)+Z1*0.3, 3.0)+QS(I,J,K)*FRAC_SICE(I,J,K)*LAMBDA_ICE
!       LAMBDAG(I,J,K)=MIN(LAMBDAG(I,J,K), 3.0)
     ENDDO
   ELSE
     LAMBDAG(I,J,0:NLEV_SOIL_HEAT(I,J))=PSIS_SAT(I,J,0)
   ENDIF
 ENDIF
 ENDDO
 ENDDO

 END SUBROUTINE THERM_CONDUCT
! ##################################################################
 SUBROUTINE QV_BARE_SOIL(QVSOIL, QSAT, QATM, KTURB, QSOIL, QSOIL_MAX, QSOIL_MIN, B, INDEX)

! Specific humidity over bare soil surface (QVSOIL) as the function
! of saturated (at soil surface temperature) speficic humidity (QSAT),
! specific humidity at the lowset atmospheric level (QATM),
! turbulent exchange coefficient at 1 m over the surface (KTURB),
! soil water content at athe soil top level (QSOIL) and its' maximum (QSOIL_MAX),
! and phisical parametr (B) that characteris the texture at the soil top level,
! INDEX indicates soil (/=1) or not soil (==1) surface (for example, water or glaciers)  

 IMPLICIT NONE

 REAL :: QVSOIL, QSAT, QATM, KTURB, KTURB_LOCAL, QSOIL, QSOIL_MAX, QSOIL_MIN, &
 B, ALFA, KMOL=2.29E-5, QREL, QVREL_SURF=1., FUNC, F1, F2, F3, Z1, ALFA_MIN
 INTEGER :: INDEX

! KTURB_LOCAL=KTURB/3. ! Normalised coefficient
! KTURB_LOCAL=MIN(KTURB_LOCAL, 1.)

 IF (INDEX /= 1) THEN ! land surface

   QREL=(QSOIL-QSOIL_MIN)/(QSOIL_MAX-QSOIL_MIN)

   IF (QREL <= 0.4) THEN ! 1)
!   IF (QREL <= 0.3) THEN ! 2)
     QVREL_SURF=QREL+0.4  ! 1)
!     QVREL_SURF=QREL*2.+0.2 ! 2)
   ELSE
     QVREL_SURF=QREL*0.33333333+0.66666667 ! 1)
!     QVREL_SURF=QREL*0.286+0.714  ! 2)
   ENDIF
   QVREL_SURF=MIN(QVREL_SURF, 1.)

! ALFA_MIN parameter is introduced for the cases of great KTURB,
! when ALFA tends to 0.,
! this is hyperbolic tangent of QREL: th(x)=exp(2*x-1)/exp(2*x+1),
! here th is modified to obtain func=[0,1] with x=[0,1]:
! func={[exp(2*b*(x+a)-1)/exp(2*b*(x+a)+1)]+1}*0.5, a=-0.5 always, b=8 and
! determs rate of function tendency to 0 and 1, if b is greater, then
! this tendency is rapider;
! 0.15 parameter is the ALFA_MIN bigest value when qrel=1., it can be changed.

   Z1=2.*8.*(QREL-0.5)
   ALFA_MIN=(EXP(Z1-1.)/EXP(Z1+1.)+1.)*0.5*0.15

! F1=C+D*(1.-QREL)**F3:
! C - parameter of KTURB influence, if C is greater, then ALFA approachs to zero quickly with KTURB rising,
! D - parameter of QREL influence, if D is greater, then ALFA approachs to zero quickly with QREL decreasing
! F3 - parameter of soil texture influence, if B is greater, then ALFA decreasing if slowly with QREL decreasing 

   F3=0.2+0.05*B

! For KTURB at 1 m
! F1=a+b*(1-qrel)**f3:
! F1 determs the "gaussian radius", if F1 is greater, then radius is smaller 
! a influences at great soil humidity only (greater a --> drier)
! b influences at all soil humidity spectrum (greater b ---> drier)
!   F1=7.*(1.2+3.0*(1.-QREL)**F3) ! 20210819
   F1=7.*(2.0+3.0*(1.-QREL)**F3) ! 20211210, very, very dry with KTURB>1.
!   F1=6.*(2.0+3.0*(1.-QREL)**F3) ! 20220512

   F2=1.0-0.8*(1.-QREL)**F3 ! dry at low qrel

!!!   ALFA=ALFA_MIN+(1.-ALFA_MIN)*2.*F2/(EXP(KTURB*F1)+EXP(-KTURB*F1)) ! for KTURB at 1 m
   ALFA=MAX( 2.*F2/(EXP(KTURB*F1)+EXP(-KTURB*F1)), ALFA_MIN) ! for KTURB at 1 m

 ELSE ! water or ice surface

   ALFA=MAX( 2./(EXP(KTURB*2.)+EXP(-KTURB*2.)), 0.5) ! for KTURB at 1 m

 ENDIF

! QVSOIL=ALFA*QSAT+(1.-ALFA)*QATM
 QVSOIL=ALFA*QSAT*QVREL_SURF+(1.-ALFA)*QATM

 END SUBROUTINE QV_BARE_SOIL
! ##################################################################
 SUBROUTINE QV_LEAF(QVLEAF, ALFA_LEAF, QSAT, QATM, KTURB, RAD_VIS_FLUX, LAI, LAI_MAX,&
 ZQS_REL_VEG, DZ_DZSUM_VEG, ZQV_BETA, NLEV)

! Specific humidity over dry traspirating plant leaf (QVLEAF) as the function
! of saturated (at leaf surface temperature) spesific humidity (QSAT),
! spesific humidity at the lowset atmospheric level (QATM),
! turbulent exchange coefficient at 1m over the surface (KTURB),
! flux of visible radiation at the terrain surface (RAD_VIS_FLUX),
! LeafAreaIndex (LAI), its maximum value permited for used vegetation dataset (LAI_MAX)
! and soil water availability (QS_REL_VEG, DZ_DZSUM_VEG)

! ALFA_LEAF is the factor used for QVLEAF definition (see bellow)
! ZQS_REL_VEG(K)=(QS_REL(K)-QS_REL_VEGWITL(K))/(QS_REL_VEGREF(K)-QS_REL_VEGWITL(K))
! DZ_DZSUM_VEG=DZ(K)/(DZ(1)+DZ(2)+...)

 IMPLICIT NONE

 REAL :: QVLEAF, ALFA, ALFA_LEAF, QSAT, QATM, KTURB, KTURB_LOCAL, KMOL=2.29E-5, RAD_VIS_FLUX, LAI, LAI_MAX,&
 FUNC_LAI, FUNC_VEG, ZQV_BETA, Z1, F1, F2
 REAL, DIMENSION(0:NLEV-1) :: ZQS_REL_VEG, DZ_DZSUM_VEG
 INTEGER :: NLEV, K

! KTURB_LOCAL=KTURB/3. ! Normalised coefficient
! KTURB_LOCAL=MIN(KTURB_LOCAL, 1.)

 ZQV_BETA=0.
 DO K=0,NLEV-1
    ZQV_BETA=ZQV_BETA+ZQS_REL_VEG(K)*DZ_DZSUM_VEG(K)
 ENDDO

! For KTURB at 1 m over the surface
 F1=30.*(-1.9*LAI/LAI_MAX+2.0) ! very dry with KTURB>1.
! F1=20.*(-1.9*LAI/LAI_MAX+2.0) ! 20220512

 Z1=MIN(RAD_VIS_FLUX**0.3/6.815, 1.) ! 600**0.3=6.815
 F2=MIN( Z1, ((LAI/LAI_MAX)**0.2) )

 ALFA_LEAF=2.*F2/(EXP(KTURB*F1)+EXP(-KTURB*F1)) ! for KTURB at 1 m
 ALFA=ALFA_LEAF*ZQV_BETA

 QVLEAF=ALFA*QSAT+(1.-ALFA)*QATM

 END SUBROUTINE QV_LEAF
! ##################################################################
 SUBROUTINE SNOW_DENSITY(I, J, NLEV)

! Snow density RHOSNOW (kg/m^3) as the function of total snow age SNOW_AGE (day) and
! period (day) when snow was be exposed to melting SNOW_MELT_AGE

 IMPLICIT NONE
 INTEGER :: NLEV, I, J, K
 REAL :: ZRHO, ZRHO_AGE, ZRHO_FIRN, ZRHO_NEW, Z01, Z02, Z03, Z04

 Z01=365.**0.3
 Z02=395.**0.3
 Z03=0.1*DTSTEP/86400. ! Snow density may change for 10% during 1 day

 DO K=0,NLEV
   ZRHO=RHOSNOW(I,J,K)
   ZRHO_AGE=MAX( MIN((SNOW_AGE(I,J,K)**0.3/Z01), 1.)*RHOSNOW_OLD, RHOSNOW_FRESH)
   ZRHO_FIRN=MAX( MIN(((SNOW_MELT_AGE(I,J,K)+30.)**0.3/Z02), 1.)*RHOSNOW_FIRN, ZRHO_AGE)
   Z04=MIN( 1.+0.5*(LEV_SNOW(I,J,K)/100.)**0.5, 1.5) 
   ZRHO_NEW=MIN(ZRHO_FIRN*Z04, RHOSNOW_FIRN)
   IF (ZRHO_NEW >= ZRHO) THEN
     RHOSNOW(I,J,K)=MIN( ZRHO_NEW, ZRHO*(1.+Z03))
   ELSE
     RHOSNOW(I,J,K)=MAX( ZRHO_NEW, ZRHO*(1.-Z03))
   ENDIF
 ENDDO

 END SUBROUTINE SNOW_DENSITY
! ##################################################################
 SUBROUTINE SNOW_ALBEDO_DEF_NO_SNOW_DIRT(SNALBED, I, J)

! Snow albedo (SNALBED, proportion) as the function of total snow age SNOW_AGE (day),
! period (day) when snow was be exposed to melting SNOW_MELT_AGE, and
! fraction of ice phase (FICE_SNOW) in upper snow layer

 IMPLICIT NONE
 REAL :: SNALBED
 INTEGER :: I, J

 IF (SNOW_AGE(I,J,0) > 0.) THEN
!!!   SNALBED=0.8/(SNOW_AGE(I,J,0)**0.15) ! Run4
   SNALBED=0.75/(SNOW_AGE(I,J,0)**0.2)
   SNALBED=MIN( MAX(SNALBED, 0.3), 0.75)
   IF (SNOW_MELT_AGE(I,J,0) > 0.) THEN
     SNALBED=SNALBED/(SNOW_MELT_AGE(I,J,0)**0.5)
     SNALBED=MIN( MAX(SNALBED, 0.3), 0.5)
   ENDIF
   IF (FICE_SNOW(I,J,0) < 0.99.AND.FICE_SNOW(I,J,0) > 0.) THEN
     SNALBED=MIN(SNALBED, 0.4)
   ENDIF
 ELSE
   SNALBED=0.7
 ENDIF

 END SUBROUTINE SNOW_ALBEDO_DEF_NO_SNOW_DIRT
! ##################################################################
 SUBROUTINE SNOW_ALBEDO_DEF(SNALBED, I, J)

! Snow albedo (SNALBED, proportion) as the function of total snow age SNOW_AGE (day),
! period (day) when snow was be exposed to melting SNOW_MELT_AGE,
! fraction of ice phase (FICE_SNOW) in upper snow layer, and
! SNOW_DIRT that is the coefficient (proportion 0-1) of dirt growth
! at snow surface due to vegetation waste, aerosol deposition, etc.

 IMPLICIT NONE
 REAL, PARAMETER :: ALB_SNOW_FRESH=0.95, ALB_SNOW_DIRT=0.35, ALB_FIRN=0.55
 REAL :: SNALBED, EXPCOEFF
 INTEGER :: I, J

 EXPCOEFF=0.17*SNOW_DIRT(I,J)+0.03

 IF (SNOW_AGE(I,J,0) > 0.) THEN
!!!   SNALBED=0.8/(SNOW_AGE(I,J,0)**0.15) ! Run4
! Before 08.06.2020 --->
!!!   SNALBED=0.75/(SNOW_AGE(I,J,0)**0.2)
!!!   SNALBED=MIN( MAX(SNALBED, 0.3), 0.75)
!!!   IF (SNOW_MELT_AGE(I,J,0) > 0.) THEN
!!!     SNALBED=SNALBED/(SNOW_MELT_AGE(I,J,0)**0.5)
!!!     SNALBED=MIN( MAX(SNALBED, 0.3), 0.5)
!!!   ENDIF
! <---
   SNALBED=ALB_SNOW_FRESH/((SNOW_AGE(I,J,0)+1.)**EXPCOEFF)
   SNALBED=MAX(SNALBED, ALB_SNOW_DIRT)
   IF (SNOW_MELT_AGE(I,J,0) > 0.) THEN
     SNALBED=SNALBED/((SNOW_MELT_AGE(I,J,0)+1.)**EXPCOEFF)
     SNALBED=MAX(SNALBED, ALB_SNOW_DIRT)
   ENDIF
   IF (FICE_SNOW(I,J,0) < 0.99.AND.FICE_SNOW(I,J,0) > 0.) THEN
     SNALBED=MIN(SNALBED, ALB_FIRN)
   ENDIF
 ELSE
   SNALBED=0.7
 ENDIF

 END SUBROUTINE SNOW_ALBEDO_DEF
! ##################################################################
 SUBROUTINE INTERP_CONSERV_1D(N_INP, X_INP, F_INP, N_OUT, X_OUT, F_OUT)

! 1D-Interpolation function F along axis X with conservation
! integral of F.

! Attention: extremes of X-point in input and in output must be the same!

! N_INP is input point number; 
! N_OUT is output point number; 
! X_INP is array with input values of point along axis X;
! X_OUT is array with output values of point along axis X;
! F_INP id array with input values of function F;
! F_OUT id array with output values of function F.

 IMPLICIT NONE
 INTEGER :: N_INP, N_OUT
 REAL, DIMENSION(N_INP) :: X_INP, F_INP
 REAL, DIMENSION(N_OUT) :: X_OUT, F_OUT
 REAL, DIMENSION(N_INP+1) :: X_INP_STAR
 REAL, DIMENSION(N_OUT+1) :: X_OUT_STAR
 REAL, DIMENSION(N_INP,N_OUT) :: WEIGHT
 INTEGER :: I, J, LAST_INP

 X_INP_STAR(1) = X_INP(1)
 DO I=2,N_INP
  X_INP_STAR(I) = (X_INP(I)+X_INP(I-1))*0.5
 ENDDO
 X_INP_STAR(N_INP+1) = X_INP(N_INP)

 X_OUT_STAR(1) = X_OUT(1)
 DO I=2,N_OUT
  X_OUT_STAR(I) = (X_OUT(I)+X_OUT(I-1))*0.5
 ENDDO
 X_OUT_STAR(N_OUT+1) = X_OUT(N_OUT)

! Full algorithm without an optimization

! DO J=1,N_OUT
!   DO I=1,N_INP
!     WEIGHT(I,J) = MAX( MIN( X_OUT_STAR(J+1), X_INP_STAR(I+1) )-MAX( X_OUT_STAR(J), X_INP_STAR(I) ), 0. ) / &
! (X_OUT_STAR(J+1) - X_OUT_STAR(J))
!   ENDDO
! ENDDO
!
! F_OUT(:)=0.
! DO J=1,N_OUT
!   DO I=1,N_INP
!     F_OUT(J)=F_OUT(J)+WEIGHT(I,J)*F_INP(I)
!   ENDDO
! ENDDO

! Optimized algorithm

 F_OUT(:)=0.
 LAST_INP=1
 DO J=1,N_OUT
  DO I=LAST_INP,N_INP
    IF (X_INP_STAR(I) > X_OUT_STAR(J+1)) THEN
      LAST_INP=MAX( I-1, 1) ! ricordo dove ripartire alla prossima iterazione
      EXIT ! ho esaurito i livelli INP che influenzano il livello J OUT
    ENDIF
    F_OUT(J)=F_OUT(J)+F_INP(I)*(MAX ( &
 MIN(X_OUT_STAR(J+1), X_INP_STAR(I+1) ) - MAX( X_OUT_STAR(J), X_INP_STAR(I)), 0.))/(X_OUT_STAR(J+1) - X_OUT_STAR(J))
   ENDDO
 ENDDO

 END SUBROUTINE INTERP_CONSERV_1D
! ##################################################################
 SUBROUTINE INTERP_1D(N_INP, X_INP, F_INP, N_OUT, X_OUT, F_OUT)

! 1D-Interpolation function F along axis X by linear algorithm

! N_INP is input point number; 
! N_OUT is output point number; 
! X_INP is array with input values of point along axis X;
! X_OUT is array with output values of point along axis X;
! F_INP id array with input values of function F;
! F_OUT id array with output values of function F.

 IMPLICIT NONE
 INTEGER :: N_INP, N_OUT
 REAL, DIMENSION(N_INP) :: X_INP, F_INP
 REAL, DIMENSION(N_OUT) :: X_OUT, F_OUT
 REAL :: WEIGHT
 INTEGER :: J, I, I1, I2

 DO J=1,N_OUT

  IF (X_OUT(J) >= X_INP(N_INP)) THEN

    F_OUT(J)=F_INP(N_INP)

  ELSEIF (X_OUT(J) <= X_INP(1)) THEN

    F_OUT(J)=F_INP(1)

  ELSE

    DO I=1,N_INP
      IF (X_OUT(J) >= X_INP(I)) EXIT
    ENDDO
    I1=I
    I2=I1+1
    WEIGHT=(X_OUT(J)-X_INP(I1))/(X_INP(I2)-X_INP(I1))
    F_OUT(J)=F_INP(I1)*(1.-WEIGHT)+F_INP(I2)*WEIGHT

  ENDIF

 ENDDO

 END SUBROUTINE INTERP_1D
! ##################################################################

END MODULE MODULE_POCHVA
! --------------------------------------------------------------------------------------------------------------------------------
