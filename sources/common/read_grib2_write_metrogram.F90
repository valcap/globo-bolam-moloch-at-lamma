!                 DATA_GRIB2_READ_ELABOR

! Procedure that reads input data on regular grib2 in grib2 format,
! and coordinates of geographical sporadic points in ascii format;
! interpolates data from regilar grib to the sporadic points,
! and writes result in convention ascii "meteogram" format.  
!-----------------------------------------------------------------------------
MODULE DATA_GRID_INPUT

! ------  INPUT DATA GRID PARAMETERS ------

 INTEGER :: ICENTRE_INP, ISUBCENTRE_INP, IMODEL_INP

! NX_INP, NY_INP : input data grid dimensions
! NLEV_ATM_INP : input data atmospheric levels number
! NLEV_SOIL_INP : input data soil levels number
! DX_INP, DY_INP: input data grid steps
! ALON_INP, ALAT0_INP: coord. in deg of the SW corner of the input data grid
! X0_INP, Y0_INP: coord. in deg of the rotation of the input data grid

  INTEGER :: NX_INP, NY_INP, NLEV_ATM_INP, NLEV_ATM_INP_MAX=70, NLEV_SOIL_INP, NLEV_SOIL_INP_MAX=20, &
             LEVEL_TYPE_INP, FLAG_ROT_INP, FLAG_CUT_PASTE_INP
  REAL :: DX_INP, DY_INP, X0_INP, Y0_INP, ALON0_INP, ALAT0_INP

! Time parameters
! IDATE0: date and time (year, month, day, hours, minutes) of initial instant of forecast
! IDATEC: validity time of forecast
! IPERIOD: time intervals (days, hours, minutes)
! IPERIOD_ACCUM: time interval for accumulated precipitation

 INTEGER :: IDATE0(5), IDATEC(5), IPERIOD_INP(4), IPERIOD(3), IPERIOD_ACCUM(3)

! Atmospheric and Soil level (m) in input data

 REAL, DIMENSION(:), ALLOCATABLE :: LEV_LIST_ATM_INP, LEV_LIST_SOIL_INP

! Input fields and related parameters:

 INTEGER, PARAMETER :: NPAR3D=20, NPAR2D=50, NPAR3D_SOIL=10

 REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: FIELD3D, FIELD3D_SOIL
 REAL, DIMENSION(:,:,:), ALLOCATABLE   :: FIELD2D
 REAL, DIMENSION(:), ALLOCATABLE       :: ALEV, BLEV

 REAL, DIMENSION(:,:), ALLOCATABLE :: ALON_INP, ALAT_INP, LSM_INP, OROG_INP, LAPSE_RATE_INP

 REAL :: VAL_MISSING=-9999.

 CHARACTER(LEN=20), DIMENSION(NPAR2D) :: NAME_FIELD2D
 CHARACTER(LEN=20), DIMENSION(NPAR3D) :: NAME_FIELD3D_ATM
 CHARACTER(LEN=20), DIMENSION(NPAR3D_SOIL) :: NAME_FIELD3D_SOIL

END MODULE DATA_GRID_INPUT
!-----------------------------------------------------------------------------
MODULE DATA_POINT_INPUT

 REAL :: GAUSS_SEMI_RADIUS, GAUSS_PAR

 INTEGER, PARAMETER :: NPAR3D_ATM_REQ0=10, NPAR3D_SOIL_REQ0=10, NPAR2D_REQ0=10,&
                       NLEV_ATM_REQ=10, NLEV_SOIL_REQ=10

 INTEGER :: INTERP_GAUSS=1, NPAR3D_ATM_REQ, NPAR3D_SOIL_REQ, NPAR2D_REQ

 INTEGER, DIMENSION(NPAR2D_REQ0) :: INDEX_FIELD2D
 INTEGER, DIMENSION(NPAR3D_ATM_REQ0) :: INDEX_FIELD3D
 INTEGER, DIMENSION(NPAR3D_SOIL_REQ0) :: INDEX_FIELD3D_SOIL
 REAL, DIMENSION(NLEV_ATM_REQ,NPAR3D_ATM_REQ0) :: LEV_ATM_REQ
 REAL, DIMENSION(NLEV_SOIL_REQ,NPAR3D_SOIL_REQ0) :: LEV_SOIL_REQ

 REAL, DIMENSION(:), ALLOCATABLE :: POINT_LAT, POINT_LON
 INTEGER, DIMENSION(:), ALLOCATABLE :: POINT_ALT
 CHARACTER(LEN=15), DIMENSION(:), ALLOCATABLE :: POINT_NAME
 INTEGER :: NPOINT

END MODULE DATA_POINT_INPUT
!-----------------------------------------------------------------------------
MODULE DATA_POINT_OUTPUT

 REAL, DIMENSION(:,:), ALLOCATABLE :: METEOPAR2D_POINT
 REAL, DIMENSION(:,:,:), ALLOCATABLE :: METEOPAR3D_ATM_POINT, METEOPAR3D_SOIL_POINT
 CHARACTER(LEN=50) :: OUTPUTNAME
 REAL :: LAPSE_RATE_CLIMATE

END MODULE DATA_POINT_OUTPUT
!-----------------------------------------------------------------------------

PROGRAM READ_GRIB2_WRITE_METEOGRAM

USE DATA_GRID_INPUT
USE DATA_POINT_INPUT
USE DATA_POINT_OUTPUT

IMPLICIT NONE

! Strings for input files names

CHARACTER*80 :: FILE_GRID_INP='input.grib2', FILE_POINT_INP='geopoint_meteogram.txt'

INTEGER :: IFILE, IGRIB, I, J, II, JJ, K, NK, NDAY, NDM, IIMON, ACCUM, IFORECAST_H, IACCUM_H,&
 NTOT, NX_BOL, NY_BOL, NFILT, NF, UNIT_IN=11, UNIT_OUT=12, IPOINT, NNN
REAL :: PI=3.141592654, ZDAY, ZDAYCOEFF

INTEGER, DIMENSION(12) ::IMON=(/31,28,31,30,31,30,31,31,30,31,30,31/)

REAL, DIMENSION(100) :: LEV_ATM_REQ0, LEV_SOIL_REQ0

NAMELIST /PROCES_PARAM/ INTERP_GAUSS, GAUSS_SEMI_RADIUS, GAUSS_PAR, &
 INDEX_FIELD2D, INDEX_FIELD3D, INDEX_FIELD3D_SOIL, LEV_ATM_REQ0, LEV_SOIL_REQ0

!--------------------------------------------------------------------

 OPEN (11, FILE='meteogram_grib2.inp', STATUS='OLD')
 READ (11, PROCES_PARAM)
 CLOSE (11)

 LEV_ATM_REQ(1:NLEV_ATM_REQ,1:NPAR3D_ATM_REQ0)=RESHAPE(LEV_ATM_REQ0,(/NLEV_ATM_REQ,NPAR3D_ATM_REQ0/))
 LEV_SOIL_REQ(1:NLEV_ATM_REQ,1:NPAR3D_ATM_REQ0)=RESHAPE(LEV_SOIL_REQ0,(/NLEV_ATM_REQ,NPAR3D_ATM_REQ0/))

! Reading of input data grid parameters only in input files

 CALL READ_GRIB2_DATA(1, FILE_GRID_INP, 1, .TRUE., FLAG_CUT_PASTE_INP, ICENTRE_INP, ISUBCENTRE_INP, IMODEL_INP, &
  NX_INP, NY_INP, NLEV_ATM_INP, NLEV_ATM_INP_MAX, NLEV_SOIL_INP, NLEV_SOIL_INP_MAX, &
  X0_INP, Y0_INP, ALON0_INP, ALAT0_INP, DX_INP, DY_INP, IDATE0, IPERIOD_INP, &
  LEV_LIST_ATM_INP, LEV_LIST_SOIL_INP, LEVEL_TYPE_INP, &
  NPAR3D, NPAR2D, NPAR3D_SOIL, FIELD3D, FIELD2D, FIELD3D_SOIL, ALEV, BLEV, VAL_MISSING)

 FLAG_ROT_INP=0
 IF (ABS(X0_INP) > 1.E-2 .OR. ABS(Y0_INP) > 1.E-2) FLAG_ROT_INP=1

 IF (ALON0_INP > 180.) ALON0_INP=ALON0_INP-360.

 PRINT *,'Grid parameters of input data on the regular grid:'
 PRINT *,'Center ',ICENTRE_INP,'(98 - ECMWF, 7 - NCEP, 80 - Rome RSMC)'
 IF (ICENTRE_INP == 80) PRINT *,'   SubCenter ',ISUBCENTRE_INP,' (102 - CNR-ISAC)'
 IF (ICENTRE_INP == 80.AND.ISUBCENTRE_INP == 102) &
 PRINT *,'Input data generated by model number ',IMODEL_INP,'(1 - Bolam, 2 - Moloch, 3 - Globo)'
 PRINT *,'Grid parameters:'
 PRINT *,'NX=',NX_INP,' NY=',NY_INP
 PRINT *,'ALON0=',ALON0_INP,'ALAT0=',ALAT0_INP,' DX=',DX_INP,'DY=',DY_INP
 PRINT *,'Rotation: ',FLAG_ROT_INP,'  X0=',X0_INP,'Y0=',Y0_INP
 PRINT *

!---------------------------------------------------------------------------

 DO I=1,NPAR2D_REQ0 
   IF (INDEX_FIELD2D(I) == INT(VAL_MISSING)) EXIT
 ENDDO
 NPAR2D_REQ=I-1

! Dew point temperature --->

 IF (NPAR2D_REQ > 0) THEN

   II=0
   DO I=1,NPAR2D_REQ
     IF (INDEX_FIELD2D(I) == 20) II=I ! Dew point temperature at 2 m
   ENDDO
!! Test --->
!   IF (II /= 0.AND.ICENTRE_INP == 80.AND.ISUBCENTRE_INP == 102) THEN
!     DO I=1,NPAR2D_REQ
!       IF (INDEX_FIELD2D(I) == 19) EXIT ! Temperature at 2 m
!     ENDDO
!     IF (I > NPAR2D_REQ) THEN
!       NPAR2D_REQ=NPAR2D_REQ+1
!       INDEX_FIELD2D(NPAR2D_REQ)=19
!     ENDIF
!     DO I=1,NPAR2D_REQ
!       IF (INDEX_FIELD2D(I) == 30) EXIT ! Relative humidity at 2 m
!     ENDDO
!     IF (I > NPAR2D_REQ) THEN
!       NPAR2D_REQ=NPAR2D_REQ+1
!       INDEX_FIELD2D(NPAR2D_REQ)=30
!     ENDIF
!     DO I=1,NPAR2D_REQ
!       IF (INDEX_FIELD2D(I) == 15) EXIT ! Pressure at M.S.L
!     ENDDO
!     IF (I > NPAR2D_REQ) THEN
!       NPAR2D_REQ=NPAR2D_REQ+1
!       INDEX_FIELD2D(NPAR2D_REQ)=15
!     ENDIF
!   ENDIF
!! <---

 ENDIF
! <---

 DO I=1,NPAR3D_ATM_REQ0 
   IF (INDEX_FIELD3D(I) == INT(VAL_MISSING)) EXIT
 ENDDO
 NPAR3D_ATM_REQ=I-1

 DO I=1,NPAR3D_SOIL_REQ0 
   IF (INDEX_FIELD3D_SOIL(I) == INT(VAL_MISSING)) EXIT
 ENDDO
 NPAR3D_SOIL_REQ=I-1

!---------------------------------------------------------------------------

! Read first file input data

! ALLOCATION OF MATRIXES

 ALLOCATE(FIELD3D(NX_INP,NY_INP,NLEV_ATM_INP_MAX,NPAR3D))
 ALLOCATE(FIELD2D(NX_INP,NY_INP,NPAR2D))
 ALLOCATE(FIELD3D_SOIL(NX_INP,NY_INP,NLEV_SOIL_INP_MAX,NPAR3D_SOIL))
 ALLOCATE(LEV_LIST_ATM_INP(NLEV_ATM_INP_MAX))
 ALLOCATE(LEV_LIST_SOIL_INP(NLEV_SOIL_INP_MAX))
 ALLOCATE(ALEV(NLEV_ATM_INP_MAX+1))
 ALLOCATE(BLEV(NLEV_ATM_INP_MAX+1))
 ALLOCATE(ALON_INP(NX_INP,NY_INP))
 ALLOCATE(ALAT_INP(NX_INP,NY_INP))
 ALLOCATE(LSM_INP(NX_INP,NY_INP))
 ALLOCATE(OROG_INP(NX_INP,NY_INP))
 ALLOCATE(LAPSE_RATE_INP(NX_INP,NY_INP))

 LSM_INP(:,:)=VAL_MISSING
 OROG_INP(:,:)=VAL_MISSING
 LAPSE_RATE_INP(:,:)=VAL_MISSING

 CALL READ_GRIB2_DATA(1, FILE_GRID_INP, 1, .FALSE., FLAG_CUT_PASTE_INP, ICENTRE_INP, ISUBCENTRE_INP, IMODEL_INP, &
  NX_INP, NY_INP, NLEV_ATM_INP, NLEV_ATM_INP_MAX, NLEV_SOIL_INP, NLEV_SOIL_INP_MAX, &
  X0_INP, Y0_INP, ALON0_INP, ALAT0_INP, DX_INP, DY_INP, IDATE0, IPERIOD_INP, &
  LEV_LIST_ATM_INP, LEV_LIST_SOIL_INP, LEVEL_TYPE_INP, &
  NPAR3D, NPAR2D, NPAR3D_SOIL, FIELD3D, FIELD2D, FIELD3D_SOIL, ALEV, BLEV, VAL_MISSING)

 PRINT *,' '
 IF (LEVEL_TYPE_INP == 1) THEN
   PRINT *,' Input data incudes atmospheric parameters at isobaric levels'
 ELSEIF (LEVEL_TYPE_INP == 2) THEN
   PRINT *,' Input data incudes atmospheric parameters at hybrid or sigma levels'
 ENDIF

! Forecast time and accumulation period

 IF (IPERIOD_INP(1)==1) THEN        ! Time unit is hour
   IPERIOD(1)=IPERIOD_INP(2)/24                        ! Days
   IPERIOD(2)=IPERIOD_INP(2)-IPERIOD(1)*24             ! Hours
   IPERIOD(3)=0                                        ! Minutes
   IPERIOD_ACCUM(1)=IPERIOD_INP(3)/24                  ! Days
   IPERIOD_ACCUM(2)=IPERIOD_INP(3)-IPERIOD_ACCUM(1)*24 ! Hours
   IPERIOD_ACCUM(3)=0                                  ! Minutes
 ELSEIF (IPERIOD_INP(1)==0) THEN     ! Time unit is minute
   IPERIOD(1)=IPERIOD_INP(2)/24/60                     ! Days
   IPERIOD(2)=(IPERIOD_INP(2)-IPERIOD(1)*24*60)/60     ! Hours
   IPERIOD(3)=IPERIOD_INP(2)-IPERIOD(1)*24*60-IPERIOD(2)*60 ! Minutes
   IPERIOD_ACCUM(1)=IPERIOD_INP(3)/24/60               ! Days
   IPERIOD_ACCUM(2)=(IPERIOD_INP(3)-IPERIOD_ACCUM(1)*24*60)/60 ! Hours
   IPERIOD_ACCUM(3)=IPERIOD_INP(3)-IPERIOD_ACCUM(1)*24*60-IPERIOD_ACCUM(2)*60 ! Minutes
 ENDIF

! Date and time for current forecast instant

 IDATEC(:)=IDATE0(:)
 IDATEC(3)=IDATEC(3)+IPERIOD(1)
 IDATEC(4)=IDATEC(4)+IPERIOD(2)
 IDATEC(5)=IDATEC(5)+IPERIOD(3)
 IF (MOD(IDATEC(1),4)==0) THEN
   IMON(2)=29
 ELSE
   IMON(2)=28
 ENDIF
 NDAY=IPERIOD(1)
 NDM=0
 DO IIMON=1,50
   NDAY=NDAY-IMON(IDATEC(2))
   IF (NDAY<0) THEN
     EXIT
   ELSE
     NDM=NDM+IMON(IDATEC(2))
     IDATEC(2)=IDATEC(2)+1
     IF (IDATEC(2)>12) THEN
       IDATEC(2)=IDATEC(2)-12
       IDATEC(1)=IDATEC(1)+1
       IF (MOD(IDATEC(1),4)==0) THEN
         IMON(2)=29
       ELSE
         IMON(2)=28
       ENDIF
     ENDIF
   ENDIF
 ENDDO
 IDATEC(3)=IDATEC(3)-NDM
 IF (IDATEC(5)>=60) THEN
   IDATEC(5)=IDATEC(5)-60
   IDATEC(4)=IDATEC(4)+1
 ENDIF
 IF (IDATEC(4)>=24) THEN
   IDATEC(4)=IDATEC(4)-24
   IDATEC(3)=IDATEC(3)+1
 ENDIF
 IF (IDATEC(3)>IMON(IDATEC(2))) THEN
   IDATEC(3)=IDATEC(3)-IMON(IDATEC(2))
   IDATEC(2)=IDATEC(2)+1
   IF (IDATEC(2)>12) THEN
     IDATEC(2)=IDATEC(2)-12
     IDATEC(1)=IDATEC(1)+1
   ENDIF
 ENDIF

 PRINT *,'Date and time parameters of input data on the regular grid:'
 PRINT *,'IDATE0         ',IDATE0(1:5)
 PRINT *,'Forecast period (day hour min) ',IPERIOD(1:3)
 PRINT *,'IDATE_CURRENT ',IDATEC(1:5)
 PRINT *,'Accumulation period (day hour min) ',IPERIOD_ACCUM(1:3)
 PRINT *,'  '

 PRINT *,' In input data found 2d surface field with following indxes:'
 DO I=1,NPAR2D
   IF (ALL(INT(FIELD2D(:,:,I)) /= INT(VAL_MISSING))) THEN
     PRINT '(A,I3)','    ',I
     IF (I == 3) LSM_INP(:,:)=FIELD2D(:,:,I)
     IF (I == 1) OROG_INP(:,:)=FIELD2D(:,:,I)/9.8
     IF (I == 9) OROG_INP(:,:)=FIELD2D(:,:,I)
     IF (I == 27) LAPSE_RATE_INP(:,:)=FIELD2D(:,:,I) ! Lapse Rate when there is one (not two) lapse rate field
!     IF (I == 27) THEN
!       LAPSE_RATE_INP(:,:)=FIELD2D(:,:,I) ! Lapse rate 1
!       PRINT*
!       PRINT *,"  Lapse rate climatic"
!       PRINT*
!     ENDIF
!     IF (I == 28) THEN
!       LAPSE_RATE_INP(:,:)=FIELD2D(:,:,I) ! Lapse rate 2
!       PRINT*
!       PRINT *,"  Lapse rate point by point"
!       PRINT*
!     ENDIF
   ENDIF
 ENDDO
 PRINT *,'  '

 PRINT *,' In input data found 3d atmospheric field with following indxes at followng levels:'
 DO I=1,NPAR3D
   DO K=1,NLEV_ATM_INP
     IF (ALL(INT(FIELD3D(:,:,K,I)) /= INT(VAL_MISSING))) THEN
       PRINT '(A,I3,A,E11.3)','    Index',I,' at level',LEV_LIST_ATM_INP(K)
     ENDIF
   ENDDO
 ENDDO
 PRINT *,'  '

 PRINT *,' In input data found 3d "soil" field with following indxes at followng levels:'
 DO I=1,NPAR3D_SOIL
   DO K=1,NLEV_SOIL_INP
     IF (ALL(INT(FIELD3D_SOIL(:,:,K,I)) /= INT(VAL_MISSING))) THEN
       PRINT '(A,I3,A,E11.3)','    Index',I,' at level',LEV_LIST_SOIL_INP(K)
     ENDIF
   ENDDO
 ENDDO
 PRINT *,'  '

! Geographical coordinate of regular grid points

 IF (FLAG_ROT_INP /= 0 ) THEN

   CALL ROT_GRID(X0_INP, Y0_INP, ALON0_INP, ALAT0_INP, DX_INP, DY_INP, ALON_INP, ALAT_INP, NX_INP, NY_INP)

! The following redefinitions might not be always correct!

 ELSE

   DO I=1,NX_INP
     ALON_INP(I,:)=ALON0_INP+FLOAT(I-1)*DX_INP
   ENDDO
   DO J=1,NY_INP
     ALAT_INP(:,J)=ALAT0_INP+FLOAT(J-1)*DY_INP
   ENDDO

 ENDIF

!---------------------------------------------------------------------------
! Climatological lapse rate

 NNN=IDATEC(3)
 DO I=1,IDATEC(2)-1
   NNN=NNN+IMON(I)
 ENDDO
 ZDAY=FLOAT(NNN)+FLOAT(IDATEC(4))/24.
 IF (ALAT_INP(NX_INP/2,NY_INP/2) > 30.) THEN
   ZDAYCOEFF = 0.5*(1.-COS(ZDAY*2.*PI/365.))   ! N.H.
 ELSEIF (ALAT_INP(NX_INP/2,NY_INP/2) < -30.) THEN
   ZDAYCOEFF = 0.5*(1.+COS(ZDAY*2.*PI/365.))   ! S.H.
 ELSE
   ZDAYCOEFF = 0.5                             ! tropics
 ENDIF
 LAPSE_RATE_CLIMATE = -(4.5e-3 + 2.8e-3*ZDAYCOEFF)
!---------------------------------------------------------------------------

 DO J=1,NY_INP
 DO I=1,NX_INP
   IF (ALON_INP(I,J) < -180.) ALON_INP(I,J)=ALON_INP(I,J) + 360.
 ENDDO
 ENDDO

!---------------------------------------------------------------------------

! IF (ICENTRE_INP /= 80.OR.ISUBCENTRE_INP /= 102) THEN
!   PRINT *,'Input data not generated by ISAC models, STOP'
!   STOP
! ENDIF

!---------------------------------------------------------------------------

! Reading of sporadic point coordinate

 OPEN (UNIT_IN, FILE=FILE_POINT_INP, STATUS='old', ERR=10)
 DO IPOINT=1,10000
   READ (UNIT_IN,*,END=101)
 ENDDO
101 REWIND (UNIT_IN)
 NPOINT=IPOINT-1

 ALLOCATE(POINT_NAME(NPOINT))
 ALLOCATE(POINT_LAT(NPOINT))
 ALLOCATE(POINT_LON(NPOINT))
 ALLOCATE(POINT_ALT(NPOINT))
 ALLOCATE(METEOPAR2D_POINT(NPOINT,NPAR2D_REQ+1))
 ALLOCATE(METEOPAR3D_ATM_POINT(NPOINT,NLEV_ATM_REQ,NPAR3D_ATM_REQ))
 ALLOCATE(METEOPAR3D_SOIL_POINT(NPOINT,NLEV_SOIL_REQ,NPAR3D_SOIL_REQ))

 DO IPOINT=1,NPOINT
   READ (UNIT_IN,*) POINT_NAME(IPOINT), POINT_LAT(IPOINT), POINT_LON(IPOINT), POINT_ALT(IPOINT)
 ENDDO

 CLOSE (UNIT_IN)

 CALL METEO_POINT

! Output writing

 CALL PAR_NAME

! Surface parameters

 DO I=1,NPAR2D_REQ

   II=INDEX_FIELD2D(I)

   WRITE (OUTPUTNAME,'(A,I4.4,4I2.2,3A)') "meteogram_inst_",IDATEC(1:5),"_surf_",TRIM(NAME_FIELD2D(II)),".dat"
   OPEN (UNIT_OUT, FILE=OUTPUTNAME, STATUS='unknown')

   WRITE (UNIT_OUT,'(I4.4,3(1X,I2.2),2X,3(1X,I2.2),2X,I4.4,3(1X,I2.2))') IDATE0(1:4),IPERIOD(1:3),IDATEC(1:4)
   DO IPOINT=1,NPOINT
     WRITE (UNIT_OUT,'(A,2(1X,F7.2),1X,I4,20(F11.1))') &
 POINT_NAME(IPOINT), POINT_LAT(IPOINT), POINT_LON(IPOINT), POINT_ALT(IPOINT),&
 METEOPAR2D_POINT(IPOINT,I)
   ENDDO

   FLUSH (UNIT_OUT)
   CLOSE (UNIT_OUT)

   PRINT *,' '
   PRINT *,'Ouput data written in ',OUTPUTNAME

 ENDDO

! Atmospheric parameters

 DO I=1,NPAR3D_ATM_REQ

   DO K=1,NLEV_ATM_REQ
     IF (INT(LEV_ATM_REQ(K,I)) == INT(VAL_MISSING)) EXIT
   ENDDO
   NK=K-1

   IF (NK == 0) CYCLE

   II=INDEX_FIELD3D(I)

   WRITE (OUTPUTNAME,'(A,I4.4,4I2.2,3A)') "meteogram_inst_",IDATEC(1:5),"_atm_",TRIM(NAME_FIELD3D_ATM(II)),".dat"
   OPEN (UNIT_OUT, FILE=OUTPUTNAME, STATUS='unknown')

   WRITE (UNIT_OUT,'(I4.4,3(1X,I2.2),2X,3(1X,I2.2),2X,I4.4,3(1X,I2.2))') IDATE0(1:4),IPERIOD(1:3),IDATEC(1:4)
   DO IPOINT=1,NPOINT
     WRITE (UNIT_OUT,'(A,2(1X,F7.2),1X,I4,20(F11.1))') &
 POINT_NAME(IPOINT), POINT_LAT(IPOINT), POINT_LON(IPOINT), POINT_ALT(IPOINT),&
 METEOPAR3D_ATM_POINT(IPOINT,1:NK,I)
   ENDDO

   FLUSH (UNIT_OUT)
   CLOSE (UNIT_OUT)

   PRINT *,' '
   PRINT *,'Ouput data written in ',OUTPUTNAME

 ENDDO

! Soil parameters

 DO I=1,NPAR3D_SOIL_REQ

   DO K=1,NLEV_SOIL_REQ
     IF (INT(LEV_SOIL_REQ(K,I)) == INT(VAL_MISSING)) EXIT
   ENDDO
   NK=K-1

   IF (NK == 0) CYCLE

   II=INDEX_FIELD3D_SOIL(I)

   WRITE (OUTPUTNAME,'(A,I4.4,4I2.2,3A)') "meteogram_inst_",IDATEC(1:5),"_soil_",TRIM(NAME_FIELD3D_SOIL(II)),".dat"
   OPEN (UNIT_OUT, FILE=OUTPUTNAME, STATUS='unknown')

   WRITE (UNIT_OUT,'(I4.4,3(1X,I2.2),2X,3(1X,I2.2),2X,I4.4,3(1X,I2.2))') IDATE0(1:4),IPERIOD(1:3),IDATEC(1:4)
   DO IPOINT=1,NPOINT
     WRITE (UNIT_OUT,'(A,2(1X,F7.2),1X,I4,20(F11.1))') &
 POINT_NAME(IPOINT), POINT_LAT(IPOINT), POINT_LON(IPOINT), POINT_ALT(IPOINT),&
 METEOPAR3D_SOIL_POINT(IPOINT,1:NK,I)
   ENDDO

   FLUSH (UNIT_OUT)
   CLOSE (UNIT_OUT)

   PRINT *,' '
   PRINT *,'Ouput data written in ',OUTPUTNAME

 ENDDO

 DEALLOCATE(FIELD3D)
 DEALLOCATE(FIELD2D)
 DEALLOCATE(FIELD3D_SOIL)

10 CONTINUE

!---------------------------------------------------------------------------

STOP
END

!------------------------------------------------------------------------------
SUBROUTINE METEO_POINT
! ----------------------------------------------------------------------------
! Procedure to define values of meteorological parameters
! for given geographical point using input data on regular grid
! ----------------------------------------------------------------------------

USE DATA_GRID_INPUT, ONLY : NLON=>NX_INP, NLAT=>NY_INP, NLEV_ATM=>NLEV_ATM_INP, NLEV_SOIL=>NLEV_SOIL_INP, &
                            ALON0=>ALON0_INP, ALAT0=>ALAT0_INP, DLON=>DX_INP, DLAT=>DY_INP,&
                            X0=>X0_INP, Y0=>Y0_INP, ALON=>ALON_INP, ALAT=>ALAT_INP,&
                            FIELD2D, FIELD3D, FIELD3D_SOIL, LSM=>LSM_INP, OROG=>OROG_INP, LAPSE_RATE=>LAPSE_RATE_INP, &
                            LEV_LIST_ATM=>LEV_LIST_ATM_INP, LEV_LIST_SOIL=>LEV_LIST_SOIL_INP, &
                            IDATE0, IDATEC, IPERIOD, ICENTRE=>ICENTRE_INP, ISUBCENTRE=>ISUBCENTRE_INP, VAL_MISSING
USE DATA_POINT_INPUT
USE DATA_POINT_OUTPUT 

IMPLICIT NONE

!                Work :
REAL, DIMENSION(NLON) :: ALON_ROT
REAL, DIMENSION(NLAT) :: ALAT_ROT
REAL, DIMENSION(NPOINT) :: POINT_LAT_ROT, POINT_LON_ROT, POINT_LAT_ROT_WORK, POINT_LON_ROT_WORK, &
 POINT_LAT_WORK, POINT_LON_WORK, OROG_POINT, LAPSE_RATE_POINT
INTEGER, DIMENSION(NPOINT) :: POINT_ALT_WORK
INTEGER, DIMENSION(NLON,NLAT) :: MASK_FIELD
CHARACTER(LEN=15), DIMENSION(NPOINT) :: POINT_NAME_WORK
REAL, DIMENSION(NPOINT) :: U_WORK, V_WORK, SPEED_WORK, DIR_WORK
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: WORK4D
REAL, DIMENSION(:,:,:), ALLOCATABLE :: WORK3D
INTEGER :: NPOINT_VALID, IPOINT, I, J, K, NK, NP, NNN, IPAR, IFORECAST_H, II, KK, KKK, &
 FLAG_ROT, FLAG_LSM=0, FLAG_OROG=0, FLAG_LAPSE_RATE=0, FLAG,&
 IPAR_T2, IPAR_TD2, IPAR_RH2, IPAR_PMSL
REAL :: WX, WY, VALX1, VALX2, ALFA, ZZZ, &
 RD=287.05, RV=461.51, G=9.807, EPS, T2, RH2, PMSL, PS, HOROG, ZLAPSE_RATE, QS, QSW, QSI, ES, ESW, ESI, Q2, TD2
INTEGER, DIMENSION(12) ::IMON=(/31,28,31,30,31,30,31,31,30,31,30,31/)

 FLAG_ROT=1
 IF (ABS(X0) < 1.E-2.AND.ABS(Y0) < 1.E-2) FLAG_ROT=0

! Definition of rotated coordinate for given geographical points

 ALLOCATE(WORK3D(NPOINT,1,4))
 WORK3D(1:NPOINT,1:1,1)=RESHAPE(POINT_LAT,(/NPOINT,1/)) 
 WORK3D(1:NPOINT,1:1,2)=RESHAPE(POINT_LON,(/NPOINT,1/)) 
 WORK3D(1:NPOINT,1:1,3)=RESHAPE(POINT_LAT_ROT,(/NPOINT,1/)) 
 WORK3D(1:NPOINT,1:1,4)=RESHAPE(POINT_LON_ROT,(/NPOINT,1/)) 

 CALL ANTI_ROT_GRID(X0, Y0, WORK3D(:,:,2), WORK3D(:,:,1), WORK3D(:,:,4), WORK3D(:,:,3), NPOINT, 1)

 POINT_LAT_ROT(1:NPOINT)=WORK3D(1:NPOINT,1,3)
 POINT_LON_ROT(1:NPOINT)=WORK3D(1:NPOINT,1,4)

 DEALLOCATE(WORK3D)

! Linear interpolation for 2D field for "simple case"
! (there are not "not valid" points in input 2d field)

 DO I=1,NLON
   ALON_ROT(I)=ALON0+DLON*FLOAT(I-1)
 ENDDO
 DO J=1,NLAT
   ALAT_ROT(J)=ALAT0+DLAT*FLOAT(J-1)
 ENDDO

! Selected points model domain inside

 NPOINT_VALID=0
 DO IPOINT=1,NPOINT
   IF (FLAG_ROT == 0) THEN ! Not rotation grid
     IF (ALON_ROT(1) < 0..AND.ALON_ROT(NLON) > 0..AND.POINT_LON_ROT(IPOINT) >=180.) &
 POINT_LON_ROT(IPOINT)=POINT_LON_ROT(IPOINT)-360.
     IF (ALON(1,1) < 0..AND.ALON(NLON,1) > 0..AND.POINT_LON(IPOINT) >= 180.) &
 POINT_LON(IPOINT)=POINT_LON(IPOINT)-360.
   ENDIF
   IF (POINT_LAT_ROT(IPOINT)>ALAT_ROT(1).AND.POINT_LAT_ROT(IPOINT)<ALAT_ROT(NLAT).AND.&
POINT_LON_ROT(IPOINT)>ALON_ROT(1).AND.POINT_LON_ROT(IPOINT)<ALON_ROT(NLON)) THEN
     NPOINT_VALID=NPOINT_VALID+1
     POINT_LON_ROT_WORK(NPOINT_VALID)=POINT_LON_ROT(IPOINT)
     POINT_LAT_ROT_WORK(NPOINT_VALID)=POINT_LAT_ROT(IPOINT)
     POINT_LON_WORK(NPOINT_VALID)=POINT_LON(IPOINT)
     POINT_LAT_WORK(NPOINT_VALID)=POINT_LAT(IPOINT)
     POINT_ALT_WORK(NPOINT_VALID)=POINT_ALT(IPOINT)
     POINT_NAME_WORK(NPOINT_VALID)=POINT_NAME(IPOINT)
   ENDIF
 ENDDO

 NPOINT=NPOINT_VALID

 IF (NPOINT == 0) THEN
   PRINT *,' '
   PRINT *,' All requested meteo. point locate outside the input grid domain'
   PRINT *,' '
   STOP
 ENDIF

 POINT_LON_ROT(1:NPOINT)=POINT_LON_ROT_WORK(1:NPOINT)
 POINT_LAT_ROT(1:NPOINT)=POINT_LAT_ROT_WORK(1:NPOINT)
 POINT_LON(1:NPOINT)=POINT_LON_WORK(1:NPOINT)
 POINT_LAT(1:NPOINT)=POINT_LAT_WORK(1:NPOINT)
 POINT_ALT(1:NPOINT)=POINT_ALT_WORK(1:NPOINT)
 POINT_NAME(1:NPOINT)=POINT_NAME_WORK(1:NPOINT)

 FLAG_LSM=0
 IF (ALL(INT(LSM(:,:)) /= INT(VAL_MISSING))) FLAG_LSM=1
  
 FLAG_OROG=0
 IF (ALL(INT(OROG(:,:)) /= INT(VAL_MISSING))) FLAG_OROG=1
  
 FLAG_LAPSE_RATE=0
 IF (ALL(INT(LAPSE_RATE(:,:)) /= INT(VAL_MISSING))) FLAG_LAPSE_RATE=1

 PRINT*
 IF (INTERP_GAUSS == 1) THEN
   PRINT *,'  Horizontal interpolation for requested geographical points with Gaussian algorithm,'
   PRINT *,' semi-radius=',GAUSS_SEMI_RADIUS,' km'
 ELSE
   PRINT *,'  Horizontal interpolation for requested geographical points with bilinear algorithm'
 ENDIF

! Meteorological perameters "interpolation" from grid to points 

 MASK_FIELD(:,:)=1

 FLAG_LSM=0
 IF (ALL(INT(LSM(:,:)) /= INT(VAL_MISSING))) FLAG_LSM=1

! Special filelds
    
 NP=1
 NK=1
 ALLOCATE (WORK4D(NLON,NLAT,NK,NP))
 ALLOCATE (WORK3D(NPOINT,NK,NP))

 IF (INTERP_GAUSS == 1) THEN

   MASK_FIELD(:,:)=1
  
   IF (FLAG_OROG == 1) THEN
  
     WORK4D(1:NLON,1:NLAT,1,1) = MAX(OROG(1:NLON,1:NLAT), 0.)
     MASK_FIELD(:,:)=1
  
! For model verification only --->
!   IF (FLAG_LSM == 1) THEN
!     DO J=1,NLAT
!     DO I=1,NLON
!       IF (LSM(I,J) < 0.9) MASK_FIELD(I,J)=0 ! Sea
!     ENDDO
!     ENDDO
!   ENDIF
! <---
  
     CALL INTERP_GAUSS_INPUT_GRID (WORK4D(1:NLON,1:NLAT,1:NK,1:NP), WORK3D(1:NPOINT,1:NK,1:NP), &
   ALON, ALAT, POINT_LON(1:NPOINT), POINT_LAT(1:NPOINT), &
   NLON, NLAT, NPOINT, NK, NP, MASK_FIELD, GAUSS_SEMI_RADIUS, GAUSS_PAR, VAL_MISSING)
  
     OROG_POINT(1:NPOINT) = WORK3D(1:NPOINT,1,1)
  
   ENDIF
  
   IF (FLAG_LAPSE_RATE == 1) THEN
  
     WORK4D(1:NLON,1:NLAT,1,1) = LAPSE_RATE(1:NLON,1:NLAT)
     MASK_FIELD(:,:)=1
  
! For model verification only --->
!     IF (FLAG_LSM == 1) THEN
!       DO J=1,NLAT
!       DO I=1,NLON
!         IF (LSM(I,J) < 0.5.OR.LAPSE_RATE(I,J) > 0.003) MASK_FIELD(I,J)=0 ! Sea or force inversion
!       ENDDO
!       ENDDO
!     ENDIF
! <---
  
     CALL INTERP_GAUSS_INPUT_GRID (WORK4D(1:NLON,1:NLAT,1:NK,1:NP), WORK3D(1:NPOINT,1:NK,1:NP), &
   ALON, ALAT, POINT_LON(1:NPOINT), POINT_LAT(1:NPOINT), &
   NLON, NLAT, NPOINT, NK, NP, MASK_FIELD, GAUSS_SEMI_RADIUS, GAUSS_PAR, VAL_MISSING)
  
     LAPSE_RATE_POINT(1:NPOINT) = WORK3D(1:NPOINT,1,1)
  
   ENDIF

 ELSE ! Bilinear interpolation

   IF (FLAG_OROG == 1) THEN
     WORK4D(1:NLON,1:NLAT,1,1) = MAX(OROG(1:NLON,1:NLAT), 0.)
     CALL INTERP_SPLINE_2D(WORK4D(1:NLON,1:NLAT,1,1), NLON, NLAT, ALON_ROT, ALAT_ROT, &
 POINT_LON_ROT(1:NPOINT), POINT_LAT_ROT(1:NPOINT), NPOINT, OROG_POINT(1:NPOINT), 1.)
   ENDIF

   IF (FLAG_LAPSE_RATE == 1) THEN
     WORK4D(1:NLON,1:NLAT,1,1) = LAPSE_RATE(1:NLON,1:NLAT)
     CALL INTERP_SPLINE_2D(WORK4D(1:NLON,1:NLAT,1,1), NLON, NLAT, ALON_ROT, ALAT_ROT, &
 POINT_LON_ROT(1:NPOINT), POINT_LAT_ROT(1:NPOINT), NPOINT, LAPSE_RATE_POINT(1:NPOINT), 1.)
   ENDIF
! <---

 ENDIF ! Interpolation mode
  
 DEALLOCATE (WORK4D)
 DEALLOCATE (WORK3D)

! 2d surface fields

 IF (NPAR2D_REQ > 0) THEN

! Dew point temperature --->

   IPAR_TD2=0
   DO IPAR=1,NPAR2D_REQ
     IF (INDEX_FIELD2D(IPAR) == 20) IPAR_TD2=IPAR ! Dew point temperature at 2 m
   ENDDO

! <---

   MASK_FIELD(:,:)=1

   NK=1
   NP=NPAR2D_REQ

   ALLOCATE (WORK4D(NLON,NLAT,NK,NP))
   ALLOCATE (WORK3D(NPOINT,NK,NP))

   IPAR_T2=0
   IPAR_TD2=0
!   IPAR_RH2=0
!   IPAR_PMSL=0
   DO IPAR=1,NPAR2D_REQ
     II=INDEX_FIELD2D(IPAR)
     IF (ALL(INT(FIELD2D(:,:,II)) /= INT(VAL_MISSING))) THEN
       WORK4D(1:NLON,1:NLAT,1,IPAR) = FIELD2D(1:NLON,1:NLAT,II)
     ELSE
       PRINT *," "
       PRINT *,"Requested 2d surface field with index ",II," not found in input data ---> stop"
       STOP
     ENDIF
     IF (II == 19) IPAR_T2=IPAR
     IF (II == 20) IPAR_TD2=IPAR
!     IF (II == 30) IPAR_RH2=IPAR
!     IF (II == 15) IPAR_PMSL=IPAR
   ENDDO

   IF (INTERP_GAUSS == 1) THEN ! Gaussian algorithm

     CALL INTERP_GAUSS_INPUT_GRID (WORK4D(1:NLON,1:NLAT,1:NK,1:NP), WORK3D(1:NPOINT,1:NK,1:NP), &
 ALON, ALAT, POINT_LON(1:NPOINT), POINT_LAT(1:NPOINT), &
 NLON, NLAT, NPOINT, NK, NP, MASK_FIELD, GAUSS_SEMI_RADIUS, GAUSS_PAR, VAL_MISSING)

     DO IPAR=1,NPAR2D_REQ
       METEOPAR2D_POINT(1:NPOINT,IPAR) = WORK3D(1:NPOINT,1,IPAR)
     ENDDO

     IF (IPAR_T2 /= 0.AND.FLAG_LSM /= 0) THEN

       MASK_FIELD(:,:)=1
! For model verification only --->
!       DO J=1,NLAT
!       DO I=1,NLON
!         IF (LSM(I,J) < 0.5) MASK_FIELD(I,J)=0 ! Sea
!       ENDDO
!       ENDDO
! <---

       NP=1
       WORK4D(1:NLON,1:NLAT,1,1) = FIELD2D(1:NLON,1:NLAT,INDEX_FIELD2D(IPAR_T2))

       CALL INTERP_GAUSS_INPUT_GRID (WORK4D(1:NLON,1:NLAT,1:NK,1:NP), WORK3D(1:NPOINT,1:NK,1:NP), &
 ALON, ALAT, POINT_LON(1:NPOINT), POINT_LAT(1:NPOINT), &
 NLON, NLAT, NPOINT, NK, NP, MASK_FIELD, GAUSS_SEMI_RADIUS, GAUSS_PAR, VAL_MISSING)

       METEOPAR2D_POINT(1:NPOINT,IPAR_T2) = WORK3D(1:NPOINT,1,1)

     ENDIF

     IF (IPAR_TD2 /= 0.AND.FLAG_LSM /= 0) THEN

       MASK_FIELD(:,:)=1
! For model verification only --->
!       DO J=1,NLAT
!       DO I=1,NLON
!         IF (LSM(I,J) < 0.5) MASK_FIELD(I,J)=0 ! Sea
!       ENDDO
!       ENDDO
! <---

       NP=1
       WORK4D(1:NLON,1:NLAT,1,1) = FIELD2D(1:NLON,1:NLAT,INDEX_FIELD2D(IPAR_TD2))

       CALL INTERP_GAUSS_INPUT_GRID (WORK4D(1:NLON,1:NLAT,1:NK,1:NP), WORK3D(1:NPOINT,1:NK,1:NP), &
 ALON, ALAT, POINT_LON(1:NPOINT), POINT_LAT(1:NPOINT), &
 NLON, NLAT, NPOINT, NK, NP, MASK_FIELD, GAUSS_SEMI_RADIUS, GAUSS_PAR, VAL_MISSING)

       METEOPAR2D_POINT(1:NPOINT,IPAR_TD2) = WORK3D(1:NPOINT,1,1)

     ENDIF

   ELSE ! Bilinear algorithm 

     DO IPAR=1,NPAR2D_REQ

       CALL INTERP_SPLINE_2D(WORK4D(1:NLON,1:NLAT,1,IPAR), NLON, NLAT, ALON_ROT, ALAT_ROT, &
 POINT_LON_ROT(1:NPOINT), POINT_LAT_ROT(1:NPOINT), NPOINT, METEOPAR2D_POINT(1:NPOINT,IPAR), 1.)

     ENDDO 

   ENDIF ! Interpolation mode

   DEALLOCATE (WORK4D)
   DEALLOCATE (WORK3D)

! Wind at 10 m: speed and direction in output

   DO I=1,NPAR2D_REQ
     IF (INDEX_FIELD2D(I) == 17 ) THEN
       IF (INDEX_FIELD2D(I+1) /= 18 ) THEN
         PRINT *,' '
         PRINT *,'Wind at 10 m:, &
 You choose U-component, but not choose V-component, impossible to create output with wind speed and direction'
         PRINT *,' '
         EXIT
       ENDIF

       NK=1
       U_WORK(1:NPOINT)=METEOPAR2D_POINT(1:NPOINT,I)
       V_WORK(1:NPOINT)=METEOPAR2D_POINT(1:NPOINT,I+1)
       CALL WIND_SPEED_DIR (U_WORK, V_WORK, SPEED_WORK, DIR_WORK, &
 NPOINT, X0, Y0, POINT_LON, POINT_LAT, POINT_LON_ROT, POINT_LAT_ROT, FLAG_ROT)
       METEOPAR2D_POINT(1:NPOINT,I)=SPEED_WORK(1:NPOINT)
       METEOPAR2D_POINT(1:NPOINT,I+1)=DIR_WORK(1:NPOINT)
       EXIT
     ENDIF
   ENDDO

! Temperature at 2 m: correction with lapse rate for orography difference between model grid points and observation point

   IF (ICENTRE == 80.AND.ISUBCENTRE == 102) THEN ! NWP Models of CNR-ISAC

     DO IPAR=1,NPAR2D_REQ
       IF (INDEX_FIELD2D(IPAR) == 19.OR.INDEX_FIELD2D(IPAR) == 20 ) THEN
         IF (FLAG_OROG == 0.OR.FLAG_LAPSE_RATE == 0) THEN
           PRINT *,' '
           PRINT *,' It is impossible to correct temperature at 2 m orography &
   for difference between model grid points and observation point, there are not data about model orography or/and lapse rate ',&
   FLAG_OROG, FLAG_LAPSE_RATE
           PRINT *,' '
         ELSE
           DO I=1,NPOINT
!             IF (ABS(OROG_POINT(I)-FLOAT(POINT_ALT(I))) <= 400..AND.INT(LAPSE_RATE_POINT(I)) /= INT(VAL_MISSING)) THEN ! For model verification only
             IF (INT(LAPSE_RATE_POINT(I)) /= INT(VAL_MISSING)) THEN
!write (34,*) 'Lapse rate correction'
!write (34,*) 'Before ',METEOPAR2D_POINT(I,IPAR)-273.15
               IF (ABS(OROG_POINT(I)-FLOAT(POINT_ALT(I))) <= 500.) THEN
                 ZLAPSE_RATE=LAPSE_RATE_POINT(I)
               ELSE
                 ZZZ=-6.66666667E-4*ABS(FLOAT(POINT_ALT(I))-OROG_POINT(I))+1.33333333
                 ZZZ=MAX (ZZZ, 0.)
                 ZLAPSE_RATE=LAPSE_RATE_POINT(I)*ZZZ+LAPSE_RATE_CLIMATE*(1.-ZZZ)
               ENDIF
               METEOPAR2D_POINT(I,IPAR)=METEOPAR2D_POINT(I,IPAR)+(FLOAT(POINT_ALT(I))-OROG_POINT(I))*ZLAPSE_RATE
!write (34,*) 'After ',METEOPAR2D_POINT(I,IPAR)-273.15
!write (34,*) 'Lapse rate ',ZLAPSE_RATE,'Model lapse rate ',LAPSE_RATE_POINT(I),' Clim. lapse rate ',LAPSE_RATE_CLIMATE,&
! ' Model orog. ',OROG_POINT(I),' Station orog. ',FLOAT(POINT_ALT(I)),&
! 'Orog difference ',FLOAT(POINT_ALT(I))-OROG_POINT(I)
             ELSE
               METEOPAR2D_POINT(I,IPAR)=VAL_MISSING
             ENDIF
           ENDDO
           PRINT *,' '
           IF (INDEX_FIELD2D(IPAR) == 19) &
 PRINT *,' Temperature at 2 m corrected with lapse rate for difference between model grid points and observation point'
           IF (INDEX_FIELD2D(IPAR) == 20) &
 PRINT *,' Dew point temperature at 2 m corrected with lapse rate for difference between model grid points and observation point'
           IF (FLAG_LAPSE_RATE /= 0) THEN
             PRINT *," Lapse rate is simulated by the model"
           ELSE
             PRINT *," Lapse rate is climatological"
           ENDIF
           PRINT *,' '
         ENDIF
       ENDIF
     ENDDO
  
!! Test --->
!! Dew point temperature at 2 m: calculation using corrected (by lapse rate)
!! temperature at 2 m, relative humidity at 2 m, and surface pressure
!  
!     IF (IPAR_TD2 /= 0) THEN
!  
!       EPS=RD/RV
!  
!       DO IPOINT=1,NPOINT
!         HOROG=OROG_POINT(IPOINT)
!         T2=METEOPAR2D_POINT(IPOINT,IPAR_T2)
!         RH2=METEOPAR2D_POINT(IPOINT,IPAR_RH2)
!         PMSL=METEOPAR2D_POINT(IPOINT,IPAR_PMSL)
!         IF (INT(HOROG) /= INT(VAL_MISSING).AND.INT(T2) /= INT(VAL_MISSING).AND.&
!  INT(RH2) /= INT(VAL_MISSING).AND.INT(PMSL) /= INT(VAL_MISSING)) THEN
!           PS=PMSL*EXP(-(G*HOROG)/(T2*RD))
!           CALL QSAT_TETENS(T2, PS, EPS, QS, QSW, QSI, ES, ESW, ESI)
!           Q2=QSW*RH2*1.E-2
!           CALL TD_TETENS(T2, PS, Q2, EPS, TD2)
!           METEOPAR2D_POINT(IPOINT,IPAR_TD2)=TD2
!         ELSE
!           METEOPAR2D_POINT(IPOINT,IPAR_TD2)=VAL_MISSING
!         ENDIF
!       ENDDO
!  
!     ENDIF
!! <---

   ELSE ! Other NWP Models

     DO IPAR=1,NPAR2D_REQ
       IF (INDEX_FIELD2D(IPAR) == 19.OR.INDEX_FIELD2D(IPAR) == 20) THEN
         IF (FLAG_OROG == 0) THEN
           PRINT *,' '
           PRINT *,' It is impossible to correct temperature or dew point temperature at 2 m orography &
   for difference between model grid points and observation point, there are not data about model orography ',FLAG_OROG
           PRINT *,' '
         ELSE
           DO I=1,NPOINT
!             IF (ABS(OROG_POINT(I)-FLOAT(POINT_ALT(I))) <= 400.) THEN ! For model verification only
             IF (ABS(OROG_POINT(I)-FLOAT(POINT_ALT(I))) <= 10000.) THEN
               METEOPAR2D_POINT(I,IPAR)=METEOPAR2D_POINT(I,IPAR)+(FLOAT(POINT_ALT(I))-OROG_POINT(I))*LAPSE_RATE_CLIMATE
             ELSE
               METEOPAR2D_POINT(I,IPAR)=VAL_MISSING
             ENDIF
           ENDDO
           PRINT *,' '
           IF (INDEX_FIELD2D(IPAR) == 19) &
 PRINT *,' Temperature at 2 m corrected with lapse rate for difference between model grid points and observation point'
           IF (INDEX_FIELD2D(IPAR) == 20) &
 PRINT *,' Dew point temperature at 2 m corrected with lapse rate for difference between model grid points and observation point'

           PRINT *,' '
         ENDIF
       ENDIF
     ENDDO

   ENDIF

 ENDIF ! 2d surface fields

! 3d atmospheric fields

 IF (NPAR3D_ATM_REQ > 0) THEN

   MASK_FIELD(:,:)=1

   NP=1

   DO IPAR=1,NPAR3D_ATM_REQ

     DO K=1,NLEV_ATM_REQ
       IF (INT(LEV_ATM_REQ(K,IPAR)) == INT(VAL_MISSING)) EXIT
     ENDDO
     NK=K-1

     IF (NK > 0) THEN

       ALLOCATE (WORK4D(NLON,NLAT,NK,NP))
       ALLOCATE (WORK3D(NPOINT,NK,NP))

       II=INDEX_FIELD3D(IPAR)

       DO K=1,NK
         FLAG=0
         DO KKK=1,NLEV_ATM
           IF (LEV_ATM_REQ(K,IPAR) == LEV_LIST_ATM(KKK)) THEN
             FLAG=1
             EXIT
           ENDIF
         ENDDO 
         IF (FLAG == 0) THEN
           PRINT *," "
           PRINT *, "Requested atmospheric level ",LEV_ATM_REQ(K,IPAR)," not found in input data ---> stop"
           STOP
         ENDIF
!         KK=KKK-1
         KK=KKK
         IF (ALL(INT(FIELD3D(:,:,KK,II)) /= INT(VAL_MISSING))) THEN
           WORK4D(1:NLON,1:NLAT,K,1) = FIELD3D(1:NLON,1:NLAT,KK,II)
         ELSE
           PRINT *," "
           PRINT *, &
"Requested 3d atmospheric field with index ",II," at level ",LEV_ATM_REQ(K,IPAR)," not found in input data ---> stop"
           STOP
         ENDIF
       ENDDO

       IF (INTERP_GAUSS == 1) THEN ! Gaussian algorithm

         CALL INTERP_GAUSS_INPUT_GRID (WORK4D(1:NLON,1:NLAT,1:NK,1:NP), WORK3D(1:NPOINT,1:NK,1:NP), &
 ALON, ALAT, POINT_LON(1:NPOINT), POINT_LAT(1:NPOINT), &
 NLON, NLAT, NPOINT, NK, NP, MASK_FIELD, GAUSS_SEMI_RADIUS, GAUSS_PAR, VAL_MISSING)

         METEOPAR3D_ATM_POINT(1:NPOINT,1:NK,IPAR) = WORK3D(1:NPOINT,1:NK,1)

       ELSE ! Bilinear algorithm

         DO K=1,NK
           CALL INTERP_SPLINE_2D(WORK4D(1:NLON,1:NLAT,K,1), NLON, NLAT, ALON_ROT, ALAT_ROT, &
 POINT_LON_ROT(1:NPOINT), POINT_LAT_ROT(1:NPOINT), NPOINT, METEOPAR3D_ATM_POINT(1:NPOINT,K,IPAR), 1.)
         ENDDO

       ENDIF ! Interpolation mode

       DEALLOCATE (WORK4D)
       DEALLOCATE (WORK3D)

     ENDIF ! NK > 0

   ENDDO ! IPAR

! Wind: speed and direction in output

   DO I=1,NPAR3D_ATM_REQ

     DO K=1,NLEV_ATM_REQ
       IF (INT(LEV_ATM_REQ(K,I)) == INT(VAL_MISSING)) EXIT
     ENDDO
     NK=K-1

     IF (NK == 0) EXIT

     IF (INDEX_FIELD3D(I) == 3 ) THEN
       IF (INDEX_FIELD3D(I+1) /= 4 ) THEN
         PRINT *,' '
         PRINT *,'Wind at an atmospheric level:, &
 You choose U-component, but not choose V-component, impossible to create output with wind speed and direction'
         PRINT *,' '
         EXIT
       ELSE
         FLAG=0
         DO K=1,NK
           IF (LEV_ATM_REQ(K,I) /= LEV_ATM_REQ(K,I+1)) THEN
             FLAG=FLAG+1
             PRINT *,' '
             PRINT *,'Wind at an atmospheric level: Choosen level of U-component ',&
 LEV_ATM_REQ(K,I),' not coicide with choosen level of V-component ',LEV_ATM_REQ(K,I+1),&
 ', impossible to create output with wind speed and direction'
         PRINT *,' '
           ENDIF 
         ENDDO
         IF (FLAG /= 0) EXIT
       ENDIF

       DO K=1,NK
         U_WORK(1:NPOINT)=METEOPAR3D_ATM_POINT(1:NPOINT,K,I)
         V_WORK(1:NPOINT)=METEOPAR3D_ATM_POINT(1:NPOINT,K,I+1)
         CALL WIND_SPEED_DIR (U_WORK, V_WORK, SPEED_WORK, DIR_WORK, &
 NPOINT, X0, Y0, POINT_LON, POINT_LAT, POINT_LON_ROT, POINT_LAT_ROT, FLAG_ROT)
         METEOPAR3D_ATM_POINT(1:NPOINT,K,I)=SPEED_WORK(1:NPOINT)
         METEOPAR3D_ATM_POINT(1:NPOINT,K,I+1)=DIR_WORK(1:NPOINT)
       ENDDO

       EXIT
     ENDIF
   ENDDO

 ENDIF ! 3d atmospheric fields

! 3d soil fields

 IF (NPAR3D_SOIL_REQ > 0) THEN

   MASK_FIELD(:,:)=1

   IF (FLAG_LSM == 1) THEN
     DO J=1,NLAT
     DO I=1,NLON
       IF (LSM(I,J) < 0.5) MASK_FIELD(I,J)=0 ! Sea
     ENDDO
     ENDDO
   ENDIF

   NP=1

   DO IPAR=1,NPAR3D_SOIL_REQ

     DO K=1,NLEV_SOIL_REQ
       IF (INT(LEV_SOIL_REQ(K,IPAR)) == INT(VAL_MISSING)) EXIT
     ENDDO
     NK=K-1

     IF (NK > 0) THEN

       ALLOCATE (WORK4D(NLON,NLAT,NK,NP))
       ALLOCATE (WORK3D(NPOINT,NK,NP))

       II=INDEX_FIELD3D_SOIL(IPAR)

       DO K=1,NK
         FLAG=0
         DO KKK=1,NLEV_SOIL
           IF (LEV_SOIL_REQ(K,IPAR) == LEV_LIST_SOIL(KKK)) THEN
             FLAG=1
             EXIT
           ENDIF
         ENDDO 
         IF (FLAG == 0) THEN
           PRINT *," "
           PRINT *, "Requested soil level ",LEV_SOIL_REQ(K,IPAR)," not found in input data ---> stop"
           STOP
         ENDIF
         KK=KKK
         IF (ALL(INT(FIELD3D_SOIL(:,:,KK,II)) /= INT(VAL_MISSING))) THEN
           WORK4D(1:NLON,1:NLAT,K,1) = FIELD3D_SOIL(1:NLON,1:NLAT,KK,II)
         ELSE
           PRINT *," "
           PRINT *, &
"Requested 3d soil field with index ",II," at level ",LEV_SOIL_REQ(K,IPAR)," not found in input data ---> stop"
           STOP
         ENDIF
       ENDDO

       IF (INTERP_GAUSS == 1) THEN ! Gaussian algorithm

         CALL INTERP_GAUSS_INPUT_GRID (WORK4D(1:NLON,1:NLAT,1:NK,1:NP), WORK3D(1:NPOINT,1:NK,1:NP), &
 ALON, ALAT, POINT_LON(1:NPOINT), POINT_LAT(1:NPOINT), &
 NLON, NLAT, NPOINT, NK, NP, MASK_FIELD, GAUSS_SEMI_RADIUS, GAUSS_PAR, VAL_MISSING)

         METEOPAR3D_SOIL_POINT(1:NPOINT,1:NK,IPAR) = WORK3D(1:NPOINT,1:NK,1)

       ELSE ! Bilinear algorithm

         DO K=1,NK
           CALL INTERP_SPLINE_2D(WORK4D(1:NLON,1:NLAT,K,1), NLON, NLAT, ALON_ROT, ALAT_ROT, &
 POINT_LON_ROT(1:NPOINT), POINT_LAT_ROT(1:NPOINT), NPOINT, METEOPAR3D_SOIL_POINT(1:NPOINT,K,IPAR), 1.)
         ENDDO

       ENDIF ! Interpolation mode

       DEALLOCATE (WORK4D)
       DEALLOCATE (WORK3D)

     ENDIF ! NK > 0

   ENDDO ! IPAR

 ENDIF ! 3d soil fields

RETURN
END SUBROUTINE METEO_POINT
!------------------------------------------------------------------------------

SUBROUTINE WIND_SPEED_DIR (U, V, SPEED, DIR, NPOINT, X0, Y0, LON, LAT, LON_ROT, LAT_ROT, FLAG_ROT)
! Procedure defines speed and direction of wind using u- and v-component of wind
! in input, also in rotated coordinate sistem
!------------------------------------------------------------------------------
IMPLICIT NONE

! Input:
INTEGER :: NPOINT ! No. of vertical velels, No. of geographical points
REAL :: X0, Y0 ! Rotation center (degree)
REAL, DIMENSION(NPOINT) :: U, V
REAL, DIMENSION(NPOINT) :: LON, LAT ! Geographical coordinates of points
REAL, DIMENSION(NPOINT) :: LON_ROT, LAT_ROT ! Rotated coordinates of points
INTEGER :: FLAG_ROT ! Flag of coordinate rotation

! Output:
REAL, DIMENSION(NPOINT) :: SPEED, DIR

! Work:
REAL, DIMENSION(NPOINT,1,8) :: WORK3D
INTEGER :: I
REAL :: PI, ZZZ, ALFA

  PI=ABS(ACOS(-1.))

 DO I=1,NPOINT
   SPEED(I)=SQRT(U(I)**2+V(I)**2)
 ENDDO

 IF (FLAG_ROT /= 0) THEN

   WORK3D(1:NPOINT,1:1,1)=RESHAPE(LON,(/NPOINT,1/))
   WORK3D(1:NPOINT,1:1,2)=RESHAPE(LAT,(/NPOINT,1/))
   WORK3D(1:NPOINT,1:1,3)=RESHAPE(LON_ROT,(/NPOINT,1/))
   WORK3D(1:NPOINT,1:1,4)=RESHAPE(LAT_ROT,(/NPOINT,1/))
  
   WORK3D(1:NPOINT,1:1,5)=RESHAPE(U(1:NPOINT),(/NPOINT,1/))
   WORK3D(1:NPOINT,1:1,6)=RESHAPE(V(1:NPOINT),(/NPOINT,1/))
   CALL ANTI_ROT_WIND(X0,Y0,WORK3D(:,:,1),WORK3D(:,:,2),WORK3D(:,:,3),WORK3D(:,:,4),&
 WORK3D(:,:,5),WORK3D(:,:,6),WORK3D(:,:,7),WORK3D(:,:,8),NPOINT,1)
   U(1:NPOINT)=WORK3D(1:NPOINT,1,7)
   V(1:NPOINT)=WORK3D(1:NPOINT,1,8)
  
! Correction of probabile errors of anti-rotation of wind
  
   DO I=1,NPOINT
     ZZZ=SPEED(I)**2-V(I)**2
     IF (ZZZ >= 0. ) THEN
       ZZZ=SQRT(ZZZ)
       IF (U(I) >= 0. ) THEN
         U(I)=ZZZ
       ELSE
         U(I)=-ZZZ
       ENDIF
     ELSE
       U(I)=0.
       IF (V(I) >= 0.) THEN
         V(I)=SPEED(I)
       ELSE
         V(I)=-SPEED(I)
       ENDIF
     ENDIF
   ENDDO

 ENDIF ! FLAG_ROT

 DO I=1,NPOINT
   IF (SPEED(I) /= 0. ) THEN
     ALFA=ASIN(ABS(V(I))/SPEED(I))
     ALFA=ALFA/PI*180.
     IF (V(I)*U(I) > 0.) THEN
       IF (U(I) > 0.) THEN
        ALFA=270.-ALFA
       ELSE
        ALFA=90.-ALFA
       ENDIF
     ELSE
       IF (U(I) > 0.) THEN
         ALFA=270.+ALFA
       ELSE
         ALFA=90.+ALFA
       ENDIF
     ENDIF
   ELSE
     ALFA=0.
   ENDIF
   DIR(I)=ALFA
 ENDDO

RETURN
END SUBROUTINE WIND_SPEED_DIR
!------------------------------------------------------------------------------
SUBROUTINE PAR_NAME
!------------------------------------------------------------------------------

USE DATA_GRID_INPUT, ONLY : NAME_FIELD2D, NAME_FIELD3D_ATM, NAME_FIELD3D_SOIL
IMPLICIT NONE

 NAME_FIELD2D(:)=""

 NAME_FIELD2D(01:02)="phig"
 NAME_FIELD2D(03)="lsm"
 NAME_FIELD2D(04)="ps"
 NAME_FIELD2D(05)="tsurf"
 NAME_FIELD2D(06)="prectot"
 NAME_FIELD2D(07)="precconv"
 NAME_FIELD2D(08)="precsnow"
 NAME_FIELD2D(09)="orog"
 NAME_FIELD2D(14)="snow"
 NAME_FIELD2D(15)="pmsl"
 NAME_FIELD2D(16)="clcover"
 NAME_FIELD2D(17)="windspeed10"
 NAME_FIELD2D(18)="winddir10"
 NAME_FIELD2D(19)="t2"
 NAME_FIELD2D(20)="td2"
 NAME_FIELD2D(21)="soiltype"
 NAME_FIELD2D(22)="fseaice"
 NAME_FIELD2D(23)="columnwater"
 NAME_FIELD2D(24)="columncondens"
 NAME_FIELD2D(25)="qsurf"
 NAME_FIELD2D(26)="thseaice"
 NAME_FIELD2D(27)="lapserate"
 NAME_FIELD2D(28)="lapserate2"
 NAME_FIELD2D(30)="rh2"
 NAME_FIELD2D(31)="cllcover"
 NAME_FIELD2D(32)="clmcover"
 NAME_FIELD2D(33)="clhcover"
 NAME_FIELD2D(34)="flsh"
 NAME_FIELD2D(35)="fllh"
 NAME_FIELD2D(36)="flswrad"
 NAME_FIELD2D(37)="fllwrad"
 NAME_FIELD2D(45)="qgmax"

 NAME_FIELD3D_ATM(:)=""

 NAME_FIELD3D_ATM(01)="phi"
 NAME_FIELD3D_ATM(02)="t"
 NAME_FIELD3D_ATM(03)="windspeed"
 NAME_FIELD3D_ATM(04)="winddir"
 NAME_FIELD3D_ATM(05)="q"
 NAME_FIELD3D_ATM(06)="rh"
 NAME_FIELD3D_ATM(07)="qc"
 NAME_FIELD3D_ATM(08)="qcw"
 NAME_FIELD3D_ATM(09)="qci"
 NAME_FIELD3D_ATM(10)="w"
 NAME_FIELD3D_ATM(11)="p"
 NAME_FIELD3D_ATM(12)="z"

 NAME_FIELD3D_SOIL(:)=""

 NAME_FIELD3D_SOIL(01)="ts"
 NAME_FIELD3D_SOIL(02)="qs"
 NAME_FIELD3D_SOIL(03)="qsw"
 NAME_FIELD3D_SOIL(04)="qsi"
 NAME_FIELD3D_SOIL(05)="tsea"
 NAME_FIELD3D_SOIL(06)="tseaice"

RETURN
END SUBROUTINE PAR_NAME
!------------------------------------------------------------------------------
