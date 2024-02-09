!-----------------------------------------------------------------------------
PROGRAM PREPAR_GRIB2_CODING_LAT_LON_DOMAIN

! Procedure read domain parameter from namelist file domain.inp,
! calculates not rotated (geografical) latitude (°N) and longitude (°E)
! of grid points, and prepares calculated for coding in grib2 format

! Uses module grib2_coding_data in write_grib2_data.F90

USE GRIB2_CODING_DATA

IMPLICIT NONE

INTEGER :: IMODEL, NLON, NLAT, NPOINT_JUMP, NX, NY
REAL :: DX, DY, X0, Y0, X1, Y1
INTEGER, DIMENSION(5) :: IDATE0, IDATEC
INTEGER, DIMENSION(3) :: IPERIOD_INP=0, IPERIOD_ACCUM=0
REAL, DIMENSION(:,:), ALLOCATABLE :: LAT, LON
CHARACTER (LEN=10) :: MODEL_NAME
INTEGER :: IFIELD, IX1, IX2, IY1, IY2

NAMELIST /DOMAIN_PARAM/ IMODEL, NLON, NLAT, NPOINT_JUMP, DX, DY, X0, Y0, X1, Y1, &
 IDATE0

  OPEN (11, FILE='domain.inp', STATUS='OLD')
  READ (11, DOMAIN_PARAM)
  CLOSE (11)
  PRINT *,'Parameters of domain (defined in domain.inp):'
  PRINT DOMAIN_PARAM

  IDATEC(:)=IDATE0(:)

  ALLOCATE(LAT(NLON,NLAT))
  ALLOCATE(LON(NLON,NLAT))

  CALL ROT_GRID (X0, Y0, X1, Y1, DX, DY, LON, LAT, NLON, NLAT)

! Name of ouput file

  IF (IMODEL == 1) MODEL_NAME="bolam"
  IF (IMODEL == 2) MODEL_NAME="moloch"
  IF (IMODEL == 3) MODEL_NAME="globo"
  IF (IMODEL == 12) MODEL_NAME="blended"

  WRITE (OUTPUT_FILE_NAME,'(2A)') TRIM(MODEL_NAME),"_domain_geograph_lat_lon.grib2"

  NFIELD = 2

! Cleaning of possible previous allocations

  DO IFIELD=1,NFIELD0
    IF (ALLOCATED(DATA(IFIELD) % FIELD) ) DEALLOCATE (DATA(IFIELD) % FIELD)
  ENDDO
  DO IFIELD=1,NFIELD
    DATA(IFIELD) % GRIB2_DESCRIPT(:) = IVALMISS
  ENDDO

! Description of data for grib2 coding (INTEGER, DIMENSION(50) :: GRIB2_DESCRIPT=IVALMISS)
! see module_write_grib2_data.F90

! Content of GRIB2_DESCRPT array - see bottom of file

! Model index: 1 - Bolam, 2 - Moloch, 3 - Globo

  DATA(1:NFIELD) % GRIB2_DESCRIPT(1) = IMODEL

! Status and type of data

  DATA(1:NFIELD) % GRIB2_DESCRIPT(6) = 0  ! operational products
  DATA(1:NFIELD) % GRIB2_DESCRIPT(7) = 1  ! forecast
  DATA(1:NFIELD) % GRIB2_DESCRIPT(30) = 0 ! bit-map absent

! Grid parameters

  DATA(1:NFIELD) % GRIB2_DESCRIPT(2) = 1 ! grid template index - horizontal grid
  NX = NLON/NPOINT_JUMP
  NY = NLAT/NPOINT_JUMP

  DATA(1:NFIELD) % NX = NX
  DATA(1:NFIELD) % NY = NY
  DATA(1:NFIELD) % X0 = X0
  DATA(1:NFIELD) % Y0 = Y0
  DATA(1:NFIELD) % DX = DX*FLOAT(NPOINT_JUMP)
  DATA(1:NFIELD) % DY = DY*FLOAT(NPOINT_JUMP)
  DATA(1:NFIELD) % X00 = X1
  DATA(1:NFIELD) % Y00 = Y1

! Initial date and time

  DO IFIELD=1,NFIELD

    DATA(IFIELD) % IDATE0(1:5) = IDATE0(1:5)

! Date and time of current forecast term

    DATA(IFIELD) % IDATEC(1:5) = IDATEC(1:5)

  ENDDO

! Definition of time unit, forecast and period length in defined time unit

  DATA(1:NFIELD) % GRIB2_DESCRIPT(4) = 1 ! Time unit for Forecast time step and Statistical period : 0 minute, 1 hour, 2 day
  DATA(1:NFIELD) % IFORECAST = IPERIOD_INP(1)*24+IPERIOD_INP(2)   ! Forecast time step
  DATA(1:NFIELD) % IPERIOD = IPERIOD_ACCUM(1)*24+IPERIOD_ACCUM(2) ! Statistical period

! Product template index

  DATA(1:NFIELD) % GRIB2_DESCRIPT(3) = 0 ! instant

! Level/layer parameters (for forecast satellite products parameters of satellite platform and sensor)

  DO IFIELD=1,NFIELD
    DATA(IFIELD) % GRIB2_DESCRIPT(10:15) = IVALMISS
  ENDDO

  DO IFIELD=1,NFIELD
    ALLOCATE(DATA(IFIELD) % FIELD(NX,NY))
  ENDDO

! Definition of data fields and data parameters in grib2 terms

  DATA(1:NFIELD) % GRIB2_DESCRIPT(10) = 1 ! Ground or water surface

! Geographical latitude of grid points
! (see WMO table 
! http://www.wmo.int/pages/prog/www/WMOCodes/WMO306_vI2/LatestVERSION/WMO306_vI2_GRIB2_CodeFlag_en.pdf)

  IFIELD = 1
  DATA(IFIELD) % GRIB2_DESCRIPT(20) =   0 ! Discipline: Meteorological products
  DATA(IFIELD) % GRIB2_DESCRIPT(21) = 191 ! Category:  miscellaneous
  DATA(IFIELD) % GRIB2_DESCRIPT(22) =   1 ! Parameter: Geographical latitude (°N)
  IY2=0
  DO IY1=1,NLAT,NPOINT_JUMP
  IY2=IY2+1
  IX2=0
  DO IX1=1,NLON,NPOINT_JUMP
  IX2=IX2+1
    DATA(IFIELD) % FIELD(IX2,IY2)=LAT(IX1,IY1)
  ENDDO
  ENDDO

! Geographical latitude of grid points
! (see WMO table 
! http://www.wmo.int/pages/prog/www/WMOCodes/WMO306_vI2/LatestVERSION/WMO306_vI2_GRIB2_CodeFlag_en.pdf)

  IFIELD = 2
  DATA(IFIELD) % GRIB2_DESCRIPT(20) =   0 ! Discipline: Meteorological products
  DATA(IFIELD) % GRIB2_DESCRIPT(21) = 191 ! Category:  miscellaneous
  DATA(IFIELD) % GRIB2_DESCRIPT(22) =   2 ! Parameter: Geographical longitude (°E)
  IY2=0
  DO IY1=1,NLAT,NPOINT_JUMP
  IY2=IY2+1
  IX2=0
  DO IX1=1,NLON,NPOINT_JUMP
  IX2=IX2+1
    DATA(IFIELD) % FIELD(IX2,IY2)=LON(IX1,IY1)
  ENDDO
  ENDDO

! End of preparation of output data for coding in grib2 format

! Write output data in grib2 format

  CALL WRITE_GRIB2_DATA

STOP
END

!-----------------------------------------------------------------------------

! Content of GRIB2_DESCRPT array:

!  1 - model index (1 Bolam, 2 Moloch, 3 Globo)
!  2 - grid template index (1 horizontal grid, 1000 vertical cross-section)
!  3 - product template index (0 instant, 8 statistical, 32 forecast satellite, 1 instant individual ensemble, 11 statistical individual ensemble)
!  4 - time unit index (0- minute, 1 hour, 2 day, 3 month, 4 year, 13 second)
!  5 - statistical elaboration type (for statistical products only) : 0 average, 1 accumulation, 2 maximum, 3 minimum
!  6 - production status of data (0 Operational products, 2 Research products, 3 Re-analysis products, 7 Sub-seasonal to seasonal prediction S2S)
!  7 - type of data (0 Analysis, 1 Forecast, 2 Analysis and forecast, 3 Control forecast, 4 Perturbed forecast)
!  8 - indicator of unit of time for the increment between successive fields used for statistical elaboration
! 10 - level (layer) type (for forecast satellite products code of satellite platform and sensor)
! 11 - first scaled value of level (layer)
! 12 - scale of first value of level (layer)
! 13 - second scaled value of level (layer)
! 14 - scale of second value of level (layer)
! 15 - level type of second level for a layer
! 20 - product discipline
! 21 - product category
! 22 - product parameter
! 30 - flag of bit-map presence: 0 not present, 1 present (use VALMISS)
! 31 - member number of Perturbed forecast
! 32 - member index of Perturbed forecast
