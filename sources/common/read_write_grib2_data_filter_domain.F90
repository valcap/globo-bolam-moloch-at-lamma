PROGRAM DATA_GRIB2_READ_WRITE

! Last update 05/10/2011

! AUTHOR: Oxana Drofa, ISAC-CNR, Bologna, Italy (o.drofa@isac.cnr.it)

! Procedure for reading and writing of meteorol. parameter fields in GRIB2 format
! using the ECMWF ECCODES library

!------------------------------------------------------------------------------

USE eccodes

IMPLICIT NONE

INTEGER :: NFILE_INPUT=1, FLAG_CUT_PASTE=0
CHARACTER, DIMENSION(5) :: FILENAME_INPUT*30='input_00.grib2'
INTEGER :: ICENTRE_CODE, ICENTRE_SUBCODE, IMODEL_CODE
REAL :: VAL_MISSING=-9999.

! For reading GRIB2 file

TYPE GRIB_FIELD
  INTEGER :: PARAM_DISCIPL, PARAM_CATEG, PARAM_IND, LEV_TYPE, LEV_TOP(2), LEV_BOT(2),&
  NXGRIB, NYGRIB, IDATE0(5),IDATEC(5), IPERIOD(5), N_VERT_COORD
  REAL :: X00GRIB, Y00GRIB, X0GRIB, Y0GRIB, XNGRIB, YNGRIB, DXGRIB, DYGRIB
  REAL, DIMENSION(:,:), ALLOCATABLE :: FIELD
  REAL, DIMENSION(:), ALLOCATABLE :: VERT_COORD
END TYPE GRIB_FIELD

TYPE (GRIB_FIELD), DIMENSION(1000) :: FIELD_2D

INTEGER :: IFILE, IGRIB, IIF, IERR
INTEGER :: JCOUNT, NCOUNT, NSIZE,&
             ISCAN, JSCAN, JSCANCONS, IINI, IFIN, DI, JINI, JFIN, DJ, ORDER(2)
REAL*4, DIMENSION(:), ALLOCATABLE :: VALUE
REAL :: ZZZ

! For writting grib2

INTEGER :: IGRIBSAMPLE1, IGRIBSAMPLE, IGRIBOUT, IGRIBCLONE, NMESSAGES
CHARACTER(LEN=30) :: GRIBOUTNAME='output.grib2'
CHARACTER(LEN=200) :: GRIBSAMPLENAME
INTEGER :: IGRIB_EDITION=2, IBITNUM=16, IDATEINI, ITIMEINI, IMESSAGE

NAMELIST /GRIB_SAMPLE/ GRIBSAMPLENAME

! For data processing

INTEGER :: III, I, J, NX_OUT, NY_OUT
REAL :: PI, X0_OUT, Y0_OUT, XN_OUT, YN_OUT
REAL, DIMENSION(:,:), ALLOCATABLE :: FIELD_OUT
INTEGER :: SLICE_XINI, SLICE_XFIN, SLICE_YINI, SLICE_YFIN

NAMELIST /SLICE_DOMAIN/ SLICE_XINI, SLICE_XFIN, SLICE_YINI, SLICE_YFIN

!------------------------------------------------------------------------------
! Level types, GRIB2 codes:

! Ground or water surface                  001
! Isobaric surface (Pa)                    100
! Mean sea level                           101
! Specified height level above m.sea.l.(m) 102
! Specified height level above ground (m)  103
! Sigma level (sigma value)                104
! Hybrid level                             105
! Depth below land surface (m)             106
! Depth below sea level (m)                160
!------------------------------------------------------------------------------

 PI = ABS(ACOS(-1.))

 OPEN (10, FILE='grib_sample.inp', STATUS='OLD')
 READ (10, GRIB_SAMPLE)
 CLOSE (10)

 OPEN (10, FILE='slice_domain.inp', STATUS='OLD')
 READ (10, SLICE_DOMAIN)
 CLOSE (10)

! Reading of GRIB2 files

  JCOUNT = 0

! Loop for GRIB2 files

 FILELOOP: DO IIF=1,NFILE_INPUT

  WRITE(FILENAME_INPUT(IIF)(7:8),'(I2.2)') IIF

! Open file

  CALL GRIB_OPEN_FILE(IFILE, TRIM(FILENAME_INPUT(IIF)),"r")

! Turn on support for multi fields messages */

  CALL GRIB_MULTI_SUPPORT_ON()

! Turn off support for multi fields messages */
!  CALL GRIB_MULTI_SUPPORT_OFF()

! Loop on all the messages in a opened file

50  CONTINUE

! Read next message into IGRIB

  CALL GRIB_NEW_FROM_FILE(IFILE,IGRIB)
  IF (IGRIB == -1 )  THEN
    GOTO 60
  ENDIF

! GRIB field's parameters

  JCOUNT = JCOUNT + 1

  IF (JCOUNT ==1 ) THEN
    CALL GRIB_GET(IGRIB,'originatingCentre', ICENTRE_CODE)          ! Code of Originating Centre
    CALL GRIB_GET(IGRIB,'subCentre', ICENTRE_SUBCODE)               ! Code of Originating SubCentre
    CALL GRIB_GET(IGRIB,'generatingProcessIdentifier', IMODEL_CODE) ! Code of generating Process (Model)
  ENDIF

! Set value (-9999) for missing value in reading message

  CALL GRIB_SET(IGRIB,'missingValue',VAL_MISSING)

! Parameter index:

  CALL GRIB_GET(IGRIB,'discipline',FIELD_2D(JCOUNT) % PARAM_DISCIPL)
  CALL GRIB_GET(IGRIB,'parameterCategory',FIELD_2D(JCOUNT) % PARAM_CATEG)
  CALL GRIB_GET(IGRIB,'parameterNumber',FIELD_2D(JCOUNT) % PARAM_IND)

! Level parameters:

  CALL GRIB_GET(IGRIB,'typeOfFirstFixedSurface', FIELD_2D(JCOUNT) % LEV_TYPE)
  CALL GRIB_GET(IGRIB,'scaledValueOfFirstFixedSurface', FIELD_2D(JCOUNT) % LEV_TOP(1))
  CALL GRIB_GET(IGRIB,'scaleFactorOfFirstFixedSurface', FIELD_2D(JCOUNT) % LEV_TOP(2))
  CALL GRIB_GET(IGRIB,'scaledValueOfSecondFixedSurface', FIELD_2D(JCOUNT) % LEV_BOT(1))
  CALL GRIB_GET(IGRIB,'scaleFactorOfSecondFixedSurface', FIELD_2D(JCOUNT) % LEV_BOT(2))
  CALL GRIB_GET(IGRIB,'numberOfCoordinatesValues', FIELD_2D(JCOUNT) % N_VERT_COORD)

! Analysis date and time

  CALL GRIB_GET(IGRIB,'year',   FIELD_2D(JCOUNT) % IDATE0(1) )
  CALL GRIB_GET(IGRIB,'month',  FIELD_2D(JCOUNT) % IDATE0(2) )
  CALL GRIB_GET(IGRIB,'day',    FIELD_2D(JCOUNT) % IDATE0(3) )
  CALL GRIB_GET(IGRIB,'hour',   FIELD_2D(JCOUNT) % IDATE0(4) )
  CALL GRIB_GET(IGRIB,'minute', FIELD_2D(JCOUNT) % IDATE0(5) )

! Forecast range

  CALL GRIB_GET(IGRIB,'indicatorOfUnitOfTimeRange', FIELD_2D(JCOUNT) % IPERIOD(1) )
  CALL GRIB_GET(IGRIB,'forecastTime', FIELD_2D(JCOUNT) % IPERIOD(2) )
  CALL GRIB_GET(IGRIB,'productDefinitionTemplateNumber', FIELD_2D(JCOUNT) % IPERIOD(4) )
  FIELD_2D(JCOUNT) % IPERIOD(3)=0
  FIELD_2D(JCOUNT) % IPERIOD(5)=0
  FIELD_2D(JCOUNT) % IDATEC(:)=0
  IF (FIELD_2D(JCOUNT) % IPERIOD(4)==8) THEN
    CALL GRIB_GET(IGRIB,'lengthOfTimeRange', FIELD_2D(JCOUNT) % IPERIOD(3) )
    FIELD_2D(JCOUNT) % IPERIOD(2)=FIELD_2D(JCOUNT) % IPERIOD(2)+FIELD_2D(JCOUNT) % IPERIOD(3)
    CALL GRIB_GET(IGRIB,'typeOfStatisticalProcessing', FIELD_2D(JCOUNT) % IPERIOD(5) )
    CALL GRIB_GET(IGRIB,'yearOfEndOfOverallTimeInterval', FIELD_2D(JCOUNT) % IDATEC(1) )
    CALL GRIB_GET(IGRIB,'monthOfEndOfOverallTimeInterval', FIELD_2D(JCOUNT) % IDATEC(2) )
    CALL GRIB_GET(IGRIB,'dayOfEndOfOverallTimeInterval', FIELD_2D(JCOUNT) % IDATEC(3) )
    CALL GRIB_GET(IGRIB,'hourOfEndOfOverallTimeInterval', FIELD_2D(JCOUNT) % IDATEC(4) )
    CALL GRIB_GET(IGRIB,'minuteOfEndOfOverallTimeInterval', FIELD_2D(JCOUNT) % IDATEC(5) )
  ENDIF

!print *,'IDATE0 ', FIELD_2D(JCOUNT) % IDATE(1:5)
!print *,'IPERIOD ', FIELD_2D(JCOUNT) % IPERIOD(1:4)

! Grid parameters

  CALL GRIB_GET(IGRIB,'numberOfPointsAlongAParallel',FIELD_2D(JCOUNT) % NXGRIB)
  CALL GRIB_GET(IGRIB,'numberOfPointsAlongAMeridian',FIELD_2D(JCOUNT) % NYGRIB)
  CALL GRIB_GET(IGRIB,'longitudeOfSouthernPoleInDegrees',FIELD_2D(JCOUNT) % X00GRIB,IERR)
  IF (IERR /= 0) FIELD_2D(JCOUNT) % X00GRIB=VAL_MISSING
  CALL GRIB_GET(IGRIB,'latitudeOfSouthernPoleInDegrees',FIELD_2D(JCOUNT) % Y00GRIB,IERR)
  IF (IERR /= 0) FIELD_2D(JCOUNT) % Y00GRIB=VAL_MISSING
  CALL GRIB_GET(IGRIB,'longitudeOfFirstGridPointInDegrees',FIELD_2D(JCOUNT) % X0GRIB)
  CALL GRIB_GET(IGRIB,'latitudeOfFirstGridPointInDegrees',FIELD_2D(JCOUNT) % Y0GRIB)
  CALL GRIB_GET(IGRIB,'longitudeOfLastGridPointInDegrees',FIELD_2D(JCOUNT) % XNGRIB)
  CALL GRIB_GET(IGRIB,'latitudeOfLastGridPointInDegrees',FIELD_2D(JCOUNT) % YNGRIB)
  CALL GRIB_GET(IGRIB,'iDirectionIncrementInDegrees',FIELD_2D(JCOUNT) % DXGRIB)
  CALL GRIB_GET(IGRIB,'jDirectionIncrementInDegrees',FIELD_2D(JCOUNT) % DYGRIB)
  IF (FIELD_2D(JCOUNT) % DXGRIB < 0.) &
       FIELD_2D(JCOUNT) % DXGRIB=ABS(FIELD_2D(JCOUNT) % XNGRIB-FIELD_2D(JCOUNT) % X0GRIB)/FLOAT(FIELD_2D(JCOUNT) % NXGRIB-1)
  IF (FIELD_2D(JCOUNT) % DYGRIB < 0.) &
      FIELD_2D(JCOUNT) % DYGRIB=ABS(FIELD_2D(JCOUNT) % YNGRIB-FIELD_2D(JCOUNT) % &
      Y0GRIB)/FLOAT(FIELD_2D(JCOUNT) % NYGRIB-1)

! Dinamic allocation array for reading array data

  CALL GRIB_GET_SIZE(IGRIB,'values',NSIZE)
  ALLOCATE(VALUE(NSIZE))

! Get values array

  CALL GRIB_GET(IGRIB,'values',VALUE)

! Dinamic arrays allocation

  ALLOCATE ( FIELD_2D(JCOUNT)%FIELD (FIELD_2D(JCOUNT)%NXGRIB,FIELD_2D(JCOUNT)%NYGRIB) )

! For hybrid or sigma levels

  IF (FIELD_2D(JCOUNT) % LEV_TYPE==105.OR.FIELD_2D(JCOUNT) % LEV_TYPE==104) THEN
   NSIZE=FIELD_2D(JCOUNT) % N_VERT_COORD
   ALLOCATE ( FIELD_2D(JCOUNT) % VERT_COORD(NSIZE) )
   CALL GRIB_GET(IGRIB,'pv',FIELD_2D(JCOUNT) % VERT_COORD)
  ENDIF

! Storing of read array

! Scanning mode Flags

  CALL GRIB_GET(IGRIB,'iScansNegatively',ISCAN)
  CALL GRIB_GET(IGRIB,'jScansPositively',JSCAN)
  CALL GRIB_GET(IGRIB,'jPointsAreConsecutive',JSCANCONS)

  IF (ISCAN == 0) THEN
    IINI = 1
    IFIN = FIELD_2D(JCOUNT)%NXGRIB
    DI = 1
  ELSE
    IINI = FIELD_2D(JCOUNT)%NXGRIB
    IFIN = 1
    DI = -1
    ZZZ=FIELD_2D(JCOUNT)%X0GRIB
    FIELD_2D(JCOUNT)%X0GRIB=FIELD_2D(JCOUNT)%XNGRIB
    FIELD_2D(JCOUNT)%XNGRIB=ZZZ
  ENDIF

  IF (JSCAN == 0) THEN
    JINI = FIELD_2D(JCOUNT)%NYGRIB
    JFIN = 1
    DJ = -1
    ZZZ=FIELD_2D(JCOUNT)%Y0GRIB
    FIELD_2D(JCOUNT)%Y0GRIB=FIELD_2D(JCOUNT)%YNGRIB
    FIELD_2D(JCOUNT)%YNGRIB=ZZZ
  ELSE
    JINI = 1
    JFIN = FIELD_2D(JCOUNT)%NYGRIB
    DJ = 1
  ENDIF

  IF ( JSCANCONS == 0) THEN
    ORDER = (/1,2/)
  ELSE
    ORDER = (/2,1/)
  ENDIF

  FIELD_2D(JCOUNT)%FIELD(IINI:IFIN:DI,JINI:JFIN:DJ) =  &
          RESHAPE(VALUE, (/FIELD_2D(JCOUNT)%NXGRIB,FIELD_2D(JCOUNT)%NYGRIB/), ORDER=ORDER)

! For global data grid (when input longitude of grid points range from 0 to 360
! deg. but local grid longitude ranges from -180 to 180 deg.)

!  IF (FIELD_2D(JCOUNT) % X0GRIB >= 0..AND.X1A < 0..AND. &
!  (FIELD_2D(JCOUNT) % XNGRIB+FIELD_2D(JCOUNT) % DXGRIB) >= 359.9) THEN

  IF (FLAG_CUT_PASTE /= 0) THEN
    FIELD_2D(JCOUNT) % X0GRIB = FIELD_2D(JCOUNT) % X0GRIB -180.
    FIELD_2D(JCOUNT) % XNGRIB = FIELD_2D(JCOUNT) % XNGRIB -180.
    DO J=1,FIELD_2D(JCOUNT) % NYGRIB
      VALUE(1 : FIELD_2D(JCOUNT) % NXGRIB)=FIELD_2D(JCOUNT) % FIELD(1 : FIELD_2D(JCOUNT) % NXGRIB, J)
      FIELD_2D(JCOUNT) % FIELD(1 : INT(FIELD_2D(JCOUNT) % NXGRIB/2), J)= &
              VALUE(INT(FIELD_2D(JCOUNT) % NXGRIB/2)+1 : FIELD_2D(JCOUNT) % NXGRIB)
      FIELD_2D(JCOUNT) % FIELD(INT(FIELD_2D(JCOUNT) % NXGRIB/2)+1 : FIELD_2D(JCOUNT) % NXGRIB, J)= &
              VALUE(1 : INT(FIELD_2D(JCOUNT) % NXGRIB/2))
    ENDDO
  ENDIF

! End of grib message elaboration

 DEALLOCATE(VALUE)

  GOTO 50

! End data extraction

60 CALL GRIB_CLOSE_FILE(IFILE)

 ENDDO FILELOOP

! End of reading of input grib2 files

 NMESSAGES=JCOUNT

!-------------------------------------------------------------------------------

! Data elaboration

 DO IMESSAGE=1,NMESSAGES
  NX_OUT=SLICE_XFIN-SLICE_XINI+1
  NY_OUT=SLICE_YFIN-SLICE_YINI+1
  X0_OUT=FIELD_2D(IMESSAGE) % X0GRIB+FLOAT(SLICE_XINI-1)*FIELD_2D(IMESSAGE) % DXGRIB
  XN_OUT=FIELD_2D(IMESSAGE) % X0GRIB+FLOAT(SLICE_XFIN-1)*FIELD_2D(IMESSAGE) % DXGRIB
  Y0_OUT=FIELD_2D(IMESSAGE) % Y0GRIB+FLOAT(SLICE_YINI-1)*FIELD_2D(IMESSAGE) % DYGRIB
  YN_OUT=FIELD_2D(IMESSAGE) % Y0GRIB+FLOAT(SLICE_YFIN-1)*FIELD_2D(IMESSAGE) % DYGRIB
  ALLOCATE(FIELD_OUT(NX_OUT,NY_OUT))
  FIELD_OUT(1:NX_OUT,1:NY_OUT)=FIELD_2D(IMESSAGE) % FIELD(SLICE_XINI:SLICE_XFIN,SLICE_YINI:SLICE_YFIN)
  DEALLOCATE(FIELD_2D(IMESSAGE) % FIELD)
  ALLOCATE(FIELD_2D(IMESSAGE) % FIELD(NX_OUT,NY_OUT))
  FIELD_2D(IMESSAGE) % FIELD(:,:)=FIELD_OUT(:,:)
  DEALLOCATE(FIELD_OUT)
  FIELD_2D(IMESSAGE) % NXGRIB=NX_OUT
  FIELD_2D(IMESSAGE) % NYGRIB=NY_OUT
  FIELD_2D(IMESSAGE) % X0GRIB=X0_OUT
  FIELD_2D(IMESSAGE) % Y0GRIB=Y0_OUT
  FIELD_2D(IMESSAGE) % XNGRIB=XN_OUT
  FIELD_2D(IMESSAGE) % YNGRIB=YN_OUT
 ENDDO

!-------------------------------------------------------------------------------

! Start of wrting of output grib2 files

! A new grib message is loaded from an existing sample

 CALL GRIB_OPEN_FILE(IGRIBSAMPLE1, TRIM(GRIBSAMPLENAME), 'r')

 CALL GRIB_NEW_FROM_FILE (IGRIBSAMPLE1, IGRIBSAMPLE, IERR) 

 IF (IERR/=0) THEN
   PRINT *,'Non found the sample of GRIB file: ',GRIBSAMPLENAME(1:III-1)
   PRINT *,'Procudere writing of output file in GRIB format intrerupts'
   STOP
 ENDIF

! Open output grib file

 CALL GRIB_OPEN_FILE (IGRIBOUT, TRIM(GRIBOUTNAME), 'w')

! Set of header parameters in Sample Grib Message

! General parameters

 CALL GRIB_SET (IGRIBSAMPLE, 'editionNumber', IGRIB_EDITION) ! Format Grib Number
 CALL GRIB_SET (IGRIBSAMPLE, 'originatingCentre', ICENTRE_CODE) ! Code of Originating Centre
 CALL GRIB_SET (IGRIBSAMPLE, 'subCentre', ICENTRE_SUBCODE)      ! Code of Originating SubCentre
 CALL GRIB_SET (IGRIBSAMPLE, 'generatingProcessIdentifier', IMODEL_CODE) ! Code of generating Process (Model)  
!!! CALL GRIB_SET (IGRIBSAMPLE, 'dataRepresentationTemplateNumber', 50002) ! Compression mode
 CALL GRIB_SET (IGRIBSAMPLE, 'bitsPerValue', IBITNUM) ! Precision
!!! CALL GRIB_SET (IGRIBSAMPLE, 'secondOrderFlags', 128)
!!! CALL GRIB_SET (IGRIBSAMPLE, 'packingType', 'grid_second_order') ! Compression mode
 CALL GRIB_SET (IGRIBSAMPLE, 'missingValue', VAL_MISSING)

 LOOP_FIELDS: DO IMESSAGE=1,NMESSAGES

! Create Clone grib message form Sample grib message
! Clone grib message will be used to create output grib message

   IGRIBCLONE=0
   CALL GRIB_CLONE (IGRIBSAMPLE, IGRIBCLONE) 

! Horizontal grid parameters

  IF (INT(FIELD_2D(IMESSAGE) % X00GRIB) /= INT(VAL_MISSING).AND.INT(FIELD_2D(IMESSAGE) % Y00GRIB) /= INT(VAL_MISSING)) THEN
    CALL GRIB_SET (IGRIBCLONE, 'typeOfGrid','rotated_ll') ! Grid type is "rotated lat-lon"
    CALL GRIB_SET (IGRIBCLONE, 'longitudeOfSouthernPoleInDegrees', FIELD_2D(IMESSAGE) % X00GRIB)
    CALL GRIB_SET (IGRIBCLONE, 'latitudeOfSouthernPoleInDegrees', FIELD_2D(IMESSAGE) % Y00GRIB)
  ENDIF
  CALL GRIB_SET (IGRIBCLONE, 'numberOfPointsAlongAParallel', FIELD_2D(IMESSAGE) % NXGRIB)
  CALL GRIB_SET (IGRIBCLONE, 'numberOfPointsAlongAMeridian', FIELD_2D(IMESSAGE) % NYGRIB)
  CALL GRIB_SET (IGRIBCLONE, 'longitudeOfFirstGridPointInDegrees', FIELD_2D(IMESSAGE) % X0GRIB)
  CALL GRIB_SET (IGRIBCLONE, 'latitudeOfFirstGridPointInDegrees', FIELD_2D(IMESSAGE) % Y0GRIB)
  CALL GRIB_SET (IGRIBCLONE, 'longitudeOfLastGridPointInDegrees', FIELD_2D(IMESSAGE) % XNGRIB)
  CALL GRIB_SET (IGRIBCLONE, 'latitudeOfLastGridPointInDegrees', FIELD_2D(IMESSAGE) % YNGRIB)
  CALL GRIB_SET (IGRIBCLONE, 'iDirectionIncrementInDegrees', FIELD_2D(IMESSAGE) % DXGRIB)
  CALL GRIB_SET (IGRIBCLONE, 'jDirectionIncrementInDegrees', FIELD_2D(IMESSAGE) % DYGRIB)
  CALL GRIB_SET (IGRIBCLONE, 'iScansNegatively', 0) ! Scanning for X X0GRIB < XNGRIB
  CALL GRIB_SET (IGRIBCLONE, 'jScansPositively', 1) ! Scanning for Y Y0GRIB < YNGRIB
  CALL GRIB_SET (IGRIBCLONE, 'jPointsAreConsecutive', 0) ! Scanning order: first - X, second - Y

! Set TemplateNumber
! And forecast range, accumulation period

  CALL GRIB_SET (IGRIBCLONE, 'productDefinitionTemplateNumber', FIELD_2D(IMESSAGE) % IPERIOD(4) )
  CALL GRIB_SET (IGRIBCLONE, 'indicatorOfUnitOfTimeRange', FIELD_2D(IMESSAGE) % IPERIOD(1))

  IF (FIELD_2D(IMESSAGE) % IPERIOD(4) == 0 ) THEN ! Instant variable

    CALL GRIB_SET (IGRIBCLONE, 'forecastTime', FIELD_2D(IMESSAGE) % IPERIOD(2) )

  ELSE ! Statistical variable

    CALL GRIB_SET (IGRIBCLONE, 'forecastTime', MAX( FIELD_2D(IMESSAGE) % IPERIOD(2)-FIELD_2D(IMESSAGE) % IPERIOD(3), 0 ))
    CALL GRIB_SET (IGRIBCLONE, 'yearOfEndOfOverallTimeInterval', FIELD_2D(IMESSAGE) % IDATEC(1))
    CALL GRIB_SET (IGRIBCLONE, 'monthOfEndOfOverallTimeInterval', FIELD_2D(IMESSAGE) % IDATEC(2))
    CALL GRIB_SET (IGRIBCLONE, 'dayOfEndOfOverallTimeInterval', FIELD_2D(IMESSAGE) % IDATEC(3))
    CALL GRIB_SET (IGRIBCLONE, 'hourOfEndOfOverallTimeInterval', FIELD_2D(IMESSAGE) % IDATEC(4))
    CALL GRIB_SET (IGRIBCLONE, 'minuteOfEndOfOverallTimeInterval', FIELD_2D(IMESSAGE) % IDATEC(5))
    CALL GRIB_SET (IGRIBCLONE, 'secondOfEndOfOverallTimeInterval', 0)
    CALL GRIB_SET (IGRIBCLONE, 'typeOfStatisticalProcessing', FIELD_2D(IMESSAGE) % IPERIOD(5))
    CALL GRIB_SET (IGRIBCLONE, 'typeOfTimeIncrement', 2)
    CALL GRIB_SET (IGRIBCLONE, 'indicatorOfUnitForTimeRange', FIELD_2D(IMESSAGE) % IPERIOD(1))
    IF (FIELD_2D(IMESSAGE) % IPERIOD(2)-FIELD_2D(IMESSAGE) % IPERIOD(3) >= 0) THEN
      CALL GRIB_SET (IGRIBCLONE, 'lengthOfTimeRange', FIELD_2D(IMESSAGE) % IPERIOD(3))
    ELSE
      CALL GRIB_SET (IGRIBCLONE, 'lengthOfTimeRange', 0)
    ENDIF

  ENDIF
!print *,'Instant or Statistical variable ',FIELD_2D(IMESSAGE) % IPERIOD(4)

! Analysis date and time

  IDATEINI=FIELD_2D(IMESSAGE) % IDATE0(1)*10000+FIELD_2D(IMESSAGE) % IDATE0(2)*100+FIELD_2D(IMESSAGE) % IDATE0(3)
  ITIMEINI=FIELD_2D(IMESSAGE) % IDATE0(4)*100+FIELD_2D(IMESSAGE) % IDATE0(5)
  CALL GRIB_SET (IGRIBCLONE, 'dataDate', IDATEINI)
  CALL GRIB_SET (IGRIBCLONE, 'dataTime', ITIMEINI)

! Set Parameter Index (Code Table 4.2)

   CALL GRIB_SET (IGRIBCLONE, 'discipline', FIELD_2D(IMESSAGE) % PARAM_DISCIPL)
   CALL GRIB_SET (IGRIBCLONE, 'parameterCategory', FIELD_2D(IMESSAGE) % PARAM_CATEG)
   CALL GRIB_SET (IGRIBCLONE, 'parameterNumber', FIELD_2D(IMESSAGE) % PARAM_IND)

! Set Level parameters

   CALL GRIB_SET (IGRIBCLONE, 'typeOfFirstFixedSurface', FIELD_2D(IMESSAGE) % LEV_TYPE) 
   CALL GRIB_SET (IGRIBCLONE, 'scaledValueOfFirstFixedSurface', FIELD_2D(IMESSAGE) % LEV_TOP(1))
   CALL GRIB_SET (IGRIBCLONE, 'scaleFactorOfFirstFixedSurface', FIELD_2D(IMESSAGE) % LEV_TOP(2))
   CALL GRIB_SET (IGRIBCLONE, 'scaledValueOfSecondFixedSurface', FIELD_2D(IMESSAGE) % LEV_BOT(1))
   CALL GRIB_SET (IGRIBCLONE, 'scaleFactorOfSecondFixedSurface', FIELD_2D(IMESSAGE) % LEV_BOT(2))

  IF (FIELD_2D(IMESSAGE) % LEV_TYPE == 104.OR.FIELD_2D(IMESSAGE) % LEV_TYPE == 105) THEN ! Sigma level or Hybrid level
     NSIZE=FIELD_2D(IMESSAGE) % N_VERT_COORD
     CALL GRIB_SET (IGRIBCLONE, 'PVPresent', 1)
     CALL GRIB_SET (IGRIBCLONE, 'numberOfVerticalCoordinateValues', NSIZE)
     ALLOCATE(VALUE(NSIZE))
     VALUE(:)=FIELD_2D(IMESSAGE) % VERT_COORD(:)
     CALL GRIB_SET (IGRIBCLONE, 'pv', VALUE)
     DEALLOCATE(VALUE)
  ELSE
     CALL GRIB_SET (IGRIBCLONE, 'PVPresent', 0)
     CALL GRIB_SET (IGRIBCLONE, 'numberOfVerticalCoordinateValues', 0)
  ENDIF

! Set 2D Data Field
  NSIZE=FIELD_2D(IMESSAGE) % NXGRIB*FIELD_2D(IMESSAGE) % NYGRIB
  ALLOCATE(VALUE(NSIZE))
  DO J=1,FIELD_2D(IMESSAGE) % NYGRIB
  DO I=1,FIELD_2D(IMESSAGE) % NXGRIB
    III=(J-1)*FIELD_2D(IMESSAGE) % NXGRIB+I
    VALUE(III)=FIELD_2D(IMESSAGE) % FIELD(I,J)
  ENDDO
  ENDDO
  CALL GRIB_SET (IGRIBCLONE, 'values', VALUE)
  DEALLOCATE(VALUE)

! Write grib message

  CALL GRIB_WRITE (IGRIBCLONE, IGRIBOUT)
  CALL GRIB_RELEASE (IGRIBCLONE)
 
 ENDDO LOOP_FIELDS

! Close output grib file

 CALL GRIB_CLOSE_FILE (IGRIBOUT)

STOP
END
!=======================================================================
