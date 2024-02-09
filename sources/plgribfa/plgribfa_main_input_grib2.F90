! PROCEDURE READS GRIB2 FORMAT DATA AND CREATES GRAPHICAL 
! PRESENTATION FOR 2D FIELDS
!------------------------------------------------------------------------------
PROGRAM PLGRIBFA_MAIN_INPUT_GRIB2

USE eccodes
USE MOD_GRAPH_PAR
USE MOD_GRAPH_PAR_DEF

IMPLICIT NONE

TYPE DATA
  INTEGER :: PARAM_DISCIPL, PARAM_CATEG, PARAM_IND, LEV_TYPE, LEV_TOP(2), LEV_BOT(2), IDATE(5), IPERIOD(4), &
  VERT_COORD_TYPE
  REAL, DIMENSION(:,:), ALLOCATABLE :: FIELD
  CHARACTER(LEN=60) :: NAME_PAR, NAME_LEV
END TYPE DATA
TYPE (DATA), DIMENSION(1000) :: READ_2D

REAL, PARAMETER :: TZER=273.15
CHARACTER :: PPF_TYPE*3, DATAFILE_NAME*14, GRIDTYPE*20, SYSTEM_COMMAND*60,&
 CONTROL1*1, CONTROL2*1, CONTROL3*1, DIR_PACKAGE*60, TABLE_NAME*100, LINE*100, &
 A3*3, A10*10, A1*1, A60*60, A61*60, A120*120, A121*120
INTEGER :: N_DATAFILE, IFILE, FILE_IND, IMESSAGE, MESSAGE, ICENTRE_CODE, ICENTRE, ISUBCENTRE, &
 IND_TEMPLATE_GRID, IND_TEMPLATE_PROD, &
 FLAG_CONST_FIELD, FLAG_STATISTIC, FLAG_STATISTIC_MODE, FLAG_MISSING, INDEX1, INDEX2, NFILE, &
 LEV_BOT_TYPE, III, IPAGE, NCHAR, IFIELD, IFL_CONTOUR, ICONTOUR, NVALUES, ISCAN, IEMOD,&
 I, J, IR, JR, IORDER, IDATEINI, ITIMEINI, ISTATISTIC, LSTATISTIC, LSTATISTIC2,&
 I1, I2, I3, NPAGE_TOT, NFIELD1, I_HALF, NX_CROSS=0, NY_CROSS=0
REAL :: RIS, DIS, ZMIN, ZMAX, ZZZ, ZSMOOTH=0.5, CWN, CWL, &
 LON_CROSS_INI, LON_CROSS_FIN, LAT_CROSS_INI, LAT_CROSS_FIN
REAL, DIMENSION(:), ALLOCATABLE :: VALUE, ZOROG
REAL, DIMENSION(:,:), ALLOCATABLE :: FIELD_WORK
REAL*8 :: ZFAC, ZX0, ZY0, ZLON, ZLAT, ZZLAT, ZAARG, ZARG
INTEGER, DIMENSION(10) :: IND

PI = ABS(ACOS(-1.))
X0=VALMISS
Y0=VALMISS
LON_CROSS_INI=VALMISS
LON_CROSS_FIN=VALMISS
LAT_CROSS_INI=VALMISS
LAT_CROSS_FIN=VALMISS

PRINT *,'Enter full path name of eccodes package directory'
PRINT *,'for example: /home/dinamica/BOLAM/PROGRAM/ECCODES/eccodes-1.1.0'
READ (*,'(A60)') DIR_PACKAGE
PRINT *,'Pass name of eccodes package directory: ',DIR_PACKAGE

NCHAR=0
DO I=60,1,-1
 IF (DIR_PACKAGE(I:I).EQ.' ') NCHAR=I
ENDDO
NCHAR=NCHAR-1

PRINT *,'Enter number of data files'
READ (*,*) N_DATAFILE
PRINT *,'Data files number ',N_DATAFILE

PRINT *,'Enter 3 flags:' 
PRINT *,'First data file is constant (static) data (0/1), Statistical elaboration of data, &
 Mode of statistical elaboration: &
 1 - accumulation, 2 - average, 3 - maximum, 4 – minimum'
READ (*,*) FLAG_CONST_FIELD, FLAG_STATISTIC, FLAG_STATISTIC_MODE
IF (FLAG_CONST_FIELD==0) THEN
  PRINT *,'The first data file do not contatin constant (static) field data' 
ELSE
  PRINT *,'The first data file contatins constant (static) field data' 
ENDIF
IF (FLAG_STATISTIC==0) THEN
  PRINT *,'No statistical elaboration of read data'
ELSE
  PRINT *,'Statistical elaboration of read datai will be done:'
  IF (FLAG_STATISTIC_MODE==1) PRINT *,'      accumulated value'
  IF (FLAG_STATISTIC_MODE==2) PRINT *,'      average value'
  IF (FLAG_STATISTIC_MODE==3) PRINT *,'      maximum value'
  IF (FLAG_STATISTIC_MODE==4) PRINT *,'      minimum value'
ENDIF

IF (FLAG_STATISTIC==0) THEN
  NFILE=1
ELSE
  IF (FLAG_CONST_FIELD==0) THEN
    NFILE=N_DATAFILE
  ELSE
    NFILE=N_DATAFILE-1
  ENDIF
ENDIF

PRINT *,'Enter flag to use the Stereographic geo. Projection (0/1)',ISTEREO 
READ (*,*) ISTEREO
IF (ISTEREO==0) THEN
  PRINT *,'You choose to not using stereographic geo. projection'
ELSE
  PRINT *,'You choose using stereographic geo. projection'
ENDIF

PRINT *,'Enter stereographic geo. Projection parameters: &
 origin latitude and longitude, projection radius, extreem latitude, &
 (default: 90. 0. 1. 15.)'
READ *,STEREOLAT0, STEREOLON0, STEREORAD, STEREOLATN
PRINT *,'Defined stereographic geo. projection parameters are following: ',STEREOLAT0, STEREOLON0, STEREORAD, STEREOLATN

PRINT *,'Enter plot divice type: psc=1, pscairo=2, pdfcairo=3, xwin=4, png=5, pngcairo=6, jpeg=7, gif=8'
READ (*,*) IDEV_TYPE
PRINT *,'Defined plot divise type is ',IDEV_TYPE

PRINT *,'Enter flag for plotting of geographical lines (0/1): &
 Coast line, State boundaries, Rivers, Administrative boundaries'
READ (*,*) IFL_LINE(1:4)
PRINT *,'Defined flags for ploting of geographical lines are: ',IFL_LINE(1:4)

PRINT *,'Enter relative zoom boundary (in x and in y), for full field plot: 0. 1. 0. 1.' 
READ (*,*) XZOOM_INI, XZOOM_FIN, YZOOM_INI, YZOOM_FIN
PRINT *,'Defined relative zoom boundary are: ',XZOOM_INI, XZOOM_FIN, YZOOM_INI, YZOOM_FIN
   
PRINT *,'Do You want to use automatic definition of fields graphic presentation parameters: &
 color filling, 1 field in 1 page, automatic definition color palette, color number, interval etc.? (y/n)'
READ (*,'(A)') CONTROL1
IF (CONTROL1=='y') THEN
  PRINT *,'You choos automatic definition of fields graphic presentation parameters'
ELSE
  PRINT *,'You will defined all parameter for fields graphic presentation'
ENDIF

PRINT *,'Do You want to define procedure parameters manualy? (y/n)'
READ (*,'(A)') CONTROL2

IF (CONTROL2=='n') THEN
  PRINT *,'You choos to define procedure by parameters memorized in definition.txt file'
  OPEN (100,FILE='definition.txt',STATUS='OLD',FORM='FORMATTED')
  OPEN (101,FILE='namelist_graf_def',STATUS='OLD')
  DO I=1,3
    READ (100,*)
  ENDDO
  PRINT *,'The definition.txt file is valid for all input data files? (y/n)'
  READ (*,'(A)') CONTROL3
  IF (FLAG_CONST_FIELD==1) THEN
    OPEN (200,FILE='definition_const.txt',STATUS='OLD',FORM='FORMATTED')
    OPEN (201,FILE='namelist_graf_def_const',STATUS='OLD')
    DO I=1,3
      READ (200,*)
    ENDDO
    ENDIF
  IF (CONTROL3=='y') THEN
    PRINT *,'Data parameters defined in definition.txt will be applied for all input files'
  ELSE
    PRINT *,'Data parameters will be defined differently for differen input files according definition.txt file' 
  ENDIF
ENDIF

IMESSAGE=0
NPAGE=0
NPAGE_TOT=0

DATAFILE_LOOP : DO IFILE=1,N_DATAFILE

  INDEX1=100
  INDEX2=101
  IF (IFILE==1.AND.FLAG_CONST_FIELD==1) THEN
    INDEX1=200
    INDEX2=201
  ENDIF

  IF ((FLAG_CONST_FIELD==1.AND.IFILE>1.AND.CONTROL3=='y').OR.&
      (FLAG_CONST_FIELD==0.AND.CONTROL3=='y')) THEN
    REWIND (INDEX1)
    REWIND (INDEX2)
    DO I=1,3
      READ (INDEX1,*)
    ENDDO
  ENDIF

! Open file

  WRITE (DATAFILE_NAME,'(A,I3.3)') 'file_grib2_',IFILE
  CALL GRIB_OPEN_FILE(FILE_IND,DATAFILE_NAME,"r")

! Turn on support for multi fields messages !!! Very important !!!

  CALL GRIB_MULTI_SUPPORT_ON()

! Loop on all the messages in a opened file

  !IMESSAGE=0

10  CONTINUE

! Read next message into MESSAGE

  CALL GRIB_NEW_FROM_FILE(FILE_IND,MESSAGE)
  IF (MESSAGE == -1 )  THEN
    GOTO 20
  ENDIF

  IMESSAGE=IMESSAGE+1

! Set value (-9999) for missing value in reading message

  CALL GRIB_SET(MESSAGE,'missingValue',VALMISS)

  IF (IFILE==1.AND.IMESSAGE==1) THEN

! General parameters

    CALL GRIB_GET(MESSAGE,'originatingCentre', ICENTRE) ! Code of Originating Centre
    CALL GRIB_GET(MESSAGE,'subCentre', ISUBCENTRE) ! Code of Originating SubCentre
    ICENTRE_CODE=ICENTRE
    IF (ICENTRE == 80 ) THEN ! Rome RSMC
      ICENTRE_CODE=ICENTRE*1000+ISUBCENTRE
    ENDIF
    CALL GRIB_GET(MESSAGE,'generatingProcessIdentifier', IMODEL_CODE) ! Code of generating Process (Model)
    CALL GRIB_GET(MESSAGE,'gridDefinitionTemplateNumber', IND_TEMPLATE_GRID) ! Index of template for definition grid type

    IF (IND_TEMPLATE_GRID <= 3) THEN ! Latitude-longitude grid

! Grid parameters

      CALL GRIB_GET(MESSAGE,'numberOfPointsAlongAParallel', NX)
      CALL GRIB_GET(MESSAGE,'numberOfPointsAlongAMeridian', NY)
      CALL GRIB_GET(MESSAGE,'longitudeOfFirstGridPointInDegrees', X00)
      CALL GRIB_GET(MESSAGE,'latitudeOfFirstGridPointInDegrees', Y00)
      CALL GRIB_GET(MESSAGE,'longitudeOfLastGridPointInDegrees', XNN)
      CALL GRIB_GET(MESSAGE,'latitudeOfLastGridPointInDegrees', YNN)
      CALL GRIB_GET(MESSAGE,'iDirectionIncrementInDegrees', DX)
      CALL GRIB_GET(MESSAGE,'jDirectionIncrementInDegrees', DY)
      IF (YNN < Y00) THEN
        ZZZ=Y00
        Y00=YNN
        YNN=ZZZ 
      ENDIF
      IF (X00>180.) X00=X00-360.
      IF ( DX < 0.) DX=ABS(XNN-X00)/FLOAT(NX-1)
      IF ( DY < 0.) DY=ABS(YNN-Y00)/FLOAT(NY-1)

      ALLOCATE(ALON(NX,NY))
      ALLOCATE(ALAT(NX,NY))
      ALLOCATE(LSM(NX,NY))
      ALON(:,:)=VALMISS
      ALAT(:,:)=VALMISS
      LSM(:,:)=VALMISS
      ALLOCATE(FIELD_WORK(NX,NY))  
  
      CALL GRIB_GET(MESSAGE,'typeOfGrid',GRIDTYPE)
      IF (GRIDTYPE(1:7)=='rotated') THEN
        CALL GRIB_GET(MESSAGE,'longitudeOfSouthernPoleInDegrees', X0)
        CALL GRIB_GET(MESSAGE,'latitudeOfSouthernPoleInDegrees', Y0)
        Y0 = ACOS(-SIN(Y0*PI/180.))*180./PI
        ZFAC = DABS(DACOS(-1.D0))/180.D0
        ZX0  = DBLE(X0)*ZFAC
        ZY0  = DBLE(Y0)*ZFAC
        DO J=1,NY
          ZLAT=(DBLE(Y00)+(J-1)*DBLE(DY))*ZFAC
          DO I=1,NX
            ZLON = (DBLE(X00)+(I-1)*DBLE(DX))*ZFAC
            ZZLAT= 1.D0/ZFAC*DASIN( DCOS(ZY0)*DSIN(ZLAT) +  &
                         DSIN(ZY0)*DCOS(ZLAT)*DCOS(ZLON) )
            ZARG = -DSIN(ZLAT)*DSIN(ZY0)+DCOS(ZY0)*DCOS(ZLAT)*DCOS(ZLON)
            ZAARG = ZARG/DCOS(ZFAC*ZZLAT)
            ALAT(I,J)=SNGL(ZZLAT)
            IF (ZAARG < -1.D0.AND.ZAARG > -1.00001D0) ZAARG = -1.D0
            IF (ZAARG >  1.D0.AND.ZAARG <  1.00001D0) ZAARG =  1.D0
            IF (ZLON < 0.D0) THEN
              ALON(I,J) = SNGL(1.D0/ZFAC*(ZX0-DACOS(ZAARG)))
            ELSE
              ALON(I,J) = SNGL(1.D0/ZFAC*(ZX0+DACOS(ZAARG)))
            ENDIF
          ENDDO
        ENDDO
      ELSE
        X0=0.
        Y0=0.
        DO I=1,NX
          IF ((X00+FLOAT(I-1)*DX)<=180.) THEN
            ALON(I,:)=X00+FLOAT(I-1)*DX
          ELSE
            ALON(I,:)=X00+FLOAT(I-1)*DX-360.
          ENDIF
        ENDDO
        DO J=1,NY
          ALAT(:,J)=Y00+FLOAT(J-1)*DY
        ENDDO
  
  ! Case of zoom in global area data
  
        IF (XZOOM_INI>XZOOM_FIN) THEN
          I_HALF=NX/2
          I1=I_HALF
          DO I=1,I_HALF
            I1=I1+1
            FIELD_WORK(I,:)=ALON(I1,:)
          ENDDO
          I1=0
          DO I=I_HALF+1,NX
            I1=I1+1
            FIELD_WORK(I,:)=ALON(I1,:)
          ENDDO
          ALON(:,:)=FIELD_WORK(:,:)
          X00=ALON(1,1)
        ENDIF
  
      ENDIF

      PRINT *
      PRINT *,'   Input data defined on lat-lon grid'
      PRINT *

    ELSEIF (IND_TEMPLATE_GRID == 1000) THEN ! Cross-section in space

      CALL GRIB_GET(MESSAGE,'numberOfVerticalPoints', NY)
      CALL GRIB_GET(MESSAGE,'meaningOfVerticalCoordinate', VERT_COORD_TYPE)
      ALLOCATE(VERT_COORD(NY))
      VERT_COORD(:)=VALMISS
      CALL GRIB_GET(MESSAGE,'verticalCoordinate', I)
      IF (I == 0) THEN ! Explicit coordinate values set
        CALL GRIB_GET(MESSAGE,'pv', VERT_COORD)
        !!!CALL GRIB_GET(MESSAGE,'????', VERT_COORD) ! In eccodes not defined (definitions/grib2/template.3.1000.def)
!!! This problem is avoided by using pv in section 4 (template 0)
      ENDIF 

      TABLE_NAME=DIR_PACKAGE(1:NCHAR)//'/definitions/grib2/tables/3/3.15.table'
      OPEN (13,FILE=TABLE_NAME,STATUS='OLD',FORM='FORMATTED')
      13 LINE(1:100)=' ' 
      READ (13,'(A100)',END=131) LINE
      IF (LINE(1:1)/='#') THEN
        I1=0
        I2=0
        DO I=1,100
          IF (LINE(I:I)==' ') THEN
            IF (I1==0) THEN
              I1=I
            ELSE
              I2=I
              GOTO 130
            ENDIF
          ENDIF
        ENDDO
    130 READ (LINE(1:I1-1),*,END=131) III
        IF (III == VERT_COORD_TYPE) THEN
          READ (LINE(I2+1:100),'(A60)') NAME_VERT_COORD
          GOTO 131
        ENDIF
      ENDIF
      GOTO 13
      131 CLOSE (13)

      PRINT *
      PRINT *,'   Input data defined on space cross-section grid with vertical coordinate in ', &
      NAME_VERT_COORD
      PRINT *

    ELSE

      PRINT *,'Not provided grid template index ',IND_TEMPLATE_GRID
      STOP

    ENDIF

  ENDIF ! IMASSAGE==1.AND.IFILE==1

! Grid parameters

  IF (IND_TEMPLATE_GRID == 1000) THEN ! Cross-section in space

    CALL GRIB_GET(MESSAGE,'numberOfHorizontalPoints', NX_CROSS)
    CALL GRIB_GET(MESSAGE,'longitudeOfFirstGridPoint', LON_CROSS_INI)
    CALL GRIB_GET(MESSAGE,'latitudeOfFirstGridPoint', LAT_CROSS_INI)
    CALL GRIB_GET(MESSAGE,'longitudeOfLastGridPoint', LON_CROSS_FIN)
    CALL GRIB_GET(MESSAGE,'latitudeOfLastGridPoint', LAT_CROSS_FIN)
    LON_CROSS_INI=LON_CROSS_INI*1.E-6
    LON_CROSS_FIN=LON_CROSS_FIN*1.E-6
    LAT_CROSS_INI=LAT_CROSS_INI*1.E-6
    LAT_CROSS_FIN=LAT_CROSS_FIN*1.E-6
    IF (LON_CROSS_INI > 180.) THEN
      LON_CROSS_INI=LON_CROSS_INI-360.
      IF (LAT_CROSS_FIN > 180.) LON_CROSS_FIN=LON_CROSS_INI-360.
    ENDIF
    NX=NX_CROSS

  ENDIF ! IND_TEMPLATE_GRID == 1000

! Coding template number

  CALL GRIB_GET(MESSAGE,'productDefinitionTemplateNumber',IND_TEMPLATE_PROD)

! Parameter index:

  CALL GRIB_GET(MESSAGE,'discipline',READ_2D(IMESSAGE)%PARAM_DISCIPL)
  CALL GRIB_GET(MESSAGE,'parameterCategory',READ_2D(IMESSAGE)%PARAM_CATEG)
  CALL GRIB_GET(MESSAGE,'parameterNumber',READ_2D(IMESSAGE)%PARAM_IND)

  TITLE_DEF=&
'########################################################################################################################'

! Parameter decoding

!!!  CALL GRIB_GET(MESSAGE,'name',READ_2D(IMESSAGE)%NAME_PAR)
!!!  CALL GRIB_GET(MESSAGE,'units',A10)

  READ_2D(IMESSAGE)%NAME_PAR='Unknown parameter code'
  WRITE (TABLE_NAME,'(2A,I0,A,I0,A)') DIR_PACKAGE(1:NCHAR),'/definitions/grib2/tables/2/4.2.',&
 READ_2D(IMESSAGE)%PARAM_DISCIPL,'.',READ_2D(IMESSAGE)%PARAM_CATEG,'.table'
  OPEN (11,FILE=TABLE_NAME,STATUS='OLD',FORM='FORMATTED')
  11 LINE(1:100)=' ' 
  READ (11,'(A100)',END=111) LINE
  IF (LINE(1:1)/='#') THEN
    I1=0
    I2=0
    DO I=1,100
      IF (LINE(I:I)==' ') THEN
        IF (I1==0) THEN
          I1=I
        ELSE
          I2=I
          GOTO 110
        ENDIF
      ENDIF
    ENDDO
110 READ (LINE(1:I1-1),*,END=111) III
    IF (III==READ_2D(IMESSAGE)%PARAM_IND) THEN
      READ (LINE(I2+1:100),'(A60)') READ_2D(IMESSAGE)%NAME_PAR
      GOTO 111
    ENDIF
  ENDIF
  GOTO 11
  111 CLOSE (11)

! Level decoding

  IF (IND_TEMPLATE_GRID <= 3) THEN ! Latitude-longitude grid

! Level parameters:

    IF (IND_TEMPLATE_PROD < 30) THEN
  
  ! Meteorological and radar data
  
      CALL GRIB_GET(MESSAGE,'typeOfFirstFixedSurface',READ_2D(IMESSAGE)%LEV_TYPE)
      CALL GRIB_GET(MESSAGE,'typeOfSecondFixedSurface',LEV_BOT_TYPE)
      CALL GRIB_GET(MESSAGE,'scaledValueOfFirstFixedSurface',READ_2D(IMESSAGE)%LEV_TOP(1))
      CALL GRIB_GET(MESSAGE,'scaleFactorOfFirstFixedSurface',READ_2D(IMESSAGE)%LEV_TOP(2))
      CALL GRIB_GET(MESSAGE,'scaledValueOfSecondFixedSurface',READ_2D(IMESSAGE)%LEV_BOT(1))
      CALL GRIB_GET(MESSAGE,'scaleFactorOfSecondFixedSurface',READ_2D(IMESSAGE)%LEV_BOT(2))
      CALL CODES_IS_MISSING(MESSAGE,'scaledValueOfSecondFixedSurface',FLAG_MISSING)
      IF (FLAG_MISSING==1) READ_2D(IMESSAGE)%LEV_BOT(1:2)=INT(VALMISS)
      CALL CODES_IS_MISSING(MESSAGE,'typeOfSecondFixedSurface',FLAG_MISSING)
      IF (FLAG_MISSING==1) READ_2D(IMESSAGE)%LEV_BOT(1:2)=INT(VALMISS)
      CALL CODES_IS_MISSING(MESSAGE,'scaledValueOfFirstFixedSurface',FLAG_MISSING)
      IF (FLAG_MISSING==1) READ_2D(IMESSAGE)%LEV_TOP(1:2)=INT(VALMISS)
      CALL CODES_IS_MISSING(MESSAGE,'typeOfFirstFixedSurface',FLAG_MISSING)
      IF (FLAG_MISSING==1) READ_2D(IMESSAGE)%LEV_TOP(1:2)=INT(VALMISS)
  
    ELSE
  
  ! Satelite data
  
      READ_2D(IMESSAGE)%LEV_TYPE = INT(VALMISS)
      CALL GRIB_GET(MESSAGE,'scaledValueOfCentralWaveNumber',READ_2D(IMESSAGE)%LEV_TOP(1))
      CALL GRIB_GET(MESSAGE,'scaleFactorOfCentralWaveNumber',READ_2D(IMESSAGE)%LEV_TOP(2))
  
    ENDIF

    READ_2D(IMESSAGE)%NAME_LEV='Unknown level code'
    TABLE_NAME=DIR_PACKAGE(1:NCHAR)//'/definitions/grib2/tables/2/4.5.table'
    OPEN (12,FILE=TABLE_NAME,STATUS='OLD',FORM='FORMATTED')
    12 LINE(1:100)=' ' 
    READ (12,'(A100)',END=121) LINE
    IF (LINE(1:1)/='#') THEN
      I1=0
      I2=0
      DO I=1,100
        IF (LINE(I:I)==' ') THEN
          IF (I1==0) THEN
            I1=I
          ELSE
            I2=I
            GOTO 120
          ENDIF
        ENDIF
      ENDDO
  120 READ (LINE(1:I1-1),*,END=121) III
      IF (III==READ_2D(IMESSAGE)%LEV_TYPE) THEN
        READ (LINE(I2+1:100),'(A60)') READ_2D(IMESSAGE)%NAME_LEV
        GOTO 121
      ENDIF
    ENDIF
    GOTO 12
    121 CLOSE (12)
  
    IF (READ_2D(IMESSAGE)%NAME_LEV(1:9)=='Specified') THEN
      DO I=1,60
        A61(I:I)=' '
      ENDDO
      A61(1:37)=READ_2D(IMESSAGE)%NAME_LEV(24:60)
      READ_2D(IMESSAGE)%NAME_LEV=A61
    ENDIF
    IF (READ_2D(IMESSAGE)%NAME_LEV(1:8)=='Specific') THEN
      DO I=1,60
        A61(I:I)=' '
      ENDDO
      A61(1:42)=READ_2D(IMESSAGE)%NAME_LEV(19:60)
      READ_2D(IMESSAGE)%NAME_LEV=A61
    ENDIF

  ELSEIF (IND_TEMPLATE_GRID == 1000) THEN ! Cross-cestion in space

    READ_2D(IMESSAGE)%LEV_TOP(:)=VALMISS
    READ_2D(IMESSAGE)%LEV_BOT(:)=VALMISS
    READ_2D(IMESSAGE)%NAME_LEV=' '

  ENDIF ! IND_TEMPLATE_GRID

  PRINT *
  PRINT *,'PARAMETER    ',READ_2D(IMESSAGE)%NAME_PAR
  IF (IND_TEMPLATE_GRID <= 3) THEN ! Latitude-longitude grid
    IF (READ_2D(IMESSAGE)%LEV_TOP(1) /= INT(VALMISS)) THEN
      PRINT *,'LEVEL        ',READ_2D(IMESSAGE)%NAME_LEV(1:30),' ',&
 READ_2D(IMESSAGE)%LEV_TOP(1)*0.1**READ_2D(IMESSAGE)%LEV_TOP(2)
    ELSE
      PRINT *,'LEVEL     Missing value '
    ENDIF
    IF (READ_2D(IMESSAGE)%LEV_BOT(1) /= INT(VALMISS))&
 PRINT *,' SECOND LEVEL            ',READ_2D(IMESSAGE)%LEV_BOT(1)*0.1**READ_2D(IMESSAGE)%LEV_BOT(2)
  ENDIF

  PRINT *,'Do You want to create graphic presentation for this field? (y/n)'
  IF (CONTROL2=='y') THEN
    READ (*,*) A1
  ELSE
    READ (INDEX1,*) A1
  ENDIF

  IF (A1=='n'.AND.&
(READ_2D(IMESSAGE)%PARAM_DISCIPL*100000+READ_2D(IMESSAGE)%PARAM_CATEG*1000+READ_2D(IMESSAGE)%PARAM_IND) /= 200000) GOTO 10 ! Pass to next message

! Orography for space cross-section (1-D matrix instead of "normal" 2-D matrix)

  IF (IND_TEMPLATE_GRID == 1000.AND.&
((READ_2D(IMESSAGE)%PARAM_DISCIPL*100000+READ_2D(IMESSAGE)%PARAM_CATEG*1000+READ_2D(IMESSAGE)%PARAM_IND) == 200007.OR.&
(READ_2D(IMESSAGE)%PARAM_DISCIPL*100000+READ_2D(IMESSAGE)%PARAM_CATEG*1000+READ_2D(IMESSAGE)%PARAM_IND) == 003000)) THEN
    NY_CROSS=NY
    NY=1
    READ_2D(IMESSAGE)%NAME_PAR=' ' ! No name for orography in space cross-section
  ENDIF

! Reading selected filed and definition of parameters for its presentation

  CALL GRIB_GET_SIZE(MESSAGE,'values',NVALUES)
  ALLOCATE (VALUE(NVALUES))
  CALL GRIB_GET(MESSAGE,'values',VALUE)

! Scanning mode Flags

  CALL GRIB_GET(MESSAGE,'scanningMode',ISCAN)

! If BTEST(ISCAN,7).EQ..F.   then Points scan in +i direction
! If BTEST(ISCAN,7).EQ..T.   then Points scan in -i direction
! If BTEST(ISCAN,8).EQ..F.   then Points scan in -j direction
! If BTEST(ISCAN,8).EQ..T.   then Points scan in +j direction
! If BTEST(ISCAN,9).EQ..F.   then Adjacent points in i direction are consecutive
! If BTEST(ISCAN,9).EQ..T.   then Adjacent points in j direction are consecutive
! i direction is defined as west to east along a parallel, or left to right along an X-axis.
!j direction is defined as south to north along a meridian, or bottom to top along a Y-axis.

  IF ((.NOT.BTEST(ISCAN,7)).AND.(.NOT.BTEST(ISCAN,6)).AND.(.NOT.BTEST(ISCAN,5))) IEMOD=111
  IF ((.NOT.BTEST(ISCAN,7)).AND.(BTEST(ISCAN,6)).AND.(.NOT.BTEST(ISCAN,5))) IEMOD=121

  ALLOCATE (READ_2D(IMESSAGE)%FIELD(NX,NY))

  IF (IEMOD==111) THEN
    DO J=1,NY
    DO I=1,NX
    I1=I+(J-1)*NX
    IR=I
    JR=NY-J+1
    READ_2D(IMESSAGE)%FIELD(IR,JR)=VALUE(I1)
    ENDDO
    ENDDO
  ENDIF
  IF (IEMOD==121) THEN
    DO J=1,NY
    DO I=1,NX
    I1=I+(J-1)*NX
    IR=I
    JR=J
    READ_2D(IMESSAGE)%FIELD(IR,JR)=VALUE(I1)
    ENDDO
    ENDDO
  ENDIF

  DEALLOCATE (VALUE)

! Orography for space cross-section

  IF (IND_TEMPLATE_GRID == 1000.AND.&
((READ_2D(IMESSAGE)%PARAM_DISCIPL*100000+READ_2D(IMESSAGE)%PARAM_CATEG*1000+READ_2D(IMESSAGE)%PARAM_IND) == 200007.OR.&
(READ_2D(IMESSAGE)%PARAM_DISCIPL*100000+READ_2D(IMESSAGE)%PARAM_CATEG*1000+READ_2D(IMESSAGE)%PARAM_IND) == 003000)) THEN
    NY=NY_CROSS
    IF (ALLOCATED(ZOROG)) DEALLOCATE (ZOROG)
    ALLOCATE (ZOROG(NX))
    ZOROG(:)=READ_2D(IMESSAGE)%FIELD(:,1)
    DEALLOCATE (READ_2D(IMESSAGE)%FIELD)
    ALLOCATE (READ_2D(IMESSAGE)%FIELD(NX,NY))
! If in cross-section vertical coordinate is pressure, then orography data
! is substituted by surface pressure;
! if in cross-section vertical coordinate is altitude, then ororgaphy data
! is true orography altitude.
    READ_2D(IMESSAGE)%PARAM_DISCIPL=2
    READ_2D(IMESSAGE)%PARAM_CATEG=0
    READ_2D(IMESSAGE)%PARAM_IND=7
    IF (VERT_COORD_TYPE == 100) THEN ! Pressure
      DO I=1,NX
        READ_2D(IMESSAGE)%FIELD(I,:)=0.
        DO J=1,NY
          IF (VERT_COORD(J) > ZOROG(I)) READ_2D(IMESSAGE)%FIELD(I,J)=1.
        ENDDO
      ENDDO
    ELSE
      DO I=1,NX
        READ_2D(IMESSAGE)%FIELD(I,:)=0.
        DO J=1,NY
          IF (VERT_COORD(J) < ZOROG(I)) READ_2D(IMESSAGE)%FIELD(I,J)=1.
        ENDDO
      ENDDO
    ENDIF
  ENDIF

! Case of zoom in global area data

  IF (.NOT.ALLOCATED(FIELD_WORK)) ALLOCATE(FIELD_WORK(NX,NY))  

  IF (XZOOM_INI>XZOOM_FIN) THEN
    I_HALF=NX/2
    I1=I_HALF
    DO I=1,I_HALF
      I1=I1+1
      FIELD_WORK(I,:)=READ_2D(IMESSAGE)%FIELD(I1,:)
    ENDDO
    I1=0
    DO I=I_HALF+1,NX
      I1=I1+1
      FIELD_WORK(I,:)=READ_2D(IMESSAGE)%FIELD(I1,:)
    ENDDO
    READ_2D(IMESSAGE)%FIELD(:,:)=FIELD_WORK(:,:)
  ENDIF

  IF (A1 == 'n'.AND.IND_TEMPLATE_GRID <= 3.AND.&
(READ_2D(IMESSAGE)%PARAM_DISCIPL*100000+READ_2D(IMESSAGE)%PARAM_CATEG*1000+READ_2D(IMESSAGE)%PARAM_IND) == 200000) THEN
    LSM(:,:)=MIN(MAX(1.-READ_2D(IMESSAGE)%FIELD(:,:),0.),1.)
    GOTO 10
  ENDIF

 IF (READ_2D(IMESSAGE)%NAME_PAR(1:11)=='Temperature'.OR.READ_2D(IMESSAGE)%NAME_PAR(1:16)=='Soil temperature'.OR.&
 READ_2D(IMESSAGE)%NAME_PAR(19:29)=='temperature'.OR.READ_2D(IMESSAGE)%NAME_PAR(22:32)=='temperature'.OR.&
 (READ_2D(IMESSAGE)%NAME_PAR(11:21)=='temperature'.AND.READ_2D(IMESSAGE)%NAME_PAR(1:21)/='Potential temperature')) THEN
   IF (READ_2D(IMESSAGE)%NAME_PAR(1:11)=='Temperature') READ_2D(IMESSAGE)%NAME_PAR(14:14)='C'
   IF (READ_2D(IMESSAGE)%NAME_PAR(1:16)=='Soil temperature') READ_2D(IMESSAGE)%NAME_PAR(19:19)='C'
   IF (READ_2D(IMESSAGE)%NAME_PAR(19:29)=='temperature') READ_2D(IMESSAGE)%NAME_PAR(32:32)='C'
   IF (READ_2D(IMESSAGE)%NAME_PAR(22:32)=='temperature') READ_2D(IMESSAGE)%NAME_PAR(35:35)='C'
   IF (READ_2D(IMESSAGE)%NAME_PAR(11:21)=='temperature') READ_2D(IMESSAGE)%NAME_PAR(24:24)='C' ! Dew point temperature
   READ_2D(IMESSAGE)%FIELD(:,:)=READ_2D(IMESSAGE)%FIELD(:,:)-TZER
 ENDIF
! IF (READ_2D(IMESSAGE)%NAME_PAR(9:14)=='lifted') THEN
!   READ_2D(IMESSAGE)%NAME_PAR(23:23)='C'
!   READ_2D(IMESSAGE)%FIELD(:,:)=READ_2D(IMESSAGE)%FIELD(:,:)-TZER
! ENDIF
 IF (READ_2D(IMESSAGE)%NAME_PAR(1:8)=='Pressure') THEN
   IF (READ_2D(IMESSAGE)%NAME_PAR(11:12)=='Pa') READ_2D(IMESSAGE)%NAME_PAR(11:14)='hPa)'
   IF (READ_2D(IMESSAGE)%NAME_PAR(26:27)=='Pa') READ_2D(IMESSAGE)%NAME_PAR(26:29)='hPa)'
   IF (READ_2D(IMESSAGE)%NAME_PAR(19:20)=='Pa') READ_2D(IMESSAGE)%NAME_PAR(19:22)='hPa)'
   READ_2D(IMESSAGE)%FIELD(:,:)=READ_2D(IMESSAGE)%FIELD(:,:)*1.E-2
 ENDIF
 IF (READ_2D(IMESSAGE)%NAME_PAR(1:19)=='Geopotential height') THEN
   IF (READ_2D(IMESSAGE)%NAME_PAR(22:24)=='gpm') READ_2D(IMESSAGE)%NAME_PAR(22:25)='dam)'
   IF (READ_2D(IMESSAGE)%NAME_PAR(30:32)=='gpm') READ_2D(IMESSAGE)%NAME_PAR(30:33)='dam)'
   READ_2D(IMESSAGE)%FIELD(:,:)=READ_2D(IMESSAGE)%FIELD(:,:)*1.E-1
 ENDIF
 IF (READ_2D(IMESSAGE)%PARAM_DISCIPL==0.AND.READ_2D(IMESSAGE)%PARAM_CATEG==1.AND.&
READ_2D(IMESSAGE)%PARAM_IND==29) THEN ! Total snowfall
   READ_2D(IMESSAGE)%FIELD(:,:)=READ_2D(IMESSAGE)%FIELD(:,:)*1.E3
   IF (READ_2D(IMESSAGE)%NAME_PAR(17:18)=="m)") READ_2D(IMESSAGE)%NAME_PAR(17:23)="kg m-2)"
 ENDIF
 IF (READ_2D(IMESSAGE)%PARAM_DISCIPL==0.AND.READ_2D(IMESSAGE)%PARAM_CATEG==1.AND.&
READ_2D(IMESSAGE)%PARAM_IND==11) THEN ! Snow depth
   READ_2D(IMESSAGE)%FIELD(:,:)=READ_2D(IMESSAGE)%FIELD(:,:)*1.E2 ! m -> cm
   IF (READ_2D(IMESSAGE)%NAME_PAR(13:14)=="m)") READ_2D(IMESSAGE)%NAME_PAR(13:15)="cm)"
 ENDIF

  IF (CONTROL1=='y') THEN

! Automatic definition of graphic presentation parameters

    NPAGE=NPAGE+1
    PAGE(NPAGE)%NFIELD=1
    IFIELD=1 

    IF (.NOT.ALLOCATED(PAGE(NPAGE)%PRESENT_TYPE)) ALLOCATE(PAGE(NPAGE)%PRESENT_TYPE(PAGE(NPAGE)%NFIELD))
    IF (.NOT.ALLOCATED(PAGE(NPAGE)%NCONTOUR)) ALLOCATE(PAGE(NPAGE)%NCONTOUR(PAGE(NPAGE)%NFIELD))
    IF (.NOT.ALLOCATED(PAGE(NPAGE)%IPALETTE)) ALLOCATE(PAGE(NPAGE)%IPALETTE(PAGE(NPAGE)%NFIELD))
    IF (.NOT.ALLOCATED(PAGE(NPAGE)%ILABFORMAT)) ALLOCATE(PAGE(NPAGE)%ILABFORMAT(PAGE(NPAGE)%NFIELD))
    IF (.NOT.ALLOCATED(PAGE(NPAGE)%FIELD_TYPE)) ALLOCATE(PAGE(NPAGE)%FIELD_TYPE(PAGE(NPAGE)%NFIELD))
    IF (.NOT.ALLOCATED(PAGE(NPAGE)%FIELD)) ALLOCATE(PAGE(NPAGE)%FIELD(NX,NY,PAGE(NPAGE)%NFIELD))
    IF (.NOT.ALLOCATED(PAGE(NPAGE)%CONTOUR_VAL))&
 ALLOCATE(PAGE(NPAGE)%CONTOUR_VAL(100,PAGE(NPAGE)%NFIELD))

    PAGE(NPAGE)%PRESENT_TYPE(IFIELD)=1
    PAGE(NPAGE)%NCONTOUR(IFIELD)=21
    PAGE(NPAGE)%IPALETTE(IFIELD)=1
    PAGE(NPAGE)%FIELD_TYPE(IFIELD)='S'
    PAGE(NPAGE)%NX_CROSS = NX_CROSS
    PAGE(NPAGE)%LON_CROSS_INI = LON_CROSS_INI
    PAGE(NPAGE)%LON_CROSS_FIN = LON_CROSS_FIN
    PAGE(NPAGE)%LAT_CROSS_INI = LAT_CROSS_INI
    PAGE(NPAGE)%LAT_CROSS_FIN = LAT_CROSS_FIN

    PAGE(NPAGE)%FIELD(:,:,IFIELD)=READ_2D(IMESSAGE)%FIELD(:,:)

    ZMAX=MAXVAL(ABS(PAGE(NPAGE)%FIELD(:,:,IFIELD)))
    IORDER=INT(ALOG10(ZMAX))

    IF (IORDER<-4) PAGE(NPAGE)%ILABFORMAT(IFIELD)=8
    IF (IORDER==-3) PAGE(NPAGE)%ILABFORMAT(IFIELD)=7
    IF (IORDER==-2) PAGE(NPAGE)%ILABFORMAT(IFIELD)=6
    IF (IORDER==-1) PAGE(NPAGE)%ILABFORMAT(IFIELD)=3
    IF (IORDER==0) PAGE(NPAGE)%ILABFORMAT(IFIELD)=2
    IF (IORDER>=1) PAGE(NPAGE)%ILABFORMAT(IFIELD)=1
    IF (IORDER>4) PAGE(NPAGE)%ILABFORMAT(IFIELD)=5

    ZMAX=MAXVAL(PAGE(NPAGE)%FIELD(:,:,IFIELD))
    ZMIN=MINVAL(PAGE(NPAGE)%FIELD(:,:,IFIELD))
    RIS=(ZMAX+ZMIN)*0.5
    RIS=FLOAT(NINT(RIS*10.**(-IORDER)))*(10.**IORDER)

    DIS=(ZMAX-ZMIN)/FLOAT(PAGE(NPAGE)%NCONTOUR(IFIELD)-1)
    IORDER=INT(ALOG10(ABS(DIS)))
    DIS=FLOAT(NINT(DIS*10.**(-IORDER)))*(10.**IORDER)

    PAGE(NPAGE)%CONTOUR_VAL(1,IFIELD)=RIS-DIS*FLOAT(PAGE(NPAGE)%NCONTOUR(IFIELD)-1)*0.5
    DO ICONTOUR=2,PAGE(NPAGE)%NCONTOUR(IFIELD)
      PAGE(NPAGE)%CONTOUR_VAL(ICONTOUR,IFIELD)=PAGE(NPAGE)%CONTOUR_VAL(1,IFIELD)+DIS*FLOAT(ICONTOUR-1)
    ENDDO

    IFL_CONTOUR=1
    ISMOOTH=0

  ELSE

   IF (CONTROL2=='n') READ (INDEX2,GRAPH_DEF) ! Read graphic presentation parameters from namelist file

! Manual definition of graphic presentation parameters

    IF (CONTROL2=='y') THEN
      PRINT *,'Enter following parameters for current field: page index, &
 total field number and index of this field in the page: NPAGE, NFIELD, IFIELD'
      READ (*,*) NPAGE, PAGE(NPAGE)%NFIELD, IFIELD
    ELSE
      NPAGE=INDPAGE
      PAGE(NPAGE)%NFIELD=NFL
      IFIELD=INDFL
    ENDIF
    
    IF (.NOT.ALLOCATED(PAGE(NPAGE)%PRESENT_TYPE)) ALLOCATE(PAGE(NPAGE)%PRESENT_TYPE(PAGE(NPAGE)%NFIELD))
    IF (.NOT.ALLOCATED(PAGE(NPAGE)%NCONTOUR)) ALLOCATE(PAGE(NPAGE)%NCONTOUR(PAGE(NPAGE)%NFIELD))
    IF (.NOT.ALLOCATED(PAGE(NPAGE)%IPALETTE)) ALLOCATE(PAGE(NPAGE)%IPALETTE(PAGE(NPAGE)%NFIELD))
    IF (.NOT.ALLOCATED(PAGE(NPAGE)%ILABFORMAT)) ALLOCATE(PAGE(NPAGE)%ILABFORMAT(PAGE(NPAGE)%NFIELD))
    IF (.NOT.ALLOCATED(PAGE(NPAGE)%FIELD_TYPE)) ALLOCATE(PAGE(NPAGE)%FIELD_TYPE(PAGE(NPAGE)%NFIELD))
    IF (.NOT.ALLOCATED(PAGE(NPAGE)%FIELD)) ALLOCATE(PAGE(NPAGE)%FIELD(NX,NY,PAGE(NPAGE)%NFIELD))
    IF (.NOT.ALLOCATED(PAGE(NPAGE)%CONTOUR_VAL))&
 ALLOCATE(PAGE(NPAGE)%CONTOUR_VAL(100,PAGE(NPAGE)%NFIELD))

    PAGE(NPAGE)%NX_CROSS = NX_CROSS
    PAGE(NPAGE)%LON_CROSS_INI = LON_CROSS_INI
    PAGE(NPAGE)%LON_CROSS_FIN = LON_CROSS_FIN
    PAGE(NPAGE)%LAT_CROSS_INI = LAT_CROSS_INI
    PAGE(NPAGE)%LAT_CROSS_FIN = LAT_CROSS_FIN

    IF (CONTROL2=='y') THEN
      PRINT *,'Enter current field type (S-Scalar,U-u-component of wind,V-v-component of wind)'
      READ (*,*) PAGE(NPAGE)%FIELD_TYPE(IFIELD)
    ELSE
      PAGE(NPAGE)%FIELD_TYPE(IFIELD)=ATYPE
    ENDIF

    IF (CONTROL2=='y') THEN
      PRINT *
      PRINT *,'Enter following parameters for current field:'
      PRINT *,'   Representation type (1-color filling, 2-isolines, 3-symbols, 4-wind barb, 5-arrows, 6-color pixels,'
      PRINT *,'                        7-color isolines, 8-color writting of values in defined points,'
      PRINT *,'                        11-collor filling without isolines,'
      PRINT *,'                        12-collor filling for mask, 61-color pixels for mask'
      PRINT *,'                        41 - vector arrows with color filling for values of scalar product of vectors,'
      PRINT *,'                        51 - wind barbs with color filling for values of scalar product of vectors)'
      PRINT *,'   Color/isoline interval number (even number),'
      PRINT *,'   Color palette index (1-Rainbow, 2-Inverted Rainbow, &
 3-Orography, 4-Physical values, 5-Cloud cover, 6-Satellite radiance, &
 7-Satellite bright temperature, 11-White-black, 12-Black-wite, 13-DPC-term, &
 14-DPC-precip, 15-DPC-humid, 16-DPC-Ht0, 21-Soil type pixel, &
 22-Veg.type(GLC1990) type pixel, 23-Veg.type(GLC2000) type pixel),'
      PRINT *,'   Label format index (1-F9.0, 2-F9.1, 3-F9.2, 4-F9.3, 5-E9.0, 6-E9.1, 7-E9.2, 8-E9.3),'
      PRINT *,'   Flag for color/isoline interval definition (1-linear, 2-free),'
      PRINT *,'   Field smoothing degree (0, 1, 2...) :'
      PRINT *,'   GRAFTYPE, NCONT, COLORPAL, INDFORMAT, TYPECONTDEF, ISMOOTH'
      READ (*,*) PAGE(NPAGE)%PRESENT_TYPE(IFIELD), PAGE(NPAGE)%NCONTOUR(IFIELD),&
 PAGE(NPAGE)%IPALETTE(IFIELD), PAGE(NPAGE)%ILABFORMAT(IFIELD), IFL_CONTOUR, ISMOOTH
    ELSE
      PAGE(NPAGE)%PRESENT_TYPE(IFIELD)=GRAFTYPE
      PAGE(NPAGE)%NCONTOUR(IFIELD)=NCONT
      PAGE(NPAGE)%IPALETTE(IFIELD)=COLORPAL
      PAGE(NPAGE)%ILABFORMAT(IFIELD)=INDFORMAT
      IFL_CONTOUR=TYPECONTDEF
    ENDIF

    IF (FLAG_STATISTIC==0) THEN
      PAGE(NPAGE)%FIELD(:,:,IFIELD)=READ_2D(IMESSAGE)%FIELD(:,:)
    ELSE
    ! Statistic elaboration of data
      IF ((FLAG_CONST_FIELD==0.AND.IFILE==1).OR.(FLAG_CONST_FIELD==1.AND.IFILE<=2)) THEN
        PAGE(NPAGE)%FIELD(:,:,IFIELD)=READ_2D(IMESSAGE)%FIELD(:,:)
      ELSE
        IF (FLAG_STATISTIC_MODE==1) THEN
          ! Accumulation
          PAGE(NPAGE)%FIELD(:,:,IFIELD)=PAGE(NPAGE)%FIELD(:,:,IFIELD)+READ_2D(IMESSAGE)%FIELD(:,:)
        ELSEIF (FLAG_STATISTIC_MODE==2) THEN
          ! Average
          PAGE(NPAGE)%FIELD(:,:,IFIELD)=PAGE(NPAGE)%FIELD(:,:,IFIELD)+READ_2D(IMESSAGE)%FIELD(:,:)
          IF (IFILE==N_DATAFILE) THEN
            ZZZ=1./FLOAT(NFILE)
            PAGE(NPAGE)%FIELD(:,:,IFIELD)=PAGE(NPAGE)%FIELD(:,:,IFIELD)*ZZZ
          ENDIF
        ELSEIF (FLAG_STATISTIC_MODE==3) THEN
          ! Maximum value
          PAGE(NPAGE)%FIELD(:,:,IFIELD)=MAXVAL(PAGE(NPAGE)%FIELD(:,:,IFIELD)+READ_2D(IMESSAGE)%FIELD(:,:))
        ELSEIF (FLAG_STATISTIC_MODE==4) THEN
          ! Minimum value
          PAGE(NPAGE)%FIELD(:,:,IFIELD)=MINVAL(PAGE(NPAGE)%FIELD(:,:,IFIELD)+READ_2D(IMESSAGE)%FIELD(:,:))
        ENDIF
      ENDIF
    ENDIF

    IF (IFL_CONTOUR==1) THEN

      ZMAX=MAXVAL(READ_2D(IMESSAGE)%FIELD(:,:))
      ZMIN=MINVAL(READ_2D(IMESSAGE)%FIELD(:,:))

      IF (CONTROL2=='y') THEN
        PRINT *,'Enter color/isoline parameters: central value and interval value: &
 RIS, DIS, considering minimum and maximum values are: ',ZMIN,ZMAX
        READ (*,*) RIS, DIS
      ELSE
        RIS=CONTDEF(1)
        DIS=CONTDEF(2)
      ENDIF

!      IF (PAGE(NPAGE)%FIELD_TYPE(IFIELD)=='S') THEN
        PAGE(NPAGE)%CONTOUR_VAL(1,IFIELD)=RIS-DIS*FLOAT(PAGE(NPAGE)%NCONTOUR(IFIELD)-1)*0.5
        DO ICONTOUR=2,PAGE(NPAGE)%NCONTOUR(IFIELD)
          PAGE(NPAGE)%CONTOUR_VAL(ICONTOUR,IFIELD)=PAGE(NPAGE)%CONTOUR_VAL(1,IFIELD)+DIS*FLOAT(ICONTOUR-1)
        ENDDO
!      ELSE ! Vector field case
!        PAGE(NPAGE)%CONTOUR_VAL(1:PAGE(NPAGE)%NCONTOUR(IFIELD),IFIELD)=CONTDEF(1:PAGE(NPAGE)%NCONTOUR(IFIELD))
!      ENDIF

    ELSE ! IFL_CONTOUR==2 

!CONTDEF (:) contains all determined contour values,
! for GRAFTYPE = 3, CONTDEF(1:2) values interval (minimum and maximum) to present on a graphic;
! for GRAFTYPE = 4 or 5, CONTDEF(1:2) contains wind speed for arrow/bard symbol scale and minimum wind speed to present on a graphic
! CONTDEF(3) - symbol inedx, CONTDEF(4) - symbol size,
! CONTDEF(5:7) - color definition: red, green, blue (0-255);
! for GRAFTYPE = 7, CONTDEF(28:30) isolines color definition: red, green, blue (0-255);
! for GRAFTYPE = 12 (color filleng for mask), CONTDEF(1:2) values interval (minimum and maximum) to color filling,
! CONTDEF(3:5) - color definition: red, green, blue (proportion 0-1);
! for GRAFTYPE = 61 (color pixels for mask), CONTDEF(1) value for mask,
! CONTDEF(3:5) - color definition: red, green, blue (proportion 0-1).

      IF (CONTROL2=='y') THEN
        PRINT *,'Enter values for all color/isoline contours'
        READ (*,*) (PAGE(NPAGE)%CONTOUR_VAL(ICONTOUR,IFIELD),ICONTOUR=1,PAGE(NPAGE)%NCONTOUR(IFIELD))
      ELSE
        PAGE(NPAGE)%CONTOUR_VAL(1:PAGE(NPAGE)%NCONTOUR(IFIELD),IFIELD)=CONTDEF(1:PAGE(NPAGE)%NCONTOUR(IFIELD))
      ENDIF

    ENDIF

  ENDIF

  PRINT *,'For current field have been defined following presentation parameter:'
  PRINT *,'Page number: ',NPAGE,' Filed number in the page: ',IFIELD,' of ',PAGE(NPAGE)%NFIELD,' total numer fields'
  PRINT *,'Presentation mode ',PAGE(NPAGE)%PRESENT_TYPE(IFIELD),&
 ' Contour number ',PAGE(NPAGE)%NCONTOUR(IFIELD),&
 ' Color palette index',PAGE(NPAGE)%IPALETTE(IFIELD),&
 ' Label format index',PAGE(NPAGE)%ILABFORMAT(IFIELD),&
 ' Contour difinition flag ',IFL_CONTOUR,' Filed type ',PAGE(NPAGE)%FIELD_TYPE(IFIELD)
  PRINT *,'Defined contour values:'
  PRINT *,(PAGE(NPAGE)%CONTOUR_VAL(ICONTOUR,IFIELD),ICONTOUR=1,PAGE(NPAGE)%NCONTOUR(IFIELD))
  PRINT *,'Field smoothing degree: ',ISMOOTH

! Smoothing (filtr) of data field

  IF (ISMOOTH>0) CALL SMOOTH(PAGE(NPAGE)%FIELD(:,:,IFIELD),NX,NY,ZSMOOTH,ISMOOTH)

  IF (NPAGE_TOT<NPAGE) NPAGE_TOT=NPAGE

  IF (.NOT.ALLOCATED(PAGE(NPAGE)%NAME_PAR)) ALLOCATE(PAGE(NPAGE)%NAME_PAR(PAGE(NPAGE)%NFIELD))
  IF (.NOT.ALLOCATED(PAGE(NPAGE)%NAME_LEV)) ALLOCATE(PAGE(NPAGE)%NAME_LEV(PAGE(NPAGE)%NFIELD))
  IF (.NOT.ALLOCATED(PAGE(NPAGE)%NCHAR_NAME_PAR)) ALLOCATE(PAGE(NPAGE)%NCHAR_NAME_PAR(PAGE(NPAGE)%NFIELD))
  IF (.NOT.ALLOCATED(PAGE(NPAGE)%NCHAR_NAME_LEV)) ALLOCATE(PAGE(NPAGE)%NCHAR_NAME_LEV(PAGE(NPAGE)%NFIELD))
  IF (.NOT.ALLOCATED(PAGE(NPAGE)%DISCIPL)) ALLOCATE(PAGE(NPAGE)%DISCIPL(PAGE(NPAGE)%NFIELD))
  IF (.NOT.ALLOCATED(PAGE(NPAGE)%CATEG)) ALLOCATE(PAGE(NPAGE)%CATEG(PAGE(NPAGE)%NFIELD))
  IF (.NOT.ALLOCATED(PAGE(NPAGE)%INDEX)) ALLOCATE(PAGE(NPAGE)%INDEX(PAGE(NPAGE)%NFIELD))

  PAGE(NPAGE)%DISCIPL(IFIELD)=READ_2D(IMESSAGE)%PARAM_DISCIPL
  PAGE(NPAGE)%CATEG(IFIELD)=READ_2D(IMESSAGE)%PARAM_CATEG
  PAGE(NPAGE)%INDEX(IFIELD)=READ_2D(IMESSAGE)%PARAM_IND

! Titles definition

! Level description

  A60(1:60)=' '

   IF (IND_TEMPLATE_GRID /= 1000) THEN

     DO I=60,1,-1
       IF (READ_2D(IMESSAGE)%NAME_LEV(I:I)/=' ') GOTO 140
     ENDDO
 140  I1=I

    IF (IND_TEMPLATE_PROD < 30) THEN

! Meteorological and radar data

      IF (READ_2D(IMESSAGE)%LEV_TYPE==100) THEN
        READ_2D(IMESSAGE)%NAME_LEV(I1-2:I1+1)='hPa)'
        READ_2D(IMESSAGE)%LEV_TOP(1)=READ_2D(IMESSAGE)%LEV_TOP(1)/100
        IF (READ_2D(IMESSAGE)%LEV_BOT(1)/=INT(VALMISS)) READ_2D(IMESSAGE)%LEV_BOT(1)=READ_2D(IMESSAGE)%LEV_BOT(1)/100
        I1=I1+1
      ENDIF

     IF (READ_2D(IMESSAGE)%LEV_TOP(1) == 0.OR.READ_2D(IMESSAGE)%LEV_TOP(1) == INT(VALMISS)) THEN
       WRITE (A60,'(2A)') ' at ',READ_2D(IMESSAGE)%NAME_LEV(1:I1)
     ELSE
       IF (READ_2D(IMESSAGE)%LEV_BOT(1)==INT(VALMISS)) THEN
         IF (READ_2D(IMESSAGE)%LEV_TOP(2)==0) THEN
           WRITE (A60,'(A,I0,2A)') ' at ',READ_2D(IMESSAGE)%LEV_TOP(1),' ',READ_2D(IMESSAGE)%NAME_LEV(1:I1)
         ELSE
           WRITE (A60,'(A,I0,A,I0,2A)') ' at ',READ_2D(IMESSAGE)%LEV_TOP(1),'.e',-READ_2D(IMESSAGE)%LEV_TOP(2),&
 ' ',READ_2D(IMESSAGE)%NAME_LEV(1:I1)
          ENDIF
        ELSE
!          WRITE (A60,'(A,I0,A,I0,2A)') ' between ',&
! INT(READ_2D(IMESSAGE)%LEV_TOP(1)*0.1**READ_2D(IMESSAGE)%LEV_TOP(2)),&
! ' and ',&
! INT(READ_2D(IMESSAGE)%LEV_BOT(1)*0.1**READ_2D(IMESSAGE)%LEV_BOT(2)),&
! ' ',READ_2D(IMESSAGE)%NAME_LEV(1:I1)
          IF (READ_2D(IMESSAGE)%LEV_TOP(2)==0) THEN
            WRITE (A60,'(A,I0,A,I0,2A)') ' between ',&
   READ_2D(IMESSAGE)%LEV_TOP(1),' and ',READ_2D(IMESSAGE)%LEV_BOT(1),&
   ' ',READ_2D(IMESSAGE)%NAME_LEV(1:I1)
          ELSE
            WRITE (A60,'(A,I0,A,I0,A,I0,A,I0,2A)') ' between ',&
   READ_2D(IMESSAGE)%LEV_TOP(1),'.e',-READ_2D(IMESSAGE)%LEV_TOP(2),' and ',&
   READ_2D(IMESSAGE)%LEV_BOT(1),'.e',-READ_2D(IMESSAGE)%LEV_BOT(2),&
   ' ',READ_2D(IMESSAGE)%NAME_LEV(1:I1)
          ENDIF
        ENDIF
      ENDIF

    ELSE

! Satelite data

      CWN=FLOAT(READ_2D(IMESSAGE)%LEV_TOP(1))*0.1**READ_2D(IMESSAGE)%LEV_TOP(2)
      IF (CWN /= 0.) THEN
        CWL=1./CWN ! Central WaveLength (m) for the sensor channel
        CWL=NINT(CWL*1.E7)*1.E-1
      ELSE
        CWL=0.
      ENDIF
      WRITE (A60,'(A,F6.1,A)') ' in ',CWL,' micron channel'

    ENDIF

  ENDIF

  READ_2D(IMESSAGE)%NAME_LEV=A60

! Analysis of date and time

  CALL GRIB_GET(MESSAGE,'dataDate',IDATEINI)
  CALL GRIB_GET(MESSAGE,'dataTime',ITIMEINI)

  PAGE(NPAGE)%IDATE0(1)=IDATEINI/10000
  PAGE(NPAGE)%IDATE0(2)=(IDATEINI-PAGE(NPAGE)%IDATE0(1)*10000)/100
  PAGE(NPAGE)%IDATE0(3)=IDATEINI-PAGE(NPAGE)%IDATE0(1)*10000-PAGE(NPAGE)%IDATE0(2)*100
  PAGE(NPAGE)%IDATE0(4)=ITIMEINI/100
  PAGE(NPAGE)%IDATE0(5)=ITIMEINI-PAGE(NPAGE)%IDATE0(4)*100

! Forecast range

  CALL GRIB_GET(MESSAGE,'indicatorOfUnitOfTimeRange',PAGE(NPAGE)%IPERIOD(4))
! Time unit indicator: 0-minute, 1-hour, 2-day, 3-month, 4-year ... (Code Table 4.4)

  IF (PAGE(NPAGE)%IPERIOD(4)<=1) THEN
    CALL GRIB_GET(MESSAGE,'forecastTime',III)
    IF (IND_TEMPLATE_PROD==8) THEN
      CALL GRIB_GET(MESSAGE,'lengthOfTimeRange',LSTATISTIC)
      III=III+LSTATISTIC
    ENDIF
    IF (PAGE(NPAGE)%IPERIOD(4)==1) THEN 
      PAGE(NPAGE)%IPERIOD(1)=INT(FLOAT(III)/24.)
      PAGE(NPAGE)%IPERIOD(2)=INT(FLOAT(III)-FLOAT(PAGE(NPAGE)%IPERIOD(1))*24.)
      PAGE(NPAGE)%IPERIOD(3)=0
    ELSEIF (PAGE(NPAGE)%IPERIOD(4)==0) THEN
      PAGE(NPAGE)%IPERIOD(1)=INT(FLOAT(III)/24./60.)
      PAGE(NPAGE)%IPERIOD(2)=INT((FLOAT(III)-FLOAT(PAGE(NPAGE)%IPERIOD(1))*24.*60.)/60.)
      PAGE(NPAGE)%IPERIOD(3)=INT(FLOAT(III)-FLOAT(PAGE(NPAGE)%IPERIOD(1))*24.*60.-FLOAT(PAGE(NPAGE)%IPERIOD(2))*60.)
    ENDIF
    IF (IND_TEMPLATE_PROD==8) THEN
! Statistical field
      CALL GRIB_GET(MESSAGE,'typeOfStatisticalProcessing',ISTATISTIC)
! 0 - Average, 1 - Accumulation, 2 - Maximum, 3 - Minimum (code table 4.10)
      LSTATISTIC=LSTATISTIC*NFILE
      A120(1:120)=' '
      A121(1:120)=' '
      A120(1:60)=READ_2D(IMESSAGE)%NAME_PAR
      DO I=120,1,-1
        IF (A120(I:I)/=' ') GOTO 133
      ENDDO
 133  I1=I
      IF (PAGE(NPAGE)%IPERIOD(4)==1) THEN 
        IF (ISTATISTIC==0) WRITE (A121,'(2A,I3,A)') A120(1:I1),' average for ',LSTATISTIC,' h'
        IF (ISTATISTIC==1) WRITE (A121,'(2A,I3,A)') A120(1:I1),' accumulated in ',LSTATISTIC,' h'
        IF (ISTATISTIC==2) WRITE (A121,'(2A,I3,A)') A120(1:I1),' maximum for ',LSTATISTIC,' h'
        IF (ISTATISTIC==3) WRITE (A121,'(2A,I3,A)') A120(1:I1),' minimum for ',LSTATISTIC,' h'
      ELSEIF (PAGE(NPAGE)%IPERIOD(4)==0) THEN
        LSTATISTIC2=LSTATISTIC/60
        LSTATISTIC=LSTATISTIC-LSTATISTIC2*60
        IF (LSTATISTIC==0) THEN
          IF (ISTATISTIC==0) WRITE (A121,'(2A,I3,A)') A120(1:I1),' average for ',LSTATISTIC2,' h'
          IF (ISTATISTIC==1) WRITE (A121,'(2A,I3,A)') A120(1:I1),' accumulated in ',LSTATISTIC2,' h'
          IF (ISTATISTIC==2) WRITE (A121,'(2A,I3,A)') A120(1:I1),' maximum for ',LSTATISTIC2,' h'
          IF (ISTATISTIC==3) WRITE (A121,'(2A,I3,A)') A120(1:I1),' minimum for ',LSTATISTIC2,' h'
        ELSE
          IF (ISTATISTIC==0) WRITE (A121,'(2A,I3,A,I2,A)') A120(1:I1),' average for ',LSTATISTIC2,' h ',LSTATISTIC,' min'
          IF (ISTATISTIC==1) WRITE (A121,'(2A,I3,A,I2,A)') A120(1:I1),' accumulated in ',LSTATISTIC2,' h ',LSTATISTIC,' min'
          IF (ISTATISTIC==2) WRITE (A121,'(2A,I3,A,I2,A)') A120(1:I1),' maximum for ',LSTATISTIC2,' h ',LSTATISTIC,' min'
          IF (ISTATISTIC==3) WRITE (A121,'(2A,I3,A,I2,A)') A120(1:I1),' minimum for ',LSTATISTIC2,' h ',LSTATISTIC,' min'
        ENDIF
      ENDIF 
      A120=A121
      READ_2D(IMESSAGE)%NAME_PAR=A120(1:60)
    ENDIF 
  ELSE
    PRINT *,'Time range unit has not been provided by procedure, error in forecast term title can be verified'
  ENDIF

! Definition of week day for Initial Date and
! definition of Current Date (including week day) with
! using of linux system command (date)

  WRITE (SYSTEM_COMMAND,'(A25,1X,I4.4,7(1X,I2.2))') &
 './date_def_fortran.script',PAGE(NPAGE)%IDATE0(1), PAGE(NPAGE)%IDATE0(2), &
 PAGE(NPAGE)%IDATE0(3), PAGE(NPAGE)%IDATE0(4), PAGE(NPAGE)%IDATE0(5), &
 PAGE(NPAGE)%IPERIOD(1), PAGE(NPAGE)%IPERIOD(2), PAGE(NPAGE)%IPERIOD(3)
  CALL SYSTEM(SYSTEM_COMMAND)
  OPEN (11,FILE='graph_date_work.txt',STATUS='OLD')
  READ (11,*) PAGE(NPAGE)%IDATE0(6), PAGE(NPAGE)%IDATEC(1:6)
  CLOSE (11)

! Definition of field title 

  PAGE(NPAGE)%NAME_PAR(IFIELD)=READ_2D(IMESSAGE)%NAME_PAR
  PAGE(NPAGE)%NAME_LEV(IFIELD)=READ_2D(IMESSAGE)%NAME_LEV
  DO I=60,1,-1
    IF (PAGE(NPAGE)%NAME_PAR(IFIELD)(I:I)/=' ') GOTO 141
  ENDDO
 141 PAGE(NPAGE)%NCHAR_NAME_PAR(IFIELD)=I
  DO I=60,1,-1
    IF (PAGE(NPAGE)%NAME_LEV(IFIELD)(I:I)/=' ') GOTO 142
   ENDDO
 142 PAGE(NPAGE)%NCHAR_NAME_LEV(IFIELD)=I

! 'Manual' definition of page title

  IF (IFIELD==1) PAGE(NPAGE)%ACOMMENT(1)=TITLE_DEF

! Definition of comment for originig center
 
  IF (NPAGE==1) THEN
    IF (ICENTRE_CODE==80102) THEN
      IF (IMODEL_CODE==1) WRITE (PAGE(NPAGE)%ACOMMENT(5),'(A)') 'Bolam Model, CNR-ISAC, Italy'
      IF (IMODEL_CODE==2) WRITE (PAGE(NPAGE)%ACOMMENT(5),'(A)') 'Moloch Model, CNR-ISAC, Italy'
      IF (IMODEL_CODE==3) WRITE (PAGE(NPAGE)%ACOMMENT(5),'(A)') 'Globo Model, CNR-ISAC, Italy'
      IF (IMODEL_CODE==12) WRITE (PAGE(NPAGE)%ACOMMENT(5),'(A)') 'Bolam and Moloch models - Blended product, CNR-ISAC, Italy'
      IF (IND_TEMPLATE_PROD == 32) THEN ! RTTOV product
        IF (IMODEL_CODE==1) WRITE (PAGE(NPAGE)%ACOMMENT(5),'(A)') 'Bolam Model, CNR-ISAC, Italy; RTTOV Model, EUMETSAT'
        IF (IMODEL_CODE==2) WRITE (PAGE(NPAGE)%ACOMMENT(5),'(A)') 'Moloch Model, CNR-ISAC, Italy; RTTOV Model, EUMETSAT'
        IF (IMODEL_CODE==3) WRITE (PAGE(NPAGE)%ACOMMENT(5),'(A)') 'Globo Model, CNR-ISAC, Italy; RTTOV Model, EUMETSAT'
      ENDIF
    ELSE
!!!      WRITE (TABLE_NAME,'(2A)') DIR_PACKAGE(1:NCHAR),'/definitions/grib1/0.table'
      WRITE (TABLE_NAME,'(2A)') DIR_PACKAGE(1:NCHAR),'/definitions/common/c-1.table'
      OPEN (11,FILE=TABLE_NAME,STATUS='OLD',FORM='FORMATTED')
    15 LINE(1:100)=' '
      READ (11,'(A100)',END=115) LINE
      IF (LINE(1:1)/='#') THEN
        I1=0
        I2=0
        DO I=1,100
          IF (LINE(I:I)==' ') THEN
            IF (I1==0) THEN
              I1=I
            ELSE
              I2=I
              GOTO 150
            ENDIF
          ENDIF
        ENDDO
   150 READ (LINE(1:I1-1),*) III
        IF (III==ICENTRE_CODE) THEN
          READ (LINE(I2+1:100),'(A120)') PAGE(NPAGE)%ACOMMENT(5)
          GOTO 115
        ENDIF
      ENDIF
      GOTO 15
   115  CLOSE (11)
    ENDIF

  ENDIF

  GOTO 10

20 CALL GRIB_CLOSE_FILE(FILE_IND)

ENDDO DATAFILE_LOOP

NPAGE=NPAGE_TOT

! Definition pages title

PAGE(2:NPAGE)%ACOMMENT(5)=PAGE(1)%ACOMMENT(5) 

DO IPAGE=1,NPAGE

  IF (PAGE(IPAGE)%ACOMMENT(1)(1:1)/='#') CYCLE 

  DO I=1,120
    PAGE(IPAGE)%ACOMMENT(1)(I:I)=' '
  ENDDO

!PRINT *,'*****************'
!DO IFIELD=1,PAGE(IPAGE)%NFIELD
!  PRINT *,PAGE(IPAGE)%NAME_PAR(IFIELD)(1:PAGE(IPAGE)%NCHAR_NAME_PAR(IFIELD)),&
! PAGE(IPAGE)%NAME_LEV(IFIELD)(1:PAGE(IPAGE)%NCHAR_NAME_LEV(IFIELD))
!ENDDO
!PRINT *,'*****************'

  IF (PAGE(IPAGE)%NFIELD==1) THEN
    I1=PAGE(IPAGE)%NCHAR_NAME_PAR(1)
    I2=PAGE(IPAGE)%NCHAR_NAME_LEV(1)
    IF (PAGE(IPAGE)%NAME_LEV(1)(5:11)/='Unknown') THEN
      WRITE (PAGE(IPAGE)%ACOMMENT(1),'(2A)') PAGE(IPAGE)%NAME_PAR(1)(1:I1),PAGE(IPAGE)%NAME_LEV(1)(1:I2)
    ELSE
      WRITE (PAGE(IPAGE)%ACOMMENT(1),'(1A)') PAGE(IPAGE)%NAME_PAR(1)(1:I1)
    ENDIF
  ELSE
! Integrating of single fields title into united page title
    III=0
    DO IFIELD=1,PAGE(IPAGE)%NFIELD
      IF (PAGE(IPAGE)%FIELD_TYPE(IFIELD)=='U') III=IFIELD
    ENDDO
    NFIELD1=PAGE(IPAGE)%NFIELD
    IF (III>0) THEN
! There is wind field
      IFIELD=III
      PAGE(IPAGE)%NAME_PAR(IFIELD)='Wind'
      PAGE(IPAGE)%NCHAR_NAME_PAR(IFIELD)=4
      DO IFIELD=III+1,PAGE(IPAGE)%NFIELD-1
        PAGE(IPAGE)%NAME_PAR(IFIELD)=PAGE(IPAGE)%NAME_PAR(IFIELD+1)
        PAGE(IPAGE)%NCHAR_NAME_PAR(IFIELD)=PAGE(IPAGE)%NCHAR_NAME_PAR(IFIELD+1)
        PAGE(IPAGE)%NAME_LEV(IFIELD)=PAGE(IPAGE)%NAME_LEV(IFIELD+1)
        PAGE(IPAGE)%NCHAR_NAME_LEV(IFIELD)=PAGE(IPAGE)%NCHAR_NAME_LEV(IFIELD+1)
      ENDDO
      NFIELD1=PAGE(IPAGE)%NFIELD-1
    ENDIF
    IND(:)=0
    III=0
    DO IFIELD=1,NFIELD1
      DO I=1,NFIELD1
        IF (I/=IFIELD) THEN
          IF (PAGE(IPAGE)%NAME_LEV(I)==PAGE(IPAGE)%NAME_LEV(IFIELD)) THEN
            IF (III==0) THEN
              III=III+1
              IND(III)=IFIELD
              III=III+1
              IND(III)=I
            ELSE
              !IF (ALL(I/=IND(1:III-1))) THEN
              IF (ALL(I/=IND(1:III))) THEN
                III=III+1
                IND(III)=I
              ENDIF
            ENDIF
          ENDIF  
        ENDIF
      ENDDO
    ENDDO
    IF (III==0) THEN
! All parameter refer at diferent levels
      I=1
      DO IFIELD=1,NFIELD1
        I1=PAGE(IPAGE)%NCHAR_NAME_PAR(IFIELD)
        I2=PAGE(IPAGE)%NCHAR_NAME_LEV(IFIELD)
        WRITE (PAGE(IPAGE)%ACOMMENT(1)(I:I+I1),'(A)') PAGE(IPAGE)%NAME_PAR(IFIELD)(1:I1)
        IF (PAGE(IPAGE)%NAME_LEV(IFIELD)(5:11)/='Unknown') THEN
          WRITE (PAGE(IPAGE)%ACOMMENT(1)(I+I1+1:I+I1+I2),'(A)') PAGE(IPAGE)%NAME_LEV(IFIELD)(1:I2)
        ELSE
          I2=0
        ENDIF
        I=I+I1+I2+1
      ENDDO
    ELSE
! There some different parameters at the same level
      I=1
      DO IFIELD=1,NFIELD1-III+1
      IF (IFIELD/=IND(2).AND.IFIELD/=IND(3).AND.IFIELD/=IND(4)) THEN
        IF (IFIELD==IND(1)) THEN
          DO I3=1,III
            I1=PAGE(IPAGE)%NCHAR_NAME_PAR(IND(I3))
            WRITE (PAGE(IPAGE)%ACOMMENT(1)(I:I+I1),'(A)') PAGE(IPAGE)%NAME_PAR(IND(I3))(1:I1)
            I=I+I1+1
          ENDDO
          I2=PAGE(IPAGE)%NCHAR_NAME_LEV(IND(1))
          IF (PAGE(IPAGE)%NAME_LEV(IND(1))(5:11)/='Unknown') THEN
            WRITE (PAGE(IPAGE)%ACOMMENT(1)(I:I+I2),'(A)') PAGE(IPAGE)%NAME_LEV(IND(1))(1:I2)
          ELSE
            I2=0
          ENDIF
          I=I+I2+1
        ELSE
          I1=PAGE(IPAGE)%NCHAR_NAME_PAR(IFIELD)
          I2=PAGE(IPAGE)%NCHAR_NAME_LEV(IFIELD)
          WRITE (PAGE(IPAGE)%ACOMMENT(1)(I:I+I1),'(A)') PAGE(IPAGE)%NAME_PAR(IFIELD)(1:I1)
          IF (PAGE(IPAGE)%NAME_LEV(IFIELD)(5:11)/='Unknown') THEN
            WRITE (PAGE(IPAGE)%ACOMMENT(1)(I+I1+1:I+I1+I2),'(A)') PAGE(IPAGE)%NAME_LEV(IFIELD)(1:I2)
          ELSE
            I2=0
          ENDIF
        I=I+I1+I2+1
        ENDIF
      ENDIF
      ENDDO
    ENDIF
  ENDIF
ENDDO

PRINT *,NPAGE
PRINT *,NX,NY,X0,Y0,X00,Y00,DX,DY
IF (IND_TEMPLATE_GRID /= 1000) THEN
  PRINT *,MINVAL(LSM),MAXVAL(LSM)
  PRINT *,MINVAL(ALON),MAXVAL(ALON)
  PRINT *,MINVAL(ALAT),MAXVAL(ALAT)
ENDIF
PRINT *,' '
DO III=1,NPAGE
  PRINT *,III
  PRINT *,PAGE(III)%ACOMMENT(1) 
  PRINT *,PAGE(III)%IDATE0(1:6)
  PRINT *,PAGE(III)%IDATEC(1:6)
  PRINT *,PAGE(III)%IPERIOD(1:4)
  PRINT *,PAGE(III)%NFIELD
  DO I=1,PAGE(III)%NFIELD
    PRINT *,'  ',PAGE(III)%PRESENT_TYPE(I),PAGE(III)%NCONTOUR(I),PAGE(III)%IPALETTE(I),&
 PAGE(III)%ILABFORMAT(I),PAGE(III)%FIELD_TYPE(I)
    PRINT *,'  ',(PAGE(III)%CONTOUR_VAL(J,I),J=1,PAGE(III)%NCONTOUR(I))
    PRINT *,'  ',MINVAL(PAGE(III)%FIELD(:,:,I)),MAXVAL(PAGE(III)%FIELD(:,:,I))
    PRINT *,'  '
  ENDDO
  PRINT *,'  '
ENDDO
PRINT *,PAGE(1)%ACOMMENT(5) 

100 CONTINUE

CLOSE (100)
CLOSE (101)
CLOSE (200)
CLOSE (201)

IFLAG_GEO=1
IF (IND_TEMPLATE_GRID <= 3) THEN
  IF (INT(LSM(1,1))==INT(VALMISS)) THEN
    IFLAG_GEO=0
    IFL_LINE(1)=1
  ENDIF
ENDIF

IF (IND_TEMPLATE_GRID == 1000) THEN
  NX=PAGE(1)%NX_CROSS
  IF (VERT_COORD_TYPE == 100) THEN ! Pressure
    VERT_COORD(:)=ALOG(VERT_COORD(:))
  ENDIF
ENDIF

! Case of zoom in global area data

IF (XZOOM_INI>XZOOM_FIN) THEN
  XZOOM_INI=MAX(XZOOM_INI-0.5,0.)
  XZOOM_FIN=MIN(XZOOM_FIN+0.5,1.)
ENDIF

IF (NPAGE>0) CALL PLPLOT_MAP

STOP
END PROGRAM PLGRIBFA_MAIN_INPUT_GRIB2 
!------------------------------------------------------------------------------
SUBROUTINE SMOOTH(VT,NI,NJ,W,NSMOOTH)

!  Smoothing of a 2-D matrix VT(NI,NJ), based on a 5 point filter (laplacian).
!  Output in VOUT(NI,NJ) and also in VT (ATTENTION: VT is redefined)

!  Si suppone che i passi di griglia in x e y siano circa uguali!
!  Il grado di filtraggio e' prop. al param. W (tra 0 e 1)
!  (ma e' meglio chiamare piu'volte la routine che aumentare W oltre circa 0.5)
!  (il valore 0.5 e' il solo che conserva la media sull'area nel caso si
!  "saltino" dei punti nella griglia di uscita)

IMPLICIT NONE

INTEGER :: NI, NJ, NSMOOTH, NS, I, J
REAL :: W
REAL, DIMENSION(NI,NJ) :: VT, VOUT

 DO 1000 NS=1,NSMOOTH

   DO J=2,NJ-1
   DO I=2,NI-1
     VOUT(I,J)=(1.-W)*VT(I,J)+.25*W*(VT(I+1,J)+VT(I-1,J)+VT(I,J+1)+VT(I,J-1))
   ENDDO
   ENDDO

!  Punti ai bordi (filtro unidimensionale)

   DO I=2, NI-1
     VOUT(I,1)=(1.-.7*W)*VT(I,1)+.35*W*(VT(I+1,1)+VT(I-1,1))
     VOUT(I,NJ)=(1.-.7*W)*VT(I,NJ)+.35*W*(VT(I+1,NJ)+VT(I-1,NJ))
   ENDDO
   DO J=2, NJ-1
     VOUT(1,J)=(1.-.7*W)*VT(1,J)+.35*W*(VT(1,J+1)+VT(1,J-1))
     VOUT(NI,J)=(1.-.7*W)*VT(NI,J)+.35*W*(VT(NI,J+1)+VT(NI,J-1))
   ENDDO

!  Punti agli angoli della griglia (invariati)

   VOUT(1,1)=VT(1,1)
   VOUT(1,NJ)=VT(1,NJ)
   VOUT(NI,NJ)=VT(NI,NJ)
   VOUT(NI,1)=VT(NI,1)

   VT(:,:)=VOUT(:,:)

1000  CONTINUE

RETURN
END
!------------------------------------------------------------------------------
