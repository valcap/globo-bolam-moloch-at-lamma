! PROCEDURE READS WORK (CONVENTIONAL) FORMAT DATA AND CREATES GRAPHICAL 
! PRESENTATION FOR 2D FIELDS
!------------------------------------------------------------------------------
PROGRAM PLGRIBFA_MAIN_INPUT_WORK

USE MOD_GRAPH_PAR
USE MOD_GRAPH_PAR_DEF

IMPLICIT NONE

TYPE DATA
  INTEGER :: PARAM_DISCIPL, PARAM_CATEG, PARAM_IND, LEV_TYPE, LEV_TOP(2), LEV_BOT(2), IDATE(5), IPERIOD(4)
  REAL, DIMENSION(:,:), ALLOCATABLE :: FIELD
  CHARACTER :: NAME_PAR*60, NAME_LEV*60 
END TYPE DATA
TYPE (DATA), DIMENSION(1000) :: READ_2D

REAL, PARAMETER :: TZER=273.15
CHARACTER :: PPF_TYPE*3, DATAFILE_NAME*14, GRIDTYPE*20, SYSTEM_COMMAND*60,&
 CONTROL1*1, CONTROL2*1, CONTROL3*1, DIR_PACKAGE*60, TABLE_NAME*100, LINE*100, &
 A3*3, A10*10, A1*1, A60*60, A61*60, A120*120, A121*120
INTEGER :: N_DATAFILE, IFILE, FILE_IND=10, IMESSAGE, MESSAGE, ICENTRE_CODE=0, IVALMIS_R=255,&
 FLAG_CONST_FIELD, FLAG_STATISTIC, FLAG_STATISTIC_MODE, INDEX1, INDEX2, NFILE, &
 LEV_BOT_TYPE, III, IPAGE, NCHAR, IFIELD, IFL_CONTOUR, ICONTOUR, NVALUES, ISCAN, IEMOD,&
 I, J, IR, JR, IORDER, IDATEINI, ITIMEINI, ITABLE, ISTATISTIC, LSTATISTIC, LSTATISTIC2,&
 I1, I2, I3, NPAGE_TOT, NFIELD1, I_HALF
REAL :: X00_READ, RIS, DIS, ZMIN, ZMAX, ZZZ, ZSMOOTH=0.5
REAL, DIMENSION(:), ALLOCATABLE :: VALUE
REAL, DIMENSION(:,:), ALLOCATABLE :: FIELD_WORK
REAL*8 :: ZFAC, ZX0, ZY0, ZLON, ZLAT, ZZLAT, ZAARG, ZARG
INTEGER, DIMENSION(10) :: IND
INTEGER :: IDATE0(5), IPERIOD(3)

PI = ABS(ACOS(-1.))
X0=VALMISS
Y0=VALMISS
IMODEL_CODE=0
READ_2D(:)%PARAM_DISCIPL=INT(VALMISS)
READ_2D(:)%PARAM_CATEG=INT(VALMISS)
READ_2D(:)%PARAM_IND=INT(VALMISS)
READ_2D(:)%LEV_TYPE=INT(VALMISS)
READ_2D(:)%LEV_TOP(1)=INT(VALMISS)
READ_2D(:)%LEV_TOP(2)=INT(VALMISS)
READ_2D(:)%LEV_BOT(1)=INT(VALMISS)
READ_2D(:)%LEV_BOT(2)=INT(VALMISS)

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

  WRITE (DATAFILE_NAME,'(A,I3.3)') 'file_input_',IFILE
  OPEN (FILE_IND, FILE=DATAFILE_NAME, FORM='UNFORMATTED', STATUS='OLD')

! Loop on all the messages in a opened file

10  CONTINUE

! Read next messages group

! General parameters: not compulsory data
  READ (FILE_IND, END=20) ICENTRE_CODE, IMODEL_CODE ! Terminology of grib2 code:Code of Originating Centre and Code of generating Process (Model)

  IMESSAGE=IMESSAGE+1

! Grid parameters: compulsory data (!)
  READ (FILE_IND) NX, NY ! Grid point numbers 
  READ (FILE_IND) X0, Y0, X00_READ, Y00, DX, DY ! Center of coordinate rotation (X0, Y0), coordinate of south-west coner (X00, Y00), grib step (Dx, DY)

! Parameter index (terminology of grib2 code): not compulsory parameters
  READ (FILE_IND) READ_2D(IMESSAGE)%PARAM_DISCIPL, READ_2D(IMESSAGE)%PARAM_CATEG, READ_2D(IMESSAGE)%PARAM_IND

! Level parameters terminology of grib2 code): not compulsory data
  READ (FILE_IND) READ_2D(IMESSAGE)%LEV_TYPE, READ_2D(IMESSAGE)%LEV_TOP(1), READ_2D(IMESSAGE)%LEV_TOP(2) ! Level type, scaled level value, scale of level value 

! Data and forecast term: compulsory data (!)
  READ (FILE_IND) IDATE0(1:5), IPERIOD(1:3) ! Year, month, day, hour, minute of forecast start (IDATE0), day number, hour numbers, minute numbers of forecast term (IPERIOD)

! 2d-field data: compulsory data (!)
  ALLOCATE (READ_2D(IMESSAGE)%FIELD(NX,NY))
  READ (FILE_IND) READ_2D(IMESSAGE)%FIELD(1:NX,1:NY)

! End of reading of the messages group

  IF (X00_READ>180.) X00_READ=X00_READ-360.

  IF (IFILE==1.AND.IMESSAGE==1) THEN

    X00=X00_READ

    ALLOCATE(ALON(NX,NY))
    ALLOCATE(ALAT(NX,NY))
    ALLOCATE(LSM(NX,NY))
    ALON(:,:)=VALMISS
    ALAT(:,:)=VALMISS
    LSM(:,:)=VALMISS
    ALLOCATE(FIELD_WORK(NX,NY))  

    IF (ABS(X0)+ABS(Y0) > 0.01) THEN ! Rotated geography
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
    ELSE ! Non rotated geography
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

  ENDIF ! IMASSAGE==1.AND.IFILE==1

  TITLE_DEF=&
'########################################################################################################################'

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
  PRINT *,'PARAMETER    ',READ_2D(IMESSAGE)%NAME_PAR
  PRINT *,'LEVEL        ',READ_2D(IMESSAGE)%NAME_LEV(1:30),' ',&
 READ_2D(IMESSAGE)%LEV_TOP(1)*0.1**READ_2D(IMESSAGE)%LEV_TOP(2)
 IF (READ_2D(IMESSAGE)%LEV_BOT(1)/=INT(VALMISS))&
 PRINT *,'                         ',READ_2D(IMESSAGE)%LEV_BOT(1)*0.1**READ_2D(IMESSAGE)%LEV_BOT(2)

  PRINT *,'Do You want to create graphic presentation for this field? (y/n)'
  IF (CONTROL2=='y') THEN
    READ (*,*) A1
  ELSE
    READ (INDEX1,*) A1
  ENDIF

  IF (A1=='n'.AND.&
(READ_2D(IMESSAGE)%PARAM_DISCIPL*100000+READ_2D(IMESSAGE)%PARAM_CATEG*1000+READ_2D(IMESSAGE)%PARAM_IND)/=200000)&
  GOTO 10 ! Pass to next message

! Case of zoom in global area data

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

  IF (A1=='n') THEN
    LSM(:,:)=READ_2D(IMESSAGE)%FIELD(:,:)
    GOTO 10
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

   DO I=60,1,-1
     IF (READ_2D(IMESSAGE)%NAME_LEV(I:I)/=' ') GOTO 130
   ENDDO
 130  I1=I

   IF (READ_2D(IMESSAGE)%LEV_TYPE==100) THEN
     READ_2D(IMESSAGE)%NAME_LEV(I1-2:I1+1)='hPa)'
     READ_2D(IMESSAGE)%LEV_TOP(1)=READ_2D(IMESSAGE)%LEV_TOP(1)/100
     IF (READ_2D(IMESSAGE)%LEV_BOT(1)/=INT(VALMISS)) READ_2D(IMESSAGE)%LEV_BOT(1)=READ_2D(IMESSAGE)%LEV_BOT(1)/100
     I1=I1+1
   ENDIF

  A60(1:60)=' '
  IF (READ_2D(IMESSAGE)%LEV_TOP(1)==0) THEN
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
!      WRITE (A60,'(A,I0,A,I0,2A)') ' between ',&
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
  READ_2D(IMESSAGE)%NAME_LEV=A60

! Data, timem forecast range

  PAGE(NPAGE)%IDATE0(1:5)=IDATE0(1:5)
  PAGE(NPAGE)%IPERIOD(1:3)=IPERIOD(1:3)

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
    IF (PAGE(NPAGE)%NAME_PAR(IFIELD)(I:I)/=' ') GOTO 131
  ENDDO
 131 PAGE(NPAGE)%NCHAR_NAME_PAR(IFIELD)=I
  DO I=60,1,-1
    IF (PAGE(NPAGE)%NAME_LEV(IFIELD)(I:I)/=' ') GOTO 132
   ENDDO
 132 PAGE(NPAGE)%NCHAR_NAME_LEV(IFIELD)=I

! 'Manual' definition of page title

  IF (IFIELD==1) PAGE(NPAGE)%ACOMMENT(1)=TITLE_DEF

! Definition of comment for originig center
 
  IF (NPAGE==1) THEN

    IF (ICENTRE_CODE==80102) THEN
      IF (IMODEL_CODE==1) WRITE (PAGE(NPAGE)%ACOMMENT(5),'(A)') 'Bolam Model, CNR-ISAC, Italy'
      IF (IMODEL_CODE==2) WRITE (PAGE(NPAGE)%ACOMMENT(5),'(A)') 'Moloch Model, CNR-ISAC, Italy'
      IF (IMODEL_CODE==3) WRITE (PAGE(NPAGE)%ACOMMENT(5),'(A)') 'Globo Model, CNR-ISAC, Italy'
      IF (IMODEL_CODE==12) WRITE (PAGE(NPAGE)%ACOMMENT(5),'(A)') 'Bolam and Moloch models - Blended product, CNR-ISAC, Italy'
    ELSE
      WRITE (TABLE_NAME,'(2A)') DIR_PACKAGE(1:NCHAR),'/definitions/grib1/0.table'
      OPEN (11,FILE=TABLE_NAME,STATUS='OLD',FORM='FORMATTED')
    13 LINE(1:100)=' '
      READ (11,'(A100)',END=113) LINE
      IF (LINE(1:1)/='#') THEN
        I1=0
        I2=0
        DO I=1,100
          IF (LINE(I:I)==' ') THEN
            IF (I1==0) THEN
              I1=I
            ELSE
              I2=I
              GOTO 140
            ENDIF
          ENDIF
        ENDDO
   140 READ (LINE(1:I1-1),*) III
        IF (III==ICENTRE_CODE) THEN
          READ (LINE(I2+1:100),'(A120)') PAGE(NPAGE)%ACOMMENT(5)
          GOTO 113
        ENDIF
      ENDIF
      GOTO 13
   113  CLOSE (11)
    ENDIF

  ENDIF

  GOTO 10

20 CLOSE (FILE_IND)

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
PRINT *,MINVAL(LSM),MAXVAL(LSM)
PRINT *,MINVAL(ALON),MAXVAL(ALON)
PRINT *,MINVAL(ALAT),MAXVAL(ALAT)
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
IF (INT(LSM(1,1))==INT(VALMISS)) THEN
  IFLAG_GEO=0
  IFL_LINE(1)=1
ENDIF

! Case of zoom in global area data

IF (XZOOM_INI>XZOOM_FIN) THEN
  XZOOM_INI=MAX(XZOOM_INI-0.5,0.)
  XZOOM_FIN=MIN(XZOOM_FIN+0.5,1.)
ENDIF

IF (NPAGE>0) CALL PLPLOT_MAP

STOP
END PROGRAM PLGRIBFA_MAIN_INPUT_WORK
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
