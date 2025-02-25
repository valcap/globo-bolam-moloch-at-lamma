!               PREBOLAM

! Last update 16/10/2024

!----------------------------------------------------------------------------
! Per eliminare la vegetazione scommentare le righe di codice dopo "bare soil everywhere"
!----------------------------------------------------------------------------

! Sett. 2024: 
! Cambiamento in subroutine physiographic_param (def_soil.F90) e scrittura
! di model_param_constant.bin - piu' ordine, eliminazione dei passaggi 
! inutile. 
! Nuovo formato mhf (scrittura record non 1d ma 2d)

! Dic. 2021: Correzione di tg del mare profondo (confusione myocean, etc.)

! Gen. 2020: Nuovo formato mhf bolam (versione 20) in input
!
! Ago. 2018: 
! Nuovo schema del suolo (pochva), i nuovi dataset fisiografici (nuovo def_soil.F90),
! parametri del suolo sono diversi su vari livelli,
! manto nevoso muti-stato, il numero dei livelli di neve e' 11,
! la coordinate verticale nello manto nevoso e' kg/m/m,
! quando su un livello la neve non c'e', tutti i paramtri di neve sono -9999.;
! nuovo formato mhf: mfs_atm, mhf_soil e il file con i campi statici (model_param_constant.bin).
! Nov. 2015: cambiata completamente le definizione dei parametri di suolo e di vegetazione,
! della definizione della temperatura e del contenuto idrico del suolo,
! eliminata la matrice soilvegpar, sostituita la scrittura dei parametri di suolo e vegetazione.
! Ago. 2015: corretto calcolo p sui livelli nel caso dati ECMWF con liv. ibridi non consecutivi
! Mag. 2015: riverificate (e diminuite in valore assoluto) le correz. T e q del suolo GFS
! Apr. 2015: ridotte un po' meno le nubi alte GFS
! Nov. 2014: introd. filtro su sst di MyOcean.
! Ott. 2014: eliminata opzione per utilizzare dati ECMWF in formato grib1
! (i dati in grib1 devono essere prima convertiti in grib2 con grib_filter di grib_api).
! Viene automaticamente individuata il centro di provenienza dei dati grib2 (GFS o ECMWF).
! Se i dati ECMWF hanno coordinate ruotate e se il centro di rotazione è diverso da quello
! scritto in prebolam.inp, l'esecuzione di prebolam si ferma con una stampa di avviso.
! Correz. su check file (uso di inquire).
! Ago. 2014: ridotta correz. che secca il suolo GFS (zflatq) e corr. di T GFS (zdtg)
! Ritocchi sul modo di calc. t laghi e in sst_isac.
! Lug. 2014: correzioni nella def. di workf (vett. 2d ausiliario per def. di SST).
! Aum. peso attribuito a SST MyOcean (da .8 a .9).
! Giu. 2014: versione unificata che puo' operare su dati di input IFS, GFS,
! e BOLAM MHF (per self-nesting - per questo utilizzo sostituisce il precedente nesting.F90).
! Il file MHF di input puo' essere generato anche da GLOBO (sempre in coord. lat-lon regolari)
! ma con le stesse caratteristiche di un file MHF di BOLAM - per questo il nome "BOLAM" della
! variabile input_model comprende anche il caso di GLOBO.
! Puo' utilizzare due tipi di files mhf in input: un unico file che contiene istanti multipli
! ("united file") oppure una serie di file mhf ciascuno contenente un singolo istante
! ("separated files") - nel primo caso INPUT_FILE_STD = 1, nel secondo INPUT_FILE_STD = 2.
! Feb. 2014: modif def. umid. suolo GFS.
! Feb. 2014: modif. routine sst_isac per i nuovi dati myocean (area piu' grande).
! Dic. 2013: Ridotta corr. T bassi strati atmosf. in funz. corr. T primo liv. di suolo.
! Riviste corr. su Po Valley - introdotta correz. di q e dipendenza dalla copertura nuvolosa.
! Rivista riduz. di qg dati GFS (ora funz. di lat. e stagione - ma non valida per tutto il globo!).
! Ago. 2013: introdotta interp. solo bilineare nel caso di dati def. solo su cornici (frame).
! Introdotta matrice di flag mask_frame_out(nlon,nlat): in caso di def. solo su cornici,
! definisce se punti della griglia del modello in output derivano da un'iterpolazione corretta
! all'interno delle cornici della griglia di input (mask_frame_out=1), oppure no (mask_frame_out=0).
! (Notare che analogam. mask_frame(iana,jana) e' una matrice di flag per la griglia di input).
! In tal caso, i punti interni alla cornice sono def. = val_missing.
! Introd. smoothing in stratosfera anche per liv. ibridi IFS-ECMWF (riduce onde e anomalie PV su orog.).
! Lug. 2013:
! eliminate correzioni di T del suolo nei dati IFS-ECMWF
! (ma il contenuto idrico del suolo ECMWF sembra eccessivo per BOLAM in molte aree).
! Modif. correz. di T e Q nel suolo dati GFS
! Apr. 2013:
! Versione che permette di includere il polo (Oxana, apr. 2013). Modificata anche interp_spline_2d.
! Inserito utilizzo di fice e ice thickness (qui iceth) letta nei files GFS (icetk) (non dispon. nei dati IFS).
! Dove il ghiaccio marino si suppone "thick" (qui e nel modello per fice >= 0.8) si pone un iceth minimo
! di 0.5 m (nel caso GFS in cui icetk e' disponibile) oppure uno spessore tra 0.5 e 1.5 m nel caso ECMWF.
! Si scrivono in un solo campo fice e iceth (codificati) nel MHF.
! Introdotta ridistribuzione snow in funz. dell'orog. (routine in redistr_snow common_routines).
! Dati GFS: prelevare ICEC (fraz. di ghiaccio marino) e ICETK (spessore ghiaccio marino).
! Dati IFS-ECMWF: prelevare fraz. di ghiaccio marino (param. 031) e sea ice T liv. 1 (035).
! Nei dati ECMWF sembra sia sempre nulla la snow depth sul ghiaccio marino.
! Qui la snow e' massimizzata a 0.2 m di acqua equiv. (ca. 2 m di spessore) - ma notare che nei
! dati ECMWF originali puo' raggiungere i 5 m di acqua equiv. sulle Alpi e 15 m sulle terre circum-polari
! anche in estate (es. Groenlandia); evidentemente nello spessore sono compresi i ghiacciai.
! Attenz.: la scelta dei livelli di suolo da cui derivare la T dei laghi non e' generale:
! dipende dai livelli del modello di input.
! Genn. 2013: cambiato modo di interp. T e Q ai livelli del suolo.
! Se presenti, assimilano dati SST ISAC-MYOCEAN (Mediterraneo).
! Usa 7 livelli di suolo.
! Dic. 2012: corretto errore nelle def. di PSU e PSV ai boundaries est e sud, risp.
! Giu. 2012: introdotta correzione per ridurre arbitrariamente l'umidita' del suolo GFS
! in funz. della lat. e del livello e tra giugno e settembre (solo emisfero nord).
! Mag. 2012:
! Nei dati recenti ECMWF ci sono varie possibilità di def. la T del mare:
! SST (solo IFS, simile a T del primo livello), TSURF, T dei 4 livelli.
! Per ECMWF, suggeriscono di usare la T del primo liv. di "suolo",
! qui TGE(:,:,1) (uguale a SST sul mare, ma SST non e' definita dappertutto).
! Per GFS, le TGE sono definite solo sulla base di TSURFE (v. subr. conv_gfs).
! Sia per IFS che GFS, i valori delle variabili al suolo vanno prelevati in
! conv_ifs o conv_gfs (prima delle correzioni che si fanno in quelle subr.).
! Mag. 2011: fatta verifica di consistenza nel tempo del bias di T dei livelli di suolo GFS:
! da tutto aprile fino al 9/5/2011 tale bias era di circa tra -1 e -0.5 gradi
! rispetto alle stesse date del 2010 - quindi il bias negativo sembrava aumentato -
! ma dal 10/5 è aumentata la T  dei livelli di suolo, per cui
! rispetto al 2011 i valori apparivano (livelli da 1 a 4: +0.4, +1.5, +1.1, +0.5.
! Pertanto è stata diminuita la correz. di T di fattori opportuni per tener conto di tali
! variazioni (v. ridef. di ZDTG dopo il commento: "Reduction of corrections after sudden
! increase of GFS ground T of 9 May 2011").
! Ago. 2010: introdotto lapse rate vicino al suolo dipendente anche da quello del modello di input.
! Filtro campi GFS applicato solo in stratosfera.
! Correz. soil moisture: basata in parte su differenze tra ECMWF e GFS -
! una correz che dipende solo dalla lat. e non dalla long. non si puo' applicare
! completamente perche' sui deserti l'ECMWF ha valori di soil moist. minori, sulle foreste
! maggiori che GFS - quindi le correz. di soil moist. applicate qui non possono essere estese
! con confidenza oltre l'Europa.
! Riduzione di cloud water e ice sopra 450 hpa per dati GFS.
!==============================================================================
! Defining the input grid of BOLAM
! Version with MHF containing cloud water + ice and T in place of theta
! Caution: if input data have different resolution at different instants,
! matrices PHIGE and PHIGE2D must be deallocated (see !*!)
! To change the spacing of the hybrid vertical coordinates, see "s(k)"
!----------------------------------------------------------------------------
! Definition of initial and boundary condition data for the BOLAM model grid
! Input:
!  - analysis and/or forecast data of the ECMWF global model (IFS) (GRIB or GRIB2 format)
! (if mixed, conversion to GRIB2 must be performed before preprocessing)
! - or analysis and/or forecast data of the NOAA/NCEP global model (GFS) (GRIB2 format)
! - or bolam MHF file (for self-nesting of BOLAM)
! Output:
!  BOLAM MHF files.
! Uses GRIB_API library (ECMWF) and JASPER library needed by subroutine read_grib2_data.
!==============================================================================
! IFS-ECMWF grid data of atmospheric parameters can be defined on hybrid model levels or
! on constant pressure levels.
! The number of GFS data (typically 26 or 47) can vary also among different instants.
! The program can deal with input data defined only on a frame of the grid domain.
! Input grid can be defined on regular geographical coordinates or on rotated coordinates
! (caution: if input grid is on rotated coordinates, then the BOLAM grid will be rotated
! with the same rotation center - the possibility of having different rotation centers
! for input grid and for model grid is not provided).
! Input file names: grib_001 (normally the analysis time), grib_002 etc.
! Each file must contain both surface-soil and atmospheric level data
!==============================================================================
! ********  LIST OF REQUIRED VARIABLES IN INPUT DATA *********
!
!            %%%%     GFS-NOAA/NCEP    %%%%

! At atmospheric (isobaric) levels (analysis and forecast instants):
! temperature, u and v wind components, geopotential, rel. humidity,
! cloud total (liquid + ice) water content
! (the last 2 fields normally are available only at levels below 100 hPa).

! At surface (analysis instant only):
! sea-land fraction, geopotential, surface temperature, sea ice fraction, sea ice thickness,
! equivalent water content in the snow cover, temperature and specific volumetric
! water content at 4 levels in the soil.

!            %%%%      IFS-ECMWF       %%%%

! Case of data on atmospheric isobaric levels (analysis and forecast instants):
! temperature, u and v wind components, geopotential, specific humidity,
! cloud liquid content, cloud ice content.

! Case of data on atmospheric hybrid levels (ML) (analysis and forecast instants):
! temperature, u and v wind components, specific humidity,
! cloud liquid content, cloud ice content, natural logarithm of surface pressure
! (formally ML: model level), geopotential at the surface (formally ML),
! consistent with ps (analysis instant only).

! At surface (SFC):
! sea-land fraction, geopotential (SFC), soil type, sea ice fraction, sea ice temperature (analysis only),
! surf temperature, equivalent water content in the snow cover,
! temperature and specific volumetric water content at 4 levels in the soil
! (or at top and bottom soil levels only - analysis only or analysis and forecast instants).

!            %%%%    BOLAM MHF    %%%%

! The purpose is either nesting BOLAM into GLOBO or self-nesting of BOLAM.
! In both cases the format of the input model must be a BOLAM-type MHF.

!==============================================================================
module param

! ------  MODEL GRID PARAMETERS ------
! NOTE: values of parameters in this module are read in file prebolam.inp

! NLON, NLAT, NLEV: dimensions of the BOLAM grid
! NLEVG : number of soil levels used in BOLAM plus 1

integer :: nlon, nlat, nlev, nlevg, ntot
integer, parameter :: nlevg0=20, nlevsnow=11

! DLON, DLAT: resolution in degrees
! X0D, Y0D: coordinates in deg. of the centre of rotation (normally near the centre of the
! grid itself) (point "T" of the Arakawa grid). For non rotated grid, define X0D=Y0D=0.
! ALON0, ALAT0: coord. of the SW corner (point "v") of the model grid.
! In case of non rotated grid, ALON0, ALAT0 must define the true coordinates of the SW corner
! and DLON should be set to DLAT/COS(LAT0) where LAT0 is a mean latitude of the model grid.
! ALFA: exponent defining the hybrid vertical coordinate (ALFA=1. => normal sigma)
! ALFAmax < P0/(P0-PS_min)
! SLT: Soil layers thickness (m)

real :: dlon, dlat, x0d, y0d, alon0, alat0, p0, alfa

real, dimension(nlevg0) :: soil_lev0
real, dimension(:), allocatable :: soil_lev

! ------  INPUT DATA GRID PARAMETERS (read in prebolam.inp) ------

! INPUT_FORMAT: source of input model grid data: 'grib2' for IFS and GFS;
!               'mhfb' for Bolam and Globo (Globo-generated MHF file).
! INPUT_FILE_STD: valid for INPUT_FORMAT = 'mhfb' only (Bolam o Globo):
!       type of input data - 1: a single file including many instants;
!                            2: different files, one for each instant.
! NLEV_ATM_INP_MAX: max. possible no. of atmospheric levels in inputa data (used to allocate arrays only).
! NLEV_SOIL_INP_MAX: max. possible number of soil levels in input data (used to allocate arrays only).
! NIST: number of instants in input data (analysis + forecasts).
!       In case input data are from GLOBO or BOLAM it is the last instant number present in input
!       data to be processed.
! INST_START: number of the first instant for GLOBO or BOLAM MHF (mhfb format) input data
!       and INPUT_FILE_STD = 1 only.
! NSFC: flag (0 or 1) for analysis of surface parameters (soil, snow depth etc.) of input data:
!       if surface parameters are defined at the first instant only, set NSFC = 1, otherwise NSFC = 0.
! ----------------------------------------------------------------------------------------
! Input model grid variables set automatically in the program:
! INPUT_MODEL: name of the global model ('IFS' or 'GFS' or 'BOLAM',
!       the latter including GLOBO input fields) - the value of INPUT_MODEL is determined
!       depending on INPUT_FORMAT and grib2 decoding.
! IANA, JANA, DXA, DYA: dimensions and resol. of the input ("analysis") grid.
! NLEVP: actual no. of atmospheric levels (const. pressure or hybrid - read from input data).
! NLEVG_INP: actual number of soil levels in input data.
! X0A, Y0A: coordinates in deg. of the centre of rotation of input data grid (0 for no rotation).
! X1A, Y1A: coord. in deg. of the SW corner of the input data grid grid.
! INP_FLAG_ROT: flag indicating rotation of the input data grid (0: no rotation; 1: rotation).
! ----------------------------------------------------------------------------------------

  character (len=10) :: input_model, input_format
  integer :: iana0, iana, jana, nlev_atm_inp_max, nlevp, nlevp1, nlevg_inp, nlev_soil_inp_max, &
             nauxp, input_file_std, nist, inst_start, nsfc, inp_inp_flag_rot, flag_cut_paste, res_change
  real    :: x0a=0., y0a=0., dxa, dya, x1a, y1a
  logical :: surf_elaborate

! SLTE: Soil layers thickness (m) in input data

  real, dimension(:), allocatable :: slte

! NST is soil types number, NVT is vegetation type number

! Flags and other output parameters of the MHF

 integer, dimension(50) :: nfdr
 real, dimension(200)   :: pdr

 integer :: nst=31, nvt=22

 integer, dimension(5) :: idate0
 integer, dimension(4) :: iperiod_inp
 integer, dimension(3) :: iperiod

! P: Pressure of full levels; PL: logarithm of pressure; PHL: pressure of half-levels

 real, dimension(:,:,:), allocatable :: p, pl, phl

! Log. of press. on auxiliary levels and model hybrid levels

 real, dimension(:), allocatable :: plaux, s, plsig, plsig2

! Lapse rate

 real :: gamma, gammac, gamma_inp

! No. of days of each months

 integer, dimension(12) :: imon=(/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)

! Physical constants

 real :: pi, rd=287.05, rv=461.51, cp=1004.6, eps, ep, rcp, &
         g0=9.807, t0=273.15, e0=611.0, ra=6371.e+3

 namelist/param_prebolam/ input_format, input_file_std,                 &
                          nlon, nlat, nlev,                             &
                          dlon, dlat, x0d, y0d, alon0, alat0, p0, alfa, &
                          soil_lev0,                                    &
                          nlev_atm_inp_max, nlev_soil_inp_max,          &
                          nist, inst_start, nsfc

 real*4 :: val_missing=-9999.
 real :: xxctr, yyctr
 real, dimension(1) :: xxct, yyct

end module param
!----------------------------------------------------------------------------
module gribana

! Input fields and related parameters:

 integer :: npar3d=20, npar2d=50, npar3d_soil=10, level_type=0

 real, dimension(:,:,:,:), allocatable :: field3d, field3d_soil
 real, dimension(:,:,:), allocatable   :: field2d
 real, dimension(:,:,:), allocatable :: ge, te, ue, ve, qe, qce, qcwe, qcie, phle
 real, dimension(:,:), allocatable   :: phige, phige2d, psel, snowe, fmaske, tsurfe, qvsurfe, &
                                        soile, pmsle, cctote, u10e, v10e, t2e, td2e,         &
                                        ficee, ficeee, ticee, ticeee, icethe, icethee
 real, dimension(:,:), allocatable   :: t80e, q80e, p80e, u80e, v80e          ! variables at 80 m above surf.
 real, dimension(:,:,:), allocatable :: tge, qge
 real, dimension(:,:), allocatable :: alon_inp, alat_inp

 integer, dimension(:,:), allocatable :: mask_frame
 logical :: frame

! For 3D input data fields at model hybrid or sigma levels:
! AK and BK vectors are used to define pressure on IFS model levels

 real, dimension(:), allocatable :: ak, bk, ph_input, p_input

! Matrices with coordinates of the input grid

 real, dimension(:), allocatable :: xe, xue, ye, yve, lev_list, lev_list_soil, sige, sige2, siginte
 real :: alfae, p0e

! Flag iflag_cloude indicates presence (1) or absence (0) of cloud water/ice in input data

 integer :: iflag_cloude=0
 integer :: ifl_xf

! For conversion of volumetric water soil content into relative soil content:
! minimum and maximum soil water content as a function of soil type in input data

 real, dimension(:), allocatable   :: qgmine, qgmaxe
 real, dimension(:,:), allocatable :: qgmine2d, qgmaxe2d

end module gribana
!----------------------------------------------------------------------------
module mod_fields

! Output fields at surface and hybrid levels

 real, dimension(:,:,:), allocatable :: u, v, t, q, qc
 real, dimension(:,:), allocatable   :: phig, fmask, tsurf, tgsurf, ps, qvsurf, qgsurf, &
   water_table_depth, tg_bottom, qg_rel_bottom, qg_rel_surf_approx, &
   cloud, totpre, conpre, snfall, runoff,  snow, fice, cswfl, clwfl, chflux, &
   cqflux, t2min, t2max, albedo, emismap1, emismap2, rgm, rgq, tice, iceth, &
   soil_albedo_dry, soil_albedo_wet, soil_emiss1_dry, soil_emiss1_wet, soil_emiss2_dry, soil_emiss2_wet, &
   veg_lai, veg_frac, veg_root_depth, veg_roughness, veg_albedo, veg_emiss1, veg_emiss2, snow_dirt, fice_soil_surf  
 real, dimension(:,:,:), allocatable :: tg, qg, tg_first_guess, qg_rel_first_guess, fice_soil, &
 soil_map, veg_map, &
 soil_qmax, soil_qmin, soil_c, soil_rho, soil_psi, soil_k, soil_par_b, soil_par_c, soil_qrel_wilt, soil_qrel_ref, &
 snow_lev, snow_t, snow_fice, snow_age, snow_melt_age, snow_dens
 real :: veg_lai_max
 integer, dimension(:,:), allocatable :: ind_lev_soil_h_bottom, ind_lev_soil_w_bottom

! Matrices with grid coordinates, for grid rotation

 real, dimension(:,:), allocatable :: xxt, xxu, xxv, yyt, yyu, yyv, xxtg, yytg

! Intermediate 3D fields (interpolated horizontally at isobaric levels).
! (Wind components are computed at U and V grid points, on both U and V grids,
! before wind rotation - see comment in subr. rot_wind)

 real, dimension(:,:,:), allocatable :: gp, tp, upu, upv, vpu, vpv, upup, vpvp, qp, qcp

 integer, dimension(:,:), allocatable :: mask_frame_out

! Copies
 integer :: iflag_soil_snow_frc
 real, dimension(:,:), allocatable :: tsurf_copy, iceth_frc, fice_frc, &
 tsurf_frc, tgsurf_frc, qvsurf_frc, qgsurf_frc, fice_soil_surf_frc, snow_frc
 real, dimension(:,:,:), allocatable :: tg_copy, tg_frc, qg_frc, fice_soil_frc, &
 snow_lev_frc, snow_t_frc, snow_fice_frc, snow_age_frc, snow_melt_age_frc, snow_dens_frc

end module mod_fields
!----------------------------------------------------------------------------

program prebolam

! Assimilation of global model data on non rotated grid

use param
use gribana
use mod_fields

! Strings for input files names

 character (len=8), dimension(:), allocatable :: input_file
 character (len=20), dimension(:), allocatable :: input_file_atm, input_file_soil
 character (len=80) :: file_name_work, file_name_work_atm, file_name_work_soil
 logical :: lex, lex2

! Caution: length of FILE_NAME_WORK string must be 80 in according to declaration
! in subroutine READ_GRIB2_DATA

 real, dimension(:,:), allocatable :: work, workf, sstcle

! Intermediate (work) 2D and 3D (soil) fields interpolated on the model grid.
! They are derived (some only virtually) from the analyses,
! except HTOPI that is derived from the global orography dataset

 real, dimension(:,:), allocatable   :: phigi, phigi2d, psil, htopi, tsurfi, sstcli
 real, dimension(:,:,:), allocatable :: tgi, qgi, tgii, qgii

! Intermediate 2D fields at 80 m above surface (GFS), interpolated on the model grid.

 real, dimension(:,:), allocatable :: t80i, tv80i, q80i, pl80i, u80u, u80v, v80u, v80v, u80i, v80i

! Sea surface temperature (SST) derived from analysis

 real, dimension(:,:), allocatable :: sste, ssti

! Intermediate work fields

 real, dimension(:,:,:), allocatable  :: tv, wspeed
 real, dimension(:), allocatable      :: gaux
 real, dimension(:,:), allocatable    :: tvs, tsair, diftop, twater
 integer, dimension(3)                :: ijk

! Work vectors to accomodate the 80 m level data used in vertical interpolations

real, dimension(:), allocatable :: plaux_add80

! Standard deviation of orography

 real, dimension(:,:), allocatable :: zhtopvar

! Pressure at the surface (defined at T, U, V input grid point respectively)

 real, dimension(:,:), allocatable :: pst, psu, psv

! Work matrices used for vertical interpolation

 real, dimension(:), allocatable :: uk, vk, tk, qk, qck, upk, vpk, qpk, qcpk, tvk, pl2,    &
                                    tsoile, qsoile, zsoil, zsoil_work, tsoili, qsoili, qsoilav, tsoilav

 integer, dimension(:), allocatable :: iv, ivg
 integer :: icentre_inp, isubcentre_inp, imodel_inp, nk, nlsnow, flag_myocean=0, flag_constant_fields
 
 integer, parameter :: npoint=1
 real, dimension(npoint,1) :: lon_point=reshape((/37.62/),(/npoint,1/)), &
 lat_point=reshape((/55.83/),(/npoint,1/)), lon_rot_point, lat_rot_point
 
 real :: z1

!---------------------------------------------------------------------------

  frame = .false. ! it is assumed that full fields are available, unless possibly in case of IFS fc. data
  pi  = abs(acos(-1.))
  eps = rd/rv
  ep  = rv/rd-1.
  rcp = rd/cp

  open (11, file='prebolam.inp', status='old')
  read (11, param_prebolam)
  close (11)
  print *,'Parameters of prebolam (defined in prebolam.inp):'
  print param_prebolam

  ntot = nlon*nlat

  k=0
  do while (.true.)
    k=k+1
    if (int(soil_lev0(k)) == -9999) exit
  enddo
  nlevg=k-1
  allocate (soil_lev(nlevg))
  soil_lev(1:nlevg) = soil_lev0(1:nlevg)

  print *,'Depth (cm) of soil levels from the surface of the output soil layers:'
  do k = 1,nlevg
  print '(i3,2f8.1)', k, soil_lev(k)*1.e2
  enddo
  print*

  allocate(input_file(nist))
  allocate(input_file_atm(nist))
  allocate(input_file_soil(nist))

  if (input_format == 'grib2' .or. input_format == 'GRIB2') then
    input_file(:) = "grib_000"
  elseif (input_format == 'mhfb' .or. input_format == 'MHFB' .and. input_file_std == 2) then
    input_file_atm(:) = "mhfb_000_atm"
    input_file_soil(:) = "mhfb_000_soil"
  endif

  do i = 1,nist
    write (input_file(i)(6:8),'(i3.3)') i
    write (input_file_atm(i)(6:8),'(i3.3)') i
    write (input_file_soil(i)(6:8),'(i3.3)') i
  enddo

  if (input_format == 'mhfb' .or. input_format == 'MHFB' .and. input_file_std == 1) then   ! case of BOLAM/GLOBO MHF file
    input_file_atm(:) = "mhfb_atm"
    input_file_soil(:) = "mhfb_soil"
  endif

  ifile = 1
  file_name_work = input_file(ifile)
  if (input_format == 'mhfb' .or. input_format == 'MHFB') then
    file_name_work      = input_file_atm(ifile)
    file_name_work_atm  = input_file_atm(ifile)
    file_name_work_soil = input_file_soil(ifile)
  endif

  do while (.true.)
    inquire(file=file_name_work,exist=lex)
    if (lex) exit
    print *,"Input file ",trim(file_name_work)," not available"
#ifdef oper
    inquire(file='prebolam_stop.txt',exist=lex2)
    if (lex2) then
      print *,"Found file prebolam_stop.txt, demand to program stop"
      stop
    endif
    print *,"Waiting (sleep) 10 s"
    call system ("sleep 10")
#else
    stop
#endif
  enddo

!---------------------------------------------------------------------------
! Allocation of arrays

   allocate(ak     (nlev_atm_inp_max+1))
   allocate(bk     (nlev_atm_inp_max+1))
   allocate(ph_input(nlev_atm_inp_max+1))
   allocate(p_input (nlev_atm_inp_max+1))
   allocate(lev_list(nlev_atm_inp_max))
   allocate(lev_list_soil(nlev_soil_inp_max))
   allocate(sige   (nlev_atm_inp_max+1))
   allocate(sige2  (nlev_atm_inp_max))
   allocate(siginte(nlev_atm_inp_max))

!---------------------------------------------------------------------------

! Definitions of input data provider centre (input model),
! and parameters of input data grid

  if (input_format == 'grib2' .or. input_format == 'GRIB2') then
    call read_grib2_data(ifile,file_name_work,1,.true.,flag_cut_paste,icentre_inp,isubcentre_inp, imodel_inp, &
              iana0,jana,nlevp,nlev_atm_inp_max,nlevg_inp,nlev_soil_inp_max,                  &
              x0a,y0a,x1a,y1a,dxa,dya,idate0,iperiod_inp,lev_list,lev_list_soil,level_type,   &
              npar3d,npar2d,npar3d_soil,field3d,field2d,field3d_soil,ak,bk,val_missing)
    if (icentre_inp == 98) then
      input_model = 'IFS'
    elseif (icentre_inp == 7) then
      input_model = 'GFS'
    else
      print *
      print *,'Centre of input model data not identified:', icentre_inp
      stop
    endif
  elseif (input_format == 'mhfb' .or. input_format == 'MHFB') then ! BOLAM
    icentre_inp = 80                          ! identifies Rome RSMC, Italy
    isubcentre_inp = 102                      ! identifies CNR-ISAC, Italy
    input_model = 'BOLAM'
    imodel_inp = 1
    call read_bolam_mhf_data(0,file_name_work_atm,file_name_work_soil,input_file_std,.true.,                     &
              iana0,jana,nlevp,nlev_atm_inp_max,nlevg_inp,nlev_soil_inp_max,             &
              x0a,y0a,x1a,y1a,dxa,dya,idate0,iperiod,lev_list,lev_list_soil,level_type,  &
              npar3d,npar2d,npar3d_soil,field3d,field2d,field3d_soil,ak,bk,val_missing)
  else
    print *
    print *,'Input data format not identified:', input_format
    stop
  endif

! For the case of rotated input data grid:
! conversion of South Pole coordinates of rotation into
! coordinates of intersection point of zero latitude and zero longitude

  if (abs(x0a)>0.001.or.abs(y0a)>0.001) then
    print*
    inp_inp_flag_rot = 1  ! grid rotation
    print*, 'The input grid uses rotated coordinates,'
    print*, 'with rotation centre (lon, lat):', x0a, y0a
    print*
    if( (x0d > -180..and.x0d < 360..and.x0d /= x0a).or.  &
        (y0d > -90..and.y0d < 90..and.y0d /= y0a)) then
      print*, 'Error: the input grid and the BOLAM grid rotation centres must'
      print*, 'coincide!'
      print*, 'The BOLAM rotation centre in file prebolam.inp'
      print*, 'is:', x0d, y0d
      print*, 'Please correct the rotation centre in prebolam.inp;'
      print*, 'alternatively, re-download input grid data with a different rotation'
      print*, 'centre (or on a non-rotated grid), before restarting prebolam!'
      stop
    endif
  else
    inp_inp_flag_rot = 0  ! no grid rotation
  endif

!---------------------------------------------------------------------------
! Allocation of arrays

   allocate(p      (nlon,nlat,nlev_atm_inp_max))
   allocate(pl     (nlon,nlat,nlev_atm_inp_max))
   allocate(phl    (nlon,nlat,nlev_atm_inp_max+1))
   allocate(plaux  (2*nlev_atm_inp_max-1))
   allocate(s      (nlev))
   allocate(plsig  (nlev))
   allocate(plsig2 (nlev))
   allocate(slte   (nlev_soil_inp_max))

   allocate(xxt    (nlon,nlat))
   allocate(xxu    (nlon,nlat))
   allocate(xxv    (nlon,nlat))
   allocate(xxtg   (nlon,nlat))
   allocate(yyt    (nlon,nlat))
   allocate(yyu    (nlon,nlat))
   allocate(yyv    (nlon,nlat))
   allocate(yytg   (nlon,nlat))

   allocate(gp     (nlon,nlat,nlev_atm_inp_max))
   allocate(tp     (nlon,nlat,nlev_atm_inp_max))
   allocate(upu    (nlon,nlat,nlev_atm_inp_max))
   allocate(upv    (nlon,nlat,nlev_atm_inp_max))
   allocate(vpu    (nlon,nlat,nlev_atm_inp_max))
   allocate(vpv    (nlon,nlat,nlev_atm_inp_max))
   allocate(upup   (nlon,nlat,nlev_atm_inp_max))
   allocate(vpvp   (nlon,nlat,nlev_atm_inp_max))
   allocate(qp     (nlon,nlat,nlev_atm_inp_max))
   allocate(qcp    (nlon,nlat,nlev_atm_inp_max))

   allocate(u      (nlon,nlat,nlev))
   allocate(v      (nlon,nlat,nlev))
   allocate(t      (nlon,nlat,nlev))
   allocate(q      (nlon,nlat,nlev))
   allocate(qc     (nlon,nlat,nlev))

   allocate(phig           (nlon,nlat))
   allocate(fmask          (nlon,nlat))
   allocate(ps             (nlon,nlat))
   allocate(tsurf          (nlon,nlat))
   allocate(tgsurf         (nlon,nlat))
   allocate(tsurf_copy     (nlon,nlat))
   allocate(tsurf_frc      (nlon,nlat))
   allocate(tgsurf_frc     (nlon,nlat))
   allocate(qvsurf         (nlon,nlat))
   allocate(qgsurf         (nlon,nlat))
   allocate(qvsurf_frc     (nlon,nlat))
   allocate(qgsurf_frc     (nlon,nlat))
   allocate(cloud          (nlon,nlat))
   allocate(totpre         (nlon,nlat))
   allocate(conpre         (nlon,nlat))
   allocate(snfall         (nlon,nlat))
   allocate(runoff         (nlon,nlat))
   allocate(snow           (nlon,nlat))
   allocate(fice           (nlon,nlat))
   allocate(iceth          (nlon,nlat))
   allocate(snow_frc       (nlon,nlat))
   allocate(fice_frc       (nlon,nlat))
   allocate(iceth_frc      (nlon,nlat))
   allocate(cswfl          (nlon,nlat))
   allocate(clwfl          (nlon,nlat))
   allocate(chflux         (nlon,nlat))
   allocate(cqflux         (nlon,nlat))
   allocate(t2min          (nlon,nlat))
   allocate(t2max          (nlon,nlat))
   allocate(albedo         (nlon,nlat))
   allocate(emismap1       (nlon,nlat))
   allocate(emismap2       (nlon,nlat))
   allocate(fice_soil_surf (nlon,nlat))
   allocate(fice_soil_surf_frc (nlon,nlat))
   allocate(rgm            (nlon,nlat))
   allocate(rgq            (nlon,nlat))
   allocate(tice           (nlon,nlat))

   allocate(tg                (nlon,nlat,nlevg))
   allocate(qg                (nlon,nlat,nlevg))
   allocate(tg_first_guess    (nlon,nlat,nlevg))
   allocate(qg_rel_first_guess(nlon,nlat,nlevg))
   allocate(fice_soil         (nlon,nlat,nlevg))
   allocate(tg_copy           (nlon,nlat,nlevg))
   allocate(tg_frc            (nlon,nlat,nlevg))
   allocate(qg_frc            (nlon,nlat,nlevg))
   allocate(fice_soil_frc     (nlon,nlat,nlevg))

   allocate(snow_lev       (nlon,nlat,nlevsnow))
   allocate(snow_t         (nlon,nlat,nlevsnow))
   allocate(snow_fice      (nlon,nlat,nlevsnow))
   allocate(snow_age       (nlon,nlat,nlevsnow))
   allocate(snow_melt_age  (nlon,nlat,nlevsnow))
   allocate(snow_dens      (nlon,nlat,nlevsnow))
   allocate(snow_lev_frc   (nlon,nlat,nlevsnow))
   allocate(snow_t_frc     (nlon,nlat,nlevsnow))
   allocate(snow_fice_frc  (nlon,nlat,nlevsnow))
   allocate(snow_age_frc   (nlon,nlat,nlevsnow))
   allocate(snow_melt_age_frc(nlon,nlat,nlevsnow))
   allocate(snow_dens_frc  (nlon,nlat,nlevsnow))
   allocate(snow_dirt      (nlon,nlat))

   allocate(water_table_depth  (nlon,nlat))
   allocate(tg_bottom          (nlon,nlat))
   allocate(qg_rel_bottom      (nlon,nlat))
   allocate(qg_rel_surf_approx (nlon,nlat))

   allocate(ind_lev_soil_h_bottom(nlon,nlat))
   allocate(ind_lev_soil_w_bottom(nlon,nlat))

   allocate(soil_map       (nlon,nlat,nst+1))
   allocate(veg_map        (nlon,nlat,nvt+1))

   allocate(soil_qmax      (nlon,nlat,nlevg))
   allocate(soil_qmin      (nlon,nlat,nlevg))
   allocate(soil_c         (nlon,nlat,nlevg))
   allocate(soil_rho       (nlon,nlat,nlevg))
   allocate(soil_psi       (nlon,nlat,nlevg))
   allocate(soil_k         (nlon,nlat,nlevg))
   allocate(soil_par_b     (nlon,nlat,nlevg))
   allocate(soil_par_c     (nlon,nlat,nlevg))
   allocate(soil_qrel_wilt (nlon,nlat,nlevg))
   allocate(soil_qrel_ref  (nlon,nlat,nlevg))

   allocate(soil_albedo_dry(nlon,nlat))
   allocate(soil_albedo_wet(nlon,nlat))
   allocate(soil_emiss1_dry(nlon,nlat))
   allocate(soil_emiss1_wet(nlon,nlat))
   allocate(soil_emiss2_dry(nlon,nlat))
   allocate(soil_emiss2_wet(nlon,nlat))

   allocate(veg_lai        (nlon,nlat))
   allocate(veg_frac       (nlon,nlat))
   allocate(veg_root_depth (nlon,nlat))
   allocate(veg_roughness  (nlon,nlat))
   allocate(veg_albedo     (nlon,nlat))
   allocate(veg_emiss1     (nlon,nlat))
   allocate(veg_emiss2     (nlon,nlat))

   allocate(phigi  (nlon,nlat))
   allocate(phigi2d(nlon,nlat))
   allocate(psil   (nlon,nlat))
   allocate(htopi  (nlon,nlat))
   allocate(tsurfi (nlon,nlat))
   allocate(sstcli (nlon,nlat))
   allocate(tgi    (nlon,nlat,nlevg))
   allocate(qgi    (nlon,nlat,nlevg))
   allocate(tgii   (nlon,nlat,nlev_soil_inp_max))
   allocate(qgii   (nlon,nlat,nlev_soil_inp_max))
   allocate(ssti   (nlon,nlat))

   allocate(t80i (nlon,nlat))
   allocate(tv80i(nlon,nlat))
   allocate(q80i (nlon,nlat))
   allocate(pl80i(nlon,nlat))
   allocate(u80u (nlon,nlat))
   allocate(u80v (nlon,nlat))
   allocate(v80u (nlon,nlat))
   allocate(v80v (nlon,nlat))
   allocate(u80i (nlon,nlat))
   allocate(v80i (nlon,nlat))

   allocate(tv     (nlon,nlat,2*nlev_atm_inp_max-1))
   allocate(gaux   (2*nlev_atm_inp_max-1))

   allocate(plaux_add80(2*nlev_atm_inp_max))

   allocate(tvs    (nlon,nlat))
   allocate(tsair  (nlon,nlat))
   allocate(diftop (nlon,nlat))
   allocate(twater (nlon,nlat))
   allocate(wspeed (nlon,nlat,nlev))

   allocate(zhtopvar(nlon,nlat))

   allocate(pst    (nlon,nlat))
   allocate(psu    (nlon,nlat))
   allocate(psv    (nlon,nlat))

   allocate(mask_frame_out(nlon,nlat))

   allocate(uk     (nlev))
   allocate(vk     (nlev))
   allocate(tk     (nlev))
   allocate(qk     (nlev))
   allocate(qck    (nlev))
   allocate(upk    (nlev_atm_inp_max+1))
   allocate(vpk    (nlev_atm_inp_max+1))
   allocate(qpk    (nlev_atm_inp_max+1))
   allocate(qcpk   (nlev_atm_inp_max))
   allocate(tvk    (2*nlev_atm_inp_max-1))
   allocate(pl2    (nlev_atm_inp_max+1))

   allocate(tsoile (nlev_soil_inp_max+1))
   allocate(qsoile (nlev_soil_inp_max+1))
   allocate(zsoil (nlevg))
   allocate(zsoil_work (nlevg))
   allocate(tsoili (nlevg))
   allocate(qsoili (nlevg))
   allocate(tsoilav (nlevg))
   allocate(qsoilav (nlevg))

   allocate(iv     (nlev))
   allocate(ivg    (nlevg))

!---------------------------------------------------------------------------

   nfdr(:)  = 0
   pdr(:)   = 0.

! Definition of hybrid semi-integer levels of BOLAM.
! Levels are passed to the model in vector PDR, starting from index 40 onward
! S(K) = Sigma(k+1/2)

   do k = 1,nlev
     z = float(k)/float(nlev)
     if ( input_model /= 'BOLAM' ) then
      s(k) = .78*z+1.44*z**3-1.22*z**4  ! hybrid coordinates used for GFS or ECMWF input
     else
      s(k) = .75*z+1.62*z**3-1.37*z**4  ! hybrid coordinates in case of BOLAM self-nesting
     endif
      if((39+k).gt.200) then ! 200 is the present dimension of pdr (to be updated if necessary!)
        print*, "The number of levels of the output model exceeds the upper limit"
        print*, "At present it is 200-39 = 161; 200 is the dimension of array pdr"
        print*, "stop!"
        stop
      endif
     pdr(39+k) = s(k)
   enddo

   call livelli(log(1000.e2),s,nlev,plsig,p0,alfa) ! compute log of pressure at integer levels for 1000 hPa mslp
   call livelli(log(700.e2),s,nlev,plsig2,p0,alfa) ! compute log of pressure at integer levels for 700 hPa msl

! Print of semi-integer and integer levels

    print '(a)', " jklev    sig        sigint  p(sigint)  p(sigint)"
    print '(a)', "                             p0=1000hPa p0=700hPa"
    do jklev = 1, nlev
    if (jklev==1) then
      print '(i4,1x,2e12.4,2f8.2)', jklev, s(jklev), 0.5*s(jklev), exp(plsig(jklev))/100., exp(plsig2(jklev))/100.
    else
      print '(i4,1x,2e12.4,2f8.2)', jklev, s(jklev), 0.5*(s(jklev)+s(jklev-1)), exp(plsig(jklev))/100., exp(plsig2(jklev))/100.
    endif
    enddo

! Initialization of some parameters and flags of MHF

    pdr(1) = dlat
    pdr(2) = dlon
    pdr(4) = alat0
    pdr(5) = alon0

! Depth of soil levels (m)

    pdr(6:5+nlevg) = soil_lev(1:nlevg)

    pdr(36) = p0
    pdr(37) = alfa
    pdr(38) = y0d
    pdr(39) = x0d

    nfdr(1)  = 2
    nfdr(2)  =  nlon
    nfdr(3)  =  nlat
    nfdr(4)  =  nlev
    nfdr(15) =  nlevg

! Definition of rotated coordinate grid used in the model

 if (inp_inp_flag_rot==0) then

! Case of non rotated input (GFS or ECMWF) data grid: XXT, YYT, XXTG, YYTG contain the lon-lat coord. of
! the BOLAM rotated grid

    call rot_grid(x0d,y0d,alon0,         alat0+dlat*0.5,dlon,dlat,xxt,yyt,nlon,nlat)
    call rot_grid(x0d,y0d,alon0+dlon*0.5,alat0+dlat*0.5,dlon,dlat,xxu,yyu,nlon,nlat)
    call rot_grid(x0d,y0d,alon0,         alat0         ,dlon,dlat,xxv,yyv,nlon,nlat)
    xxtg = xxt
    yytg = yyt

    do j = 1,nlat
    do i = 1,nlon
     if(xxt(i,j) <= -180.) xxt(i,j) = xxt(i,j) + 360.
     if(xxu(i,j) <= -180.) xxu(i,j) = xxu(i,j) + 360.
     if(xxv(i,j) <= -180.) xxv(i,j) = xxv(i,j) + 360.
     if(xxt(i,j) >= 360.) xxt(i,j) = xxt(i,j) - 360.
     if(xxu(i,j) >= 360.) xxu(i,j) = xxu(i,j) - 360.
     if(xxv(i,j) >= 360.) xxv(i,j) = xxv(i,j) - 360.
    enddo
    enddo

 else

! Case of rotated (ECMWF) input data grid - In this case only XXTG and YYTG
! (used for ground parameter interpolation, given in regular lat-lon), that
! contain the lon-lat coord. of the BOLAM rotated grid, require rotation to be computed

    do i = 1,nlon
      xxt(i,:) = alon0+float(i-1)*dlon
      xxu(i,:) = xxt(i,:)+dlon*0.5
      xxv(i,:) = xxt(i,:)
    enddo
    do j = 1,nlat
      yyv(:,j) = alat0+float(j-1)*dlat
      yyt(:,j) = yyv(:,j)+dlat*0.5
      yyu(:,j) = yyt(:,j)
    enddo
    call rot_grid(x0d,y0d,alon0,alat0+dlat*0.5,dlon,dlat,xxtg,yytg,nlon,nlat)

 endif

! Case of input data on a global grid (typically for historical GFS data or ERA40 data...)

   if (inp_inp_flag_rot==0.and.x1a==0..and.x1a+dxa*float(iana0-1)>357.) then
     do jlat = 1,nlat
     do jlon = 1,nlon
       if (xxt(jlon,jlat) < 0.) xxt(jlon,jlat) = xxt(jlon,jlat)+360.
       if (xxu(jlon,jlat) < 0.) xxu(jlon,jlat) = xxu(jlon,jlat)+360.
       if (xxv(jlon,jlat) < 0.) xxv(jlon,jlat) = xxv(jlon,jlat)+360.
     enddo
     enddo
   endif

! Comput. of extreme coordinates of the rotated grid (for check)

   aminxt = minval(xxt(:,:))
   aminxu = minval(xxu(:,:))
   aminxv = minval(xxv(:,:))
   amaxxt = maxval(xxt(:,:))
   amaxxu = maxval(xxu(:,:))
   amaxxv = maxval(xxv(:,:))
   aminxtu= min(aminxt,aminxu)
   aminx  = min(aminxtu,aminxv)
   amaxxtu= max(amaxxt,amaxxu)
   amaxx  = max(amaxxtu,amaxxv)

   aminyt = minval(yyt(:,:))
   aminyu = minval(yyu(:,:))
   aminyv = minval(yyv(:,:))
   amaxyt = maxval(yyt(:,:))
   amaxyu = maxval(yyu(:,:))
   amaxyv = maxval(yyv(:,:))
   aminytu= min(aminyt,aminyu)
   aminy  = min(aminytu,aminyv)
   amaxytu= max(amaxyt,amaxyu)
   amaxy  = max(amaxytu,amaxyv)
   print*
   print*, 'Extremes of the selected output model grid:'
   print*, '    xmin           xmax           ymin           ymax'
   print*, aminx, amaxx, aminy, amaxy

   xxctr = alon0 +            0.5*dlon*float(nlon-1)
   yyctr = alat0 + 0.5*dlat + 0.5*dlat*float(nlat-1)
   call rot_grid(x0d,y0d,xxctr,yyctr,0.,0.,xxct,yyct,1,1)
   write(*,'(a,2f9.4)') ' Centre of the t-grid (lon., lat.):', xxct(1), yyct(1)
   xxctr = alon0 + 0.5*dlon + 0.5*dlon*float(nlon-1)
   yyctr = alat0 + 0.5*dlat + 0.5*dlat*float(nlat-1)
   call rot_grid(x0d,y0d,xxctr,yyctr,0.,0.,xxct,yyct,1,1)
   write(*,'(a,2f9.4)') ' Centre of the u-grid (lon., lat.):', xxct(1), yyct(1)
   xxctr = alon0 + 0.5*dlon*float(nlon-1)
   yyctr = alat0 + 0.5*dlat*float(nlat-1)
   call rot_grid(x0d,y0d,xxctr,yyctr,0.,0.,xxct,yyct,1,1)
   write(*,'(a,2f9.4)') ' Centre of the v-grid (lon., lat.):', xxct(1), yyct(1)

   iana = iana0

! Case when the model domain includes the whole longitude circle (typically for circumpolar domain location)
! One degree is a suitable interval to check if the model domain covers a full circumpolar area

   if (input_model /= 'BOLAM'.and.inp_inp_flag_rot == 0) then
     if ((abs(aminx) < 1. .and. abs(amaxx-360.) < 1.).or.(abs(aminx+180.) < 1. .and. abs(amaxx-180.) < 1.)) then
       iana = iana0 + 1
     endif
   endif
   flag_cut_paste = 0

! For control:
 call anti_rot_grid(x0d, y0d, lon_point(1:npoint,1), lat_point(1:npoint,1), &
 lon_rot_point(1:npoint,1), lat_rot_point(1:npoint,1), npoint, 1)
 ipoint=nint((lon_rot_point(1,1)-alon0)/dlon)
 jpoint=nint((lat_rot_point(1,1)-alat0-dlat*0.5)/dlat)
 print *,'Point ',ipoint,jpoint
! open (31,file="point_out.dat",status="unknown",form="formatted",position="append")

!---------------------------------------
! Start of loop for processing of input data files
!---------------------------------------
! Progressive number of input data file: IFILE
! Total number of input data files: NIST

 if (input_model == 'BOLAM' .and. input_file_std == 1) nist = nist-inst_start+1

do ifile = 1, nist

  if(input_model /= 'BOLAM') print *,'Input grib2 data read from ',input_file(ifile),' file'

! Reading of input data grid parameters only in grib2 files

  file_name_work = input_file(ifile)
  if (input_model == 'BOLAM') then
    file_name_work      = input_file_atm(ifile)
    file_name_work_atm  = input_file_atm(ifile)
    file_name_work_soil = input_file_soil(ifile)
  endif

  do while (.true.)
    inquire(file=file_name_work,exist=lex)
    if (lex) exit
    print *,"Input file ",trim(file_name_work)," not available"
#ifdef oper
    print *,"Waiting (sleep) 10 s"
    call system ("sleep 10")
#else
    stop
#endif
  enddo

! One more reading of input data grid parameters in case the grid changes
! for various analysis/forecast instants
! Note: if input data is IFS or GFS, then input grid parameters may change for different instants,
! if input data is Bolam or Globo, then input grid parameters are constant,
! and must be defined only once.

  if (input_format == 'grib2' .or. input_format == 'GRIB2') then
    call read_grib2_data(ifile,file_name_work,1,.true.,flag_cut_paste,icentre_inp,isubcentre_inp,imodel_inp, &
              iana0,jana,nlevp,nlev_atm_inp_max,nlevg_inp,nlev_soil_inp_max,                  &
              x0a,y0a,x1a,y1a,dxa,dya,idate0,iperiod_inp,lev_list,lev_list_soil,level_type,   &
              npar3d,npar2d,npar3d_soil,field3d,field2d,field3d_soil,ak,bk,val_missing)
  else                                        ! BOLAM
    if (ifile == 1)  &
      call read_bolam_mhf_data(ifile,file_name_work_atm,file_name_work_soil,input_file_std,.true.,              &
              iana0,jana,nlevp,nlev_atm_inp_max,nlevg_inp,nlev_soil_inp_max,            &
              x0a,y0a,x1a,y1a,dxa,dya,idate0,iperiod,lev_list,lev_list_soil,level_type, &
              npar3d,npar2d,npar3d_soil,field3d,field2d,field3d_soil,ak,bk,val_missing)
  endif

  if (x1a > 180.) x1a = x1a-360. ! Necessary, otherwise negative long. coordinates are mistaken

  iana = iana0
  if (input_model /= 'BOLAM'.and.inp_inp_flag_rot == 0) then
    if ((abs(aminx) < 1. .and. abs(amaxx-360.) < 1.).or.(abs(aminx+180.) < 1. .and. abs(amaxx-180.) < 1.)) then
      iana = iana0 + 1
    endif
  endif

!---------------------------------------------------------------------------
! ALLOCATION OF ARRAYS

   if (.not.allocated(field3d))      allocate(field3d(iana0,jana,nlev_atm_inp_max,npar3d))
   if (.not.allocated(field2d))      allocate(field2d(iana0,jana,npar2d))
   if (.not.allocated(field3d_soil)) allocate(field3d_soil(iana0,jana,nlev_soil_inp_max,npar3d_soil))
   if (.not.allocated(ge     )) allocate(ge  (iana,jana,nlev_atm_inp_max))
   if (.not.allocated(te     )) allocate(te  (iana,jana,nlev_atm_inp_max))
   if (.not.allocated(ue     )) allocate(ue  (iana,jana,nlev_atm_inp_max))
   if (.not.allocated(ve     )) allocate(ve  (iana,jana,nlev_atm_inp_max))
   if (.not.allocated(qe     )) allocate(qe  (iana,jana,nlev_atm_inp_max))
   if (.not.allocated(qce    )) allocate(qce (iana,jana,nlev_atm_inp_max))
   if (.not.allocated(qcwe   )) allocate(qcwe(iana,jana,nlev_atm_inp_max))
   if (.not.allocated(qcie   )) allocate(qcie(iana,jana,nlev_atm_inp_max))
   if (.not.allocated(phle   )) allocate(phle(iana,jana,nlev_atm_inp_max))

   if (.not.allocated(soile  ))  allocate(soile  (iana,jana))
   if (.not.allocated(phige  ))  allocate(phige  (iana,jana))
   if (.not.allocated(phige2d))  allocate(phige2d(iana,jana))
   if (.not.allocated(fmaske ))  allocate(fmaske (iana,jana))
   if (.not.allocated(qgmine2d)) allocate(qgmine2d(iana,jana))
   if (.not.allocated(qgmaxe2d)) allocate(qgmaxe2d(iana,jana))

   if (.not.allocated(psel   )) allocate(psel   (iana,jana))
   if (.not.allocated(tsurfe )) allocate(tsurfe (iana,jana))
   if (.not.allocated(qvsurfe)) allocate(qvsurfe (iana,jana))
   if (.not.allocated(pmsle  )) allocate(pmsle  (iana,jana))
   if (.not.allocated(snowe  )) allocate(snowe  (iana,jana))
   if (.not.allocated(ficee  )) allocate(ficee  (iana,jana))
   if (.not.allocated(ficeee )) allocate(ficeee (iana,jana))
   if (.not.allocated(icethe )) allocate(icethe (iana,jana))
   if (.not.allocated(icethee)) allocate(icethee(iana,jana))
   if (.not.allocated(ticee  )) allocate(ticee  (iana,jana))
   if (.not.allocated(ticeee )) allocate(ticeee (iana,jana))
   if (.not.allocated(xe     )) allocate(xe     (iana))
   if (.not.allocated(xue    )) allocate(xue    (iana))
   if (.not.allocated(ye     )) allocate(ye     (jana))
   if (.not.allocated(yve    )) allocate(yve    (jana))
   if (.not.allocated(alon_inp)) allocate(alon_inp(iana,jana))
   if (.not.allocated(alat_inp)) allocate(alat_inp(iana,jana))
   if (.not.allocated(mask_frame)) allocate(mask_frame(iana,jana))

   if (.not.allocated(tge    )) allocate(tge    (iana,jana,nlev_soil_inp_max))
   if (.not.allocated(qge    )) allocate(qge    (iana,jana,nlev_soil_inp_max))
   if (.not.allocated(cctote )) allocate(cctote (iana,jana))
   if (.not.allocated(u10e   )) allocate(u10e   (iana,jana))
   if (.not.allocated(v10e   )) allocate(v10e   (iana,jana))
   if (.not.allocated(t2e    )) allocate(t2e    (iana,jana))
   if (.not.allocated(td2e   )) allocate(td2e   (iana,jana))
   if (.not.allocated(sste   )) allocate(sste   (iana,jana))
   if (.not.allocated(work   )) allocate(work   (iana,jana))
   if (.not.allocated(workf  )) allocate(workf  (iana,jana))
   if (.not.allocated(sstcle )) allocate(sstcle (iana,jana))

   if (.not.allocated(t80e   )) allocate(t80e  (iana,jana))
   if (.not.allocated(q80e   )) allocate(q80e  (iana,jana))
   if (.not.allocated(p80e   )) allocate(p80e  (iana,jana))
   if (.not.allocated(u80e   )) allocate(u80e  (iana,jana))
   if (.not.allocated(v80e   )) allocate(v80e  (iana,jana))
!---------------------------------------------------------------------------

   flag_cut_paste = 0

! Case of input data on a global grid (typically for historical GFS data or ERA40 data...)
! In case of global area in input data, grid points having positive longitudes but west of
! Greenwich are "cut" and "pasted" (giving them negative longitudes) in order to reconstruct
! a grid crossing the Greenwich meridian

  if (iana0 == iana) then

    if (inp_inp_flag_rot==0.and.x1a==0..and.x1a+dxa*float(iana-1)>357.) then
     zx1a_ini = x1a
     if (aminx*amaxx>=0.) then
       if (x1a==0.) ifl_xf = 0
       if (x1a==-180.) ifl_xf = 1
     else
       if (x1a==0.) ifl_xf = 1
       if (x1a==-180.) ifl_xf = 0
     endif
     if (ifl_xf==1) then
       if (x1a==0.) then
         x1a = -180.
       elseif (x1a==-180.) then
         x1a = 0.
       endif

       flag_cut_paste = 1

       print*,'Shift of the input grid longitude:',' x1a=',x1a,'zx1a_ini=',zx1a_ini
     endif
    endif

    if (inp_inp_flag_rot==0) then
      if (x1a > 0..and.minval(xxt) <= 0.) then
        if (x1a-360.+float(iana-1)*dxa >= maxval(xxt)) x1a = x1a-360.
      elseif (x1a < 0..and.minval(xxt) >= 0.) then
        if (x1a+360. <= minval(xxt)) x1a = x1a+360.
      endif
    endif

  else

    flag_cut_paste = 0

  endif

  if(ifile.eq.1 .or. input_model /= 'BOLAM') then
   write (*,*) 'Input grid parameters:'
   write(*,'(a,i4,a,i4,a,f9.4,a,f9.4)') ' nx =',iana0,', ny =',jana,',  x0 =',x1a,',  y0 =',y1a
   write(*,'(a,f6.3,a,f6.3)') ' dx =',dxa,', dy =',dya
  endif

! Definition of geographical coordinates (in deg.) of the input data grid (rotated or non rotated)

   if (input_model == "IFS" .or. input_model == "GFS") then
     call rot_grid(x0a,y0a,x1a,y1a,dxa,dya,alon_inp,alat_inp,iana,jana)
   else if (input_model == "BOLAM") then
     call rot_grid(x0a,y0a,x1a,y1a+dya*0.5,dxa,dya,alon_inp,alat_inp,iana,jana)
   endif

! Check that the rotated model grid is inside the input data grid

   do i = 1,iana
    xe(i)  = x1a + (i-1)*dxa  ! x-coord. of T and v-points of the staggered BOLAM-type grid (also all points of the IFS and GFS grids)
    xue(i) = xe(i) + dxa*0.5  ! x-coord. of u-points of the staggered BOLAM-type grid
   enddo
   xemin = minval(xe(:))
   xemax = maxval(xe(:))
   if (input_model == "BOLAM") xemax = maxval(xue(:))

   if(aminx<xemin .or. amaxx>xemax) then
     print*
     print*,'Error: output model grid not contained in the input grid:'
     print*,'xmin  =', aminx, 'xmax  =', amaxx
     print*,'xemin =', xemin, 'xemax =', xemax
      if (input_model == "IFS" .or. input_model == "GFS" .and. inp_inp_flag_rot == 0) then
       print*,'Note: if the input model data (GFS or IFS-ECMWF, non rotated) contain'
       print*,'an entire latitude circle, the prescribed longitude interval in the'
       print*,'retrieval script must be (0, 360) or (-180, 180)'
      endif
    stop ' Aborting.'
    endif

    do j = 1,jana
      if (input_model == "IFS" .or. input_model == "GFS") then
       ye(j) = y1a + (j-1)*dya      ! y-coord. of all points of the IFS and GFS grids
      else if (input_model == "BOLAM") then
       yve(j) = y1a + (j-1)*dya     ! y-coord. of v-points of the staggered BOLAM-type grid
       ye(j)  = yve(j) + dya*0.5    ! y-coord. of T and u-points of the staggered BOLAM-type grid
      endif
    enddo
    yemin = minval(ye(:))
    yemax = maxval(ye(:))
    if (input_model == "BOLAM") yemin = minval(yve(:))

    if(aminy<yemin .or. amaxy>yemax) then
      print*
      print*,'Error: output model grid not contained in the analysis grid'
      print*,'ymin=',aminy,'ymax=',amaxy,'yemin=',yemin,'yemax=',yemax
      stop 'Aborting.'
    endif

! Call of the procedure for reading and decoding input data files in grib2 format
! or in BOLAM MHF (mhfb) format

  file_name_work = input_file(ifile)

  if (input_format == 'grib2' .or. input_format == 'GRIB2') then
    call read_grib2_data(ifile,file_name_work,1,.false.,flag_cut_paste,icentre_inp,isubcentre_inp,imodel_inp, &
             iana0,jana,nlevp,nlev_atm_inp_max,nlevg_inp,nlev_soil_inp_max,                    &
             x0a,y0a,x1a,y1a,dxa,dya,idate0,iperiod_inp,lev_list,lev_list_soil,level_type,     &
             npar3d,npar2d,npar3d_soil,field3d,field2d,field3d_soil,ak,bk,val_missing)
  else                   ! BOLAM
    if (input_file_std == 1 .and. ifile == 1) then ! skip of initial instants in the united input file
      do ifile2 = 1, inst_start-1
          print *,"Start reading and skipping instant no.", ifile2
          call read_bolam_mhf_data(ifile,file_name_work_atm,file_name_work_soil,input_file_std,.false.,         &
              iana0,jana,nlevp,nlev_atm_inp_max,nlevg_inp,nlev_soil_inp_max,            &
              x0a,y0a,x1a,y1a,dxa,dya,idate0,iperiod,lev_list,lev_list_soil,level_type, &
              npar3d,npar2d,npar3d_soil,field3d,field2d,field3d_soil,ak,bk,val_missing)
      enddo
    endif
    print*, "Start reading and processing instant no.", ifile + inst_start -1
    call read_bolam_mhf_data(ifile,file_name_work_atm,file_name_work_soil,input_file_std,.false.,               &
              iana0,jana,nlevp,nlev_atm_inp_max,nlevg_inp,nlev_soil_inp_max,            &
              x0a,y0a,x1a,y1a,dxa,dya,idate0,iperiod,lev_list,lev_list_soil,level_type, &
              npar3d,npar2d,npar3d_soil,field3d,field2d,field3d_soil,ak,bk,val_missing)
  endif

  nlevp1 = nlevp+1
  nauxp  = 2*nlevp-1

  surf_elaborate=.false.
  if (ifile == 1.or.nsfc == 0) surf_elaborate=.true.

! Preliminary elaboration of read and decoded input data

  if (input_model == 'IFS') then
    call conv_ifs_data(ifile)
  elseif (input_model == 'GFS') then
    call conv_gfs_data(ifile)
  elseif (input_model == 'BOLAM') then
    call conv_bolam_data(ifile)
  endif

  nfdr(5:9)   = idate0(1:5)
  nfdr(10:12) = iperiod(1:3)

!  Definition of iflag_80m: defined =1 in case of the presence of 80 m (above surf.) data (GFS only)

    iflag_80m = 0
    if (input_model=="GFS") then
      if (t80e(1,1)>val_missing.and.q80e(1,1)>val_missing.and.p80e(1,1)>val_missing &
        .and.u80e(1,1)>val_missing.and.v80e(1,1)>val_missing) then
        print*
        print*, "Variables at the level of 80 m above surface are present in GFS data and used."
        iflag_80m = 1
      else
        print*, "Variables at the level of 80 m above surface are not present in GFS data."
        print*
      endif
    endif

! iflag_80m = 0 ! for test only: if uncommented, data at the 80 m level above ground are not used

  if (ifile == 1) then ! At the first instant only

    print *
    if (level_type == 1) then
       print *,'Type of atmospheric levels: isobaric'
    elseif (level_type == 2) then
       print *,'Type of atmospheric levels: hybrid'
    else
       print *,'Type of atmospheric levels not defined'
    endif
    print *

! Definition of day in the year (neglecting leap-year because it
! is used only to approximate seasonal variation) at the first instant

    inyar = nfdr(5)
    inmon = nfdr(6)
    inday = nfdr(7)

    ndayr = inday
    do jmon = 1,inmon-1
      ndayr = ndayr + imon(jmon)
    enddo
    zday = float(ndayr-1)+float(nfdr(10))+float(nfdr(11))/24. ! day of the year as real (with hour precision)

! Definition of climatological lapse rate depending on season (in both hemispheres),
! varying from 4.5E-3 at the end of Dec. to 7.3E-3 at the end of June in the N.H.
! (opposite in the S.H. - no seasonal change if the domain centre is in the tropics).

    if(y0d > 30.) then
    zcoeday = 0.5*(1.-cos(zday*2.*pi/365.))   ! N.H.
    elseif(y0d < -30.) then
    zcoeday = 0.5*(1.+cos(zday*2.*pi/365.))   ! S.H.
    else
    zcoeday = 0.5                             ! tropics
    endif
    gammac = 4.5e-3 + 2.8e-3*zcoeday

    write(*, '(a,f8.4)') " Climatological lapse rate:", gammac

  endif ! end of definitions at the first instant only (ifile=1)

! call outgraph(7,1,iana,jana,0.,0.,x1a,y1a,dxa,dya,&
! 2, 0,  0,  1, 1, 0,idate0(1:5),iperiod(1:3),fmaske(:,:),-1.,1.)
! call outgraph(7,1,iana,jana,0.,0.,x1a,y1a,dxa,dya,&
! 2, 0,  9,106, 5, 0,idate0(1:5),iperiod(1:3),qge(:,:,1),1.,0.)

! --------------------------------------------------
! Definition of mask frame for output (model) grid

  mask_frame_out(:,:) = 1 ! 1 means a valid value, but in case frame=false, mask_frame_out is not used

  if (frame) then

   j = jana/2
   do i = 1,iana/2
    if (mask_frame(i,j) == 0) then
    ifw = i-1
    exit
    endif
   enddo

   do i = iana/2,iana
    if (mask_frame(i,j) == 1) then
    ife = i
    exit
    endif
   enddo

   i = iana/2
   do j = 1,jana/2
    if (mask_frame(i,j) == 0) then
    jfs = j-1
    exit
    endif
   enddo

   do j = jana/2,jana
    if (mask_frame(i,j) == 1) then
    jfn = j
    exit
    endif
   enddo

   zlonfw = x1a+float(ifw-1)*dxa
   zlonfe = x1a+float(ife-1)*dxa
   zlatfs = y1a+float(jfs-1)*dya
   zlatfn = y1a+float(jfn-1)*dya

   ifw = int((zlonfw-alon0)/dlon)+1
   ife = nlon-int((alon0+float(nlon-1)*dlon+dlon*0.5-zlonfe)/dlon)
   jfs = int((zlatfs-alat0)/dlat)+1
   jfn = nlat-int((alat0+float(nlat-1)*dlat+dlat*0.5-zlatfn)/dlat)

   mask_frame_out(:,:) = 1
   mask_frame_out(ifw+1:ife-1,jfs+1:jfn-1) = 0

  endif

  if (surf_elaborate) then      ! very long if!

! Over land T is prescribed as a deep soil T (at about 50 cm)
! to be used below to define lake temperature.

!   call plotout(tsurfe,iana,jana,99)
!   call plotout(tge(:,:,1),iana,jana,99)
!   stop

  if (input_model /= 'BOLAM') then   ! GFS or ECMWF
    do j = 1,jana
    do i = 1,iana
      if (fmaske(i,j) >= .5) then
        if (int(tsurfe(1,1)) /= int(val_missing) .and. input_model == 'GFS') then  ! for GFS, SST can be defined only if TSURFE is available
          workf(i,j) = tsurfe(i,j)
        else
          workf(i,j) = tge(i,j,1)
        endif
      else
        workf(i,j) = 0.5*(tge(i,j,nlevg_inp/2)+tge(i,j,nlevg_inp/2+1)) ! guess for T of lakes
      endif
    enddo
    enddo
  else                             ! BOLAM
    do j = 1,jana
    do i = 1,iana
      if (fmaske(i,j) >= .5) then
        workf(i,j) = tge(i,j,1)
      else
        workf(i,j) = 0.5*(tge(i,j,nlevg_inp/2+1)+tge(i,j,nlevg_inp/2+2)) ! guess for T of lakes
      endif
    enddo
    enddo
  endif

! Seatemp extends T from the sea towards the land

  call seatemp(workf,fmaske,sste,iana,jana,3,1,0.9)
  call seatemp(tge(:,:,nlevg_inp),fmaske,sstcle,iana,jana,4,1,0.9) ! constant "deep" SST (unchanged in BOLAM or GLOBO runs)

! Sea ice fraction, ice temperature and surface level temperature
! (sea ice thickness not available in IFS data: set below as a function of fice)

!   print*, "iana, jana", iana,jana
!   call plotout(fmaske,iana,jana,99)
!   call plotout(ficeee,iana,jana,99)
!   call plotout(snowe,iana,jana,99)
!   call plotout(tsurfe,iana,jana,99)
!   call plotout(tge(:,:,1),iana,jana,99)
!   call plotout(workf(:,:),iana,jana,99)
!   call plotout(t80e(:,:),iana,jana,99)
!   call plotout(q80e(:,:),iana,jana,99)
!   call plotout(p80e(:,:)/100.,iana,jana,99)
!   call plotout(u80e(:,:),iana,jana,99)
!   call plotout(v80e(:,:),iana,jana,99)
!   stop

   if (input_model == 'IFS') then

     if(maxval(ficeee).gt.0.02) then
       call seatemp(ficeee,fmaske,ficee,iana,jana,3,1,0.9)
     else
       ficee = 0.
     endif

     if(minval(ticeee).gt.180..and.maxval(ficeee).gt.0.02) then  ! excluding case of ticeee not defined
       call seatemp(ticeee,ficeee,ticee,iana,jana,3,1,0.9)   ! expands tice towards the ice-free sea
     else
       ticee(:,:) = min(workf(:,:), 271.4)  ! workf contains a temperature suitable for sea and lake ice
       ticee(:,:) = max(ticee(:,:), 230.)
     endif

   elseif (input_model == 'GFS') then

! For GFS, definition of ticeee from tge

     ticeee(:,:) = tge(:,:,1)

     do j = 1, jana
     do i = 1, iana
     ticeee(i,j) = workf(i,j)  ! workf contains a temperature suitable for sea and lake ice
     if(iflag_80m.eq.1) then
       ticeee(i,j) = max(ticeee(i,j), t80e(i,j)-10.)
     else
       ticeee(i,j) = max(ticeee(i,j), te(i,j,nlevp)-14.)
     endif
     ticeee(i,j) = min(ticeee(i,j), 271.4)
     enddo
     enddo

     call seatemp(ficeee,fmaske,ficee,iana,jana,3,1,0.9)
     call seatemp(icethee,fmaske,icethe,iana,jana,3,1,0.9)
     call seatemp(ticeee,ficeee,ticee,iana,jana,3,1,0.9)   ! expands tice towards the ice-free sea

   else                           ! BOLAM

     ticeee(:,:) = tge(:,:,1)

     call seatemp(ficeee,fmaske,ficee,iana,jana,4,1,0.9)
     call seatemp(icethee,fmaske,icethe,iana,jana,4,1,0.9)
     call seatemp(ticeee,ficeee,ticee,iana,jana,4,1,0.9)   ! expands tice towards the ice-free sea

   endif

! Landtemp extends T, Q and albedo from the land towards the sea:
! these variables are redefined, so they must not be used to define quantities over the sea

  if(minval(tsurfe).gt.100.) call landtemp(tsurfe,fmaske,iana,jana,4,0,0.9)
  do k = 1,nlevg_inp
    call landtemp(tge(1,1,k),fmaske,iana,jana,4,0,0.9)
    call landtemp(qge(1,1,k),fmaske,iana,jana,4,0,0.9)
  enddo

! --------------------------------------------------
! Horizontal interpolation of surface fields
! --------------------------------------------------

  if (input_model == "BOLAM") then
    al = 0.6
  else
    al = 0.7
  endif
  if(frame) al = 1.     ! in case of data on frames only, bilinear interpolation is used

  if(minval(tsurfe).gt.100.) call interp_spline_2d(tsurfe,iana,jana,xe,ye,xxt,yyt,ntot,tsurfi,al)
  if (input_model == "BOLAM") then
    call interp_spline_2d(qvsurfe,iana,jana,xe,ye,xxt,yyt,ntot,qvsurf,al)
  endif
  call interp_spline_2d(sste,iana,jana,xe,ye,xxt,yyt,ntot,ssti,al)
  call interp_spline_2d(ticee,iana,jana,xe,ye,xxt,yyt,ntot,tice,al)
  call interp_spline_2d(ficee,iana,jana,xe,ye,xxt,yyt,ntot,fice,al)
  call interp_spline_2d(icethe,iana,jana,xe,ye,xxt,yyt,ntot,iceth,al)
  call interp_spline_2d(sstcle,iana,jana,xe,ye,xxt,yyt,ntot,sstcli,al)
  do k = 1, nlevg_inp
  call interp_spline_2d (tge(1,1,k), iana, jana, xe, ye, xxt, yyt, ntot, tgii(1,1,k), al)
  call interp_spline_2d (qge(1,1,k), iana, jana, xe, ye, xxt, yyt, ntot, qgii(1,1,k), al)
  enddo
  al = 0.9
  call interp_spline_2d(snowe,iana,jana,xe,ye,xxt,yyt,ntot,snow,al)

! Modification of SST by assimilating SST from ISAC-MyOcean product (Mediterranean,
! Black Sea and near Atlantic) - useful mainly for GFS input model (ECMWF sst is very
! similar to the MyOcean data). The input grid can be in geographical or rotated coordinates. is accepted).
! It is advised not to call sst_isac in the case of Bolam input model for self nesting:
! the SST at deep output-model levels should not be modified, since sstcli (used below to
! define deep sea levels tg) is used by Bolam as reference sst for relaxation.
! However, if "BOLAM" indicates an initial condition derived from Globo and MyOcean data
! are available, calling sst_isac modifies all levels of sea temperature (through sstcli).

  if (ifile == 1) then
    call sst_isac(ssti,xxtg,yytg,nlon,nlat,ntot,dlon,dlat,alon0,alat0,x0d,y0d,idate0,flag_myocean)
    if (flag_myocean == 1) sstcli(:,:) = ssti(:,:) ! case of ssti modified by ISAC-MyOcean product
  endif

!! call outgraph(80102,1,nlon,nlat,x0d,y0d,alon0,alat0+dlat*0.5,dlon,dlat,&
!! 2, 0,  0,  1, 1, 0,idate0(1:5),iperiod(1:3),fmask(:,:),1.,0.)
!! call outgraph(80102,1,nlon,nlat,x0d,y0d,alon0,alat0+dlat*0.5,dlon,dlat,&
!! 0, 0,  0,  1, 1, 0,idate0(1:5),iperiod(1:3),ssti(:,:),1.,-273.15)

  if (ifile == 1) then ! At the first instant only

! Definition of all physiographic parameters on the model grid:
! topography height, variance of topography height,
! land-sea fraction, soil and vegetation parameters;
! soil temperature and water content on the bottom level and also at
! all soil levels in following some climatic hypothesis and surface data

    print *

! -----------------------------------------------------------------------

! Definition of all constant (in time) model physiographical parameters

! If model_param_constant.bin exists, then reads constant model
! physiographical parameters from this file,
! else parameters are defined using def_soil.F90 procedures

! If model_param_constant.bin exists (flag_constant_fields=0), then reads constant model
! physiographical parameters from this file,
! else (flag_constant_fields/=0) parameters are defined using
! subroutine physiographic_param (in def_soil.F90)

! Physiographical parameters variable in time (LAI, vegetation frac,
! soil temperature and soil water content vertical approximated profies)
! are difined by subroutine physiographic_param (in def_soil.F90)

    call read_param_const(x0d, y0d, alon0, alat0+dlat*0.5, htopi, zhtopvar, flag_constant_fields)

    print *

    call physiographic_param (nlon, nlat, nst, nvt, nlevg, &
 alon0, alat0+dlat*0.5, dlon, dlat, x0d, y0d, xxtg, yytg, tsurfi, soil_lev, ndayr, flag_constant_fields, &
 htopi, zhtopvar, fmask, &
 water_table_depth, tg_bottom, qg_rel_bottom, qg_rel_surf_approx, &
 ind_lev_soil_h_bottom, ind_lev_soil_w_bottom, &
 tg_first_guess, qg_rel_first_guess,&
 soil_map, veg_map, &
 soil_qmax, soil_qmin, soil_c, soil_rho, soil_psi, soil_k, &
 soil_par_b, soil_par_c, soil_qrel_wilt, soil_qrel_ref, &
 soil_albedo_dry, soil_albedo_wet, &
 soil_emiss1_dry, soil_emiss1_wet, soil_emiss2_dry, soil_emiss2_wet, &
 veg_lai, veg_frac, veg_root_depth, veg_roughness, veg_albedo, veg_emiss1, veg_emiss2, veg_lai_max, &
 snow_dirt)

    phig(:,:) = htopi(:,:)*g0

! call outgraph(80102,1,nlon,nlat,x0d,y0d,alon0,alat0+dlat*0.5,dlon,dlat,&
! 2, 0,  0,  1, 1, 0,idate0(1:5),iperiod(1:3),fmask(:,:),1.,0.)
! call outgraph(80102,1,nlon,nlat,x0d,y0d,alon0,alat0+dlat*0.5,dlon,dlat,&
! 0, 0,  0,  1, 1, 0,idate0(1:5),iperiod(1:3),veg_lai(:,:),1.,0.)

  endif ! end of definitions at the first instant only

! Reset of snow, fice, iceth unfeasible values and blend of ssti inland (used for lakes)

  snow(:,:)  = max(snow(:,:),0.)
  fice(:,:)  = max(fice(:,:), 0.)
  fice(:,:)  = min(fice(:,:), 1.)
  iceth(:,:) = max(iceth(:,:), 0.)
  tice(:,:)  = min(tice(:,:), 271.4)
  do j = 1, nlat
  do i = 1, nlon
  if(fmask(i,j).lt.0.5) then
  fice(i,j)  = 0.
  iceth(i,j) = 0.
  ssti(i,j) = 0.85*ssti(i,j) + 0.15*0.5*(tgii(i,j,nlevg_inp/2)+tgii(i,j,nlevg_inp/2+1))
  endif
  enddo
  enddo

! --------------------------------------------------
! Vertical interpolation of soil parameters (temp. and water cont.)
! --------------------------------------------------

  alf = 0.
  ex1 = 1.
  ex2 = 1.
  z1=sum(lev_list_soil(1:nlevg_inp))
  z2=sum(soil_lev(1:nlevg))

  if (abs(z1-z2)>1.e-1.or.nlevg_inp<nlevg) then  ! case vertical interpolation is needed

! lev_list_soil: levels of input model soil
! soil_lev: levels and depths of output model soil layers

    zsoil(1:nlevg) = sqrt(soil_lev(1:nlevg)) ! to make the depth coordinate more linear
    zsoil_work(1:nlevg_inp) = sqrt(lev_list_soil(1:nlevg_inp)) ! to make the depth coordinate more linear

    do j = 1, nlat
    do i = 1, nlon
  
! Temperature of soil and sea water

     tsoile(1:nlevg_inp) = tgii(i,j,1:nlevg_inp)

     if (fmask(i,j) < 0.5.or.fice(i,j) >= 0.8) then

! Soil or thick sea ice

        tgi(i,j,1:nlevg) = tg_first_guess(i,j,1:nlevg)
  
        nk = ind_lev_soil_h_bottom(i,j)

        if (soil_lev(nk) >= lev_list_soil(nlevg_inp)+0.5) then

! The bottom temperature soil level is below the bottom soil level in input data:
! temperature at upper soil levels is defined by vertical interpolation of
! input data and temperature at lower soil level is defined by first guess temperature

          do jklev = 1,nk
            if (soil_lev(jklev) >= lev_list_soil(nlevg_inp)+0.5) exit
          enddo
          nk = jklev
          zsoil_work(nlevg_inp+1) = zsoil(nk)
          tsoile(nlevg_inp+1) = tg_first_guess(i,j,nk)
          call near (zsoil(1:nk), nk, zsoil_work(1:nlevg_inp+1), nlevg_inp+1, ivg(1:nk))
          call interp_spline_1d (tsoili(1:nk), zsoil(1:nk), nk, tsoile(1:nlevg_inp+1), zsoil_work(1:nlevg_inp+1), &
   nlevg_inp+1, ivg(1:nk), alf, ex1, ex2)
          tgi(i,j,1:nk) = tsoili(1:nk)

        else

! The bottom temperature soil level is upper then bottom soil level in input data:
! temperature at all soil levels is defined by vertical interpolation of input data

          call near (zsoil(1:nk), nk, zsoil_work(1:nlevg_inp), nlevg_inp, ivg(1:nk))
          call interp_spline_1d (tsoili(1:nk), zsoil(1:nk), nk, tsoile(1:nlevg_inp), zsoil_work(1:nlevg_inp), &
   nlevg_inp, ivg(1:nk), alf, ex1, ex2)
          tgi(i,j,1:nk) = tsoili(1:nk)

        endif

      else

! Sea water

        call near (zsoil(1:nlevg), nlevg, zsoil_work(1:nlevg_inp), nlevg_inp, ivg(1:nlevg))
        call interp_spline_1d (tsoili(1:nlevg), zsoil(1:nlevg), nlevg, tsoile(1:nlevg_inp), zsoil_work(1:nlevg_inp), &
   nlevg_inp, ivg(1:nlevg), alf, ex1, ex2)
        tgi(i,j,1:nlevg) = tsoili(1:nlevg)

      endif

! Relative soil water content

      qsoile(1:nlevg_inp) = qgii(i,j,1:nlevg_inp)
    
      nk = ind_lev_soil_w_bottom(i,j)

      if (nk > 0) then ! No water body, no glacier

        qgi(i,j,1:nlevg) = qg_rel_first_guess(i,j,1:nlevg)

        if (soil_lev(nk) >= lev_list_soil(nlevg_inp)+0.5) then
  
! The bottom soil water content level is below the bottom soil level in input data:
! water content at upper soil levels is defined by vertical interpolation of
! input data and water content at lower soil level is defined by first guess temperature
  
          do jklev = 1,nk
            if (soil_lev(jklev) >= lev_list_soil(nlevg_inp)+0.5) exit
          enddo
          nk = jklev
          zsoil_work(nlevg_inp+1) = zsoil(nk)
          qsoile(nlevg_inp+1) = qg_rel_first_guess(i,j,nk)
          call near (zsoil(1:nk), nk, zsoil_work(1:nlevg_inp+1), nlevg_inp+1, ivg(1:nk))
          call interp_spline_1d (qsoili(1:nk), zsoil(1:nk), nk, qsoile(1:nlevg_inp+1), zsoil_work(1:nlevg_inp+1), &
     nlevg_inp+1, ivg(1:nk), alf, ex1, ex2)
          qgi(i,j,1:nk) = qsoili(1:nk)
  
        else
  
! The bottom soil water content level is upper then bottom soil level in input data:
! water content at all soil levels is defined by vertical interpolation of input data
  
          call near (zsoil(1:nk), nk, zsoil_work(1:nlevg_inp), nlevg_inp, ivg(1:nk))
          call interp_spline_1d (qsoili(1:nk), zsoil(1:nk), nk, qsoile(1:nlevg_inp), zsoil_work(1:nlevg_inp), &
     nlevg_inp, ivg(1:nk), alf, ex1, ex2)
          qgi(i,j,1:nk) = qsoili(1:nk)
    
        endif

      else
        
        qgi(i,j,1:nlevg) = 0.

      endif ! No water body, no glacier
  
    enddo
    enddo

  else

    tgi(:,:,1:nlevg) = tgii(:,:,1:nlevg)
    qgi(:,:,1:nlevg) = qgii(:,:,1:nlevg)

  endif

  qgi(:,:,:) = min(max(qgi(:,:,:),0.),1.) ! relative value

! Redistribution of snow depending on orography

  if (input_model /= "BOLAM") then
    call redistr_snow(snow, htopi, nlon, nlat, 7)
  else
    if (maxval(snow) > 1.e-2) call redistr_snow(snow, htopi, nlon, nlat, 4)
  endif
i=474; j=391
print *
print *,'poisk fmask',fmask(i,j)
print *,'poisk htopi',phig(i,j)/g0
print *,'poisk zhtopvar',zhtopvar(i,j)
print *,'poisk tg_bottom, qg_rel_bottom, qg_rel_surf_approx, water_table_depth',&
 tg_bottom(i,j), qg_rel_bottom(i,j), qg_rel_surf_approx(i,j), water_table_depth(i,j)
print *,'poisk ind_lev_soil_h_bottom ind_lev_soil_w_bottom',ind_lev_soil_h_bottom(i,j),ind_lev_soil_w_bottom(i,j)
print *,'poisk soil_map',soil_map(i,j,:)
print *,'poisk veg_map',veg_map(i,j,:)
print *,'poisk soil_qmax',soil_qmax(i,j,:)
print *,'poisk soil_qmin',soil_qmin(i,j,:)
print *,'poisk soil_c',soil_c(i,j,:)
print *,'poisk soil_rho',soil_rho(i,j,:)
print *,'poisk soil_psi',soil_psi(i,j,:)
print *,'poisk soil_k',soil_k(i,j,:)
print *,'poisk soil_par_b',soil_par_b(i,j,:)
print *,'poisk soil_qrel_wilt',soil_qrel_wilt(i,j,:)
print *,'poisk soil_qrel_ref',soil_qrel_ref(i,j,:)
print *,'poisk soil_albedo_dry',soil_albedo_dry(i,j)
print *,'poisk soil_albedo_dry',soil_albedo_dry(i,j)
print *,'poisk soil_albedo_wet',soil_albedo_wet(i,j)
print *,'poisk soil_emiss1_dry',soil_emiss1_dry(i,j)
print *,'poisk soil_emiss1_wet',soil_emiss1_wet(i,j)
print *,'poisk soil_emiss2_dry',soil_emiss2_dry(i,j)
print *,'poisk soil_emiss2_wet',soil_emiss2_wet(i,j)
print *,'poisk veg_root_depth',veg_root_depth(i,j)
print *,'poisk veg_roughness',veg_roughness(i,j)
print *,'poisk veg_albedo',veg_albedo(i,j)
print *,'poisk veg_emiss1',veg_emiss1(i,j)
print *,'poisk veg_emiss2',veg_emiss2(i,j)
print *,'poisk snow_dirt',snow_dirt(i,j)
print *,'poisk veg_lai_max',veg_lai_max
print *,'poisk tg_first_guess',tg_first_guess(i,j,:)
print *,'poisk qg_rel_first_guess',qg_rel_first_guess(i,j,:)
print *,'poisk veg_lai',veg_lai(i,j)
print *,'poisk veg_frac',veg_frac(i,j)
print *,'poisk tgi',tgi(i,j,:)
print *,'poisk qgi',qgi(i,j,:)
print *

 endif ! surf_elaborate

! --------------------------------------------------
!  Only in the case of hybrid levels: horizontal interpolation of IFS orography
!  (for boundary files if available, but necessary in case of different resolution)
!  and horizontal interpolation of surface pressure.
!  Note: PHIGE and PHIGE2D are the same for GFS data and different for IFS data (hybrid levels).
!  If IFILE >1, PHIGE2D is not defined for IFS data.
! --------------------------------------------------

 if (level_type == 2) then
  if (input_model == "BOLAM") then
    al = 0.7
  else
    al = 0.9
  endif
  if(frame) al = 1. ! in case of data on frames only, bilinear interpolation is used

  if(minval(phige).gt.-9000.) then
    call interp_spline_2d(phige,iana,jana,xe,ye,xxt,yyt,ntot,phigi,al)
  elseif (minval(phige).lt.-9000. .and. res_change == 1) then
    print*, "Phige (3-D model-level surf. geopotential) is not available"
    print*, "in IFS-ECMWF input data at instant", ifile
    print*, "and grid resolution has changed from the previous instant."
    print*, "This implies large errors in computed geopotential: prebolam stops!"
    stop
  else
    print*, "Phige (3-D model-level surf. geopotential) is not available"
    print*, "in IFS-ECMWF input data at instant", ifile
    print*, "but grid resolution has not changed from the previous instant, therefore"
    print*, "the interpolated field (phigi) of the previous instant is used."
  endif

  if (minval(phige2d).gt.-9000.) then
    call interp_spline_2d(phige2d,iana,jana,xe,ye,xxt,yyt,ntot,phigi2d,al)
  else
    print*, "Phige2d (2-D surface geopotential) is not available in IFS-ECMWF"
    print*, "input data at instant", ifile
    print*, "therefore the interpolated field (phigi2d) of the previous instant is used."
    if(res_change == 1) then
      print*, "Since grid resolution has changed, this may imply some errors"
      print*, "in the elaboration of surface fields."
    endif
  endif
  call interp_spline_2d(psel,iana,jana,xe,ye,xxt,yyt,ntot,psil,al)

 endif ! level_type=2 (hybrid levels)

! --------------------------------------------------
!  Only in the case of IFS-ECMWF isobaric level data:
!  horizontal interpolation of ECMWF orography geopot.
! (first instant only considered, even in case of resolution change
!  of boundary files with respect to initial condition)
! --------------------------------------------------

 if(input_model == "IFS".and.level_type.eq.1.and.ifile.eq.1) then
  if(minval(phige2d).gt.-9000.) then
    al = 0.9
    call interp_spline_2d(phige2d,iana,jana,xe,ye,xxt,yyt,ntot,phigi2d,al)
  else
    print*, "Phige2d (input model surface geopotential) is not available"
    print*, "in IFS-ECMWF data at the first instant: prebolam stops!"
    stop
  endif
 endif

! --------------------------------------------------
!  Only in the case of GFS:
!  horizontal interpolation of GFS orography geopot. (for boundary files if available)
!  and horizontal interpolation of fields at the level of 80 m above surface
! --------------------------------------------------

 if(input_model == "GFS") then

  if(minval(phige2d).gt.-9000.) then
    al = 0.9
    call interp_spline_2d(phige2d,iana,jana,xe,ye,xxt,yyt,ntot,phigi2d,al)
  else
    print*, "Phige2d (input model surface geopotential) is not available at this instant"
    print*
  endif

  if (iflag_80m.eq.1) then
   al = 0.9
   work(:,:) = log(p80e(:,:))
   call interp_spline_2d(t80e,iana,jana,xe,ye,xxt,yyt,ntot,t80i,al)
   call interp_spline_2d(q80e,iana,jana,xe,ye,xxt,yyt,ntot,q80i,al)
   call interp_spline_2d(work,iana,jana,xe,ye,xxt,yyt,ntot,pl80i,al)
   call interp_spline_2d(u80e,iana,jana,xe,ye,xxu,yyu,ntot,u80u,al)
   call interp_spline_2d(u80e,iana,jana,xe,ye,xxv,yyv,ntot,u80v,al)
   call interp_spline_2d(v80e,iana,jana,xe,ye,xxu,yyu,ntot,v80u,al)
   call interp_spline_2d(v80e,iana,jana,xe,ye,xxv,yyv,ntot,v80v,al)
   q80i(:,:)  = max(1.e-7, q80i(:,:))
   tv80i(:,:) = (1.+ep*q80i(:,:))*t80i(:,:)

! Rotation of the wind components on the rotated grid

   if (inp_inp_flag_rot==0) then
     call rot_wind(u80u,v80u,u80v,v80v,xxu,xxv,yyu,yyv,u80i,v80i, &
                   nlon,nlat,x0d,y0d,alon0,alat0,dlon,dlat)
   else
     u80i(:,:) = u80u(:,:)
     v80i(:,:) = v80v(:,:)
   endif

  endif  ! iflag_80m
 endif   ! input_model=GFS

! ----------------------------------------------
! Definition of fields on output model levels
! ----------------------------------------------
! Note that the same 3D matrices P and PL are used for pressure on input model levels,
! for both cases of input pressure levels (in this case only the vertical index varies)
! and of hybrid levels (considered later)
! Note also that the no. of press. levels can be different at different instants

 if (level_type==1) then  ! Isobaric levels

! Definition of auxiliary isobaric levels (logaritmic interp.):
! P(1,1,K)      standard levels
! PL(1,1,K)     logP(K)
! PLAUX(K)      log Paux (interpolated)

  i0 = 1
  j0 = 1
  do k = 1,nlevp
    p(i0,j0,k)   = lev_list(k)
    pl(i0,j0,k)  = log(p(i0,j0,k))
    plaux(2*k-1) = pl(i0,j0,k)
  enddo
  do k = 1,nlevp-1
    plaux(2*k) = (pl(i0,j0,k)+pl(i0,j0,k+1))/2.
  enddo

! GFS only:
! recomputation of geopotential at isobaric levels below ground surface of the input model:
! sometimes these values are not consistent with the average temperature in the
! layer above, so using them for computing the temperature at auxiliary levels can
! produce unrealistic values; therefore, the geopotential is recomputed using the
! virtual temperature of the layer above.
! The correction is applied depending on the level difference.

  if (input_model == 'GFS') then

  do j = 1,jana
  do i = 1,iana

! Search of the index of the first level above ground

    ki = 0
    do k = 1,nlevp
      if(phige2d(i,j) < ge(i,j,k)) ki = k
    enddo

! Geopotential correction starting from above (the updated geopotential is
! computed starting from the first level located below ground KI+1).
! A compromise solution is applied: weighted average between the old and the
! new geopotential, corrected as a function of the distance between the GFS orography
! and the level considered (ZWF=0.: no correct.; ZWF=1.: full correct.)

    if (ki < nlevp) then
      do k = ki+1,nlevp
        zdiftp = max(phige2d(i,j)-ge(i,j,k), 0.)
        zwf = min(zdiftp/1500.,1.)
        zge = ge(i,j,k)
        ge(i,j,k) = ge(i,j,k-1) + alog(p(i0,j0,k-1)/p(i0,j0,k))*rd*.5*          &
                   ((1.+ep*qe(i,j,k))*te(i,j,k)+(1.+ep*qe(i,j,k-1))*te(i,j,k-1))
       ge(i,j,k) = zwf*ge(i,j,k)+(1.-zwf)*zge
      enddo
    endif

  enddo
  enddo

! GFS only: horizontal filtering of the upper level (atmospheric) input fields

!    print*, "Horizontal smoothing applied to GFS fields at stratospheric levels"
    wei = 0.5
    do k = 1,nlevp
    zzz = -8000.*log(p(1,1,k)/1000.e2)
    nsmooth = int((zzz/1000.)**2/45.)
    nsmooth = nsmooth/4
    if( ifile == 1.and.nsmooth > 0) &
 write(*,'(a,i4,a,f8.3,a,i3)') " GFS lev.",nlevp-k+1,",  p(1,1) hPa",p(1,1,k)/100.,", nsmooth",nsmooth
    if (nsmooth >= 1) then
    call smooth_soil(ge(1,1,k), work,iana,jana,wei,nsmooth)
    call smooth_soil(te(1,1,k), work,iana,jana,wei,nsmooth)
    call smooth_soil(ue(1,1,k), work,iana,jana,wei,nsmooth)
    call smooth_soil(ve(1,1,k), work,iana,jana,wei,nsmooth)
    call smooth_soil(qe(1,1,k), work,iana,jana,wei,nsmooth)
    call smooth_soil(qce(1,1,k),work,iana,jana,wei,nsmooth)
    endif
    enddo

  endif ! GFS

! The following definition of gamma_inp is only for the case of isobaric levels (with or without frame).
! The two lowest levels of the input model are used

   gamma_inp = 0.
   iii = 0
   do j = 1,jana
   do i = 1,iana
     if (mask_frame(i,j)==1) then
      if(xe(i)>=aminx.and.xe(i)<=amaxx.and.ye(j)>=aminy.and.ye(j)<=amaxy) then
       iii = iii+1
       gamma_inp = gamma_inp + (te(i,j,nlevp)-te(i,j,nlevp-1))/(ge(i,j,nlevp-1)-ge(i,j,nlevp))*g0
      endif
     endif
   enddo
   enddo
     if(iii.eq.0) then
     print*, "iii = 0: gamma_inp cannot be computed: stop"
     stop
     endif
   gamma_inp = gamma_inp/float(iii)
   write(*,'(a,f8.4)') " Mean lapse rate of input model near the surface:", gamma_inp
!   print*, "Mean lapse rate computed on", iii, "points"

 else  ! Hybrid levels

   if (input_model == "IFS") then

    do j = 1, nlat
    do i = 1, nlon

! Computation of P on half levels, with the expression:
! P(k+1/2) = A(k+1/2) + B(k+1/2)*PS

      do k = 1, nlev_atm_inp_max+1
        ph_input(k) = ak(k) + bk(k)*exp(psil(i,j))
      enddo

! Definition of P on full levels, with the expression:
! P(k) = 0.5*(P(k+1/2) + P(k-1/2))

      do k = 1, nlev_atm_inp_max
        p_input(k) = 0.5*(ph_input(k) + ph_input(k+1))
      enddo

      do k = 1, nlevp
        k1 = int(lev_list(k))
        phl(i,j,k) = ph_input(k1)
        p(i,j,k) = p_input(k1)
        pl(i,j,k) = log(p(i,j,k))
      enddo
      k1 = int(lev_list(nlevp)) + 1
      phl(i,j,nlevp+1) = ak(k1) + bk(k1)*exp(psil(i,j))

    enddo
    enddo

    do j = 1, jana
    do i = 1, iana

! Computation of P on half levels, with the expression:
! P(k+1/2) = A(k+1/2) + B(k+1/2)*PS

      do k = 1, nlev_atm_inp_max+1
        ph_input(k) = ak(k) + bk(k)*exp(psel(i,j))
      enddo

! Definition of P on full levels, with the expression:
!  P(k) = 0.5*(P(k+1/2) + P(k-1/2))

      do k = 1, nlevp
        k1 = int(lev_list(k))
        phle(i,j,k) = ph_input(k1)
      enddo
      k1 = int(lev_list(nlevp))+1
      phle(i,j,nlevp+1) = ak(k1) + bk(k1)*exp(psel(i,j))

    enddo
    enddo

   else          ! BOLAM (defining the hybrid BOLAM coordinates)

     p0e   = ak(1)
     alfae = ak(2)
     sige(1) = 0.
     sige(2:nlevp+1) = bk(1:nlevp)
     sige2(1:nlevp) = bk(1:nlevp)

! Computation of P on half levels

     do k = 1,nlevp
       k1 = int(lev_list(k))
       phl(:,:,k)  = p0e*sige(k1)-(p0e-exp(psil(:,:)))*sige(k1)**3*(1.+alfae*(1.-sige(k1))*(2.-sige(k1)))
       phle(:,:,k) = p0e*sige(k1)-(p0e-exp(psel(:,:)))*sige(k1)**3*(1.+alfae*(1.-sige(k1))*(2.-sige(k1)))
     enddo
     k1 = int(lev_list(nlevp))+1
     phl (:,:,nlevp+1) = p0e*sige(k1)-(p0e-exp(psil(:,:)))*sige(k1)**3*(1.+alfae*(1.-sige(k1))*(2.-sige(k1)))
     phle(:,:,nlevp+1) = p0e*sige(k1)-(p0e-exp(psel(:,:)))*sige(k1)**3*(1.+alfae*(1.-sige(k1))*(2.-sige(k1)))

! Definition of P on full levels

    do k=1,nlevp
      siginte(k) = 0.5*(sige(k)+sige(k+1))
    enddo

     do k = 1,nlevp
       p(:,:,k) = p0e*siginte(k)-(p0e-exp(psil(:,:)))*siginte(k)**3*(1.+alfae*(1.-siginte(k))*(2.-siginte(k)))
       pl(:,:,k) = log(p(:,:,k))
     enddo

   endif

  if (input_model == 'IFS' .and. .not.frame) then

! Horizontal filtering of the upper level (hybrid atmospheric levels) input fields

    print*, "Horizontal smoothing applied to IFS fields at stratospheric hybrid levels"
    wei = 0.5
    do k = 1,nlevp
    zzz = -8000.*log(p(1,1,k)/1000.e2)
    nsmooth = int((zzz/1000.)**2/45.)
    nsmooth = nsmooth/3  ! more filtering than for GFS data due to higher resol.
    if( ifile == 1.and.nsmooth > 0) print*, "IFS lev.",nlevp-k+1,"p(1,1), hPa",p(1,1,k)/100.,"nsmooth",nsmooth
    if (nsmooth >= 1) then
    call smooth_soil(te(1,1,k), work,iana,jana,wei,nsmooth)
    call smooth_soil(ue(1,1,k), work,iana,jana,wei,nsmooth)
    call smooth_soil(ve(1,1,k), work,iana,jana,wei,nsmooth)
    call smooth_soil(qe(1,1,k), work,iana,jana,wei,nsmooth)
    call smooth_soil(qce(1,1,k),work,iana,jana,wei,nsmooth)
    endif
    enddo

  endif ! end of case IFS and not frame

! gamma_inp, in the case of hybrid levels, is computed from second and fourth levels above the ground
! (to avoid extreme values near the surface and assuming that hybrid levels have high vert. res.)

  k1 = nlevp
  k2 = nlevp-1
  k3 = nlevp-2
  k4 = nlevp-3
  gamma_inp = 0.
  iii = 0

  do j = 1,jana
  do i = 1,iana
    if (mask_frame(i,j)==1) then
     if(xe(i)>=aminx.and.xe(i)<=amaxx.and.ye(j)>=aminy.and.ye(j)<=amaxy) then
      iii = iii+1
      zpl1 = log(0.5*(phle(i,j,k1) + phle(i,j,k2)))
      zpl2 = log(0.5*(phle(i,j,k3) + phle(i,j,k4)))
      gamma_inp = gamma_inp+g0/rd*(log(te(i,j,k4))-log(te(i,j,k2)))/(zpl2-zpl1)
     endif
    endif
  enddo
  enddo

  gamma_inp = gamma_inp/float(iii)
  write(*,'(a,f8.4)') " Mean lapse rate of input model near the surface", gamma_inp
!  print*, "Mean lapse rate computed on", iii, "points"

 endif    ! level type

! Comput. of (minus) spatially averaged istantaneous lapse rate of the input model near the surface

 if (input_model /= "BOLAM") then
   gamma = max(0.35*gammac + 0.65*gamma_inp, 0.)
   write(*,'(a,f8.4)') " Blended lapse rate:", gamma
 else
   gamma = gammac
   write(*,'(a,f8.4)') " Lapse rate", gamma
 endif

! ----------------------------------------------
! Horizontal interpolation of the upper level fields
! ----------------------------------------------
! AL: spline tension parameter

  al = 0.3
  if(frame) al = 1. ! in case of data on frames only, bilinear interpolation is used

  if (level_type==1) then  ! isobaric levels
    do k = 1,nlevp
    call interp_spline_2d(ge(1,1,k),iana,jana,xe,ye,xxt,yyt,ntot,gp(1,1,k),al)
    enddo
  endif

  do k = 1,nlevp
    call interp_spline_2d(te(1,1,k),iana,jana,xe,ye,xxt,yyt,ntot,tp(1,1,k),al)
  enddo

  if (input_model == "BOLAM") then
    do k = 1,nlevp
      call interp_spline_2d(ue(1,1,k),iana,jana,xue,ye,xxu,yyu,ntot,upu(1,1,k),al)
      call interp_spline_2d(ue(1,1,k),iana,jana,xue,ye,xxv,yyv,ntot,upv(1,1,k),al)
      call interp_spline_2d(ve(1,1,k),iana,jana,xe,yve,xxu,yyu,ntot,vpu(1,1,k),al)
      call interp_spline_2d(ve(1,1,k),iana,jana,xe,yve,xxv,yyv,ntot,vpv(1,1,k),al)
    enddo
  else    ! IFS or GFS
    do k = 1,nlevp
      call interp_spline_2d(ue(1,1,k),iana,jana,xe,ye,xxu,yyu,ntot,upu(1,1,k),al)
      call interp_spline_2d(ue(1,1,k),iana,jana,xe,ye,xxv,yyv,ntot,upv(1,1,k),al)
      call interp_spline_2d(ve(1,1,k),iana,jana,xe,ye,xxu,yyu,ntot,vpu(1,1,k),al)
      call interp_spline_2d(ve(1,1,k),iana,jana,xe,ye,xxv,yyv,ntot,vpv(1,1,k),al)
    enddo
  endif

  al = 0.4
  do k = 1,nlevp
    call interp_spline_2d(qe(1,1,k),iana,jana,xe,ye,xxt,yyt,ntot,qp(1,1,k),al)
    call interp_spline_2d(qce(1,1,k),iana,jana,xe,ye,xxt,yyt,ntot,qcp(1,1,k),al)
  enddo

! Elimination of possible negative values of specific humidity

  qp(:,:,:)  = max(1.e-7, qp(:,:,:))
  qcp(:,:,:) = max(0.,   qcp(:,:,:))

! Rotation of the wind components on the rotated grid

  if (inp_inp_flag_rot==0) then
    do k = 1,nlevp
      call rot_wind(upu(1,1,k),vpu(1,1,k),upv(1,1,k),vpv(1,1,k),xxu,xxv,yyu,yyv,    &
                    upup(1,1,k),vpvp(1,1,k),nlon,nlat,x0d,y0d,alon0,alat0,dlon,dlat)
    enddo
  else
    upup(:,:,:) = upu(:,:,:)
    vpvp(:,:,:) = vpv(:,:,:)
  endif

! Geopotential at the ground surface (topography) of the output model

  phig(:,:) = g0*htopi(:,:)

!----------------------------------
! Start of vertical interpolation
!----------------------------------

! Definition of TV (virtual temperature) and PS (surface pressure)

  if (level_type==1) then ! isobaric levels

! p_tr: threshold in deltap/p (applied to differences of log(p)), used to incorporate 80m level data
! (it should correspond to at least 10-15 m to avoid noise in vert. interp.)

   p_tr = g0/(rd*275.)*18.

   i0 = 1
   j0 = 1
   do j = 1,nlat
   do i = 1,nlon

! Definition of virtual T at the auxiliary pressure levels

! The layer average temperature is imposed to be consistent with the hydrostatic eq.

     tv(i,j,1) = (1+ep*qp(i,j,1))*tp(i,j,1)
     gaux(1) = gp(i,j,1)
     do k = 1,nlevp-1
     k1 = 2*k+1
     tv(i,j,k1) = (1.+ep*qp(i,j,k+1))*tp(i,j,k+1)
     gaux(k1) = gp(i,j,k+1)
     k2 = 2*k
!     tvm = -(gp(i,j,k)-gp(i,j,k+1))/(pl(k)-pl(k+1))/rd  ! standard: assumes TVM as the mean Tv in the layer

! Alternative expression for TVM: after tests on GFS data compar. 26 and 47 levels (apr. 2010),
! a definition of delp depending on p is introduced, with a transition from DP/P at low levels
! to D(ln(p)) at higher levels - this empirically reduces "noise" al low levels but maintains
! realistic corrections at higher levels

     peso = 0.5*(p(i0,j0,k)+p(i0,j0,k+1))/1.1e5          ! weight for the denomin. (p/prif)
     delp1 = pl(i0,j0,k)-pl(i0,j0,k+1)                   ! delta(log(p))
     delp2 = 2.*(p(i0,j0,k)-p(i0,j0,k+1))/(p(i0,j0,k)+p(i0,j0,k+1)) ! deltap/pmed
     delp3 = peso*delp2+(1.-peso)*delp1
     tvm = -(gp(i,j,k)-gp(i,j,k+1))/delp3/rd

!    tv(i,j,k2) = 2.*tvm-(tv(i,j,k2-1)+tv(i,j,k2+1))/2.  ! standard: TVM as average TV in the layer (OK for IFS data)
!    tv(i,j,k2) = 0.5*(tv(i,j,k2-1)+tv(i,j,k2+1))        ! with this the auxiliary information from geop. is not used
!    tv(i,j,k2) = 0.25*(tv(i,j,k2-1)+tv(i,j,k2+1))+0.5*tvm ! possible intermediate expression
!    tv(i,j,k2) = 1.2*tvm-0.1*(tv(i,j,k2-1)+tv(i,j,k2+1))  ! intermediate expression

! Below: blending depending again on p, since empirically it seems that the correction based
! on the use of GPH data is more realistic at higher levels than at lower

     tv1 = 1.6*tvm-0.3*(tv(i,j,k2-1)+tv(i,j,k2+1))
     tv2 = 0.5*(tv(i,j,k2-1)+tv(i,j,k2+1))
     peso2 = sqrt(peso)
     tv(i,j,k2) = peso2*tv2+(1.-peso2)*tv1
     gaux(k2) = gaux(k1)-rd*(plaux(k2)-plaux(k1))*(tv(i,j,k2)+tv(i,j,k1))/2.
     enddo

! Computation of surface pressure using T virt. at the auxiliary levels, as computed above.
! The hydrostatic relation is used to extrapolate (or interpolate) TVIRT using
! geopotential as independent variable (instead of log(p) to avoid a second order equation).
! The geopotential is computed first at the auxiliary intermediate levels
! (consistent with the comput. of TVIRT at the same levels, using thickness between standard levels).

    ki = 0
    do k = 1,nauxp
    if (phig(i,j).lt.gaux(k)) ki = k
    enddo
    if (ki.eq.nauxp) ki = nauxp-1

    if (ki.eq.0) then
    print *,'Error: ki=0 in computing geopotential!'
    print *,'i = ',i,'  j = ',j
    print *,'phig(i,j)=',phig(i,j),' gaux(1)=',gaux(1),' gaux(nauxp)=',gaux(nauxp)
    stop 'Aborting.'
    endif

! TVS: defined below as virtual air temperature at the surface, extrapolated from air temp.
! In case of temp. inversion, TVS is averaged between constant value and extrapolated value
! if the output model surface (topography) is below the lowest isobaric level of the input model

      tvs(i,j) = tv(i,j,ki+1)+(tv(i,j,ki+1)-tv(i,j,ki))*(phig(i,j)-gaux(ki+1))/(gaux(ki+1)-gaux(ki))
      if(tv(i,j,ki+1).lt.tv(i,j,ki).and.(phig(i,j)-gaux(ki+1)).lt.0.) tvs(i,j) = 0.4*tvs(i,j)+0.6*tv(i,j,ki+1)

! Definition of surface pressure PST (logarithm of ps on T points) - closest level to topog. chosen

      if(abs(phig(i,j)-gaux(ki+1)).lt.(abs(phig(i,j)-gaux(ki)))) then
        kin = ki+1
      else
        kin = ki
      endif
      pst(i,j) = plaux(kin) - (2./rd)*(phig(i,j)-gaux(kin))/(tvs(i,j)+tv(i,j,kin))

! Data at 80 m are used below to recompute surface pressure only in the case the 80 m lev. is below the
! lowest isobaric level above the topography (minus an offset) - geopotential at 80 m must be recomputed
! consistently with the hydrostatic eqn. to prevent noise. In case of inversion, TVS is averaged as above.

    if(iflag_80m.eq.1) then
      ki80 = 0
      do k = 1,nauxp
      if (pl80i(i,j).gt.plaux(k)) ki80 = k
      enddo
      if (ki80.eq.nauxp) ki80 = nauxp-1

      if(abs(pl80i(i,j)-plaux(ki80+1)).lt.(abs(pl80i(i,j)-plaux(ki80)))) then
        ki80n = ki80+1
      else
        ki80n = ki80
      endif
      delphi = rd*0.5*(tv(i,j,ki80n)+tv80i(i,j))*(pl80i(i,j)-plaux(ki80n))
      g80i = gaux(ki80n) - delphi

      if(gaux(ki+1).gt.phig(i,j).and.g80i.lt.(gaux(ki+1)-15.*g0)) then
        pldo = pl80i(i,j)
        tvdo = tv80i(i,j)
        tvs(i,j) = tvdo + (tvdo-tv(i,j,ki+1))*(phig(i,j)-g80i)/(g80i-gaux(ki+1))
        if(tvdo.lt.tv(i,j,ki+1).and.(phig(i,j)-g80i).lt.0.) tvs(i,j) = 0.4*tvs(i,j)+0.6*tvdo
        pst(i,j) = pldo - (2./rd)*(phig(i,j)-g80i)/(tvs(i,j)+tvdo)
      endif
    endif ! case iflag_80m = 1

    ps(i,j) = exp(pst(i,j))

! End of loop on the grid points

  enddo
  enddo

  else  ! from here case of hybrid levels

! Definition of geopotential at input model levels, not available in input data.
! It is computed using PHIGI and pressure data.
! Geopotential is used to compute the surface pressure only (not for temperature definition).
! Fields, which are the result of horizontal interpolation, are used.

! Definition of virtual temperature

   do k = 1,nlevp
     tv(:,:,k) = (1.+ep*qp(:,:,k))*tp(:,:,k)
   enddo

! Definition of geopotential at internal levels (GP)

  gp(:,:,nlevp) = phigi(:,:)-rd*tv(:,:,nlevp)*(pl(:,:,nlevp)-psil(:,:))

  kf = nlevp-1
  do k = kf,1,-1
    gp(:,:,k) = gp(:,:,k+1)-rd*0.5*(tv(:,:,k)+tv(:,:,k+1))*(pl(:,:,k)-pl(:,:,k+1))
  enddo

  do j = 1,nlat
  do i = 1,nlon

! Definition of surface pressure PS
! Search of the index of the first level above ground

    ki = 0
    do k = 1,nlevp
      if (phig(i,j) < gp(i,j,k)) ki = k
    enddo

    if (ki==0) then
      print *,'Error: ki=0 in computing geopotential!'
      print *,'i = ',i,'  j = ',j
      print *,'phig(i,j)=',phig(i,j),' gp(i,j,1)=',gp(i,j,1)
      stop 'Aborting.'
    endif

! TVS: defined below as virtual air temperature at the surface, extrapolated from air temp.
! If the grid point is under the orography, TVS is computed as a weighted average
! between the result of the extrapolation with lapse rate GAMMA and
! the result of the linear extrapolation from higher levels.

   if (ki==nlevp) then

! In this case the ground surface is below the lowest level of the input model:
! in cases of inversion extrapol. is problematic, therefore only the lapse rate GAMMA
! is used; otherwise a weighted average of the lapse rate GAMMA
! and linear extrapolation from levels above is applied

    if ((tv(i,j,ki)-tv(i,j,ki-1)) > 0.) then
      tvs(i,j) = tv(i,j,ki)+0.7*(gamma*(gp(i,j,ki)-phig(i,j))/g0)                       &
         +0.3*((tv(i,j,ki)-tv(i,j,ki-1))*(phig(i,j)-gp(i,j,ki))/(gp(i,j,ki)-gp(i,j,ki-1)))
    else
      tvs(i,j) = tv(i,j,ki)+gamma*(gp(i,j,ki)-phig(i,j))/g0
    endif

! Extrapolation of surface pressure PS at T points
! PST contains the logarithm of PS on T points

    pst(i,j) = pl(i,j,ki)+(gp(i,j,ki)-phig(i,j))/(rd*0.5*(tv(i,j,ki)+tvs(i,j)))
    ps(i,j) = exp(pst(i,j))

   else

! In this case the ground surface is above the lowest level of the input model:
! interpolation is applied

    tvs(i,j) = tv(i,j,ki+1)+(tv(i,j,ki+1)-tv(i,j,ki))*(phig(i,j)-gp(i,j,ki+1))/(gp(i,j,ki+1)-gp(i,j,ki))

! Extrapolation of surface pressure PS at T points
! PST contains the logarithm of PS on T points

    pst(i,j) = pl(i,j,ki+1)+(gp(i,j,ki+1)-phig(i,j))/(rd*0.5*(tvs(i,j)+tv(i,j,ki+1)))
    ps(i,j) = exp(pst(i,j))

   endif

  enddo
  enddo

 endif ! end of distinction between level types

! Verification of PS minimum to check the value of ALFA that defines the
! hybrid coordinates is not consistent (too small for orography)

  zpsmin = minval(ps(:,:))
  call livelli(log(zpsmin),s,nlev,plsig,p0,alfa)
  write(*,'(a,f8.3,a,f8.3)') " Min(ps) =", zpsmin/100.," - corresp. p of the lowest model lev. =", exp(plsig(nlev))/100.
  if(zpsmin <= (exp(plsig(nlev)) + 40.)) then ! min. model level height above ground of about 4 m
   print*, "alfa =", alfa
   print*, "Error: alfa is too small for the orography in the domain."
   print*, "The value of alfa in prebolam.inp must be increased."
   stop ' Stop.'
  endif

! Definition of logarithm of pressure at the ground surface on u, v points

  do j = 1,nlat
    do i = 1,nlon-1
      psu(i,j) = (pst(i,j)+pst(i+1,j))/2.
    enddo
    psu(nlon,j) = pst(nlon,j)
  enddo

  do i = 1,nlon
    psv(i,1) = pst(i,1)
    do j = 2,nlat
      psv(i,j) = (pst(i,j)+pst(i,j-1))/2.
    enddo
  enddo

! ----------------------------------------------------------
! Plots for check:
! Intermediate working fields after the horizontal iterpolation
! from the input model grid into the output model grid
! ----------------------------------------------------------

!  if (ifile == 1) then
!    do k=1,nlevp
!      call plotout(tp(1,1,k),nlon,nlat,99)
!      call plotout(gp(1,1,k),nlon,nlat,99)
!      call plotout(p(1,1,k),nlon,nlat,99)
!      call plotout(qp(1,1,k),nlon,nlat,99)
!      call plotout(upup(1,1,k),nlon,nlat,99)
!      call plotout(vpvp(1,1,k),nlon,nlat,99)
!    enddo

!    call plotout(tp(1,1,nlevp),nlon,nlat,99)
!    call plotout(tp(1,1,nlevp-1),nlon,nlat,99)
!    call plotout(gp(1,1,nlevp),nlon,nlat,99)
!    call plotout(p(1,1,nlevp),nlon,nlat,99)
!    call plotout(qp(1,1,nlevp),nlon,nlat,99)
!    call plotout(upup(1,1,nlevp),nlon,nlat,99)
!    call plotout(vpvp(1,1,nlevp),nlon,nlat,99)
!    call plotout(phigi,nlon,nlat,99)
!    call plotout(phigi2d,nlon,nlat,99)
!    call plotout(psil,nlon,nlat,99)
!    call plotout(tsurfi,nlon,nlat,99)
!    call plotout(ssti,nlon,nlat,99)
!    call plotout(tgi(1,1,1),nlon,nlat,99)
!    call plotout(tgi(1,1,2),nlon,nlat,99)
!    call plotout(tgi(1,1,3),nlon,nlat,99)
!    call plotout(tgi(1,1,4),nlon,nlat,99)
!    call plotout(tgi(1,1,5),nlon,nlat,99)
!    call plotout(tgi(1,1,6),nlon,nlat,99)
!    call plotout(qgi(1,1,1),nlon,nlat,99)
!    call plotout(qgi(1,1,2),nlon,nlat,99)
!    call plotout(qgi(1,1,3),nlon,nlat,99)
!    call plotout(qgi(1,1,4),nlon,nlat,99)
!    call plotout(qgi(1,1,5),nlon,nlat,99)
!    call plotout(qgi(1,1,6),nlon,nlat,99)
!    call plotout(snow,nlon,nlat,99)
!    call plotout(htopi,nlon,nlat,99)
!    call plotout(fmask,nlon,nlat,99)
!    call plotout(fice,nlon,nlat,99)
!    call plotout(iceth,nlon,nlat,99)
!    call plotout(tice,nlon,nlat,99)
!    call plotout(t80i,nlon,nlat,99)
!    call plotout(q80i,nlon,nlat,99)
!    call plotout(exp(pl80i)/100.,nlon,nlat,99)
!    call plotout(u80i,nlon,nlat,99)
!    call plotout(v80i,nlon,nlat,99)
!  endif
!  stop

! ----------------------------------------------------------
! Vertical interpolation of fields at the atmospheric levels
! Start of the main loop on grid points
! ----------------------------------------------------------

! Interpolation parameters (spline tension, top and bottom extrap., in the order, for each variable)

  alf_t  = 0.6
  ex1_t  = 0.6
  ex2_t  = 0.95

  alf_q  = 0.7
  ex1_q  = 0.5
  ex2_q  = 0.85

  alf_qc = 0.7
  ex1_qc = 0.5
  ex2_qc = 0.8

  alf_uv = 0.6
  ex1_uv = 0.2
  ex2_uv = 0.75

  do j = 1,nlat
  do i = 1,nlon

! Computation of log of pressure at T points on the output model hybrid levels

    call livelli(pst(i,j),s,nlev,plsig,p0,alfa)

! Definition of working vectors (virt. temp. tvk is defined later below)

    do k = 1,nlevp
     upk(k) = upup(i,j,k)
     vpk(k) = vpvp(i,j,k)
     qpk(k) = qp (i,j,k)
     qcpk(k)= qcp(i,j,k)
    enddo

! Case in which 80 m level data are available.
! Note that, in vectors used to accomodate the 80 m data, the last (lowest) element is not
! defined in the case 80 m level data are discarded (in case they are close than p_tr to isobaric
! levels). K80 is the index of the inserted 80 m level - if k80=0, no 80 m level inserted.

    if(iflag_80m == 1) then
     i0 = 1
     j0 = 1
     k80 = 0
     do k = 2, nlevp
      if(pl80i(i,j).lt.pl(i0,j0,k).and.                                             &
       (pl(i0,j0,k)-pl80i(i,j)).gt.p_tr.and.(pl80i(i,j)-pl(i0,j0,k-1)).gt.p_tr) then
       k80 = k
       exit
      endif
     enddo

     if(k80.eq.0.and.(pl80i(i,j)-pl(i0,j0,nlevp)).gt.p_tr) k80 = nlevp + 1 ! case 80m level below lowest isobaric level

     if(k80.gt.0) then !  80 m level data used, vectors are filled from 1 to nlevp+1
      do k = 1, k80-1
       upk(k) = upup(i,j,k)
       vpk(k) = vpvp(i,j,k)
       qpk(k) = qp  (i,j,k)
      enddo
      upk(k80) = u80i(i,j)
      vpk(k80) = v80i(i,j)
      qpk(k80) = q80i(i,j)

      if(k80.lt.nlevp+1) then
       do k = k80+1, nlevp+1
        upk(k) = upup(i,j,k-1)
        vpk(k) = vpvp(i,j,k-1)
        qpk(k) = qp  (i,j,k-1)
       enddo
      endif
     endif  ! note that in the case k80 = 0 vectors are already defined
    endif  ! case iflag_80m == 1

    if (level_type==1) then  ! isobaric levels

     i0 = 1
     j0 = 1
     if(iflag_80m == 0 .or. k80 == 0) then
      pl2(1:nlevp) = pl(i0,j0,1:nlevp)
     else             ! the value of k80 defined above is used here
      do k = 1, k80-1
       pl2(k) = pl(i0,j0,k)
      enddo
      pl2(k80) = pl80i(i,j)

      if(k80.lt.nlevp+1) then
       do k = k80+1, nlevp+1
        pl2(k) = pl(i0,j0,k-1)
       enddo
      endif
     endif ! case iflag_80m=0 or k80=0

! Definition of working vector tvk (virtual t)

     if(iflag_80m == 0) then ! case of isobaric level data, no 80 m level data
      tvk(1:nauxp) = tv(i,j,1:nauxp)
     elseif(iflag_80m == 1) then  ! case of GFS model with 80 m level data
      k80a = 0
      do k = 1, nauxp
       if(pl80i(i,j) .le. plaux(k)) then
        k80a = k
        exit
       endif
      enddo

! Elimination of the 80 m inserted level if too close (in pressure) to an isobaric level

      if(k80a.gt.0) then
       if((plaux(k80a)-pl80i(i,j)).le.p_tr.or.(pl80i(i,j)-plaux(k80a-1)).le.p_tr) k80a = 0
      endif

! Case of 80m level below the lowest aux. isobaric level: 80 m level added

      if(k80a.eq.0.and.(pl80i(i,j)-plaux(nauxp)).gt.p_tr) k80a = nauxp + 1

      if(k80a.gt.0) then !  80 m level data used, vectors are filled from 1 to nauxp + 1
       do k = 1, k80a-1
        plaux_add80(k) = plaux(k)
        tvk(k) = tv(i,j,k)
       enddo
       plaux_add80(k80a) = pl80i(i,j)
       tvk(k80a) = tv80i(i,j)
       if(k80a.lt.nauxp+1) then
        do k = k80a+1, nauxp+1
         plaux_add80(k) = plaux(k-1)
         tvk(k) = tv(i,j,k-1)
        enddo
       endif
      else              !  k80a = 0: no 80 m level used, vectors are filled from 1 to nauxp
       do k = 1, nauxp
        plaux_add80(k) = plaux(k)
        tvk(k) = tv(i,j,k)
       enddo
      endif
     endif  ! case iflag_80m=1

! Interpolation of Tv
! In case of temp. inversion, Tv is averaged between constant value and extrapolated value
! (but in this case Tv is recomputed later below)

     if(iflag_80m == 0.or.k80a == 0) then
      call near(plsig,nlev,plaux,nauxp,iv)
      call interp_spline_1d(tk,plsig,nlev,tvk,plaux,nauxp,iv,alf_t,ex1_t,ex2_t)
      do k = 1,nlev
       if(plsig(k) > plaux(nauxp).and.tvk(nauxp) < tvk(nauxp-1)) tk(k) = 0.4*tk(k)+0.6*tvk(nauxp)
      enddo
     elseif(iflag_80m == 1.and.k80a > 0) then
      call near(plsig,nlev,plaux_add80,nauxp+1,iv)
      call interp_spline_1d(tk,plsig,nlev,tvk,plaux_add80,nauxp+1,iv,alf_t,ex1_t,ex2_t)
      do k = 1,nlev
       if(plsig(k) > plaux_add80(nauxp+1).and.tvk(nauxp+1) < tvk(nauxp)) tk(k) = 0.4*tk(k)+0.6*tvk(nauxp+1)
      enddo
     endif

    else   ! hybrid levels

     do k = 1,nlevp
      tvk(k) = tp(i,j,k)
      pl2(k) = pl(i,j,k)
     enddo

! Interpolation of T
! T is kept constant in case of extrapolation toward the ground surface with T inversion
! In case of temp. inversion, Tv is averaged between constant value and extrapolated value
! (but in this case T is recomputed below)

     call near(plsig,nlev,pl2,nlevp,iv)
     call interp_spline_1d(tk,plsig,nlev,tvk,pl2,nlevp,iv,alf_t,ex1_t,ex2_t)
     do k = 1,nlev
       if(plsig(k) > pl2(nlevp).and.tvk(nlevp) < tvk(nlevp-1)) tk(k) = 0.4*tk(k)+0.6*tvk(nlevp)
     enddo

    endif   ! end of condition on level types

! From here: for all level types (but distinction in case of 80 m level data)
! Interpolation of specific humidity and cloud content

    if(iflag_80m == 0.or.k80 == 0) then
     call near(plsig,nlev,pl2,nlevp,iv)
     call interp_spline_1d(qk,plsig,nlev,qpk,pl2,nlevp,iv,alf_q,ex1_q,ex2_q)
    elseif(iflag_80m == 1.and.k80 > 0) then
     call near(plsig,nlev,pl2,nlevp+1,iv)
     call interp_spline_1d(qk,plsig,nlev,qpk,pl2,nlevp+1,iv,alf_q,ex1_q,ex2_q)
    endif
    q(i,j,1:nlev) = max(1.e-7, qk(1:nlev))

    call interp_spline_1d(qck,plsig,nlev,qcpk,pl2,nlevp,iv,alf_qc,ex1_qc,ex2_qc)
    qc(i,j,1:nlev) = max(0.,qck(1:nlev))

    if (level_type == 1) then  ! case of isobaric levels input model

! Conversion of TK, that in case of isob. levels contains virt. temp., into temperature

     t(i,j,1:nlev) = tk(1:nlev)/(1+ep*q(i,j,1:nlev))

!  Only in case of IFS-ECMWF isobaric data: relative humidity limited and extrapolated
!  constant if the lowest input model level (1000 hPa) is above the surface.
!  Find output model level lbot located immediately above the lowest input model level

     if (input_model == "IFS") then
     p1 = lev_list(nlevp)
     do k = nlev, 1, -1
     if (exp(plsig(k)) <= p1) then
     lbot = k
     exit
     endif
     enddo
     if (lbot < nlev) then
     call qsat_tetens (t(i,j,lbot+1), exp(plsig(lbot+1)), eps, qsat, qsw, qsi, es, esw, esi)
     rh1 = min (q(i,j,lbot+1)/qsat, .96)
      do k = lbot+1, nlev
      pks = exp(plsig(k))
      call qsat_tetens (t(i,j,k), pks, eps, qsat, qsw, qsi, es, esw, esi)
      q(i,j,k) = rh1*qsat
      if (q(i,j,k)/qsat<0.95) qc(i,j,k) = 0.
      t(i,j,k) = tk(k)/(1+ep*q(i,j,k))          !?
      enddo
     endif
     endif  !  IFS (isobaric) data

!  Only in case of NOAA-GFS (isobaric) data: relative humidity limited and extrapolated
!  constant if the lowest input model level (1000 hPa or 80 m level) is above the surface.
!  Find output model level lbot located immediately above the lowest input model level

     if (input_model == "GFS") then

      if(iflag_80m == 0.or.k80 == 0) then
      p1 = exp(pl2(nlevp))
      elseif(iflag_80m == 1.and.k80 > 0) then
      p1 = exp(pl2(nlevp+1))
      endif

     do k = nlev, 1, -1
     if (exp(plsig(k)) <= p1) then
     lbot = k
     exit
     endif
     enddo
     if (lbot < nlev) then
     call qsat_tetens (t(i,j,lbot+1), exp(plsig(lbot+1)), eps, qsat, qsw, qsi, es, esw, esi)
     rh1 = min (q(i,j,lbot+1)/qsat, .99)
      do k = lbot+1, nlev
      pks = exp(plsig(k))
      call qsat_tetens (t(i,j,k), pks, eps, qsat, qsw, qsi, es, esw, esi)
      q(i,j,k) = rh1*qsat
      if (q(i,j,k)/qsat<0.92) qc(i,j,k) = 0.
      enddo
     endif
     endif  !  GFS (isobaric) data

    else  ! case of hybrid levels (terrain-following) input model

!  Find output model level lbot located immediately above the lowest input model level

     do k = nlev, 1, -1
       if (plsig(k) <= pl(i,j,nlevp)) then
        lbot = k
        exit
       endif
     enddo

! Redefinition of variables in the case when an output model level is located
! below the lowest level of the input model

     if(lbot < nlev) then

!  Under LBOT extrapolation with lapse rate GAMMA for T,
!  blended with the extrapolated values (computed above with interp_spline_1d).
!  zlaps: local lapse rate (comp. over the 3 lowest layers)
!  w1: weight given to the modified profiles of t and q defined below (function of local lapse rate)

      delzf = -0.5*rd/g0*(pl(i,j,nlevp-3) - pl(i,j,nlevp))*(tp(i,j,nlevp-3) + tp(i,j,nlevp))
      zlaps = (tp(i,j,nlevp) - tp(i,j,nlevp-3))/delzf
      w1 = min(1., -1./9.*zlaps*1000. + 1.)  ! w1=1. for zlaps<=0.deg/km; w1=0. for zlaps>=9.deg/km
      w1 = max(w1, .1)
!      w1 = 0.  ! pure extrap.
!      w1 = 1. ! pure redefinitions below

      omega = 0.5*rd/g0*(0.8*gamma+0.2*zlaps)
      beta = omega*(plsig(nlev)-pl(i,j,nlevp))
      tk(nlev) = w1*(tp(i,j,nlevp)*(1.+beta)/(1.-beta)) + (1.-w1)*tk(nlev)
      do k = lbot, nlev-1
        zk1 = pl(i,j,nlevp) - (plsig(nlev)-plsig(k)) ! log(p) used as vert. coordinate
        do kb = 1, nlevp-1
          if(pl(i,j,kb).gt.zk1) exit   ! kb: index of the input model level immediately below the k level
        enddo
        ztk = (tp(i,j,kb)*(zk1-pl(i,j,kb-1)) + tp(i,j,kb-1)*(pl(i,j,kb) - zk1))/(pl(i,j,kb) - pl(i,j,kb-1))
        zdo = omega*(pl(i,j,kb) - zk1)
        tk(k) = w1*(ztk*(1.+zdo)/(1.-zdo)) + (1.-w1)*tk(k)
      enddo

     endif

     t(i,j,1:nlev) = tk(1:nlev)

     if(lbot < nlev) then

!  Q is extrapolated assuming constant rel. humidity if lapse rate is not large,
!  otherwise values are mainly result of extrapolation with interp_spline_1d.

      p1 = exp(pl(i,j,nlevp))
      call qsat_tetens(tp(i,j,nlevp), p1, eps, qsat, qsw, qsi, es, esw, esi)
      rh1 = min(qp(i,j,nlevp)/qsat, 1.)
      do k = lbot, nlev
        kf = k - nlev + nlevp
        pkf = exp(pl(i,j,kf))
        call qsat_tetens(tp(i,j,kf), pkf, eps, qsat, qsw, qsi, es, esw, esi)
        rh = min(qp(i,j,kf)/qsat, 1.)
        rh = 0.85*rh + 0.15*rh1
        pks = exp(plsig(k))
        call qsat_tetens(t(i,j,k), pks, eps, qsat, qsw, qsi, es, esw, esi)
        q(i,j,k) = w1*rh*qsat + (1.-w1)*q(i,j,k)
        q(i,j,k) = min(q(i,j,k), qsat)
        if(q(i,j,k)/qsat<0.95) qc(i,j,k) = 0.
      enddo
     endif
    endif  ! end of condition on level types

! From here: for all level types

! CONTROLQ is called only in the case in which cloud water/ice is absent in input;
! it checks that Q does not exceed (significantly) saturation and
! increases slightly q values close to saturation

    do k = 1,nlev
      psig = exp(plsig(k))
      if (iflag_cloude == 0 .and. input_model /= 'BOLAM') then
        call controlq (q(i,j,k),t(i,j,k),psig,eps)
      endif
    enddo

! Interpolation of wind
! In case of BOLAM input model, the staggering of the input mod. grid must be taken into account

! Interpolation of u-wind component

    if(i == nlon) then
      psilu = psil(i,j)
    else
      psilu = 0.5*(psil(i,j) + psil(i+1,j))
    endif

    call livelli(psu(i,j),s,nlev,plsig,p0,alfa) ! computes log(p) of output model levels at u-points
    if(input_model == 'BOLAM') call livelli(psilu,sige2,nlevp,pl2,p0e,alfae)

    if(iflag_80m==0.or.k80==0) then
     call near(plsig,nlev,pl2,nlevp,iv)
     call interp_spline_1d(uk,plsig,nlev,upk,pl2,nlevp,iv,alf_uv,ex1_uv,ex2_uv)
    elseif(iflag_80m==1.and.k80 > 0) then
      call near(plsig,nlev,pl2,nlevp+1,iv)
      call interp_spline_1d(uk,plsig,nlev,upk,pl2,nlevp+1,iv,alf_uv,ex1_uv,ex2_uv)
    endif
    u(i,j,1:nlev) = uk(1:nlev)

! Find output model level lbot located immediately above the lowest input model level

     do k = nlev, 1, -1
       if (plsig(k) <= pl2(nlevp)) then
        lbot = k
        exit
       endif
     enddo

! Redefinition of u below the lowest input model level for hybrid coordinate
! input models (for which lbot is defined above).
! The wind of lowest level of the output model is re-defined equal to the lowest level of the input model;
! linear interpolation between lbot and the lowest level is then applied.

    if (level_type == 2) then
      if (lbot < nlev) then
        u(i,j,nlev) = upup(i,j,nlevp)
        do k = lbot+1, nlev-1
        u(i,j,k) = u(i,j,nlev) + (plsig(k)-plsig(nlev))*(u(i,j,lbot)-u(i,j,nlev))/(plsig(lbot)-plsig(nlev))
        enddo
      elseif (lbot == nlev) then
        u(i,j,nlev) = 0.8*u(i,j,nlev) + 0.2*upup(i,j,nlevp) ! should reduce wind at lowest level when interpolated
      endif
    endif

! Interpolation of v-wind component

    if(j == 1) then
      psilv = psil(i,j)
    else
      psilv = 0.5*(psil(i,j-1) + psil(i,j))
    endif

    call livelli(psv(i,j),s,nlev,plsig,p0,alfa) ! computes log(p) of output model levels at v-points
    if(input_model == 'BOLAM') call livelli(psilv,sige2,nlevp,pl2,p0e,alfae)

    if(iflag_80m==0.or.k80==0) then
     call near(plsig,nlev,pl2,nlevp,iv)
     call interp_spline_1d(vk,plsig,nlev,vpk,pl2,nlevp,iv,alf_uv,ex1_uv,ex2_uv)
    elseif(iflag_80m==1.and.k80 > 0) then
      call near(plsig,nlev,pl2,nlevp+1,iv)
      call interp_spline_1d(vk,plsig,nlev,vpk,pl2,nlevp+1,iv,alf_uv,ex1_uv,ex2_uv)
    endif
    v(i,j,1:nlev) = vk(1:nlev)

! Find output model level lbot located immediately above the lowest input model level

     do k = nlev, 1, -1
       if (plsig(k) <= pl2(nlevp)) then
        lbot = k
        exit
       endif
     enddo

! Redefinition of v below the lowest input model level for hybrid coordinate
! input models (for which lbot is defined above).
! The wind of lowest level of the output model is re-defined equal to the lowest level of the input model;
! linear interpolation between lbot and the lowest level is then applied.

    if (level_type == 2 ) then
      if (lbot < nlev) then
        v(i,j,nlev) = vpvp(i,j,nlevp)
        do k = lbot+1, nlev-1
        v(i,j,k) = v(i,j,nlev) + (plsig(k)-plsig(nlev))*(v(i,j,lbot)-v(i,j,nlev))/(plsig(lbot)-plsig(nlev))
        enddo
      elseif (lbot == nlev) then
        v(i,j,nlev) = 0.8*v(i,j,nlev) + 0.2*vpvp(i,j,nlevp) ! should reduce wind at lowest level when interpolated
      endif
    endif

  enddo
  enddo  ! End of vertical interpolations (loop on grid points)

! Comput. of max. wind speed (only for printout)

  wspeed = sqrt(u**2+v**2)
  wm  = maxval(wspeed)
  ijk = maxloc(wspeed)
  write (*,'(a,f8.2,a,3i5)') " Max. wind speed in BOLAM grid:", wm," at grid point",ijk(1),ijk(2),ijk(3)

!---------------------------------------------------------
! Start of definition of fields at the ground surface and in the soil
!---------------------------------------------------------

 if (surf_elaborate) then

  zzh1 = 20. + 70.*(dlat/0.12)**.5
  zzh2 = 2.*zzh1

  if (nlevg_inp < 6) then
    ng = 1
  else
    ng = 2
  endif

  if (input_model /= "BOLAM") then
    zdiffer = 450.
    zdifcoeff = 0.6
  else
    zdiffer = 600.
    zdifcoeff = 1.00
  endif

! Start of the loop on the grid points

  do j=1,nlat
  do i=1,nlon

! Computation of a weighting factor depending on the difference between the
! input model topography interpolated on the model grid (PHIGI2D) and topography
! defined in output grid (HTOPI). DIFTOP > 0 if PHIGI2D > HTOPI.
! In case of "valley" with respect to the input model orography,
! and in winter only, the lapse rate is further reduced

    diftop(i,j) = (phigi2d(i,j)/g0)-htopi(i,j)
    if (diftop(i,j) > 0..and.zcoeday < .5) then
      zlapse = gamma*(.5+zcoeday)
    else
      zlapse = gamma
    endif

! Definition of twater: temperature of surface water (sea, lakes, rivers, including sea ice and lake ice).
! The value is selected between ssti (that is correct for open sea, but not for lakes)
! and temp. at a deep soil level of input model (reasonable for lakes)
! as a function of fmask and topography (by a linear interpolation).
! Twater is not corrected using diftop at an altitude < zzh1 to avoid contamination at sea near the
! coast where the model sea level can be at some height above sea level.

    ztlake = 0.5*(tgii(i,j,nlevg_inp-ng-1)+tgii(i,j,nlevg_inp-ng)) + 0.65*gammac*diftop(i,j)

    if (htopi(i,j) < zzh1) then
      twater(i,j) = ssti(i,j)
    elseif (htopi(i,j) > zzh2) then
      twater(i,j) = ztlake
    else
      twater(i,j) = ssti(i,j)+(htopi(i,j)-zzh1)*(ztlake-ssti(i,j))/(zzh2-zzh1)
    endif

! Re-definition of fice (fraction of sea/lake ice with respect to the sole part of sea/lake water in the grid box)

    if (twater(i,j)<271.4.and.fmask(i,j)>=0.5) fice(i,j) = max(fice(i,j), (271.4-twater(i,j))/9.) ! over sea (and lakes)
    if (twater(i,j)<273.0.and.fmask(i,j)>=0.5.and.htopi(i,j)>zzh1) &
                                                 fice(i,j) = max(fice(i,j),(273.-twater(i,j))/6.)  ! over higher lakes
    fice(i,j) = min(fice(i,j), 1.)
    fice(i,j) = max(fice(i,j), 0.)
    if(fmask(i,j).lt.0.5) fice(i,j) = 0.

    if(fice(i,j).gt.0.01) then
    tice(i,j) = min(tice(i,j), twater(i,j))
    endif

! Definition of iceth in case of IFS data, as a function of fice

    if (input_model == 'IFS') then
    iceth(i,j) = -3.5 + 5.*fice(i,j)     ! sets 0.5<iceth<1.5 for 0.8<fice<1.
    endif

    if(fice(i,j).lt.0.8) then
      iceth(i,j) = 0.                   ! treshold must be consistent with the BOLAM model
    else
      iceth(i,j) = max(iceth(i,j), 0.5) ! treshold and value must be consistent with the BOLAM model
    endif
    iceth(i,j) = min(iceth(i,j), 3.0)

! Snow redefinition (in meters of equivalent water)

    fhard = 1. - fmask(i,j) + fice(i,j)*fmask(i,j)
    snow(i,j) = max(snow(i,j)*fhard, 0.)
    if (fmask(i,j).gt.0.5.and.fice(i,j).lt.0.5) snow(i,j) = 0.
    if (snow(i,j) < 0.0005) snow(i,j) = 0.
    snow(i,j) = min(snow(i,j), 0.08 + htopi(i,j)*(0.2-0.08)/3.4e3) ! 80 cm at sea level, 2 m at 3400 m
    fsnow = min((snow(i,j)/.020)**.67, max(0.9, fice(i,j))) ! over compact sea ice, fsnow can reach 1.

! Definition of surface temperature:
! over land, the value is intermediate between TSURF from the input model and the
! temperature extrapolated from the lower atmosphere TSAIR, using a weighting factor WF,
! depending on difference between the topography of the input and output models.
! TSAIR is abs. temp., defined from the virt. temp. tvs using the humidity of the lowest layer

    wf = zdifcoeff*zdiffer/(abs(diftop(i,j))+zdiffer)
    tsair(i,j) = tvs(i,j)/(1.+ep*q(i,j,nlev))

! TLAND is defined as as weighted average between TSAIR and TSKINI.
! Check that TSURF is present in the input data. If TSURF not available, TLAND=TSAIR

    if (tsurfi(1,1) > 190..and.tsurfi(1,1) < 340.) then
      tland = wf*(tsurfi(i,j)+zlapse*diftop(i,j)) + (1.-wf)*tsair(i,j)
    else
      if(i == 1.and.j == 1) print*, "Caution: tsurf is not available in input data!"
      tland = tsair(i,j)
    endif

! Definition of TSURF

    if(fmask(i,j) > .5) then
      tsurf(i,j) = twater(i,j)
    else
      tsurf(i,j) = tland
    endif

! Case of presence of snow over land or sea/lake ice

    if(fsnow > 0.05) then
    tsurf(i,j) = min(tsurf(i,j), (1.-fsnow)*tsurf(i,j) + fsnow*t0)
    endif

! Over glacier

    if (int(soil_map(i,j,1)) == nst-1) tsurf(i,j) = min(tsurf(i,j), 273.)

! Definition of TG.
! For definition of orographic corrections of soil temperature at various
! levels, different lapse rates are applied, with seasonal variations decreasing
! from near the surface to deeper levels.
! In case of Bolam input model (self-nesting), the difference between surface
! sea temp. and deep sea temp. is preserved.

    if (fmask(i,j) > 0.5) then        ! water surface
      tg(i,j,1) = tsurf(i,j)
      if (phigi2d(i,j)/g0 < zzh1) then
        tg(i,j,2:nlevg) = sstcli(i,j) ! water surface is sea surface
      else
        tg(i,j,2:nlevg) = tsurf(i,j)  ! water surface is lakes surface
      endif
    else                              ! land surface
      tg(i,j,1) = tgi(i,j,1) + .5*(zlapse + gammac)*diftop(i,j)
      tg(i,j,2) = tgi(i,j,2) + 1./3.*(zlapse + gammac + 6.e-3)*diftop(i,j)
      tg(i,j,3) = tgi(i,j,3) + .5*(gammac + 6.e-3)*diftop(i,j)
      tg(i,j,4:nlevg) = tgi(i,j,4:nlevg) + 6.e-3*diftop(i,j)
    endif

! Redefinition of TG in case of sea/lake ice
! (important: used in the model as reference temp.)

    tg(i,j,:) = tg(i,j,:)*(1.-fice(i,j)) + tice(i,j)*fice(i,j)

! Redefinition of TG (deep levels) in case of sea/lake ice
! Case of thin ice: deep soil TG used as relaxation temp.

    if(fice(i,j) > 0.1 .and. fice(i,j) < 0.8 .and. fmask(i,j) > 0.5 .and. htopi(i,j) < zzh1) then
     tg(i,j,2:nlevg) = min(tg(i,j,2:nlevg), (1.-fice(i,j))*271.4 + fice(i,j)*tg(i,j,2:nlevg))
     tg(i,j,2:nlevg) = max(tg(i,j,2:nlevg), (1.-fice(i,j))*271.4 + fice(i,j)*247.)
    elseif(fice(i,j) > 0.1 .and. fmask(i,j) > 0.5 .and. htopi(i,j) >= zzh1) then
     tg(i,j,2:nlevg) = min(tg(i,j,2:nlevg), (1.-fice(i,j))*273.1 + fice(i,j)*tg(i,j,2:nlevg))
     tg(i,j,2:nlevg) = max(tg(i,j,2:nlevg), (1.-fice(i,j))*273.1 + fice(i,j)*247.)
    endif

! Case of thick ice: all soil TG defined as in soil/glacier scheme

    if(fice(i,j) >= 0.8 .and. fmask(i,j) > 0.5 .and. htopi(i,j) < zzh1) then
     do jklev = 2, nlevg
     tg(i,j,jklev) = tg(i,j,1) + (soil_lev(jklev)-soil_lev(1))*(271.4-tg(i,j,1))/(soil_lev(nlevg)-soil_lev(1))
     enddo
    elseif(fice(i,j) >= 0.8 .and. fmask(i,j) > 0.5 .and. htopi(i,j) >= zzh1) then
     do jklev = 2, nlevg
     tg(i,j,jklev) = tg(i,j,1) + (soil_lev(jklev)-soil_lev(1))*(273.1-tg(i,j,1))/(soil_lev(nlevg)-soil_lev(1))
     enddo
    endif

! Glaciers and thick sea ice (temporary sea glacier): 
! T must remain below 0°C, snow must be defind, also all soil parameters

    if (fice(i,j) >= 0.8.or.int(soil_map(i,j,1)) == nst-1) then

      snow(i,j) = max (0.05, snow(i,j))

      if (fice(i,j) >= 0.8) then ! sea ice
        soil_psi(i,j,:) = 2.2 ! thermal conductivity of ice
        do jklev = 1, nlevg
          if (soil_lev(jklev) > iceth(i,j)) exit
        enddo
        ind_lev_soil_h_bottom(i,j) = max(jklev-1, 3)
        ind_lev_soil_w_bottom(i,j) = 0
        soil_albedo_dry(i,j) = 0.70
        soil_albedo_wet(i,j) = 0.70
        soil_emiss1_dry(i,j) = 0.98
        soil_emiss2_dry(i,j) = 0.98
        soil_emiss1_wet(i,j) = 0.98
        soil_emiss2_wet(i,j) = 0.98
        albedo(i,j) = 0.7
        tsurf(i,j) = min(tsurf(i,j), 273.0)
        do jklev = 1, ind_lev_soil_h_bottom(i,j)
          tg(i,j,jklev) = min(tg(i,j,jklev), 273.0)
        enddo
      else ! land glacier
        tsurf(i,j) = min(tsurf(i,j), 273.1)
        tg(i,j,1) = min(tg(i,j,1), 272.5)
        tg(i,j,2) = min(tg(i,j,2), 272.0)
        tg(i,j,3) = min(tg(i,j,3), 271.5)
        tg(i,j,4:nlevg) = min(tg(i,j,4:nlevg), 271.0)
      endif

    endif

! Further snow reduction over land as a function of TG(1) and TSAIR
! and redefinition of TSURF

    ttme = 0.8*tg(i,j,1) + 0.2*tsair(i,j)
    if (ttme.gt.278..and.fmask(i,j).le.0.5) then
    snow(i,j) = 0.
    tsurf(i,j) = tland
    endif

! Defintion of soil temperature at the soil top

    if (snow(i,j) < 0.01) then
      tgsurf(i,j)=tsurf(i,j)
    else
      tgsurf(i,j)=min(tg(i,j,1), 273.) 
    endif

! Definition of ice fraction in soil water content

    fice_soil(i,j,:)=0.
    if (fmask(i,j) < 0.5.or.fice(i,j) >= 0.8) then
      if (fice(i,j) >= 0.8.or.int(soil_map(i,j,1)) == nst-1) then ! Glaciers and thick sea ice (sea glacier)
        fice_soil(i,j,:)=1.
        fice_soil_surf(i,j)=1.
      else ! Soil
        call frac_soil_ice_point(tgsurf(i,j), fice_soil_surf(i,j), soil_par_b(i,j,1))
        do jklev = 1, ind_lev_soil_h_bottom(i,j)
          call frac_soil_ice_point(tg(i,j,jklev), fice_soil(i,j,jklev), soil_par_b(i,j,jklev))
        enddo
      endif
    endif

! Conversion of interpolated fields of relative soil water content
! into volumetric soil content (m**3/m**3)
! (soil water is extrapolated from land to sea across coastlines)

    qg(i,j,1:nlevg) = qgi(i,j,1:nlevg)*(soil_qmax(i,j,1:nlevg)-soil_qmin(i,j,1:nlevg)) + soil_qmin(i,j,1:nlevg)

    qgsurf(i,j) = qg(i,j,1)

! Denition of snow variables at snow levels
! 1-st level is snow top, bottom level is soil surface

    snow(i,j)=snow(i,j)*1.e3
    call snow_lev_def(snow(i,j), nlsnow, snow_lev(i,j,:))
    snow(i,j)=snow(i,j)*1.e-3

! Definition of snow variables

    snow_t(i,j,:) = val_missing
    snow_fice(i,j,:) = val_missing
    snow_age(i,j,:) = val_missing
    snow_melt_age(i,j,:) = val_missing
    snow_dens(i,j,:) = val_missing
  
    if (nlsnow > 0) then
      snow_t(i,j,1) = min(tsurf(i,j), 273.)
      if (nlsnow < 2) then
        snow_t(i,j,1) = min(tg(i,j,1), 273.)
      else
        snow_t(i,j,nlsnow+1) = min(tg(i,j,1), 273.)
        z1=1./float(nlsnow+1)
        do jklev=2,nlsnow
          snow_t(i,j,jklev) = snow_t(i,j,1)*(1.-z1*float(jklev-1)) + snow_t(i,j,nlsnow+1)*z1*float(jklev-1)
        enddo
      endif
      snow_fice(i,j,1:nlsnow+1) = 1.
      snow_age(i,j,1:nlsnow+1) = 10. ! days
      snow_melt_age(i,j,1:nlsnow+1) = 0. ! days
      if (fmask(i,j) < 0.5.and.int(soil_map(i,j,1)) == nst-1) then ! glacier
        snow_age(i,j,1:nlsnow+1) = 180. ! days
        snow_melt_age(i,j,1:nlsnow+1) = 30. ! days
      endif
      if (fmask(i,j) >= 0.5.and.fice(i,j) >= 0.8) then ! thick sea ice
        snow_age(i,j,1:nlsnow+1) = 180. ! days
        snow_melt_age(i,j,1:nlsnow+1) = 0. ! days
      endif
    endif

!  Definition of roughness length (over sea it will be redefined in the model)

    rgm(i,j) = max( (1.-fmask(i,j))*veg_roughness(i,j) + 1.e-3*(fmask(i,j)-fice(i,j)) + 5.e-2*fice(i,j), 1.e-3)
    rgq(i,j) = rgm(i,j)

! Soil water content definition for water bodies: sea, lakes, rivers, glaciers

    if (fmask(i,j) > 0.5.or.nint(soil_map(i,j,1)) >= nst-1) qg(i,j,:) = 1.

! Other 2D fields included in the MHF (it is not necessary to initialise them here)

    cloud (i,j) = 0.
    totpre(i,j) = 0.
    conpre(i,j) = 0.
    snfall(i,j) = 0.
    runoff(i,j) = 0.

! Final redefinition of qvsurf

    if (input_model /= "BOLAM") then
      qvsurf (i,j) = q(i,j,nlev)
    else
      qvsurf (i,j) = q(i,j,nlev)*0.5 + qvsurf(i,j)*0.5
    endif

    cswfl (i,j) = 0.
    clwfl (i,j) = 0.
    chflux(i,j) = 0.
    cqflux(i,j) = 0.
    t2min (i,j) = 0.
    t2max (i,j) = 0.

  enddo
  enddo  ! End of the loop on grid points

!  call plotout(twater-273.15,nlon,nlat,99)

  call ccloud ! computation of total clouds (in "cloud" - used for Po Valley corrections)

  goto 567  ! no corrections for the Po Valley

!--------------------------------------------------------------------------------------------
! For case of NOAA-GFS input data only:
! Increase of humidity and cloud water at low levels over the Po Valley to simulate fog
! Corrections are applied from Nov. to mid Feb. and decrease with cloud cover

    zday1 =  45. ! day of the year for the end of applic. of q correction (mid of Feb.)
    zday2 = 304. ! day of the year for the beginning of applic. of q correction (beginning of Nov.)
    icc = 0  ! flag to avoid printout if the Po Valley area is not in the domain
    if (input_model == 'GFS'.and.(zday < zday1.or.zday > zday2)) then
    call ccloud
    do j = 1,nlat
    do i = 1,nlon
     zlon = xxtg(i,j)
     zlat = yytg(i,j)
     if(zlat > 44.35.and.zlat < 45.75.and.zlon > 7.9.and.zlon < 12.6.and.fmask(i,j) < 0.5) then
     icc = 1
     call livelli(pst(i,j),s,nlev,plsig,p0,alfa)
     do jklev = nlev/2, nlev
      zp = exp(plsig(jklev))
      zaltrel = (ps(i,j)-zp)/11.9      ! approx. height of model levels above ground
      zaltabs = phig(i,j)/g0 + zaltrel ! approx. height of model levels above sea level
      if(zaltabs < 200. .and. zaltrel < 150.) then
      call qsat_tetens(t(i,j,jklev), zp, eps, zqs, zqsw, zqsi, zes, zesw, zesi)
      q(i,j,jklev)  = max(0.91*zqs*(1.-cloud(i,j)**2), q(i,j,jklev)) ! specific humidity set closer to saturation (mixed)
      qc(i,j,jklev) = max(0.1e-3*(1.-cloud(i,j)**2),  qc(i,j,jklev)) ! force cloud water (fog)
      endif
     enddo
     qvsurf(i,j) = q(i,j,nlev)
     endif
    enddo
    enddo
    if(icc==1) print*, "Humidity and cloud water increased near surface over the Po Valley"
    endif

 567 continue

!--------------------------------------------------------------------------------------------

! Using of forecast soil/snow data in place of analysis data 

  call read_forecast_mhf_soil(iflag_soil_snow_frc)
!  call read_forecast_mhf_ascii_soil(iflag_soil_snow_frc)

  write (*,*)
  if (iflag_soil_snow_frc == 1) then

    do j=1,nlat
    do i=1,nlon

      tsurf_copy(i,j)=tsurf(i,j)
      tg_copy(i,j,1:nlevg)=tg(i,j,1:nlevg)

      if (fmask(i,j) < 0.5 ) then ! land or temperary sea glacier

!if (i==278.and.j==209) then
!print *,' Bologna'
!print *,'tsurf ',tsurf(i,j),tsurf_frc(i,j)
!print *,'tgsurf ',tgsurf(i,j),tgsurf_frc(i,j)
!do jklev=1,nlevg
!print *,jklev,tg(i,j,jklev),tg_frc(i,j,jklev)
!enddo
!print *,'qvsurf ',qvsurf(i,j),qvsurf_frc(i,j)
!print *,'qgsurf ',qgsurf(i,j),qgsurf_frc(i,j)
!do jklev=1,nlevg
!print *,jklev,qg(i,j,jklev),qg_frc(i,j,jklev)
!enddo
!print *,' fice_soil_surf ',fice_soil_surf(i,j),fice_soil_surf_frc(i,j)
!do jklev=1,nlevg
!print *,jklev,fice_soil(i,j,jklev),fice_soil_frc(i,j,jklev)
!enddo
!print *,'snow ',snow(i,j),snow_frc(i,j)
!do jklev=1,nlevsnow
!print *,jklev,snow_lev(i,j,jklev),snow_lev_frc(i,j,jklev)
!enddo
!do jklev=1,nlevsnow
!print *,jklev,snow_t(i,j,jklev),snow_t_frc(i,j,jklev)
!enddo
!do jklev=1,nlevsnow
!print *,jklev,snow_fice(i,j,jklev),snow_fice_frc(i,j,jklev)
!enddo
!do jklev=1,nlevsnow
!print *,jklev,snow_age(i,j,jklev),snow_age_frc(i,j,jklev)
!enddo
!do jklev=1,nlevsnow
!print *,jklev,snow_melt_age(i,j,jklev),snow_melt_age_frc(i,j,jklev)
!enddo
!do jklev=1,nlevsnow
!print *,jklev,snow_dens(i,j,jklev),snow_dens_frc(i,j,jklev)
!enddo
!endif
!
!if (i==286.and.j==299) then
!print *,' Leipzig'
!print *,'tsurf ',tsurf(i,j),tsurf_frc(i,j)
!print *,'tgsurf ',tgsurf(i,j),tgsurf_frc(i,j)
!do jklev=1,nlevg
!print *,jklev,tg(i,j,jklev),tg_frc(i,j,jklev)
!enddo
!print *,'qvsurf ',qvsurf(i,j),qvsurf_frc(i,j)
!print *,'qgsurf ',qgsurf(i,j),qgsurf_frc(i,j)
!do jklev=1,nlevg
!print *,jklev,qg(i,j,jklev),qg_frc(i,j,jklev)
!enddo
!print *,' fice_soil_surf ',fice_soil_surf(i,j),fice_soil_surf_frc(i,j)
!do jklev=1,nlevg
!print *,jklev,fice_soil(i,j,jklev),fice_soil_frc(i,j,jklev)
!enddo
!print *,'snow ',snow(i,j),snow_frc(i,j)
!do jklev=1,nlevsnow
!print *,jklev,snow_lev(i,j,jklev),snow_lev_frc(i,j,jklev)
!enddo
!do jklev=1,nlevsnow
!print *,jklev,snow_t(i,j,jklev),snow_t_frc(i,j,jklev)
!enddo
!do jklev=1,nlevsnow
!print *,jklev,snow_fice(i,j,jklev),snow_fice_frc(i,j,jklev)
!enddo
!do jklev=1,nlevsnow
!print *,jklev,snow_age(i,j,jklev),snow_age_frc(i,j,jklev)
!enddo
!do jklev=1,nlevsnow
!print *,jklev,snow_melt_age(i,j,jklev),snow_melt_age_frc(i,j,jklev)
!enddo
!do jklev=1,nlevsnow
!print *,jklev,snow_dens(i,j,jklev),snow_dens_frc(i,j,jklev)
!enddo
!endif
!
!if (i==244.and.j==228) then
!print *,' Plato Rosa'
!print *,'tsurf ',tsurf(i,j),tsurf_frc(i,j)
!print *,'tgsurf ',tgsurf(i,j),tgsurf_frc(i,j)
!do jklev=1,nlevg
!print *,jklev,tg(i,j,jklev),tg_frc(i,j,jklev)
!enddo
!print *,'qvsurf ',qvsurf(i,j),qvsurf_frc(i,j)
!print *,'qgsurf ',qgsurf(i,j),qgsurf_frc(i,j)
!do jklev=1,nlevg
!print *,jklev,qg(i,j,jklev),qg_frc(i,j,jklev)
!enddo
!print *,' fice_soil_surf ',fice_soil_surf(i,j),fice_soil_surf_frc(i,j)
!do jklev=1,nlevg
!print *,jklev,fice_soil(i,j,jklev),fice_soil_frc(i,j,jklev)
!enddo
!print *,'snow ',snow(i,j),snow_frc(i,j)
!do jklev=1,nlevsnow
!print *,jklev,snow_lev(i,j,jklev),snow_lev_frc(i,j,jklev)
!enddo
!do jklev=1,nlevsnow
!print *,jklev,snow_t(i,j,jklev),snow_t_frc(i,j,jklev)
!enddo
!do jklev=1,nlevsnow
!print *,jklev,snow_fice(i,j,jklev),snow_fice_frc(i,j,jklev)
!enddo
!do jklev=1,nlevsnow
!print *,jklev,snow_age(i,j,jklev),snow_age_frc(i,j,jklev)
!enddo
!do jklev=1,nlevsnow
!print *,jklev,snow_melt_age(i,j,jklev),snow_melt_age_frc(i,j,jklev)
!enddo
!do jklev=1,nlevsnow
!print *,jklev,snow_dens(i,j,jklev),snow_dens_frc(i,j,jklev)
!enddo
!endif
!
!if (i==237.and.j==227) then
!print *,' Mont Blanc'
!print *,'tsurf ',tsurf(i,j),tsurf_frc(i,j)
!print *,'tgsurf ',tgsurf(i,j),tgsurf_frc(i,j)
!do jklev=1,nlevg
!print *,jklev,tg(i,j,jklev),tg_frc(i,j,jklev)
!enddo
!print *,'qvsurf ',qvsurf(i,j),qvsurf_frc(i,j)
!print *,'qgsurf ',qgsurf(i,j),qgsurf_frc(i,j)
!do jklev=1,nlevg
!print *,jklev,qg(i,j,jklev),qg_frc(i,j,jklev)
!enddo
!print *,' fice_soil_surf ',fice_soil_surf(i,j),fice_soil_surf_frc(i,j)
!do jklev=1,nlevg
!print *,jklev,fice_soil(i,j,jklev),fice_soil_frc(i,j,jklev)
!enddo
!print *,'snow ',snow(i,j),snow_frc(i,j)
!do jklev=1,nlevsnow
!print *,jklev,snow_lev(i,j,jklev),snow_lev_frc(i,j,jklev)
!enddo
!do jklev=1,nlevsnow
!print *,jklev,snow_t(i,j,jklev),snow_t_frc(i,j,jklev)
!enddo
!do jklev=1,nlevsnow
!print *,jklev,snow_fice(i,j,jklev),snow_fice_frc(i,j,jklev)
!enddo
!do jklev=1,nlevsnow
!print *,jklev,snow_age(i,j,jklev),snow_age_frc(i,j,jklev)
!enddo
!do jklev=1,nlevsnow
!print *,jklev,snow_melt_age(i,j,jklev),snow_melt_age_frc(i,j,jklev)
!enddo
!do jklev=1,nlevsnow
!print *,jklev,snow_dens(i,j,jklev),snow_dens_frc(i,j,jklev)
!enddo
!endif

        tsurf(i,j) = tsurf_frc(i,j)
        tgsurf(i,j) = tgsurf_frc(i,j)
        tg(i,j,1:nlevg) = tg_frc(i,j,1:nlevg)
        qvsurf(i,j) = qvsurf_frc(i,j)
        qgsurf(i,j) = qgsurf_frc(i,j)
        qg(i,j,1:nlevg) = qg_frc(i,j,1:nlevg)
        fice_soil_surf(i,j) = fice_soil_surf_frc(i,j)
        fice_soil(i,j,1:nlevg) = fice_soil_frc(i,j,1:nlevg)
        snow_lev(i,j,1:nlevsnow) = val_missing
        snow_t(i,j,1:nlevsnow) = val_missing
        snow_fice(i,j,1:nlevsnow) = val_missing
        snow_age(i,j,1:nlevsnow) = val_missing
        snow_melt_age(i,j,1:nlevsnow) = val_missing
        snow_dens(i,j,1:nlevsnow) = val_missing
        do k=1,nlevsnow
          if (int(snow_t_frc(i,j,k)) == int(val_missing)) exit
        enddo
        nk=k-1
        if (nk >= 1) then
          snow_lev(i,j,1:nk) = snow_lev_frc(i,j,1:nk)
          snow_t(i,j,1:nk) = snow_t_frc(i,j,1:nk)
          snow_fice(i,j,1:nk) = snow_fice_frc(i,j,1:nk)
          snow_age(i,j,1:nk) = snow_age_frc(i,j,1:nk)
          snow_melt_age(i,j,1:nk) = snow_melt_age_frc(i,j,1:nk)
          snow_dens(i,j,1:nk) = snow_dens_frc(i,j,1:nk)
          snow(i,j) = snow_lev(i,j,nk)*1.e-3
        else
          snow(i,j) = 0.
        endif

      endif ! land

! Case of change grid point status:
! "sea glacier" turn to sea point status, or
! sea icea becames "sea glacier", transfers to land point status

      if (fmask(i,j) >= 0.5) then ! sea

        if (fice_frc(i,j) < 0.8.and.fice(i,j) >= 0.8) then ! sea point became temporary sea glacier

          soil_psi(i,j,:) = 2.2 ! thermal conductivity of ice
          soil_albedo_dry(i,j) = 0.70
          soil_albedo_wet(i,j) = 0.70
          soil_emiss1_dry(i,j) = 0.98
          soil_emiss2_dry(i,j) = 0.98
          soil_emiss1_wet(i,j) = 0.98
          soil_emiss2_wet(i,j) = 0.98

          do jklev = 1, nlevg
            if (soil_lev(jklev) > iceth(i,j)) exit
          enddo
          ind_lev_soil_h_bottom(i,j) = max(jklev-1, 3)
          ind_lev_soil_w_bottom(i,j) = 0

          tsurf(i,j) = min(tsurf(i,j), 273.0)
          do jklev = 1, ind_lev_soil_h_bottom(i,j)
            tg(i,j,jklev) = min(tg(i,j,jklev), 273.0)
          enddo

          fice_soil(i,j,:) = 1.
          fice_soil_surf(i,j) = 1.

          snow(i,j) = max (0.05, snow(i,j))
          snow(i,j)=snow(i,j)*1.e3
          call snow_lev_def(snow(i,j), nlsnow, snow_lev(i,j,:))
          snow(i,j)=snow(i,j)*1.e-3
          snow_t(i,j,:) = val_missing
          snow_fice(i,j,:) = val_missing
          snow_age(i,j,:) = val_missing
          snow_melt_age(i,j,:) = val_missing
          snow_dens(i,j,:) = val_missing
          if (nlsnow > 0) then
            snow_t(i,j,1) = min(tsurf(i,j), 273.)
            if (nlsnow < 2) then
              snow_t(i,j,1) = min(tg(i,j,1), 273.)
            else
              snow_t(i,j,nlsnow+1) = min(tg(i,j,1), 273.)
              z1=1./float(nlsnow+1)
              do jklev=2,nlsnow
                snow_t(i,j,jklev) = snow_t(i,j,1)*(1.-z1*float(jklev-1)) + snow_t(i,j,nlsnow+1)*z1*float(jklev-1)
              enddo
            endif
            snow_fice(i,j,1:nlsnow+1) = 1.
            snow_age(i,j,1:nlsnow+1) = 180. ! days
            snow_melt_age(i,j,1:nlsnow+1) = 0. ! days
          endif

        endif ! sea point became temporary sea glacier

        if (fice_frc(i,j) >= 0.8.and.fice(i,j) < 0.8) then ! temporary sea glacier turn to sea point status

          soil_psi(i,j,:) = 0.57 ! thermal conductivity of water
          soil_albedo_dry(i,j) = 0.07
          soil_albedo_wet(i,j) = 0.07
          soil_emiss1_dry(i,j) = 0.97
          soil_emiss2_dry(i,j) = 0.97
          soil_emiss1_wet(i,j) = 0.97
          soil_emiss2_wet(i,j) = 0.97

          ind_lev_soil_h_bottom(i,j) = 0
          ind_lev_soil_w_bottom(i,j) = 0

          tsurf(i,j) = tsurf_copy(i,j)
          tg(i,j,1:nlevg) = tg(i,j,1:nlevg)

          fice_soil(i,j,:) = 0.
          fice_soil_surf(i,j) = 0.

          snow(i,j) = 0.
          snow_lev(i,j,:) = val_missing
          snow_t(i,j,:) = val_missing
          snow_fice(i,j,:) = val_missing
          snow_age(i,j,:) = val_missing
          snow_melt_age(i,j,:) = val_missing
          snow_dens(i,j,:) = val_missing

        endif ! temporary sea glacier turn to sea point status

      endif ! sea

    enddo
    enddo
    
    write (*,*) "Soil and snow prognostic parameters (temperature, water content, ect.) "
    write (*,*) "are initialisated with forecast data generated by a preceded model run"

  else ! iflag_soil_snow_frc

    write (*,*) "Soil and snow prognostic parameters (temperature, water content, ect.) are initialisated with analysis data"
    write (*,*) "    Cold soil and snow start" 

  endif ! iflag_soil_snow_frc
  write (*,*)

! Definition of the radiation parameters at the surface (not needed at all instants)

   call surfradpar(fmask, fice, qg(:,:,1), veg_map, soil_qmax(:,:,1), soil_qmin(:,:,1), &
 soil_albedo_dry, soil_albedo_wet, soil_emiss1_dry, soil_emiss1_wet, soil_emiss2_dry, soil_emiss2_wet, &
 veg_frac, veg_albedo, veg_emiss1, veg_emiss2, &
 nlon, nlat, nlevg, nvt, 1, albedo, emismap1, emismap2)

!--------------------------------------------------------------------------------------------

! Computation of average soil temperature and water content (only for diagnostics)

  tsoilav(:) = 0.
  qsoilav(:) = 0.

  icnt = 0
  do j=1,nlat
  do i=1,nlon
  if(fmask(i,j).lt.0.5) then
    icnt = icnt+1
    tsoilav(1:nlevg) = tsoilav(1:nlevg) + tg(i,j,1:nlevg)
    qsoilav(1:nlevg) = qsoilav(1:nlevg) + qg(i,j,1:nlevg)
  endif
  enddo
  enddo

  if (icnt > 0) then
    tsoilav(1:nlevg) = tsoilav(1:nlevg)/float(icnt)
    qsoilav(1:nlevg) = qsoilav(1:nlevg)/float(icnt)
    print*
    print*, "Averages of soil temperature (deg.C) at various soil levels, BOLAM grid:"
    print '(30f9.2)',tsoilav(1:nlevg)-273.15
    print*
    print*, "Averages of soil water content (m**3/m**3) at various soil levels, BOLAM grid:"
    print '(30f9.4)',qsoilav(1:nlevg)
  endif

 endif ! surf_elaborate

! Correction of the wind at the lowest level as a function of roughness

 do j = 1,nlat
  do i = 1,nlon
    zrgm = min(rgm(i,j),8.)
    fatt =  1./(1.+0.12*zrgm)
    u(i,j,nlev) = u(i,j,nlev)*fatt
    v(i,j,nlev) = v(i,j,nlev)*fatt
  enddo
 enddo

! -------------------
!  Plots for check: fields on the output model grid
! -------------------

! if (ifile == 1) then

!    do k=1,nlev
!      call plotout (t(1,1,k),nlon,nlat,99)
!      call plotout (q(1,1,k),nlon,nlat,99)
!      call plotout (u(1,1,k),nlon,nlat,99)
!      call plotout (v(1,1,k),nlon,nlat,99)
!    enddo

!    call plotout (t(1,1,1),nlon,nlat,99)
!    call plotout (q(1,1,1),nlon,nlat,99)
!    call plotout (u(1,1,1),nlon,nlat,99)
!    call plotout (v(1,1,1),nlon,nlat,99)
!    call plotout (t(1,1,nlev),nlon,nlat,99)
!    call plotout (q(1,1,nlev),nlon,nlat,99)
!    call plotout (u(1,1,nlev),nlon,nlat,99)
!    call plotout (v(1,1,nlev),nlon,nlat,99)

!    call plotout (htopi,nlon,nlat,99)
!    call plotout (phigi2d/g0,nlon,nlat,99)
!    call plotout (htopi-phigi2d/g0,nlon,nlat,99)
!    call plotout (fmask,nlon,nlat,99)
!    call plotout (albedo,nlon,nlat,99)
!    call plotout (emismap1,nlon,nlat,99)
!    call plotout (emismap2,nlon,nlat,99)
!    call plotout (tsair,nlon,nlat,99)
!    call plotout (tsurf,nlon,nlat,99)
!    call plotout (tg(1,1,1),nlon,nlat,99)
!    call plotout (tg(1,1,2),nlon,nlat,99)
!    call plotout (tg(1,1,3),nlon,nlat,99)
!    call plotout (tg(1,1,4),nlon,nlat,99)
!    call plotout (tg(1,1,5),nlon,nlat,99)
!    call plotout (tg(1,1,6),nlon,nlat,99)
!    call plotout (qg(1,1,1),nlon,nlat,99)
!    call plotout (qg(1,1,2),nlon,nlat,99)
!    call plotout (qg(1,1,3),nlon,nlat,99)
!    call plotout (qg(1,1,4),nlon,nlat,99)
!    call plotout (qg(1,1,5),nlon,nlat,99)
!    call plotout (qg(1,1,6),nlon,nlat,99)
!    call plotout (rgm,nlon,nlat,99)
!    call plotout (rgq,nlon,nlat,99)
!    call plotout (snow,nlon,nlat,99)
!    call plotout (fice,nlon,nlat,99)
!    call plotout (diftop,nlon,nlat,99)
!    stop
! endif

! Writing of processed model fields in MHF format:

  if (frame) then

    do j = 1,nlat
    do i = 1,nlon
      if (mask_frame_out(i,j) == 0) then
        ps(i,j) = val_missing
      endif
    enddo
    enddo
    do jklev = 1,nlev
      do j = 1,nlat
       do i = 1,nlon
        if (mask_frame_out(i,j) == 0) then
          u (i,j,jklev) = val_missing
          v (i,j,jklev) = val_missing
          t (i,j,jklev) = val_missing
          q (i,j,jklev) = val_missing
          qc(i,j,jklev) = val_missing
        endif
       enddo
      enddo
    enddo

  endif

!i=ipoint; j=jpoint
!write (31,*) ' ------------------------'
!write (31,*) ' p '
!write (31,*) p(i,j,1:nlev)
!write (31,*) ' u '
!write (31,*) u(i,j,1:nlev)
!write (31,*) ' v '
!write (31,*) v(i,j,1:nlev)
!write (31,*) ' w '
!write (31,*) w(i,j,1:nlev+1)
!write (31,*) ' t '
!write (31,*) t(i,j,1:nlev)
!write (31,*) ' q '
!write (31,*) q(i,j,1:nlev)
!write (31,*) ' qcw '
!write (31,*) qcw(i,j,1:nlev)
!write (31,*) ' qci '
!write (31,*) qci(i,j,1:nlev)
!write (31,*) ' veg_lai '
!write (31,*) veg_lai(i,j)
!write (31,*) ' veg_frac '
!write (31,*) veg_frac(i,j)
!write (31,*) ' rgm '
!write (31,*) rgm(i,j)
!write (31,*) ' rgq '
!write (31,*) rgq(i,j)
!write (31,*) ' iceth '
!write (31,*) iceth(i,j)
!write (31,*) ' fice '
!write (31,*) fice(i,j)
!write (31,*) ' albedo '
!write (31,*) albedo(i,j)
!write (31,*) ' emismap1 '
!write (31,*) emismap1(i,j)
!write (31,*) ' emismap2 '
!write (31,*) emismap2(i,j)
!write (31,*) ' tsurf '
!write (31,*) tsurf(i,j)
!write (31,*) ' tgsurf '
!write (31,*) tgsurf(i,j)
!write (31,*) ' tg '
!write (31,*) tg(i,j,1:nlevg)
!write (31,*) ' qvsurf '
!write (31,*) qvsurf(i,j)
!write (31,*) ' qgsurf '
!write (31,*) qgsurf(i,j)
!write (31,*) ' qg '
!write (31,*) qg(i,j,1:nlevg)
!write (31,*) ' fice_soil_surf '
!write (31,*) fice_soil_surf(i,j)
!write (31,*) ' fice_soil '
!write (31,*) fice_soil(i,j,1:nlevg)
!write (31,*) ' snow '
!write (31,*) snow(i,j)
!write (31,*) ' snow_lev '
!write (31,*) snow_lev(i,j,1:nlevsnow)
!write (31,*) ' snow_t '
!write (31,*) snow_t(i,j,1:nlevsnow)
!write (31,*) ' snow_fice '
!write (31,*) snow_fice(i,j,1:nlevsnow)
!write (31,*) ' snow_age '
!write (31,*) snow_age(i,j,1:nlevsnow)
!write (31,*) ' snow_melt_age '
!write (31,*) snow_melt_age(i,j,1:nlevsnow)
!write (31,*) ' snow_dens '
!write (31,*) snow_dens(i,j,1:nlevsnow)
!write (31,*) ' snow_albedo '
!write (31,*) 0.71
!write (31,*) ' htop'
!write (31,*) phig(i,j)/g0
!write (31,*) ' htopvar'
!write (31,*) htopvar(i,j)
!write (31,*) ' soil_map'
!write (31,*) soil_map(i,j,1:nst+1)
!write (31,*) ' veg_map'
!write (31,*) veg_map(i,j,1:nvt+1)
!write (31,*) ' soil_qmax'
!write (31,*) soil_qmax(i,j,1:nlevg)
!write (31,*) ' soil_qmin'
!write (31,*) soil_qmin(i,j,1:nlevg)
!write (31,*) ' soil_c'
!write (31,*) soil_c(i,j,1:nlevg)
!write (31,*) ' soil_rho'
!write (31,*) soil_rho(i,j,1:nlevg)
!write (31,*) ' soil_psi'
!write (31,*) soil_psi(i,j,1:nlevg)
!write (31,*) ' soil_k'
!write (31,*) soil_k(i,j,1:nlevg)
!write (31,*) ' soil_par_b'
!write (31,*) soil_par_b(i,j,1:nlevg)
!write (31,*) ' soil_par_c'
!write (31,*) soil_par_c(i,j,1:nlevg)
!write (31,*) ' soil_qrel_wilt'
!write (31,*) soil_qrel_wilt(i,j,1:nlevg)
!write (31,*) ' soil_qrel_ref'
!write (31,*) soil_qrel_ref(i,j,1:nlevg)
!write (31,*) ' soil_albedo_dry'
!write (31,*) soil_albedo_dry(i,j)
!write (31,*) ' soil_albedo_wet'
!write (31,*) soil_albedo_wet(i,j)
!write (31,*) ' soil_emiss1_dry'
!write (31,*) soil_emiss1_dry(i,j)
!write (31,*) ' soil_emiss1_wet'
!write (31,*) soil_emiss1_wet(i,j)
!write (31,*) ' soil_emiss2_dry'
!write (31,*) soil_emiss2_dry(i,j)
!write (31,*) ' soil_emiss2_wet'
!write (31,*) soil_emiss2_wet(i,j)
!write (31,*) ' veg_root_depth'
!write (31,*) veg_root_depth(i,j)
!write (31,*) ' veg_roughness'
!write (31,*) veg_roughness(i,j)
!write (31,*) ' veg_albedo'
!write (31,*) veg_albedo(i,j)
!write (31,*) ' veg_emiss1'
!write (31,*) veg_emiss1(i,j)
!write (31,*) ' veg_emiss2'
!write (31,*) veg_emiss2(i,j)
!write (31,*) ' ind_lev_soil_h_bottom'
!write (31,*) ind_lev_soil_h_bottom(i,j)
!write (31,*) ' ind_lev_soil_w_bottom'
!write (31,*) ind_lev_soil_w_bottom(i,j)
!write (31,*) ' veg_lai_max'
!write (31,*) veg_lai_max

  if (ifile == 1.and.flag_constant_fields == 1) then
    call write_param_const(x0d, y0d, alon0, alat0+dlat*0.5, htopi, zhtopvar)
  endif

  call write_mhf_atm(ifile)
  call write_mhf_soil(ifile)

! call outgraph(80102,1,nlon,nlat,x0d,y0d,alon0,alat0+dlat*0.5,dlon,dlat,&
! 2, 0,  0,  1, 1, 0,idate0(1:5),iperiod(1:3),fmask(:,:),1.,0.)
! call outgraph(80102,1,nlon,nlat,x0d,y0d,alon0,alat0+dlat*0.5,dlon,dlat,&
! 0, 0,  0,  1, 1, 0,idate0(1:5),iperiod(1:3),((qg(:,:,11)-soil_qmin(:,:,11))/(soil_qmax(:,:,11)-soil_qmin(:,:,11))),1.,0.)

!---------------------------------------------------------------------------
! DEALLOCATION OF ARRAYS
! (note that suolo_inp, vegeta_inp, soilvegpar_inp are used only in case of BOLAM as input model,
! so they are defined at the first instant only and do not need to be deallocated).

   deallocate(field3d)
   deallocate(field2d)
   deallocate(field3d_soil)
   deallocate(ge)
   deallocate(te)
   deallocate(ue)
   deallocate(ve)
   deallocate(qe)
   deallocate(qce)
   deallocate(qcwe)
   deallocate(qcie)
   deallocate(phle)
   deallocate(phige)
   deallocate(phige2d)
   deallocate(fmaske)

   deallocate(psel)
   deallocate(tsurfe)
   deallocate(qvsurfe)
   deallocate(pmsle)
   deallocate(snowe)
   deallocate(ficee)
   deallocate(ficeee)
   deallocate(icethe)
   deallocate(icethee)
   deallocate(ticee)
   deallocate(ticeee)
   deallocate(xe)
   deallocate(xue)
   deallocate(ye)
   deallocate(yve)
   deallocate(alon_inp)
   deallocate(alat_inp)
   deallocate(mask_frame)

   deallocate(tge)
   deallocate(qge)
   deallocate(cctote)
   deallocate(u10e)
   deallocate(v10e)
   deallocate(t2e)
   deallocate(td2e)
   deallocate(sste)
   deallocate(sstcle)
   deallocate(work)
   deallocate(workf)

   deallocate(t80e)
   deallocate(q80e)
   deallocate(p80e)
   deallocate(u80e)
   deallocate(v80e)

   if (input_model == 'GFS') then
     deallocate(soile)
     deallocate(qgmine2d)
     deallocate(qgmaxe2d)
   endif

!---------------------------------------------------------------------------
! End of loop IFILE on all input files (for initial and boundary conditions)

 enddo

 stop
 end
!=======================================================================
subroutine wrec2 (kunit, nlon, nlat, vect)

implicit none

integer :: kunit, nlon, nlat
real, dimension(nlon,nlat) :: vect

 write(kunit) vect(1:nlon,1:nlat)

return
end subroutine wrec2
!=======================================================================
subroutine wrec2_int (kunit, nlon, nlat, ivect)

implicit none

integer :: kunit, nlon, nlat
integer, dimension(nlon,nlat) :: ivect

 write(kunit) ivect(1:nlon,1:nlat)

return
end subroutine wrec2_int
!=======================================================================
subroutine write_mhf_atm(nf)

! Writes a MHF file of BOLAM with atmospheric variables

use param
use mod_fields

implicit none

real, dimension(nlon,nlat) :: field2d_add
integer :: nf, iunit = 60, i, j, k
character(len=30) :: file_output

 write (file_output,'(a,i2.2,a)') 'input_atm_',nf,'.mhf'
 
 open (iunit,file=file_output,form='unformatted')
 
 write(iunit) nfdr
 write(iunit) pdr

 call wrec2 (iunit, nlon, nlat, ps(1:nlon,1:nlat))

 do k = 1,nlev
   call wrec2 (iunit, nlon, nlat, u(1:nlon,1:nlat,k))
 enddo

 do k = 1,nlev
   call wrec2 (iunit, nlon, nlat, v(1:nlon,1:nlat,k))
 enddo

 do k = 1,nlev
   call wrec2 (iunit, nlon, nlat, t(1:nlon,1:nlat,k))
 enddo

 do k = 1,nlev
   call wrec2 (iunit, nlon, nlat, q(1:nlon,1:nlat,k))
 enddo

 do k = 1,nlev
   call wrec2 (iunit, nlon, nlat, qc(1:nlon,1:nlat,k))
 enddo

 close (iunit)

 write (*,*)
 write (*,*) 'Output mhf file ',trim(file_output),' written' 
 write (*,*)
 
#ifdef oper
!!! call system("sync")
!!! call system("ls -l -L "//file_output)
!!! call system("date")
 open  (iunit, file=trim(file_output)//'.txt', status='unknown')
 write (iunit,'(2a)') trim(file_output),' is full and closed'
 close (iunit)
!!! call system("sync")
#endif

return
end
!=======================================================================
subroutine write_mhf_soil(nf)

! Writes a MHF file of BOLAM with surface, sea, soil variables

use param
use mod_fields

implicit none

real, dimension(nlon,nlat) :: field2d_add, snow_albedo, runoff_tot
integer :: nf, iunit = 60, i, j, k, iwr
character(len=30) :: file_output

 snow_albedo(:,:) = 0.71
 runoff_tot(:,:) = 0.

 write (file_output,'(a,i2.2,a)') 'input_soil_',nf,'.mhf'
 
 open (iunit,file=file_output,form='unformatted')
 
 write(iunit) nfdr
 write(iunit) pdr

 call wrec2 (iunit, nlon, nlat, veg_lai(1:nlon,1:nlat))

 call wrec2 (iunit, nlon, nlat, veg_frac(1:nlon,1:nlat))

 call wrec2 (iunit, nlon, nlat, rgm(1:nlon,1:nlat))

 call wrec2 (iunit, nlon, nlat, rgq(1:nlon,1:nlat))

 call wrec2 (iunit, nlon, nlat, iceth(1:nlon,1:nlat))

 call wrec2 (iunit, nlon, nlat, fice(1:nlon,1:nlat))

 call wrec2 (iunit, nlon, nlat, albedo(1:nlon,1:nlat))

 call wrec2 (iunit, nlon, nlat, emismap1(1:nlon,1:nlat))

 call wrec2 (iunit, nlon, nlat, emismap2(1:nlon,1:nlat))

 call wrec2 (iunit, nlon, nlat, cloud(1:nlon,1:nlat))

 call wrec2 (iunit, nlon, nlat, totpre(1:nlon,1:nlat))

 call wrec2 (iunit, nlon, nlat, conpre(1:nlon,1:nlat))

 call wrec2 (iunit, nlon, nlat, snfall(1:nlon,1:nlat))

 call wrec2 (iunit, nlon, nlat, tsurf(1:nlon,1:nlat))

 call wrec2 (iunit, nlon, nlat, tgsurf(1:nlon,1:nlat))

 do k = 1,nlevg
   call wrec2 (iunit, nlon, nlat, tg(1:nlon,1:nlat,k))
 enddo

 call wrec2 (iunit, nlon, nlat, qvsurf(1:nlon,1:nlat))

 call wrec2 (iunit, nlon, nlat, qgsurf(1:nlon,1:nlat))

 do k = 1,nlevg
   call wrec2 (iunit, nlon, nlat, qg(1:nlon,1:nlat,k))
 enddo

 call wrec2 (iunit, nlon, nlat, fice_soil_surf(1:nlon,1:nlat))

 do k = 1,nlevg
   call wrec2 (iunit, nlon, nlat, fice_soil(1:nlon,1:nlat,k))
 enddo

 call wrec2 (iunit, nlon, nlat, snow(1:nlon,1:nlat))

 do k = 1,nlevsnow
   call wrec2 (iunit, nlon, nlat, snow_lev(1:nlon,1:nlat,k))
 enddo

 do k = 1,nlevsnow
   call wrec2 (iunit, nlon, nlat, snow_t(1:nlon,1:nlat,k))
 enddo

 do k = 1,nlevsnow
   call wrec2 (iunit, nlon, nlat, snow_fice(1:nlon,1:nlat,k))
 enddo

 do k = 1,nlevsnow
   call wrec2 (iunit, nlon, nlat, snow_age(1:nlon,1:nlat,k))
 enddo

 do k = 1,nlevsnow
   call wrec2 (iunit, nlon, nlat, snow_melt_age(1:nlon,1:nlat,k))
 enddo

 do k = 1,nlevsnow
   call wrec2 (iunit, nlon, nlat, snow_dens(1:nlon,1:nlat,k))
 enddo

 call wrec2 (iunit, nlon, nlat, snow_albedo(1:nlon,1:nlat))

 call wrec2 (iunit, nlon, nlat, cswfl(1:nlon,1:nlat))

 call wrec2 (iunit, nlon, nlat, clwfl(1:nlon,1:nlat))

 call wrec2 (iunit, nlon, nlat, chflux(1:nlon,1:nlat))

 call wrec2 (iunit, nlon, nlat, cqflux(1:nlon,1:nlat))

 call wrec2 (iunit, nlon, nlat, t2min(1:nlon,1:nlat))

 call wrec2 (iunit, nlon, nlat, t2max(1:nlon,1:nlat))

 call wrec2 (iunit, nlon, nlat, runoff(1:nlon,1:nlat))

 call wrec2 (iunit, nlon, nlat, runoff_tot(1:nlon,1:nlat))

! Writing of additional 2D fields

 field2d_add(:,:) = 0.
 do iwr = 1,3
   call wrec2 (iunit, nlon, nlat, field2d_add(1:nlon,1:nlat))
 enddo

 close (iunit)

 write (*,*)
 write (*,*) 'Output mhf file ',trim(file_output),' written' 
 write (*,*)
 
#ifdef oper
!!! call system("sync")
!!! call system("ls -l -L "//file_output)
!!! call system("date")
 open  (iunit, file=trim(file_output)//'.txt', status='unknown')
 write (iunit,'(2a)') trim(file_output),' is full and closed'
 close (iunit)
!!! call system("sync")
#endif

return
end
!=======================================================================
subroutine write_param_const(x0, y0, alon0, alat0, htop, htopvar)

! Writes an additional output file with all constant (in time) 
! model physiographical parameters

use param, only : nlon, nlat, nlevg, dlon, dlat, nst, nvt
use mod_fields, only : fmask, rgm, rgq, soil_map, veg_map, &
 soil_qmax, soil_qmin, soil_c, soil_rho, soil_psi, soil_k, soil_par_b, soil_par_c, soil_qrel_wilt, soil_qrel_ref, &
 soil_albedo_dry, soil_albedo_wet, soil_emiss1_dry, soil_emiss1_wet, soil_emiss2_dry, soil_emiss2_wet, &
 water_table_depth, tg_bottom, qg_rel_bottom, qg_rel_surf_approx, &
 veg_lai, veg_frac, veg_root_depth, veg_roughness, veg_albedo, veg_emiss1, veg_emiss2, &
 snow_dirt, ind_lev_soil_h_bottom, ind_lev_soil_w_bottom, veg_lai_max

implicit none
integer :: iunit=60, k
real :: x0, y0, alon0, alat0
real, dimension(nlon,nlat) :: htop, htopvar
character(len=30) :: file_output="model_param_constant.bin"

 open (iunit, file=file_output, form='unformatted', status='unknown')

 write (iunit) nlon, nlat, nlevg, dlon, dlat, x0, y0, alon0, alat0, nst, nvt

 call wrec2 (iunit, nlon, nlat, fmask(1:nlon,1:nlat))

 call wrec2 (iunit, nlon, nlat, htop(1:nlon,1:nlat))

 call wrec2 (iunit, nlon, nlat, htopvar(1:nlon,1:nlat))

 do k = 1,nst+1
   call wrec2 (iunit, nlon, nlat, soil_map(1:nlon,1:nlat,k))
 enddo

 do k = 1,nvt+1
   call wrec2 (iunit, nlon, nlat, veg_map(1:nlon,1:nlat,k))
 enddo

 do k = 1,nlevg
   call wrec2 (iunit, nlon, nlat, soil_qmax(1:nlon,1:nlat,k))
 enddo

 do k = 1,nlevg
   call wrec2 (iunit, nlon, nlat, soil_qmin(1:nlon,1:nlat,k))
 enddo

 do k = 1,nlevg
   call wrec2 (iunit, nlon, nlat, soil_c(1:nlon,1:nlat,k))
 enddo

 do k = 1,nlevg
   call wrec2 (iunit, nlon, nlat, soil_rho(1:nlon,1:nlat,k))
 enddo

 do k = 1,nlevg
   call wrec2 (iunit, nlon, nlat, soil_psi(1:nlon,1:nlat,k))
 enddo

 do k = 1,nlevg
   call wrec2 (iunit, nlon, nlat, soil_k(1:nlon,1:nlat,k))
 enddo

 do k = 1,nlevg
   call wrec2 (iunit, nlon, nlat, soil_par_b(1:nlon,1:nlat,k))
 enddo

 do k = 1,nlevg
   call wrec2 (iunit, nlon, nlat, soil_par_c(1:nlon,1:nlat,k))
 enddo

 do k = 1,nlevg
   call wrec2 (iunit, nlon, nlat, soil_qrel_wilt(1:nlon,1:nlat,k))
 enddo

 do k = 1,nlevg
   call wrec2 (iunit, nlon, nlat, soil_qrel_ref(1:nlon,1:nlat,k))
 enddo

 call wrec2 (iunit, nlon, nlat, soil_albedo_dry(1:nlon,1:nlat))

 call wrec2 (iunit, nlon, nlat, soil_albedo_wet(1:nlon,1:nlat))

 call wrec2 (iunit, nlon, nlat, soil_emiss1_dry(1:nlon,1:nlat))

 call wrec2 (iunit, nlon, nlat, soil_emiss1_wet(1:nlon,1:nlat))

 call wrec2 (iunit, nlon, nlat, soil_emiss2_dry(1:nlon,1:nlat))

 call wrec2 (iunit, nlon, nlat, soil_emiss2_wet(1:nlon,1:nlat))

 call wrec2 (iunit, nlon, nlat, veg_root_depth(1:nlon,1:nlat))

 call wrec2 (iunit, nlon, nlat, veg_roughness(1:nlon,1:nlat))

 call wrec2 (iunit, nlon, nlat, veg_albedo(1:nlon,1:nlat))

 call wrec2 (iunit, nlon, nlat, veg_emiss1(1:nlon,1:nlat))

 call wrec2 (iunit, nlon, nlat, veg_emiss2(1:nlon,1:nlat))

 call wrec2 (iunit, nlon, nlat, water_table_depth(1:nlon,1:nlat))

 call wrec2 (iunit, nlon, nlat, tg_bottom(1:nlon,1:nlat))

 call wrec2 (iunit, nlon, nlat, qg_rel_bottom(1:nlon,1:nlat))

 call wrec2 (iunit, nlon, nlat, qg_rel_surf_approx(1:nlon,1:nlat))

 call wrec2_int (iunit, nlon, nlat, ind_lev_soil_h_bottom(1:nlon,1:nlat))

 call wrec2_int (iunit, nlon, nlat, ind_lev_soil_w_bottom(1:nlon,1:nlat))

 call wrec2 (iunit, nlon, nlat, snow_dirt(1:nlon,1:nlat))

 write(iunit) veg_lai_max

! flush (iunit)
 close (iunit)

 write (*,*)
 write (*,*) 'Output file ',trim(file_output),' written' 
 write (*,*)

#ifdef oper
 open  (iunit, file=trim(file_output)//'.txt', status='unknown')
 write (iunit,'(2a)') trim(file_output),' is full and closed'
 close (iunit)
#endif

return
end
!=======================================================================
subroutine read_param_const(x0, y0, alon0, alat0, htop, htopvar, flag_data)

! Reads input file with all constant (in time) 
! model physiographical parameters

use param, only : nlon, nlat, nlevg, dlon, dlat, nst, nvt
use mod_fields, only : fmask, rgm, rgq, soil_map, veg_map, &
 water_table_depth, tg_bottom, qg_rel_bottom, qg_rel_surf_approx, &
 soil_qmax, soil_qmin, soil_c, soil_rho, soil_psi, soil_k, soil_par_b, soil_par_c, soil_qrel_wilt, soil_qrel_ref, &
 soil_albedo_dry, soil_albedo_wet, soil_emiss1_dry, soil_emiss1_wet, soil_emiss2_dry, soil_emiss2_wet, &
 veg_lai, veg_frac, veg_root_depth, veg_roughness, veg_albedo, veg_emiss1, veg_emiss2, &
 snow_dirt, ind_lev_soil_h_bottom, ind_lev_soil_w_bottom, veg_lai_max

implicit none
integer :: iunit=60, flag_data, &
 nlon_read, nlat_read, nlevg_read, nst_read, nvt_read, ierr_open, flag, k
real :: x0, y0, x0_read, y0_read, alon0, alat0, alon0_read, alat0_read, dlon_read, dlat_read
real, dimension(nlon,nlat) :: htop, htopvar
character(len=30) :: file_data="model_param_constant.bin"

 flag_data=0

 open (iunit, file=file_data, form='unformatted', status='old', iostat=ierr_open)

 if (ierr_open /= 0) then
   flag_data=1
   return
 endif

 read (iunit) nlon_read, nlat_read, nlevg_read, dlon_read, dlat_read, x0_read, y0_read, alon0_read, alat0_read, nst_read, nvt_read

 flag=1
 if (nlon_read /= nlon) flag=0
 if (nlat_read /= nlat) flag=0
 if (nlevg_read /= nlevg) flag=0
 if (dlon_read /= dlon) flag=0
 if (dlat_read /= dlat) flag=0
 if (x0_read /= x0) flag=0
 if (y0_read /= y0) flag=0
 if (alon0_read /= alon0) flag=0
 if (alat0_read /= alat0) flag=0

 if (flag == 0) then
   write (*,*)
   write (*,*) "Error in header parameters in input file ,",trim(file_data),", not coincident with defined parameters"
   write (*,*) "Model nlon, nlat, nlevg, dlon, dlat, x0, y0, alon0, alat0 :", &
 nlon, nlat, nlevg, dlon, dlat, x0, y0, alon0, alat0
   write (*,*)"Read nlon, nlat, nlevg, dlon, dlat, x0, y0, alon0, alat0 :", &
 nlon_read, nlat_read, nlevg_read, dlon_read, dlat_read, x0_read, y0_read, alon0_read, alat0_read
   write (*,*) "File with constant model parameter fields is erroneous and not will be used"
   write (*,*)
   flag_data=1
   return
 endif

 call rrec2 (iunit, nlon, nlat, fmask(1:nlon,1:nlat))

 call rrec2 (iunit, nlon, nlat, htop(1:nlon,1:nlat))

 call rrec2 (iunit, nlon, nlat, htopvar(1:nlon,1:nlat))

 do k = 1,nst+1
   call rrec2 (iunit, nlon, nlat, soil_map(1:nlon,1:nlat,k))
 enddo

 do k = 1,nvt+1
   call rrec2 (iunit, nlon, nlat, veg_map(1:nlon,1:nlat,k))
 enddo

 do k = 1,nlevg
   call rrec2 (iunit, nlon, nlat, soil_qmax(1:nlon,1:nlat,k))
 enddo

 do k = 1,nlevg
   call rrec2 (iunit, nlon, nlat, soil_qmin(1:nlon,1:nlat,k))
 enddo

 do k = 1,nlevg
   call rrec2 (iunit, nlon, nlat, soil_c(1:nlon,1:nlat,k))
 enddo

 do k = 1,nlevg
   call rrec2 (iunit, nlon, nlat, soil_rho(1:nlon,1:nlat,k))
 enddo

 do k = 1,nlevg
   call rrec2 (iunit, nlon, nlat, soil_psi(1:nlon,1:nlat,k))
 enddo

 do k = 1,nlevg
   call rrec2 (iunit, nlon, nlat, soil_k(1:nlon,1:nlat,k))
 enddo

 do k = 1,nlevg
   call rrec2 (iunit, nlon, nlat, soil_par_b(1:nlon,1:nlat,k))
 enddo

 do k = 1,nlevg
   call rrec2 (iunit, nlon, nlat, soil_par_c(1:nlon,1:nlat,k))
 enddo

 do k = 1,nlevg
   call rrec2 (iunit, nlon, nlat, soil_qrel_wilt(1:nlon,1:nlat,k))
 enddo

 do k = 1,nlevg
   call rrec2 (iunit, nlon, nlat, soil_qrel_ref(1:nlon,1:nlat,k))
 enddo

 call rrec2 (iunit, nlon, nlat, soil_albedo_dry(1:nlon,1:nlat))

 call rrec2 (iunit, nlon, nlat, soil_albedo_wet(1:nlon,1:nlat))

 call rrec2 (iunit, nlon, nlat, soil_emiss1_dry(1:nlon,1:nlat))

 call rrec2 (iunit, nlon, nlat, soil_emiss1_wet(1:nlon,1:nlat))

 call rrec2 (iunit, nlon, nlat, soil_emiss2_dry(1:nlon,1:nlat))

 call rrec2 (iunit, nlon, nlat, soil_emiss2_wet(1:nlon,1:nlat))

 call rrec2 (iunit, nlon, nlat, veg_root_depth(1:nlon,1:nlat))

 call rrec2 (iunit, nlon, nlat, veg_roughness(1:nlon,1:nlat))

 call rrec2 (iunit, nlon, nlat, veg_albedo(1:nlon,1:nlat))

 call rrec2 (iunit, nlon, nlat, veg_emiss1(1:nlon,1:nlat))

 call rrec2 (iunit, nlon, nlat, veg_emiss2(1:nlon,1:nlat))

 call rrec2 (iunit, nlon, nlat, water_table_depth(1:nlon,1:nlat))

 call rrec2 (iunit, nlon, nlat, tg_bottom(1:nlon,1:nlat))

 call rrec2 (iunit, nlon, nlat, qg_rel_bottom(1:nlon,1:nlat))

 call rrec2 (iunit, nlon, nlat, qg_rel_surf_approx(1:nlon,1:nlat))

 call rrec2_int (iunit, nlon, nlat, ind_lev_soil_h_bottom(1:nlon,1:nlat))

 call rrec2_int (iunit, nlon, nlat, ind_lev_soil_w_bottom(1:nlon,1:nlat))

 call rrec2 (iunit, nlon, nlat, snow_dirt(1:nlon,1:nlat))

 read (iunit) veg_lai_max

 close (iunit)

 write (*,*)
 write (*,*) 'File with model constant parameter fields ',trim(file_data),' read'
 write (*,*)

return
end
!=======================================================================
subroutine rrec2 (kunit, nlon, nlat, vect)

implicit none

integer :: kunit, nlon, nlat
real, dimension(nlon,nlat) :: vect

 read(kunit) vect(1:nlon,1:nlat)

return
end subroutine rrec2
!=======================================================================
subroutine rrec2_int (kunit, nlon, nlat, ivect)

implicit none

integer :: kunit, nlon, nlat
integer, dimension(nlon,nlat) :: ivect

 read(kunit) ivect(1:nlon,1:nlat)

return
end subroutine rrec2_int
!=======================================================================
subroutine read_forecast_mhf_soil(iflag)

! Reads a MHF file of BOLAM with surface and soil variables 
! simulated by a forecast run with siutable validity date and hour

use param
use mod_fields, only : iceth, fice, tsurf, tgsurf, tg, qvsurf, qgsurf, qg, fice_soil_surf, fice_soil, &
 snow, snow_lev, snow_t, snow_fice, snow_age, snow_melt_age, snow_dens, &
 iceth_frc, fice_frc, tsurf_frc, tgsurf_frc, tg_frc, qvsurf_frc, qgsurf_frc, qg_frc, fice_soil_surf_frc, fice_soil_frc, &
 snow_frc, snow_lev_frc, snow_t_frc, snow_fice_frc, snow_age_frc, snow_melt_age_frc, snow_dens_frc

implicit none

integer :: iflag, iunit=11, iopen_err=0, ird, ierr=0, i, j, k, ndayr, &
 year_frc, month_frc, day_frc, hour_frc, minute_frc
character(len=30) :: file_inp="bolam_forecast_soil.mhf",str_date0,str_frc_term,str_date_frc
integer, dimension(50) :: nfdr_frc
real, dimension(200)   :: pdr_frc
real, dimension(nlon,nlat) :: field2d_add, snow_albedo_rfc, runoff_tot_frc

 iflag=1

 open (iunit,file=file_inp,form='unformatted',status='old',iostat=iopen_err)

 if (iopen_err /= 0) then
   write (*,*)
   write (*,*) '   Not found file ',trim(file_inp),' with forecast surface, soil and snow data' 
   write (*,*)
   iflag=0
   return
 else
   write (*,*)
   write (*,*) '   Reading on unit ',iunit,' file ',trim(file_inp),&
 ' with forecast surface, soil and snow data, variables redefinition' 
   write (*,*)
 endif

 read(iunit) nfdr_frc
 read(iunit) pdr_frc

 ierr=0
 if (nfdr_frc(2) /= nfdr(2)) ierr=ierr+1 ! nlon
 if (nfdr_frc(3) /= nfdr(3)) ierr=ierr+1 ! nlat
 if (abs(nfdr_frc(15)-nfdr(15)) > 1.e-4) ierr=ierr+1 ! nlevg
 if (abs(pdr_frc(2)-pdr(2)) > 1.e-4) ierr=ierr+1 ! dlon
 if (abs(pdr_frc(1)-pdr(1)) > 1.e-4) ierr=ierr+1 ! dlat
 if (abs(pdr_frc(39)-pdr(39)) > 1.e-4) ierr=ierr+1 ! x0
 if (abs(pdr_frc(38)-pdr(38)) > 1.e-4) ierr=ierr+1 ! y0
 if (abs(pdr_frc(5)-pdr(5)) > 1.e-4) ierr=ierr+1 ! alon0
 if (abs(pdr_frc(4)-pdr(4)) > 1.e-4) ierr=ierr+1 ! alat0

 if (ierr /= 0) then
   write (*,*)
   write (*,*) "Error in header parameters in input file ,",trim(file_inp),", not coincident with defined parameters"
   write (*,*) "Model nlon, nlat, nlevg, dlon, dlat, x0, y0, alon0, alat0 :", &
 nfdr(2), nfdr(3), nfdr(15), pdr(2), pdr(1), pdr(39), pdr(38), pdr(5), pdr(4)
   write (*,*)"Read nlon, nlat, nlevg, dlon, dlat, x0, y0, alon0, alat0 :", &
 nfdr_frc(2), nfdr_frc(3), nfdr_frc(15), pdr_frc(2), pdr_frc(1), pdr_frc(39), pdr_frc(38), pdr_frc(5), pdr_frc(4)
   write (*,*)
   iflag=0
   return
 endif

 call calendar (nfdr_frc(5), nfdr_frc(6), nfdr_frc(7), nfdr_frc(8), nfdr_frc(9), &
 nfdr_frc(10), nfdr_frc(11), nfdr_frc(12), &
 year_frc, month_frc, day_frc, hour_frc, minute_frc, ndayr)

 ierr=0
 if (year_frc /= nfdr(5)) ierr=ierr+1
 if (month_frc /= nfdr(6)) ierr=ierr+1
 if (day_frc /= nfdr(7)) ierr=ierr+1
 if (hour_frc /= nfdr(8)) ierr=ierr+1
 if (minute_frc /= nfdr(9)) ierr=ierr+1

 if (ierr /= 0) then
   write (*,*)
   write (*,*) "Error in date and time parameters in input file ",trim(file_inp)
   write (*,*) "Initialisation date and time (year, month, day, hour, minute): ", nfdr(5:9)
   write (*,*) "Read data valid at followind date and time (year, month, day, hour, minute): ", & 
 year_frc, month_frc, day_frc, hour_frc, minute_frc
   write (*,*) "Read initial date and time (year, month, day, hour, minute): ", nfdr_frc(5:9)
   print *,"Read forecast term (day, hour, minute): ",nfdr_frc(10:12)
   write (*,*)
   iflag=0
   return
 endif

 write (str_date0,'(4(i2.2,a),i4.4)') nfdr_frc(8),':',nfdr_frc(9),' ',nfdr_frc(7),'/',nfdr_frc(6),'/',nfdr_frc(5)
 write (str_frc_term,'(a,i3,a,2(i2,a))') '+',nfdr_frc(10),'day ',nfdr_frc(11),'hour ',nfdr_frc(12),'minute'
 write (str_date_frc,'(4(i2.2,a),i4.4)') hour_frc,':',minute_frc,' ',day_frc,'/',month_frc,'/',year_frc

 do ird=1,4 ! veg_lai, veg_frac, rgm, rgq
   call rrec2 (iunit, nlon, nlat, field2d_add(1:nlon,1:nlat))
 enddo

 call rrec2 (iunit, nlon, nlat, iceth_frc(1:nlon,1:nlat))

 call rrec2 (iunit, nlon, nlat, fice_frc(1:nlon,1:nlat))

 do ird=1,7 ! albedo, emismap1, emismap2, cloud, totpre, conpre, snfall
   call rrec2 (iunit, nlon, nlat, field2d_add(1:nlon,1:nlat))
 enddo

 call rrec2 (iunit, nlon, nlat, tsurf_frc(1:nlon,1:nlat))

 call rrec2 (iunit, nlon, nlat, tgsurf_frc(1:nlon,1:nlat))

 do k = 1,nlevg
   call rrec2 (iunit, nlon, nlat, tg_frc(1:nlon,1:nlat,k))
 enddo

 call rrec2 (iunit, nlon, nlat, qvsurf_frc(1:nlon,1:nlat))

 call rrec2 (iunit, nlon, nlat, qgsurf_frc(1:nlon,1:nlat))

 do k = 1,nlevg
   call rrec2 (iunit, nlon, nlat, qg_frc(1:nlon,1:nlat,k))
 enddo

 call rrec2 (iunit, nlon, nlat, fice_soil_surf_frc(1:nlon,1:nlat))

 do k = 1,nlevg
   call rrec2 (iunit, nlon, nlat, fice_soil_frc(1:nlon,1:nlat,k))
 enddo

 call rrec2 (iunit, nlon, nlat, snow_frc(1:nlon,1:nlat))

 do k = 1,nlevsnow
   call rrec2 (iunit, nlon, nlat, snow_lev_frc(1:nlon,1:nlat,k))
 enddo

 do k = 1,nlevsnow
   call rrec2 (iunit, nlon, nlat, snow_t_frc(1:nlon,1:nlat,k))
 enddo

 do k = 1,nlevsnow
   call rrec2 (iunit, nlon, nlat, snow_fice_frc(1:nlon,1:nlat,k))
 enddo

 do k = 1,nlevsnow
   call rrec2 (iunit, nlon, nlat, snow_age_frc(1:nlon,1:nlat,k))
 enddo

 do k = 1,nlevsnow
   call rrec2 (iunit, nlon, nlat, snow_melt_age_frc(1:nlon,1:nlat,k))
 enddo

 do k = 1,nlevsnow
   call rrec2 (iunit, nlon, nlat, snow_dens_frc(1:nlon,1:nlat,k))
 enddo

! do ird=1,12 ! snow_albedo, cswfl, clwfl, chflux, cqflux, t2min, t2max, runoff, runoff_tot, 3 additional 2D fields 
!   call rrec2 (iunit, nlon, nlat, field2d_add(1:nlon,1:nlat))
! enddo

 close (iunit)

 write (*,*)
 write (*,*) "Soil and snow prognostic parameters (temperature, water content, ect.) values"
 write (*,*) "are read from forecast output data of run start"
 write (*,*) trim(str_date0)," forcast term ",trim(str_frc_term)," valid ",trim(str_date_frc)

return
end
!=======================================================================
subroutine read_forecast_mhf_soil_old(iflag)

! Reads a MHF file of BOLAM with surface and soil variables 
! simulated by a forecast run with siutable validity date and hour

use param
use mod_fields, only : iceth, fice, tsurf, tgsurf, tg, qvsurf, qgsurf, qg, fice_soil_surf, fice_soil, &
 snow, snow_lev, snow_t, snow_fice, snow_age, snow_melt_age, snow_dens, &
 iceth_frc, fice_frc, tsurf_frc, tgsurf_frc, tg_frc, qvsurf_frc, qgsurf_frc, qg_frc, fice_soil_surf_frc, fice_soil_frc, &
 snow_frc, snow_lev_frc, snow_t_frc, snow_fice_frc, snow_age_frc, snow_melt_age_frc, snow_dens_frc

implicit none

integer :: iflag, iunit=11, iopen_err=0, ird, ierr=0, i, j, k, ndayr, &
 year_frc, month_frc, day_frc, hour_frc, minute_frc
character(len=30) :: file_inp="bolam_forecast_soil.mhf",str_date0,str_frc_term,str_date_frc
integer, dimension(50) :: nfdr_frc
real, dimension(200)   :: pdr_frc
real, dimension(nlon,nlat) :: field2d_add, snow_albedo_rfc, runoff_tot_frc

 iflag=1

 open (iunit,file=file_inp,form='unformatted',status='old',iostat=iopen_err)

 if (iopen_err /= 0) then
   write (*,*)
   write (*,*) '   Not found file ',trim(file_inp),' with forecast surface, soil and snow data' 
   write (*,*)
   iflag=0
   return
 else
   write (*,*)
   write (*,*) '   Reading on unit ',iunit,' file ',trim(file_inp),&
 ' with forecast surface, soil and snow data, variables redefinition' 
   write (*,*)
 endif

 read(iunit) nfdr_frc
 read(iunit) pdr_frc

 ierr=0
 if (nfdr_frc(2) /= nfdr(2)) ierr=ierr+1 ! nlon
 if (nfdr_frc(3) /= nfdr(3)) ierr=ierr+1 ! nlat
 if (abs(nfdr_frc(15)-nfdr(15)) > 1.e-4) ierr=ierr+1 ! nlevg
 if (abs(pdr_frc(2)-pdr(2)) > 1.e-4) ierr=ierr+1 ! dlon
 if (abs(pdr_frc(1)-pdr(1)) > 1.e-4) ierr=ierr+1 ! dlat
 if (abs(pdr_frc(39)-pdr(39)) > 1.e-4) ierr=ierr+1 ! x0
 if (abs(pdr_frc(38)-pdr(38)) > 1.e-4) ierr=ierr+1 ! y0
 if (abs(pdr_frc(5)-pdr(5)) > 1.e-4) ierr=ierr+1 ! alon0
 if (abs(pdr_frc(4)-pdr(4)) > 1.e-4) ierr=ierr+1 ! alat0

 if (ierr /= 0) then
   write (*,*)
   write (*,*) "Error in header parameters in input file ,",trim(file_inp),", not coincident with defined parameters"
   write (*,*) "Model nlon, nlat, nlevg, dlon, dlat, x0, y0, alon0, alat0 :", &
 nfdr(2), nfdr(3), nfdr(15), pdr(2), pdr(1), pdr(39), pdr(38), pdr(5), pdr(4)
   write (*,*)"Read nlon, nlat, nlevg, dlon, dlat, x0, y0, alon0, alat0 :", &
 nfdr_frc(2), nfdr_frc(3), nfdr_frc(15), pdr_frc(2), pdr_frc(1), pdr_frc(39), pdr_frc(38), pdr_frc(5), pdr_frc(4)
   write (*,*)
   iflag=0
   return
 endif

 call calendar (nfdr_frc(5), nfdr_frc(6), nfdr_frc(7), nfdr_frc(8), nfdr_frc(9), &
 nfdr_frc(10), nfdr_frc(11), nfdr_frc(12), &
 year_frc, month_frc, day_frc, hour_frc, minute_frc, ndayr)

 ierr=0
 if (year_frc /= nfdr(5)) ierr=ierr+1
 if (month_frc /= nfdr(6)) ierr=ierr+1
 if (day_frc /= nfdr(7)) ierr=ierr+1
 if (hour_frc /= nfdr(8)) ierr=ierr+1
 if (minute_frc /= nfdr(9)) ierr=ierr+1

 if (ierr /= 0) then
   write (*,*)
   write (*,*) "Error in date and time parameters in input file ",trim(file_inp)
   write (*,*) "Initialisation date and time (year, month, day, hour, minute): ", nfdr(5:9)
   write (*,*) "Read data valid at followind date and time (year, month, day, hour, minute): ", & 
 year_frc, month_frc, day_frc, hour_frc, minute_frc
   write (*,*) "Read initial date and time (year, month, day, hour, minute): ", nfdr_frc(5:9)
   print *,"Read forecast term (day, hour, minute): ",nfdr_frc(10:12)
   write (*,*)
   iflag=0
   return
 endif

 write (str_date0,'(4(i2.2,a),i4.4)') nfdr_frc(8),':',nfdr_frc(9),' ',nfdr_frc(7),'/',nfdr_frc(6),'/',nfdr_frc(5)
 write (str_frc_term,'(a,i3,a,2(i2,a))') '+',nfdr_frc(10),'day ',nfdr_frc(11),'hour ',nfdr_frc(12),'minute'
 write (str_date_frc,'(4(i2.2,a),i4.4)') hour_frc,':',minute_frc,' ',day_frc,'/',month_frc,'/',year_frc

 do ird=1,4 ! veg_lai, veg_frac, rgm, rgq
   do j = 1,nlat
     read(iunit) (field2d_add(i,j),i=1,nlon)
   enddo
 enddo

 do j = 1,nlat
   read(iunit) (iceth_frc(i,j),i=1,nlon)
 enddo
 do j = 1,nlat
   read(iunit) (fice_frc(i,j),i=1,nlon)
 enddo

 do ird=1,7 ! albedo, emismap1, emismap2, cloud, totpre, conpre, snfall
   do j = 1,nlat
     read(iunit) (field2d_add(i,j),i=1,nlon)
   enddo
 enddo

 do j = 1,nlat
   read(iunit) (tsurf_frc(i,j),i=1,nlon)
 enddo
 do j = 1,nlat
   read(iunit) (tgsurf_frc(i,j),i=1,nlon)
 enddo
 do k = 1,nlevg
   do j = 1,nlat
     read(iunit) (tg_frc(i,j,k),i=1,nlon)
   enddo
 enddo
 do j = 1,nlat
   read(iunit) (qvsurf_frc(i,j),i=1,nlon)
 enddo
 do j = 1,nlat
   read(iunit) (qgsurf_frc(i,j),i=1,nlon)
 enddo
 do k = 1,nlevg
   do j = 1,nlat
     read(iunit) (qg_frc(i,j,k),i=1,nlon)
   enddo
 enddo
 do j = 1,nlat
   read(iunit) (fice_soil_surf_frc(i,j),i=1,nlon)
 enddo
 do k = 1,nlevg
   do j = 1,nlat
     read(iunit) (fice_soil_frc(i,j,k),i=1,nlon)
   enddo
 enddo
 do j = 1,nlat
   read(iunit) (snow_frc(i,j),i=1,nlon)
 enddo
 do k = 1,nlevsnow
   do j = 1,nlat
     read(iunit) (snow_lev_frc(i,j,k),i=1,nlon)
   enddo
 enddo
 do k = 1,nlevsnow
   do j = 1,nlat
     read(iunit) (snow_t_frc(i,j,k),i=1,nlon)
   enddo
 enddo
 do k = 1,nlevsnow
   do j = 1,nlat
     read(iunit) (snow_fice_frc(i,j,k),i=1,nlon)
   enddo
 enddo
 do k = 1,nlevsnow
   do j = 1,nlat
     read(iunit) (snow_age_frc(i,j,k),i=1,nlon)
   enddo
 enddo
 do k = 1,nlevsnow
   do j = 1,nlat
     read(iunit) (snow_melt_age_frc(i,j,k),i=1,nlon)
   enddo
 enddo
 do k = 1,nlevsnow
   do j = 1,nlat
     read(iunit) (snow_dens_frc(i,j,k),i=1,nlon)
   enddo
 enddo

! do ird=1,12 ! snow_albedo, cswfl, clwfl, chflux, cqflux, t2min, t2max, runoff, runoff_tot, 3 additional 2D fields 
!   do j = 1,nlat
!     read(iunit) (field2d_add(i,j),i=1,nlon)
!   enddo
! enddo

 close (iunit)

 write (*,*)
 write (*,*) "Soil and snow prognostic parameters (temperature, water content, ect.) values"
 write (*,*) "are read from forecast output data of run start"
 write (*,*) trim(str_date0)," forcast term ",trim(str_frc_term)," valid ",trim(str_date_frc)

return
end
!=======================================================================
subroutine read_forecast_mhf_soil_access_direct(iflag)

! Reads a MHF file of BOLAM, written in access direct form, 
! with surface and soil variables 
! simulated by a forecast run with siutable validity date and hour

use param
use mod_fields, only : iceth, fice, tsurf, tgsurf, tg, qvsurf, qgsurf, qg, fice_soil_surf, fice_soil, &
 snow, snow_lev, snow_t, snow_fice, snow_age, snow_melt_age, snow_dens, &
 iceth_frc, fice_frc, tsurf_frc, tgsurf_frc, tg_frc, qvsurf_frc, qgsurf_frc, qg_frc, fice_soil_surf_frc, fice_soil_frc, &
 snow_frc, snow_lev_frc, snow_t_frc, snow_fice_frc, snow_age_frc, snow_melt_age_frc, snow_dens_frc

implicit none

integer :: iflag, iunit=11, iopen_err=0, rec_ind, rec_lon, ird, ierr=0, i, j, k, ndayr, &
 year_frc, month_frc, day_frc, hour_frc, minute_frc
character(len=30) :: file_inp="bolam_forecast_soil.bin",str_date0,str_frc_term,str_date_frc
integer, dimension(50) :: nfdr_frc
real, dimension(200)   :: pdr_frc
real, dimension(nlon,nlat) :: field2d_add, snow_albedo_rfc, runoff_tot_frc
real, dimension (nlon) :: array_1d

 rec_lon = nlon
 rec_ind = 0

 iflag=1

 open (iunit,file=file_inp,form='unformatted',status='old',access='direct', recl=rec_lon*4,iostat=iopen_err)

 if (iopen_err /= 0) then
   write (*,*)
   write (*,*) '   Not found file ',trim(file_inp),' with forecast surface, soil and snow data' 
   write (*,*)
   iflag=0
   return
 else
   write (*,*)
   write (*,*) '   Reading on unit ',iunit,' file ',trim(file_inp),&
 ' with forecast surface, soil and snow data, variables redefinition' 
   write (*,*)
 endif

 rec_ind = rec_ind+1
 read(iunit, rec=rec_ind) array_1d
 nfdr_frc(1:50) = int(array_1d(1:50))
 rec_ind = rec_ind+1
 read(iunit, rec=rec_ind) array_1d
 pdr_frc(1:200) = array_1d(1:200)

 ierr=0
 if (nfdr_frc(2) /= nfdr(2)) ierr=ierr+1 ! nlon
 if (nfdr_frc(3) /= nfdr(3)) ierr=ierr+1 ! nlat
 if (nfdr_frc(15) /= nfdr(15)) ierr=ierr+1 ! nlevg
 if (pdr_frc(2) /= pdr(2)) ierr=ierr+1 ! dlon
 if (pdr_frc(1) /= pdr(1)) ierr=ierr+1 ! dlat
 if (pdr_frc(39) /= pdr(39)) ierr=ierr+1 ! x0
 if (pdr_frc(38) /= pdr(38)) ierr=ierr+1 ! y0
 if (pdr_frc(5) /= pdr(5)) ierr=ierr+1 ! alon0
 if (pdr_frc(4) /= pdr(4)) ierr=ierr+1 ! alat0

 if (ierr /= 0) then
   write (*,*)
   write (*,*) "Error in header parameters in input file ,",trim(file_inp),", not coincident with defined parameters"
   write (*,*) "Model nlon, nlat, nlevg, dlon, dlat, x0, y0, alon0, alat0 :", &
 nfdr(2), nfdr(3), nfdr(15), pdr(2), pdr(1), pdr(39), pdr(38), pdr(5), pdr(4)
   write (*,*)"Read nlon, nlat, nlevg, dlon, dlat, x0, y0, alon0, alat0 :", &
 nfdr_frc(2), nfdr_frc(3), nfdr_frc(15), pdr_frc(2), pdr_frc(1), pdr_frc(39), pdr_frc(38), pdr_frc(5), pdr_frc(4)
   write (*,*)
   iflag=0
   return
 endif

 call calendar (nfdr_frc(5), nfdr_frc(6), nfdr_frc(7), nfdr_frc(8), nfdr_frc(9), &
 nfdr_frc(10), nfdr_frc(11), nfdr_frc(12), &
 year_frc, month_frc, day_frc, hour_frc, minute_frc, ndayr)

 ierr=0
 if (year_frc /= nfdr(5)) ierr=ierr+1
 if (month_frc /= nfdr(6)) ierr=ierr+1
 if (day_frc /= nfdr(7)) ierr=ierr+1
 if (hour_frc /= nfdr(8)) ierr=ierr+1
 if (minute_frc /= nfdr(9)) ierr=ierr+1

 if (ierr /= 0) then
   write (*,*)
   write (*,*) "Error in date and time parameters in input file ",trim(file_inp)
   write (*,*) "Initialisation date and time (year, month, day, hour, minute): ", nfdr(5:9)
   write (*,*) "Read data valid at followind date and time (year, month, day, hour, minute): ", & 
 year_frc, month_frc, day_frc, hour_frc, minute_frc
   write (*,*) "Read initial date and time (year, month, day, hour, minute): ", nfdr_frc(5:9)
   print *,"Read forecast term (day, hour, minute): ",nfdr_frc(10:12)
   write (*,*)
   iflag=0
   return
 endif

 write (str_date0,'(4(i2.2,a),i4.4)') nfdr_frc(8),':',nfdr_frc(9),' ',nfdr_frc(7),'/',nfdr_frc(6),'/',nfdr_frc(5)
 write (str_frc_term,'(a,i3,a,2(i2,a))') '+',nfdr_frc(10),'day ',nfdr_frc(11),'hour ',nfdr_frc(12),'minute'
 write (str_date_frc,'(4(i2.2,a),i4.4)') hour_frc,':',minute_frc,' ',day_frc,'/',month_frc,'/',year_frc

 do ird=1,4 ! veg_lai, veg_frac, rgm, rgq
   do j = 1,nlat
     rec_ind = rec_ind+1
     read(iunit, rec=rec_ind) array_1d
   enddo
 enddo

 do j = 1,nlat
   rec_ind = rec_ind+1
   read(iunit, rec=rec_ind) array_1d
   iceth_frc(1:nlon,j) = array_1d(1:rec_lon)
 enddo
 do j = 1,nlat
   rec_ind = rec_ind+1
   read(iunit, rec=rec_ind) array_1d
   fice_frc(1:nlon,j) = array_1d(1:rec_lon)
 enddo

 do ird=1,7 ! albedo, emismap1, emismap2, cloud, totpre, conpre, snfall
   do j = 1,nlat
     rec_ind = rec_ind+1
     read(iunit, rec=rec_ind) array_1d
   enddo
 enddo

 do j = 1,nlat
   rec_ind = rec_ind+1
   read(iunit, rec=rec_ind) array_1d
   tsurf_frc(1:nlon,j) = array_1d(1:rec_lon)
 enddo
 do j = 1,nlat
   rec_ind = rec_ind+1
   read(iunit, rec=rec_ind) array_1d
   tgsurf_frc(1:nlon,j) = array_1d(1:rec_lon)
 enddo
 do k = 1,nlevg
   do j = 1,nlat
     rec_ind = rec_ind+1
     read(iunit, rec=rec_ind) array_1d
     tg_frc(1:nlon,j,k) = array_1d(1:rec_lon)
   enddo
 enddo
 do j = 1,nlat
   rec_ind = rec_ind+1
   read(iunit, rec=rec_ind) array_1d
   qvsurf_frc(1:nlon,j) = array_1d(1:rec_lon)
 enddo
 do j = 1,nlat
   rec_ind = rec_ind+1
   read(iunit, rec=rec_ind) array_1d
   qgsurf_frc(1:nlon,j) = array_1d(1:rec_lon)
 enddo
 do k = 1,nlevg
   do j = 1,nlat
     rec_ind = rec_ind+1
     read(iunit, rec=rec_ind) array_1d
     qg_frc(1:nlon,j,k) = array_1d(1:rec_lon)
   enddo
 enddo
 do j = 1,nlat
   rec_ind = rec_ind+1
   read(iunit, rec=rec_ind) array_1d
   fice_soil_surf_frc(1:nlon,j) = array_1d(1:rec_lon)
 enddo
 do k = 1,nlevg
   do j = 1,nlat
     rec_ind = rec_ind+1
     read(iunit, rec=rec_ind) array_1d
     fice_soil_frc(1:nlon,j,k) = array_1d(1:rec_lon)
   enddo
 enddo
 do j = 1,nlat
   rec_ind = rec_ind+1
   read(iunit, rec=rec_ind) array_1d
   snow_frc(1:nlon,j) = array_1d(1:rec_lon)
 enddo
 do k = 1,nlevsnow
   do j = 1,nlat
     rec_ind = rec_ind+1
     read(iunit, rec=rec_ind) array_1d
     snow_lev_frc(1:nlon,j,k) = array_1d(1:rec_lon)
   enddo
 enddo
 do k = 1,nlevsnow
   do j = 1,nlat
     rec_ind = rec_ind+1
     read(iunit, rec=rec_ind) array_1d
     snow_t_frc(1:nlon,j,k) = array_1d(1:rec_lon)
   enddo
 enddo
 do k = 1,nlevsnow
   do j = 1,nlat
     rec_ind = rec_ind+1
     read(iunit, rec=rec_ind) array_1d
     snow_fice_frc(1:nlon,j,k) = array_1d(1:rec_lon)
   enddo
 enddo
 do k = 1,nlevsnow
   do j = 1,nlat
     rec_ind = rec_ind+1
     read(iunit, rec=rec_ind) array_1d
     snow_age_frc(1:nlon,j,k) = array_1d(1:rec_lon)
   enddo
 enddo
 do k = 1,nlevsnow
   do j = 1,nlat
     rec_ind = rec_ind+1
     read(iunit, rec=rec_ind) array_1d
     snow_melt_age_frc(1:nlon,j,k) = array_1d(1:rec_lon)
   enddo
 enddo
 do k = 1,nlevsnow
   do j = 1,nlat
     rec_ind = rec_ind+1
     read(iunit, rec=rec_ind) array_1d
     snow_dens_frc(1:nlon,j,k) = array_1d(1:rec_lon)
   enddo
 enddo

! do ird=1,12 ! snow_albedo, cswfl, clwfl, chflux, cqflux, t2min, t2max, runoff, runoff_tot, 3 additional 2D fields 
!   do j = 1,nlat
!     rec_ind = rec_ind+1
!     read(iunit, rec=rec_ind) array_1d
!   enddo
! enddo

 close (iunit)

 write (*,*)
 write (*,*) "Soil and snow prognostic parameters (temperature, water content, ect.) values"
 write (*,*) "are read from forecast output data of run start"
 write (*,*) trim(str_date0)," forcast term ",trim(str_frc_term)," valid ",trim(str_date_frc)

return
end
!=======================================================================
subroutine read_forecast_mhf_ascii_soil(iflag)

! Reads a MHF file of BOLAM with surface and soil variables 
! simulated by forecast (+24) run of previous day 

use param
use mod_fields, only : iceth, fice, tsurf, tgsurf, tg, qvsurf, qgsurf, qg, fice_soil_surf, fice_soil, &
 snow, snow_lev, snow_t, snow_fice, snow_age, snow_melt_age, snow_dens, &
 iceth_frc, fice_frc, tsurf_frc, tgsurf_frc, tg_frc, qvsurf_frc, qgsurf_frc, qg_frc, fice_soil_surf_frc, fice_soil_frc, &
 snow_frc, snow_lev_frc, snow_t_frc, snow_fice_frc, snow_age_frc, snow_melt_age_frc, snow_dens_frc

implicit none

integer :: iflag, iunit=11, iopen_err=0, ird, ierr=0, i, j, k, ndayr, &
 year_frc, month_frc, day_frc, hour_frc, minute_frc
character(len=30) :: file_inp="bolam_forecast_soil.dat",str_date0,str_frc_term,str_date_frc
integer, dimension(50) :: nfdr_frc
real, dimension(200)   :: pdr_frc
real, dimension(nlon,nlat) :: field2d_add, snow_albedo_rfc, runoff_tot_frc

 iflag=1

 open (iunit,file=file_inp,form='formatted',status='old',iostat=iopen_err)

 if (iopen_err /= 0) then
   write (*,*)
   write (*,*) '   Not found file ',trim(file_inp),' with forecast surface, soil and snow data' 
   write (*,*)
   iflag=0
   return
 else
   write (*,*)
   write (*,*) '   Reading on unit ',iunit,' file ',trim(file_inp),&
 ' with forecast surface, soil and snow data, variables redefinition' 
   write (*,*)
 endif

 read(iunit,1002) nfdr_frc
 read(iunit,1001) pdr_frc

 ierr=0
 if (nfdr_frc(2) /= nfdr(2)) ierr=ierr+1 ! nlon
 if (nfdr_frc(3) /= nfdr(3)) ierr=ierr+1 ! nlat
 if (nfdr_frc(15) /= nfdr(15)) ierr=ierr+1 ! nlevg
 if (pdr_frc(2) /= pdr(2)) ierr=ierr+1 ! dlon
 if (pdr_frc(1) /= pdr(1)) ierr=ierr+1 ! dlat
 if (pdr_frc(39) /= pdr(39)) ierr=ierr+1 ! x0
 if (pdr_frc(38) /= pdr(38)) ierr=ierr+1 ! y0
 if (pdr_frc(5) /= pdr(5)) ierr=ierr+1 ! alon0
 if (pdr_frc(4) /= pdr(4)) ierr=ierr+1 ! alat0

 if (ierr /= 0) then
   write (*,*)
   write (*,*) "Error in header parameters in input file ,",trim(file_inp),", not coincident with defined parameters"
   write (*,*) "Model nlon, nlat, nlevg, dlon, dlat, x0, y0, alon0, alat0 :", &
 nfdr(2), nfdr(3), nfdr(15), pdr(2), pdr(1), pdr(39), pdr(38), pdr(5), pdr(4)
   write (*,*)"Read nlon, nlat, nlevg, dlon, dlat, x0, y0, alon0, alat0 :", &
 nfdr_frc(2), nfdr_frc(3), nfdr_frc(15), pdr_frc(2), pdr_frc(1), pdr_frc(39), pdr_frc(38), pdr_frc(5), pdr_frc(4)
   write (*,*)
   iflag=0
   return
 endif

! nfdr(5:9) - initial year, month, day, hor, minute
! nfdr(10:12) - days, hours, minutes of forecast validity

 call calendar (nfdr_frc(5), nfdr_frc(6), nfdr_frc(7), nfdr_frc(8), nfdr_frc(9),  &
 nfdr_frc(10), nfdr_frc(11), nfdr_frc(12), &
 year_frc, month_frc, day_frc, hour_frc, minute_frc, ndayr)

 ierr=0
 if (year_frc /= nfdr(5)) ierr=ierr+1
 if (month_frc /= nfdr(6)) ierr=ierr+1
 if (day_frc /= nfdr(7)) ierr=ierr+1
 if (hour_frc /= nfdr(8)) ierr=ierr+1
 if (minute_frc /= nfdr(9)) ierr=ierr+1

 if (ierr /= 0) then
   write (*,*)
   write (*,*) "Error in date and time parameters in input file ,",trim(file_inp)
   write (*,*) "Initialisation date and time (year, month, day, hour, minute): ", nfdr(5:9)
   write (*,*) "Read data valid at followind date and time (year, month, day, hour, minute): ", & 
 year_frc, month_frc, day_frc, hour_frc, minute_frc
   write (*,*) "Read initial date and time (year, month, day, hour, minute): ", nfdr_frc(5:9)
   print *,"Read forecast term (day, hour, minute): ",nfdr_frc(10:12)
   write (*,*)
   iflag=0
   return
 endif

 write (str_date0,'(4(i2.2,a),i4.4)') nfdr_frc(8),':',nfdr_frc(9),' ',nfdr_frc(7),'/',nfdr_frc(6),'/',nfdr_frc(5)
 write (str_frc_term,'(a,i3,a,2(i2,a))') '+',nfdr_frc(10),'day ',nfdr_frc(11),'hour ',nfdr_frc(12),'minute'
 write (str_date_frc,'(4(i2.2,a),i4.4)') hour_frc,':',minute_frc,' ',day_frc,'/',month_frc,'/',year_frc

 do ird=1,2 ! veg_lai, veg_frac
   do j = 1,nlat
     read(iunit,1001) (field2d_add(i,j),i=1,nlon)
   enddo
 enddo

 do j = 1,nlat
   read(iunit,1001) (iceth_frc(i,j),i=1,nlon)
 enddo
 do j = 1,nlat
   read(iunit,1001) (fice_frc(i,j),i=1,nlon)
 enddo

 do ird=1,7 ! albedo, emismap1, emismap2, cloud, totpre, conpre, snfall
   do j = 1,nlat
     read(iunit,1001) (field2d_add(i,j),i=1,nlon)
   enddo
 enddo

 do j = 1,nlat
   read(iunit,1001) (tsurf_frc(i,j),i=1,nlon)
 enddo
 do j = 1,nlat
   read(iunit,1001) (tgsurf_frc(i,j),i=1,nlon)
 enddo
 do k = 1,nlevg
   do j = 1,nlat
     read(iunit,1001) (tg_frc(i,j,k),i=1,nlon)
   enddo
 enddo
 do j = 1,nlat
   read(iunit,1001) (qvsurf_frc(i,j),i=1,nlon)
 enddo
 do j = 1,nlat
   read(iunit,1001) (qgsurf_frc(i,j),i=1,nlon)
 enddo
 do k = 1,nlevg
   do j = 1,nlat
     read(iunit,1001) (qg_frc(i,j,k),i=1,nlon)
   enddo
 enddo
 do j = 1,nlat
   read(iunit,1001) (fice_soil_surf_frc(i,j),i=1,nlon)
 enddo
 do k = 1,nlevg
   do j = 1,nlat
     read(iunit,1001) (fice_soil_frc(i,j,k),i=1,nlon)
   enddo
 enddo
 do j = 1,nlat
   read(iunit,1001) (snow_frc(i,j),i=1,nlon)
 enddo
 do k = 1,nlevsnow
   do j = 1,nlat
     read(iunit,1001) (snow_lev_frc(i,j,k),i=1,nlon)
   enddo
 enddo
 do k = 1,nlevsnow
   do j = 1,nlat
     read(iunit,1001) (snow_t_frc(i,j,k),i=1,nlon)
   enddo
 enddo
 do k = 1,nlevsnow
   do j = 1,nlat
     read(iunit,1001) (snow_fice_frc(i,j,k),i=1,nlon)
   enddo
 enddo
 do k = 1,nlevsnow
   do j = 1,nlat
     read(iunit,1001) (snow_age_frc(i,j,k),i=1,nlon)
   enddo
 enddo
 do k = 1,nlevsnow
   do j = 1,nlat
     read(iunit,1001) (snow_melt_age_frc(i,j,k),i=1,nlon)
   enddo
 enddo
 do k = 1,nlevsnow
   do j = 1,nlat
     read(iunit,1001) (snow_dens_frc(i,j,k),i=1,nlon)
   enddo
 enddo

! do ird=1,12 ! snow_albedo, cswfl, clwfl, chflux, cqflux, t2min, t2max, runoff, runoff_tot, 3 additional 2D fields 
!   do j = 1,nlat
!     read(iunit,1001) (field2d_add(i,j),i=1,nlon)
!   enddo
! enddo

 close (iunit)

 write (*,*)
 write (*,*) "Soil and snow prognostic parameters (temperature, water content, ect.) values"
 write (*,*) "are read from forecast output data of run start"
 write (*,*) trim(str_date0)," forcast term ",trim(str_frc_term)," valid ",trim(str_date_frc)

1001 format (100e16.8)
1002 format (100i10)

return
end
!=======================================================================
      subroutine livelli(psl,s,nlev,plsig,p0,alfa)

! Given the surface pressure psl=log(ps) and hybrid semi-integer levels s,
! computes plsig=log(p) on hybrid integer levels

       real s(nlev),plsig(nlev)

      zps = exp(psl)
      zsigint = .5*s(1)
      plsig(1) = log (p0*zsigint-(p0-zps)*zsigint**3*(1.+alfa*(1.-zsigint)*(2.-zsigint)))
      do k = 2, nlev
      zsigint = .5*(s(k)+s(k-1))
      plsig(k) = log (p0*zsigint-(p0-zps)*zsigint**3*(1.+alfa*(1.-zsigint)*(2.-zsigint)))
      enddo

       return
       end
! =============================================================
 subroutine conv_ifs_data(ist)

! Procedure to convert meteorological fields derived from input grib2 files
! For IFS-ECMWF data
!------------------------------------------------------------------------

use param, only : iana0, iana, jana, nlev_atm_inp_max, nlev_soil_inp_max, nlevp, nlevg_inp, &
 val_missing, pi, idate0, iperiod, iperiod_inp, slte, surf_elaborate, res_change
use gribana

implicit none

integer :: ist, ind_field, i, j, k, iii, jjj, ifl, np, iste, npt
real :: zqcmin, qcrit, zzz, zflatt, day1, day2
real, dimension(iana,jana,nlev_atm_inp_max) :: rh
integer, dimension(iana,jana) :: imaske
real, dimension(iana,jana,nlev_soil_inp_max) :: worki, workf
real, dimension(nlev_soil_inp_max) :: zdtg
real, save :: day
integer, save :: iday, nday, nste,iana_old, jana_old
integer, dimension(12) :: imon=(/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
real :: depth, depthsum
real :: tgav1,tgav2,tgav3,tgav4,qgav1,qgav2,qgav3,qgav4

 res_change = 0
 if(ist > 1.and.(iana /= iana_old.or.jana /= jana_old)) then
   print*
   print*, "The IFS-ECMWF grid resolution has changed at instant", ist
   res_change = 1
 endif
 iana_old = iana
 jana_old = jana

 mask_frame(:,:) = 1
 frame = .false.

 if (iperiod_inp(1)==1) then         ! Time unit is hour
   iperiod(1) = iperiod_inp(2)/24                             ! Days
   iperiod(2) = iperiod_inp(2)-iperiod(1)*24                  ! Hours
   iperiod(3) = 0                                             ! Minutes
 elseif (iperiod_inp(1)==0) then     ! Time unit is minute
   iperiod(1) = iperiod_inp(2)/24/60                          ! Days
   iperiod(2) = (iperiod_inp(2)-iperiod(1)*24*60)/60          ! Hours
   iperiod(3) = iperiod_inp(2)-iperiod(1)*24*60-iperiod(2)*60 ! Minutes
 endif

! For the case in which input data are defined only on frames:

 if (ist > 1) then
   ind_field = 2 ! Temperature at atm. levels in input data
   k = nlevp
   if (any(int(field3d(:,:,k,ind_field)) == int(val_missing))) then
     frame = .true.
     print *,'Input data fields defined only on frames'
     do j = 1,jana
     do i = 1,iana
       if (int(field3d(i,j,k,ind_field)) == int(val_missing)) then
         mask_frame(i,j) = 0
         field3d(i,j,:,:) = field3d(1,1,:,:)
         field3d_soil(i,j,:,:) = field3d_soil(1,1,:,:)
         field2d(i,j,:) = field2d(1,1,:)
       endif
     enddo
     enddo
   else
     print *,'Input data fields defined on full area (not on frames only)'
   endif
 endif

! Definition of meteorological fields declared in main program
! List to be checked - see subroutine READ_GRIB2_DATA

 ge(1:iana0,:,1:nlevp) = field3d(1:iana0,:,1:nlevp,1) ! Geopotential in m**2 s**-2
 te(1:iana0,:,1:nlevp) = field3d(1:iana0,:,1:nlevp,2)
 ue(1:iana0,:,1:nlevp) = field3d(1:iana0,:,1:nlevp,3)
 ve(1:iana0,:,1:nlevp) = field3d(1:iana0,:,1:nlevp,4)
 qe(1:iana0,:,1:nlevp) = field3d(1:iana0,:,1:nlevp,5)
 rh(1:iana0,:,1:nlevp) = field3d(1:iana0,:,1:nlevp,6)

 if (iana0 /= iana) then
   ge(iana,:,1:nlevp) = ge(1,:,1:nlevp)
   te(iana,:,1:nlevp) = te(1,:,1:nlevp)
   ue(iana,:,1:nlevp) = ue(1,:,1:nlevp)
   ve(iana,:,1:nlevp) = ve(1,:,1:nlevp)
   qe(iana,:,1:nlevp) = qe(1,:,1:nlevp)
   rh(iana,:,1:nlevp) = rh(1,:,1:nlevp)
 endif

! Note that in the following the last index denotes:
! 7: total cloud water content (sum of water and ice - parameter normally not available from IFS-ECMWF);
! 8: cloud liquid water; 9: cloud ice.

 zqcmin = 1.e-7
 iflag_cloude = 0
 qce (:,:,:) = 0.
 qcwe(:,:,:) = 0.
 qcie(:,:,:) = 0.

! Check availability of total cloud water (liquid + ice) in input data

 if (any(int(field3d(:,:,1:nlevp,7)) == int(val_missing)).or.maxval(field3d(:,:,1:nlevp,7)) < zqcmin) then
   if (ist<=2) print*, "Total cloud water (liquid + ice) is not available in input"
   iflag_cloude = 0
 else
   iflag_cloude = 1
   if (ist<=2) print*, "Total cloud water (liquid + ice) is available in input"
   qce(1:iana0,:,1:nlevp) = max(0., field3d(1:iana0,:,1:nlevp,7))
   if (iana0 /= iana) qce(iana,:,1:nlevp) = qce(1,:,1:nlevp)
 endif

! Check availability of cloud liquid water in input data

 if (any(int(field3d(:,:,1:nlevp,8)) == int(val_missing)).or.maxval(field3d(:,:,1:nlevp,8)) < zqcmin) then
   if (ist<=2) print*, "Cloud liquid water is not available in input"
   iflag_cloude = max(iflag_cloude, 0)
 else
   iflag_cloude = 1
   if (ist<=2) print*, "Cloud liquid water is available in input"
   qcwe(1:iana0,:,1:nlevp) = max(0., field3d(1:iana0,:,1:nlevp,8))
   if (iana0 /= iana) qcwe(iana,:,1:nlevp) = qcwe(1,:,1:nlevp)
 endif

! Check availability of cloud ice in input data

 if (any(int(field3d(:,:,1:nlevp,9)) == int(val_missing)).or.maxval(field3d(:,:,1:nlevp,9)) < zqcmin) then
   if (ist<=2) print*, "Cloud ice is not available in input"
   iflag_cloude = max(iflag_cloude, 0)
 else
   iflag_cloude = 1
   if (ist<=2) print*, "Cloud ice is available in input"
   qcie(1:iana0,:,1:nlevp) = max(0., field3d(1:iana0,:,1:nlevp,9))
   if (iana0 /= iana) qcie(iana,:,1:nlevp) = qcie(1,:,1:nlevp)
 endif

 if (iflag_cloude==1) qce(:,:,1:nlevp) = max(qce(:,:,1:nlevp), qcwe(:,:,1:nlevp)+qcie(:,:,1:nlevp))

 if (ist==1) then

! If available, the 2-D geopotential at the ground surface is used, otherwise
! the "3-D" geop. (model level) at the ground surface is used to define phige2d.
! The viceversa is not possible in case of hybrid level data: if the model level
! geopotential is missing, the program stops with a message.
! In case of input data on isobaric levels, the 3-D geopotential is not required.

   if ((any(int(field2d(:,:,2))==int(val_missing))).and.(any(int(field2d(:,:,1))==int(val_missing)))) then
     print*, "Neither the 3-D model level geopotential nor the 2-D geopotential"
     print*, "at the surface are available in IFS-ECMWF input grib data: prebolam stops!"
     stop
   endif

   if ((any(int(field2d(:,:,2))==int(val_missing))).and.(level_type == 0)) then ! case of hybrid level data
     print*, "The 3-D model level geopotential is missing in IFS-ECMWF input grib data."
     print*, "But this field is essential for a correct processing of geopotential"
     print*, "when input data are on hybrid levels, therefore prebolam stops!"
     print*, "Check the downloading of grib/grib2 fields!"
     stop
   endif

   if (any(int(field2d(:,:,1))==int(val_missing))) then
     print*, "Warning: 2-D geopotential at the surface (gaussian grid topography)"
     print*, "is missing in IFS-ECMWF input grib data."
     print*, "Therefore it is substituted with the 3-D model level geopotential."
     print*, "This may imply some errors in interpolating surface variables"
     print*, "from the IFS-ECMWF ground surface to the BOLAM ground surface."
     phige2d(1:iana0,:) = field2d(1:iana0,:,2)
   else
     phige2d(1:iana0,:) = field2d(1:iana0,:,1)
   endif

   if (any(int(field2d(:,:,3))==int(val_missing))) then
     print*, "The land-sea mask field is missing in IFS-ECMWF input grib data,"
     print*, "therefore prebolam stops!"
     stop
   else
     fmaske(1:iana0,:) = field2d(1:iana0,:,3)
   endif

   soile(1:iana0,:) = field2d(1:iana0,:,21) ! Check of missing value for soile is done below

   if (iana0 /= iana) then
     phige2d(iana,:) = phige2d(1,:)
     fmaske(iana,:)  = fmaske(1,:)
     soile(iana,:)   = soile(1,:)
   endif

! Conversion of FMASK: from 0-sea, 1-land to 1-sea, 0-land

   fmaske(:,:) = min(max(1.-fmaske(:,:),0.),1.)

 endif  ! ist=1

 phige(1:iana0,:)  = field2d(1:iana0,:,2)
 psel(1:iana0,:)   = field2d(1:iana0,:,4)
 pmsle(1:iana0,:)  = field2d(1:iana0,:,15)
 cctote(1:iana0,:) = field2d(1:iana0,:,16)
 u10e(1:iana0,:)   = field2d(1:iana0,:,17)
 v10e(1:iana0,:)   = field2d(1:iana0,:,18)
 t2e(1:iana0,:)    = field2d(1:iana0,:,19)
 td2e(1:iana0,:)   = field2d(1:iana0,:,20)

 if (iana0 /= iana) then
   phige(iana,:)  = phige(1,:)
   psel(iana,:)   = psel(1,:)
   pmsle(iana,:)  = pmsle(1,:)
   cctote(iana,:) = cctote(1,:)
   u10e(iana,:)   = u10e(1,:)
   v10e(iana,:)   = v10e(1,:)
   t2e(iana,:)    = t2e(1,:)
   td2e(iana,:)   = td2e(1,:)
 endif

 if (ist==1) then

   if (level_type == 1) then
     print *, 'Input data on isobaric levels with the following values:'
     do k = 1,nlevp
       print '(i3,a2,f5.0,a4)',k,'  ',lev_list(k)*0.01,'e+02'
     enddo
     print *
   else
     print *,'Input data on hybrid levels:'
     print*, "actual levels available (left col.) and ECMWF level index (right col.)"
     do k = 1,nlevp
       print '(i3,a,i3)',k,'   ',int(lev_list(k))
     enddo
     print*
   endif

   print *,'Thickness and mid-level depth (cm) from the surface of the input soil layers:'
   depthsum = 0.
   do k = 1,nlevg_inp
     depth = 2.*(lev_list_soil(k) - depthsum)
     depthsum = depthsum + depth
     print '(i3,2f8.1)', k, depth*1.e2, lev_list_soil(k)*1.e2
   enddo
   print*

 endif ! ist=1

 qe(:,:,:)  = max(qe(:,:,:), zqcmin)
 qce(:,:,:) = max(qce(:,:,:), 0.)

! Definition of maximum and minumum volumetric soil content
! for input parameters (on the input data grid)
! Conversion of volumetric water soil content into relative soil content

 if (ist==1) then

! Parameters of soil in input data:
! NSTE - Number of soil texture types
! QGMINE - minimum volumetric water content (m**3/m**3)
! QGMAXE - maximum (saturation, porosity) volumetric water content (m**3/m**3)

! IFS-ECMWF global model

  iday = idate0(1)*10000+idate0(2)*100+idate0(3)

  if (iday >= 20070606) then

! After 06/06/2007:
! Soil texture types: 1 - Coarse, 2 - Medium, 3 - Medium-Fine,
!                     4 - Fine, 5 - Very Fine, 6 - Organic
!                     7 - unknown!

    nste = 7
    allocate(qgmine(nste))
    allocate(qgmaxe(nste))

    qgmine(1) = 0.000
    qgmine(2) = 0.000
    qgmine(3) = 0.000
    qgmine(4) = 0.000
    qgmine(5) = 0.000
    qgmine(6) = 0.000
    qgmine(7) = 0.000
    qgmaxe(1) = 0.403
    qgmaxe(2) = 0.439
    qgmaxe(3) = 0.430
    qgmaxe(4) = 0.520
    qgmaxe(5) = 0.614
    qgmaxe(6) = 0.766
    qgmaxe(7) = 0.590

    if (any(int(soile(:,:))==val_missing)) then
       write (*,*) 'Error - IFS_ECMWF data:'
       write (*,*) 'after 06/06/2007 parameter 043 (soil type) is required in input'
       stop 'Aborting.'
    endif

    do j = 1,jana
    do i = 1,iana
      iste = int(soile(i,j))
      if (iste>=1.and.iste<=nste) then
        qgmine2d(i,j) = qgmine(iste)
        qgmaxe2d(i,j) = qgmaxe(iste)
      else
        qgmine2d(i,j) = 0.
        qgmaxe2d(i,j) = 0.
        if (iste > nste) write (*,*)                                                       &
        "Caution: soil texture type in input file is outside the allowed range at point:", &
        i, j, iste
      endif
    enddo
    enddo

  else

! Before 06/06/2007 (DATE=2007-06-06):
! a unique soil texture type is used in IFS-ECMWF model: 1 - Loamy

    soile(:,:) = 0.
    nste = 1
    allocate(qgmine(nste))
    allocate(qgmaxe(nste))
    qgmine(1) = 0.000
    qgmaxe(1) = 0.472
    qgmine2d(:,:) = qgmine(1)
    qgmaxe2d(:,:) = qgmaxe(1)

  endif ! iday >= 20070606

 endif ! ist = 1

 if (surf_elaborate) then

   if (any(int(field2d(:,:,5))==int(val_missing))) then
     print*, "Caution: the tsurf field is missing in IFS-ECMWF input grib data"
     print*, "at instant =", ist
     print*, "But since it is not a mandatory field, prebolam does not stop"
     print*, "If instant > 1, try with NSFC=1 in prebolam.inp, which"
     print*, "avoids the analysis of surface variables after the first instant."
   else
     tsurfe(1:iana0,:) = field2d(1:iana0,:,5)
   endif

   ficeee(1:iana0,:) = min(max(field2d(1:iana0,:,22), 0.), 1.) ! Sets ficeee=0. in case of missing data
   ticeee(1:iana0,:) = field3d_soil(1:iana0,:,1,6)

   if (any(int(field2d(:,:,14))==int(val_missing))) then
      print*, "The snow cover field is missing in IFS-ECMWF input grib data,"
      print*, "therefore prebolam stops!"
      print*, "Instant =", ist
      print*, "If instant > 1, try with NSFC=1 in prebolam.inp, which"
      print*, "avoids the analysis of surface variables after the first instant."
      stop
   else
     snowe(1:iana0,:) = max(0., field2d(1:iana0,:,14)*1.e-3) ! Conversion from kg m**-2 (mm) into m of water
   endif

   if (any(int(field3d_soil(:,:,1:nlevg_inp,1))==int(val_missing))) then
      print*, "The soil temperature field is missing in IFS-ECMWF input grib data,"
      print*, "therefore prebolam stops!"
      print*, "Instant =", ist
      print*, "If instant > 1, try with NSFC=1 in prebolam.inp, which"
      print*, "avoids the analysis of surface variables after the first instant."
      stop
   else
      tge(1:iana0,:,1:nlevg_inp) = field3d_soil(1:iana0,:,1:nlevg_inp,1)
   endif

   if (any(int(field3d_soil(:,:,1:nlevg_inp,2))==int(val_missing))) then
      print*, "The soil moisture field is missing in IFS-ECMWF input grib data,"
      print*, "therefore prebolam is stops!"
      print*, "Instant =", ist
      print*, "If instant > 1, try with NSFC=1 in prebolam.inp, which"
      print*, "avoids the analysis of surface variables after the first instant."
      stop
   else
      qge(1:iana0,:,1:nlevg_inp) = max(0., field3d_soil(1:iana0,:,1:nlevg_inp,2))
   endif

   if (iana0 /= iana) then
     tge(iana,:,1:nlevg_inp) = tge(1,:,1:nlevg_inp)
     qge(iana,:,1:nlevg_inp) = qge(1,:,1:nlevg_inp)
   endif

   if (iana0 /= iana) then
     tsurfe(iana,:) = tsurfe(1,:)
     ficeee(iana,:) = ficeee(1,:)
     ticeee(iana,:) = ticeee(1,:)
     snowe(iana,:) = snowe(1,:)
   endif

   do j = 1,jana
   do i = 1,iana
     if (int(soile(i,j))>=1.and.int(soile(i,j))<=nste) then
       qge(i,j,1:nlevg_inp) = min(max((qge(i,j,1:nlevg_inp)-qgmine2d(i,j))/(qgmaxe2d(i,j)-qgmine2d(i,j)), 0.), 1.)
     else
       qge(i,j,1:nlevg_inp) = 0.
     endif
   enddo
   enddo

! QGE redefined over the sea to an average value to be used to defined soil moisture in small islands

   do j = 1,jana
   do i = 1,iana
    if(fmaske(i,j) >= 0.5) qge(i,j,:) = 0.33
   enddo
   enddo

! Ad hoc corrections of some fields because of problems found in input data

   day1 = day+365.*0.5
   if (day1>365.) day1 = day1-365.

   print*, "No correction applied to IFS-ECMWF soil data" ! No correction applied to IFS after july 2013!
   goto 567

! Correction of systematic error of soil temperature in ECMWF input data
! Verified only for European area on the basis of bias study in Climagri project
! Extended to all longitudes at N.H. latitudes only

! CAUTION: the following correction should be verified in time and
! for different geographical areas
! Last update July 2008 (Oxana Drofa), after 1 year verification with Central Europe stations
! Corrections based on analytical sinusoidal functions
! (describing annual and semi-annul periods - curve fitting of gnuplot used)

  if (ist==1) then
    nday = idate0(3)
    do j = 1,idate0(2)-1
      nday = nday + imon(j)
    enddo
    day = float(nday-1)+float(idate0(4))/24. ! day of the year as real (with hour precision)
  endif

  day1 = day+365.*0.5
  if (day1>365.) day1 = day1-365.

! Averages of T and Q over land before applying corrections

  tgav1 = 0.
  tgav2 = 0.
  tgav3 = 0.
  tgav4 = 0.
  qgav1 = 0.
  qgav2 = 0.
  qgav3 = 0.
  qgav4 = 0.
  npt = 0

  do j = 1,jana
  do i = 1,iana
    if(fmaske(i,j).lt.0.5) then
      npt = npt+1
      tgav1 = tgav1+tge(i,j,1)
      tgav2 = tgav2+tge(i,j,2)
      tgav3 = tgav3+tge(i,j,3)
      tgav4 = tgav4+tge(i,j,4)
      qgav1 = qgav1+qge(i,j,1)
      qgav2 = qgav2+qge(i,j,2)
      qgav3 = qgav3+qge(i,j,3)
      qgav4 = qgav4+qge(i,j,4)
    endif
  enddo
  enddo

  tgav1 = tgav1/float(npt)
  tgav2 = tgav2/float(npt)
  tgav3 = tgav3/float(npt)
  tgav4 = tgav4/float(npt)
  qgav1 = qgav1/float(npt)
  qgav2 = qgav2/float(npt)
  qgav3 = qgav3/float(npt)
  qgav4 = qgav4/float(npt)
  print*, "No. of points for averages", npt
  print*, "Averages of soil temper. (lev. 1 to 4) before corrections, input data:", tgav1, tgav2, tgav3, tgav4
  print*, "Averages of soil wetness (lev. 1 to 4) before corrections, input data:", qgav1, qgav2, qgav3, qgav4

  do j = 1,jana
  do i = 1,iana

    if (alat_inp(i,j) >= 0.) then
      day2 = day    ! N. Hemis.
    else
      day2 = day1   ! S. Hemis.: time shift of half year
    endif

    zflatt = exp(-0.002*(alat_inp(i,j)-50.)**2)

    zdtg(1) = 0.47-1.24*cos((day2-17.4)/365.*2.*pi)+0.28*cos((day2+29.9)/365.*4.*pi)
    zdtg(2) = 1.02-1.19*cos((day2- 3.9)/365.*2.*pi)+0.32*cos((day2+36.9)/365.*4.*pi)
    zdtg(3) = 1.04-0.80*cos((day2+16.3)/365.*2.*pi)+0.26*cos((day2+41.2)/365.*4.*pi)
    zdtg(4) = 1.01-0.62*cos((day2+58.1)/365.*2.*pi)+0.30*cos((day2+44.4)/365.*4.*pi)
    zdtg(1:4) = zdtg(1:4)*zflatt

    if(fmaske(i,j) < 0.5) tge(i,j,1:4) = tge(i,j,1:4) + zdtg(1:4)

    if (i.eq.iana/2.and.j.eq.jana/2) then
      print*, "Corrections for T of soil levels at a mid domain latitude:"
      do k=1,4
        print*, "   Lev.", k,"    dt =",zdtg(k)
      enddo
    endif

  enddo
  enddo
  print*, "Caution: T data of soil layers corrected!"

! Averages of T and Q over land after applying all corrections

   tgav1 = 0.
   tgav2 = 0.
   tgav3 = 0.
   tgav4 = 0.
   qgav1 = 0.
   qgav2 = 0.
   qgav3 = 0.
   qgav4 = 0.
   npt = 0
   do j = 1,jana
   do i = 1,iana
     if(fmaske(i,j).lt.0.5) then
       npt = npt+1
       tgav1 = tgav1+tge(i,j,1)
       tgav2 = tgav2+tge(i,j,2)
       tgav3 = tgav3+tge(i,j,3)
       tgav4 = tgav4+tge(i,j,4)
       qgav1 = qgav1+qge(i,j,1)
       qgav2 = qgav2+qge(i,j,2)
       qgav3 = qgav3+qge(i,j,3)
       qgav4 = qgav4+qge(i,j,4)
     endif
   enddo
   enddo
   tgav1 = tgav1/float(npt)
   tgav2 = tgav2/float(npt)
   tgav3 = tgav3/float(npt)
   tgav4 = tgav4/float(npt)
   qgav1 = qgav1/float(npt)
   qgav2 = qgav2/float(npt)
   qgav3 = qgav3/float(npt)
   qgav4 = qgav4/float(npt)
   print*, "Averages of soil temper. (lev. 1 to 4) after all corrections, input grid:", tgav1, tgav2, tgav3, tgav4
   print*, "Averages of soil wetness (lev. 1 to 4) after all corrections, input grid:", qgav1, qgav2, qgav3, qgav4

567 continue

 endif ! surf_elaborate

 return
 end
!=======================================================================
 subroutine conv_gfs_data(ist)

! Procedure to convert meteorological fields derived from input data files
! For GFS data
!------------------------------------------------------------------------

use param, only : iana0, iana, jana, nlev_atm_inp_max, nlev_soil_inp_max, nlevp, nlevg_inp, &
 eps, val_missing, pi, idate0, iperiod, iperiod_inp, surf_elaborate, res_change, g0
use gribana

implicit none

integer :: ist
integer :: i, j, k, khumini, nday, npt, idayi, idayf, iint
integer, dimension(3) :: iii
real, dimension(iana,jana,nlev_atm_inp_max) :: rh, wspeed
real :: wm, qsat, qsatw, qsati, esat, esatw, esati, eee, rhcwi, day1, day2, &
 zflatt, zflatq
real :: tgav1,tgav2,tgav3,tgav4,qgav1,qgav2,qgav3,qgav4,zcq1,zcq2,zcq3,zcq4
integer, save :: nste, iana_old, jana_old
real, save :: day
real, dimension(nlev_atm_inp_max) :: p
real, dimension(nlev_soil_inp_max) :: zdtg
integer, dimension(12) :: imon=(/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
real :: depth, depthsum, zf1, zf2

 res_change = 0
 if(ist > 1.and.(iana /= iana_old.or.jana /= jana_old)) then
   print*
   print*, "The GFS grid resolution has changed at instant: ist =", ist
   res_change = 1
 endif
 iana_old = iana
 jana_old = jana

 mask_frame(:,:) = 1

! Note: the definitions below are used only for diagnostics
! Date, time etc.

 if (iperiod_inp(1)==1) then         ! Time unit is hour
   iperiod(1) = iperiod_inp(2)/24                             ! Days
   iperiod(2) = iperiod_inp(2)-iperiod(1)*24                  ! Hours
   iperiod(3) = 0                                             ! Minutes
 elseif (iperiod_inp(1)==0) then     ! Time unit is minute
   iperiod(1) = iperiod_inp(2)/24/60                          ! Days
   iperiod(2) = (iperiod_inp(2)-iperiod(1)*24*60)/60          ! Hours
   iperiod(3) = iperiod_inp(2)-iperiod(1)*24*60-iperiod(2)*60 ! Minutes
 endif

! Definition of meteorological fields declared in main program
! LIST TO BE CHECKED ! see SUBROUTINE READ_GRIB2_DATA

 ge(1:iana0,:,:) = field3d(1:iana0,:,:,1)
 te(1:iana0,:,:) = field3d(1:iana0,:,:,2)
 ue(1:iana0,:,:) = field3d(1:iana0,:,:,3)
 ve(1:iana0,:,:) = field3d(1:iana0,:,:,4)
! qe(1:iana0,:,:) = field3d(1:iana0,:,:,5)
 rh(1:iana0,:,1:nlevp) = field3d(1:iana0,:,1:nlevp,6)

 if (iana0 /= iana) then
   ge(iana,:,:) = ge(1,:,:)
   te(iana,:,:) = te(1,:,:)
   ue(iana,:,:) = ue(1,:,:)
   ve(iana,:,:) = ve(1,:,:)
!   qe(iana,:,:) = qe(1,:,:)
   rh(iana,:,1:nlevp) = rh(1,:,1:nlevp)
 endif

 do k = 1,nlevp
   if (lev_list(k)>15000.) goto 100 ! for pressure levels > 150 hpa no modification is introduced
 enddo
100 khumini = k

 if (any(int(field3d(:,:,khumini:nlevp,7))==-9999).or.maxval(field3d(:,:,khumini:nlevp,7))<1.e-7) then
   iflag_cloude = 0
   qce(:,:,:) = 0.
   if(ist<=2) print*, "Cloud total water (liquid+ice) is not available in input"
 else
   iflag_cloude = 1
   qce(1:iana0,:,khumini:nlevp) = field3d(1:iana0,:,khumini:nlevp,7)
   if (iana0 /= iana) qce(iana,:,khumini:nlevp) = qce(1,:,khumini:nlevp)
   qce(:,:,1:khumini) = 0.   ! Definition of missing values of cloud water+ice at high levels
   if(ist<=2) print*, "Cloud total water (liquid+ice) is available in input"
 endif

! List of isobaric level:

   p(1:nlevp) = lev_list(1:nlevp)

  if (ist.le.2) then
   write(*,'(a,i3)') " Number of isobaric levels nlevp =", nlevp
   print *,'Isobaric levels for which T is available in input data:'
   do k = 1,nlevp
   print '(i4,a6,f5.0,a2)',k,'   p =',p(k)*0.01,'e2'
   enddo
   print*
  endif

! Auxiliary 2d data from 80 m level above ground (t, p, q, u, v)

  t80e(1:iana0,:)  = field2d(1:iana0,:,38)
  q80e(1:iana0,:)  = field2d(1:iana0,:,39)
  p80e(1:iana0,:)  = field2d(1:iana0,:,40)
  u80e(1:iana0,:)  = field2d(1:iana0,:,41)
  v80e(1:iana0,:)  = field2d(1:iana0,:,42)

  if (iana0 /= iana) then
   t80e(iana,:) = t80e(1,:)
   q80e(iana,:) = q80e(1,:)
   p80e(iana,:) = p80e(1,:)
   u80e(iana,:) = u80e(1,:)
   v80e(iana,:) = v80e(1,:)
  endif

 if (ist==1) then
   print *,'Thickness and mid-level depth (cm) from the surface of the input soil layers:'
   depthsum = 0.
   do k = 1,nlevg_inp
   depth = 2.*(lev_list_soil(k) - depthsum)
   depthsum = depthsum + depth
   print '(i3,2f8.1)', k, depth*1.e2, lev_list_soil(k)*1.e2
   enddo
   print*
 endif

 if (surf_elaborate) then

! Topography geopotential.
! For the GFS model, PHIGE2D and PHIGE are the same

   if (any(int(field2d(:,:,1))==int(val_missing))) then
     print*, "The topography field is missing in GFS input grib data,"
     print*, "but surface analysis has been required, therefore prebolam stops!"
     print*, "Instant =", ist
     print*, "If instant > 1, try with NSFC=1 in prebolam.inp, which"
     print*, "avoids the analysis of surface variables after the first instant."
     stop
   else
     phige2d(1:iana0,:) = field2d(1:iana0,:,1)
     phige(1:iana0,:)   = phige2d(1:iana0,:)
     if (iana0 /= iana) then
       phige2d(iana,:) = phige2d(1,:)
       phige(iana,:)   = phige(1,:)
     endif
   endif

 elseif (surf_elaborate.eqv..false..and.res_change == 1) then

   if (any(int(field2d(:,:,1))==int(val_missing))) then
     print*, "The topography field is missing in GFS input grib data"
     print*, "and GFS grid resolution has changed at instant =", ist
     print*, "In this case the (interpolated) input model topography of the previous instant"
     print*, "is used, which may imply some errors in surface or near surface fields."
   endif
   phige2d(1:iana0,:) = field2d(1:iana0,:,1) ! this sets values in case topography is missing
   phige(1:iana0,:)   = phige2d(1:iana0,:)   ! this sets values in case topography is missing
   if (iana0 /= iana) then
     phige2d(iana,:) = phige2d(1,:)
     phige(iana,:)   = phige(1,:)
   endif

 endif  ! surf_elaborate

 if (surf_elaborate) then

   if (any(int(field2d(:,:,3))==int(val_missing)).and.any(int(field2d(:,:,46))==int(val_missing))) then
     print*, "The land-sea mask field is missing in GFS input grib data,"
     print*, "but surface analysis has been required, therefore prebolam stops!"
     print*, "Instant =", ist
     print*, "If instant > 1, try with NSFC=1 in prebolam.inp, which"
     print*, "avoids the analysis of surface variables after the first instant."
     print*, "In case of change of GFS grid resolution with respect to the first"
     print*, "instant (analysis), the (interpolated) fmask of the analysis is used,"
     print*, "which may imply some errors in surface fields."
     stop
   endif

! Land-sea mask LANDN (coding discipline=2, category=0, parameter=218, index=46 here), valid
! after 19.07.2017, takes precedence over the "old" LAND (discipline=2, category=0,
! parameter=0, index=3 here), which is the standard land-sea mask to be used before 19.07.2017.

   if (any(int(field2d(:,:,3))/=int(val_missing))) then
     fmaske(1:iana0,:) = field2d(1:iana0,:,3)
   endif

   if (any(int(field2d(:,:,46))/=int(val_missing))) then
     fmaske(1:iana0,:) = field2d(1:iana0,:,46)
   endif

   if (iana0 /= iana) fmaske(iana,:) = fmaske(1,:)
   fmaske(:,:) = min(max(1.-fmaske(:,:),0.),1.) ! Conversion of FMASK: 0-sea, 1-land to 1-sea, 0-land

! Definition of tge and qge - check on missing data is made below, because it must depend on fmaske

   tge(1:iana0,:,1:nlevg_inp) = field3d_soil(1:iana0,:,1:nlevg_inp,1)
   qge(1:iana0,:,1:nlevg_inp) = field3d_soil(1:iana0,:,1:nlevg_inp,2)

! Check on the availability of soil temperature data over land (always missing over the sea)

   do j=1,jana
   do i=1,iana0
     if (fmaske(i,j).lt.1.e-7.and.tge(i,j,1).lt.0.) then
       print*, "The soil temperature field is missing in GFS input grib data"
       print*, "over land, therefore prebolam stops!"
       print*, "Instant =", ist
       print*, "If instant > 1, try with NSFC=1 in prebolam.inp, which"
       print*, "avoids the analysis of surface variables after the first instant."
!       stop
     endif
   enddo
   enddo

! Check on the availability of soil moisture data over land (always missing over the sea)

   do j=1,jana
   do i=1,iana0
     if (fmaske(i,j).lt.1.e-7.and.qge(i,j,1).lt.0.) then
       print*, "The soil moisture field is missing in GFS input grib data"
       print*, "over land, therefore prebolam stops!"
       print*, "Instant =", ist
       print*, "If instant > 1, try with NSFC=1 in prebolam.inp, which"
       print*, "avoids the analysis of surface variables after the first instant."
!       stop
     endif
   enddo
   enddo

! Check on the availability of surf temperature

   if (any(int(field2d(:,:,5))==int(val_missing))) then
     print*, "The tsurf field is missing in GFS input grib data."
     print*, "Tsurf is mandatory at the first instant for defining SST,"
     print*, "therefore prebolam stops!"
     print*, "Instant =", ist
     print*, "If instant > 1, try with NSFC=1 in prebolam.inp, which"
     print*, "avoids the analysis of surface variables after the first instant."
     stop
   else
     tsurfe(1:iana0,:)  = field2d(1:iana0,:,5)
   endif

! Check on the availability of snow depth over land (missing over non frozen sea)

   do j=1,jana
   do i=1,iana0
     if (fmaske(i,j).lt.1.e-7.and.(int(field2d(i,j,14))==int(val_missing))) then
       print*, "The snow cover field is missing in GFS input grib data,"
       print*, "therefore prebolam stops!"
       print*, "Instant =", ist
       print*, "If instant > 1, try with NSFC=1 in prebolam.inp, which"
       print*, "avoids the analysis of surface variables after the first instant."
!       stop
       snowe(i,j) = 0.
     else
       snowe(i,j) = max(0., field2d(i,j,14)*1.e-3) ! Conversion from kg m**-2 (mm) into m of water
     endif
   enddo
   enddo

   ficeee(1:iana0,:)  = min(max(field2d(1:iana0,:,22), 0.), 1.) ! Sets ficeee=0. in case of missing data
   icethee(1:iana0,:) = max(field2d(1:iana0,:,26), 0.)

   if (iana0 /= iana) then
     tge(iana,:,1:nlevg_inp) = tge(1,:,1:nlevg_inp)
     qge(iana,:,1:nlevg_inp) = qge(1,:,1:nlevg_inp)
     tsurfe(iana,:)  = tsurfe(1,:)
     snowe(iana,:)   = snowe(1,:)
     ficeee(iana,:)  = ficeee(1,:)
     icethee(iana,:) = icethee(1,:)
   endif

! Definition of temperature and soil water content in the case of missing data (-9999) over the sea
! (TGE and QGE are not defined over sea in GFS data)

   do k = 1,nlevg_inp
   do j = 1,jana
   do i = 1,iana
     if (int(tge(i,j,k))==int(val_missing)) tge(i,j,k) = tsurfe(i,j)
     if (int(qge(i,j,k))==int(val_missing)) qge(i,j,k) = 0.33 ! to be used to define qg on small islands
   enddo
   enddo
   enddo

! Definition of maximum and minimum volumetric soil content
! for input parameters (on the input data grid)
! Conversion of volumetric water soil content into relative soil content

! Parameters of soil in input data:
! NSTE - Number of soil texture types
! QGMINE - minimum volumetric water content (m**3/m**3)
! QGMAXE - maximum (saturation, porosity) volumetric water conent (m**3/m**3)

! GFS (NOAA-NCEP) global model
! 0.470 is the max. water soil content ((m**3/m**3), this value has been
! empirically verified by data during period 01.06.2008-30.06.2009

   nste = 1
   if(.not.allocated(qgmine)) allocate(qgmine(nste))
   if(.not.allocated(qgmaxe)) allocate(qgmaxe(nste))
   qgmine(1) = 0.000
   qgmaxe(1) = 0.470

   qgmine2d(:,:) = qgmine(1)
   qgmaxe2d(:,:) = qgmaxe(1)

   do k = 1,nlevg_inp
     qge(:,:,k) = (qge(:,:,k)-qgmine2d(:,:))/(qgmaxe2d(:,:)-qgmine2d(:,:))
   enddo

 endif ! surf_elaborate

! Relative humidity is correctly defined only up to 100-150 hPa only

 do k = 1,nlevp
! if(p(k).lt.150.e2) rh(:,:,k) = 4.
 if(p(k).le.100.e2) rh(:,:,k) = 4.
 if(p(k).le.70.e2)  rh(:,:,k) = 2.
 if(p(k).le.30.e2)  rh(:,:,k) = 1.
 if(p(k).le.20.e2)  rh(:,:,k) = .5
 if(p(k).le.10.e2)  rh(:,:,k) = .05
 if(p(k).le. 3.e2)  rh(:,:,k) = .004
 if(p(k).le. 1.e2)  rh(:,:,k) = .0003
 enddo

! Conversion of relative humidity into specific humidity

 do k = 1,nlevp
 do j = 1,jana
 do i = 1,iana
 rh(i,j,k) = max(rh(i,j,k),0.)
 rh(i,j,k) = min(rh(i,j,k),102.)
 call qsat_tetens(te(i,j,k), p(k), eps, qsat, qsatw, qsati, esat, esatw, esati)
 eee = esat*rh(i,j,k)*1.e-2
 qe(i,j,k) = eps*eee/(eps*eee-eee+p(k))
 qe(i,j,k) = max(qe(i,j,k), 1.e-7)

! Reduction of cloud water+ice

 if(p(k).lt.280.e2) then
 qce(i,j,k) = qce(i,j,k)*(p(k)/280.e2)**2 ! Reduction of cloud water/ice at high levels
 endif
 qce(i,j,k) = max(qce(i,j,k), 0.)
 qce(i,j,k) = min(qce(i,j,k),1.5e-3)   ! to avoid excessive condensate (further reduction in DEFCLI in BOLAM)
 if(rh(i,j,k) > 0..and.rh(i,j,k) <= 60.) qce(i,j,k) = 0.  ! to avoid condensate at low relative humidity
 enddo
 enddo
 enddo

! --- Ad hoc corrections of some fields because of problems found in input data ----

 if (ist==1) then
   nday = idate0(3)
   do j = 1,idate0(2)-1
     nday = nday + imon(j)
   enddo
   day = float(nday-1)+float(idate0(4))/24. ! day of the year as real (with hour precision)
 endif

 day1 = day+365.*0.5
 if (day1>365.) day1 = day1-365.

 if (surf_elaborate) then

! Averages of T and Q over land before applying corrections

    tgav1 = 0.
    tgav2 = 0.
    tgav3 = 0.
    tgav4 = 0.
    qgav1 = 0.
    qgav2 = 0.
    qgav3 = 0.
    qgav4 = 0.
    npt = 0
    do j = 1,jana
    do i = 1,iana
    if(fmaske(i,j).lt.0.5) then
    npt = npt+1
    tgav1 = tgav1+tge(i,j,1)
    tgav2 = tgav2+tge(i,j,2)
    tgav3 = tgav3+tge(i,j,3)
    tgav4 = tgav4+tge(i,j,4)
    qgav1 = qgav1+qge(i,j,1)
    qgav2 = qgav2+qge(i,j,2)
    qgav3 = qgav3+qge(i,j,3)
    qgav4 = qgav4+qge(i,j,4)
    endif
    enddo
    enddo
    tgav1 = tgav1/float(npt)
    tgav2 = tgav2/float(npt)
    tgav3 = tgav3/float(npt)
    tgav4 = tgav4/float(npt)
    qgav1 = qgav1/float(npt)
    qgav2 = qgav2/float(npt)
    qgav3 = qgav3/float(npt)
    qgav4 = qgav4/float(npt)
    write(*,'(a,i6)') " No. of points for averages:", npt
    print*, "Averages of soil temper. (lev. 1 to 4) before corrections, input data:", tgav1, tgav2, tgav3, tgav4
    print*, "Averages of soil wetness (lev. 1 to 4) before corrections, input data:", qgav1, qgav2, qgav3, qgav4

    do j = 1,jana
    do i = 1,iana

! Correction of apparent overestimation of water soil content in some areas:
! Calabria, Sicilia, Sardegna

! Sardegna

      if (alon_inp(i,j)>8..and.alon_inp(i,j)<10..and.alat_inp(i,j)>38.8.and.alat_inp(i,j)<41.15) then
        qge(i,j,1)   = qge(i,j,1)*.90*(1.-fmaske(i,j))
        qge(i,j,2:4) = qge(i,j,2:4)*.80*(1.-fmaske(i,j))
      endif

! Sicilia e Calabria

      if (alon_inp(i,j)>12..and.alon_inp(i,j)<17.2.and.alat_inp(i,j)>36.5.and.alat_inp(i,j)<39.1) then
        qge(i,j,1)   = qge(i,j,1)*.90*(1.-fmaske(i,j))
        qge(i,j,2:4) = qge(i,j,2:4)*.80*(1.-fmaske(i,j))
      endif

! Correction of the bias of soil T and Qrel (levels 1 to 4)

! CAUTION: the following corrections should be verified in time
! Last update Aug. 2009 (Oxana Drofa), after 1 year verification of soil T using
! Central Europe stations (July 2008) and comparison between GFS and ECMWF values
! (August 2009)
! Corrections for soil T are based on differences between GFS and observed values
! over Central Europe, and extended to other areas based on differences between GFS
! and ECMWF values (using data in the period 01.06.2008-30.06.2009)
! Corrections for soil moisture are based only on differences between GFS and ECMWF
! values (using data 01.06.2008-30.06.2009), assuming that ECMWF mean values are correct
! Corrections for T and water have been generalized to all latitudes (incl. SH) after
! bias was verified at different geographical zones
! Corrections in time based on analytical sinusoidal functions
! (describing annual and semi-annual periods - curve fitting of gnuplot used) -
! the correction amplitudes are the same for NH and SH, with time dependency anti-symmetric
! From comparison of soil T between ECMWF and GFS data in 2012 over Europe and Med., it seems that at
! least in the NH the sign of correction might even change over north Africa - therefore reduced
! correction (only a constant part - see below zdtg) is applied between about 35 deg. N and S

      if (alat_inp(i,j) >= 0.) then
        day2 = day    ! N. Hemis.
      else
        day2 = day1   ! S. Hemis.: time shift of half year
      endif

! Soil temperature

!      zflatt  = 1.-exp(-(alat_inp(i,j)/30.)**2.)
      zflatt  = 1.-exp(-(alat_inp(i,j)/32.)**2.)  ! July 2016

      if(abs(alat_inp(i,j)).lt.39.) zflatt = zflatt*max((abs(alat_inp(i,j))-34.)/5., 0.)

      zdtg(1) = 0.5+zflatt*(0.25-0.63*cos((day2-11.0)/365.*2.*pi)+0.67*cos((day2+18.3)/365.*4.*pi))
      zdtg(2) = 1.3+zflatt*(1.19-3.77*cos((day2+ 8.1)/365.*2.*pi)+1.08*cos((day2+ 9.6)/365.*4.*pi))
      zdtg(3) = 0.8+zflatt*(1.63-4.11*cos((day2+ 3.3)/365.*2.*pi)+0.84*cos((day2+ 0.7)/365.*4.*pi))
      zdtg(4) = 0.8+zflatt*(1.51-3.97*cos((day2- 2.7)/365.*2.*pi)+0.61*cos((day2-10.5)/365.*4.*pi))

! Revision of T corrections against ECMWF soil (July 2013)

!      zdtg(1) = 1.0* zdtg(1)
!      zdtg(2) = 0.96*zdtg(2)
!      zdtg(3) = 0.94*zdtg(3)
!      zdtg(4) = 0.92*zdtg(4)

! Revision of T corrections after comparing BOLAM equil. values vs GFS (Aug. 2014)

!      zdtg(1) = 0.7 *zdtg(1)
!      zdtg(2) = 0.95*zdtg(2)
!      zdtg(3) = 0.93*zdtg(3)
!      zdtg(4) = 0.91*zdtg(4)

! Revision of T corrections after comparing GFS soil T in 2015 vs 2014 and 2012 (May. 2015 - in addition to the one above))

!      zdtg(2) = 0.925*zdtg(2)
!      zdtg(3) = 0.94 *zdtg(3)
!      zdtg(4) = 0.95 *zdtg(4)

! Revision of T corrections after comparing BOLAM equil. values vs GFS (July 2016 - cumulative over the two above)
! Correction of atmospheric layers temp. as a function of zdtg eliminated (July 2016)

      zdtg(1) = 0.50*zdtg(1)
      zdtg(2) = 0.88*zdtg(2)
      zdtg(3) = 0.87*zdtg(3)
      zdtg(4) = 0.86*zdtg(4)
      if(alat_inp(i,j).lt.34.) zdtg = 0.  ! July 2016: GFS soil T seems too high over Sahara

      tge(i,j,1:4) = tge(i,j,1:4)+zdtg(1:4)*(1.-fmaske(i,j))

      if (i.eq.iana.and.j.eq.jana) then
      print*, "Corrections for T of soil levels at northernmost latitude of the GFS domain:"
        do k = 1,4
          write(*,'(a,i4,a,f10.5)') "   Lev.", k,"    dt =",zdtg(k)
        enddo
      endif

      if (i.eq.iana/2.and.j.eq.jana/2) then
      print*, "Corrections for T of soil levels at mid latitude of the GFS domain:"
        do k = 1,4
          write(*,'(a,i4,a,f10.5)') "   Lev.", k,"    dt =",zdtg(k)
        enddo
      endif

      if (i.eq.1.and.j.eq.1) then
      print*, "Corrections for T of soil levels at southernmost latitude of the GFS domain:"
        do k = 1,4
          write(*,'(a,i4,a,f10.5)') "   Lev.", k,"    dt =",zdtg(k)
        enddo
      endif

! Correction of soil moisture depend. on latitude and season (larger reduction in summer, with also
! a latitudinal seasonal shift - but intended here only for the North Africa-Europe area)

      zf1 = cos((day2-35.0)/365.*2.*pi)
      zf2 = cos((day -35.0)/365.*2.*pi)

! reduction of drying correction, larger at mid-high lat. (aug. 2014)

      zflatq = 0.5*(1.-0.65*(0.9-0.1*zf1)*exp(-((abs(alat_inp(i,j))-27.+2.5*(1.+zf2))/ 5.)**2.))  &
            +  0.5*(1.-0.65*(0.9-0.1*zf1)*exp(-((abs(alat_inp(i,j))-29.+2.5*(1.+zf2))/12.)**2.))
      zflatq = 0.5*zflatq + 0.5

      if(alat_inp(i,j).gt.0.) then  ! if uncommented, excludes moisture correction in the South. Hemis.
        qge(i,j,1) = zflatq*qge(i,j,1)
        qge(i,j,2) = zflatq*qge(i,j,2)
        qge(i,j,3) = zflatq*qge(i,j,3)
        qge(i,j,4) = zflatq*qge(i,j,4)
      endif

    enddo
    enddo

    qge(:,:,:) = max(min(qge(:,:,:),1.),0.) ! reset to min-max of relative values

    print*, "T of soil layers and lowest atmospheric layers over land corrected"
    print*, "Q of soil layers corrected"

! Averages of T and Q over land after applying all corrections

    tgav1 = 0.
    tgav2 = 0.
    tgav3 = 0.
    tgav4 = 0.
    qgav1 = 0.
    qgav2 = 0.
    qgav3 = 0.
    qgav4 = 0.
    npt = 0

    do j = 1,jana
    do i = 1,iana
      if(fmaske(i,j).lt.0.5) then
        npt = npt+1
        tgav1 = tgav1+tge(i,j,1)
        tgav2 = tgav2+tge(i,j,2)
        tgav3 = tgav3+tge(i,j,3)
        tgav4 = tgav4+tge(i,j,4)
        qgav1 = qgav1+qge(i,j,1)
        qgav2 = qgav2+qge(i,j,2)
        qgav3 = qgav3+qge(i,j,3)
        qgav4 = qgav4+qge(i,j,4)
      endif
    enddo
    enddo

    tgav1 = tgav1/float(npt)
    tgav2 = tgav2/float(npt)
    tgav3 = tgav3/float(npt)
    tgav4 = tgav4/float(npt)
    qgav1 = qgav1/float(npt)
    qgav2 = qgav2/float(npt)
    qgav3 = qgav3/float(npt)
    qgav4 = qgav4/float(npt)
    print*, "Averages of soil temper. (lev. 1 to 4) after all corrections, input grid:", tgav1, tgav2, tgav3, tgav4
    print*, "Averages of soil wetness (lev. 1 to 4) after all corrections, input grid:", qgav1, qgav2, qgav3, qgav4

 endif  ! surf_elaborate

! Comput. of max. wind speed (only for printout)

! wspeed(:,:,1:nlevp)=sqrt(ue(:,:,1:nlevp)**2+ve(:,:,1:nlevp)**2)
! wm=maxval(wspeed)
! iii=maxloc(wspeed)
! print*, "max wind speed in gfs domain", wm
! print*, "at grid point",iii(1),iii(2),iii(3)

return
end

! =============================================================
 subroutine conv_bolam_data(ist)

! Procedure to convert meteorological fields derived from input files of MHF-BOLAM type
! for BOLAM or GLOBO models
!------------------------------------------------------------------------

use param, only : iana0, iana, jana, nlevp, nlevg, nlevg_inp, &
 dxa, dya, x0a, y0a, x1a, y1a, nst, nvt, surf_elaborate, idate0
use gribana

implicit none

integer :: ist, i, j, k
real, dimension(iana,jana) :: fmaske_work
real :: depth, depthsum
real,dimension(iana, jana, nlevg_inp) :: soil_qmax_inp, soil_qmin_inp

 mask_frame(:,:) = 1
 frame = .false.

! Definition of meteorological fields declared in main program
! LIST TO BE CHECKED ! see subroutine read_grib2_data

 te (1:iana0,:,1:nlevp) = field3d(1:iana0,:,1:nlevp,2)
 ue (1:iana0,:,1:nlevp) = field3d(1:iana0,:,1:nlevp,3)
 ve (1:iana0,:,1:nlevp) = field3d(1:iana0,:,1:nlevp,4)
 qe (1:iana0,:,1:nlevp) = field3d(1:iana0,:,1:nlevp,5)
 qce(1:iana0,:,1:nlevp) = field3d(1:iana0,:,1:nlevp,7)

 if (iana0 /= iana) then
   te (iana,:,1:nlevp) = te (1,:,1:nlevp)
   ue (iana,:,1:nlevp) = ue (1,:,1:nlevp)
   ve (iana,:,1:nlevp) = ve (1,:,1:nlevp)
   qe (iana,:,1:nlevp) = qe (1,:,1:nlevp)
   qce(iana,:,1:nlevp) = qce(1,:,1:nlevp)
 endif

 qe(:,:,:)  = max(qe(:,:,:), 0.)
 qce(:,:,:) = max(qce(:,:,:), 0.)

 psel(1:iana0,:)    = log(field2d(1:iana0,:,4))

 if (iana0 /= iana) then
   psel(iana,:)    = psel(1,:)
 endif

 if (ist==1) then

   phige2d(1:iana0,:) = field2d(1:iana0,:,2)
   phige(1:iana0,:)   = field2d(1:iana0,:,2)
   fmaske(1:iana0,:)  = field2d(1:iana0,:,3) ! 1-sea, 0-land

   if (iana0 /= iana) then
     phige2d(iana,:) = phige2d(1,:)
     phige(iana,:)   = phige(1,:)
     fmaske(iana,:)  = fmaske(1,:)
   endif

   if (level_type == 1) then
     print *, 'Input data on isobaric levels with the following values:'
     do k=1,nlevp
       print '(i3,a2,f5.0,a4)',k,'  ',lev_list(k)*0.01,'e+02'
     enddo
     print *
   else
     print *,'Input data on hybrid levels:'
     print*, "actual levels available (left col.) and GLOBO/BOLAM level index (right col.)"
     do k=1,nlevp
       print '(i3,a,i3)',k,'   ',int(lev_list(k))
     enddo
     print*
   endif

   print *,'Depth (cm) of soil levels from the surface of the input soil layers:'
   depthsum = 0.
   do k=1,nlevg_inp
     depth = 2.*(lev_list_soil(k) - depthsum)
     depthsum = depthsum + depth
     print '(i3,2f8.1)', k, lev_list_soil(k)*1.e2
   enddo
   print*

 endif

 if (surf_elaborate) then

   snowe(1:iana0,:)   = field2d(1:iana0,:,14)
   tsurfe(1:iana0,:)  = field2d(1:iana0,:,5)
   qvsurfe(1:iana0,:)  = field2d(1:iana0,:,25)
   ficeee(1:iana0,:)  = min(max( field2d(1:iana0,:,22), 0.), 1.)
   icethee(1:iana0,:) = max(field2d(1:iana0,:,26), 0.)

   tge(1:iana0,:,1:nlevg_inp) = field3d_soil(1:iana0,:,1:nlevg_inp,1)
   qge(1:iana0,:,1:nlevg_inp) = field3d_soil(1:iana0,:,1:nlevg_inp,2)
   soil_qmax_inp(1:iana0,:,1:nlevg_inp) = field3d_soil(1:iana0,:,1:nlevg_inp,7)
   soil_qmin_inp(1:iana0,:,1:nlevg_inp) = field3d_soil(1:iana0,:,1:nlevg_inp,8)

   if (iana0 /= iana) then

     snowe(iana,:)   = snowe(1,:)
     tsurfe(iana,:)  = tsurfe(1,:)
     qvsurfe(iana,:)  = qvsurfe(1,:)
     ficeee(iana,:)  = ficeee(1,:)
     icethee(iana,:) = icethee(1,:)

     tge(iana,:,1:nlevg_inp) = tge(1,:,1:nlevg_inp)
     qge(iana,:,1:nlevg_inp) = qge(1,:,1:nlevg_inp)
     soil_qmax_inp(iana,:,1:nlevg_inp) = soil_qmax_inp(1,:,1:nlevg_inp)
     soil_qmin_inp(iana,:,1:nlevg_inp) = soil_qmin_inp(1,:,1:nlevg_inp)

   endif

! In GLOBO or BOLAM in water body and glaciers qg is 1., not 0.

   do j = 1,jana
   do i = 1,iana

     if (qge(i,j,1) > 0.99) qge(i,j,1:nlevg_inp)=0.

     if (fmaske(i,j) < 0.5) then
       qge(i,j,1:nlevg_inp) = min(max((qge(i,j,1:nlevg_inp)-soil_qmin_inp(i,j,1:nlevg_inp)) &
                           /(soil_qmax_inp(i,j,1:nlevg_inp)-soil_qmin_inp(i,j,1:nlevg_inp)), 0.), 1.)
      else
       qge(i,j,1:nlevg_inp) = 0.
     endif

! Redefinition of qg over input model sea to a value that defines
! deep soil moisture for small islands appearing in open seas

     if (fmaske(i,j) >= 0.5.and.qge(i,j,1) < 0.01) qge(i,j,1:nlevg_inp)=0.50

   enddo
   enddo

 endif ! surf_elaborate

return
end
!=======================================================================
 subroutine sst_isac(ssti,xxtg,yytg,nlon,nlat,ntot,dlon,dlat,alon0,alat0,x0d,y0d,idate0,ifl)

! Defines SST on the Bolam grid, using matrix of data SST ISAC-MYOCEAN on the Mediterranean,
! Black Sea and a portion of near Atlantic, combining it with the ECMWF or GFS SST analysis.
! SST data are in deg. C /100, 2 byte integers. Missing values: < -999
! (including the ext. framework).
! Lat.-lon. regular grid, step 1/16 deg (0.0625 deg).
! Grid dimensions: 1441x721 (not read in the sst file - must be written explicitely below).
! Grid extremes: lat 20.0, 65.0; lon -35.0, 55.0 (not read in the sst file: min. lat and lon
! are written explicitely below).

parameter(im=1441,jm=721)    ! dimensions of the input grid
integer*2 sst_med(im,jm)
real, dimension(im,jm) :: sst_ext, fmask_ext, workf, fmask_w
real xe_ext(im), ye_ext(jm)
real, dimension(nlon,nlat) :: xxtg, yytg, ssti, fmaski, work, work2
integer, dimension(5) :: idate0
character*12 filesst
character*30 str
integer, dimension(12) :: imon=(/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/) ! odd years not considered!
integer today, yesterday, day_bef_yesterday, flag
logical ex

ifl = 0

! Reading of SST data file (file name format: yyyymmdd.dat - if yyyymmdd_l.dat, it must be renamed)

filesst(9:12) = ".dat"

! SST files are taken into account to define the sst only if the data is today or yesterday or
! day before yesterday, otherwise the sst of the input model is not modified

! idate0(3): day of month; idate0(2): month; idate0(1): year (4 digits)

if(mod(idate0(1),4) == 0) imon(2)=29 ! to take into account (most) odd years
today = idate0(3) + idate0(2)*100 + idate0(1)*10000
!print*, "today", today
write (filesst(1:8), '(i8.8)') today
inquire(file=filesst, exist=ex)
if(ex.eqv..true.) then           ! case of sst file of today
open (10, file=filesst, status='old', form='unformatted', err=105)
else
yesterday = idate0(3)-1 + idate0(2)*100 + idate0(1)*10000
if(idate0(2) == 1.and.idate0(3) == 1) yesterday = 31 + 12*100 + (idate0(1)-1)*10000 ! case 1st Jan
do iii = 2,12
if(idate0(2) == iii.and.idate0(3) == 1) yesterday = imon(iii-1) + (iii-1)*100 + idate0(1)*10000 ! feb. to dec.
enddo
!print*, "yesterday", yesterday
write (filesst(1:8), '(i8.8)') yesterday
inquire(file=filesst, exist=ex)
if(ex.eqv..true.) then           ! case of sst file of yesterday
open (10, file=filesst, status='old', form='unformatted', err=105)
else
day_bef_yesterday = idate0(3)-2 + idate0(2)*100 + idate0(1)*10000
if(idate0(2) == 1.and.idate0(3) == 1) day_bef_yesterday = 30 + 12*100 + (idate0(1)-1)*10000 ! case 1st Jan
if(idate0(2) == 1.and.idate0(3) == 2) day_bef_yesterday = 31 + 12*100 + (idate0(1)-1)*10000 ! case 2nd Jan
do iii = 2,12
if(idate0(2) == iii.and.idate0(3) == 1) day_bef_yesterday = imon(iii-1)-1 + (iii-1)*100 + idate0(1)*10000 ! feb. to dec.
if(idate0(2) == iii.and.idate0(3) == 2) day_bef_yesterday = imon(iii-1) + (iii-1)*100 + idate0(1)*10000 ! feb. to dec.
enddo
!print*, "day_bef_yesterday", day_bef_yesterday
write (filesst(1:8), '(i8.8)') day_bef_yesterday
inquire(file=filesst, exist=ex)
if(ex.eqv..true.) then           ! case of sst file of day before yesterday
open (10, file=filesst, status='old', form='unformatted', err=105)
else
print*
print*, "File with SST ISAC-MYOCEAN not found or date not suitable."
print*

return
endif
endif
endif

read (10) sst_med
print*
print*, "File ", filesst, " with the ISAC-MYOCEAN SST read, to be used to define SST."
print*

do j = 1,jm
do i = 1,im
if(sst_med(i,j) < -50) then
sst_med(i,j) = -9999
fmask_ext(i,j) = 0.
else
fmask_ext(i,j) = 1.
endif
enddo
enddo

! Comput. of average sst in valid points, to be attributed to invalid points
! (land or boundaries or not covered area in north Atlantic)

t_aver = 0
icount = 0
do j = 1,jm
do i = 1,im
if(fmask_ext(i,j) >= 0.5) then
t_aver = t_aver + sst_med(i,j)
icount = icount + 1
endif
enddo
enddo
t_aver = t_aver/float(icount)
print*, "Average SST in file ", filesst, " in C:", t_aver/100.
print*

do j = 1,jm
do i = 1,im
if(fmask_ext(i,j) < 0.5) sst_med(i,j) = t_aver
enddo
enddo

! T in degrees Kelvin

sst_ext = float(sst_med)/100. + 273.15
tav = t_aver/100. + 273.15

! Extension of SST values towards the land or to invalid areas

workf = sst_ext
call seatemp(workf,fmask_ext,sst_ext,im,jm,6,2,0.6)

workf = fmask_ext
call seatemp(workf,fmask_ext,fmask_w,im,jm,6,1,1.)  ! define fmask_w as an extended fmask towards land or invalid areas

! Interpolation on the Bolam grid

do i = 1,im
xe_ext(i) = -35.+(i-1)*0.0625
enddo

do j = 1,jm
ye_ext(j) = 20.+(j-1)*0.0625
enddo

! Check of geographical coordinates: the sst extended grid must cover the full Bolam domain,
! otherwise it is not used (the check is based on the above chosen coordinates of the extended domain)

if(minval(xxtg) < minval(xe_ext).or.maxval(xxtg) > maxval(xe_ext).or.   &
   minval(yytg) < minval(ye_ext).or.maxval(yytg) > maxval(ye_ext)) then
print*, "Data of the ISAC-MYOCEAN SST file do not cover the Bolam grid"
print*, " and therefore are not used to define SST"
print*, "Min. Bolam long.", minval(xxtg), "Max. Bolam long.", maxval(xxtg)
print*, "Min. Bolam lat.", minval(yytg), "Max. Bolam lat.", maxval(yytg)
print*, "Min. SST data long.", minval(xe_ext), "Max. SST data long.", maxval(xe_ext)
print*, "Min. SST data lat.", minval(ye_ext), "Max. SST data lat.", maxval(ye_ext)
print*
return
endif

call interp_spline_2d(sst_ext,im,jm,xe_ext,ye_ext,xxtg,yytg,ntot,work,0.8)
call interp_spline_2d(fmask_w,im,jm,xe_ext,ye_ext,xxtg,yytg,ntot,fmaski,0.8) ! fmaski computed from interpol. of the extended fmask_w
fmaski = min(1., fmaski)
fmaski = max(0., fmaski)

!str = 'work1'
!call plotout(work-273.15,fmaski,nlon,nlat,str,,99)

call smooth_soil(work,work2,nlon,nlat,0.5,2)

!str = 'work2'
!call plotout(work-273.15,fmaski,nlon,nlat,str,99)
!str = 'ssti'
!call plotout(ssti-273.15,fmaski,nlon,nlat,str,99)

! Combination of ssti computed with GFS analysis data and of sst computed with observation data
! wei is the weight to be attributed to the MY-OCEAN product.
! Note that matrix work does not contain values useful to define temp. of lakes inland.
! Note also that fmaski is different from fmask (so 0.7 is used in place of 0.5).

wei = 0.88
do jlat = 1,nlat
do jlon = 1,nlon
if(fmaski(jlon,jlat) > 0.7) ssti(jlon,jlat) = (1.-wei)*ssti(jlon,jlat) + wei*work(jlon,jlat)
enddo
enddo
!str = 'ssti'
!call plotout(ssti-273.15,fmaski,nlon,nlat,str,99)
close (10)
ifl = 1
return

105 print*, "File with SST ISAC-MYOCEAN existing but unreadable."
return
end subroutine sst_isac
! ======================================================================
      subroutine ccloud

! Definition of cloud cover parameterized as a function of
! relative humidity and cloud ice+cloud water
! Weight of zrh tuned on GFS analyses) - no reduction as a function of static stability

  use param
  use mod_fields

      real zneb(nlev),zb1(nlev),zsigint(nlev),huc(nlev)

      yliv = 2834170.5
      yliw = 333560.5
      ylwv = yliv-yliw
      zeps = 1.e-6
      qccrit = 2.2e-4

! Definition of critical rel. hum. for cloud formation
! depending on level (do not set huc = 1. to avoid overflow)

      hlow = 0.88
      hmid = 0.87
      hhig = 0.92
      slow = 0.900
      smid = 0.550
      shig = 0.250
      zalf = (hlow-hmid)/(1.-slow)
      zbet = hlow-zalf
      zalf2 = (hmid-hhig)/(smid-shig)
      zbet2 = hmid-zalf2*smid

      zsigint(1) = .5*s(1)
      do k = 2, nlev
      zsigint(k) = .5*(s(k)+s(k-1))
      enddo

      do jklev = 1, nlev
       if(zsigint(jklev).lt.slow.and.zsigint(jklev).gt.smid) then
       huc(jklev) = hmid
       elseif(zsigint(jklev).ge.slow) then
       huc(jklev) = zalf * zsigint(jklev) + zbet
       elseif(zsigint(jklev).le.smid.and.zsigint(jklev).gt.shig) then
       huc(jklev) = zalf2 * zsigint(jklev) + zbet2
       else
       huc(jklev) = hhig
       endif
      enddo

      do jlat = 1, nlat
      do jlon = 1, nlon

      do jklev = nlev,2,-1
      zppp = p0*zsigint(jklev)-(p0-ps(jlon,jlat))*zsigint(jklev)**3*(1.+alfa*(1.-zsigint(jklev))*(2.-zsigint(jklev)))
      ztmp = t(jlon,jlat,jklev)
      call comp_esk(zesat,zqs1,ztmp,zppp,2)  ! blended saturation
      call comp_esk(zesat,zqs2,ztmp,zppp,3)  ! sat. to water below 0
      zqs = 0.70*zqs1+0.30*zqs2              ! ad hoc to limit contrib. to high clouds
      zrh1 = min(1., q(jlon,jlat,jklev)/zqs)
      zrh = (zrh1-huc(jklev))/(1.-huc(jklev))
      zrh = min(1.-zeps, zrh)
      zrh = max(zeps,    zrh)
!      zneb(jklev) = (0.07*zrh + qc(jlon,jlat,jklev)/qccrit)**.7
      zneb(jklev) = (0.5*zrh + qc(jlon,jlat,jklev)/qccrit)**.7
      zneb(jklev) = max(0., min(1.,zneb(jklev)))
      enddo

      zb1(2) = 1.-zneb(2)
      do jklev = 3,nlev
      znebmax = max(zneb(jklev),zneb(jklev-1))
      if (1.-zneb(jklev-1) /= 0.) then
      zb1(jklev) = (1.-znebmax)/(1.-zneb(jklev-1))
      else
      zb1(jklev) = 0.
      endif
      enddo

      zprod = 1.
      do jklev = 2, nlev
      zprod = zprod*zb1(jklev)
      enddo

      cloud (jlon,jlat) = max(0., 1.-zprod)
      cloud (jlon,jlat) = min(cloud (jlon,jlat),1.-zeps)

      enddo
      enddo

      return
      end
! ======================================================================
subroutine read_bolam_mhf_data(ifl,filename_atm,filename_soil,data_mode,ini_flag,                  &
    nlon_inp,nlat_inp,nlev_atm_inp,nlev_atm_inp_max,nlev_soil_inp,nlev_soil_inp_max,     &
    x0a,y0a,x1a,y1a,dxa,dya,idate,iperiod,lev_list,lev_list_soil,level_type,     &
    npar3d,npar2d,npar3d_soil,field3d,field2d,field3d_soil,alev,blev,val_missing)

! Reads a BOLAM-type MHF

implicit none

! Input

integer :: ifl
character (len=80) :: filename_atm, filename_soil
integer :: data_mode ! 1 - united file, 2 - separated files
logical :: ini_flag

! Output

integer :: level_type, nlon_inp0, nlon_inp, nlat_inp, nlev_atm_inp, nlev_atm_inp_max, nlev_soil_inp, nlev_soil_inp_max, &
           idate(5), iperiod(3), npar3d, npar3d_soil, npar2d
real :: x0a, y0a, x1a, y1a, dxa, dya, lev_list(nlev_atm_inp_max), lev_list_soil(nlev_soil_inp_max), &
        alev(nlev_atm_inp_max+1), blev(nlev_atm_inp_max+1), val_missing
real, dimension(nlon_inp,nlat_inp,nlev_atm_inp_max,npar3d) :: field3d
real, dimension(nlon_inp,nlat_inp,npar2d) :: field2d
real, dimension(nlon_inp,nlat_inp,nlev_soil_inp_max,npar3d_soil) :: field3d_soil

! For reading

integer :: unit_atm = 21, unit_soil=22, unit_const=23
integer, dimension(50) :: nfdr, nfdr_soil
real, dimension(200)   :: pdr, pdr_soil
character (len=80) :: filename_const = "model_param_constant_inp.bin"
integer, parameter :: nlevsnow=11

! For working

integer :: ierr_open_atm=0, ierr_open_soil=0, ierr_open_const=0, ierr_read=0, &
 ipar2d, ipar3d, ipar3d_soil, jlon, jlat, jklev, jklev_soil, jklev_snow, ird, &
 ierr, nlon_const, nlat_const, nlevg_const, nst_inp, nvt_inp
real :: dlon_const, dlat_const, x0_const, y0_const, alon0_const, alat0_const
real, parameter :: g0=9.807
real, dimension(nlon_inp,nlat_inp) :: read2d_work

 if (data_mode==2) then
   open (unit_atm, file=trim(filename_atm), status='old', form='unformatted', iostat=ierr_open_atm)
   open (unit_soil, file=trim(filename_soil), status='old', form='unformatted', iostat=ierr_open_soil)
 else
   if (ifl == 0) then
     open (unit_atm, file=trim(filename_atm), status='old', form='unformatted', iostat=ierr_open_atm)
     open (unit_soil, file=trim(filename_soil), status='old', form='unformatted', iostat=ierr_open_soil)
   endif
 endif

 if (ierr_open_atm /= 0) then
   print *,"Error in opening of input file ", trim(filename_atm)
   stop
 endif

 if (ierr_open_soil /= 0) then
   print *,"Error in opening of input file ", trim(filename_soil)
   stop
 endif

!  1.  Read descriptor records

 read(unit_soil, iostat=ierr_read) nfdr_soil
 read(unit_atm, iostat=ierr_read) nfdr
 if (ierr_read /= 0) then
   print *,"End of file encountered in reading input file ", trim(filename_atm)
   stop
 endif
 read(unit_atm) pdr
 read(unit_soil) pdr_soil

 if (ini_flag) then

   if (data_mode == 1) then
     rewind (unit_atm)
     rewind (unit_soil)
   else
     close (unit_atm)
     close (unit_soil)
   endif

   if (any(nfdr(:) /= nfdr_soil(:))) then
     print *,"Not concidence of nfdr parameters in ",trim(filename_atm)," and ",trim(filename_soil)," input files, stop"
     stop
   endif

   if (any(pdr(:) /= pdr_soil(:))) then
     print *,"Not concidence of pdr parameters in ",trim(filename_atm)," and ",trim(filename_soil)," input files, stop"
     stop
   endif

   level_type    = 2       ! hybrid or sigma levels
   nlon_inp0     = nfdr(2)
   nlon_inp      = nlon_inp0
   nlat_inp      = nfdr(3)
   nlev_atm_inp  = nfdr(4)
   nlev_soil_inp = nfdr(15)
   x0a           = pdr(39)
   y0a           = pdr(38)
   dxa           = pdr(2)
   dya           = pdr(1)
   x1a           = pdr(5)
   y1a           = pdr(4)
   alev(1)       = pdr(36)
   alev(2)       = pdr(37)
   blev(1:nlev_atm_inp) = pdr(39+1:39+nlev_atm_inp)

   do jklev = 1,nlev_atm_inp
     lev_list(jklev)=float(jklev)
   enddo

   lev_list_soil(1:nlev_soil_inp) = pdr(6:6+nlev_soil_inp-1)

   if (any(nfdr(5:9) /= nfdr_soil(5:9))) then
     print *,"Not coincidence of initial data read from ",trim(filename_atm)," and ",trim(filename_soil),":"
     print *,trim(filename_atm)," ",nfdr(5:9)
     print *,trim(filename_soil)," ",nfdr_soil(5:9)
     print *," stop"
     stop
   endif

   idate(1:5) = nfdr(5:9)

   print*
   print *,'Date (YYYY MM DD) and time (HH MM, analysis time) of the input fields:'
   print *,idate(1:3),'            ',idate(4:5)

   return

 endif

 if (ifl == 1) then

! Reading of constant physiographical parameters of input model

   open (unit_const, file=trim(filename_const), status='old', form='unformatted', iostat=ierr_open_const)

   if (ierr_open_const /= 0) then
     print *,"Error in opening of input file ", trim(filename_const),", stop"
     stop
   endif

   read (unit_const) nlon_const, nlat_const, nlevg_const, dlon_const, dlat_const, &
 x0_const, y0_const, alon0_const, alat0_const, nst_inp, nvt_inp

   ierr=0
   if (nlon_inp /= nlon_const) ierr=ierr+1
   if (nlat_inp /= nlat_const) ierr=ierr+1
   if (nlev_soil_inp /= nlevg_const) ierr=ierr+1
   if (dxa /= dlon_const) ierr=ierr+1
   if (dya /= dlat_const) ierr=ierr+1
   if (x0a /= x0_const) ierr=ierr+1
   if (y0a /= y0_const) ierr=ierr+1
   if (x1a /= alon0_const) ierr=ierr+1
   if (y1a /= alat0_const-dlat_const*0.5) ierr=ierr+1

   if (ierr /= 0) then
     write (*,*)
     write (*,*) "Error in header parameters in input file, ",trim(filename_const),", not coincident with parameters read from ",&
 trim(filename_atm)
     write (*,*) trim(filename_atm)," nlon, nlat, nlevg, dlon, dlat, x0, y0, alon0, alat0 :", &
 nlon_inp, nlat_inp, nlev_soil_inp, dxa, dya, x0a, y0a, x1a, y1a
     write (*,*) trim(filename_const)," nlon, nlat, nlevg, dlon, dlat, x0, y0, alon0, alat0 :", &
 nlon_const, nlat_const, nlevg_const, dlon_const, dlat_const, x0_const, y0_const, alon0_const, alat0_const-dlat_const*0.5
     write (*,*) "  stop"
     stop
   endif

! Sea-land fraction

   ipar2d = 3        ! fmask
   call rrec2 (unit_const, nlon_inp, nlat_inp, field2d(1:nlon_inp,1:nlat_inp,ipar2d))

! Topography hight

   ipar2d = 2        ! Topography
   call rrec2 (unit_const, nlon_inp, nlat_inp, field2d(1:nlon_inp,1:nlat_inp,ipar2d))
   field2d(:,:,ipar2d) = field2d(:,:,ipar2d)*g0

   do ird = 1,1 ! Topography variation
     call rrec2 (unit_const, nlon_inp, nlat_inp, read2d_work(1:nlon_inp,1:nlat_inp))
   enddo

   do ird = 1,nst_inp+1 ! Soil types
     call rrec2 (unit_const, nlon_inp, nlat_inp, read2d_work(1:nlon_inp,1:nlat_inp))
   enddo

   do ird = 1,nvt_inp+1 ! Vegetation types
     call rrec2 (unit_const, nlon_inp, nlat_inp, read2d_work(1:nlon_inp,1:nlat_inp))
   enddo

   ipar3d_soil = 7     ! qgmax
   do jklev_soil = 1,nlev_soil_inp
     call rrec2 (unit_const, nlon_inp, nlat_inp, field3d_soil(1:nlon_inp,1:nlat_inp,jklev,ipar3d_soil))
   enddo

   ipar3d_soil = 8     ! qgmin
   do jklev_soil = 1,nlev_soil_inp
     call rrec2 (unit_const, nlon_inp, nlat_inp, field3d_soil(1:nlon_inp,1:nlat_inp,jklev,ipar3d_soil))
   enddo

   close (unit_const)

 endif ! ifl == 1

 if (any(nfdr(10:12) /= nfdr_soil(10:12))) then
   print *,"Not coincidence of validity time read from ",trim(filename_atm)," and ",trim(filename_soil),":"
   print *,trim(filename_atm)," ",nfdr(10:12)
   print *,trim(filename_soil)," ",nfdr_soil(10:12)
   print *," stop"
   stop
 endif

 iperiod(1:3) = nfdr(10:12)

 print *,'Validity time if forecast (00 if analysis) (days, hours, minutes):'
 print *,iperiod(1:3)

! Atmospheric variables

 ipar2d = 4        ! Surface pressure
 call rrec2 (unit_atm, nlon_inp, nlat_inp, field2d(1:nlon_inp,1:nlat_inp,ipar2d))

 ipar3d = 3        ! u
 do jklev = 1,nlev_atm_inp
   call rrec2 (unit_atm, nlon_inp, nlat_inp, field3d(1:nlon_inp,1:nlat_inp,jklev,ipar3d))
 enddo

 ipar3d = 4        ! v
 do jklev = 1,nlev_atm_inp
   call rrec2 (unit_atm, nlon_inp, nlat_inp, field3d(1:nlon_inp,1:nlat_inp,jklev,ipar3d))
 enddo

 ipar3d = 2        ! t
 do jklev = 1,nlev_atm_inp
   call rrec2 (unit_atm, nlon_inp, nlat_inp, field3d(1:nlon_inp,1:nlat_inp,jklev,ipar3d))
 enddo

 ipar3d = 5        ! q
 do jklev = 1,nlev_atm_inp
   call rrec2 (unit_atm, nlon_inp, nlat_inp, field3d(1:nlon_inp,1:nlat_inp,jklev,ipar3d))
 enddo

 ipar3d = 7        ! qc
 do jklev = 1,nlev_atm_inp
   call rrec2 (unit_atm, nlon_inp, nlat_inp, field3d(1:nlon_inp,1:nlat_inp,jklev,ipar3d))
 enddo

! Physiographical parameters changing in time

 do ird = 1,4 ! Vegetation LAI, Vegetation fraction, RGM, RGQ
   call rrec2 (unit_soil, nlon_inp, nlat_inp, read2d_work(1:nlon_inp,1:nlat_inp))
 enddo

 ipar2d = 22       ! Sea ice fraction
 call rrec2 (unit_soil, nlon_inp, nlat_inp, field2d(1:nlon_inp,1:nlat_inp,ipar2d))

 ipar2d = 26       ! Sea ice thickness (m)
 call rrec2 (unit_soil, nlon_inp, nlat_inp, field2d(1:nlon_inp,1:nlat_inp,ipar2d))

! Surface, soil, sea variables

 do ird = 1,7 ! Albedo, 2 emisivities, Tot. cloudness, Tot. prec., Conv. prec., Snow prec.
   call rrec2 (unit_soil, nlon_inp, nlat_inp, read2d_work(1:nlon_inp,1:nlat_inp))
 enddo

 ipar2d = 5        ! tsurf
 call rrec2 (unit_soil, nlon_inp, nlat_inp, field2d(1:nlon_inp,1:nlat_inp,ipar2d))

 do ird = 1,1 ! tg surface
   call rrec2 (unit_soil, nlon_inp, nlat_inp, read2d_work(1:nlon_inp,1:nlat_inp))
 enddo

 ipar3d_soil = 1   ! tg
 do jklev = 1,nlev_soil_inp
   call rrec2 (unit_soil, nlon_inp, nlat_inp, field3d_soil(1:nlon_inp,1:nlat_inp,jklev,ipar3d_soil))
 enddo

 ipar2d = 25       ! qvsurf
 call rrec2 (unit_soil, nlon_inp, nlat_inp, field2d(1:nlon_inp,1:nlat_inp,ipar2d))

 do ird = 1,1 ! qg surface
   call rrec2 (unit_soil, nlon_inp, nlat_inp, read2d_work(1:nlon_inp,1:nlat_inp))
 enddo

 ipar3d_soil = 2   ! qg
 do jklev = 1,nlev_soil_inp
   call rrec2 (unit_soil, nlon_inp, nlat_inp, field3d_soil(1:nlon_inp,1:nlat_inp,jklev,ipar3d_soil))
 enddo

 do ird = 1,nlev_soil_inp+1 ! Soil ice fraction at surface and at soil levels
   call rrec2 (unit_soil, nlon_inp, nlat_inp, read2d_work(1:nlon_inp,1:nlat_inp))
 enddo

 ipar2d = 14       ! snow
 call rrec2 (unit_soil, nlon_inp, nlat_inp, field2d(1:nlon_inp,1:nlat_inp,ipar2d))

 do ird = 1,6 ! Snow mass, temperature, ice fraction, age, melting age, density at snow levels
   do jklev_snow = 1,nlevsnow
     call rrec2 (unit_soil, nlon_inp, nlat_inp, read2d_work(1:nlon_inp,1:nlat_inp))
   enddo
 enddo

 do ird = 1, 12   ! Snow albedo, CSWFl, CLWFl, CHFlux, CQFlux, T2Min, T2Max, Runoff, Runoff total, 3 empty
   call rrec2 (unit_soil, nlon_inp, nlat_inp, read2d_work(1:nlon_inp,1:nlat_inp))
 enddo

 if (data_mode == 2) then
   close (unit_atm)
   close (unit_soil)
 endif

return
end

!=======================================================================
subroutine read_bolam_mhf_data_old(ifl,filename_atm,filename_soil,data_mode,ini_flag,                  &
    nlon_inp,nlat_inp,nlev_atm_inp,nlev_atm_inp_max,nlev_soil_inp,nlev_soil_inp_max,     &
    x0a,y0a,x1a,y1a,dxa,dya,idate,iperiod,lev_list,lev_list_soil,level_type,     &
    npar3d,npar2d,npar3d_soil,field3d,field2d,field3d_soil,alev,blev,val_missing)

! Reads a BOLAM-type MHF

implicit none

! Input

integer :: ifl
character (len=80) :: filename_atm, filename_soil
integer :: data_mode ! 1 - united file, 2 - separated files
logical :: ini_flag

! Output

integer :: level_type, nlon_inp0, nlon_inp, nlat_inp, nlev_atm_inp, nlev_atm_inp_max, nlev_soil_inp, nlev_soil_inp_max, &
           idate(5), iperiod(3), npar3d, npar3d_soil, npar2d
real :: x0a, y0a, x1a, y1a, dxa, dya, lev_list(nlev_atm_inp_max), lev_list_soil(nlev_soil_inp_max), &
        alev(nlev_atm_inp_max+1), blev(nlev_atm_inp_max+1), val_missing
real, dimension(nlon_inp,nlat_inp,nlev_atm_inp_max,npar3d) :: field3d
real, dimension(nlon_inp,nlat_inp,npar2d) :: field2d
real, dimension(nlon_inp,nlat_inp,nlev_soil_inp_max,npar3d_soil) :: field3d_soil

! For reading

integer :: unit_atm = 21, unit_soil=22, unit_const=23
integer, dimension(50) :: nfdr, nfdr_soil
real, dimension(200)   :: pdr, pdr_soil
character (len=80) :: filename_const = "model_param_constant_inp.bin"
integer, parameter :: nlevsnow=11

! For working

integer :: ierr_open_atm=0, ierr_open_soil=0, ierr_open_const=0, ierr_read=0, &
 ipar2d, ipar3d, ipar3d_soil, jlon, jlat, jklev, jklev_soil, jklev_snow, ird, &
 ierr, nlon_const, nlat_const, nlevg_const, nst_inp, nvt_inp
real :: dlon_const, dlat_const, x0_const, y0_const, alon0_const, alat0_const
real, parameter :: g0=9.807
real, dimension(nlon_inp,nlat_inp) :: read2d_work

 if (data_mode==2) then
   open (unit_atm, file=trim(filename_atm), status='old', form='unformatted', iostat=ierr_open_atm)
   open (unit_soil, file=trim(filename_soil), status='old', form='unformatted', iostat=ierr_open_soil)
 else
   if (ifl == 0) then
     open (unit_atm, file=trim(filename_atm), status='old', form='unformatted', iostat=ierr_open_atm)
     open (unit_soil, file=trim(filename_soil), status='old', form='unformatted', iostat=ierr_open_soil)
   endif
 endif

 if (ierr_open_atm /= 0) then
   print *,"Error in opening of input file ", trim(filename_atm)
   stop
 endif

 if (ierr_open_soil /= 0) then
   print *,"Error in opening of input file ", trim(filename_soil)
   stop
 endif

!  1.  Read descriptor records

 read(unit_soil, iostat=ierr_read) nfdr_soil
 read(unit_atm, iostat=ierr_read) nfdr
 if (ierr_read /= 0) then
   print *,"End of file encountered in reading input file ", trim(filename_atm)
   stop
 endif
 read(unit_atm) pdr
 read(unit_soil) pdr_soil

 if (ini_flag) then

   if (data_mode == 1) then
     rewind (unit_atm)
     rewind (unit_soil)
   else
     close (unit_atm)
     close (unit_soil)
   endif

   if (any(nfdr(:) /= nfdr_soil(:))) then
     print *,"Not concidence of nfdr parameters in ",trim(filename_atm)," and ",trim(filename_soil)," input files, stop"
     stop
   endif

   if (any(pdr(:) /= pdr_soil(:))) then
     print *,"Not concidence of pdr parameters in ",trim(filename_atm)," and ",trim(filename_soil)," input files, stop"
     stop
   endif

   level_type    = 2       ! hybrid or sigma levels
   nlon_inp0     = nfdr(2)
   nlon_inp      = nlon_inp0
   nlat_inp      = nfdr(3)
   nlev_atm_inp  = nfdr(4)
   nlev_soil_inp = nfdr(15)
   x0a           = pdr(39)
   y0a           = pdr(38)
   dxa           = pdr(2)
   dya           = pdr(1)
   x1a           = pdr(5)
   y1a           = pdr(4)
   alev(1)       = pdr(36)
   alev(2)       = pdr(37)
   blev(1:nlev_atm_inp) = pdr(39+1:39+nlev_atm_inp)

   do jklev = 1,nlev_atm_inp
     lev_list(jklev)=float(jklev)
   enddo

   lev_list_soil(1:nlev_soil_inp) = pdr(6:6+nlev_soil_inp-1)

   if (any(nfdr(5:9) /= nfdr_soil(5:9))) then
     print *,"Not coincidence of initial data read from ",trim(filename_atm)," and ",trim(filename_soil),":"
     print *,trim(filename_atm)," ",nfdr(5:9)
     print *,trim(filename_soil)," ",nfdr_soil(5:9)
     print *," stop"
     stop
   endif

   idate(1:5) = nfdr(5:9)

   print*
   print *,'Date (YYYY MM DD) and time (HH MM, analysis time) of the input fields:'
   print *,idate(1:3),'            ',idate(4:5)

   return

 endif

 if (ifl == 1) then

! Reading of constant physiographical parameters of input model

   open (unit_const, file=trim(filename_const), status='old', form='unformatted', iostat=ierr_open_const)

   if (ierr_open_const /= 0) then
     print *,"Error in opening of input file ", trim(filename_const),", stop"
     stop
   endif

   read (unit_const) nlon_const, nlat_const, nlevg_const, dlon_const, dlat_const, &
 x0_const, y0_const, alon0_const, alat0_const, nst_inp, nvt_inp

   ierr=0
   if (nlon_inp /= nlon_const) ierr=ierr+1
   if (nlat_inp /= nlat_const) ierr=ierr+1
   if (nlev_soil_inp /= nlevg_const) ierr=ierr+1
   if (dxa /= dlon_const) ierr=ierr+1
   if (dya /= dlat_const) ierr=ierr+1
   if (x0a /= x0_const) ierr=ierr+1
   if (y0a /= y0_const) ierr=ierr+1
   if (x1a /= alon0_const) ierr=ierr+1
   if (y1a /= alat0_const-dlat_const*0.5) ierr=ierr+1

   if (ierr /= 0) then
     write (*,*)
     write (*,*) "Error in header parameters in input file, ",trim(filename_const),", not coincident with parameters read from ",&
 trim(filename_atm)
     write (*,*) trim(filename_atm)," nlon, nlat, nlevg, dlon, dlat, x0, y0, alon0, alat0 :", &
 nlon_inp, nlat_inp, nlev_soil_inp, dxa, dya, x0a, y0a, x1a, y1a
     write (*,*) trim(filename_const)," nlon, nlat, nlevg, dlon, dlat, x0, y0, alon0, alat0 :", &
 nlon_const, nlat_const, nlevg_const, dlon_const, dlat_const, x0_const, y0_const, alon0_const, alat0_const-dlat_const*0.5
     write (*,*) "  stop"
     stop
   endif

! Sea-land fraction

   ipar2d = 3        ! fmask
   do jlat = 1, nlat_inp
     read (unit_const) (field2d(jlon,jlat,ipar2d), jlon = 1, nlon_inp)
   enddo

! Topography hight

   ipar2d = 2        ! Topography
   do jlat = 1, nlat_inp
     read (unit_const) (field2d(jlon,jlat,ipar2d), jlon = 1, nlon_inp)
   enddo
   field2d(:,:,ipar2d) = field2d(:,:,ipar2d)*g0

   do ird = 1,1 ! Topography variation
     do jlat = 1, nlat_inp
       read (unit_const) (read2d_work(jlon,jlat), jlon = 1, nlon_inp)
     enddo
   enddo

   do ird = 1,nst_inp+1 ! Soil types
     do jlat = 1, nlat_inp
       read (unit_const) (read2d_work(jlon,jlat), jlon = 1, nlon_inp)
     enddo
   enddo

   do ird = 1,nvt_inp+1 ! Vegetation types
     do jlat = 1, nlat_inp
       read (unit_const) (read2d_work(jlon,jlat), jlon = 1, nlon_inp)
     enddo
   enddo

   ipar3d_soil = 7     ! qgmax
   do jklev_soil = 1,nlev_soil_inp
     do jlat = 1, nlat_inp
       read (unit_const) (field3d_soil(jlon,jlat,jklev_soil,ipar3d_soil), jlon = 1, nlon_inp)
     enddo
   enddo

   ipar3d_soil = 8     ! qgmin
   do jklev_soil = 1,nlev_soil_inp
     do jlat = 1, nlat_inp
       read (unit_const) (field3d_soil(jlon,jlat,jklev_soil,ipar3d_soil), jlon = 1, nlon_inp)
     enddo
   enddo

   close (unit_const)

 endif ! ifl == 1

 if (any(nfdr(10:12) /= nfdr_soil(10:12))) then
   print *,"Not coincidence of validity time read from ",trim(filename_atm)," and ",trim(filename_soil),":"
   print *,trim(filename_atm)," ",nfdr(10:12)
   print *,trim(filename_soil)," ",nfdr_soil(10:12)
   print *," stop"
   stop
 endif

 iperiod(1:3) = nfdr(10:12)

 print *,'Validity time if forecast (00 if analysis) (days, hours, minutes):'
 print *,iperiod(1:3)

! Atmospheric variables

 ipar2d = 4        ! Surface pressure
 do jlat = 1, nlat_inp
   read (unit_atm) (field2d(jlon,jlat,ipar2d), jlon = 1, nlon_inp)
 enddo

 ipar3d = 3        ! u
 do jklev = 1,nlev_atm_inp
   do jlat = 1, nlat_inp
     read (unit_atm) (field3d(jlon,jlat,jklev,ipar3d), jlon = 1, nlon_inp)
   enddo
 enddo

 ipar3d = 4        ! v
 do jklev = 1,nlev_atm_inp
   do jlat = 1, nlat_inp
     read (unit_atm) (field3d(jlon,jlat,jklev,ipar3d), jlon = 1, nlon_inp)
   enddo
 enddo

 ipar3d = 2        ! t
 do jklev = 1,nlev_atm_inp
   do jlat = 1, nlat_inp
     read (unit_atm) (field3d(jlon,jlat,jklev,ipar3d), jlon = 1, nlon_inp)
   enddo
 enddo

 ipar3d = 5        ! q
 do jklev = 1,nlev_atm_inp
   do jlat = 1, nlat_inp
     read (unit_atm) (field3d(jlon,jlat,jklev,ipar3d), jlon = 1, nlon_inp)
   enddo
 enddo

 ipar3d = 7        ! qc
 do jklev = 1,nlev_atm_inp
   do jlat = 1, nlat_inp
     read (unit_atm) (field3d(jlon,jlat,jklev,ipar3d), jlon = 1, nlon_inp)
   enddo
 enddo

! Physiographical parameters changing in time

 do ird = 1,4 ! Vegetation LAI, Vegetation fraction, RGM, RGQ
   do jlat = 1, nlat_inp
     read (unit_soil) (read2d_work(jlon,jlat), jlon = 1, nlon_inp)
   enddo
 enddo

 ipar2d = 22       ! Sea ice fraction
 do jlat = 1, nlat_inp
   read (unit_soil) (field2d(jlon,jlat,ipar2d), jlon = 1, nlon_inp)
 enddo

 ipar2d = 26       ! Sea ice thickness (m)
 do jlat = 1, nlat_inp
   read (unit_soil) (field2d(jlon,jlat,ipar2d), jlon = 1, nlon_inp)
 enddo

! Surface, soil, sea variables

 do ird = 1,7 ! Albedo, 2 emisivities, Tot. cloudness, Tot. prec., Conv. prec., Snow prec.
   do jlat = 1, nlat_inp
     read (unit_soil) (read2d_work(jlon,jlat), jlon = 1, nlon_inp)
   enddo
 enddo

 ipar2d = 5        ! tsurf
 do jlat = 1, nlat_inp
   read (unit_soil) (field2d(jlon,jlat,ipar2d), jlon = 1, nlon_inp)
 enddo

 do ird = 1,1 ! tg surface
   do jlat = 1, nlat_inp
     read (unit_soil) (read2d_work(jlon,jlat), jlon = 1, nlon_inp)
   enddo
 enddo

 ipar3d_soil = 1   ! tg
 do jklev = 1,nlev_soil_inp
   do jlat = 1, nlat_inp
     read (unit_soil) (field3d_soil(jlon,jlat,jklev,ipar3d_soil), jlon = 1, nlon_inp)
   enddo
 enddo

 ipar2d = 25       ! qvsurf
 do jlat = 1, nlat_inp
   read (unit_soil) (field2d(jlon,jlat,ipar2d), jlon = 1, nlon_inp)
 enddo

 do ird = 1,1 ! qg surface
   do jlat = 1, nlat_inp
     read (unit_soil) (read2d_work(jlon,jlat), jlon = 1, nlon_inp)
   enddo
 enddo

 ipar3d_soil = 2   ! qg
 do jklev = 1,nlev_soil_inp
   do jlat = 1, nlat_inp
     read (unit_soil) (field3d_soil(jlon,jlat,jklev,ipar3d_soil), jlon = 1, nlon_inp)
   enddo
 enddo

 do ird = 1,nlev_soil_inp+1 ! Soil ice fraction at surface and at soil levels
   do jlat = 1, nlat_inp
     read (unit_soil) (read2d_work(jlon,jlat), jlon = 1, nlon_inp)
   enddo
 enddo

 ipar2d = 14       ! snow
 do jlat = 1, nlat_inp
   read (unit_soil) (field2d(jlon,jlat,ipar2d), jlon = 1, nlon_inp)
 enddo

 do ird = 1,6 ! Snow mass, temperature, ice fraction, age, melting age, density at snow levels
   do jklev_snow = 1,nlevsnow
     do jlat = 1, nlat_inp
       read (unit_soil) (read2d_work(jlon,jlat), jlon = 1, nlon_inp)
     enddo
   enddo
 enddo

 do ird = 1, 12   ! Snow albedo, CSWFl, CLWFl, CHFlux, CQFlux, T2Min, T2Max, Runoff, Runoff total, 3 empty
   do jlat = 1, nlat_inp
     read (unit_soil) (read2d_work(jlon,jlat), jlon = 1, nlon_inp)
   enddo
 enddo

 if (data_mode == 2) then
   close (unit_atm)
   close (unit_soil)
 endif

return
end

!=======================================================================
