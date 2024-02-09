!
!     program postmoloch

! Last update 04/11/2022

! Nov. 2022: prodotti del manto nevoso: profondita', temperatura, densita'
! ed eta della neve alla superficie, in fondo e alla meta' del manto.
! Nov. 2022: altidutine (sopra il livello di mare) della nevicata,
! 2 campi, 1 campo calcola sulla base della temperatura del bulbo bagnato (+1,5 dec.C),
! e 1 campo corretto con la frazione degli idrometeori ghiacciati. 
! Nov. 2022: PS (Surface Pressure) in output in grib2 (per MeteoAosta)

! Dic. 2020: correzzione dell'errore della griglia diradata in input (mhfr);
! correzzione dell'iterpolazione della temperatura sotto l'orografia (sottocolle).

! Nov. 2020: scrittura dei parametri della coordinata verticale in grib2;
! PV sui livelli isobarici;
! Altezza dello zero termico.

! Ott. 2019: nuova espress. coord. vert. (introd. param. a0).
! In output due wind gust:
! 1. gust istantanea basata sul calcolo dell'altezza del PBL,
!    nella matrice gust, campo 57 in wrpost;
! 2. gust max. nell'interv. tra due istanti di postprocessing (posson comprendere
!    piu' intervalli tra mhf se njump>1), calcolata in moloch in funz. di TKE,
!    nella matrice ws10max, campo 49 in wrpost.

! Ago. 2018: nuovo formato mhf: mfs_atm, mhf_soil e il file con i campi statici (model_param_constant.bin).

! Ott. 2017: inserito calcolo di CIN, altezza top del PBL e raffiche del vento a 10 m (gust)
! Inserito fix calcolo quantita' al suolo su montagne rough
! Sett. 2017: rivista cclouds (riduz. nubi).
! Ago. 2017: roughness per calcolo quantita' al suolo in vtsurf posta a 0.05
! Mag. 2017:
! Definita nuova routine cclouds che usa algoritmo (Maurizio) per calcolo nubi alte,
! medie, basse e tot.
! Definito il numero di Richardson (usato in cloud_post per il calcolo delle nubi)
! al primo semi-livello sopra il suolo come zri (definito in vtsurf).
! Introdotta matrice richs (cloud_post si deve chiamare dopo vtsurf)
! Nuova static stability per ridurre cloud fraction (si usa Richardson calc. con Durran&Klemp
! per la parte umida); nuovo algoritmo (Maurizio) per calcolo nubi alte, medie, basse e tot.
! Mar. 2017: velocizzato evitando calcoli inutili nel caso njump>1.

! Versione che puo' plottare livelli a Z costante - attenz.: qui commentata WRITE di IZFLAG
! Se si lasciano coord. standard per le cross-section, va commentata la parte di codice
! del DO sotto il commento "geographycal points with prescribed coordinates"

!-----------------------------------------------------------------------
 module mod_postmoloch

! npar: no. of variables on pressure level in output
! nlevpo: no. of constant pressure levels in output
! ncrossx, ncrossy: no. of cross-sections in long. and lat.
! nzlev: no. of constant height levels in output

   integer, parameter :: npar=8, nlevpo=20, ncrossx=3, ncrossy=3, nzlev=20, npoint_cross_max=1000
   integer :: njump, nlivz, cross_number, cross_interp_gauss
   real :: zcbot, zctop        ! top and bottom zeta for cross-sections
   logical :: ncrossec, output_format_ppf, output_format_grib2
   integer, dimension(ncrossy)  :: loncr
   integer, dimension(ncrossx)  :: latcr
   integer, dimension(80)       :: isflag
   integer, dimension(npar,nlevpo) :: ipflag
   integer, dimension(npar*nlevpo) :: ipflag0
   integer, dimension(nzlev)    :: izflag
   real, dimension(10) :: lonini_cross, latini_cross, lonfin_cross, latfin_cross
   integer, dimension(10) :: npoint_cross
   real :: valmiss=-9999.
   integer :: nunic1, nunic2, ierr_read1=0, ierr_read2=0

   namelist/postp/ njump,zcbot,zctop,nlivz,ncrossec,output_format_ppf,output_format_grib2,              &
                   loncr,latcr,                                                                         &
                   lonini_cross,latini_cross,lonfin_cross,latfin_cross,npoint_cross,cross_interp_gauss, &
                   isflag,ipflag0,izflag

      CONTAINS

      function gzita (zita,h,a0)    ! Decay function
      gzita = 1. -a0*(zita/h)-(3.-2.*a0)*(zita/h)**2 +(2.-a0)*(zita/h)**3
      end function gzita

      function gzitap (zita,h,a0)   ! Derivative of decay function
      gzitap = (-a0-(6.-4.*a0)*(zita/h) +(6.-3.*a0)*(zita/h)**2)/h
      end function gzitap

      function bzita (zita,h,b0)    ! Stretching function
      bzita = b0 + (1.-b0)*(zita/h)
      end function bzita

      function bzitap (zita,h,b0)   ! Derivative of stretching function
      bzitap = (1.-b0)/h
      end function bzitap

 end module mod_postmoloch
!--------------------------------------------------------------------------------------------------
 module data_observ_aux
   real, dimension(:,:), allocatable :: t05, t005, q05rel, tg005, tg010, tg020, tg050, tg100
 end module data_observ_aux
!--------------------------------------------------------------------------------------------------

! ---- For RTTOV simulation ---

#ifdef rttov
 module mod_rttov
   integer, parameter :: nchan_max=20
   integer :: nchan, rttov_sat_series, rttov_sat_id, rttov_sat_sensor, flag_visible
   integer, dimension (nchan_max) :: channel_elabor_list
   character(len=256) :: path_rttov_emis_atlas, path_rttov_brdf_atlas
   namelist /rttov_par/ rttov_sat_series, rttov_sat_id, rttov_sat_sensor, channel_elabor_list, flag_visible, &
 path_rttov_emis_atlas, path_rttov_brdf_atlas
   integer, dimension(5) :: date_time
   real, dimension(:,:,:), allocatable :: p_invert, t_invert, q_invert, qcw_invert, qci_invert, clfrac_invert
   real, dimension(:,:), allocatable :: fseaice
   integer, dimension(:,:), allocatable :: flag_cloud_conv
   real :: w_crit_conv=1.5, fracw, z1, z2, z3, z4, z5, z6
   integer, dimension(:), allocatable :: &
 sensor_chan_id  ! List of "true" channels index, not RTTOV channels list
   real*8, dimension(:), allocatable :: &
 sensor_chan_cw  ! Central Wave Number (m^-1) of elaborated channels
   real, dimension(:,:,:), allocatable :: radiance, radiance_bt, radiance_refl, &
 radiance_clear, radiance_clear_bt, radiance_clear_refl, emis_rttov
   character(len=50) :: name_rttov_out_grib2
   real :: zoom_xini=0.00, zoom_xfin=1.00, zoom_yini=0.00, zoom_yfin=1.00, alon1, alat1, val_missing=-9999.
   integer :: npoint_jump, iini, ifin, jini, jfin
 end module mod_rttov
#endif
!--------------------------------------------------------------------------------------------------
      program postmoloch

! Postprocessing of MOLOCH with NLEVG soil levels

 use mod_postmoloch
 use data_observ_aux

! ---- For RTTOV simulation ---

#ifdef rttov
 use mod_rttov
#endif

! -----------------------------

      real, parameter                     :: tzer=273.15, pzer=1.e5, ezer=611., cw=4186.8, ci=2093.4
      real, parameter                     :: yliv=2834170.5, ylwv=2500610.
      integer, dimension(50)              :: nfdr
      real, dimension(100)                :: pdr
      integer, parameter                  :: nlev_snow=11, nlev_snow_prod=3
      real, dimension(:), allocatable     :: zlcr
      real, dimension(nlevpo)             :: plevo, zalp
      real, dimension(nzlev)              :: zlev
      real, dimension(1)                  :: zsq

      integer, dimension(:,:),allocatable :: izout, izoutx, izouty
      real, dimension(:),     allocatable :: fmyu, fz, fzh, zaux, zauxp, zlnp, zlnph

      real, dimension(:,:),   allocatable :: phig, tskin, qskin, cloudt, prectot, precconv, precsolid, runoff, runoff_tot, &
                                             snow, fsnow, albedo, &
                                             fice, iceth, rgm, rgq, fmask, emis1, emis2, hx, hy,                       &
                                             alont, alatt, alonu, alatu, alonv, alatv, alonz, alatz, fcorio,           &
                                             zorogr, ps, zpz0, zus, zvs, zts, zqs, zthes, ztskin,                      &
                                             ztgr1, ztgr2, ztgr3, ztgr4, zqgr1, zqgr2, zqgr3, zqgr4,                   &
                                             zwork, zwor2, zwor3, zwor4, zwor5, hflux, qflux, t2, td2, q2, q2rel,      &
                                             u10, v10, gust, zlev0,                                                    &
                                             cloudl, cloudm, cloudh, cape, capek, cswfl, clwfl, chflux, cqflux, t2min, &
                                             t2max, ws10max, zlift, ztvirtup, orog_cross, ztarr, zqarr,                &
                                             lapse_rate_1, lapse_rate_2, iwv, richs, cin, cink, pbl,                   &
                                             tgsurf, qgsurf, fice_soil_surf, snow_albedo,                              &
                                             cwvflux, csoilh_bott_flux, csoilw_bott_flux, snowfall_height, snowfall_height0
      integer, dimension(:,:), allocatable :: snowfall_level

      real, dimension(:,:,:), allocatable :: tg, qg, qgmax, qgmin, p, u, v, w, t, q, qcw, qci, s, zthee, fmz, fmzh, potvor,   &
                                             zrelh, tvirt, zph, zt, zq, zu, zv, zw, zqcw, zqci, zpv, teta,        &
                                             zclwi, zrh, zthetae, tvlift, ztvlift, zparr, zeta, tg_post, qg_post, &
                                             ztcrx, ztecrx, zucrx, zvcrx, zwcrx, zrhcrx, zcwcrx, zcicrx,          &
                                             ztcry, ztecry, zucry, zvcry, zwcry, zrhcry, zcwcry, zcicry,          &
                                             t_cross, theta_cross, thetae_cross, u_cross, v_cross, u_tang_cross,  &
                                             v_norm_cross, w_cross, rh_cross, qcw_cross, qci_cross,               &
                                             fice_soil, snow_lev, snow_t, snow_fice, snow_age, snow_melt_age, snow_dens, &
                                             snow_lev_depth, snow_lev_depth_prod, snow_t_prod, snow_age_prod, snow_dens_prod

      real, dimension(:), allocatable     :: zzzg, zztg
      real, dimension(nlev_snow)          :: snow_dens_h

      data zlev/ 0.5e3, 1.0e3, 1.5e3, 2.0e3, 2.5e3, 3.0e3, 3.5e3, 4.0e3, 4.5e3, 5.0e3, &
                 5.5e3, 6.0e3, 6.5e3, 7.0e3, 7.5e3, 8.0e3, 8.5e3, 9.0e3, 9.5e3,10.0e3/

   real, dimension(:,:), allocatable :: prectotr, precconvr, precsolidr, runoffr, runoff_tot_r, &
                                        cswflr, clwflr, chfluxr, cqfluxr, t2minr, t2maxr, ws10maxr, &
                                        lon_crossy, lat_crossy, lon_crossx, lat_crossx, orog_crossx, orog_crossy, &
                                        cwvfluxr, csoilh_bott_fluxr, csoilw_bott_fluxr

   integer, dimension(5) :: idate0, idatec
   integer, dimension(3) :: iperiod=0, iperiod_accum=0
   integer, dimension(12):: imon=(/31,28,31,30,31,30,31,31,30,31,30,31/)
   integer :: nday, ndm, iimon, ierr_open1, ierr_open2
   real :: zlat, day, coeday, lapse_rate, lapse_rate_clim, zfrac
   integer :: kbott, kmedium

!  Definition of postprocessing procedure parameters

   open (11,file='postmoloch.inp',status='old')
   read(11,postp)
!   write(*,postp)
   close (11)

  ipflag(1:npar,1:nlevpo)=reshape(ipflag0,(/npar,nlevpo/))

!  Definition of values of pressure levels in Pa

      plevo(1) = 1000.e2
      plevo(2) =  950.e2
      plevo(3) =  900.e2
      plevo(4) =  850.e2
      plevo(5) =  800.e2
      plevo(6) =  750.e2
      plevo(7) =  700.e2
      plevo(8) =  650.e2
      plevo(9) =  600.e2
      plevo(10)=  550.e2
      plevo(11)=  500.e2
      plevo(12)=  450.e2
      plevo(13)=  400.e2
      plevo(14)=  350.e2
      plevo(15)=  300.e2
      plevo(16)=  250.e2
      plevo(17)=  200.e2
      plevo(18)=  150.e2
      plevo(19)=  100.e2
      plevo(20)=   50.e2

      do jklev = 1, nlevpo
      zalp(jklev) = alog(plevo(jklev))
      enddo

! ---- For RTTOV simulation ---

#ifdef rttov
 open (11, file="rttov_interface.inp", status="old")
 read (11, rttov_par)
 close (11)
 do ix=1,nchan_max
   if (channel_elabor_list(ix) <= 0 ) exit
 enddo
 nchan=ix-1
#endif

! -----------------------------

      nunic1 = 21
      nunic2 = 22

      open (nunic1,file='moloch_atm.mhf',status='old',form='unformatted',iostat=ierr_open1)
      print *
      if (ierr_open1 /= 0) then
        print *,'Not found input file moloch_atm.mhf'
        stop
      else
        print *,'Input file moloch_atm.mhf opened on the unit ',nunic1
      endif
      open (nunic2,file='moloch_soil.mhf',status='old',form='unformatted',iostat=ierr_open2)
      print *
      if (ierr_open2 /= 0) then
        print *,'Not found input file moloch_soil.mhf'
        stop
      else
        print *,'Input file moloch_soil.mhf opened on the unit ',nunic2
      endif
      print *

      read (nunic1) nfdr
      read (nunic1) pdr

      rewind nunic1

      nlon   = nfdr(2)
      nlat   = nfdr(3)
      nlev   = nfdr(4)
      nlevg  = nfdr(15)
      mhfr   = nfdr(20)
      if (mhfr /= 2) mhfr = 1
      nlevp1 = nlev+1
      nlonm1 = nlon-1
      nlatm1 = nlat-1

      write(*, '(a,4i5)'), " nlon, nlat, nlev, nlevg:", nlon, nlat, nlev, nlevg
      write(*, '(a,i5)'), " mhfr=", mhfr

!------------------------------------------------------------
! Definition of longitudes and latitudes for cross sections
! Caution: redefines values defined in postmoloch.inp file
!------------------------------------------------------------

      if(loncr(3).gt.nlon) then
      loncr(1) = nlon/4
      loncr(2) = nlon/2
      loncr(3) = 3*nlon/4
      print*, "Longitudes of cross-sections redefined because some values defined in file"
      print*, "postmoloch.inp are out of the model domain"
      endif
      if(latcr(3).gt.nlat) then
      latcr(1) = nlat/4
      latcr(2) = nlat/2
      latcr(3) = 3*nlat/4
      print*, "Latitudes of cross-sections redefined because some values defined in file"
      print*, "postmoloch.inp are out of the model domain"
      endif

! -----------------------------

! Definition of the number of "free" cross-sections

  cross_number = 0
  do icross = 1,10
    if (lonini_cross(icross) == valmiss.or.latini_cross(icross) == valmiss.or. &
lonfin_cross(icross) == valmiss.and.latfin_cross(icross) == valmiss) exit
    cross_number = cross_number+1
    if (npoint_cross(icross) == 0) then
      dist_deg=sqrt((lonini_cross(icross)-lonfin_cross(icross))**2+ &
 (latini_cross(icross)-latfin_cross(icross))**2)
      npoint_cross(icross)=nint(dist_deg/pdr(1))+1
    endif
    npoint_cross(icross) = min( npoint_cross(icross), npoint_cross_max )

! Cross-section "vector" must be directed from west to east: it is important for interpretation
! of direct. of tangential and normal wind components relative to cross-section vector

    if (lonini_cross(icross) > lonfin_cross(icross)) then
      zzz=lonfin_cross(icross)
      lonfin_cross(icross)=lonini_cross(icross)
      lonini_cross(icross)=zzz
      zzz=latfin_cross(icross)
      latfin_cross(icross)=latini_cross(icross)
      latini_cross(icross)=zzz
    endif

  enddo

! -----------------------------
!---------------------------------------------------------------------------------------------------------------
!  Array dynamic allocation
!---------------------------------------------------------------------------------------------------------------

      allocate(   fmyu (nlat)  , stat=ierr)
      allocate(     fz (nlev)  , stat=ierr)
      allocate(  zlnph (nlev)  , stat=ierr)
      allocate(    fzh (nlevp1), stat=ierr)
      allocate(   zaux (nlevp1), stat=ierr)
      allocate(   zauxp(nlevp1), stat=ierr)
      allocate(   zlnp (nlevp1), stat=ierr)
      allocate(   zzzg (nlevg ), stat=ierr)
      allocate(   zztg (nlevg ), stat=ierr)
      allocate(   zlcr (nlivz ), stat=ierr)

      allocate(  izout (nlon, nlat ), stat=ierr)
      allocate( izoutx (nlon, nlivz), stat=ierr)
      allocate( izouty (nlat, nlivz), stat=ierr)

      allocate(  ztcrx (nlon, nlivz, ncrossx), stat=ierr)
      allocate( ztecrx (nlon, nlivz, ncrossx), stat=ierr)
      allocate(  zucrx (nlon, nlivz, ncrossx), stat=ierr)
      allocate(  zvcrx (nlon, nlivz, ncrossx), stat=ierr)
      allocate(  zwcrx (nlon, nlivz, ncrossx), stat=ierr)
      allocate( zrhcrx (nlon, nlivz, ncrossx), stat=ierr)
      allocate( zcwcrx (nlon, nlivz, ncrossx), stat=ierr)
      allocate( zcicrx (nlon, nlivz, ncrossx), stat=ierr)

      allocate(  ztcry (nlat, nlivz, ncrossy), stat=ierr)
      allocate( ztecry (nlat, nlivz, ncrossy), stat=ierr)
      allocate(  zucry (nlat, nlivz, ncrossy), stat=ierr)
      allocate(  zvcry (nlat, nlivz, ncrossy), stat=ierr)
      allocate(  zwcry (nlat, nlivz, ncrossy), stat=ierr)
      allocate( zrhcry (nlat, nlivz, ncrossy), stat=ierr)
      allocate( zcwcry (nlat, nlivz, ncrossy), stat=ierr)
      allocate( zcicry (nlat, nlivz, ncrossy), stat=ierr)

      allocate(   phig    (nlon, nlat), stat=ierr)
      allocate(  tskin    (nlon, nlat), stat=ierr)
      allocate(  qskin    (nlon, nlat), stat=ierr)
      allocate( cloudt    (nlon, nlat), stat=ierr)
      allocate( cloudl    (nlon, nlat), stat=ierr)
      allocate( cloudm    (nlon, nlat), stat=ierr)
      allocate( cloudh    (nlon, nlat), stat=ierr)
      allocate( prectot   (nlon, nlat), stat=ierr)
      allocate( precconv  (nlon, nlat), stat=ierr)
      allocate( precsolid (nlon, nlat), stat=ierr)
      allocate( runoff    (nlon, nlat), stat=ierr)
      allocate( runoff_tot(nlon, nlat), stat=ierr)
      allocate(   fice    (nlon, nlat), stat=ierr)
      allocate(  iceth    (nlon, nlat), stat=ierr)
      allocate(   snow    (nlon, nlat), stat=ierr)
      allocate(  fsnow    (nlon, nlat), stat=ierr)
      allocate( albedo    (nlon, nlat), stat=ierr)
      allocate(    rgm    (nlon, nlat), stat=ierr)
      allocate(    rgq    (nlon, nlat), stat=ierr)
      allocate(  fmask    (nlon, nlat), stat=ierr)
      allocate(  emis1    (nlon, nlat), stat=ierr)
      allocate(  emis2    (nlon, nlat), stat=ierr)
      allocate(     hx    (nlon, nlat), stat=ierr)
      allocate(     hy    (nlon, nlat), stat=ierr)
      allocate(  alont    (nlon, nlat), stat=ierr)
      allocate(  alatt    (nlon, nlat), stat=ierr)
      allocate(  alonu    (nlon, nlat), stat=ierr)
      allocate(  alatu    (nlon, nlat), stat=ierr)
      allocate(  alonv    (nlon, nlat), stat=ierr)
      allocate(  alatv    (nlon, nlat), stat=ierr)
      allocate(  alonz    (nlon, nlat), stat=ierr)
      allocate(  alatz    (nlon, nlat), stat=ierr)
      allocate( fcorio    (nlon, nlat), stat=ierr)
      allocate( zorogr    (nlon, nlat), stat=ierr)
      allocate(     ps    (nlon, nlat), stat=ierr)
      allocate(   zpz0    (nlon, nlat), stat=ierr)
      allocate(    zus    (nlon, nlat), stat=ierr)
      allocate(    zvs    (nlon, nlat), stat=ierr)
      allocate(    zts    (nlon, nlat), stat=ierr)
      allocate(    zqs    (nlon, nlat), stat=ierr)
      allocate(  zthes    (nlon, nlat), stat=ierr)
      allocate( ztskin    (nlon, nlat), stat=ierr)
      allocate(  ztgr1    (nlon, nlat), stat=ierr)
      allocate(  ztgr2    (nlon, nlat), stat=ierr)
      allocate(  ztgr3    (nlon, nlat), stat=ierr)
      allocate(  ztgr4    (nlon, nlat), stat=ierr)
      allocate(  zqgr1    (nlon, nlat), stat=ierr)
      allocate(  zqgr2    (nlon, nlat), stat=ierr)
      allocate(  zqgr3    (nlon, nlat), stat=ierr)
      allocate(  zqgr4    (nlon, nlat), stat=ierr)
      allocate(  zwork    (nlon, nlat), stat=ierr)
      allocate(  zwor2    (nlon, nlat), stat=ierr)
      allocate(  zwor3    (nlon, nlat), stat=ierr)
      allocate(  zwor4    (nlon, nlat), stat=ierr)
      allocate(  zwor5    (nlon, nlat), stat=ierr)
      allocate(  hflux    (nlon, nlat), stat=ierr)
      allocate(  qflux    (nlon, nlat), stat=ierr)
      allocate(     t2    (nlon, nlat), stat=ierr)
      allocate(    td2    (nlon, nlat), stat=ierr)
      allocate(     q2    (nlon, nlat), stat=ierr)
      allocate(  q2rel    (nlon, nlat), stat=ierr)
      allocate(    u10    (nlon, nlat), stat=ierr)
      allocate(    v10    (nlon, nlat), stat=ierr)
      allocate(   cape    (nlon, nlat), stat=ierr)
      allocate(    cin    (nlon, nlat), stat=ierr)
      allocate(   cink    (nlon, nlat), stat=ierr)
      allocate(  capek    (nlon, nlat), stat=ierr)
      allocate(  zlift    (nlon, nlat), stat=ierr)
      allocate(    t05    (nlon, nlat), stat=ierr)
      allocate(   t005    (nlon, nlat), stat=ierr)
      allocate( q05rel    (nlon, nlat), stat=ierr)
      allocate(  tg005    (nlon, nlat), stat=ierr)
      allocate(  tg010    (nlon, nlat), stat=ierr)
      allocate(  tg020    (nlon, nlat), stat=ierr)
      allocate(  tg050    (nlon, nlat), stat=ierr)
      allocate(  tg100    (nlon, nlat), stat=ierr)
      allocate(  cswfl    (nlon, nlat), stat=ierr)
      allocate(  clwfl    (nlon, nlat), stat=ierr)
      allocate( chflux    (nlon, nlat), stat=ierr)
      allocate( cqflux    (nlon, nlat), stat=ierr)
      allocate(  t2min    (nlon, nlat), stat=ierr)
      allocate(  t2max    (nlon, nlat), stat=ierr)
      allocate(ws10max    (nlon, nlat), stat=ierr)
      allocate(ztvirtup   (nlon, nlat), stat=ierr)
      allocate(  ztarr    (nlon, nlat), stat=ierr)
      allocate(  zqarr    (nlon, nlat), stat=ierr)
      allocate(    iwv    (nlon, nlat), stat=ierr)
      allocate(    pbl    (nlon, nlat), stat=ierr)
      allocate(   gust    (nlon, nlat), stat=ierr)
      allocate(  zlev0    (nlon, nlat), stat=ierr)
      allocate(  richs    (nlon,nlat),  stat=ierr)
      allocate(lapse_rate_1(nlon, nlat), stat=ierr)
      allocate(lapse_rate_2(nlon, nlat), stat=ierr)
      allocate( tgsurf (nlon,nlat),  stat=ierr)
      allocate( qgsurf (nlon,nlat),  stat=ierr)
      allocate(fice_soil_surf(nlon,nlat),  stat=ierr)
      allocate( snow_albedo (nlon,nlat),  stat=ierr)
      allocate( cwvflux (nlon,nlat),  stat=ierr)
      allocate( csoilh_bott_flux (nlon,nlat),  stat=ierr)
      allocate( csoilw_bott_flux (nlon,nlat),  stat=ierr)
      allocate( snowfall_level   (nlon,nlat),  stat=ierr)
      allocate( snowfall_height0 (nlon,nlat),  stat=ierr)
      allocate( snowfall_height  (nlon,nlat),  stat=ierr)

      allocate( tg_post(nlon, nlat, nlevg), stat=ierr)
      allocate( qg_post(nlon, nlat, nlevg), stat=ierr)
      allocate( qgmax  (nlon, nlat, nlevg), stat=ierr)
      allocate( qgmin  (nlon, nlat, nlevg), stat=ierr)

      allocate(     u (nlon, nlat, nlev), stat=ierr)
      allocate(     v (nlon, nlat, nlev), stat=ierr)
      allocate(     t (nlon, nlat, nlev), stat=ierr)
      allocate(     p (nlon, nlat, nlev), stat=ierr)
      allocate(     q (nlon, nlat, nlev), stat=ierr)
      allocate(   qcw (nlon, nlat, nlev), stat=ierr)
      allocate(   qci (nlon, nlat, nlev), stat=ierr)
      allocate( zthee (nlon, nlat, nlev), stat=ierr)
      allocate( zrelh (nlon, nlat, nlev), stat=ierr)
      allocate( tvirt (nlon, nlat, nlev), stat=ierr)
      allocate(tvlift (nlon, nlat, nlev), stat=ierr)
      allocate(  zeta (nlon, nlat, nlev), stat=ierr)
      allocate(  teta (nlon, nlat, nlev), stat=ierr)
      allocate(   fmz (nlon, nlat, nlev), stat=ierr)

      allocate(    tg (nlon, nlat, nlevg), stat=ierr)
      allocate(    qg (nlon, nlat, nlevg), stat=ierr)
      allocate( zparr (nlon, nlat, 2    ), stat=ierr)
      allocate(ztvlift(nlon, nlat, 2    ), stat=ierr)
      allocate(fice_soil (nlon, nlat, nlevg), stat=ierr)
      allocate(     snow_lev (nlon, nlat, nlev_snow), stat=ierr)
      allocate(       snow_t (nlon, nlat, nlev_snow), stat=ierr)
      allocate(    snow_fice (nlon, nlat, nlev_snow), stat=ierr)
      allocate(     snow_age (nlon, nlat, nlev_snow), stat=ierr)
      allocate(snow_melt_age (nlon, nlat, nlev_snow), stat=ierr)
      allocate(    snow_dens (nlon, nlat, nlev_snow), stat=ierr)
      allocate(snow_lev_depth(nlon, nlat, nlev_snow), stat=ierr)
      allocate(snow_lev_depth_prod (nlon, nlat, nlev_snow_prod), stat=ierr)
      allocate(        snow_t_prod (nlon, nlat, nlev_snow_prod), stat=ierr)
      allocate(      snow_age_prod (nlon, nlat, nlev_snow_prod), stat=ierr)
      allocate(     snow_dens_prod (nlon, nlat, nlev_snow_prod), stat=ierr)

      allocate(   w (nlon, nlat, nlevp1), stat=ierr)
      allocate( fmzh (nlon, nlat, nlevp1), stat=ierr)
      allocate(potvor(nlon, nlat, nlevp1), stat=ierr)

      allocate(   zph (nlon, nlat, nlevpo), stat=ierr)
      allocate(    zu (nlon, nlat, nlevpo), stat=ierr)
      allocate(    zv (nlon, nlat, nlevpo), stat=ierr)
      allocate(    zt (nlon, nlat, nlevpo), stat=ierr)
      allocate(    zq (nlon, nlat, nlevpo), stat=ierr)
      allocate(    zw (nlon, nlat, nlevpo), stat=ierr)
      allocate(  zqcw (nlon, nlat, nlevpo), stat=ierr)
      allocate(  zqci (nlon, nlat, nlevpo), stat=ierr)
      allocate( zclwi (nlon, nlat, nlevpo), stat=ierr)
      allocate(   zrh (nlon, nlat, nlevpo), stat=ierr)
      allocate(   zpv (nlon, nlat, nlevpo), stat=ierr)
      allocate(zthetae(nlon, nlat, nlevpo), stat=ierr)

      allocate(prectotr   (nlon, nlat), stat=ierr)
      allocate(precconvr  (nlon, nlat), stat=ierr)
      allocate(precsolidr (nlon, nlat), stat=ierr)
      allocate(runoffr (nlon, nlat), stat=ierr)
      allocate(runoff_tot_r (nlon, nlat), stat=ierr)
      allocate( cswflr (nlon, nlat), stat=ierr)
      allocate( clwflr (nlon, nlat), stat=ierr)
      allocate(chfluxr (nlon, nlat), stat=ierr)
      allocate(cqfluxr (nlon, nlat), stat=ierr)
      allocate( t2minr (nlon, nlat), stat=ierr)
      allocate( t2maxr (nlon, nlat), stat=ierr)
      allocate(ws10maxr(nlon, nlat), stat=ierr)
      allocate(cwvfluxr(nlon, nlat), stat=ierr)
      allocate(csoilh_bott_fluxr(nlon, nlat), stat=ierr)
      allocate(csoilw_bott_fluxr(nlon, nlat), stat=ierr)

      allocate(lon_crossy(ncrossy,2), stat=ierr)
      allocate(lat_crossy(ncrossy,2), stat=ierr)
      allocate(lon_crossx(ncrossx,2), stat=ierr)
      allocate(lat_crossx(ncrossx,2), stat=ierr)
      allocate(orog_crossy(nlat,ncrossx), stat=ierr)
      allocate(orog_crossx(nlon,ncrossx), stat=ierr)

      allocate(orog_cross(npoint_cross_max,cross_number), stat=ierr)

      allocate(t_cross      (npoint_cross_max,nlivz,cross_number), stat=ierr)
      allocate(theta_cross  (npoint_cross_max,nlivz,cross_number), stat=ierr)
      allocate(thetae_cross (npoint_cross_max,nlivz,cross_number), stat=ierr)
      allocate(u_cross      (npoint_cross_max,nlivz,cross_number), stat=ierr)
      allocate(v_cross      (npoint_cross_max,nlivz,cross_number), stat=ierr)
      allocate(u_tang_cross (npoint_cross_max,nlivz,cross_number), stat=ierr)
      allocate(v_norm_cross (npoint_cross_max,nlivz,cross_number), stat=ierr)
      allocate(w_cross      (npoint_cross_max,nlivz,cross_number), stat=ierr)
      allocate(rh_cross     (npoint_cross_max,nlivz,cross_number), stat=ierr)
      allocate(qcw_cross    (npoint_cross_max,nlivz,cross_number), stat=ierr)
      allocate(qci_cross    (npoint_cross_max,nlivz,cross_number), stat=ierr)

! ---- For RTTOV simulation ---

#ifdef rttov
 allocate (p_invert(nlon,nlat,nlev))
 allocate (t_invert(nlon,nlat,nlev))
 allocate (q_invert(nlon,nlat,nlev))
 allocate (qcw_invert(nlon,nlat,nlev))
 allocate (qci_invert(nlon,nlat,nlev))
 allocate (clfrac_invert(nlon,nlat,nlev))
 allocate (fseaice(nlon,nlat))
 allocate (flag_cloud_conv(nlon,nlat))
 allocate (sensor_chan_id(nchan))
 sensor_chan_id(:) = int(val_missing)
 allocate (sensor_chan_cw(nchan))
 sensor_chan_cw(:) = val_missing
 allocate (radiance(nlon,nlat,nchan))
 radiance(:,:,:) = val_missing
 allocate (radiance_bt(nlon,nlat,nchan))
 radiance_bt(:,:,:) = val_missing
 allocate (radiance_refl(nlon,nlat,nchan))
 radiance_refl(:,:,:) = val_missing
 allocate (radiance_clear(nlon,nlat,nchan))
 radiance_clear(:,:,:) = val_missing
 allocate (radiance_clear_bt(nlon,nlat,nchan))
 radiance_clear_bt(:,:,:) = val_missing
 allocate (radiance_clear_refl(nlon,nlat,nchan))
 radiance_clear_refl(:,:,:) = val_missing
 allocate (emis_rttov(nlon,nlat,nchan))
 emis_rttov(:,:,:) = val_missing
#endif

! -----------------------------
!---------------------------------------------------------------------------------------------------------------

      iist = 0
      iist2= 0

      prectot(:,:) = 0.
      precconv(:,:) = 0.
      precsolid(:,:) = 0.
      runoff(:,:) = 0.
      runoff_tot(:,:) = 0.
      cswfl(:,:)  = 0.
      clwfl(:,:)  = 0.
      chflux(:,:) = 0.
      cqflux(:,:) = 0.
      t2min(:,:)  = 999.
      t2max(:,:)  = 0.
      ws10max(:,:)= 0.
      cwvflux(:,:)= 0.
      csoilh_bott_flux(:,:)= 0.
      csoilw_bott_flux(:,:)= 0.

!  Define model and physical parameters

      dlat   = pdr(1)
      dlon   = pdr(2)
      dtstep = pdr(3)
      alat0  = pdr(4)
      alon0  = pdr(5)
      x0     = pdr(39)
      y0     = pdr(38)
      h      = pdr(40)
      qccrit = pdr(41)
      b0     = pdr(42)
      a0     = pdr(43)

      if(qccrit.lt.1.e-5.or.qccrit.gt.1.e-3) then
      qccrit = 2.2e-4 ! for compatibility with old mhf, missing definition of qccrit
      print*, "qccrit changed, defined in postmoloch"
      endif

      hsnc   = .02
      rd     = 287.05
      rv     = 461.51
      eps    = rd/rv
      ep     = 1./eps-1.
      cpd    = 1004.6
      cvd    = cpd-rd
      cpv    = 1869.46
      gamma  = cpd/cvd
      rdrcp  = rd/cpd
      g      = 9.807
      a      = 6371.e+3
      pi     = abs(acos(-1.))
      ome    = 7.292e-5
      dx     = a*dlon*pi/180.
      dy     = a*dlat*pi/180.
      dz     = h/float(nlev)

! Readinf of constant (in time) model physiographical parameters
! (phig, fmask)

      call rd_param_const(nlon, nlat, nlevg, x0, y0, alon0, alat0+dlat*0.5, dlon, dlat, phig, fmask, qgmax, qgmin)

 do while (.true.)

! Read postprocessing file

      print *
      print *, 'Read ',iist+1,' instant from moloch_atm.mhf and moloch_soil.mhf'

      call rdmhf_atm(nunic1, nlon, nlat, nlev, nfdr, pdr, p, u, v, w, t, q, qcw, qci)
      if (ierr_read1 /= 0) exit

      call rdmhf_soil(nunic2, nlon, nlat, nlevg, nlev_snow, nfdr, pdr, &
                      rgm, rgq, fice, iceth, albedo, emis1, emis2, &
                      cloudt, prectotr, precconvr, precsolidr, &
                      tskin, tgsurf, tg, qskin, qgsurf, qg, fice_soil_surf, fice_soil, &
                      snow, snow_lev, snow_t, snow_fice, snow_age, snow_melt_age, snow_dens, snow_albedo, &
                      cswflr, clwflr, chfluxr, cqfluxr, t2minr, t2maxr, ws10maxr, runoffr, runoff_tot_r, &
                      cwvfluxr, csoilh_bott_fluxr, csoilw_bott_fluxr, snowfall_level)
      if (ierr_read2 /= 0) exit

      ivalt = nfdr(10)*10000+nfdr(11)*100+nfdr(12)
      nhist=nfdr(17)

      if(iist == 0) then
       if(ivalt == 0) print*, "The 1st processed instant is an initial condition"
       if(njump /= 1) print '(A8,I2,A53)'," njump =",njump,": note that the 1st forecast instant is not processed"
      endif

      q   = max(q, 0.)
      qcw = max(qcw, 0.)
      qci = max(qci, 0.)

! Computation of accumulated precipitation

      iist=iist+1

      do jlat=1,nlat
      do jlon=1,nlon
      prectot(jlon,jlat)=prectot(jlon,jlat)+prectotr(jlon,jlat)
      precconv(jlon,jlat)=precconv(jlon,jlat)+precconvr(jlon,jlat)
      precsolid(jlon,jlat)=precsolid(jlon,jlat)+precsolidr(jlon,jlat)
      runoff(jlon,jlat)=runoff(jlon,jlat)+runoffr(jlon,jlat)
      runoff_tot(jlon,jlat)=runoff_tot(jlon,jlat)+runoff_tot_r(jlon,jlat)
      cswfl(jlon,jlat)=cswfl(jlon,jlat)+cswflr(jlon,jlat)
      clwfl(jlon,jlat)=clwfl(jlon,jlat)+clwflr(jlon,jlat)
      chflux(jlon,jlat)=chflux(jlon,jlat)+chfluxr(jlon,jlat)
      cqflux(jlon,jlat)=cqflux(jlon,jlat)+cqfluxr(jlon,jlat)
      if (t2minr(jlon,jlat)<t2min(jlon,jlat)) t2min(jlon,jlat)=t2minr(jlon,jlat)
      if (t2maxr(jlon,jlat)>t2max(jlon,jlat)) t2max(jlon,jlat)=t2maxr(jlon,jlat)
      if (ws10maxr(jlon,jlat)>ws10max(jlon,jlat)) ws10max(jlon,jlat)=ws10maxr(jlon,jlat)
      cwvflux(jlon,jlat)=cwvflux(jlon,jlat)+cwvfluxr(jlon,jlat)
      csoilh_bott_flux(jlon,jlat)=csoilh_bott_flux(jlon,jlat)+csoilh_bott_fluxr(jlon,jlat)
      csoilw_bott_flux(jlon,jlat)=csoilw_bott_flux(jlon,jlat)+csoilw_bott_fluxr(jlon,jlat)
      enddo
      enddo

      if (iist==1) then
        if (ivalt==0) then
          iist0=-1 ! validity time of first field=0: analysis
        else
          iist0=0
        endif
      endif

      if(mod(iist+iist0,njump) == 0) then
      iist2=iist2+1

!---------------------------------------------------------------------------------------------------------------
!  Definition of orography

      zorogr = phig/g

!  Pressure check&fix (to avoid problems due to pressure increasing or nearly constant with height)

    do jlat = 1, nlat
    do jlon = 1, nlon
     do jklev = nlev-2, 1, -1
     if (p(jlon,jlat,jklev).le.p(jlon,jlat,jklev+1)+1.e2) then
     p(jlon,jlat,jklev) = p(jlon,jlat,jklev+1) + 1.e2
     endif
     enddo
    enddo
    enddo

!---------------------------------------------------------------------------------------------------------------
!  fmy = 1./cos(latitude)
!  fz = 1.-zita/h
!  alatt, alont = lat. and lon. of t-points
!  alatz = lat. of vorticity points
!  fcorio = Coriolis param. (on vorticity points only)
!---------------------------------------------------------------------------------------------------------------

      call rot_grid(x0,y0,alon0,                     alat0+dlat*0.5/float(mhfr),dlon,dlat,alont,alatt,nlon,nlat)
      call rot_grid(x0,y0,alon0+dlon*0.5/float(mhfr),alat0+dlat*0.5/float(mhfr),dlon,dlat,alonu,alatu,nlon,nlat)
      call rot_grid(x0,y0,alon0,                     alat0,                     dlon,dlat,alonv,alatv,nlon,nlat)
      call rot_grid(x0,y0,alon0+dlon*0.5/float(mhfr),alat0,                     dlon,dlat,alonz,alatz,nlon,nlat)

      zfac = pi/180.
      do jlat = 1, nlat
      zlatt = (alat0 + dlat*.5/float(mhfr) +(jlat-1)*dlat)*zfac
      fmyu(jlat) = 1./cos(zlatt)
      do jlon = 1, nlon
      fcorio(jlon,jlat) = 2.*ome*sin(pi*alatz(jlon,jlat)/180.)
      enddo
      enddo

      do jklev=1,nlev
      zita  =(jklev-1)*dz+dz/2.
      zitah =(jklev-1)*dz
      fz (jklev) =1.-zita /h
      fzh(jklev) =1.-zitah/h
      do jlat=1,nlat
      do jlon=1,nlon
      zeta(jlon,jlat,jklev) = zorogr(jlon,jlat)*gzita(zita,h,a0) -h*bzita(zita,h,b0)*log(fz(jklev))
      fmz (jlon,jlat,jklev) = fz(jklev)/(bzita(zita,h,b0)+zorogr(jlon,jlat)*fz(jklev)*gzitap(zita,h,a0) &
                          - h*fz(jklev)*log(fz(jklev))*bzitap(zita,h,b0))
      fmzh(jlon,jlat,jklev) = fzh(jklev)/(bzita(zitah,h,b0)+zorogr(jlon,jlat)*fzh(jklev)*gzitap(zitah,h,a0) &
                              - h*fzh(jklev)*log(fzh(jklev))*bzitap(zitah,h,b0))
      enddo
      enddo
      enddo
      fzh(nlevp1) = 0.
      fmzh(:,:,nlevp1) = 0.

      do jlat=1,nlat
      jlatm1=max(jlat-1,1)
      do jlon=1,nlon
      jlonp1=min(jlon+1,nlon)
      hx(jlon,jlat)=(phig(jlonp1,jlat)-phig(jlon,jlat))/g/dx*fmyu(jlat)
      hy(jlon,jlat)=(phig(jlon,jlat)-phig(jlon,jlatm1))/g/dy
      enddo
      enddo

!  Definition of levels of cross-sections (regular z-coordinates)

      do jklev=1,nlivz
      zlcr(jklev) = (zctop-zcbot)/float(nlivz-1)*(jklev-1)
      enddo

!---------------------------------------------------------------------------------------------------------------
!  Definition of teta, thetae and relative humidity relh
!---------------------------------------------------------------------------------------------------------------
!  Constants used to compute partial pressure at saturation

      zcw1 = (cpv-cw)/rv
      zci1 = (cpv-ci)/rv
      zcw2 = ylwv/tzer/rv-(cpv-cw)/rv
      zci2 = yliv/tzer/rv-(cpv-ci)/rv

      do jklev =1, nlev
      do jlat = 1, nlat
      do jlon = 1, nlon

      zf = max(1.e-8, q(jlon,jlat,jklev))
      call comp_esk(zzpvs,zfs,t(jlon,jlat,jklev),p(jlon,jlat,jklev),1) ! Comp. saturation to water and ice separately

!  3-d relative humidity for cross-sections

      zzpv=(zf*p(jlon,jlat,jklev))/(eps + zf -eps*zf)
      zrelh(jlon,jlat,jklev) = zzpv/zzpvs*100.

!  3-d equivalent pot. temperature (Bolton, MWR, 1980)
!  ztl: temp. at condens. lev. (computed with expr. by Bolton)

      zmr = zf/(1.-zf)
      ze = p(jlon,jlat,jklev)/100. * zmr/(eps+zmr)
      ztl = 55. + 2840./(3.5*alog(t(jlon,jlat,jklev)) - alog(ze) - 4.805)
      zespon = (3.376/ztl - 0.00254)*zmr*1000*(1.+.81*zmr)
      teta (jlon,jlat,jklev) = t(jlon,jlat,jklev)*(pzer/p(jlon,jlat,jklev))**rdrcp
      zthee(jlon,jlat,jklev) = teta(jlon,jlat,jklev)*exp(zespon)

      enddo
      enddo
      enddo

!---------------------------------------------------------------------------------------------------------------
!  Tvirt and CAPE
!---------------------------------------------------------------------------------------------------------------

     tvirt = t * (1.+ep*q-qcw-qci)

     if(isflag(30).eq.1) then

      cape = 0.
      cin = 0.

      jk1 = nlev-7     ! top level for CAPE computation

      do jk0 = 1, 2, 1   ! loop over starting level

      call lift_parcel(q(:,:,jk0), p, t(:,:,jk0), tvlift, nlon, nlat, nlev, jk0, jk1, 0)

      do jlat = 1, nlat
      do jlon = 1, nlon
      capek(jlon,jlat) = 0.
      cink(jlon,jlat) = 0.
      do jklev = jk0, jk1-1
      zcap=.5*rd*( (tvlift(jlon,jlat,jklev  )-tvirt(jlon,jlat,jklev  ))/p(jlon,jlat,jklev  ) +  &
                   (tvlift(jlon,jlat,jklev+1)-tvirt(jlon,jlat,jklev+1))/p(jlon,jlat,jklev+1) )* &
                   (p(jlon,jlat,jklev)-p(jlon,jlat,jklev+1))
      capek(jlon,jlat)=capek(jlon,jlat)+max(zcap,0.)  ! sum over positive buoyancy only (no convective inhibition)
      if (capek(jlon,jlat).lt.20.) then   ! important - replaces condition on level of free convection
      cink(jlon,jlat) = cink(jlon,jlat) + min(zcap,0.)
      endif
      enddo

      enddo
      enddo
      cape = max (cape, capek)  ! save only maximum CAPE over starting level
      cin  = min (cin, cink)    ! save only minimum CIN  over starting level
      enddo   ! end loop over starting level

      do jlat = 1, nlat
      do jlon = 1, nlon
      if (cape(jlon,jlat).lt.  50.) cin(jlon,jlat) = -999.
      if (cin (jlon,jlat).lt.-400.) cin(jlon,jlat) = -999.
      enddo
      enddo

      endif

!---------------------------------------------------------------------------------------------------------------
! Level of zero (above sea level)
!---------------------------------------------------------------------------------------------------------------

      do jlat = 1, nlat
      do jlon = 1, nlon

        zlev0(jlon,jlat) = -99999.

        if (t(jlon,jlat,nlev-nlev/10) > tzer) then
          zlev0(jlon,jlat) = 10000.
        else
          do jklev = nlev-nlev/10, 2, -1
            if (t(jlon,jlat,jklev).le.tzer.and.t(jlon,jlat,jklev-1).gt.tzer) then  ! linear extrapolation down
              zlev0(jlon,jlat) = zeta(jlon,jlat,jklev-1) + (zeta(jlon,jlat,jklev)-zeta(jlon,jlat,jklev-1)) * &
                         (t(jlon,jlat,jklev-1)-tzer)/(t(jlon,jlat,jklev-1)-t(jlon,jlat,jklev))
              exit
            endif
          enddo
        endif

! Case in which the level of 0 is below ground: extrap. from second level with standard lapse rate

        if (zlev0(jlon,jlat).eq.-99999.) then
          zlev0(jlon,jlat) = zeta(jlon,jlat,2) - (tzer-t(jlon,jlat,2))/6.5e-3
        endif

      enddo
      enddo

      zlev0 = max (zlev0, -50.)

!---------------------------------------------------------------------------------------------------------------
! Integrated Water Vapour IWV
!---------------------------------------------------------------------------------------------------------------

      do jlat=1,nlat
      do jlon=1,nlon
      iwv (jlon,jlat) = 0.
      do jk=1,nlev
      iwv(jlon,jlat)=iwv(jlon,jlat) + q(jlon,jlat,jk)*p(jlon,jlat,jk)/(rd*tvirt(jlon,jlat,jk)*fmz(jlon,jlat,jk))*dz
      enddo
      enddo
      enddo

!---------------------------------------------------------------------------------------------------------------
!  Lifted index
!---------------------------------------------------------------------------------------------------------------

! Definition of lifted index: virtual temp. difference between a parcel near about 500 hPa (average
! between two model levels, one above and one below 500 hPa) and a parcel
! lifted pseudo-adiabatically, starting in a layer of depth of zpdepth above the surface.
! Since ps is variable over orography, the 500 hPa level is actually arbitrarily defined as
! a pressure 500.e2 - (1005.e2-p(1))/4. - zpdepth also depends on the surface pressure.

      if(isflag(50).eq.1) then
      zlift(:,:)=999.
      ztarr(:,:)=0.
      zqarr(:,:)=0.
      ztvlift(:,:,:)=0.
      zparr(:,:,:)=0.

! Average of pot. temp. and specif. humidity in the pressure layer above the first model level
! (regardless of model layer depths, that in MOLOCH are roughly proportional to mass)

      do jlat=1,nlat
      do jlon=1,nlon

      jklol=-1
      jkupl=-1
      do jk=1,nlev
      zp0p=p(jlon,jlat,jk)
      zpdepth=85.e2*sqrt(p(jlon,jlat,1)/pzer) ! decreases the press. depth of mixing lower layer over topography

      if(zp0p.lt.p(jlon,jlat,1)-zpdepth.and.jklol.lt.0) jklol=jk-1 ! index of the top of the lower layer
       if(zp0p.lt.(500.e2-(pzer-p(jlon,jlat,1))/4.).and.jkupl.lt.0) then
       jkupl=jk ! index of the first layer above the upper (near 500 hPa) level
       goto 7890
       endif
      enddo

! Average of the pot. temp., spec. hum. and press. in the lower layer

7890  zpavlo=0.
      ztavlo=0.
      zqavlo=0.
      do jk=1,jklol
      zpavlo=zpavlo+p(jlon,jlat,jk)
      ztavlo=ztavlo+teta(jlon,jlat,jk)
      zqavlo=zqavlo+q(jlon,jlat,jk)
      enddo
      zdz=float(jklol)
      zpavlo=zpavlo/zdz
      ztavlo=ztavlo/zdz
      zqavlo=zqavlo/zdz
      ztempavlo=ztavlo*(zpavlo/pzer)**rdrcp  ! from theta to temp.

! Average of virtual temperature and press. in the upper layer

      ztvirtup(jlon,jlat)=0.5*(tvirt(jlon,jlat,jkupl)+tvirt(jlon,jlat,jkupl-1))
      zpavup=0.5*(p(jlon,jlat,jkupl)+p(jlon,jlat,jkupl-1))
      zqarr(jlon,jlat)=zqavlo
      ztarr(jlon,jlat)=ztempavlo
      zparr(jlon,jlat,1)=zpavlo
      zparr(jlon,jlat,2)=zpavup
      enddo
      enddo

      call lift_parcel(zqarr, zparr, ztarr, ztvlift, nlon, nlat, 2, 1, 2, 0)
      zlift(:,:)=ztvlift(:,:,2)-ztvirtup(:,:)

      endif ! iflag 50

!---------------------------------------------------------------------------------------------------------------
!  Surface fields
!---------------------------------------------------------------------------------------------------------------

      do jlat = 1, nlat-1
      do jlon = 2, nlon
      zus (jlon,jlat) = ((2.*float(mhfr)-1.)*u(jlon,jlat,1)+u(jlon-1,jlat,1))/(2.*float(mhfr))  ! destaggering
      zvs (jlon,jlat) = (v(jlon,jlat+1,1)+(2.*float(mhfr)-1.)*v(jlon,jlat,1))/(2.*float(mhfr))  ! destaggering
      enddo
      enddo

      zgam=6.0e-3
      do jlat=1,nlat
      jlatp1=min(jlat+1,nlat)
      do jlon=1,nlon
      jlonm1=max(jlon-1,1)

      zqtot=q(jlon,jlat,1)+qcw(jlon,jlat,1)+qci(jlon,jlat,1)
      tvirt(jlon,jlat,1)=t(jlon,jlat,1)*(1.-zqtot+q(jlon,jlat,1)/eps)

      zeta(jlon,jlat,1) = zorogr(jlon,jlat)*gzita(dz*.5,h,a0) -h*bzita(dz*.5,h,b0)*log(fz(1))
      ztz0=tvirt(jlon,jlat,1)+zgam*zeta(jlon,jlat,1)
      ztbar=.5*(ztz0+tvirt(jlon,jlat,1))
      zpz0(jlon,jlat)=p(jlon,jlat,1)/100.*exp(g*zeta(jlon,jlat,1)/rd/ztbar)

      zdz=zeta(jlon,jlat,1)-zorogr(jlon,jlat)
      ztbar=tvirt(jlon,jlat,1)
      ps(jlon,jlat)=p(jlon,jlat,1)*exp(g*zdz/rd/ztbar)

      zts (jlon,jlat) = t(jlon,jlat,1)-tzer
      zthes(jlon,jlat)= zthee(jlon,jlat,1)
      zqs(jlon,jlat)  = q(jlon,jlat,1)
      ztskin(jlon,jlat)= tskin(jlon,jlat) - tzer
      ztgr1(jlon,jlat) = tg(jlon,jlat,1) - tzer
!      ztgr2(jlon,jlat) = tg(jlon,jlat,2) - tzer
!      ztgr3(jlon,jlat) = tg(jlon,jlat,3) - tzer
!      ztgr4(jlon,jlat) = tg(jlon,jlat,4) - tzer
      ztgr2(jlon,jlat) = tg(jlon,jlat,3) - tzer
      ztgr3(jlon,jlat) = tg(jlon,jlat,5) - tzer
      ztgr4(jlon,jlat) = tg(jlon,jlat,7) - tzer

      zqgr1(jlon,jlat) = qg(jlon,jlat,1)
!      zqgr2(jlon,jlat) = qg(jlon,jlat,2)
!      zqgr3(jlon,jlat) = qg(jlon,jlat,3)
!      zqgr4(jlon,jlat) = qg(jlon,jlat,4)
      zqgr2(jlon,jlat) = qg(jlon,jlat,3)
      zqgr3(jlon,jlat) = qg(jlon,jlat,5)
      zqgr4(jlon,jlat) = qg(jlon,jlat,7)
      enddo
      enddo

      tg_post(:,:,:) = tg(:,:,:) - tzer
      qg_post(:,:,:) = qg(:,:,:)

!  Def.of ground temperature at standard observation levels

      zzzg(1:nlevg)= pdr(6:5+nlevg)

      do jlat = 1, nlat
      do jlon = 1, nlon
      do ig=1,nlevg
      zztg(ig)=tg(jlon,jlat,ig)
      enddo
      zsq(1)=sqrt(0.05)
      call interp(0.,1.,1.,nlevg,sqrt(zzzg),zztg,zsq(1),tg005(jlon,jlat),1)
      zsq(1)=sqrt(0.10)
      call interp(0.,1.,1.,nlevg,sqrt(zzzg),zztg,zsq(1),tg010(jlon,jlat),1)
      zsq(1)=sqrt(0.20)
      call interp(0.,1.,1.,nlevg,sqrt(zzzg),zztg,zsq(1),tg020(jlon,jlat),1)
      zsq(1)=sqrt(0.50)
      call interp(0.,1.,1.,nlevg,sqrt(zzzg),zztg,zsq(1),tg050(jlon,jlat),1)
      zsq(1)=sqrt(1.00)
      call interp(0.,1.,1.,nlevg,sqrt(zzzg),zztg,zsq(1),tg100(jlon,jlat),1)
      enddo
      enddo

      do jlat = 1, nlat
      do jlon = 1, nlon
      tg005(jlon,jlat) = tg005(jlon,jlat)-tzer
      tg010(jlon,jlat) = tg010(jlon,jlat)-tzer
      tg020(jlon,jlat) = tg020(jlon,jlat)-tzer
      tg050(jlon,jlat) = tg050(jlon,jlat)-tzer
      tg100(jlon,jlat) = tg100(jlon,jlat)-tzer
      enddo
      enddo

!---------------------------------------------------------------------------------------------------------------
!  Snow cover products:
!  at 3 levels of snow products:
! 1 - snow surface, 3 - snow bottom, 2 - level of medium of geometrical snow depth
! variables: level depth (m), temperature (K), age (days), density (kg m^-3)
! grib2 codes:
! snow age 0-1-17, snow albedo 0-19-19, snow cover 0-1-42,
! snow depth 0-1-11, snow depth water equivalent 0-1-60,
! snow temperature (at snow top) 0-0-18, snow temperature 192-201-203 or 0-1-208
!  Code table 4.5 - Fixed surface types and units: 114 Snow level (Numeric)
!---------------------------------------------------------------------------------------------------------------

      snow_lev_depth(:,:,:)      = valmiss
      snow_lev_depth_prod(:,:,:) = valmiss
      snow_t_prod(:,:,:)         = valmiss
      snow_age_prod(:,:,:)       = valmiss
      snow_dens_prod(:,:,:)      = valmiss
      do jlat = 1, nlat
      do jlon = 1, nlon

        kbott=int(valmiss)
        do jklev=1,nlev_snow
          if (int(snow_lev(jlon,jlat,jklev)) == int(valmiss)) exit
        enddo
        kbott=jklev-1 

        if (kbott > 1 ) then ! snow cover exists

          snow_lev_depth_prod(jlon,jlat,:) = 0.
          snow_lev_depth(jlon,jlat,1) = 0.

          if (kbott > 2) then ! There are some snow layers

            do jklev=1,kbott-1
              snow_dens_h(jklev)=(snow_dens(jlon,jlat,jklev+1)+snow_dens(jlon,jlat,jklev))*0.5
              zdepth=(snow_lev(jlon,jlat,jklev+1)-snow_lev(jlon,jlat,jklev))/snow_dens_h(jklev)
              snow_lev_depth_prod(jlon,jlat,3)=snow_lev_depth_prod(jlon,jlat,3)+zdepth
              snow_lev_depth(jlon,jlat,jklev+1)=snow_lev_depth(jlon,jlat,jklev)+zdepth
            enddo
            
            snow_lev_depth_prod(jlon,jlat,2) = snow_lev_depth_prod(jlon,jlat,3)*0.5

            do jklev=1,kbott
              if (snow_lev_depth(jlon,jlat,jklev) >= snow_lev_depth_prod(jlon,jlat,2)) exit
            enddo
            kmedium=jklev

            zfrac=(snow_lev_depth_prod(jlon,jlat,2)-snow_lev_depth(jlon,jlat,kmedium-1))/ &
 (snow_lev_depth(jlon,jlat,kmedium)-snow_lev_depth(jlon,jlat,kmedium-1))

          else ! There is one snow layer only

            snow_dens_h(1)=(snow_dens(jlon,jlat,2)+snow_dens(jlon,jlat,1))*0.5
            snow_lev_depth_prod(jlon,jlat,3) = (snow_lev(jlon,jlat,2)-snow_lev(jlon,jlat,1))/snow_dens_h(1)
            snow_lev_depth_prod(jlon,jlat,2) = snow_lev_depth_prod(jlon,jlat,3)*0.5

            kmedium=2
            zfrac=0.5

          endif

          snow_t_prod(jlon,jlat,1) = snow_t(jlon,jlat,1)
          snow_t_prod(jlon,jlat,3) = snow_t(jlon,jlat,kbott)
          snow_t_prod(jlon,jlat,2) = snow_t(jlon,jlat,kmedium-1)*(1.-zfrac)+snow_t(jlon,jlat,kmedium)*zfrac

          snow_age_prod(jlon,jlat,1) = snow_age(jlon,jlat,1)
          snow_age_prod(jlon,jlat,3) = snow_age(jlon,jlat,kbott)
          snow_age_prod(jlon,jlat,2) = snow_age(jlon,jlat,kmedium-1)*(1.-zfrac)+snow_age(jlon,jlat,kmedium)*zfrac

          snow_dens_prod(jlon,jlat,1) = snow_dens(jlon,jlat,1)
          snow_dens_prod(jlon,jlat,3) = snow_dens(jlon,jlat,kbott)
          snow_dens_prod(jlon,jlat,2) = snow_dens(jlon,jlat,kmedium-1)*(1.-zfrac)+snow_dens(jlon,jlat,kmedium)*zfrac

        endif ! snow cover exists

      enddo
      enddo

!---------------------------------------------------------------------------------------------------------------
! Height of snow fall limit
     call snowfall_height_def (q, p, t, snowfall_level, zeta, zorogr, nlon, nlat, nlev, ylwv, cpd, &
 snowfall_height, snowfall_height0)

!---------------------------------------------------------------------------------------------------------------

!  Smoothing of zpz0
!  Filtering of msl pressure as a function of orographic height
!  (zfp is a filtering factor, empirically set)

      do 105 iter = 1, 25

      do jlat = 1, nlat
      jlatp1=min(jlat+1,nlat)
      jlatm1=max(jlat-1,1)
      do jlon = 1, nlon
      jlonp1 = min(jlon+1,nlon)
      jlonm1 = max(jlon-1,1)
      zfp = min (zorogr(jlon,jlat)/3500.,.8)
      zfp = max (zfp, 0.)
      zwork(jlon,jlat) = (1.-zfp) * zpz0(jlon,jlat) + zfp/4.*  &
           (zpz0(jlonm1,jlat) + zpz0(jlonp1,jlat) + zpz0(jlon,jlatm1) + zpz0(jlon,jlatp1))
      enddo
      enddo

      zpz0 = zwork

  105 continue

      do kk=1,2
      do jlat = 1, nlat
      do jlon = 2, nlon-1
      zwork(jlon,jlat)=.25*(zpz0(jlon-1,jlat)+zpz0(jlon+1,jlat))+.5*zpz0(jlon,jlat)
      enddo
      enddo
      do jlat = 2, nlat-1
      do jlon = 2, nlon-1
      zpz0(jlon,jlat)=.25*(zwork(jlon,jlat+1)+zwork(jlon,jlat-1))+.5*zwork(jlon,jlat)
      enddo
      enddo
      enddo

      do jlon=1,nlon
      do jlat=1,nlat
      fsnow(jlon,jlat)=min((snow(jlon,jlat)/hsnc)**.67,.9)
      enddo
      enddo
      call vtsurf(u10, v10, t2, q2, hflux, qflux, nlon, nlat, nlev, fmask, rgm, rgq, &
                  richs, phig, ps, fsnow, tskin, qskin, u, v, t, p, q, h, dz, a0, b0, mhfr)

!  Fix over high mountains.....

      do jlat = 1, nlat
      do jlon = 1, nlon
      zglin = max (0., min((phig(jlon,jlat)-11000.)/15000., 1.))
      u10(jlon,jlat) = (1.-zglin)*u10(jlon,jlat) + zglin*u(jlon,jlat,1)
      v10(jlon,jlat) = (1.-zglin)*v10(jlon,jlat) + zglin*v(jlon,jlat,1)
      if (t2(jlon,jlat)+tzer.lt.t(jlon,jlat,1)) then
      zglin = max (0., min((phig(jlon,jlat)-10000.)/25000., 1.))
      t2(jlon,jlat) = (1.-zglin)*t2(jlon,jlat) + zglin*(t(jlon,jlat,1)-tzer)
      q2(jlon,jlat) = (1.-zglin)*q2(jlon,jlat) + zglin*q(jlon,jlat,1)
      endif

!!  FIXFIX low wind
!
!    if (fmask(jlon,jlat).lt.0.1) then
!      if (tskin(jlon,jlat)-tzer.gt.t2(jlon,jlat)+.1) then
!      zglin = u10(jlon,jlat)**2 +v10(jlon,jlat)**2
!      zglin = max ((25.-zglin)/25., 0.)
!      zglin = min (zglin, .8)
!      t2(jlon,jlat) = t2(jlon,jlat) + min(zglin*(tskin(jlon,jlat)-tzer-t2(jlon,jlat)), 2.0)
!      endif
!    endif

!  Relative humidity at 2 m (with respect to water! - WMO definition)

      call comp_esk (zesk,zqsat,t2(jlon,jlat)+tzer,ps(jlon,jlat),3)  ! partial pressure over water
      q2p = min(q2(jlon,jlat), zqsat*1.01) ! ensures that q2rel does not exceed 101% (for graphic smoothness)
      eee = ps(jlon,jlat)*q2p/(eps*(1.-q2p)+q2p)
      q2rel(jlon,jlat) = eee/zesk*100.

! Dew point temperature at 2 m

      call td_tetens (t2(jlon,jlat)+tzer, ps(jlon,jlat), q2(jlon,jlat), eps, td2(jlon,jlat))
      td2(jlon,jlat) = td2(jlon,jlat)-tzer
      enddo
      enddo

!  Definition of clouds

      call cclouds (nlon, nlat, nlev, h, dz, a0, b0, qccrit, phig, richs, zeta, u, v, teta,    &
                    cloudt, cloudh, cloudm, cloudl, t, p, q, qcw, qci, tvirt, pbl, zorogr)

!  Definition of Potential Vorticity (PV)

      call pv (fmzh, u, v, w, teta, hx, hy, fmyu, dx, dy, dz, h, a0, p, rd, tvirt, &
               fcorio, nlon, nlat, nlev, nlevp1, potvor)

! Wind gust computation as a function of the PBL top height

      do jlat = 1, nlat
      do jlon = 1, nlon
      gust(jlon,jlat) = sqrt(u10(jlon,jlat)**2+v10(jlon,jlat)**2)
       do jklev = 1, 8
       jkup = jklev
       if (pbl(jlon,jlat).lt.zeta(jlon,jlat,jklev)-zorogr(jlon,jlat)) exit
       enddo
      if (jkup.gt.1) then
       do jklev = 1, jkup-1
       zwink = sqrt(u(jlon,jlat,jklev)**2+v(jlon,jlat,jklev)**2)
       gust(jlon,jlat) = max (gust(jlon,jlat), zwink)
       enddo
      endif
      enddo
      enddo

! Integrated Vapour Transport (care: overwritten on surface wind components zus, zvs)

!!      print*, "The integrated vapour transport vector is computed and written"
!!      print*, "in place of the wind vector at the lowest model level."
!!      do jlat=1,nlat
!!      jlatp1=min(jlat+1,nlat)
!!      do jlon=1,nlon
!!      jlonm1=max(jlon-1,1)
!!      zus(jlon,jlat) = 0.
!!      zvs(jlon,jlat) = 0.
!!      do jk=1,nlev
!!      if(p(jlon,jlat,jk).gt.30000.) then                         ! selection in terms of pressure
!!!      if(zeta(jlon,jlat,jk)-zorogr(jlon,jlat).lt.3001.) then   ! selection with height over orography
!!      zus(jlon,jlat)=zus(jlon,jlat) + q(jlon,jlat,jk)*.5*(u(jlon,jlat,jk)+u(jlonm1,jlat,jk))* &
!!                    p(jlon,jlat,jk)/(rd*tvirt(jlon,jlat,jk)*fmz(jlon,jlat,jk))*dz
!!      zvs(jlon,jlat)=zvs(jlon,jlat) + q(jlon,jlat,jk)*.5*(v(jlon,jlat,jk)+v(jlon,jlatp1,jk))* &
!!                    p(jlon,jlat,jk)/(rd*tvirt(jlon,jlat,jk)*fmz(jlon,jlat,jk))*dz
!!      endif
!!      enddo
!!      enddo
!!      enddo

!---------------------------------------------------------------------------------------------------------------
!  Interpolation of variables on constant pressure surfaces
!---------------------------------------------------------------------------------------------------------------

      do jlat=1,nlat
      jlatp1=min(jlat+1,nlat)
      do jlon=1,nlon
      jlonm1=max(jlon-1,1)

      zalf = .5
      do jklev=1,nlev
      zlnp (nlevp1-jklev)=alog(p(jlon,jlat,jklev))
      enddo
      zlnp(nlevp1)=alog(zpz0(jlon,jlat)*100.)

      do jklev = 1, nlev-1
      zlnph(nlevp1-jklev) = log(.5*(p(jlon,jlat,jklev)+p(jlon,jlat,jklev+1)))
      enddo
      zlnph(1) = log(.5*p(jlon,jlat,nlev))

!  Here and below, zauxp is used in input of interp due to possible changes done on this vector
!  in subr. interp in the case pressure is not monotonic;
!  when zauxp is changed in interp, also the dependent variables to be interpolated must be changed,
!  so zauxp must be refreshed every time interp is called

!  Interpolation of zeta (converted to geopotential meters)

      ze1 = 1.
      ze2 = 1.
      zauxp=zlnp
      do jklev= 1, nlev
      zaux(nlevp1-jklev)=zeta(jlon,jlat,jklev)*g/9.8
      enddo
      zaux(nlevp1)=0.

! In case zlnp(nlevp1) <= zlnp(nlev) (lowest model level pressure >= mslp), the auxiliary level
! at sea level pressure cannot be used - this concerns vert. interpolation of t and phi only

      if(zlnp(nlevp1).le.zlnp(nlev)) then
      nlevint=nlev
      else
      nlevint=nlevp1
      endif

      call interp(zalf, ze1, ze2, nlevint, zauxp, zaux, zalp, zph(jlon,jlat,1:nlevpo), nlevpo)

!  Interpolation of t (converted to Celsius)

      ze1 = .8
      ze2 = .8
      zauxp=zlnp
      do jklev= 2, nlev
      zaux(nlevp1-jklev)=t(jlon,jlat,jklev)-tzer
      enddo
      zaux(nlev  ) = t(jlon,jlat,2)-tzer + zgam*(zeta(jlon,jlat,2)-zeta(jlon,jlat,1)) ! sottocolle
      zaux(nlevp1) = t(jlon,jlat,2)-tzer + zgam*zeta(jlon,jlat,2)                     ! sottocolle
      call interp(zalf, ze1, ze2, nlevint, zauxp, zaux, zalp, zt(jlon,jlat,1:nlevpo), nlevpo)

!  Interpolation of q

      ze1 = .7
      ze2 = .7
      zauxp=zlnp
      do jklev= 1, nlev
      q(jlon,jlat,jklev)=max(q(jlon,jlat,jklev),1.e-9)
      zaux(nlevp1-jklev)=sqrt(q(jlon,jlat,jklev))
      enddo
      call interp(zalf, ze1, ze2, nlev, zauxp, zaux, zalp, zq(jlon,jlat,1:nlevpo), nlevpo)
      do jklev=1, nlevpo
      zq(jlon,jlat,jklev)=max(zq(jlon,jlat,jklev), 0.)
      zq(jlon,jlat,jklev)=zq(jlon,jlat,jklev)**2
      enddo

!  Interpolation of qcw

      zauxp=zlnp
      do jklev= 1, nlev
      qcw(jlon,jlat,jklev)=max(qcw(jlon,jlat,jklev),1.e-12)
      zaux(nlevp1-jklev)=qcw(jlon,jlat,jklev)
      enddo
      call interp(zalf, ze1, ze2, nlev, zauxp, zaux, zalp, zqcw(jlon,jlat,1:nlevpo), nlevpo)

!  Interpolation of qci

      zauxp=zlnp
      do jklev= 1, nlev
      qci(jlon,jlat,jklev)=max(qci(jlon,jlat,jklev), 1.e-12)
      zaux(nlevp1-jklev)=qci(jlon,jlat,jklev)
      enddo
      call interp(zalf, ze1, ze2, nlev, zauxp, zaux, zalp, zqci(jlon,jlat,1:nlevpo), nlevpo)

!  Interpolation of u (averaged on t points)

      ze1 = .4
      ze2 = .4
      if (jlon > 1) then
        zauxp = zlnp
        do jklev = 1, nlev
        zaux(nlevp1-jklev) = ((2.*float(mhfr)-1.)*u(jlon,jlat,jklev)+u(jlon-1,jlat,jklev))/(2.*float(mhfr))  ! destaggering
        enddo
        call interp(zalf, ze1, ze2, nlev, zauxp, zaux, zalp, zu(jlon,jlat,1:nlevpo), nlevpo)
      endif

!  Interpolation of v (averaged on t points)

      if (jlat < nlat) then
        zauxp = zlnp
        do jklev = 1, nlev
        zaux(nlevp1-jklev) = (v(jlon,jlat+1,jklev)+(2.*float(mhfr)-1.)*v(jlon,jlat,jklev))/(2.*float(mhfr))  ! destaggering
        enddo
        call interp (zalf, ze1, ze2, nlev, zauxp, zaux, zalp, zv(jlon,jlat,1:nlevpo), nlevpo)
      endif

!  Interpolation of w (averaged vertically)

      ze1 = .5
      ze2 = .5
      zauxp(1:nlev) = zlnph(1:nlev)
      do jklev = 1, nlev
      zaux(nlevp1-jklev) = w(jlon,jlat,jklev+1)
      enddo
      call interp (zalf, ze1, ze2, nlev, zauxp, zaux, zalp, zw(jlon,jlat,1:nlevpo), nlevpo)

!  Interpolation of potvor

      if (jlon.ne.1.and.jlat.ne.nlat) then
        ze1 = .5
        ze2 = .5
        zauxp(1:nlev) = zlnph(1:nlev)
        do jklev = 1, nlev
          zaux(nlevp1-jklev) = .25*(potvor(jlon,jlat  ,jklev+1)+potvor(jlon-1,jlat  ,jklev+1) + &
                                potvor(jlon,jlat+1,jklev+1)+potvor(jlon-1,jlat+1,jklev+1))*1.e6 ! conversion to PV units
        enddo
        call interp (zalf, ze1, ze2, nlev, zauxp, zaux, zalp, zpv(jlon,jlat,1:nlevpo), nlevpo)
      endif

      do jklev = 1, nlevpo
      zq  (jlon,jlat,jklev)=max(zq  (jlon,jlat,jklev), 1.e-9)
      zqcw(jlon,jlat,jklev)=max(zqcw(jlon,jlat,jklev), 0.)
      zqci(jlon,jlat,jklev)=max(zqci(jlon,jlat,jklev), 0.)
      if(zalp(jklev).gt.alog(ps(jlon,jlat))) then
      zq  (jlon,jlat,jklev)=0.
      zqcw(jlon,jlat,jklev)=0.
      zqci(jlon,jlat,jklev)=0.
! Wind under orography
      zu  (jlon,jlat,jklev)=0. !!u(jlon,jlat,1)
      zv  (jlon,jlat,jklev)=0. !!v(jlon,jlat,1)
      zw  (jlon,jlat,jklev)=0.
      zpv (jlon,jlat,jklev) = 0.
      endif
      enddo

      enddo
      enddo

!---------------------------------------------------------------------------------------------------------------
!  Comput. of total cloud water and ice
!  Defin. of relative hum. and thetae on pressure levels
!  Defin. of specific hum. below ground
!---------------------------------------------------------------------------------------------------------------

      zclwi = zqcw + zqci

      do jlat=1,nlat
      do jlon=1,nlon

!  First level above ground

      jlev1=nlevpo
      do jklev = 1, nlevpo
        if(plevo(jklev).lt.ps(jlon,jlat)) jlev1=min(jlev1,jklev)
      enddo

!  Above ground only

      do jklev=jlev1,nlevpo
      zf = min(0.1, zq(jlon,jlat,jklev))
      ztdan = zt(jlon,jlat,jklev)+tzer

      call comp_esk(zzpvs,zfs,ztdan,plevo(jklev),1) ! Computes saturation to water and ice separately
      zzpv=(zf*plevo(jklev))/(eps + zf -eps*zf)
      if(jklev.eq.jlev1) zrelh1=min(zzpv/zzpvs*100., 100.) ! Used to extrapolate below ground
      zrh(jlon,jlat,jklev) = min(zzpv/zzpvs*100., 102.)

!  Equivalent potential temperature (Bolton, MWR, 1980)
!  ztl: temp. at condens. lev. (computed with expr. by Bolton)

      zmr = zf/(1.-zf)
      ze = plevo(jklev)/100. * zmr/(eps+zmr)
      ztl = 55. + 2840./(3.5*alog(ztdan) - alog(ze) - 4.805)
      zespon = (3.376/ztl - 0.00254)*zmr*1000*(1.+.81*zmr)
      zthetae(jlon,jlat,jklev) = ztdan*(1.e5/plevo(jklev))**rdrcp*exp(zespon)
      enddo

!  Definitions at levels below ground (some quantities are then smoothed)

      do jklev=1,jlev1-1
      zq(jlon,jlat,jklev) = zq(jlon,jlat,jlev1)
      zrh(jlon,jlat,jklev) = zrelh1
      zthetae(jlon,jlat,jklev)=zthetae(jlon,jlat,jlev1)
      enddo

      enddo
      enddo

!---------------------------------------------------------------------------------------------------------------
!  Filtering of temp., geop., thetae, relh and q at pressure levels below ground
!  Lateral boundary points are not filtered, so extrap. is applied later
!  Only vertical levels below 400e2 Pa are considered
!---------------------------------------------------------------------------------------------------------------

      nlevmx = 1
      do jklev = 1,nlevpo
      if(plevo(jklev).gt.400.e2) nlevmx = jklev
      enddo

      do jklev = 1, nlevmx
      do iter = 1,18

      do 150 jlat = 2,nlatm1
      do 150 jlon = 2,nlonm1
      if(plevo(jklev).gt.p(jlon,jlat,1)) then
      zwork(jlon,jlat) = .4*zt(jlon,jlat,jklev)  + .6/4.*(zt(jlon-1,jlat,jklev)+zt(jlon+1,jlat,jklev)+   &
                            zt(jlon,jlat-1,jklev)+zt(jlon,jlat+1,jklev))
      zwor2(jlon,jlat) = .4*zph(jlon,jlat,jklev) + .6/4.*(zph(jlon-1,jlat,jklev)+zph(jlon+1,jlat,jklev)+ &
                            zph(jlon,jlat-1,jklev)+zph(jlon,jlat+1,jklev))
      zwor3(jlon,jlat) = .4*zthetae(jlon,jlat,jklev) + .6/4.*                                            &
                           (zthetae(jlon-1,jlat,jklev)+zthetae(jlon+1,jlat,jklev)+                       &
                            zthetae(jlon,jlat-1,jklev)+zthetae(jlon,jlat+1,jklev))
      zwor4(jlon,jlat) = .4*zrh(jlon,jlat,jklev) + .6/4.*(zrh(jlon-1,jlat,jklev)+zrh(jlon+1,jlat,jklev)+ &
                            zrh(jlon,jlat-1,jklev)+zrh(jlon,jlat+1,jklev))
      zwor5(jlon,jlat) = .4*zq(jlon,jlat,jklev)  + .6/4.*(zq(jlon-1,jlat,jklev)+zq(jlon+1,jlat,jklev)+   &
                            zq(jlon,jlat-1,jklev)+zq(jlon,jlat+1,jklev))
      else
      zwork(jlon,jlat) = zt (jlon,jlat,jklev)
      zwor2(jlon,jlat) = zph(jlon,jlat,jklev)
      zwor3(jlon,jlat) = zthetae(jlon,jlat,jklev)
      zwor4(jlon,jlat) = zrh(jlon,jlat,jklev)
      zwor5(jlon,jlat) = zq(jlon,jlat,jklev)
      endif
 150  continue

      do 152 jlat = 2,nlatm1
      do 152 jlon = 2,nlonm1
      zt (jlon,jlat,jklev)     = zwork(jlon,jlat)
      zph(jlon,jlat,jklev)     = zwor2(jlon,jlat)
      zthetae(jlon,jlat,jklev) = zwor3(jlon,jlat)
      zrh(jlon,jlat,jklev)     = zwor4(jlon,jlat)
      zq(jlon,jlat,jklev)      = zwor5(jlon,jlat)
 152  continue

      enddo
      enddo

!  Extrapolation of some surface fields at lateral boundaries

      zeex = .5
      call bextr(hflux ,nlon,nlat,zeex,.false.)
      call bextr(qflux ,nlon,nlat,zeex,.false.)
      call bextr(prectot,nlon,nlat,zeex,.false.)
      call bextr(precconv,nlon,nlat,zeex,.false.)
      call bextr(precsolid,nlon,nlat,zeex,.false.)
      call bextr(snow,  nlon,nlat,zeex,.false.)
      call bextr(runoff,nlon,nlat,zeex,.false.)
      call bextr(runoff_tot,nlon,nlat,zeex,.false.)
      call bextr(cloudt,nlon,nlat,zeex,.false.)
      call bextr(cloudh,nlon,nlat,zeex,.false.)
      call bextr(cloudm,nlon,nlat,zeex,.false.)
      call bextr(cloudl,nlon,nlat,zeex,.false.)
      call bextr(t2,    nlon,nlat,zeex,.false.)
      call bextr(td2,   nlon,nlat,zeex,.false.)
      call bextr(q2,    nlon,nlat,zeex,.false.)
      call bextr(q2rel, nlon,nlat,zeex,.false.)
      call bextr(ztskin,nlon,nlat,zeex,.false.)
      call bextr(ztgr1, nlon,nlat,zeex,.false.)
      call bextr(ztgr2, nlon,nlat,zeex,.false.)
      call bextr(ztgr3, nlon,nlat,zeex,.false.)
      call bextr(ztgr4, nlon,nlat,zeex,.false.)
      call bextr(zqgr1, nlon,nlat,zeex,.false.)
      call bextr(zqgr2, nlon,nlat,zeex,.false.)
      call bextr(zqgr3, nlon,nlat,zeex,.false.)
      call bextr(zqgr4, nlon,nlat,zeex,.false.)
      call bextr(qskin, nlon,nlat,zeex,.false.)
      call bextr(t05,   nlon,nlat,zeex,.false.)
      call bextr(q05rel,nlon,nlat,zeex,.false.)
      call bextr(t005,  nlon,nlat,zeex,.false.)
      call bextr(tg005, nlon,nlat,zeex,.false.)
      call bextr(tg010, nlon,nlat,zeex,.false.)
      call bextr(tg020, nlon,nlat,zeex,.false.)
      call bextr(tg050, nlon,nlat,zeex,.false.)
      call bextr(tg100, nlon,nlat,zeex,.false.)
      call bextr(cswfl, nlon,nlat,zeex,.false.)
      call bextr(clwfl, nlon,nlat,zeex,.false.)
      call bextr(chflux,nlon,nlat,zeex,.false.)
      call bextr(cqflux,nlon,nlat,zeex,.false.)
      call bextr(t2min, nlon,nlat,zeex,.false.)
      call bextr(t2max, nlon,nlat,zeex,.false.)
      call bextr(ws10max,nlon,nlat,zeex,.false.)
      call bextr(pbl,   nlon,nlat,zeex,.false.)

      do k=1,nlevg
        call bextr(tg_post(:,:,k), nlon,nlat,zeex,.false.)
      enddo

!  Resetting of unphysical values (due only to extrapolation at boundaries)

      do jlat = 1, nlat
      do jlon = 1, nlon
      prectot(jlon,jlat) = max(prectot(jlon,jlat), 0.)
      precconv(jlon,jlat) = max(precconv(jlon,jlat), 0.)
      precsolid(jlon,jlat) = max(precsolid(jlon,jlat), 0.)
      snow  (jlon,jlat) = max(snow  (jlon,jlat), 0.)
      runoff(jlon,jlat) = max(runoff(jlon,jlat), 0.)
      runoff_tot(jlon,jlat) = max(runoff_tot(jlon,jlat), 0.)
      q2    (jlon,jlat) = max(q2    (jlon,jlat), 0.)
      q2rel (jlon,jlat) = max(q2rel (jlon,jlat), 0.)
      q2rel (jlon,jlat) = min(q2rel (jlon,jlat), 102.)
      q05rel(jlon,jlat) = max(q05rel(jlon,jlat), 0.)
      q05rel(jlon,jlat) = min(q05rel(jlon,jlat), 102.)
      cloudt(jlon,jlat) = max(cloudt(jlon,jlat), 0.)
      cloudt(jlon,jlat) = min(cloudt(jlon,jlat), 1.)
      cloudl(jlon,jlat) = max(cloudl(jlon,jlat), 0.)
      cloudl(jlon,jlat) = min(cloudl(jlon,jlat), 1.)
      cloudm(jlon,jlat) = max(cloudm(jlon,jlat), 0.)
      cloudm(jlon,jlat) = min(cloudm(jlon,jlat), 1.)
      cloudh(jlon,jlat) = max(cloudh(jlon,jlat), 0.)
      cloudh(jlon,jlat) = min(cloudh(jlon,jlat), 1.)
      enddo
      enddo

! Accumulated fluxes in kjoule/m2

      cswfl(:,:)=cswfl(:,:)*1.e-3
      clwfl(:,:)=clwfl(:,:)*1.e-3
      chflux(:,:)=chflux(:,:)*1.e-3
      cqflux(:,:)=cqflux(:,:)*1.e-3

! Extreme temperatures in C degrees

      t2min(:,:)=t2min(:,:)-tzer
      t2max(:,:)=t2max(:,:)-tzer

!---------------------------------------------------------------------------------------------------------------
! Calculation and writing of forecast data

 idate0(1:5)=nfdr(5:9)
 iperiod(1:3)=nfdr(10:12)

 iforecast_h=iperiod(1)*24+iperiod(2)
 zaccum_min=nint((float(nhist*njump)*dtstep)/60.)
 iperiod_accum(1)=int(zaccum_min/60./24.)
 iperiod_accum(2)=int((zaccum_min-float(iperiod_accum(1))*60.*24.)/60.)
 iperiod_accum(3)=int(zaccum_min-float(iperiod_accum(1))*60.*24.-float(iperiod_accum(2))*60.)

! Date and time for current forecast instant

 idatec(:)=idate0(:)
 idatec(3)=idatec(3)+iperiod(1)
 idatec(4)=idatec(4)+iperiod(2)
 idatec(5)=idatec(5)+iperiod(3)
 if (mod(idatec(1),4)==0) then
   imon(2)=29
 else
   imon(2)=28
 endif
 nday=iperiod(1)
 ndm=0
 do iimon=1,50
   nday=nday-imon(idatec(2))
   if (nday<0) then
     exit
   else
     ndm=ndm+imon(idatec(2))
     idatec(2)=idatec(2)+1
     if (idatec(2)>12) then
       idatec(2)=idatec(2)-12
       idatec(1)=idatec(1)+1
       if (mod(idatec(1),4)==0) then
         imon(2)=29
       else
         imon(2)=28
       endif
     endif
   endif
 enddo
 idatec(3)=idatec(3)-ndm
 if (idatec(5)>=60) then
   idatec(5)=idatec(5)-60
   idatec(4)=idatec(4)+1
 endif
 if (idatec(4)>=24) then
   idatec(4)=idatec(4)-24
   idatec(3)=idatec(3)+1
 endif
 if (idatec(3)>imon(idatec(2))) then
   idatec(3)=idatec(3)-imon(idatec(2))
   idatec(2)=idatec(2)+1
   if (idatec(2)>12) then
     idatec(2)=idatec(2)-12
     idatec(1)=idatec(1)+1
   endif
 endif

!=======================================================================

! Lapse rate

   day=0.
   do i=1,idatec(2)-1
     day=day+float(imon(i))
   enddo
   day=day+float(idatec(3))+float(idatec(4))/24.

   zlat=alatt(nlon/2,nlat/2)
   if(zlat > 30.) then
   coeday = 0.5*(1.-cos(day*2.*pi/365.))   ! N.H.
   elseif(zlat < -30.) then
   coeday = 0.5*(1.+cos(day*2.*pi/365.))   ! S.H.
   else
   coeday = 0.5                            ! tropics
   endif

! lapse_rate_clim is a climatological prescribed lapse rate, changing between 4.5e-3 (end of dec.) and 7.4e-3 (end of june)

   lapse_rate_clim = 4.5e-3 + 2.8e-3 * coeday

! computation of average lapse rate of input model between lowest levels

   lapse_rate=0.
   nnn=0
   do jlat=1,nlat
   do jlon=1,nlon
     if (fmask(jlon,jlat) < 0.5) then
      lapse_rate=lapse_rate+(t(jlon,jlat,5)-t(jlon,jlat,1))/(zeta(jlon,jlat,5)-zeta(jlon,jlat,1))
      nnn=nnn+1
     endif
   enddo
   enddo
   lapse_rate=lapse_rate/float(nnn)
   lapse_rate=min (0.70*lapse_rate - 0.30*lapse_rate_clim, 0.) ! limited to avoid inversions

   lapse_rate_1(:,:)=lapse_rate

   do jlat=1,nlat
   do jlon=1,nlon
     jklev1=1
     jklev2=2
     do jklev=1,nlev
       if (zeta(jlon,jlat,jklev)-zorogr(jlon,jlat) >= 200.) exit
     enddo
     jklev1=jklev
     do jklev=nlev,2,-1
       if (zeta(jlon,jlat,jklev)-zorogr(jlon,jlat) <= 1500.) exit
     enddo
     jklev2=max(jklev, jklev1+1)
     lapse_rate=(t(jlon,jlat,jklev2)-t(jlon,jlat,jklev1))/(zeta(jlon,jlat,jklev2)-zeta(jlon,jlat,jklev1))
     lapse_rate_2(jlon,jlat)=lapse_rate
   enddo
   enddo

!=======================================================================

 do jklev = 1, nlevpo
   isobarlev = int(plevo(jklev)/100.)
   if (isobarlev==850) lev850=jklev
   if (isobarlev==700) lev700=jklev
 enddo

! ---- For RTTOV simulation ---

#ifdef rttov

 do jklev=1,nlev
   jklev1=nlev-jklev+1
   p_invert(:,:,jklev1)=p(:,:,jklev)
   t_invert(:,:,jklev1)=t(:,:,jklev)
   q_invert(:,:,jklev1)=q(:,:,jklev)
   qcw_invert(:,:,jklev1)=qcw(:,:,jklev)
   qci_invert(:,:,jklev1)=qci(:,:,jklev)
 enddo

 fseaice = fice

 do jlat=1,nlat
 do jlon=1,nlon
   flag_cloud_conv(jlon,jlat) = 0
   do jklev=6,nlev
     if (w(jlon,jlat,jklev) > w_crit_conv) flag_cloud_conv(jlon,jlat) = 1
   enddo
 enddo
 enddo

 date_time(1:5) = idatec(1:5)
 iini=max(nint(nlon*zoom_xini), 1)
 ifin=min(nint(nlon*zoom_xfin), nlon)
 jini=max(nint(nlat*zoom_yini), 1)
 jfin=min(nint(nlat*zoom_yfin), nlat)
 npoint_jump = max( (ifin-iini+1)*(jfin-jini+1)/300000, 1)
 alon1=pdr(5)+pdr(2)*(iini-1)
 alat1=pdr(4)+pdr(1)*0.5/float(mhfr)+pdr(1)*(jini-1)

 print *
 print *,"------------ Call RTTOV_fwd model ----------------"

 call rttov11_fwd_interface &
 (nchan, rttov_sat_series, rttov_sat_id, rttov_sat_sensor, channel_elabor_list(1:nchan), flag_visible, &
 path_rttov_emis_atlas,  path_rttov_brdf_atlas, &
 ifin-iini+1, jfin-jini+1, npoint_jump, nlev, &
 alatt(iini:ifin,jini:jfin), alont(iini:ifin,jini:jfin), fmask(iini:ifin,jini:jfin), &
 phig(iini:ifin,jini:jfin)/g, emis1(iini:ifin,jini:jfin), &
 ps(iini:ifin,jini:jfin), t2(iini:ifin,jini:jfin)+tzer, q2(iini:ifin,jini:jfin), &
 u10(iini:ifin,jini:jfin), v10(iini:ifin,jini:jfin), tskin(iini:ifin,jini:jfin), &
 fsnow(iini:ifin,jini:jfin), qg(iini:ifin,jini:jfin,1), &
 fseaice(iini:ifin,jini:jfin), flag_cloud_conv(iini:ifin,jini:jfin), &
 p_invert(iini:ifin,jini:jfin,:), t_invert(iini:ifin,jini:jfin,:), q_invert(iini:ifin,jini:jfin,:), &
 qcw_invert(iini:ifin,jini:jfin,:), qci_invert(iini:ifin,jini:jfin,:), clfrac_invert(iini:ifin,jini:jfin,:), &
 date_time, sensor_chan_id, sensor_chan_cw, &
 radiance(iini:ifin,jini:jfin,:), radiance_bt(iini:ifin,jini:jfin,:), radiance_refl(iini:ifin,jini:jfin,:), &
 radiance_clear(iini:ifin,jini:jfin,:), radiance_clear_bt(iini:ifin,jini:jfin,:), radiance_clear_refl(iini:ifin,jini:jfin,:), &
 emis_rttov(iini:ifin,jini:jfin,:))

 print *
 print *,"------------ RTTOV-11 forward model simulated data for ----------------"
 print *,"--- Satellite code ",rttov_sat_series," number ",rttov_sat_id," instrument (sensor) code ",rttov_sat_sensor," ----"
 print *,"---- Elaborated Channel list ----"
 print *,sensor_chan_id(:)
 print *,"---- Elaborated Channel Central Wavenumber ----"
 print *,nint(sensor_chan_cw(:))
 write (*,"(4(a,i4),a,i1,a)") &
 "---- fiels defined at subdomain (",iini,":",ifin,",",jini,":",jfin,") with point jump ",npoint_jump," ----"

 call write_grib2_horizontal_grid_rttov_data (2, ifin-iini+1, jfin-jini+1, npoint_jump, &
 pdr(39), pdr(38), alon1, alat1, pdr(2), pdr(1),                                        &
 idate0, iperiod, iperiod_accum, idatec,                                                &
 nchan, rttov_sat_series, rttov_sat_id, rttov_sat_sensor,                               &
 sensor_chan_id, sensor_chan_cw,                                                        &
 radiance(iini:ifin,jini:jfin,:), radiance_bt(iini:ifin,jini:jfin,:),                   &
 radiance_clear(iini:ifin,jini:jfin,:), radiance_clear_bt(iini:ifin,jini:jfin,:),       &
 emis_rttov(iini:ifin,jini:jfin,:))

 print *
 print *,"------------------ RTTOV-11 products created ---------------------"
 print *
#endif

!---------------------------------------------------------------------------------------------------------------

 if (ncrossec.and.cross_number > 0) then

! Free cross section definition

   call free_cross_section(nlon, nlat, nlev, nlivz, cross_number, npoint_cross_max, npoint_cross,   &
 lonini_cross, latini_cross, lonfin_cross, latfin_cross, zlcr,                                      &
 x0, y0, dlon, dlat, alon0, alat0+dlat*0.5, dz, h, a0, b0, fz, fzh, rdrcp, valmiss, cross_interp_gauss, &
 alont, alatt, alonu, alatu, alonv, alatv,                                                          &
 zorogr, p, t, zthee, u, v, w, zrelh, qcw, qci,                                                     &
 orog_cross, t_cross, theta_cross, thetae_cross, u_cross, v_cross, u_tang_cross, v_norm_cross,      &
 w_cross, rh_cross, qcw_cross, qci_cross)

 endif

!---------------------------------------------------------------------------------------------------------------
!  Output definitions and write
!  fmask, longitude and latitude of rotated grid
!---------------------------------------------------------------------------------------------------------------

if (output_format_ppf) then

      if(iist2.eq.1) then
      open(80,file='moloch.ppf',status='unknown',form='formatted')
      write(80,1111) isflag(1:80)
      write(80,1112) ipflag
!      write(80,'(i1)') izflag
 1111 format(80i1)
 1112 format(8i1)
      call wrpost(fmask, izout, nlon, nlat, nfdr, pdr)
      call wrpost(alont, izout, nlon, nlat, nfdr, pdr)
      call wrpost(alatt, izout, nlon, nlat, nfdr, pdr)
      endif

!  orography (in m)

      if(isflag(1).eq.1) call wrpost(zorogr, izout, nlon, nlat, nfdr, pdr)

!  pz0 (mslp, in hPa)

      if(isflag(2).eq.1) call wrpost(zpz0,   izout, nlon, nlat, nfdr, pdr)

!  prectot (tot. accum. precip. in kg/m2 eq. to mm)

      if(isflag(3).eq.1) call wrpost(prectot, izout, nlon, nlat, nfdr, pdr)

!  precsolid (tot. accum. snowfall in kg/m2 eq. to mm of equivalent water)

      if(isflag(4).eq.1) call wrpost(precsolid, izout, nlon, nlat, nfdr, pdr)

!  u10 (wind comp. at 10 m)

      if(isflag(5).eq.1) call wrpost(u10,    izout, nlon, nlat, nfdr, pdr)

!  v10 (wind comp. at 10 m)

      if(isflag(5).eq.1) call wrpost(v10,    izout, nlon, nlat, nfdr, pdr)

!  t2 (air temperature at 2 m)

      if(isflag(6).eq.1) call wrpost(t2,     izout, nlon, nlat, nfdr, pdr)

!  q2 (specific humidity at 2 m)

      if(isflag(7).eq.1) call wrpost(q2,     izout, nlon, nlat, nfdr, pdr)

!  q2rel (relative humidity at 2 m)

      if(isflag(8).eq.1) call wrpost(q2rel,  izout, nlon, nlat, nfdr, pdr)

!  zqgr1 (ground water lev. 1, in meters)

      if(isflag(9).eq.1) call wrpost(zqgr1,  izout, nlon, nlat, nfdr, pdr)

!  zqgr2 (ground water lev. 2, in meters)

      if(isflag(10).eq.1) call wrpost(zqgr2, izout, nlon, nlat, nfdr, pdr)

!  zqgr3 (ground water lev. 3, in meters)

      if(isflag(11).eq.1) call wrpost(zqgr3,  izout, nlon, nlat, nfdr, pdr)

!  zqgr4 (ground water lev. 4, in meters)

      if(isflag(12).eq.1) call wrpost(zqgr4, izout, nlon, nlat, nfdr, pdr)

!  ztskin (skin temperature)

      if(isflag(13).eq.1) call wrpost(ztskin,izout, nlon, nlat, nfdr, pdr)

!  ztgr1 (ground temp. lev 1)

      if(isflag(14).eq.1) call wrpost(ztgr1, izout, nlon, nlat, nfdr, pdr)

!  ztgr2 (ground temp. lev 2)

      if(isflag(15).eq.1) call wrpost(ztgr2, izout, nlon, nlat, nfdr, pdr)

!  ztgr3 (ground temp. lev 3)

      if(isflag(16).eq.1) call wrpost(ztgr3, izout, nlon, nlat, nfdr, pdr)

!  ztgr4 (ground temp. lev 4)

      if(isflag(17).eq.1) call wrpost(ztgr4, izout, nlon, nlat, nfdr, pdr)

!  hfl (upward flux of sens. heat in watt/m**2)

      if(isflag(18).eq.1) call wrpost(hflux, izout, nlon, nlat, nfdr, pdr)

!  qfl (upward flux of latent heat in watt/m**2)

      if(isflag(19).eq.1) call wrpost(qflux, izout, nlon, nlat, nfdr, pdr)

!  temperature at lev. 1

      if(isflag(20).eq.1) call wrpost(zts,   izout, nlon, nlat, nfdr, pdr)

!  equiv. pot. temp. at lev.1

      if(isflag(21).eq.1) call wrpost(zthes, izout, nlon, nlat, nfdr, pdr)

!  us (wind vel. at lev.1)

      if(isflag(22).eq.1) call wrpost(zus,   izout, nlon, nlat, nfdr, pdr)

!  vs (wind vel. at lev.1)

      if(isflag(22).eq.1) call wrpost(zvs,   izout, nlon, nlat, nfdr, pdr)

!  zqs (spec. humid. at lev.1)

      if(isflag(23).eq.1) call wrpost(zqs,   izout, nlon, nlat, nfdr, pdr)

!  total, high, mid and low cloud cover

      if(isflag(24).eq.1) call wrpost(cloudt,izout, nlon, nlat, nfdr, pdr)
      if(isflag(25).eq.1) call wrpost(cloudh,izout, nlon, nlat, nfdr, pdr)
      if(isflag(26).eq.1) call wrpost(cloudm,izout, nlon, nlat, nfdr, pdr)
      if(isflag(27).eq.1) call wrpost(cloudl,izout, nlon, nlat, nfdr, pdr)

!  snow height (m of equivalent water)

      if(isflag(28).eq.1) call wrpost(snow,  izout, nlon, nlat, nfdr, pdr)

!  runoff (tot. accum. run-off in kg/m2 eq. to mm/m2)

      if(isflag(29).eq.1) call wrpost(runoff, izout, nlon, nlat, nfdr, pdr)

!  cape (j/kg)

      if(isflag(30).eq.1) call wrpost(cape,  izout, nlon, nlat, nfdr, pdr)

! qskin (land surface specific humidity)

      if(isflag(31).eq.1) call wrpost(qskin,  izout, nlon, nlat, nfdr, pdr)

! t05 (air temperature at 0.5 m)

      if(isflag(32).eq.1) call wrpost(t05,  izout, nlon, nlat, nfdr, pdr)

! q05rel (air relative humidity at 0.5 m, %)

      if(isflag(33).eq.1) call wrpost(q05rel,  izout, nlon, nlat, nfdr, pdr)

! t005 (air temperature at 0.05 m)

      if(isflag(34).eq.1) call wrpost(t005,  izout, nlon, nlat, nfdr, pdr)

! tg005 (ground temperature at 5 cm depth)

      if(isflag(35).eq.1) call wrpost(tg005,  izout, nlon, nlat, nfdr, pdr)

! tg010 (ground temperature at 10 cm depth)

      if(isflag(36).eq.1) call wrpost(tg010,  izout, nlon, nlat, nfdr, pdr)

! tg020 (ground temperature at 20 cm depth)

      if(isflag(37).eq.1) call wrpost(tg020,  izout, nlon, nlat, nfdr, pdr)

! tg050 (ground temperature at 50 cm depth)

      if(isflag(38).eq.1) call wrpost(tg050,  izout, nlon, nlat, nfdr, pdr)

! tg100 (ground temperature at 100 cm depth)

      if(isflag(39).eq.1) call wrpost(tg100,  izout, nlon, nlat, nfdr, pdr)

! albedo

      if(isflag(40).eq.1) call wrpost(albedo,  izout, nlon, nlat, nfdr, pdr)

! emissivity in broadband window and in 8-12 micron window

      if(isflag(41).eq.1) call wrpost(emis1,  izout, nlon, nlat, nfdr, pdr)
      if(isflag(42).eq.1) call wrpost(emis2,  izout, nlon, nlat, nfdr, pdr)

! accumulated surface fluxes (upward flux of sens. heat in watt/m**2)

! short wave radiation

      if(isflag(43).eq.1) call wrpost(cswfl,  izout, nlon, nlat, nfdr, pdr)

! long wave radiation

      if(isflag(44).eq.1) call wrpost(clwfl,  izout, nlon, nlat, nfdr, pdr)

! sensible heat

      if(isflag(45).eq.1) call wrpost(chflux, izout, nlon, nlat, nfdr, pdr)

! latent heat

      if(isflag(46).eq.1) call wrpost(cqflux, izout, nlon, nlat, nfdr, pdr)

! extreme values at standard observation levels

! t 2m min

      if(isflag(47).eq.1) call wrpost(t2min,  izout, nlon, nlat, nfdr, pdr)

! t 2m max

      if(isflag(48).eq.1) call wrpost(t2max,  izout, nlon, nlat, nfdr, pdr)

! wind speed 10m max

      if(isflag(49).eq.1) call wrpost(ws10max,izout, nlon, nlat, nfdr, pdr)

! lifted index (K)

      if(isflag(50).eq.1) call wrpost(zlift,  izout, nlon, nlat, nfdr, pdr)

! sea ice fraction

      if(isflag(51).eq.1) call wrpost(fice,  izout, nlon, nlat, nfdr, pdr)

! sea ice thickness (m)

      if(isflag(52).eq.1) call wrpost(iceth, izout, nlon, nlat, nfdr, pdr)

! dew point temperature at 2 m (C)

      if(isflag(53).eq.1) call wrpost(td2, izout, nlon, nlat, nfdr, pdr)

! integrated water vapour (kg/m2)

      if(isflag(54).eq.1) call wrpost(iwv, izout, nlon, nlat, nfdr, pdr)

! CIN

      if(isflag(55).eq.1) call wrpost(cin, izout, nlon, nlat, nfdr, pdr)

! PBL top height

      if(isflag(56).eq.1) call wrpost(pbl, izout, nlon, nlat, nfdr, pdr)

! gust

      if(isflag(57).eq.1) call wrpost(gust, izout, nlon, nlat, nfdr, pdr)

!  loop for fields defined on pressure levels

      do jklev = 1, nlevpo
      isobarlev = nint(plevo(jklev)/100.)

!  smoothing of geopotential

      zrd2 = dx**2/dy**2
      zw0  = -2.-2.*zrd2
      do jlat = 2, nlatm1
      do jlon = 2, nlonm1
      zwork(jlon,jlat) = zph(jlon,jlat,jklev) +0.15*                                      &
         ( (zph(jlon,jlat-1,jklev)+zph(jlon,jlat+1,jklev))*zrd2 + zph(jlon-1,jlat,jklev)+ &
         zph(jlon,jlat,jklev)*zw0 + zph(jlon+1,jlat,jklev) )
      enddo
      enddo
      do jlat = 2, nlatm1
      do jlon = 2, nlonm1
      zph(jlon,jlat,jklev) = zwork(jlon,jlat)
      enddo
      enddo

!  geopotential height

      if(ipflag(1, jklev).eq.1) call wrpost(zph    (1,1,jklev), izout, nlon, nlat, nfdr, pdr)

!  t (temperature, deg. Celsius)

      if(ipflag(2, jklev).eq.1) call wrpost(zt     (1,1,jklev), izout, nlon, nlat, nfdr, pdr)

!  u (wind comp., m/s)

      if(ipflag(3, jklev).eq.1) call wrpost(zu     (1,1,jklev), izout, nlon, nlat, nfdr, pdr)

!  v (wind comp., m/s)

      if(ipflag(3, jklev).eq.1) call wrpost(zv     (1,1,jklev), izout, nlon, nlat, nfdr, pdr)

!  q (specific humidity)

      if(ipflag(4, jklev).eq.1) call wrpost(zq     (1,1,jklev), izout, nlon, nlat, nfdr, pdr)

!  rh (relative humidity, % from 0 to 100)

      if(ipflag(5, jklev).eq.1) call wrpost(zrh    (1,1,jklev), izout, nlon, nlat, nfdr, pdr)

!  w (wind comp., m/s)

      if(ipflag(6, jklev).eq.1) call wrpost(zw     (1,1,jklev), izout, nlon, nlat, nfdr, pdr)

!  equivalent potential temperature

      if(ipflag(7, jklev).eq.1) call wrpost(zthetae(1,1,jklev), izout, nlon, nlat, nfdr, pdr)

!  total cloud water + ice (kg/kg) or Potential porticity

      if (ipflag(8, jklev).eq.1) then
        if (plevo(jklev).le.300.e2) then ! PV
          call wrpost(zpv    (1,1,jklev), izout, nlon, nlat, nfdr, pdr)
        else ! total cloud water + ice
          call wrpost(zclwi  (1,1,jklev), izout, nlon, nlat, nfdr, pdr)
        endif
      endif

      enddo ! jklev

!  loop for fields defined on z(m) levels

!      do jklev=1,nzlev
!      if(izflag(1, jklev).eq.1) call wrpost(<zfield>(1,1,jklev), izout, nlon, nlat, nfdr, pdr)
!      enddo

!---------------------------------------------------------------------------------------------------------------
!  North-south cross-sections
!---------------------------------------------------------------------------------------------------------------

      zalf= .5
      ze1 = .7
      ze2 = .7

      if (ncrossec) then

      if(iist2.eq.1) then
      open(81,file='moloch_crossy.ppf',status='unknown',form='formatted')
      open(82,file='moloch_crossx.ppf',status='unknown',form='formatted')
      endif

      do 430 j=1,ncrossy
      jlon=loncr(j)
      jlonm1=max(jlon-1,1)

      do 420 jlat=1,nlat
      jlatp1=min(jlat+1,nlat)

!  zlnp contains hereafter the height of moloch levels

      do jklev= 1, nlev
      zita = (jklev-1)*dz + .5*dz
      zlnp(jklev) = - h*bzita(zita,h,b0)*log(fz(jklev)) + zorogr(jlon,jlat)*gzita(zita,h,a0)
      enddo

!  theta y-sections

      do jklev= 1, nlev
      zaux(jklev) = t(jlon,jlat,jklev)*(1.e5/p(jlon,jlat,jklev))**rdrcp
      enddo
      call interp(zalf, ze1, ze2, nlev, zlnp, zaux, zlcr, ztcry(jlat,1:nlivz,j), nlivz)

!  thetae y-sections

      do jklev= 1, nlev
      zaux(jklev)=zthee(jlon,jlat,jklev)
      enddo
      call interp(zalf, ze1, ze2, nlev, zlnp, zaux, zlcr, ztecry(jlat,1:nlivz,j), nlivz)

!  u-wind y-sections

      do jklev= 1, nlev
      zaux(jklev)=.5*(u(jlon,jlat,jklev)+u(jlonm1,jlat,jklev))
      enddo
      call interp(zalf, ze1, ze2, nlev, zlnp, zaux, zlcr, zucry(jlat,1:nlivz,j), nlivz)

!  v-wind y-sections

      do jklev= 1, nlev
      zaux(jklev)=.5*(v(jlon,jlat,jklev)+v(jlon,jlatp1,jklev))
      enddo
      call interp(zalf, ze1, ze2, nlev, zlnp, zaux, zlcr, zvcry(jlat,1:nlivz,j), nlivz)

!  w  y-sections

      do jklev= 1, nlev
      zaux(jklev)=.5*(w(jlon,jlat,jklev)+w(jlon,jlat,jklev+1))
      enddo
      call interp(zalf, ze1, ze2, nlev, zlnp, zaux, zlcr, zwcry(jlat,1:nlivz,j), nlivz)

!  relh y-sections

      do jklev= 1, nlev
      zaux(jklev)=zrelh(jlon,jlat,jklev)
      enddo
      call interp(zalf, ze1, ze2, nlev, zlnp, zaux, zlcr, zrhcry(jlat,1:nlivz,j), nlivz)

!  cloud water y-sections

      do jklev= 1, nlev
      zaux(jklev)=qcw(jlon,jlat,jklev)
      enddo
      call interp(zalf, ze1, ze2, nlev, zlnp, zaux, zlcr, zcwcry(jlat,1:nlivz,j), nlivz)

!  cloud ice y-sections

      do jklev= 1, nlev
      zaux(jklev)=qci(jlon,jlat,jklev)
      enddo
      call interp(zalf, ze1, ze2, nlev, zlnp, zaux, zlcr, zcicry(jlat,1:nlivz,j), nlivz)

 420  continue

!  Reset of negative values of humidity and water due to interpolation

      do jklev= 1, nlivz
      do jlat= 1, nlat
      if(zrhcry(jlat,jklev,j).lt.0) zrhcry(jlat,jklev,j)=0.
      if(zcwcry(jlat,jklev,j).lt.0) zcwcry(jlat,jklev,j)=0.
      if(zcicry(jlat,jklev,j).lt.0) zcicry(jlat,jklev,j)=0.
      enddo
      enddo

      call wrposy( ztcry(1:nlat,1:nlivz,j), izouty, zorogr(jlon,1:nlat), zlcr, jlon, nlon, nlat, nlivz, nfdr, pdr)
      call wrposy( zucry(1:nlat,1:nlivz,j), izouty, zorogr(jlon,1:nlat), zlcr, jlon, nlon, nlat, nlivz, nfdr, pdr)
      call wrposy(ztecry(1:nlat,1:nlivz,j), izouty, zorogr(jlon,1:nlat), zlcr, jlon, nlon, nlat, nlivz, nfdr, pdr)
      call wrposy( zvcry(1:nlat,1:nlivz,j), izouty, zorogr(jlon,1:nlat), zlcr, jlon, nlon, nlat, nlivz, nfdr, pdr)
      call wrposy( zwcry(1:nlat,1:nlivz,j), izouty, zorogr(jlon,1:nlat), zlcr, jlon, nlon, nlat, nlivz, nfdr, pdr)
      call wrposy(zrhcry(1:nlat,1:nlivz,j), izouty, zorogr(jlon,1:nlat), zlcr, jlon, nlon, nlat, nlivz, nfdr, pdr)
      call wrposy(zcwcry(1:nlat,1:nlivz,j), izouty, zorogr(jlon,1:nlat), zlcr, jlon, nlon, nlat, nlivz, nfdr, pdr)
      call wrposy(zcicry(1:nlat,1:nlivz,j), izouty, zorogr(jlon,1:nlat), zlcr, jlon, nlon, nlat, nlivz, nfdr, pdr)

 430  continue

!---------------------------------------------------------------------------------------------------------------
!  Same as above but for west-east cross-sections
!---------------------------------------------------------------------------------------------------------------

      do 530 j = 1, ncrossx
      jlat=latcr(j)
      jlatp1=min(jlat+1,nlat)

      do 520 jlon = 1, nlon
      jlonm1=max(jlon-1,1)

      do jklev= 1, nlev
      zita = (jklev-1)*dz + .5*dz
      zlnp(jklev) = -h*bzita(zita,h,b0)*log(fz(jklev)) + zorogr(jlon,jlat)*gzita(zita,h,a0)
      enddo

!  theta x-sections

      do jklev= 1,nlev
      zaux(jklev)=t(jlon,jlat,jklev)*(1.e5/p(jlon,jlat,jklev))**rdrcp
      enddo
      call interp(zalf, ze1, ze2, nlev, zlnp, zaux, zlcr, ztcrx(jlon,1:nlivz,j), nlivz)

!  thetae x-sections

      do jklev= 1,nlev
      zaux(jklev)=zthee(jlon,jlat,jklev)
      enddo
      call interp(zalf, ze1, ze2, nlev, zlnp, zaux, zlcr, ztecrx(jlon,1:nlivz,j), nlivz)

!  u-wind x-sections

      do jklev= 1,nlev
      zaux(jklev)=.5*(u(jlon,jlat,jklev)+u(jlonm1,jlat,jklev))
      enddo
      call interp(zalf, ze1, ze2, nlev, zlnp, zaux, zlcr, zucrx(jlon,1:nlivz,j), nlivz)

!  v-wind x-sections

      do jklev= 1,nlev
      zaux(jklev)=.5*(v(jlon,jlat,jklev)+v(jlon,jlatp1,jklev))
      enddo
      call interp(zalf, ze1, ze2, nlev, zlnp, zaux, zlcr, zvcrx(jlon,1:nlivz,j), nlivz)

!  w  x-sections

      do jklev= 1,nlev
      zaux(jklev)=.5*(w(jlon,jlat,jklev)+w(jlon,jlat,jklev+1))
      enddo
      call interp(zalf, ze1, ze2, nlev, zlnp, zaux, zlcr, zwcrx(jlon,1:nlivz,j), nlivz)

!  relh x-sections

      do jklev= 1, nlev
      zaux(jklev)=zrelh(jlon,jlat,jklev)
      enddo
      call interp(zalf, ze1, ze2, nlev, zlnp, zaux, zlcr, zrhcrx(jlon,1:nlivz,j), nlivz)

!  cloud water x-sections

      do jklev= 1, nlev
      zaux(jklev)=qcw(jlon,jlat,jklev)
      enddo
      call interp(zalf, ze1, ze2, nlev, zlnp, zaux, zlcr, zcwcrx(jlon,1:nlivz,j), nlivz)

!  cloud ice x-sections

      do jklev= 1, nlev
      zaux(jklev)=qci(jlon,jlat,jklev)
      enddo
      call interp(zalf, ze1, ze2, nlev, zlnp, zaux, zlcr, zcicrx(jlon,1:nlivz,j), nlivz)

 520  continue

!  Reset of negative values of humidity and water due to interpolation

      do jklev= 1, nlivz
      do jlon= 1, nlon
      if(zrhcrx(jlon,jklev,j).lt.0) zrhcrx(jlon,jklev,j)=0.
      if(zcwcrx(jlon,jklev,j).lt.0) zcwcrx(jlon,jklev,j)=0.
      if(zcicrx(jlon,jklev,j).lt.0) zcicrx(jlon,jklev,j)=0.
      enddo
      enddo

      call wrposx( ztcrx(1:nlon,1:nlivz,j), izoutx, zorogr(1,jlat), zlcr, jlat, nlon, nlat, nlivz, nfdr, pdr)
      call wrposx( zvcrx(1:nlon,1:nlivz,j), izoutx, zorogr(1,jlat), zlcr, jlat, nlon, nlat, nlivz, nfdr, pdr)
      call wrposx(ztecrx(1:nlon,1:nlivz,j), izoutx, zorogr(1,jlat), zlcr, jlat, nlon, nlat, nlivz, nfdr, pdr)
      call wrposx( zucrx(1:nlon,1:nlivz,j), izoutx, zorogr(1,jlat), zlcr, jlat, nlon, nlat, nlivz, nfdr, pdr)
      call wrposx( zwcrx(1:nlon,1:nlivz,j), izoutx, zorogr(1,jlat), zlcr, jlat, nlon, nlat, nlivz, nfdr, pdr)
      call wrposx(zrhcrx(1:nlon,1:nlivz,j), izoutx, zorogr(1,jlat), zlcr, jlat, nlon, nlat, nlivz, nfdr, pdr)
      call wrposx(zcwcrx(1:nlon,1:nlivz,j), izoutx, zorogr(1,jlat), zlcr, jlat, nlon, nlat, nlivz, nfdr, pdr)
      call wrposx(zcicrx(1:nlon,1:nlivz,j), izoutx, zorogr(1,jlat), zlcr, jlat, nlon, nlat, nlivz, nfdr, pdr)

 530  continue

      endif ! ncrossec

 endif ! output_format_ppf

 if (output_format_grib2) then

! Coding-writing output data in grib2 format

 call write_grib2_horizontal_grid (nlon,nlat,nlev,nlevg,nlev_snow_prod,     &
      mhfr,idate0,iperiod,iperiod_accum,idatec,pdr,plevo,                   &
      tzer,alont,alatt,fmask,zlev,                                          &
      zph,zt,zu,zv,zq,zrh,zw,zthetae,zclwi,zpv,                             &
      tg_post,qg_post,qgmax,qgmin,                                          &
      zorogr,zpz0,prectot,precsolid,                                        &
      u10,v10,t2,td2,q2,q2rel,ztskin,hflux,qflux,                           &
      zts,zthes,zus,zvs,zqs,cloudt,cloudh,cloudm,cloudl,snow,runoff,        &
      cape,qskin,albedo,emis1,emis2,cswfl,clwfl,chflux,cqflux,              &
      t2min,t2max,ws10max,zlift,fice,iceth,                                 &
      lapse_rate_1,lapse_rate_2, iwv, cin, pbl, gust, zlev0,                &
      snow_lev_depth_prod, snow_t_prod, snow_age_prod, snow_dens_prod,      &
      snowfall_height, snowfall_height0, ps)

! For constant (static) field writing (fmask, orography)

  if (iist2==1) then
#ifndef rttov
    call write_grib2_horizontal_grid_static_data(2, nlon, nlat,          &
         pdr(2), pdr(1), pdr(39), pdr(38), pdr(5), pdr(4)+pdr(1)*0.5, 1, &
         idate0, iperiod, iperiod_accum, idatec, zorogr, fmask, 0)
#else
    call write_grib2_horizontal_grid_static_data(2, nlon, nlat,          &
         pdr(2), pdr(1), pdr(39), pdr(38), pdr(5), pdr(4)+pdr(1)*0.5, npoint_jump, &
         idate0, iperiod, iperiod_accum, idatec, zorogr, fmask, 0)
#endif
  endif

    if (ncrossec) then ! Space cross-sections

      do j=1,ncrossy
        lon_crossy(j,1) = alont(loncr(j),1)
        lon_crossy(j,2) = alont(loncr(j),nlat)
        lat_crossy(j,1) = alatt(loncr(j),1)
        lat_crossy(j,2) = alatt(loncr(j),nlat)
        do jlat=1,nlat
          orog_crossy(jlat,j)=zorogr(loncr(j),jlat)
        enddo
      enddo

      do j=1,ncrossx
        lon_crossx(j,1) = alont(1,latcr(j))
        lon_crossx(j,2) = alont(nlon,latcr(j))
        lat_crossx(j,1) = alatt(1,latcr(j))
        lat_crossx(j,2) = alatt(nlon,latcr(j))
        do jlon=1,nlon
          orog_crossy(jlon,j)=zorogr(jlon,latcr(j))
        enddo
      enddo

#ifdef oper

!!      call write_grib2_cross_space_i_j (idate0, iperiod, iperiod_accum, idatec,             &
!! nlon, nlat, nlivz, ncrossy, ncrossx, zlcr, lon_crossy, lat_crossy, lon_crossx, lat_crossx, &
!! x0, y0, orog_crossy, orog_crossx,                                                          &
!! ztcry, ztecry, zvcry, zucry, zwcry, zrhcry, zcwcry, zcicry,                                &
!! ztcrx, ztecrx, zucrx, zvcrx, zwcrx, zrhcrx, zcwcrx, zcicrx)

      if (cross_number > 0) then

        call write_grib2_cross_space (idate0, iperiod, iperiod_accum, idatec,                     &
 cross_number, nlivz, npoint_cross_max, npoint_cross(1:cross_number), zlcr,                       &
 lonini_cross(1:cross_number), latini_cross(1:cross_number),                                      &
 lonfin_cross(1:cross_number), latfin_cross(1:cross_number),                                      &
 x0, y0, orog_cross(:,1:cross_number),                                                            &
 t_cross(:,:,1:cross_number), theta_cross(:,:,1:cross_number), thetae_cross(:,:,1:cross_number),  &
 u_tang_cross(:,:,1:cross_number), v_norm_cross(:,:,1:cross_number), w_cross(:,:,1:cross_number), &
 rh_cross(:,:,1:cross_number), qcw_cross(:,:,1:cross_number), qci_cross(:,:,1:cross_number))

      else

        print *,'No free cross-section is defined'

      endif

#endif

  endif ! Space cross-sections

endif ! output_format_grib2

      prectot(:,:)=0.
      precconv(:,:)=0.
      precsolid(:,:)=0.
      runoff(:,:)=0.
      runoff_tot(:,:)=0.
      cswfl(:,:)=0.
      clwfl(:,:)=0.
      chflux(:,:)=0.
      cqflux(:,:)=0.
      t2min(:,:)=999.
      t2max(:,:)=0.
      ws10max(:,:)=0.
      cwvflux(:,:)=0.
      csoilh_bott_flux(:,:)=0.
      csoilw_bott_flux(:,:)=0.

      write(*, '(a,i4,a)') ' Instant number',iist,' post-processed'

      endif ! mod(iist+iist0,njump) == 0

 enddo ! while (.true.)

stop
end
!###############################################################################################################
      subroutine rdmhf_atm(kunit, nlon, nlat, nlev, nfdr, pdr, p, u, v, w, t, q, qcw, qci)

!   Read Model History File with atmospheric variables from kunit

      use mod_postmoloch, only : ierr_read1

      implicit none

      integer                           :: kunit, nlon, nlat, nlev, nlevg, jlon, jlat, jklev
      integer, dimension(50)            :: nfdr
      real, dimension(100)              :: pdr
      real, dimension(nlon,nlat,nlev)   :: p, u, v, t, q, qcw, qci
      real, dimension(nlon,nlat,nlev+1) :: w

      read(kunit, iostat=ierr_read1) nfdr

      if (ierr_read1 /= 0) then
        close(kunit)
        print*,'EOF reached on file unit ', kunit
        return
      endif

      read(kunit) pdr

      do jklev=1,nlev
      do jlat=1,nlat
        read (kunit) (p(jlon,jlat,jklev),jlon=1,nlon)
      enddo
      enddo
      do jklev=1,nlev
      do jlat=1,nlat
        read(kunit) (u(jlon,jlat,jklev),jlon=1,nlon)
      enddo
      enddo
      do jklev=1,nlev
      do jlat=1,nlat
        read(kunit) (v(jlon,jlat,jklev),jlon=1,nlon)
      enddo
      enddo
      do jklev=1,nlev+1
      do jlat=1,nlat
        read(kunit) (w(jlon,jlat,jklev),jlon=1,nlon)
      enddo
      enddo
      do jklev=1,nlev
      do jlat=1,nlat
        read(kunit) (t(jlon,jlat,jklev),jlon=1,nlon)
      enddo
      enddo
      do jklev=1,nlev
      do jlat=1,nlat
        read(kunit) (q(jlon,jlat,jklev),jlon=1,nlon)
      enddo
      enddo
      do jklev=1,nlev
      do jlat=1,nlat
        read(kunit) (qcw(jlon,jlat,jklev),jlon=1,nlon)
      enddo
      enddo
      do jklev=1,nlev
      do jlat=1,nlat
        read(kunit) (qci(jlon,jlat,jklev),jlon=1,nlon)
      enddo
      enddo

      return
      end
!###############################################################################################################
      subroutine rdmhf_soil(kunit, nlon, nlat, nlevg, nlev_snow, nfdr, pdr, &
                     rgm, rgq, fice, iceth, albedo, emis1, emis2, &
                     cloudt, prectotr, precconvr, precsolidr, &
                     tskin, tgsurf, tg, qskin, qgsurf, qg, fice_soil_surf, fice_soil, &
                     snow, snow_lev, snow_t, snow_fice, snow_age, snow_melt_age, snow_dens, snow_albedo, &
                     cswflr, clwflr, chfluxr, cqfluxr, t2minr, t2maxr, ws10maxr, runoffr, runoff_tot_r, &
                     cwvfluxr, csoilh_bott_fluxr, csoilw_bott_fluxr, snowfall_level)

!   Read Model History File with surface/sea/soil variables from kunit

      use mod_postmoloch, only : ierr_read2

      implicit none

      integer                           :: kunit, nlon, nlat, nlevg, nlev_snow,  jlon, jlat, jklev, ird
      integer, dimension(50)            :: nfdr, nfdr_local
      real, dimension(100)              :: pdr, pdr_local
      real, dimension(nlon,nlat)        :: rgm, rgq, fice, iceth, albedo, emis1, emis2, &
                                           cloudt, prectotr, precconvr, precsolidr, &
                                           tskin, tgsurf, qskin, qgsurf, fice_soil_surf, &
                                           snow, snow_albedo, cswflr, clwflr, chfluxr, cqfluxr, t2minr, t2maxr, ws10maxr, &
                                           runoffr, runoff_tot_r, cwvfluxr, csoilh_bott_fluxr, csoilw_bott_fluxr
      real, dimension(nlon,nlat,nlevg)  :: tg, qg, fice_soil
      real, dimension(nlon,nlat,nlev_snow) :: snow_lev, snow_t, snow_fice, snow_age, snow_melt_age, snow_dens
      integer, dimension(nlon,nlat)     :: snowfall_level
      real, dimension(nlon,nlat)        :: field2d_add

      read(kunit, iostat=ierr_read2) nfdr_local

      if (ierr_read2 /= 0) then
        close(kunit)
        print*,'EOF reached on file unit ', kunit
        return
      endif

      read(kunit) pdr_local

      if (any(nfdr_local(:) /= nfdr(:)).or.any(pdr_local(:) /= pdr(:))) then
          print *, ' Header parameter (nfdr, pdr) read from unit ',kunit, &
 ' not coincide with defined parameters'
          print *,"   STOP"
          stop
      endif

!  Physiographical parameters changing in time

      do jlat=1,nlat
        read (kunit) (field2d_add(jlon,jlat),jlon=1,nlon) ! lai
      enddo

      do jlat=1,nlat
        read (kunit) (field2d_add(jlon,jlat),jlon=1,nlon) ! fveg
      enddo

      do jlat=1,nlat
        read (kunit) (rgm(jlon,jlat),jlon=1,nlon)
      enddo

      do jlat=1,nlat
        read (kunit) (rgq(jlon,jlat),jlon=1,nlon)
      enddo

      do jlat=1,nlat
        read (kunit) (iceth(jlon,jlat),jlon=1,nlon)
      enddo

      do jlat=1,nlat
        read (kunit) (fice(jlon,jlat),jlon=1,nlon)
      enddo

      do jlat=1,nlat
        read (kunit) (albedo(jlon,jlat),jlon=1,nlon)
      enddo

      do jlat=1,nlat
        read (kunit) (emis1(jlon,jlat),jlon=1,nlon)
      enddo

      do jlat=1,nlat
        read (kunit) (emis2(jlon,jlat),jlon=1,nlon)
      enddo

! Prognostic cloud and precipitation variables at the surface

      do jlat=1,nlat
        read(kunit) (cloudt(jlon,jlat),jlon=1,nlon)
      enddo

      do jlat=1,nlat
        read(kunit) (prectotr(jlon,jlat),jlon=1,nlon)
      enddo

!      do jlat=1,nlat
!        read(kunit) (precconvr(jlon,jlat),jlon=1,nlon)
!      enddo

      do jlat=1,nlat
        read(kunit) (precsolidr(jlon,jlat),jlon=1,nlon)
      enddo

! Prognostic surface and soil/sea fields

      do jlat=1,nlat
        read(kunit) (tskin(jlon,jlat),jlon=1,nlon)
      enddo

      do jlat=1,nlat
        read (kunit) (tgsurf(jlon,jlat),jlon=1,nlon)
      enddo

      do jklev=1,nlevg
      do jlat=1,nlat
        read(kunit) (tg(jlon,jlat,jklev),jlon=1,nlon)
      enddo
      enddo

      do jlat=1,nlat
        read(kunit) (qskin(jlon,jlat),jlon=1,nlon)
      enddo

      do jlat=1,nlat
        read (kunit) (qgsurf(jlon,jlat),jlon=1,nlon)
      enddo

      do jklev=1,nlevg
      do jlat=1,nlat
        read(kunit) (qg(jlon,jlat,jklev),jlon=1,nlon)
      enddo
      enddo

      do jlat=1,nlat
        read (kunit) (fice_soil_surf(jlon,jlat),jlon=1,nlon)
      enddo

      do jklev=1,nlevg
      do jlat=1,nlat
        read(kunit) (fice_soil(jlon,jlat,jklev),jlon=1,nlon)
      enddo
      enddo

      do jlat=1,nlat
        read(kunit) (snow(jlon,jlat),jlon=1,nlon)
      enddo

      do jklev=1,nlev_snow
      do jlat=1,nlat
        read(kunit) (snow_lev(jlon,jlat,jklev),jlon=1,nlon)
      enddo
      enddo

      do jklev=1,nlev_snow
      do jlat=1,nlat
        read(kunit) (snow_t(jlon,jlat,jklev),jlon=1,nlon)
      enddo
      enddo

      do jklev=1,nlev_snow
      do jlat=1,nlat
        read(kunit) (snow_fice(jlon,jlat,jklev),jlon=1,nlon)
      enddo
      enddo

      do jklev=1,nlev_snow
      do jlat=1,nlat
        read(kunit) (snow_age(jlon,jlat,jklev),jlon=1,nlon)
      enddo
      enddo

      do jklev=1,nlev_snow
      do jlat=1,nlat
        read(kunit) (snow_melt_age(jlon,jlat,jklev),jlon=1,nlon)
      enddo
      enddo

      do jklev=1,nlev_snow
      do jlat=1,nlat
        read(kunit) (snow_dens(jlon,jlat,jklev),jlon=1,nlon)
      enddo
      enddo

      do jlat=1,nlat
        read(kunit) (snow_albedo(jlon,jlat),jlon=1,nlon)
      enddo

      do jlat=1,nlat
        read(kunit) (cswflr(jlon,jlat),jlon=1,nlon)
      enddo

      do jlat=1,nlat
        read(kunit) (clwflr(jlon,jlat),jlon=1,nlon)
      enddo

      do jlat=1,nlat
        read(kunit) (chfluxr(jlon,jlat),jlon=1,nlon)
      enddo

      do jlat=1,nlat
        read(kunit) (cqfluxr(jlon,jlat),jlon=1,nlon)
      enddo

      do jlat=1,nlat
        read(kunit) (t2minr(jlon,jlat),jlon=1,nlon)
      enddo

      do jlat=1,nlat
        read(kunit) (t2maxr(jlon,jlat),jlon=1,nlon)
      enddo

      do jlat=1,nlat
        read(kunit) (ws10maxr(jlon,jlat),jlon=1,nlon)
      enddo

      do jlat=1,nlat
        read(kunit) (runoffr(jlon,jlat),jlon=1,nlon)
      enddo

      do jlat=1,nlat
        read(kunit) (runoff_tot_r(jlon,jlat),jlon=1,nlon)
      enddo

      do jlat=1,nlat
        read(kunit) (cwvfluxr(jlon,jlat),jlon=1,nlon)
      enddo

      do jlat=1,nlat
        read(kunit) (csoilh_bott_fluxr(jlon,jlat),jlon=1,nlon)
      enddo

      do jlat=1,nlat
        read(kunit) (csoilw_bott_fluxr(jlon,jlat),jlon=1,nlon)
      enddo

      do jlat=1,nlat
        read(kunit) (field2d_add(jlon,jlat),jlon=1,nlon)
      enddo
      snowfall_level(:,:) = int(field2d_add(:,:))

! Reading of additional 2d fields

     do ird=1,1
      do jlat=1,nlat
      read(kunit) (field2d_add(jlon,jlat),jlon=1,nlon)
      enddo
     enddo

      return
      end
!###############################################################################################################
      subroutine rd_param_const(nlon, nlat, nlevg, x0, y0, alon0, alat0, dlon, dlat, phig, fmask, qgmax, qgmin)

! Reads from additional input file
! all constant (in time) model physiographical parameters

      implicit none

      integer                           :: nlon, nlat, nlev, nlevg, &
                                           nlon_local, nlat_local, nlev_local, nlevg_local, &
                                           iunit=11, ierr_open, ierr, nst_local, nvt_local, &
 jlon, jlat, ird, jklev
      real                              :: x0, y0, alon0, alat0, dlon, dlat, &
                                           x0_local, y0_local, alon0_local, alat0_local, dlon_local, dlat_local
      real, parameter                   :: g=9.807
      real, dimension(nlon,nlat)        :: phig, fmask, zread
      real, dimension(nlon,nlat,nlevg)  :: qgmax, qgmin
      character(len=30) :: filerd="model_param_constant.bin"

     open (iunit,file=trim(filerd),status='old',form='unformatted',iostat=ierr_open)
     if (ierr_open /= 0) then
        print *
        print *,'Not found input ',trim(filerd)
        print *,"   stop"
        stop
      endif

      read (iunit) nlon_local, nlat_local, nlevg_local, dlon_local, dlat_local, x0_local, y0_local, alon0_local, alat0_local, &
 nst_local, nvt_local

      if (nlon_local /= nlon) ierr=ierr+1
      if (nlat_local /= nlat) ierr=ierr+1
      if (nlevg_local /= nlevg) ierr=ierr+1
      if (dlon_local /= dlon) ierr=ierr+1
      if (dlat_local /= dlat) ierr=ierr+1
      if (x0_local /= x0) ierr=ierr+1
      if (y0_local /= y0) ierr=ierr+1
      if (alon0_local /= alon0) ierr=ierr+1
      if (alat0_local /= alat0) ierr=ierr+1

      if (ierr /= 0) then
        print *,"Error in header parameters in input file ,",trim(filerd),", not coincident with defined parameters"
        print *,"Model nlon, nlat, nlevg, dlon, dlat, x0, y0, alon0, alat0 :", &
 nlon, nlat, nlevg, dlon, dlat, x0, y0, alon0, alat0
        print *,"Read nlon, nlat, nlevg, dlon, dlat, x0, y0, alon0, alat0 :", &
 nlon_local, nlat_local, nlevg_local, dlon_local, dlat_local, x0_local, y0_local, alon0_local, alat0_local
        print *,"   stop"
        stop
      endif

!  land-sea fraction and orography

      do jlat=1,nlat
        read(iunit) (fmask(jlon,jlat),jlon=1,nlon)
      enddo

      do jlat=1,nlat
        read(iunit) (phig(jlon,jlat),jlon=1,nlon)
      enddo
      phig(:,:) = phig(:,:)*g

      do jlat=1,nlat
        read(iunit) (zread(jlon,jlat),jlon=1,nlon) ! htopvar
      enddo
      do ird=1,nst_local+1
        do jlat=1,nlat
          read(iunit) (zread(jlon,jlat),jlon=1,nlon) ! soil_map
        enddo
      enddo
      do ird=1,nvt_local+1
        do jlat=1,nlat
          read(iunit) (zread(jlon,jlat),jlon=1,nlon) ! veg_map
        enddo
      enddo

      do jklev=1,nlevg
        do jlat=1,nlat
          read(iunit) (qgmax(jlon,jlat,jklev),jlon=1,nlon)
        enddo
      enddo

      do jklev=1,nlevg
        do jlat=1,nlat
          read(iunit) (qgmin(jlon,jlat,jklev),jlon=1,nlon)
        enddo
      enddo

      close(iunit)

      return
      end
!###############################################################################################################
      subroutine wrpost(a, izout, nlon, nlat, nfdr, pdr)

!  Writes horizontal fields (different variables) on unit 80
!  Fields are written as positive integers (4 digits) after finding the best rescaling
!  (zfac: scaling factor; zoffs: offset)

 use mod_postmoloch

      implicit none
      integer                        :: nlon, nlat, mhfr, jlon, jlat, ix1, ix2, iy1, iy2, nhist, iminu, ihour
      real                           :: dlat, dlon, alat0, alon0, dtstep, zmin, zmax, zoffs, zfac, zy0, zx0
      integer, dimension(nlon,nlat)  :: izout
      integer, dimension(50)         :: nfdr
      real,    dimension(nlon,nlat)  :: a
      real,    dimension(100)        :: pdr
      real :: aihour

!  Def. of geographical parameters for output (t-points)

      nhist  = nfdr(17)
      mhfr   = nfdr(20)
      dlat   = pdr(1)
      dlon   = pdr(2)
      dtstep = pdr(3)
      alat0  = pdr(4) + 0.5*dlat/float(mhfr)
      alon0  = pdr(5)
      zy0    = pdr(38)
      zx0    = pdr(39)

!  Comp. of accumulation time (hours and minutes) for precipitation

      iminu = nint((nhist*dtstep)/60.)
      aihour = iminu/60
      ihour = int(aihour)
      iminu = iminu - ihour*60

! The following to accumulate precip. on multiple intervals

      ihour = ihour*njump
      iminu = iminu*njump

      aihour = iminu/60.
      ihour = ihour + int(aihour)
      iminu = iminu - int(aihour)*60

!  Zoom (i1 starting point, i2 ending point)

      ix1=1           ! longitude
      ix2=nlon
      iy1=1           ! latitude
      iy2=nlat

!  Redefinition of very small values of za (may be necessary)

      do jlat=iy1, iy2
      do jlon=ix1, ix2
      if(abs(a(jlon,jlat)).lt.1.e-20) a(jlon,jlat)=0.
      enddo
      enddo

!  Definition of offset and scaling factors

      zmin = 1.e35
      zmax =-1.e35
      do jlat = iy1, iy2
      do jlon = ix1, ix2
      zmax = max(zmax,a(jlon,jlat))
      zmin = min(zmin,a(jlon,jlat))
      enddo
      enddo
      zoffs = zmin

      if ((zmax-zmin).ge.1.e-33) then
      zfac = 9999./(zmax-zmin)
      elseif (abs(zmax).gt.1.e-32) then
      zfac = 9999./zmax
      else
      zfac = 1.e35
      endif

      do jlat = iy1, iy2
      do jlon = ix1, ix2
      izout(jlon,jlat) = nint((a(jlon,jlat)-zoffs)*zfac)
      enddo
      enddo

      write(80,9000) nfdr(5:12), ihour, iminu, zoffs, zfac
      write(80,9200) iy2-iy1+1, ix2-ix1+1, alat0, alon0, dlat, dlon, zy0, zx0
      write(80,9100) ((izout(jlon,jlat), jlon=ix1,ix2), jlat=iy1,iy2)

 9000 format(10i10,2(1x,e12.5))
 9100 format(20i4)
 9200 format(2i5,6f11.5)

      return
      end
!###############################################################################################################
      subroutine wrposx(a, izout, zpz, zlcr, ncr, nlon, nlat, nlivz, nfdr, pdr)

!  Writes zonal cross-sections at latitude ncr on unit 82
!  a is the matrix containing the main output field on cross-sect.
!  zlcr is the vector of cross-section levels
!  zpz contains orography along the cross-sect.
!  Fields are written as positive integers (4 digits) after finding the best rescaling
!  (zfac: scaling factor; zoffs: offset)

      implicit none
      integer                        :: ncr, nlon, nlat, nlivz, jlon, jlev, i1, i2
      real                           :: dlat, dlon, alat0, alon0, zalat, zmin, zmax, zoffs, zfac
      integer, dimension(50)         :: nfdr
      integer, dimension(nlon,nlivz) :: izout
      real,    dimension(100)        :: pdr
      real,    dimension(nlon,nlivz) :: a
      real,    dimension(nlon)       :: zpz
      real,    dimension(nlivz)      :: zlcr

!  Def. of geographical parameters for output (t-points)

      dlat  = pdr(1)
      dlon  = pdr(2)
      alat0 = pdr(4) + dlat/2.
      alon0 = pdr(5)

!  Comp. of latitude of cross-section (in model coordinates)

      zalat = alat0 + float(ncr-1)*dlat

!  Zoom (i1 starting point, i2 ending point)

      i1=1
      i2=nlon

! Redefinition of very small values of za (may be necessary)

      do jlev=1, nlivz
      do jlon=i1, i2
      if(abs(a(jlon,jlev)).lt.1.e-20) a(jlon,jlev)=0.
      enddo
      enddo

!  Definition of offset and scaling factors

      zmin = 1.e35
      zmax = -1.e35
      do jlev = 1, nlivz
      do jlon= i1, i2
      zmax = max(zmax,a(jlon,jlev))
      zmin = min(zmin,a(jlon,jlev))
      enddo
      enddo
      zoffs = zmin

      if ((zmax-zmin).ge.1.e-33) then
      zfac = 9999./(zmax-zmin)
      elseif (abs(zmax).gt.1.e-32) then
      zfac = 9999./zmax
      else
      zfac = 1.e35
      endif

      do jlev = 1, nlivz
      do jlon = i1, i2
      izout(jlon,jlev)=nint((a(jlon,jlev)-zoffs)*zfac)
      enddo
      enddo

      write(82,9000) nfdr(5:12), zoffs, zfac
      write(82,9200) i2-i1+1, nlat, nlivz, alon0, alat0, dlon, dlat, zalat
      write(82,9300) (zpz(jlon), jlon=i1,i2)
      write(82,9300) zlcr
      write(82,9100) ((izout(jlon,jlev), jlon=i1,i2), jlev=1,nlivz)

 9000 format(8i10,2(1x,e12.5))
 9100 format(20i4)
 9200 format(3i5,5f11.5)
 9300 format(8f11.5)

      return
      end
!###############################################################################################################
      subroutine wrposy(a, izout, zpz, zlcr, ncr, nlon, nlat, nlivz, nfdr, pdr)

!  Writes meridional cross-sections at longitude ncr on unit 81
!  a is the matrix containing the cross-sect.
!  zlcr is the vector of cross-section levels
!  zpz contains orography along the cross-sect.
!  Fields are written as positive integers (4 digits) after finding the best rescaling
!  (zfac: scaling factor; zoffs: offset)

      implicit none
      integer                        :: ncr, nlon, nlat, nlivz, jlat, jlev, i1, i2
      real                           :: dlat, dlon, alat0, alon0, zalon, zmin, zmax, zoffs, zfac
      integer, dimension(50)         :: nfdr
      integer, dimension(nlat,nlivz) :: izout
      real,    dimension(100)        :: pdr
      real,    dimension(nlat,nlivz) :: a
      real,    dimension(nlat)       :: zpz, zaux
      real,    dimension(nlivz)      :: zlcr

!  Def. of geographical parameters for output (t-points)

      dlat  = pdr(1)
      dlon  = pdr(2)
      alat0 = pdr(4) + dlat/2.
      alon0 = pdr(5)

!  Comp. of longitude of cross-section (in model coordinates)

      zalon = alon0 + float(ncr-1)*dlon

!  Zoom (i1 starting point, i2 ending point)

      i1=1
      i2=nlat

! Redefinition of very small values of za (may be necessary)

      do jlev=1, nlivz
      do jlat=i1, i2
      if(abs(a(jlat,jlev)).lt.1.e-20) a(jlat,jlev)=0.
      enddo
      enddo

!  Definition of offset and scaling factors

      zmin = 1.e35
      zmax = -1.e35
      do jlev = 1, nlivz
      do jlat = i1, i2
      zmax = max(zmax,a(jlat,jlev))
      zmin = min(zmin,a(jlat,jlev))
      enddo
      enddo
      zoffs = zmin

      if ((zmax-zmin).ge.1.e-33) then
      zfac = 9999./(zmax-zmin)
      elseif (abs(zmax).gt.1.e-32) then
      zfac = 9999./zmax
      else
      zfac = 1.e35
      endif

      do jlev = 1, nlivz
      do jlat = i1, i2
      izout(i2-jlat+i1,jlev)=nint((a(jlat,jlev)-zoffs)*zfac)
      enddo
      enddo

      do jlat=i1,i2
      zaux(i2-jlat+i1)=zpz(jlat)
      enddo

      write(81,9000) nfdr(5:12), zoffs, zfac
      write(81,9200) i2-i1+1, nlon, nlivz, alat0, alon0, dlat, dlon, zalon
      write(81,9300) (zaux(jlat), jlat=i1,i2)
      write(81,9300) zlcr
      write(81,9100) ((izout(jlat,jlev), jlat=i1,i2), jlev=1,nlivz)

 9000 format(8i10,2(1x,e12.5))
 9100 format(20i4)
 9200 format(3i5,5f11.5)
 9300 format(8f11.5)

      return
      end
!###############################################################################################################
      subroutine vtsurf (u10, v10, t2, q2, hflux, qflux, nlon, nlat, nlev, fmask, rgm, rgq, &
                         richs, phig, ps, fsnow, tskin, qskin, u, v, t, p, q, h, dz, a0, b0, mhfr)

      use data_observ_aux
      use mod_postmoloch, only : gzita, bzita

!  Interpolates wind at 10 m, temperature and spec. humid. at 2 m.

      real, dimension(nlon,nlat,nlev) :: u, v, t, p, q
      real, dimension(nlon,nlat)      :: u10, v10, t2, td2, q2, q2rel, hflux, qflux, richs, &
                                         fmask, phig, rgm, rgq, ps, fsnow, tskin, qskin
      integer :: mhfr

      real, parameter :: zaholt=1., zbholt=2./3., zcholt=5., zdholt=0.35,                                &
                         ak=.4, zbet=5., zgam=16., yliv=2834170.5, yliw=333560.5, ylwv=yliv-yliw,        &
                         tzer=273.15, pzer=1.e5, ezer= 611., rd=287.05, rv=461.51, eps=rd/rv,            &
                         cpd=1004.6, cvd=cpd-rd, cpv=1869.46, cw=4186.8, ci=cw/2., rdrcp=rd/cpd, g=9.807

!  Constants for computing saturation partial pressure

      real, parameter :: ccw1=(cpv-cw)/rv, ccw2=ylwv/tzer/rv-ccw1, cci1=(cpv-ci)/rv, cci2=yliv/tzer/rv-cci1

!  Businger functions when Ri<0

      psium(zx1,zx2)=alog((1.+zx1)**2*(1.+zx1**2)/((1.+zx2)**2*(1.+zx2**2)))-2.*(atan(zx1)-atan(zx2))
      psiuh(zy1,zy2)=2.*alog((1.+zy1)/(1.+zy2))

!  Holtslag functions when Ri>0

      psism(zx1,zx2)=-zaholt*zx1-zbholt*(zx1-zcholt/zdholt)*exp(-zdholt*zx1) &
                     +zaholt*zx2+zbholt*(zx2-zcholt/zdholt)*exp(-zdholt*zx2)
      psish(zy1,zy2)=-(1.+2./3.*zaholt*zy1)**1.5-zbholt*(zy1-zcholt/zdholt)*exp(-zdholt*zy1) &
                     +(1.+2./3.*zaholt*zy2)**1.5+zbholt*(zy2-zcholt/zdholt)*exp(-zdholt*zy2)

!  Turbulent fluxes at the ground (positive upward)

      zep=1./eps-1.

      do 100 jlat=1,nlat
      jlatp1=min(jlat+1,nlat)
      do 100 jlon=1,nlon
      jlonm1=max(jlon-1,1)

      zphi  = phig(jlon,jlat)*gzita(.5*dz,h,a0)-g*h*bzita(.5*dz,h,b0)*log(1.-.5*dz/h)
      za    = (zphi-phig(jlon,jlat))/g
      zua   = 0.5*(u(jlon,jlat,1)+u(jlonm1,jlat,1))
      zva   = 0.5*(v(jlon,jlatp1,1)+v(jlon,jlat,1))
      zua   = ((2.*float(mhfr)-1.)*u(jlon,jlat,1)+u(jlonm1,jlat,1))/(2.*float(mhfr))  ! destaggering
      zva   = (v(jlon,jlatp1,1)+(2.*float(mhfr)-1.)*v(jlon,jlat,1))/(2.*float(mhfr))  ! destaggering
      zmod2  = zua**2 + zva**2 + 0.07
      zmod  = sqrt(zmod2)
      zros  = ps(jlon,jlat)/( rd*tskin(jlon,jlat)*(1.+zep*qskin(jlon,jlat)) )

!  Virtual potential temperature computed with skin temperature

      zconvg = (pzer/ps(jlon,jlat) )**rdrcp*(1.+zep*qskin(jlon,jlat))
      ztevg  = tskin(jlon,jlat)*zconvg

      zconv1 = (pzer/p(jlon,jlat,1))**rdrcp*(1.+zep*q(jlon,jlat,1)  )
      ztetav = t(jlon,jlat,1)*zconv1

!  Richardson number above the ground

      ztebar = .5*(ztevg+ztetav)
      zri    = za*g*(ztetav-ztevg)/(ztebar*zmod2)
      richs(jlon,jlat) = zri ! to be used for the definition of clouds

!  Definition of roughness zrgm, zrgt, zrgq

      if(fmask(jlon,jlat).gt..5) then

!  Computation of Charnok roughness

      zchar = 5.e-4
      zcoch1=.0185*ak**2/g*zmod2
        do jiter = 1, 5
        zcoch2= za/zchar
        zchar = zcoch1/alog(zcoch2)**2
        enddo
      zrgm=zchar

        if(zri.ge.0.25) then
        zrgt = 2.2e-9
        zrgq = zrgt
        elseif(zri.lt.0.) then
        zrgt = 5.e-5
        zrgq = 9.e-5
        else
        zrgt = 3.80227e-6*exp(-zri*29.819387)
        zrgq = zrgt
        endif

        zrgmd=zrgm
        zrgtd=zrgt
        zrgqd=zrgq
      else
      zsea=5.e-5

      zrgm=fmask(jlon,jlat)*zsea+(1.-fmask(jlon,jlat))*rgm(jlon,jlat)
      zrgt=fmask(jlon,jlat)*zsea+(1.-fmask(jlon,jlat))*rgq(jlon,jlat)
      zrgq=zrgt

! Local roughness set to a few cm if larger (for interpolated
! diagnostics at 2 m (and other near surface levels) and 10 m only,
! not for computing fluxes)

      zrgmd=fmask(jlon,jlat)*zsea+(1.-fmask(jlon,jlat))*min(rgm(jlon,jlat), .05)
      zrgtd=fmask(jlon,jlat)*zsea+(1.-fmask(jlon,jlat))*min(rgq(jlon,jlat), .05)
      zrgqd=zrgtd
      endif

!  Computation of za/l

      zalzam=alog((za+zrgm)/zrgm)
      zalzat=alog((za+zrgt)/zrgt)
      zalzaq=alog((za+zrgq)/zrgq)

      zalzamd=alog((za+zrgmd)/zrgmd)
      zalzatd=alog((za+zrgtd)/zrgtd)
      zalzaqd=alog((za+zrgqd)/zrgqd)

      if(zri.ge.0.) then      ! Holtslag functions

      zpsim=0.
      zpsimd=0.
      zpsih=0.
      zpsihd=0.

      do jiter=1,4
      zal=zri*(zalzam-zpsim)**2/(zalzat-zpsih)
      zald=zri*(zalzamd-zpsimd)**2/(zalzatd-zpsihd)
      zx1=zal*(za+zrgm)/za
      zx1d=zald*(za+zrgmd)/za
      zx2=zal*zrgm/za
      zx2d=zald*zrgmd/za
      zy1=zal*(za+zrgt)/za
      zy1d=zald*(za+zrgtd)/za
      zy2=zal*zrgt/za
      zy2d=zald*zrgtd/za
      zpsim=psism(zx1,zx2)
      zpsimd=psism(zx1d,zx2d)
      zpsih=psish(zy1,zy2)
      zpsihd=psish(zy1d,zy2d)
      enddo
      zphim=1.+zaholt*zal+zbholt*zal*(1.+zcholt-zdholt*zal)*exp(-zdholt*zal)
      zphimd=1.+zaholt*zald+zbholt*zald*(1.+zcholt-zdholt*zald)*exp(-zdholt*zald)
      zpsiq=zpsih
      zpsiqd  = zpsihd
      zpsim10= zpsimd*(10./za)
      zpsih2 = zpsihd*(2. /za)
      zpsiq2 = zpsih2
!--------------------------------------------------
          zpsih05 = zpsih*(0.5 /za)
          zpsih005 = zpsih*(0.05 /za)
          zpsiq05 = zpsih05
!--------------------------------------------------

      else                    ! Businger functions

      zpsim=0.
      zpsimd=0.
      zpsih=0.
      zpsihd=0.
      do jiter=1,4
      zal=zri*(zalzam-zpsim)**2/(zalzat-zpsih)
      zald=zri*(zalzamd-zpsimd)**2/(zalzatd-zpsihd)
      zx1=(1.-zgam*zal*(za+zrgm)/za)**.25
      zx1d=(1.-zgam*zald*(za+zrgmd)/za)**.25
      zx2=(1.-zgam*zal*zrgm/za)**.25
      zx2d=(1.-zgam*zald*zrgmd/za)**.25
      zy1=sqrt(1.-zgam*zal*(za+zrgt)/za)
      zy1d=sqrt(1.-zgam*zald*(za+zrgtd)/za)
      zy2=sqrt(1.-zgam*zal*zrgt/za)
      zy2d=sqrt(1.-zgam*zald*zrgtd/za)
      zpsim=psium(zx1,zx2)
      zpsimd=psium(zx1d,zx2d)
      zpsih=psiuh(zy1,zy2)
      zpsihd=psiuh(zy1d,zy2d)
      enddo
      zz1=sqrt(1.-zgam*zal*(za+zrgq)/za)
      zz1d=sqrt(1.-zgam*zald*(za+zrgqd)/za)
      zz2=sqrt(1.-zgam*zal*zrgq/za)
      zz2d=sqrt(1.-zgam*zald*zrgqd/za)
      zpsiq=psiuh(zz1,zz2)
      zpsiqd=psiuh(zz1d,zz2d)
      zx1d=(1.-zgam*(10.+zrgmd)/za*zald)**0.25
      zy1d=sqrt(1.-zgam*zald*(2.+zrgtd)/za)
      zz1d=sqrt(1.-zgam*zald*(2.+zrgqd)/za)
      zpsim10=psium(zx1d,zx2d)
      zpsih2 =psiuh(zy1d,zy2d)
      zpsiq2 =psiuh(zz1d,zz2d)
!--------------------------------------------------
          zyy1=sqrt(1.-zgam*zal*(0.5+zrgt)/za)
          zyyy1=sqrt(1.-zgam*zal*(0.05+zrgt)/za)
          zzz1=sqrt(1.-zgam*zal*(0.5+zrgq)/za)
          zpsih05 =psiuh(zyy1,zy2)
          zpsih005 =psiuh(zyyy1,zy2)
          zpsiq05 =psiuh(zzz1,zz2)
!--------------------------------------------------

      endif

      zustar = ak*zmod/(zalzam-zpsim)
      ztstar = ak*(ztetav-ztevg)/(zalzat-zpsih)
      zqstar = ak*(q(jlon,jlat,1)-qskin(jlon,jlat))/(zalzaq-zpsiq)

      zustard = ak*zmod/(zalzamd-zpsimd)
      ztstard = ak*(ztetav-ztevg)/(zalzatd-zpsihd)
      zqstard = ak*(q(jlon,jlat,1)-qskin(jlon,jlat))/(zalzaqd-zpsiqd)

!  Surface fluxes of sensible and latent heat (positive upward)

      zcdt = zustar*ak/(zalzat-zpsih)
      zperflut = -zros*zcdt*cpd/zconvg
      hflux(jlon,jlat)=zperflut*(ztetav-ztevg)

      zcdq = zustar*ak/(zalzaq-zpsiq)
      zlate=ylwv-2360.*(tskin(jlon,jlat)-tzer)
      if(tskin(jlon,jlat).lt.tzer) zlate=zlate+fsnow(jlon,jlat)*(yliw+2255.*(tskin(jlon,jlat)-tzer))
      zperfluq = -zros*zcdq*zlate
      qflux(jlon,jlat)=zperfluq*(q(jlon,jlat,1)-qskin(jlon,jlat))
!-----------------------------------------------------------------------
      zv10=zustard/ak*(alog((10.+zrgmd)/zrgmd)-zpsim10)
      zt2 =ztevg+ztstard/ak*(alog((2. +zrgtd)/zrgtd)-zpsih2)
      zq2 =qskin(jlon,jlat)+zqstard/ak*(alog((2. +zrgqd)/zrgqd)-zpsiq2)
!-----------------------------------------------------------------------
      zt05 =ztevg+ztstar/ak*(alog((0.5 +zrgt)/zrgt)-zpsih05)
      zt005 =ztevg+ztstar/ak*(alog((0.05 +zrgt)/zrgt)-zpsih005)
      zq05 =qskin(jlon,jlat)+zqstar/ak*(alog((0.5 +zrgq)/zrgq)-zpsiq05)
!--------------------------------------------------

!  Wind components

      u10(jlon,jlat) = zua/zmod*zv10
      v10(jlon,jlat) = zva/zmod*zv10

!  From potential temperature to temperature at 2 m in Celsius

      zconv2 = (pzer/ps(jlon,jlat) )**rdrcp*(1.+zep*zq2)
      zt2=zt2/zconv2
      t2(jlon,jlat) = zt2-tzer
      q2(jlon,jlat) = zq2

!--------------------------------------------------
! Temperature at 0.5 m and 0.05 m and relative humidity at 0.5 m
! (with respect to water! - WMO, not standard)

      t05(jlon,jlat)=zt05/zconv2-tzer   ! approx. using zconv2
      t005(jlon,jlat)=zt005/zconvg-tzer ! approx. using zconvg
      zt05=zt05/zconv2
      call comp_esk(zesk, zqsat, zt05, ps(jlon,jlat), 3)  ! partial pressure over water
      q05p=min(zq05,zqsat*1.01) ! ensures that q05p does not exceed 101% (1% mor for graphic smoothness)
      eee=ps(jlon,jlat)*q05p/(eps*(1.-q05p)+q05p)
      q05rel(jlon,jlat) = eee/zesk*100.
!      q05rel(jlon,jlat)=min(zq05/zqsat*100., 102.)
!--------------------------------------------------

 100  continue

      return
      end
!###############################################################################################################
      subroutine cclouds (nlon, nlat, nlev, h, dz, a0, b0, qccrit, phig, richs, zeta, u, v, teta,   &
                          cloudt, cloudh, cloudm, cloudl, t, p, q, qcw, qci, tvirt, pbl, zorogr)

! Computes total and high, middle and low cloud fraction as a function of
! cloud water, cloud ice and relative humidity.
! Reduction of cloud fraction computed as a function of the Richardson number,
! depending on moist static stability as Durran & Klemp (1982).
! Cloud cover algorithm as revised from Geleyn's (Maurizio).
! Low, middle, high clouds as WMO definition: < 2000 m, 2000-6000 m, > 6000 m.

! ---- For RTTOV simulation ---

#ifdef rttov
 use mod_rttov, only : clfrac_invert
#endif

      real, dimension(nlon,nlat)     :: phig, richs, cloudt, cloudh, cloudm, cloudl, pbl, zorogr
      real, dimension(nlon,nlat,nlev):: t, p, q, tvirt, u, v, qcw, qci, zeta, teta
      real, dimension(nlon,nlev)     :: fcloud, zqs, zteta, ztetav, zrich
      real, dimension(nlev)          :: zcol

      real, parameter :: tzer=273.15, pzer=1.e5, ezer= 611., rd=287.05, rv=461.51, eps=rd/rv, &
                         cpd=1004.6, cpv=1869.46, rdrcp=rd/cpd, cw = 4186.8,                  &
                         yliv=2834170.5, yliw=333560.5, ylwv=yliv-yliw,                       &
                         ccw1 = (cpv-cw)/rv, ccw2 = ylwv/tzer/rv-ccw1, g = 9.807
      ntop = nlev - 3
      zqcrit = qccrit
      huc = 0.89
      pbl = 1.

! Computation of dry and moist static stability and of effective static stability
! as a function of relative humidity.
! Computation of the Richardson number.

      do jlat = 2, nlat-1

! Comput. of virtual theta

      do jklev = 1, ntop+1
      do jlon = 2, nlon-1
      zt0t = tzer/t(jlon,jlat,jklev)
      zesk = ezer*exp(-ccw1*log(zt0t)+ccw2*(1.-zt0t))     ! partial pressure over water
      zqs(jlon,jklev) = zesk*eps/(p(jlon,jlat,jklev)+zesk*(eps-1.))
      zteta (jlon,jklev) = teta(jlon,jlat,jklev)
      ztetav(jlon,jklev) = tvirt(jlon,jlat,jklev)*teta(jlon,jlat,jklev)/t(jlon,jlat,jklev)
      enddo
      enddo

!  Computation of the Richardson no. depending on dry and moist (saturated) static stability
!  as Durran & Klemp (1982). (Specific humidity is used in place of mixing ratio).

      do jklev = 2, ntop+1   ! loop on half levels
      do jlon = 2, nlon-1
      dthd   = 2.*(ztetav(jlon,jklev)-ztetav(jlon,jklev-1))/(ztetav(jlon,jklev)+ztetav(jlon,jklev-1)) ! dry
      r_up = zqs(jlon,jklev  )
      r_do = zqs(jlon,jklev-1)
      r_av = 0.5*(r_up + r_do)
      t_av = 0.5*(t(jlon,jlat,jklev-1) + t(jlon,jlat,jklev))
      theta_up = zteta(jlon,jklev)
      theta_do = zteta(jlon,jklev-1)
      zaa = (1. + ylwv*r_av/(rd*t_av))/(1. + eps*ylwv**2*r_av/(cpd*rd*t_av**2))
      dthm = zaa*((theta_up-theta_do)*2./(theta_up + theta_do) + ylwv/(cpd*t_av)*(r_up - r_do)) & ! Durran&Klemp
             - q(jlon,jlat,jklev  ) - qcw(jlon,jlat,jklev-1) -  qci(jlon,jlat,jklev  )          &
             + q(jlon,jlat,jklev-1) + qcw(jlon,jlat,jklev-1) +  qci(jlon,jlat,jklev-1)

! Average relative humidity computed giving some more weight to the layer below

      zrhm = 0.55*q(jlon,jlat,jklev-1)/zqs(jlon,jklev-1) + 0.45*q(jlon,jlat,jklev)/zqs(jlon,jklev)
      zcof = max(-24.5 + 25.*zrhm, 0.)    ! zcof=0. for rh<=0.98, zcof=0.5 for rh=1
      zcof = min(zcof, .85)
      dlogthe = zcof*dthm + (1.-zcof)*dthd ! effective stability near saturation
      zrdz  = 1./(zeta(jlon,jlat,jklev)-zeta(jlon,jlat,jklev-1))
      zdtdz = dlogthe*zrdz

      jlonm1 = max(jlon-1,1)
      jlatp1 = min(jlat+1,nlat)
      zdudz = .5*(u(jlon,jlat,jklev)+u(jlonm1,jlat,jklev)-u(jlon,jlat,jklev-1)-u(jlonm1,jlat,jklev-1))*zrdz
      zdvdz = .5*(v(jlon,jlat,jklev)+v(jlon,jlatp1,jklev)-v(jlon,jlat,jklev-1)-v(jlon,jlatp1,jklev-1))*zrdz
      zshear = zdudz**2 + zdvdz**2
      zbuoy = g*dlogthe*zrdz
      zrich(jlon,jklev) = min (zbuoy/(zshear + 1.e-6), 500.)
      enddo
      enddo   ! end of loop on half levels
      zrich(:,1) = richs(:,jlat)  ! Richardson no. at the surface (computed in vtsurf)

      do jklev = 1, ntop
      do jlon = 2, nlon-1
      call comp_esk(zesat, zqs1, t(jlon,jlat,jklev), p(jlon,jlat,jklev), 2) ! blended saturation
      call comp_esk(zesat, zqs2, t(jlon,jlat,jklev), p(jlon,jlat,jklev), 3) ! sat. to water below 0
      zzqs = 0.70*zqs1 + 0.30*zqs2                   ! ad hoc to limit contrib. to high clouds
      zrh1  =  min(1., q(jlon,jlat,jklev)/zzqs)
      zrh = (zrh1-huc)/(1.-huc)
      zrh = min(1., zrh)
      zrh = max(0., zrh)

      fcloud(jlon,jklev) = (max(0.15*zrh, (qcw(jlon,jlat,jklev) + 0.77*qci(jlon,jlat,jklev))/qccrit))**.7

! Reduction of clouds in unstable layers as a function of the local Richardson no.

      zstabr = 0.5*zrich(jlon,jklev) + 0.5*zrich(jlon,jklev+1)

      fc1 = 1.
      if(zstabr.lt..25.and.zstabr.gt.0) fc1 = .875 + 0.5*zstabr
      if(zstabr.gt.0..and.zrich(jlon,jklev).lt.0.) fc1 = 0.82
      if(zstabr.lt.0.) fc1 = 0.72
      fcloud(jlon,jklev) = fc1*fcloud(jlon,jklev)
      fcloud(jlon,jklev) = max(0., min(1., fcloud(jlon,jklev)))
      enddo
      enddo

      do 10 jlon = 2, nlon-1
      nlclo = 0
      nmclo = 0
      do jklev = 1, nlev
      if(zeta(jlon,jlat,jklev).le.2000.) nlclo = nlclo + 1 ! no. of levels above 2000 m
      if(zeta(jlon,jlat,jklev).le.6000.) nmclo = nmclo + 1 ! no. of levels above 6000 m
      enddo

 ! High clouds

      n1 = ntop
      n2 = nmclo + 1
      ntot = n1 - n2 + 1
      zcol(1:ntot) = fcloud(jlon, n1:n2:-1)
      call ccolumn (zcol, ntot, zprod)
      cloudh(jlon,jlat) = max(0., 1.-zprod)

! Middle clouds

      n1 = nmclo
      n2 = nlclo + 1
      ntot = n1 - n2 + 1
      if(ntot.ge.1) then
      zcol(1:ntot) = fcloud(jlon, n1:n2:-1)
      call ccolumn (zcol, ntot, zprod)
      cloudm(jlon,jlat) = max(0., 1.-zprod)
      endif

! Low clouds

      n1 = nlclo
      n2 = 1
      ntot = n1 - n2 + 1
      if(ntot.ge.1) then
      zcol(1:ntot) = fcloud(jlon, n1:n2:-1)
      call ccolumn (zcol, ntot, zprod)
      cloudl(jlon,jlat) = max(0., 1.-zprod)
      endif

! Total clouds: redefines cloudt in output from MOLOCH

      n1 = ntop
      n2 = 1
      ntot = n1 - n2 + 1
      if(ntot.ge.1) then
      zcol(1:ntot) = fcloud(jlon, n1:n2:-1)
      call ccolumn (zcol, ntot, zprod)
      cloudt(jlon,jlat) = max(0., 1.-zprod)
      endif

! ---- For RTTOV simulation ---

#ifdef rttov
     clfrac_invert(jlon,jlat,:) = 0.
     do jklev=1,nlev
       jklev1=nlev-jklev+1
       clfrac_invert(jlon,jlat,jklev1)=min(max(fcloud(jlon,jklev),0.),1.)
     enddo
#endif

 10   continue

! Compute height (over surface) of top of pbl

      do jlon = 2, nlon-1
      jkpbl = 0
       do jklev = 1, nlev/4
       zriav = .5*(zrich(jlon,jklev+1)+zrich(jlon,jklev))
       if (zriav.gt..25) exit
       jkpbl = jklev
       enddo
      if (jkpbl.eq.0.and.zrich(jlon,1).lt..25) jkpbl = 1
      if (jkpbl.gt.0) pbl(jlon,jlat) = zeta(jlon,jlat,jkpbl)-zorogr(jlon,jlat)
      pbl(jlon,jlat) = min (pbl(jlon,jlat), 2500.)
      enddo

      enddo ! jlat

      return
      end
!###############################################################################################################
   subroutine ccolumn (a, n, prod)

! Computes total cloud cover in a column with an algorithm modified from Geleyn's
! It is assumed that clouds are vertically correlated
! Calculation must proceed from top to bottom

   implicit none
   integer n, kmax(n), kmin(n), k, kma, kmi
   real a(n), prod, sig, sigm, pnum, pden

   if (n.eq.1) then
   prod = (1.0-a(1))
   return
   endif

   a(n) = .999*a(n)  ! to avoid "cloud holes" when cloud cover is constant at bottom

! kmax: index of relative maxima

   kmax = 0
   kmin = 0
   kma = 0
   kmi = 0
   sigm = 1.

   do k = 2, n
   if (a(k).eq.a(k-1)) cycle     ! avoids equal values
   sig = sign(1., a(k)-a(k-1))   ! 1 if a(k)>=a(k-1), else -1
   if (sig*sigm.eq.-1.) then     ! if opposite signs, it is an extreme
   sigm = sig
     if (sig.eq.-1.) then        ! if the second is < 0 ...
     kma = kma+1
     kmax(kma) = k-1             ! ... then k-1 was a max. ...
     elseif (sig.eq.1.) then
     kmi = kmi+1
     kmin(kmi) = k-1             ! ... else it was a min.
     endif
   endif
   enddo

   if (a(n).gt.a(n-1)) then
   kma = kma+1
   kmax(kma) = n                 ! also the bottom level can be a max.
   endif

! product of probabilities of maxima at numerator

   pnum = 1.
   do k = 1, kma
   pnum = pnum*(1.-a(kmax(k)))
   enddo

! product of probabilities of minima at denominator

   pden = 1.
   do k = 1, kmi
   pden = pden*(1.-a(kmin(k)))
   enddo

! ratio of the two

   if (pden.ne.0.) prod = pnum/pden

   return
   end subroutine ccolumn
!###############################################################################################################
      subroutine lift_parcel (qmix, p, tmix, tvlift, nlon, nlat, nlev, jk0, jk1, iwl)

! Computes virtual temperature of a moist air parcel lifted from level jk0
! to a generic level jklev such that: jk0 <= jklev <= jk1
! qmix and tmix are the parcel properties at level jk0
! Results are saved in tvlift
! Liquid water may be removed during lifting (case iwl=0)

      real, dimension(nlon,nlat)      :: qmix, tmix
      real, dimension(nlon,nlat,nlev) :: p, tvlift
      real, parameter :: yliv=2834170.5, yliw=333560.5, ylwv=yliv-yliw, tzer=273.15, pzer=1.e5, ezer=611.0, &
                         rd=287.05, rv=461.51, eps=rd/rv, cpd=1004.6, cpv=1869.46, cw=4186.8,               &
                         ccw1=(cpv-cw)/rv, ccw2=ylwv/tzer/rv-ccw1

      tvlift(:,:,jk0) =  tmix(:,:)*(1.+(1./eps-1.)*qmix(:,:))

      do 100 jlat = 1, nlat
      do 100 jlon = 1, nlon

      jkconl = nlev

! Calculation of parcel entropy before lifting (zs0)

      zesk  = qmix(jlon,jlat)*p(jlon,jlat,jk0)/( eps+(1.-eps)*qmix(jlon,jlat) )
      zctot = (1.-qmix(jlon,jlat))*cpd+qmix(jlon,jlat)*cpv
      zsa   = -(1.-qmix(jlon,jlat))*rd*log((p(jlon,jlat,jk0)-zesk)/pzer)
      zsb   = -qmix(jlon,jlat)*rv*log(zesk/ezer+1.e-18)
      zsc   =  qmix(jlon,jlat)*yliv/tzer
      zs0   = zctot*log(tmix(jlon,jlat)/tzer) + zsa + zsb + zsc

      do jklev = jk0+1, jk1

! Lifted parcel temperature by conservation of entropy in the hypothesis of no condensation

      zesk  =  qmix(jlon,jlat)*p(jlon,jlat,jklev)/( eps+(1.-eps)*qmix(jlon,jlat) )
      zctot =  (1.-qmix(jlon,jlat))*cpd + qmix(jlon,jlat)*cpv
      zsa   = -(1.-qmix(jlon,jlat))*rd*log((p(jlon,jlat,jklev)-zesk)/pzer)
      zsb   = -qmix(jlon,jlat)*rv*log(zesk/ezer+1.e-18)
      zsc   =  qmix(jlon,jlat)*yliv/tzer
      tlift =  tzer*exp((zs0-zsa-zsb-zsc)/zctot)

! Saturation with respect to water

      zt0t = tzer/tlift
      zesk = ezer*exp( -ccw1*log(zt0t) + ccw2*(1.-zt0t) )  ! saturated pressure over water (in pa)
      zqsw = zesk*eps/(p(jlon,jlat,jklev)+zesk*(eps-1.))

      if (qmix(jlon,jlat).le.zqsw) then   ! no condensation occurs
      tvlift(jlon,jlat,jklev) =  tlift*(1.+(1./eps-1.)*qmix(jlon,jlat))
      else
      jkconl = jklev    ! index of the model level just above the lifting condensation level
      go to 9
      endif
      enddo  ! close loop over unsaturated parcel
 9    continue

! Iterative solution in case of condensation with Newton steps

      zt1   = tvlift(jlon,jlat,jkconl-1)
      zq0   = qmix(jlon,jlat)

      do jklev = jkconl, jk1

! Entropy value at zt1 to start the Newton step (zs1)

      zt0t  =  tzer/zt1
      zesk  =  ezer*exp( -ccw1*log(zt0t) + ccw2*(1.-zt0t) )  ! saturated pressure over water (in Pa)
      zqsw  =  zesk*eps/(p(jlon,jlat,jklev)+zesk*(eps-1.))
      zwat  =  zq0 - zqsw
      zctot =  (1.-zq0)*cpd + zqsw*cpv + zwat*cw
      zsa   = -(1.-zq0)*rd*log((p(jlon,jlat,jklev)-zesk)/pzer)
      zsb   = -zqsw*rv*log(zesk/ezer)
      zsc   =  (zqsw*yliv + zwat*yliw)/tzer
      zs1   = -zctot*log(zt0t) + zsa + zsb + zsc

      tlift =  zt1 - 1.   ! first guess of lifted parcel temperature
      icount=0

 10   continue

! Entropy of saturated air with condensed water/ice

      zt0t  =  tzer/tlift
      zesk  =  ezer*exp( -ccw1*log(zt0t) + ccw2*(1.-zt0t) )   ! saturated pressure over water (in pa)
      zqsw  =  zesk*eps/(p(jlon,jlat,jklev)+zesk*(eps-1.))
      zwat  =  zq0 - zqsw
      zctot =  (1.-zq0)*cpd + zqsw*cpv + zwat*cw
      zsa   = -(1.-zq0)*rd*log((p(jlon,jlat,jklev)-zesk)/pzer)
      zsb   = -zqsw*rv*log(zesk/ezer)
      zsc   =  (zqsw*yliv + zwat*yliw)/tzer
      zs    = -zctot*log(zt0t) + zsa + zsb + zsc

! Newton step

      if (abs(tlift-zt1).gt.0.1 .and. icount.le.5) then

      zrds  = (zt1-tlift)/(zs1-zs)
      zt1 = tlift
      zs1 = zs
      tlift = tlift + (zs0-zs)*zrds
      icount = icount+1
      go to 10

      else

      if (iwl.eq.0) then
      tvlift(jlon,jlat,jklev) = tlift*(1.+(1./eps-1.)*zqsw)        ! virtual temperature with no water loading
      elseif (iwl.eq.1) then
      tvlift(jlon,jlat,jklev) = tlift*(1.+(1./eps-1.)*zqsw - zwat) ! virtual temperature with water loading
      endif

      endif

! Removal of liquid water from the lifted parcel

      if (iwl.eq.0) then
      zq0 = zqsw          ! new total water content of lifted parcel
      zs0 = zs0 - zwat*cw*log(tlift/tzer) - zwat*yliw/tzer
      endif

      enddo   ! close loop over model levels where parcel is saturated

 100  continue

      return
      end
!###############################################################################################################
    subroutine pv (fmzh, u, v, w, teta, hx, hy, fmyu, dx, dy, dz, h, a0, p, rd, tv, &
                   fcorio, nlon, nlat, nlev, nlevp1, potvor)

! computes potential vorticity potvor on Z points and half-levels

    use mod_postmoloch, only : gzita

    real, dimension(nlon,nlat,nlev)   :: u, v, teta, p, tv, rom1
    real, dimension(nlon,nlat,nlevp1) :: w, potvor, fmzh
    real, dimension(nlon,nlat)        :: hx, hy, fcorio
    real, dimension(nlat)             :: fmyu

    potvor = 0.
    rom1 = rd*tv/p

    do k = 2, nlev
    zitah = (k-1)*dz
    gzh = gzita (zitah, h, a0)
    do j = 2, nlat
    do i = 1, nlon-1
    zhx = .5*(hx(i,j)+hx(i,j-1))
    zhy = .5*(hy(i,j)+hy(i+1,j))

    tz = .25/dz*( fmzh(i  ,j  ,k)*(teta(i  ,j  ,k)-teta(i  ,j  ,k-1)) + &
                  fmzh(i+1,j,  k)*(teta(i+1,j  ,k)-teta(i+1,j  ,k-1)) + &
                  fmzh(i  ,j-1,k)*(teta(i  ,j-1,k)-teta(i  ,j-1,k-1)) + &
                  fmzh(i+1,j-1,k)*(teta(i+1,j-1,k)-teta(i+1,j-1,k-1)) )

    wz = .125/dz*( fmzh(i  ,j  ,k)*(w(i  ,j  ,k+1)-w(i  ,j  ,k-1)) + &
                   fmzh(i+1,j,  k)*(w(i+1,j  ,k+1)-w(i+1,j  ,k-1)) + &
                   fmzh(i  ,j-1,k)*(w(i  ,j-1,k+1)-w(i  ,j-1,k-1)) + &
                   fmzh(i+1,j-1,k)*(w(i+1,j-1,k+1)-w(i+1,j-1,k-1)) )

    uz = .25/dz*( (fmzh(i,j  ,k)+fmzh(i+1,j  ,k))*(u(i,j  ,k)-u(i,j  ,k-1)) + &
                  (fmzh(i,j-1,k)+fmzh(i+1,j-1,k))*(u(i,j-1,k)-u(i,j-1,k-1)) )

    vz = .25/dz*( (fmzh(i  ,j,k)+fmzh(i  ,j-1,k))*(v(i  ,j,k)-v(i  ,j,k-1)) + &
                  (fmzh(i+1,j,k)+fmzh(i+1,j-1,k))*(v(i+1,j,k)-v(i+1,j,k-1)) )

    tx = .25/dx*( (teta(i+1,j  ,k)-teta(i,j  ,k)+teta(i+1,j  ,k-1)-teta(i,j  ,k-1))*fmyu(j  ) + &
                  (teta(i+1,j-1,k)-teta(i,j-1,k)+teta(i+1,j-1,k-1)-teta(i,j-1,k-1))*fmyu(j-1) ) - gzh*zhx*tz
                  
    ty = .25/dy*( teta(i  ,j,k)-teta(i  ,j-1,k)+teta(i  ,j,k-1)-teta(i  ,j-1,k-1) + &
                  teta(i+1,j,k)-teta(i+1,j-1,k)+teta(i+1,j,k-1)-teta(i+1,j-1,k-1) )   - gzh*zhy*tz

    wx = .5/dx* ((w(i+1,j,k)-w(i,j,k))*fmyu(j) + (w(i+1,j-1,k)-w(i,j-1,k))*fmyu(j-1)) - gzh*zhx*wz
    wy = .5/dy* ( w(i,j,k)-w(i,j-1,k) + w(i+1,j,k)-w(i+1,j-1,k) )                     - gzh*zhy*wz
    vx = .5/dx* ( v(i+1,j,k)-v(i,j,k) + v(i+1,j,k-1)-v(i,j,k-1) )*fmyu(j)             - gzh*zhx*vz
    uy = .5/dy* ( u(i,j,k)-u(i,j-1,k) + u(i,j,k-1)-u(i,j-1,k-1) )                     - gzh*zhy*uz

    z1 = wy - vz
    z2 = uz - wx
    z3 = vx - uy + fcorio(i,j)

    zrom1 = .125*( rom1(i,j,k  )+rom1(i+1,j,k  )+rom1(i,j-1,k  )+rom1(i+1,j-1,k  ) + &
                   rom1(i,j,k-1)+rom1(i+1,j,k-1)+rom1(i,j-1,k-1)+rom1(i+1,j-1,k-1) )

    potvor(i,j,k) = (z1*tx + z2*ty + z3*tz)*zrom1
    enddo
    enddo
    enddo

    return
    end
!###############################################################################################################
subroutine snowfall_height_def (qv, p, t, snowfall_level, zeta, horog, nlon, nlat, nlev, ylwv, cpd, &
 snowfall_height, snowfall_height0)

! Computes snowfall height (limit) that is altitude over sea level
! when ice precipitation predominate water precipitation.
! there are 2 aoutput field: the first is calculate using wet bulb temperature
! (equal +1,3 deg C) only, and the second is corrected by
! the lowest level of ice hydrometeores predomination.

implicit none

real, dimension(nlon,nlat), intent(out) :: snowfall_height, snowfall_height0

integer :: nlon, nlat, nlev
real, dimension(nlon,nlat,nlev), intent(in)  :: qv, t, p, zeta
real, dimension(nlon,nlat), intent(in) :: horog
integer, dimension(nlon,nlat), intent(in) :: snowfall_level
real, intent(in) :: ylwv, cpd

real, parameter :: tzer=273.15, twb_ref = 1.3+tzer

integer :: jlon, jlat, jklev, ktop, kt, ksnowfall
real :: zp, zt, zqv, zqs, zes, zrh, ztwb1, ztwb2, zfr
real, dimension(nlev) :: ztwb

! Searching of index of level lower then 8000 m (above the surface)

 do jlat = 1, nlat
 do jlon = 1, nlon

   ktop = nlev
   do jklev = 1, nlev
     if (zeta(jlon,jlat,jklev)-horog(jlon,jlat) > 8000.) then
       ktop = jklev-1
       exit
     endif
   enddo

   do jklev = 1,ktop

     zp  = p(jlon,jlat,jklev)
     zt  = t(jlon,jlat,jklev)
     zqv = qv(jlon,jlat,jklev)

! Specific humidity of water vapour saturation respect to water also for t<tzer

     call comp_esk(zes, zqs, zt, zp, 3)

! Relative humidity

     zrh = zqv/zqs

! Wet bulb temperature 1st guess by empirical formula     

     ztwb1 = zt + (0.04*(zt-276.) + 1.)*(8.4 + 5.6e-5*(800.e2 - zp))*(zrh - 1.)

! Wet bulb temperature for water by analytical formula

     call comp_esk(zes, zqs, ztwb1, zp, 3)
     ztwb2 = zt - ylwv/cpd*(zqs-zqv)     

! Wet bulb temperature by empirical approximation

     ztwb(jklev) = 0.5*ztwb1 + 0.5*ztwb2

   enddo

! Definition of altitude of reference value of wet bulb temperature
! that is the fist variant of snowfall height

   snowfall_height0(jlon,jlat) = 6000.
   kt = 0
   do jklev = 1,10
     if (ztwb(jklev) < twb_ref) kt = kt+1
   enddo

   if (kt == 10) then

     snowfall_height0(jlon,jlat) = horog(jlon,jlat)

   else

     do jklev = ktop-1, 1, -1
       if (ztwb(jklev) >= twb_ref) then
         zfr = (ztwb(jklev)-twb_ref)/(ztwb(jklev)-ztwb(jklev+1))
         zfr = max( min( zfr, 1.), 0.)
         snowfall_height0(jlon,jlat) = (1.-zfr)*zeta(jlon,jlat,jklev) + zfr*zeta(jlon,jlat,jklev+1)
         exit
       endif
     enddo

   endif

! Definition of second variant of snowfall height using information
! about height of ice hydrometeors predomination

  snowfall_height(jlon,jlat) = snowfall_height0(jlon,jlat)

  ksnowfall=snowfall_level(jlon,jlat)

  if (ksnowfall > 0.and.ksnowfall <= ktop) then
    if (ksnowfall == 1) then
      snowfall_height(jlon,jlat) = horog(jlon,jlat)
    else
      snowfall_height(jlon,jlat) = zeta(jlon,jlat,ksnowfall)
    endif
  endif

! Case of very small precipitation then snowfall_height is very high:

  if (snowfall_height(jlon,jlat)-snowfall_height0(jlon,jlat) > 1000.) snowfall_height(jlon,jlat)=snowfall_height0(jlon,jlat)

 enddo
 enddo

return
end subroutine snowfall_height_def
!###############################################################################################################
      subroutine outgraf(a,nx,ny,alon0r,alat0r,dlonr,dlatr,idd1,idd2,idd3,zf1,zf2,ifl1,ifl2)
      dimension a(nx,ny),za(nx,ny)

      open (11,file='output_graf_md.dat',status='unknown',position='append')

      do j=1,ny
      do i=1,nx
      za(i,j)=a(i,j)*zf1+zf2
      enddo
      enddo

      if (ifl1.eq.1) write (11,*) ifl2
      write (11,*) idd1,idd2,idd3,0,0
      write (11,*) ny,nx,alat0r,alon0r,dlatr,dlonr
      write (11,*) ((za(i,j),i=1,nx),j=1,ny)

      close (11)

      return
      end
! ======================================================================
subroutine free_cross_section(nlon, nlat, nlev, nlev_cross, cross_number, npoint_cross_max, npoint_cross, &
 lonini_cross, latini_cross, lonfin_cross, latfin_cross, h_lev_cross,                                     &
 x0, y0, dlon, dlat, alon0, alat0, dz, h, a0, b0, fz, fzh, rdrcp, valmiss, interp_gauss,                      &
 alont, alatt, alonu, alatu, alonv, alatv,                                                                &
 orogr, p, t, thetae, u, v, w, rh, qcw, qci,                                                              &
 orog_cross, t_cross, theta_cross, thetae_cross, u_cross, v_cross, u_tang_cross, v_norm_cross,            &
 w_cross, rh_cross, qcw_cross, qci_cross)

use mod_postmoloch, only : gzita, bzita

implicit none

! Input:

integer :: nlon, nlat, nlev, nlev_cross, cross_number, npoint_cross_max, interp_gauss
integer, dimension(cross_number) :: npoint_cross
real, dimension(cross_number) :: lonini_cross, latini_cross, lonfin_cross, latfin_cross
real, dimension(nlev_cross) :: h_lev_cross
real :: x0, y0, dlat, dlon, alon0, alat0, dz, h, a0, b0, rdrcp, valmiss
real, dimension(nlev) :: fz
real, dimension(nlev+1) :: fzh
real, dimension(nlon,nlat) :: alont, alatt, alonu, alatu, alonv, alatv, orogr
real, dimension(nlon,nlat,nlev) :: p, t, thetae, u, v, rh, qcw, qci
real, dimension(nlon,nlat,nlev+1) :: w

! Output :

real, dimension(npoint_cross_max,cross_number) :: orog_cross
real, dimension(npoint_cross_max,nlev_cross,cross_number) :: t_cross, theta_cross, thetae_cross, &
 u_cross, v_cross, u_tang_cross, v_norm_cross, w_cross, rh_cross, qcw_cross, qci_cross

! Work :

real, dimension(npoint_cross_max,cross_number) :: lon_cross, lat_cross, lon_cross_rot, lat_cross_rot
real, dimension(nlon,nlat,nlev) :: theta
real, dimension(npoint_cross_max,nlev,cross_number) :: t_cross_lev, theta_cross_lev, thetae_cross_lev, &
 u_cross_lev, v_cross_lev, w_cross_lev, rh_cross_lev, qcw_cross_lev, qci_cross_lev
real, dimension(nlon) :: lon_rot, lon_u_rot
real, dimension(nlat) :: lat_rot, lat_v_rot
real, dimension(nlev) :: h_lev, h_lev_1
integer, dimension(nlev_cross) :: iv
integer, dimension(nlon,nlat) :: mask_field
real, dimension(nlon,nlat) :: work2d
integer, parameter :: npar=9
real, dimension(nlon,nlat,nlev,npar) :: work4d
real, dimension(npoint_cross_max,nlev,npar) :: work3d

integer :: icross, i, jlon, jlat, jklev, np, nlev_inp, im1, ip1, nsmooth=20
real :: zdlon, zdlat, zita, &
 semi_radius=6., gauss_par=1., wei_smooth=0.5, alfa_hinterp=0.6, &
 alfa_vinterp=0.5, ex_bott=0.6, ex_top=0.6, pi, angle, zcos, zsin
character(len=20) :: file_out_work

  pi = abs(acos(-1.))

  mask_field(:,:) = 1

! Definiton of coordinates (lon and lat) of free cross-section

  do icross = 1,cross_number
    write (file_out_work, '(a,I2.2,a)') 'cross_lon_lat_',icross,'.txt'
    open (11, file=file_out_work, form="formatted", status="unknown")
    zdlon=(lonfin_cross(icross)-lonini_cross(icross))/float(npoint_cross(icross)-1)
    zdlat=(latfin_cross(icross)-latini_cross(icross))/float(npoint_cross(icross)-1)
    do i=1,npoint_cross(icross)
      lon_cross(i,icross)=lonini_cross(icross)+float(i-1)*zdlon
      lat_cross(i,icross)=latini_cross(icross)+float(i-1)*zdlat
      write (11, '(2f11.6)') lon_cross(i,icross),lat_cross(i,icross)
    enddo
    close (11)
  enddo

! Definition of cross-section parameter on model grid and levels

  do jklev = 1, nlev
  do jlat = 1, nlat
  do jlon = 1, nlon

    theta(jlon,jlat,jklev) = t(jlon,jlat,jklev)*(1.e5/p(jlon,jlat,jklev))**rdrcp

  enddo
  enddo
  enddo

! Horizontal interpolation to cross-section points

  do jlon=1,nlon
    lon_rot(jlon)=alon0+float(jlon-1)*dlon
    lon_u_rot(jlon)=alon0+float(jlon-1)*dlon+dlon*0.5
  enddo
  do jlat=1,nlat
    lat_rot(jlat)=alat0+float(jlat-1)*dlat
    lat_v_rot(jlat)=alat0+float(jlat-1)*dlat-dlat*0.5
  enddo

  do icross=1,cross_number

    np = npoint_cross(icross)

    call anti_rot_grid(x0, y0, lon_cross(1:np,icross), lat_cross(1:np,icross), &
 lon_cross_rot(1:np,icross), lat_cross_rot(1:np,icross), np, 1)

    if (maxval(lon_cross_rot(1:np,icross)) > maxval(lon_rot)) then
      print*
      print *,'Vertical cross-section coordinate is outside if model domain:'
      print *,'Cross-section number ',icross
      print *,'in rotation coordinates: max longitude of cross-section: ', &
 maxval(lon_cross_rot(1:np,icross)),' of domain: ',maxval(lon_rot)
      stop
    endif

    if (minval(lon_cross_rot(1:np,icross)) < minval(lon_rot)) then
      print*
      print *,'Vertical cross-section coordinate is outside if model domain:'
      print *,'Cross-section number ',icross
      print *,'in rotation coordinates: min longitude of cross-section: ', &
 minval(lon_cross_rot(1:np,icross)),' of domain: ',minval(lon_rot)
      stop
    endif

    if (maxval(lat_cross_rot(1:np,icross)) > maxval(lat_rot)) then
      print*
      print *,'Vertical cross-section coordinate is outside if model domain:'
      print *,'Cross-section number ',icross
      print *,'in rotation coordinates: max latitude of cross-section: ', &
 maxval(lat_cross_rot(1:np,icross)),' of domain: ',maxval(lat_rot)
      stop
    endif

    if (minval(lat_cross_rot(1:np,icross)) < minval(lat_rot)) then
      print*
      print *,'Vertical cross-section coordinate is outside if model domain:'
      print *,'Cross-section number ',icross
      print *,'in rotation coordinates: min latitude of cross-section: ', &
 minval(lat_cross_rot(1:np,icross)),' of domain: ',minval(lat_rot)
      stop
    endif

  enddo

  print*
  if (interp_gauss == 1 ) then

    print *,'Horizontal interpolation for free cross-section with Gaussian algorithm,'
    print *,'semi-radius=',semi_radius,' km'

! Interpolation with gaussian 2-d function

    work4d(1:nlon,1:nlat,1,1) = orogr (1:nlon,1:nlat)
    do icross=1,cross_number
      np = npoint_cross(icross)
      call interp_gauss_input_grid (work4d(1:nlon,1:nlat,1:1,1:1), work3d(1:np,1:1,1:1), &
 alont, alatt, lon_cross(1:np,icross), lat_cross(1:np,icross),                           &
 nlon, nlat, np, 1, 1, mask_field, semi_radius, gauss_par, valmiss)
      orog_cross(1:np,icross) = work3d(1:np,1,1)
    enddo

    work4d(1:nlon,1:nlat,1:nlev,1) = t     (1:nlon,1:nlat,1:nlev)
    work4d(1:nlon,1:nlat,1:nlev,2) = theta (1:nlon,1:nlat,1:nlev)
    work4d(1:nlon,1:nlat,1:nlev,3) = thetae(1:nlon,1:nlat,1:nlev)
    work4d(1:nlon,1:nlat,1:nlev,4) = w     (1:nlon,1:nlat,1:nlev)
    work4d(1:nlon,1:nlat,1:nlev,5) = rh    (1:nlon,1:nlat,1:nlev)
    work4d(1:nlon,1:nlat,1:nlev,6) = qcw   (1:nlon,1:nlat,1:nlev)
    work4d(1:nlon,1:nlat,1:nlev,7) = qci   (1:nlon,1:nlat,1:nlev)
    do icross=1,cross_number
      np = npoint_cross(icross)
      call interp_gauss_input_grid (work4d(1:nlon,1:nlat,1:nlev,1:7), work3d(1:np,1:nlev,1:7), &
 alont, alatt, lon_cross(1:np,icross), lat_cross(1:np,icross), &
 nlon, nlat, np, nlev, 7, mask_field, semi_radius, gauss_par, valmiss)
      t_cross_lev(1:np,1:nlev,icross)      = work3d(1:np,1:nlev,1)
      theta_cross_lev(1:np,1:nlev,icross)  = work3d(1:np,1:nlev,2)
      thetae_cross_lev(1:np,1:nlev,icross) = work3d(1:np,1:nlev,3)
      w_cross_lev(1:np,1:nlev,icross)      = work3d(1:np,1:nlev,4)
      rh_cross_lev(1:np,1:nlev,icross)     = work3d(1:np,1:nlev,5)
      qcw_cross_lev(1:np,1:nlev,icross)    = work3d(1:np,1:nlev,6)
      qci_cross_lev(1:np,1:nlev,icross)    = work3d(1:np,1:nlev,7)
    enddo

    work4d(1:nlon,1:nlat,1:nlev,1) = u(1:nlon,1:nlat,1:nlev)
    do icross=1,cross_number
      np = npoint_cross(icross)
      call interp_gauss_input_grid (work4d(1:nlon,1:nlat,1:nlev,1:1), work3d(1:np,1:nlev,1:1), &
 alonu, alatu, lon_cross(1:np,icross), lat_cross(1:np,icross),                                 &
 nlon, nlat, np, nlev, 1, mask_field, semi_radius, gauss_par, valmiss)
      u_cross_lev(1:np,1:nlev,icross) = work3d(1:np,1:nlev,1)
    enddo

    work4d(1:nlon,1:nlat,1:nlev,1) = v(1:nlon,1:nlat,1:nlev)
    do icross=1,cross_number
      np = npoint_cross(icross)
      call interp_gauss_input_grid (work4d(1:nlon,1:nlat,1:nlev,1:1), work3d(1:np,1:nlev,1:1), &
 alonv, alatv, lon_cross(1:np,icross), lat_cross(1:np,icross),                                 &
 nlon, nlat, np, nlev, 1, mask_field, semi_radius, gauss_par, valmiss)
      v_cross_lev(1:np,1:nlev,icross) = work3d(1:np,1:nlev,1)
    enddo

  else

! Interpolation with spline algorithm

    print *,'Horizontal interpolation for free cross-sections with spline algorithm and smoothed input fields,'
    print *,'parameters of spatial smoothing: weight ',wei_smooth,' numer ',nsmooth

! All 2d fields are smoothed before horizontal interpolation

    call smooth(orogr, work2d, nlon, nlat, wei_smooth, nsmooth)
    do icross=1,cross_number
      np = npoint_cross(icross)
      call interp_spline_2d(work2d, nlon, nlat, lon_rot, lat_rot, &
 lon_cross_rot(1:np,icross), lat_cross_rot(1:np,icross),          &
 np, orog_cross(1:np,icross), alfa_hinterp)
    enddo

    do jklev = 1,nlev

      call smooth(t(:,:,jklev), work2d, nlon, nlat, wei_smooth, nsmooth)
      do icross=1,cross_number
        np = npoint_cross(icross)
        call interp_spline_2d(work2d, nlon, nlat, lon_rot, lat_rot, &
 lon_cross_rot(1:np,icross), lat_cross_rot(1:np,icross),            &
 np, t_cross_lev(1:np,jklev,icross), alfa_hinterp)
      enddo

      call smooth(theta(:,:,jklev), work2d, nlon, nlat, wei_smooth, nsmooth)
      do icross=1,cross_number
        np = npoint_cross(icross)
        call interp_spline_2d(work2d, nlon, nlat, lon_rot, lat_rot, &
 lon_cross_rot(1:np,icross), lat_cross_rot(1:np,icross),            &
 np, theta_cross_lev(1:np,jklev,icross), alfa_hinterp)
      enddo

      call smooth(thetae(:,:,jklev), work2d, nlon, nlat, wei_smooth, nsmooth)
      do icross=1,cross_number
        np = npoint_cross(icross)
        call interp_spline_2d(work2d, nlon, nlat, lon_rot, lat_rot, &
 lon_cross_rot(1:np,icross), lat_cross_rot(1:np,icross),            &
 np, thetae_cross_lev(1:np,jklev,icross), alfa_hinterp)
      enddo

      call smooth(u(:,:,jklev), work2d, nlon, nlat, wei_smooth, nsmooth)
      do icross=1,cross_number
        np = npoint_cross(icross)
        call interp_spline_2d(work2d, nlon, nlat, lon_u_rot, lat_rot, &
 lon_cross_rot(1:np,icross), lat_cross_rot(1:np,icross),              &
 np, u_cross_lev(1:np,jklev,icross), alfa_hinterp)
      enddo

      call smooth(v(:,:,jklev), work2d, nlon, nlat, wei_smooth, nsmooth)
      do icross=1,cross_number
        np = npoint_cross(icross)
        call interp_spline_2d(work2d, nlon, nlat, lon_rot, lat_v_rot, &
 lon_cross_rot(1:np,icross), lat_cross_rot(1:np,icross),              &
 np, v_cross_lev(1:np,jklev,icross), alfa_hinterp)
      enddo

      call smooth(w(:,:,jklev), work2d, nlon, nlat, wei_smooth, nsmooth)
      do icross=1,cross_number
        np = npoint_cross(icross)
        call interp_spline_2d(work2d, nlon, nlat, lon_rot, lat_rot,  &
 lon_cross_rot(1:np,icross), lat_cross_rot(1:np,icross),             &
 np, w_cross_lev(1:np,jklev,icross), alfa_hinterp)
      enddo

      call smooth(rh(:,:,jklev), work2d, nlon, nlat, wei_smooth, nsmooth)
      do icross=1,cross_number
        np = npoint_cross(icross)
        call interp_spline_2d(work2d, nlon, nlat, lon_rot, lat_rot, &
 lon_cross_rot(1:np,icross), lat_cross_rot(1:np,icross),            &
 np, rh_cross_lev(1:np,jklev,icross), alfa_hinterp)
      enddo

      call smooth(qcw(:,:,jklev), work2d, nlon, nlat, wei_smooth, nsmooth)
      do icross=1,cross_number
        np = npoint_cross(icross)
        call interp_spline_2d(work2d, nlon, nlat, lon_rot, lat_rot, &
 lon_cross_rot(1:np,icross), lat_cross_rot(1:np,icross),            &
 np, qcw_cross_lev(1:np,jklev,icross), alfa_hinterp)
      enddo

      call smooth(qci(:,:,jklev), work2d, nlon, nlat, wei_smooth, nsmooth)
      do icross=1,cross_number
        np = npoint_cross(icross)
        call interp_spline_2d(work2d, nlon, nlat, lon_rot, lat_rot, &
 lon_cross_rot(1:np,icross), lat_cross_rot(1:np,icross),            &
 np, qci_cross_lev(1:np,jklev,icross), alfa_hinterp)
      enddo

    enddo

  endif ! Horizontal iterpolation type
  print*

! Vertical interpolation to cross-section z-levels

  do icross=1,cross_number

    np = npoint_cross(icross)

    do i=1,np

! Definition of vertical axis in input: geometric height of model levels

      do jklev=1,nlev
        zita = (jklev-1)*dz + .5*dz
        h_lev(jklev)= -h*bzita(zita,h,b0)*log(fz(jklev)) + orog_cross(i,icross)*gzita(zita,h,a0)
      enddo

      do jklev=1,nlev
        zita = (jklev-1)*dz
        h_lev_1(jklev)= -h*bzita(zita,h,b0)*log(fzh(jklev)) + orog_cross(i,icross)*gzita(zita,h,a0)
      enddo

      call near(h_lev_cross, nlev_cross, h_lev, nlev, iv)

      call interp_spline_1d(t_cross(i,1:nlev_cross,icross), h_lev_cross, nlev_cross, &
 t_cross_lev(i,1:nlev,icross), h_lev, nlev, iv, alfa_vinterp, ex_bott, ex_top)

      call interp_spline_1d(theta_cross(i,1:nlev_cross,icross), h_lev_cross, nlev_cross, &
 theta_cross_lev(i,1:nlev,icross), h_lev, nlev, iv, alfa_vinterp, ex_bott, ex_top)

      call interp_spline_1d(thetae_cross(i,1:nlev_cross,icross), h_lev_cross, nlev_cross, &
 thetae_cross_lev(i,1:nlev,icross), h_lev, nlev, iv, alfa_vinterp, ex_bott, ex_top)

      call interp_spline_1d(u_cross(i,1:nlev_cross,icross), h_lev_cross, nlev_cross, &
 u_cross_lev(i,1:nlev,icross), h_lev, nlev, iv, alfa_vinterp, ex_bott, ex_top)

      call interp_spline_1d(v_cross(i,1:nlev_cross,icross), h_lev_cross, nlev_cross, &
 v_cross_lev(i,1:nlev,icross), h_lev, nlev, iv, alfa_vinterp, ex_bott, ex_top)

      call interp_spline_1d(w_cross(i,1:nlev_cross,icross), h_lev_cross, nlev_cross, &
 w_cross_lev(i,1:nlev,icross), h_lev_1, nlev, iv, alfa_vinterp, ex_bott, ex_top)

      call interp_spline_1d(rh_cross(i,1:nlev_cross,icross), h_lev_cross, nlev_cross, &
 rh_cross_lev(i,1:nlev,icross), h_lev, nlev, iv, alfa_vinterp, ex_bott, ex_top)

      call interp_spline_1d(qcw_cross(i,1:nlev_cross,icross), h_lev_cross, nlev_cross, &
 qcw_cross_lev(i,1:nlev,icross), h_lev, nlev, iv, alfa_vinterp, ex_bott, ex_top)

      call interp_spline_1d(qci_cross(i,1:nlev_cross,icross), h_lev_cross, nlev_cross, &
 qci_cross_lev(i,1:nlev,icross), h_lev, nlev, iv, alfa_vinterp, ex_bott, ex_top)

    enddo

! Definition of tangential (u_tang) and normal (v_norm) wind components
! relatively to cross-section "vector"

    do i=1,np

! Definition of angle between model grid x-axis and cross-section vector

      im1=max(i-1,1)
      ip1=min(i+1,np)

      zdlat=lat_cross_rot(ip1,icross)-lat_cross_rot(im1,icross)
      zdlon=lon_cross_rot(ip1,icross)-lon_cross_rot(im1,icross)

      if (zdlat > 0.) then
        if (zdlon > 0.) then
          angle=atan(zdlat/zdlon)
        else
          angle=pi*0.5
        endif
      elseif (zdlat < 0.) then
        if (zdlon > 0.) then
          angle=2.*pi-atan(-zdlat/zdlon)
        else
          angle=pi*1.5
        endif
      else
        angle=0.
      endif

      zcos=cos(angle)
      zsin=sin(angle)

! Definition of wind component tangential to cross-section vector

      u_tang_cross(i,:,icross)=u_cross(i,:,icross)*zcos+v_cross(i,:,icross)*zsin

! Definition of wind component normal to cross-section vector

      v_norm_cross(i,:,icross)=v_cross(i,:,icross)*zcos-u_cross(i,:,icross)*zsin

    enddo

!  Reset of negative values of humidity and water due to interpolation

    do jklev=1,nlev_cross
      do i=1,np
        rh_cross(i,jklev,icross)=max( rh_cross(i,jklev,icross), 0. )
        qcw_cross(i,jklev,icross)=max( qcw_cross(i,jklev,icross), 0. )
        qci_cross(i,jklev,icross)=max( qci_cross(i,jklev,icross), 0. )
        if (h_lev_cross(jklev) < orog_cross(i,icross)) then
          u_tang_cross(i,jklev,icross)=0.
          v_norm_cross(i,jklev,icross)=0.
          w_cross(i,jklev,icross)=0.
        endif
      enddo
    enddo

  enddo

return
end

!-----------------------------------------------------------------------------

!----------- Subroutines for coding ouput data in grib2 format ---------------

SUBROUTINE WRITE_GRIB2_HORIZONTAL_GRID(NLON,NLAT,NLEV,NLEVG,NLEV_SNOW_PROD,   &
 MHFR,IDATE0,IPERIOD_INP,IPERIOD_ACCUM,IDATEC,PDR,PLEVO,                      &
 TTR,ALONT,ALATT,FMASK,ZLEV,                                                  &
 ZPH,ZT,ZU,ZV,ZQ,ZRH,ZW,ZTHETAE,ZCLWI,ZPV,                                    &
 TG,QG,QGMAX,QGMIN,                                                           &
 ZOROGR,ZPZ0,PRECTOT,PRECSOLID,U10,V10,T2,TD2,Q2,Q2REL,ZTSKIN,HFLUX,QFLUX,        &
 ZTS,ZTHES,ZUS,ZVS,ZQS,CLOUDT,CLOUDH,CLOUDM,CLOUDL,SNOW,RUNOFF,               &
 CAPE,QSKIN,ALBEDO,EMIS1,EMIS2,CSWFL,CLWFL,CHFLUX,CQFLUX,                     &
 T2MIN,T2MAX,WS10MAX,ZLIFT,SEAICE_F,SEAICE_TH,LAPSE_RATE_1,LAPSE_RATE_2, IWV, &
 CIN, PBL, GUST, ZLEV0,                                                       &
 SNOW_LEV_DEPTH_PROD, SNOW_T_PROD, SNOW_AGE_PROD, SNOW_DENS_PROD,             &
 SNOWFALL_HEIGHT, SNOWFALL_HEIGHT0, PS)

! Procedure that prepares fields of post-processed model data (on horizontal grid of various
! level types) for coding in grib2 format

! Uses module grib2_coding_data in write_grib2_data.F90

USE GRIB2_CODING_DATA
USE MOD_POSTMOLOCH, ONLY : NJUMP, NLEVPO, NZLEV, ISFLAG, IPFLAG, IZFLAG, BZITA, GZITA
USE DATA_OBSERV_AUX

IMPLICIT NONE

INTEGER :: NLON, NLAT, NLEV, NLEVG, NLEV_SNOW_PROD, MHFR
INTEGER, DIMENSION(5) :: IDATE0, IDATEC
INTEGER, DIMENSION(3) :: IPERIOD_INP, IPERIOD_ACCUM
REAL, DIMENSION(100) :: PDR
REAL, DIMENSION(NLEVPO) :: PLEVO
REAL :: TTR
REAL, DIMENSION(NLON,NLAT) :: ALONT, ALATT, FMASK
REAL, DIMENSION(NZLEV) :: ZLEV
REAL, DIMENSION(NLON,NLAT,NLEVPO) :: ZPH, ZT, ZU, ZV, ZQ, ZRH, ZW, ZTHETAE, ZCLWI, ZPV
REAL, DIMENSION(NLON,NLAT,NLEVG) :: TG, QG, QGMAX, QGMIN
REAL, DIMENSION(NLON,NLAT,NLEV_SNOW_PROD) :: SNOW_LEV_DEPTH_PROD, SNOW_T_PROD, SNOW_AGE_PROD, SNOW_DENS_PROD
REAL, DIMENSION(NLON,NLAT) :: ZOROGR, ZPZ0, PRECTOT, PRECSOLID, U10, V10, T2, TD2, Q2, Q2REL, ZTSKIN,  &
 HFLUX, QFLUX, ZTS, ZTHES, ZUS, ZVS, ZQS, CLOUDT, CLOUDH, CLOUDM, CLOUDL, SNOW, RUNOFF,            &
 CAPE, QSKIN, ALBEDO, EMIS1, EMIS2, CSWFL, CLWFL, CHFLUX, CQFLUX,                                  &
 T2MIN, T2MAX, WS10MAX, ZLIFT, SEAICE_F, SEAICE_TH, LAPSE_RATE_1, LAPSE_RATE_2, IWV, CIN, PBL, GUST, ZLEV0, &
 SNOWFALL_HEIGHT, SNOWFALL_HEIGHT0, PS

INTEGER :: IFIELD, I, J, K, FLAG_QG, FLAG_TG, NPAR_ISOBAR_LEV, NPAR_ZETTA_LEV, NPAR_SURF
REAL, DIMENSION(NLEVG) :: LEV_SOIL
REAL :: PERIOD_SEC, H, A0, B0, DZ, ZITA, AK, BK
INTEGER, DIMENSION(1) :: ARRAY_SHAPE_1D
INTEGER, DIMENSION(2) :: ARRAY_SHAPE_2D

  ARRAY_SHAPE_2D = SHAPE(IPFLAG)
  NPAR_ISOBAR_LEV = ARRAY_SHAPE_2D(1)

!!  ARRAY_SHAPE_2D = SHAPE(IZFLAG)
!!  NPAR_ZETTA_LEV = ARRAY_SHAPE_2D(1)
  NPAR_ZETTA_LEV = 1

  ARRAY_SHAPE_1D = SHAPE(ISFLAG)
  NPAR_SURF = ARRAY_SHAPE_1D(1)

! Name of ouput file

  WRITE (OUTPUT_FILE_NAME,'(A,I4.4,3I2.2,A,I3.3,2I2.2,A)') "moloch_",IDATE0(1:4),"_",IPERIOD_INP(1:3),".grib2"

! Depth of soil layers

  LEV_SOIL(1:NLEVG)=PDR(6:5+NLEVG)

  FLAG_QG=0
  IF (ANY(ISFLAG(09:12) == 1)) FLAG_QG=1
  FLAG_TG=0
  IF (ANY(ISFLAG(14:17) == 1)) FLAG_TG=1

! Statistical period in second

  PERIOD_SEC = IPERIOD_ACCUM(1)*86400 + IPERIOD_ACCUM(2)*3600 + IPERIOD_ACCUM(3)*60
  PERIOD_SEC=MAX(PERIOD_SEC, 1.)

! Cleaning of possible previous allocations

  DO IFIELD=1,NFIELD0
    IF (ALLOCATED(DATA(IFIELD) % FIELD) ) DEALLOCATE (DATA(IFIELD) % FIELD)
    IF (ALLOCATED(DATA(IFIELD) % VERT_COORD_PAR) ) DEALLOCATE (DATA(IFIELD) % VERT_COORD_PAR)
    DATA(IFIELD) % GRIB2_DESCRIPT(:) = IVALMISS
    DATA(IFIELD) % GRIB2_DESCRIPT(11:14) = 0
  ENDDO

! Description of data for grib2 coding (INTEGER, DIMENSION(50) :: GRIB2_DESCRIPT=IVALMISS
! see module_write_grib2_data.F90)

! Content of GRIB2_DESCRPT array - see bottom of file

! Model index: 1 - Bolam, 2 - Moloch, 3 - Globo

  DATA(1) % GRIB2_DESCRIPT(1) = 2

! Status and type of data

  DATA(1) % GRIB2_DESCRIPT(6) = 0  ! operational products
  DATA(1) % GRIB2_DESCRIPT(7) = 1  ! forecast
  DATA(1) % GRIB2_DESCRIPT(30) = 0 ! bit-map absent

! Model vertical coordinate parameters

  DATA(1) % N_VERT_COORD_PAR = NLEV*2
  ALLOCATE (DATA(1) % VERT_COORD_PAR (DATA(1) % N_VERT_COORD_PAR) )
  H  = PDR(40)
  B0 = PDR(42)
  A0 = PDR(43)
  DZ = H/FLOAT(NLEV)
  DO K=1,NLEV
    ZITA = (K-1)*DZ + DZ*0.5
    AK = -BZITA(ZITA,H,B0)*H*LOG(1.-ZITA/H)
    BK = GZITA(ZITA,H,A0)
    DATA(1) % VERT_COORD_PAR(K) = AK
    DATA(1) % VERT_COORD_PAR(K+NLEV) = BK ! zeta(i,j,k) = horog(i,j)*bk+ak
  ENDDO

! Grid parameters

  DATA(1) % GRIB2_DESCRIPT(2) = 1 ! grid template index - horizontal grid
  DATA(1) % NX = NLON
  DATA(1) % NY = NLAT
  DATA(1) % X0 = PDR(39)
  DATA(1) % Y0 = PDR(38)
  DATA(1) % DX = PDR(2)
  DATA(1) % DY = PDR(1)
  DATA(1) % X00 = PDR(5)
  DATA(1) % Y00 = PDR(4) + DATA(1) % DY*0.5/FLOAT(MHFR)

! Initial date and time

  DATA(1) % IDATE0(1:5) = IDATE0(1:5)

! Date and time of current forecast term

  DATA(1) % IDATEC(1:5) = IDATEC(1:5)

! Definition of time unit, forecast and period length in defined time unit

  DATA(1) % GRIB2_DESCRIPT(4) = 1 ! Time unit for Forecast time step and Statistical period : 0 minute, 1 hour, 2 day
  DATA(1) % IFORECAST = IPERIOD_INP(1)*24+IPERIOD_INP(2)   ! Forecast time step
  DATA(1) % IPERIOD = IPERIOD_ACCUM(1)*24+IPERIOD_ACCUM(2) ! Statistical period

! Product template index

  DATA(1) % GRIB2_DESCRIPT(3) = 0 ! instant

! Level/layer parameters (for forecast satellite products parameters of satellite platform and sensor)

  DATA(1) % GRIB2_DESCRIPT(10:15) = IVALMISS

! Calculation of number of data fields and copy of common parameters to all data description arrays

   IFIELD=0

! Data at isobaric levels

  DO K=1,NLEVPO
    DO I=1,NPAR_ISOBAR_LEV ! Number of parameters defined at isobaric levels
      IF (IPFLAG(I,K) == 1) THEN
        IFIELD = IFIELD+1
        IF (I == 3) IFIELD = IFIELD+1 ! 1 field for second wind component
      ENDIF
    ENDDO
  ENDDO

! Data at zetta levels

  DO K=1,NZLEV
!!    DO I=1,NPAR_ZETTA_LEV ! Number of parameters defined at zetta (geometric hight) levels
      IF (IZFLAG(K) == 1) THEN
        IFIELD = IFIELD+1
      ENDIF
!!    ENDDO
  ENDDO

! Data at the surface

  DO I=1,NPAR_SURF ! Number of parameters defined at the surface
    IF (ISFLAG(I) == 1) THEN
      IFIELD = IFIELD+1
      IF (I == 5.OR.I == 22) IFIELD = IFIELD+1 ! 1 field for second wind component
      IF (I == 59) IFIELD = IFIELD-1+NLEV_SNOW_PROD*4 ! snow cover products: 4 variables at NLEV_SNOW_PROD levels
    ENDIF
  ENDDO

! For QG and TG at all soil level

  IF (FLAG_QG == 1) THEN
    DO I=9,12
      IF (ISFLAG(I) == 1) IFIELD = IFIELD-1
    ENDDO
    DO K=1,NLEVG
      IFIELD = IFIELD+1
    ENDDO
  ENDIF

  IF (FLAG_TG == 1) THEN
    DO I=14,17
      IF (ISFLAG(I) == 1) IFIELD = IFIELD-1
    ENDDO
    DO K=1,NLEVG
      IFIELD = IFIELD+1
    ENDDO
  ENDIF

  NFIELD = IFIELD

  DO IFIELD=1,NFIELD
    ALLOCATE(DATA(IFIELD) % FIELD(NLON,NLAT))
  ENDDO

  DO IFIELD=2,NFIELD

! Grid parameters

    DATA(IFIELD) % NX = DATA(1) % NX
    DATA(IFIELD) % NY = DATA(1) % NY
    DATA(IFIELD) % X0 = DATA(1) % X0
    DATA(IFIELD) % Y0 = DATA(1) % Y0
    DATA(IFIELD) % DX = DATA(1) % DX
    DATA(IFIELD) % DY = DATA(1) % DY
    DATA(IFIELD) % X00 = DATA(1) % X00
    DATA(IFIELD) % Y00 = DATA(1) % Y00

! Initial date and time

    DATA(IFIELD) % IDATE0(:) = DATA(1) % IDATE0(:)

! Length of forecast term and period (for statistical data fields)

    DATA(IFIELD) % IFORECAST = DATA(1) % IFORECAST
    DATA(IFIELD) % IPERIOD = DATA(1) % IPERIOD

! Date and time for current forecast term

    DATA(IFIELD) % IDATEC(:) = DATA(1) % IDATEC(:)

! Description array

    DATA(IFIELD) % GRIB2_DESCRIPT(:) = DATA(1) % GRIB2_DESCRIPT(:)

! Parameters of model vertical coordinate

    DATA(IFIELD) % N_VERT_COORD_PAR = DATA(1) % N_VERT_COORD_PAR
    ALLOCATE (DATA(IFIELD) % VERT_COORD_PAR (DATA(IFIELD) % N_VERT_COORD_PAR) )
    DATA(IFIELD) % VERT_COORD_PAR (:) = DATA(1) % VERT_COORD_PAR (:)

  ENDDO

! Definition of data fields and data parameters in grib2 terms

  IFIELD=0

! Data on isobaric levels

  DO K=1,NLEVPO
    DO I=1,NPAR_ISOBAR_LEV ! Number of parameters defined at isobaric levels

      IF (IPFLAG(I,K)==1) THEN

        IFIELD = IFIELD+1

        DATA(IFIELD) % GRIB2_DESCRIPT(10) = 100 ! Isobaric surface  (Pa)
        DATA(IFIELD) % GRIB2_DESCRIPT(11) = INT(PLEVO(K))
        DATA(IFIELD) % GRIB2_DESCRIPT(12) = 0
        DATA(IFIELD) % GRIB2_DESCRIPT(13:14) = IVALMISS
        DATA(IFIELD) % GRIB2_DESCRIPT(20) = 0   ! Discipline: Meteorological products

        IF (I == 1) THEN ! Geopotential
          DATA(IFIELD) % GRIB2_DESCRIPT(21) = 3 ! Category: Mass
          DATA(IFIELD) % GRIB2_DESCRIPT(22) = 5 ! Parameter: Geopotential height (gpm)
          DATA(IFIELD) % FIELD(:,:) = ZPH(:,:,K)
        ENDIF

        IF (I == 2) THEN ! Temperature
          DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0 ! Category: Temperature
          DATA(IFIELD) % GRIB2_DESCRIPT(22) = 0 ! Parameter: Temperature (K)
          DATA(IFIELD) % FIELD(:,:) = ZT(:,:,K)+TTR
        ENDIF

        IF (I == 3) THEN ! Wind components
          DATA(IFIELD) % GRIB2_DESCRIPT(21) = 2 ! Category: Momentum
          DATA(IFIELD) % GRIB2_DESCRIPT(22) = 2 ! Parameter: u-component of wind (m s-1)
          DATA(IFIELD) % FIELD(:,:) = ZU(:,:,K)
          IFIELD = IFIELD+1
          DATA(IFIELD) % GRIB2_DESCRIPT(10) = 100 ! Isobaric surface  (Pa)
          DATA(IFIELD) % GRIB2_DESCRIPT(11) = NINT(PLEVO(K))
          DATA(IFIELD) % GRIB2_DESCRIPT(12) = 0
          DATA(IFIELD) % GRIB2_DESCRIPT(13:14) = IVALMISS
          DATA(IFIELD) % GRIB2_DESCRIPT(20) = 0 ! Discipline: Meteorological products
          DATA(IFIELD) % GRIB2_DESCRIPT(21) = 2 ! Category: Momentum
          DATA(IFIELD) % GRIB2_DESCRIPT(22) = 3 ! Parameter: v-component of wind (m s-1)
          DATA(IFIELD) % FIELD(:,:) = ZV(:,:,K)
        ENDIF

        IF (I == 4) THEN ! Specific humidity
          DATA(IFIELD) % GRIB2_DESCRIPT(21) = 1 ! Category: Moisture
          DATA(IFIELD) % GRIB2_DESCRIPT(22) = 0 ! Parameter: Specific humidity (kg kg-1)
          DATA(IFIELD) % FIELD(:,:) = ZQ(:,:,K)
        ENDIF

        IF (I == 5) THEN ! Relative humidity
          DATA(IFIELD) % GRIB2_DESCRIPT(21) = 1 ! Category: Moisture
          DATA(IFIELD) % GRIB2_DESCRIPT(22) = 1 ! Parameter: Relative humidity (%)
          DATA(IFIELD) % FIELD(:,:) = ZRH(:,:,K)
        ENDIF

        IF (I == 6) THEN ! W (Vertical velocity)
          DATA(IFIELD) % GRIB2_DESCRIPT(21) = 2 ! Category: Momentum
          DATA(IFIELD) % GRIB2_DESCRIPT(22) = 9 ! Parameter: Vertical velocity [geometric] (m s-1)
          DATA(IFIELD) % FIELD(:,:) = ZW(:,:,K)
        ENDIF

        IF (I == 7) THEN ! THETAE
          DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0 ! Category: Temperature
          DATA(IFIELD) % GRIB2_DESCRIPT(22) = 3 ! Parameter: Pseudo-adiabatic potential temperature or equivalent potential temperature (K)
          DATA(IFIELD) % FIELD(:,:) = ZTHETAE(:,:,K)
        ENDIF

        IF (I == 8) THEN ! Total cloud water content or Potential vorticity
          IF (PLEVO(K) > 300.E2) THEN
            DATA(IFIELD) % GRIB2_DESCRIPT(21) = 6  ! Category: Cloud
            DATA(IFIELD) % GRIB2_DESCRIPT(22) = 17 ! Parameter: Total condensate (kg kg-1)
            DATA(IFIELD) % FIELD(:,:) = ZCLWI(:,:,K)
          ELSE
            DATA(IFIELD) % GRIB2_DESCRIPT(21) = 2  ! Category: Momentum
            DATA(IFIELD) % GRIB2_DESCRIPT(22) = 14 ! Parameter: Potential vorticity (K m2 kg-1 s-1)
            DATA(IFIELD) % FIELD(:,:) = ZPV(:,:,K)*1.E-6 ! conversion from PV units
          ENDIF
        ENDIF

      ENDIF

    ENDDO
  ENDDO

! Data at zeta (geometric height) levels

!  DO K=1,NZLEV
!
!    IF (ITFLAG(K) == 1) THEN
!
!      IFIELD = IFIELD+1
!
!      DATA(IFIELD) % GRIB2_DESCRIPT(10) = 103 ! Specified height level above ground (m)
!      DATA(IFIELD) % GRIB2_DESCRIPT(11) = NINT(ZLEV(K))
!      DATA(IFIELD) % GRIB2_DESCRIPT(12) = 0
!      DATA(IFIELD) % GRIB2_DESCRIPT(13:14) = IVALMISS
!      DATA(IFIELD) % GRIB2_DESCRIPT(20) =
!      DATA(IFIELD) % GRIB2_DESCRIPT(21) =
!      DATA(IFIELD) % GRIB2_DESCRIPT(22) =
!      DATA(IFIELD) % FIELD(:,:) =
!
!    ENDIF
!
!  ENDDO

! Data at the model soil levels

  IF (FLAG_QG == 1) THEN

    DO K=1,NLEVG

      IFIELD = IFIELD+1

      DATA(IFIELD) % GRIB2_DESCRIPT(10) = 106 ! Depth below land surface (m)
      DATA(IFIELD) % GRIB2_DESCRIPT(11) = NINT(LEV_SOIL(K)*1.E3) ! First scaled value of level (layer)
      DATA(IFIELD) % GRIB2_DESCRIPT(12) = 3   ! Scale of first value of level (layer)
      DATA(IFIELD) % GRIB2_DESCRIPT(13:14) = IVALMISS
      DATA(IFIELD) % GRIB2_DESCRIPT(20) = 2   ! Discipline: Land surface products
      DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0   ! Category: Vegetation/Biomass
      DATA(IFIELD) % GRIB2_DESCRIPT(22) = 9   ! Parameter: Volumetric soil moisture content (Proportion)
!      DATA(IFIELD) % FIELD(:,:) = QG(:,:,K)
      DATA(IFIELD) % FIELD(:,:) = (QG(:,:,K)-QGMIN(:,:,K))/(QGMAX(:,:,K)-QGMIN(:,:,K)) ! relative value

    ENDDO

  ENDIF

  IF (FLAG_TG == 1) THEN

    DO K=1,NLEVG

      IFIELD = IFIELD+1

      DATA(IFIELD) % GRIB2_DESCRIPT(10) = 106 ! Depth below land surface (m)
      DATA(IFIELD) % GRIB2_DESCRIPT(11) = NINT(LEV_SOIL(K)*1.E3) ! First scaled value of level (layer)
      DATA(IFIELD) % GRIB2_DESCRIPT(12) = 3   ! Scale of first value of level (layer)
      DATA(IFIELD) % GRIB2_DESCRIPT(13:14) = IVALMISS
      DATA(IFIELD) % GRIB2_DESCRIPT(20) = 2   ! Discipline: Land surface products
      DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0   ! Category: Vegetation/Biomass
      DATA(IFIELD) % GRIB2_DESCRIPT(22) = 2   ! Parameter: Soil temperature (K)
      DATA(IFIELD) % FIELD(:,:) = TG(:,:,K)+TTR

    ENDDO

  ENDIF

! Data at the surface

  DO I=1,NPAR_SURF ! Number of parameters defined at the surface

    IF (I >=  9.AND.I <= 12) CYCLE ! QG at model soil levels
    IF (I >= 14.AND.I <= 17) CYCLE ! TG at model soil levels

    IF (ISFLAG(I) == 1) THEN

      IFIELD = IFIELD+1

      DATA(IFIELD) % GRIB2_DESCRIPT(10) = 1         ! Ground or water surface
      DATA(IFIELD) % GRIB2_DESCRIPT(11) = IVALMISS  ! First scaled value of level (layer)
      DATA(IFIELD) % GRIB2_DESCRIPT(12) = IVALMISS  ! Scale of first value of level (layer)
      DATA(IFIELD) % GRIB2_DESCRIPT(13) = IVALMISS  ! Second scaled value of level (layer)
      DATA(IFIELD) % GRIB2_DESCRIPT(14) = IVALMISS  ! Scale of second value of level (layer)
      DATA(IFIELD) % GRIB2_DESCRIPT(20) = 0         ! Discipline: Meteorological products

      IF (I == 1) THEN ! Orography
        DATA(IFIELD) % GRIB2_DESCRIPT(20) = 2 ! Discipline: Land surface products
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0 ! Category: Vegetation/Biomass
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 7 ! Parameter: Model terrain height (m)
        DATA(IFIELD) % FIELD(:,:) = ZOROGR(:,:)
      ENDIF

      IF (I == 2) THEN ! M.s.l. pressure
        DATA(IFIELD) % GRIB2_DESCRIPT(10) = 101 ! Mean sea level
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 3   ! Category: Mass
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 1   ! Parameter: Pressure reduced to MSL (Pa)
        DATA(IFIELD) % FIELD(:,:) = ZPZ0(:,:)*1.E2
      ENDIF

      IF (I == 3) THEN ! Total precipitation
        DATA(IFIELD) % GRIB2_DESCRIPT(3) = 8  ! Product template: Statistical
        DATA(IFIELD) % GRIB2_DESCRIPT(5) = 1  ! Statistical elaboration type: accumulation
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 1 ! Category: Moisture
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 8 ! Parameter: Total precipitation (kg m-2)
        DATA(IFIELD) % FIELD(:,:) = PRECTOT(:,:)
      ENDIF

      IF (I == 4) THEN ! Solid phase precipitation
        DATA(IFIELD) % GRIB2_DESCRIPT(3) = 8   ! Product template: Statistical
        DATA(IFIELD) % GRIB2_DESCRIPT(5) = 1   ! Statistical elaboration type: accumulation
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 1  ! Category: Moisture
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 29 ! Parameter: Total snowfall (m)
        DATA(IFIELD) % FIELD(:,:) = PRECSOLID(:,:)*1.E-3
      ENDIF

      IF (I == 5) THEN ! Wind at 10 m
        DATA(IFIELD) % GRIB2_DESCRIPT(10) = 103 ! Specified height level above ground (m)
        DATA(IFIELD) % GRIB2_DESCRIPT(11) = 10  ! First scaled value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(12) = 0   ! Scale of first value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 2   ! Category: Momentum
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 2   ! Parameter: u-component of wind (m s-1)
        DATA(IFIELD) % FIELD(:,:) = U10(:,:)
        IFIELD = IFIELD+1
        DATA(IFIELD) % GRIB2_DESCRIPT(10) = 103 ! Specified height level above ground (m)
        DATA(IFIELD) % GRIB2_DESCRIPT(11) = 10  ! First scaled value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(12) = 0   ! Scale of first value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(13:14) = IVALMISS
        DATA(IFIELD) % GRIB2_DESCRIPT(20) = 0   ! Discipline: Meteorological products
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 2   ! Category: Momentum
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 3   ! Parameter: v-component of wind (m s-1)
        DATA(IFIELD) % FIELD(:,:) = V10(:,:)
      ENDIF

      IF (I == 6) THEN ! Temper. at 2m
        DATA(IFIELD) % GRIB2_DESCRIPT(10) = 103 ! Specified height level above ground (m)
        DATA(IFIELD) % GRIB2_DESCRIPT(11) = 2   ! First scaled value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(12) = 0   ! Scale of first value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0   ! Category: Temperature
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 0   ! Parameter: Temperature (K)
        DATA(IFIELD) % FIELD(:,:) = T2(:,:)+TTR
      ENDIF

      IF (I == 7) THEN ! Specific hum. at 2 m
        DATA(IFIELD) % GRIB2_DESCRIPT(10) = 103 ! Specified height level above ground (m)
        DATA(IFIELD) % GRIB2_DESCRIPT(11) = 2   ! First scaled value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(12) = 0   ! Scale of first value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 1   ! Category: Moisture
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 0   ! Parameter: Specific humidity (kg kg-1)
        DATA(IFIELD) % FIELD(:,:) = Q2(:,:)
      ENDIF

      IF (I == 8) THEN ! Relative hum. at 2 m
        DATA(IFIELD) % GRIB2_DESCRIPT(10) = 103 ! Specified height level above ground (m)
        DATA(IFIELD) % GRIB2_DESCRIPT(11) = 2   ! First scaled value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(12) = 0   ! Scale of first value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 1   ! Category: Moisture
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 1   ! Parameter: Relative humidity (%)
        DATA(IFIELD) % FIELD(:,:) = Q2REL(:,:)
      ENDIF

!      IF (I == 9) THEN ! Ground water lev. 1
!        K=1
!        DATA(IFIELD) % GRIB2_DESCRIPT(10) = 106 ! Depth below land surface (m)
!        DATA(IFIELD) % GRIB2_DESCRIPT(11) = NINT(LEV_SOIL(K)*1.E3) ! First scaled value of level (layer)
!        DATA(IFIELD) % GRIB2_DESCRIPT(12) = 3   ! Scale of first value of level (layer)
!        DATA(IFIELD) % GRIB2_DESCRIPT(20) = 2   ! Discipline: Land surface products
!        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0   ! Category: Vegetation/Biomass
!        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 9   ! Parameter: Volumetric soil moisture content (Proportion)
!        DATA(IFIELD) % FIELD(:,:) = QG(:,:,K)
!      ENDIF
!
!      IF (I == 10) THEN ! Ground water lev. 3
!        K=3
!        DATA(IFIELD) % GRIB2_DESCRIPT(10) = 106 ! Depth below land surface (m)
!        DATA(IFIELD) % GRIB2_DESCRIPT(11) = NINT(LEV_SOIL(K)*1.E3) ! First scaled value of level (layer)
!        DATA(IFIELD) % GRIB2_DESCRIPT(12) = 3   ! Scale of first value of level (layer)
!        DATA(IFIELD) % GRIB2_DESCRIPT(20) = 2   ! Discipline: Land surface products
!        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0   ! Category: Vegetation/Biomass
!        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 9   ! Parameter: Volumetric soil moisture content (Proportion)
!        DATA(IFIELD) % FIELD(:,:) = QG(:,:,K)
!      ENDIF
!
!      IF (I == 11) THEN ! Ground water lev. 5
!        K=5
!        DATA(IFIELD) % GRIB2_DESCRIPT(10) = 106 ! Depth below land surface (m)
!        DATA(IFIELD) % GRIB2_DESCRIPT(11) = NINT(LEV_SOIL(K)*1.E3) ! First scaled value of level (layer)
!        DATA(IFIELD) % GRIB2_DESCRIPT(12) = 3   ! Scale of first value of level (layer)
!        DATA(IFIELD) % GRIB2_DESCRIPT(20) = 2   ! Discipline: Land surface products
!        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0   ! Category: Vegetation/Biomass
!        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 9   ! Parameter: Volumetric soil moisture content (Proportion)
!        DATA(IFIELD) % FIELD(:,:) = QG(:,:,K)
!      ENDIF
!
!      IF (I == 12) THEN ! Ground water lev. 7
!        K=7
!        DATA(IFIELD) % GRIB2_DESCRIPT(10) = 106 ! Depth below land surface (m)
!        DATA(IFIELD) % GRIB2_DESCRIPT(11) = NINT(LEV_SOIL(K)*1.E3) ! First scaled value of level (layer)
!        DATA(IFIELD) % GRIB2_DESCRIPT(12) = 3   ! Scale of first value of level (layer)
!        DATA(IFIELD) % GRIB2_DESCRIPT(20) = 2   ! Discipline: Land surface products
!        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0   ! Category: Vegetation/Biomass
!        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 9   ! Parameter: Volumetric soil moisture content (Proportion)
!        DATA(IFIELD) % FIELD(:,:) = QG(:,:,K)
!      ENDIF

      IF (I == 13) THEN ! Skin temperature
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0 ! Category: Temperature
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 0 ! Parameter: Temperature (K)
        DATA(IFIELD) % FIELD(:,:) = ZTSKIN(:,:)+TTR
      ENDIF

!      IF (I == 14) THEN ! Ground temp. lev. 1
!        K=1
!        DATA(IFIELD) % GRIB2_DESCRIPT(10) = 106 ! Depth below land surface (m)
!        DATA(IFIELD) % GRIB2_DESCRIPT(11) = NINT(LEV_SOIL(K)*1.E3) ! First scaled value of level (layer)
!        DATA(IFIELD) % GRIB2_DESCRIPT(12) = 3   ! Scale of first value of level (layer)
!        DATA(IFIELD) % GRIB2_DESCRIPT(20) = 2   ! Discipline: Land surface products
!        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0   ! Category: Vegetation/Biomass
!        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 2   ! Parameter: Soil temperature (K)
!        DATA(IFIELD) % FIELD(:,:) = TG(:,:,K)+TTR
!      ENDIF
!
!      IF (I == 15) THEN ! Ground temp. lev. 3
!        K=3
!        DATA(IFIELD) % GRIB2_DESCRIPT(10) = 106 ! Depth below land surface (m)
!        DATA(IFIELD) % GRIB2_DESCRIPT(11) = NINT(LEV_SOIL(K)*1.E3) ! First scaled value of level (layer)
!        DATA(IFIELD) % GRIB2_DESCRIPT(12) = 3   ! Scale of first value of level (layer)
!        DATA(IFIELD) % GRIB2_DESCRIPT(20) = 2   ! Discipline: Land surface products
!        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0   ! Category: Vegetation/Biomass
!        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 2   ! Parameter: Soil temperature (K)
!        DATA(IFIELD) % FIELD(:,:) = TG(:,:,K)+TTR
!      ENDIF
!
!      IF (I == 16) THEN ! Ground temp. lev. 5
!        K=5
!        DATA(IFIELD) % GRIB2_DESCRIPT(10) = 106 ! Depth below land surface (m)
!        DATA(IFIELD) % GRIB2_DESCRIPT(11) = NINT(LEV_SOIL(K)*1.E3) ! First scaled value of level (layer)
!        DATA(IFIELD) % GRIB2_DESCRIPT(12) = 3   ! Scale of first value of level (layer)
!        DATA(IFIELD) % GRIB2_DESCRIPT(20) = 2   ! Discipline: Land surface products
!        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0   ! Category: Vegetation/Biomass
!        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 2   ! Parameter: Soil temperature (K)
!        DATA(IFIELD) % FIELD(:,:) = TG(:,:,K)+TTR
!      ENDIF
!
!      IF (I == 17) THEN ! Ground temp. lev. 7
!        K=7
!        DATA(IFIELD) % GRIB2_DESCRIPT(10) = 106 ! Depth below land surface (m)
!        DATA(IFIELD) % GRIB2_DESCRIPT(11) = NINT(LEV_SOIL(K)*1.E3) ! First scaled value of level (layer)
!        DATA(IFIELD) % GRIB2_DESCRIPT(12) = 3   ! Scale of first value of level (layer)
!        DATA(IFIELD) % GRIB2_DESCRIPT(20) = 2   ! Discipline: Land surface products
!        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0   ! Category: Vegetation/Biomass
!        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 2   ! Parameter: Soil temperature (K)
!        DATA(IFIELD) % FIELD(:,:) = TG(:,:,K)+TTR
!      ENDIF

      IF (I == 18) THEN ! Flux of sensible heat
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0  ! Category: Temperature
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 11 ! Parameter: Sensible heat net flux (W m-2)
        DATA(IFIELD) % FIELD(:,:) = HFLUX(:,:)
      ENDIF

      IF (I == 19) THEN ! Flux of latent heat
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0  ! Category: Temperature
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 10 ! Parameter: Latent heat net flux (W m-2)
        DATA(IFIELD) % FIELD(:,:) = QFLUX(:,:)
      ENDIF

      IF (I == 20) THEN ! Temper. at lowest atm. level
        DATA(IFIELD) % GRIB2_DESCRIPT(10) = 105 ! Hybrid level
        DATA(IFIELD) % GRIB2_DESCRIPT(11) = 1   ! First scaled value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(12) = 0   ! Scale of first value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0   ! Category: Temperature
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 0   ! Parameter: Temperature (K)
        DATA(IFIELD) % FIELD(:,:) = ZTS(:,:)+TTR
      ENDIF

      IF (I == 21) THEN ! Equiv. pot. temper. at lowest atm. level
        DATA(IFIELD) % GRIB2_DESCRIPT(10) = 105 ! Hybrid level
        DATA(IFIELD) % GRIB2_DESCRIPT(11) = 1   ! First scaled value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(12) = 0   ! Scale of first value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0   ! Category: Temperature
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 3   ! Parameter: Pseudo-adiabatic potential temperature or equivalent pot. temp. (K)
        DATA(IFIELD) % FIELD(:,:) = ZTHES(:,:)+TTR
      ENDIF

      IF (I == 22) THEN ! Wind at lowest atm. level
        DATA(IFIELD) % GRIB2_DESCRIPT(10) = 105 ! Hybrid level
        DATA(IFIELD) % GRIB2_DESCRIPT(11) = 1   ! First scaled value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(12) = 0   ! Scale of first value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 2   ! Category: Momentum
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 2   ! Parameter: u-component of wind (m s-1)
        DATA(IFIELD) % FIELD(:,:) = ZUS(:,:)
        IFIELD = IFIELD+1
        DATA(IFIELD) % GRIB2_DESCRIPT(10) = 105 ! Hybrid level
        DATA(IFIELD) % GRIB2_DESCRIPT(11) = 1   ! First scaled value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(12) = 0   ! Scale of first value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(13:14) = IVALMISS
        DATA(IFIELD) % GRIB2_DESCRIPT(20) = 0   ! Discipline: Meteorological products
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 2   ! Category: Momentum
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 3   ! Parameter: v-component of wind (m s-1)
        DATA(IFIELD) % FIELD(:,:) = ZVS(:,:)
      ENDIF

      IF (I == 23) THEN ! Spec. hum. at lowest atm. level
        DATA(IFIELD) % GRIB2_DESCRIPT(10) = 105 ! Hybrid level
        DATA(IFIELD) % GRIB2_DESCRIPT(11) = 1   ! First scaled value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(12) = 0   ! Scale of first value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 1   ! Category: Moisture
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 0   ! Parameter: Specific humidity (kg kg-1)
        DATA(IFIELD) % FIELD(:,:) = ZQS(:,:)
      ENDIF

      IF (I == 24) THEN ! Total cloud cover post proc. or model
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 6 ! Category: Cloud
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 1 ! Parameter: Total cloud cover (%)
        DATA(IFIELD) % FIELD(:,:) = CLOUDT(:,:)*1.E2
      ENDIF

      IF (I == 25) THEN ! High cloud cover
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 6 ! Category: Cloud
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 5 ! Parameter: High cloud cover (%)
        DATA(IFIELD) % FIELD(:,:) = CLOUDH(:,:)*1.E2
      ENDIF

      IF (I == 26) THEN ! Middle cloud cover
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 6 ! Category: Cloud
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 4 ! Parameter: Medium cloud cover (%)
        DATA(IFIELD) % FIELD(:,:) = CLOUDM(:,:)*1.E2
      ENDIF

      IF (I == 27) THEN ! Low cloud cover
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 6 ! Category: Cloud
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 3 ! Parameter: Low cloud cover (%)
        DATA(IFIELD) % FIELD(:,:) = CLOUDL(:,:)*1.E2
      ENDIF

      IF (I == 28) THEN ! Snow cover
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 1  ! Category: Moisture
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 13 ! Parameter: Water equivalent of accumulated snow depth (kg m-2)
        DATA(IFIELD) % FIELD(:,:) = SNOW(:,:)*1.E3
      ENDIF

      IF (I == 29) THEN ! Runoff
        DATA(IFIELD) % GRIB2_DESCRIPT(20) = 2 ! Discipline: Land surface products
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0 ! Category: Vegetation/Biomass
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 5 ! Parameter: Water runoff (kg m-2
        DATA(IFIELD) % FIELD(:,:) = RUNOFF(:,:)
      ENDIF

      IF (I == 30) THEN ! CAPE
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 7 ! Category: Thermodynamic Stability Indices
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 6 ! Parameter: Convective available potential energy (J kg-1)
        DATA(IFIELD) % FIELD(:,:) = CAPE(:,:)
      ENDIF

      IF (I == 31) THEN ! Skin Specific hum.
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 1 ! Category: Moisture
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 0 ! Parameter: Specific humidity (kg kg-1)
        DATA(IFIELD) % FIELD(:,:) = QSKIN(:,:)
      ENDIF

      IF (I == 32) THEN ! Temper. at 0.5 m above surf.
        DATA(IFIELD) % GRIB2_DESCRIPT(10) = 103 ! Specified height level above ground (m)
        DATA(IFIELD) % GRIB2_DESCRIPT(11) = 50  ! First scaled value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(12) = 2   ! Scale of first value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0   ! Category: Temperature
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 0   ! Parameter: Temperature (K)
        DATA(IFIELD) % FIELD(:,:) = T05(:,:)+TTR
      ENDIF

      IF (I == 33) THEN ! Relative hum. at 0.5 m above surf.
        DATA(IFIELD) % GRIB2_DESCRIPT(10) = 103 ! Specified height level above ground (m)
        DATA(IFIELD) % GRIB2_DESCRIPT(11) = 50  ! First scaled value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(12) = 2   ! Scale of first value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 1   ! Category: Moisture
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 1   ! Parameter: Relative humidity (%)
        DATA(IFIELD) % FIELD(:,:) = Q05REL(:,:)
      ENDIF

      IF (I == 34) THEN ! Temper. at 0.05 m above surf.
        DATA(IFIELD) % GRIB2_DESCRIPT(10) = 103 ! Specified height level above ground (m)
        DATA(IFIELD) % GRIB2_DESCRIPT(11) = 5   ! First scaled value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(12) = 2   ! Scale of first value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0   ! Category: Temperature
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 0   ! Parameter: Temperature (K)
        DATA(IFIELD) % FIELD(:,:) = T005(:,:)+TTR
      ENDIF

      IF (I == 35) THEN ! Ground temp. at 5 cm depth
        DATA(IFIELD) % GRIB2_DESCRIPT(10) = 106 ! Depth below land surface (m)
        DATA(IFIELD) % GRIB2_DESCRIPT(11) = 5   ! First scaled value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(12) = 2   ! Scale of first value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(20) = 2   ! Discipline: Land surface products
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0   ! Category: Vegetation/Biomass
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 2   ! Parameter: Soil temperature (K)
        DATA(IFIELD) % FIELD(:,:) = TG005(:,:)+TTR
      ENDIF

      IF (I == 36) THEN ! Ground temp. at 10 cm depth
        DATA(IFIELD) % GRIB2_DESCRIPT(10) = 106 ! Depth below land surface (m)
        DATA(IFIELD) % GRIB2_DESCRIPT(11) = 10  ! First scaled value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(12) = 2   ! Scale of first value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(20) = 2   ! Discipline: Land surface products
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0   ! Category: Vegetation/Biomass
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 2   ! Parameter: Soil temperature (K)
        DATA(IFIELD) % FIELD(:,:) = TG010(:,:)+TTR
      ENDIF

      IF (I == 37) THEN ! Ground temp. at 20 cm depth
        DATA(IFIELD) % GRIB2_DESCRIPT(10) = 106 ! Depth below land surface (m)
        DATA(IFIELD) % GRIB2_DESCRIPT(11) = 20  ! First scaled value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(12) = 2   ! Scale of first value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(20) = 2   ! Discipline: Land surface products
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0   ! Category: Vegetation/Biomass
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 2   ! Parameter: Soil temperature (K)
        DATA(IFIELD) % FIELD(:,:) = TG020(:,:)+TTR
      ENDIF

      IF (I == 38) THEN ! Ground temp. at 50 cm depth
        DATA(IFIELD) % GRIB2_DESCRIPT(10) = 106 ! Depth below land surface (m)
        DATA(IFIELD) % GRIB2_DESCRIPT(11) = 50  ! First scaled value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(12) = 2   ! Scale of first value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(20) = 2   ! Discipline: Land surface products
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0   ! Category: Vegetation/Biomass
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 2   ! Parameter: Soil temperature (K)
        DATA(IFIELD) % FIELD(:,:) = TG050(:,:)+TTR
      ENDIF

      IF (I == 39) THEN ! Ground temp. at 100 cm depth
        DATA(IFIELD) % GRIB2_DESCRIPT(10) = 106 ! Depth below land surface (m)
        DATA(IFIELD) % GRIB2_DESCRIPT(11) = 100 ! First scaled value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(12) = 2   ! Scale of first value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(20) = 2   ! Discipline: Land surface products
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0   ! Category: Vegetation/Biomass
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 2   ! Parameter: Soil temperature (K)
        DATA(IFIELD) % FIELD(:,:) = TG100(:,:)+TTR
      ENDIF

      IF (I == 40) THEN ! Albedo
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 19 ! Category: Physical atmospheric properties
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 1  ! Parameter: Albedo (%)
        DATA(IFIELD) % FIELD(:,:) = ALBEDO(:,:)*1.E2
      ENDIF

      IF (I == 41) THEN ! Emissivity (broadband)
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 19  ! Category: Physical atmospheric properties
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 201 ! Parameter: local use Emissivity (%)
        DATA(IFIELD) % FIELD(:,:) = EMIS1(:,:)*1.E2
      ENDIF

      IF (I == 42) THEN ! Emissivity (window)
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 19  ! Category: Physical atmospheric properties
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 202 ! Parameter: local use Emissivity (%)
        DATA(IFIELD) % FIELD(:,:) = EMIS2(:,:)*1.E2
      ENDIF

      IF (I == 43) THEN ! Cumulated short wave radiative flux
        DATA(IFIELD) % GRIB2_DESCRIPT(3) = 8  ! Product template: Statistical
        DATA(IFIELD) % GRIB2_DESCRIPT(5) = 0  ! Statistical elaboration type: average
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 4 ! Category: Short-wave Radiation
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 0 ! Parameter: Net short-wave radiation flux (surface) (W m-2)
        DATA(IFIELD) % FIELD(:,:) = CSWFL(:,:)*1.E3/PERIOD_SEC
      ENDIF

      IF (I == 44) THEN ! Cumulated long wave radiative flux
        DATA(IFIELD) % GRIB2_DESCRIPT(3) = 8  ! Product template: Statistical
        DATA(IFIELD) % GRIB2_DESCRIPT(5) = 0  ! Statistical elaboration type: average
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 5 ! Category: Long-wave Radiation
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 0 ! Parameter: Net long wave radiation flux (surface) (W m-2)
        DATA(IFIELD) % FIELD(:,:) = CLWFL(:,:)*1.E3/PERIOD_SEC
      ENDIF

      IF (I == 45) THEN ! Cumulated sensible heat flux
        DATA(IFIELD) % GRIB2_DESCRIPT(3) = 8   ! Product template: Statistical
        DATA(IFIELD) % GRIB2_DESCRIPT(5) = 0   ! Statistical elaboration type: average
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0  ! Category: Temperature
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 11 ! Parameter: Sensible heat net flux (W m-2)
        DATA(IFIELD) % FIELD(:,:) = CHFLUX(:,:)*1.E3/PERIOD_SEC
      ENDIF

      IF (I == 46) THEN ! Cumulated latent heat flux
        DATA(IFIELD) % GRIB2_DESCRIPT(3) = 8   ! Product template: Statistical
        DATA(IFIELD) % GRIB2_DESCRIPT(5) = 0   ! Statistical elaboration type: average
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0  ! Category: Temperature
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 10 ! Parameter: Latent heat net flux (W m-2)
        DATA(IFIELD) % FIELD(:,:) = CQFLUX(:,:)*1.E3/PERIOD_SEC
      ENDIF

      IF (I == 47) THEN ! T min. 2m
        DATA(IFIELD) % GRIB2_DESCRIPT(3) = 8    ! Product template: Statistical
        DATA(IFIELD) % GRIB2_DESCRIPT(5) = 3    ! Statistical elaboration type: minimum
        DATA(IFIELD) % GRIB2_DESCRIPT(10) = 103 ! Specified height level above ground (m)
        DATA(IFIELD) % GRIB2_DESCRIPT(11) = 2   ! First scaled value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(12) = 0   ! Scale of first value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0   ! Category: Temperature
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 0   ! Parameter: Temperature (K)
        DATA(IFIELD) % FIELD(:,:) = T2MIN(:,:)+TTR
      ENDIF

      IF (I == 48) THEN ! T max. 2m
        DATA(IFIELD) % GRIB2_DESCRIPT(3) = 8    ! Product template: Statistical
        DATA(IFIELD) % GRIB2_DESCRIPT(5) = 2    ! Statistical elaboration type: maximum
        DATA(IFIELD) % GRIB2_DESCRIPT(10) = 103 ! Specified height level above ground (m)
        DATA(IFIELD) % GRIB2_DESCRIPT(11) = 2   ! First scaled value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(12) = 0   ! Scale of first value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0   ! Category: Temperature
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 0   ! Parameter: Temperature (K)
        DATA(IFIELD) % FIELD(:,:) = T2MAX(:,:)+TTR
      ENDIF

      IF (I == 49) THEN ! Wind speed max. 10m
        DATA(IFIELD) % GRIB2_DESCRIPT(3) = 8    ! Product template: Statistical
        DATA(IFIELD) % GRIB2_DESCRIPT(5) = 2    ! Statistical elaboration type: maximum
        DATA(IFIELD) % GRIB2_DESCRIPT(10) = 103 ! Specified height level above ground (m)
        DATA(IFIELD) % GRIB2_DESCRIPT(11) = 10  ! First scaled value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(12) = 0   ! Scale of first value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0   ! Category: Momentum
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 1   ! Parameter: Wind speed (m s-1)
        DATA(IFIELD) % FIELD(:,:) = WS10MAX(:,:)
      ENDIF

      IF (I == 50) THEN ! Lifted index
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 7  ! Category: Thermodynamic Stability Indices
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 10 ! Parameter: Surface lifted index (K)
        DATA(IFIELD) % FIELD(:,:) = ZLIFT(:,:)
      ENDIF

      IF (I == 51) THEN ! Sea ice fraction
        DATA(IFIELD) % GRIB2_DESCRIPT(20) = 10 ! Discipline: Oceanographic products
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 2  ! Category: Ice
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 0  ! Parameter: Ice cover (Proportion)
        DATA(IFIELD) % FIELD(:,:) = SEAICE_F(:,:)
      ENDIF

      IF (I == 52) THEN ! Sea ice thickness (m)
        DATA(IFIELD) % GRIB2_DESCRIPT(20) = 10 ! Discipline: Oceanographic products
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 2  ! Category: Ice
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 1  ! Parameter: Ice thickness (m)
        DATA(IFIELD) % FIELD(:,:) = SEAICE_TH(:,:)
      ENDIF

      IF (I == 53) THEN ! Dew point temperature at 2 m
        DATA(IFIELD) % GRIB2_DESCRIPT(10) = 103 ! Specified height level above ground (m)
        DATA(IFIELD) % GRIB2_DESCRIPT(11) = 2   ! First scaled value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(12) = 0   ! Scale of first value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0   ! Category: Temperature
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 6   ! Parameter: Dew point temperature (K)
        DATA(IFIELD) % FIELD(:,:) = TD2(:,:)+TTR
      ENDIF

      IF (I == 54) THEN ! Integrated Water Vapour (kg/m2)
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 1  ! Category: Moisture
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 64 ! Parameter: Total column integrated water vapour (kg m-2)
        DATA(IFIELD) % FIELD(:,:) = IWV(:,:)
      ENDIF

      IF (I == 55) THEN ! Convective inibition (J/kg)
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 7 ! Category: Thermodynamic Stability Indices
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 7 ! Parameter: Convective inhibition (J kg-1)
        DATA(IFIELD) % FIELD(:,:) = CIN(:,:)
      ENDIF

!      IF (I == 56) THEN ! PBL top height (above surface)
!        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 3  ! Category: Mass
!        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 18 ! Parameter: Planetary boundary layer height (M)
!        DATA(IFIELD) % FIELD(:,:) = PBL(:,:)
!      ENDIF

      IF (I == 57) THEN ! Wind Gust
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 2  ! Category: Momentum
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 22 ! Parameter: Wind speed (gust)  (M sec-1)
        DATA(IFIELD) % FIELD(:,:) = GUST(:,:)
      ENDIF

      IF (I == 58) THEN ! Hight (m) of 0 deg. C isotherm
        DATA(IFIELD) % GRIB2_DESCRIPT(10) = 4   ! Level of 0 C isotherm
        DATA(IFIELD) % GRIB2_DESCRIPT(11) = 0   ! First scaled value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(12) = 0   ! Scale of first value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 3   ! Category: Mass
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 6   ! Parameter: Geometric height (m)
        DATA(IFIELD) % FIELD(:,:) = ZLEV0(:,:)
      ENDIF

      IF (I == 59) THEN ! Snow cover products

        IFIELD = IFIELD-1

        DO K=1,NLEV_SNOW_PROD
! Data at snow product levels: 
! 1 - snow surface, 3 - snow bottom, 2 - level of medium of geometrical snow depth
          DO J=1,4 ! variables: level depth (m), temperature (K), age (days), density (kg m^-3)

            IFIELD = IFIELD+1

            DATA(IFIELD) % GRIB2_DESCRIPT(10) = 114 ! Snow level (Numeric)
            DATA(IFIELD) % GRIB2_DESCRIPT(11) = K
            DATA(IFIELD) % GRIB2_DESCRIPT(12) = 0
            DATA(IFIELD) % GRIB2_DESCRIPT(13:14) = IVALMISS
            DATA(IFIELD) % GRIB2_DESCRIPT(20) = 0   ! Discipline: Meteorological products

            IF (J == 1) THEN ! Geometric depth (from top to bottom)
              DATA(IFIELD) % GRIB2_DESCRIPT(21) =  1 ! Category: Moisture
              DATA(IFIELD) % GRIB2_DESCRIPT(22) = 11 ! Parameter: Snow depth (m)
              DATA(IFIELD) % FIELD(:,:) = SNOW_LEV_DEPTH_PROD(:,:,K)
            ENDIF

            IF (J == 2) THEN ! Temperature
              DATA(IFIELD) % GRIB2_DESCRIPT(21) =  0  ! Category: Temperature
              DATA(IFIELD) % GRIB2_DESCRIPT(22) =  0  ! Parameter: Temperature (K)
              DATA(IFIELD) % FIELD(:,:) = SNOW_T_PROD(:,:,K)
            ENDIF

            IF (J == 3) THEN ! Density (kg m-3)
              DATA(IFIELD) % GRIB2_DESCRIPT(21) =  1 ! Category: Moisture
              DATA(IFIELD) % GRIB2_DESCRIPT(22) = 61 ! Parameter: Snow density (kg m-3)
              DATA(IFIELD) % FIELD(:,:) = SNOW_DENS_PROD(:,:,K)
            ENDIF

            IF (J == 4) THEN ! Age (days)
              DATA(IFIELD) % GRIB2_DESCRIPT(21) =  1 ! Category: Moisture
              DATA(IFIELD) % GRIB2_DESCRIPT(22) = 17 ! Parameter: Snow age (d)
              DATA(IFIELD) % FIELD(:,:) = SNOW_AGE_PROD(:,:,K)
            ENDIF

          ENDDO
        ENDDO

      ENDIF

      IF (I == 60) THEN ! Lapse rate in lower troposphere
        DATA(IFIELD) % GRIB2_DESCRIPT(11) = 0   ! First scaled value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(10) = 1 ! Ground or water surface
        DATA(IFIELD) % GRIB2_DESCRIPT(20) = 0 ! Discipline: Meteorological products
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0 ! Category: Temperature
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 8 ! Parameter: Lapse rate (K m-1)
        DATA(IFIELD) % FIELD(:,:) = LAPSE_RATE_2(:,:) ! Local value
!        IFIELD = IFIELD+1
!        DATA(IFIELD) % GRIB2_DESCRIPT(10) = 1 ! Ground or water surface
!        DATA(IFIELD) % GRIB2_DESCRIPT(20) = 0 ! Discipline: Meteorological products
!        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0 ! Category: Temperature
!        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 8 ! Parameter: Lapse rate (K m-1)
!        DATA(IFIELD) % FIELD(:,:) = LAPSE_RATE_1(:,:) ! Climatological
      ENDIF

      IF (I == 61) THEN ! Height of snow fall limit (m) above sea level
        DATA(IFIELD) % GRIB2_DESCRIPT(10) = 101 ! Mean sea level
        DATA(IFIELD) % GRIB2_DESCRIPT(11) = 0   ! First scaled value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(20) = 0   ! Discipline: Meteorological products
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 1   ! Category: Moisture
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 204 ! Parameter: Height of snow fall limit (m) (ECMWF)
        DATA(IFIELD) % FIELD(:,:) = SNOWFALL_HEIGHT(:,:)
      ENDIF

      IF (I == 62) THEN ! First guess of height of snow fall limit (m) above sea level
        DATA(IFIELD) % GRIB2_DESCRIPT(10) = 101 ! Mean sea level
        DATA(IFIELD) % GRIB2_DESCRIPT(11) = 0   ! First scaled value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(20) = 0   ! Discipline: Meteorological products
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 1   ! Category: Moisture
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 205 ! Parameter: Height of snow fall limit (m) (ISAC)
        DATA(IFIELD) % FIELD(:,:) = SNOWFALL_HEIGHT0(:,:)
      ENDIF

      IF (I == 63) THEN ! Surface Pressure
        DATA(IFIELD) % GRIB2_DESCRIPT(10) = 1 ! Ground or water surface
        DATA(IFIELD) % GRIB2_DESCRIPT(20) = 0 ! Discipline: Meteorological products
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 3 ! Category: Mass
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 0 ! Parameter: Pressure (Pa)
        DATA(IFIELD) % FIELD(:,:) = PS(:,:)
      ENDIF

    ENDIF

  ENDDO

! End of preparation of output data for coding in grib2 format

! Write output data in grib2 format

  CALL WRITE_GRIB2_DATA

RETURN
END

!-----------------------------------------------------------------------------

SUBROUTINE WRITE_GRIB2_CROSS_SPACE(IDATE0, IPERIOD_INP, IPERIOD_ACCUM, IDATEC,                       &
 NCROSS, NLEVEL, NPOINT_MAX, NPOINT, HLEVEL, LONINI_CROSS, LATINI_CROSS, LONFIN_CROSS, LATFIN_CROSS, &
 X0, Y0, OROG_CROSS,                                                                                 &
 T_CROSS, THETA_CROSS, THETAE_CROSS, U_TANG_CROSS, V_NORM_CROSS, W_CROSS,                            &
 RH_CROSS, QCW_CROSS, QCI_CROSS)

! Procedure that prepares fields of post-processed model data (on various space-crossing grids)
! for coding in grib2 format

! Uses module grib2_coding_data in write_grib2_data.F90

USE GRIB2_CODING_DATA

IMPLICIT NONE

INTEGER, DIMENSION(5) :: IDATE0, IDATEC
INTEGER, DIMENSION(3) :: IPERIOD_INP, IPERIOD_ACCUM
INTEGER :: NCROSS, NLEVEL, NPOINT_MAX
INTEGER :: NPARAM = 9 ! (prognostic parameters + orography)
INTEGER, DIMENSION(NCROSS) :: NPOINT
REAL, DIMENSION(NLEVEL) :: HLEVEL
REAL, DIMENSION(NCROSS) :: LONINI_CROSS, LATINI_CROSS, LONFIN_CROSS, LATFIN_CROSS
REAL :: X0, Y0

REAL, DIMENSION(NPOINT_MAX,NCROSS) :: OROG_CROSS

REAL, DIMENSION(NPOINT_MAX,NLEVEL,NCROSS) :: T_CROSS, THETA_CROSS, THETAE_CROSS, &
 U_TANG_CROSS, V_NORM_CROSS, W_CROSS, RH_CROSS, QCW_CROSS, QCI_CROSS

INTEGER :: IFIELD, ICROSS, IPARAM, NX, NZ

! Name of ouput file

  WRITE (OUTPUT_FILE_NAME,'(A,I4.4,3I2.2,A,I3.3,2I2.2,A)') "moloch_cross_",IDATE0(1:4),"_",IPERIOD_INP(1:3),".grib2"

! Cleaning of possible previous allocations

  DO IFIELD=1,NFIELD0
    IF (ALLOCATED(DATA(IFIELD) % FIELD) ) DEALLOCATE (DATA(IFIELD) % FIELD)
    IF (ALLOCATED(DATA(IFIELD) % VERT_COORD_PAR) ) DEALLOCATE (DATA(IFIELD) % VERT_COORD_PAR)
    DATA(IFIELD) % GRIB2_DESCRIPT(:) = IVALMISS
    DATA(IFIELD) % GRIB2_DESCRIPT(11:14) = 0
  ENDDO

! Description of data for grib2 coding (INTEGER, DIMENSION(50) :: GRIB2_DESCRIPT=IVALMISS
! see module_write_grib2_data.F90)

! Content of GRIB2_DESCRPT array see bottom of file

! Model index: 1 - Bolam, 2 - Moloch, 3 - Globo

  DATA(1) % GRIB2_DESCRIPT(1) = 2

! Status and type of data

  DATA(1) % GRIB2_DESCRIPT(6) = 0 ! operational products
  DATA(1) % GRIB2_DESCRIPT(7) = 1 ! forecast
  DATA(1) % GRIB2_DESCRIPT(30) = 0 ! bit-map absent

! Grid parameters

  DATA(1) % GRIB2_DESCRIPT(2) = 1000 ! grid template index - vertical cross-section
  DATA(1) % X0 = X0
  DATA(1) % Y0 = Y0

! Initial date and time

  DATA(1) % IDATE0(1:5) = IDATE0(1:5)

! Date and time of current forecast term

  DATA(1) % IDATEC(1:5) = IDATEC(1:5)

! Definition of time unit, forecast and period length in defined time unit

  DATA(1) % GRIB2_DESCRIPT(4) = 1 ! Time unit for Forecast time step and Statistical period : 0 - minute, 1 -hour, 2 - day
  DATA(1) % IFORECAST = IPERIOD_INP(1)*24+IPERIOD_INP(2) ! Forecast time step
  DATA(1) % IPERIOD = IPERIOD_ACCUM(1)*24+IPERIOD_ACCUM(2) ! Statistical period

! Product template index

  DATA(1) % GRIB2_DESCRIPT(3) = 0 ! instant

! Level/layer parameters (for forecast satellite products parameters of satellite platform and sensor)

  DATA(1) % GRIB2_DESCRIPT(10:15) = IVALMISS

! Calculation of number of data fields and copy of common parameters to all data description arrays

  NFIELD = NCROSS*NPARAM

  DO IFIELD=2,NFIELD

! Grid parameters

    DATA(IFIELD) % X0 = DATA(1) % X0
    DATA(IFIELD) % Y0 = DATA(1) % Y0

! Initial date and time

    DATA(IFIELD) % IDATE0(:) = DATA(1) % IDATE0(:)

! Length of forecast term and period (for statistica data fields)

    DATA(IFIELD) % IFORECAST = DATA(1) % IFORECAST
    DATA(IFIELD) % IPERIOD = DATA(1) % IPERIOD

! Date and time of current forecast term

    DATA(IFIELD) % IDATEC(:) = DATA(1) % IDATEC(:)

! Description array

    DATA(IFIELD) % GRIB2_DESCRIPT(:) = DATA(1) % GRIB2_DESCRIPT(:)

  ENDDO

! Definition of data fields and data parameters in grib2 terms

  IFIELD=0

! Data on space-crossing grids

  DO ICROSS=1,NCROSS

    DO IPARAM=1,NPARAM

      IFIELD = IFIELD+1

! Grid parameters

      NX = NPOINT(ICROSS)
      DATA(IFIELD) % NX = NX
      DATA(IFIELD) % LON_INI = LONINI_CROSS(ICROSS)
      DATA(IFIELD) % LAT_INI = LATINI_CROSS(ICROSS)
      DATA(IFIELD) % LON_FIN = LONFIN_CROSS(ICROSS)
      DATA(IFIELD) % LAT_FIN = LATFIN_CROSS(ICROSS)

      IF (IPARAM /= NPARAM ) THEN ! Prognostic parameters
        NZ = NLEVEL
        DATA(IFIELD) % NY = NZ
        ALLOCATE(DATA(IFIELD) % VERT_COORD_PAR(NZ))
        DATA(IFIELD) % VERT_COORD_PAR(:) = HLEVEL(:)
      ELSE ! Orography
        NZ = 1
        DATA(IFIELD) % NY = NZ
        ALLOCATE(DATA(IFIELD) % VERT_COORD_PAR(NZ))
        DATA(IFIELD) % VERT_COORD_PAR(:) = VALMISS
      ENDIF
      DATA(IFIELD) % IND_VERT_COORD = 102 ! Altitude above mean sea level (m) 3.15.table

! Definition of data fields

      ALLOCATE(DATA(IFIELD) % FIELD(NX,NZ))

      DATA(IFIELD) % GRIB2_DESCRIPT(20) = 0 ! Discipline: Meteorological products

      IF (IPARAM == 1) THEN ! Temperature (K)
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0 ! Category: Temperature
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 0 ! Parameter: Temperature (K)
        DATA(IFIELD) % FIELD(:,:) = T_CROSS(1:NX,1:NZ,ICROSS)
      ENDIF

      IF (IPARAM == 2) THEN ! Potential temperature (K)
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0 ! Category: Temperature
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 2 ! Parameter: Potential temperature (K)
        DATA(IFIELD) % FIELD(:,:) = THETA_CROSS(1:NX,1:NZ,ICROSS)
      ENDIF

      IF (IPARAM == 3) THEN ! Equivalent potential temperature
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0 ! Category: Temperature
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 3 ! Parameter: Potential temperature (K)
        DATA(IFIELD) % FIELD(:,:) = THETAE_CROSS(1:NX,1:NZ,ICROSS)
      ENDIF

      IF (IPARAM == 4) THEN ! U-component of wind
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 2 ! Category: Momentum
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 2 ! Parameter: u-component of wind (m s-1)
        DATA(IFIELD) % FIELD(:,:) = U_TANG_CROSS(1:NX,1:NZ,ICROSS)
      ENDIF

      IF (IPARAM == 5) THEN ! V-component of wind
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 2 ! Category: Momentum
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 3 ! Parameter: v-component of wind (m s-1)
        DATA(IFIELD) % FIELD(:,:) = V_NORM_CROSS(1:NX,1:NZ,ICROSS)
      ENDIF

      IF (IPARAM == 6) THEN ! Vertical velocity
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 2 ! Category: Momentum
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 9 ! Parameter: Vertical velocity [geometric] (m s-1)
        DATA(IFIELD) % FIELD(:,:) = W_CROSS(1:NX,1:NZ,ICROSS)
      ENDIF

      IF (IPARAM == 7) THEN ! Relative humidity
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 1 ! Category: Moisture
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 1 ! Parameter: Relative humidity (%)
        DATA(IFIELD) % FIELD(:,:) = RH_CROSS(1:NX,1:NZ,ICROSS)
      ENDIF

      IF (IPARAM == 8) THEN ! Total water content of cloud
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 6  ! Category: Cloud
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 17 ! Parameter: Total condensate (kg kg-1)
        DATA(IFIELD) % FIELD(:,:) = QCW_CROSS(1:NX,1:NZ,ICROSS) + QCI_CROSS(1:NX,1:NZ,ICROSS)
      ENDIF

      IF (IPARAM == NPARAM) THEN ! Orography height
        DATA(IFIELD) % GRIB2_DESCRIPT(20) = 2 ! Discipline: Land surface products
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0 ! Category: Vegetation/Biomass
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 7 ! Parameter: Model terrain height (m)
        DATA(IFIELD) % FIELD(:,1) = OROG_CROSS(1:NX,ICROSS)
      ENDIF

    ENDDO

  ENDDO

! End of preparation of output data for coding in grib2 format

! Write output data in grib2 format

  CALL WRITE_GRIB2_DATA

RETURN
END

!-----------------------------------------------------------------------------

SUBROUTINE WRITE_GRIB2_CROSS_SPACE_I_J(IDATE0, IPERIOD_INP, IPERIOD_ACCUM, IDATEC,              &
 NLON, NLAT, NLEVEL,  NCROSSY, NCROSSX, HLEVEL, LON_CROSSY, LAT_CROSSY, LON_CROSSX, LAT_CROSSX, &
 X0, Y0, OROG_CROSSY, OROG_CROSSX,                                                              &
 THETA_CROSSY, THETAE_CROSSY, U_CROSSY, V_CROSSY, W_CROSSY, RH_CROSSY, QCW_CROSSY, QCI_CROSSY,  &
 THETA_CROSSX, THETAE_CROSSX, U_CROSSX, V_CROSSX, W_CROSSX, RH_CROSSX, QCW_CROSSX, QCI_CROSSX)

! Procedure that prepares fields of postprocesses model data on space-crossing grids
! for coding in grib2 format.
! Space-crossing grids are oriented along model horizontal grid axes:
! some cross-sections along x-axis and some cross-sections along y-axis.

! Uses module grib2_coding_data in write_grib2_data.F90

USE GRIB2_CODING_DATA

IMPLICIT NONE

INTEGER, DIMENSION(5) :: IDATE0, IDATEC
INTEGER, DIMENSION(3) :: IPERIOD_INP, IPERIOD_ACCUM
INTEGER :: NLON, NLAT, NLEVEL, NCROSSY, NCROSSX
INTEGER :: NPARAM = 9 ! (prognostic parameters + orography)
REAL, DIMENSION(NLEVEL) :: HLEVEL
REAL, DIMENSION(NCROSSY,2) :: LON_CROSSY, LAT_CROSSY
REAL, DIMENSION(NCROSSX,2) :: LON_CROSSX, LAT_CROSSX
REAL :: X0, Y0

REAL, DIMENSION(NLAT,1,NCROSSY) :: OROG_CROSSY
REAL, DIMENSION(NLON,1,NCROSSX) :: OROG_CROSSX

REAL, DIMENSION(NLAT,NLEVEL,NCROSSY) :: THETA_CROSSY, THETAE_CROSSY, U_CROSSY, V_CROSSY, W_CROSSY, &
 RH_CROSSY, QCW_CROSSY, QCI_CROSSY
REAL, DIMENSION(NLON,NLEVEL,NCROSSX) :: THETA_CROSSX, THETAE_CROSSX, U_CROSSX, V_CROSSX, W_CROSSX, &
 RH_CROSSX, QCW_CROSSX, QCI_CROSSX

INTEGER :: IFIELD, IAXIS, NCROSS, ICROSS, IPARAM, NX, NZ

! Name of ouput file

  WRITE (OUTPUT_FILE_NAME,'(A,I4.4,3I2.2,A,I3.3,2I2.2,A)') "moloch_cross_i_j_",IDATE0(1:4),"_",IPERIOD_INP(1:3),".grib2"

! Cleaning of possible previous allocations

  DO IFIELD=1,NFIELD0
    IF (ALLOCATED(DATA(IFIELD) % FIELD) ) DEALLOCATE (DATA(IFIELD) % FIELD)
    IF (ALLOCATED(DATA(IFIELD) % VERT_COORD_PAR) ) DEALLOCATE (DATA(IFIELD) % VERT_COORD_PAR)
    DATA(IFIELD) % GRIB2_DESCRIPT(:) = IVALMISS
  ENDDO

! Description of data for grib2 coding (INTEGER, DIMENSION(50) :: GRIB2_DESCRIPT=IVALMISS, see module_write_grib2_data.F90)

! Content of GRIB2_DESCRPT array see bottom of file

! Model index: 1 - Bolam, 2 - Moloch, 3 - Globo

  DATA(1) % GRIB2_DESCRIPT(1) = 2

! Status and type of data

  DATA(1) % GRIB2_DESCRIPT(6) = 0 ! operational products
  DATA(1) % GRIB2_DESCRIPT(7) = 1 ! forecast
  DATA(1) % GRIB2_DESCRIPT(30) = 0 ! bit-map absent

! Grid parameters

  DATA(1) % GRIB2_DESCRIPT(2) = 1000 ! grid template index - vertical cross-section
  DATA(1) % X0 = X0
  DATA(1) % Y0 = Y0

! Initial date and time

  DATA(1) % IDATE0(1:5) = IDATE0(1:5)

! Date and time for current forecast term

  DATA(1) % IDATEC(1:5) = IDATEC(1:5)

! Definition of time unit, forecast and period length in defined time unit

  DATA(1) % GRIB2_DESCRIPT(4) = 1 ! Time unit for Forecast time step and Statistical period : 0 - minute, 1 -hour, 2 - day
  DATA(1) % IFORECAST = IPERIOD_INP(1)*24+IPERIOD_INP(2)   ! Forecast time step
  DATA(1) % IPERIOD = IPERIOD_ACCUM(1)*24+IPERIOD_ACCUM(2) ! Statistical period

! Product template index

  DATA(1) % GRIB2_DESCRIPT(3) = 0 ! instant

! Level/layer parameters (for forecast satellite products parameters of satellite platform and sensor)

  DATA(1) % GRIB2_DESCRIPT(10:15) = IVALMISS

! Calculation of number of data fields and copy of common parameters to all data description arrays

  NFIELD = NCROSSX*NPARAM + NCROSSY*NPARAM

  DO IFIELD=2,NFIELD

! Grid parameters

    DATA(IFIELD) % X0 = DATA(1) % X0
    DATA(IFIELD) % Y0 = DATA(1) % Y0

! Initial date and time

    DATA(IFIELD) % IDATE0(:) = DATA(1) % IDATE0(:)

! Length of forecast term and period (for statistica data fields)

    DATA(IFIELD) % IFORECAST = DATA(1) % IFORECAST
    DATA(IFIELD) % IPERIOD = DATA(1) % IPERIOD

! Date and time for current forecast term

    DATA(IFIELD) % IDATEC(:) = DATA(1) % IDATEC(:)

! Description array

    DATA(IFIELD) % GRIB2_DESCRIPT(:) = DATA(1) % GRIB2_DESCRIPT(:)

  ENDDO

! Definition of data fields and data parameters in grib2 terms

  IFIELD=0

! Data on space-crossing grids

  DO IAXIS=1,2 ! 1 - meridian (y) axis, 2 - parallel (x) axis

    DO ICROSS=1,NCROSSY

      DO IPARAM=1,NPARAM

        IFIELD = IFIELD+1

! Grid parameters

        IF (IAXIS == 1) THEN ! y-axis
          NX = NLAT
          DATA(IFIELD) % NX = NX
          DATA(IFIELD) % LON_INI = LON_CROSSY(ICROSS,1)
          DATA(IFIELD) % LON_FIN = LON_CROSSY(ICROSS,2)
          DATA(IFIELD) % LAT_INI = LAT_CROSSY(ICROSS,1)
          DATA(IFIELD) % LAT_FIN = LAT_CROSSY(ICROSS,2)
        ELSE ! x-axis
          NX = NLON
          DATA(IFIELD) % NX = NX
          DATA(IFIELD) % LON_INI = LON_CROSSX(ICROSS,1)
          DATA(IFIELD) % LON_FIN = LON_CROSSX(ICROSS,2)
          DATA(IFIELD) % LAT_INI = LAT_CROSSX(ICROSS,1)
          DATA(IFIELD) % LAT_FIN = LAT_CROSSX(ICROSS,2)
        ENDIF

        IF (IPARAM /= NPARAM ) THEN ! Prognostic parameters
          NZ = NLEVEL
          DATA(IFIELD) % NY = NZ
          ALLOCATE(DATA(IFIELD) % VERT_COORD_PAR(NZ))
          DATA(IFIELD) % VERT_COORD_PAR(:) = HLEVEL(:)
        ELSE ! Orography
          NZ = 1
          DATA(IFIELD) % NY = NZ
          ALLOCATE(DATA(IFIELD) % VERT_COORD_PAR(NZ))
          DATA(IFIELD) % VERT_COORD_PAR(:) = VALMISS
        ENDIF
        DATA(IFIELD) % IND_VERT_COORD = 102 ! Altitude above mean sea level (m) 3.15.table

! Definition of data fields

        ALLOCATE(DATA(IFIELD) % FIELD(NX,NZ))

        DATA(IFIELD) % GRIB2_DESCRIPT(20) = 0 ! Discipline: Meteorological products

        IF (IPARAM == 1) THEN ! Potential temperature (K)
          DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0 ! Category: Temperature
          DATA(IFIELD) % GRIB2_DESCRIPT(22) = 2 ! Parameter: Potential temperature (K)
          IF (IAXIS == 1) THEN ! y-axis
            DATA(IFIELD) % FIELD(:,:) = THETA_CROSSY(1:NX,1:NZ,ICROSS)
          ELSE ! x-axis
            DATA(IFIELD) % FIELD(:,:) = THETA_CROSSX(1:NX,1:NZ,ICROSS)
          ENDIF
        ENDIF

        IF (IPARAM == 2) THEN ! Equivalent potential temperature
          DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0 ! Category: Temperature
          DATA(IFIELD) % GRIB2_DESCRIPT(22) = 3 ! Parameter: Potential temperature (K)
          IF (IAXIS == 1) THEN ! y-axis
            DATA(IFIELD) % FIELD(:,:) = THETA_CROSSY(1:NX,1:NZ,ICROSS)
          ELSE ! x-axis
            DATA(IFIELD) % FIELD(:,:) = THETA_CROSSX(1:NX,1:NZ,ICROSS)
          ENDIF
        ENDIF

        IF (IPARAM == 3) THEN ! U-component of wind
          DATA(IFIELD) % GRIB2_DESCRIPT(21) = 2 ! Category: Momentum
          DATA(IFIELD) % GRIB2_DESCRIPT(22) = 2 ! Parameter: u-component of wind (m s-1)
          IF (IAXIS == 1) THEN ! y-axis
            DATA(IFIELD) % FIELD(:,:) = U_CROSSY(1:NX,1:NZ,ICROSS)
          ELSE ! x-axis
            DATA(IFIELD) % FIELD(:,:) = U_CROSSX(1:NX,1:NZ,ICROSS)
          ENDIF
        ENDIF

        IF (IPARAM == 4) THEN ! V-component of wind
          DATA(IFIELD) % GRIB2_DESCRIPT(21) = 2 ! Category: Momentum
          DATA(IFIELD) % GRIB2_DESCRIPT(22) = 3 ! Parameter: v-component of wind (m s-1)
          IF (IAXIS == 1) THEN ! y-axis
            DATA(IFIELD) % FIELD(:,:) = V_CROSSY(1:NX,1:NZ,ICROSS)
          ELSE ! x-axis
            DATA(IFIELD) % FIELD(:,:) = V_CROSSX(1:NX,1:NZ,ICROSS)
          ENDIF
        ENDIF

        IF (IPARAM == 5) THEN ! Vertical velocity
          DATA(IFIELD) % GRIB2_DESCRIPT(21) = 2 ! Category: Momentum
          DATA(IFIELD) % GRIB2_DESCRIPT(22) = 9 ! Parameter: Vertical velocity [geometric] (m s-1)
          IF (IAXIS == 1) THEN ! y-axis
            DATA(IFIELD) % FIELD(:,:) = W_CROSSY(1:NX,1:NZ,ICROSS)
          ELSE ! x-axis
            DATA(IFIELD) % FIELD(:,:) = W_CROSSX(1:NX,1:NZ,ICROSS)
          ENDIF
        ENDIF

        IF (IPARAM == 6) THEN ! Relative humidity
          DATA(IFIELD) % GRIB2_DESCRIPT(21) = 1 ! Category: Moisture
          DATA(IFIELD) % GRIB2_DESCRIPT(22) = 1 ! Parameter: Relative humidity (%)
          IF (IAXIS == 1) THEN ! y-axis
            DATA(IFIELD) % FIELD(:,:) = RH_CROSSY(1:NX,1:NZ,ICROSS)
          ELSE ! x-axis
            DATA(IFIELD) % FIELD(:,:) = RH_CROSSX(1:NX,1:NZ,ICROSS)
          ENDIF
        ENDIF

        IF (IPARAM == 7) THEN ! QCW
          DATA(IFIELD) % GRIB2_DESCRIPT(21) = 1  ! Category: Moisture
          DATA(IFIELD) % GRIB2_DESCRIPT(22) = 83 ! Parameter: Specific cloud liquid water content (kg kg-1)
          IF (IAXIS == 1) THEN ! y-axis
            DATA(IFIELD) % FIELD(:,:) = QCW_CROSSY(1:NX,1:NZ,ICROSS)
          ELSE ! x-axis
            DATA(IFIELD) % FIELD(:,:) = QCW_CROSSX(1:NX,1:NZ,ICROSS)
          ENDIF
        ENDIF

        IF (IPARAM == 8) THEN ! QCI
          DATA(IFIELD) % GRIB2_DESCRIPT(21) = 1  ! Category: Moisture
          DATA(IFIELD) % GRIB2_DESCRIPT(22) = 84 ! Parameter: Specific cloud ice water content (kg kg-1)
          IF (IAXIS == 1) THEN ! y-axis
            DATA(IFIELD) % FIELD(:,:) = QCI_CROSSY(1:NX,1:NZ,ICROSS)
          ELSE ! x-axis
            DATA(IFIELD) % FIELD(:,:) = QCI_CROSSX(1:NX,1:NZ,ICROSS)
          ENDIF
        ENDIF

        IF (IPARAM == NPARAM) THEN ! Orography height
          DATA(IFIELD) % GRIB2_DESCRIPT(20) = 2 ! Discipline: Land surface products
          DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0 ! Category: Vegetation/Biomass
          DATA(IFIELD) % GRIB2_DESCRIPT(22) = 7 ! Parameter: Model terrain height (m)
          IF (IAXIS == 1) THEN ! y-axis
            DATA(IFIELD) % FIELD(:,1) = OROG_CROSSY(1:NX,1,ICROSS)
          ELSE ! x-axis
            DATA(IFIELD) % FIELD(:,1) = OROG_CROSSX(1:NX,1,ICROSS)
          ENDIF
        ENDIF

      ENDDO

    ENDDO

  ENDDO

! End of preparation of output data for coding in grib2 format

! Write output data in grib2 format

  CALL WRITE_GRIB2_DATA

RETURN
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

!-----------------------------------------------------------------------------
! ======================================================================
