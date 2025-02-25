!      program postbolam

! Last update 16/10/2024

! Sett. 2024: Cambiata la scrittura di model_param_constant.bin - piu' grandezze

! Ago. 2024: Nuovo formato mhf (scrittura record non 1d ma 2d)

! Mag. 2023: PS (Surface Pressure) in output in grib2 (per MeteoAosta)

! Corretto il formato del head di ppf

! Cambiato valore di zrc (riduz. nubi) sul mare (da 0.52 a 0.40).
! Calcolo di gust basato sull'altezza del PBL.
! Calcolo opzionale del flusso (trasporto) di vapore integrato verticalmente:
! le compon. del vettore vengono scritte al posto di zus e zvs.

! Ago. 2018: nuovo formato mhf: mfs_atm, mhf_soil e il file con i campi statici (model_param_constant.bin).

! Sett. 2017: modif. calcolo CAPE e lifted index (con nuova subr. lift_parcel_bol)
! Lug. 2017: calcolo rid. nubi includendo il metodo Geleyn
! ma cercando di separare cumuli da stratocumuli, per ridurre solo i secondi (specie sul mare)
! Mag. 2017: definita nuova routine ccloud (Maurizio) per calcolo nubi.
! Definito il numero di Richardson (usato in cloud_post per il calcolo delle nubi)
! al primo semi-livello sopra il suolo come zri (definito in vtsurflux).
! Introdotta matrice richs -nuova static stability per ridurre cloud fraction (si usa
! Richardson calc. con Durran&Klemp per la parte umida).
! Apr. 2017: versione con calcolo diverse quantita' alla tropopausa e
! al livello di vento massimo (output previsto solo in grib2).
! Richiede cambiamenti nel file postbolam.inp, inclusa scrittura dei liv. p in questo file.
! Per mantenere la compatibilita' con programmi NCARG offline (boldis.def) senza
! cambiare dimensioni, introdotta matrice di flag ipflag20(npar x 20) (al posto di
! ipflag(npar x nlevpo0)) - ipflag20 viene scritta nel file bolam.ppf
! (in alternativa, se si scrive ipflag nel file bolam.ppf, occorre settare le dimensioni
! corrispondenti nei programmi di grafica NCARG).
! Mar. 2017: velocizzato evitando calcoli inutili nel caso njump>1.
! Introdotto calcolo del flusso (trasporto) di vapore integrato verticalmente:
! il calcolo e' opzionale, e le compon. del vettore vengono scritte al posto di zus e zvs.
! Occorre scommentare codice (v. "Integrated Vapour Transport") e modificare boldis.def:
! sostituire le 2 righe dove si def. "wind at lowest level" con:
! 0.,    80.,   'INT. VAPOUR FLUX (KG/M/S)', 'V$', 15,  0, 0, 3, 3, 1
! 'INT. VAPOUR FLUX (KG/M/S)'
! Apr. 2015: aggiunte subroutines per scrittura in grib2 - messo a 2m il valore di roughness del
! momento sopra il quale si usa il secondo livello sopra il suolo per il calcolo delle variab. a 2m e 10m
! Nov. 2014: Versione con interfaccia per RTTOV-11 (output di RTTOV in formato grib2)
! (usa direttive #ifdef rttov)
! Feb. 2014: inserito calcolo di Integrated Water Vapour
! Gen. 2014: calcola lapse-rate in 2 modi e la temperatura di rugiada a 2 m (per la verifica)
! Ago. 2013: calcola cross-sections lungo direzioni arbitrarie (Oxana - solo per output grib2)
! Ridefinisce i missing values (def. in prebolam: -9999.) come valori nel punto (1,1) per poter
! post-processare (e plottare) mhf contenenti campi di bound. files definiti solo su cornici
! Lug. 2013: sono scritti nel file .ppf 4 campi di TG e QG ai livelli 1,3,5,7 invece che 1,2,3,4
! (v. def. di zqgr1 - ztgr4 e ztgr1 - ztgr4).
! Apr. 2013: aggiunta lettura e plot sea ice thickness (iceth), scritto assieme a fice
! Inserite formule rotaz. griglia in doppia precisione (necess. per griglie che includono il polo)
! Mar. 2013: aggiunto sea ice fraction fice (in lettura e plot)
! Genn. 2013: ritoccati param. di interp. e generalizz. per NLEVG livelli suolo (NLEVG in nfdr(17))
! Genn. 2012: cambiato modo di interp. T ai livelli std. del suolo
! Dic. 2012: corretto errore di interp. u10 e v10
! Mag. 2012:
! Versione che legge cloud water + ice nel mhf di bolam e calcola nubi di conseguenza
! Prevede inoltre che nel mhf sia scritta T e non theta
! nuovo calcolo nubi (ccloud) -
! messa roughness per T uguale a quella per q sul mare nel caso instabile (in vtsurf e vtsurflux)
! Valori di roughness per T e q uguali a variabili per 0<Ri<0.25 (come in BOLAM).
! Non si usa cloud water&ice come variabile da plottare, ma solo per calcolare nubi e T virt.
! Nov. 2011: introd. PZER variabile (definita da prebolam), che def. livelli ibridi
! Introd. matrice TEMP
! Giu. 2011: messa sat. risp. ad acqua nella tvlift (def. convenzionali di CAPE e lifted index)
! Genn. 2011: nuova CCLOUD
! RH ai livelli P non ha sovra-saturazione (eccetto 2% per smoothing grafico)
! Per calcolo THETAE, si usa la saturaz. mista rispetto ad acqua e ghiaccio

! Def. RH a 2m come umidità relativa rispetto all'acqua (WMO standard)

! Attenz.: si usa in vari punti NLBPER, ma non viene mai def., per cui
! e' FALSE - si potrebbe implementare, per il canale periodico,
! definendola in postbolam.inp (variabile POSTP)

! Roughness length piccola per calcoli diagnostici di variabili a 2m e 10 m, ma non per flussi
!-----------------------------------------------------------------------
 module mod_dimensions
   integer :: nlon, nlat, nlev, nlevg, nlonm1, nlonm2, nlonm3, nlonp1, nlevm1, nlevp1, nlatm1, nlatm2
 end module mod_dimensions

 module mod_postbolam
   use mod_dimensions
   integer, parameter :: nswtch = 15, nvarm = 7

!  dates and descriptor records

   integer :: nfdr(50)
   real    :: pdr(200)
   integer :: nday0,nhou0,nmin0,ndayr,nyrin,nmonin,ndayin,nhouin,nminin

!  model parameters

   integer, parameter :: nlev_snow=11
   integer :: ntsrc, ntslm
   real :: dx, dy, diffcv, diffct, diffcq, areat, areav, tddamp, diffdv, dbdiv
   real, dimension(:), allocatable :: hxv, hxt, tangu, tangv, spco, sig, dsig, sigint, dsgint, siginl, bndrel, sigalf, dsigalf
   real, dimension(:,:), allocatable :: fv, alatt, alont, alatz, alonz, alatu, alonu, alatv, alonv

!  namelist variables

   real :: dtstep
   integer :: nunic1, nunic2, nunhi, nunbc, nundg, nadj, nstep, ntsbou, nbc, nhist, ndrunt, nx1, nx2, ny1, ny2, nz1, nz2, &
              ierr_read1=0, ierr_read2=0
   real :: hordc, divdc, pzer, alfa

!  physical constants

   real :: rd, ra, cp, rcp, omd, eps, ep, g, we, rv, r00, stef, alp0, pi, ttr

!  auxiliary variables

   real, dimension(:,:), allocatable :: phig, alps, pst, transts, lapse_rate_1, lapse_rate_2
   real, dimension(:,:,:), allocatable :: phi, temp, tvirt, omeg, zdi1, zdi2

! main prognostic variables

   real, dimension(:,:), allocatable :: ps
   real, dimension(:,:,:), allocatable :: u, v, t, q, qc

! 'physics'

   integer :: nbolay,nltop,nlclo,nmclo,ncnlev
   real ::  rscs, hsnc, d1, d2, d3, d4, aksn, hdifs, qg1max, qg2max, qgdifs, ak, tadjt, tadjq, alsn, emsn, emis
   real, dimension(:), allocatable :: huc, uvml, hqml, vfac
   real, dimension(:,:), allocatable :: totpre, snow, conpre, albedo, snow_albedo, &
                                        rgm, rgq, fmask, snfall, runoff, runoff_tot, tskin, tgsurf, &
                                        cswfl, clwfl, cswdfl, chflux, cqflux, &
                                        cwvflux, csoilh_bott_flux, csoilw_bott_flux,       &
                                        cswflr, clwflr, cswdflr, chfluxr, cqfluxr, &
                                        cwvfluxr, csoilh_bott_fluxr, csoilw_bott_fluxr,&
                                        t2min, t2max, emismap1, emismap2, fice, iceth,            &
                                        cdm, cdt, cdq, hflux, qflux, qskin, qgsurf, fsnow, fice_soil_surf,               &
                                        cloudt, cloudh, cloudm, cloudl, &
                                        frvis, frirr, cloud, u10, v10, t2, q2, q2rel, td2, &
                                        precon, precls, snocon, snowls, iwv, richs, pbl, gust
   real, dimension(:,:,:), allocatable :: tg, qg, qgmax, qgmin, fice_soil, &
 snow_lev, snow_t, snow_fice, snow_age, snow_melt_age, snow_dens
   real, dimension(:), allocatable :: zztg, zzzg, zd

!  control switches

   logical :: nlana, nlbper, nlbfix, nlmpg, nlddmp, nldiff, nlorogd, nlrad, nlvdff, nlsurf, nlrain, nlcadj, &
              nlecrad, nlwaf, nlspong

! npar: no. of variables on constant pressure surfaces in output
! nlevpo: no. of constant pressure surfaces in output
! nlevto: no. of constant theta surfaces in output

  integer, parameter :: npar = 12, nlevpo0 = 50, nlevto = 4
  integer :: nlevpo

  integer, dimension(80) :: isflag           ! Surface parameters
  integer, dimension(npar,nlevpo0) :: ipflag ! Parameters at isobaric levels
  integer, dimension(npar,20) :: ipflag20    ! Parameters at 20 isobaric levels for compatibility with NCAR Graphics of ppf file
  integer, dimension(nlevto) :: itflag       ! Parameters at constant theta levels
  integer, dimension(nlevpo0) :: iplevo      ! Integer value (hPa) of isobaric levels
  integer, dimension(npar) :: itropflag, imwflag ! Parameters at tropopause and maximum wind levels

  logical :: output_format_ppf, output_format_grib2
  real(kind=8) :: zfac, zx0, zy0, zlatt, zlatv
  real :: valmiss = -9999.
  integer :: ivalmiss = -9999

! Variables at tropopause level
  real, dimension(:,:), allocatable :: p_tropopause, h_tropopause, t_tropopause, u_tropopause, v_tropopause, &
        q_tropopause, rh_tropopause, om_tropopause, rv_tropopause, pv_tropopause, the_tropopause

! Variables at maximum winf level
  real, dimension(:,:), allocatable :: p_maxwind, h_maxwind,  t_maxwind, u_maxwind, v_maxwind, &
        q_maxwind, rh_maxwind, om_maxwind, rv_maxwind, pv_maxwind, the_maxwind
  integer, dimension(:,:), allocatable :: klev_maxwind

 end module mod_postbolam
!-----------------------------------------------------------------------

!  Module containing defined values of pressure levels in Pa

 module mod_rplev
   use mod_postbolam, only : nlevpo0, nlevto
   real plevo(nlevpo0)
 end module mod_rplev

 module mod_jump
   integer :: njump
 end module mod_jump

 module mod_rwprer
   real, dimension(:,:), allocatable :: totprer, conprer, snfallr, runoffr, runoff_tot_r
 end module mod_rwprer

 module mod_surfadd
   real, dimension(:,:), allocatable :: t05, t005, q05, q05rel, tg005, tg010, tg020, tg050, tg100
 end module mod_surfadd

 module mod_cross
   logical :: ncrossec
   integer :: loncr1, loncr2, loncr3, latcr1, latcr2, latcr3, cross_number=0, npoint_cross_max=1000, cross_interp_gauss
   real, dimension(10) :: lonini_cross, latini_cross, lonfin_cross, latfin_cross
   integer, dimension(10) :: npoint_cross
 end module mod_cross

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
   real, dimension(:,:,:), allocatable :: p, qcw, qci, clfrac
   real, dimension(:,:), allocatable :: fseaice
   integer, dimension(:,:), allocatable :: flag_cloud_conv
   real :: w_crit_conv=0.5, fracw, z1, z2, z3, z4, z5, z6
   integer, dimension(:), allocatable :: &
 sensor_chan_id  ! List of "true" channels index, not RTTOV channels list
   real*8, dimension(:), allocatable :: &
 sensor_chan_cw  ! Central Wave Number (m^-1) of elaborated channels
   real, dimension(:,:,:), allocatable :: radiance, radiance_bt, radiance_refl, &
 radiance_clear, radiance_clear_bt, radiance_clear_refl, emis_rttov
   real :: zoom_xini=0.00, zoom_xfin=1.00, zoom_yini=0.00, zoom_yfin=1.00, alon1, alat1, val_missing=-9999.
   integer :: npoint_jump, iini, ifin, jini, jfin
 end module mod_rttov
#endif

! -----------------------------
!########################################################################################################
 program postbolam

 use mod_postbolam
 use mod_rplev
 use mod_jump
 use mod_rwprer
 use mod_surfadd
 use mod_cross

! ---- For RTTOV simulation ---

#ifdef rttov
 use mod_rttov
#endif

! -----------------------------

 logical nlwafp
 integer, dimension(nlevpo0*npar) :: ipflag0
 integer :: ierr_open1, ierr_open2

 integer, dimension(5) :: idate0=(/2015, 12, 1, 0, 0/)
 integer, dimension(3) :: iperiod=(/0, 0, 0/)

 integer, parameter :: npoint=1
 integer, dimension(npoint) :: ipoint, jpoint
 real, dimension(npoint,1) :: lon_point=reshape((/11.34/),(/npoint,1/)), &
 lat_point=reshape((/44.49/),(/npoint,1/)), lon_rot_point, lat_rot_point

 namelist/postp/njump,nlwafp,prtop,iterf,ncrossec,output_format_ppf,output_format_grib2,             &
                loncr1,loncr2,loncr3,latcr1,latcr2,latcr3,                                           &
                lonini_cross,latini_cross,lonfin_cross,latfin_cross,npoint_cross,cross_interp_gauss, &
                isflag,ipflag0,itflag,iplevo,itropflag,imwflag

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
      open (nunic1,file='bolam_atm.mhf',status='old',form='unformatted',iostat=ierr_open1)
      print *
      if (ierr_open1 /= 0) then
        print *,'Not found input file bolam_atm.mhf'
        stop
      else
        print *,'Input file bolam_atm.mhf opened on the unit ',nunic1
      endif
      open (nunic2,file='bolam_soil.mhf',status='old',form='unformatted',iostat=ierr_open2)
      print *
      if (ierr_open2 /= 0) then
        print *,'Not found input file bolam_soil.mhf'
        stop
      else
        print *,'Input file bolam_soil.mhf opened on the unit ',nunic2
      endif
      print *

      read(nunic1) nfdr
      read(nunic1) pdr

      rewind nunic1

      nlon =nfdr(2)
      nlat =nfdr(3)
      nlev =nfdr(4)
      nlevg=nfdr(15)
      nlonm1=nlon-1
      nlonm2=nlon-2
      nlonm3=nlon-3
      nlonp1=nlon+1
      nlevm1=nlev-1
      nlevp1=nlev+1
      nlatm1=nlat-1
      nlatm2=nlat-2

! Dinamic array allocation:

 allocate(fv(nlon,nlat), stat=ierr)
 allocate(hxv(nlat), stat=ierr)
 allocate(hxt(nlat), stat=ierr)
 allocate(alatt(nlon,nlat), stat=ierr)
 allocate(alont(nlon,nlat), stat=ierr)
 allocate(alatz(nlon,nlat), stat=ierr)
 allocate(alonz(nlon,nlat), stat=ierr)
 allocate(alatu(nlon,nlat), stat=ierr)
 allocate(alonu(nlon,nlat), stat=ierr)
 allocate(alatv(nlon,nlat), stat=ierr)
 allocate(alonv(nlon,nlat), stat=ierr)
 allocate(tangu(nlat), stat=ierr)
 allocate(tangv(nlat), stat=ierr)
 allocate(spco(nlev), stat=ierr)
 allocate(sig(nlev), stat=ierr)
 allocate(dsig(nlev), stat=ierr)
 allocate(sigint(nlev), stat=ierr)
 allocate(dsgint(nlev), stat=ierr)
 allocate(siginl(nlev), stat=ierr)
 allocate(sigalf(nlev), stat=ierr)
 allocate(dsigalf(nlev),stat=ierr)
 allocate(zzzg(nlevg), stat=ierr)
 allocate(zztg(nlevg), stat=ierr)
 allocate(zd  (nlevg), stat=ierr)
 allocate(phig(nlon,nlat), stat=ierr)
 allocate(phi(nlon,nlat,nlev), stat=ierr)
 allocate(temp(nlon,nlat,nlev), stat=ierr)
 allocate(tvirt(nlon,nlat,nlev), stat=ierr)
 allocate(alps(1:nlonp1,0:nlat), stat=ierr)
 allocate(pst(nlon,nlat), stat=ierr)
 allocate(omeg(nlon,nlat,nlev), stat=ierr)
 allocate(transts(nlon,nlat), stat=ierr)
 allocate(lapse_rate_1(nlon,nlat), stat=ierr)
 allocate(lapse_rate_2(nlon,nlat), stat=ierr)
 allocate(zdi1(nlon,nlat,nlev), stat=ierr)
 allocate(zdi2(nlon,nlat,nlev), stat=ierr)
 allocate(u(nlon,nlat,nlev), stat=ierr)
 allocate(v(nlon,nlat,nlev), stat=ierr)
 allocate(t(nlon,nlat,nlev), stat=ierr)
 allocate(q(nlon,nlat,nlev), stat=ierr)
 allocate(qc(nlon,nlat,nlev), stat=ierr)
 allocate(ps(1:nlonp1,0:nlat), stat=ierr)

 allocate(tg(nlon,nlat,nlevg), stat=ierr)
 allocate(qg(nlon,nlat,nlevg), stat=ierr)
 allocate(qgmax(nlon,nlat,nlevg), stat=ierr)
 allocate(qgmin(nlon,nlat,nlevg), stat=ierr)
 allocate(fice_soil(nlon,nlat,nlevg), stat=ierr)
 allocate(snow_lev(nlon,nlat,nlev_snow), stat=ierr)
 allocate(snow_t(nlon,nlat,nlev_snow), stat=ierr)
 allocate(snow_fice(nlon,nlat,nlev_snow), stat=ierr)
 allocate(snow_age(nlon,nlat,nlev_snow), stat=ierr)
 allocate(snow_melt_age(nlon,nlat,nlev_snow), stat=ierr)
 allocate(snow_dens(nlon,nlat,nlev_snow), stat=ierr)

 allocate(totpre(nlon,nlat), stat=ierr)
 allocate(snow(nlon,nlat), stat=ierr)
 allocate(runoff(nlon,nlat), stat=ierr)
 allocate(runoff_tot(nlon,nlat), stat=ierr)
 allocate(conpre(nlon,nlat), stat=ierr)
 allocate(albedo(nlon,nlat), stat=ierr)
 allocate(snow_albedo(nlon,nlat), stat=ierr)
 allocate(rgm(nlon,nlat), stat=ierr)
 allocate(rgq(nlon,nlat), stat=ierr)
 allocate(fmask(nlon,nlat), stat=ierr)
 allocate(snfall(nlon,nlat), stat=ierr)
 allocate(tskin(nlon,nlat), stat=ierr)
 allocate(tgsurf(nlon,nlat), stat=ierr)
 allocate(cswfl(nlon,nlat), stat=ierr)
 allocate(clwfl(nlon,nlat), stat=ierr)
 allocate(cswdfl(nlon,nlat), stat=ierr)
 allocate(chflux(nlon,nlat), stat=ierr)
 allocate(cqflux(nlon,nlat), stat=ierr)
 allocate(cwvflux(nlon,nlat), stat=ierr)
 allocate(csoilh_bott_flux(nlon,nlat), stat=ierr)
 allocate(csoilw_bott_flux(nlon,nlat), stat=ierr)
 allocate(t2min(nlon,nlat), stat=ierr)
 allocate(t2max(nlon,nlat), stat=ierr)
 allocate(emismap1(nlon,nlat), stat=ierr)
 allocate(emismap2(nlon,nlat), stat=ierr)
 allocate(fice(nlon,nlat), stat=ierr)
 allocate(iceth(nlon,nlat), stat=ierr)
 allocate(huc(nlev), stat=ierr)
 allocate(uvml(nlev), stat=ierr)
 allocate(hqml(nlev), stat=ierr)
 allocate(vfac(nlev), stat=ierr)
 allocate(cdm(nlon,nlat), stat=ierr)
 allocate(cdt(nlon,nlat), stat=ierr)
 allocate(cdq(nlon,nlat), stat=ierr)
 allocate(hflux (nlon,nlat), stat=ierr)
 allocate(qflux (nlon,nlat), stat=ierr)
 allocate(qskin(nlon,nlat), stat=ierr)
 allocate(qgsurf(nlon,nlat), stat=ierr)
 allocate(fice_soil_surf(nlon,nlat), stat=ierr)
 allocate(fsnow(nlon,nlat), stat=ierr)
 allocate(cloudt(nlon,nlat), stat=ierr)
 allocate(cloudh(nlon,nlat), stat=ierr)
 allocate(cloudm(nlon,nlat), stat=ierr)
 allocate(cloudl(nlon,nlat), stat=ierr)
 allocate(frvis(nlon,nlat), stat=ierr)
 allocate(frirr(nlon,nlat), stat=ierr)
 allocate(cloud(nlon,nlat), stat=ierr)
 allocate(u10(nlon,nlat), stat=ierr)
 allocate(v10(nlon,nlat), stat=ierr)
 allocate(t2(nlon,nlat), stat=ierr)
 allocate(q2(nlon,nlat), stat=ierr)
 allocate(q2rel(nlon,nlat), stat=ierr)
 allocate(td2(nlon,nlat), stat=ierr)
 allocate(precon(nlon,nlat), stat=ierr)
 allocate(precls(nlon,nlat), stat=ierr)
 allocate(snocon(nlon,nlat), stat=ierr)
 allocate(snowls(nlon,nlat), stat=ierr)
 allocate(totprer(nlon,nlat), stat=ierr)
 allocate(conprer(nlon,nlat), stat=ierr)
 allocate(snfallr(nlon,nlat), stat=ierr)
 allocate(runoffr(nlon,nlat), stat=ierr)
 allocate(runoff_tot_r(nlon,nlat), stat=ierr)
 allocate(cswflr(nlon,nlat), stat=ierr)
 allocate(clwflr(nlon,nlat), stat=ierr)
 allocate(cswdflr(nlon,nlat), stat=ierr)
 allocate(chfluxr(nlon,nlat), stat=ierr)
 allocate(cqfluxr(nlon,nlat), stat=ierr)
 allocate(cwvfluxr(nlon,nlat), stat=ierr)
 allocate(csoilh_bott_fluxr(nlon,nlat), stat=ierr)
 allocate(csoilw_bott_fluxr(nlon,nlat), stat=ierr)
 allocate(t05(nlon,nlat), stat=ierr)
 allocate(t005(nlon,nlat), stat=ierr)
 allocate(q05(nlon,nlat), stat=ierr)
 allocate(q05rel(nlon,nlat), stat=ierr)
 allocate(tg005(nlon,nlat), stat=ierr)
 allocate(tg010(nlon,nlat), stat=ierr)
 allocate(tg020(nlon,nlat), stat=ierr)
 allocate(tg050(nlon,nlat), stat=ierr)
 allocate(tg100(nlon,nlat), stat=ierr)
 allocate(iwv(nlon,nlat), stat=ierr)
 allocate(p_tropopause(nlon,nlat), stat=ierr)
 allocate(h_tropopause(nlon,nlat), stat=ierr)
 allocate(t_tropopause(nlon,nlat), stat=ierr)
 allocate(u_tropopause(nlon,nlat), stat=ierr)
 allocate(v_tropopause(nlon,nlat), stat=ierr)
 allocate(q_tropopause(nlon,nlat), stat=ierr)
 allocate(rh_tropopause(nlon,nlat), stat=ierr)
 allocate(om_tropopause(nlon,nlat), stat=ierr)
 allocate(rv_tropopause(nlon,nlat), stat=ierr)
 allocate(pv_tropopause(nlon,nlat), stat=ierr)
 allocate(the_tropopause(nlon,nlat), stat=ierr)
 allocate(klev_maxwind(nlon,nlat), stat=ierr)
 allocate(p_maxwind(nlon,nlat), stat=ierr)
 allocate(h_maxwind(nlon,nlat), stat=ierr)
 allocate(t_maxwind(nlon,nlat), stat=ierr)
 allocate(u_maxwind(nlon,nlat), stat=ierr)
 allocate(v_maxwind(nlon,nlat), stat=ierr)
 allocate(q_maxwind(nlon,nlat), stat=ierr)
 allocate(rh_maxwind(nlon,nlat), stat=ierr)
 allocate(om_maxwind(nlon,nlat), stat=ierr)
 allocate(rv_maxwind(nlon,nlat), stat=ierr)
 allocate(pv_maxwind(nlon,nlat), stat=ierr)
 allocate(the_maxwind(nlon,nlat), stat=ierr)
 allocate(richs(nlon,nlat), stat=ierr)
 allocate(pbl(nlon,nlat), stat=ierr)
 allocate(gust(nlon,nlat), stat=ierr)

! ---- For RTTOV simulation ---

#ifdef rttov
 allocate (p(nlon,nlat,nlev))
 allocate (qcw(nlon,nlat,nlev))
 allocate (qci(nlon,nlat,nlev))
 allocate (clfrac(nlon,nlat,nlev))
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

      open (11,file='postbolam.inp',status='old')
      rewind 11
      read(11,postp)
#ifdef oper
      write(*,postp)
#endif
      close (11)

      ipflag(1:npar,1:nlevpo0)=reshape(ipflag0,(/npar,nlevpo0/))

!  Definition of values of pressure levels in Pa on which fields are computed

      do jlev = 1,nlevpo0
        if (iplevo(jlev) /= ivalmiss) then
          plevo(jlev) = float(iplevo(jlev))*1.e2
        else
          exit
        endif
      enddo
      nlevpo = jlev - 1

      iist = 0
      iist2= 0

      totpre(:,:)=0.
      conpre(:,:)=0.
      snfall(:,:)=0.
      runoff(:,:)=0.
      runoff_tot(:,:)=0.
      cswfl(:,:)=0.
      clwfl(:,:)=0.
      chflux(:,:)=0.
      cqflux(:,:)=0.
      cwvflux(:,:)=0.
      csoilh_bott_flux(:,:)=0.
      csoilw_bott_flux(:,:)=0.

! Define physical parameters

       rd     = 287.05
       ra     = 6371.e+3
       cp     = 1004.6
       rcp    = rd/cp
       zcpv   = 1869.46
       zdelta = zcpv/cp
       omd    = zdelta-1.
       rv     = 461.51
       eps    = rd/rv
       ep     = rv/rd -1.
       g      = 9.807
       we     = 7.292e-5
       r00    = 1373.
       stef   = 5.67e-8
       pi     = abs(acos(-1.))
       alp0   = alog(1.e5)
       ttr    = 273.15
       ak     = 0.4
       hsnc   =  .02

! Readinf of constant (in time) model physiographical parameters
! (phig, fmask, rgm, rgq)

      call rd_param_const

! call outgraph(80102,1,nfdr(2),nfdr(3),pdr(39),pdr(38),pdr(5),pdr(4)+0.5*pdr(1),pdr(2),pdr(1),&
! 2, 0,  0,  1, 1, 0,idate0(1:5),iperiod(1:3),fmask(:,:),1.,0.)
! call outgraph(80102,1,nfdr(2),nfdr(3),pdr(39),pdr(38),pdr(5),pdr(4)+0.5*pdr(1),pdr(2),pdr(1),&
! 0, 0,193,  1, 1, 0,idate0(1:5),iperiod(1:3),phig(:,:)/g,1.,0.)
! call outgraph(80102,1,nfdr(2),nfdr(3),pdr(39),pdr(38),pdr(5),pdr(4)+0.5*pdr(1),pdr(2),pdr(1),&
! 0, 0,193,  1, 1, 0,idate0(1:5),iperiod(1:3),rgm(:,:),1.,0.)
! call outgraph(80102,1,nfdr(2),nfdr(3),pdr(39),pdr(38),pdr(5),pdr(4)+0.5*pdr(1),pdr(2),pdr(1),&
! 0, 0,193,  1, 1, 0,idate0(1:5),iperiod(1:3),rgq(:,:),1.,0.)

 do while (.true.)

! Read postprocessing file

      print*, 'Read ',iist+1,' instant from bolam_atm.mhf and bolam_soil.mhf'

      call rdmhf_atm(nunic1)
      if (ierr_read1 /= 0) exit
      call rdmhf_soil(nunic2)
      if (ierr_read2 /= 0) exit

!! For control:
!print *,'***************'
! call anti_rot_grid(pdr(39), pdr(38), lon_point(npoint,1), lat_point(npoint,1), &
! lon_rot_point(npoint,1), lat_rot_point(npoint,1), npoint, 1)
! do i=1,npoint
!   ipoint(i)=nint((lon_rot_point(i,1)-pdr(5))/pdr(2))
!   jpoint(i)=nint((lat_rot_point(i,1)-pdr(4))/pdr(1))
!   print *,'Point ',i,ipoint(i),jpoint(i)
! enddo

do jlat = 1,nlat
do jlon = 1,nlon
  if (runoffr(jlon,jlat)>1.E3) print *,'runoff *** ',runoffr(jlon,jlat),jlon,jlat
enddo
enddo

 print *,'poisk 1',totprer(358,164),maxval(totprer(:,:)),maxloc(totprer(:,:))

! Redefinition of missing values (case of use of frames for boundary cond. files)

      do jlat = 1,nlat
      do jlon = 1,nlon
       if (ps(jlon,jlat) == valmiss)  ps(jlon,jlat) = ps(1,1)
      enddo
      enddo

      do jklev = 1,nlev
      do jlat = 1,nlat
      do jlon = 1,nlon
      if(u(jlon,jlat,jklev)    == valmiss) u   (jlon,jlat,jklev) = u   (1,1,jklev)
      if(v(jlon,jlat,jklev)    == valmiss) v   (jlon,jlat,jklev) = v   (1,1,jklev)
      if(temp(jlon,jlat,jklev) == valmiss) temp(jlon,jlat,jklev) = temp(1,1,jklev)
      if(q(jlon,jlat,jklev)    == valmiss) q   (jlon,jlat,jklev) = q   (1,1,jklev)
      if(qc(jlon,jlat,jklev)   == valmiss) qc  (jlon,jlat,jklev) = qc  (1,1,jklev)
      enddo
      enddo
      enddo

      q  = max(q, 0.)
      qc = max(qc, 0.)

! Define model parameters

      zdlat  = pdr(1)
      zdlon  = pdr(2)
      dtstep = pdr(3)
      ivalt  = nfdr(10)*10000+nfdr(11)*100+nfdr(12)
      nhist  = nfdr(17)

      if(iist==0) then
       if(ivalt==0) print*, "The 1st processed instant is an initial condition"
       if(njump.ne.1) print '(A8,I2,A53)'," njump =",njump,": note that the 1st forecast instant is not processed"
      endif

! zalat0 is the southernmost latitude (v-point)
! zalon0 is the westernmost longitude (t-point)

      zalat0 = pdr(4)
      zalon0 = pdr(5)
      if (zalon0 > 180.) zalon0=zalon0-360.
      dx     = ra*pi/180.*zdlon
      dy     = ra*pi/180.*zdlat

! Definition of real grid point latitudes and longitudes (in degrees) using double precision
! Definition of geographical factors hxv, hxt

    call rot_grid(pdr(39),pdr(38),zalon0,          zalat0+zdlat*0.5,zdlon,zdlat,alont,alatt,nlon,nlat)
    call rot_grid(pdr(39),pdr(38),zalon0+zdlon*0.5,zalat0+zdlat*0.5,zdlon,zdlat,alonu,alatu,nlon,nlat)
    call rot_grid(pdr(39),pdr(38),zalon0,          zalat0         ,zdlon,zdlat,alonv,alatv,nlon,nlat)
    call rot_grid(pdr(39),pdr(38),zalon0+zdlon*0.5,zalat0         ,zdlon,zdlat,alonz,alatz,nlon,nlat)

    zfac = dabs(dacos(-1.d0))/180.d0
    zx0  = dble(pdr(39))*zfac
    zy0  = dble(pdr(38))*zfac
    do jlat = 1, nlat
    zlatv=(dble(zalat0)+                (jlat-1)*dble(zdlat))*zfac
    zlatt=(dble(zalat0)+dble(zdlat)/2.d0+(jlat-1)*dble(zdlat))*zfac
    hxt(jlat)    = dcos(zlatt)
    hxv(jlat)    = dcos(zlatv)

      do jlon = 1, nlon

!  Coriolis parameter fv(jlon,jlat), defined on vorticity points only

      fv(jlon,jlat) = 2.*we*sin(pi*alatz(jlon,jlat)/180.)
      enddo
    enddo

! Definition of model levels (read from pdr(40)...) and of depending variables alfk and betk
! Note that sig(1) corresponds to level 1+1/2

      pzer=pdr(36)
      alfa=pdr(37)
      do jklev=1,nlev
      sig(jklev) = pdr(39+jklev)
      enddo
      dsig(1)  = sig(1)
      sigint(1)=.5*sig(1)
      siginl(1)=alog(sigint(1))

      do jklev = 2,nlev
      dsig(jklev) = sig(jklev)-sig(jklev-1)
      sigint(jklev)=.5*(sig(jklev)+sig(jklev-1))
      siginl(jklev)=alog(sigint(jklev))
      enddo

      sigalf(:)  = sigint(:)**3*(1.+alfa*(1.-sigint(:))*(2.-sigint(:)))
      dsigalf(:) = dsig(:)*sigint(:)**2*(3.+alfa*(6.-12.*sigint(:)+5.*sigint(:)**2))

      do jklev=1,nlev-1
      dsgint(jklev) = sigint(jklev+1) - sigint(jklev)
      enddo

! Definition of potential temp. (put in t), geopotential and vertical velocity

      call alpstv
      call phicomp

! Different computations of omega, depending on type of advection used in the model (waf or fbas)
! In case of waf, the model time step is used.
! If nlwafp true, but dtstep=0 (typically analysis to be plotted) the fbas option is used

      if(nlwafp.and.dtstep.lt.1.e-5) then
      print*, "Attempting option WAF for computing OMEGA"
      print*, "from non model-output fields: switched to option FBAS"
      print*, 'dtstep =', dtstep
      nlwafp=.false.
      endif

      call diverg(nlwafp)

! Definition of surface parameters v10, t2, q2, q2rel, td2

      call vtsurf
      call vtsurflux

!  Fix over high mountains.....

      do jlat = 1, nlat
      do jlon = 1, nlon
      zglin = max (.0, min((phig(jlon,jlat)-10000.)/15000., 1.))
      u10(jlon,jlat) = (1.-zglin)*u10(jlon,jlat) + zglin*u(jlon,jlat,nlev)
      v10(jlon,jlat) = (1.-zglin)*v10(jlon,jlat) + zglin*v(jlon,jlat,nlev)
      if (t2(jlon,jlat).lt.temp(jlon,jlat,nlev)-tzer) then
      zglin = max (.0, min((phig(jlon,jlat)-8000.)/20000., 1.))
      t2(jlon,jlat) = (1.-zglin)*t2(jlon,jlat) + zglin*(temp(jlon,jlat,nlev)-tzer)
      q2(jlon,jlat) = (1.-zglin)*q2(jlon,jlat) + zglin*q(jlon,jlat,nlev)
      call td_tetens (t2(jlon,jlat), ps(jlon,jlat), q2(jlon,jlat), eps, td2(jlon,jlat))
      endif
      enddo
      enddo

! Definition of cloud cover

      call ccloud

! Wind gust computation as a function of the PBL top height

      do jlat = 1, nlat
      do jlon = 1, nlon
        gust(jlon,jlat) = sqrt(u10(jlon,jlat)**2+v10(jlon,jlat)**2)
        do jklev = nlev, nlev-8, -1
          jkup = jklev
          if (pbl(jlon,jlat)*g < phi(jlon,jlat,jklev)-phig(jlon,jlat)-1.) exit
        enddo
        if (jkup < nlev) then
          do jklev = nlev, jkup+1, -1
            zwink = sqrt(u(jlon,jlat,jklev)**2+v(jlon,jlat,jklev)**2)
            gust(jlon,jlat) = max (gust(jlon,jlat), zwink)
          enddo
        endif
      enddo
      enddo

!***********************************************************************
      iist=iist+1
      do jlat=1,nlat
      do jlon=1,nlon
      totpre(jlon,jlat)=totpre(jlon,jlat)+totprer(jlon,jlat)
      conpre(jlon,jlat)=conpre(jlon,jlat)+conprer(jlon,jlat)
      snfall(jlon,jlat)=snfall(jlon,jlat)+snfallr(jlon,jlat)
      runoff(jlon,jlat)=runoff(jlon,jlat)+runoffr(jlon,jlat)
      runoff_tot(jlon,jlat)=runoff_tot(jlon,jlat)+runoff_tot_r(jlon,jlat)
      cswfl(jlon,jlat)=cswfl(jlon,jlat)+cswflr(jlon,jlat)
      clwfl(jlon,jlat)=clwfl(jlon,jlat)+clwflr(jlon,jlat)
      chflux(jlon,jlat)=chflux(jlon,jlat)+chfluxr(jlon,jlat)
      cqflux(jlon,jlat)=cqflux(jlon,jlat)+cqfluxr(jlon,jlat)
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
      if(mod(iist+iist0,njump).eq.0) then
      iist2=iist2+1
      call postpcpm(iist2,prtop,iterf)
!***********************************************************************
      print*, 'Instant number ',iist,' post-processed'
      totpre(:,:)=0.
      conpre(:,:)=0.
      snfall(:,:)=0.
      runoff(:,:)=0.
      runoff_tot(:,:)=0.
      cswfl(:,:)=0.
      clwfl(:,:)=0.
      chflux(:,:)=0.
      cqflux(:,:)=0.
      cwvflux(:,:)=0.
      csoilh_bott_flux(:,:)=0.
      csoilw_bott_flux(:,:)=0.
      endif

 enddo ! while (.true.)

stop
end
! ######################################################################
      subroutine rdmhf_atm(kunit)

!   Read Model History File with atmospheric variables from kunit

 use mod_postbolam

!  Read descriptor records

      read(kunit, iostat=ierr_read1) nfdr

      if (ierr_read1 /= 0) then
        close(kunit)
        print*,'EOF reached on file unit ', kunit
        return
      endif

      read(kunit) pdr

! Read atmospheric variables

      call rrec2 (kunit, nlon, nlat, ps(1:nlon,1:nlat))

      do jklev=1,nlev
      call rrec2 (kunit, nlon, nlat, u(1:nlon,1:nlat,jklev))
      enddo
      do jklev=1,nlev
      call rrec2 (kunit, nlon, nlat, v(1:nlon,1:nlat,jklev))
      enddo
      do jklev=1,nlev
      call rrec2 (kunit, nlon, nlat, temp(1:nlon,1:nlat,jklev))
      enddo
      do jklev=1,nlev
      call rrec2 (kunit, nlon, nlat, q(1:nlon,1:nlat,jklev))
      enddo
      do jklev=1,nlev
      call rrec2 (kunit, nlon, nlat, qc(1:nlon,1:nlat,jklev))
      enddo

      return
      end
!#######################################################################
      subroutine rdmhf_soil(kunit)

!   Read Model History File with surface/sea/soil variables from kunit

 use mod_postbolam
 use mod_rwprer
 integer :: nfdr_local(50)
 real    :: pdr_local(200)
 real, dimension(nlon,nlat) :: field2d_add

!  Read descriptor records

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

      call rrec2 (kunit, nlon, nlat, field2d_add) ! lai
      call rrec2 (kunit, nlon, nlat, field2d_add) ! fveg
      call rrec2 (kunit, nlon, nlat, rgm)
      call rrec2 (kunit, nlon, nlat, rgq)
      call rrec2 (kunit, nlon, nlat, iceth)
      call rrec2 (kunit, nlon, nlat, fice)
      call rrec2 (kunit, nlon, nlat, albedo)
      call rrec2 (kunit, nlon, nlat, emismap1)
      call rrec2 (kunit, nlon, nlat, emismap2)

! Prognostic cloud and precipitation variables at the surface

      call rrec2 (kunit, nlon, nlat, cloudt)
      call rrec2 (kunit, nlon, nlat, totprer)
      call rrec2 (kunit, nlon, nlat, conprer)
      call rrec2 (kunit, nlon, nlat, snfallr)

! Prognostic surface and soil/sea fields

      call rrec2 (kunit, nlon, nlat, tskin)
      call rrec2 (kunit, nlon, nlat, tgsurf)
      do jklev = 1, nlevg
      call rrec2 (kunit, nlon, nlat, tg(1,1,jklev))
      enddo

      call rrec2 (kunit, nlon, nlat, qskin)
      call rrec2 (kunit, nlon, nlat, qgsurf)
      do jklev = 1, nlevg
      call rrec2 (kunit, nlon, nlat, qg(1,1,jklev))
      enddo

      call rrec2 (kunit, nlon, nlat, fice_soil_surf)
      do jklev = 1, nlevg
      call rrec2 (kunit, nlon, nlat, fice_soil(1,1,jklev))
      enddo

      call rrec2 (kunit, nlon, nlat, snow)
      do jklev = 1, nlev_snow
      call rrec2 (kunit, nlon, nlat, snow_lev(1,1,jklev))
      enddo
      do jklev = 1, nlev_snow
      call rrec2 (kunit, nlon, nlat, snow_t(1,1,jklev))
      enddo
      do jklev = 1, nlev_snow
      call rrec2 (kunit, nlon, nlat, snow_fice(1,1,jklev))
      enddo
      do jklev = 1, nlev_snow
      call rrec2 (kunit, nlon, nlat, snow_age(1,1,jklev))
      enddo
      do jklev = 1, nlev_snow
      call rrec2 (kunit, nlon, nlat, snow_melt_age(1,1,jklev))
      enddo
      do jklev = 1, nlev_snow
      call rrec2 (kunit, nlon, nlat, snow_dens(1,1,jklev))
      enddo
      call rrec2 (kunit, nlon, nlat, snow_albedo)

      call rrec2 (kunit, nlon, nlat, cswflr)
      call rrec2 (kunit, nlon, nlat, clwflr)
      call rrec2 (kunit, nlon, nlat, chfluxr)
      call rrec2 (kunit, nlon, nlat, cqfluxr)
      call rrec2 (kunit, nlon, nlat, t2min)
      call rrec2 (kunit, nlon, nlat, t2max)
      call rrec2 (kunit, nlon, nlat, runoffr)
      call rrec2 (kunit, nlon, nlat, runoff_tot_r)

! Additional 2d fields

      do iwr = 1, 3
      call rrec2 (kunit, nlon, nlat, field2d_add)
      enddo
!      call rrec2 (kunit, nlon, nlat, cwvfluxr)
!      call rrec2 (kunit, nlon, nlat, csoilh_bott_fluxr)
!      call rrec2 (kunit, nlon, nlat, csoilw_bott_fluxr)

      return
      end
!#######################################################################
subroutine rd_param_const

! Reads from additional input file
! constant (in time) model physiographical parameters:
! land-sea fraction and orography*g

use mod_postbolam

character(len=30) :: filerd="model_param_constant.bin"
integer :: iunit=11, ierr_open, nlon_local, nlat_local, nlevg_local, nst_local, nvt_local, ierr=0, ird, jklev
real :: dlon_local, dlat_local, x0_local, y0_local, alon0_local, alat0_local
real, dimension(nlon,nlat) :: zread

 open (iunit,file=trim(filerd),status='old',form='unformatted',iostat=ierr_open)
 if (ierr_open /= 0) then
   print *
   print *,'Not found input ',trim(filerd)
   print *,"   STOP"
   stop
 endif

 read (iunit) nlon_local, nlat_local, nlevg_local, dlon_local, dlat_local, &
 x0_local, y0_local, alon0_local, alat0_local, nst_local, nvt_local

 alat0_local=alat0_local-dlat_local*0.5

 if (nlon_local /= nfdr(2)) ierr=ierr+1
 if (nlat_local /= nfdr(3)) ierr=ierr+1
 if (nlevg_local /= nfdr(15)) ierr=ierr+1
 if (dlon_local /= pdr(2)) ierr=ierr+1
 if (dlat_local /= pdr(1)) ierr=ierr+1
 if (x0_local /= pdr(39)) ierr=ierr+1
 if (y0_local /= pdr(38)) ierr=ierr+1
 if (alon0_local /= pdr(5)) ierr=ierr+1
 if (alat0_local /= pdr(4)) ierr=ierr+1

 if (ierr /= 0) then
   print *,"Error in header parameters in input file ,",trim(filerd),", not coincident with defined parameters"
   print *,"Model nlon, nlat, nlevg, dlon, dlat, x0, y0, alon0, alat0 :", &
 nfdr(2), nfdr(3), nfdr(15), pdr(2), pdr(1), pdr(39), pdr(38), pdr(5), pdr(4)
   print *,"Read nlon, nlat, nlevg, dlon, dlat, x0, y0, alon0, alat0 :", &
 nlon_local, nlat_local, nlevg_local, dlon_local, dlat_local, x0_local, y0_local, alon0_local, alat0_local
   print *,"   STOP"
   stop
 endif

! land sea mask and orography

 call rrec2 (iunit, nlon, nlat, fmask)
 call rrec2 (iunit, nlon, nlat, phig)
 phig(:,:) = phig(:,:)*g

 call rrec2 (iunit, nlon, nlat, zread) ! htopvar
 do ird=1,nst_local+1
   call rrec2 (iunit, nlon, nlat, zread) ! soil_map
 enddo
 do ird=1,nvt_local+1
   call rrec2 (iunit, nlon, nlat, zread) ! veg_map
 enddo

 do jklev=1,nlevg
   call rrec2 (iunit, nlon, nlat, qgmax(:,:,jklev))
 enddo

 do jklev=1,nlevg
   call rrec2 (iunit, nlon, nlat, qgmin(:,:,jklev))
 enddo

 close (iunit)

return
end
!#######################################################################
subroutine rrec2 (kunit, nlon, nlat, vect)
implicit none

 integer :: kunit, nlon, nlat
 real(4), dimension(nlon,nlat) :: vect

 read(kunit) vect(1:nlon, 1:nlat)

return
end
!#######################################################################
      subroutine rrec2_old (kunit, nlon, nlat, vect)
      real vect (nlon,nlat)

      do jlat = 1, nlat
      read (kunit) (vect(jlon,jlat), jlon = 1, nlon)
      enddo

      return
      end
!#######################################################################
      subroutine alpstv

! Postprocessing version

 use mod_postbolam

!  Extra ps values at southern and eastern boundary by linear extrapolation.
!  Ps is dimensioned as ps(1:nlonp1,0:nlat)
!  These values are needed to compute divergence on such boundaries

      do jlon = 1, nlon
      ps(jlon,0) = 2.*ps(jlon,1)-ps(jlon,2)
      enddo

      do jlat = 0, nlat
      if(nlbper) then

!  Case of periodic b.c.

      ps(nlonp1,jlat) = ps(3,jlat)

      else

!  Linear extrapolation

      ps(nlonp1,jlat) = 2.*ps(nlon,jlat)-ps(nlonm1,jlat)
      endif
      enddo

      do jlat=0,nlat
      do jlon=1,nlonp1
      alps (jlon,jlat)=alog(ps(jlon,jlat))
      enddo
      enddo

!  Definition of pot. temp. (in t) and virtual temp.

      do jklev = 1, nlev
      do jlat = 1, nlat
      do jlon = 1, nlon
      zzpp=pzer*sigint(jklev)-(pzer-ps(jlon,jlat))*sigalf(jklev)
      t(jlon,jlat,jklev) = temp(jlon,jlat,jklev)*(1.e5/zzpp)**rcp
      zwat = max(0., qc(jlon,jlat,jklev))
      zq   = max(0., q(jlon,jlat,jklev))
      tvirt(jlon,jlat,jklev) = temp(jlon,jlat,jklev)*(1.+ep*zq-zwat)
      enddo
      enddo
      enddo

      return
      end
!#######################################################################
      subroutine phicomp

! Computes geopotential

      use mod_postbolam
      implicit none
      integer jlon, jlat, jklev
      real(4) zc1, zsigmed, zdsgalf, zsgalf

!  Phi is geop. at integer levels

      zsigmed=.5*(1.+sigint(nlev))
      zsgalf  = zsigmed**3 * (1.+alfa*(1.-zsigmed)*(2.-zsigmed))
      zdsgalf = zsigmed**2 * (3.+alfa*(6.-12.*zsigmed+5.*zsigmed**2))
      do jlat = 1, nlat
      do jlon = 1, nlon
      zc1 = ( pzer-(pzer-ps(jlon,jlat))*zdsgalf )/( pzer*zsigmed-(pzer-ps(jlon,jlat))*zsgalf )
      phi(jlon,jlat,nlev) = phig(jlon,jlat) + rd*(1.-sigint(nlev))*zc1*tvirt(jlon,jlat,nlev)
      enddo
      enddo

      do jklev = nlevm1, 1, -1
      zsigmed = .5*(sigint(jklev)+sigint(jklev+1))
      zsgalf  = zsigmed**3 * (1.+alfa*(1.-zsigmed)*(2.-zsigmed))
      zdsgalf = zsigmed**2 * (3.+alfa*(6.-12.*zsigmed+5.*zsigmed**2))
      do jlat = 1, nlat
      do jlon = 1, nlon
      zc1 = ( pzer-(pzer-ps(jlon,jlat))*zdsgalf )/( pzer*zsigmed-(pzer-ps(jlon,jlat))*zsgalf )
      phi(jlon,jlat,jklev) = phi(jlon,jlat,jklev+1) + .5*rd*(sigint(jklev+1)-sigint(jklev))*zc1* &
                             (tvirt(jlon,jlat,jklev)+tvirt(jlon,jlat,jklev+1))
      enddo
      enddo
      enddo

      return
      end subroutine phicomp
!#######################################################################
      subroutine diverg (nlwafp)

      use mod_postbolam

      real zdiv(nlon,nlev), zom1(nlon,nlev), zdvint(nlon,nlevp1)
      real zpbx(nlon), zpby(nlon), zpbyn(nlon)
      real zwork(nlonp1,nlat)

      logical nlwafp

      do jlat=1,nlat
      do jlon=1,nlon
      pst(jlon,jlat)=0.
      enddo
      enddo

      do 15 jlat = 2, nlatm1
      zrdx   = 1.         /(hxt(jlat)*dx)
      zhxvt  = hxv(jlat  )/(hxt(jlat)*dy)
      zhxvtn = hxv(jlat+1)/(hxt(jlat)*dy)

      do jlon = 1, nlonm1
      zpbx (jlon) = .5*zrdx  *(ps(jlon,jlat  )+ps(jlon+1,jlat  ))
      zpby (jlon) = .5*zhxvt *(ps(jlon,jlat  )+ps(jlon  ,jlat-1))
      zpbyn(jlon) = .5*zhxvtn*(ps(jlon,jlat+1)+ps(jlon  ,jlat  ))
      enddo

      do jklev = 1, nlev
      do jlon = 2, nlonm1
      zdi1(jlon,jlat,jklev)=(u(jlon,jlat,jklev)-u(jlon-1,jlat,jklev))*      &
                            zrdx+zhxvtn*v(jlon,jlat+1,jklev)-zhxvt*v(jlon,jlat,jklev)
      zdi2(jlon,jlat,jklev)=u(jlon,jlat,jklev)*zpbx(jlon)-u(jlon-1,jlat,jklev)*zpbx(jlon-1)+ &
                            v(jlon,jlat+1,jklev)*zpbyn(jlon)-v(jlon,jlat,jklev)*zpby(jlon)
      enddo
      enddo
 15   continue

! If WAF, zdi2 is redefined

      if(nlwafp) call wafps(ps,u,v,zdi2,nlon,nlat,nlev,nlonp1,hxt,hxv,dx,dy,dtstep,zwork)

! The following for smoothing vertical velocity

!      call divdam (zdi1, nlon, nlat, nlev, zwork)
!      call divdam (zdi2, nlon, nlat, nlev, zwork)

!  Definition of tendency of ps and of (-) integral from top to
!  a given half-level of horizontal divergence (zdvint)

      omeg(:,:,:)=0.
      do 20 jlat = 2, nlatm1

      do jklev = 1, nlev
      zdsgalf = sigint(jklev)**2 * (3.+alfa*(6.-12.*sigint(jklev)+5.*sigint(jklev)**2))
      do jlon = 2, nlonm1
      zdiv(jlon,jklev) = pzer*(1.-zdsgalf)*zdi1(jlon,jlat,jklev) + zdsgalf*zdi2(jlon,jlat,jklev)
      zom1(jlon,jklev) = zdi2(jlon,jlat,jklev)-ps(jlon,jlat)*zdi1(jlon,jlat,jklev)   ! to compute omega
      pst(jlon,jlat)=pst(jlon,jlat)-zdiv(jlon,jklev)*dsig(jklev)
      zdvint(jlon,jklev)=pst(jlon,jlat)
      enddo
      enddo

!  Omega (omeg) at t-points and integer levels

      do jklev = 1, nlev
      do jlon = 2, nlonm1
      omeg(jlon,jlat,jklev)= sigalf(jklev)*zom1(jlon,jklev) +             &
                zdvint(jlon,jklev) + zdiv(jlon,jklev)*(sig(jklev)-sigint(jklev))
      enddo
      enddo

 20   continue

      return
      end
!#######################################################################
 module mod_rwcross
   integer, parameter :: ncrossx = 3, ncrossy = 3, npress = 100
   integer :: lncr(ncrossy), ltcr(ncrossx)
   real :: zplcr(npress)
   real, dimension(:,:), allocatable :: zspcr, zspcrx
   real, dimension(:,:,:), allocatable :: ztcr, zrhcr, zucr, zvcr, ztecr, zwcr, zpvcr, ztcrx, zrhcrx, zucrx, &
                                          zvcrx, ztecrx, zwcrx, zpvcrx
   real, dimension(:,:), allocatable :: lon_cross, lat_cross, ps_cross
   real, dimension(:,:,:), allocatable :: theta_cross, thetae_cross, u_cross, v_cross, u_tang_cross, v_norm_cross, &
 w_cross, pv_cross, rh_cross
 end module mod_rwcross

 module mod_rtlev
   use mod_postbolam, only : nlevpo, nlevto
   real tlevo(nlevto) ! per sup. theta
 end module mod_rtlev

!  Module containing vectors (for surface fields) or vertical slabs
!  (for fields at pre-selected pressure levels, defined in plevo)

 module mod_rwoutp
   use mod_postbolam, only : nlevpo, nlevto
   real, dimension(:), allocatable :: zvtz0, zvpz0, zvts, zvus, zvvs, zvqs
   real, dimension(:,:), allocatable :: zvt, zvph, zvu, zvv, zvq, zvom, zvpv, zvthetae
 end module mod_rwoutp

 module mod_rwoutt
   real, dimension(:,:), allocatable :: zvpvt ! for theta surfaces
 end module mod_rwoutt

!  Potential vorticity

 module mod_ptvor
   real, dimension(:,:,:), allocatable :: potv
   real, dimension(:,:), allocatable :: potvn
 end module mod_ptvor
!#######################################################################
      subroutine postpcpm(kist,prtop,iterf)

!  Prepares 2-d and 3-d fields (on horiz. surfaces and cross-sections)
!  Called routines specific for postprocessing are:
!  defout, interp, wrpost, wrposx, wrposy, bextr, potvor.
!  Nlevpo is the number of output standard pressure levels where model
!  fields are interpolated (surface values apart).
!  Specific values of pressure for these levels are defined below.
!  Levels and fields are selected using the flag matrix ipflag.

 use mod_postbolam
 use mod_surfadd
 use mod_cross
 use mod_rwcross
 use mod_rplev
 use mod_rtlev
 use mod_rwoutp
 use mod_rwoutt
 use mod_ptvor
 use mod_jump

! ---- For RTTOV simulation ---

#ifdef rttov
 use mod_rttov
#endif

! -----------------------------
!=======================================================================
!  Cross-sections:
!  ncrossy: number of north-south cross-sections
!  ncrossx: number of west-east cross-sections
!  Cross_sections can be only along meridians and parallels.
!  Npress: number of pressure level for cross-sections (arbitrary)
!  (to be changed here and in subr. defout)
!  Loncr1, loncr2 ...: grid points defining longitudes where N-S cross sections are defined.
!  The number of these param. must be equal to ncrossy - same for longitude.

!=======================================================================

!  Horiz. lon-lat matrixes containing output surface fields

  real, dimension(:,:), allocatable :: ztz0, zpz0, ztgrs, ztmin, ztmax, zqgrs, ztgr1, ztgr2, ztgr3, ztgr4, zqgr1, &
                                       zqgr2, zqgr3, zqgr4, zts, zus, zvs, zqs
  real, dimension(:,:,:), allocatable :: tg_post, qg_post

!  3-d matrixes containing output fields on selected pressure levels

  real, dimension(:,:,:), allocatable :: zt, zph, zu, zv, zq, zom, zpv, zrv, zrh, zthetae
  real, dimension(:,:,:), allocatable :: zpvt   ! for theta surfaces

!  Surface pressure

  real, dimension(:,:), allocatable :: zsp

! Orography and CAPE

  real, dimension(:,:), allocatable :: zorogr, zcape, zcin, zlift, zlev0
  real, dimension(:), allocatable :: ztvp, zpout

!  Auxiliary vectors and matrixes used for interpolation on staggered grid

  real, dimension(:), allocatable :: zaxus, zaxvs
  real, dimension(:,:), allocatable :: zaxu, zaxv
  real, dimension(:,:,:), allocatable :: zpva, zrva
  real, dimension(:,:,:), allocatable :: zpvat ! per sup. theta
  real, dimension(:,:), allocatable :: zwork, zwor1, zwor2

  integer, dimension(5) :: idate0, idatec
  integer, dimension(3) :: iperiod=0, iperiod_accum=0
  integer, dimension(12) ::imon=(/31,28,31,30,31,30,31,31,30,31,30,31/)
  integer :: nday, ndm, iimon
  real, dimension(:,:), allocatable :: lon_crossy, lat_crossy, lon_crossx, lat_crossx
  real, dimension(1) :: zsq
  real, dimension(2) :: zpavup, ztvlift
  real :: alon0, alat0, dlon, dlat, zlat, coeday, lapse_rate, lapse_rate_clim
  real, dimension(nlev) :: zvgrt, zh, zspeed
  integer, dimension(1) :: kmax

! Definition of the number of "free" cross-sections

  if (ncrossec) then

    do icross = 1,10
      if (lonini_cross(icross) == valmiss.or.latini_cross(icross) == valmiss.or.     &
           lonfin_cross(icross) == valmiss.and.latfin_cross(icross) == valmiss) exit
      cross_number = cross_number+1
      if (npoint_cross(icross) == 0) then
        dist_deg=sqrt((lonini_cross(icross)-lonfin_cross(icross))**2 + &
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

  endif

! Allocation of dinamic arrays

  allocate(zspcr(nlat,ncrossy), stat=ierr)
  allocate(zspcrx(nlon,ncrossx), stat=ierr)
  allocate(ztcr(nlat,npress,ncrossy), stat=ierr)
  allocate(zrhcr(nlat,npress,ncrossy), stat=ierr)
  allocate(zucr(nlat,npress,ncrossy), stat=ierr)
  allocate(zvcr(nlat,npress,ncrossy), stat=ierr)
  allocate(ztecr(nlat,npress,ncrossy), stat=ierr)
  allocate(zwcr(nlat,npress,ncrossy), stat=ierr)
  allocate(zpvcr(nlat,npress,ncrossy), stat=ierr)
  allocate(ztcrx(nlon,npress,ncrossx), stat=ierr)
  allocate(zrhcrx(nlon,npress,ncrossx), stat=ierr)
  allocate(zucrx(nlon,npress,ncrossx), stat=ierr)
  allocate(zvcrx(nlon,npress,ncrossx), stat=ierr)
  allocate(ztecrx(nlon,npress,ncrossx), stat=ierr)
  allocate(zwcrx(nlon,npress,ncrossx), stat=ierr)
  allocate(zpvcrx(nlon,npress,ncrossx), stat=ierr)

  allocate(lon_cross(npoint_cross_max,cross_number), stat=ierr)
  allocate(lat_cross(npoint_cross_max,cross_number), stat=ierr)
  allocate(ps_cross(npoint_cross_max,cross_number), stat=ierr)
  allocate(theta_cross(npoint_cross_max,npress,cross_number), stat=ierr)
  allocate(thetae_cross(npoint_cross_max,npress,cross_number), stat=ierr)
  allocate(u_cross(npoint_cross_max,npress,cross_number), stat=ierr)
  allocate(v_cross(npoint_cross_max,npress,cross_number), stat=ierr)
  allocate(u_tang_cross(npoint_cross_max,npress,cross_number), stat=ierr)
  allocate(v_norm_cross(npoint_cross_max,npress,cross_number), stat=ierr)
  allocate(w_cross(npoint_cross_max,npress,cross_number), stat=ierr)
  allocate(pv_cross(npoint_cross_max,npress,cross_number), stat=ierr)
  allocate(rh_cross(npoint_cross_max,npress,cross_number), stat=ierr)

  allocate(zvtz0(nlon), stat=ierr)
  allocate(zvpz0(nlon), stat=ierr)
  allocate(zvts(nlon), stat=ierr)
  allocate(zvus(nlon), stat=ierr)
  allocate(zvvs(nlon), stat=ierr)
  allocate(zvqs(nlon), stat=ierr)
  allocate(zvt(nlon,nlevpo), stat=ierr)
  allocate(zvph(nlon,nlevpo), stat=ierr)
  allocate(zvu(nlon,nlevpo), stat=ierr)
  allocate(zvv(nlon,nlevpo), stat=ierr)
  allocate(zvq(nlon,nlevpo), stat=ierr)
  allocate(zvom(nlon,nlevpo), stat=ierr)
  allocate(zvpv(nlon,nlevpo), stat=ierr)
  allocate(zvthetae(nlon,nlevpo), stat=ierr)
  allocate(zvpvt(nlon,nlevto), stat=ierr)

  allocate(potv(nlon,nlat,nlev), stat=ierr)
  allocate(potvn(nlon,nlev), stat=ierr)

  allocate(ztz0(nlon,nlat), stat=ierr)
  allocate(zpz0(nlon,nlat), stat=ierr)
  allocate(ztgrs(nlon,nlat), stat=ierr)
  allocate(ztmin(nlon,nlat), stat=ierr)
  allocate(ztmax(nlon,nlat), stat=ierr)
  allocate(zqgrs(nlon,nlat), stat=ierr)
  allocate(ztgr1(nlon,nlat), stat=ierr)
  allocate(ztgr2(nlon,nlat), stat=ierr)
  allocate(ztgr3(nlon,nlat), stat=ierr)
  allocate(ztgr4(nlon,nlat), stat=ierr)
  allocate(zqgr1(nlon,nlat), stat=ierr)
  allocate(zqgr2(nlon,nlat), stat=ierr)
  allocate(zqgr3(nlon,nlat), stat=ierr)
  allocate(zqgr4(nlon,nlat), stat=ierr)
  allocate(zts(nlon,nlat), stat=ierr)
  allocate(zus(nlon,nlat), stat=ierr)
  allocate(zvs(nlon,nlat), stat=ierr)
  allocate(zqs(nlon,nlat), stat=ierr)
  allocate(tg_post(nlon,nlat,nlevg), stat=ierr)
  allocate(qg_post(nlon,nlat,nlevg), stat=ierr)

  allocate(zt(nlon,nlat,nlevpo), stat=ierr)
  allocate(zph(nlon,nlat,nlevpo), stat=ierr)
  allocate(zu(nlon,nlat,nlevpo), stat=ierr)
  allocate(zv(nlon,nlat,nlevpo), stat=ierr)
  allocate(zq(nlon,nlat,nlevpo), stat=ierr)
  allocate(zom(nlon,nlat,nlevpo), stat=ierr)
  allocate(zpv(nlon,nlat,nlevpo), stat=ierr)
  allocate(zrv(nlon,nlat,nlevpo), stat=ierr)
  allocate(zrh(nlon,nlat,nlevpo), stat=ierr)
  allocate(zthetae(nlon,nlat,nlevpo), stat=ierr)

  allocate(zpvt(nlon,nlat,nlevpo), stat=ierr)

  allocate(zsp(nlon,nlat), stat=ierr)

  allocate(zorogr(nlon,nlat), stat=ierr)
  allocate(zcape(nlon,nlat), stat=ierr)
  allocate(zcin(nlon,nlat), stat=ierr)
  allocate(zlift(nlon,nlat), stat=ierr)
  allocate(zlev0(nlon,nlat), stat=ierr)
  allocate(ztvp(nlev), stat=ierr)
  allocate(zpout(nlev), stat=ierr)

  allocate(zaxus(nlat), stat=ierr)
  allocate(zaxu(nlat,nlevpo), stat=ierr)
  allocate(zaxvs(nlon), stat=ierr)
  allocate(zaxv(nlon,nlevpo), stat=ierr)
  allocate(zpva(nlon,nlat,nlevpo), stat=ierr)
  allocate(zrva(nlon,nlat,nlevpo), stat=ierr)
  allocate(zwork(nlon,nlat), stat=ierr)
  allocate(zwor1(nlon,nlat), stat=ierr)
  allocate(zwor2(nlon,nlat), stat=ierr)
  allocate(zpvat(nlon,nlat,nlevto), stat=ierr)

  allocate(lon_crossy(ncrossy,2), stat=ierr)
  allocate(lat_crossy(ncrossy,2), stat=ierr)
  allocate(lon_crossx(ncrossx,2), stat=ierr)
  allocate(lat_crossx(ncrossx,2), stat=ierr)
!-----------------------------------------------------------------------------
!     Initialisation of tlevo vector ! for theta surfaces

      tlevo(1) = 315
      tlevo(2) = 320
      tlevo(3) = 325
      tlevo(4) = 330

!=======================================================================
!  Definition of longitude and latitude for cross sections

      lncr(1) = loncr1
      lncr(2) = loncr2
      lncr(3) = loncr3
      ltcr(1) = latcr1
      ltcr(2) = latcr2
      ltcr(3) = latcr3

! Redefinition of longitudes and latitudes for cross sections
! Caution: redefines values defined in postbolam.inp file

      if(lncr(3).gt.nlon) then
      lncr(1) = nlon/4
      lncr(2) = nlon/2
      lncr(3) = 3*nlon/4
      print*, "Longitudes of cross-sections redefined because some values defined in file"
      print*, "postbolam.inp are out of the model domain"
      endif
      if(ltcr(3).gt.nlat) then
      ltcr(1) = nlat/4
      ltcr(2) = nlat/2
      ltcr(3) = 3*nlat/4
      print*, "Latitudes of cross-sections redefined because some values defined in file"
      print*, "postbolam.inp are out of the model domain"
      endif

!  Definition of values of log(p) for levels of cross-sections -
!  pl1 and pl2 are bottom and top of cross-sections

      zpl1 = alog(104000.)
      zpl2 = alog(prtop*100.)

      zdpl = (zpl2-zpl1)/float(npress-1)
      do jk = 1, npress
      zplcr(jk) = zpl1 + zdpl*float(jk-1)
      enddo
!=======================================================================

!  Resetting of negative values of q2

      do jlat = 1, nlat
      do jlon = 1, nlon
      q2(jlon,jlat) = max(q2(jlon,jlat), 1.e-8)
      enddo
      enddo

!  Definition of orography

      do jlat = 1, nlat
      do jlon = 1, nlon
      zorogr(jlon,jlat) = phig(jlon,jlat)/g
      enddo
      enddo

!  Definition of CAPE: it is assumed that air parcel does not contain liquid water at start

      if (isflag(24).eq.1) then
      do jlat = 1, nlat
      do jlon = 1, nlon

      do jklev = nlev, 2, -1
      zpout(jklev) = pzer*sigint(jklev)-(pzer-ps(jlon,jlat))*sigalf(jklev)
      enddo
      zcape(jlon,jlat) = 0.
      zcin(jlon,jlat) = 0.

      do jk0 = nlev, nlev-2, -1
      zq0p = q(jlon,jlat,jk0)
      zt0p = temp(jlon,jlat,jk0)
      call lift_parcel_bol (zq0p, zt0p, zpout, ztvp, nlev, jk0, 7, 0)

      zcap = 0.
      do jklev = jk0, 8, -1
      zcap1 = .5*rd*(zpout(jklev)-zpout(jklev-1))*                        &
              ((ztvp(jklev  )-tvirt(jlon,jlat,jklev  ))/zpout(jklev  )+   &
               (ztvp(jklev-1)-tvirt(jlon,jlat,jklev-1))/zpout(jklev-1))
      zcap = zcap + max(zcap1,0.)  ! elimina il contributo se < 0 (cin)
      if (zcap < 20.) then   ! important - replaces condition on level of free convection
        zci = zci + min(zcap1,0.)
      endif
      enddo

      zcape(jlon,jlat) = max(zcape(jlon,jlat),zcap)
      zcin(jlon,jlat)  = min(zcin(jlon,jlat), zci)
      enddo

      if (zcape(jlon,jlat) <   50.) zcin(jlon,jlat) = -9999.
      if (zcin (jlon,jlat) < -400.) zcin(jlon,jlat) = -9999.

      enddo
      enddo
      endif

! Integrated Water Vapour IWV

      do jlat = 1, nlat
      do jlon = 1, nlon
       iwv (jlon,jlat) = 0.
       do jk=1,nlev
       iwv(jlon,jlat)=iwv(jlon,jlat) + q(jlon,jlat,jk)/g*( pzer*dsig(jk)-(pzer-ps(jlon,jlat))*dsigalf(jk) )
       enddo
      enddo
      enddo

! Definition of lifted index: virtual temp. difference between a parcel near about 500 hPa
! (averaged between two model levels, one above and one below 500 hPa) and a parcel
! lifted pseudo-adiabatically, starting in a layer of depth of ZPDEPTH above the surface.
! Since ps is variable over orography, the 500 hPa level is actually arbitrarily defined as
! a pressure 500.E2 - (1.E5-PS)/4.  ZPDEPTH also depends on the surface pressure.

      if(isflag(50).eq.1) then
      do jlat = 1, nlat
      do jlon = 1, nlon
      zlift(jlon,jlat)=999.

! Average of pot. temp. and specif. humidity in the pressure layer above the ground
! (weighting with pressure thickness of model layers)

     jklol=-1
     jkupl=-1
     do jk=nlev,1,-1
     zp0p=pzer*sigint(jk)-(pzer-ps(jlon,jlat))*sigalf(jk)
     zpdepth=90.e2*sqrt(ps(jlon,jlat)/1.e5)  ! decreases the press. depth of mixing lower layer over topography
     if(zp0p.lt.ps(jlon,jlat)-zpdepth.and.jklol.lt.0) jklol=jk+1 ! index of the top of the lower layer
     if(zp0p.lt.(500.e2-(1.e5-ps(jlon,jlat))/4.).and.jkupl.lt.0) then
       jkupl=jk ! index of the first layer above the upper (near 500 hpa) level
       goto 7890
       endif
     enddo

! Average of temp., spec. hum. and press. in the lower layer

7890 zpavlo=0.
     ztempavlo=0.
     zqavlo=0.
     zdz=0.
     do jk=jklol,nlev
     zmss=pzer*dsig(jk)-(pzer-ps(jlon,jlat))*dsigalf(jk)
     zpavlo=zpavlo+(pzer*sigint(jk)-(pzer-ps(jlon,jlat))*sigalf(jk))*zmss
     ztempavlo=ztempavlo+temp(jlon,jlat,jk)*zmss
     zqavlo=zqavlo+q(jlon,jlat,jk)*zmss
     zdz=zdz+zmss
     enddo
     zpavlo=zpavlo/zdz
     ztempavlo=ztempavlo/zdz
     zqavlo=zqavlo/zdz

! Average of virtual temperature and press. in the upper layer

     ztvirtup = 0.5*(tvirt(jlon,jlat,jkupl)+tvirt(jlon,jlat,jkupl+1))
     zpavup(2) = zpavlo
     zpavup(1) = 0.5*(pzer*sigint(jkupl)-(pzer-ps(jlon,jlat))*sigalf(jkupl)+  &
                 pzer*sigint(jkupl+1)-(pzer-ps(jlon,jlat))*sigalf(jkupl+1))
      call lift_parcel_bol (zqavlo, ztempavlo, zpavup, ztvlift, 2, 2, 1, 0)
      zlift(jlon,jlat)=ztvlift(1)-ztvirtup
      enddo
      enddo
      endif

!  End of computation of lifted index

! Definition of level of 0°C isotherm - computed from above, so in cases
! of multiple levels, the highest is considered (because defines melting level)

      if(isflag(51).eq.1) then
      do jlat = 1, nlat
      do jlon = 1, nlon
      zlev0(jlon,jlat)=-99999.
      ztempab=-99999.
      do jklev=10,nlev
      if (jklev>10) ztempab=temp(jlon,jlat,jklev-1)
      if(ztempab.le.ttr.and.temp(jlon,jlat,jklev).gt.ttr) then  ! linear extrapolation down
      zlev0(jlon,jlat)=(phi(jlon,jlat,jklev)-phi(jlon,jlat,jklev-1))/g/(temp(jlon,jlat,jklev)-ztempab)* &
                      (ttr-ztempab) + phi(jlon,jlat,jklev-1)/g
      goto 7910
      endif
      enddo
7910  continue

! Case in which the level of 0 is below ground: extrap. from second level with standard lapse rate

      if(zlev0(jlon,jlat).eq.-99999.) then
      zlev0(jlon,jlat)=phi(jlon,jlat,nlev-1)/g - (ttr-temp(jlon,jlat,nlev-1))/6.5e-3
      endif
      enddo
      enddo
      endif

      call divdam(zlev0, nlon, nlat, 1, zwork) ! Filter to smooth discontinuities at fronts or intersection with ground
      zlev0 = max(zlev0, -50.)  ! avoids plotting values << 0.

!  Loop in latitude

      do jlat = 1, nlat

!  Def. of pot. vorticity in the north slab (to compute pv on t-points)
!  Pot. vorticity is not defined on eastern and southern boundary

      if(jlat.eq.1) then
        potv(:,jlat,:)=0.
        potvn(:,:)=0.
        call potvor(2,potv(:,jlat,:))
      else
        do jlon=1,nlon
        do jklev=1,nlev
        potv(jlon,jlat,jklev)=potvn(jlon,jklev)
        enddo
        enddo
      endif

      if(jlat.lt.nlat) call potvor(jlat+1,potvn)

      zvt(:,:)=0.
      zvph(:,:)=0.
      zvu(:,:)=0.
      zvv(:,:)=0.
      zvq(:,:)=0.
      zvom(:,:)=0.
      zvpv(:,:)=0.
      zvthetae(:,:)=0.
      zvpvt(:,:)=0.
      call defout(jlat)
      if (jlat.gt.1) call thetdef(jlat)

!  Def. of horiz. matrixes from vectors of output variables

      do jlon = 1, nlon
      ztz0(jlon,jlat)  = zvtz0(jlon)
      zpz0(jlon,jlat)  = zvpz0(jlon)
      zts (jlon,jlat)  = zvts (jlon)
      zus (jlon,jlat)  = zvus (jlon)
      zvs (jlon,jlat)  = zvvs (jlon)
      zqs (jlon,jlat)  = zvqs (jlon)
      zsp(jlon,jlat)   = ps(jlon,jlat)

        do jlev = 1, nlevpo
        zt (jlon,jlat,jlev) = zvt (jlon,jlev)
        zph(jlon,jlat,jlev) = zvph(jlon,jlev)
        zu (jlon,jlat,jlev) = zvu (jlon,jlev)
        zv (jlon,jlat,jlev) = zvv (jlon,jlev)
        zq (jlon,jlat,jlev) = zvq (jlon,jlev)
        zom(jlon,jlat,jlev) = zvom(jlon,jlev)
        zpv(jlon,jlat,jlev) = zvpv(jlon,jlev)
        zthetae(jlon,jlat,jlev) = zvthetae(jlon,jlev)
        enddo
        do jlev = 1, nlevto                    ! for theta surf.
        zpvt(jlon,jlat,jlev)= zvpvt(jlon,jlev)
        enddo
      enddo
      enddo ! loop in lat.

!  Computation of relative vorticity (zrv) at z points
!  (undefined along south and east boundaries)
!  zrv is set to zero if any of the surrounding t points is
!  below ground, to avoid contamination from winds below ground

      do jlon = 1, nlonm1
      do jlat = 2, nlat
      do jlev = 1, nlevpo
      if (plevo(jlev).gt.zsp(jlon,jlat).or.plevo(jlev).gt.zsp(jlon+1,jlat).or.   &
        plevo(jlev).gt.zsp(jlon,jlat-1).or.plevo(jlev).gt.zsp(jlon+1,jlat-1)) then
        zrv(jlon,jlat,jlev)=0.
      else
        zrv (jlon,jlat,jlev) = (-(zu(jlon,jlat,jlev)*hxt(jlat)-zu(jlon,jlat-1,jlev)*hxt(jlat-1))/dy+ &
        (zv(jlon+1,jlat,jlev)-zv(jlon,jlat,jlev))/dx)/(hxt(jlat)+hxt(jlat-1))*2.
      endif
      enddo
      enddo
      enddo

!  Lin. extrap. at lateral bound. of missing variables (see grid struct.)
!  for variables missing only at some of the bound. lines: pv, rv
!  East

      do jlat = 2, nlat
      do jlev = 1, nlevpo
      zpv(nlon,jlat,jlev) = 2.*zpv(nlonm1,jlat,jlev) - zpv(nlon-2,jlat,jlev)
      zrv(nlon,jlat,jlev) = 2.*zrv(nlonm1,jlat,jlev) - zrv(nlon-2,jlat,jlev)
      enddo
       do jlev = 1, nlevto         ! per sup. theta
        zpvt(nlon,jlat,jlev) = 2.*zpvt(nlonm1,jlat,jlev) - zpvt(nlon-2,jlat,jlev)
       enddo
      enddo

      do jlon = 1, nlon

!  South

      do jlev = 1, nlevpo
      zpv(jlon,1,jlev) = 2.*zpv(jlon,2,jlev) - zpv(jlon,3,jlev)
      zrv(jlon,1,jlev) = 2.*zrv(jlon,2,jlev) - zrv(jlon,3,jlev)
      enddo

      do jlev = 1, nlevto                  ! per sup. theta
      zpvt(jlon,1,jlev) = 2.*zpvt(jlon,2,jlev)-zpvt(jlon,3,jlev)
      enddo

      enddo

!  Extrapolation of omega at all lateral boundaries

      zeex = .5
      do jlev = 1, nlevpo
      call bextr(zom(1,1,jlev),nlon,nlat,zeex,nlbper)
      enddo

!=======================================================================
!  Boundary extrapolation of vertical velocity

      do k=1,npress
        do n=1,ncrossy
        zwcr(1   ,k,n) = 1.5*zwcr(2     ,k,n)-0.5*zwcr(3     ,k,n)
        zwcr(nlat,k,n) = 1.5*zwcr(nlat-1,k,n)-0.5*zwcr(nlat-2,k,n)
        enddo
        do n=1,ncrossx
        zwcrx(1   ,k,n) = 1.5*zwcrx(2     ,k,n)-0.5*zwcrx(3     ,k,n)
        zwcrx(nlon,k,n) = 1.5*zwcrx(nlon-1,k,n)-0.5*zwcrx(nlon-2,k,n)
        enddo
      enddo

! Definition of the level of Tropopause

      if (any(itropflag(:) /= 0)) then

      do jlat = 1, nlat
      do jlon = 1, nlon

      do iteration = 1,2

        if (iteration == 1) then
          hlim_down_1 = 1000.
          hlim_up_1 = 8000.
          hlim_down_2 = 15000.
          hlim_up_2 = 25000.
        else
          hlim_down_1 = 1000.
          hlim_up_1 = max(h_tropopause(jlon,jlat)-500., 4000.)
          hlim_down_2 = h_tropopause(jlon,jlat)+500.
          hlim_up_2 = h_tropopause(jlon,jlat)+10000.
        endif

        do jklev = 1, nlev
          zh(jklev) = phi(jlon,jlat,jklev)/g
        enddo

        do jklev = 1,nlev-1
          zvgrt(jklev) = (temp(jlon,jlat,jklev)-temp(jlon,jlat,jklev+1))/(zh(jklev)-zh(jklev+1))
        enddo
        zvgrt(nlev) = -0.01

! Method of regressions

        n = 0
        sum1 = 0.
        sum2 = 0.
        sum3 = 0.
        sum4 = 0.
        do jklev = nlev, 1, -1
          if ( (phi(jlon,jlat,jklev)-phig(jlon,jlat))/g > hlim_down_1.and.zh(jklev) < hlim_up_1) then
            n = n+1
            sum1 = sum1+zh(jklev)
            sum2 = sum2+temp(jlon,jlat,jklev)
            sum3 = sum3+zh(jklev)**2
            sum4 = sum4+zh(jklev)*temp(jlon,jlat,jklev)
          endif
        enddo
        if (n > 0) then
          a1 = (sum4-sum1*sum2/n)/(sum3-sum1**2/n)
          b1 = (sum2-a1*sum1)/n
        else
          a1 = float(ivalmiss)
          b1 = float(ivalmiss)
        endif

        n = 0
        sum1 = 0.
        sum2 = 0.
        sum3 = 0.
        sum4 = 0.
        do jklev = nlev, 1, -1
          if ( zh(jklev) > hlim_down_2.and.zh(jklev) < hlim_up_2) then
            n = n+1
            sum1 = sum1+zh(jklev)
            sum2 = sum2+temp(jlon,jlat,jklev)
            sum3 = sum3+zh(jklev)**2
            sum4 = sum4+zh(jklev)*temp(jlon,jlat,jklev)
          endif
        enddo
        if (n > 0) then
          a2 = (sum4-sum1*sum2/n)/(sum3-sum1**2/n)
          b2 = (sum2-a2*sum1)/n
        else
          a2 = float(ivalmiss)
          b2 = float(ivalmiss)
        endif

        iflag = 0
        if (int(a1) /= ivalmiss.and.int(a2) /= ivalmiss) then

          h_tropopause(jlon,jlat) = (b2-b1)/(a1-a2)

        else
          iflag = 1
        endif

        if (h_tropopause(jlon,jlat) < 7000.or.h_tropopause(jlon,jlat) > 20000.) iflag = 1

        if (iflag == 1) exit

      enddo ! iteration

iflag=1 !!!!

        if (iflag == 1) then

! Method of vertical gradient of temperature

          do jklev = 1,nlev
            if (zh(jklev) < 20000.) exit
          enddo
          lev_ini=max(jklev-1, 1)
          do jklev = nlev, 1, -1
            if (zh(jklev) >  7000.) exit
          enddo
          lev_fin=min(jklev+1, nlev)

          h_tropopause(jlon,jlat) = 10000.
          do jklev = lev_fin, lev_ini, -1
            if (zvgrt(jklev) > -2.e-3) then
              h_tropopause(jlon,jlat) = zh(jklev)
              exit
            endif
          enddo

        endif

      enddo
      enddo

      call smooth(h_tropopause, zwork, nlon, nlat, 0.5, 5)
      h_tropopause(:,:) = zwork(:,:)

! Vertical interpolation of parameters to tropopause level

      call defout_tropopause

      endif ! Tropopause definition

! Definition of the level of Maximum wind

      if (any(imwflag(:) /= 0)) then

      do jlat = 1, nlat
      do jlon = 1, nlon

        zspeed(:) = 0.
        do jklev = 3, nlev
          if (jlon /= 1.and.jlat /= nlat) then
            zspeed(jklev) = (0.5*(u(jlon,jlat,jklev)+u(jlon-1,jlat,jklev)))**2+(0.5*(v(jlon,jlat,jklev)+v(jlon,jlat+1,jklev)))**2
          elseif (jlon == 1.and.jlat /= nlat) then
            zspeed(jklev) = u(jlon,jlat,jklev)**2+(0.5*(v(jlon,jlat,jklev)+v(jlon,jlat+1,jklev)))**2
          elseif (jlon /= 1.and.jlat == nlat) then
            zspeed(jklev) = (0.5*(u(jlon,jlat,jklev)+u(jlon-1,jlat,jklev)))**2+v(jlon,jlat,jklev)**2
          else
            zspeed(jklev) = u(jlon,jlat,jklev)**2+v(jlon,jlat,jklev)**2
          endif
        enddo

        kmax = maxloc(zspeed)
        klev_maxwind(jlon,jlat) = kmax(1)

      enddo
      enddo

! Definition of parameters at the level of maximum wind

      call defout_maxwind

      endif ! Maximum wind level definition

!  Final rescaling of rel. hum. (from 0 to 100) and pv (to pv units)
!  Reset of unphysical values of theta, thetae (may be due to
!  extrapolation below ground) and of relative humidity.

      do k=1,npress
        do n=1,ncrossy
          do jlat=1,nlat
          zrhcr(jlat,k,n) = zrhcr(jlat,k,n)*100.
          zrhcr(jlat,k,n) = max(zrhcr(jlat,k,n),0.)
          zpvcr(jlat,k,n) = zpvcr(jlat,k,n)*1.e6
          ztcr (jlat,k,n) = max(ztcr (jlat,k,n),150.)
          ztecr(jlat,k,n) = max(ztecr(jlat,k,n),150.)
          enddo
        enddo
        do n=1,ncrossx
          do jlon=1,nlon
          zrhcrx(jlon,k,n) = zrhcrx(jlon,k,n)*100.
          zrhcrx(jlon,k,n) = max(zrhcrx(jlon,k,n),0.)
          zpvcrx(jlon,k,n) = zpvcrx(jlon,k,n)*1.e6
          ztcrx (jlon,k,n) = max(ztcrx (jlon,k,n),150.)
          ztecrx(jlon,k,n) = max(ztecrx(jlon,k,n),150.)
          enddo
        enddo
      enddo

!=======================================================================

   if (cross_number > 0) call free_cross_section

!=======================================================================

!  Def.of ground temperatures, tmin and tmax

      do jlat = 1, nlat
      do jlon = 1, nlon
      ztgrs(jlon,jlat) = tskin(jlon,jlat)- ttr
      ztgr1(jlon,jlat) = tg(jlon,jlat,1) - ttr
!      ztgr2(jlon,jlat) = tg(jlon,jlat,2) - ttr
!      ztgr3(jlon,jlat) = tg(jlon,jlat,3) - ttr
!      ztgr4(jlon,jlat) = tg(jlon,jlat,4) - ttr
      ztgr2(jlon,jlat) = tg(jlon,jlat,3) - ttr
      ztgr3(jlon,jlat) = tg(jlon,jlat,5) - ttr
      ztgr4(jlon,jlat) = tg(jlon,jlat,7) - ttr
      ztmin(jlon,jlat) = t2min(jlon,jlat)- ttr
      ztmax(jlon,jlat) = t2max(jlon,jlat)- ttr
      zts(jlon,jlat)   = zts(jlon,jlat)  - ttr
      t05(jlon,jlat)   = t05(jlon,jlat)  - ttr
      t005(jlon,jlat)  = t005(jlon,jlat) - ttr
      enddo
      enddo

      tg_post(:,:,:) = tg(:,:,:) - ttr
      qg_post(:,:,:) = qg(:,:,:)

!  Defin. of ground water

      do jlat = 1, nlat
      do jlon = 1, nlon
      zqgrs(jlon,jlat) = qskin(jlon,jlat)
      zqgr1(jlon,jlat) = qg(jlon,jlat,1)
!      zqgr2(jlon,jlat) = qg(jlon,jlat,2)
!      zqgr3(jlon,jlat) = qg(jlon,jlat,3)
!      zqgr4(jlon,jlat) = qg(jlon,jlat,4)
      zqgr2(jlon,jlat) = qg(jlon,jlat,3)
      zqgr3(jlon,jlat) = qg(jlon,jlat,5)
      zqgr4(jlon,jlat) = qg(jlon,jlat,7)
      enddo
      enddo

!  Ground temperature

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
      do ig=1,nlevg
      tg(jlon,jlat,ig)=tg(jlon,jlat,ig)-ttr
      enddo
      tg005(jlon,jlat)=tg005(jlon,jlat)-ttr
      tg010(jlon,jlat)=tg010(jlon,jlat)-ttr
      tg020(jlon,jlat)=tg020(jlon,jlat)-ttr
      tg050(jlon,jlat)=tg050(jlon,jlat)-ttr
      tg100(jlon,jlat)=tg100(jlon,jlat)-ttr
      enddo
      enddo

!  Surface fluxes of sensib. and latent heat (w/m**2), accumul. prec., clouds, 2 m temp., humidity etc.

      zeex = .5
      call bextr(hflux ,nlon,nlat,zeex,nlbper)
      call bextr(qflux ,nlon,nlat,zeex,nlbper)
      call bextr(totpre,nlon,nlat,zeex,nlbper)
      call bextr(conpre,nlon,nlat,zeex,nlbper)
      call bextr(snfall,nlon,nlat,zeex,nlbper)
      call bextr(snow,  nlon,nlat,zeex,nlbper)
      call bextr(runoff,nlon,nlat,zeex,nlbper)
      call bextr(runoff_tot,nlon,nlat,zeex,nlbper)
      call bextr(cloud ,nlon,nlat,zeex,nlbper)
      call bextr(cloudt,nlon,nlat,zeex,nlbper)
      call bextr(cloudh,nlon,nlat,zeex,nlbper)
      call bextr(cloudm,nlon,nlat,zeex,nlbper)
      call bextr(cloudl,nlon,nlat,zeex,nlbper)
      call bextr(t2,    nlon,nlat,zeex,nlbper)
      call bextr(q2,    nlon,nlat,zeex,nlbper)
      call bextr(q2rel, nlon,nlat,zeex,nlbper)
      call bextr(td2,   nlon,nlat,zeex,nlbper)
      call bextr(zqgrs, nlon,nlat,zeex,nlbper)
      call bextr(ztgrs, nlon,nlat,zeex,nlbper)
      call bextr(ztgr1, nlon,nlat,zeex,nlbper)
      call bextr(ztgr2, nlon,nlat,zeex,nlbper)
      call bextr(ztgr3, nlon,nlat,zeex,nlbper)
      call bextr(ztgr4, nlon,nlat,zeex,nlbper)
      call bextr(ztmin, nlon,nlat,zeex,nlbper)
      call bextr(ztmax, nlon,nlat,zeex,nlbper)
      call bextr(zts,   nlon,nlat,zeex,nlbper)
      call bextr(t05,   nlon,nlat,zeex,nlbper)
      call bextr(q05rel,nlon,nlat,zeex,nlbper)
      call bextr(t005,  nlon,nlat,zeex,nlbper)
      call bextr(tg005, nlon,nlat,zeex,nlbper)
      call bextr(tg010, nlon,nlat,zeex,nlbper)
      call bextr(tg020, nlon,nlat,zeex,nlbper)
      call bextr(tg050, nlon,nlat,zeex,nlbper)
      call bextr(tg100, nlon,nlat,zeex,nlbper)
      call bextr(cswfl, nlon,nlat,zeex,nlbper)
      call bextr(clwfl, nlon,nlat,zeex,nlbper)
      call bextr(chflux,nlon,nlat,zeex,nlbper)
      call bextr(cqflux,nlon,nlat,zeex,nlbper)
      call bextr(cwvflux,nlon,nlat,zeex,nlbper)
      call bextr(csoilh_bott_flux,nlon,nlat,zeex,nlbper)
      call bextr(csoilw_bott_flux,nlon,nlat,zeex,nlbper)
      call bextr(iwv   ,nlon,nlat,zeex,nlbper)

      do k=1,nlevg
        call bextr(tg_post(:,:,k), nlon,nlat,zeex,nlbper)
      enddo

!  Resetting of unphysical values (due only to extrap. at bound.)

      do  jlat = 1, nlat
      do  jlon = 1, nlon
      totpre (jlon,jlat) = max(totpre(jlon,jlat), 0.)
      conpre (jlon,jlat) = max(conpre(jlon,jlat), 0.)
      snfall (jlon,jlat) = max(snfall(jlon,jlat), 0.)
      snow   (jlon,jlat) = max(snow  (jlon,jlat), 0.)
      runoff (jlon,jlat) = max(runoff(jlon,jlat), 0.)
      runoff_tot (jlon,jlat) = max(runoff_tot(jlon,jlat), 0.)
      cloud  (jlon,jlat) = max(cloud (jlon,jlat), 0.)
      cloud  (jlon,jlat) = min(cloud (jlon,jlat), 1.)
      cloudt (jlon,jlat) = max(cloudt(jlon,jlat), 0.)
      cloudt (jlon,jlat) = min(cloudt(jlon,jlat), 1.)
      cloudh (jlon,jlat) = max(cloudh(jlon,jlat), 0.)
      cloudh (jlon,jlat) = min(cloudh(jlon,jlat), 1.)
      cloudm (jlon,jlat) = max(cloudm(jlon,jlat), 0.)
      cloudm (jlon,jlat) = min(cloudm(jlon,jlat), 1.)
      cloudl (jlon,jlat) = max(cloudl(jlon,jlat), 0.)
      cloudl (jlon,jlat) = min(cloudl(jlon,jlat), 1.)
      q2     (jlon,jlat) = max(q2    (jlon,jlat), 0.)
      q2rel  (jlon,jlat) = max(q2rel (jlon,jlat), 0.)
      q2rel  (jlon,jlat) = min(q2rel (jlon,jlat), 102.)
      zqgrs  (jlon,jlat) = max(zqgrs (jlon,jlat), 0.)
      q05rel (jlon,jlat) = max(q05rel(jlon,jlat), 0.)
      q05rel (jlon,jlat) = min(q05rel(jlon,jlat), 102.)
      iwv    (jlon,jlat) = max(iwv   (jlon,jlat), 0.)
      enddo
      enddo

!  Rescaling of mslp (in hPa) and def. of geop. heights

      do jlat = 1, nlat
      do jlon = 1, nlon
      zpz0(jlon,jlat) = zpz0(jlon,jlat)/100.
      enddo
      enddo

!  Horizontal smoothing of fields when extrapolated below ground;
!  filtering of msl pressure as a function of orographic height
!  (zfp is a filtering factor, empirically depending on zfpg)

      zfpg=1./(g*3500.)
      do iter = 1,iterf
      do jlat = 2,nlat-1
      do jlon = 2,nlon-1
      zfp = phig(jlon,jlat)*zfpg
      zfp=min(zfp,0.8)
      zfp=max(zfp,0.)
      zwork(jlon,jlat) = (1.-zfp) * zpz0(jlon,jlat) + zfp/4.* (zpz0(jlon-1,jlat)+zpz0(jlon+1,jlat)+  &
                          zpz0(jlon,jlat-1)+zpz0(jlon,jlat+1))
      enddo
      enddo

!  Filtering of mslp at lateral boundaries

      do jlat = 2,nlatm1
      jlon=1
      zfp = phig(jlon,jlat)*zfpg
      zfp=min(zfp,0.8)
      zfp=max(zfp,0.)
      zwork(jlon,jlat) = (1.-zfp) * zpz0(jlon,jlat) + zfp/2.* (zpz0(jlon,jlat-1)+zpz0(jlon,jlat+1))
      jlon=nlon
      zfp = phig(jlon,jlat)*zfpg
      zfp=min(zfp,0.8)
      zfp=max(zfp,0.)
      zwork(jlon,jlat) = (1.-zfp) * zpz0(jlon,jlat) + zfp/2.* (zpz0(jlon,jlat-1)+zpz0(jlon,jlat+1))
      enddo

      do jlon = 2,nlon-1
      jlat=1
      zfp = phig(jlon,jlat)*zfpg
      zfp=min(zfp,0.8)
      zfp=max(zfp,0.)
      zwork(jlon,jlat) = (1.-zfp) * zpz0(jlon,jlat) + zfp/2.* (zpz0(jlon-1,jlat)+zpz0(jlon+1,jlat))
      jlat=nlat
      zfp = phig(jlon,jlat)*zfpg
      zfp=min(zfp,0.8)
      zfp=max(zfp,0.)
      zwork(jlon,jlat) = (1.-zfp) * zpz0(jlon,jlat) + zfp/2.* (zpz0(jlon-1,jlat)+zpz0(jlon+1,jlat))
      enddo

      zwork(1,1) = zpz0(1,1)
      zwork(1,nlat) = zpz0(1,nlat)
      zwork(nlon,1) = zpz0(nlon,1)
      zwork(nlon,nlat) = zpz0(nlon,nlat)

      do jlat = 1,nlat
      do jlon = 1,nlon
      zpz0(jlon,jlat) = zwork(jlon,jlat)
      enddo
      enddo
      enddo

!  Filtering of temp. and geop. at pressure levels where below ground
!  Lateral boundary points are not filtered, so extrap. is applied later;
!  only levels below 400e2 Pa are considered

      nlevmx = 1
      do jlev = 1,nlevpo
      if(plevo(jlev).gt.400.e2) nlevmx = jlev
      enddo

! Values used to interpolate below ground are those from the second model level above ground
! (to assure smoothing)

      do jlev = 1, nlevmx
      do iter = 1,iterf
      do jlat = 2,nlat-1
      do jlon = 2,nlon-1
      zzpp=pzer*sigint(nlev-1)-(pzer-ps(jlon,jlat))*sigalf(nlev-1)
      if(plevo(jlev).gt.zzpp) then
      zwork(jlon,jlat) = .3*zt(jlon,jlat,jlev) + .7/4.*(zt(jlon-1,jlat,jlev)+zt(jlon+1,jlat,jlev)+  &
            zt(jlon,jlat-1,jlev)+zt(jlon,jlat+1,jlev))
      zwor1(jlon,jlat) = .3*zthetae(jlon,jlat,jlev) + .7/4.* (zthetae(jlon-1,jlat,jlev)+       &
            zthetae(jlon+1,jlat,jlev)+zthetae(jlon,jlat-1,jlev)+zthetae(jlon,jlat+1,jlev))
      zwor2(jlon,jlat) = .3*zph(jlon,jlat,jlev) + .7/4.* (zph(jlon-1,jlat,jlev)+zph(jlon+1,jlat,jlev)+ &
            zph(jlon,jlat-1,jlev)+zph(jlon,jlat+1,jlev))
      else
      zwork(jlon,jlat) = zt(jlon,jlat,jlev)
      zwor1(jlon,jlat) = zthetae(jlon,jlat,jlev)
      zwor2(jlon,jlat) = zph(jlon,jlat,jlev)
      endif
      enddo
      enddo

      do jlat = 2,nlatm1
      do jlon = 2,nlonm1
      zt (jlon,jlat,jlev) = zwork(jlon,jlat)
      zthetae(jlon,jlat,jlev) = zwor1(jlon,jlat)
      zph(jlon,jlat,jlev) = zwor2(jlon,jlat)
      enddo
      enddo
      enddo ! end of loop in ITERF

!  Extrap. to lateral boundaries where filtering was not applied

      zeex = .7
      call bextr(zt (1,1,jlev),nlon,nlat,zeex,nlbper)
      call bextr(zph(1,1,jlev),nlon,nlat,zeex,nlbper)

      enddo ! loop in JLEV

!  End of horizontal smoothing

      do jlat = 1, nlat
      do jlon = 1, nlon

!  First level above topography

      jlev1=nlevpo
      do jlev = 1, nlevpo
      if(plevo(jlev).lt.zsp(jlon,jlat)) jlev1=min(jlev1,jlev)
      enddo

      do jlev = 1, nlevpo

!  Relative humidity on pressure levels: values are limited to saturation
!  (slight over-saturation only for plot smoothing)
!  (zrh is not used in cross-sections, where supersaturation is allowed)
!  RH is computed with respect to water and ice separately

      zf = zq(jlon,jlat,jlev)
      ztdan = zt(jlon,jlat,jlev)
      call comp_esk(zzpvs,zzqs,ztdan,plevo(jlev),1)! computes saturation to water and ice separately
      zzpv=(zf*plevo(jlev))/(.622 + zf -.622*zf)
      if(jlev.eq.jlev1) zrelh1=min(zzpv/zzpvs, 1.) ! used to extrapolate below ground
      zrh(jlon,jlat,jlev) = min(zzpv/zzpvs*100., 102.)
      enddo

! Extrapolation of Q below ground assuming constant RH

      do jlev=1,jlev1-1
      ztdan = zt(jlon,jlat,jlev)
      call comp_esk(zzpvs,zzqs,ztdan,plevo(jlev),1)! computes saturation to water and ice separately
      zfs =.622*zzpvs/(plevo(jlev)-zzpvs)
      zfs = zfs/(1.+zfs)
      zrelh1=min(0.5*(zrelh1+0.01*q2rel(jlon,jlat)),0.99) ! media tra rh al primo liv. e rh a 2m
      zq(jlon,jlat,jlev) =  zrelh1*zfs
      zrh(jlon,jlat,jlev) = zrelh1*100.
      enddo

      enddo
      enddo

!  Interpolation of various fields defined at u, v and z-points onto t-points over the staggered grid
!  Interp. of u and v at western bound. lin. extrap. for dlon/2

      do jlat = 1, nlat
      zaxus(jlat) = (3.*zus(1,jlat) - zus(2,jlat))/2.
      do jlev = 1, nlevpo
      zaxu(jlat,jlev) = (3.*zu(1,jlat,jlev) - zu(2,jlat,jlev))/2.
      enddo
      enddo
      do jlat = 1, nlat
      do jlon = nlon, 2, -1
      zus(jlon,jlat) = (zus(jlon,jlat) + zus(jlon-1,jlat))/2.
      do jlev = 1, nlevpo
      zu(jlon,jlat,jlev) = (zu(jlon,jlat,jlev) + zu(jlon-1,jlat,jlev))/2.
      enddo
      enddo
      enddo
      do jlat = 1, nlat
      zus(1,jlat) = zaxus(jlat)
      do jlev = 1, nlevpo
      zu(1,jlat,jlev) = zaxu(jlat,jlev)
      enddo
      enddo

!  At northern bound. lin. extrap. for dlat/2

      do jlon = 1, nlon
      zaxvs(jlon)   = (3.* zvs(jlon,nlat) -  zvs(jlon,nlat-1))/2.
      do jlev = 1, nlevpo
      zaxv(jlon,jlev) = (3.*zv(jlon,nlat,jlev) - zv(jlon,nlat-1,jlev))/2.
      enddo
      enddo

      do jlat = 1, nlatm1
      do jlon = 1, nlon
      zvs(jlon,jlat)  = ( zvs(jlon,jlat) +  zvs(jlon,jlat+1))/2.
      do jlev = 1, nlevpo
      zv(jlon,jlat,jlev) = (zv(jlon,jlat,jlev) + zv(jlon,jlat+1,jlev))/2.
      enddo
      enddo
      enddo
      do jlon = 1, nlon
      zvs(jlon,nlat)  =zaxvs(jlon)
      do jlev = 1, nlevpo
      zv(jlon,nlat,jlev) =zaxv(jlon,jlev)
      enddo
      enddo

! Integrated Vapour Transport (put in zus,zvs)

!!      do jlat = 1, nlat
!!      jlatp1 = min(jlat+1,nlat)
!!      do jlon = 1, nlon
!!      jlonm1 = max(jlon-1,1)
!!       zus (jlon,jlat) = 0.
!!       zvs (jlon,jlat) = 0.
!!       do jk = 1,nlev
!!        zzpp = pzer*sigint(jk)-(pzer-ps(jlon,jlat))*sigalf(jk)
!!        if (zzpp.gt.30000.) then
!!        zzdp = pzer*dsig(jk)-(pzer-ps(jlon,jlat))*dsigalf(jk)
!!        zus(jlon,jlat) = zus(jlon,jlat)+q(jlon,jlat,jk)*.5*(u(jlon,jlat,jk)+u(jlonm1,jlat,jk))/g*zzdp
!!        zvs(jlon,jlat) = zvs(jlon,jlat)+q(jlon,jlat,jk)*.5*(v(jlon,jlat,jk)+v(jlon,jlatp1,jk))/g*zzdp
!!        endif
!!       enddo
!!      enddo
!!      enddo

!  Interp. of pot. vort. and rel. vorticity, comp. at z points using the four surrounding z points
!  (unlike for velocity comp., there is no extrapol. on half grid)

      do jlat = 1, nlatm1
      do jlon = 2, nlon
      do jlev = 1, nlevpo
      zpva(jlon,jlat,jlev) = ((zpv(jlon,jlat,jlev) + zpv(jlon-1,jlat,jlev))+ zpv(jlon-1,jlat+1,jlev) + &
                               zpv(jlon,jlat+1,jlev))/4.
      zrva(jlon,jlat,jlev) = ((zrv(jlon,jlat,jlev) + zrv(jlon-1,jlat,jlev))+ zrv(jlon-1,jlat+1,jlev) + &
                               zrv(jlon,jlat+1,jlev))/4.
      enddo

      do jlev = 1, nlevto    ! for theta surf.
      zpvat(jlon,jlat,jlev) = ((zpvt(jlon,jlat,jlev) + zpvt(jlon-1,jlat,jlev))+ zpvt(jlon-1,jlat+1,jlev) &
                               + zpvt(jlon,jlat+1,jlev))/4.
      enddo
      enddo
      enddo

!  At western bound. lin. extrap. using already defined interior values

      do jlat = 1, nlat-1
        do jlev = 1, nlevpo
         zpva(1,jlat,jlev) = 2.*zpva(2,jlat,jlev) - zpva(3,jlat,jlev)
         zrva(1,jlat,jlev) = 2.*zrva(2,jlat,jlev) - zrva(3,jlat,jlev)
        enddo
        do jlev = 1, nlevto ! per sup. theta
         zpvat(1,jlat,jlev) = 2.*zpvat(2,jlat,jlev)-zpvat(3,jlat,jlev)
        enddo
      enddo

!  At northern bound. lin. extrap. using already defined interior values

      do jlon = 1, nlon
         do jlev = 1, nlevpo
         zpva(jlon,nlat,jlev) = 2.*zpva(jlon,nlat-1,jlev) - zpva(jlon,nlat-2,jlev)
         zrva(jlon,nlat,jlev) = 2.*zrva(jlon,nlat-1,jlev) - zrva(jlon,nlat-2,jlev)
         enddo

         do jlev = 1, nlevto ! per sup. theta
         zpvat(jlon,nlat,jlev) = 2.*zpvat(jlon,nlat-1,jlev) - zpvat(jlon,nlat-2,jlev)
         enddo
       enddo

!  Final conversions and resettings
!  Rescaling of pot. vort. to get it in "pv units"
!  Rescaling of rel. vort. with 1.e-4 (mid-lat. Coriolis)
!  Conversion from geopotential to geopotential height
!  Conversion from Kelvin to Celsius

      do jlon = 1, nlon
      do jlat = 1, nlat
      t2(jlon,jlat)   = t2(jlon,jlat)-ttr
      td2(jlon,jlat)  = td2(jlon,jlat)-ttr
      do jlev = 1, nlevpo
      zpva(jlon,jlat,jlev) = zpva(jlon,jlat,jlev)*1.e6
      zrva(jlon,jlat,jlev) = zrva(jlon,jlat,jlev)*1.e4
      zph(jlon,jlat,jlev)  = zph(jlon,jlat,jlev)/9.8
      zt(jlon,jlat,jlev)   = zt(jlon,jlat,jlev)-ttr
      enddo

      do jlev = 1, nlevto ! for theta surfaces
      zpvat(jlon,jlat,jlev) = zpvat(jlon,jlat,jlev)*1.e6
      enddo

! Resetting of underground values of u, v, pv after the horiz. interpol on t points
! (note that for u, v, pv, resetting to 0 had already been made on non-t points after vert. interp. -
! this has the effect of making wind speed smaller at some points near intersection with the ground

      do jlev = 1, nlevmx
         if(plevo(jlev).gt.zsp(jlon,jlat)) then
         zu(jlon,jlat,jlev) = 0.
         zv(jlon,jlat,jlev) = 0.
         endif
      enddo

! Vorticity and pv zeroed if under or very close to the ground (20 hPa)

       do jlev = 1, nlevmx
         if(plevo(jlev).gt.(zsp(jlon,jlat)-2000.)) then
         zpva(jlon,jlat,jlev) = 0.
         zrva(jlon,jlat,jlev) = 0.
         endif
      enddo

! Downward short-wave radiation flux (W m-2)

      zfsnow=min( ((snow(jlon,jlat)/5.e-3)**0.5), 1.0)
      zalbedo=albedo(jlon,jlat)*(1.-zfsnow)+snow_albedo(jlon,jlat)*zfsnow
      cswdfl(jlon,jlat)=cswfl(jlon,jlat)/(1.-zalbedo)

! Cumulated fluxes in kjoule/m2

      cswfl(jlon,jlat) =cswfl(jlon,jlat)*1.e-3
      clwfl(jlon,jlat) =clwfl(jlon,jlat)*1.e-3
      chflux(jlon,jlat)=chflux(jlon,jlat)*1.e-3
      cqflux(jlon,jlat)=cqflux(jlon,jlat)*1.e-3

      enddo
      enddo

! Calculation and writing of forecast data

  idate0(1:5)  = nfdr(5:9)
  iperiod(1:3) = nfdr(10:12)

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
! call outgraph(80102,1,nlon,nlat,pdr(39),pdr(38),pdr(5),pdr(4)+pdr(1)*0.5,pdr(2),pdr(1), &
! 2, 0,  0,  1, 0, 0,idate0(1:5),iperiod(1:3),fmask,-1.,1.)
! call outgraph(80102,1,nlon,nlat,pdr(39),pdr(38),pdr(5),pdr(4)+pdr(1)*0.5,pdr(2),pdr(1), &
! 0, 3,  6,  7, 0, 0,idate0(1:5),iperiod(1:3),h_tropopause,1.,0.)
! call outgraph(80102,1,nlon,nlat,pdr(39),pdr(38),pdr(5),pdr(4)+pdr(1)*0.5,pdr(2),pdr(1), &
! 0, 3,  0,  7, 0, 0,idate0(1:5),iperiod(1:3),p_tropopause,1.e-2,0.)
! call outgraph(80102,1,nlon,nlat,pdr(39),pdr(38),pdr(5),pdr(4)+pdr(1)*0.5,pdr(2),pdr(1), &
! 0, 0,  0,  7, 0, 0,idate0(1:5),iperiod(1:3),t_tropopause,1.,-273.15)
! call outgraph(80102,1,nlon,nlat,pdr(39),pdr(38),pdr(5),pdr(4)+pdr(1)*0.5,pdr(2),pdr(1), &
! 0, 2,  2,  7, 0, 0,idate0(1:5),iperiod(1:3),u_tropopause,1.,0.)
! call outgraph(80102,1,nlon,nlat,pdr(39),pdr(38),pdr(5),pdr(4)+pdr(1)*0.5,pdr(2),pdr(1), &
! 0, 2,  3,  7, 0, 0,idate0(1:5),iperiod(1:3),v_tropopause,1.,0.)

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

! Lapse_rate_clim is a climatological prescribed lapse rate, changing
! between 4.5e-3 (end of dec.) and 7.4e-3 (end of june)

   lapse_rate_clim = 4.5e-3 + 2.8e-3 * coeday

! computation of average lapse rate of input model between lower levels

   lapse_rate=0.
   nnn=0
   do jlat=1,nlat
   do jlon=1,nlon
     if (fmask(jlon,jlat) < 0.5) then
      lapse_rate=lapse_rate+(temp(jlon,jlat,nlev-4)-temp(jlon,jlat,nlev))*g/(phi(jlon,jlat,nlev-4)-phi(jlon,jlat,nlev))
      nnn=nnn+1
     endif
   enddo
   enddo
   lapse_rate=lapse_rate/float(nnn)
   lapse_rate=min (0.70*lapse_rate - 0.30*lapse_rate_clim, 0.) ! limited to avoid inversions

   lapse_rate_1(:,:)=lapse_rate

   do jlat=1,nlat
   do jlon=1,nlon
     jklev1=nlev
     jklev2=nlev-1
     do jklev=nlev,1,-1
       if ((phi(jlon,jlat,jklev)-phig(jlon,jlat))/g >= 200.) exit
     enddo
     jklev1=jklev
     do jklev=1,nlev-1
       if ((phi(jlon,jlat,jklev)-phig(jlon,jlat))/g <= 500.) exit
     enddo
     jklev2=min(jklev, jklev1-1)
     lapse_rate=(temp(jlon,jlat,jklev2)-temp(jlon,jlat,jklev1))*g/(phi(jlon,jlat,jklev2)-phi(jlon,jlat,jklev1))
     lapse_rate_2(jlon,jlat)=lapse_rate
   enddo
   enddo

!=======================================================================

 do jklev = 1, nlevpo
   isobarlev = int(plevo(jklev)/100.)
   if (isobarlev==850) lev850=jklev
   if (isobarlev==700) lev700=jklev
 enddo

 x0    = pdr(39)
 y0    = pdr(38)
 alat0 = pdr(4)
 alon0 = pdr(5)
 if (alon0 > 180.) alon0=alon0-360.
 dlat  = pdr(1)
 dlon  = pdr(2)

! ---- For RTTOV simulation ---

#ifdef rttov

 fseaice=fice

 do jklev=1,nlev
 do jlat=1,nlat
 do jlon=1,nlon
   p(jlon,jlat,jklev)=pzer*sigint(jklev)-(pzer-ps(jlon,jlat))*sigalf(jklev)
   if (qc(jlon,jlat,jklev) > 1.e-8) then
     call qsat_entropy(temp(jlon,jlat,jklev), p(jlon,jlat,jklev), z1, z2, z3, z4, z5, z6, fracw)
     qcw(jlon,jlat,jklev)=qc(jlon,jlat,jklev)*fracw
     qci(jlon,jlat,jklev)=qc(jlon,jlat,jklev)*(1.-fracw)
   else
     qcw(jlon,jlat,jklev)=0.
     qci(jlon,jlat,jklev)=0.
   endif
 enddo
 enddo
 enddo

 do jlat=1,nlat
 do jlon=1,nlon
   if (totpre(jlon,jlat) > 0.) then
     if (conpre(jlon,jlat)/totpre(jlon,jlat) > 0.5) then
       flag_cloud_conv(jlon,jlat) = 1
     else
       flag_cloud_conv(jlon,jlat) = 0
     endif
   else
     flag_cloud_conv(jlon,jlat) = 0
     do jklev=1,nlev-5
       zrho = p(jlon,jlat,jklev)/rd/tvirt(jlon,jlat,jklev)
       if (omeg(jlon,jlat,jklev)/(zrho*g) > w_crit_conv) flag_cloud_conv(jlon,jlat) = 1
     enddo
   endif
 enddo
 enddo

 date_time(1:5) = idatec(1:5)
 iini=max(nint(nlon*zoom_xini), 1)
 ifin=min(nint(nlon*zoom_xfin), nlon)
 jini=max(nint(nlat*zoom_yini), 1)
 jfin=min(nint(nlat*zoom_yfin), nlat)
 npoint_jump = max( (ifin-iini+1)*(jfin-jini+1)/300000, 1)
 alon1=pdr(5)+pdr(2)*(iini-1)
 alat1=pdr(4)+pdr(1)*0.5+pdr(1)*(jini-1)

 print *
 print *,"------------ Call RTTOV_fdw model Version 11 ----------------"

 call rttov11_fwd_interface &
 (nchan, rttov_sat_series, rttov_sat_id, rttov_sat_sensor, channel_elabor_list(1:nchan), flag_visible, &
 path_rttov_emis_atlas, path_rttov_brdf_atlas,                                                         &
 ifin-iini+1, jfin-jini+1, npoint_jump, nlev,                                                          &
 alatt(iini:ifin,jini:jfin), alont(iini:ifin,jini:jfin), fmask(iini:ifin,jini:jfin),                   &
 phig(iini:ifin,jini:jfin)/g, emismap1(iini:ifin,jini:jfin),                                           &
 ps(iini:ifin,jini:jfin), t2(iini:ifin,jini:jfin)+ttr, q2(iini:ifin,jini:jfin),                        &
 u10(iini:ifin,jini:jfin), v10(iini:ifin,jini:jfin), tskin(iini:ifin,jini:jfin),                       &
 fsnow(iini:ifin,jini:jfin), qg(iini:ifin,jini:jfin,1),                                                &
 fseaice(iini:ifin,jini:jfin), flag_cloud_conv(iini:ifin,jini:jfin),                                   &
 p(iini:ifin,jini:jfin,:), temp(iini:ifin,jini:jfin,:), q(iini:ifin,jini:jfin,:),                      &
 qcw(iini:ifin,jini:jfin,:), qci(iini:ifin,jini:jfin,:), clfrac(iini:ifin,jini:jfin,:), date_time,     &
 sensor_chan_id, sensor_chan_cw, &
 radiance(iini:ifin,jini:jfin,:), radiance_bt(iini:ifin,jini:jfin,:), radiance_refl(iini:ifin,jini:jfin,:), &
 radiance_clear(iini:ifin,jini:jfin,:), radiance_clear_bt(iini:ifin,jini:jfin,:),                           &
 radiance_clear_refl(iini:ifin,jini:jfin,:), emis_rttov(iini:ifin,jini:jfin,:))

 print *
 print *,"------------ RTTOV-11 forward model simulated data for ------------"
 print *,"---- Satellite code ",rttov_sat_series," number ",rttov_sat_id," instrument (sensor) code ",rttov_sat_sensor," ----"
 print *,"---- Elaborated Channel list ----"
 print *,sensor_chan_id(:)
 print *,"---- Elaborated Channel Central WaveNumber ----"
 print *,nint(sensor_chan_cw(:))
 write (*,"(4(a,i4),a,i1,a)") &
 "---- fiels defined at subdomain (",iini,":",ifin,",",jini,":",jfin,") with point jump ",npoint_jump," ----"

! call outgraph(80102,1,ifin-iini+1,jfin-jini+1,pdr(39),pdr(38),alon1,alat1,pdr(2),pdr(1), &
! 2, 0,  0,  1, 1, 0,idate0(1:5),iperiod(1:3),fmask(iini:ifin,jini:jfin),-1.,1.)
! call outgraph(80102,1,ifin-iini+1,jfin-jini+1,pdr(39),pdr(38),alon1,alat1,pdr(2),pdr(1), &
! 3, 1, 16,  1, 1, 0,idate0(1:5),iperiod(1:3),radiance(iini:ifin,jini:jfin,1),1.,0.)
! call outgraph(80102,1,ifin-iini+1,jfin-jini+1,pdr(39),pdr(38),alon1,alat1,pdr(2),pdr(1), &
! 3, 1, 14,  1, 1, 0,idate0(1:5),iperiod(1:3),radiance_clear(iini:ifin,jini:jfin,1),1.,0.)
! call outgraph(80102,1,ifin-iini+1,jfin-jini+1,pdr(39),pdr(38),alon1,alat1,pdr(2),pdr(1), &
! 3, 1, 17,  1, 1, 0,idate0(1:5),iperiod(1:3),radiance_bt(iini:ifin,jini:jfin,1),1.,-273.15)
! call outgraph(80102,1,ifin-iini+1,jfin-jini+1,pdr(39),pdr(38),alon1,alat1,pdr(2),pdr(1), &
! 3, 1, 15,  1, 1, 0,idate0(1:5),iperiod(1:3),radiance_clear_bt(iini:ifin,jini:jfin,1),1.,-273.15)
! call outgraph(80102,1,ifin-iini+1,jfin-jini+1,pdr(39),pdr(38),alon1,alat1,pdr(2),pdr(1), &
! 3, 1,201,  1, 1, 0,idate0(1:5),iperiod(1:3),emis_rttov(iini:ifin,jini:jfin,1),1.,0.)

 call write_grib2_horizontal_grid_rttov_data (1, ifin-iini+1, jfin-jini+1, npoint_jump, &
 pdr(39), pdr(38), alon1, alat1, pdr(2), pdr(1),                                        &
 idate0, iperiod, iperiod_accum, idatec,                                                &
 nchan, rttov_sat_series, rttov_sat_id, rttov_sat_sensor,                               &
 sensor_chan_id, sensor_chan_cw,                                                        &
 radiance(iini:ifin,jini:jfin,:), radiance_bt(iini:ifin,jini:jfin,:),                   &
 radiance_clear(iini:ifin,jini:jfin,:), radiance_clear_bt(iini:ifin,jini:jfin,:), emis_rttov(iini:ifin,jini:jfin,:) )

 print *
 print *,"------------------ RTTOV-11 products created ---------------------"
 print *

#endif

! -----------------------------

!******************* Output definitions and write **********************

      if (output_format_ppf) then

! Land-sea mask, longitude and latitude or rotated grid

      if(kist.eq.1) then
      open(80,file='bolam.ppf',status='unknown',form='formatted')
      write(80,1111) (isflag(i),i=1,80)

!      write(80,1112) ipflag
      ipflag20(:,:) = ipflag(:, 1:20) ! provisional, for compatibility with NCARG plotting of ppf file
      write(80,1112) ipflag20         ! provisional, for compatibility with NCARG plotting of ppf file

      write(80,1112) (itflag(i),i=1,nlevto)
 1111 format(80i1)
 1112 format(12i1)

      call wrpost(fmask)
      call wrpost(alont)
      call wrpost(alatt)
      endif

! Orography (in m)

      if(isflag(1).eq.1) call wrpost(zorogr)

! PZ0 (mslp, in hPa)

      if(isflag(2).eq.1) call wrpost(zpz0)

! TOTPRE (tot. accum. precip. in kg/m2 eq. to mm)

      if(isflag(3).eq.1) call wrpost(totpre)

! CONPRE (convective accum. precip.)

      if(isflag(4).eq.1) call wrpost(conpre)

! Accumulated snow fall (mm of equivalent water)

      if(isflag(5).eq.1) call wrpost(snfall)

! U10 (wind comp. at 10 m)

      if(isflag(6).eq.1) call wrpost(u10)

! V10 (wind comp. at 10 m)

      if(isflag(6).eq.1) call wrpost(v10)

! T2 (air temperature at 2 m)

      if(isflag(7).eq.1) call wrpost(t2)

! Q2 (specific humidity at 2 m)

      if(isflag(8).eq.1) call wrpost(q2)

! Q2REL (relative humidity at 2 m)

      if(isflag(9).eq.1) call wrpost(q2rel)

! ZQGRS (skin specific humidity)

      if(isflag(10).eq.1) call wrpost(zqgrs)

! ZTGRS (skin temperature)

      if(isflag(11).eq.1) call wrpost(ztgrs)

! HFL (upward flux of sens. heat in watt/m**2)

      if(isflag(12).eq.1) call wrpost(hflux)

! QFL (upward flux of latent heat in watt/m**2)

      if(isflag(13).eq.1) call wrpost(qflux)

! US (wind vel. at nlev: lowest sigma)

      if(isflag(14).eq.1) call wrpost(zus)

! VS (wind vel. at nlev: lowest sigma)

      if(isflag(14).eq.1) call wrpost(zvs)

! Temperature at nlev: lowest sigma

      if(isflag(15).eq.1) call wrpost(zts)

! ZQS (spec. humid. at nlev: lowest sigma)

      if(isflag(16).eq.1) call wrpost(zqs)

! Total cloud cover (comp. post mortem)

      if(isflag(17).eq.1) call wrpost(cloud)

! Cloud cover (3 layers - comp. post mortem)

      if(isflag(18).eq.1) call wrpost(cloudh)
      if(isflag(19).eq.1) call wrpost(cloudm)
      if(isflag(20).eq.1) call wrpost(cloudl)

! Total cloud cover (from BOLAM MHF)

      if(isflag(21).eq.1) call wrpost(cloudt)

! Snow height (m of equivalent water)

      if(isflag(22).eq.1) call wrpost(snow)

! Accumulated run-off (kg/m2)

      if(isflag(23).eq.1) call wrpost(runoff)

! CAPE (j/kg)

      if(isflag(24).eq.1) call wrpost(zcape)

! ZQGR: water of ground layers

      if(isflag(25).eq.1) call wrpost(zqgr1)
      if(isflag(26).eq.1) call wrpost(zqgr2)
      if(isflag(27).eq.1) call wrpost(zqgr3)
      if(isflag(28).eq.1) call wrpost(zqgr4)

! ZTGR: temp. of ground layers

      if(isflag(29).eq.1) call wrpost(ztgr1)
      if(isflag(30).eq.1) call wrpost(ztgr2)
      if(isflag(31).eq.1) call wrpost(ztgr3)
      if(isflag(32).eq.1) call wrpost(ztgr4)

! Ground temperature at observation levels (5, 10, 20, 50 and 100 cm)

      if(isflag(33).eq.1) call wrpost(tg005)
      if(isflag(34).eq.1) call wrpost(tg010)
      if(isflag(35).eq.1) call wrpost(tg020)
      if(isflag(36).eq.1) call wrpost(tg050)
      if(isflag(37).eq.1) call wrpost(tg100)

! Min. and max. temperature at 2 m

      if(isflag(38).eq.1) call wrpost(ztmin)
      if(isflag(39).eq.1) call wrpost(ztmax)

! Temperature and relative humidity at 0.5 m

      if(isflag(40).eq.1) call wrpost(t05)
      if(isflag(41).eq.1) call wrpost(q05rel)

! Temperature at 0.05 m

      if(isflag(42).eq.1) call wrpost(t005)

! Cumulated fluxes: short wave radiation, long wave radiation,
! sensible heat, latent heat (in kjoule/m2)

      if(isflag(43).eq.1) call wrpost(cswfl)
      if(isflag(44).eq.1) call wrpost(clwfl)
      if(isflag(45).eq.1) call wrpost(chflux)
      if(isflag(46).eq.1) call wrpost(cqflux)

! Radiative parameters of land surface (without considering snow at ground)
! Emissivity: emismap1 broadband; emismap2 in 8-12 micron window

      if(isflag(47).eq.1) call wrpost(albedo)
      if(isflag(48).eq.1) call wrpost(emismap1)
      if(isflag(49).eq.1) call wrpost(emismap2)

! Lifted index (K)

      if(isflag(50).eq.1) call wrpost(zlift)

! Level of 0°c (m)

      if(isflag(51).eq.1) call wrpost(zlev0)

! Sea ice fraction (0-1)

      if(isflag(52).eq.1) call wrpost(fice)

! Sea ice thickness (m)

      if(isflag(53).eq.1) call wrpost(iceth)

! TD2 - dew point temperature at 2 m

      if(isflag(54).eq.1) call wrpost(td2)

! IWV - integrated water vapour (kg/m2)

      if(isflag(55).eq.1) call wrpost(iwv)

!-----------------------------------------------------------------------
!--------------- Loop for fields defined on pressure levels  -----------
!-----------------------------------------------------------------------

      do jlev = 1, nlevpo
      isobarlev = nint(plevo(jlev)/100.)

! PHI (here is geopotential height)

      if(ipflag(1, jlev).eq.1) call wrpost(zph(1,1,jlev))

! T (temperature, Celsius)

      if(ipflag(2, jlev).eq.1) call wrpost(zt(1,1,jlev))

! U (wind comp., m/s)

      if(ipflag(3, jlev).eq.1) call wrpost(zu(1,1,jlev))

! V (wind comp., m/s)

      if(ipflag(3, jlev).eq.1) call wrpost(zv(1,1,jlev))

! Q (specific humidity)

      if(ipflag(4, jlev).eq.1) call wrpost(zq(1,1,jlev))

! RH (relative humidity, % from 0 to 100)

      if(ipflag(5, jlev).eq.1) call wrpost(zrh(1,1,jlev))

! Omega (vert. vel., Pa/s)

      if(ipflag(6, jlev).eq.1) call wrpost(zom(1,1,jlev))

! Relative vorticity

      if(ipflag(7, jlev).eq.1) call wrpost(zrva(1,1,jlev))

! Ertel potential vorticity

      if(ipflag(8, jlev).eq.1) call wrpost(zpva(1,1,jlev))

! Equiv. potent. temp. thetae

      if(ipflag(9, jlev).eq.1) call wrpost(zthetae(1,1,jlev))

      enddo

! Pot. vort. on selected isentropic surfaces

      do jlev = 1, nlevto ! per sup. theta
       if(itflag(jlev).eq.1) call wrpost(zpvat(1,1,jlev))
      enddo

!=======================================================================
! Cross-section variables are written on a separate file
! Note that the order is by each cross-section line
! The order of variables depends on superpositions made in graphics

      if (ncrossec) then

      if(kist.eq.1) then
      open(81,file='bolam_crossy.ppf',status='unknown',form='formatted')
      open(82,file='bolam_crossx.ppf',status='unknown',form='formatted')
      endif

      do jnc = 1, ncrossy

! Potential temperature

      call wrposy(ztcr(1,1,jnc) ,zspcr(1,jnc),zplcr,lncr(jnc),npress)

! Zonal velocity component

      call wrposy(zucr(1,1,jnc) ,zspcr(1,jnc),zplcr,lncr(jnc),npress)

! Equivalent potential temperature

      call wrposy(ztecr(1,1,jnc),zspcr(1,jnc),zplcr,lncr(jnc),npress)

! Merid. velocity component

      call wrposy(zvcr(1,1,jnc) ,zspcr(1,jnc),zplcr,lncr(jnc),npress)

! Vertical. velocity component w, computed approx. from omega

      call wrposy(zwcr(1,1,jnc) ,zspcr(1,jnc),zplcr,lncr(jnc),npress)

! Potential vorticity

      call wrposy(zpvcr(1,1,jnc),zspcr(1,jnc),zplcr,lncr(jnc),npress)

! Relative humidity

      call wrposy(zrhcr(1,1,jnc),zspcr(1,jnc),zplcr,lncr(jnc),npress)
      enddo

! Same for west-east cross sections

      do jnc = 1, ncrossx
      call wrposx(ztcrx(1,1,jnc) ,zspcrx(1,jnc),zplcr,ltcr(jnc),npress)
      call wrposx(zvcrx(1,1,jnc) ,zspcrx(1,jnc),zplcr,ltcr(jnc),npress)
      call wrposx(ztecrx(1,1,jnc),zspcrx(1,jnc),zplcr,ltcr(jnc),npress)
      call wrposx(zucrx(1,1,jnc) ,zspcrx(1,jnc),zplcr,ltcr(jnc),npress)
      call wrposx(zwcrx(1,1,jnc), zspcrx(1,jnc),zplcr,ltcr(jnc),npress)
      call wrposx(zpvcrx(1,1,jnc),zspcrx(1,jnc),zplcr,ltcr(jnc),npress)
      call wrposx(zrhcrx(1,1,jnc),zspcrx(1,jnc),zplcr,ltcr(jnc),npress)
      enddo

      endif ! if ncrossec
!=======================================================================

  endif ! output_format_ppf

  if (output_format_grib2) then

! Coding-writing output data in grib2 format

    call write_grib2_horizontal_grid (idate0,iperiod,iperiod_accum,idatec, &
      ipflag,itflag,isflag,                                                &
      zph,zt,zu,zv,zq,zrh,zom,zrva,zpva,zthetae,zpvat,                     &
      tg_post, qg_post,                                                    &
      zorogr,zpz0,ztgrs,zqgrs,zus,zvs,zts,zqs,zcape,ztmin,ztmax,zlift,     &
      zlev0, fice, iceth, lapse_rate_1, lapse_rate_2, iwv, zcin, pbl, gust)

! For writing constant (static) fields (fmask, orography)

!    if (kist == 1) then
#ifndef rttov
      call write_grib2_horizontal_grid_static_data(1, nlon, nlat,     &
      pdr(2), pdr(1), pdr(39), pdr(38), pdr(5), pdr(4)+pdr(1)*0.5, 1, &
      idate0, iperiod, iperiod_accum, idatec, zorogr, fmask, 0)
#else
      call write_grib2_horizontal_grid_static_data(1, nlon, nlat,     &
      pdr(2), pdr(1), pdr(39), pdr(38), pdr(5), pdr(4)+pdr(1)*0.5, npoint_jump, &
      idate0, iperiod, iperiod_accum, idatec, zorogr, fmask, 0)
#endif
!    endif

    if (ncrossec) then !  space cross-sections

      do jnc=1,ncrossy
        lon_crossy(jnc,1) = alont(lncr(jnc),1)
        lon_crossy(jnc,2) = alont(lncr(jnc),nlat)
        lat_crossy(jnc,1) = alatt(lncr(jnc),1)
        lat_crossy(jnc,2) = alatt(lncr(jnc),nlat)
      enddo

      do jnc=1,ncrossx
        lon_crossx(jnc,1) = alont(1,ltcr(jnc))
        lon_crossx(jnc,2) = alont(nlon,ltcr(jnc))
        lat_crossx(jnc,1) = alatt(1,ltcr(jnc))
        lat_crossx(jnc,2) = alatt(nlon,ltcr(jnc))
      enddo

#ifdef oper

!!      call write_grib2_cross_space_i_j (idate0, iperiod, iperiod_accum, idatec,                          &
!!       nlon, nlat, npress, ncrossy, ncrossx, exp(zplcr), lon_crossy, lat_crossy, lon_crossx, lat_crossx, &
!!       pdr(39), pdr(38), zspcr, zspcrx,                                                                  &
!!       ztcr, ztecr, zvcr, zucr, zwcr, zpvcr, zrhcr,                                                      &
!!       ztcrx, ztecrx, zucrx, zvcrx, zwcrx, zpvcrx, zrhcrx)

      if (cross_number > 0) then

        call write_grib2_cross_space (idate0, iperiod, iperiod_accum, idatec, cross_number,     &
         npress, npoint_cross_max, npoint_cross(1:cross_number), exp(zplcr),                     &
         lonini_cross(1:cross_number), latini_cross(1:cross_number),                             &
         lonfin_cross(1:cross_number), latfin_cross(1:cross_number),                             &
         pdr(39), pdr(38), ps_cross(:,1:cross_number),                                           &
         theta_cross(:,:,1:cross_number), thetae_cross(:,:,1:cross_number),                      &
         u_tang_cross(:,:,1:cross_number), v_norm_cross(:,:,1:cross_number),                     &
         w_cross(:,:,1:cross_number), pv_cross(:,:,1:cross_number), rh_cross(:,:,1:cross_number))

      else

        print *,'No free cross-section is defined'

      endif

#endif

    endif ! Space cross-sections

    endif ! output_format_grib2

    return
    end
! ######################################################################
      subroutine wrpost(za)

!  Writes horizontal fields (different variables)
!  Data are written on file 80
!  Fields are written as positive integers (4 digits) after finding the best rescaling
!  (zfac: scaling factor; zoffs: offset)

 use mod_postbolam
 use mod_jump

      dimension za(nlon,nlat), izout(nlon,nlat)

!  Def. of geographical parameters for output (t-points)

      zdlat  = pdr(1)
      zdlon  = pdr(2)
      zalat0 = pdr(4) + zdlat/2.
      zalon0 = pdr(5)
      zy0    = pdr(38)
      zx0    = pdr(39)

!  Comp. of accumulation time (hours and minutes) for rain

      iminu = nint((nhist * dtstep)/60.)
      aihour = iminu/60.                ! must be real
      ihour = int(aihour)
      iminu = iminu - ihour*60

!  The following for accumulating on multiple intervals

      ihour = ihour*njump
      iminu = iminu*njump
      aihour = iminu/60.
      ihour = ihour + int(aihour)
      iminu = iminu - int(aihour)*60

!  Redefinition of very small values of za (may be necessary)

      do jlat=1,nlat
      do jlon=1,nlon
      if(abs(za(jlon,jlat)).lt.1.e-20) za(jlon,jlat)=0.
      enddo
      enddo

!  Definition of offset and scaling factors

      zmin = 1.e35
      zmax =-1.e35
      do jlat = 1, nlat
      do jlon = 1, nlon
      zmax = max(zmax,za(jlon,jlat))
      zmin = min(zmin,za(jlon,jlat))
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

      do jlat = 1, nlat
      do jlon = 1, nlon
      izout(jlon,jlat) = nint((za(jlon,jlat)-zoffs)*zfac)
      enddo
      enddo

      write(80,9000) nfdr(5:12),ihour,iminu,zoffs,zfac
      write(80,9200) nlat,nlon,zalat0,zalon0,zdlat,zdlon,zy0,zx0
      write(80,9100) izout

 9000 format(10i10,2(1x,e12.5))
 9100 format(20i4)
 9200 format(2i5,6f11.5)
      return
      end
!#######################################################################
      subroutine wrposx(za,zpz,zplcr,kncr,npress)

!  Writes vert. west-east cross-sections.
!  za is the matrix containing the main output field on cross-sect.;
!  zplcr is the vector ln(p) of cross-section levels;
!  zpz contains surface pressure along the cross-sect.
!  Data are written on file 82
!  Fields are written as positive integers (4 digits) after finding the best rescaling (zfac: scaling factor; zoffs: offset)

 use mod_postbolam

      dimension za(nlon,npress), izout(nlon,200),zpz(nlon),zaux(nlon),zplcr(npress)

!  Def. of geographical parameters for output (t-points)

      zdlat  = pdr(1)
      zdlon  = pdr(2)
      zalat0 = pdr(4) + zdlat/2.
      zalon0 = pdr(5)

!  Comp. of latitude of cross-section (rotated coordinates)

      zalat = zalat0 + float(kncr-1)*zdlat

!  Redefinition of very small values of za (may be necessary)

      do jlev=1,npress
      do jlon=1,nlon
      if(abs(za(jlon,jlev)).lt.1.e-20) za(jlon,jlev)=0.
      enddo
      enddo

!  Definition of offset and scaling factors

      zmin = 1.e35
      zmax = -1.e35
      do jlev = 1, npress
      do jlon= 1, nlon
      zmax = max(zmax,za(jlon,jlev))
      zmin = min(zmin,za(jlon,jlev))
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

      do jlev = 1, npress
      do jlon = 1, nlon
      izout(jlon,jlev)=nint((za(jlon,jlev)-zoffs)*zfac)
      enddo
      enddo

!  Transf. of surf. pressure into ln of surf. press.

      do jlon = 1, nlon
      zaux(jlon) = alog(zpz(jlon))
      enddo

      write(82,9000) nfdr(5:12),zoffs,zfac
      write(82,9200) nlon,nlat,npress,zalon0,zalat0,zdlon,zdlat,zalat
      write(82,9300) zaux
      write(82,9300) zplcr
      write(82,9100) ((izout(j1,j2),j1=1,nlon),j2=1,npress)

 9000 format(8i10,2(1x,e12.5))
 9100 format(20i4)
 9200 format(3i5,5f11.5)
 9300 format(8f9.5)
      return
      end
!#######################################################################
      subroutine wrposy(za,zpz,zplcr,kncr,npress)

!  Writes vert. north-south cross-sections.
!  za is the matrix containing the main output field on cross-sect.;
!  zplcr is the vector ln(p) of cross-section levels;
!  zpz contains surface pressure along the cross-sect.
!  Data are written on file 81
!  Fields are written as positive integers (4 digits) after finding the best rescaling (zfac: scaling factor; zoffs: offset)

 use mod_postbolam

      dimension za(nlat,npress), izout(nlat,200),zpz(nlat),zaux(nlat),zplcr(npress)

!  Def. of geographical parameters for output (t-points)

      zdlat  = pdr(1)
      zdlon  = pdr(2)
      zalat0 = pdr(4) + zdlat/2.
      zalon0 = pdr(5)

!  Comp. of longitude of cross-section (rotated coordinates)

      zalon = zalon0 + float(kncr-1)*zdlon

! Redefinition of very small values of za (may be necessary)

      do jlev=1,npress
      do jlat=1,nlat
      if(abs(za(jlat,jlev)).lt.1.e-20) za(jlat,jlev)=0.
      enddo
      enddo

!  Definition of offset and scaling factors

      zmin = 1.e35
      zmax = -1.e35
      do jlev = 1, npress
      do jlat= 1, nlat
      zmax = max(zmax,za(jlat,jlev))
      zmin = min(zmin,za(jlat,jlev))
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

      do jlev = 1, npress
      do jlat = 1, nlat
      izout(nlat-jlat+1,jlev)=nint((za(jlat,jlev)-zoffs)*zfac)
      enddo
      enddo

!  Transf. of surf. pressure into ln of surf. press.

      do jlat = 1, nlat
      zaux(nlat-jlat+1) = alog(zpz(jlat))
      enddo

      write(81,9000) nfdr(5:12),zoffs,zfac
      write(81,9200) nlat,nlon,npress,zalat0,zalon0,zdlat,zdlon,zalon
      write(81,9300) zaux
      write(81,9300) zplcr
      write(81,9100) ((izout(j1,j2),j1=1,nlat),j2=1,npress)

 9000 format(8i10,2(1x,e12.5))
 9100 format(20i4)
 9200 format(3i5,5f11.5)
 9300 format(8f9.5)
      return
      end
!#######################################################################
      subroutine defout(klat)

!  COMPUTATION OF VARIABLES AT THE GROUND OR AT Z = 0 AND AT
!  STANDARD PRESSURE LEVELS

 use mod_postbolam
 use mod_rplev
 use mod_rwoutp
 use mod_ptvor
 use mod_rwcross

 dimension zalp(nlevpo)
 real, dimension(nlev) :: ztcra, zrhcra, zucra, zvcra, ztecra,zwcra,zpvcra

!  Pressure levels (log(p)) for cross-sections

 real, dimension(nlevp1) :: zaux,zlnp,zlnpu,zlnpv,zlnppv,zlnppvc

      do jlev = 1, nlevpo
      zalp(jlev) = alog(plevo(jlev))
      enddo

      do 600 jlon = 1, nlon
      jlonm1=max(jlon-1,1)

! zpsu is used to compute u - averaged between jlon and jlon+1 and extrap. at east bound.
! zpsv is used to compute v - averaged between klat and klat-1
! zpspv is used to compute pv on z-points and is averaged from 4 surrounding t-points

      zpsu = .5*(ps(jlon,klat) + ps(jlon+1,klat))
      zpspv= .25*(ps(jlon,klat) + ps(jlon+1,klat) + ps(jlon,klat-1) + ps(jlon+1,klat-1))
      zpsv = .5*(ps(jlon,klat) + ps(jlon,klat-1))
      zalps = alog(ps(jlon,klat))
      zalpsu= alog(zpsu)
      zalpsv= alog(zpsv)

      do jklev = 1, nlev
      zlnp(jklev) = alog(pzer*sigint(jklev)-(pzer-ps(jlon,klat))*sigalf(jklev))
      zlnpu(jklev)= alog(pzer*sigint(jklev)-(pzer-zpsu)*sigalf(jklev))
      zlnpv(jklev)= alog(pzer*sigint(jklev)-(pzer-zpsv)*sigalf(jklev))
      zsgalf = sig(jklev)**3*(1.+alfa*(1.-sig(jklev))*(2.-sig(jklev)))
      zlnppv(jklev) =alog(pzer*sig(jklev)-(pzer-zpspv        )*zsgalf)
      zlnppvc(jklev)=alog(pzer*sig(jklev)-(pzer-ps(jlon,klat))*zsgalf)
      enddo

!  Computation of mslp
!  For extrapolation of pressure to the msl, the temperature at level nlev-1
!  is used to avoid problems due to "extreme" values of t near the ground

      zvts(jlon) = temp(jlon,klat,nlev)
      zvtsm = temp(jlon,klat,nlev-1)

!  Average of computed stability with climatological value 6.5 k/km

      zss = 6.5e-03

!  Computation of t at sea level and at ground level with a lapse rate of zss

      zvtz0(jlon) = zvtsm + zss*(phi(jlon,klat,nlev-1)/g) ! T AT Z=0
      zvtg=zvtsm+zss*(zalps-zlnp(nlev-1))*rd/g*zvts(jlon)*(1.+ep*q(jlon,klat,nlev)) ! T AT GROUND

!  Computation of mean sea level pressure (hydrostatically) using approx. virtual temperature,
!  computed with spec. hum. at lowest model level

      ztbarv = (zvtz0(jlon)+zvtg)/2.*(1.+ep*q(jlon,klat,nlev))
      zvpz0(jlon) = ps(jlon,klat)*exp(phig(jlon,klat)/ztbarv/rd)

!  u component at sigma = nlev (lowest sigma)

      zvus(jlon) = u(jlon,klat,nlev)

!  v component at sigma = nlev (lowest sigma)

      zvvs(jlon) = v(jlon,klat,nlev)

!  q at sigma = nlev (lowest sigma)

      zvqs(jlon) = q(jlon,klat,nlev)

!  Vertical interpol. of t, q, u, v, phi, omega, pot. vort.
!  where available, auxiliary levels at ground (or below) are used.
!  In this case nlevp1 levels are considered for input values.
!  The z=0 level is used to interp. t and phi below ground (for consistency between t and phi).
!  For wind and spec. hum., the values at 10m and 2m are used.
!  pv is extrapolated const. below ground, then set to zero;
!  q, u, v and omega are set to zero below ground.

      zlnp(nlevp1) = alog(zvpz0(jlon))

! In case zlnp(nlevp1) <= zlnp(nlev) (lowest model level below sea level) the auxiliary level
! at sea level pressure cannot be used - this concerns vert. interpolation of t and phi only

      if(zlnp(nlevp1).le.zlnp(nlev)) then
      nlevint=nlev
      else
      nlevint=nlevp1
      endif

!  Interpolation of t

      do jklev= 1, nlev
      zaux(jklev) = t(jlon,klat,jklev)
      enddo

!  From t to theta at z=0

      zaux(nlevp1) = zvtz0(jlon)*exp(-rcp*(zlnp(nlevp1) - alp0))

!  Def. parameters for vert. interpolation routine

      zalf= .5
      ze1 = .8
      ze2 = .8

      call interp(zalf,ze1,ze2,nlevint,zlnp,zaux,zalp,zvt(jlon,1:nlevpo),nlevpo)

! Purely linear interpolation below lowest sigma level, with a standard lapse rate for theta,
! to limit problems due to spline interpolation in case of strong inversions

      do jlev = 1, nlevpo
      if(zalp(jlev).gt.zlnp(nlev)) then
      zvt(jlon,jlev)=zaux(nlev)-(3./0.105)*(zalp(jlev)-zlnp(nlev))
      zvt(jlon,jlev)=0.5*zvt(jlon,jlev)+ 0.5*(zaux(nlev-1)-(3./0.105)*(zalp(jlev)-zlnp(nlev-1)))
      endif

!  Transf. from theta to t

      zvt(jlon,jlev) = zvt(jlon,jlev)*exp(rcp*(zalp(jlev) - alp0))
      enddo

!  Interpolation of q. The sqrt of q is used for interp.

      do jklev= 1, nlev
      if (q(jlon,klat,jklev).lt.1.e-12) then
      q(jlon,klat,jklev) = 1.e-12
      endif
      zaux(jklev) = sqrt(q(jlon,klat,jklev))
      enddo
      zaux(nlevp1) = sqrt(q2(jlon,klat))
      zlnp(nlevp1) = zalps

      zalf= .5
      ze1 = .7
      ze2 = .7
      call interp(zalf,ze1,ze2,nlevp1,zlnp,zaux,zalp,zvq(jlon,1:nlevpo),nlevpo)

!  Retransf. to q and reset to 0 below ground

      do jlev = 1, nlevpo
      zvq(jlon,jlev) = max(zvq(jlon,jlev), 0.)
      zvq(jlon,jlev) = zvq(jlon,jlev)**2
      if (zalp(jlev) .gt. zalps) zvq(jlon,jlev) = 0.
      enddo

!  Interpolation of u

      do jklev= 1, nlev
      zaux(jklev) = u(jlon,klat,jklev)
      enddo
      zaux(nlevp1) = 0.5*(u10(jlonm1,klat)+u10(jlon,klat)) ! u10 is defined on T points (in vtsurf)
      zlnpu(nlevp1) = zalpsu

      zalf= .5
      ze1 = .4
      ze2 = .4
      call interp(zalf,ze1,ze2,nlevp1,zlnpu,zaux,zalp,zvu(jlon,1:nlevpo),nlevpo)

!  Reset to 0 below ground

      do jlev = 1, nlevpo
      if (zalp(jlev) .gt. zalpsu) zvu(jlon,jlev) = 0.
      enddo

!  Interpolation of v

      do jklev= 1, nlev
      zaux(jklev) = v(jlon,klat,jklev)
      enddo
      klatp1=min(klat+1,nlat)
      zaux(nlevp1) = 0.5*(v10(jlon,klat)+v10(jlon,klatp1)) ! v10 is defined on T points (in vtsurf)
      zlnpv(nlevp1) = zalpsv
      call interp(zalf,ze1,ze2,nlevp1,zlnpv,zaux,zalp,zvv(jlon,1:nlevpo),nlevpo)

!  Reset to 0 below ground

      do jlev = 1, nlevpo
      if (zalp(jlev) .gt. zalpsv) zvv(jlon,jlev) = 0.
      enddo

!  Interpolation of phi

      do jklev= 1, nlev
      zaux(jklev) = phi(jlon,klat,jklev)
      enddo
      zaux(nlevp1) = 0.
      zlnp(nlevp1) = alog(zvpz0(jlon))

      zalf= .5
      ze1 = 1.
      ze2 = 1.
      call interp(zalf,ze1,ze2,nlevint,zlnp,zaux,zalp,zvph(jlon,1:nlevpo),nlevpo)

!  Skip of omega at west and north bound.; skip also pot. vort. (and cross-sect.) at east bound.

      if(klat .eq. nlat) go to 590
      if(jlon.eq.1) goto 590
      if(jlon .eq. nlon)  go to 595

!  Interpolation of omega

      do jklev= 1, nlev
      zaux(jklev) = omeg(jlon,klat,jklev)
      enddo

      zalf= .5
      ze1 = .5
      ze2 = .5
      call interp(zalf,ze1,ze2,nlev,zlnp,zaux,zalp,zvom(jlon,1:nlevpo),nlevpo)

!  Reset to 0 below ground

      do jlev = 1, nlevpo
      if (zalp(jlev) .gt. zalps) zvom(jlon,jlev) = 0.
      enddo

 590  continue

!  Interpolation of pot. vort

      do jklev= 1, nlev
      zaux(jklev) = potv(jlon,klat,jklev)
      enddo

      zalf= .5
      ze1 = .2
      ze2 = .2
      call interp(zalf,ze1,ze2,nlev,zlnppv,zaux,zalp,zvpv(jlon,1:nlevpo),nlevpo)

!  Reset to 0 below last model level

      do jlev = 1, nlevpo
      if (zalp(jlev).gt.zlnppv(nlev)) zvpv(jlon,jlev)=0.
      enddo

 595  continue

!=======================================================================
!  Definition of vertical vectors for cross-section variables.
!  Transf. from specific to relative humidity and computation of equivalent potential temperature
!  are done here on sigma levels, before vertical interpolation.
!  Note that for cross sections the horizontal interp. is done before vertical interp. for
!   staggered variables (u, v, potv), while the opposite holds for fields on pressure levels

      do 400 jklev = 1, nlev
      ztcra(jklev) = t(jlon,klat,jklev)
      ztmpr = temp(jlon,klat,jklev)

!  Quantities to compute rel. hum. and equiv. pot. temp.
!  zpres: air pressure; ze: vapour pressure; zmr: mix. ratio
!  ztl: temp. at condens. lev. (computed with expr. by Bolton)

      zpres = pzer*sigint(jklev)-(pzer-ps(jlon,klat))*sigalf(jklev)
      zqq = q(jlon,klat,jklev)
      zmr = zqq/(1.-zqq)
      ze = zpres/100. * zmr/(eps+zmr)
      if(ze.lt.1.e-13) ze = 1.e-13

! qtorh does not limit zrh to saturation because it is used in cross-sections to extract cloud water+ice

      call qtorh(zqq,ztmpr,zpres,zrh)
      zrhcra(jklev)= zrh
      if(jlon.gt.1) then
      zucra(jklev) = .5*(u(jlon,klat,jklev)+u(jlon-1,klat,jklev))
      else
      zucra(jklev) = .5*(3.*u(1,klat,jklev)-u(2,klat,jklev))
      endif

      if(klat.lt.nlat) then
      zvcra(jklev) = .5*(v(jlon,klat,jklev)+v(jlon,klat+1,jklev))
      else
      zvcra(jklev) = .5*(3.*v(jlon,nlat,jklev)-v(jlon,nlatm1,jklev))
      endif

!  In the following the vert. vel. is comp. as -omega/p

      zwcra(jklev) = -omeg(jlon,klat,jklev)/zpres

!  Definition of pot. vorticity (interpolated on t-points)

      if(jlon.eq.1) then
      zpvcra(jklev) = .5*(potv(1,klat,jklev)+potvn(1,jklev))
      elseif(jlon.eq.nlon) then
      zpvcra(jklev) = .5*(potv(nlon-1,klat,jklev)+potvn(nlon-1,jklev))
      else
      zpvcra(jklev) = .25*(potv (jlon,klat,jklev)+potv (jlon-1,klat,jklev)+potvn(jlon,jklev)+potvn(jlon-1,jklev))
      endif

! Limitation of mixing ratio to subsaturation for comput. of equivalent pot. temp.
! Saturation is computed using the blended curve below freezing

      call comp_esk(zzpvs,zzqs,ztmpr,zpres,2)

      zzpv = (zqq*zpres)/(0.622+zqq-0.622*zqq)
      if(zzpv.gt.zzpvs) then
       zzpv = zzpvs
       zmr = 0.622*zzpv/(zpres-zzpv)
       ze = zpres/100.*zmr/(eps+zmr)
      endif

!  Definition of equivalent pot. temp. for cross-sections
!  and for const. pressure maps.
!  The definition follows Bolton (MWR,1980)
!  ztl: temp. at condens. lev. (computed with expr. by Bolton)

      ztl = 55. + 2840./(3.5*alog(ztmpr) - alog(ze) - 4.805)
      zespon = (3.376/ztl - 0.00254)*zmr*1000*(1.+.81*zmr)
      ztecra(jklev) = ztcra(jklev)*exp(zespon)
      zaux(jklev)=ztecra(jklev)
 400  continue

! Vertical interp. of thetae for const. press. maps

      zalf= .5
      ze1 = .0
      ze2 = .0
      call interp(zalf,ze1,ze2,nlev,zlnp,zaux,zalp,zvthetae(jlon,1:nlevpo),nlevpo)

!---------------------------------------
!   Definition of cross section matrices
!---------------------------------------

      do jnc=1,ncrossy
      if(jlon.eq.lncr(jnc)) then
      zspcr(klat,jnc) = ps(jlon,klat)

!   Interpolation of vertical vectors for cross-sect. variables from sigma levels to pressure levels
!   Note that only model levels are used (not values at the ground)

      zalf = .5
      call interp(zalf,.6,.6,nlev,zlnp,ztcra,    zplcr,ztcr (klat,1:npress,jnc),npress)
      call interp(zalf,.6,.6,nlev,zlnp,zrhcra,   zplcr,zrhcr(klat,1:npress,jnc),npress)
      call interp(zalf,.6,.6,nlev,zlnp,zucra,    zplcr,zucr (klat,1:npress,jnc),npress)
      call interp(zalf,.6,.6,nlev,zlnp,zvcra,    zplcr,zvcr (klat,1:npress,jnc),npress)
      call interp(zalf,.6,.6,nlev,zlnp,ztecra,   zplcr,ztecr(klat,1:npress,jnc),npress)
      call interp(zalf,.6,.6,nlev,zlnp,zwcra,    zplcr,zwcr (klat,1:npress,jnc),npress)
      call interp(zalf,.6,.6,nlev,zlnppvc,zpvcra,zplcr,zpvcr(klat,1:npress,jnc),npress)
      endif
      enddo

      do jnc=1,ncrossx
      if(klat.eq.ltcr(jnc)) then
      zspcrx(jlon,jnc) = ps(jlon,klat)
      zalf = .5
      call interp(zalf,.6,.6,nlev,zlnp,ztcra,    zplcr,ztcrx (jlon,1:npress,jnc),npress)
      call interp(zalf,.6,.6,nlev,zlnp,zrhcra,   zplcr,zrhcrx(jlon,1:npress,jnc),npress)
      call interp(zalf,.6,.6,nlev,zlnp,zucra,    zplcr,zucrx (jlon,1:npress,jnc),npress)
      call interp(zalf,.6,.6,nlev,zlnp,zvcra,    zplcr,zvcrx (jlon,1:npress,jnc),npress)
      call interp(zalf,.6,.6,nlev,zlnp,ztecra,   zplcr,ztecrx(jlon,1:npress,jnc),npress)
      call interp(zalf,.6,.6,nlev,zlnp,zwcra,    zplcr,zwcrx (jlon,1:npress,jnc),npress)
      call interp(zalf,.6,.6,nlev,zlnppvc,zpvcra,zplcr,zpvcrx(jlon,1:npress,jnc),npress)
      endif
      enddo
!=======================================================================

 600  continue

      return
      end
!#######################################################################
subroutine defout_tropopause

! Definition of variables at the Tropopause level

 use mod_postbolam, only : nlon, nlat, nlev, phi, g, pzer, sigint, sigalf, ps, temp, u, v, q, omeg, t, &
 dy, dx, hxt, eps, itropflag, &
 h_tropopause, p_tropopause, t_tropopause, u_tropopause, v_tropopause, q_tropopause,&
 rh_tropopause, om_tropopause, rv_tropopause, pv_tropopause, the_tropopause
 use mod_ptvor, only : potv

 implicit none

 real, dimension(nlev) :: zh, zprofile
 real, dimension(1) :: zout
 integer, dimension(1) :: iv
 integer :: jlon, jlat, jklev, jklev1
 real :: zalf=0.5, ze1=0.2, ze2=0.2, zqs, zpvs, zpv, zu, zu1, zv, zv1, zp, ze, zmr, ztl, zespon

 if (itropflag(10) == 1.or.itropflag(5) == 1) then ! Pressure

   do jlat = 1, nlat
   do jlon = 1, nlon

     do jklev = 1, nlev
       jklev1 = nlev-jklev+1
       zh(jklev) = phi(jlon,jlat,jklev1)/g
       zprofile(jklev) = pzer*sigint(jklev1)-(pzer-ps(jlon,jlat))*sigalf(jklev1)
     enddo

     call near (h_tropopause(jlon,jlat), 1, zh, nlev, iv(1:1))
     call interp_spline_1d(zout(1:1), h_tropopause(jlon,jlat), 1, zprofile, zh, nlev, iv(1:1), zalf, ze1, ze2)
     p_tropopause(jlon,jlat) = zout(1)

   enddo
   enddo

 endif

 if (itropflag(2) == 1.or.itropflag(5) == 1) then ! Temperature

   do jlat = 1, nlat
   do jlon = 1, nlon

     do jklev = 1, nlev
       jklev1 = nlev-jklev+1
       zh(jklev) = phi(jlon,jlat,jklev1)/g
       zprofile(jklev) = temp(jlon,jlat,jklev1)
     enddo

     call near (h_tropopause(jlon,jlat), 1, zh, nlev, iv(1:1))
     call interp_spline_1d(zout(1:1), h_tropopause(jlon,jlat), 1, zprofile, zh, nlev, iv(1:1), zalf, ze1, ze2)
     t_tropopause(jlon,jlat) = zout(1)

   enddo
   enddo

 endif

 if (itropflag(3) == 1.or.itropflag(7) == 1) then ! Wind components

   do jlat = 1, nlat
   do jlon = 1, nlon

     do jklev = 1, nlev
       jklev1 = nlev-jklev+1
       zh(jklev) = phi(jlon,jlat,jklev1)/g
       if (jlon /= 1) then
         zprofile(jklev) = (u(jlon-1,jlat,jklev1)+u(jlon,jlat,jklev))*0.5
       else
         zprofile(jklev) = u(jlon,jlat,jklev)
       endif
     enddo

     call near (h_tropopause(jlon,jlat), 1, zh, nlev, iv(1:1))
     call interp_spline_1d(zout(1:1), h_tropopause(jlon,jlat), 1, zprofile, zh, nlev, iv(1:1), zalf, ze1, ze2)
     u_tropopause(jlon,jlat) = zout(1)

     do jklev = 1, nlev
       jklev1 = nlev-jklev+1
       if (jlat /= nlat) then
         zprofile(jklev) = (v(jlon,jlat,jklev1)+v(jlon,jlat+1,jklev))*0.5
       else
         zprofile(jklev) = v(jlon,jlat,jklev)
       endif
     enddo

     call near (h_tropopause(jlon,jlat), 1, zh, nlev, iv(1:1))
     call interp_spline_1d(zout(1:1), h_tropopause(jlon,jlat), 1, zprofile, zh, nlev, iv(1:1), zalf, ze1, ze2)
     v_tropopause(jlon,jlat) = zout(1)

   enddo
   enddo

 endif

 if (itropflag(4) == 1.or.itropflag(5) == 1) then ! Specific humidity

   do jlat = 1, nlat
   do jlon = 1, nlon

     do jklev = 1, nlev
       jklev1 = nlev-jklev+1
       zh(jklev) = phi(jlon,jlat,jklev1)/g
       zprofile(jklev) = sqrt(q(jlon,jlat,jklev1))
     enddo

     call near (h_tropopause(jlon,jlat), 1, zh, nlev, iv(1:1))
     call interp_spline_1d(zout(1:1), h_tropopause(jlon,jlat), 1, zprofile, zh, nlev, iv(1:1), zalf, ze1, ze2)
     q_tropopause(jlon,jlat) = zout(1)**2

   enddo
   enddo

 endif

 if (itropflag(5) == 1) then ! Relative humidity

   do jlat = 1, nlat
   do jlon = 1, nlon

     call comp_esk(zpvs,zqs,t_tropopause(jlon,jlat),p_tropopause(jlon,jlat),1) ! computes saturation to water and ice separately
     zpv=(q_tropopause(jlon,jlat)*p_tropopause(jlon,jlat))/(0.622+q_tropopause(jlon,jlat)-0.622*q_tropopause(jlon,jlat))
     rh_tropopause(jlon,jlat) = min(zpv/zpvs*100., 102.)

   enddo
   enddo

 endif

 if (itropflag(6) == 1) then ! Vertical velocity

   do jlat = 1, nlat
   do jlon = 1, nlon

     do jklev = 1, nlev
       jklev1 = nlev-jklev+1
       zh(jklev) = phi(jlon,jlat,jklev1)/g
       zprofile(jklev) = omeg(jlon,jlat,jklev1)
     enddo

     call near (h_tropopause(jlon,jlat), 1, zh, nlev, iv(1:1))
     call interp_spline_1d(zout(1:1), h_tropopause(jlon,jlat), 1, zprofile, zh, nlev, iv(1:1), zalf, ze1, ze2)
     om_tropopause(jlon,jlat) = zout(1)

   enddo
   enddo

 endif

 if (itropflag(7) == 1) then ! Relative vorticity

   rv_tropopause(jlon,jlat) = 0.

   do jlat = 2, nlat
   do jlon = 1, nlon-1

     zu = 0.5*(u_tropopause(jlon,jlat)+u_tropopause(jlon+1,jlat))
     zu1 = 0.5*(u_tropopause(jlon,jlat-1)+u_tropopause(jlon+1,jlat-1))
     zv = 0.5*(v_tropopause(jlon,jlat)+v_tropopause(jlon,jlat-1))
     zv1 = 0.5*(v_tropopause(jlon+1,jlat)+v_tropopause(jlon+1,jlat-1))
     rv_tropopause(jlon,jlat) = (-(zu*hxt(jlat)-zu1*hxt(jlat-1))/dy+ &
 (zv1-zv)/dx)/(hxt(jlat)+hxt(jlat-1))*2.

   enddo
   enddo

 endif

 if (itropflag(8) == 1) then ! Potential vorticity

   do jlat = 1, nlat
   do jlon = 1, nlon

     do jklev = 1, nlev-1
       jklev1 = nlev-jklev+1
       zh(jklev1) = 0.5*(phi(jlon,jlat,jklev)+phi(jlon,jlat,jklev+1))/g
     enddo
     zh(1) = phi(jlon,jlat,nlev)/g*0.5

     do jklev = 1, nlev
       jklev1 = nlev-jklev+1
       if (jlon /= 1.and.jlat /= nlat) then
         zprofile(jklev) = 0.25* &
 (potv(jlon,jlat,jklev1)+potv(jlon-1,jlat,jklev1)+potv(jlon,jlat+1,jklev1)+potv(jlon-1,jlat+1,jklev1))
       elseif (jlon == 1.and.jlat /= nlat) then
         zprofile(jklev) = 0.50*(potv(jlon,jlat,jklev1)+potv(jlon,jlat+1,jklev1))
       elseif (jlon /= 1.and.jlat == nlat) then
         zprofile(jklev) = 0.50*(potv(jlon,jlat,jklev1)+potv(jlon-1,jlat,jklev1))
       else
         zprofile(jklev) = potv(jlon,jlat,jklev1)
       endif
     enddo

     call near (h_tropopause(jlon,jlat), 1, zh, nlev, iv(1:1))
     call interp_spline_1d(zout(1:1), h_tropopause(jlon,jlat), 1, zprofile, zh, nlev, iv(1:1), zalf, ze1, ze2)
     pv_tropopause(jlon,jlat) = zout(1)

   enddo
   enddo

 endif

 if (itropflag(9) == 1) then ! Equivalent Potential Temperature

   do jlat = 1, nlat
   do jlon = 1, nlon

     do jklev = 1, nlev
       jklev1 = nlev-jklev+1
       zh(jklev1) = phi(jlon,jlat,jklev)/g
       zp = pzer*sigint(jklev)-(pzer-ps(jlon,jlat))*sigalf(jklev)
       zmr = q(jlon,jlat,jklev)/(1.-q(jlon,jlat,jklev))
       ze = max( zp/100. * zmr/(eps+zmr), 1.e-13)
       call comp_esk(zpvs,zqs,temp(jlon,jlat,jklev),zp,2)
       zpv = (q(jlon,jlat,jklev)*zp)/(0.622+q(jlon,jlat,jklev)-0.622*q(jlon,jlat,jklev))
       if (zpv > zpvs) then
        zpv = zpvs
        zmr = 0.622*zpv/(zp-zpv)
        zpv = zp/100.*zmr/(eps+zmr)
       endif
       ztl = 55. + 2840./(3.5*alog(temp(jlon,jlat,jklev)) - alog(ze) - 4.805)
       zespon = (3.376/ztl - 0.00254)*zmr*1000*(1.+.81*zmr)
       zprofile(jklev1) = t(jlon,jlat,jklev)*exp(zespon)
     enddo

     call near (h_tropopause(jlon,jlat), 1, zh, nlev, iv(1:1))
     call interp_spline_1d(zout(1:1), h_tropopause(jlon,jlat), 1, zprofile, zh, nlev, iv(1:1), zalf, ze1, ze2)
     the_tropopause(jlon,jlat) = zout(1)

   enddo
   enddo

 endif

return
end
!#######################################################################
subroutine defout_maxwind

! Definition of variables at the level of maximum wind speed

 use mod_postbolam, only : nlon, nlat, nlev, phi, g, pzer, sigint, sigalf, ps, temp, u, v, q, omeg, t, &
 dy, dx, hxt, eps, imwflag, &
 klev_maxwind, h_maxwind, p_maxwind, t_maxwind, u_maxwind, v_maxwind, q_maxwind,&
 rh_maxwind, om_maxwind, rv_maxwind, pv_maxwind, the_maxwind
 use mod_ptvor, only : potv

 implicit none

 real, dimension(nlon,nlat) :: zwork
 integer :: nsmooth=1, jlon, jlat, jklev
 real :: wei=0.5, zqs, zpvs, zpv, zu, zu1, zv, zv1, pv_up, pv_down, zp, ze, zmr, ztl, zespon

 if (imwflag(1) == 1) then ! Altitude

   do jlat = 1, nlat
   do jlon = 1, nlon

     jklev = klev_maxwind(jlon,jlat)
     if (jklev /= 1) then
       h_maxwind(jlon,jlat) = 0.5*(phi(jlon,jlat,jklev)+phi(jlon,jlat,jklev-1))/g
     else
       h_maxwind(jlon,jlat) = phi(jlon,jlat,jklev)/g
     endif

   enddo
   enddo

   call smooth(h_maxwind, zwork, nlon, nlat, wei, nsmooth)
   h_maxwind(:,:) = zwork(:,:)

 endif

 if (imwflag(10) == 1.or.imwflag(5) == 1) then ! Pressure

   do jlat = 1, nlat
   do jlon = 1, nlon

     jklev = klev_maxwind(jlon,jlat)
     p_maxwind(jlon,jlat) = pzer*sigint(jklev)-(pzer-ps(jlon,jlat))*sigalf(jklev)

   enddo
   enddo

   call smooth(p_maxwind, zwork, nlon, nlat, wei, nsmooth)
   p_maxwind(:,:) = zwork(:,:)

 endif

 if (imwflag(2) == 1.or.imwflag(5) == 1) then ! Temperature

   do jlat = 1, nlat
   do jlon = 1, nlon

     jklev = klev_maxwind(jlon,jlat)
     t_maxwind(jlon,jlat) = temp(jlon,jlat,jklev)

   enddo
   enddo

   call smooth(t_maxwind, zwork, nlon, nlat, wei, nsmooth)
   t_maxwind(:,:) = zwork(:,:)

 endif

 if (imwflag(3) == 1.or.imwflag(7) == 1) then ! Wind components

   do jlat = 1, nlat
   do jlon = 1, nlon

     jklev = klev_maxwind(jlon,jlat)
     if (jlon /= 1) then
       u_maxwind(jlon,jlat) = (u(jlon-1,jlat,jklev)+u(jlon,jlat,jklev))*0.5
     else
       u_maxwind(jlon,jlat) = u(jlon,jlat,jklev)
     endif
     if (jlat /= nlat) then
       v_maxwind(jlon,jlat) = (v(jlon,jlat,jklev)+v(jlon,jlat+1,jklev))*0.5
     else
       v_maxwind(jlon,jlat) = v(jlon,jlat,jklev)
     endif

   enddo
   enddo

   call smooth(u_maxwind, zwork, nlon, nlat, wei, nsmooth)
   u_maxwind(:,:) = zwork(:,:)
   call smooth(v_maxwind, zwork, nlon, nlat, wei, nsmooth)
   v_maxwind(:,:) = zwork(:,:)

 endif

 if (imwflag(4) == 1.or.imwflag(5) == 1) then ! Specific humidity

   do jlat = 1, nlat
   do jlon = 1, nlon

     jklev = klev_maxwind(jlon,jlat)
     q_maxwind(jlon,jlat) = q(jlon,jlat,jklev)

   enddo
   enddo

   call smooth(q_maxwind, zwork, nlon, nlat, wei, nsmooth)
   q_maxwind(:,:) = zwork(:,:)

 endif

 if (imwflag(5) == 1) then ! Relative humidity

   do jlat = 1, nlat
   do jlon = 1, nlon

     call comp_esk(zpvs,zqs,t_maxwind(jlon,jlat),p_maxwind(jlon,jlat),1) ! computes saturation to water and ice separately
     zpv=(q_maxwind(jlon,jlat)*p_maxwind(jlon,jlat))/(0.622+q_maxwind(jlon,jlat)-0.622*q_maxwind(jlon,jlat))
     rh_maxwind(jlon,jlat) = min(zpv/zpvs*100., 102.)

   enddo
   enddo

 endif

 if (imwflag(6) == 1) then ! Vertical velocity

   do jlat = 1, nlat
   do jlon = 1, nlon

     jklev = klev_maxwind(jlon,jlat)
     om_maxwind(jlon,jlat) = omeg(jlon,jlat,jklev)

   enddo
   enddo

   call smooth(om_maxwind, zwork, nlon, nlat, wei, nsmooth)
   om_maxwind(:,:) = zwork(:,:)

 endif

 if (imwflag(7) == 1) then ! Relative vorticity

   rv_maxwind(jlon,jlat) = 0.

   do jlat = 2, nlat
   do jlon = 1, nlon-1

     zu = 0.5*(u_maxwind(jlon,jlat)+u_maxwind(jlon+1,jlat))
     zu1 = 0.5*(u_maxwind(jlon,jlat-1)+u_maxwind(jlon+1,jlat-1))
     zv = 0.5*(v_maxwind(jlon,jlat)+v_maxwind(jlon,jlat-1))
     zv1 = 0.5*(v_maxwind(jlon+1,jlat)+v_maxwind(jlon+1,jlat-1))
     rv_maxwind(jlon,jlat) = (-(zu*hxt(jlat)-zu1*hxt(jlat-1))/dy+ &
 (zv1-zv)/dx)/(hxt(jlat)+hxt(jlat-1))*2.

   enddo
   enddo

 endif

 if (imwflag(8) == 1) then ! Potential vorticity

   do jlat = 1, nlat
   do jlon = 1, nlon

     jklev = klev_maxwind(jlon,jlat)

     if (jlon /= 1.and.jlat /= nlat) then
       pv_down = 0.25*(potv(jlon,jlat,jklev)+potv(jlon-1,jlat,jklev)+potv(jlon,jlat+1,jklev)+potv(jlon-1,jlat+1,jklev))
       if (jklev > 1 ) &
 pv_up = 0.25*(potv(jlon,jlat,jklev-1)+potv(jlon-1,jlat,jklev-1)+potv(jlon,jlat+1,jklev-1)+potv(jlon-1,jlat+1,jklev-1))
     elseif (jlon == 1.and.jlat /= nlat) then
       pv_down = 0.50*(potv(jlon,jlat,jklev)+potv(jlon,jlat+1,jklev))
       if (jklev > 1 ) pv_up = 0.50*(potv(jlon,jlat,jklev-1)+potv(jlon,jlat+1,jklev-1))
     elseif (jlon /= 1.and.jlat == nlat) then
       pv_down = 0.50*(potv(jlon,jlat,jklev)+potv(jlon-1,jlat,jklev))
       if (jklev > 1 ) pv_up = 0.50*(potv(jlon,jlat,jklev-1)+potv(jlon-1,jlat,jklev-1))
     else
       pv_down = potv(jlon,jlat,jklev)
       if (jklev > 1 ) pv_up = potv(jlon,jlat,jklev-1)
     endif
     if (jklev == 1) pv_up = 0.

     pv_maxwind(jlon,jlat) = 0.5*(pv_up+pv_down)

   enddo
   enddo

   call smooth(pv_maxwind, zwork, nlon, nlat, wei, nsmooth)
   pv_maxwind(:,:) = zwork(:,:)

 endif

 if (imwflag(9) == 1) then ! Equivalent Potential Temperature

   do jlat = 1, nlat
   do jlon = 1, nlon

     jklev = klev_maxwind(jlon,jlat)

     zp = pzer*sigint(jklev)-(pzer-ps(jlon,jlat))*sigalf(jklev)
     zmr = q(jlon,jlat,jklev)/(1.-q(jlon,jlat,jklev))
     ze = max( zp/100. * zmr/(eps+zmr), 1.e-13)
     call comp_esk(zpvs,zqs,temp(jlon,jlat,jklev),zp,2)
     zpv = (q(jlon,jlat,jklev)*zp)/(0.622+q(jlon,jlat,jklev)-0.622*q(jlon,jlat,jklev))
     if (zpv > zpvs) then
      zpv = zpvs
      zmr = 0.622*zpv/(zp-zpv)
      zpv = zp/100.*zmr/(eps+zmr)
     endif
     ztl = 55. + 2840./(3.5*alog(temp(jlon,jlat,jklev)) - alog(ze) - 4.805)
     zespon = (3.376/ztl - 0.00254)*zmr*1000*(1.+.81*zmr)
     the_maxwind(jlon,jlat) = t(jlon,jlat,jklev)*exp(zespon)

   enddo
   enddo

 endif

return
end
!#######################################################################
      subroutine potvor(klat,potv)

!  Computes Ertel potential vorticity for the postprocessing.
!  Pot. vort. is computed at sigma semi-integer levels (from level 1+1/2 to level nlev-1/2)
!  and at z-points of the horizontal grid.
!  Note that pot. vort. at nlev (corresp. to nlev+1/2, that is the ground surface) cannot be
!  computed and is set equal to level nlev-1  (corresp. to nlev - 1/2).
!  Pot. vort. is not computed along southern and eastern boundaries.

 use mod_postbolam

      dimension potv(nlon,nlev)
      dimension zdthetx(nlon,nlev), zdthety(nlon,nlev)
      dimension zdthxh(nlon,nlevm1), zdthyh(nlon,nlevm1)
      dimension zita(nlon,nlev), zitah(nlon,nlevm1), zplnt(nlon)

!  Definition of absolute vorticity zita (divided by dpk as in the model) on integer sigma levels

      zhxt  = hxt(klat)
      zhxv  = hxv(klat)
      zhxts = hxt(klat-1)
      zhxvs = hxv(klat-1)

      do jlon = 1, nlonm1
      zplnt(jlon) = fv(jlon,klat)*.5*(zhxt+zhxts)
      enddo

      do jklev = 1, nlev
      do jlon = 1, nlonm1
      zita(jlon,jklev) = -(u(jlon,klat,jklev)*zhxt-u(jlon,klat-1,jklev)*zhxts)              &
                          /dy+(v(jlon+1,klat,jklev)-v(jlon,klat,jklev))/dx + zplnt(jlon)
      enddo
      enddo

!  Comp. of dtheta/dx and dtheta/dy (averaged horizontally)

      do jklev = 1, nlev
      do jlon = 1, nlonm1
      zdthetx(jlon,jklev) = ((t(jlon+1,klat,jklev) + t(jlon+1,klat-1,jklev))/2.              &
                           - (t(jlon,klat,jklev) + t(jlon,klat-1,jklev))/2.)/dx/zhxv
      zdthety(jlon,jklev) = ((t(jlon+1,klat,jklev) + t(jlon,klat,jklev))/2.                  &
                           - (t(jlon+1,klat-1,jklev) + t(jlon,klat-1,jklev))/2.)/dy
      enddo
      enddo

!  Averaging at half levels the quantities previously computed over
!  integer levels (a simple arithmetic average is computed)

      do jklev = 1, nlevm1
      do jlon = 1, nlonm1
      zitah(jlon,jklev) = (zita(jlon,jklev)+zita(jlon,jklev+1))/2.
      zdthxh(jlon,jklev) = (zdthetx(jlon,jklev)+zdthetx(jlon,jklev+1))/2.
      zdthyh(jlon,jklev) = (zdthety(jlon,jklev)+zdthety(jlon,jklev+1))/2.
      enddo
      enddo

!  Computation of potential vorticity

      do jklev = 1,nlevm1
      zsgal=sig(jklev)**2*(3.+alfa*(6.-12.*sig(jklev)+5.*sig(jklev)**2))
      do jlon = 1, nlonm1
      zdthds = ((t(jlon,klat,jklev+1)+t(jlon+1,klat,jklev+1)+t(jlon,klat-1,jklev+1)+t(jlon+1,klat-1,jklev+1))/4.- &
               (t(jlon,klat,jklev)+t(jlon+1,klat,jklev)+t(jlon,klat-1,jklev)+t(jlon+1,klat-1,jklev))/4.)/  &
                dsgint(jklev)
      zdhxuds = ((zhxt*(u(jlon,klat,jklev+1)-u(jlon,klat,jklev))) +                            &
                zhxts*(u(jlon,klat-1,jklev+1)-u(jlon,klat-1,jklev))) / (2.* dsgint(jklev))
      zdvds= ((v(jlon+1,klat,jklev+1)+v(jlon,klat,jklev+1))/2. -                               &
             (v(jlon+1,klat,jklev)+v(jlon,klat,jklev))/2.)/dsgint(jklev)

      zpp1=pzer-(pzer-ps(jlon+1,klat  ))*zsgal
      zpp2=pzer-(pzer-ps(jlon  ,klat  ))*zsgal
      zpp3=pzer-(pzer-ps(jlon+1,klat-1))*zsgal
      zpp4=pzer-(pzer-ps(jlon  ,klat-1))*zsgal
      zfac=4./(zhxt*(zpp1+zpp2)+zhxts*(zpp3+zpp4))

      potv(jlon,jklev) = -g*zfac*(zitah(jlon,jklev)*zdthds+zdhxuds*zdthyh(jlon,jklev)-zdvds*zdthxh(jlon,jklev))
      enddo
      enddo

!  Defin. of bottom layer

      do jlon = 1, nlonm1
      potv(jlon,nlev) = potv(jlon,nlevm1)
      enddo

      return
      end
!#######################################################################
       subroutine qtorh (q,t,p,rh)

!  Transforms specif. hum. (q in input) into rel. humidity (0 < 1, but allowing oversaturation)
!  t and p (input) are input temp. and press.
!  Satur. press. is computed with respect to water and ice separately

       parameter (eps=0.622)
       call comp_esk(esat,zqsat,t,p,1)
       e = p*q/(eps-eps*q+q)
!       rh = min(e/esat, 1.01)
       rh = e/esat
       return
       end
!#######################################################################
    subroutine vtsurf

! Postprocessing routine - interpolates wind at 10 m, temperature and spec. humid. at 2 m.
! A reduced roughness over land is used to better simulate measurement conditions
! The second level above ground is used for momentum roughness > 2 m

    use mod_postbolam
    use mod_surfadd

    real, parameter :: zgam=16., zaholt=1., zbholt=.667, zcholt=5., zdholt=.35
    real :: zwork(nlon,nlat)

!  Businger functions

    psium(zz1,zz2) = log( (1.+zz1)**2*(1.+zz1**2)/((1.+zz2)**2*(1.+zz2**2)) ) -2.*(atan(zz1)-atan(zz2))
    psiuh(zz1,zz2) = 2.*log((1.+zz1**2)/(1.+zz2**2))

!  Holtslag functions

    psism(zz1,zz2) = -zaholt*zz1-zbholt*(zz1-zcholt/zdholt)*exp(-zdholt*zz1)  &
                     +zaholt*zz2+zbholt*(zz2-zcholt/zdholt)*exp(-zdholt*zz2)
    psish(zz1,zz2) = -(1.+2./3.*zaholt*zz1)**1.5-zbholt*(zz1-zcholt/zdholt)*exp(-zdholt*zz1) &
                     +(1.+2./3.*zaholt*zz2)**1.5+zbholt*(zz2-zcholt/zdholt)*exp(-zdholt*zz2)

! Loop over horizontal grid

  do jlat = 1, nlat
    jlatp1 = min(jlat+1,nlat)

    do jlon = 1, nlon
      jlonm1= max(1,jlon-1)

! Definition of potential temp. corresponding to the ground temperature

      zconv = exp( rcp*(alps(jlon,jlat)-alp0) )
      ztetg = tskin(jlon,jlat)/zconv
      zthvg = ztetg*(1.+ep*qskin(jlon,jlat))

      if (rgm(jlon,jlat).gt.1.8) then
      nlast = nlev-1
      else
      nlast = nlev
      endif
      zph   = phi(jlon,jlat,nlast)-phig(jlon,jlat)
      zua   = .5*(u(jlon,jlat  ,nlast)+u(jlonm1,jlat,nlast))
      zva   = .5*(v(jlon,jlatp1,nlast)+v(jlon  ,jlat,nlast))
      zdq   = q(jlon,jlat,nlast)-qskin(jlon,jlat)
      zdte  = t(jlon,jlat,nlast)-ztetg
      zzpp  = pzer*sigint(nlast)-(pzer-ps(jlon,jlat))*sigalf(nlast)
      zthvn = tvirt(jlon,jlat,nlast)*(1.e5/zzpp)**rcp
      zza   = zph/g + rgm(jlon,jlat)
      zmod2  = zua**2 + zva**2 + 0.07
      zmod  = sqrt(zmod2)

! Bulk Richardson number

      ztbarv = 0.5*(zthvg+zthvn)
      zri = g*zza*(zthvn-zthvg)/(ztbarv*zmod2)

! Definition of roughness zrgm, zrgt, zrgq

        if (fmask(jlon,jlat).gt..5.and.fice(jlon,jlat).lt.0.1) then ! Charnok roughness
        zchar = 5.e-4
        zcoch1=.0185*ak**2/g*zmod2
        do jiter = 1, 5
        zcoch2= zza/zchar
        zchar = zcoch1/log(zcoch2)**2
        enddo
        zrgm = zchar

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

        else  ! Local roughness set to a few cm if larger
        zsea = 5.e-5
        zrgm = fmask(jlon,jlat)*zsea+(1.-fmask(jlon,jlat))*min(rgm(jlon,jlat), .05)
        zrgt = fmask(jlon,jlat)*zsea+(1.-fmask(jlon,jlat))*min(rgq(jlon,jlat), .05)
        zrgq = zrgt
        endif

      zalzam = log(zza/zrgm)
      zalzat = log(zza/zrgt)
      zalzaq = log(zza/zrgq)

      if(zri.gt.0.) then         ! Holtslag functions
        zal = .2 + 4.*zri
        zpsim = psism(zal, zal*zrgm/zza)
        zpsih = psish(zal, zal*zrgt/zza)
          do jiter = 1, 10
          zalm = zal
          zal = zri*(zalzam-zpsim)**2/(zalzat-zpsih)
          error = abs(zal-zalm)/zal
          if (error.lt.1.e-2) go to 1
          zpsim = psism(zal, zal*zrgm/zza)
          zpsih = psish(zal, zal*zrgt/zza)
          enddo
 1      continue
        zpsiq   = zpsih          ! because rgt=rgq in the stable case
        zpsim10 = psism(zal*(10.+zrgm)/zza, zal*zrgm/zza)
        zpsih2  = psish(zal*(2. +zrgt)/zza, zal*zrgt/zza)
        zpsiq2  = zpsih2
        zpsih05 = psish(zal*(.5 +zrgt)/zza, zal*zrgt/zza)
        zpsih005= psish(zal*(.05+zrgt)/zza, zal*zrgt/zza)
        zpsiq05 = zpsih05
      else                       ! Businger functions
        zpsim = 0.
        zpsih = 0.
          do jiter = 1, 4
          zal = zri*(zalzam-zpsim)**2/(zalzat-zpsih)   ! za/l from eq (6.47) of Businger
          x1 = (1.-zgam*zal)**.25
          x2 = (1.-zgam*zal*zrgm/zza)**.25
          y2 = (1.-zgam*zal*zrgt/zza)**.25
          zpsim = psium(x1,x2)
          zpsih = psiuh(x1,y2)
          enddo
        z2 = (1.-zgam*zal*zrgq/zza)**.25
        zpsiq   = psiuh(x1, z2)
        x1 = (1.-zgam*zal*(10.+zrgm)/zza)**.25
        zpsim10 = psium(x1, x2)
        x1 = (1.-zgam*zal*(2. +zrgt)/zza)**.25
        zpsih2  = psiuh(x1, y2)
        x1 = (1.-zgam*zal*(2. +zrgq)/zza)**.25
        zpsiq2  = psiuh(x1, z2)
        x1 = (1.-zgam*zal*(.5 +zrgt)/zza)**.25
        zpsih05 = psiuh(x1, y2)
        x1 = (1.-zgam*zal*(.05+zrgt)/zza)**.25
        zpsih005= psiuh(x1, y2)
        x1 = (1.-zgam*zal*(.5 +zrgq)/zza)**.25
        zpsiq05 = psiuh(x1, z2)
      endif

! Turbulent fluxes of momentum, temperature, and specific humidity

      zustar = ak/(zalzam-zpsim)*zmod
      ztstar = ak/(zalzat-zpsih)*zdte
      zqstar = ak/(zalzaq-zpsiq)*zdq

      zv10  =                    zustar/ak*(log((10.  +zrgm)/zrgm)-zpsim10 )
      zt2   = ztetg            + ztstar/ak*(log((2.   +zrgt)/zrgt)-zpsih2  )
      zt05  = ztetg            + ztstar/ak*(log((0.5  +zrgt)/zrgt)-zpsih05 )
      zt005 = ztetg            + ztstar/ak*(log((0.05 +zrgt)/zrgt)-zpsih005)
      zq2   = qskin(jlon,jlat) + zqstar/ak*(log((2.   +zrgq)/zrgq)-zpsiq2  )
      zq05  = qskin(jlon,jlat) + zqstar/ak*(log((0.5  +zrgq)/zrgq)-zpsiq05 )

! Wind components

      u10(jlon,jlat) = zua/zmod*zv10
      v10(jlon,jlat) = zva/zmod*zv10

! From potential temperature to temperature

      t2  (jlon,jlat) = zt2  *zconv
      t05 (jlon,jlat) = zt05 *zconv
      t005(jlon,jlat) = zt005*zconv

! Specific humidity near surface must not exceed saturation with respect to water and ice
! Relative humidity computed with respect to water (WMO standard)
! Dew point temperature

      call comp_esk(zeskl,zqsat,t2(jlon,jlat),ps(jlon,jlat),2) ! Computes blended saturation below freezing
      q2(jlon,jlat) = min(zq2,zqsat)

      call comp_esk(zeskl,zqsat,t2(jlon,jlat),ps(jlon,jlat),3) ! Computes saturation to water also below freezing
      q2p=min(zq2,zqsat*1.01) ! assures that Q2REL does not exceed 101% (1% more for smoothing graphics)
      eee=ps(jlon,jlat)*q2p/(0.622*(1.-q2p)+q2p)
      q2rel(jlon,jlat) = eee/zeskl*100.

      call comp_esk(zeskl,zqsat,t05(jlon,jlat),ps(jlon,jlat),3) ! Computes saturation to water also below freezing
      q05p=min(zq05,zqsat*1.01) ! assures that Q2REL does not exceed 101% (1% more for smoothing graphics)
      eee=ps(jlon,jlat)*q05p/(0.622*(1.-q05p)+q05p)
      q05rel(jlon,jlat) = eee/zeskl*100.

      call td_tetens(t2(jlon,jlat), ps(jlon,jlat), q2(jlon,jlat), eps, td2(jlon,jlat))

    enddo

  enddo

!      call divdam (t2, nlon, nlat, 1, zwork) ! filter
!      call divdam (q2rel, nlon, nlat, 1, zwork) ! filter

      return
      end
!#######################################################################
    subroutine vtsurflux

! Computes surface fluxes of sensible and latent heat consistently with BOLAM a.b.l.

    use mod_postbolam

    real, parameter :: zgam=16., zaholt=1., zbholt=.667, zcholt=5., zdholt=.35

! Businger functions

    psium(zz1,zz2) = log( (1.+zz1)**2*(1.+zz1**2)/((1.+zz2)**2*(1.+zz2**2)) ) -2.*(atan(zz1)-atan(zz2))
    psiuh(zz1,zz2) = 2.*log((1.+zz1**2)/(1.+zz2**2))

! Holtslag functions

    psism(zz1,zz2)=-zaholt*zz1-zbholt*(zz1-zcholt/zdholt)*exp(-zdholt*zz1)  &
                   +zaholt*zz2+zbholt*(zz2-zcholt/zdholt)*exp(-zdholt*zz2)
    psish(zz1,zz2)=-(1.+2./3.*zaholt*zz1)**1.5-zbholt*(zz1-zcholt/zdholt)*exp(-zdholt*zz1) &
                   +(1.+2./3.*zaholt*zz2)**1.5+zbholt*(zz2-zcholt/zdholt)*exp(-zdholt*zz2)

! Latent heat

    healat(z) = 2.5008e6 - 2.36e3 * (z-ttr)

! Loop over horizontal grid

  do jlat = 1, nlat
    jlatp1 = min(jlat+1,nlat)

    do jlon = 1, nlon
      jlonm1= max(1,jlon-1)

! Definition of potential temp. corresponding to the ground temperature

      zconv = exp( rcp*(alps(jlon,jlat)-alp0) )
      ztetg = tskin(jlon,jlat)/zconv
      zthvg = ztetg*(1.+ep*qskin(jlon,jlat))

      zph   = phi(jlon,jlat,nlev)-phig(jlon,jlat)
      zua   = .5*(u(jlon,jlat  ,nlev)+u(jlonm1,jlat,nlev))
      zva   = .5*(v(jlon,jlatp1,nlev)+v(jlon  ,jlat,nlev))
      zdq   = q(jlon,jlat,nlev)-qskin(jlon,jlat)
      zdte  = t(jlon,jlat,nlev)-ztetg
      zzpp  = pzer*sigint(nlev)-(pzer-ps(jlon,jlat))*sigalf(nlev)
      zthvn = tvirt(jlon,jlat,nlev)*(1.e5/zzpp)**rcp
      zza   = zph/g + rgm(jlon,jlat)
      zmod2  = zua**2 + zva**2 + 0.07
      zmod  = sqrt(zmod2)
      zro   = ps(jlon,jlat)/rd/tvirt(jlon,jlat,nlev)

! Bulk Richardson number

      ztbarv = 0.5*(zthvg+zthvn)
      zri = g*zza*(zthvn-zthvg)/(ztbarv*zmod2)
      richs(jlon,jlat) = zri ! to be used for the definition of clouds

! Definition of roughness zrgm, zrgt, zrgq

        if (fmask(jlon,jlat).gt..5) then ! Computation of Charnok roughness
        zchar = 5.e-4
        zcoch1=.0185*ak**2/g*zmod2
        do jiter = 1, 5
        zcoch2= zza/zchar
        zchar = zcoch1/log(zcoch2)**2
        enddo
        zrgm = zchar

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

        else
        zsea = 5.e-5
        zrgm = fmask(jlon,jlat)*zsea+(1.-fmask(jlon,jlat))*rgm(jlon,jlat)
        zrgt = fmask(jlon,jlat)*zsea+(1.-fmask(jlon,jlat))*rgq(jlon,jlat)
        zrgq = zrgt
        endif

      zalzam = log(zza/zrgm)
      zalzat = log(zza/zrgt)
      zalzaq = log(zza/zrgq)

      if(zri.gt.0.) then         ! Holtslag functions
        zal = .2 + 4.*zri
        zpsim = psism(zal, zal*zrgm/zza)
        zpsih = psish(zal, zal*zrgt/zza)
          do jiter = 1, 10
          zalm = zal
          zal = zri*(zalzam-zpsim)**2/(zalzat-zpsih)
          error = abs(zal-zalm)/zal
          if (error.lt.1.e-2) go to 1
          zpsim = psism(zal, zal*zrgm/zza)
          zpsih = psish(zal, zal*zrgt/zza)
          enddo
 1      continue
        zpsiq   = zpsih       ! because rgt=rgq in the stable case
      else                    ! Businger functions
        zpsim = 0.
        zpsih = 0.
          do jiter = 1, 4
          zal = zri*(zalzam-zpsim)**2/(zalzat-zpsih)   ! za/l from eq (6.47) of Businger
          x1 = (1.-zgam*zal)**.25
          x2 = (1.-zgam*zal*zrgm/zza)**.25
          y2 = (1.-zgam*zal*zrgt/zza)**.25
          zpsim = psium(x1,x2)
          zpsih = psiuh(x1,y2)
          enddo
        z2 = (1.-zgam*zal*zrgq/zza)**.25
        zpsiq   = psiuh(x1, z2)
      endif

! Turbulent fluxes of momentum, temperature, and specific humidity

      zustar = ak/(zalzam-zpsim)*zmod
      ztstar = ak/(zalzat-zpsih)*zdte
      zqstar = ak/(zalzaq-zpsiq)*zdq

! Surface fluxes of sensible and latent heat (positive upward)

      hflux(jlon,jlat) =-zro*cp*zustar*ztstar*zconv
      qflux(jlon,jlat) =-zro*healat(tskin(jlon,jlat))*zustar*zqstar

    enddo

  enddo

      return
      end
!#######################################################################
      subroutine ccloud

! Computes total and high, middle and low cloud fraction as a function of
! cloud condensate and relative humidity.
! Reduction of cloud fraction computed as a function of the Richardson number,
! depending on moist static stability as Durran & Klemp (1982).
! Cloud cover algorithm as revised from Geleyn's (Maurizio).
! Low, middle, high clouds as WMO definition: < 2000 m, 2000-6000 m, > 6000 m.

 use mod_postbolam

! ---- For RTTOV simulation ---

#ifdef rttov
 use mod_rttov, only : clfrac
#endif

! -----------------------------
    real zqs(nlon,nlev), ztetav(nlon,nlev), fcloud(nlon,nlev), rich(nlon,nlevp1)
    real zstabg(nlon), zindn(nlon), zcol(nlev)

    ntop = 4
    yliv = 2834170.5
    yliw = 333560.5
    ylwv = yliv - yliw
    ezer = 611.
    cpv  = 1869.46
    cw   = 4186.8
    ccw1 = (cpv-cw)/rv
    ccw2 = ylwv/ttr/rv-ccw1
    zomd = cpv/cp-1.

    qccrit = 2.2e-4

! Definition of critical rel. hum. for cloud formation depending on level
! (do not set huc = 1. to avoid overflow)

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

    do jklev = 1, nlev
     if(sigint(jklev).lt.slow.and.sigint(jklev).gt.smid) then
     huc(jklev) = hmid
     elseif(sigint(jklev).ge.slow) then
     huc(jklev) = zalf * sigint(jklev) + zbet
     elseif(sigint(jklev).le.smid.and.sigint(jklev).gt.shig) then
     huc(jklev) = zalf2 * sigint(jklev) + zbet2
     else
     huc(jklev) = hhig
     endif
    enddo

! Computation of dry and moist static stability and of effective static stability
! as a function of relative humidity.
! Computation of the Richardson number.

    do jlat = 1, nlat

! Comput. of virtual theta

    do jklev = ntop-1, nlev
    do jlon = 1, nlon
     zzpp = pzer*sigint(jklev) - (pzer-ps(jlon,jlat))*sigalf(jklev)
     zconv = exp(rcp*(alp0-log(zzpp)))
     zt0t = ttr/temp(jlon,jlat,jklev)
     zesk = ezer*exp(-ccw1*log(zt0t)+ccw2*(1.-zt0t))     ! partial pressure over water
     zqs(jlon,jklev) = zesk*eps/(zzpp+zesk*(eps-1.))
     ztetav(jlon,jklev) = tvirt(jlon,jlat,jklev)*zconv
    enddo
    enddo

!  Computation of the Richardson no. depending on dry and moist (saturated) static stability
!  as Durran & Klemp (1982). (Specific humidity is used in place of mixing ratio).

    do jklev = ntop, nlev   ! loop on half levels
    do jlon = 1, nlon
     dthd   = 2.*(ztetav(jlon,jklev-1)-ztetav(jlon,jklev))/(ztetav(jlon,jklev-1)+ztetav(jlon,jklev)) ! dry
     r_up = zqs(jlon,jklev-1)
     r_do = zqs(jlon,jklev  )
     r_av = 0.5*(r_up + r_do)
     t_av = 0.5*(temp(jlon,jlat,jklev) + temp(jlon,jlat,jklev-1))
     theta_up = t(jlon,jlat,jklev-1)
     theta_do = t(jlon,jlat,jklev  )
     zaa = (1. + ylwv*r_av/(rd*t_av))/(1. + eps*ylwv**2*r_av/(cp*rd*t_av**2))
     dthm = zaa*((theta_up-theta_do)*2./(theta_up + theta_do) + ylwv/(cp*t_av)*(r_up - r_do)) & ! Durran&Klemp
            - q(jlon,jlat,jklev-1) - qc(jlon,jlat,jklev-1)                                    &
            + q(jlon,jlat,jklev  ) + qc(jlon,jlat,jklev  )

! Average relative humidity computed giving some more weight to the layer below

     zrhm = 0.55*q(jlon,jlat,jklev)/zqs(jlon,jklev) + 0.45*q(jlon,jlat,jklev-1)/zqs(jlon,jklev-1)
     zcof = max(-24.5 + 25.*zrhm, 0.)    ! zcof=0. for rh<=0.98, zcof=0.5 for rh=1
     zcof = min(zcof, .85)
     dlogthe = zcof*dthm + (1.-zcof)*dthd ! effective stability near saturation
     zrdz  = g/(phi(jlon,jlat,jklev-1)-phi(jlon,jlat,jklev))
     zdtdz = dlogthe*zrdz

     jlonm1 = max(jlon-1,1)
     jlatp1 = min(jlat+1,nlat)
     zdudz = .5*(u(jlon,jlat,jklev-1)+u(jlonm1,jlat,jklev-1)-u(jlon,jlat,jklev)-u(jlonm1,jlat,jklev))*zrdz
     zdvdz = .5*(v(jlon,jlat,jklev-1)+v(jlon,jlatp1,jklev-1)-v(jlon,jlat,jklev)-v(jlon,jlatp1,jklev))*zrdz
     zshear = zdudz**2 + zdvdz**2
     zbuoy = g*dlogthe*zrdz
     rich(jlon,jklev) = min (zbuoy/(zshear + 1.e-6), 500.)
    enddo
    enddo
    rich(:,nlev+1) = richs(:,jlat)  ! Richardson no. at the surface (computed in vtsurflux)

    do jklev = ntop, nlev
    do jlon = 1, nlon
     zppp = pzer*sigint(jklev)-(pzer-ps(jlon,jlat))*sigalf(jklev)
     ztmp = temp(jlon,jlat,jklev)
     call comp_esk(zesat,zqs1,ztmp,zppp,2)  ! blended saturation
     call comp_esk(zesat,zqs2,ztmp,zppp,3)  ! sat. to water below 0
     zzqs = 0.70*zqs1 + 0.30*zqs2           ! ad hoc to limit contrib. to high clouds
     zrh1 = min(1., q(jlon,jlat,jklev)/zzqs)
     zrh = (zrh1-huc(jklev))/(1.-huc(jklev))
     zrh = min(1., zrh)
     zrh = max(0., zrh)

! Approximate separation between water and ice clouds

     if(ztmp.ge.ttr) then
      zratio = 1.
     elseif (ztmp.lt.ttr-28.) then
      zratio = 0.
     else
      zratio = 1.04979*(0.5 + 0.5*tanh((ztmp-ttr+9.)/6.))
     endif
     zqcw = zratio*qc(jlon,jlat,jklev)
     zqci = qc(jlon,jlat,jklev) - zqcw
     fcloud(jlon,jklev) = (max(0.15*zrh, (1.*zqcw + 0.77*zqci)/qccrit))**.7
    enddo
    enddo

! Initialization for computing Geleyn stability:
! weighted average of static energy between the surface and the lowest model level.
! More weight is given to surface t and q (with respect to air) over sea than over land,
! to reduce clouds more over sea than over land in unstable pbl.

    do jlon = 1, nlon
     zwe = 0.1 + fmask(jlon,jlat)*0.6
     zstabg(jlon) = zwe*(cp*(1. + zomd*qskin(jlon,jlat))*tskin(jlon,jlat) + phig(jlon,jlat)) +           &
                   (1.-zwe)*(cp*(1. + zomd*q(jlon,jlat,nlev))*temp(jlon,jlat,nlev) + phi(jlon,jlat,nlev))
     zindn(jlon) = 1.
    enddo

! Reduction of clouds in unstable layers

    fcrit = 0.3
    do jklev = nlev, ntop, -1
    do jlon = 1, nlon

! Geleyn stability, based on (dry, virtual) static energy profile (non-local, for pbl only)

     zdstabg = zstabg(jlon) - cp*(1. + zomd*q(jlon,jlat,jklev))*temp(jlon,jlat,jklev) - phi(jlon,jlat,jklev)
     zindn(jlon) = zindn(jlon)*max(0., sign(1., zdstabg))

! Local stability derived from the Richardson no.

     zstabr = 0.5*rich(jlon,jklev+1) + 0.5*rich(jlon,jklev)

     fc1 = 1.
     fc2 = 1.
     if(zstabr.lt..25.and.zstabr.gt.0.) fc1 = 0.875 + 0.5*zstabr
     if(zstabr.gt.0..and.rich(jlon,jklev+1).lt.0.) fc1 = 0.8
     if(zstabr.le.0.) fc1 = 0.7

! Over the sea: computation of cloud depth and cloud top to distinguish between cumuli
! and stratocumuli, the latter to be more strongly reduced by Geleyn method

     if(fmask(jlon,jlat).gt.0.5) then ! sea
     cloudd = 0.
     ktop = nlev
     ifl = 0
     if(zindn(jlon).gt..9) then
      do jk = jklev-1, nlev/2, -1
      if(fcloud(jlon,jk).gt.fcrit) then
      cloudd = cloudd + 0.5*(phi(jlon,jlat,jk-1)-phi(jlon,jlat,jk+1))/g
      ktop = jk
      ifl = 1
      else
      exit
      endif
      enddo
      if(ifl.eq.1) then
       ritop = rich(jlon,ktop)  ! Ri at half lev. above cloud top
       cltop = (phi(jlon,jlat,ktop)-phig(jlon,jlat))/g  ! height of cloud top
       if(cloudd.gt.480..or.cltop.gt.1550..or.ritop.lt.0.5) zindn(jlon) = 0.55 ! cumuli, not stratocumuli
      endif
     endif
     zrc = 0.40
     else  ! land
     zrc = 0.22
     endif !  sea or land

     fc2 = 1. - zrc*zindn(jlon)
     fcloud(jlon,jklev) = min(fc1, fc2)*fcloud(jlon,jklev)
     fcloud(jlon,jklev) = max(0., min(1., fcloud(jlon,jklev)))
    enddo
    enddo

    do jlon = 1, nlon
     nlclo = 0
     nmclo = 0

     do jklev = 1, nlev
      if(phi(jlon,jlat,jklev)/g.gt.2000.) nlclo = nlclo + 1 ! no. of levels above 2000 m
      if(phi(jlon,jlat,jklev)/g.gt.6000.) nmclo = nmclo + 1 ! no. of levels above 6000 m
     enddo

! High clouds

     n1 = ntop
     n2 = nmclo
     ntot = n2 - n1 + 1
     zcol(1:ntot) = fcloud(jlon, n1:n2)
     call ccolumn (zcol, ntot, zprod)
     cloudh(jlon,jlat) = max(0., 1.-zprod)

! Middle clouds

     n1 = nmclo + 1
     n2 = nlclo
     ntot = n2 - n1 + 1
     if(ntot.ge.1) then
     zcol(1:ntot) = fcloud(jlon, n1:n2)
     call ccolumn (zcol, ntot, zprod)
     cloudm(jlon,jlat) = max(0., 1.-zprod)
     endif

! Low clouds

     n1 = nlclo + 1
     n2 = nlev
     ntot = n2 - n1 + 1
     if(ntot.ge.1) then
     zcol(1:ntot) = fcloud(jlon, n1:n2)
     call ccolumn (zcol, ntot, zprod)
     cloudl(jlon,jlat) = max(0., 1.-zprod)
     endif

! Total clouds

     n1 = ntop
     n2 = nlev
     ntot = n2 - n1 + 1
     if(ntot.ge.1) then
     zcol(1:ntot) = fcloud(jlon, n1:n2)
     call ccolumn (zcol, ntot, zprod)
     cloud(jlon,jlat) = max(0., 1.-zprod)
     endif

! ---- For RTTOV simulation ---

#ifdef rttov
     clfrac(jlon,jlat,:) = 0.
     do jklev=1,nlev
       clfrac(jlon,jlat,jklev) = min(max(fcloud(jlon,jklev),0.),1.)
     enddo
#endif

    enddo ! jlon

! -----------------------------

! Compute height (over surface) of top of pbl

    do jlon = 1, nlon
    jkpbl = 0
     do jklev = nlev, nlev/2, -1
     zriav = .5*(rich(jlon,jklev+1)+rich(jlon,jklev))
     if (zriav > 0.25) exit
     jkpbl = jklev
     enddo
    if (jkpbl == 0.and.rich(jlon,nlev+1) < 0.25) jkpbl = nlev
    if (jkpbl > 0) pbl(jlon,jlat) = (phi(jlon,jlat,jkpbl)-phig(jlon,jlat))/g
    pbl(jlon,jlat) = min (pbl(jlon,jlat), 2500.)
    enddo

    enddo ! jlat

    return
    end
!#######################################################################
   subroutine ccolumn (a,n,prod)

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
!#######################################################################
      subroutine thetdef(klat)

! Computes pv on specific isentropic surfaces (by interpolation from pv computed on sigma surfaces)

 use mod_postbolam
 use mod_rtlev
 use mod_rwoutt
 use mod_ptvor

       dimension ztpv(nlev),ztpv2(nlev),zaux(nlev)

       do 600 jlon=1, nlon-1

!  Interpolation of pot. vort and interpolation of theta on zita (pv) points and on
!  semi-integer levels where pv is defined.
!  Note that zaux and ztpv are inverted in order because ztpv must be increasing for routine interp

       do jlev= 1, nlev
       zaux(nlev-jlev+1) = potv(jlon,klat,jlev)
       enddo

       do jlev= 1, nlev-1
       ztpv(nlev-jlev+1)= 0.5*(.25*t(jlon,klat,jlev)+.25*t(jlon+1,klat,jlev)+.25*t(jlon,klat-1,jlev)+     &
                               .25*t(jlon+1,klat-1,jlev) +.25*t(jlon,klat,jlev+1)+.25*t(jlon+1,klat,jlev+1) + &
                               .25*t(jlon,klat-1,jlev+1)+.25*t(jlon+1,klat-1,jlev+1))
       enddo

! Extrapolation to define ztpv at nlev+1/2 (surface) (provided not super-adiabatic)

       ztpv(1) = min(2.*ztpv(2) - ztpv(3), ztpv(2))

      do jlev= 1, nlev
      ztpv2(jlev) = ztpv(jlev)
      enddo

! Check and removal of superadiabatic layers to have monotonic theta
! starting from the top, then from the bottom and averaging

      do jlev= nlev,2,-1
      if (ztpv(jlev-1).ge.ztpv(jlev)) ztpv(jlev-1)=ztpv(jlev)-0.05
      enddo

      do jlev= 1, nlev-1
      if (ztpv2(jlev+1).le.ztpv2(jlev)) ztpv2(jlev+1)=ztpv2(jlev)+0.05
      enddo

      do jlev= 1, nlev
      ztpv(jlev) = 0.5*(ztpv(jlev)+ztpv2(jlev))
      enddo

! Interpolation of pv on theta

      zalf = .5
      ze1 = .7
      ze2 = .7
      call interp(zalf,ze1,ze2,nlev,ztpv,zaux,tlevo,zvpvt(jlon,1:nlevto),nlevto)

! Reset to zero of pv for theta surfaces below ground

      do jlev=1,nlevto
      if (tlevo(jlev).lt.ztpv(1)) zvpvt(jlon,jlev)=0.
      enddo

 600  continue

      return
      end
!#######################################################################
      subroutine lift_parcel_bol (qmix, tmix, p, tvlift, nlev, jk0, jk1, iwl)

! Computes virtual temperature of a moist air parcel lifted from level jk0
! to a generic level jklev such that: jk0 >= jklev >= jk1
! qmix and tmix are the parcel properties at level jk0
! Results are saved in tvlift
! Liquid water may be removed during lifting (case iwl=0)

      real, dimension(nlev) :: p, tvlift
      real, parameter :: yliv=2834170.5, yliw=333560.5, ylwv=yliv-yliw, tzer=273.15, pzer=1.e5, ezer=611.0, &
                         rd=287.05, rv=461.51, eps=rd/rv, cpd=1004.6, cpv=1869.46, cw=4186.8,               &
                         ccw1=(cpv-cw)/rv, ccw2=ylwv/tzer/rv-ccw1

      tvlift(jk0) =  tmix*(1.+(1./eps-1.)*qmix)

      jkconl = 1

! Calculation of parcel entropy before lifting (zs0)

      zesk  = qmix*p(jk0)/( eps+(1.-eps)*qmix )
      zctot = (1.-qmix)*cpd+qmix*cpv
      zsa   = -(1.-qmix)*rd*log((p(jk0)-zesk)/pzer)
      zsb   = -qmix*rv*log(zesk/ezer+1.e-18)
      zsc   =  qmix*yliv/tzer
      zs0   = zctot*log(tmix/tzer) + zsa + zsb + zsc

      do jklev = jk0-1, jk1, -1

! Lifted parcel temperature by conservation of entropy in the hypothesis of no condensation

      zesk  =  qmix*p(jklev)/( eps+(1.-eps)*qmix )
      zctot =  (1.-qmix)*cpd + qmix*cpv
      zsa   = -(1.-qmix)*rd*log((p(jklev)-zesk)/pzer)
      zsb   = -qmix*rv*log(zesk/ezer+1.e-18)
      zsc   =  qmix*yliv/tzer
      tlift =  tzer*exp((zs0-zsa-zsb-zsc)/zctot)

! Saturation with respect to water

      zt0t = tzer/tlift
      zesk = ezer*exp( -ccw1*log(zt0t) + ccw2*(1.-zt0t) )  ! saturated pressure over water (in pa)
      zqsw = zesk*eps/(p(jklev)+zesk*(eps-1.))

      if (qmix.le.zqsw) then   ! no condensation occurs
      tvlift(jklev) =  tlift*(1.+(1./eps-1.)*qmix)
      else
      jkconl = jklev    ! index of the model level just above the lifting condensation level
      go to 9
      endif
      enddo  ! close loop over unsaturated parcel
 9    continue

! Iterative solution in case of condensation with Newton steps

      zt1   = tvlift(jkconl+1)
      zq0   = qmix

      do jklev = jkconl, jk1, -1

! Entropy value at zt1 to start the Newton step (zs1)

      zt0t  =  tzer/zt1
      zesk  =  ezer*exp( -ccw1*log(zt0t) + ccw2*(1.-zt0t) )  ! saturated pressure over water (in Pa)
      zqsw  =  zesk*eps/(p(jklev)+zesk*(eps-1.))
      zwat  =  zq0 - zqsw
      zctot =  (1.-zq0)*cpd + zqsw*cpv + zwat*cw
      zsa   = -(1.-zq0)*rd*log((p(jklev)-zesk)/pzer)
      zsb   = -zqsw*rv*log(zesk/ezer)
      zsc   =  (zqsw*yliv + zwat*yliw)/tzer
      zs1   = -zctot*log(zt0t) + zsa + zsb + zsc

      tlift =  zt1 - 1.   ! first guess of lifted parcel temperature
      icount = 0

 10   continue

! Entropy of saturated air with condensed water/ice

      zt0t  =  tzer/tlift
      zesk  =  ezer*exp( -ccw1*log(zt0t) + ccw2*(1.-zt0t) )   ! saturated pressure over water (in pa)
      zqsw  =  zesk*eps/(p(jklev)+zesk*(eps-1.))
      zwat  =  zq0 - zqsw
      zctot =  (1.-zq0)*cpd + zqsw*cpv + zwat*cw
      zsa   = -(1.-zq0)*rd*log((p(jklev)-zesk)/pzer)
      zsb   = -zqsw*rv*log(zesk/ezer)
      zsc   =  (zqsw*yliv + zwat*yliw)/tzer
      zs    = -zctot*log(zt0t) + zsa + zsb + zsc

! Newton step

      if (abs(tlift-zt1).gt.0.1 .and. icount.le.5) then

      zrds  = (zt1-tlift)/(zs1-zs)
      zt1 = tlift
      zs1 = zs
      tlift = tlift + (zs0-zs)*zrds
      icount = icount + 1
      go to 10

      else

      if (iwl.eq.0) then
      tvlift(jklev) = tlift*(1.+(1./eps-1.)*zqsw)        ! virtual temperature with no water loading
      elseif (iwl.eq.1) then
      tvlift(jklev) = tlift*(1.+(1./eps-1.)*zqsw - zwat) ! virtual temperature with water loading
      endif

      endif

! Removal of liquid water from the lifted parcel

      if (iwl.eq.0) then
      zq0 = zqsw          ! new total water content of lifted parcel
      zs0 = zs0 - zwat*cw*log(tlift/tzer) - zwat*yliw/tzer
      endif

      enddo   ! close loop over model levels where parcel is saturated

      return
      end
!#######################################################################
      subroutine wafps(ps,u,v,zdiv,nlon,nlat,nlev,nlonp1,hxt,hxv,dx,dy,dt,waflux)

!  Divergence on t points of the 2-d mass flux with WAF scheme

      real u(nlon,nlat,nlev),v(nlon,nlat,nlev)
      real ps(1:nlonp1,0:nlat),hxt(nlat),hxv(nlat)
      real waflux(nlonp1,nlat),zdiv(nlon,nlat,nlev)
      real denrx(nlonp1,nlat), denry(nlonp1,nlat)

      rdt=1./dt
      do jlat=1,nlat-1
      do jlon=2,nlonp1
      denr=ps(jlon,jlat)-ps(jlon-1,jlat)
        if(abs(denr).lt.1.e-15) then
        if(denr.ge.0.) then
        denr=1.e-15
        else
        denr=-1.e-15
        endif
        endif
      denrx(jlon,jlat)=denr
      enddo
      enddo

      do jlat=1,nlat
      do jlon=2,nlon
      denr=ps(jlon,jlat)-ps(jlon,jlat-1)
        if(abs(denr).lt.1.e-15) then
        if(denr.ge.0.) then
        denr=1.e-15
        else
        denr=-1.e-15
        endif
        endif
      denry(jlon,jlat)=denr
      enddo
      enddo

      do 1000 jklev=1,nlev

      do jlat=1,nlat-1
      zcost=dt/(dx*hxt(jlat))
      do jlon=2,nlonp1
      amu=u(jlon-1,jlat,jklev)*zcost
        if(amu.ge.0.) then
        is=1
        j1=jlon-1
        j1m1=j1-1
        if(j1m1.lt.1) j1m1=1
        else
        is=-1
        j1=jlon+1
        if(j1.gt.nlonp1) j1=nlonp1
        j1m1=j1-1
        endif
       r=(ps(j1,jlat)-ps(j1m1,jlat))/denrx(jlon,jlat)

       if(r.gt.0.) then
       if(r.gt.2.) then
       b=2.
       elseif(r.gt.1.) then
       b=r
       elseif(r.gt.0.5) then
       b=1.
       else
       b=2.*r
       endif
       else
       b=0.
       endif

      phi=is+amu*b-is*b
      waflux(jlon,jlat)=amu*((1.+phi)*ps(jlon-1,jlat)+(1.-phi)*ps(jlon  ,jlat))
      enddo
      do jlon=2,nlon
      zdiv(jlon,jlat,jklev)=.5*waflux(jlon+1,jlat)-.5*waflux(jlon,jlat)
      enddo
      enddo

      zcost=dt/dy
      do jlat=1,nlat
      do jlon=2,nlon
      amu=v(jlon,jlat,jklev)*zcost
        if(amu.ge.0.) then
        is=1
        j1=jlat-1
        j1m1=j1-1
        if(j1m1.lt.0) j1m1=0
        else
        is=-1
        j1=jlat+1
        if(j1.gt.nlat) j1=nlat
        j1m1=j1-1
        endif
       r=(ps(jlon,j1)-ps(jlon,j1m1))/denry(jlon,jlat)

       if(r.gt.0.) then
       if(r.gt.2.) then
       b=2.
       elseif(r.gt.1.) then
       b=r
       elseif(r.gt.0.5) then
       b=1.
       else
       b=2.*r
       endif
       else
       b=0.
       endif

      phi=is+amu*b-is*b
      waflux(jlon,jlat)=0.5*amu*((1.+phi)*ps(jlon,jlat-1)+(1.-phi)*ps(jlon,jlat))
      enddo
      enddo
      do jlat=1,nlat-1
      rhxt=1./hxt(jlat)
      areap=hxv(jlat+1)*rhxt
      aream=hxv(jlat  )*rhxt
      do jlon=2,nlon
      zdiv(jlon,jlat,jklev)=(zdiv(jlon,jlat,jklev)+waflux(jlon,jlat+1)*areap-waflux(jlon,jlat)*aream)*rdt
      enddo
      enddo

 1000 continue

      return
      end
!######################################################################
      subroutine divdam (div, nlon, nlat, nlev, p2)

!  2-grid interval wave filter over the whole domain

      real div(nlon,nlat,nlev), p2(nlon,nlat)

      do jklev = 1, nlev
      do jlat = 2, nlat-1
      do jlon = 2, nlon-1
      p2(jlon,jlat)=.25*(div(jlon-1,jlat,jklev)+div(jlon+1,jlat,jklev)) +.5*div(jlon,jlat,jklev)
      enddo
      enddo

      do jlat = 3, nlat-2
      do jlon = 3, nlon-2
      div(jlon,jlat,jklev)=.25*(p2(jlon,jlat+1)+p2(jlon,jlat-1)) +.5*p2(jlon,jlat)
      enddo
      enddo
      enddo

      return
      end
!#######################################################################
subroutine free_cross_section

use mod_cross
use mod_rwcross
use mod_postbolam, only : nlon, nlat, nlev, pdr, sigint, sigalf, sig, alfa, pzer, eps, pi, &
 alont, alatt, alonu, alatu, alonv, alatv, alonz, alatz,                                   &
 ps, t, temp, q, u, v, omeg, valmiss
use mod_ptvor
implicit none

real, dimension(npoint_cross_max,cross_number) :: lon_cross_rot, lat_cross_rot
real, dimension(nlon,nlat,nlev) :: theta_loc, thetae_loc, w_loc, rh_loc
real, dimension(npoint_cross_max,nlev,cross_number) :: theta_cross_lev, thetae_cross_lev, &
 u_cross_lev, v_cross_lev, w_cross_lev, pv_cross_lev, rh_cross_lev
real, dimension(nlon) :: lon_rot, lon_u_rot
real, dimension(nlat) :: lat_rot, lat_v_rot
real, dimension(npress) :: lnp_out, work_1d
real, dimension(nlev) :: lnp_inp, lnp_inp_1
integer, dimension(npress) :: iv
integer, dimension(nlon,nlat) :: mask_field
real, dimension(nlon,nlat) :: work2d
integer, parameter :: npar=7
real, dimension(nlon,nlat,nlev,npar) :: work4d
real, dimension(npoint_cross_max,nlev,npar) :: work3d

integer :: icross, i, jlon, jlat, jklev, jklev1, np, nlev_inp, im1, ip1, nsmooth=5
real :: zdlon, zdlat, ztmpr, zpres, zqq, ze, zmr, zrh, zzqs, zzpvs, zzpv, ztl, zespon, &
 x0, y0, dlon, dlat, alon0, alat0, semi_radius=20., gauss_par=1., wei_smooth=0.5,      &
 alfa_hinterp=0.6, alfa_vinterp=0.5, ex_bott=0.6, ex_top=0.6,                          &
 angle, zcos, zsin, zsgalf
character(len=20) :: file_out_work

real :: zss = 6.5e-03 ! Average of computed stability with climatological value 6.5 K/km

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

    theta_loc(jlon,jlat,jklev) = t(jlon,jlat,jklev)
    ztmpr = temp(jlon,jlat,jklev)

!  Quantities to compute rel. hum. and equiv. pot. temp.
!  zpres: air pressure; ze: vapour pressure; zmr: mix. ratio
!  ztl: temp. at condens. lev. (computed with expr. by Bolton)

    zpres = pzer*sigint(jklev)-(pzer-ps(jlon,jlat))*sigalf(jklev)
    zqq = q(jlon,jlat,jklev)
    zmr = zqq/(1.-zqq)
    ze = zpres/100. * zmr/(eps+zmr)
    if (ze < 1.e-13) ze = 1.e-13

! QTORH does not limit ZRH to saturation because it is used in cross-sections to extract
! cloud water+ice

    call qtorh(zqq, ztmpr, zpres, zrh)
    rh_loc(jlon,jlat,jklev)= zrh

!  In the following the vert. vel. is comp. as -omega/p

    w_loc(jlon,jlat,jklev) = -omeg(jlon,jlat,jklev)/zpres

! Limitation of mixing ratio to subsaturation for comput. of equivalent pot. temp.
! Saturation is computed using the blended curve below freezing

    call comp_esk(zzpvs, zzqs, ztmpr, zpres, 2)

    zzpv = (zqq*zpres)/(0.622+zqq-0.622*zqq)
    if (zzpv > zzpvs) then
      zzpv = zzpvs
      zmr = 0.622*zzpv/(zpres-zzpv)
      ze = zpres/100.*zmr/(eps+zmr)
    endif

!  Definition of equivalent pot. temp. for cross-sections and for const. pressure maps.
!  The definition follows Bolton (MWR,1980).
!  ztl: temp. at condens. lev. (computed with expr. by Bolton).

    ztl = 55. + 2840./(3.5*alog(ztmpr) - alog(ze) - 4.805)
    zespon = (3.376/ztl - 0.00254)*zmr*1000*(1.+.81*zmr)
    thetae_loc(jlon,jlat,jklev) = theta_loc(jlon,jlat,jklev)*exp(zespon)

  enddo
  enddo
  enddo

! Horizontal interpolation to cross-section points

! Interpolation with spline algorithm

  x0 = pdr(39)
  y0 = pdr(38)
  dlon = pdr(2)
  dlat = pdr(1)
  alon0 = pdr(5)
  alat0 = pdr(4)+dlat*0.5

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
      print *,'in rotated coordinates: max longitude of cross-section: ', &
 maxval(lon_cross_rot(1:np,icross)),' of domain: ',maxval(lon_rot)
      stop
    endif

    if (minval(lon_cross_rot(1:np,icross)) < minval(lon_rot)) then
      print*
      print *,'Vertical cross-section coordinate is outside if model domain:'
      print *,'Cross-section number ',icross
      print *,'in rotated coordinates: min longitude of cross-section: ', &
 minval(lon_cross_rot(1:np,icross)),' of domain: ',minval(lon_rot)
      stop
    endif

    if (maxval(lat_cross_rot(1:np,icross)) > maxval(lat_rot)) then
      print*
      print *,'Vertical cross-section coordinate is outside if model domain:'
      print *,'Cross-section number ',icross
      print *,'in rotated coordinates: max latitude of cross-section: ', &
 maxval(lat_cross_rot(1:np,icross)),' of domain: ',maxval(lat_rot)
      stop
    endif

    if (minval(lat_cross_rot(1:np,icross)) < minval(lat_rot)) then
      print*
      print *,'Vertical cross-section coordinate is outside if model domain:'
      print *,'Cross-section number ',icross
      print *,'in rotated coordinates: min latitude of cross-section: ', &
 minval(lat_cross_rot(1:np,icross)),' of domain: ',minval(lat_rot)
      stop
    endif

  enddo

  print*
  if (cross_interp_gauss == 1 ) then

    print *,'Horizontal interpolation for free cross-sections with Gaussian algorithm,'
    print *,'semi-radius=',semi_radius,' km'

! Interpolation with gaussian 2-d function

    work4d(1:nlon,1:nlat,1,1) = ps (1:nlon,1:nlat)
    do icross=1,cross_number
      np = npoint_cross(icross)
      call interp_gauss_input_grid (work4d(1:nlon,1:nlat,1:1,1:1), work3d(1:np,1:1,1:1), &
 alont, alatt, lon_cross(1:np,icross), lat_cross(1:np,icross),                           &
 nlon, nlat, np, 1, 1, mask_field, semi_radius, gauss_par, valmiss)
      ps_cross(1:np,icross) = work3d(1:np,1,1)
    enddo

    work4d(1:nlon,1:nlat,1:nlev,1) = theta_loc  (1:nlon,1:nlat,1:nlev)
    work4d(1:nlon,1:nlat,1:nlev,2) = thetae_loc (1:nlon,1:nlat,1:nlev)
    work4d(1:nlon,1:nlat,1:nlev,3) = w_loc      (1:nlon,1:nlat,1:nlev)
    work4d(1:nlon,1:nlat,1:nlev,4) = rh_loc     (1:nlon,1:nlat,1:nlev)
    do icross=1,cross_number
      np = npoint_cross(icross)
      call interp_gauss_input_grid (work4d(1:nlon,1:nlat,1:nlev,1:4), work3d(1:np,1:nlev,1:4), &
 alont, alatt, lon_cross(1:np,icross), lat_cross(1:np,icross),                                 &
 nlon, nlat, np, nlev, 4, mask_field, semi_radius, gauss_par, valmiss)
      theta_cross_lev(1:np,1:nlev,icross)  = work3d(1:np,1:nlev,1)
      thetae_cross_lev(1:np,1:nlev,icross) = work3d(1:np,1:nlev,2)
      w_cross_lev(1:np,1:nlev,icross)      = work3d(1:np,1:nlev,3)
      rh_cross_lev(1:np,1:nlev,icross)     = work3d(1:np,1:nlev,4)
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

    work4d(1:nlon,1:nlat,1:nlev,1) = potv(1:nlon,1:nlat,1:nlev)
    do icross=1,cross_number
      np = npoint_cross(icross)
      call interp_gauss_input_grid (work4d(1:nlon,1:nlat,1:nlev,1:1), work3d(1:np,1:nlev,1:1), &
 alonz, alatz, lon_cross(1:np,icross), lat_cross(1:np,icross),                                 &
 nlon, nlat, np, nlev, 1, mask_field, semi_radius, gauss_par, valmiss)
      pv_cross_lev(1:np,1:nlev,icross) = work3d(1:np,1:nlev,1)
    enddo

  else

! Interpolation with spline algorithm

    print *,'Horizontal interpolation for free cross-sections with spline algorithm and smoothed input fields,'
    print *,'parameters of spatial smoothing: weight ',wei_smooth,' numer ',nsmooth

! All 2d field is smoothed before horizontal interpolation

    work2d=ps(1:nlon,1:nlat)
    call smooth(ps(1:nlon,1:nlat), work2d, nlon, nlat, wei_smooth, nsmooth)
    do icross=1,cross_number
      np = npoint_cross(icross)
      call interp_spline_2d(work2d, nlon, nlat, lon_rot, lat_rot, &
 lon_cross_rot(1:np,icross), lat_cross_rot(1:np,icross),          &
 np, ps_cross(1:np,icross), 1.)
    enddo

    do jklev = 1,nlev

      work2d=theta_loc(:,:,jklev)
      call smooth(theta_loc(:,:,jklev), work2d, nlon, nlat, wei_smooth, nsmooth)
      do icross=1,cross_number
        np = npoint_cross(icross)
        call interp_spline_2d(work2d, nlon, nlat, lon_rot, lat_rot, &
 lon_cross_rot(1:np,icross), lat_cross_rot(1:np,icross),            &
 np, theta_cross_lev(1:np,jklev,icross), alfa_hinterp)
      enddo

      work2d=thetae_loc(:,:,jklev)
      call smooth(thetae_loc(:,:,jklev), work2d, nlon, nlat, wei_smooth, nsmooth)
      do icross=1,cross_number
        np = npoint_cross(icross)
        call interp_spline_2d(work2d, nlon, nlat, lon_rot, lat_rot, &
 lon_cross_rot(1:np,icross), lat_cross_rot(1:np,icross),            &
 np, thetae_cross_lev(1:np,jklev,icross), alfa_hinterp)
      enddo

      work2d=u(:,:,jklev)
      call smooth(u(:,:,jklev), work2d, nlon, nlat, wei_smooth, nsmooth)
      do icross=1,cross_number
        call interp_spline_2d(work2d, nlon, nlat, lon_u_rot, lat_rot, &
 lon_cross_rot(1:np,icross), lat_cross_rot(1:np,icross),              &
 np, u_cross_lev(1:np,jklev,icross), alfa_hinterp)
      enddo

      work2d=v(:,:,jklev)
      call smooth(v(:,:,jklev), work2d, nlon, nlat, wei_smooth, nsmooth)
      do icross=1,cross_number
        call interp_spline_2d(work2d, nlon, nlat, lon_rot, lat_v_rot, &
 lon_cross_rot(1:np,icross), lat_cross_rot(1:np,icross),              &
 np, v_cross_lev(1:np,jklev,icross), alfa_hinterp)
      enddo

      work2d=w_loc(:,:,jklev)
      call smooth(w_loc(:,:,jklev), work2d, nlon, nlat, wei_smooth, nsmooth)
      do icross=1,cross_number
        call interp_spline_2d(work2d, nlon, nlat, lon_rot, lat_rot, &
 lon_cross_rot(1:np,icross), lat_cross_rot(1:np,icross),            &
 np, w_cross_lev(1:np,jklev,icross), alfa_hinterp)
      enddo

      work2d=potv(:,:,jklev)
      call smooth(potv(:,:,jklev), work2d, nlon, nlat, wei_smooth, nsmooth)
      do icross=1,cross_number
        call interp_spline_2d(work2d, nlon, nlat, lon_u_rot, lat_v_rot, &
 lon_cross_rot(1:np,icross), lat_cross_rot(1:np,icross),                &
 np, pv_cross_lev(1:np,jklev,icross), alfa_hinterp)
      enddo

      work2d=rh_loc(:,:,jklev)
      call smooth(rh_loc(:,:,jklev), work2d, nlon, nlat, wei_smooth, nsmooth)
      do icross=1,cross_number
        call interp_spline_2d(work2d, nlon, nlat, lon_rot, lat_rot, &
 lon_cross_rot(1:np,icross), lat_cross_rot(1:np,icross),            &
 np, rh_cross_lev(1:np,jklev,icross), alfa_hinterp)
      enddo

    enddo

 endif  ! Horizontal iterpolation type
 print*

! Vertical interpolation to cross-section p-levels

  do jklev=1,npress
    jklev1=npress-jklev+1
    lnp_out(jklev1)=zplcr(jklev)
  enddo

  do icross=1,cross_number

    np = npoint_cross(icross)

    do i=1,np

! Definition of vertical axis in input: log(p)

      do jklev=1,nlev
        zsgalf = sig(jklev)**3*(1.+alfa*(1.-sig(jklev))*(2.-sig(jklev)))
        lnp_inp  (jklev) = alog (pzer*sigint(jklev)-(pzer-ps_cross(i,icross))*sigalf(jklev))
        lnp_inp_1(jklev) = alog (pzer*sig   (jklev)-(pzer-ps_cross(i,icross))*zsgalf       )
      enddo

      call near(lnp_out, npress, lnp_inp, nlev, iv)

      call interp_spline_1d(theta_cross(i,1:npress,icross), lnp_out, npress, &
 theta_cross_lev(i,1:nlev,icross), lnp_inp, nlev, iv, alfa_vinterp, ex_bott, ex_top)

      work_1d(1:npress)=theta_cross(i,1:npress,icross)
      do jklev=1,npress
        jklev1=npress-jklev+1
        theta_cross(i,jklev,icross)=work_1d(jklev1)
      enddo

      call interp_spline_1d(thetae_cross(i,1:npress,icross), lnp_out, npress, &
 thetae_cross_lev(i,1:nlev,icross), lnp_inp, nlev, iv, alfa_vinterp, ex_bott, ex_top)

      work_1d(1:npress)=thetae_cross(i,1:npress,icross)
      do jklev=1,npress
        jklev1=npress-jklev+1
        thetae_cross(i,jklev,icross)=work_1d(jklev1)
      enddo

      call interp_spline_1d(u_cross(i,1:npress,icross), lnp_out, npress, &
 u_cross_lev(i,1:nlev,icross), lnp_inp, nlev, iv, alfa_vinterp, ex_bott, ex_top)

      work_1d(1:npress)=u_cross(i,1:npress,icross)
      do jklev=1,npress
        jklev1=npress-jklev+1
        u_cross(i,jklev,icross)=work_1d(jklev1)
      enddo

      call interp_spline_1d(v_cross(i,1:npress,icross), lnp_out, npress, &
 v_cross_lev(i,1:nlev,icross), lnp_inp, nlev, iv, alfa_vinterp, ex_bott, ex_top)

      work_1d(1:npress)=v_cross(i,1:npress,icross)
      do jklev=1,npress
        jklev1=npress-jklev+1
        v_cross(i,jklev,icross)=work_1d(jklev1)
      enddo

      call interp_spline_1d(w_cross(i,1:npress,icross), lnp_out, npress, &
 w_cross_lev(i,1:nlev,icross), lnp_inp, nlev, iv, alfa_vinterp, ex_bott, ex_top)

      work_1d(1:npress)=w_cross(i,1:npress,icross)
      do jklev=1,npress
        jklev1=npress-jklev+1
        w_cross(i,jklev,icross)=work_1d(jklev1)
      enddo

      call interp_spline_1d(pv_cross(i,1:npress,icross), lnp_out, npress, &
 pv_cross_lev(i,1:nlev,icross), lnp_inp_1, nlev, iv, alfa_vinterp, ex_bott, ex_top)

      work_1d(1:npress)=pv_cross(i,1:npress,icross)
      do jklev=1,npress
        jklev1=npress-jklev+1
        pv_cross(i,jklev,icross)=work_1d(jklev1)
      enddo

      call interp_spline_1d(rh_cross(i,1:npress,icross), lnp_out, npress, &
 rh_cross_lev(i,1:nlev,icross), lnp_inp, nlev, iv, alfa_vinterp, ex_bott, ex_top)

      work_1d(1:npress)=rh_cross(i,1:npress,icross)
      do jklev=1,npress
        jklev1=npress-jklev+1
        rh_cross(i,jklev,icross)=work_1d(jklev1)
      enddo

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

! Final rescaling of rel. hum. (from 0 to 100) and pv (to pv units).
! Reset of uphysical values of theta, thetae (may be due to extrapolation below ground)
! and of relative humidity.

    do jklev=1,npress
      jklev1=npress-jklev+1
      do i=1,np
        rh_cross(i,jklev,icross)=max( min( rh_cross(i,jklev,icross)*100., 100.), 0.)
        pv_cross(i,jklev,icross)=pv_cross(i,jklev,icross)*1.e6
        theta_cross(i,jklev,icross)=max( theta_cross(i,jklev,icross), 150. )
        thetae_cross(i,jklev,icross)=max( thetae_cross(i,jklev,icross), 150. )
        if (lnp_out(jklev1) > alog(ps_cross(i,icross))) then
          u_tang_cross(i,jklev,icross)=0.
          v_norm_cross(i,jklev,icross)=0.
          w_cross(i,jklev,icross)=0.
        endif
      enddo
    enddo

  enddo

return
end
!#######################################################################
!----------- Subroutines for coding ouput data in grib2 format ---------------

SUBROUTINE WRITE_GRIB2_HORIZONTAL_GRID(IDATE0,IPERIOD_INP,IPERIOD_ACCUM,IDATEC,IPFLAG,ITFLAG,ISFLAG, &
 ZPH,ZT,ZU,ZV,ZQ,ZRH,ZOM,ZRVA,ZPVA,ZTHETAE,ZPVAT, &
 TG,QG, &
 ZOROGR,ZPZ0,ZTGRS,ZQGRS,ZUS,ZVS,ZTS,ZQS,ZCAPE,ZTMIN,ZTMAX,ZLIFT,ZLEV0, &
 SEAICE_F,SEAICE_TH, LAPSE_RATE_1, LAPSE_RATE_2, IWV, ZCIN, PBL, GUST)

! Procedure that prepares fields of post-processed model data (on horizontal grid of various
! level types) for coding in grib2 format

! Uses module grib2_coding_data in write_grib2_data.F90

USE GRIB2_CODING_DATA
USE MOD_DIMENSIONS, ONLY : NLON, NLAT, NLEV, NLEVG
USE MOD_POSTBOLAM, ONLY : TTR, PDR, &
 PS, TOTPRE, CONPRE, SNFALL, U10, V10, T2, Q2, Q2REL, TD2, HFLUX, QFLUX, &
 CLOUD, CLOUDH, CLOUDM, CLOUDL, SNOW, RUNOFF, RUNOFF_TOT, CSWFL, CLWFL, CSWDFL, CHFLUX, CQFLUX, &
 CWVFLUX, CSOILH_BOTT_FLUX, CSOILW_BOTT_FLUX, &
 ALBEDO, EMISMAP1, EMISMAP2, NLEVPO, NLEVTO, &
 ITROPFLAG, IMWFLAG, &
 H_TROPOPAUSE, P_TROPOPAUSE, T_TROPOPAUSE, U_TROPOPAUSE, V_TROPOPAUSE, Q_TROPOPAUSE, &
 RH_TROPOPAUSE, OM_TROPOPAUSE, RV_TROPOPAUSE, PV_TROPOPAUSE, THE_TROPOPAUSE, P_TROPOPAUSE, &
 H_MAXWIND, P_MAXWIND, T_MAXWIND, U_MAXWIND, V_MAXWIND, Q_MAXWIND, &
 RH_MAXWIND, OM_MAXWIND, RV_MAXWIND, PV_MAXWIND, THE_MAXWIND, P_MAXWIND, &
 QGMIN, QGMAX
USE MOD_RPLEV
USE MOD_RTLEV
USE MOD_SURFADD

IMPLICIT NONE

INTEGER, DIMENSION(5) :: IDATE0, IDATEC
INTEGER, DIMENSION(3) :: IPERIOD_INP, IPERIOD_ACCUM
INTEGER, DIMENSION(12,NLEVPO) :: IPFLAG
INTEGER, DIMENSION(NLEVTO) :: ITFLAG
INTEGER, DIMENSION(80) :: ISFLAG
REAL, DIMENSION(NLON,NLAT,NLEVPO) :: ZPH, ZT, ZU, ZV, ZQ, ZRH, ZOM, ZRVA, ZPVA, ZTHETAE
REAL, DIMENSION(NLON,NLAT,NLEVTO) :: ZPVAT
REAL, DIMENSION(NLON,NLAT,NLEVG) :: TG, QG
REAL, DIMENSION(NLON,NLAT) :: ZOROGR, ZPZ0, ZTGRS, ZQGRS, ZUS, ZVS, ZTS, ZQS, ZCAPE, &
 ZTMIN, ZTMAX, ZLIFT, ZLEV0, SEAICE_F, SEAICE_TH, LAPSE_RATE_1, LAPSE_RATE_2, IWV, &
 ZCIN, PBL, GUST

INTEGER :: NPAR=12, IFIELD, I, K, FLAG_QG, FLAG_TG
REAL, DIMENSION(NLEVG) :: LEV_SOIL
REAL :: PERIOD_SEC

! Name of ouput file

  WRITE (OUTPUT_FILE_NAME,'(A,I4.4,3I2.2,A,I3.3,2I2.2,A)') "bolam_",IDATE0(1:4),"_",IPERIOD_INP(1:3),".grib2"

#ifdef climate
  WRITE (OUTPUT_FILE_NAME,'(A,I4.4,3I2.2,A)') "bolam_",IDATEC(1:4),".grib2"
#endif

! Depth of soil levels

!  LEV_SOIL(1)=PDR(6)*0.5
!  DO K=2,NLEVG
!    LEV_SOIL(K) = LEV_SOIL(K-1)+(PDR(K-1+5)+PDR(K+5))*0.5
!  ENDDO

  LEV_SOIL(1:NLEVG)= PDR(6:5+NLEVG)

  FLAG_QG=0
  IF (ANY(ISFLAG(25:28) == 1)) FLAG_QG=1
  FLAG_TG=0
  IF (ANY(ISFLAG(29:32) == 1)) FLAG_TG=1

! Statistical period in second

  PERIOD_SEC = IPERIOD_ACCUM(1)*86400 + IPERIOD_ACCUM(2)*3600 + IPERIOD_ACCUM(3)*60

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

  DATA(1) % GRIB2_DESCRIPT(1) = 1

! Status and type of data

  DATA(1) % GRIB2_DESCRIPT(6) = 0  ! operational products
  DATA(1) % GRIB2_DESCRIPT(7) = 1  ! forecast
  DATA(1) % GRIB2_DESCRIPT(30) = 0 ! bit-map absent

! Model vertical coordinate parameters

  DATA(1) % N_VERT_COORD_PAR = NLEV+2
  ALLOCATE (DATA(1) % VERT_COORD_PAR (DATA(1) % N_VERT_COORD_PAR) )
  DATA(1) % VERT_COORD_PAR(1)=PDR(40)*0.5
  DO K=2,NLEV
    DATA(1) % VERT_COORD_PAR(K)=(PDR(39+K)+PDR(39+K-1))*0.5 ! sigma of integer levels
  ENDDO
  DATA(1) % VERT_COORD_PAR(NLEV+1) = PDR(37) ! alfa
  DATA(1) % VERT_COORD_PAR(NLEV+2) = PDR(36) ! p0

! Grid parameters

  DATA(1) % GRIB2_DESCRIPT(2) = 1 ! grid template index - horizontal grid
  DATA(1) % NX = NLON
  DATA(1) % NY = NLAT
  DATA(1) % X0 = PDR(39)
  DATA(1) % Y0 = PDR(38)
  DATA(1) % DX = PDR(2)
  DATA(1) % DY = PDR(1)
  DATA(1) % X00 = PDR(5)
  DATA(1) % Y00 = PDR(4) + DATA(1) % DY*0.5

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
    DO I=1,NPAR ! Number of parameters defined at isobaric levels
      IF (IPFLAG(I,K) == 1) THEN
        IFIELD = IFIELD+1
        IF (I == 3) IFIELD = IFIELD+1 ! 1 field for second wind component
      ENDIF
    ENDDO
  ENDDO

! Data at Isentropic (theta) levels (Potential Temperature)

  DO K=1,NLEVTO
    IF (ITFLAG(K) == 1) THEN
      IFIELD = IFIELD+1
    ENDIF
  ENDDO

! Data at Tropopause level

  DO I=1,NPAR ! Number of parameters defined at Tropopause level
    IF (ITROPFLAG(I) == 1) THEN
      IFIELD = IFIELD+1
      IF (I == 3) IFIELD = IFIELD+1 ! 1 field for second wind component
    ENDIF
  ENDDO

! Data at Maximum Wind level

  DO I=1,NPAR ! Number of parameters defined at Maximum Wind level
    IF (IMWFLAG(I) == 1) THEN
      IFIELD = IFIELD+1
      IF (I == 3) IFIELD = IFIELD+1 ! 1 field for second wind component
    ENDIF
  ENDDO

! Data at the surface

  DO I=1,80 ! Number of parameters defined at the surface
    IF (ISFLAG(I) == 1) THEN
      IFIELD = IFIELD+1
      IF (I == 6.OR.I == 14) IFIELD = IFIELD+1 ! 1 field for second wind component
    ENDIF
  ENDDO

! For QG and TG at all soil levels

  IF (FLAG_QG == 1) THEN
    DO I=25,28
      IF (ISFLAG(I) == 1) IFIELD = IFIELD-1
    ENDDO
    DO K=1,NLEVG
      IFIELD = IFIELD+1
    ENDDO
  ENDIF

  IF (FLAG_TG == 1) THEN
    DO I=29,32
      IF (ISFLAG(I) == 1) IFIELD = IFIELD-1
    ENDDO
    DO K=1,NLEVG
      IFIELD = IFIELD+1
    ENDDO
  ENDIF

! For (eventual) cumulated fluxes at soil bottom level and water flux at soil surface, and surface runoff

!!!  IF (ISFLAG(23) == 1) IFIELD = IFIELD+4 ! 23 - Runoff

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
    DO I=1,NPAR ! Number of parameters defined at isobaric levels

      IF (IPFLAG(I,K)==1) THEN

        IFIELD = IFIELD+1

        DATA(IFIELD) % GRIB2_DESCRIPT(10) = 100 ! Isobaric surface  (Pa)
        DATA(IFIELD) % GRIB2_DESCRIPT(11) = NINT(PLEVO(K))
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

        IF (I == 6) THEN ! Omega
          DATA(IFIELD) % GRIB2_DESCRIPT(21) = 2 ! Category: Momentum
          DATA(IFIELD) % GRIB2_DESCRIPT(22) = 8 ! Parameter: Vertical velocity [pressure] (Pa s-1)
          DATA(IFIELD) % FIELD(:,:) = ZOM(:,:,K)
        ENDIF

        IF (I == 7) THEN ! Relative vorticity
          DATA(IFIELD) % GRIB2_DESCRIPT(21) = 2  ! Category: Momentum
          DATA(IFIELD) % GRIB2_DESCRIPT(22) = 12 ! Parameter: Relative vorticity (s-1)
          DATA(IFIELD) % FIELD(:,:) = ZRVA(:,:,K)
        ENDIF

        IF (I == 8) THEN ! Potential vorticity
          DATA(IFIELD) % GRIB2_DESCRIPT(21) = 2  ! Category: Momentum
          DATA(IFIELD) % GRIB2_DESCRIPT(22) = 14 ! Parameter: Potential vorticity (K m2 kg-1 s-1)
          DATA(IFIELD) % FIELD(:,:) = ZPVA(:,:,K)
        ENDIF

        IF (I == 9) THEN ! THETAE
          DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0 ! Category: Temperature
          DATA(IFIELD) % GRIB2_DESCRIPT(22) = 3 ! Parameter: Pseudo-adiabatic potential temperature or equivalent potential temperature (K)
          DATA(IFIELD) % FIELD(:,:) = ZTHETAE(:,:,K)
        ENDIF

      ENDIF

    ENDDO
  ENDDO

! Data at Isentropic (theta) levels (Potential Temperature)

  DO K=1,NLEVTO

    IF (ITFLAG(K) == 1) THEN

      IFIELD = IFIELD+1

      DATA(IFIELD) % GRIB2_DESCRIPT(10) = 107 ! pt Isentropic (theta) level Potential Temp.(K)
      DATA(IFIELD) % GRIB2_DESCRIPT(11) = NINT(TLEVO(K))
      DATA(IFIELD) % GRIB2_DESCRIPT(12) = 0
      DATA(IFIELD) % GRIB2_DESCRIPT(13:14) = IVALMISS
      DATA(IFIELD) % GRIB2_DESCRIPT(20) = 0  ! Discipline: Meteorological products
      DATA(IFIELD) % GRIB2_DESCRIPT(21) = 2  ! Category: Momentum
      DATA(IFIELD) % GRIB2_DESCRIPT(22) = 14 ! Parameter: Potential vorticity (K m2 kg-1 s-1)
      DATA(IFIELD) % FIELD(:,:) = ZPVAT(:,:,K)

    ENDIF

  ENDDO

! Data at the Tropopause level

  DO I=1,NPAR ! Number of parameters defined at isobaric levels

    IF (ITROPFLAG(I)==1) THEN

      IFIELD = IFIELD+1

      DATA(IFIELD) % GRIB2_DESCRIPT(10) = 7 ! Tropopause
      DATA(IFIELD) % GRIB2_DESCRIPT(11:14) = IVALMISS
      DATA(IFIELD) % GRIB2_DESCRIPT(20) = 0   ! Discipline: Meteorological products

      IF (I == 1) THEN ! Geometric Height (m)
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 3 ! Category: Mass
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 6 ! Parameter: Geopotential height (m)
        DATA(IFIELD) % FIELD(:,:) = H_TROPOPAUSE(:,:)
      ENDIF

      IF (I == 2) THEN ! Temperature
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0 ! Category: Temperature
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 0 ! Parameter: Temperature (K)
        DATA(IFIELD) % FIELD(:,:) = T_TROPOPAUSE(:,:)
      ENDIF

      IF (I == 3) THEN ! Wind components
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 2 ! Category: Momentum
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 2 ! Parameter: u-component of wind (m s-1)
        DATA(IFIELD) % FIELD(:,:) = U_TROPOPAUSE(:,:)
        IFIELD = IFIELD+1
        DATA(IFIELD) % GRIB2_DESCRIPT(10) = 7 ! Tropopause
        DATA(IFIELD) % GRIB2_DESCRIPT(11:14) = IVALMISS
        DATA(IFIELD) % GRIB2_DESCRIPT(20) = 0 ! Discipline: Meteorological products
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 2 ! Category: Momentum
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 3 ! Parameter: v-component of wind (m s-1)
        DATA(IFIELD) % FIELD(:,:) = V_TROPOPAUSE(:,:)
      ENDIF

      IF (I == 4) THEN ! Specific humidity
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 1 ! Category: Moisture
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 0 ! Parameter: Specific humidity (kg kg-1)
        DATA(IFIELD) % FIELD(:,:) = Q_TROPOPAUSE(:,:)
      ENDIF

      IF (I == 5) THEN ! Relative humidity
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 1 ! Category: Moisture
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 1 ! Parameter: Relative humidity (%)
        DATA(IFIELD) % FIELD(:,:) = RH_TROPOPAUSE(:,:)
      ENDIF

      IF (I == 6) THEN ! Omega
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 2 ! Category: Momentum
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 8 ! Parameter: Vertical velocity [pressure] (Pa s-1)
        DATA(IFIELD) % FIELD(:,:) = OM_TROPOPAUSE(:,:)
      ENDIF

      IF (I == 7) THEN ! Relative vorticity
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 2  ! Category: Momentum
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 12 ! Parameter: Relative vorticity (s-1)
        DATA(IFIELD) % FIELD(:,:) = RV_TROPOPAUSE(:,:)
      ENDIF

      IF (I == 8) THEN ! Potential vorticity
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 2  ! Category: Momentum
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 14 ! Parameter: Potential vorticity (K m2 kg-1 s-1)
        DATA(IFIELD) % FIELD(:,:) = PV_TROPOPAUSE(:,:)
      ENDIF

      IF (I == 9) THEN ! THETAE
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0 ! Category: Temperature
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 3 ! Parameter: Pseudo-adiabatic potential temperature or equivalent potential temperature (K)
        DATA(IFIELD) % FIELD(:,:) = THE_TROPOPAUSE(:,:)
      ENDIF

      IF (I == 10) THEN ! P
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 3 ! Category: Mass
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 0 ! Parameter: Pressure (Pa)
        DATA(IFIELD) % FIELD(:,:) = P_TROPOPAUSE(:,:)*1.E6 ! PV units
      ENDIF

    ENDIF

  ENDDO

! Data at the Maximum Wind level

  DO I=1,NPAR ! Number of parameters defined at isobaric levels

    IF (IMWFLAG(I)==1) THEN

      IFIELD = IFIELD+1

      DATA(IFIELD) % GRIB2_DESCRIPT(10) = 6 ! Maximum wind level
      DATA(IFIELD) % GRIB2_DESCRIPT(11:14) = IVALMISS
      DATA(IFIELD) % GRIB2_DESCRIPT(20) = 0   ! Discipline: Meteorological products

      IF (I == 1) THEN ! Geometric Height (m)
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 3 ! Category: Mass
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 6 ! Parameter: Geopotential height (m)
        DATA(IFIELD) % FIELD(:,:) = H_MAXWIND(:,:)
      ENDIF

      IF (I == 2) THEN ! Temperature
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0 ! Category: Temperature
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 0 ! Parameter: Temperature (K)
        DATA(IFIELD) % FIELD(:,:) = T_MAXWIND(:,:)
      ENDIF

      IF (I == 3) THEN ! Wind components
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 2 ! Category: Momentum
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 2 ! Parameter: u-component of wind (m s-1)
        DATA(IFIELD) % FIELD(:,:) = U_MAXWIND(:,:)
        IFIELD = IFIELD+1
        DATA(IFIELD) % GRIB2_DESCRIPT(10) = 6 ! Maximum wind level
        DATA(IFIELD) % GRIB2_DESCRIPT(11:14) = IVALMISS
        DATA(IFIELD) % GRIB2_DESCRIPT(20) = 0 ! Discipline: Meteorological products
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 2 ! Category: Momentum
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 3 ! Parameter: v-component of wind (m s-1)
        DATA(IFIELD) % FIELD(:,:) = V_MAXWIND(:,:)
      ENDIF

      IF (I == 4) THEN ! Specific humidity
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 1 ! Category: Moisture
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 0 ! Parameter: Specific humidity (kg kg-1)
        DATA(IFIELD) % FIELD(:,:) = Q_MAXWIND(:,:)
      ENDIF

      IF (I == 5) THEN ! Relative humidity
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 1 ! Category: Moisture
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 1 ! Parameter: Relative humidity (%)
        DATA(IFIELD) % FIELD(:,:) = RH_MAXWIND(:,:)
      ENDIF

      IF (I == 6) THEN ! Omega
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 2 ! Category: Momentum
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 8 ! Parameter: Vertical velocity [pressure] (Pa s-1)
        DATA(IFIELD) % FIELD(:,:) = OM_MAXWIND(:,:)
      ENDIF

      IF (I == 7) THEN ! Relative vorticity
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 2  ! Category: Momentum
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 12 ! Parameter: Relative vorticity (s-1)
        DATA(IFIELD) % FIELD(:,:) = RV_MAXWIND(:,:)
      ENDIF

      IF (I == 8) THEN ! Potential vorticity
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 2  ! Category: Momentum
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 14 ! Parameter: Potential vorticity (K m2 kg-1 s-1)
        DATA(IFIELD) % FIELD(:,:) = PV_MAXWIND(:,:)*1.E6 ! PV units
      ENDIF

      IF (I == 9) THEN ! THETAE
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0 ! Category: Temperature
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 3 ! Parameter: Pseudo-adiabatic potential temperature or equivalent potential temperature (K)
        DATA(IFIELD) % FIELD(:,:) = THE_MAXWIND(:,:)
      ENDIF

      IF (I == 10) THEN ! P
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 3 ! Category: Mass
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 0 ! Parameter: Pressure (Pa)
        DATA(IFIELD) % FIELD(:,:) = P_MAXWIND(:,:)
      ENDIF

    ENDIF

  ENDDO

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

  DO I=1,80 ! Number of parameters defined at the surface


    IF (I >= 25.AND.I <= 32) CYCLE ! QG and TG at model soil levels

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
        DATA(IFIELD) % FIELD(:,:) = TOTPRE(:,:)
      ENDIF

      IF (I == 4) THEN ! Convective precipitation
        DATA(IFIELD) % GRIB2_DESCRIPT(3) = 8   ! Product template: Statistical
        DATA(IFIELD) % GRIB2_DESCRIPT(5) = 1   ! Statistical elaboration type: accumulation
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 1  ! Category: Moisture
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 10 ! Parameter: Convective precipitation (kg m-2)
        DATA(IFIELD) % FIELD(:,:) = CONPRE(:,:)
      ENDIF

      IF (I == 5) THEN ! Solid phase precipitation
        DATA(IFIELD) % GRIB2_DESCRIPT(3) = 8   ! Product template: Statistical
        DATA(IFIELD) % GRIB2_DESCRIPT(5) = 1   ! Statistical elaboration type: accumulation
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 1  ! Category: Moisture
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 29 ! Parameter: Total snowfall (m)
        DATA(IFIELD) % FIELD(:,:) = SNFALL(:,:)*1.E-3
      ENDIF

      IF (I == 6) THEN ! Wind at 10 m
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

      IF (I == 7) THEN ! Temper. at 2m
        DATA(IFIELD) % GRIB2_DESCRIPT(10) = 103 ! Specified height level above ground (m)
        DATA(IFIELD) % GRIB2_DESCRIPT(11) = 2   ! First scaled value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(12) = 0   ! Scale of first value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0   ! Category: Temperature
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 0   ! Parameter: Temperature (K)
        DATA(IFIELD) % FIELD(:,:) = T2(:,:)+TTR
      ENDIF

      IF (I == 8) THEN ! Specific hum. at 2 m
        DATA(IFIELD) % GRIB2_DESCRIPT(10) = 103 ! Specified height level above ground (m)
        DATA(IFIELD) % GRIB2_DESCRIPT(11) = 2   ! First scaled value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(12) = 0   ! Scale of first value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 1   ! Category: Moisture
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 0   ! Parameter: Specific humidity (kg kg-1)
        DATA(IFIELD) % FIELD(:,:) = Q2(:,:)
      ENDIF

      IF (I == 9) THEN ! Relative hum. at 2 m
        DATA(IFIELD) % GRIB2_DESCRIPT(10) = 103 ! Specified height level above ground (m)
        DATA(IFIELD) % GRIB2_DESCRIPT(11) = 2   ! First scaled value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(12) = 0   ! Scale of first value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 1   ! Category: Moisture
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 1   ! Parameter: Relative humidity (%)
        DATA(IFIELD) % FIELD(:,:) = Q2REL(:,:)
      ENDIF

      IF (I == 10) THEN ! Skin Specific hum.
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 1 ! Category: Moisture
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 0 ! Parameter: Specific humidity (kg kg-1)
        DATA(IFIELD) % FIELD(:,:) = ZQGRS(:,:)
      ENDIF

      IF (I == 11) THEN ! Skin temperature
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0 ! Category: Temperature
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 0 ! Parameter: Temperature (K)
        DATA(IFIELD) % FIELD(:,:) = ZTGRS(:,:)+TTR
      ENDIF

      IF (I == 12) THEN ! Flux of sensible heat
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0  ! Category: Temperature
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 11 ! Parameter: Sensible heat net flux (W m-2)
        DATA(IFIELD) % FIELD(:,:) = HFLUX(:,:)
      ENDIF

      IF (I == 13) THEN ! Flux of latent heat
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0  ! Category: Temperature
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 10 ! Parameter: Latent heat net flux (W m-2)
        DATA(IFIELD) % FIELD(:,:) = QFLUX(:,:)
      ENDIF

      IF (I == 14) THEN ! Wind at lowest atm. level
        DATA(IFIELD) % GRIB2_DESCRIPT(10) = 105  ! Hybrid level
        DATA(IFIELD) % GRIB2_DESCRIPT(11) = NLEV ! First scaled value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(12) = 0    ! Scale of first value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 2    ! Category: Momentum
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 2    ! Parameter: u-component of wind (m s-1)
        DATA(IFIELD) % FIELD(:,:) = ZUS(:,:)
        IFIELD = IFIELD+1
        DATA(IFIELD) % GRIB2_DESCRIPT(10) = 105  ! Hybrid level
        DATA(IFIELD) % GRIB2_DESCRIPT(11) = NLEV ! First scaled value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(12) = 0 ! Scale of first value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(13:14) = IVALMISS
        DATA(IFIELD) % GRIB2_DESCRIPT(20) = 0    ! Discipline: Meteorological products
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 2    ! Category: Momentum
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 3    ! Parameter: v-component of wind (m s-1)
        DATA(IFIELD) % FIELD(:,:) = ZVS(:,:)
      ENDIF

      IF (I == 15) THEN ! Temper. at lowest atm. level
        DATA(IFIELD) % GRIB2_DESCRIPT(10) = 105  ! Hybrid level
        DATA(IFIELD) % GRIB2_DESCRIPT(11) = NLEV ! First scaled value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(12) = 0    ! Scale of first value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0    ! Category: Temperature
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 0    ! Parameter: Temperature (K)
        DATA(IFIELD) % FIELD(:,:) = ZTS(:,:)+TTR
      ENDIF

      IF (I == 16) THEN ! Spec. hum. at lowest atm. level
        DATA(IFIELD) % GRIB2_DESCRIPT(10) = 105  ! Hybrid level
        DATA(IFIELD) % GRIB2_DESCRIPT(11) = NLEV ! First scaled value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(12) = 0    ! Scale of first value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 1    ! Category: Moisture
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 0    ! Parameter: Specific humidity (kg kg-1)
        DATA(IFIELD) % FIELD(:,:) = ZQS(:,:)
      ENDIF

      IF (I == 17.OR.I == 21) THEN ! Total cloud cover post proc. or model
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 6 ! Category: Cloud
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 1 ! Parameter: Total cloud cover (%)
        DATA(IFIELD) % FIELD(:,:) = CLOUD(:,:)*1.E2
      ENDIF

      IF (I == 18) THEN ! High cloud cover
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 6 ! Category: Cloud
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 5 ! Parameter: High cloud cover (%)
        DATA(IFIELD) % FIELD(:,:) = CLOUDH(:,:)*1.E2
      ENDIF

      IF (I == 19) THEN ! Middle cloud cover
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 6 ! Category: Cloud
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 4 ! Parameter: Medium cloud cover (%)
        DATA(IFIELD) % FIELD(:,:) = CLOUDM(:,:)*1.E2
      ENDIF

      IF (I == 20) THEN ! Low cloud cover
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 6 ! Category: Cloud
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 3 ! Parameter: Low cloud cover (%)
        DATA(IFIELD) % FIELD(:,:) = CLOUDL(:,:)*1.E2
      ENDIF

      IF (I == 22) THEN ! Snow cover
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 1  ! Category: Moisture
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 13 ! Parameter: Water equivalent of accumulated snow depth (kg m-2)
        DATA(IFIELD) % FIELD(:,:) = SNOW(:,:)*1.E3
      ENDIF

      IF (I == 23) THEN ! Runoff
        DATA(IFIELD) % GRIB2_DESCRIPT(20) = 2 ! Discipline: Land surface products
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0 ! Category: Vegetation/Biomass
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 5 ! Parameter: Surface water runoff (kg m-2
!!!        DATA(IFIELD) % FIELD(:,:) = RUNOFF_TOT(:,:)
        DATA(IFIELD) % FIELD(:,:) = RUNOFF(:,:)

!!!        IFIELD=IFIELD+1
!!!        DATA(IFIELD) % GRIB2_DESCRIPT(10) = 1         ! Ground or water surface
!!!        DATA(IFIELD) % GRIB2_DESCRIPT(11:14) = IVALMISS  ! Parameters of first and second levels (layer)
!!!        DATA(IFIELD) % GRIB2_DESCRIPT(20) = 2 ! Discipline: Land surface products
!!!        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0 ! Category: Vegetation/Biomass
!!!        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 34 ! Parameter: Surface water runoff (kg m-2
!!!        DATA(IFIELD) % FIELD(:,:) = RUNOFF(:,:)
      ENDIF

      IF (I == 24) THEN ! CAPE
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 7 ! Category: Thermodynamic Stability Indices
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 6 ! Parameter: Convective available potential energy (J kg-1)
        DATA(IFIELD) % FIELD(:,:) = ZCAPE(:,:)
      ENDIF

!      IF (I == 25) THEN ! Ground water lev. 1
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
!      IF (I == 26) THEN ! Ground water lev. 3
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
!      IF (I == 27) THEN ! Ground water lev. 5
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
!      IF (I == 28) THEN ! Ground water lev. 7
!        K=7
!        DATA(IFIELD) % GRIB2_DESCRIPT(10) = 106 ! Depth below land surface (m)
!        DATA(IFIELD) % GRIB2_DESCRIPT(11) = NINT(LEV_SOIL(K)*1.E3) ! First scaled value of level (layer)
!        DATA(IFIELD) % GRIB2_DESCRIPT(12) = 3   ! Scale of first value of level (layer)
!        DATA(IFIELD) % GRIB2_DESCRIPT(20) = 2   ! Discipline: Land surface products
!        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0   ! Category: Vegetation/Biomass
!        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 9   ! Parameter: Volumetric soil moisture content (Proportion)
!        DATA(IFIELD) % FIELD(:,:) = QG(:,:,K)
!      ENDIF
!
!      IF (I == 29) THEN ! Ground temp. lev. 1
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
!      IF (I == 30) THEN ! Ground temp. lev. 3
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
!      IF (I == 31) THEN ! Ground temp. lev. 5
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
!      IF (I == 32) THEN ! Ground temp. lev. 7
!        K=7
!        DATA(IFIELD) % GRIB2_DESCRIPT(10) = 106 ! Depth below land surface (m)
!        DATA(IFIELD) % GRIB2_DESCRIPT(11) = NINT(LEV_SOIL(K)*1.E3) ! First scaled value of level (layer)
!        DATA(IFIELD) % GRIB2_DESCRIPT(12) = 3   ! Scale of first value of level (layer)
!        DATA(IFIELD) % GRIB2_DESCRIPT(20) = 2   ! Discipline: Land surface products
!        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0   ! Category: Vegetation/Biomass
!        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 2   ! Parameter: Soil temperature (K)
!        DATA(IFIELD) % FIELD(:,:) = TG(:,:,K)+TTR
!      ENDIF

      IF (I == 33) THEN ! Ground temp. at 5 cm depth
        DATA(IFIELD) % GRIB2_DESCRIPT(10) = 106 ! Depth below land surface (m)
        DATA(IFIELD) % GRIB2_DESCRIPT(11) = 5   ! First scaled value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(12) = 2   ! Scale of first value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(20) = 2   ! Discipline: Land surface products
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0   ! Category: Vegetation/Biomass
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 2   ! Parameter: Soil temperature (K)
        DATA(IFIELD) % FIELD(:,:) = TG005(:,:)+TTR
      ENDIF

      IF (I == 34) THEN ! Ground temp. at 10 cm depth
        DATA(IFIELD) % GRIB2_DESCRIPT(10) = 106 ! Depth below land surface (m)
        DATA(IFIELD) % GRIB2_DESCRIPT(11) = 10  ! First scaled value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(12) = 2   ! Scale of first value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(20) = 2   ! Discipline: Land surface products
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0   ! Category: Vegetation/Biomass
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 2   ! Parameter: Soil temperature (K)
        DATA(IFIELD) % FIELD(:,:) = TG010(:,:)+TTR
      ENDIF

      IF (I == 35) THEN ! Ground temp. at 20 cm depth
        DATA(IFIELD) % GRIB2_DESCRIPT(10) = 106 ! Depth below land surface (m)
        DATA(IFIELD) % GRIB2_DESCRIPT(11) = 20  ! First scaled value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(12) = 2   ! Scale of first value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(20) = 2   ! Discipline: Land surface products
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0   ! Category: Vegetation/Biomass
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 2   ! Parameter: Soil temperature (K)
        DATA(IFIELD) % FIELD(:,:) = TG020(:,:)+TTR
      ENDIF

      IF (I == 36) THEN ! Ground temp. at 50 cm depth
        DATA(IFIELD) % GRIB2_DESCRIPT(10) = 106 ! Depth below land surface (m)
        DATA(IFIELD) % GRIB2_DESCRIPT(11) = 50  ! First scaled value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(12) = 2   ! Scale of first value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(20) = 2   ! Discipline: Land surface products
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0   ! Category: Vegetation/Biomass
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 2   ! Parameter: Soil temperature (K)
        DATA(IFIELD) % FIELD(:,:) = TG050(:,:)+TTR
      ENDIF

      IF (I == 37) THEN ! Ground temp. at 100 cm depth
        DATA(IFIELD) % GRIB2_DESCRIPT(10) = 106 ! Depth below land surface (m)
        DATA(IFIELD) % GRIB2_DESCRIPT(11) = 100 ! First scaled value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(12) = 2   ! Scale of first value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(20) = 2   ! Discipline: Land surface products
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0   ! Category: Vegetation/Biomass
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 2   ! Parameter: Soil temperature (K)
        DATA(IFIELD) % FIELD(:,:) = TG100(:,:)+TTR
      ENDIF

      IF (I == 38) THEN ! T min. 2m
        DATA(IFIELD) % GRIB2_DESCRIPT(3) = 8    ! Product template: Statistical
        DATA(IFIELD) % GRIB2_DESCRIPT(5) = 3    ! Statistical elaboration type: minimum
        DATA(IFIELD) % GRIB2_DESCRIPT(10) = 103 ! Specified height level above ground (m)
        DATA(IFIELD) % GRIB2_DESCRIPT(11) = 2   ! First scaled value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(12) = 0   ! Scale of first value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0   ! Category: Temperature
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 0   ! Parameter: Temperature (K)
        DATA(IFIELD) % FIELD(:,:) = ZTMIN(:,:)+TTR
      ENDIF

      IF (I == 39) THEN ! T max. 2m
        DATA(IFIELD) % GRIB2_DESCRIPT(3) = 8    ! Product template: Statistical
        DATA(IFIELD) % GRIB2_DESCRIPT(5) = 2    ! Statistical elaboration type: maximum
        DATA(IFIELD) % GRIB2_DESCRIPT(10) = 103 ! Specified height level above ground (m)
        DATA(IFIELD) % GRIB2_DESCRIPT(11) = 2   ! First scaled value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(12) = 0   ! Scale of first value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0   ! Category: Temperature
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 0   ! Parameter: Temperature (K)
        DATA(IFIELD) % FIELD(:,:) = ZTMAX(:,:)+TTR
      ENDIF

      IF (I == 40) THEN ! Temper. at 0.5 m above surf.
        DATA(IFIELD) % GRIB2_DESCRIPT(10) = 103 ! Specified height level above ground (m)
        DATA(IFIELD) % GRIB2_DESCRIPT(11) = 50  ! First scaled value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(12) = 2   ! Scale of first value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0   ! Category: Temperature
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 0   ! Parameter: Temperature (K)
        DATA(IFIELD) % FIELD(:,:) = T05(:,:)+TTR
      ENDIF

      IF (I == 41) THEN ! Relative hum. at 0.5 m above surf.
        DATA(IFIELD) % GRIB2_DESCRIPT(10) = 103 ! Specified height level above ground (m)
        DATA(IFIELD) % GRIB2_DESCRIPT(11) = 50  ! First scaled value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(12) = 2   ! Scale of first value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 1   ! Category: Moisture
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 1   ! Parameter: Relative humidity (%)
        DATA(IFIELD) % FIELD(:,:) = Q05REL(:,:)
      ENDIF

      IF (I == 42) THEN ! Temper. at 0.05 m above surf.
        DATA(IFIELD) % GRIB2_DESCRIPT(10) = 103 ! Specified height level above ground (m)
        DATA(IFIELD) % GRIB2_DESCRIPT(11) = 5   ! First scaled value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(12) = 2   ! Scale of first value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0   ! Category: Temperature
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 0   ! Parameter: Temperature (K)
        DATA(IFIELD) % FIELD(:,:) = T005(:,:)+TTR
      ENDIF

      IF (I == 43) THEN ! Cumulated short wave radiative flux
        DATA(IFIELD) % GRIB2_DESCRIPT(3) = 8  ! Product template: Statistical
        DATA(IFIELD) % GRIB2_DESCRIPT(5) = 0  ! Statistical elaboration type: average
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 4 ! Category: Short-wave Radiation
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 0 ! Parameter: Net short-wave radiation flux (surface) (W m-2)
        IF (PERIOD_SEC > 0) THEN
          DATA(IFIELD) % FIELD(:,:) = CSWFL(:,:)*1.E3/PERIOD_SEC
        ELSE
          DATA(IFIELD) % FIELD(:,:) = CSWFL(:,:)*1.E3
        ENDIF
      ENDIF

      IF (I == 44) THEN ! Cumulated long wave radiative flux
        DATA(IFIELD) % GRIB2_DESCRIPT(3) = 8  ! Product template: Statistical
        DATA(IFIELD) % GRIB2_DESCRIPT(5) = 0  ! Statistical elaboration type: average
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 5 ! Category: Long-wave Radiation
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 0 ! Parameter: Net long wave radiation flux (surface) (W m-2)
        IF (PERIOD_SEC > 0) THEN
          DATA(IFIELD) % FIELD(:,:) = CLWFL(:,:)*1.E3/PERIOD_SEC
        ELSE
          DATA(IFIELD) % FIELD(:,:) = CLWFL(:,:)*1.E3
        ENDIF
      ENDIF

      IF (I == 45) THEN ! Cumulated sensible heat flux
        DATA(IFIELD) % GRIB2_DESCRIPT(3) = 8   ! Product template: Statistical
        DATA(IFIELD) % GRIB2_DESCRIPT(5) = 0   ! Statistical elaboration type: average
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0  ! Category: Temperature
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 11 ! Parameter: Sensible heat net flux (W m-2)
        IF (PERIOD_SEC > 0) THEN
          DATA(IFIELD) % FIELD(:,:) = CHFLUX(:,:)*1.E3/PERIOD_SEC
        ELSE
          DATA(IFIELD) % FIELD(:,:) = CHFLUX(:,:)*1.E3
        ENDIF
      ENDIF

      IF (I == 46) THEN ! Cumulated latent heat flux
        DATA(IFIELD) % GRIB2_DESCRIPT(3) = 8   ! Product template: Statistical
        DATA(IFIELD) % GRIB2_DESCRIPT(5) = 0   ! Statistical elaboration type: average
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0  ! Category: Temperature
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 10 ! Parameter: Latent heat net flux (W m-2)
        IF (PERIOD_SEC > 0) THEN
          DATA(IFIELD) % FIELD(:,:) = CQFLUX(:,:)*1.E3/PERIOD_SEC
        ELSE
          DATA(IFIELD) % FIELD(:,:) = CQFLUX(:,:)*1.E3
        ENDIF
      ENDIF

      IF (I == 47) THEN ! Albedo
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 19 ! Category: Physical atmospheric properties
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 1  ! Parameter: Albedo (%)
        DATA(IFIELD) % FIELD(:,:) = ALBEDO(:,:)*1.E2
      ENDIF

      IF (I == 48) THEN ! Emissivity (broadband)
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 19  ! Category: Physical atmospheric properties
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 201 ! Parameter: local use Emissivity (%)
        DATA(IFIELD) % FIELD(:,:) = EMISMAP1(:,:)*1.E2
      ENDIF

      IF (I == 49) THEN ! Emissivity (window)
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 19  ! Category: Physical atmospheric properties
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 202 ! Parameter: local use Emissivity (%)
        DATA(IFIELD) % FIELD(:,:) = EMISMAP2(:,:)*1.E2
      ENDIF

      IF (I == 50) THEN ! Lifted index
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 7  ! Category: Thermodynamic Stability Indices
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 10 ! Parameter: Surface lifted index (K)
        DATA(IFIELD) % FIELD(:,:) = ZLIFT(:,:)
      ENDIF

      IF (I == 51) THEN ! Level of 0 C
        DATA(IFIELD) % GRIB2_DESCRIPT(10) = 4 ! Level of 0o C isotherm
        DATA(IFIELD) % GRIB2_DESCRIPT(11) = 0 ! First scaled value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(12) = 0 ! Scale of first value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 3 ! Category: Mass
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 6 ! Parameter: Geometric height (m)
        DATA(IFIELD) % FIELD(:,:) = ZLEV0(:,:)
      ENDIF

      IF (I == 52) THEN ! Sea ice fraction
        DATA(IFIELD) % GRIB2_DESCRIPT(20) = 10 ! Discipline: Oceanographic products
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 2  ! Category: Ice
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 0  ! Parameter: Ice cover (Proportion)
        DATA(IFIELD) % FIELD(:,:) = SEAICE_F(:,:)
      ENDIF

      IF (I == 53) THEN ! Sea ice thickness (m)
        DATA(IFIELD) % GRIB2_DESCRIPT(20) = 10 ! Discipline: Oceanographic products
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 2  ! Category: Ice
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 1  ! Parameter: Ice thickness (m)
        DATA(IFIELD) % FIELD(:,:) = SEAICE_TH(:,:)
      ENDIF

      IF (I == 54) THEN ! Dew point temperature at 2 m
        DATA(IFIELD) % GRIB2_DESCRIPT(10) = 103 ! Specified height level above ground (m)
        DATA(IFIELD) % GRIB2_DESCRIPT(11) = 2   ! First scaled value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(12) = 0   ! Scale of first value of level (layer)
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0   ! Category: Temperature
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 6   ! Parameter: Dew point temperature (K)
        DATA(IFIELD) % FIELD(:,:) = TD2(:,:)+TTR
      ENDIF

      IF (I == 55) THEN ! Integrated Water Vapour (kg/m2)
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 1  ! Category: Moisture
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 64 ! Parameter: Total column integrated water vapour (kg m-2)
        DATA(IFIELD) % FIELD(:,:) = IWV(:,:)
      ENDIF

      IF (I == 56) THEN ! Convective inibition (J/kg)
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 7 ! Category: Thermodynamic Stability Indices
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 7 ! Parameter: Convective inhibition (J kg-1)
        DATA(IFIELD) % FIELD(:,:) = ZCIN(:,:)
      ENDIF

      IF (I == 57) THEN ! PBL top height (above surface)
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 3  ! Category: Mass
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 18 ! Parameter: Planetary boundary layer height (M)
        DATA(IFIELD) % FIELD(:,:) = PBL(:,:)
      ENDIF

      IF (I == 58) THEN ! Wind Gust
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 2  ! Category: Momentum
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 22 ! Parameter: Wind speed (gust)  (M sec-1)
        DATA(IFIELD) % FIELD(:,:) = GUST(:,:)
      ENDIF

      IF (I == 59) THEN ! Average short wave radiative downward flux
        DATA(IFIELD) % GRIB2_DESCRIPT(3) = 8  ! Product template: Statistical
        DATA(IFIELD) % GRIB2_DESCRIPT(5) = 0  ! Statistical elaboration type: average
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 4 ! Category: Short-wave Radiation
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 7 ! Parameter: Downward short-wave radiation flux (surface) (W m-2)
        IF (PERIOD_SEC > 0) THEN
          DATA(IFIELD) % FIELD(:,:) = CSWDFL(:,:)/PERIOD_SEC
        ELSE
          DATA(IFIELD) % FIELD(:,:) = CSWDFL(:,:)
        ENDIF
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

      IF (I == 61) THEN ! Surface Pressure (Pa)
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 3  ! Category: Mass
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 0  ! Parameter: Pressure (Pa)
        DATA(IFIELD) % FIELD(1:NLON,1:NLAT) = PS(1:NLON,1:NLAT)
      ENDIF

    ENDIF

  ENDDO

!!!  IF (ISFLAG(23) == 1) THEN ! 23 - Runoff
!!!
!!!    IFIELD = IFIELD+1
!!!    DATA(IFIELD) % GRIB2_DESCRIPT(10) = 1 ! Ground or water surface
!!!    DATA(IFIELD) % GRIB2_DESCRIPT(20) = 2   ! Discipline: Land surface products
!!!    DATA(IFIELD) % GRIB2_DESCRIPT(21) = 3   ! Category: Soil Products
!!!    DATA(IFIELD) % GRIB2_DESCRIPT(22) = 201 ! Parameter: Water vapour flux at soil surface (kg m-2 s-1)
!!!    IF (PERIOD_SEC > 0) THEN
!!!      DATA(IFIELD) % FIELD(:,:) = CWVFLUX(:,:)/PERIOD_SEC
!!!    ELSE
!!!      DATA(IFIELD) % FIELD(:,:) = CWVFLUX(:,:)
!!!    ENDIF
!!!
!!!    IFIELD = IFIELD+1
!!!    DATA(IFIELD) % GRIB2_DESCRIPT(10) = 106 ! Depth below land surface (m)
!!!    DATA(IFIELD) % GRIB2_DESCRIPT(11) = NINT(LEV_SOIL(NLEVG)*1.E3) ! First scaled value of level (layer)
!!!    DATA(IFIELD) % GRIB2_DESCRIPT(12) = 3   ! Scale of first value of level (layer)
!!!    DATA(IFIELD) % GRIB2_DESCRIPT(20) = 2   ! Discipline: Land surface products
!!!    DATA(IFIELD) % GRIB2_DESCRIPT(21) = 3   ! Category: Soil Products
!!!    DATA(IFIELD) % GRIB2_DESCRIPT(22) = 202 ! Parameter: Heat flux at soil bottom level (W m-2)
!!!    IF (PERIOD_SEC > 0) THEN
!!!      DATA(IFIELD) % FIELD(:,:) = CSOILH_BOTT_FLUX(:,:)/PERIOD_SEC
!!!    ELSE
!!!      DATA(IFIELD) % FIELD(:,:) = CSOILH_BOTT_FLUX(:,:)
!!!    ENDIF
!!!
!!!    IFIELD = IFIELD+1
!!!    DATA(IFIELD) % GRIB2_DESCRIPT(10) = 106 ! Depth below land surface (m)
!!!    DATA(IFIELD) % GRIB2_DESCRIPT(11) = NINT(LEV_SOIL(NLEVG)*1.E3) ! First scaled value of level (layer)
!!!    DATA(IFIELD) % GRIB2_DESCRIPT(12) = 3   ! Scale of first value of level (layer)
!!!    DATA(IFIELD) % GRIB2_DESCRIPT(20) = 2   ! Discipline: Land surface products
!!!    DATA(IFIELD) % GRIB2_DESCRIPT(21) = 3   ! Category: Soil Products
!!!    DATA(IFIELD) % GRIB2_DESCRIPT(22) = 203 ! Parameter: Water flux at soil bottom level (kg m-2 s-1)
!!!    IF (PERIOD_SEC > 0) THEN
!!!      DATA(IFIELD) % FIELD(:,:) = CSOILW_BOTT_FLUX(:,:)/PERIOD_SEC
!!!    ELSE
!!!      DATA(IFIELD) % FIELD(:,:) = CSOILW_BOTT_FLUX(:,:)
!!!    ENDIF
!!!
!!!  ENDIF

! End of preparation of output data for coding in grib2 format

! Write output data in grib2 format

  CALL WRITE_GRIB2_DATA

RETURN
END
!#######################################################################
SUBROUTINE WRITE_GRIB2_CROSS_SPACE(IDATE0, IPERIOD_INP, IPERIOD_ACCUM, IDATEC, &
 NCROSS, NLEVEL, NPOINT_MAX, NPOINT, PLEVEL, LONINI_CROSS, LATINI_CROSS, LONFIN_CROSS, LATFIN_CROSS, &
 X0, Y0, PS_CROSS, &
 THETA_CROSS, THETAE_CROSS, U_TANG_CROSS, V_NORM_CROSS, W_CROSS, PV_CROSS, RH_CROSS)

! Procedure that prepares fields of post-processed model data (on various space-crossing grids)
! for coding in grib2 format

! Uses module grib2_coding_data in write_grib2_data.F90

USE GRIB2_CODING_DATA

IMPLICIT NONE

INTEGER, DIMENSION(5) :: IDATE0, IDATEC
INTEGER, DIMENSION(3) :: IPERIOD_INP, IPERIOD_ACCUM
INTEGER :: NCROSS, NLEVEL, NPOINT_MAX
INTEGER :: NPARAM = 8 ! (prognostic parameters + "orography")
INTEGER, DIMENSION(NCROSS) :: NPOINT
REAL, DIMENSION(NLEVEL) :: PLEVEL
REAL, DIMENSION(NCROSS) :: LONINI_CROSS, LATINI_CROSS, LONFIN_CROSS, LATFIN_CROSS
REAL :: X0, Y0

REAL, DIMENSION(NPOINT_MAX,NCROSS) :: PS_CROSS

REAL, DIMENSION(NPOINT_MAX,NLEVEL,NCROSS) :: THETA_CROSS, THETAE_CROSS, &
 U_TANG_CROSS, V_NORM_CROSS, W_CROSS, PV_CROSS, RH_CROSS

INTEGER :: IFIELD, ICROSS, IPARAM, NX, NZ

! Name of ouput file

  WRITE (OUTPUT_FILE_NAME,'(A,I4.4,3I2.2,A,I3.3,2I2.2,A)') "bolam_cross_",IDATE0(1:4),"_",IPERIOD_INP(1:3),".grib2"

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

  DATA(1) % GRIB2_DESCRIPT(1) = 1

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
        DATA(IFIELD) % VERT_COORD_PAR(:) = PLEVEL(:)
      ELSE ! "Orography"
        NZ = 1
        DATA(IFIELD) % NY = NZ
        ALLOCATE(DATA(IFIELD) % VERT_COORD_PAR(NZ))
        DATA(IFIELD) % VERT_COORD_PAR(:) = VALMISS
      ENDIF
      DATA(IFIELD) % IND_VERT_COORD = 100 ! Vertical coordinate type: Pressure (Pa)

! Definition of data fields

      ALLOCATE(DATA(IFIELD) % FIELD(NX,NZ))

      DATA(IFIELD) % GRIB2_DESCRIPT(20) = 0 ! Discipline: Meteorological products

      IF (IPARAM == 1) THEN ! Potential temperature (K)
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0 ! Category: Temperature
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 2 ! Parameter: Potential temperature (K)
        DATA(IFIELD) % FIELD(:,:) = THETA_CROSS(1:NX,1:NZ,ICROSS)
      ENDIF

      IF (IPARAM == 2) THEN ! Equivalent potential temperature
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 0 ! Category: Temperature
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 3 ! Parameter: Potential temperature (K)
        DATA(IFIELD) % FIELD(:,:) = THETAE_CROSS(1:NX,1:NZ,ICROSS)
      ENDIF

      IF (IPARAM == 3) THEN ! U-component of wind
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 2 ! Category: Momentum
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 2 ! Parameter: u-component of wind (m s-1)
        DATA(IFIELD) % FIELD(:,:) = U_TANG_CROSS(1:NX,1:NZ,ICROSS)
      ENDIF

      IF (IPARAM == 4) THEN ! V-component of wind
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 2 ! Category: Momentum
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 3 ! Parameter: v-component of wind (m s-1)
        DATA(IFIELD) % FIELD(:,:) = V_NORM_CROSS(1:NX,1:NZ,ICROSS)
      ENDIF

      IF (IPARAM == 5) THEN ! Vertical velocity
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 2 ! Category: Momentum
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 9 ! Parameter: Vertical velocity [geometric] (m s-1)
        DATA(IFIELD) % FIELD(:,:) = W_CROSS(1:NX,1:NZ,ICROSS)
      ENDIF

      IF (IPARAM == 6) THEN ! Potential vorticity
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 2  ! Category: Momentum
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 14 ! Parameter: Potential vorticity (K m2 kg-1 s-1)
        DATA(IFIELD) % FIELD(:,:) = PV_CROSS(1:NX,1:NZ,ICROSS)
      ENDIF

      IF (IPARAM == 7) THEN ! Relative humidity
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 1 ! Category: Moisture
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 1 ! Parameter: Relative humidity (%)
        DATA(IFIELD) % FIELD(:,:) = RH_CROSS(1:NX,1:NZ,ICROSS)
      ENDIF

      IF (IPARAM == NPARAM) THEN ! "Orography" (Surface pressure)
        DATA(IFIELD) % GRIB2_DESCRIPT(21) = 3 ! Category: Mass
        DATA(IFIELD) % GRIB2_DESCRIPT(22) = 0 ! Parameter: Pressure (Pa)
        DATA(IFIELD) % FIELD(:,1) = PS_CROSS(1:NX,ICROSS)
      ENDIF

    ENDDO

  ENDDO

! End of preparation of output data for coding in grib2 format

! Write output data in grib2 format

  CALL WRITE_GRIB2_DATA

RETURN
END
!#######################################################################
SUBROUTINE WRITE_GRIB2_CROSS_SPACE_I_J(IDATE0, IPERIOD_INP, IPERIOD_ACCUM, IDATEC, &
 NLON, NLAT, NLEVEL,  NCROSSY, NCROSSX, PLEVEL, LON_CROSSY, LAT_CROSSY, LON_CROSSX, LAT_CROSSX, &
 X0, Y0, PS_CROSSY, PS_CROSSX, &
 THETA_CROSSY, THETAE_CROSSY, U_CROSSY, V_CROSSY, W_CROSSY, PV_CROSSY, RH_CROSSY, &
 THETA_CROSSX, THETAE_CROSSX, U_CROSSX, V_CROSSX, W_CROSSX, PV_CROSSX, RH_CROSSX)

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
INTEGER :: NPARAM = 8 ! (prognostic parameters + "orography")
REAL, DIMENSION(NLEVEL) :: PLEVEL
REAL, DIMENSION(NCROSSY,2) :: LON_CROSSY, LAT_CROSSY
REAL, DIMENSION(NCROSSX,2) :: LON_CROSSX, LAT_CROSSX
REAL :: X0, Y0

REAL, DIMENSION(NLAT,1,NCROSSY) :: PS_CROSSY
REAL, DIMENSION(NLON,1,NCROSSX) :: PS_CROSSX

REAL, DIMENSION(NLAT,NLEVEL,NCROSSY) :: THETA_CROSSY, THETAE_CROSSY, U_CROSSY, V_CROSSY, W_CROSSY, &
 PV_CROSSY, RH_CROSSY
REAL, DIMENSION(NLON,NLEVEL,NCROSSX) :: THETA_CROSSX, THETAE_CROSSX, U_CROSSX, V_CROSSX, W_CROSSX, &
 PV_CROSSX, RH_CROSSX

INTEGER :: IFIELD, IAXIS, NCROSS, ICROSS, IPARAM, NX, NZ

! Name of ouput file

  WRITE (OUTPUT_FILE_NAME,'(A,I4.4,3I2.2,A,I3.3,2I2.2,A)') "bolam_cross_i_j_",IDATE0(1:4),"_",IPERIOD_INP(1:3),".grib2"

! Cleaning of possible previous allocations

  DO IFIELD=1,NFIELD0
    IF (ALLOCATED(DATA(IFIELD) % FIELD) ) DEALLOCATE (DATA(IFIELD) % FIELD)
    IF (ALLOCATED(DATA(IFIELD) % VERT_COORD_PAR) ) DEALLOCATE (DATA(IFIELD) % VERT_COORD_PAR)
    DATA(IFIELD) % GRIB2_DESCRIPT(:) = IVALMISS
    DATA(IFIELD) % GRIB2_DESCRIPT(11:14) = 0
  ENDDO

! Description of data for grib2 coding (INTEGER, DIMENSION(50) :: GRIB2_DESCRIPT=IVALMISS,
! see module_write_grib2_data.F90)

! Content of GRIB2_DESCRPT array - see bottom of file

! Model index: 1 - Bolam, 2 - Moloch, 3 - Globo

  DATA(1) % GRIB2_DESCRIPT(1) = 1

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
          DATA(IFIELD) % VERT_COORD_PAR(:) = PLEVEL(:)
        ELSE ! "Orography"
          NZ = 1
          DATA(IFIELD) % NY = NZ
          ALLOCATE(DATA(IFIELD) % VERT_COORD_PAR(NZ))
          DATA(IFIELD) % VERT_COORD_PAR(:) = VALMISS
        ENDIF
        DATA(IFIELD) % IND_VERT_COORD = 100 ! Vertical coordinate type: Pressure (Pa)

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

        IF (IPARAM == 6) THEN ! Potential vorticity
          DATA(IFIELD) % GRIB2_DESCRIPT(21) = 2  ! Category: Momentum
          DATA(IFIELD) % GRIB2_DESCRIPT(22) = 14 ! Parameter: Potential vorticity (K m2 kg-1 s-1)
          IF (IAXIS == 1) THEN ! y-axis
            DATA(IFIELD) % FIELD(:,:) = PV_CROSSY(1:NX,1:NZ,ICROSS)
          ELSE ! x-axis
            DATA(IFIELD) % FIELD(:,:) = PV_CROSSX(1:NX,1:NZ,ICROSS)
          ENDIF
        ENDIF

        IF (IPARAM == 7) THEN ! Relative humidity
          DATA(IFIELD) % GRIB2_DESCRIPT(21) = 1 ! Category: Moisture
          DATA(IFIELD) % GRIB2_DESCRIPT(22) = 1 ! Parameter: Relative humidity (%)
          IF (IAXIS == 1) THEN ! y-axis
            DATA(IFIELD) % FIELD(:,:) = RH_CROSSY(1:NX,1:NZ,ICROSS)
          ELSE ! x-axis
            DATA(IFIELD) % FIELD(:,:) = RH_CROSSX(1:NX,1:NZ,ICROSS)
          ENDIF
        ENDIF

        IF (IPARAM == NPARAM) THEN ! "Orography" (Surface pressure)
          DATA(IFIELD) % GRIB2_DESCRIPT(21) = 3 ! Category: Mass
          DATA(IFIELD) % GRIB2_DESCRIPT(22) = 0 ! Parameter: Pressure (Pa)
          IF (IAXIS == 1) THEN ! y-axis
            DATA(IFIELD) % FIELD(:,1) = PS_CROSSY(1:NX,1,ICROSS)
          ELSE ! x-axis
            DATA(IFIELD) % FIELD(:,1) = PS_CROSSX(1:NX,1,ICROSS)
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
!#######################################################################
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
