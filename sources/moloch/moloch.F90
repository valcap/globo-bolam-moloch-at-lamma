! Last update 16/10/2024

! Version 24.1.1

! Sett. 2024: Cambiata la scrittura di model_param_constant.bin - piu' grandezze

! Giu. 2024: Nuovo formato mhf (scrittura record non 1d ma 2d)

! Ott. 2023: Option of compilation without ECMWF radiation scheme pachage,
! astronomical variables definition and aerosols and ozone content definitions
! are extrated from that scheme.
  
! Feb. 2023: Restart option

! Nov. 2022: Introduction in soil_mhf output 1 field 2d with atmospheric level
! number when ice precipitation predominate liquid precipitation
! (for snowfall altitude definition in postprocessing elaboration)

! May 2022: Correction of vegetation fraction in subroutine surfradpar

! Dec. 2021: Input/output units have been reordered, for each type of data
! a specific unit number has been assigned.

! Dec. 2021: Correction of error in pochva scheme, the definition of qv_surf
! (air humidity) at the surface (mainly at bare soil surface) as the function
! of turbulent coefficient in the surface layer and top soil water content,
! error was due to turbulent coefficient that was not normalized to
! 1 m over surface.

! Nov. 2021: Substitution !$ symbols for mpi library implementation by ifdef mpi command

! Oct. 2021: nudging of T, U, V (from Version 19)

! Dic. 2020: divergence damping (subroutine divdamp)

! Set. 2020: parametrizzazione della convezione (subroutine convection_kf)

! Ago. 2019: vertical advection (wafone) may be (optionaly) executed
! in half time step version (nltwice)

! Gen. 2019: turb. orograf. drag (subr. tofd)

! Ago. 2018: 
! Nuovo schema del suolo (pochva), i nuovi dataset fisiografici,
! parametri del suolo sono diversi su vari livelli,
! manto nevoso muti-stato, il numero dei livelli di neve e' 11,
! la coordinate verticale nello manto nevoso e' kg/m/m,
! quando su un livello la neve non c'e', tutti i paramtri di neve sono -9999.;
! nuovo formato mhf: mfs_atm, mhf_soil e il file con i campi statici (model_param_constant.bin).

! Version theta-pai (v 17)

! Apr. 2018: aggiunto effetto della slope (slopeff) sul flusso radiativo solare diretto.
! Corretto ist. di chiamata della radiaz. ECMWF.
! Aggiunta riga per aggiorn. parallelismo nell'avv. di u e v.
! Diminuita diff. orizz. (elim. diff. su incrementi)
! Mar. 2018: corretto problema diff. vert. qcw e qci (lower b. c.)
! Dic. 2017: corretto errore calcolo cloudt
! Nov. 2017: rivista traspiraz. vegetaz. (messa in funz. della rad. netta) e stominr
! (problema eccessiva RH a 2m al tramonto).
! Stabilizz. codice posponendo il calcolo della mix. length nel caso stabile
! e minimizzando la tend. di T (problema nel caso di eccessiva evaporazione dalle foglie).
! Ott. 2017: modif. calcolo t2m e vento 10 m su montagne "rough" (v. "fix over").
! Aumentata conduc. termica suolo su montagne > 1000 m (e rough)
! Per problemi di stabilita':
! lasciato come nel moloch precedente il ritardo di roscdm e roscdt;
! resa implicita la diff. vert. di TKE e w (commentata la vers. esplicita);
! la stabilita' dipende anche dalla def. di mixing length mlz nel caso stabile:
! qui lasciata Deardorff (sembra migliore per la precip.), ma Blackadar e' piu' stabile.
! Commentate chiamate a hdifft e hdiffu (turbolenza).
! Sett. 2017: rivista cloudfr (riduz. nubi).
! Giu. 2017: introd. irrigazione e risaie (in soil - v. irrig e ricefields)
! Corretti vari parametri che definiscono le variabili su suolo urbano
! (ma da rivedere la def. di qg in funz. di qgwilt su suolo urbano).
! Tolta subroutine advect_waf_12 (su pc forse piu' veloce ma su cluster piu' lenta)
! Mag. 2017:
! Il param. mhfr nella namelist (moloch.inp):
! mhfr = 1: full mhf horiz. resolution; 2: half mhf horiz. resolution.
! Nuova vdiff: scorporato calcolo surface layer (ore in surflayer), vdiff usa static stability
! Durran&Klemp, nuove mixing length; nuovo calcolo nubi.
! Calcolo diff. orizzontale turbolenta per T e U (hdiffT e hdiffU).
! Gen. 2017: introdotti nuovi campi diagnostici superficiali in shf:
! average latent heat flux, average sensible heat flux, average evaporation flux
! shf_accum, lhf_accum, qf_accum (su richiesta CNR-IAMC).
! Rivista def. snow (snow cover) in funz. della neve prevista a 24 ore

! Direttive ifdef: "oper" per la gestione operativa.

! -----------------------------------------------------------------------------------------------------
    module mod_moloch

    include 'dimensions.inc'

!    integer, parameter :: gnlon=1154, gnlat=1154, nlev=60, nlevg=7
!    integer, parameter :: nprocsx=12, nprocsy=16

! pochva --->
    integer, parameter :: nlev_snow=11
! <---

    real    :: dtstep=30., dtstepa, ddamp=0.15, hrun=24., hist=3., hist_full_res=24., &
               hbound=1., hdiag=0.25, srad=480., htop=0.92, mswshf=60., mslfilter=0.5, &
               hrst=12., hback=-0.1, hspray=-0.5, xsorg=10., ysorg=45., nudging_time=14400.
    integer :: nbl=8, nbc, nradm=2, nstep, nstep0, nadv=3, nsound=4, nhist, ntsbou, ndrunt, &
               ntsrad, ntop, ntswshf, mhfr=1, nstep_sl_filter, &
               nhist_full_res, ntnudg, nrst, nxspray=70, nyspray=70, nzspray=50
    logical :: nlbfix=.false., nlmic2=.false., nlana=.true., nlradar=.false., nlclimate=.false., nlrst=.false., &
 nltwice=.false., nlhdiff=.false., nlconv=.false., nlsnudg=.false., nlord=.false.
    namelist /model/ dtstep, nadv, nsound, nbl, nradm, hrun, hist, hist_full_res, hbound, hdiag, &
                     mhfr, srad, mslfilter, htop, &
                     nlmic2, nlbfix, nlana, &
                     mswshf, ddamp, nlbfix, nlana, nlord, &
                     nlconv, nsnudg, nlradar, nltwice, nlhdiff, nlclimate, nlrst, &
                     hrst, hback, hspray, nxspray, nyspray, nzspray, xsorg, ysorg

    integer, parameter :: nlon = (gnlon-2)/nprocsx+2, nlat = (gnlat-2)/nprocsy+2, &
                          km=gnlon/40, lm=gnlat/40
    integer, parameter :: nlonm1 = nlon-1, nlatm1 = nlat-1, nlevp1 = nlev+1, nbuffer = 2*max(nlon,nlat)*nlev
    integer, parameter :: ntype_aerosol = 6 ! defined in subroutine aerdef, must be eqaul to kaer variable in subroutine radintec,
! or in subroutine radiat_init (if ecmwf radiation scheme is not available)
!    integer, parameter :: nst_old = 15, nvt_old = 14  !  nst # of soil types; nvt # of vegetation types.
    integer, parameter :: nst = 31, nvt = 22  !  nst # of soil types; nvt # of vegetation types.
    integer            :: nyrin, nmonin, ndayin, nhouin, nminin
    integer, parameter :: nlonr  = nlon/2-1
    real               :: dlon, dlat, alon0, alat0, x0, y0, h, a0, b0, dx, dy, dz, dt, co2ppm, rdecli, reqtim
    real,    parameter :: yliv = 2834170.5, yliw = 333560.5, ylwv = yliv-yliw, tzer = 273.15, pzer = 1.e5, ezer = 611.,  &
                          rd = 287.05, rv = 461.51, eps = rd/rv, cpd = 1004.6, cvd = cpd-rd, cpv = 1869.46, cw = 4186.8, &
                          ci = cw/2., gamma = cpd/cvd, rdrcp = rd/cpd, a = 6371.e3, g = 9.807, omega = 7.292e-5,         &
                          pi = 3.14159265, tkemin = 1.e-9, qccrit = 2.2e-4, cvv = cpv-rv, rdrcv = rd/cvd
    real,    parameter :: bcldwat = 0.40e-3, bcldice = 0.08e-3  ! limits of cloud water and ice at lateral boundaries
!!!    real,    parameter :: alsn    = .71    ! mean snow albedo (except over glaciers)
    real,    parameter :: radarmval = 0.   ! radar reflectivity (dbz) min. value (from -200 to 0)
    integer, dimension(50)        :: nfdr, nfdr0, nfdrr
    real, dimension(100)          :: pdr, pdr0, pdrr
    real, dimension(256)          :: bndrel
    real, dimension(nlev)         :: ffilt
    real, dimension(nlevg)        :: d, lev_soil
    real, dimension(nlat)         :: fmyu, fmyv, clv, hxt, tang
    real, dimension(nlon,nlat)    :: alont, alatt, fcorio, hx, hy, ps, fmask, phig, orogvar, roscdm, &
                                     roscdt, rgm, rgmd, rgq, albedo, emisg1, emisg2, alsn, &
                                     prectot, precconv, precsolid, raini, snowi, raicon, snocon, &
                                     snow, fsnow, bvf, tskin, qskin, &
                                     cloudt, hflux, qflux, frirr, frvis, swsdtf, swsddr, &
                                     fveg, lai, proot, veg_rough, &
                                     deriv, qveg, fwetl, fvegs, cw1, cw2, cw3, ustar, tstar, qstar, cla1,  &
                                     runoff, runoff_tot, dfrvis, dfrirr, slopeff, snowfor,                      &
                                     fice, fice_rd, iceth, iceth_rd, zdtg, lapse_rate,                          &
                                     hsnc, fsnowmax, tetavs, &
                                     fices, pbl, gust, &
                                     albedo_soil_dry, albedo_soil_wet, albedo_veg,                              &
                                     emiss1_soil_dry, emiss1_soil_wet, emiss1_veg,                              &
                                     emiss2_soil_dry, emiss2_soil_wet, emiss2_veg
    real, dimension(nlon,nlat,nst+1)  :: suolo
    real, dimension(nlon,nlat,nvt+1)  :: vegeta
    real, dimension(gnlon,gnlat)      :: gfield, gfield1
    real, dimension(nlon,nlat,nlevg)  :: tg, qg, qgice
    real, dimension(nlon,nlat,nlev)   :: tvirt, fmz, u, v, t, p, q, qcw, qci, qpw, qpi1, qpi2, div2,          &
                                         ub1, vb1, tb1, pb1, qb1, qcwb1, qcib1, ub2, vb2, tb2, pb2, pf,       &
                                         qb2, qcwb2, qcib2, rradar, dtdt, ut, vt, wt, tket, ncw, nci, fcloud, &
                                         zeta, pai, tetav, qsatw, chm, dqdt, dqcwdt, dqcidt, dqpwdt, dqpi1dt
    real, dimension(nlon,nlat,nlevp1) :: w, s, wb1, wb2, fmzh, wwkw, tke, cvm, mlz, prandtl, deltaw, rich, tetavh
    integer, dimension(nlon,nlat) :: level_snowfall

! pochva --->

    integer, dimension(nlon,nlat) :: mask_soil, ind_lev_soil_h_bottom, ind_lev_soil_w_bottom
    real, dimension(nlon,nlat)    :: kturb_surf_m, kturb_surf_h, kturb_surf_q,  flh_specif, flh_latent, fl_wv, &
                                     h_lev_bottom, fl_rain, fl_snow, fl_rad_tot, &
                                     fl_heat_soil_bottom, fl_water_soil_bottom, fl_runoff, fl_runoff_tot, &
                                     cwvflux, cfl_heat_soil_bottom, cfl_water_soil_bottom, &
                                     tg_surf, qg_surf, fice_soil_surf, snow_dirt, &
                                     tgclim, qgclim, &
                                     water_table_depth, tg_bottom, qg_rel_bottom, qg_rel_surf_approx
    real, dimension(nlon,nlat,20)  :: kturb_surf_m_mem, kturb_surf_h_mem, kturb_surf_q_mem
    real, dimension(nlon,nlat,nlev_snow) :: snow_lev, snow_t, snow_fice, snow_age, snow_melt_age, snow_dens
    real(4), dimension(nlon,nlat,nlevg) :: qsoil_max, qsoil_min, c_soil, rho_soil, &
                                           psi_soil, kw_soil, par_b_soil, par_c_soil, &
                                           qsrel_wilt, qsrel_ref, &
                                           fice_soil
    real, parameter :: valmiss=-9999.
    real :: lai_max
    integer :: flag_pochva_initial=1, flag_pochva_par_change=0, flag_surf_flux=0

! if flag_surf_flux = 1, then surface turbulent heat and water vapor fluxes are defined by surface layer scheme,
! if flag_surf_flux = 0, then surface turbulent heat and water vapor fluxes are defined by soil scheme.
! Suggesting: use flag_surf_flux=0 for time step over 45 s for numerical stability

! <---

! Surface cumulated fields

    real, dimension(nlon,nlat) :: cswfl, clwfl, cshflux, clhflux, t2min, t2max, ws10max, totsky, soldir, &
                                  shf_accum, lhf_accum, qf_accum

! Additional surface fields in shf

    real, dimension(nlon,nlat) :: qprectot=0., qprecsolid=0.

! For "interpolation" of wind, temp. and humidity at prescribed levels of geometric height (m)
! above and/or below the earth surface

    integer, parameter :: n_std_lev_atm = 5, n_std_lev_soil = 0
    integer, dimension(nlon,nlat) :: n_std_lev_sl
    real(4), dimension(n_std_lev_atm) :: std_lev_atm=(/2., 10., 50., 80., 100./)
    real(4), dimension(n_std_lev_soil) :: std_lev_soil
    real(4), dimension(nlon,nlat,n_std_lev_atm) :: t_std_lev, u_std_lev, v_std_lev, q_std_lev, rh_std_lev, td_std_lev
    real(4), dimension(nlon,nlat,n_std_lev_soil) :: tg_std_lev, qg_std_lev

! costants used in computing partial pressures at saturation

    real, parameter :: ccw1=(cpv-cw)/rv, ccw2=ylwv/tzer/rv-ccw1, cci1=(cpv-ci)/rv, cci2=yliv/tzer/rv-cci1

! mpi veriables

    integer :: myid=0, iprocs=1, infx=1, infy=1, supx, supy, &
               ip_e=-1, ip_n=-1, ip_s=-1, ip_w=-1, ip_null=-1, &
               comm_row=0, comm_col=0, row_color=0, col_color=0

! radiation schemes

    integer, parameter                :: mcica = 0
    real, dimension(nlon,nlat,nlev,ntype_aerosol) :: aerosol
    real, dimension(nlon,nlat)        :: corvis, corirr, gelvis, gelirr
    real, dimension(nlon,nlat,nlev)   :: corrdt, geldt, ozon, aerotot
    integer :: imhf=0, ishf=0, irf=0, inst_spray=0

    contains

    function rdeno (t1, t2, t3, t4)   ! in WAF scheme
    real rdeno, t1, t2, t3, t4, zzden
    zzden = t3-t4
    rdeno = (t1-t2)/sign(max(abs(zzden),1.e-15),zzden)
    end function rdeno

    function gzita (zita)    ! Decay function
    gzita = 1. -a0*(zita/h)-(3.-2.*a0)*(zita/h)**2 +(2.-a0)*(zita/h)**3
    end function gzita

    function gzitap (zita)   ! Derivative of decay function
    gzitap = (-a0-(6.-4.*a0)*(zita/h) +(6.-3.*a0)*(zita/h)**2)/h
    end function gzitap

    function bzita (zita)    ! Stretching function
    bzita = b0 + (1.-b0)*(zita/h)
    end function bzita

    function bzitap (zita)   ! Derivative of stretching function
    bzitap = (1.-b0)/h
    end function bzitap

    end module mod_moloch
!###############################################################################################################
    program moloch

    use mod_moloch
    use pochva_scheme, only : pochva

    implicit none

#ifdef mpi
      include 'mpif.h'
#endif

    integer      :: ierr, comm ! mpi
    integer      :: iunit_inp=10, ndim, i, j, k, n, jstep, jsound, interv, &
                    jklevp1, jklevm1, jlatm1, jlonp1, jlo1, jlo2, &
                    jla1, jlon, jlat, jklev, iday, ihou, imin, &
                    jlatp1, jkl1, ndayr, nyrc, nmonc,  ndayc, nhouc, nminc, nyrc1, &
                    nmonc1, ndayc1, nhouc1, nminc1, ndayr1, jsday, jadv
    integer      :: iday0, ihou0, imin0
    real         :: zt0, zeps, zita, zitah, zfz, zfzh, time, ztime, zmax, zdummy, zfilt, zgamq,    &
                    zdtrdx, zdtrdy, zdtrdz, zcs2, zrrcv, zuh, zvh, zrom1w,                         &
                    zwexpl, zp, zm, zrapp, zdiv, zcx, zcu, zcy, zcyp, zcv, zcw, zrom1u, zrom1v,    &
                    zucor, zvcor, zu, zv, zw, zum, zup, zvm, zvp, zwm, zwp, zrdx, zrdy, zsig,      &
                    zz1, zdgz, z1, z2, zub, zvb, zwb, zpb, ztb, zqb, zqcwb, zqcib, zqgleq, zws10,  &
                    zalsn, zrrcv1, zhea, zrfmzu, zrfmzv, zrfmzum, zrfmzvp, zrid, snwe
    real         :: zqs, zdth, zt0t, zesk
    real*8       :: zfac, zx0, zy0, zlont, zlatv, zlatt, zaarg, zargt, zzlatt
    integer*4    :: stime1, stime2, countrate, countmax
    integer imonth(12)
!    character(len=15) :: filesoil
    character(len=40) :: command
    character(len=30) :: str
    real         :: pd1, pd2, pd4, pd5, pd38, pd39
    integer      :: gnlonr, gnlatr, ntback, ntspray, i1spray, j1spray
! pochva --->
    real         :: zdlev, dlev_base=10., dlev_min=5.
    integer :: n_outpoint=0
    integer, dimension(10) :: &
 i_outpoint=    (/ 70,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1/), &
 j_outpoint=    (/ 48,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1/), &
 iproc_outpoint=(/280,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1/)
    character (len=30) :: file_outpoint
! <---

    data imonth/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/

    call system_clock (stime1, countrate, countmax)

!--------------------------------------------------------------------------------------------------------
! initialize the scalar case
!--------------------------------------------------------------------------------------------------------

    nstep0   = 1
    supx    = nlon
    supy    = nlat

!--------------------------------------------------------------------------------------------------------
! initialize MPI environment
!--------------------------------------------------------------------------------------------------------

#ifdef mpi
      comm = mpi_comm_world
      call mpi_init(ierr)
      if (ierr .ne. 0) then
      print *,'Error starting MPI program. Terminating.'
      call mpi_abort(comm, ierr)
      endif

      call mpi_comm_size(comm, iprocs, ierr)  ! find out number of processors
#endif

    if (iprocs.ne.nprocsx*nprocsy) then
    print *,'Error in the definition of number of processes:'
    print *,'check values of nprocsx and nprocsy.',iprocs,nprocsx,nprocsy
    print *,'Terminating.'
#ifdef mpi
      call mpi_abort(comm, ierr)
#endif
    stop
    endif

#ifdef mpi
      if (gnlat-2.ne.(nlat-2)*nprocsy.or.gnlon-2.ne.(nlon-2)*nprocsx) then
        print *,'Errors in domain decomposition: terminating.'
        call mpi_abort(comm, ierr)
      endif
#endif

    if (mod(nlon,2).ne.0.or.mod(nlat,2).ne.0) then
    print *,'nlon or nlat are not even numbers: terminating.'  ! for ECMWF radiation computation
#ifdef mpi
      call mpi_abort(comm, ierr)
#endif
    endif

#ifdef mpi

      call mpi_comm_rank(comm, myid, ierr) ! find out my process id
      if (myid.eq.0) then
        print *
        print*,'#### MOLOCH parallel version ####'
        write(*,6000) nprocsx, nprocsy, gnlon, gnlat
 6000 format(/,' Number of processors =',i2,' x',i2,/,' Global x dim. =',i5,/,' Global y dim. =',i5)
      endif

! set exterior limits of the local domain

      infx = 1+(myid/nprocsy)*(nlon-2)
      supx = infx+nlon-1
      infy = 1+(myid-(myid/nprocsy)*nprocsy)*(nlat-2)
      supy = infy+nlat-1

      write (*,6001) myid, nlon, nlat, nlev, infx, supx, infy, supy
 6001 format(' myid=',i3,'     nlon,nlat,nlev=',3i4,'     position in global domain= ',4i4)

! define id's of neighbour processes

      ip_e = myid+nprocsy
      ip_w = myid-nprocsy
      if (ip_e.gt.iprocs-1) ip_e = mpi_proc_null
      if (ip_w.lt.0       ) ip_w = mpi_proc_null
      ip_s = myid-1
      ip_n = myid+1
      if (mod(myid,nprocsy).eq.0        ) ip_s = mpi_proc_null
      if (mod(myid,nprocsy).eq.nprocsy-1) ip_n = mpi_proc_null
      ip_null = mpi_proc_null

! define communicators of rows and columns of processes

      row_color = myid - (myid/nprocsy)*nprocsy
      col_color = myid/nprocsy
      call mpi_comm_split (comm, row_color, iprocs, comm_row, ierr)
      call mpi_comm_split (comm, col_color, iprocs, comm_col, ierr)

!      call system('hostname')
#endif

!--------------------------------------------------------------------------------------------------------
    if (myid == 0) then
      print *
      print *,'    --- Moloch model Version 23.0.1 ---'
      print *
    endif
!--------------------------------------------------------------------------------------------------------
!  Read parameters from namelist file
!--------------------------------------------------------------------------------------------------------

    open (iunit_inp, file='moloch.inp', status='old')
    read (iunit_inp, model)
    close (iunit_inp)

#ifndef rad_ecmwf
    nradm = min(nradm, 1) ! No ECMWF radiation scheme
#endif

    if (myid.eq.0) print *,'Parameters of moloch (moloch.inp):'
    if (myid.eq.0) print model

    nstep   = hrun *3600./dtstep+.5
    nhist   = hist *3600./dtstep+.5
    nhist_full_res   = hist_full_res *3600./dtstep+.5
    ntsbou  = hbound *3600./dtstep+.5
    nbc     = (nstep-1)/ntsbou+2
    ndrunt  = max( int(hdiag*3600./dtstep+.5), 1)
    ntsrad  = srad/dtstep+.5
    ntnudg  = ntsrad
    ntop    = htop*nlev
    ntswshf = int(mswshf*60./dtstep+.5) ! if negative, shf not written
    dtstepa = dtstep/float(nadv)
    nstep_sl_filter = int(mslfilter*60./dtstep)
    nrst    = hrst*3600./dtstep+.5
    ntback  = hback *3600./dtstep+.5  ! if negative, traj. file not written
    ntspray = hspray*3600./dtstep+.5  ! if negative, SPRAY file not written

    if(mhfr.eq.0.or.mhfr.gt.2) then
    if (myid.eq.0) write(*,'(a,i2)') " Parameter mhfr invalid or not defined in moloch.inp: stop!"
#ifdef mpi
      call mpi_abort(comm, ierr)
#endif
    stop
    endif

    if (myid.eq.0) then
      write(*,'(a,i4)') " Number of levels =", nlev
      print*
      if (mhfr.eq.1) then
        write(*,'(a,i2)') " A full mhf file will be written, mhfr =", mhfr
      elseif (mhfr.eq.2) then
        write(*,'(a,i2)') " A reduced mhf file will be written, mhfr =", mhfr
      endif
      if (nlconv) then
        print*
        print*, '****** Parameterization of convection activated ******'
        print*
      endif
    endif

!--------------------------------------------------------------------------------------------------------
!  Read initial condition
!--------------------------------------------------------------------------------------------------------

    if (nlrst) then
      call rdrf (1) ! for nstep0 definition only
      interv = nstep0/ntsbou + 1
    else
      interv = 1
    endif

    call rdmhf_atm (interv, 1)
    call rdmhf_soil (interv, 1)
    call rd_param_const
    if (mhfr /= 1) call wr_param_const

    if (nlrst) then
      call rdrf (1) ! for nfdr0, pdr0, nfdr, pdr definition only
    endif

!--------------------------------------------------------------------------------------------------------
! read soil and vegetation parameters
!--------------------------------------------------------------------------------------------------------

! Soil, vegetation, snow version 19 and less --->
!    call rdgeo
! <---

!--------------------------------------------------------------------------------------------------------
!  Constants
!--------------------------------------------------------------------------------------------------------

    dlat   = pdr0(1)
    dlon   = pdr0(2)
    alat0  = pdr0(4) + (infy-1)*dlat  ! latit.  of the 1st v-point of the local domain
    alon0  = pdr0(5) + (infx-1)*dlon  ! longit. of the 1st t-point of the local domain
    lev_soil(1:nlevg) = pdr0(6:5+nlevg)
    d(1) = lev_soil(1)*2.
    do jklev = 2, nlevg
      d(jklev)=2.*(lev_soil(jklev)-lev_soil(jklev-1))-d(jklev-1)
    enddo
    y0     = pdr0(38)
    x0     = pdr0(39)
    h      = pdr0(40)      ! h = rd*t0/g
    b0     = pdr0(42)
    a0     = pdr0(43)
    zt0    = g*h/rd
    dx     = a*dlon*pi/180.
    dy     = a*dlat*pi/180.
    dz     = h/nlev
    zeps   = 1./eps

    if (ntspray.gt.0) then            ! punto iniziale del ritaglio per SPRAY
!    nxspray=330; nyspray=150; nzspray=45; i1spray=245+1; j1spray=780+1  ! pimpa su operativo
    call raglio (x0,y0, pdr0(5),pdr0(4), dlon,dlat, xsorg,ysorg, gnlon,gnlat, nxspray,nyspray, i1spray,j1spray)
    if (myid.eq.0) print*, 'SW origin of SPRAY file:', i1spray, j1spray
    endif

!--------------------------------------------------------------------------------------------------------
!  Initial date
!--------------------------------------------------------------------------------------------------------

! save old validity date

    iday0  = nfdr0(10)
    ihou0  = nfdr0(11)
    imin0  = nfdr0(12)

! decode year/month/day/hour/min. of initial condition

    nyrin = nfdr0(5)
    if (mod(nyrin,4).eq.0) imonth(2)=29
    nmonin = nfdr0(6)
    ndayin = nfdr0(7)
    nhouin = nfdr0(8)
    nminin = nfdr0(9)

! update initial condition date (if initial condition is a forecast)

      nminin = nminin + imin0
      if (nminin.ge.60) then
      nminin = nminin - 60
      nhouin = nhouin + 1
      endif
      nhouin = nhouin + ihou0
      if (nhouin.ge.24) then
      nhouin = nhouin - 24
      ndayin = ndayin + 1
      endif
      ndayin = ndayin + iday0
      if (ndayin.gt.imonth(nmonin)) then
      ndayin=ndayin-imonth(nmonin)
      nmonin=nmonin+1
      endif
      if (nmonin.eq.13) then
      nmonin=1
      nyrin=nyrin+1
      endif
      if (nyrin.eq.100) nyrin=0


    if (myid.eq.0) write(*,6003) nyrin, nmonin, ndayin, nhouin, nminin
 6003 format(/,' Date of initial condition  -- YYYY/MM/DD/HH/MM --', I6.2,4I4.2)

! define day of the year ( 1 < ndayr < 366 )

    ndayr = ndayin
    do j=1,nmonin-1
    ndayr = ndayr + imonth(j)
    enddo

! CO2 concentration defined as a function of the year, assuming linear trend around 2015

    co2ppm = 400. +2.*(nyrin - 2015)
    co2ppm = max (co2ppm, 280.)
    if(myid.eq.0) write(*, '(a, f7.2)') ' CO2 concentration in ppm =', co2ppm

!--------------------------------------------------------------------------------------------------------
!  Metrics - Rotated grid - Coriolis - Slope
!--------------------------------------------------------------------------------------------------------

    do jlat=1,nlat
    do jlon=1,nlon
    do jklev=1,nlev
    zita =(jklev-1)*dz+dz/2.
    zitah=(jklev-1)*dz
    zfz  = 1.-zita/h
    fmz (jlon,jlat,jklev) = zfz/(bzita(zita )+phig(jlon,jlat)/g*zfz*gzitap(zita ) &
                          - h*zfz*log(zfz)*bzitap(zita ))
    zeta(jlon,jlat,jklev) = max(phig(jlon,jlat)/g*gzita(zita) -h*bzita(zita)*log(zfz) &
                          - phig(jlon,jlat)/g, 0.) ! height above orography
    zfz  = 1.-zitah/h
    fmzh(jlon,jlat,jklev) = zfz/(bzita(zitah)+phig(jlon,jlat)/g*zfz*gzitap(zitah) &
                          - h*zfz*log(zfz)*bzitap(zitah))
    enddo
    fmzh(jlon,jlat,nlevp1) = 0.

    do jklev = 1, n_std_lev_atm
    if (std_lev_atm(jklev) > zeta(jlon,jlat,1)) then
    n_std_lev_sl(jlon,jlat) = jklev-1
    exit
    endif
    enddo

    enddo
    enddo

    zfac = dabs(dacos(-1.d0))/180.d0
    zx0  = dble(x0)*zfac
    zy0  = dble(y0)*zfac
    do jlat=1,nlat
    zlatv=(dble(alat0)+                (jlat-1)*dble(dlat))*zfac
    zlatt=(dble(alat0)+dble(dlat)/2.d0+(jlat-1)*dble(dlat))*zfac
    hxt (jlat) = dcos(zlatt)
    clv (jlat) = dcos(zlatv)
    fmyv(jlat) = 1.d0/dcos(zlatv)
    fmyu(jlat) = 1.d0/dcos(zlatt)
    tang(jlat) = sin(zlatt)*fmyu(jlat)/a
    do jlon=1,nlon
    zlont = (dble(alon0)+(jlon-1)*dble(dlon))*zfac
    zzlatt= 1.d0/zfac*dasin( dcos(zy0)*dsin(zlatt) + dsin(zy0)*dcos(zlatt)*dcos(zlont) )
    zargt = -dsin(zlatt)*dsin(zy0)+dcos(zy0)*dcos(zlatt)*dcos(zlont)
    zaarg = zargt/dcos(zfac*zzlatt)
    alatt(jlon,jlat)=zzlatt
    if(zaarg.lt.-1..and.zaarg.gt.-1.00001) zaarg = -1.d0
    if(zaarg.gt. 1..and.zaarg.lt. 1.00001) zaarg =  1.d0
      if (zlont.lt.0.d0) then
      alont(jlon,jlat) = 1.d0/zfac*(zx0-dacos(zaarg))
      else
      alont(jlon,jlat) = 1.d0/zfac*(zx0+dacos(zaarg))
      endif
    fcorio(jlon,jlat)=2.*omega*sin(pi*alatt(jlon,jlat)/180.)
    enddo
    enddo

!  orographic slope

    do jlat = 2, nlat
    do jlon = 1, nlonm1
    hx(jlon,jlat)=(phig(jlon+1,jlat)-phig(jlon,jlat))*fmyu(jlat)/g/dx
    hy(jlon,jlat)=(phig(jlon,jlat)-phig(jlon,jlat-1))/g/dy
    enddo
    enddo

!--------------------------------------------------------------------------------------------------------
!  Check Courant numbers
!--------------------------------------------------------------------------------------------------------

  if (myid.eq.0) then
    print*
    write(*,'(a, f8.2, a, f8.2, a, f7.2)') ' dx =', dx*cos(pi/180.*(alat0+dlat*gnlat/2)), ', dy =', dy, ', dzita =', dz
    print*
    write(*,'(a, f8.2)') ' Density scale H and zita top =', h
    print*
    write(*,'(a, f5.2)') ' Vertical coord. stretching factor b0 =', b0
    print*
    write(*,'(a, f5.2)') ' Parameter in vert. coord. decay function: a0 =', a0
    print*
    write(*,'(a, f7.4)') ' Courant number of horizontal sound waves =', sqrt(2.)*sqrt(gamma*rd*300.)*dtstepa/nsound/dx
    print*
  endif

    zmax=0.
    do jklev = 1, nlev
    do jlat = 2, nlat
    do jlon = 1, nlon
    zmax = max (zmax, sqrt(u(jlon,jlat,jklev)**2+v(jlon,jlat,jklev)**2))
    enddo
    enddo
    enddo

#ifdef mpi
      call mpi_reduce (zmax, zdummy, 1, mpi_real, mpi_max, 0, comm, ierr)    ! calculate global maximum
      zmax = zdummy
#endif
    if (myid.eq.0) then
      write(*,'(a, f7.4)') ' Max. Courant number for horizontal advection =', sqrt(2.)*zmax*dtstepa/dx
      print*
    endif

  if (myid.eq.0) then
    write(*,'(a,i4)') ' Top level of microphysics and turbulence ntop =', ntop
    print*
    print*,'               zeta       zetah'
    do jklev = nlev, 1, -1
    zitah=(jklev-1)*dz
    zita =(jklev-1)*dz+dz/2.
    if(jklev.le.9) write(*,6007) jklev, -h*bzita(zita)*log(1.-zita/h), -h*bzita(zitah)*log(1.-zitah/h)
    enddo
 6007 format(' Lev =', i3, 1x, 2f11.2)
    print*
    write(*,'(a,i2,a)') ' Radiation selection nradm =', nradm, ' (0: none; 1: Geleyn; 2: ECMWF)'
  endif

!--------------------------------------------------------------------------------------------------------
!  Initialization of soil and vegetation scheme
!--------------------------------------------------------------------------------------------------------

  if (myid.eq.0) then
    print *
    print *,'Depth of soil levels (m):'
    print '(30f6.2)', lev_soil(1:nlevg)
  endif

  if (.not.nlrst) then

    nstep0 = 1

    do jlat = 1, nlat
    do jlon = 1, nlon

    alsn(jlon,jlat) = .71      ! mean snow albedo (exept glaciers) will be changed by soil processes scheme

    fice(jlon,jlat) = max(fice(jlon,jlat), 0.)
    fice(jlon,jlat) = min(fice(jlon,jlat), 1.)
    if (fmask(jlon,jlat).lt.0.5) fice(jlon,jlat) = 0.

    qg(jlon,jlat,1:nlevg) = max( min( qg(jlon,jlat,1:nlevg), qsoil_max(jlon,jlat,1:nlevg)), qsoil_min(jlon,jlat,1:nlevg))

    rgmd(jlon,jlat)  = rgm(jlon,jlat)

! pochva --->

    if (fmask(jlon,jlat) < 0.5.or.fice(jlon,jlat) >= 0.8) then
      mask_soil(jlon,jlat) = 1
    else
      mask_soil(jlon,jlat) = 0
    endif

! Constant parameters that must be modified in the case of thick sea ice
! (all variable parameters of thick sea ice have been correctly defied by
! preprocessing procedure and read from input_soil_01.mhf file)

    if (fice(jlon,jlat) >= 0.8) then ! Thick sea ice
      psi_soil(jlon,jlat,:) = 2.2 ! thermal conductivity of ice
      do jklev = 1, nlevg
        if (lev_soil(jklev) > iceth(jlon,jlat)) exit
      enddo
      ind_lev_soil_h_bottom(jlon,jlat) = max(jklev-1, 3)
      ind_lev_soil_w_bottom(jlon,jlat) = 0
    endif

! <---

    enddo
    enddo

! ---> pochva
    mask_soil(1,   1:nlat) = 0
    mask_soil(nlon,1:nlat) = 0
    mask_soil(1:nlon,   1) = 0
    mask_soil(1:nlon,nlat) = 0
! <---

!  Albedo and emissivity

 if (myid == 0) then
   print *
   print *,"Renewal of radiation parameters of the surface"
 endif
 call surfradpar (fmask, fice, qg_surf, vegeta, qsoil_max, qsoil_min, &
 albedo_soil_dry, albedo_soil_wet, emiss1_soil_dry, emiss1_soil_wet, emiss2_soil_dry, emiss2_soil_wet, &
 fveg, albedo_veg, emiss1_veg, emiss2_veg, &
 nlon, nlat, nlevg, nvt, 1, albedo, emisg1, emisg2)

!--------------------------------------------------------------------------------------------------------
!  Dynamics workspace
!--------------------------------------------------------------------------------------------------------

    w(:,:,nlevp1) =0.
    wwkw(:,:,1)   =0.       ! tridiagonal inversion
    s             =0.       ! generalized vertical velocity
    deltaw        =0.       ! nonhydrostatic term in pressure gradient force
    div2          =0.

  endif ! .not.nlrst

    dt     = dtstepa/float(nsound)
    zdtrdx = dt/dx
    zdtrdy = dt/dy
    zdtrdz = dt/dz
    zcs2   = zdtrdz**2*rdrcv

!--------------------------------------------------------------------------------------------------------
!  Sponge layer at the top of the atmosphere
!--------------------------------------------------------------------------------------------------------

    do jklev = 2, nlev
      if (jklev.le.ntop-2) then  !  filter of w at top
      ffilt(jklev) = 0.
      else
      zfilt = (ntop-3)*dz
      zz1 = (dz*(jklev-1)-zfilt)/(h-zfilt)
      ffilt(jklev) = 1.*sin(0.5*pi*zz1)**2
      endif
    enddo

!--------------------------------------------------------------------------------------------------------
!  Initialization of physical parameters
!--------------------------------------------------------------------------------------------------------

    dtdt      = 0.
    geldt     = 0.
    corrdt    = 0.
    frvis     = 0.
    frirr     = 0.
    gelvis    = 0.
    gelirr    = 0.
    corvis    = 0.
    corirr    = 0.
    fcloud    = 1.e-10     ! box fraction covered by cloud
    qpw       = 0.
    qpi1      = 0.
    qpi2      = 0.
    ncw       = 3.e7       ! density numbers of water clouds (used for 2 moment microph.)
    nci       = 6.5e7      ! density numbers of ice clouds (used for 2 moment microph.)
    tke       = tkemin     ! turbulent kinetic energy
    mlz       = 1.         ! mixing length
    prectot   = 0.         ! sum of liquid and soil precipitation (kg/m2) accum. between 2 checkpoints
    precconv  = 0.         ! convective rain+snow (kg/m2) accum. between 2 checkpoints
    precsolid = 0.         ! solid precipitation (kg/m2 of equiv. water) acc. between 2 checkpoints
    runoff    = 0.         ! runoff of upper soil layer (kg/m2) between 2 checkpoints
    runoff_tot= 0.         ! total runoff (kg/m2) between 2 checkpoints
    raini     = 0.         ! rain in one timestep (kg/m2)
    snowi     = 0.         ! snow in one timestep (kg/m2 of equiv. water)
    hflux     = 0.
    qflux     = 0.
    rradar    = radarmval  ! background radar reflectivity (dbz)
    roscdm    = 1.e-2
    roscdt    = 1.e-2      ! 'drag coefficient' at the surface (also used in soil scheme)
    cswfl     = 0.
    clwfl     = 0.
    cshflux   = 0.
    clhflux   = 0.
    t2min     = 999.
    t2max     = 0.
    ws10max   = 0.
    totsky    = 0.
    soldir    = 0.
    rich      = 10.
    shf_accum = 0.
    lhf_accum = 0.
    qf_accum  = 0.
    slopeff   = 0.
    cwvflux   = 0.
    cfl_heat_soil_bottom  = 0.
    cfl_water_soil_bottom = 0.

!--------------------------------------------------------------------------------------------------------
! Initialization of basic pressure (for dynamic) variable (pai)
!--------------------------------------------------------------------------------------------------------

    call paidef (p, t, q, qcw, qci, pai)

!--------------------------------------------------------------------------------------------------------
!  Initialization of boundary conditions
!  Definition of boundary relaxation coefficients  - relax requires nbl to be a power of 2
!--------------------------------------------------------------------------------------------------------

    call relax (nbl, .01, 1., bndrel)
    do jlon = nbl,2,-1
    bndrel(jlon) = bndrel(jlon-1)
    enddo
    bndrel(1) = 1.

! Boundary condition at time t1 (start in time interval) read from mhf files

    pb1   = pai
    ub1   = u
    vb1   = v
    wb1   = w
    tb1   = t
    qb1   = q
    qcwb1 = qcw
    qcib1 = qci

    if (.not.nlbfix) then

! Boundary condition at time t2 (end of time interval) read from mhf files

      if (nlrst) then
        interv = nstep0/ntsbou + 2
      else
        interv = 2
      endif
      call rdmhf_atm (interv, 0)
      call paidef (pb2, tb2, qb2, qcwb2, qcib2, pf)
      pb2 = pf

    endif

!--------------------------------------------------------------------------------------------------------
!  Initialization of all prognostic and diagnostic variables using restart data files
!--------------------------------------------------------------------------------------------------------

  if (nlrst) then

! Reading of restart file

    call rdrf (0)
    call rdrf_pochva(nstep0, nprocsx, nprocsy, myid, gfield, gnlon, gnlat, dtstep)
    flag_pochva_initial = 0

  endif ! nlrst

!--------------------------------------------------------------------------------------------------------
!  Initialization of descriptor records
!--------------------------------------------------------------------------------------------------------

  if (.not.nlrst) then

!  define descriptor records

    nfdr  = nfdr0
    pdr   = pdr0

!  nfdr

    nfdr(5) = nyrin
    nfdr(6) = nmonin
    nfdr(7) = ndayin
    nfdr(8) = nhouin
    nfdr(9) = nminin
    nfdr(10:12) = 0
    nfdr(16)= nstep
    nfdr(14)= nbl
    nfdr(13)= nsound
    nfdr(17)= nhist
    nfdr(18)= ntswshf
    nfdr(20)= 1
    if (mhfr == 2) nfdr(20) = mhfr

! pdr

    pdr(3)  = dtstep
!    do jklev = 1, nlevg
!    pdr(5+jklev) = lev_soil(jklev)
!    enddo
    pdr(41) = qccrit ! threshold for cloud definition

  endif

! nfdrr and pdrr for output mhfr only (case of mhfr equal 1 or 2)

    nfdrr = nfdr
    nfdrr(2) = nfdr(2)/float(mhfr)
    nfdrr(3) = nfdr(3)/float(mhfr)
    pdrr = pdr
    pdrr(1) = float(mhfr)*pdr(1)
    pdrr(2) = float(mhfr)*pdr(2)

!--------------------------------------------------------------------------------------------------------
!  Initial condition on MHF
!--------------------------------------------------------------------------------------------------------

    if (.not.nlrst) then

      call cloudfr
      call ccloudt
      call snowfall_level

      if (nlana.and.myid == 0) imhf=imhf+1

      if (nlana) then

        call wrmhf_atm
        call wrmhf_soil
        if (myid == 0) write(*,'(a)') ' Initial condition saved in the moloch mhf_atm and mhf_soil files'

      endif

      call runout (0)

    else ! nlrst

      if (myid == 0) print *,'----------------- Restart --------------------'
      call runout (nstep0)
      if (myid == 0) print *,'----------------------------------------------'

      imhf = (nstep0-1)/nhist
      if (nlana) imhf = imhf+1
      ishf = (nstep0-1)/ntswshf+1
      irf = nstep0/nrst

    endif ! nlrst

!********************************************************************************************************
!  START INTEGRATION
!********************************************************************************************************

    call system_clock (stime1, countrate)  !.....elapsed time

    do 2000 jstep = nstep0, nstep

!  compute forecast time in days, hours, and minutes (iday, ihou, imin)

    iday = int( jstep*dtstep/86400. )
    ihou = int((jstep*dtstep-iday*86400.)/3600.)
    imin = int((jstep*dtstep-iday*86400.-ihou*3600.)/60. +.5)
    if(imin.ge.60) then
    imin = imin - 60
    ihou = ihou + 1
    endif
    if(ihou.ge.24) then
    ihou = ihou -24
    iday = iday +1
    endif

!  compute date of current time step (nyrc, nmonc, ndayc, nhouc, nminc) and day-of-the-year (ndayr)

    call calendar (nyrin, nmonin, ndayin, nhouin, nminin, iday, ihou, imin, &
                   nyrc, nmonc, ndayc, nhouc, nminc, ndayr)

    tetav = t*(1. -(q+qcw+qci+qpw+qpi1+qpi2) +zeps*q)/pai

    do 501 jadv = 1, nadv
    pf = pai

#ifdef mpi
      call u_ghost (tetav(:,nlatm1,:), ip_n, tetav(:   ,1,:), ip_s, nlon*nlev)
      call u_ghost (tetav(2,:     ,:), ip_w, tetav(nlon,:,:), ip_e, nlat*nlev)
#endif

    do jklev = 2, nlev
    do jlat = 2, nlatm1
    do jlon = 2, nlonm1
    tetavh(jlon,jlat,jklev) = .5*(tetav(jlon,jlat,jklev)+tetav(jlon,jlat,jklev-1))
    enddo
    enddo
    enddo

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    do 500 jsound = 1, nsound   ! sound waves
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! partial definition of the generalized vertical velocity

#ifdef mpi
      call u_ghost (v(:,2,:), ip_s, v(:,nlat,:), ip_n, nlon*nlev)
      call u_ghost (u(nlonm1,:,:), ip_e, u(1,:,:), ip_w, nlat*nlev)
#endif

    do jlat = 2, nlatm1
    do jlon = 2, nlonm1
    zuh = u(jlon,jlat,1)*hx(jlon,jlat)+u(jlon-1,jlat,1)*hx(jlon-1,jlat)
    zvh = v(jlon,jlat,1)*hy(jlon,jlat)+v(jlon,jlat+1,1)*hy(jlon,jlat+1)
    w(jlon,jlat,1) = .5*(zuh+zvh)
    s(jlon,jlat,1) = -w(jlon,jlat,1)
    enddo
    enddo

    do jklev = 2, nlev
    zitah = (jklev-1)*dz
    do jlat = 2, nlatm1
    do jlon = 2, nlonm1
    zuh = (u(jlon,jlat,jklev-1)+u(jlon,jlat,jklev))*hx(jlon,jlat)+   &
          (u(jlon-1,jlat,jklev-1)+u(jlon-1,jlat,jklev))*hx(jlon-1,jlat)
    zvh = (v(jlon,jlat,jklev-1)+v(jlon,jlat,jklev))*hy(jlon,jlat)+   &
          (v(jlon,jlat+1,jklev-1)+v(jlon,jlat+1,jklev))*hy(jlon,jlat+1)
    s(jlon,jlat,jklev) = -.25*(zuh+zvh)*gzita(zitah)
    enddo
    enddo
    enddo

! part of divergence (except w contribution) put in div2

    do jklev = 1, nlev
    do jlat = 2, nlatm1
    zcx = zdtrdx*fmyu(jlat)
    zcyp= zdtrdy*fmyu(jlat)*clv(jlat+1)
    zcy = zdtrdy*fmyu(jlat)*clv(jlat  )
    do jlon = 2, nlonm1
    zrfmzu  = 2./(fmz(jlon,jlat,jklev)+fmz(jlon+1,jlat,jklev))
    zrfmzum = 2./(fmz(jlon,jlat,jklev)+fmz(jlon-1,jlat,jklev))
    zrfmzv  = 2./(fmz(jlon,jlat,jklev)+fmz(jlon,jlat-1,jklev))
    zrfmzvp = 2./(fmz(jlon,jlat,jklev)+fmz(jlon,jlat+1,jklev))
    zup = u(jlon  ,jlat,jklev)*zrfmzu
    zum = u(jlon-1,jlat,jklev)*zrfmzum
    zvp = v(jlon,jlat+1,jklev)*zrfmzvp
    zvm = v(jlon,jlat  ,jklev)*zrfmzv
    div2(jlon,jlat,jklev) = ((zup-zum)*zcx +zvp*zcyp -zvm*zcy) * fmz(jlon,jlat,jklev)
    enddo
    enddo
    enddo

    call divdamp (ddamp)  ! divergence damping

    do jklev = 1, nlev
    do jlat = 2, nlatm1
    do jlon = 2, nlonm1
    div2(jlon,jlat,jklev) = div2(jlon,jlat,jklev)+zdtrdz*fmz(jlon,jlat,jklev)*(s(jlon,jlat,jklev+1)-s(jlon,jlat,jklev))
    enddo
    enddo
    enddo

! new w (implicit scheme)

    do jklev = 2, nlev
    do jlat = 2, nlatm1
    do jlon = 2, nlonm1

    if (.not.nlconv) deltaw(jlon,jlat,jklev) = -w(jlon,jlat,jklev)

! explicit w: it must be consistent with the initialization of pai

    tetavh(jlon,jlat,jklev) = tetavh(jlon,jlat,jklev)-w(jlon,jlat,jklev)*fmzh(jlon,jlat,jklev)*zdtrdz &
                            * (tetav(jlon,jlat,jklev)-tetav(jlon,jlat,jklev-1))
    zrom1w = cpd*tetavh(jlon,jlat,jklev)*fmzh(jlon,jlat,jklev)

    zwexpl = w(jlon,jlat,jklev) - zrom1w*zdtrdz*(pai(jlon,jlat,jklev)-pai(jlon,jlat,jklev-1)) - g*dt
    zwexpl = zwexpl + rdrcv*zrom1w*zdtrdz* &
           (pai(jlon,jlat,jklev)*div2(jlon,jlat,jklev)-pai(jlon,jlat,jklev-1)*div2(jlon,jlat,jklev-1))

! computation of the tridiagonal matrix coefficients
! - zp*w(k+1) + (1+zp+zm)*w(k) - zm*w(k-1) = zwexpl

    zp = zcs2*fmz(jlon,jlat,jklev  )*zrom1w*pai(jlon,jlat,jklev  ) + ffilt(jklev)
    zm = zcs2*fmz(jlon,jlat,jklev-1)*zrom1w*pai(jlon,jlat,jklev-1) + ffilt(jklev)

! 1st loop for the tridiagonal inversion

    zrapp = 1./(1.+zm+zp-zm*wwkw(jlon,jlat,jklev-1))
    w(jlon,jlat,jklev)    = zrapp*(zwexpl+zm*w(jlon,jlat,jklev-1))
    wwkw(jlon,jlat,jklev) = zrapp*zp
    enddo
    enddo
    enddo

! 2nd loop for the tridiagonal inversion

    do jklev = nlev, 2, -1
    do jlat = 2, nlatm1
    do jlon = 2, nlonm1
    w(jlon,jlat,jklev) = w(jlon,jlat,jklev) + wwkw(jlon,jlat,jklev)*w(jlon,jlat,jklev+1)
    if (.not.nlconv) deltaw(jlon,jlat,jklev) = deltaw(jlon,jlat,jklev) + w(jlon,jlat,jklev)
    enddo
    enddo
    enddo

! new Exner function

    do jklev = 1, nlev
    do jlat = 2, nlatm1
    do jlon = 2, nlonm1
    zdiv = div2(jlon,jlat,jklev)+zdtrdz*fmz(jlon,jlat,jklev)*(w(jlon,jlat,jklev+1)-w(jlon,jlat,jklev))

! no approximation ...
!    zdiv =zdiv*(1.+.86*q(jlon,jlat,jklev)+3.2*qcw(jlon,jlat,jklev))/(1.+.96*q(jlon,jlat,jklev)+4.8*qcw(jlon,jlat,jklev))
!    tetav(jlon,jlat,jklev) = tetav(jlon,jlat,jklev) *(1. +rdrcv*zdiv*(.25*q(jlon,jlat,jklev) +4.2*  &
!                             (qcw(jlon,jlat,jklev)+qpw(jlon,jlat,jklev)) +2.1*qci(jlon,jlat,jklev)) )

    pai(jlon,jlat,jklev) = pai(jlon,jlat,jklev)*(1.-rdrcv*zdiv)
    enddo
    enddo
    enddo

!#ifdef mpi
!      call u_ghost (tetav(:,nlatm1,:), ip_n, tetav(:   ,1,:), ip_s, nlon*nlev)
!      call u_ghost (tetav(2,:     ,:), ip_w, tetav(nlon,:,:), ip_e, nlat*nlev)
!#endif

! horizontal momentum equations

#ifdef mpi
      call u_ghost (pai(:,nlatm1,:), ip_n, pai(:   ,1,:), ip_s, nlon*nlev)
      call u_ghost (pai(2,:     ,:), ip_w, pai(nlon,:,:), ip_e, nlat*nlev)
#endif
    if (.not.nlconv) then
#ifdef mpi
        call u_ghost (deltaw(:,nlatm1,:), ip_n, deltaw(:   ,1,:), ip_s, nlon*nlevp1)
        call u_ghost (deltaw(2,:     ,:), ip_w, deltaw(nlon,:,:), ip_e, nlat*nlevp1)
#endif
    endif

    do jklev = 1, nlev
    zita = (jklev-1)*dz+dz/2.
    do jlat = 2, nlatm1
    zcx = zdtrdx*fmyu(jlat)
    do jlon = 2, nlonm1

    if (.not.nlconv) then
    zfz = .25*(deltaw(jlon  ,jlat,jklev+1) + deltaw(jlon  ,jlat,jklev)  &
             + deltaw(jlon+1,jlat,jklev+1) + deltaw(jlon+1,jlat,jklev)) + g*dt
    else
    zfz = g*dt
    endif

    zrom1u = .5*cpd*(tetav(jlon,jlat,jklev)+tetav(jlon+1,jlat,jklev))
    u(jlon,jlat,jklev) = u(jlon,jlat,jklev) - zrom1u*zcx*(pai(jlon+1,jlat,jklev)-pai(jlon,jlat,jklev)) &
                       - zfz*hx(jlon,jlat)*gzita(zita) + fcorio(jlon,jlat)*v(jlon,jlat,jklev)*dt

    if (.not.nlconv) then
    zfz = .25*(deltaw(jlon,jlat  ,jklev+1) + deltaw(jlon,jlat  ,jklev)  &
             + deltaw(jlon,jlat-1,jklev+1) + deltaw(jlon,jlat-1,jklev)) + g*dt
    else
    zfz = g*dt
    endif

    zrom1v = .5*cpd*(tetav(jlon,jlat,jklev)+tetav(jlon,jlat-1,jklev))
    v(jlon,jlat,jklev) = v(jlon,jlat,jklev) - zrom1v*zdtrdy*(pai(jlon,jlat,jklev)-pai(jlon,jlat-1,jklev)) &
                       - zfz*hy(jlon,jlat)*gzita(zita) - fcorio(jlon,jlat)*u(jlon,jlat,jklev)*dt
    enddo
    enddo
    enddo

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 500 continue  !  end of time step for sound waves
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 ! Complete computation of generalized vertical velocity

    do jklev = 2, nlev
    do jlat = 2, nlatm1
    do jlon = 2, nlonm1
    s(jlon,jlat,jklev) = (w(jlon,jlat,jklev) + s(jlon,jlat,jklev))*fmzh(jlon,jlat,jklev)
    enddo
    enddo
    enddo
    s(:,:,1) = 0.
    s(:,:,nlevp1) = 0.
    w(:,:,nlevp1) = w(:,:,nlev) !!!!!!!!!!!!!!

!--------------------------------------------------------------------------------------------------------
!  WAF advection of all variables
!--------------------------------------------------------------------------------------------------------

    call advect_waf (dtstepa)

!    pai = pai - pf
!    call filt3d (pai, .05)
!    pai = pai + pf

 501 continue

! reset of residual cloud and precip. up to NTOP and above NTOP

    do jklev = 1, ntop
    do jlat = 2, nlatm1
    do jlon = 2, nlonm1
    if (qcw(jlon,jlat,jklev).lt.1.e-9)  qcw (jlon,jlat,jklev)=0.
    if (qci(jlon,jlat,jklev).lt.1.e-9)  qci (jlon,jlat,jklev)=0. ! saving some pristine cristals...
    if (qpw(jlon,jlat,jklev).lt.1.e-8)  qpw (jlon,jlat,jklev)=0.
    if (qpi1(jlon,jlat,jklev).lt.1.e-8) qpi1(jlon,jlat,jklev)=0.
    if (qpi2(jlon,jlat,jklev).lt.1.e-8) qpi2(jlon,jlat,jklev)=0.
    enddo
    enddo
    enddo
    do jklev = ntop+1, nlev
    do jlat = 2, nlatm1
    do jlon = 2, nlonm1
    qcw (jlon,jlat,jklev) = 0.
    qci (jlon,jlat,jklev) = 0.
    qpw (jlon,jlat,jklev) = 0.
    qpi1(jlon,jlat,jklev) = 0.
    qpi2(jlon,jlat,jklev) = 0.
    enddo
    enddo
    enddo

!do jlat=1,nlat
!do jlon=1,nlon
!  if (abs(tskin(jlon,jlat))>400..or.abs(tskin(jlon,jlat))<100..or.isnan(tskin(jlon,jlat))) &
! print '(a,i5,e16.8,3x,3i4.3)',"*3 oshibka * tskin ",jstep,tskin(jlon,jlat),jlon,jlat,myid
!  if (abs(qskin(jlon,jlat))>10.) print *,"*3 oshibka * qskin ",jstep,qskin(jlon,jlat),jlon,jlat,myid
!  if (abs(snow(jlon,jlat))>1000.) print *,"*3* snow ",jstep,snow(jlon,jlat),jlon,jlat,myid
!  if (abs(fsnow(jlon,jlat))>10.) print *,"*3* fsnow ",jstep,fsnow(jlon,jlat),jlon,jlat,myid
!  if (isnan(tskin(jlon,jlat))) print *,"*3* oshibka tskin ",jstep,tskin(jlon,jlat),jlon,jlat,myid
!  if (isnan(qskin(jlon,jlat))) print *,"*3* oshibka qskin ",jstep,qskin(jlon,jlat),jlon,jlat,myid
!  if (isnan(snow(jlon,jlat))) print *,"*3* oshibka snow ",jstep,snow(jlon,jlat),jlon,jlat,myid
!  if (isnan(fsnow(jlon,jlat))) print *,"*3* oshibka fsnow ",jstep,fsnow(jlon,jlat),jlon,jlat,myid
!  if (isnan(v(jlon,jlat,1))) print *,"*3* oshibka v ",jstep,v(jlon,jlat,1),jlon,jlat,myid
!  if (isnan(t(jlon,jlat,1))) print *,"*3* oshibka t ",jstep,t(jlon,jlat,1),jlon,jlat,myid
!  if (isnan(q(jlon,jlat,1))) print *,"*3* oshibka q ",jstep,q(jlon,jlat,1),jlon,jlat,myid
!  if (isnan(tke(jlon,jlat,1))) print *,"*3* oshibka tke ",jstep,tke(jlon,jlat,1),jlon,jlat,myid
!enddo
!enddo

!--------------------------------------------------------------------------------------------------------
!  Turbulent diffusion
!--------------------------------------------------------------------------------------------------------

    p     = pzer*pai**(cpd/rd)
    tvirt = tetav*pai
    t     = tvirt/(1. -(q+qcw+qci+qpw+qpi1+qpi2) +zeps*q)
    zz1 = -g*h*bzita(.5*dz)*log(1.-.5*dz/h)
    do jlat = 1, nlat
    do jlon = 1, nlon
    zdgz = phig(jlon,jlat)*(gzita(.5*dz)-1.) + zz1
    ps(jlon,jlat) = p(jlon,jlat,1)*exp(zdgz/(rd*tvirt(jlon,jlat,1)))
    enddo
    enddo

    tke   = max(tke, tkemin)

    do jklev = 1, nlev
    do jlat = 1, nlat
    do jlon = 1, nlon
    zt0t = tzer/t(jlon,jlat,jklev)
    zesk = ezer*exp(-ccw1*alog(zt0t)+ccw2*(1.-zt0t))    ! partial pressure over water
    qsatw(jlon,jlat,jklev) = zesk*eps/(p(jlon,jlat,jklev)+zesk*(eps-1.)) ! over water
    enddo
    enddo
    enddo

    call tofd

#ifdef mpi
      call u_ghost (v(:     ,2,:), ip_s, v(:,nlat,:), ip_n, nlon*nlev)
      call u_ghost (u(nlonm1,:,:), ip_e, u(1,:   ,:), ip_w, nlat*nlev)
#endif

    call surflayer(jstep)

    if (flag_surf_flux == 0 ) then
! Specific and latent heat fluxes and water vapour fluxes defined by pochva
! scheme
      do jlat=1,nlat
      do jlon=1,nlon
        if (mask_soil(jlon,jlat) == 1) then
          hflux(jlon,jlat)=-(flh_specif(jlon,jlat)+flh_latent(jlon,jlat))
          qflux(jlon,jlat)=-fl_wv(jlon,jlat)
        endif
      enddo
      enddo
    endif

    call vdiff(jstep)

!  orographic blocking and gravity wave drag

    if (nlord) call orogdrag (jstep)

!do jlat=1,nlat
!do jlon=1,nlon
!  if (abs(tskin(jlon,jlat))>400..or.abs(tskin(jlon,jlat))<100..or.isnan(tskin(jlon,jlat))) &
! print '(a,i5,e16.8,3x,3i4.3)',"*4* oshibka * tskin ",jstep,tskin(jlon,jlat),jlon,jlat,myid
!  if (abs(qskin(jlon,jlat))>10.) print *,"*4* oshibka * qskin ",jstep,qskin(jlon,jlat),jlon,jlat,myid
!  if (abs(snow(jlon,jlat))>1000.) print *,"*4* snow ",jstep,snow(jlon,jlat),jlon,jlat,myid
!  if (abs(fsnow(jlon,jlat))>10.) print *,"*4* fsnow ",jstep,fsnow(jlon,jlat),jlon,jlat,myid
!  if (isnan(tskin(jlon,jlat))) print *,"*4* oshibka tskin ",jstep,tskin(jlon,jlat),jlon,jlat,myid
!  if (isnan(qskin(jlon,jlat))) print *,"*4* oshibka qskin ",jstep,qskin(jlon,jlat),jlon,jlat,myid
!  if (isnan(snow(jlon,jlat))) print *,"*4* oshibka snow ",jstep,snow(jlon,jlat),jlon,jlat,myid
!  if (isnan(fsnow(jlon,jlat))) print *,"*4* oshibka fsnow ",jstep,fsnow(jlon,jlat),jlon,jlat,myid
!  if (isnan(v(jlon,jlat,1))) print *,"*4* oshibka v ",jstep,v(jlon,jlat,1),jlon,jlat,myid
!  if (isnan(t(jlon,jlat,1))) print *,"*4* oshibka t ",jstep,t(jlon,jlat,1),jlon,jlat,myid
!  if (isnan(q(jlon,jlat,1))) print *,"*4* oshibka q ",jstep,q(jlon,jlat,1),jlon,jlat,myid
!  if (isnan(tke(jlon,jlat,1))) print *,"*4* oshibka tke ",jstep,tke(jlon,jlat,1),jlon,jlat,myid
!enddo
!enddo

    if (nlhdiff) then
    call hdiffu   ! horizontal turbulent diffusion of u and v ! in progress
    call hdifft   ! horizontal turbulent diffusion of t and q ! in progress
    endif

    tvirt = tetav*pai
    t     = tvirt/(1. -(q+qcw+qci+qpw+qpi1+qpi2) +zeps*q)

!--------------------------------------------------------------------------------------------------------
!  Radiation
!--------------------------------------------------------------------------------------------------------

!  Update albedo and emissivity

 if (mod(jstep,10*ntsrad) == 0) call surfradpar (fmask, fice, qg_surf, vegeta, qsoil_max, qsoil_min, &
 albedo_soil_dry, albedo_soil_wet, emiss1_soil_dry, emiss1_soil_wet, emiss2_soil_dry, emiss2_soil_wet, &
 fveg, albedo_veg, emiss1_veg, emiss2_veg, &
 nlon, nlat, nlevg, nvt, 1, albedo, emisg1, emisg2)

    if (jstep == nstep0.or.mod(jstep-1,ntsrad) == 0) then

    if (nradm > 0) then ! radiation

! definition of aerosol, ozone, solar time reqtim and solar declin. rdecli,
! used in radint for Geleyn radiation - the last two quantites vary slowly so can be kept
! constant for a few days, as in the case radintec is not called)

    if (mod(jstep-1, ntsrad*4*80).eq.0) then
#ifdef rad_ecmwf
      call aerdef (nyrc, nmonc, ndayc, nhouc, nminc)
#else
      call radiat_init(rdecli, reqtim, ozon, aerosol, aerotot, &
   nyrc, nmonc, ndayc, nhouc, nminc, &
   nlon, nlat, nlev, ntype_aerosol, dlat, dlon, &
   alatt, alont, p, ps, t, tskin, phig(:,:)/g, myid)
#endif
    endif

! time lapse of ntsrad*dtstep sec. to compute tendencies of surface radiation fluxes (useful only for solar rad.)

    call calendar (nyrc , nmonc , ndayc , nhouc , nminc , 0, 0, nint(ntsrad*dtstep/60.),  &
                   nyrc1, nmonc1, ndayc1, nhouc1, nminc1, ndayr1)

    call cloudfr                                ! cloud fraction
    call ccloudt                                ! total cloud cover
    call radint (ndayr1, nhouc1, nminc1, jstep) ! Geleyn radiation

#ifdef rad_ecmwf
! ------------ ECMWF radiation scheme --------------

      if (nradm.eq.2 .and. mod(jstep-1,ntsrad*4).eq.0) then  ! correction with ECMWF radiation
      if (myid.eq.0.and.jstep.eq.1) write(*,'(a,i2)') ' ECMWF radiation with McICA =', mcica

      corrdt = 0.
      corvis = 0.
      corirr = 0.

      do 30 jlat = 2, nlat-2, 2
      call radintec (jlat, 2, nyrc, nmonc, ndayc, nhouc, nminc)  ! ECMWF radiation

! definition of correction terms (ECMWF minus Geleyn)

      do jklev = 1, nlev
      do jlon = 2, nlon-2, 2
      corrdt(jlon,jlat,jklev) = corrdt(jlon,jlat,jklev) - geldt(jlon,jlat,jklev)
      enddo
      enddo
      do jlon = 2, nlon-2, 2
      corvis(jlon,jlat)   = corvis(jlon,jlat) - gelvis(jlon,jlat)
      corirr(jlon,jlat)   = corirr(jlon,jlat) - gelirr(jlon,jlat)
      enddo

!  longitudinal interpolation

#ifdef mpi
        call u_ghost (corrdt(2,jlat,:), ip_w, corrdt(nlon,jlat,:), ip_e, nlev)
        call u_ghost (corvis(2,jlat), ip_w, corvis(nlon,jlat), ip_e, 1)
        call u_ghost (corirr(2,jlat), ip_w, corirr(nlon,jlat), ip_e, 1)
        call u_ghost (swsdtf(2,jlat), ip_w, swsdtf(nlon,jlat), ip_e, 1)
        call u_ghost (swsddr(2,jlat), ip_w, swsddr(nlon,jlat), ip_e, 1)
#endif

      do jklev = 1, nlev
      do jlon = 3, nlonm1, 2
      corrdt(jlon,jlat,jklev) = 0.5*(corrdt(jlon-1,jlat,jklev)+corrdt(jlon+1,jlat,jklev))
      enddo
      enddo
      do jlon = 3, nlonm1, 2
      corvis(jlon,jlat) = 0.5*(corvis(jlon-1,jlat)+corvis(jlon+1,jlat))
      corirr(jlon,jlat) = 0.5*(corirr(jlon-1,jlat)+corirr(jlon+1,jlat))
      swsdtf(jlon,jlat) = 0.5*(swsdtf(jlon-1,jlat)+swsdtf(jlon+1,jlat))
      swsddr(jlon,jlat) = 0.5*(swsddr(jlon-1,jlat)+swsddr(jlon+1,jlat))
      enddo
 30   continue

!  latitudinal interpolation

#ifdef mpi
        call u_ghost (corrdt(2:nlonm1,2,:), ip_s, corrdt(2:nlonm1,nlat,:), ip_n, (nlon-2)*nlev)
        call u_ghost (corvis(2:nlonm1,2), ip_s, corvis(2:nlonm1,nlat), ip_n, nlon-2)
        call u_ghost (corirr(2:nlonm1,2), ip_s, corirr(2:nlonm1,nlat), ip_n, nlon-2)
        call u_ghost (swsdtf(2:nlonm1,2), ip_s, swsdtf(2:nlonm1,nlat), ip_n, nlon-2)
        call u_ghost (swsddr(2:nlonm1,2), ip_s, swsddr(2:nlonm1,nlat), ip_n, nlon-2)
#endif

      do jklev = 1, nlev
      do jlat = 3, nlatm1, 2
      do jlon = 2, nlonm1
      corrdt(jlon,jlat,jklev) = 0.5*(corrdt(jlon,jlat-1,jklev)+corrdt(jlon,jlat+1,jklev))
      enddo
      enddo
      enddo
      do jlat = 3, nlatm1, 2
      do jlon = 2, nlonm1
      corvis(jlon,jlat) = 0.5*(corvis(jlon,jlat-1)+corvis(jlon,jlat+1))
      corirr(jlon,jlat) = 0.5*(corirr(jlon,jlat-1)+corirr(jlon,jlat+1))
      swsdtf(jlon,jlat) = 0.5*(swsdtf(jlon,jlat-1)+swsdtf(jlon,jlat+1))
      swsddr(jlon,jlat) = 0.5*(swsddr(jlon,jlat-1)+swsddr(jlon,jlat+1))
      enddo
      enddo

      call filt3d (corrdt, 1.)   !  Smoothing of RADINTEC correction
      call filt2d (corvis, 1.)   !  Smoothing of RADINTEC correction
      call filt2d (corirr, 1.)   !  Smoothing of RADINTEC correction

      endif ! condition on ECMWF radiation
! ------------ End of ECMWF radiation scheme --------------
#endif

      endif ! end of radiation

! Kain-Fritsh parameterization of convection

      raicon = 0.
      snocon = 0.
      dtdt = 0.
      dqdt = 0.
      dqcwdt = 0.
      dqcidt = 0.
      dqpwdt = 0.
      dqpi1dt = 0.

      if (nlconv) then ! convection

        call convection_kf

        call filt3d (dtdt  , .8)   !  Smoothing of convection
        call filt3d (dqdt  , .8)   !  Smoothing of convection
        call filt2d (raicon, .8)   !  Smoothing of convection
        call filt2d (snocon, .8)   !  Smoothing of convection

        raicon = raicon/float(ntsrad)  ! convective rain per timestep
        snocon = snocon/float(ntsrad)

      endif                          ! nlconv end of convection

    dtdt  = dtdt + geldt  + corrdt
    if (jstep.eq.1) frvis = max(0., gelvis + corvis)
!    if (jstep.eq.1) frirr = gelirr + corirr
    frirr = gelirr + corirr
    dfrvis = (max(0., gelvis + corvis)-frvis)/float(ntsrad)
!    dfrirr = (gelirr+corirr           -frirr)/float(ntsrad)

    endif ! condition on ntsrad

! update fields with radiation and convection contributions

    if (nradm.ne.0) then

    do jklev = 1, nlev
    do jlat = 2, nlatm1
    do jlon = 2, nlonm1
    t(jlon,jlat,jklev) = t(jlon,jlat,jklev) + dtdt(jlon,jlat,jklev)*dtstep
    q(jlon,jlat,jklev)   = max (0., q  (jlon,jlat,jklev) + dqdt  (jlon,jlat,jklev)*dtstep)
    qcw(jlon,jlat,jklev) = max (0., qcw(jlon,jlat,jklev) + dqcwdt(jlon,jlat,jklev)*dtstep)
    qci(jlon,jlat,jklev) = max (0., qci(jlon,jlat,jklev) + dqcidt(jlon,jlat,jklev)*dtstep)
    qpw(jlon,jlat,jklev) = max (0., qpw(jlon,jlat,jklev) + dqpwdt(jlon,jlat,jklev)*dtstep)
    qpi1(jlon,jlat,jklev)= max (0.,qpi1(jlon,jlat,jklev) +dqpi1dt(jlon,jlat,jklev)*dtstep)
    enddo
    enddo
    enddo

    frvis = frvis + dfrvis
!    frirr = frirr + dfrirr

    prectot    = prectot + raicon + snocon
    precconv   = precconv + raicon + snocon
    precsolid  = precsolid  + snocon
    qprectot   = qprectot + raicon + snocon
    qprecsolid = qprecsolid + snocon

    endif

!--------------------------------------------------------------------------------------------------------
!  Soil and sea
!--------------------------------------------------------------------------------------------------------

    call sea_surface

! pochva --->

    h_lev_bottom(:,:) = phig(:,:)*(gzita(.5*dz)-1.)/g -h*bzita(.5*dz)*log(1.-.5*dz/h)
    fl_rain(:,:) = raini(:,:)/dtstep
    fl_snow(:,:) = snowi(:,:)/dtstep
! City Heat Island, urban vegetation index is 20 
    fl_rad_tot(:,:) = frvis(:,:)+frirr(:,:)*(1.-0.1*vegeta(:,:,21))
    snow(:,:) = snow(:,:)*1.e3 ! m of equiv. water ---> kg/m/m

!  For option of using heat fluxes as input fileds for Pochva scheme

    if (flag_surf_flux == 1 ) then
! Specific and latent heat fluxes and water vapour fluxes defined by
! surface layer scheme
      flh_specif(:,:) = -hflux(:,:)
      flh_latent(:,:) = -qflux(:,:) * (ylwv + (cpv-cw)*(tskin(:,:)-tzer))
      fl_wv(:,:) = -qflux(:,:)
    endif

if (jstep ==1) then
do k=1,n_outpoint
  if (myid == iproc_outpoint(k)) then
    i=i_outpoint(k)
    j=j_outpoint(k)
    write (file_outpoint,'(a,i2.2,a)') "qsoil_max_",K,".dat"
    open (31, file=file_outpoint, form="formatted", position="append")
    write (31,'(100e16.8)') qsoil_max(i,j,:)
    close (31)
    write (file_outpoint,'(a,i2.2,a)') "qsoil_min_",K,".dat"
    open (31, file=file_outpoint, form="formatted", position="append")
    write (31,'(100e16.8)') qsoil_min(i,j,:)
    close (31)
    write (file_outpoint,'(a,i2.2,a)') "c_soil_",K,".dat"
    open (31, file=file_outpoint, form="formatted", position="append")
    write (31,'(100e16.8)') c_soil(i,j,:)
    close (31)
    write (file_outpoint,'(a,i2.2,a)') "rho_soil_",K,".dat"
    open (31, file=file_outpoint, form="formatted", position="append")
    write (31,'(100e16.8)') rho_soil(i,j,:)
    close (31)
    write (file_outpoint,'(a,i2.2,a)') "psi_soil_",K,".dat"
    open (31, file=file_outpoint, form="formatted", position="append")
    write (31,'(100e16.8)') psi_soil(i,j,:)
    close (31)
    write (file_outpoint,'(a,i2.2,a)') "kw_soil_",K,".dat"
    open (31, file=file_outpoint, form="formatted", position="append")
    write (31,'(100e16.8)') kw_soil(i,j,:)
    close (31)
    write (file_outpoint,'(a,i2.2,a)') "par_b_soil_",K,".dat"
    open (31, file=file_outpoint, form="formatted", position="append")
    write (31,'(100e16.8)') par_b_soil(i,j,:)
    close (31)
    write (file_outpoint,'(a,i2.2,a)') "qsrel_wilt_",K,".dat"
    open (31, file=file_outpoint, form="formatted", position="append")
    write (31,'(100e16.8)') qsrel_wilt(i,j,:)
    close (31)
    write (file_outpoint,'(a,i2.2,a)') "qsrel_ref_",K,".dat"
    open (31, file=file_outpoint, form="formatted", position="append")
    write (31,'(100e16.8)') qsrel_ref(i,j,:)
    close (31)
    write (file_outpoint,'(a,i2.2,a)') "pochva_par_",K,".dat"
    open (31, file=file_outpoint, form="formatted", position="append")
    write (31,'(4f16.3,2i3)') lai(i,j), lai_max, fveg(i,j), proot(i,j),ind_lev_soil_h_bottom(i,j), ind_lev_soil_w_bottom(i,j)
    close (31)
  endif
enddo
endif

    call pochva &
                (lev_soil, tg, tskin, qg, qskin, &
                 snow, fsnow, alsn, &
                 fl_heat_soil_bottom, fl_water_soil_bottom, fl_runoff, fl_runoff_tot, &
                 mask_soil, h_lev_bottom, ps, p(:,:,1), t(:,:,1), q(:,:,1), &
                 fl_rain, fl_snow, fl_rad_tot, frvis, flh_specif, flh_latent, fl_wv, kturb_surf_h, kturb_surf_q, &
                 qsoil_max, qsoil_min, c_soil,  rho_soil, psi_soil, kw_soil,  par_b_soil, par_c_soil, &
                 qsrel_wilt, qsrel_ref, lai, lai_max, fveg, proot, ind_lev_soil_h_bottom, ind_lev_soil_w_bottom, &
                 valmiss, dtstep, jstep, flag_surf_flux, flag_pochva_initial, flag_pochva_par_change, myid, &
                 snow_dirt_ext=snow_dirt, ts_surf_ext=tg_surf, qs_surf_ext=qg_surf, &
                 fice_soil_ext=fice_soil, fice_soil_surf_ext=fice_soil_surf, &
                 flag_snow_init=1, lev_snow_ext=snow_lev, tsnow_ext=snow_t, fice_snow_ext=snow_fice, &
                 snow_age_ext=snow_age, snow_melt_age_ext=snow_melt_age, rho_snow_ext=snow_dens)
    flag_pochva_initial = 0
    flag_pochva_par_change = 0

    snow(:,:) = snow(:,:)*1.e-3 ! kg/m/m ---> m of equiv. water
! <---

!do jlat=1,nlat
!do jlon=1,nlon
!  if (abs(tskin(jlon,jlat))>400..or.abs(tskin(jlon,jlat))<100..or.isnan(tskin(jlon,jlat))) &
! print '(a,i5,e16.8,3x,3i4.3)',"*5* oshibka tskin ",jstep,tskin(jlon,jlat),jlon,jlat,myid
!  if (abs(qskin(jlon,jlat))>10.) print *,"*5* oshibka qskin ",jstep,qskin(jlon,jlat),jlon,jlat,myid
!  if (abs(snow(jlon,jlat))>1000.) print *,"*5* snow ",jstep,snow(jlon,jlat),jlon,jlat,myid
!  if (abs(fsnow(jlon,jlat))>10.) print *,"*5* fsnow ",jstep,fsnow(jlon,jlat),jlon,jlat,myid
!  if (isnan(tskin(jlon,jlat))) print *,"*5* oshibka tskin ",jstep,tskin(jlon,jlat),jlon,jlat,myid
!  if (isnan(qskin(jlon,jlat))) print *,"*5* oshibka qskin ",jstep,qskin(jlon,jlat),jlon,jlat,myid
!  if (isnan(snow(jlon,jlat))) print *,"*5* oshibka snow ",jstep,snow(jlon,jlat),jlon,jlat,myid
!  if (isnan(fsnow(jlon,jlat))) print *,"*5* oshibka fsnow ",jstep,fsnow(jlon,jlat),jlon,jlat,myid
!  if (isnan(v(jlon,jlat,1))) print *,"*5* oshibka v ",jstep,v(jlon,jlat,1),jlon,jlat,myid
!  if (isnan(t(jlon,jlat,1))) print *,"*5* oshibka t ",jstep,t(jlon,jlat,1),jlon,jlat,myid
!  if (isnan(q(jlon,jlat,1))) print *,"*5* oshibka q ",jstep,q(jlon,jlat,1),jlon,jlat,myid
!  if (isnan(tke(jlon,jlat,1))) print *,"*5* oshibka tke ",jstep,tke(jlon,jlat,1),jlon,jlat,myid
!enddo
!enddo

do k=1,n_outpoint
  if (myid == iproc_outpoint(k)) then
    i=i_outpoint(k)
    j=j_outpoint(k)
    close (31)
    write (file_outpoint,'(a,i2.2,a)') "atm_",K,".dat"
    open (31, file=file_outpoint, form="formatted", position="append")
    write (31,'(100e16.8)') float(jstep)*dtstep/3600., &
 p(i,j,1), &               !  2 
 p(i,j,5), &               !  3 
 p(i,j,10), &              !  4 
 p(i,j,30), &              !  5 
 t(i,j,1)-273.15, &        !  6 
 t(i,j,5)-273.15, &        !  7 
 t(i,j,10)-273.15, &       !  8 
 t(i,j,30)-273.15, &       !  9 
 u(i,j,1), &               !  10
 u(i,j,5), &               !  11
 u(i,j,10), &              !  12
 u(i,j,30), &              !  13
 q(i,j,1), &               !  14
 q(i,j,5), &               !  15
 q(i,j,10), &              !  16
 q(i,j,30)                 !  17
    close (31)
    write (file_outpoint,'(a,i2.2,a)') "surf_",K,".dat"
    open (31, file=file_outpoint, form="formatted", position="append")
    write (31,'(100e16.8)') float(jstep)*dtstep/3600., &
 ps(i,j), &                !  2 
 p(i,j,1), &               !  3 
 t(i,j,1)-273.15, &        !  4 
 q(i,j,1), &               !  5 
 h_lev_bottom(i,j), &      !  6
 fl_rain(i,j),&            !  7
 fl_snow(i,j),&            !  8
 fl_rad_tot(i,j),&         !  9
 frvis(i,j),&              ! 10
 flh_specif(i,j),&         ! 11
 flh_latent(i,j),&         ! 12
 fl_wv(i,j),&              ! 13
 kturb_surf_h(i,j),&       ! 14
 kturb_surf_q(i,j),&       ! 15
 albedo(i,j),&             ! 16
 emisg1(i,j),&             ! 17
 emisg2(i,j)               ! 18
    close (31)
    write (file_outpoint,'(a,i2.2,a)') "tsoil_",K,".dat"
    open (31, file=file_outpoint, form="formatted", position="append")
    write (31,'(100e16.8)') float(jstep)*dtstep/3600., &
 tskin(i,j)-273.15, &      !  2 
 tg_surf(i,j)-273.15, &    !  3 
 tg(i,j,1:nlevg)-273.15    !  4:12 
    close (31)
    write (file_outpoint,'(a,i2.2,a)') "qsoil_",K,".dat"
    open (31, file=file_outpoint, form="formatted", position="append")
    write (31,'(100e16.8)') float(jstep)*dtstep/3600., &
 qskin(i,j), &             !  2 
 qg_surf(i,j), &           !  3
 qg(i,j,1:nlevg)           !  4:12
    close (31)
    write (file_outpoint,'(a,i2.2,a)') "qsoil_rel_",K,".dat"
    open (31, file=file_outpoint, form="formatted", position="append")
    write (31,'(100e16.8)') float(jstep)*dtstep/3600., &
 qskin(i,j), &            !  2 
 (qg_surf(i,j)-qsoil_min(i,j,1))/(qsoil_max(i,j,1)-qsoil_min(i,j,1)), & !  3 
 (qg(i,j,1:nlevg)-qsoil_min(i,j,1:nlevg))/(qsoil_max(i,j,1:nlevg)-qsoil_min(i,j,1:nlevg)) !  4:12 
    close (31)
    write (file_outpoint,'(a,i2.2,a)') "snow_lev_",K,".dat"
    open (31, file=file_outpoint, form="formatted", position="append")
    write (31,'(100e16.8)') float(jstep)*dtstep/3600., &
 snow(i,j), &             !  2 
 fsnow(i,j), &            !  3 
 alsn(i,j),&              !  4 
 snow_lev(i,j,1:nlev_snow) !  5:17 
    close (31)
    write (file_outpoint,'(a,i2.2,a)') "fice_soil_",K,".dat"
    open (31, file=file_outpoint, form="formatted", position="append")
    write (31,'(100e16.8)') float(jstep)*dtstep/3600., &
 fice_soil_surf(i,j), &            !  2 
 fice_soil(i,j,1:nlevg)            !  3:11 
    close (31)
    write (file_outpoint,'(a,i2.2,a)') "snow_t_",K,".dat"
    open (31, file=file_outpoint, form="formatted", position="append")
    write (31,'(100e16.8)') float(jstep)*dtstep/3600., &
 snow_t(i,j,1:nlev_snow)           !  2:14
    close (31)
    write (file_outpoint,'(a,i2.2,a)') "snow_dens_",K,".dat"
    open (31, file=file_outpoint, form="formatted", position="append")
    write (31,'(100e16.8)') float(jstep)*dtstep/3600., &
 snow_dens(i,j,1:nlev_snow)        !  2:14
    close (31)
  endif
enddo
!#ifdef mpi
!        call mpi_barrier(comm, ierr)
!#endif
!do jlat=1,nlat
!do jlon=1,nlon
!  if (abs(tskin(jlon,jlat))>400..or.abs(tskin(jlon,jlat))<100..or.isnan(tskin(jlon,jlat))) then
!#ifdef mpi
!      call mpi_abort(comm, ierr)
!#endif
!   stop
!  endif
!enddo
!enddo

!--------------------------------------------------------------------------------------------------------
!  Microphysics -- update prectot and snow fall
!--------------------------------------------------------------------------------------------------------

    call micro2m (jstep)

    prectot    = prectot + raini + snowi
    precsolid  = precsolid + snowi
    qprectot   = qprectot + raini + snowi
    qprecsolid = qprecsolid +snowi

!--------------------------------------------------------------------------------------------------------
! computation of accumulated, averaged and min/max quantities
!--------------------------------------------------------------------------------------------------------

    do jlat = 2, nlatm1
    do jlon = 2, nlonm1
    totsky(jlon,jlat) = totsky(jlon,jlat) + swsdtf(jlon,jlat)
    soldir(jlon,jlat) = soldir(jlon,jlat) + swsddr(jlon,jlat)
    cswfl (jlon,jlat) = cswfl (jlon,jlat) + frvis(jlon,jlat)*dtstep
    clwfl (jlon,jlat) = clwfl (jlon,jlat) + frirr(jlon,jlat)*dtstep
    cshflux(jlon,jlat) = cshflux(jlon,jlat) - hflux(jlon,jlat)*dtstep
    shf_accum(jlon,jlat) = shf_accum(jlon,jlat) - hflux(jlon,jlat)
    zhea = ylwv + (cpv-cw)*(tskin(jlon,jlat)-tzer)
    clhflux(jlon,jlat) = clhflux(jlon,jlat) - qflux(jlon,jlat)*dtstep*zhea
    lhf_accum(jlon,jlat) = lhf_accum(jlon,jlat) - qflux(jlon,jlat)*zhea
    qf_accum(jlon,jlat) = qf_accum(jlon,jlat) - qflux(jlon,jlat)*dtstep
    cwvflux(jlon,jlat) = cwvflux(jlon,jlat) - qflux(jlon,jlat)*dtstep
    t2min(jlon,jlat) = min (t2min(jlon,jlat), t_std_lev(jlon,jlat,1))
    t2max(jlon,jlat) = max (t2max(jlon,jlat), t_std_lev(jlon,jlat,1))
    ws10max(jlon,jlat) = max (ws10max(jlon,jlat), sqrt(u_std_lev(jlon,jlat,2)**2+v_std_lev(jlon,jlat,2)**2))
    cfl_heat_soil_bottom(jlon,jlat) = cfl_heat_soil_bottom(jlon,jlat) + fl_heat_soil_bottom(jlon,jlat)*dtstep
    cfl_water_soil_bottom(jlon,jlat) = cfl_water_soil_bottom(jlon,jlat) + fl_water_soil_bottom(jlon,jlat)*dtstep
    runoff(jlon,jlat) = runoff(jlon,jlat) + fl_runoff(jlon,jlat)*dtstep
    runoff_tot(jlon,jlat) = runoff_tot(jlon,jlat) + fl_runoff_tot(jlon,jlat)*dtstep
    enddo
    enddo

!-----------------------------------------------------------------------
! lateral boundary conditions
!-----------------------------------------------------------------------

      if (nlbfix) then

        ztime = 0.

      else

!  defines the boundary conditions at time jstep*dtstep by linear interpolation in time.

        ztime = float(mod(jstep,ntsbou))/float(ntsbou)
        if (jstep.ge.(nbc-1)*ntsbou) ztime = 1.

        if (mod(jstep,ntsbou).eq.0.and.jstep.lt.(nbc-1)*ntsbou) then    !  read new boundary file

          ub1   = ub2
          vb1   = vb2
          wb1   = wb2
          pb1   = pb2
          tb1   = tb2
          qb1   = qb2
          qcwb1 = qcwb2
          qcib1 = qcib2

          interv = jstep/ntsbou + 1

          call rdmhf_atm (interv+1, 0)
          call paidef (pb2, tb2, qb2, qcwb2, qcib2, pf)
          pb2 = pf

          call rdmhf_soil (interv, 0)

          if (nlclimate) then

!!            tg(:,:,nlevg) = tgclim(:,:)  ! uncomment here for long runs (with variable soil T)
!!            qg(:,:,nlevg) = qgclim(:,:)  ! uncomment here for long runs (with variable soil T)
            do jlat = 2, nlatm1
            do jlon = 2, nlonm1
              snow(jlon,jlat) = snow(jlon,jlat)*1.e3 ! m of equiv. water ---> kg/m/m
              if (fmask(jlon,jlat) >= 0.5) then ! sea
                if (fice(jlon,jlat) < 0.8.and.fice_rd(jlon,jlat) >= 0.8) then ! sea point becames "soil" because of sea ice increasing until the threshold value
                  mask_soil(jlon,jlat) = 1
                  psi_soil(jlon,jlat,:) = 2.2 ! thermal conductivity of ice
                  albedo_soil_dry(jlon,jlat)=0.70
                  albedo_soil_wet(jlon,jlat)=0.70
                  emiss1_soil_dry(jlon,jlat)=0.98
                  emiss1_soil_wet(jlon,jlat)=0.98
                  emiss2_soil_dry(jlon,jlat)=0.98
                  emiss2_soil_wet(jlon,jlat)=0.98
                  snow_lev(jlon,jlat,:) = valmiss
                  snow_t(jlon,jlat,:) = valmiss
                  snow_fice(jlon,jlat,:) = valmiss
                  snow_age(jlon,jlat,:) = valmiss
                  snow_melt_age(jlon,jlat,:) = valmiss
                  snow_dens(jlon,jlat,:) = valmiss
                  fl_heat_soil_bottom(jlon,jlat) = 0.
                  fl_water_soil_bottom(jlon,jlat) = 0.
                  fl_runoff(jlon,jlat) = 0.
                  fl_runoff_tot(jlon,jlat) = 0.
                  tskin(jlon,jlat) = min(tskin(jlon,jlat), tzer-0.1)
                  zdlev=dlev_base
                  do while (.true.)
                    if (zdlev*float(nlev_snow-1) >= snow(jlon,jlat)) exit
                    zdlev=zdlev+dlev_base
                  enddo
                  snow_lev(jlon,jlat,1)=0.
                  if (snow(jlon,jlat) >= zdlev+dlev_min) then
                    n=int(snow(jlon,jlat)/zdlev)+1
                    if (snow(jlon,jlat)-zdlev*float(n-1) >= dlev_min) n=n+1
                    do jklev=2,n
                      snow_lev(jlon,jlat,jklev)=snow(jlon,jlat)-float(n-jklev)*zdlev
                    enddo
                    n=n-1
                  else
                    n=1
                    snow_lev(jlon,jlat,2)=snow(jlon,jlat)
                  endif
                  snow_t(jlon,jlat,1) = tskin(jlon,jlat)
                  if (n < 2) then
                    snow_t(jlon,jlat,2) = min(tg(jlon,jlat,1), tzer-0.1)
                  else
                    snow_t(jlon,jlat,n+1) = min(tg(jlon,jlat,1), tzer-0.1)
                    z1=1./float(n)
                    do jklev=2,n
                      snow_t(jlon,jlat,jklev) = snow_t(jlon,jlat,1)*(1.-z1*float(jklev-1)) + &
 snow_t(jlon,jlat,n+1)*z1*float(jklev-1)
                    enddo
                  endif
                  snow_fice(jlon,jlat,1:n+1) = 1.
                  snow_age(jlon,jlat,1:n+1) = 180.
                  snow_melt_age(jlon,jlat,1:n+1) = 0.
                  do jklev = 1, nlevg
                    if (lev_soil(jklev) > iceth(jlon,jlat)) exit
                  enddo
                  ind_lev_soil_h_bottom(jlon,jlat) = max(jklev-1, 3)
                  ind_lev_soil_w_bottom(jlon,jlat) = 0
                  do jklev = 1, ind_lev_soil_h_bottom(jlon,jlat)
                    tg(jlon,jlat,jklev) = min(tg(jlon,jlat,jklev), tzer-0.1)
                  enddo
                endif
                if (fice(jlon,jlat) >= 0.8.and.fice_rd(jlon,jlat) < 0.8) then ! sea point returns to be "sea" because of sea ice decreasing under the threshold value
                  mask_soil(jlon,jlat) = 0
                  psi_soil(jlon,jlat,:) = 0.57 ! thermal conductivity of water
                  albedo_soil_dry(jlon,jlat)=0.07
                  albedo_soil_wet(jlon,jlat)=0.07
                  emiss1_soil_dry(jlon,jlat)=0.97
                  emiss1_soil_wet(jlon,jlat)=0.97
                  emiss2_soil_dry(jlon,jlat)=0.97
                  emiss2_soil_wet(jlon,jlat)=0.97
                  tskin(jlon,jlat) = tgclim(jlon,jlat)
                  tg(jlon,jlat,:) = tgclim(jlon,jlat)
                  snow(jlon,jlat) = 0.
                  snow_lev(jlon,jlat,:) = valmiss
                  snow_t(jlon,jlat,:) = valmiss
                  snow_fice(jlon,jlat,:) = valmiss
                  snow_age(jlon,jlat,:) = valmiss
                  snow_melt_age(jlon,jlat,:) = valmiss
                  snow_dens(jlon,jlat,:) = valmiss
                  fl_heat_soil_bottom(jlon,jlat) = 0.
                  fl_water_soil_bottom(jlon,jlat) = 0.
                  fl_runoff(jlon,jlat) = 0.
                  fl_runoff_tot(jlon,jlat) = 0.
                endif
              endif
              if (mask_soil(jlon,jlat) == 0) then  ! over sea (excluding thick sea ice)
                tg(jlon,jlat,2:nlevg) = tgclim(jlon,jlat)
              endif
              snow(jlon,jlat) = snow(jlon,jlat)*1.e-3 ! kg/m/m ---> m of equiv. water
            enddo
            enddo
            fice(:,:) = fice_rd(:,:)
            iceth(:,:) = iceth_rd(:,:)
            call collect (fice, gfield)
            if (myid == 0) print *,'Updated fice, ',maxval(gfield),maxloc(gfield)
          endif ! nlclimate

        endif ! read new boundary file

      endif ! .not.nlbfix

! ------------------------------------------------------
! Large Scale Nudging of T, U, V
! (only for the cases of large grid spacing at a large domain)

  if (nlconv.and.nlsnudg) then
    if (jstep == 1.or.mod(jstep,ntnudg) == 0) then
      call nudging_tuv(jstep, ztime)
    endif
  endif
! ------------------------------------------------------

!  linear interpolation in time between b1 and b2

      if (ip_w.eq.ip_null) then
      do jklev = 1, nlev
      do jlat = 1, nlat
         jlo2 = nbl
         if (myid == 0) jlo2 = min (jlat,nbl)
         if (myid.eq.nprocsy-1) jlo2 = min (nlat-jlat+1,nbl)
      do jlon = 1, jlo2
      zub   = ub1  (jlon,jlat,jklev) + ztime*(ub2  (jlon,jlat,jklev)-ub1  (jlon,jlat,jklev))
      zvb   = vb1  (jlon,jlat,jklev) + ztime*(vb2  (jlon,jlat,jklev)-vb1  (jlon,jlat,jklev))
      zwb   = wb1  (jlon,jlat,jklev) + ztime*(wb2  (jlon,jlat,jklev)-wb1  (jlon,jlat,jklev))
      zpb   = pb1  (jlon,jlat,jklev) + ztime*(pb2  (jlon,jlat,jklev)-pb1  (jlon,jlat,jklev))
      ztb   = tb1  (jlon,jlat,jklev) + ztime*(tb2  (jlon,jlat,jklev)-tb1  (jlon,jlat,jklev))
      zqb   = qb1  (jlon,jlat,jklev) + ztime*(qb2  (jlon,jlat,jklev)-qb1  (jlon,jlat,jklev))
      zqcwb = qcwb1(jlon,jlat,jklev) + ztime*(qcwb2(jlon,jlat,jklev)-qcwb1(jlon,jlat,jklev))
      zqcib = qcib1(jlon,jlat,jklev) + ztime*(qcib2(jlon,jlat,jklev)-qcib1(jlon,jlat,jklev))

! limitation of cloud water and ice to prevent excessive precip. at the boundaries

      if (.not.nlmic2) zqcwb = min (zqcwb, bcldwat)
      if (.not.nlmic2) zqcib = min (zqcib, bcldice)

      z1 = bndrel(jlon)
      z2 = 1.-z1
      u  (jlon,jlat,jklev) = z1*zub   + z2*u  (jlon,jlat,jklev)
      v  (jlon,jlat,jklev) = z1*zvb   + z2*v  (jlon,jlat,jklev)
      w  (jlon,jlat,jklev) = z1*zwb   + z2*w  (jlon,jlat,jklev)
      pai(jlon,jlat,jklev) = z1*zpb   + z2*pai(jlon,jlat,jklev)
      t  (jlon,jlat,jklev) = z1*ztb   + z2*t  (jlon,jlat,jklev)
      q  (jlon,jlat,jklev) = z1*zqb   + z2*q  (jlon,jlat,jklev)
      qcw(jlon,jlat,jklev) = z1*zqcwb + z2*qcw(jlon,jlat,jklev)
      qci(jlon,jlat,jklev) = z1*zqcib + z2*qci(jlon,jlat,jklev)
      if (nlmic2) ncw(jlon,jlat,jklev) = z1*1.e10 + z2*ncw(jlon,jlat,jklev)
      if (nlmic2) nci(jlon,jlat,jklev) = z1*1.e10 + z2*nci(jlon,jlat,jklev)
      enddo
      enddo
      enddo
      do jlat = 1, nlat
      qskin(1,jlat) = qskin(2,jlat)
      tskin(1,jlat) = tskin(2,jlat)
      enddo
      endif

      if (ip_e.eq.ip_null) then
      do jklev = 1, nlev
      do jlat = 1, nlat
         jlo2 = 1
         if (myid.eq.iprocs-nprocsy) jlo2 = max(nbl-jlat+1   ,1)
         if (myid.eq.iprocs-1      ) jlo2 = max(nbl+jlat-nlat,1)
      do jlon = jlo2, nbl
      jlo1 = nlon-nbl+jlon
      zub   = ub1  (jlo1,jlat,jklev) + ztime*(ub2  (jlo1,jlat,jklev)-ub1  (jlo1,jlat,jklev))
      zvb   = vb1  (jlo1,jlat,jklev) + ztime*(vb2  (jlo1,jlat,jklev)-vb1  (jlo1,jlat,jklev))
      zwb   = wb1  (jlo1,jlat,jklev) + ztime*(wb2  (jlo1,jlat,jklev)-wb1  (jlo1,jlat,jklev))
      zpb   = pb1  (jlo1,jlat,jklev) + ztime*(pb2  (jlo1,jlat,jklev)-pb1  (jlo1,jlat,jklev))
      ztb   = tb1  (jlo1,jlat,jklev) + ztime*(tb2  (jlo1,jlat,jklev)-tb1  (jlo1,jlat,jklev))
      zqb   = qb1  (jlo1,jlat,jklev) + ztime*(qb2  (jlo1,jlat,jklev)-qb1  (jlo1,jlat,jklev))
      zqcwb = qcwb1(jlo1,jlat,jklev) + ztime*(qcwb2(jlo1,jlat,jklev)-qcwb1(jlo1,jlat,jklev))
      zqcib = qcib1(jlo1,jlat,jklev) + ztime*(qcib2(jlo1,jlat,jklev)-qcib1(jlo1,jlat,jklev))
      if (.not.nlmic2) zqcwb = min (zqcwb, bcldwat)
      if (.not.nlmic2) zqcib = min (zqcib, bcldice)

      z1 = bndrel(nbl-jlon+1)
      z2 = 1.-z1
      u  (jlo1,jlat,jklev) = z1*zub   + z2*u  (jlo1,jlat,jklev)
      v  (jlo1,jlat,jklev) = z1*zvb   + z2*v  (jlo1,jlat,jklev)
      w  (jlo1,jlat,jklev) = z1*zwb   + z2*w  (jlo1,jlat,jklev)
      pai(jlo1,jlat,jklev) = z1*zpb   + z2*pai(jlo1,jlat,jklev)
      t  (jlo1,jlat,jklev) = z1*ztb   + z2*t  (jlo1,jlat,jklev)
      q  (jlo1,jlat,jklev) = z1*zqb   + z2*q  (jlo1,jlat,jklev)
      qcw(jlo1,jlat,jklev) = z1*zqcwb + z2*qcw(jlo1,jlat,jklev)
      qci(jlo1,jlat,jklev) = z1*zqcib + z2*qci(jlo1,jlat,jklev)
      if (nlmic2) ncw(jlo1,jlat,jklev) = z1*1.e10 + z2*ncw(jlo1,jlat,jklev)
      if (nlmic2) nci(jlo1,jlat,jklev) = z1*1.e10 + z2*nci(jlo1,jlat,jklev)
      enddo
      enddo
      enddo
      do jlat = 1, nlat
      qskin(nlon,jlat) = qskin(nlonm1,jlat)
      tskin(nlon,jlat) = tskin(nlonm1,jlat)
      enddo
      endif

      if (ip_s.eq.ip_null) then
      do jklev = 1, nlev
      do jlat = 1, nbl
         jlo1 = 1
         if (myid == 0) jlo1 = jlat+1
         jlo2 = nlon
         if (myid.eq.iprocs-nprocsy) jlo2 = nlon-jlat
      z1 = bndrel(jlat)
      z2 = 1.-z1
      do jlon = jlo1, jlo2
      zub   = ub1  (jlon,jlat,jklev) + ztime*(ub2  (jlon,jlat,jklev)-ub1  (jlon,jlat,jklev))
      zvb   = vb1  (jlon,jlat,jklev) + ztime*(vb2  (jlon,jlat,jklev)-vb1  (jlon,jlat,jklev))
      zwb   = wb1  (jlon,jlat,jklev) + ztime*(wb2  (jlon,jlat,jklev)-wb1  (jlon,jlat,jklev))
      zpb   = pb1  (jlon,jlat,jklev) + ztime*(pb2  (jlon,jlat,jklev)-pb1  (jlon,jlat,jklev))
      ztb   = tb1  (jlon,jlat,jklev) + ztime*(tb2  (jlon,jlat,jklev)-tb1  (jlon,jlat,jklev))
      zqb   = qb1  (jlon,jlat,jklev) + ztime*(qb2  (jlon,jlat,jklev)-qb1  (jlon,jlat,jklev))
      zqcwb = qcwb1(jlon,jlat,jklev) + ztime*(qcwb2(jlon,jlat,jklev)-qcwb1(jlon,jlat,jklev))
      zqcib = qcib1(jlon,jlat,jklev) + ztime*(qcib2(jlon,jlat,jklev)-qcib1(jlon,jlat,jklev))
      if (.not.nlmic2) zqcwb = min (zqcwb, bcldwat)
      if (.not.nlmic2) zqcib = min (zqcib, bcldice)

      u  (jlon,jlat,jklev) = z1*zub   + z2*u  (jlon,jlat,jklev)
      v  (jlon,jlat,jklev) = z1*zvb   + z2*v  (jlon,jlat,jklev)
      w  (jlon,jlat,jklev) = z1*zwb   + z2*w  (jlon,jlat,jklev)
      pai(jlon,jlat,jklev) = z1*zpb   + z2*pai(jlon,jlat,jklev)
      t  (jlon,jlat,jklev) = z1*ztb   + z2*t  (jlon,jlat,jklev)
      q  (jlon,jlat,jklev) = z1*zqb   + z2*q  (jlon,jlat,jklev)
      qcw(jlon,jlat,jklev) = z1*zqcwb + z2*qcw(jlon,jlat,jklev)
      qci(jlon,jlat,jklev) = z1*zqcib + z2*qci(jlon,jlat,jklev)
      if (nlmic2) ncw(jlon,jlat,jklev) = z1*1.e10 + z2*ncw(jlon,jlat,jklev)
      if (nlmic2) nci(jlon,jlat,jklev) = z1*1.e10 + z2*nci(jlon,jlat,jklev)
      enddo
      enddo
      enddo
      do jlon = 1, nlon
      qskin(jlon,1) = qskin(jlon,2)
      tskin(jlon,1) = tskin(jlon,2)
      enddo
      endif

      if (ip_n.eq.ip_null) then
      do jklev = 1, nlev
      do jlat = 1, nbl
      jla1 = nlat-nbl+jlat
         jlo1 = 1
         if (myid.eq.nprocsy-1) jlo1 = nbl-jlat+2
         jlo2 = nlon
         if (myid.eq.iprocs-1) jlo2 = nlonm1-nbl+jlat
      z1 = bndrel(nbl-jlat+1)
      z2 = 1.-z1
      do jlon = jlo1, jlo2
      zub   = ub1  (jlon,jla1,jklev) + ztime*(ub2  (jlon,jla1,jklev)-ub1  (jlon,jla1,jklev))
      zvb   = vb1  (jlon,jla1,jklev) + ztime*(vb2  (jlon,jla1,jklev)-vb1  (jlon,jla1,jklev))
      zwb   = wb1  (jlon,jla1,jklev) + ztime*(wb2  (jlon,jla1,jklev)-wb1  (jlon,jla1,jklev))
      zpb   = pb1  (jlon,jla1,jklev) + ztime*(pb2  (jlon,jla1,jklev)-pb1  (jlon,jla1,jklev))
      ztb   = tb1  (jlon,jla1,jklev) + ztime*(tb2  (jlon,jla1,jklev)-tb1  (jlon,jla1,jklev))
      zqb   = qb1  (jlon,jla1,jklev) + ztime*(qb2  (jlon,jla1,jklev)-qb1  (jlon,jla1,jklev))
      zqcwb = qcwb1(jlon,jla1,jklev) + ztime*(qcwb2(jlon,jla1,jklev)-qcwb1(jlon,jla1,jklev))
      zqcib = qcib1(jlon,jla1,jklev) + ztime*(qcib2(jlon,jla1,jklev)-qcib1(jlon,jla1,jklev))
      if (.not.nlmic2) zqcwb = min (zqcwb, bcldwat)
      if (.not.nlmic2) zqcib = min (zqcib, bcldice)

      u  (jlon,jla1,jklev) = z1*zub   + z2*u  (jlon,jla1,jklev)
      v  (jlon,jla1,jklev) = z1*zvb   + z2*v  (jlon,jla1,jklev)
      w  (jlon,jla1,jklev) = z1*zwb   + z2*w  (jlon,jla1,jklev)
      pai(jlon,jla1,jklev) = z1*zpb   + z2*pai(jlon,jla1,jklev)
      t  (jlon,jla1,jklev) = z1*ztb   + z2*t  (jlon,jla1,jklev)
      q  (jlon,jla1,jklev) = z1*zqb   + z2*q  (jlon,jla1,jklev)
      qcw(jlon,jla1,jklev) = z1*zqcwb + z2*qcw(jlon,jla1,jklev)
      qci(jlon,jla1,jklev) = z1*zqcib + z2*qci(jlon,jla1,jklev)
      if (nlmic2) ncw(jlon,jla1,jklev) = z1*1.e10 + z2*ncw(jlon,jla1,jklev)
      if (nlmic2) nci(jlon,jla1,jklev) = z1*1.e10 + z2*nci(jlon,jla1,jklev)
      enddo
      enddo
      enddo
      do jlon = 1, nlon
      qskin(jlon,nlat) = qskin(jlon,nlatm1)
      tskin(jlon,nlat) = tskin(jlon,nlatm1)
      enddo
      endif

!--------------------------------------------------------------------------------------------------------
! Check point
!--------------------------------------------------------------------------------------------------------

    call pbl_height

    iday =  int(jstep*dtstep/86400.)
    ihou =  int((jstep*dtstep-iday*86400.)/3600.)
    imin =  int((jstep*dtstep-iday*86400.-ihou*3600.)/60. +.5)
      if(imin.ge.60) then
      imin = imin - 60
      ihou = ihou + 1
      endif
      if(ihou.ge.24) then
      ihou = ihou -24
      iday = iday +1
      endif
    nfdr(10) = iday
    nfdr(11) = ihou
    nfdr(12) = imin
    nfdrr(10:12) = nfdr(10:12)

    if (mod(jstep,nhist).eq.0) then

      call cloudfr
      call ccloudt
      call snowfall_level
      if (mhfr.eq.2) then ! precipitation averaged conserving area integral
        call aver (prectot)
        call aver (precsolid)
      endif
      if (myid == 0) then
        write (*,'(a)') ' Check-point: writing on moloch mhf_atm and mhf_soil files'
        imhf = imhf + 1
      endif
      call wrmhf_atm
      call wrmhf_soil
      if (nlradar) call wrshf_radar(jstep)

! Reset cumulated quantities written in the MHF

      prectot    = 0.
      precconv   = 0.
      precsolid  = 0.
      runoff     = 0.
      runoff_tot = 0.
      cswfl      = 0.
      clwfl      = 0.
      cshflux    = 0.
      clhflux    = 0.
      t2min      = 999.
      t2max      = 0.
      ws10max    = 0.
      cwvflux    = 0.
      cfl_heat_soil_bottom = 0.
      cfl_water_soil_bottom = 0.

    endif ! mod(jstep,nhist).eq.0

    if (mhfr == 2.and.mod(jstep,nhist_full_res) == 0) then
      if (myid == 0) then
!!!        write (*,'(a)') ' Check-point: writing on moloch mhf_atm and mhf_soil files with full resoltion '
        write (*,'(a)') ' Check-point: writing on moloch mhf_soil files with full resoltion '
      endif
!!!      call wrmhf_atm_full_res
      call wrmhf_soil_full_res
    endif

    if (mod(jstep,ndrunt).eq.0) call runout (jstep)  ! run time diagnostic

    if (ntswshf > 0.and.(mod(jstep,ntswshf).eq.0.or.jstep.eq.1)) then
    call cloudfr
    call ccloudt
    call wrshf (jstep) ! write SHF
    endif

    if (ntback.gt.0.and.(mod(jstep,ntback).eq.0)) then
    call wrback                      ! write trajectory file
    endif
    if (ntspray.gt.0.and.(mod(jstep,ntspray).eq.0)) then
      inst_spray = inst_spray+1
      call wrspray (i1spray, j1spray,inst_spray)  ! write SPRAY file
    endif

! Writting of restart file
    if (mod(jstep,nrst) == 0.and.jstep > nstep0) then
      call wrrf(jstep)
      call wrrf_pochva(jstep, myid, gfield, gnlon, gnlat, irf)
    endif

!********************************************************************************************************
2000 continue  ! End of timestep
!********************************************************************************************************

!    call collect (fmask, gfield1) ! examples for 2-D plotting
!    call collect (phig, gfield)
!    str = 'phig'
!    if (myid == 0) call plotout (gfield, gfield1, gnlon, gnlat, str, 99)
!    str = 'tg1'
!    call collect (tg(1:nlon,1:nlat,1), gfield)
!    if (myid == 0) call plotout (gfield, gfield1, gnlon, gnlat, str, 99)
!#ifdef mpi
!      call mpi_finalize (ierr)
!#endif
!    stop

    call system_clock (stime2, countrate, countmax)
    if (myid == 0) then
      if (stime2.ge.stime1) then
      print*,'System time (sec)', float(stime2-stime1)/countrate
      else
      print*,'System time (sec)', float(countmax-stime1+stime2)/countrate
      endif
    endif

#ifdef mpi
      call mpi_finalize (ierr)          ! clean up MPI environment
#endif

    end program moloch
!###############################################################################################################
 subroutine rdmhf_atm (nun, init_flag)

!  Reads the Model History File (MHF) with prognostic atmospheric variables

    use mod_moloch
    implicit none

    integer :: nun, init_flag, iunit=12, jklev, isleep, istop, no_input, ierr_open
    character(len=30) :: filerd, anun
    integer :: comm, error

#ifdef mpi
      include 'mpif.h'

      comm = mpi_comm_world
#endif

    istop = 2
    no_input = 0
    isleep = 0

    if (myid == 0) then

      print *

      write (anun,'(i2.2)') nun

      filerd="input_atm_"//trim(anun)//".mhf"

      do while (.true.)
        open (iunit, file=filerd, form='unformatted', status='old', iostat=ierr_open)
        if (ierr_open == 0) then
#ifdef oper
          call system("sync")
          call system("ls -l -L "//filerd)
          call system("date")
          call system("sleep 1")
#endif
          exit
        else
          print *,"Input file ",filerd," is not ready: wait - sleep 60 s..."
          isleep = isleep + 1
#ifndef oper
          if (isleep == istop) then
            print*, 'Input file no.', nun, ' not found. Program stopped!'
            no_input = 1
            exit
          endif
#endif
          call system ("sleep 60")
        endif
      enddo

    endif

#ifdef mpi
      call mpi_bcast (no_input, 1, mpi_integer, 0, comm, error)
#endif
    if (no_input == 1) then
#ifdef mpi
        call mpi_finalize (error)
#endif
      stop
    endif

!  read descriptor records and broadcast to all processes

    if (myid == 0) then
      if (init_flag == 1) then
        read(iunit) nfdr0
        read(iunit) pdr0
      else
        read(iunit)
        read(iunit)
      endif
    endif

#ifdef mpi
      if (init_flag == 1) call mpi_bcast (nfdr0 , 50 , mpi_integer, 0, comm, error)
      if (init_flag == 1) call mpi_bcast (pdr0  , 100, mpi_real   , 0, comm, error)
#endif

    if (nfdr0(2) /= gnlon.or.nfdr0(3) /= gnlat.or.nfdr0(4) /= nlev.or.nfdr0(15) /= nlevg) then  ! check on MHF length
      if (myid == 0) then
        print*, ' ***** MHF atm error: check grid dimensions *****'
        print *,' Defined in dimension.inc gnlon, gnlat, nlev, nlevg: ',gnlon, gnlat, nlev, nlevg
        print *,' Read from input mhf:                                ',nfdr0(2), nfdr0(3), nfdr0(4), nfdr0(15)
      endif
#ifdef mpi
        call mpi_finalize (error)        ! clean up MPI environment
#endif
      stop
    endif

    do jklev = 1, nlev
      if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
      if (init_flag == 1) call disperse (gfield, p  (1,1,jklev))
      if (init_flag /= 1) call disperse (gfield, pb2(1,1,jklev))
    enddo
    do jklev = 1, nlev
      if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
      if (init_flag == 1) call disperse (gfield, u  (1,1,jklev))
      if (init_flag /= 1) call disperse (gfield, ub2(1,1,jklev))
    enddo
    do jklev = 1, nlev
      if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
      if (init_flag == 1) call disperse (gfield, v  (1,1,jklev))
      if (init_flag /= 1) call disperse (gfield, vb2(1,1,jklev))
    enddo
    do jklev = 1, nlevp1
      if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
      if (init_flag == 1) call disperse (gfield, w  (1,1,jklev))
      if (init_flag /= 1) call disperse (gfield, wb2(1,1,jklev))
    enddo
    do jklev = 1, nlev
      if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
      if (init_flag == 1) call disperse (gfield, t  (1,1,jklev))
      if (init_flag /= 1) call disperse (gfield, tb2(1,1,jklev))
    enddo
    do jklev = 1, nlev
      if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
      if (init_flag == 1) call disperse (gfield, q  (1,1,jklev))
      if (init_flag /= 1) call disperse (gfield, qb2(1,1,jklev))
    enddo
    do jklev = 1, nlev
      if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
      if (init_flag == 1) call disperse (gfield, qcw  (1,1,jklev))
      if (init_flag /= 1) call disperse (gfield, qcwb2(1,1,jklev))
    enddo
    do jklev = 1, nlev
      if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
      if (init_flag == 1) call disperse (gfield, qci  (1,1,jklev))
      if (init_flag /= 1) call disperse (gfield, qcib2(1,1,jklev))
    enddo

    if (myid == 0) then
      close (iunit)
!     call system("date")
      write (*,'(2a)') ' Read file ',trim(filerd)
      write (*,*)
    endif

 return
 end subroutine rdmhf_atm
!###############################################################################################################
 subroutine rdmhf_soil (nun, init_flag)

!  Reads the Model History File (MHF) with prognostic surface, sea and soil variables

    use mod_moloch
    implicit none

    integer :: nun, init_flag, iunit=13, k, isleep, istop, no_input, ierr_open, ierr_read
    character(len=30) :: filerd, anun
    integer, dimension(50) :: nfdr_local
    real, dimension(100) :: pdr_local
    real, dimension(nlon,nlat) :: field2d_add
    integer :: comm, error

#ifdef mpi
      include 'mpif.h'
      comm = mpi_comm_world
#endif

    istop = 2
    no_input = 0
    isleep = 0

    if (myid == 0) then

      print *

      write (anun,'(i2.2)') nun

      filerd="input_soil_"//trim(anun)//".mhf"

      do while (.true.)
        open (iunit, file=filerd, form='unformatted', status='old', iostat=ierr_open)
        if (ierr_open == 0) then
#ifdef oper   
          call system("sync")
          call system("ls -l -L "//filerd)
          call system("date")
          call system("sleep 1")
#endif
          exit
        else
          print *,"Input file ",filerd," is not ready: wait - sleep 60 s..."
          isleep = isleep + 1
#ifndef oper
          if (isleep == istop) then
            print*, 'Input file no.', nun, ' not found. Program stopped!'
            no_input = 1
            exit
          endif
#endif
          call system ("sleep 60")
        endif
      enddo

    endif

#ifdef mpi
      call mpi_bcast (no_input, 1, mpi_integer, 0, comm, error)
#endif
    if (no_input == 1) then
#ifdef mpi
        call mpi_finalize (error)
#endif
      stop
    endif

!  read descriptor records

    if (myid == 0) then
      if (init_flag == 1) then
        read(iunit) nfdr_local
        read(iunit) pdr_local
        if (any(nfdr_local(:) /= nfdr0(:)).or.any(pdr_local(:) /= pdr0(:))) then
          ierr_read = 1
          print *, ' Header parameter (nfdr, pdr) of ',trim(filerd), &
 ' input file not coincide with defined parameters.   STOP the program'
#ifdef mpi
          call mpi_abort(comm, error)
#endif
          stop
        endif
      else
        read(iunit)
        read(iunit)
      endif
    endif

#ifdef mpi
      call mpi_bcast (ierr_read, 1, mpi_integer, 0, comm, error)
#endif
    if (ierr_read == 1) then
#ifdef mpi
        call mpi_finalize (error)
#endif
      stop
    endif

!  Physiographical parameters changing in time

    if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
    if (init_flag == 1.or.nlclimate) call disperse (gfield, lai)

    if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
    if (init_flag == 1.or.nlclimate) call disperse (gfield, fveg)

    if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
    if (init_flag == 1.or.nlclimate) call disperse (gfield, rgm)

    if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
    if (init_flag == 1.or.nlclimate) call disperse (gfield, rgq)

    if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
    if (init_flag == 1.or.nlclimate) call disperse (gfield, fice_rd)
    if (init_flag == 1) fice(:,:)=fice_rd(:,:)

    if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
    if (init_flag == 1.or.nlclimate) call disperse (gfield, iceth_rd)
    if (init_flag == 1) iceth(:,:)=iceth_rd(:,:)

    if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
    if (init_flag == 1) call disperse (gfield, albedo)

    if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
    if (init_flag == 1) call disperse (gfield, emisg1)

    if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
    if (init_flag == 1) call disperse (gfield, emisg2)

! Prognostic cloud and precipitation variables at the surface

    if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
    if (init_flag == 1) call disperse (gfield, cloudt)

    if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
    if (init_flag == 1) call disperse (gfield, prectot)

!    if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
!    if (init_flag == 1) call disperse (gfield, precconv)

    if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
    if (init_flag == 1) call disperse (gfield, precsolid)

! Prognostic surface and soil/sea fields

    if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
    if (init_flag == 1) call disperse (gfield, tskin)

    if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
    if (init_flag == 1) call disperse (gfield, tg_surf)

    if (init_flag == 1) then
      do k = 1, nlevg
        if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
        call disperse (gfield, tg(1,1,k))
      enddo
    else
      do k = 1, nlevg-1
        if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
      enddo
      if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
      call disperse (gfield, tgclim)
    endif

    if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
    if (init_flag == 1) then
      call disperse (gfield, qskin)
    endif

    if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
    if (init_flag == 1) call disperse (gfield, qg_surf)

    if (init_flag == 1) then
      do k = 1, nlevg
        if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
        call disperse (gfield, qg(1,1,k))
      enddo
    else
      do k = 1, nlevg-1
        if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
      enddo
      if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
      call disperse (gfield, qgclim)
    endif

    if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
    if (init_flag == 1) call disperse (gfield, fice_soil_surf)

    do k = 1, nlevg
      if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
      call disperse (gfield, fice_soil(1,1,k))
    enddo

    if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
    if (init_flag == 1) call disperse (gfield, snow)

    do k = 1, nlev_snow
      if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
      call disperse (gfield, snow_lev(1,1,k))
    enddo

    do k = 1, nlev_snow
      if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
      call disperse (gfield, snow_t(1,1,k))
    enddo

    do k = 1, nlev_snow
      if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
      call disperse (gfield, snow_fice(1,1,k))
    enddo

    do k = 1, nlev_snow
      if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
      call disperse (gfield, snow_age(1,1,k))
    enddo

    do k = 1, nlev_snow
      if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
      call disperse (gfield, snow_melt_age(1,1,k))
    enddo

    do k = 1, nlev_snow
      if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
      call disperse (gfield, snow_dens(1,1,k))
    enddo

    if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
    call disperse (gfield, alsn)

    if (myid == 0) then
      close (iunit)
!     call system("date")
      write (*,'(2a)') ' Read file ',trim(filerd)
      write (*,*)
    endif

    if (nun >= 3) flag_pochva_par_change=1

 return
 end subroutine rdmhf_soil
!###############################################################################################################
 subroutine rd_param_const

! Reads from additional input file
! all constant (in time) model physiographical parameters

    use mod_moloch
    implicit none

    integer :: iunit=11, ierr_open, nlon_local, nlat_local, nlevg_local, nst_local, nvt_local, ierr=0, ird, jklev
    character (len=30) :: filerd="model_param_constant.bin"
    real :: dlon_local, dlat_local, x0_local, y0_local, alon0_local, alat0_local
    integer, dimension(gnlon,gnlat) :: igfield
    real, dimension(1) :: z1
    integer :: comm, error

#ifdef mpi
      include 'mpif.h'
      comm = mpi_comm_world
#endif

 if (myid == 0) then

   write(*,*)

   open (iunit, file=trim(filerd), status='old', form='unformatted', iostat=ierr_open)
   if (ierr_open /= 0) then
     print *,"Error in the open input data file,",trim(filerd)," constant (in time) physiographycal parameters: terminating."
#ifdef mpi
       call mpi_abort(comm, error)
#endif
     stop
   endif

   read (iunit) nlon_local, nlat_local, nlevg_local, dlon_local, dlat_local, x0_local, y0_local, alon0_local, alat0_local, &
 nst_local, nst_local

   alat0_local=alat0_local-dlat_local*0.5

   if (nlon_local /= gnlon) ierr=ierr+1
   if (nlat_local /= gnlat) ierr=ierr+1
   if (nlevg_local /= nlevg) ierr=ierr+1
   if (dlon_local /= pdr0(2)) ierr=ierr+1
   if (dlat_local /= pdr0(1)) ierr=ierr+1
   if (x0_local /= pdr0(39)) ierr=ierr+1
   if (y0_local /= pdr0(38)) ierr=ierr+1
   if (alon0_local /= pdr0(5)) ierr=ierr+1
   if (alat0_local /= pdr0(4)) ierr=ierr+1

   if (ierr /= 0) then
     print *,"Error in header parameters in input file ,",trim(filerd),", not coincident with defined parameters"
     print *,"Model nlon, nlat, nlevg, dlon, dlat, x0, y0, alon0, alat0 :", &
 gnlon, gnlat, nlevg, pdr0(2), pdr0(1), pdr0(39), pdr0(38), pdr0(5), pdr0(4)
     print *,"Read nlon, nlat, nlevg, dlon, dlat, x0, y0, alon0, alat0 :", &
 nlon_local, nlat_local, nlevg_local, dlon_local, dlat_local, x0_local, y0_local, alon0_local, alat0_local
     print *,"   STOP"
#ifdef mpi
       call mpi_abort(comm, error)
#endif
     stop
   endif

 endif ! myid == 0

!  orography*g and land sea mask

 if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, fmask)

 if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, phig)
 phig(:,:) = phig(:,:)*g

 if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield) ! Orography variance
 call disperse (gfield, orogvar)

!  soil types

 do ird = 1, nst+1
   if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, suolo(1,1,ird))
 enddo

!  vegetation types

 do ird = 1, nvt+1
   if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, vegeta(1,1,ird))
 enddo

!  physical soil parameters at soil levels

 do jklev = 1, nlevg
   if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, qsoil_max(1,1,jklev))
 enddo
 do jklev = 1, nlevg
   if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, qsoil_min(1,1,jklev))
 enddo
 do jklev = 1, nlevg
   if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, c_soil(1,1,jklev))
 enddo
 do jklev = 1, nlevg
   if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, rho_soil(1,1,jklev))
 enddo
 do jklev = 1, nlevg
   if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, psi_soil(1,1,jklev))
 enddo
 do jklev = 1, nlevg
   if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, kw_soil(1,1,jklev))
 enddo
 do jklev = 1, nlevg
   if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, par_b_soil(1,1,jklev))
 enddo
do jklev = 1, nlevg
   if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, par_c_soil(1,1,jklev))
 enddo
 do jklev = 1, nlevg
   if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, qsrel_wilt(1,1,jklev))
 enddo
 do jklev = 1, nlevg
   if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, qsrel_ref(1,1,jklev))
 enddo

!  radiation parameters of soil surface

 if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, albedo_soil_dry)
 if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, albedo_soil_wet)
 if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, emiss1_soil_dry)
 if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, emiss1_soil_wet)
 if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, emiss2_soil_dry)
 if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, emiss2_soil_wet)

!  physical parameters of vegetation

 if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, proot)
 if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield) ! veg_roughness

!  radiation parameters of vegetation surface

 if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, albedo_veg)
 if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, emiss1_veg)
 if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, emiss2_veg)

! top and bottom soil parameters

 if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, water_table_depth)
 if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, tg_bottom)
 if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, qg_rel_bottom)
 if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, qg_rel_surf_approx)

!  soil level index of bottom contidion for heat transport and water transport schemes

 if (myid == 0) call rrec2_int (iunit, gnlon, gnlat, igfield)
 call disperse_int (igfield, ind_lev_soil_h_bottom)
 if (myid == 0) call rrec2_int (iunit, gnlon, gnlat, igfield)
 call disperse_int (igfield, ind_lev_soil_w_bottom)

! snow_dirt: snow "dirtibility", that is weight (proportion 0-1) of dirt
! growth of snow surface due to vegetation waste, aerosol deposition, etc., 
! it is used in snow albedo definition.

 if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, snow_dirt)

!  maximum values of LAI in used dataset (for a normalization)

 if (myid == 0) read (iunit) lai_max
 if (myid == 0) z1(1)=lai_max
#ifdef mpi
   call mpi_bcast (z1, 1, mpi_real   , 0, comm, error)
#endif
 if (myid /= 0) lai_max=z1(1)

 if (myid == 0) then
   close (iunit)
   write(*,*) " Read file ",trim(filerd)
   write(*,*)
 endif

 return
 end subroutine rd_param_const
!###############################################################################################################
      subroutine disperse (gfield, lfield)

!  root distributes local domains to other processes
!  other processes receive their sub-domains from root

      use mod_moloch, only : nprocsx, nprocsy, myid, nlon, nlat, gnlon, gnlat
      implicit none
      integer                          :: jpr, count, comm, error, tag=0
      integer                          :: sx, ex, sy, ey
      real(4), dimension(nlon,nlat)    :: lfield
      real(4), dimension(gnlon,gnlat)  :: gfield

#ifdef mpi

        include 'mpif.h'
        integer                        :: status(mpi_status_size)

        count = nlon*nlat
        comm = mpi_comm_world
        call mpi_barrier(comm, error)

        if (myid == 0) then
          lfield(1:nlon,1:nlat) = gfield(1:nlon,1:nlat)
          do jpr = 1, nprocsx*nprocsy-1
            sx = 1+(jpr/nprocsy)*(nlon-2)
            ex = sx+nlon-1
            sy = 1+(jpr-(jpr/nprocsy)*nprocsy)*(nlat-2)
            ey = sy+nlat-1
            call mpi_send (gfield(sx:ex,sy:ey), count, mpi_real, jpr, tag, comm, error)
          enddo
        else
          call mpi_recv (lfield, count, mpi_real, 0, tag, comm, status, error)
        endif

        call mpi_barrier(comm, error)

#else

        lfield(1:nlon,1:nlat) = gfield(1:nlon,1:nlat)

#endif

      return
      end subroutine disperse
!###############################################################################################################
      subroutine disperse_int (igfield, ilfield)

!  root distributes local inegeter domains to other processes
!  other processes receive their sub-domains from root

      use mod_moloch, only : nprocsx, nprocsy, myid, nlon, nlat, gnlon, gnlat
      implicit none
      integer                          :: jpr, count, comm, error, tag=0
      integer                          :: sx, ex, sy, ey
      integer, dimension(nlon,nlat)    :: ilfield
      integer, dimension(gnlon,gnlat)  :: igfield

#ifdef mpi

        include 'mpif.h'
        integer                          :: status(mpi_status_size)

        count = nlon*nlat
        comm = mpi_comm_world
        call mpi_barrier(comm, error)

        if (myid == 0) then
          ilfield(1:nlon,1:nlat) = igfield(1:nlon,1:nlat)
          do jpr = 1, nprocsx*nprocsy-1
            sx = 1+(jpr/nprocsy)*(nlon-2)
            ex = sx+nlon-1
            sy = 1+(jpr-(jpr/nprocsy)*nprocsy)*(nlat-2)
            ey = sy+nlat-1
            call mpi_send (igfield(sx:ex,sy:ey), count, mpi_real, jpr, tag, comm, error)
          enddo
        else
          call mpi_recv (ilfield, count, mpi_real, 0, tag, comm, status, error)
        endif

        call mpi_barrier(comm, error)

#else

        ilfield(1:nlon,1:nlat) = igfield(1:nlon,1:nlat)

#endif

      return
      end subroutine disperse_int
!###############################################################################################################
    subroutine rrec2 (kunit, nlon, nlat, vect)

    implicit none

    integer :: kunit, nlon, nlat
    real, dimension(nlon, nlat) :: vect

    read(kunit) vect(1:nlon,1:nlat)

    return
    end subroutine rrec2
!###############################################################################################################
    subroutine rrec2_int (kunit, nlon, nlat, ivect)

    implicit none

    integer :: kunit, nlon, nlat
    integer, dimension(nlon, nlat) :: ivect

    read(kunit) ivect(1:nlon,1:nlat)

    return
    end subroutine rrec2_int
!###############################################################################################################
    subroutine rrec2_old (kunit, nlon, nlat, vect)

    implicit none

    integer :: kunit, nlon, nlat, jlat
    real, dimension(nlon, nlat) :: vect

    do jlat=1,nlat
      read(kunit) vect(1:nlon,jlat)
    enddo

    return
    end subroutine rrec2_old
!###############################################################################################################
    subroutine rrec2_int_old (kunit, nlon, nlat, ivect)

    implicit none

    integer :: kunit, nlon, nlat, jlat
    integer, dimension(nlon, nlat) :: ivect

    do jlat=1,nlat
      read(kunit) ivect(1:nlon,jlat)
    enddo

    return
    end subroutine rrec2_int_old
!###############################################################################################################
    subroutine wrec2 (kunit, nlon, nlat, vect)

    implicit none

    integer :: kunit, nlon, nlat
    real, dimension(nlon, nlat) :: vect

    write(kunit) vect(1:nlon,1:nlat)

!    call flush (kunit)

    return
    end subroutine wrec2
!###############################################################################################################
    subroutine wrec2r (kunit, nlon, nlat, vect)

    use mod_moloch, only: mhfr

    implicit none

! mhfr - write every mhfr points

    integer :: kunit, nlon, nlat
    real, dimension(nlon, nlat) :: vect

    write(kunit) vect(1:nlon:mhfr,1:nlat:mhfr)

!    call flush (kunit)

    return
    end subroutine wrec2r
!###############################################################################################################
    subroutine wrec2r_old (kunit, nlon, nlat, vect)

    use mod_moloch, only: mhfr

    implicit none

! mhfr - write every mhfr points

    integer :: kunit, nlon, nlat, jlat
    real, dimension(nlon, nlat) :: vect

    do jlat=1,nlat,mhfr
      write(kunit) vect(1:nlon:mhfr, jlat)
    enddo

!    call flush (kunit)

    return
    end subroutine wrec2r_old
!###############################################################################################################
    subroutine wrec2_int (kunit, nlon, nlat, ivect)

    implicit none

    integer :: kunit, nlon, nlat
    integer, dimension(nlon, nlat) :: ivect

    write(kunit) ivect(1:nlon,1:nlat)

!    call flush (kunit)

    return
    end subroutine wrec2_int
!###############################################################################################################
    subroutine wrec2r_int (kunit, nlon, nlat, ivect)

    use mod_moloch, only: mhfr

    implicit none

! mhfr - write every mhfr points

    integer :: kunit, nlon, nlat
    integer, dimension(nlon, nlat) :: ivect

    write(kunit) ivect(1:nlon:mhfr,1:nlat:mhfr)

!    call flush (kunit)

    return
    end subroutine wrec2r_int
!###############################################################################################################
    subroutine wrec2r_int_old (kunit, nlon, nlat, ivect)

    use mod_moloch, only: mhfr

    implicit none

! mhfr - write every mhfr points

    integer :: kunit, nlon, nlat, jlat
    integer, dimension(nlon, nlat) :: ivect

    do jlat=1,nlat,mhfr
      write(kunit) ivect(1:nlon:mhfr, jlat)
    enddo

!    call flush (kunit)

    return
    end subroutine wrec2r_int_old
!###############################################################################################################
    subroutine filt2d (p, anu2)

    use mod_moloch, only: nlon, nlat, nlonm1, nlatm1, ip_e, ip_n, ip_s, ip_w
    implicit none
    integer :: ierr, comm, tag=1, isend, irecv
    integer jlon, jlat
    real p(nlon,nlat), p2(nlon,nlat), anu2

#ifdef mpi

      include 'mpif.h'
      integer :: status(mpi_status_size)

      comm  = mpi_comm_world

!------------------------
!  Update all ghostlines
!------------------------

      call mpi_isend (p(2:nlonm1,nlatm1), nlon-2, mpi_real, ip_n, tag, comm, isend, ierr)
      call mpi_irecv (p(2:nlonm1,1     ), nlon-2, mpi_real, ip_s, tag, comm, irecv, ierr)
      call mpi_wait  (isend, status, ierr)
      call mpi_wait  (irecv, status, ierr)
      call mpi_isend (p(2:nlonm1,2     ), nlon-2, mpi_real, ip_s, tag, comm, isend, ierr)
      call mpi_irecv (p(2:nlonm1,nlat  ), nlon-2, mpi_real, ip_n, tag, comm, irecv, ierr)
      call mpi_wait  (isend, status, ierr)
      call mpi_wait  (irecv, status, ierr)
!      call u_ghost (p(2:nlonm1,nlatm1), ip_n, p(2:nlonm1,1   ), ip_s, nlon-2)
!      call u_ghost (p(2:nlonm1,2     ), ip_s, p(2:nlonm1,nlat), ip_n, nlon-2)
      call u_ghost (p(nlonm1,:), ip_e, p(1   ,:), ip_w, nlat)
      call u_ghost (p(2     ,:), ip_w, p(nlon,:), ip_e, nlat)

#endif

!-----------------------
!  Horizontal diffusion
!-----------------------

    do jlat = 2, nlatm1
    do jlon = 2, nlonm1
    p2(jlon,jlat) = .125 *(p(jlon,jlat-1)+p(jlon-1,jlat)+p(jlon+1,jlat)+p(jlon,jlat+1))-.5*p(jlon,jlat)
    enddo
    enddo

    do jlat = 2, nlatm1
    do jlon = 2, nlonm1
    p(jlon,jlat) = p(jlon,jlat) + anu2*p2(jlon,jlat)
    enddo
    enddo

    return
    end subroutine filt2d
!###############################################################################################################
    subroutine aver (p)

    use mod_moloch, only: nlon, nlat, nlonm1, nlatm1, ip_e, ip_n, ip_s, ip_w
    implicit none
    integer :: ierr, comm, tag=1, isend, irecv
    integer jlon, jlat
    real p(nlon,nlat), p2(nlon,nlat)

#ifdef mpi

      include 'mpif.h'
      integer :: status(mpi_status_size)

      comm  = mpi_comm_world

!------------------------
!  Update all ghostlines
!------------------------

      call mpi_isend (p(2:nlonm1,nlatm1), nlon-2, mpi_real, ip_n, tag, comm, isend, ierr)
      call mpi_irecv (p(2:nlonm1,1     ), nlon-2, mpi_real, ip_s, tag, comm, irecv, ierr)
      call mpi_wait  (isend, status, ierr)
      call mpi_wait  (irecv, status, ierr)
      call mpi_isend (p(2:nlonm1,2     ), nlon-2, mpi_real, ip_s, tag, comm, isend, ierr)
      call mpi_irecv (p(2:nlonm1,nlat  ), nlon-2, mpi_real, ip_n, tag, comm, irecv, ierr)
      call mpi_wait  (isend, status, ierr)
      call mpi_wait  (irecv, status, ierr)
!      call u_ghost (p(2:nlonm1,nlatm1), ip_n, p(2:nlonm1,1   ), ip_s, nlon-2)
!      call u_ghost (p(2:nlonm1,2     ), ip_s, p(2:nlonm1,nlat), ip_n, nlon-2)
      call u_ghost (p(nlonm1,:), ip_e, p(1   ,:), ip_w, nlat)
      call u_ghost (p(2     ,:), ip_w, p(nlon,:), ip_e, nlat)

#endif

!-------------------------------------------
!  Average preserving area-integrated values
!-------------------------------------------

    do jlat = 2, nlatm1
    do jlon = 2, nlonm1
    p2(jlon,jlat) = .125*(p(jlon,jlat-1)+p(jlon-1,jlat)+p(jlon+1,jlat)+p(jlon,jlat+1)) +.25*p(jlon,jlat) + &
                    .0625*(p(jlon-1,jlat-1)+p(jlon-1,jlat+1)+p(jlon+1,jlat-1)+p(jlon+1,jlat+1))
    enddo
    enddo

    do jlat = 2, nlatm1
    do jlon = 2, nlonm1
    p(jlon,jlat) = p2(jlon,jlat)
    enddo
    enddo

    return
    end subroutine aver
!###############################################################################################################
    subroutine filt3d (p, anu2)

    use mod_moloch, only: nlon, nlat, nlev, nlonm1, nlatm1, ip_e, ip_n, ip_s, ip_w
    implicit none
    integer :: ierr, comm, tag=1, isend, irecv
    integer jlon, jlat, k
    real p(nlon,nlat,nlev), p2(nlon,nlat), anu2

#ifdef mpi

      include 'mpif.h'
      integer :: status(mpi_status_size)

      comm  = mpi_comm_world

!------------------------
!  Update all ghostlines
!------------------------

      do k = 1, nlev
        call mpi_isend (p(2:nlonm1,nlatm1,k), nlon-2, mpi_real, ip_n, tag, comm, isend, ierr)
        call mpi_irecv (p(2:nlonm1,1     ,k), nlon-2, mpi_real, ip_s, tag, comm, irecv, ierr)
        call mpi_wait  (isend, status, ierr)
        call mpi_wait  (irecv, status, ierr)
        call mpi_isend (p(2:nlonm1,2     ,k), nlon-2, mpi_real, ip_s, tag, comm, isend, ierr)
        call mpi_irecv (p(2:nlonm1,nlat  ,k), nlon-2, mpi_real, ip_n, tag, comm, irecv, ierr)
        call mpi_wait  (isend, status, ierr)
        call mpi_wait  (irecv, status, ierr)
!        call u_ghost (p(2:nlonm1,nlatm1,k), ip_n, p(2:nlonm1,1   ,k), ip_s, nlon-2)
!        call u_ghost (p(2:nlonm1,2     ,k), ip_s, p(2:nlonm1,nlat,k), ip_n, nlon-2)
        call u_ghost (p(nlonm1,:,k), ip_e, p(1   ,:,k), ip_w, nlat)
        call u_ghost (p(2     ,:,k), ip_w, p(nlon,:,k), ip_e, nlat)
      enddo

#endif

!-----------------------
!  Horizontal diffusion
!-----------------------

    do 50 k = 1, nlev
    do jlat = 2, nlatm1
    do jlon = 2, nlonm1
    p2(jlon,jlat) = .125 *(p(jlon,jlat-1,k)+p(jlon-1,jlat,k)+p(jlon+1,jlat,k)+p(jlon,jlat+1,k))-.5*p(jlon,jlat,k)
    enddo
    enddo

    do jlat = 2, nlatm1
    do jlon = 2, nlonm1
    p(jlon,jlat,k) = p(jlon,jlat,k) + anu2*p2(jlon,jlat)
    enddo
    enddo
50  continue

    return
    end subroutine filt3d
!###############################################################################################################
      subroutine relax (is, gammin, gammax, alpha)

!  Computes optimal relaxation coefficients for lateral
!  boundary conditions (Lehmann, MAP, 1993,1-14)
!  See the paper for more comments
!  Input:  is       width of boundary relaxation zone (power of 2)
!          gammin   minimal Courant number (c*dt/dx)
!          gammax   maximal Courant number
!  Output: alpha()  weight of externally specified values in the boundary
!                   zone (corresponding to optimal relax. coefficients)

      parameter (ismax=16, ismax2=32)
      dimension alpha(is)
      dimension p (0:ismax2), q (0:ismax2)
      dimension pp(0:ismax2), qq(0:ismax2)
      real my,kk,kdt2

      n = 1
      p(0) = 0.
      p(1) = 1.
      q(0) = 1.
      q(1) = 0.
      my = sqrt(gammax/gammin)
 1000 my = sqrt((my+1./my)/2.)

      do i=0,n+n
      pp(i) = 0.
      qq(i) = 0.
      enddo

      do i=0,n
      do j=0,n
      pp(i+j) = pp(i+j)+p(i)*p(j)+q(i)*q(j)
      qq(i+j) = qq(i+j)+2.*my*p(i)*q(j)
      enddo
      enddo

      do i=0,n+n
      p(i) = pp(i)
      q(i) = qq(i)
      enddo

      n = 2*n
      if(n.lt.is) goto 1000
        if(n.ne.is.and.is.ne.1) then
        write (*,*) "Caution: nbl is not a power of 2"
        write (*,*) "Stop in subr. relax"
        stop
        endif

      do i=n,1,-1
      kk=p(i)/q(i-1)

      do j=i,1,-1
      xxx=q(j)
      q(j)=p(j)-kk*q(j-1)
      p(j) = xxx
      enddo

      xxx = q(0)
      q(0) = p(0)
      p(0) = xxx
      kdt2 = kk*sqrt(gammin*gammax)
      alpha(i) = kdt2/(1.+kdt2)
      enddo

!  Remark: this alpha corresponds to the leapfrog scheme,
!  whereas kdt2 is independent of the integration scheme

      return
      end subroutine relax
!###############################################################################################################
 subroutine wrmhf_atm

! Writes the Model History File (MHF) with prognostic atmospheric variables

 use mod_moloch
 implicit none

 integer :: iunit=22, iunit_work=29, jklev
 character(len=30) file_out, amhf

   if (myid == 0) then

#ifdef oper
     write (amhf,'(i3.3)') imhf
     file_out='moloch_atm_'//trim(amhf)//'.mhf' 
     open (iunit, file=trim(file_out), form='unformatted', status='unknown')
#else
     open (iunit, file='moloch_atm.mhf', form='unformatted', status='unknown', position='append')
#endif

!  write descriptor records

     write(iunit) nfdrr
     write(iunit) pdrr

   endif

   do jklev=1,nlev
     call collect (p(1,1,jklev), gfield)
     if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)
   enddo
   do jklev=1,nlev
     call collect (u(1,1,jklev), gfield)
     if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)
   enddo
   do jklev=1,nlev
     call collect (v(1,1,jklev), gfield)
     if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)
   enddo
   do jklev=1,nlevp1
     call collect (w(1,1,jklev), gfield)
     if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)
   enddo
   do jklev=1,nlev
     call collect (t(1,1,jklev), gfield)
     if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)
   enddo
   do jklev=1,nlev
     call collect (q(1,1,jklev), gfield)
     if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)
   enddo
   do jklev=1,nlev
     call collect (qcw(1,1,jklev), gfield)
     if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)
   enddo
   do jklev=1,nlev
     call collect (qci(1,1,jklev), gfield)
     if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)
   enddo

   if (myid == 0) then
     flush (iunit)
     close (iunit)

#ifdef oper
     call system("sync")
!!!     call system("ls -l -L "//file_out)
!!!     call system("date")
     open (iunit_work, file=trim(file_out)//'.txt', status='unknown')
     write (iunit_work,'(2a)') trim(file_out),' is full and closed'
     close (iunit_work)
     call system("sync")
     print *,'Output written on file ', trim(file_out)
#else
     print *, 'File moloch_atm.mhf written'
#endif

   endif

 return
 end subroutine wrmhf_atm
!###############################################################################################################
 subroutine wrmhf_soil

! Writes the Model History File (MHF) with prognostic surface, sea and soil variables

 use mod_moloch
 implicit none

 integer :: iunit=23, iunit_work=29, jklev
 character(len=30) file_out, amhf
 real, dimension(nlon,nlat) :: zwork

   if (myid == 0) then

#ifdef oper
     write (amhf,'(i3.3)') imhf
     file_out='moloch_soil_'//trim(amhf)//'.mhf' 
     open (iunit, file=trim(file_out), form='unformatted', status='unknown')
#else
     open (iunit, file='moloch_soil.mhf', form='unformatted', status='unknown', position='append')
#endif

!  write descriptor records

     write(iunit) nfdrr
     write(iunit) pdrr

   endif

!  Physiographical parameters changing in time

   call collect (lai, gfield)
   if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)

   call collect (fveg, gfield)
   if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)

   call collect (rgm, gfield)
   if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)

   call collect (rgq, gfield)
   if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)

   call collect (fice, gfield)
   if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)

   call collect (iceth, gfield)
   if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)

   call collect (albedo, gfield)
   if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)

   call collect (emisg1, gfield)
   if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)

   call collect (emisg2, gfield)
   if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)

! Prognostic cloud and precipitation variables at the surface

   call collect (cloudt, gfield)
   if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)

   call collect (prectot, gfield)
   if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)

!   call collect (precconv, gfield)
!   if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)

   call collect (precsolid, gfield)
   if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)

! Prognostic surface and soil/sea fields

   call collect (tskin, gfield)
   if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)

   call collect (tg_surf, gfield)
   if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)

   do jklev=1,nlevg
     call collect (tg(1,1,jklev), gfield)
     if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)
   enddo

   call collect (qskin, gfield)
   if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)

   call collect (qg_surf, gfield)
   if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)

   do jklev=1,nlevg
     call collect (qg(1,1,jklev), gfield)
     if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)
   enddo

   call collect (fice_soil_surf, gfield)
   if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)

   do jklev=1,nlevg
     call collect (fice_soil(1,1,jklev), gfield)
     if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)
   enddo

   call collect (snow, gfield)
   if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)

   do jklev=1,nlev_snow
     call collect (snow_lev(1,1,jklev), gfield)
     if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)
   enddo

   do jklev=1,nlev_snow
     call collect (snow_t(1,1,jklev), gfield)
     if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)
   enddo

   do jklev=1,nlev_snow
     call collect (snow_fice(1,1,jklev), gfield)
     if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)
   enddo

   do jklev=1,nlev_snow
     call collect (snow_age(1,1,jklev), gfield)
     if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)
   enddo

   do jklev=1,nlev_snow
     call collect (snow_melt_age(1,1,jklev), gfield)
     if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)
   enddo

   do jklev=1,nlev_snow
     call collect (snow_dens(1,1,jklev), gfield)
     if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)
   enddo

   call collect (alsn, gfield)
   if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)

! Statistical (cumulated, average, max, min) surface variables

   call collect (cswfl, gfield)
   if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)

   call collect (clwfl, gfield)
   if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)

   call collect (cshflux, gfield)
   if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)

   call collect (clhflux, gfield)
   if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)

   call collect (t2min, gfield)
   if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)

   call collect (t2max, gfield)
   if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)

   call collect (ws10max, gfield)
   if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)

   call collect (runoff, gfield)
   if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)

   call collect (runoff_tot, gfield)
   if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)

   call collect (cwvflux, gfield)
   if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)

   call collect (cfl_heat_soil_bottom, gfield)
   if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)

   call collect (cfl_water_soil_bottom, gfield)
   if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)

   call collect (float(level_snowfall), gfield)
   if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)

! Writting of additional 2D fields

   zwork = 0.
   call collect (zwork, gfield)
   if(myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)

   if (myid == 0) then

     flush (iunit)
     close (iunit)

#ifdef oper
     call system("sync")
!!!     call system("ls -l -L "//file_out)
!!!     call system("date")
     open (iunit_work, file=trim(file_out)//'.txt', status='unknown')
     write (iunit_work,'(2a)') trim(file_out),' is full and closed'
     close (iunit_work)
     call system("sync")
     print *,'Output written on file ', trim(file_out)
#else
     print *, 'File moloch_soil.mhf written'
#endif

   endif

 return
 end subroutine wrmhf_soil
!###############################################################################################################
 subroutine wr_param_const

! Writes file with all constant (in time) model physiographical parameters in not full resolution
! It is used only in case of writting of mhf outputs in half resolution (mhfr=2)

    use mod_moloch
    implicit none

    integer :: iunit=21, iwr, jklev
    character (len=50) :: filewr="model_param_constant_not_full_res.bin"
    integer, dimension(gnlon,gnlat) :: igfield
    integer :: comm, error

#ifdef mpi
      include 'mpif.h'

      comm = mpi_comm_world
#endif

    if (myid == 0) then

      write(*,*)

      open (iunit, file=trim(filewr), form='unformatted', status='unknown')

! gnlon, gnlat, nlevg, dlon, dlat, x0, y0, alon0, alat0, nst, nvt
      write (iunit) nfdr0(2)/mhfr, nfdr0(3)/mhfr, nfdr0(15), pdr0(2)*float(mhfr), pdr0(1)*float(mhfr), &
 pdr0(39), pdr0(38), pdr0(5), pdr0(4)+pdr0(1)*0.5*float(mhfr), nst, nvt

    endif ! myid == 0

!  orography*g and land sea mask

 call collect (fmask, gfield)
 if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)

 call collect (phig/g, gfield)
 if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)

 call collect (orogvar, gfield)
 if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield) ! Orography variance

!  soil types

 do iwr = 1, nst+1
   call collect (suolo(1,1,iwr), gfield)
   if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)
 enddo

!  vegetation types

 do iwr = 1, nvt+1
   call collect (vegeta(1,1,iwr), gfield)
   if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)
 enddo

!  physical soil parameters at soil levels

 do jklev = 1, nlevg
   call collect (qsoil_max(1,1,jklev), gfield)
   if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)
 enddo
 do jklev = 1, nlevg
   call collect (qsoil_min(1,1,jklev), gfield)
   if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)
 enddo
 do jklev = 1, nlevg
   call collect (c_soil(1,1,jklev), gfield)
   if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)
 enddo
 do jklev = 1, nlevg
   call collect (rho_soil(1,1,jklev), gfield)
   if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)
 enddo
 do jklev = 1, nlevg
   call collect (psi_soil(1,1,jklev), gfield)
   if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)
 enddo
 do jklev = 1, nlevg
   call collect (kw_soil(1,1,jklev), gfield)
   if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)
 enddo
 do jklev = 1, nlevg
   call collect (par_b_soil(1,1,jklev), gfield)
   if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)
 enddo
 do jklev = 1, nlevg
   call collect (par_c_soil(1,1,jklev), gfield)
   if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)
 enddo
 do jklev = 1, nlevg
   call collect (qsrel_wilt(1,1,jklev), gfield)
   if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)
 enddo
 do jklev = 1, nlevg
   call collect (qsrel_ref(1,1,jklev), gfield)
   if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)
 enddo

!  radiation parameters of soil surface

 call collect (albedo_soil_dry, gfield)
 if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)
 call collect (albedo_soil_wet, gfield)
 if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)
 call collect (emiss1_soil_dry, gfield)
 if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)
 call collect (emiss1_soil_wet, gfield)
 if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)
 call collect (emiss2_soil_dry, gfield)
 if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)
 call collect (emiss2_soil_wet, gfield)
 if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)

!  physical parameters of vegetation

 call collect (proot, gfield)
 if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)
 call collect (veg_rough, gfield)
 if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield) ! veg_roughness

!  radiation parameters of vegetation surface

 call collect (albedo_veg, gfield)
 if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)
 call collect (emiss1_veg, gfield)
 if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)
 call collect (emiss2_veg, gfield)
 if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)

! soil top and bottom parameters

 call collect (water_table_depth, gfield)
 if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)
 call collect (tg_bottom, gfield)
 if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)
 call collect (qg_rel_bottom, gfield)
 if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)
 call collect (qg_rel_surf_approx, gfield)
 if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)

!  soil level index of bottom contidion for heat transport and water transport schemes

 call collect_int (ind_lev_soil_h_bottom, igfield)
 if (myid == 0) call wrec2r_int_old (iunit, gnlon, gnlat, igfield)
 call collect_int (ind_lev_soil_w_bottom, igfield)
 if (myid == 0) call wrec2r_int_old (iunit, gnlon, gnlat, igfield)

! snow_dirt: snow "dirtibility", that is weight (proportion 0-1) of dirt
! growth of snow surface due to vegetation waste, aerosol deposition, etc., 
! it is used in snow albedo definition.

 call collect (snow_dirt, gfield)
 if (myid == 0) call wrec2r (iunit, gnlon, gnlat, gfield)

!  maximum values of LAI in used dataset (for a normalization)

 if (myid == 0) write (iunit) lai_max

 if (myid == 0) then
   flush (iunit)
   close (iunit)
   print *, 'File ',trim(filewr),' written'
   print *
 endif

 return
 end subroutine wr_param_const
!###############################################################################################################
 subroutine wrmhf_atm_full_res

! Writes the Model History File (MHF) with prognostic atmospheric variables
! with full horizontal resolution

 use mod_moloch
 implicit none

 integer :: iunit=24, iunit_work=29, jklev
 character(len=50) file_out, aini, aterm

 if (myid == 0) then

   write (aini,'(i4.4,3i2.2)') nfdr(5:8)
   write (aterm,'(i3.3,2i2.2)') nfdr(10:12)
   file_out='moloch_full_res_atm_'//trim(aini)//'_'//trim(aterm)//'.mhf' 
   open (iunit, file=trim(file_out), form='unformatted', status='unknown')

!  write descriptor records

   write(iunit) nfdr
   write(iunit) pdr

 endif

   do jklev=1,nlev
     call collect (p(1,1,jklev), gfield)
     if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)
   enddo
   do jklev=1,nlev
     call collect (u(1,1,jklev), gfield)
     if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)
   enddo
   do jklev=1,nlev
     call collect (v(1,1,jklev), gfield)
     if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)
   enddo
   do jklev=1,nlevp1
     call collect (w(1,1,jklev), gfield)
     if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)
   enddo
   do jklev=1,nlev
     call collect (t(1,1,jklev), gfield)
     if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)
   enddo
   do jklev=1,nlev
     call collect (q(1,1,jklev), gfield)
     if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)
   enddo
   do jklev=1,nlev
     call collect (qcw(1,1,jklev), gfield)
     if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)
   enddo
   do jklev=1,nlev
     call collect (qci(1,1,jklev), gfield)
     if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)
   enddo

   if (myid == 0) then

     flush (iunit)
     close (iunit)

#ifdef oper
     call system("sync")
!!!     call system("ls -l -L "//file_out)
!!!     call system("date")
     open (iunit_work, file=trim(file_out)//'.txt', status='unknown')
     write (iunit_work,'(2a)') trim(file_out),' is full and closed'
     close (iunit_work)
     call system("sync")
     print *,'Output written on file ', trim(file_out)
#else
     print *, 'File moloch_atm.mhf written'
#endif

   endif

 return
 end subroutine wrmhf_atm_full_res
!###############################################################################################################
 subroutine wrmhf_soil_full_res

! Writes the Model History File (MHF) with prognostic surface, sea and soil variables
! with full horizontal resoltuion

! Attention: all cumulated, miximum, minimum and average fields are false!
! (prectot, precconv, precsolid, cswfl, clwfl, cshflux, clhflux, t2min, t2max,ws10max, runoff, runoff_tot, cwvflux, cfl_heat_soil_bottom, cfl_water_soil_bottom)

 use mod_moloch
 implicit none

 integer :: iunit=25, iunit_work=29, jklev
 character(len=50) file_out, aini, aterm
 real, dimension(nlon,nlat) :: zwork

 if (myid == 0) then

   write (aini,'(i4.4,3i2.2)') nfdr(5:8)
   write (aterm,'(i3.3,2i2.2)') nfdr(10:12)
   file_out='moloch_full_res_soil_'//trim(aini)//'_'//trim(aterm)//'.mhf' 
   open (iunit, file=trim(file_out), form='unformatted', status='unknown')

!  write descriptor records

     write(iunit) nfdr
     write(iunit) pdr

 endif

!  Physiographical parameters changing in time

   call collect (lai, gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

   call collect (fveg, gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

   call collect (rgm, gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

   call collect (rgq, gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

   call collect (fice, gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

   call collect (iceth, gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

   call collect (albedo, gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

   call collect (emisg1, gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

   call collect (emisg2, gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

! Prognostic cloud and precipitation variables at the surface

   call collect (cloudt, gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

   call collect (prectot, gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

!   call collect (precconv, gfield)
!   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

   call collect (precsolid, gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

! Prognostic surface and soil/sea fields

   call collect (tskin, gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

   call collect (tg_surf, gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

   do jklev=1,nlevg
     call collect (tg(1,1,jklev), gfield)
     if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)
   enddo

   call collect (qskin, gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

   call collect (qg_surf, gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

   do jklev=1,nlevg
     call collect (qg(1,1,jklev), gfield)
     if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)
   enddo

   call collect (fice_soil_surf, gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

   do jklev=1,nlevg
     call collect (fice_soil(1,1,jklev), gfield)
     if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)
   enddo

   call collect (snow, gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

   do jklev=1,nlev_snow
     call collect (snow_lev(1,1,jklev), gfield)
     if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)
   enddo

   do jklev=1,nlev_snow
     call collect (snow_t(1,1,jklev), gfield)
     if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)
   enddo

   do jklev=1,nlev_snow
     call collect (snow_fice(1,1,jklev), gfield)
     if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)
   enddo

   do jklev=1,nlev_snow
     call collect (snow_age(1,1,jklev), gfield)
     if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)
   enddo

   do jklev=1,nlev_snow
     call collect (snow_melt_age(1,1,jklev), gfield)
     if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)
   enddo

   do jklev=1,nlev_snow
     call collect (snow_dens(1,1,jklev), gfield)
     if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)
   enddo

   call collect (alsn, gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

! Statistical (cumulated, average, max, min) surface variables

   call collect (cswfl, gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

   call collect (clwfl, gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

   call collect (cshflux, gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

   call collect (clhflux, gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

   call collect (t2min, gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

   call collect (t2max, gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

   call collect (ws10max, gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

   call collect (runoff, gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

   call collect (runoff_tot, gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

   call collect (cwvflux, gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

   call collect (cfl_heat_soil_bottom, gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

   call collect (cfl_water_soil_bottom, gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

   call collect (float(level_snowfall), gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

! Writting of additional 2D fields

   zwork = 0.
   call collect (zwork, gfield)
   if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

   if (myid == 0) then

     flush (iunit)
     close (iunit)

#ifdef oper
     call system("sync")
!!!     call system("ls -l -L "//file_out)
!!!     call system("date")
     open (iunit_work, file=trim(file_out)//'.txt', status='unknown')
     write (iunit_work,'(2a)') trim(file_out),' is full and closed'
     close (iunit_work)
     call system("sync")
     print *,'Output written on file ', trim(file_out)
#else
     print *, 'File moloch_soil.mhf written'
#endif

   endif

 return
 end subroutine wrmhf_soil_full_res
!###############################################################################################################
      subroutine collect (lfield, gfield)

!  root process receives all sub-domains to reconstruct the global domain
!  other processes send theyr sub-domains to root

    use mod_moloch, only : nprocsx, nprocsy, myid, nlon, nlat, gnlon, gnlat
    implicit none
    integer                         :: comm, error, tag=0
    integer                         :: count, sx, ex, sy, ey, jpr
    real(4), dimension(nlon,nlat)   :: lfield
    real(4), dimension(gnlon,gnlat) :: gfield

#ifdef mpi

      include 'mpif.h'
      integer                         :: status(mpi_status_size)

      comm = mpi_comm_world

      if (myid == 0) then
        gfield(1:nlon,1:nlat) = lfield(1:nlon,1:nlat)
        do jpr = 1, nprocsx*nprocsy-1
          sx = 2+(jpr/nprocsy)*(nlon-2)
          ex = sx+nlon-3
          sy = 2+(jpr-(jpr/nprocsy)*nprocsy)*(nlat-2)
          ey = sy+nlat-3
          if (jpr.lt.nprocsy               ) sx = sx -1
          if (jpr.ge.(nprocsx-1)*nprocsy   ) ex = ex +1
          if (mod(jpr,nprocsy).eq.0        ) sy = sy -1
          if (mod(jpr,nprocsy).eq.nprocsy-1) ey = ey +1
          count = (ex-sx+1)*(ey-sy+1)
          call mpi_recv (gfield(sx:ex,sy:ey), count, mpi_real, jpr, jpr, comm, status, error)
        enddo
      else
        sx = 2
        ex = nlon-1
        sy = 2
        ey = nlat-1
        if (myid.lt.nprocsy               ) sx = sx -1
        if (myid.ge.(nprocsx-1)*nprocsy   ) ex = ex +1
        if (mod(myid,nprocsy).eq.0        ) sy = sy -1
        if (mod(myid,nprocsy).eq.nprocsy-1) ey = ey +1
        count = (ex-sx+1)*(ey-sy+1)
        call mpi_send (lfield(sx:ex,sy:ey), count, mpi_real, 0, myid, comm, error)
      endif

      call mpi_barrier(comm, error)

#else

      gfield(1:nlon,1:nlat) = lfield(1:nlon,1:nlat)

#endif

    return
    end subroutine collect
!###############################################################################################################
    subroutine collect_int (ilfield, igfield)

!  root process receives all sub-domains to reconstruct the global domain
!  other processes send theyr sub-domains to root

    use mod_moloch, only : nprocsx, nprocsy, myid, nlon, nlat, gnlon, gnlat
    implicit none
    integer                         :: comm, error, tag=0
    integer                         :: count, sx, ex, sy, ey, jpr
    integer, dimension(nlon,nlat)   :: ilfield
    integer, dimension(gnlon,gnlat) :: igfield

#ifdef mpi
      include 'mpif.h'
      integer                       :: status(mpi_status_size)

      comm = mpi_comm_world

      if (myid == 0) then
        igfield(1:nlon,1:nlat) = ilfield(1:nlon,1:nlat)
        do jpr = 1, nprocsx*nprocsy-1
          sx = 2+(jpr/nprocsy)*(nlon-2)
          ex = sx+nlon-3
          sy = 2+(jpr-(jpr/nprocsy)*nprocsy)*(nlat-2)
          ey = sy+nlat-3
          if (jpr.lt.nprocsy               ) sx = sx -1
          if (jpr.ge.(nprocsx-1)*nprocsy   ) ex = ex +1
          if (mod(jpr,nprocsy).eq.0        ) sy = sy -1
          if (mod(jpr,nprocsy).eq.nprocsy-1) ey = ey +1
          count = (ex-sx+1)*(ey-sy+1)
          call mpi_recv (igfield(sx:ex,sy:ey), count, mpi_integer, jpr, jpr, comm, status, error)
        enddo
      else
        sx = 2
        ex = nlon-1
        sy = 2
        ey = nlat-1
        if (myid.lt.nprocsy               ) sx = sx -1
        if (myid.ge.(nprocsx-1)*nprocsy   ) ex = ex +1
        if (mod(myid,nprocsy).eq.0        ) sy = sy -1
        if (mod(myid,nprocsy).eq.nprocsy-1) ey = ey +1
        count = (ex-sx+1)*(ey-sy+1)
        call mpi_send (ilfield(sx:ex,sy:ey), count, mpi_integer, 0, myid, comm, error)
      endif

      call mpi_barrier(comm, error)

#else

      igfield(1:nlon,1:nlat) = ilfield(1:nlon,1:nlat)

#endif

    return
    end subroutine collect_int
!###############################################################################################################
      subroutine runout (jstep)

      use mod_moloch, only : nlon, nlat, nlev, nprocsx, nprocsy, u, v, w, p, tskin, dtstep, ip_e, ip_n, ip_s, ip_w, myid
    integer :: ierr, comm

#ifdef mpi
      include 'mpif.h'
      real, dimension(12, 0:(nprocsx*nprocsy-1)) :: buffer
      real, dimension(12, 0:(nprocsx*nprocsy-1)) :: buffer_2d
      integer :: status(mpi_status_size)

      comm = mpi_comm_world
#endif

      umax=0.
      vmax=0.
      wmax=0.
      zpmin=1.e8
      zpmax=-1.e8
      ztsmin=1.e8
      ztsmax=-1.e8

      j1 = 2
      j2 = nlat-1
      if (ip_s.lt.0) j1 = 1
      if (ip_n.lt.0) j2 = nlat
      i1 = 2
      i2 = nlon-1
      if (ip_w.lt.0) i1 = 1
      if (ip_e.lt.0) i2 = nlon

      do jlat = j1, j2
        do jlon = i1, i2
          if (p(jlon,jlat,1).lt.zpmin) then
            zpmin = p(jlon,jlat,1)
            ipmin = jlon + (myid/nprocsy)*(nlon-2)
            jpmin = jlat + (myid-(myid/nprocsy)*nprocsy)*(nlat-2)
          endif
          if (p(jlon,jlat,1).gt.zpmax) then
            zpmax = p(jlon,jlat,1)
            ipmax = jlon + (myid/nprocsy)*(nlon-2)
            jpmax = jlat + (myid-(myid/nprocsy)*nprocsy)*(nlat-2)
          endif
          if (tskin(jlon,jlat).lt.ztsmin) then
            ztsmin = tskin(jlon,jlat)
            itsmin = jlon + (myid/nprocsy)*(nlon-2)
            jtsmin = jlat + (myid-(myid/nprocsy)*nprocsy)*(nlat-2)
          endif
          if (tskin(jlon,jlat).gt.ztsmax) then
            ztsmax = tskin(jlon,jlat)
            itsmax = jlon + (myid/nprocsy)*(nlon-2)
            jtsmax = jlat + (myid-(myid/nprocsy)*nprocsy)*(nlat-2)
          endif
          do jklev = 1, nlev
            if (abs(u(jlon,jlat,jklev)).gt.umax) then
            umax = abs(u(jlon,jlat,jklev))
            iumax = jlon + (myid/nprocsy)*(nlon-2)
            jumax = jlat + (myid-(myid/nprocsy)*nprocsy)*(nlat-2)
            kumax = jklev
            endif
            if (abs(v(jlon,jlat,jklev)).gt.vmax) then
            vmax = abs(v(jlon,jlat,jklev))
            ivmax= jlon + (myid/nprocsy)*(nlon-2)
            jvmax= jlat + (myid-(myid/nprocsy)*nprocsy)*(nlat-2)
            kvmax= jklev
            endif
            if (abs(w(jlon,jlat,jklev)).gt.wmax) then
            wmax = abs(w(jlon,jlat,jklev))
            iwmax= jlon + (myid/nprocsy)*(nlon-2)
            jwmax= jlat + (myid-(myid/nprocsy)*nprocsy)*(nlat-2)
            kwmax= jklev
            endif
          enddo
        enddo
      enddo

#ifdef mpi

        buffer( 1,myid) = umax
        buffer( 2,myid) = iumax
        buffer( 3,myid) = jumax
        buffer( 4,myid) = kumax
        buffer( 5,myid) = vmax
        buffer( 6,myid) = ivmax
        buffer( 7,myid) = jvmax
        buffer( 8,myid) = kvmax
        buffer( 9,myid) = wmax
        buffer(10,myid) = iwmax
        buffer(11,myid) = jwmax
        buffer(12,myid) = kwmax
    
        buffer_2d( 1,myid) = zpmax
        buffer_2d( 2,myid) = ipmax
        buffer_2d( 3,myid) = jpmax
        buffer_2d( 4,myid) = zpmin
        buffer_2d( 5,myid) = ipmin
        buffer_2d( 6,myid) = jpmin
        buffer_2d( 7,myid) = ztsmax
        buffer_2d( 8,myid) = itsmax
        buffer_2d( 9,myid) = jtsmax
        buffer_2d(10,myid) = ztsmin
        buffer_2d(11,myid) = itsmin
        buffer_2d(12,myid) = jtsmin

!  gather info on subdomain local maxima to root processor

        if (myid.gt.0) then
          call mpi_send (buffer(:,myid), 12, mpi_real, 0, myid, comm, ierr)
          call mpi_send (buffer_2d(:,myid), 12, mpi_real, 0, myid, comm, ierr)
        else
          do ip = 1, nprocsx*nprocsy-1
            call mpi_recv (buffer(:,ip), 12, mpi_real, ip, ip, comm, status, ierr)
            call mpi_recv (buffer_2d(:,ip), 12, mpi_real, ip, ip, comm, status, ierr)
          enddo
        endif

!  calculate global minimum

!        call mpi_reduce (zpmin, zdummy, 1, mpi_real, mpi_min, 0, comm, ierr)
!        zpmin = zdummy

        if (myid == 0) then

          zpmin = 1.e8
          ztsmin = 1.e8
          do i = 0, nprocsx*nprocsy-1
            if (buffer_2d( 4,i).lt.zpmin) then
              zpmin = buffer_2d( 4,i)
              ipmin = INT(buffer_2d( 5,i))
              jwmin = INT(buffer_2d( 6,i))
            endif
            if (buffer_2d(10,i).lt.ztsmin) then
              ztsmin = buffer_2d(10,i)
              itsmin = INT(buffer_2d(11,i))
              jtsmin = INT(buffer_2d(12,i))
            endif
          enddo

!  compute global maximum

          umax = 0.
          vmax = 0.
          wmax = 0.
          zpmax = -1.e8
          ztsmax = -1.e8
          do i = 0, nprocsx*nprocsy-1
            if (buffer(1,i).gt.umax) then
              umax  = buffer(1,i)
              iumax = INT(buffer(2,i))
              jumax = INT(buffer(3,i))
              kumax = INT(buffer(4,i))
            endif
            if (buffer(5,i).gt.vmax) then
              vmax  = buffer(5,i)
              ivmax = INT(buffer(6,i))
              jvmax = INT(buffer(7,i))
              kvmax = INT(buffer(8,i))
            endif
            if (buffer(9,i).gt.wmax) then
              wmax  = buffer(9,i)
              iwmax = INT(buffer(10,i))
              jwmax = INT(buffer(11,i))
              kwmax = INT(buffer(12,i))
            endif
            if (buffer_2d( 1,i).gt.zpmax) then
              zpmax = buffer_2d( 1,i)
              ipmax = INT(buffer_2d( 2,i))
              jwmax = INT(buffer_2d( 3,i))
            endif
            if (buffer_2d( 7,i).gt.ztsmax) then
              ztsmax = buffer_2d( 7,i)
              itsmax = INT(buffer_2d( 8,i))
              jtsmax = INT(buffer_2d( 9,i))
            endif
          enddo

        endif

#endif

    if (myid == 0) then

      print*
      write(*,'(a,i5,a,f8.3)') ' jstep',jstep, "   -   integration time (hours):", jstep*dtstep/3600.
      print 1010, umax, vmax, wmax, zpmin/100.
      print 1011, iumax, jumax, kumax, ivmax, jvmax, kvmax, iwmax, jwmax, kwmax
      write (*,'(2(a,f7.2))') ' Max tskin =    ',ztsmax-273.15,'   Min tskin = ',ztsmin-273.15
      write (*,'(2(a,2i4))') ' at gridpoints:  ',itsmax,jtsmax,'              ',itsmin,jtsmin

    endif

 1010 format(' Max abs u,v,w =',3f13.3,'      min ps =',f8.2)
 1011 format(' at gridpoints:    ',3i4,2x,3i4,2x,3i4)
      return
      end subroutine runout
!###############################################################################################################
      subroutine surflayer(kstep)

!  Vertical diffusion in the atmospheric surface layer

      use mod_moloch, only : nlon, nlat, nlonm1, nlatm1, dz, h, g, phig, gzita, bzita, tetavs, fsnow, fmask,  &
                             tskin, qskin, ps, u, v, tetav, rd, rv, cpd, rdrcp, roscdm, roscdt, qstar, tstar, &
                             ustar, rich, kturb_surf_m, kturb_surf_h, kturb_surf_q, &
                             hflux, qflux, rgm, rgmd, rgq, std_lev_atm, n_std_lev_atm, u_std_lev,       &
                             v_std_lev, t_std_lev, q_std_lev, n_std_lev_sl, t, q, cvm, &
                             nstep_sl_filter, kturb_surf_m_mem, kturb_surf_h_mem, kturb_surf_q_mem, myid

      integer :: kstep
      real, parameter :: zak=.4, zgam=16., zaholt=1., zbholt=2./3., zcholt=5., zdholt=0.35
      real, dimension(n_std_lev_atm) :: psim_sl, psih_sl, psiq_sl
      real :: zlev_bottom
      real, dimension(nlon,nlat) :: kturb_surf_m_filter, kturb_surf_h_filter, kturb_surf_q_filter

!  Businger functions

      psium(zz1,zz2) = log( (1.+zz1)**2*(1.+zz1**2)/((1.+zz2)**2*(1.+zz2**2)) ) -2.*(atan(zz1)-atan(zz2))
      psiuh(zz1,zz2) = 2.*log((1.+zz1**2)/(1.+zz2**2))

!  Holtslag functions

      psism(zz1,zz2) = -zaholt*zz1-zbholt*(zz1-zcholt/zdholt)*exp(-zdholt*zz1)  &
                       +zaholt*zz2+zbholt*(zz2-zcholt/zdholt)*exp(-zdholt*zz2)
      psish(zz1,zz2) = -(1.+2./3.*zaholt*zz1)**1.5-zbholt*(zz1-zcholt/zdholt)*exp(-zdholt*zz1) &
                       +(1.+2./3.*zaholt*zz2)**1.5+zbholt*(zz2-zcholt/zdholt)*exp(-zdholt*zz2)

      zalp0 = log(1.e5)
      zza0  = -h*bzita(.5*dz)*log(1.-.5*dz/h)
      zep   = rv/rd-1.

      do 100 jlat = 2, nlatm1
      do 100 jlon = 2, nlonm1

!--------------------------------------------------------
!  Turbulent fluxes at the ground (positive upward)
!--------------------------------------------------------

      zza   = phig(jlon,jlat)*(gzita(.5*dz)-1.)/g + zza0 + rgm(jlon,jlat)   ! zita=0 is located at z=rgm
      zua   = 0.5*(u(jlon-1,jlat,1)+u(jlon,jlat,1))
      zva   = 0.5*(v(jlon,jlat,1)+v(jlon,jlat+1,1))
      zmod2 = zua**2 + zva**2 + 0.07
      zmod  = sqrt(zmod2)

!  Virtual potential temperature computed with skin temperature

      zconvg = exp(rdrcp*(zalp0-log(ps(jlon,jlat))))*(1.+zep*qskin(jlon,jlat))
      tetavs(jlon,jlat) = tskin(jlon,jlat)*zconvg

!  Bulk Richardson number

      zri = 2.*zza*g*(tetav(jlon,jlat,1)-tetavs(jlon,jlat))/((tetav(jlon,jlat,1)+tetavs(jlon,jlat))*zmod2)
      rich(jlon,jlat,1) = min (zri, 500.)

      if (fmask(jlon,jlat).ge..5) then

!  Computation of Charnok roughness

        zchar = 5.e-4
        zcoch1 = .0185*zak**2/g*zmod2
        do jiter = 1, 5
        zcoch2 = zza/zchar
        zchar  = zcoch1/alog(zcoch2)**2
        enddo
        zrgm = zchar

!  Roughness lengths over sea interpolated logarithmically between 0.<Ri<0.25 with Large&Pond values

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
        zrgm = fmask(jlon,jlat)*zsea+(1.-fmask(jlon,jlat))*rgm(jlon,jlat)*  &
               (1.-0.5*fsnow(jlon,jlat)*0.4/(0.4+rgm(jlon,jlat)))
        zrgt = fmask(jlon,jlat)*zsea+(1.-fmask(jlon,jlat))*rgq(jlon,jlat)*  &
               (1.-0.5*fsnow(jlon,jlat)*0.4/(0.4+rgq(jlon,jlat)))
        zrgq = zrgt
      endif

      rgmd(jlon,jlat) = zrgm

      zalzam = alog (zza/zrgm)
      zalzat = alog (zza/zrgt)
      zalzaq = alog (zza/zrgq)

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
        zpsiq = zpsih  ! because rgt=rgq in the stable case
        ustar(jlon,jlat) = zak*zmod/(zalzam-zpsim)
        tstar(jlon,jlat) = zak/(zalzat-zpsih)*(tetav(jlon,jlat,1)-tetavs(jlon,jlat))
        qstar(jlon,jlat) = zak/(zalzaq-zpsiq)*(q    (jlon,jlat,1)-qskin (jlon,jlat))
        zcdm = ustar(jlon,jlat)**2/zmod
        zcdt = ustar(jlon,jlat)*zak/(zalzat-zpsih)
        zcdq = zcdt

           do k = 1, n_std_lev_sl(jlon,jlat)
           psim_sl(k) = psism(zal*(std_lev_atm(k)+zrgm)/zza, zal*zrgm/zza)
           psih_sl(k) = psish(zal*(std_lev_atm(k)+zrgt)/zza, zal*zrgt/zza)
           psiq_sl(k) = psih_sl(k) ! because rgt=rgq in the stable case
           enddo

      else                       ! Businger functions

        zpsim = 0.
        zpsih = 0.
        do jiter = 1, 4
        zal = zri*(zalzam-zpsim)**2/(zalzat-zpsih)   ! za/L from eq (6.47) of Businger
        zz1 = (1.-zgam*zal)**.25
        zx2 = (1.-zgam*zal*zrgm/zza)**.25
        zy2 = (1.-zgam*zal*zrgt/zza)**.25
        zpsim = psium(zz1,zx2)
        zpsih = psiuh(zz1,zy2)
        enddo
        zz2 = (1.-zgam*zal*zrgq/zza)**.25
        zpsiq = psiuh(zz1,zz2)
        ustar(jlon,jlat) = zak*zmod/(zalzam-zpsim)
        tstar(jlon,jlat) = zak/(zalzat-zpsih)*(tetav(jlon,jlat,1)-tetavs(jlon,jlat))
        qstar(jlon,jlat) = zak/(zalzaq-zpsiq)*(q    (jlon,jlat,1)-qskin (jlon,jlat))
        zcdm = ustar(jlon,jlat)**2/zmod
        zcdt = ustar(jlon,jlat)*zak/(zalzat-zpsih)
        zcdq = ustar(jlon,jlat)*zak/(zalzaq-zpsiq)

           do k = 1, n_std_lev_sl(jlon,jlat)
           zz1 = (1.-zgam*zal*(std_lev_atm(k)+zrgm)/zza)**.25
           psim_sl(k) = psium(zz1, zx2)
           zz1 = (1.-zgam*zal*(std_lev_atm(k)+zrgt)/zza)**.25
           psih_sl(k) = psiuh(zz1, zy2)
           zz1 = (1.-zgam*zal*(std_lev_atm(k)+zrgq)/zza)**.25
           psiq_sl(k) = psiuh(zz1, zz2)
           enddo

      endif

!--------------------------------------------------------------------------------------------------
! Output variables
!--------------------------------------------------------------------------------------------------

    zlev_bottom = phig(jlon,jlat)*(gzita(.5*dz)-1.)/g + zza0

    kturb_surf_m(jlon,jlat) = abs( -zcdm*zlev_bottom )
    kturb_surf_h(jlon,jlat) = abs( (-zcdt/zconvg)*zlev_bottom )
    kturb_surf_q(jlon,jlat) = abs( -zcdq*zlev_bottom)

!--------------------------------------------------------------------------------------------------

! "Interpolation" at standard levels (m) above the surface, inside the surface layer, for output purposes

    do k = 1, n_std_lev_sl(jlon,jlat)
    zuv = ustar(jlon,jlat)/zak*(log((std_lev_atm(k)+zrgm)/zrgm)-psim_sl(k))
    u_std_lev(jlon,jlat,k) = zuv*zua/zmod
    v_std_lev(jlon,jlat,k) = zuv*zva/zmod
    q_std_lev(jlon,jlat,k) = qskin(jlon,jlat) +qstar(jlon,jlat)/zak*(log((std_lev_atm(k)+zrgq)/zrgq)-psiq_sl(k))
    zconv2 = zconvg*(1.+zep*q_std_lev(jlon,jlat,k))/(1.+zep*qskin(jlon,jlat))
    t_std_lev(jlon,jlat,k) = (tetavs(jlon,jlat)+tstar(jlon,jlat)/zak*(log((std_lev_atm(k)+zrgt)/zrgt)-psih_sl(k)))/zconv2

!  Fix over high mountains...

    zglin = max (0., min((phig(jlon,jlat)-11000.)/15000., 1.))
    u_std_lev(jlon,jlat,k) = (1.-zglin)*u_std_lev(jlon,jlat,k) + zglin*u(jlon,jlat,1)
    v_std_lev(jlon,jlat,k) = (1.-zglin)*v_std_lev(jlon,jlat,k) + zglin*v(jlon,jlat,1)
    if (t_std_lev(jlon,jlat,k).lt.t(jlon,jlat,1)) then
    zglin = max (0., min((phig(jlon,jlat)-10000.)/25000., 1.))
    t_std_lev(jlon,jlat,k) = (1.-zglin)*t_std_lev(jlon,jlat,k) + zglin*t(jlon,jlat,1)
    q_std_lev(jlon,jlat,k) = (1.-zglin)*q_std_lev(jlon,jlat,k) + zglin*q(jlon,jlat,1)
    endif
    enddo

!!  FIXFIX low wind
!
!    if (fmask(jlon,jlat).lt.0.1) then
!      if (tskin(jlon,jlat).gt.t_std_lev(jlon,jlat,1)+.1) then
!      zzz = u_std_lev(jlon,jlat,2)**2 + v_std_lev(jlon,jlat,2)**2
!      zzz = max ((25.-zzz)/25., 0.)
!      zzz = min (zzz, .8)
!      t_std_lev(jlon,jlat,1) = t_std_lev(jlon,jlat,1) + min (zzz*(tskin(jlon,jlat)-t_std_lev(jlon,jlat,1)), 2.0)
!      endif
!    endif

!--------------------------------------------------------
!  Exchange coefficients at the surface
!--------------------------------------------------------

!      zros = ps(jlon,jlat)/(rd*tskin(jlon,jlat)*(1.+zep*qskin(jlon,jlat)))
!      roscdm(jlon,jlat) = .75*roscdm(jlon,jlat) +.25*zros*zcdm
!      roscdt(jlon,jlat) = .75*roscdt(jlon,jlat) +.25*zros*zcdt ! time average to stabilize tskin computation
!!      roscdm(jlon,jlat) = .5*roscdm(jlon,jlat) +.5*zros*zcdm
!!      roscdt(jlon,jlat) = .5*roscdt(jlon,jlat) +.5*zros*zcdt  ! time average to stabilize tskin computation
!      hflux (jlon,jlat) = -roscdt(jlon,jlat)*cpd/zconvg       ! sensible heat flux coefficient
!      cvm(jlon,jlat,1)  = roscdm(jlon,jlat)*zza/zros

 100  continue

!--------------------------------------------------------
!  Exchange coefficients at the surface
!--------------------------------------------------------

! Time filter for turbulent coefficient (heat and water vapour) in the surface layer

 if (nstep_sl_filter > 1) then

   if (kstep > nstep_sl_filter) then
     do n = 1,nstep_sl_filter-1
       kturb_surf_m_mem(:,:,n) = kturb_surf_m_mem(:,:,n+1)
       kturb_surf_h_mem(:,:,n) = kturb_surf_h_mem(:,:,n+1)
       kturb_surf_q_mem(:,:,n) = kturb_surf_q_mem(:,:,n+1)
     enddo
     kturb_surf_m_mem(:,:,nstep_sl_filter) = kturb_surf_m(:,:)
     kturb_surf_h_mem(:,:,nstep_sl_filter) = kturb_surf_h(:,:)
     kturb_surf_q_mem(:,:,nstep_sl_filter) = kturb_surf_q(:,:)
     kturb_surf_m_filter(:,:) = 0.
     kturb_surf_h_filter(:,:) = 0.
     kturb_surf_q_filter(:,:) = 0.
     do n = 1,nstep_sl_filter
       kturb_surf_m_filter(:,:) = kturb_surf_m_filter(:,:) + kturb_surf_m_mem(:,:,n)
       kturb_surf_h_filter(:,:) = kturb_surf_h_filter(:,:) + kturb_surf_h_mem(:,:,n)
       kturb_surf_q_filter(:,:) = kturb_surf_q_filter(:,:) + kturb_surf_q_mem(:,:,n)
     enddo
     kturb_surf_m_filter(:,:) = kturb_surf_m_filter(:,:) / float(nstep_sl_filter)
     kturb_surf_h_filter(:,:) = kturb_surf_h_filter(:,:) / float(nstep_sl_filter)
     kturb_surf_q_filter(:,:) = kturb_surf_q_filter(:,:) / float(nstep_sl_filter)
     kturb_surf_m(:,:) = kturb_surf_m_filter(:,:)
     kturb_surf_h(:,:) = kturb_surf_h_filter(:,:)
     kturb_surf_q(:,:) = kturb_surf_q_filter(:,:)
   else
     kturb_surf_m_mem(:,:,kstep) = kturb_surf_m(:,:)
     kturb_surf_h_mem(:,:,kstep) = kturb_surf_h(:,:)
     kturb_surf_q_mem(:,:,kstep) = kturb_surf_q(:,:)
   endif

 endif

! Surface fluxes

 do jlat = 2, nlatm1
 do jlon = 2, nlonm1
   zros = ps(jlon,jlat)/(rd*tskin(jlon,jlat)*(1.+zep*qskin(jlon,jlat)))
   zlev_bottom = phig(jlon,jlat)*(gzita(.5*dz)-1.)/g + zza0
   zconvg = exp(rdrcp*(zalp0-log(ps(jlon,jlat))))*(1.+zep*qskin(jlon,jlat))
   roscdm(jlon,jlat) = zros*kturb_surf_m(jlon,jlat)/zlev_bottom
   roscdt(jlon,jlat) = zros*kturb_surf_h(jlon,jlat)/zlev_bottom*zconvg
   cvm(jlon,jlat,1) = kturb_surf_m(jlon,jlat)
   hflux(jlon,jlat) = kturb_surf_h(jlon,jlat)*(zros*cpd)/zlev_bottom*(tetavs(jlon,jlat)-tetav(jlon,jlat,1))
   qflux(jlon,jlat) = kturb_surf_q(jlon,jlat)*zros/zlev_bottom*(qskin(jlon,jlat)-q(jlon,jlat,1))
 enddo
 enddo

      return
      end subroutine surflayer
!###############################################################################################################
      subroutine tofd

      use mod_moloch, only : nlon, nlat, ntop, nlonm1, nlatm1, u, v, dtstep, orogvar, zeta

      ctofd  = 5.e-9        ! turbulent orographic form drag coefficient

      do jklev = 1, ntop
      do jlat = 2, nlatm1
      do jlon = 2, nlonm1
      zilev = zeta(jlon,jlat,jklev)+.01
      if (zilev.lt.3000.) then  ! Turbulent orographic form drag - Beljaars et al, QJRMS, 2004
      zvar = orogvar(jlon,jlat)**2
      ztofd = ctofd*sqrt(u(jlon,jlat,jklev)**2+v(jlon,jlat,jklev)**2)*zvar/zilev**1.2
      if (zilev.gt.500.) ztofd = ztofd*exp(-(zilev/1500.)**1.5)
      ztofd = 1./(1.+dtstep*ztofd)
      u(jlon,jlat,jklev) = u(jlon,jlat,jklev)*ztofd
      v(jlon,jlat,jklev) = v(jlon,jlat,jklev)*ztofd
      endif
      enddo
      enddo
      enddo

      return
      end subroutine tofd
!###############################################################################################################
      subroutine orogdrag (jstep)

! This parametrization has showed bad results

! Computes orographic gravity wave and block drag

      use mod_moloch, only: nlon, nlat, nlev, nlevp1, nlonm1, nlatm1, phig, u, v, p, tvirt, hxt, hx, hy, zeta, rd, &
                            g, dx, dy, dz, dlon, dlat, dtstep, fmz, ip_e, ip_n, ip_s, ip_w, bvf, rich, orogvar
      implicit none
      real, dimension(nlon,nlat,nlevp1) :: fmx, fmy
      real, dimension(nlon,nlat)        :: zorxx, zoryy
      real, dimension(nlev)             :: zwe
      real, dimension(nlevp1)           :: fup, fdw
      integer, dimension(nlon,nlat)     :: iormask, lcrit
      integer jstep, jlon, jlat, jklev
      real zuav, zvav, zn, zro, zsprod, zrich, zut, zvt, hra, znfa, zinx, ziny, zdc, ordragx, ordragy
      real zorx, zory, zuinc, zvinc, bldragx, bldragy, zorogstd, zdz, zridz

! zdragc is a general amplitude coefficient;
! cext and zrichden determine the extinction coefficient zwe
! cext is the inverse of a typical atmospheric depth of absorption

      real, parameter :: zdragc=0.03, cext=2.e-3, zrichden=1.25, zrichcrit=-.5
      real, parameter :: zbldragc=.5e-4  ! general amplitude coefficient for block drag

      save zorxx, zoryy, iormask

! znfa: normaliz. factor to keep the average drag in x nearly independent of model resolution
! 0.40+0.28 is dlon+dlat for globo 31 km res. - expon. 2./3. found empirically

      znfa = ((0.40+0.28)/(dlon+dlat))**(2./3.)

!-------------------------------------
!  1- Initialization (first call only)
!-------------------------------------

      if (jstep.eq.1) then
      iormask = 0

! Second derivatives of orography
! Only crests, not valleys, are considered - a const. is subtracted to avoid "noise" over sea or hills
! A function of subgrid orography variance is added
! Threshold values for drag application (using iormask)

      do jlat = 2, nlatm1
      do jlon = 2, nlonm1
      hra = max(hxt(jlat)*dlon/dlat, 0.1) ! to limit the singularity at poles
      zorxx(jlon,jlat) = (phig(jlon-1,jlat)+phig(jlon+1,jlat)-2.*phig(jlon,jlat))/g/sqrt(hra)
      zoryy(jlon,jlat) = (phig(jlon,jlat-1)+phig(jlon,jlat+1)-2.*phig(jlon,jlat))/g
      zorxx(jlon,jlat) = min(0., zorxx(jlon,jlat)+10.)
      zoryy(jlon,jlat) = min(0., zoryy(jlon,jlat)+10.)
      zorxx(jlon,jlat) = znfa*(zorxx(jlon,jlat)-0.35*orogvar(jlon,jlat))
      zoryy(jlon,jlat) = znfa*(zoryy(jlon,jlat)-0.35*orogvar(jlon,jlat))
        if(zorxx(jlon,jlat).lt.-50..or.zoryy(jlon,jlat).lt.-50.) then
        iormask(jlon,jlat) = 1
        else
        iormask(jlon,jlat) = 0
        endif
      enddo
      enddo
      endif  ! end of computations at the first instant only

#ifdef mpi
        call u_ghost (u(nlonm1,:,:), ip_e, u(1,:,:), ip_w, nlat*nlev)
        call u_ghost (v(:,2,:), ip_s, v(:,nlat,:), ip_n, nlon*nlev)
#endif

      fmx = 0.
      fmy = 0.
      lcrit = nlev

      do jlat = 2, nlatm1
      do jlon = 2, nlonm1
      if (iormask(jlon,jlat).eq.1.and.bvf(jlon,jlat).gt.1.e-6) then

!-------------------------------------------------------------------
!  2- Computation of momentum flux at the surface (ordragx, ordragy)
!-------------------------------------------------------------------

! All computations on T points
! Comput. of u and v as averages over the three lowest levels
! Brunt-Vaisala freq. (squared) bvf is computed in vdiff over the 3 lowest layers)

      zuav = (u(jlon-1,jlat,1)+u(jlon,jlat,1)+u(jlon-1,jlat,2)+u(jlon,jlat,2)+u(jlon-1,jlat,3)+u(jlon,jlat,3))/6.
      zvav = (v(jlon,jlat,1)+v(jlon,jlat+1,1)+v(jlon,jlat,2)+v(jlon,jlat+1,2)+v(jlon,jlat,3)+v(jlon,jlat+1,3))/6.
      zn   = sqrt(bvf(jlon,jlat))         ! bvf is the square of Brunt-Vaisala frequency near surface
      if (zn.gt.1.2e-2) zn = 1.2e-2       ! to limit g.w. drag for high static stability
      zro = p(jlon,jlat,1)/(rd*tvirt(jlon,jlat,1)) ! density at first level
      ordragx = zdragc*zro*zn*zuav*zorxx(jlon,jlat)
      ordragy = zdragc*zro*zn*zvav*zoryy(jlon,jlat)

! Computation of the critical level, defined as the level where the wind speed becomes zero
! or the wind vector becomes perpendicular to the orographic drag (that can have a different
! direction with respect to the surface (3 level-mean) wind vector).
! The critical level is defined where the scalar product between the wind and the drag
! becomes (for the first time going upward) null or positive.
! lcrit contains the index of the half-level below the level where the above condition is verified.
! In case of no critical level, lcrit = nlev

      do jklev = 4, nlev
      zut = .5*(u(jlon-1,jlat,jklev) + u(jlon,jlat,  jklev))
      zvt = .5*(v(jlon,  jlat,jklev) + v(jlon,jlat+1,jklev))
      zsprod = zut*ordragx + zvt*ordragy
        if (zsprod.ge.0..or.rich(jlon,jlat,jklev).lt.0.) then
        lcrit(jlon,jlat) = min (lcrit(jlon,jlat), jklev)
        exit
        endif
      enddo

!-----------------------------------------------------------
!  3- Computation of momentum flux (fmx, fmy) on half-levels
!-----------------------------------------------------------

! Gravity wave drag is applied up to the critical level 
! The moist Richardson number computed in vdiff is used
! Extintion coefficient (zwe) for upward and downward momentum fluxes on levels, 
! depending on Richardson number ( f(k+1)-f(k)=-zwe*.5*(f(k+1)+f(k)) )

      do jklev = 1, lcrit(jlon,jlat)-1
      zrich = .5*(rich(jlon,jlat,jklev+1)+rich(jlon,jlat,jklev))
      zrich = min (zrich, 5.)
      zrich = max (zrich, zrichcrit)
      zwe(jklev) = cext*dz/(fmz(jlon,jlat,jklev)*(2.*zrich + zrichden))
      zwe(jklev) = min (zwe(jklev), 1.) ! zwe = 2 implica assorbimento totale in un layer
      enddo
      if (rich(jlon,jlat,lcrit(jlon,jlat)).lt.0.) zwe(lcrit(jlon,jlat)-1) = 1.

! Upward and downward momentum flux profiles on half levels

      fup = 0.
      fdw = 0.
      fup(1) = 1.
      do jklev = 2, lcrit(jlon,jlat)
      fup(jklev) = fup(jklev-1)*(1.-.5*zwe(jklev-1))/(1.+.5*zwe(jklev-1))
      enddo
      fdw(lcrit(jlon,jlat)) = fup(lcrit(jlon,jlat))  ! reflection at critical level
      do jklev = lcrit(jlon,jlat)-1, 1, -1
      fdw(jklev) = fdw(jklev+1)*(1.-.5*zwe(jklev))/(1.+.5*zwe(jklev))
      enddo

      do jklev = 1, lcrit(jlon,jlat)   ! x and y total momentum flux (on halh-levels)
      fmx(jlon,jlat,jklev) = ordragx*(fup(jklev)-fdw(jklev))
      fmy(jlon,jlat,jklev) = ordragy*(fup(jklev)-fdw(jklev))
      enddo

      endif ! cond. on iormask and bvf
      enddo
      enddo

!--------------------------------------------------------
!  4- Divergence of vertical momentum flux and u,v update
!--------------------------------------------------------

#ifdef mpi
        call u_ghost (fmx(2,:     ,:), ip_w, fmx(nlon,:,:), ip_e, nlat*nlev)
        call u_ghost (fmy(:,nlatm1,:), ip_n, fmy(:   ,1,:), ip_s, nlon*nlev)
#endif

      do jlat = 2, nlatm1
      do jlon = 2, nlonm1
      if (iormask(jlon,jlat).eq.1.and.bvf(jlon,jlat).gt.1.e-6) then

!  Fluxes averaged on u,v points

      do jklev = 1, lcrit(jlon,jlat)-1
      zdc = .5*dtstep*rd*tvirt(jlon,jlat,jklev)*fmz(jlon,jlat,jklev)/(p(jlon,jlat,jklev)*dz) !.5*dt/dz/ro
      zinx = zdc*((fmx(jlon,jlat,jklev+1)+fmx(jlon+1,jlat,jklev+1))-  &
                  (fmx(jlon,jlat,jklev  )+fmx(jlon+1,jlat,jklev  )))
      ziny = zdc*((fmy(jlon,jlat-1,jklev+1)+fmy(jlon,jlat,jklev+1))-  &
                  (fmy(jlon,jlat-1,jklev  )+fmy(jlon,jlat,jklev  )))
      u(jlon,jlat,jklev) = u(jlon,jlat,jklev) - zinx
      v(jlon,jlat,jklev) = v(jlon,jlat,jklev) - ziny
      enddo

      endif ! cond. on iormask and bvf
      enddo
      enddo

!----------------------------------------
!  5- Computation of drag due to blocking
!----------------------------------------

! Comput. of u e v as averages over the three lowest levels
! Brunt-Vaisala freq. squared bvf is computed in vdiff over the lowest 3 layers

      do jlat = 2, nlatm1
      do jlon = 2, nlonm1
      if (bvf(jlon,jlat).gt.1.e-5) then

      zorogstd = orogvar(jlon,jlat)+.01
      zn = sqrt(bvf(jlon,jlat)) ! bvf is the square of Brunt-Vaisala freq.
      zro = p(jlon,jlat,1)/(rd*tvirt(jlon,jlat,1))

! U-component

      zuav = (u(jlon,jlat,1) + u(jlon,jlat,2) + u(jlon,jlat,3))/3.
      zuav = sign (sqrt(abs(zuav)), zuav)     ! empirico...per diminuire l'effetto per alte wind speed
      zorx = znfa*hx(jlon,jlat)*dx*hxt(jlat)
      bldragx = .35*zuav*zorogstd                                    ! st. dev. contribution
      if (zuav*zorx.gt.0..and.zn*abs(zorx).gt.1.) bldragx = bldragx +sign(zuav*zorx, zuav) ! applied only at upslope flows
      bldragx = zbldragc*zn*zro*bldragx

! V-component

      zvav = (v(jlon,jlat,1) + v(jlon,jlat,2) + v(jlon,jlat,3))/3.
      zvav = sign (sqrt(abs(zvav)), zvav)
      zory = znfa*hy(jlon,jlat)*dy
      bldragy = .35*zvav*zorogstd
      if (zvav*zory.gt.0..and.zn*abs(zory).gt.1.) bldragy = bldragy +sign(zvav*zory, zvav)
      bldragy = zbldragc*zn*zro*bldragy

! Application of drag at lower levels, decaying with height above ground

      if (abs(bldragx).gt.1.e-6.or.abs(bldragy).gt.1.e-6) then

      do jklev = 1, nlev/2
      zdz = (zeta(jlon,jlat,jklev)-zeta(jlon,jlat,1))/zorogstd
        if (zdz.lt.1.0) then  ! linear decay at the lowest levels
        zridz = 1.-.3*zdz 
        else
        zridz = .7*exp(-13.*bvf(jlon,jlat)*(zdz-1.))
        endif
      if (u(jlon,jlat,jklev)*bldragx.gt.0.) then
      zuinc = dtstep*bldragx*zridz
      zuinc = sign (min(.7*abs(u(jlon,jlat,jklev)), abs(zuinc)), zuinc) ! limitation of increment
      u(jlon,jlat,jklev) = u(jlon,jlat,jklev) - zuinc
      endif
      if (v(jlon,jlat,jklev)*bldragy.gt.0.) then
      zvinc = dtstep*bldragy*zridz
      zvinc = sign (min(.7*abs(v(jlon,jlat,jklev)), abs(zvinc)), zvinc)
      v(jlon,jlat,jklev) = v(jlon,jlat,jklev) - zvinc
      endif
      enddo

      endif

      endif  ! condition on bvf
      enddo
      enddo

      return
      end subroutine orogdrag
!###############################################################################################################
      subroutine vdiff(kstep)  ! Theta-Pai version

!  E-L closure - vertical derivatives only
!  vertical diffusion of momentum, specific humidity and virtual pot. temperature
!  computes mixing length, TKE, heat and specific humidity fluxes at the ground (fluxes are positive upward)

      use mod_moloch, only : nlon, nlat, nlev, nlevp1, ntop, nlonm1, nlatm1, dz, h, g, phig, gzita, bzita,  &
                             u, v, w, tetav, pai, p, q, tke, mlz, rd, cpd, ylwv, fmz, fmzh, tkemin, dx, dy, &
                             qsatw, ustar, tstar, roscdm, roscdt, hflux, qflux, tetavs, dtstep, tvirt, t,   &
                             qcw, qci, qskin, chm, cvm, prandtl, eps, rich, bvf, myid

      real, dimension(nlevp1)      :: zita0
      real, dimension(nlon)        :: za1m, za1t, za1q, za1mm1
      real, dimension(nlon,nlevp1) :: zhlev, zam, zamm1, zcm, zcmm1, zah, zch, dlogthe
      real, dimension(nlon,nlev)   :: zeu, zev, zet, zeq, zfe, zfw
      real, parameter              :: zak=.4, zumltun=1.0
      integer                      :: kstep

      chm = 0.
      cvm(:,:,2:nlevp1) = 0.

      zce = .17             ! ustar**2 = zce*tke
      zmlmax = 100.         ! Blackadar mixing length maximum value
      zmlcut = sqrt(dx*dy)  ! mixing length maximum value

!  Array used to compute half-level heights

      do jk = 1, nlev
      zitah = (jk-1)*dz
      zita0(jk) = -h*bzita(zitah)*log(1.-zitah/h)
      enddo

      zzc0 = -dtstep/dz**2

!--------------------------------------------
!  Loop in latitude to the end of the routine
!--------------------------------------------

      do 1000 jlat = 2, nlatm1

!  height of half-levels above orography, eq. potential and virtual temperature

      do jklev = 1, ntop+1
      zitah = (jklev-1)*dz
      do jlon = 2, nlonm1
      zhlev(jlon,jklev) = phig(jlon,jlat)*(gzita(zitah)-1.)/g + zita0(jklev)
      enddo
      enddo

!  Computation of Dtheta/theta (dlogthe).
!  If air is saturated, moist saturates stability (Durran&Klemp, 1982 ) is considered

      do jk = 2, ntop               ! loop on half levels
      do jlon = 2, nlonm1
      r_up = qsatw(jlon,jlat,jk  )
      r_do = qsatw(jlon,jlat,jk-1)
      dthd = 2.*(tetav(jlon,jlat,jk)-tetav(jlon,jlat,jk-1))/(tetav(jlon,jlat,jk)+tetav(jlon,jlat,jk-1)) ! dry
      r_av = 0.5*(r_up + r_do)
      t_av = 0.5*(t(jlon,jlat,jk) + t(jlon,jlat,jk-1))
      theta_up = t(jlon,jlat,jk  )/pai(jlon,jlat,jk  )
      theta_do = t(jlon,jlat,jk-1)/pai(jlon,jlat,jk-1)
      zaa = (1. + ylwv*r_av/(rd*t_av))/(1. + eps*ylwv**2*r_av/(cpd*rd*t_av**2))
      dthm = zaa*((theta_up-theta_do)*2./(theta_up + theta_do) + ylwv/(cpd*t_av)*(r_up - r_do)) &
             - q(jlon,jlat,jk  ) - qcw(jlon,jlat,jk  ) - qci(jlon,jlat,jk  )                    &
             + q(jlon,jlat,jk-1) + qcw(jlon,jlat,jk-1) + qci(jlon,jlat,jk-1)

      zrhm = 0.45*(q(jlon,jlat,jk)/qsatw(jlon,jlat,jk)) + 0.55*(q(jlon,jlat,jk-1)/qsatw(jlon,jlat,jk-1))
      zcof = max(-24.5 + 25.*zrhm, 0.)    ! zcof=0. for rh<=0.98, zcof=0.5 for rh=1
      zcof = min(zcof, .85)
      dlogthe(jlon,jk) = zcof*dthm + (1.-zcof)*dthd  ! effective stability near saturation
      enddo
      enddo

!--------------------------------------------------------
!  TKE lower boundary condition
!--------------------------------------------------------

      do jlon = 2, nlonm1

      tke(jlon,jlat,1) = ustar(jlon,jlat)**2/zce+tkemin

      zwstar = 0.
      jkubl  = 1  ! top of unstable boundary layer
      ztflux = -ustar(jlon,jlat)*tstar(jlon,jlat)

      if (ztflux.gt.0.) then  ! in case of unstable surface layer
        do jk = 2, ntop   ! loop on half levels
        if (dlogthe(jlon,jk).gt.0..or.jk.eq.ntop) then
        jkubl = jk
        exit
        endif
        enddo
      zwstar = (g*zhlev(jlon,jkubl)*ztflux/tetavs(jlon,jlat))**(1./3.)
      tke(jlon,jlat,1) = tke(jlon,jlat,1) + .4*zwstar**2
      endif

      enddo

!---------------------------------------------------------------------------
!  TKE equation
!  Mixing length and Prandtl number on W points
!---------------------------------------------------------------------------

      do jklev = 2, ntop
      do jlon = 2, nlonm1
      zrdz  = fmzh(jlon,jlat,jklev)/dz
      zdudz = .5*(u(jlon,jlat,jklev)+u(jlon-1,jlat,jklev)-u(jlon,jlat,jklev-1)-u(jlon-1,jlat,jklev-1))*zrdz
      zdvdz = .5*(v(jlon,jlat,jklev)+v(jlon,jlat+1,jklev)-v(jlon,jlat,jklev-1)-v(jlon,jlat+1,jklev-1))*zrdz
      zdwdz = .5*(w(jlon,jlat,jklev+1)-w(jlon,jlat,jklev-1))*zrdz
      zbuoy = g*dlogthe(jlon,jklev)*zrdz
      zrich = min (zbuoy/(zdudz**2+zdvdz**2+1.e-6), 500.)
      rich(jlon,jlat,jklev) = zrich

! Brunt-Vaisala averaged near surface

      if(jklev==2) bvf(jlon,jlat) = 0.2*zbuoy
      if(jklev==3) bvf(jlon,jlat) = bvf(jlon,jlat) + 0.35*zbuoy
      if(jklev==4) bvf(jlon,jlat) = bvf(jlon,jlat) + 0.45*zbuoy

      if (zbuoy > 0.) then
      prandtl(jlon,jlat,jklev) = .75 + 1.5*zrich + .5*zrich*sqrt(zrich) ! average (after Anderson, Bound. L. Met., 2009)
      else
      prandtl(jlon,jlat,jklev) = .75
      endif
      zshear = zdudz**2 + zdvdz**2 + 2.*zdwdz**2
      zp = mlz(jlon,jlat,jklev)*(zshear-zbuoy/prandtl(jlon,jlat,jklev))
      zd = zce*tke(jlon,jlat,jklev)/(mlz(jlon,jlat,jklev)+.01)
      zdtke = sqrt(zce*tke(jlon,jlat,jklev))*dtstep*(zp-zd)
      zdtke = max (zdtke, -.8*tke(jlon,jlat,jklev))    ! remove no more than 80% of TKE in a single timestep
      tke(jlon,jlat,jklev) = tke(jlon,jlat,jklev) + zdtke
      enddo
      enddo

!  Frictional heating

      do jklev = 1, ntop-1
      do jlon = 2, nlonm1
      zmlzm = .5*(mlz(jlon,jlat,jklev)+mlz(jlon,jlat,jklev+1))
      ztkem = .5*zce*(tke(jlon,jlat,jklev)+tke(jlon,jlat,jklev+1))
      zdtfric = ztkem*sqrt(ztkem)/(zmlzm+.01)
      tetav(jlon,jlat,jklev) = tetav(jlon,jlat,jklev) + zdtfric*dtstep/(cpd*pai(jlon,jlat,jklev))
      enddo
      enddo

!  Mixing length

      do jklev = 2, ntop
      do jlon = 2, nlonm1
      zblack = zak*zhlev(jlon,jklev)*zmlmax/(zak*zhlev(jlon,jklev)+zmlmax)
      zrdz  = fmzh(jlon,jlat,jklev)/dz
      zbuoy = g*dlogthe(jlon,jklev)*zrdz
      if (zbuoy > 0.) then
      mlz(jlon,jlat,jklev) = min (zblack, .52*sqrt(zce*tke(jlon,jlat,jklev)/(zbuoy+1.e-6))) ! Deardorff (1980)

! The following alternative

!      mlz(jlon,jlat,jklev) = zblack/(1. + 12.*rich(jlon,jlat,jklev))
      else
      mlz(jlon,jlat,jklev) = zumltun*max (1./zrdz, zblack) ! unstable - layer thickness
      mlz(jlon,jlat,jklev) = min (mlz(jlon,jlat,jklev), zmlcut)
      endif
      enddo
      enddo

!  Vertical diffusion of TKE and W on half levels (tke(ntop+1)=tkemin)
!  Flux limited explicit scheme

      zlim = .05*dx**2/dtstep

      do jklev = 1, ntop
      do jlon = 2, nlonm1
      zkkk = sqrt(.5*zce*(tke(jlon,jlat,jklev+1)+tke(jlon,jlat,jklev)))
      chm(jlon,jlat,jklev) = min (zmlcut*zkkk, zlim) ! coefficient of hor. flux of momentum (T points)
      zkkk = .5*(mlz(jlon,jlat,jklev+1)+mlz(jlon,jlat,jklev))*zkkk
      zch(jlon,jklev)   = zzc0*fmz(jlon,jlat,jklev)**2*zkkk  ! approx.
      zah(jlon,jklev+1) = zch(jlon,jklev)                    ! approx.
      zkkk = min (zkkk, .2*dz**2/dtstep)     !! flux limited
      zdzh = zhlev(jlon,jklev+1)-zhlev(jlon,jklev)
!      zfe(jlon,jklev) = -zkkk*(tke(jlon,jlat,jklev+1)-tke(jlon,jlat,jklev))/zdzh    ! vertical flux of TKE (for explic. diff.)
!      zfw(jlon,jklev) = -2.*zkkk*(w(jlon,jlat,jklev+1)-w(jlon,jlat,jklev))/zdzh + & ! vertical flux of w (for explic. diff.)
!                        (tke(jlon,jlat,jklev+1)+tke(jlon,jlat,jklev))/3.
!      zfw(jlon,jklev) = (tke(jlon,jlat,jklev+1)+tke(jlon,jlat,jklev))/3.    ! partial vert. flux of w
      enddo
      enddo

! Explicit w and TKE vertical diffusion

!      do jk = 2, ntop
!      do jlon = 2, nlonm1
!      w(jlon,jlat,jk) = w(jlon,jlat,jk) - (zfw(jlon,jk)-zfw(jlon,jk-1))*dtstep*fmzh(jlon,jlat,jk)/dz
!      enddo
!      enddo

! Implicit w and TKE vertical diffusion (more stable, especially for TKE)

      do jlon = 2, nlonm1
      zrb = 1./(1. -zah(jlon,2) -zch(jlon,2))
      zet(jlon,2)      = -zch(jlon,2)*zrb
!      w  (jlon,jlat,2) = (w  (jlon,jlat,2)-zah(jlon,2)*w  (jlon,jlat,1))*zrb
      tke(jlon,jlat,2) = (tke(jlon,jlat,2)-zah(jlon,2)*tke(jlon,jlat,1))*zrb
      enddo
      do jk = 3, ntop
      do jlon = 2, nlonm1
      zb    = 1. - zah(jlon,jk) - zch(jlon,jk)
      zrden = 1./(zb +zah(jlon,jk)*zet(jlon,jk-1))
      zet(jlon,jk)      = -zch(jlon,jk)*zrden
!      w  (jlon,jlat,jk) = (w  (jlon,jlat,jk)-zah(jlon,jk)*w  (jlon,jlat,jk-1))*zrden
      tke(jlon,jlat,jk) = (tke(jlon,jlat,jk)-zah(jlon,jk)*tke(jlon,jlat,jk-1))*zrden
      enddo
      enddo
      do jk = ntop-1, 2, -1
      do jlon = 2, nlonm1
!      w  (jlon,jlat,jk) = zet(jlon,jk)*w  (jlon,jlat,jk+1) + w  (jlon,jlat,jk)
      tke(jlon,jlat,jk) = zet(jlon,jk)*tke(jlon,jlat,jk+1) + tke(jlon,jlat,jk)
      if (tke(jlon,jlat,jk).lt.tkemin) tke(jlon,jlat,jk) = tkemin
      enddo
      enddo

!-----------------------------------------------
!  Vertical diffusion of u, v, teta, q, qcw, qci
!  coefficients c(k), k=1,ntop  -  c(ntop)=0.
!  coefficients a(k), k=2,ntop
!-----------------------------------------------

      do jlon = 2, nlonm1
      zzc = -dtstep*fmz(jlon,jlat,1)*rd*tvirt(jlon,jlat,1)/(p(jlon,jlat,1)*dz)
      za1m(jlon) = zzc*roscdm(jlon,jlat)
      za1t(jlon) = zzc*roscdt(jlon,jlat)
      za1q(jlon) = zzc*roscdt(jlon,jlat)
      enddo

      do jklev = 1, ntop-1
      do jlon = 2, nlonm1
      zkkm  = mlz(jlon,jlat,jklev+1)*sqrt(zce*tke(jlon,jlat,jklev+1))
      cvm(jlon,jlat,jklev+1) = zkkm
      zroka = p(jlon,jlat,jklev)*tvirt(jlon,jlat,jklev+1)/(p(jlon,jlat,jklev+1)*tvirt(jlon,jlat,jklev))
      zrokc = .5*(1.+1./zroka)
      zcm(jlon,jklev)   = zzc0*zrokc*fmzh(jlon,jlat,jklev+1)*fmz(jlon,jlat,jklev)*zkkm
      zam(jlon,jklev+1) = zcm(jlon,jklev)*zroka*fmz(jlon,jlat,jklev+1)/fmz(jlon,jlat,jklev)
      zch(jlon,jklev)   = zcm(jlon,jklev  )/prandtl(jlon,jlat,jklev+1)
      zah(jlon,jklev+1) = zam(jlon,jklev+1)/prandtl(jlon,jlat,jklev+1)
      enddo
      enddo

      zcm(:,ntop) = 0.
      zch(:,ntop) = 0.

      if (jlat.eq.2) then   ! non proprio parallelo...
      za1mm1(:)  = za1m(:)
      zamm1(:,:) = zam(:,:)
      zcmm1(:,:) = zcm(:,:)
      endif

!-------------------------------------------------
!  tridiagonal matrix inversion:                 -
!  a(k)*psi(k-1)+b(k)*psi(k)+c(k)*psi(k+1)=r(k)  -
!  b(k) = 1-a(k)-c(k)                            -
!  r(k) = psim(k),   k>1                         -
!  r(1) = -a(1)*psisurf + psim(1)                -
!  where psi(k) = u,v,teta,q  at new time level  -
!  psim(k) = u,v,teta,q  at old time level       -
!  and psisurf their values at the surface.      -
!  resol. formula: psi(k)=e(k)*psi(k+1)+f(k)     -
!  a(1)=0 for zero surface flux
!-------------------------------------------------

      do jlon = 2, nlonm1
      jlonp1 = min (jlon+1, nlonm1)
      zrbu = 1./(1. -.5*(za1m(jlon)+za1m(jlonp1)) -.5*(zcm(jlon,1)+zcm(jlonp1,1)))
      zrbv = 1./(1. -.5*(za1m(jlon)+za1mm1(jlon)) -.5*(zcm(jlon,1)+zcmm1(jlon,1)))
      zrbt = 1./(1. -za1t(jlon) -zch(jlon,1))
      zrbq = 1./(1. -za1q(jlon) -zch(jlon,1))
      zeu(jlon,1) = -.5*(zcm(jlon,1)+zcm(jlonp1,1))*zrbu
      zev(jlon,1) = -.5*(zcm(jlon,1)+zcmm1(jlon,1))*zrbv
      zet(jlon,1) = -zch(jlon,1)*zrbt
      zeq(jlon,1) = -zch(jlon,1)*zrbq
      u(jlon,jlat,1) = u(jlon,jlat,1)*zrbu
      v(jlon,jlat,1) = v(jlon,jlat,1)*zrbv
      tetav(jlon,jlat,1) = (tetav(jlon,jlat,1)-za1t(jlon)*tetavs(jlon,jlat))*zrbt
      q(jlon,jlat,1) = (q(jlon,jlat,1)-za1q(jlon)*qskin(jlon,jlat))*zrbq
      qcw(jlon,jlat,1) = qcw(jlon,jlat,1)*zrbq
      qci(jlon,jlat,1) = qci(jlon,jlat,1)*zrbq
      enddo

      do jk = 2, ntop    ! psi(ntop)=f(ntop)
      do jlon = 2, nlonm1
      jlonp1 = min (jlon+1, nlonm1)
      zbu = 1.-.5*(zam(jlon,jk)+zam(jlonp1,jk)) -.5*(zcm(jlon,jk)+zcm(jlonp1,jk))
      zbv = 1.-.5*(zam(jlon,jk)+zamm1(jlon,jk)) -.5*(zcm(jlon,jk)+zcmm1(jlon,jk))
      zbh = 1.-zah(jlon,jk)-zch(jlon,jk)
      zrdenu = 1./(zbu +.5*(zam(jlon,jk)+zam(jlonp1,jk))*zeu(jlon,jk-1))
      zrdenv = 1./(zbv +.5*(zam(jlon,jk)+zamm1(jlon,jk))*zev(jlon,jk-1))
      zrdent = 1./(zbh +zah(jlon,jk)*zet(jlon,jk-1))
      zrdenq = 1./(zbh +zah(jlon,jk)*zeq(jlon,jk-1))
      zeu(jlon,jk) = -.5*(zcm(jlon,jk)+zcm(jlonp1,jk))*zrdenu
      zev(jlon,jk) = -.5*(zcm(jlon,jk)+zcmm1(jlon,jk))*zrdenv
      zet(jlon,jk) = -zch(jlon,jk)*zrdent
      zeq(jlon,jk) = -zch(jlon,jk)*zrdenq
      u    (jlon,jlat,jk) = (u    (jlon,jlat,jk) -.5*(zam(jlon,jk)+zam(jlonp1,jk))*u    (jlon,jlat,jk-1))*zrdenu
      v    (jlon,jlat,jk) = (v    (jlon,jlat,jk) -.5*(zam(jlon,jk)+zamm1(jlon,jk))*v    (jlon,jlat,jk-1))*zrdenv
      tetav(jlon,jlat,jk) = (tetav(jlon,jlat,jk) -    zah(jlon,jk)                *tetav(jlon,jlat,jk-1))*zrdent
      q  (jlon,jlat,jk) = (q  (jlon,jlat,jk) -zah(jlon,jk)*q  (jlon,jlat,jk-1))*zrdenq
      qcw(jlon,jlat,jk) = (qcw(jlon,jlat,jk) -zah(jlon,jk)*qcw(jlon,jlat,jk-1))*zrdenq
      qci(jlon,jlat,jk) = (qci(jlon,jlat,jk) -zah(jlon,jk)*qci(jlon,jlat,jk-1))*zrdenq
      enddo
      enddo

      do jk = ntop-1, 1, -1
      do jlon = 2, nlonm1
      u    (jlon,jlat,jk) = zeu(jlon,jk)*u    (jlon,jlat,jk+1) + u    (jlon,jlat,jk)
      v    (jlon,jlat,jk) = zev(jlon,jk)*v    (jlon,jlat,jk+1) + v    (jlon,jlat,jk)
      tetav(jlon,jlat,jk) = zet(jlon,jk)*tetav(jlon,jlat,jk+1) + tetav(jlon,jlat,jk)
      q    (jlon,jlat,jk) = zeq(jlon,jk)*q    (jlon,jlat,jk+1) + q    (jlon,jlat,jk)
      qcw  (jlon,jlat,jk) = zeq(jlon,jk)*qcw  (jlon,jlat,jk+1) + qcw  (jlon,jlat,jk)
      qci  (jlon,jlat,jk) = zeq(jlon,jk)*qci  (jlon,jlat,jk+1) + qci  (jlon,jlat,jk)
      enddo
      enddo

      za1mm1(:)  = za1m(:)
      zamm1(:,:) = zam(:,:)
      zcmm1(:,:) = zcm(:,:)

 1000 continue !  end of loop in latitude

      return
      end subroutine vdiff
!###############################################################################################################
      subroutine cloudfr

!  Defines fcloud: fraction of cloudy sky, as a function of
!  explicit cloud water, cloud ice, snow and relative humidity

      use mod_moloch, only : nlon, nlat, rich, p, t, q, qcw, qci, fcloud, qpi1, qccrit, ntop
      implicit none

      real zesat, zqs, zqs1, zqs2, zrh, zrh1, huc, zepneb, zrich
      integer jlat, jlon, jklev
      parameter (zepneb = 1.e-10, huc = 0.89)

      do jklev = 1, ntop
      do jlat = 2, nlat-1
      do jlon = 2, nlon-1
      call comp_esk (zesat, zqs1, t(jlon,jlat,jklev), p(jlon,jlat,jklev), 2)  ! blended saturation
      call comp_esk (zesat, zqs2, t(jlon,jlat,jklev), p(jlon,jlat,jklev), 3)  ! sat. to water below 0
      zqs = 0.70*zqs1 + 0.30*zqs2                    ! ad hoc to limit contrib. to high clouds
      zrh1 = min(1., q(jlon,jlat,jklev)/zqs)
      zrh = (zrh1-huc)/(1.-huc)
      zrh = min(1., zrh)
      zrh = max(0., zrh)

      fcloud(jlon,jlat,jklev) = (max(0.15*zrh, (qcw(jlon,jlat,jklev) + 0.77*qci(jlon,jlat,jklev) + &
                                0.2*qpi1(jlon,jlat,jklev))/qccrit))**0.7

      zrich = 0.5*(rich(jlon,jlat,jklev+1) + rich(jlon,jlat,jklev))
      if(zrich.lt..25.and.zrich.ge.0) fcloud(jlon,jlat,jklev) = (.875 + 0.5*zrich)*fcloud(jlon,jlat,jklev)
      if(zrich.gt.0..and.rich(jlon,jlat,jklev).lt.0.) fcloud(jlon,jlat,jklev) = 0.82*fcloud(jlon,jlat,jklev)
      if(zrich.lt.0.) fcloud(jlon,jlat,jklev) = 0.72*fcloud(jlon,jlat,jklev)
      fcloud(jlon,jlat,jklev) = max(zepneb, min(.999995, fcloud(jlon,jlat,jklev)))
      enddo
      enddo
      enddo

      return
      end subroutine cloudfr
!###############################################################################################################
      subroutine ccloudt

! Definition of total cloud cover cloudt
! Computes total cloud cover in a column with an algorithm modified from Geleyn's
! It is assumed that clouds are vertically correlated
! Calculation must proceed from top to bottom

      use mod_moloch, only : nlon, nlat, nlev, fcloud, cloudt, ntop
      implicit none

      real zprod, zcol(nlev), pnum, pden, sig, sigm
      integer jlat, jlon, jklev, k, kma, kmi, kmax(nlev), kmin(nlev)

      fcloud(:,:,1) = .999*fcloud(:,:,1)  ! to avoid "cloud holes" when cloud cover is constant at bottom

      do jlat = 2, nlat-1
      do jlon = 2, nlon-1
      zcol(1:ntop) = fcloud(jlon,jlat,ntop:1:-1)

      kmax = 0
      kmin = 0
      kma = 0
      kmi = 0
      sigm = 1.

      do k = 2, ntop
      if (zcol(k).eq.zcol(k-1)) cycle    ! avoids equal values
      sig = sign(1., zcol(k)-zcol(k-1))  ! 1 if zcol(k)>=zcol(k-1), else -1
      if (sig*sigm.eq.-1.) then          ! if opposite signs, it is an extreme
      sigm = sig
        if (sig.eq.-1.) then             ! if the second is < 0 ...
        kma = kma+1
        kmax(kma) = k-1                  ! ... then k-1 was a max. ...
        elseif (sig.eq.1.) then
        kmi = kmi+1
        kmin(kmi) = k-1                  ! ... else it was a min.
        endif
      endif
      enddo

      if (zcol(ntop).gt.zcol(ntop-1)) then
      kma = kma+1
      kmax(kma) = ntop        ! also the bottom level can be a max.
      endif

! product of probabilities of maxima at numerator

      pnum = 1.
      do k = 1, kma
      pnum = pnum*(1.-zcol(kmax(k)))
      enddo

! product of probabilities of minima at denominator

      pden = 1.
      do k = 1, kmi
      pden = pden*(1.-zcol(kmin(k)))
      enddo

! ratio of the two

      if (pden.ne.0.) zprod = pnum/pden

      cloudt(jlon,jlat) = max (0., 1.-zprod)
      cloudt(jlon,jlat) = min (cloudt(jlon,jlat), 1.)
      enddo
      enddo

      return
      end subroutine ccloudt
!###############################################################################################################
#ifndef rad_ecmwf
subroutine radiat_init(rds, ret, ozon, aerosol, aerotot, &
   nyrc, nmonc, ndayc, nhouc, nminc, &
   nlon, nlat, nlev, ntaer, dlat, dlon, &
   alat, alon, p, ps, t, tsurf, htop, proc_id)

! Initialization of astronomical parameters and
! total aerosols content for Ritter-Geleyn radiation scheme;
! ozone content and content of single type aerosol are defined also.
!
! Version for Moloch
!

use module_radiat_param
implicit none

! Input and output variables

integer, intent(in) :: nyrc, nmonc, ndayc, nhouc, nminc, nlon, nlat, nlev, ntaer, proc_id
real, intent(in) :: dlat, dlon
real, dimension(nlon,nlat), intent(in) :: alat, alon, ps, tsurf, htop
real, dimension(nlon,nlat,nlev), intent(in) :: p, t

real, intent(out) :: rds, ret
real, dimension(nlon,nlat,nlev), intent(out) :: ozon, aerotot
real, dimension(nlon,nlat,nlev,ntaer), intent(out) :: aerosol

! Internal variables

integer :: nsss, nzzaa, nzzmm, iiyr, nindat, iminut, jlon, jlat, jklev, jf, jk
real :: julian_day, rtime, rteta, rel, rem, rlls, rllls, zangozc

real, dimension(nlon,kaer,nlev) :: paer5
real, dimension(nlon,nlev) :: pap5, pt5, pozon5
real, dimension(nlon,nlev+1) :: paph5, pth5
real, dimension(nlon) :: plon5, plat5, pgelam5, pgemu5, pclon5, pslon5
real, dimension(nlev) :: zeta
real, dimension(nlev+1) :: zetah

! ---------------------------------------
  if (proc_id == 0) print *, 'Definition of astronomical variables, aerosol and ozone for radiation scheme'
! ---------------------------------------

! ---------------------------------------
! Astronomical variables:
! RDS is the declination of the Earth
! RET is the equation of time
! -------------------------------------

! nindat:  YYYYMMDD date of the simulation, for example 19870701

  iiyr = nyrc
  nindat =  iiyr*10000 + nmonc*100 + ndayc

! nsss: no. of seconds in a day (rounded to one minute)

  nsss = nhouc*3600 + nminc*60 

  iminut=int(float(nsss)/60.)

  nzzaa = nyrc - ( (1-SIGN(1,nmonc-3))/2 )
  nzzmm = nmonc + 6*(1-SIGN(1,nmonc-3))

  julian_day=1720994.5 + REAL(2-nzzaa/100 + (nzzaa/100)/4 + &
                         INT(365.25*REAL(nzzaa)) + &
                         INT(30.601*REAL(nzzmm+1)) + ndayc)

! Time in sec. after 12:00 gmt of 1 jan 2000

  rtime=(julian_day-2451545.)*rday + REAL(nsss)

  rteta=rtime/(rday*365.25)

!!!PTETA=rteta(rtime)

! Orbital parameterts of the Earth

  rel   = 1.7535 + 6.283076*rteta
  rem   = 6.240075 + 6.283020*rteta

! Relative movement Sun/Earth

  rlls  = 4.8951 + 6.283076*rteta
  rllls = 4.8952 + 6.283320*rteta - 0.0075*SIN(rel) - &
                  0.0326*COS(rel) - 0.0003*SIN(2.0*rel) + &
                  0.0002*COS(2.0*rel)

! Declination of the Earth

  rds = ASIN(SIN(repsm)*SIN(rllls))

! Equation of time

  ret = 591.8*SIN(2.0*rlls) - 459.4*SIN(rem) + &
             39.5*SIN(rem)*COS(2.0*rlls) - &
             12.7*SIN(4.*rlls) - 4.8*SIN(2.0*rem)

! ---------------------------------------

! Ozon and Aerosol definifion

  zch4 = 1.72e-06*zch4mwg/zairmwg
  zn2o = 310.e-09*zn2omwg/zairmwg
  zno2 = 500.e-13*zno2mwg/zairmwg
  zo3  =   1.e-06*zo3mwg /zairmwg
  zc11 = 280.e-12*zc11mwg/zairmwg
  zc12 = 484.e-12*zc12mwg/zairmwg
  zc22 =   1.e-12*zc22mwg/zairmwg
  zcl4 =   1.e-12*zcl4mwg/zairmwg

  call surdi
  call suaerh

!- Fortuin-Langematz O3 climatology

  call suecozc ( nindat , iminut )

!- ECMWF Geleyn O3 climatology

  zangozc=rel-1.7535
  call suecozo ( zangozc )

! climatological aerosol
! TEGEN/GISS AEROSOL CLIMATOLOGY

  if (naer == 1) then
    call suecaebc                  ! black carbon aerosol (fires)
    call suecaeor                  ! organic type aerosol
    call suecaesd                  ! soil-dust aerosol
    call suecaess                  ! sea-salt aerosol
    call suecaesu                  ! sulfate-type aerosol
    call suecaec (nindat, iminut)  ! climatol. distrib. of volcanic aerosol
  endif

! Preparation of atmospheric variable for aerosol and ozon definition (subroutine radaca)

! zeta and zetah: sigma-like coordinate for the routines that set parameters for aerosol

  do jk = 1, nlev
    zeta(jk)  = (jk-0.5)/float(nlev)
    zetah(jk) = (jk-1.0)/float(nlev)
  enddo
  zetah(nlev+1) = 1.

  loop_latitude_index: do jlat =1,nlat

! Geographical coordinate parameters

    do jlon = 1, nlon
      plat5(jlon)  = alat(jlon,jlat) * zdegrad
      plon5(jlon)  = alon(jlon,jlat) * zdegrad  ! long. < 0 not accepted
      if (plon5(jlon) < 0.) plon5(jlon) = plon5(jlon) + 2.*rpi
      pgelam5(jlon)= plon5(jlon)
      pgemu5(jlon) = sin(plat5(jlon))
      pclon5(jlon) = cos(plon5(jlon))
      pslon5(jlon) = sin(plon5(jlon))
    enddo

! Pressure at full levels (filtered in the vertical)

    do jlon = 1, nlon
      pap5(jlon,nlev) = p(jlon,jlat,1   )
      pap5(jlon,1   ) = p(jlon,jlat,nlev)
    enddo
    do jklev = 2, nlev-1
      jk = nlev+1-jklev
      do jlon = 1, nlon
          pap5(jlon,jk) = .25*(p(jlon,jlat,jklev-1)+p(jlon,jlat,jklev+1))+.5*p(jlon,jlat,jklev)
      enddo
    enddo

! Pressure at half levels

    do jlon = 1, nlon
      paph5(jlon,1)      = 0.
      paph5(jlon,nlev+1) = ps(jlon,jlat)
    enddo
    do jk = 2, nlev
      do jlon = 1, nlon
        paph5(jlon,jk) = 0.5*(pap5(jlon,jk-1)+pap5(jlon,jk))
      enddo
    enddo

! Temperature at full levels

    do jklev = 1, nlev
      jk = nlev+1-jklev
      do jlon = 1, nlon
        pt5(jlon,jk) = t(jlon,jlat,jklev)
      enddo
    enddo

! Temperature at half-levels: at top defined as lev. 1
! at bottom as tskin, elsewhere is the arithmetic mean

    do jlon = 1, nlon
      pth5(jlon,1)      = pt5(jlon,1)
      pth5(jlon,nlev+1) = tsurf(jlon,jlat)
    enddo

    do jk = 2, nlev
      do jlon = 1, nlon
        pth5(jlon,jk) = 0.5*(pt5(jlon,jk-1)+pt5(jlon,jk))
      enddo
    enddo

! Call to other set-up routines for various coefficients. These are dependent on the vertical resolution

    call suaerv (nlev, zetah, cvdaes, cvdael, cvdaeu, cvdaed, rctrbga, rcvobga, rcstbga, rcaeops, rcaeopl, &
                 rcaeopu, rcaeopd, rctrpt, rcaeadk, rcaeadm, rcaeros)

! Aerosol and ozon content definition

    call radaca (1, nlon, nlon, 1, nlev, &
                 paph5, pgelam5, pgemu5, pclon5, pslon5, pth5, paer5, pozon5)

! The computation of ozone must be done separately at each point in longitude, because radocz
! assumes (incorrectly for rotated grid) that latitude is the same for all vectors in longitude

    do jlon= 1, nlon
      call radozc (1, 1, 1, 1, nlev, 1, 1, 0, paph5(jlon,:), pgemu5(jlon), pozon5(jlon,:))
    enddo

    if (naer == 0) paer5(:,:,:) = zepaer

    do jklev = 1,nlev
      do jlon = 1,nlon
        ozon(jlon,jlat,jklev) = pozon5(jlon,  jklev)     ! ozone
        aerosol(jlon,jlat,jklev,1) = paer5(jlon,1,jklev) ! land (organic + sulfate) aerosol
        aerosol(jlon,jlat,jklev,2) = paer5(jlon,2,jklev) ! sea salt aerosol
        aerosol(jlon,jlat,jklev,3) = paer5(jlon,3,jklev) ! desert dust aerosol
        aerosol(jlon,jlat,jklev,4) = paer5(jlon,4,jklev) ! urban + black carbon aerosol
        aerosol(jlon,jlat,jklev,5) = paer5(jlon,5,jklev) ! volcanic aerosol
        aerosol(jlon,jlat,jklev,6) = paer5(jlon,6,jklev) ! stratospheric background aerosol
      enddo
    enddo

! Ad hoc increase of some aerosol (urban and sulfate) over the Po Valley in the lower troposphere

    do jklev = 1,nlev
      do jlon = 1,nlon
        if (alat(jlon,jlat) > 44.3.and.alat(jlon,jlat) < 46.2.and.   &
 alon(jlon,jlat) > 7.0.and.alon(jlon,jlat) < 13.4.and.   &
 htop(jlon,jlat) < 500..and.pap5(jlon,jklev) > 800.e2) then
      aerosol(jlon,jlat,jklev,1) = 1.10*aerosol(jlon,jlat,jklev,1)
      aerosol(jlon,jlat,jklev,4) = 1.20*aerosol(jlon,jlat,jklev,4)
        endif
      enddo
    enddo

  enddo loop_latitude_index

  do jf = 1, 4
    call filt3d (ozon, 1.)
  enddo
  ozon = max (ozon, 0.)
  do jf = 1, 20
    do jk = 1, 6
      call filt3d (aerosol(1:nlon,1:nlat,1:nlev,jk), 1.)
    enddo
  enddo
  do jf = 1, 90
    call filt3d (aerosol(1:nlon,1:nlat,1:nlev, 2), 1.)
  enddo
  aerosol = max (aerosol, zepaer)


  aerotot(:,:,:) = aerosol(:,:,:,1) + aerosol(:,:,:,2) + aerosol(:,:,:,3) +  &
                   aerosol(:,:,:,4) + aerosol(:,:,:,5) + aerosol(:,:,:,6)

return
end subroutine radiat_init
#endif
!###############################################################################################################
      subroutine radint (ndayr, nhouc, nminc, jstep)

!  Interface to radial

      use mod_moloch, only : nlon, nlat, nlev, nlonm1, nlatm1, nlevp1, hx, hy, h, dz, ps, p, t, q, qcw,     &
                             qci, qpi1, emisg1, fsnow, albedo, alsn, fmask, tskin, qskin, rd, g, cpd, cpv,  &
                             pi, geldt, gelvis, gelirr, alont, alatt, myid, aerotot, fcloud, slopeff, fice, &
                             cloudt, reqtim, rdecli

      implicit none

#ifdef mpi
        include 'mpif.h'
#endif

      common /radiazio/ stefan, ro, xts, gg, sdec, cdec, zalpn(nlon),                        &
              qfs(nlon,nlev), tfs(nlon,nlev), zdp(nlon,nlev), zpm(nlon,nlev), zp(nlon,nlev), &
              zqli(nlon,nlev), zqice(nlon,nlev), zneb(nlon,nlev), zcp(nlon,nlevp1),          &
              daer(nlon,nlev), tsf(nlon), alb(nlon), emis(nlon), zmu0(nlon),                 &
              dtfr(nlon,nlev), rat(nlon), rg(nlon), rirs(nlon), rvis(nlon)
      real(kind=8) stefan, ro, xts, gg, sdec, cdec, zalpn, qfs, tfs, zdp, zpm, zp, &
                   zqli, zqice, zneb, zcp, daer, tsf, alb, emis, zmu0,             &
                   dtfr, rat, rg, rirs, rvis
      real    ::  zz1, zgmt, zqhl, zclat, zslat(nlon), zhsv(nlon)
      real    ::  amux, amuy, amuz, zhx, zhy, znorm
      real    ::  aer(nlev), zarg, zalsn, zm, zhard
      integer ::  jstep, jlon, jlat, jklev, jk, ndayr, nhouc, nminc, ierr, comm

#ifdef mpi
      comm = mpi_comm_world
#endif

!  Definition of radiation parameters

      ro = 1365.       ! ECMWF rad. value (Geleyn orig.: ro = 1373.)
      ro = ro*0.982    ! further tuning vs. ECMWF vis. radiation at surface
      stefan = 5.67e-8
      gg     = g

!  Seasonal variation of solar constant (factor - Pielke, p. 211)

      zarg = 2.*pi*float(ndayr-1)/365.
      xts  = 1.00011 + .034221*cos(zarg) + .00128*sin(zarg) + .000719*cos(2.*zarg) + .000077*sin(2.*zarg)

!  Solar declination (rdecli defined in aerdef (first step) and radintec)

      sdec = sin(rdecli)
      cdec = cos(rdecli)

!  Vertical profile of aerosol content (aer is a vertical integral; daer is aerosol on integer levels)
!  (case in which ECMWF aerosol is not used)

!      do jklev = 1, nlev
!      jk = nlev + 1 - jklev
!      zz1 = 1. - (jklev-1)*dz/h
!      aer(jk) = zz1*(.28 - .43*zz1**2 + .44*zz1**4)
!      enddo

!  Loop over all latitudes

      do 1000 jlat = 2, nlatm1

!  (Case in which ECMWF aerosol is not used)

!      do jklev = 1, nlev-1
!      do jlon = 2, nlonm1
!      daer(jlon,jklev) = aer(jklev+1) - aer(jklev)
!      enddo
!      enddo
!      daer(:,nlev) = daer(:,nlev-1)

!  Case in which the sum of ECMWF aerosol is used

      do jklev = 1, nlev
      do jlon = 2, nlonm1
      daer(jlon,jklev) = 1.08*aerotot(jlon,jlat,jklev) ! tuning vs. ECMWF vis. radiation at surf.
      enddo
      enddo

!  Temp., spec. humidity, water and ice clouds
!  Fraction of cloud cover fcloud put in zneb

      do jklev = 1, nlev
      jk = nlevp1 - jklev
      do jlon = 2, nlonm1
      tfs(jlon,jk)  = t  (jlon,jlat,jklev)
      qfs(jlon,jk)  = q  (jlon,jlat,jklev)
      zqli(jlon,jk) = qcw(jlon,jlat,jklev)
      zqice(jlon,jk)= qci(jlon,jlat,jklev) + 0.2*qpi1(jlon,jlat,jklev)
      if(zqli (jlon,jk).lt.1.d-12) zqli (jlon,jk)=1.d-12
      if(zqice(jlon,jk).lt.1.d-12) zqice(jlon,jk)=1.d-12
      zneb(jlon,jk) = fcloud(jlon,jlat,jklev)
      enddo
      enddo

!  Specific heat of moist air at half-levels

      do jk = 2, nlev
      do jlon = 2, nlonm1
      zqhl = .5*(qfs(jlon,jk) + qfs(jlon,jk-1))
      zcp(jlon,jk) = cpd*(1.-zqhl) + cpv*zqhl
      enddo
      enddo
      do jlon = 2, nlonm1
      zcp(jlon,1)      = cpd
      zcp(jlon,nlevp1) = cpd*(1.-qskin(jlon,jlat)) + cpv*qskin(jlon,jlat)
      enddo

!  Smoothing of pressure profile

      do jlon = 2, nlonm1
      zpm(jlon,nlev) = p(jlon,jlat,1)
      zpm(jlon,1   ) = p(jlon,jlat,nlev)
      enddo
      do jklev = 2, nlev-1
      jk = nlevp1-jklev
      do jlon = 2, nlonm1
      zpm(jlon,jk) = .25*(p(jlon,jlat,jklev-1)+p(jlon,jlat,jklev+1))+.5*p(jlon,jlat,jklev)
      enddo
      enddo

!  Pressure on half levels. zp(nlev) is the surface pressure
!  Pressure at top is zp(1) minus the thickness of the last integer layer zdp(1)
!  Pressure increments between half-levels

      do jk = 1, nlev-1
      do jlon = 2, nlonm1
      zp(jlon,jk) = 0.5*(zpm(jlon,jk+1) + zpm(jlon,jk))
      enddo
      enddo
      do jlon = 2, nlonm1
      zp(jlon,nlev) = ps(jlon,jlat)
      zdp(jlon,1)   = zp(jlon,1)
      enddo
      do jk = 2, nlev
      jklev = nlev-jk+1
      do jlon = 2, nlonm1
      zdp(jlon,jk) = zp(jlon,jk)-zp(jlon,jk-1)
        if (zdp(jlon,jk).le.0.) then
        print*,' Negative pressure thickness in radial --- jlon / jlat / jklev ='
        print*,'                               ',jlon,jlat,jklev
          if (jklev.gt.nlev/2) then
          print*,' ######## Aborted ######## at jstep =',jstep
#ifdef mpi
            call mpi_abort(comm, ierr)
#endif
          stop
          else
          zdp(jlon,jk) = 10.
          endif
        endif
      enddo
      enddo

!  Surface temperature, albedo and emissivity
!  Change of snow albedo (alb) in a range depending on skin temperature (as in radintec)

      do jlon = 2, nlonm1
      tsf(jlon) = tskin(jlon,jlat)
!      if (tsf(jlon).lt.277.) then  ! snow rejuvenation
!      zalsn = alsn -.1 + .2*(277. - tsf(jlon))/14.
!      zalsn = min (zalsn, alsn + .1)
!      else
!      zalsn = alsn-.1
!      endif
!      alb(jlon)  = albedo(jlon,jlat)*(1.-fsnow(jlon,jlat)) + zalsn*fsnow(jlon,jlat)
      alb(jlon)  = albedo(jlon,jlat)*(1.-fsnow(jlon,jlat)) + alsn(jlon,jlat)*fsnow(jlon,jlat)
      emis(jlon) = emisg1(jlon,jlat)*(1.-fsnow(jlon,jlat)) + 0.99* fsnow(jlon,jlat)
      enddo

!  True solar hour (units: hours and fractions of hours)
!  Reqtim (in sec) accounts for astron. time; it is computed in aerdef and radintec;
!  it is required to assure syncronization with ECMWF solar radiation
!  Albedo: for parallel (direct) radiation, dependency on solar zenith angle is computed
!  over land (Yang et al) and over sea (Taylor et al, as ECMWF)

      zgmt = float(nhouc) + float(nminc)/60. + reqtim/3600.
      do jlon = 2, nlonm1
      zhsv(jlon)  = zgmt + alont(jlon,jlat)*24./360.
      zslat(jlon) = sin(alatt(jlon,jlat)*pi/180.)
      zmu0(jlon)  = zslat(jlon)*sdec-sqrt(1.-zslat(jlon)**2)*cdec*cos(.2618993878d0*zhsv(jlon))
      zm          = max(zmu0(jlon), 1.e-12)
      if (fmask(jlon,jlat).lt.0.5) then
      zalpn(jlon) = min(1., alb(jlon)*(1. + 0.26)/(1. + 2.*0.26*zm))
      else
      zhard = 1. - fmask(jlon,jlat) + fice(jlon,jlat)*fmask(jlon,jlat)
      zalpn(jlon) = (1.-zhard)*0.037/(1.1*zm**1.4 + 0.15) + alb(jlon)*zhard
      endif
      enddo

!  Call of radiation subroutine

      call radial

!  Output variables

      do jklev = 1, nlev
      jk = nlevp1 - jklev
      do jlon = 2, nlonm1
      geldt(jlon,jlat,jklev) = sngl(dtfr(jlon,jk))
      enddo
      enddo

!  Surface fluxes of visible and infrared radiation (positive downward)

      do jlon = 2, nlonm1

!  True solar declination (to define slopeff: effective slope for solar radiation)

      zclat = cos(alatt(jlon,jlat)*pi/180.)
      amuz  = sdec*zslat(jlon) + cdec*zclat*cos(.2618*(zhsv(jlon)-12.))
      if (amuz.gt.1.e-5) then
      amuy = sdec*zclat - cdec*zslat(jlon)*cos(.2618*(zhsv(jlon)-12.))
      amux = -cdec*sin(.2618*(zhsv(jlon)-12.))
      zhx  = .5*(hx(jlon,jlat) + hx(jlon-1,jlat))
      zhy  = .5*(hy(jlon,jlat) + hy(jlon,jlat+1))
      znorm= sqrt(1. + zhx**2+ zhy**2)
      slopeff(jlon,jlat) = (amuz - amux*zhx - amuy*zhy)/(amuz*znorm)
      slopeff(jlon,jlat) = min (slopeff(jlon,jlat), 5.) ! 05/04/2024
      else
      slopeff(jlon,jlat) = 1.
      endif

!      gelvis(jlon,jlat) = sngl(rg(jlon))
      gelvis(jlon,jlat) = (0.3 + 0.7*slopeff(jlon,jlat))*sngl(rg(jlon)) ! assuming constant ratio of direct to tot. sky sw rad.
!      gelirr(jlon,jlat) = max (sngl(rat(jlon)), sngl(rat(jlon))*0.855) ! tuning vs. ECMWF ir radiation at surface
      gelirr(jlon,jlat) = max (sngl(rat(jlon)), sngl(rat(jlon))*0.855*0.92*(1.-0.05*cloudt(jlon,jlat))) ! more tuning!

      enddo

1000  continue

      return
      end subroutine radint
!###############################################################################################################
      subroutine radial

! Ritter -Geleyn radiation scheme (modified)

      use mod_moloch, only: nlon, nlonm1, nlev, nlevp1, co2ppm
      implicit integer (i-n), real(kind=8) (a-h, o-z)
      parameter (kideb = 2, kiut = nlonm1, nk = nlev)

      common /radiazio/ stefan, ro, xts, gg, sdec, cdec, zalpn(nlon),                     &
              qfs(nlon,nlev), tfs(nlon,nlev), dp(nlon,nlev), pm(nlon,nlev), p(nlon,nlev), &
              zqli(nlon,nlev), zqice(nlon,nlev), zneb(nlon,nlev), cph(nlon,nlevp1),       &
              daer(nlon,nlev), tsf(nlon), alb(nlon), emis(nlon), zmu0(nlon),              &
              dtfr(nlon,nlev), rat(nlon), rg(nlon), rirs(nlon), rvis(nlon)

!  ZONE DE TRAVAIL DE YRADNEB

       common /zotraf/                                                   &
             zfpc(nlon,nlevp1),zfpn(nlon,nlevp1),zfdc(nlon,nlevp1),      &
             zfdn(nlon,nlevp1),zfmc(nlon,nlevp1),zfmn(nlon,nlevp1),      &
             ztu1(nlon,nlev),ztu2(nlon,nlev),ztu3(nlon,nlev),            &
             ztu4(nlon,nlev),ztu5(nlon,nlev),ztu6(nlon,nlev),            &
             ztu7(nlon,nlev),ztu8(nlon,nlev),ztu9(nlon,nlev),            &
             zbb(nlon,nlevp1),zzfdc(nlon,nlevp1),zzfdn(nlon,nlevp1),     &
             zzfmc(nlon,nlevp1),zzfmn(nlon,nlevp1),za1c(nlon),           &
             za2c(nlon),za3c(nlon),za4c(nlon),za5c(nlon),za1n(nlon),     &
             za2n(nlon),za3n(nlon),za4n(nlon),za5n(nlon),zal1(nlon),     &
             zal2(nlon),zbe1(nlon),zbe2(nlon),zga1(nlon),zga2(nlon),     &
             zde1(nlon),zde2(nlon),znebl(nlon),znebcl(nlon),             &
             zeoss(nlon),zeots(nlon),                                    &
             zmu0i(nlon),zusa(nlon),zusn(nlon),zusi(nlon)
      dimension zueogs(nlon,nlev),zdeogs(nlon,nlev),zdeogt(nlon,nlev),   &
             zb1(nlon,nlev),zb2(nlon,nlev),zb3(nlon,nlev),               &
             zb4(nlon,nlev),zueogt(nlon,nlev)
      equivalence (zueogs(nlon,nlev),zfpc(nlon,nlevp1)),                 &
                  (zdeogs(nlon,nlev),zfpn(nlon,nlevp1)),                 &
                  (zdeogt(nlon,nlev),zfdc(nlon,nlevp1)),                 &
                  (zueogt(nlon,nlev),zzfdc(nlon,nlevp1)),                &
                  (zb1(1,1),ztu1(1,1)),(zb2(1,1),ztu2(1,1)),             &
                  (zb3(1,1), ztu3(1,1)),(zb4(1,1),ztu4(1,1))
      dimension znmaxb(nlon),znminb(nlon),znmaxh(nlon),znminh(nlon),     &
                zstab(nlon),ztrans(nlon),ztranb(nlon),                   &
                zfnet(nlon),zztrans(nlon),zztranb(nlon),znsh(nlon),      &
                znsc(nlon),znso(nlon),znth(nlon),zntc(nlon),znto(nlon),  &
                zrsh(nlon),zrsc(nlon),zrso(nlon),zrth(nlon),zcth(nlon),  &
                zrtc(nlon),zrto(nlon),zeos(nlon),zeot(nlon),zeosh(nlon), &
                zeoth(nlon),zeosb(nlon),zeotb(nlon),znsor(nlon),         &
                zdmu0i(nlon),zueots(nlon)
      equivalence (znmaxb(1),ztrans(1),znsh(1),za1c(1)),(znminb(1),      &
                   ztranb(1),znsc(1),za2c(1)),(znmaxh(1),zfnet(1),       &
                   znso(1),za3c(1)),(znminh(1),znth(1),za4c(1)),         &
                  (zstab(1),zntc(1),za5c(1)),(zztrans(1),                &
                   znto(1),za1n(1)),(zztranb(1),zrsh(1),za2n(1)),        &
                  (zrsc(1),za3n(1)),(zrso(1),za4n(1)),(zrth(1),za5n(1)), &
                  (zcth(1),zal1(1)),(zrtc(1),zal2(1)),(zrto(1),zbe1(1)), &
                  (zeos(1),zbe2(1)),(zeot(1),zga1(1)),(zeosh(1),         &
                   zga2(1)),(zeoth(1),zde1(1)),(zeosb(1),zde2(1)),       &
                  (zeotb(1),znebl(1)),(znsor(1),znebcl(1)),              &
                  (zdmu0i(1),zmu0i(1)),(zueots(1),zusa(1))

!-----------------------------------------------------------------------
!     O - DEFINITION DES DEGRES DE LIBERTE (CONSTANTES OPTIQUES
!         ESSENTIELLEMENT).

      integer iideb(100), iiut(100)
      logical lsideb,lsi
      dimension zga(6),zgb(6),zgc(6),znp(5,6),zdp(5,6),zgas2b(6),zg4b(6)
      data                                                              &
       zga/.8041d-01,.1456d+00,.4787d+01,.2102d+04,.1334d+01,.1551d-04/ &
      ,zgb/.8968d+07,.2413d+10,.3548d+05,.6370d+10,.8499d+11,.1012d+06/ &
      ,zgc/.5925d-10,.1842d-10,.2532d-07,.1953d+07,.1734d-11,.1225d-16/ &
      ,zgd4/.3608d-69/,zge4/.7563d+04/
      data znp/.69926d+01,.63915d+00,.28896d-03,.10050d-08,.99037d-16,  &
               .31105d+01,.14225d-01,.69355d-06,.36087d-12,.44113d-20,  &
               .23659d+01,.11139d-02,.10618d-05,0.,0.,                  &
               .19810d+03,.46954d+04,.22512d+04,.52461d+02,.11645d+00,  &
               .46348d+02,.35630d+02,.33005d+01,.18045d-01,.88667d-05,  &
               .47413d+01,.16334d+01,.48164d+00,.56140d-02,.67790d-04/
      data zdp/.21868d+02,.17453d+02,.68918d+00,.23456d-03,.22317d-09,  &
               .48401d+02,.18648d+02,.35199d-01,.63691d-06,.86395d-13,  &
               .76948d+02,.41056d+01,.42667d-03,0.,0.,                  &
               .27241d+03,.57091d+04,.14393d+05,.29879d+04,.25382d+02,  &
               .86942d+02,.32186d+03,.10775d+03,.21261d+01,.40003d-02,  &
               .24408d+02,.81919d+01,.72193d+01,.56230d+00,.25384d-02/

!  SOME CONSTANT VALUES MODIF. AFTER GELEYN'S SUGGEST., SEPT. 92

      data zeoray/.8606d-06/
      data zeoata/.5037d-01/,zeodta/.1766d-01/,zbsfta/.3471d+00/
      data zeoasa/.6727d-01/,zeodsa/.3665d+00/,zbsfsa/.3490d+00/
      data zusaa/-.3020d+00/,zusba/0./
      data zeoatn/.9087d+01/,zeodtn/.4329d+01/,zbsftn/.3584d+00/
      data zeoasn/.3953d-01/,zeodsn/.6558d+01/,zbsfsn/.3300d+00/
      data zusan/-.3400d+00/,zusbn/0./
      data zeoati/.2994d+01/,zeodti/.1066d+01/,zbsfti/.3268d+00/
      data zeoasi/.1592d-01/,zeodsi/.1419d+01/,zbsfsi/.3246d+00/
      data zusai/-.3508d+00/,zusbi/0./

      do 1 jga=1,6
      zgas2b(jga)=zga(jga)/(2.*zgb(jga))
      zg4b(jga)=4.*zgb(jga)
    1 continue
      zeo2ta=2.*zbsfta*zeodta
      zeo1ta=zeo2ta+2.*zeoata
      zeo2tn=2.*zbsftn*zeodtn
      zeo1tn=zeo2tn+2.*zeoatn
      zeo2ti=2.*zbsfti*zeodti
      zeo1ti=zeo2ti+2.*zeoati
      zeo2sa=2.*zbsfsa*zeodsa
      zeo1sa=zeo2sa+2.*zeoasa
      zeo2sn=2.*zbsfsn*zeodsn
      zeo1sn=zeo2sn+2.*zeoasn
      zeo2si=2.*zbsfsi*zeodsi
      zeo1si=zeo2si+2.*zeoasi
      zeosa=zeodsa+zeoasa
      zeosn=zeodsn+zeoasn
      zeosi=zeodsi+zeoasi

      zrae=0.001324d+00
      zzrae=zrae*(zrae+2.)
      zirae=1./zrae
      zroxts=ro*xts
      zpeps=1.d-03
      zqeps=1.d-09
      zepres=1.d-12
      zargli=-250.d+00

!     O.1   CALCUL DU COSINUS DE L'ANGLE SOLAIRE ZENITHAL ET DES
!           LIMITES JOUR/NUIT

      do 10 ji=kideb,kiut
      zdmu0i(ji)=0.5*(sqrt(zmu0(ji)*zmu0(ji)+zzrae)-zmu0(ji))*zirae
      rg(ji)=0.
      rvis(ji)=0.
   10 continue
      iaucr=0
      lsideb=zmu0(kideb).gt.0.
      if(.not.lsideb) go to 11
      iaucr=1
      iideb(1)=kideb
   11 continue
      do 14 ji=kideb,kiut
      lsi=zmu0(ji).gt.0.
      if(.not.(lsi.neqv.lsideb)) go to 14
      if(.not.lsi) go to 12
      iaucr=iaucr+1
      iideb(iaucr)=ji
      go to 13
   12 continue
      iiut(iaucr)=ji-1
   13 continue
      lsideb=lsi
   14 continue
      if(.not.lsi) go to 15
      iiut(iaucr)=kiut
   15 continue

!     O.3   FLUX DE CORPS NOIR.

      do 31 jk=1,nk
      do 30 ji=kideb,kiut
      zbb(ji,jk)=stefan*(tfs(ji,jk)*tfs(ji,jk))*(tfs(ji,jk)*tfs(ji,jk))
   30 continue
   31 continue
      do 32 ji=kideb,kiut
      zbb(ji,nk+1)=stefan*(tsf(ji)*tsf(ji))*(tsf(ji)*tsf(ji))
   32 continue

!-----------------------------------------------------------------------
!     I - CALCUL DES EPAISSEURS OPTIQUES GASEUSES, DE LA NEBULOSITE ET
!         DES CONDITIONS AUX LIMITES AU SOL.

      zmd=0.5*exp(0.5d+00)
      zco2=dble(co2ppm*1.e-6*44./29.)
      zo31=0.6012d-01
      zo32=0.3166d+04
      zo32=sqrt(zo32*zo32*zo32)

!     I.1   SOMMET DE L'ATMOSPHERE DU MODELE (OU VALEUR ARBITRAIREMENT
!           FAIBLE).

      do 110 ji=kideb,kiut
      zp=max(zpeps,p(ji,1)-dp(ji,1))
      zq=max(zqeps,qfs(ji,1))
      zrt=sqrt(tfs(ji,1))
      zti=1./tfs(ji,1)
      zt2=tfs(ji,1)*tfs(ji,1)
      zt2i=zti*zti
      zt4i=zt2i*zt2i
      ztdeta=zrt**49*exp(zge4/tfs(ji,1))
      zt3=zt2*tfs(ji,1)
      znsh(ji)=(2.*zp)*zq
      znsc(ji)=(2.*zp)*zco2
      znso(ji)=(2.*zo31)/(1.+zo32/sqrt(zp*zp*zp))
      znth(ji)=znsh(ji)*zti
      zntc(ji)=znsc(ji)
      znto(ji)=znso(ji)*zt2
      zrsh(ji)=znsh(ji)*(zmd*0.5*zp)*zrt
      zrsc(ji)=znsc(ji)*(zmd*0.5*zp)*tfs(ji,1)
      zrso(ji)=znso(ji)*(zmd*0.5*zp)*zrt
      zrth(ji)=zrsh(ji)*zt2i*zrt
      zcth(ji)=zrth(ji)*zt4i*(1.+zgd4*zq*ztdeta)
      zrtc(ji)=zrsc(ji)*tfs(ji,1)
      zrto(ji)=zrso(ji)*zt3*zrt
      znsor(ji)=znso(ji)
      znsh(ji)=znsh(ji)*zdmu0i(ji)
      znsc(ji)=znsc(ji)*zdmu0i(ji)
      znso(ji)=znso(ji)*zdmu0i(ji)
      zrsh(ji)=zrsh(ji)*(zdmu0i(ji)/zmd)
      zrsc(ji)=zrsc(ji)*(zdmu0i(ji)/zmd)
      zrso(ji)=zrso(ji)*(zdmu0i(ji)/zmd)
  110 continue
      do 111 jn=1,iaucr
      do 111 ji=iideb(jn),iiut(jn)
      zwsh=(zgas2b(1)*(zrsh(ji)/znsh(ji)))*(sqrt(1.+zg4b(1)*znsh(ji) &
          /(zrsh(ji)/znsh(ji)))-1.)+zgc(1)*zrsh(ji)
      zwsc=(zgas2b(2)*(zrsc(ji)/znsc(ji)))*(sqrt(1.+zg4b(2)*znsc(ji) &
          /(zrsc(ji)/znsc(ji)))-1.)+zgc(2)*zrsc(ji)
      zwso=(zgas2b(3)*(zrso(ji)/znso(ji)))*(sqrt(1.+zg4b(3)*znso(ji) &
          /(zrso(ji)/znso(ji)))-1.)+zgc(3)*zrso(ji)
      zeos(ji)=(zwsh*(1.+zwsh*(znp(1,1)+zwsh*(znp(2,1)+zwsh*(znp(3,1)          &
              +zwsh*(znp(4,1)+zwsh*znp(5,1)))))))/(1.+zwsh*(zdp(1,1)           &
              +zwsh*(zdp(2,1)+zwsh*(zdp(3,1)+zwsh*(zdp(4,1)+zwsh*zdp(5,1)))))) &
              +(zwsc*(1.+zwsc*(znp(1,2)+zwsc*(znp(2,2)+zwsc*(znp(3,2)          &
              +zwsc*(znp(4,2)+zwsc*znp(5,2)))))))/(1.+zwsc*(zdp(1,2)           &
              +zwsc*(zdp(2,2)+zwsc*(zdp(3,2)+zwsc*(zdp(4,2)+zwsc*zdp(5,2)))))) &
              +(zwso*(1.+zwso*(znp(1,3)+zwso*(znp(2,3)+zwso*znp(3,3)))))       &
              /(1.+zwso*(zdp(1,3)+zwso*(zdp(2,3)+zwso*zdp(3,3))))
      zeoss(ji)=zeos(ji)/zdmu0i(ji)
  111 continue
      do 112 ji=kideb,kiut
      zwth=(zgas2b(4)*(zrth(ji)/znth(ji)))*(sqrt(1.+zg4b(4)*znth(ji)          &
          /(zrth(ji)/znth(ji)))-1.)+zgc(4)*zcth(ji)
      zwtc=(zgas2b(5)*(zrtc(ji)/zntc(ji)))*(sqrt(1.+zg4b(5)*zntc(ji)          &
          /(zrtc(ji)/zntc(ji)))-1.)+zgc(5)*zrtc(ji)
      zwto=(zgas2b(6)*(zrto(ji)/znto(ji)))*(sqrt(1.+zg4b(6)*znto(ji)          &
          /(zrto(ji)/znto(ji)))-1.)+zgc(6)*zrto(ji)
      zeot(ji)=(zwth*(1.+zwth*(znp(1,4)+zwth*(znp(2,4)+zwth*(znp(3,4)          &
              +zwth*(znp(4,4)+zwth*znp(5,4)))))))/(1.+zwth*(zdp(1,4)           &
              +zwth*(zdp(2,4)+zwth*(zdp(3,4)+zwth*(zdp(4,4)+zwth*zdp(5,4)))))) &
              +(zwtc*(1.+zwtc*(znp(1,5)+zwtc*(znp(2,5)+zwtc*(znp(3,5)          &
              +zwtc*(znp(4,5)+zwtc*znp(5,5)))))))/(1.+zwtc*(zdp(1,5)           &
              +zwtc*(zdp(2,5)+zwtc*(zdp(3,5)+zwtc*(zdp(4,5)+zwtc*zdp(5,5)))))) &
              +(zwto*(1.+zwto*(znp(1,6)+zwto*(znp(2,6)+zwto*(znp(3,6)          &
              +zwto*(znp(4,6)+zwto*znp(5,6)))))))/(1.+zwto*(zdp(1,6)           &
              +zwto*(zdp(2,6)+zwto*(zdp(3,6)+zwto*(zdp(4,6)+zwto*zdp(5,6))))))
      zeots(ji)=zeot(ji)
  112 continue

!     I.2   BOUCLE SUR LES NIVEAUX VERTICAUX.

      do 126 jk=1,nk
      do 120 jn=1,iaucr
      do 120 ji=iideb(jn),iiut(jn)
      zeosh(ji)=zeos(ji)
  120 continue
      do 121 ji=kideb,kiut
      zeoth(ji)=zeot(ji)
  121 continue
      do 122 ji=kideb,kiut
      zp=p(ji,jk)
      zq=max(zqeps,qfs(ji,jk))
      zrt=sqrt(tfs(ji,jk))
      zti=1./tfs(ji,jk)
      zt2=tfs(ji,jk)*tfs(ji,jk)
      zt2i=zti*zti
      zt4i=zt2i*zt2i
      ztdeta=zrt**49*exp(zge4/tfs(ji,jk))
      zt3=zt2*tfs(ji,jk)
      zduh=(2.*dp(ji,jk))*zq
      zduc=(2.*dp(ji,jk))*zco2
      zduo=(2.*zo31)/(1.+zo32/sqrt(zp*zp*zp))-znsor(ji)
      znsor(ji)=znsor(ji)+zduo
      znsh(ji)=znsh(ji)+zduh*zdmu0i(ji)
      znsc(ji)=znsc(ji)+zduc*zdmu0i(ji)
      znso(ji)=znso(ji)+zduo*zdmu0i(ji)
      znth(ji)=znth(ji)+zduh*zti
      zntc(ji)=zntc(ji)+zduc
      znto(ji)=znto(ji)+zduo*zt2
      zduh=zduh*(zmd*pm(ji,jk))*zrt
      zduc=zduc*(zmd*pm(ji,jk))*tfs(ji,jk)
      zduo=zduo*(zmd*pm(ji,jk))*zrt
      zrsh(ji)=zrsh(ji)+zduh*(zdmu0i(ji)/zmd)
      zrsc(ji)=zrsc(ji)+zduc*(zdmu0i(ji)/zmd)
      zrso(ji)=zrso(ji)+zduo*(zdmu0i(ji)/zmd)
      zrth(ji)=zrth(ji)+(zduh*zt2i*zrt)
      zcth(ji)=zcth(ji)+(zduh*zt2i*zrt)*zt4i*(1.+zgd4*zq*ztdeta)
      zrtc(ji)=zrtc(ji)+zduc*tfs(ji,jk)
      zrto(ji)=zrto(ji)+zduo*zt3*zrt
  122 continue
      do 123 jn=1,iaucr
      do 123 ji=iideb(jn),iiut(jn)
      zwsh=(zgas2b(1)*(zrsh(ji)/znsh(ji)))*(sqrt(1.+zg4b(1)*znsh(ji) &
          /(zrsh(ji)/znsh(ji)))-1.)+zgc(1)*zrsh(ji)
      zwsc=(zgas2b(2)*(zrsc(ji)/znsc(ji)))*(sqrt(1.+zg4b(2)*znsc(ji) &
          /(zrsc(ji)/znsc(ji)))-1.)+zgc(2)*zrsc(ji)
      zwso=(zgas2b(3)*(zrso(ji)/znso(ji)))*(sqrt(1.+zg4b(3)*znso(ji) &
          /(zrso(ji)/znso(ji)))-1.)+zgc(3)*zrso(ji)
      zeos(ji)=(zwsh*(1.+zwsh*(znp(1,1)+zwsh*(znp(2,1)+zwsh*(znp(3,1)          &
              +zwsh*(znp(4,1)+zwsh*znp(5,1)))))))/(1.+zwsh*(zdp(1,1)           &
              +zwsh*(zdp(2,1)+zwsh*(zdp(3,1)+zwsh*(zdp(4,1)+zwsh*zdp(5,1)))))) &
              +(zwsc*(1.+zwsc*(znp(1,2)+zwsc*(znp(2,2)+zwsc*(znp(3,2)          &
              +zwsc*(znp(4,2)+zwsc*znp(5,2)))))))/(1.+zwsc*(zdp(1,2)           &
              +zwsc*(zdp(2,2)+zwsc*(zdp(3,2)+zwsc*(zdp(4,2)+zwsc*zdp(5,2)))))) &
              +(zwso*(1.+zwso*(znp(1,3)+zwso*(znp(2,3)+zwso*znp(3,3)))))       &
              /(1.+zwso*(zdp(1,3)+zwso*(zdp(2,3)+zwso*zdp(3,3))))
      zdeogs(ji,jk)=(zeos(ji)-zeosh(ji))/zdmu0i(ji)
  123 continue
      do 124 ji=kideb,kiut
      zwth=(zgas2b(4)*(zrth(ji)/znth(ji)))*(sqrt(1.+zg4b(4)*znth(ji) &
          /(zrth(ji)/znth(ji)))-1.)+zgc(4)*zcth(ji)
      zwtc=(zgas2b(5)*(zrtc(ji)/zntc(ji)))*(sqrt(1.+zg4b(5)*zntc(ji) &
          /(zrtc(ji)/zntc(ji)))-1.)+zgc(5)*zrtc(ji)
      zwto=(zgas2b(6)*(zrto(ji)/znto(ji)))*(sqrt(1.+zg4b(6)*znto(ji) &
          /(zrto(ji)/znto(ji)))-1.)+zgc(6)*zrto(ji)
      zeot(ji)=(zwth*(1.+zwth*(znp(1,4)+zwth*(znp(2,4)+zwth*(znp(3,4)          &
              +zwth*(znp(4,4)+zwth*znp(5,4)))))))/(1.+zwth*(zdp(1,4)           &
              +zwth*(zdp(2,4)+zwth*(zdp(3,4)+zwth*(zdp(4,4)+zwth*zdp(5,4)))))) &
              +(zwtc*(1.+zwtc*(znp(1,5)+zwtc*(znp(2,5)+zwtc*(znp(3,5)          &
              +zwtc*(znp(4,5)+zwtc*znp(5,5)))))))/(1.+zwtc*(zdp(1,5)           &
              +zwtc*(zdp(2,5)+zwtc*(zdp(3,5)+zwtc*(zdp(4,5)+zwtc*zdp(5,5)))))) &
              +(zwto*(1.+zwto*(znp(1,6)+zwto*(znp(2,6)+zwto*(znp(3,6)          &
              +zwto*(znp(4,6)+zwto*znp(5,6)))))))/(1.+zwto*(zdp(1,6)           &
              +zwto*(zdp(2,6)+zwto*(zdp(3,6)+zwto*(zdp(4,6)+zwto*zdp(5,6))))))
      zdeogt(ji,jk)=zeot(ji)-zeoth(ji)
  124 continue
  126 continue
      do 125 ji=kideb,kiut
      zeot(ji)=0.
      znth(ji)=0.
      zntc(ji)=0.
      znto(ji)=0.
      zrth(ji)=0.
      zcth(ji)=0.
      zrtc(ji)=0.
      zrto(ji)=0.
  125 continue
      do 129 jk=nk,1,-1
      do 127 jn=1,iaucr
      do 127 ji=iideb(jn),iiut(jn)
      zeosb(ji)=zeos(ji)
  127 continue
      do 131 ji=kideb,kiut
      zeotb(ji)=zeot(ji)
  131 continue
      do 128 ji=kideb,kiut
      zp=max(zpeps,p(ji,jk)-dp(ji,jk))
      zq=max(zqeps,qfs(ji,jk))
      zrt=sqrt(tfs(ji,jk))
      zti=1./tfs(ji,jk)
      zt2=tfs(ji,jk)*tfs(ji,jk)
      zt2i=zti*zti
      zt4i=zt2i*zt2i
      ztdeta=zrt**49*exp(zge4/tfs(ji,jk))
      zt3=zt2*tfs(ji,jk)
      zduh=(2.*dp(ji,jk))*zq
      zduc=(2.*dp(ji,jk))*zco2
      zduo=znsor(ji)-(2.*zo31)/(1.+zo32/sqrt(zp*zp*zp))
      znsor(ji)=znsor(ji)-zduo
      znsh(ji)=znsh(ji)+zduh
      znsc(ji)=znsc(ji)+zduc
      znso(ji)=znso(ji)+zduo
      znth(ji)=znth(ji)+zduh*zti
      zntc(ji)=zntc(ji)+zduc
      znto(ji)=znto(ji)+zduo*zt2
      zduh=zduh*(zmd*pm(ji,jk))*zrt
      zduc=zduc*(zmd*pm(ji,jk))*tfs(ji,jk)
      zduo=zduo*(zmd*pm(ji,jk))*zrt
      zrsh(ji)=zrsh(ji)+zduh
      zrsc(ji)=zrsc(ji)+zduc
      zrso(ji)=zrso(ji)+zduo
      zrth(ji)=zrth(ji)+(zduh*zt2i*zrt)
      zcth(ji)=zcth(ji)+(zduh*zt2i*zrt)*zt4i*(1.+zgd4*zq*ztdeta)
      zrtc(ji)=zrtc(ji)+zduc*tfs(ji,jk)
      zrto(ji)=zrto(ji)+zduo*zt3*zrt
  128 continue
      do 132 jn=1,iaucr
      do 132 ji=iideb(jn),iiut(jn)
      zwsh=(zgas2b(1)*(zrsh(ji)/znsh(ji)))*(sqrt(1.+zg4b(1)*znsh(ji)  &
          /(zrsh(ji)/znsh(ji)))-1.)+zgc(1)*zrsh(ji)
      zwsc=(zgas2b(2)*(zrsc(ji)/znsc(ji)))*(sqrt(1.+zg4b(2)*znsc(ji)  &
          /(zrsc(ji)/znsc(ji)))-1.)+zgc(2)*zrsc(ji)
      zwso=(zgas2b(3)*(zrso(ji)/znso(ji)))*(sqrt(1.+zg4b(3)*znso(ji)  &
          /(zrso(ji)/znso(ji)))-1.)+zgc(3)*zrso(ji)
      zeos(ji)=(zwsh*(1.+zwsh*(znp(1,1)+zwsh*(znp(2,1)+zwsh*(znp(3,1)          &
              +zwsh*(znp(4,1)+zwsh*znp(5,1)))))))/(1.+zwsh*(zdp(1,1)           &
              +zwsh*(zdp(2,1)+zwsh*(zdp(3,1)+zwsh*(zdp(4,1)+zwsh*zdp(5,1)))))) &
              +(zwsc*(1.+zwsc*(znp(1,2)+zwsc*(znp(2,2)+zwsc*(znp(3,2)          &
              +zwsc*(znp(4,2)+zwsc*znp(5,2)))))))/(1.+zwsc*(zdp(1,2)           &
              +zwsc*(zdp(2,2)+zwsc*(zdp(3,2)+zwsc*(zdp(4,2)+zwsc*zdp(5,2)))))) &
              +(zwso*(1.+zwso*(znp(1,3)+zwso*(znp(2,3)+zwso*znp(3,3)))))       &
              /(1.+zwso*(zdp(1,3)+zwso*(zdp(2,3)+zwso*zdp(3,3))))
      zueogs(ji,jk)=(zeos(ji)-zeosb(ji))
  132 continue
      do 133 ji=kideb,kiut
      zwth=(zgas2b(4)*(zrth(ji)/znth(ji)))*(sqrt(1.+zg4b(4)*znth(ji)  &
          /(zrth(ji)/znth(ji)))-1.)+zgc(4)*zcth(ji)
      zwtc=(zgas2b(5)*(zrtc(ji)/zntc(ji)))*(sqrt(1.+zg4b(5)*zntc(ji)  &
          /(zrtc(ji)/zntc(ji)))-1.)+zgc(5)*zrtc(ji)
      zwto=(zgas2b(6)*(zrto(ji)/znto(ji)))*(sqrt(1.+zg4b(6)*znto(ji)  &
          /(zrto(ji)/znto(ji)))-1.)+zgc(6)*zrto(ji)
      zeot(ji)=(zwth*(1.+zwth*(znp(1,4)+zwth*(znp(2,4)+zwth*(znp(3,4)          &
              +zwth*(znp(4,4)+zwth*znp(5,4)))))))/(1.+zwth*(zdp(1,4)           &
              +zwth*(zdp(2,4)+zwth*(zdp(3,4)+zwth*(zdp(4,4)+zwth*zdp(5,4)))))) &
              +(zwtc*(1.+zwtc*(znp(1,5)+zwtc*(znp(2,5)+zwtc*(znp(3,5)          &
              +zwtc*(znp(4,5)+zwtc*znp(5,5)))))))/(1.+zwtc*(zdp(1,5)           &
              +zwtc*(zdp(2,5)+zwtc*(zdp(3,5)+zwtc*(zdp(4,5)+zwtc*zdp(5,5)))))) &
              +(zwto*(1.+zwto*(znp(1,6)+zwto*(znp(2,6)+zwto*(znp(3,6)          &
              +zwto*(znp(4,6)+zwto*znp(5,6)))))))/(1.+zwto*(zdp(1,6)           &
              +zwto*(zdp(2,6)+zwto*(zdp(3,6)+zwto*(zdp(4,6)+zwto*zdp(5,6))))))
      zueogt(ji,jk)=min(zeot(ji)-zeotb(ji),zdeogt(ji,jk))
  133 continue
  129 continue
      do 134 ji=kideb,kiut
      zeotb(ji)=zeot(ji)
  134 continue
      do 135 ji=kideb,kiut
      zp=max(zpeps,p(ji,1)-dp(ji,1))
      zq=max(zqeps,qfs(ji,1))
      zrt=sqrt(tfs(ji,1))
      zti=1./tfs(ji,1)
      zt2=tfs(ji,1)*tfs(ji,1)
      zt2i=zti*zti
      zt4i=zt2i*zt2i
      ztdeta=zrt**49*exp(zge4/tfs(ji,1))
      zt3=zt2*tfs(ji,1)
      zduh=(2.*zp)*zq
      zduc=(2.*zp)*zco2
      zduo=znsor(ji)
      znth(ji)=znth(ji)+zduh*zti
      zntc(ji)=zntc(ji)+zduc
      znto(ji)=znto(ji)+zduo*zt2
      zduh=zduh*(zmd*0.5*zp)*zrt
      zduc=zduc*(zmd*0.5*zp)*tfs(ji,1)
      zduo=zduo*(zmd*0.5*zp)*zrt
      zrth(ji)=zrth(ji)+(zduh*zt2i*zrt)
      zcth(ji)=zcth(ji)+(zduh*zt2i*zrt)*zt4i*(1.+zgd4*zq*ztdeta)
      zrtc(ji)=zrtc(ji)+zduc*tfs(ji,1)
      zrto(ji)=zrto(ji)+zduo*zt3*zrt
      zwth=(zgas2b(4)*(zrth(ji)/znth(ji)))*(sqrt(1.+zg4b(4)*znth(ji) &
          /(zrth(ji)/znth(ji)))-1.)+zgc(4)*zcth(ji)
      zwtc=(zgas2b(5)*(zrtc(ji)/zntc(ji)))*(sqrt(1.+zg4b(5)*zntc(ji) &
          /(zrtc(ji)/zntc(ji)))-1.)+zgc(5)*zrtc(ji)
      zwto=(zgas2b(6)*(zrto(ji)/znto(ji)))*(sqrt(1.+zg4b(6)*znto(ji) &
          /(zrto(ji)/znto(ji)))-1.)+zgc(6)*zrto(ji)
      zeot(ji)=(zwth*(1.+zwth*(znp(1,4)+zwth*(znp(2,4)+zwth*(znp(3,4)          &
              +zwth*(znp(4,4)+zwth*znp(5,4)))))))/(1.+zwth*(zdp(1,4)           &
              +zwth*(zdp(2,4)+zwth*(zdp(3,4)+zwth*(zdp(4,4)+zwth*zdp(5,4)))))) &
              +(zwtc*(1.+zwtc*(znp(1,5)+zwtc*(znp(2,5)+zwtc*(znp(3,5)          &
              +zwtc*(znp(4,5)+zwtc*znp(5,5)))))))/(1.+zwtc*(zdp(1,5)           &
              +zwtc*(zdp(2,5)+zwtc*(zdp(3,5)+zwtc*(zdp(4,5)+zwtc*zdp(5,5)))))) &
              +(zwto*(1.+zwto*(znp(1,6)+zwto*(znp(2,6)+zwto*(znp(3,6)          &
              +zwto*(znp(4,6)+zwto*znp(5,6)))))))/(1.+zwto*(zdp(1,6)           &
              +zwto*(zdp(2,6)+zwto*(zdp(3,6)+zwto*(zdp(4,6)+zwto*zdp(5,6))))))
      zueots(ji)=min(zeot(ji)-zeotb(ji),zeots(ji))
  135 continue

!-----------------------------------------------------------------------
!     II - CALCUL DES ELEMENTS DE GEOMETRIE NUAGEUSE POUR LE THERMIQUE.

      ikm1=nk-1
      do 200 ji=kideb,kiut
      zb1(ji,1)=1.-zneb(ji,1)
      zb3(ji,1)=1.
        if(zneb(ji,1).gt.zneb(ji,2)) then
        znmaxb(ji)=zneb(ji,1)
        znminb(ji)=zneb(ji,2)
        zb1(ji,2)=1.
        zb3(ji,2)=zneb(ji,2)/zneb(ji,1)
        else
        znmaxb(ji)=zneb(ji,2)
        znminb(ji)=zneb(ji,1)
        zb1(ji,2)=(1.-zneb(ji,2))/(1.-zneb(ji,1))
        zb3(ji,2)=1.
        endif
  200 continue
      do 203 jk=2,ikm1
      do 201 ji=kideb,kiut
      znmaxh(ji)=znmaxb(ji)
      znminh(ji)=znminb(ji)
  201 continue
      do 202 ji=kideb,kiut
      znmaxb(ji)=max(zneb(ji,jk),zneb(ji,jk+1))
      znminb(ji)=min(zneb(ji,jk),zneb(ji,jk+1))
      zb1(ji,jk+1)=(1.-znmaxb(ji))*(1./(1.-zneb(ji,jk)))
      zb2(ji,jk-1)=(1.-znmaxh(ji))*(1./(1.-zneb(ji,jk)))
      zb3(ji,jk+1)=znminb(ji)*(1./zneb(ji,jk))
      zb4(ji,jk-1)=znminh(ji)*(1./zneb(ji,jk))
  202 continue
  203 continue
      do 204 ji=kideb,kiut
      znmaxh(ji)=znmaxb(ji)
      znminh(ji)=znminb(ji)
      zb2(ji,nk-1)=(1.-znmaxh(ji))/(1.-zneb(ji,nk))
      zb4(ji,nk-1)=znminh(ji)/zneb(ji,nk)
      zb2(ji,nk)=1.
      zb4(ji,nk)=1.
  204 continue

!-----------------------------------------------------------------------
!     III - CALCUL DES COMPLEMENTS A UN DES FLUX THERMIQUES EN
!           ATMOSPHERE ISOTHERME NORMALISEE PAR ELIMINATION/SUBSTITUTION
!           PUIS UTILISATION DES RESULTATS POUR LES CALCULS DE FLUX ET
!           DE DIVERGENCES.

!     III.1   SOMMET DE L'ATMOSPHERE.

      do 310 ji=kideb,kiut
      zfdc(ji,1)=exp(max(-zeots(ji),zargli))
      zfdn(ji,1)=0.
  310 continue

!     III.2   PREMIERE COUCHE: CALCUL DES EPAISSEURS OPTIQUES, DES
!             TRANSMISSIVITES/REFLECTIVITES ET DES COEFFICIENTS LOCAUX
!             DE GEOMETRIE NUAGEUSE.

      do 320 ji=kideb,kiut
      zeo1=zdeogt(ji,1)+zeo1ta*daer(ji,1)
      zeo2=zeo2ta*daer(ji,1)
      zeps=sqrt(zeo1*zeo1-zeo2*zeo2)
      ztau=exp(max(-zeps,zargli))
      zrho=zeo2/(zeo1+zeps)
      za4c(ji)=ztau*(1.-(zrho*zrho))*(1./(1.-(zrho*zrho)*(ztau*ztau)))
      za5c(ji)=zrho*(1.-(ztau*ztau))*(1./(1.-(zrho*zrho)*(ztau*ztau)))
      zeo1 = zeo1 + zeo1tn*dp(ji,1)*zqli(ji,1) + zeo1ti*dp(ji,1)*zqice(ji,1)
      zeo2 = zeo2 + zeo2tn*dp(ji,1)*zqli(ji,1) + zeo2ti*dp(ji,1)*zqice(ji,1)
      zeps=sqrt(zeo1*zeo1-zeo2*zeo2)
      ztau=exp(max(-zeps,zargli))
      zrho=zeo2/(zeo1+zeps)
      za4n(ji)=ztau*(1.-(zrho*zrho))*(1./(1.-(zrho*zrho)*(ztau*ztau)))
      za5n(ji)=zrho*(1.-(ztau*ztau))*(1./(1.-(zrho*zrho)*(ztau*ztau)))
      zal1(ji)=zb2(ji,1)
      zal2(ji)=zb1(ji,1)
      zbe1(ji)=1.-zb2(ji,1)
      zbe2(ji)=1.-zb1(ji,1)
      zga1(ji)=1.-zb4(ji,1)
      zga2(ji)=1.-zb3(ji,1)
      zde1(ji)=zb4(ji,1)
      zde2(ji)=zb3(ji,1)
  320 continue

!     III.3   PREMIERE COUCHE, ELIMINATION (SANS PROBLEMES).

      do 330 ji=kideb,kiut
      zfmc(ji,1)=za5c(ji)*(zal2(ji)*zfdc(ji,1))
      zfdc(ji,2)=za4c(ji)*(zal2(ji)*zfdc(ji,1))
      zfmn(ji,1)=za5n(ji)*(zbe2(ji)*zfdc(ji,1))
      zfdn(ji,2)=za4n(ji)*(zbe2(ji)*zfdc(ji,1))
      ztu1(ji,1)=0.
      ztu2(ji,1)=zal1(ji)*za4c(ji)
      ztu3(ji,1)=zga1(ji)*za4c(ji)
      ztu4(ji,1)=zbe1(ji)*za4n(ji)
      ztu5(ji,1)=zde1(ji)*za4n(ji)
      ztu6(ji,1)=zal1(ji)*za5c(ji)
      ztu7(ji,1)=zga1(ji)*za5c(ji)
      ztu8(ji,1)=zbe1(ji)*za5n(ji)
      ztu9(ji,1)=zde1(ji)*za5n(ji)
  330 continue

!     III.4   BOUCLE SUR LES NIVEAUX, CALCULS PRELIMINAIRES PUIS
!             ELIMINATION.

      do 343 jk=2,nk
      do 340 ji=kideb,kiut
      zeo1=zdeogt(ji,jk)+zeo1ta*daer(ji,jk)
      zeo2=zeo2ta*daer(ji,jk)
      zeps=sqrt(zeo1*zeo1-zeo2*zeo2)
      ztau=exp(max(-zeps,zargli))
      zrho=zeo2/(zeo1+zeps)
      za4c(ji)=ztau*(1.-(zrho*zrho))*(1./(1.-(zrho*zrho)*(ztau*ztau)))
      za5c(ji)=zrho*(1.-(ztau*ztau))*(1./(1.-(zrho*zrho)*(ztau*ztau)))
      zeo1 = zeo1 + zeo1tn*dp(ji,jk)*zqli(ji,jk) + zeo1ti*dp(ji,jk)*zqice(ji,jk)
      zeo2 = zeo2 + zeo2tn*dp(ji,jk)*zqli(ji,jk) + zeo2ti*dp(ji,jk)*zqice(ji,jk)
      zeps=sqrt(zeo1*zeo1-zeo2*zeo2)
      ztau=exp(max(-zeps,zargli))
      zrho=zeo2/(zeo1+zeps)
      za4n(ji)=ztau*(1.-(zrho*zrho))*(1./(1.-(zrho*zrho)*(ztau*ztau)))
      za5n(ji)=zrho*(1.-(ztau*ztau))*(1./(1.-(zrho*zrho)*(ztau*ztau)))
      zal1(ji)=zb2(ji,jk)
      zal2(ji)=zb1(ji,jk)
      zbe1(ji)=1.-zb2(ji,jk)
      zbe2(ji)=1.-zb1(ji,jk)
      zga1(ji)=1.-zb4(ji,jk)
      zga2(ji)=1.-zb3(ji,jk)
      zde1(ji)=zb4(ji,jk)
      zde2(ji)=zb3(ji,jk)
  340 continue
      do 341 ji=kideb,kiut
      ztd1=1./(1.-za5c(ji)*(zal2(ji)*ztu6(ji,jk-1)+zga2(ji)*ztu8(ji,jk-1)))
      zfmc(ji,jk)=(ztd1*za5c(ji))*(zal2(ji)*zfdc(ji,jk)+zga2(ji)*zfdn(ji,jk))
      ztu1(ji,jk)=(ztd1*za5c(ji))*(zal2(ji)*ztu7(ji,jk-1)+zga2(ji)*ztu9(ji,jk-1))
      ztu2(ji,jk)=(ztd1*za4c(ji))*zal1(ji)
      ztu3(ji,jk)=(ztd1*za4c(ji))*zga1(ji)
      ztd2=za5n(ji)*(zbe2(ji)*ztu6(ji,jk-1)+zde2(ji)*ztu8(ji,jk-1))
      ztd3=1./(1.-za5n(ji)*(zbe2(ji)*ztu7(ji,jk-1)+zde2(ji)*ztu9(ji,jk-1))-ztd2*ztu1(ji,jk))
      zfmn(ji,jk)=ztd3*(za5n(ji)*(zbe2(ji)*zfdc(ji,jk)+zde2(ji)*zfdn(ji,jk))+ztd2*zfmc(ji,jk))
      ztu4(ji,jk)=ztd3*(za4n(ji)*zbe1(ji)+ztd2*ztu2(ji,jk))
      ztu5(ji,jk)=ztd3*(za4n(ji)*zde1(ji)+ztd2*ztu3(ji,jk))
  341 continue
      do 342 ji=kideb,kiut
      ztd4=za4c(ji)*(zal2(ji)*ztu6(ji,jk-1)+zga2(ji)*ztu8(ji,jk-1))
      ztd5=za4c(ji)*(zal2(ji)*ztu7(ji,jk-1)+zga2(ji)*ztu9(ji,jk-1))
      zfdc(ji,jk+1)=za4c(ji)*(zal2(ji)*zfdc(ji,jk)+zga2(ji)*zfdn(ji,jk)) &
                   +ztd4*zfmc(ji,jk)+ztd5*zfmn(ji,jk)
      ztu6(ji,jk)=za5c(ji)*zal1(ji)+ztd4*ztu2(ji,jk)+ztd5*ztu4(ji,jk)
      ztu7(ji,jk)=za5c(ji)*zga1(ji)+ztd4*ztu3(ji,jk)+ztd5*ztu5(ji,jk)
      ztd6=za4n(ji)*(zbe2(ji)*ztu6(ji,jk-1)+zde2(ji)*ztu8(ji,jk-1))
      ztd7=za4n(ji)*(zbe2(ji)*ztu7(ji,jk-1)+zde2(ji)*ztu9(ji,jk-1))
      zfdn(ji,jk+1)=za4n(ji)*(zbe2(ji)*zfdc(ji,jk)+zde2(ji)*zfdn(ji,jk)) &
                   +ztd6*zfmc(ji,jk)+ztd7*zfmn(ji,jk)
      ztu8(ji,jk)=za5n(ji)*zbe1(ji)+ztd6*ztu2(ji,jk)+ztd7*ztu4(ji,jk)
      ztu9(ji,jk)=za5n(ji)*zde1(ji)+ztd6*ztu3(ji,jk)+ztd7*ztu5(ji,jk)
  342 continue
  343 continue

!     III.5   TRAITEMENT EN SURFACE.

      do 350 ji=kideb,kiut
      zal=1.-emis(ji)
      ztds1=1./(1.-zal*ztu6(ji,nk))
      zfmc(ji,nk+1)=(ztds1*zal)*zfdc(ji,nk+1)
      ztus1=(ztds1*zal)*ztu7(ji,nk)
      ztds2=zal*ztu8(ji,nk)
      ztds3=1./(1.-zal*ztu9(ji,nk)-ztds2*ztus1)
      zfmn(ji,nk+1)=ztds3*(zal*zfdn(ji,nk+1)+ztds2*zfmc(ji,nk+1))
      zfmc(ji,nk+1)=zfmc(ji,nk+1)+ztus1*zfmn(ji,nk+1)
  350 continue

!     III.6   SUBSTITUTION NIVEAU PAR NIVEAU.

      do 361 jk=nk,1,-1
      do 360 ji=kideb,kiut
      zfdn(ji,jk+1)=zfdn(ji,jk+1)+ztu8(ji,jk)*zfmc(ji,jk+1)+ztu9(ji,jk)*zfmn(ji,jk+1)
      zfdc(ji,jk+1)=zfdc(ji,jk+1)+ztu6(ji,jk)*zfmc(ji,jk+1)+ztu7(ji,jk)*zfmn(ji,jk+1)
      zfmn(ji,jk)=zfmn(ji,jk)+ztu4(ji,jk)*zfmc(ji,jk+1)+ztu5(ji,jk)*zfmn(ji,jk+1)
      zfmc(ji,jk)=zfmc(ji,jk)+ztu2(ji,jk)*zfmc(ji,jk+1)+ztu3(ji,jk)*zfmn(ji,jk+1)+ztu1(ji,jk)*zfmn(ji,jk)
  360 continue
  361 continue

!     III.7   CALCUL DES FLUX ET DES DELTA-FLUX AVEC L'HYPOTHESE DU
!             COOLING TO SPACE.

      do 370 ji=kideb,kiut
      ztrans(ji)=zfdc(ji,nk+1)-zfmc(ji,nk+1)
      ztrans(ji)=ztrans(ji)+zfdn(ji,nk+1)-zfmn(ji,nk+1)
      zfnet(ji)=-zbb(ji,nk+1)*ztrans(ji)
      rat(ji)=zfnet(ji)/emis(ji)+zbb(ji,nk+1)
  370 continue
      do 373 jk=nk,1,-1
      do 371 ji=kideb,kiut
      ztranb(ji)=ztrans(ji)
  371 continue
      do 372 ji=kideb,kiut
      ztrans(ji)=zfdc(ji,jk)-zfmc(ji,jk)
      ztrans(ji)=ztrans(ji)+zfdn(ji,jk)-zfmn(ji,jk)
      dtfr(ji,jk)=zbb(ji,jk)*(ztranb(ji)-ztrans(ji))
      zfnet(ji)=zfnet(ji)+dtfr(ji,jk)
  372 continue
  373 continue
      do 374 ji=kideb,kiut
      rirs(ji)=zfnet(ji)
  374 continue

!-----------------------------------------------------------------------
!     XII - RECALCUL DE LA GEOMETRIE NUAGEUSE POUR LE THERMIQUE.

      ikm1=nk-1
      do 1200 ji=kideb,kiut
      zb1(ji,1)=1.-zneb(ji,1)
      zb3(ji,1)=1.
      znmaxb(ji)=max(zneb(ji,1),zneb(ji,2))
      znminb(ji)=min(zneb(ji,1),zneb(ji,2))
      zb1(ji,2)=(1.-znmaxb(ji))/(1.-zneb(ji,1))
      zb3(ji,2)=znminb(ji)/zneb(ji,1)
 1200 continue
      do 1203 jk=2,ikm1
      do 1201 ji=kideb,kiut
      znmaxh(ji)=znmaxb(ji)
      znminh(ji)=znminb(ji)
 1201 continue
      do 1202 ji=kideb,kiut
      znmaxb(ji)=max(zneb(ji,jk),zneb(ji,jk+1))
      znminb(ji)=min(zneb(ji,jk),zneb(ji,jk+1))
      zb1(ji,jk+1)=(1.-znmaxb(ji))*(1./(1.-zneb(ji,jk)))
      zb2(ji,jk-1)=(1.-znmaxh(ji))*(1./(1.-zneb(ji,jk)))
      zb3(ji,jk+1)=znminb(ji)*(1./zneb(ji,jk))
      zb4(ji,jk-1)=znminh(ji)*(1./zneb(ji,jk))
 1202 continue
 1203 continue
      do 1204 ji=kideb,kiut
      znmaxh(ji)=znmaxb(ji)
      znminh(ji)=znminb(ji)
      zb2(ji,nk-1)=(1.-znmaxh(ji))/(1.-zneb(ji,nk))
      zb4(ji,nk-1)=znminh(ji)/zneb(ji,nk)
      zb2(ji,nk)=1.
      zb4(ji,nk)=1.
 1204 continue

!-----------------------------------------------------------------------
!     XIII - CALCUL DES COMPLEMENTS AU FLUX DU CORPS NOIR DES FLUX
!            THERMIQUES 'ANTI-SURESTIMATION' POUR LE PROFIL DE
!            TEMPERATURE REEL ET LE PROFIL ISOTHERME NORMALISE AFIN
!            D'OBTENIR PAR DIFFERENCE LES TERMES D'ECHANGE QUI
!            AJOUTES AU RESULTATS DU COOLING TO SPACE DONNERONT LE
!            CHAMP DE RAYONNEMENT FINAL.

!     XIII.1   SOMMET DE L'ATMOSPHERE.

      do 1310 ji=kideb,kiut
      zfdc(ji,1)=exp(max(-zueots(ji),zargli))
      zzfdc(ji,1)=zbb(ji,1)*zfdc(ji,1)
      zfdn(ji,1)=0.
      zzfdn(ji,1)=0.
 1310 continue

!     XIII.2   PREMIERE COUCHE: CALCUL DES EPAISSEURS OPTIQUES, DES
!              TRANSMISSIVITES/REFLECTIVITES ET DES COEFFICIENTS LOCAUX
!              DE GEOMETRIE NUAGEUSE.

      do 1320 ji=kideb,kiut
      znebl(ji)=zneb(ji,1)
      znebcl(ji)=1.-zneb(ji,1)
      zeo1=zueogt(ji,1)+zeo1ta*daer(ji,1)
      zeo2=zeo2ta*daer(ji,1)
      zeps=sqrt(zeo1*zeo1-zeo2*zeo2)
      ztau=exp(max(-zeps,zargli))
      zrho=zeo2/(zeo1+zeps)
      za4c(ji)=ztau*(1.-(zrho*zrho))*(1./(1.-(zrho*zrho)*(ztau*ztau)))
      za5c(ji)=zrho*(1.-(ztau*ztau))*(1./(1.-(zrho*zrho)*(ztau*ztau)))
      zeo1 = zeo1 + zeo1tn*dp(ji,1)*zqli(ji,1) + zeo1ti*dp(ji,1)*zqice(ji,1)
      zeo2 = zeo2 + zeo2tn*dp(ji,1)*zqli(ji,1) + zeo2ti*dp(ji,1)*zqice(ji,1)
      zeps=sqrt(zeo1*zeo1-zeo2*zeo2)
      ztau=exp(max(-zeps,zargli))
      zrho=zeo2/(zeo1+zeps)
      za4n(ji)=ztau*(1.-(zrho*zrho))*(1./(1.-(zrho*zrho)*(ztau*ztau)))
      za5n(ji)=zrho*(1.-(ztau*ztau))*(1./(1.-(zrho*zrho)*(ztau*ztau)))
      zal1(ji)=zb2(ji,1)
      zal2(ji)=zb1(ji,1)
      zbe1(ji)=1.-zb2(ji,1)
      zbe2(ji)=1.-zb1(ji,1)
      zga1(ji)=1.-zb4(ji,1)
      zga2(ji)=1.-zb3(ji,1)
      zde1(ji)=zb4(ji,1)
      zde2(ji)=zb3(ji,1)
 1320 continue

!     XIII.3   PREMIERE COUCHE, ELIMINATION (SANS PROBLEMES).

      do 1330 ji=kideb,kiut
      zfmc(ji,1)=za5c(ji)*(zal2(ji)*zfdc(ji,1))
      zzfmc(ji,1)=za5c(ji)*(zal2(ji)*zzfdc(ji,1))+za4c(ji)*znebcl(ji)*(zbb(ji,1)-zbb(ji,2))
      zfdc(ji,2)=za4c(ji)*(zal2(ji)*zfdc(ji,1))
      zzfdc(ji,2)=za4c(ji)*(zal2(ji)*zzfdc(ji,1))+za5c(ji)*znebcl(ji)*(zbb(ji,1)-zbb(ji,2))
      zfmn(ji,1)=za5n(ji)*(zbe2(ji)*zfdc(ji,1))
      zzfmn(ji,1)=za5n(ji)*(zbe2(ji)*zzfdc(ji,1))+za4n(ji)*znebl (ji)*(zbb(ji,1)-zbb(ji,2))
      zfdn(ji,2)=za4n(ji)*(zbe2(ji)*zfdc(ji,1))
      zzfdn(ji,2)=za4n(ji)*(zbe2(ji)*zzfdc(ji,1))+za5n(ji)*znebl (ji)*(zbb(ji,1)-zbb(ji,2))
      ztu1(ji,1)=0.
      ztu2(ji,1)=zal1(ji)*za4c(ji)
      ztu3(ji,1)=zga1(ji)*za4c(ji)
      ztu4(ji,1)=zbe1(ji)*za4n(ji)
      ztu5(ji,1)=zde1(ji)*za4n(ji)
      ztu6(ji,1)=zal1(ji)*za5c(ji)
      ztu7(ji,1)=zga1(ji)*za5c(ji)
      ztu8(ji,1)=zbe1(ji)*za5n(ji)
      ztu9(ji,1)=zde1(ji)*za5n(ji)
 1330 continue

!     XIII.4   BOUCLE SUR LES NIVEAUX, CALCULS PRELIMINAIRES PUIS ELIMINATION.

      do 1343 jk=2,nk
      do 1340 ji=kideb,kiut
      znebl(ji)=zneb(ji,jk)
      znebcl(ji)=1.-zneb(ji,jk)
      zeo1=zueogt(ji,jk)+zeo1ta*daer(ji,jk)
      zeo2=zeo2ta*daer(ji,jk)
      zeps=sqrt(zeo1*zeo1-zeo2*zeo2)
      ztau=exp(max(-zeps,zargli))
      zrho=zeo2/(zeo1+zeps)
      za4c(ji)=ztau*(1.-(zrho*zrho))*(1./(1.-(zrho*zrho)*(ztau*ztau)))
      za5c(ji)=zrho*(1.-(ztau*ztau))*(1./(1.-(zrho*zrho)*(ztau*ztau)))
      zeo1 = zeo1 + zeo1tn*dp(ji,jk)*zqli(ji,jk) + zeo1ti*dp(ji,jk)*zqice(ji,jk)
      zeo2 = zeo2 + zeo2tn*dp(ji,jk)*zqli(ji,jk) + zeo2ti*dp(ji,jk)*zqice(ji,jk)
      zeps=sqrt(zeo1*zeo1-zeo2*zeo2)
      ztau=exp(max(-zeps,zargli))
      zrho=zeo2/(zeo1+zeps)
      za4n(ji)=ztau*(1.-(zrho*zrho))*(1./(1.-(zrho*zrho)*(ztau*ztau)))
      za5n(ji)=zrho*(1.-(ztau*ztau))*(1./(1.-(zrho*zrho)*(ztau*ztau)))
      zal1(ji)=zb2(ji,jk)
      zal2(ji)=zb1(ji,jk)
      zbe1(ji)=1.-zb2(ji,jk)
      zbe2(ji)=1.-zb1(ji,jk)
      zga1(ji)=1.-zb4(ji,jk)
      zga2(ji)=1.-zb3(ji,jk)
      zde1(ji)=zb4(ji,jk)
      zde2(ji)=zb3(ji,jk)
 1340 continue
      do 1341 ji=kideb,kiut
      ztd1=1./(1.-za5c(ji)*(zal2(ji)*ztu6(ji,jk-1)+zga2(ji)*ztu8(ji,jk-1)))
      zfmc(ji,jk)=(ztd1*za5c(ji))*(zal2(ji)*zfdc(ji,jk)+zga2(ji)*zfdn(ji,jk))
      zzfmc(ji,jk)=(ztd1*za5c(ji))*(zal2(ji)*zzfdc(ji,jk)+zga2(ji) &
                  *zzfdn(ji,jk))+ztd1*znebcl(ji)*(za4c(ji)         &
                  *(zbb(ji,jk)-zbb(ji,jk+1))+za5c(ji)*(zbb(ji,jk)-zbb(ji,jk-1)))
      ztu1(ji,jk)=(ztd1*za5c(ji))*(zal2(ji)*ztu7(ji,jk-1)+zga2(ji)*ztu9(ji,jk-1))
      ztu2(ji,jk)=(ztd1*za4c(ji))*zal1(ji)
      ztu3(ji,jk)=(ztd1*za4c(ji))*zga1(ji)
      ztd2=za5n(ji)*(zbe2(ji)*ztu6(ji,jk-1)+zde2(ji)*ztu8(ji,jk-1))
      ztd3=1./(1.-za5n(ji)*(zbe2(ji)*ztu7(ji,jk-1)+zde2(ji)*ztu9(ji,jk-1))-ztd2*ztu1(ji,jk))
      zfmn(ji,jk)=ztd3*(za5n(ji)*(zbe2(ji)*zfdc(ji,jk)+zde2(ji)*zfdn(ji,jk))+ztd2*zfmc(ji,jk))
      zzfmn(ji,jk)=ztd3*(za5n(ji)*(zbe2(ji)*zzfdc(ji,jk)+zde2(ji)  &
                  *zzfdn(ji,jk))+ztd2*zzfmc(ji,jk))+ztd3*znebl(ji) &
                  *(za4n(ji)*(zbb(ji,jk)-zbb(ji,jk+1))+za5n(ji)*(zbb(ji,jk)-zbb(ji,jk-1)))
      ztu4(ji,jk)=ztd3*(za4n(ji)*zbe1(ji)+ztd2*ztu2(ji,jk))
      ztu5(ji,jk)=ztd3*(za4n(ji)*zde1(ji)+ztd2*ztu3(ji,jk))
 1341 continue
      do 1342 ji=kideb,kiut
      ztd4=za4c(ji)*(zal2(ji)*ztu6(ji,jk-1)+zga2(ji)*ztu8(ji,jk-1))
      ztd5=za4c(ji)*(zal2(ji)*ztu7(ji,jk-1)+zga2(ji)*ztu9(ji,jk-1))
      zfdc(ji,jk+1)=za4c(ji)*(zal2(ji)*zfdc(ji,jk)+zga2(ji)*zfdn(ji,jk))+ztd4*zfmc(ji,jk)+ztd5*zfmn(ji,jk)
      zzfdc(ji,jk+1)=za4c(ji)*(zal2(ji)*zzfdc(ji,jk)+zga2(ji)            &
                    *zzfdn(ji,jk))+ztd4*zzfmc(ji,jk)+ztd5*zzfmn(ji,jk)   &
                    +znebcl(ji)*(za5c(ji)*(zbb(ji,jk)                    &
                    -zbb(ji,jk+1))+za4c(ji)*(zbb(ji,jk)-zbb(ji,jk-1)))
      ztu6(ji,jk)=za5c(ji)*zal1(ji)+ztd4*ztu2(ji,jk)+ztd5*ztu4(ji,jk)
      ztu7(ji,jk)=za5c(ji)*zga1(ji)+ztd4*ztu3(ji,jk)+ztd5*ztu5(ji,jk)
      ztd6=za4n(ji)*(zbe2(ji)*ztu6(ji,jk-1)+zde2(ji)*ztu8(ji,jk-1))
      ztd7=za4n(ji)*(zbe2(ji)*ztu7(ji,jk-1)+zde2(ji)*ztu9(ji,jk-1))
      zfdn(ji,jk+1)=za4n(ji)*(zbe2(ji)*zfdc(ji,jk)+zde2(ji)*zfdn(ji,jk))+ztd6*zfmc(ji,jk)+ztd7*zfmn(ji,jk)
      zzfdn(ji,jk+1)=za4n(ji)*(zbe2(ji)*zzfdc(ji,jk)+zde2(ji)            &
                    *zzfdn(ji,jk))+ztd6*zzfmc(ji,jk)+ztd7*zzfmn(ji,jk)   &
                    +znebl(ji)*(za5n(ji)*(zbb(ji,jk)-zbb(ji,jk+1))       &
                    +za4n(ji)*(zbb(ji,jk)-zbb(ji,jk-1)))
      ztu8(ji,jk)=za5n(ji)*zbe1(ji)+ztd6*ztu2(ji,jk)+ztd7*ztu4(ji,jk)
      ztu9(ji,jk)=za5n(ji)*zde1(ji)+ztd6*ztu3(ji,jk)+ztd7*ztu5(ji,jk)
 1342 continue
 1343 continue

!     XIII.5   TRAITEMENT EN SURFACE.

      do 1350 ji=kideb,kiut
      zal=1.-emis(ji)
      ztds1=1./(1.-zal*ztu6(ji,nk))
      zfmc(ji,nk+1)=(ztds1*zal)*zfdc(ji,nk+1)
      zzfmc(ji,nk+1)=(ztds1*zal)*zzfdc(ji,nk+1)+(ztds1*zal)*(zbb(ji,nk+1)-zbb(ji,nk))
      zzfmc(ji,nk+1)=zzfmc(ji,nk+1)-(ztds1*zal)*znebl(ji)  *(zbb(ji,nk+1)-zbb(ji,nk))
      ztus1=(ztds1*zal)*ztu7(ji,nk)
      ztds2=zal*ztu8(ji,nk)
      ztds3=1./(1.-zal*ztu9(ji,nk)-ztds2*ztus1)
      zfmn(ji,nk+1)=ztds3*(zal*zfdn(ji,nk+1)+ztds2*zfmc(ji,nk+1))
      zzfmn(ji,nk+1)=ztds3*(zal*zzfdn(ji,nk+1)+ztds2*zzfmc(ji,nk+1)) &
                    +ztds3*zal*znebl(ji)*(zbb(ji,nk+1)-zbb(ji,nk))
      zfmc(ji,nk+1)=zfmc(ji,nk+1)+ztus1*zfmn(ji,nk+1)
      zzfmc(ji,nk+1)=zzfmc(ji,nk+1)+ztus1*zzfmn(ji,nk+1)
 1350 continue

!     XIII.6   SUBSTITUTION NIVEAU PAR NIVEAU.

      do 1361 jk=nk,1,-1
      do 1360 ji=kideb,kiut
      zfdn(ji,jk+1)=zfdn(ji,jk+1)+ztu8(ji,jk)*zfmc(ji,jk+1)+ztu9(ji,jk)*zfmn(ji,jk+1)
      zzfdn(ji,jk+1)=zzfdn(ji,jk+1)+ztu8(ji,jk)*zzfmc(ji,jk+1)+ztu9(ji,jk)*zzfmn(ji,jk+1)
      zfdc(ji,jk+1)=zfdc(ji,jk+1)+ztu6(ji,jk)*zfmc(ji,jk+1)+ztu7(ji,jk)*zfmn(ji,jk+1)
      zzfdc(ji,jk+1)=zzfdc(ji,jk+1)+ztu6(ji,jk)*zzfmc(ji,jk+1)+ztu7(ji,jk)*zzfmn(ji,jk+1)
      zfmn(ji,jk)=zfmn(ji,jk)+ztu4(ji,jk)*zfmc(ji,jk+1)+ztu5(ji,jk)*zfmn(ji,jk+1)
      zzfmn(ji,jk)=zzfmn(ji,jk)+ztu4(ji,jk)*zzfmc(ji,jk+1)+ztu5(ji,jk)*zzfmn(ji,jk+1)
      zfmc(ji,jk)=zfmc(ji,jk)+ztu2(ji,jk)*zfmc(ji,jk+1)+ztu3(ji,jk)     &
                 *zfmn(ji,jk+1)+ztu1(ji,jk)*zfmn(ji,jk)
      zzfmc(ji,jk)=zzfmc(ji,jk)+ztu2(ji,jk)*zzfmc(ji,jk+1)+ztu3(ji,jk)  &
                  *zzfmn(ji,jk+1)+ztu1(ji,jk)*zzfmn(ji,jk)
 1360 continue
 1361 continue

!     XIII.7   CALCUL DES FLUX ET DES DELTA-FLUX DEFINITIFS.

      do 1370 ji=kideb,kiut
      ztrans(ji)=zfdc(ji,nk+1)-zfmc(ji,nk+1)
      zztrans(ji)=zzfdc(ji,nk+1)-zzfmc(ji,nk+1)+zbb(ji,nk+1)-zbb(ji,nk)
      ztrans(ji)=ztrans(ji)+zfdn(ji,nk+1)-zfmn(ji,nk+1)
      zztrans(ji)=zztrans(ji)+zzfdn(ji,nk+1)-zzfmn(ji,nk+1)
      zfnet(ji)=(rat(ji)-zbb(ji,nk+1))*emis(ji)+zbb(ji,nk+1)*ztrans(ji)-zztrans(ji)

!  ******** OUR CORRECTION **********
!  (ORIGINAL IN THE FIRST LINE)

!        rat(ji)=zfnet(ji)/emis(ji)+zbb(ji,nk+1)
         rat(ji)=zfnet(ji)
 1370 continue
      do 1373 jk=nk,2,-1
      do 1371 ji=kideb,kiut
      ztranb(ji)=ztrans(ji)
      zztranb(ji)=zztrans(ji)
 1371 continue
      do 1372 ji=kideb,kiut
      ztrans(ji)=zfdc(ji,jk)-zfmc(ji,jk)
      zztrans(ji)=zzfdc(ji,jk)-zzfmc(ji,jk)+zbb(ji,jk)-zbb(ji,jk-1)
      ztrans(ji)=ztrans(ji)+zfdn(ji,jk)-zfmn(ji,jk)
      zztrans(ji)=zztrans(ji)+zzfdn(ji,jk)-zzfmn(ji,jk)
      dtfr(ji,jk)=dtfr(ji,jk)-zbb(ji,jk)*(ztranb(ji)-ztrans(ji))+(zztranb(ji)-zztrans(ji))
      zfnet(ji)=zfnet(ji)+dtfr(ji,jk)
      dtfr(ji,jk)=dtfr(ji,jk)*gg/(dp(ji,jk)*cph(ji,jk))
 1372 continue
 1373 continue
      do 1374 ji=kideb,kiut
      ztranb(ji)=ztrans(ji)
      zztranb(ji)=zztrans(ji)
 1374 continue
      do 1375 ji=kideb,kiut
      ztrans(ji)=zfdc(ji,1)-zfmc(ji,1)
      zztrans(ji)=zzfdc(ji,1)-zzfmc(ji,1)
      ztrans(ji)=ztrans(ji)+zfdn(ji,1)-zfmn(ji,1)
      zztrans(ji)=zztrans(ji)+zzfdn(ji,1)-zzfmn(ji,1)
      dtfr(ji,1)=dtfr(ji,1)-zbb(ji,1)*(ztranb(ji)-ztrans(ji))+(zztranb(ji)-zztrans(ji))
      zfnet(ji)=zfnet(ji)+dtfr(ji,1)
      dtfr(ji,1)=dtfr(ji,1)*gg/(dp(ji,1)*cph(ji,1))
 1375 continue
      do 1376 ji=kideb,kiut
      rirs(ji)=zfnet(ji)
 1376 continue

!-----------------------------------------------------------------------
!     IV - CALCUL DES TABLEAUX DEPENDANT DE L'ANGLE SOLAIRE.

      do 400 jn=1,iaucr
      do 400 ji=iideb(jn),iiut(jn)
      zmu0i(ji)=(sqrt(zmu0(ji)*zmu0(ji)+zzrae)-zmu0(ji))*zirae
      zusa(ji)=(0.5+zusaa*zmu0(ji))/(1.+zusba*zmu0(ji))
      zusn(ji)=(0.5+zusan*zmu0(ji))/(1.+zusbn*zmu0(ji))
      zusi(ji)=(0.5+zusai*zmu0(ji))/(1.+zusbi*zmu0(ji))
  400 continue

!-----------------------------------------------------------------------
!     V - RECALCUL DE LA GEOMETRIE NUAGEUSE POUR LE SOLAIRE.

      ikm1=nk-1
      do 500 jn=1,iaucr
      do 500 ji=iideb(jn),iiut(jn)
      zb1(ji,1)=1.-zneb(ji,1)
      zb3(ji,1)=1.
      znmaxb(ji)=max(zneb(ji,1),zneb(ji,2))
      znminb(ji)=min(zneb(ji,1),zneb(ji,2))
      zb1(ji,2)=(1.-znmaxb(ji))/(1.-zneb(ji,1))
      zb3(ji,2)=znminb(ji)/zneb(ji,1)
  500 continue
      do 503 jk=2,ikm1
      do 501 jn=1,iaucr
      do 501 ji=iideb(jn),iiut(jn)
      znmaxh(ji)=znmaxb(ji)
      znminh(ji)=znminb(ji)
  501 continue
      do 502 jn=1,iaucr
      do 502 ji=iideb(jn),iiut(jn)
      znmaxb(ji)=max(zneb(ji,jk),zneb(ji,jk+1))
      znminb(ji)=min(zneb(ji,jk),zneb(ji,jk+1))
      zb1(ji,jk+1)=(1.-znmaxb(ji))*(1./(1.-zneb(ji,jk)))
      zb2(ji,jk-1)=(1.-znmaxh(ji))*(1./(1.-zneb(ji,jk)))
      zb3(ji,jk+1)=znminb(ji)*(1./zneb(ji,jk))
      zb4(ji,jk-1)=znminh(ji)*(1./zneb(ji,jk))
  502 continue
  503 continue
      do 504 jn=1,iaucr
      do 504 ji=iideb(jn),iiut(jn)
      znmaxh(ji)=znmaxb(ji)
      znminh(ji)=znminb(ji)
      zb2(ji,nk-1)=(1.-znmaxh(ji))/(1.-zneb(ji,nk))
      zb4(ji,nk-1)=znminh(ji)/zneb(ji,nk)
      zb2(ji,nk)=1.
      zb4(ji,nk)=1.
  504 continue

!-----------------------------------------------------------------------
!     VI - CALCUL DES FLUX SOLAIRES PAR ELIMINATION/SUBSTITUTION PUIS
!          UTILISATION DES RESULTATS POUR LES CALCULS DE DIVERGENCES.

!     VI.1   SOMMET DE L'ATMOSPHERE.

      do 610 jn=1,iaucr
      do 610 ji=iideb(jn),iiut(jn)
      zfpc(ji,1)=zroxts*zmu0(ji)*exp(max(-0.5*zmu0i(ji)*zeoss(ji),zargli))
      zfdc(ji,1)=0.
      zfpn(ji,1)=0.
      zfdn(ji,1)=0.
  610 continue

!     VI.2   PREMIERE COUCHE: CALCUL DES EPAISSEURS OPTIQUES, DES
!            TRANSMISSIVITES/REFLECTIVITES ET DES COEFFICIENTS LOCAUX
!            DE GEOMETRIE NUAGEUSE.

      do 620 jn=1,iaucr
      do 620 ji=iideb(jn),iiut(jn)
      zeo1=zueogs(ji,1)+zeo1sa*daer(ji,1)+(zeoray*dp(ji,1))
      zeo2=zeo2sa*daer(ji,1)+(zeoray*dp(ji,1))
      zeps=sqrt(zeo1*zeo1-zeo2*zeo2)
      ztau=exp(max(-zeps,zargli))
      zrho=zeo2/(zeo1+zeps)
      za4c(ji)=ztau*(1.-(zrho*zrho))*(1./(1.-(zrho*zrho)*(ztau*ztau)))
      za5c(ji)=zrho*(1.-(ztau*ztau))*(1./(1.-(zrho*zrho)*(ztau*ztau)))
      zeo=0.5*zdeogs(ji,1)+zeosa*daer(ji,1)+(zeoray*dp(ji,1))
      zmu0ic=(zeps/zeo)+dsign(max(dabs(zmu0i(ji)-(zeps/zeo)),zepres),zmu0i(ji)-(zeps/zeo))
      zeo3=((zeodsa*daer(ji,1))*zusa(ji)+(0.5*(zeoray*dp(ji,1))))*zmu0ic
      zeo4=((zeodsa*daer(ji,1))*(1.-zusa(ji))+(0.5*(zeoray*dp(ji,1))))*zmu0ic
      zeo5=zeo*zmu0ic
      za1c(ji)=exp(max(-zeo5,zargli))
      zg1=(zeo3*(zeo5-zeo1)-zeo2*zeo4)*(1./(zeo5*zeo5-zeps*zeps))
      zg2=-(zeo4*(zeo5+zeo1)+zeo2*zeo3)*(1./(zeo5*zeo5-zeps*zeps))
      za2c(ji)=zg2*(za1c(ji)-za4c(ji))-zg1*za5c(ji)*za1c(ji)
      za3c(ji)=zg1*(1.-za4c(ji)*za1c(ji))-zg2*za5c(ji)
      zeo1 = zeo1 + zeo1sn*dp(ji,1)*zqli(ji,1) + zeo1si*dp(ji,1)*zqice(ji,1)
      zeo2 = zeo2 + zeo2sn*dp(ji,1)*zqli(ji,1) + zeo2si*dp(ji,1)*zqice(ji,1)
      zeps=sqrt(zeo1*zeo1-zeo2*zeo2)
      ztau=exp(max(-zeps,zargli))
      zrho=zeo2/(zeo1+zeps)
      za4n(ji)=ztau*(1.-(zrho*zrho))*(1./(1.-(zrho*zrho)*(ztau*ztau)))
      za5n(ji)=zrho*(1.-(ztau*ztau))*(1./(1.-(zrho*zrho)*(ztau*ztau)))
      zeo = zeo + zeosn*dp(ji,1)*zqli(ji,1) + zeosi*dp(ji,1)*zqice(ji,1)
      zmu0in=(zeps/zeo)+dsign(max(dabs(zmu0i(ji)-(zeps/zeo)),zepres),zmu0i(ji)-(zeps/zeo))
      zeo3=(zeo3/zmu0ic+zusn(ji)*((dp(ji,1)*zqli(ji,1))*zeodsn)        &
                       +zusi(ji)*((dp(ji,1)*zqice(ji,1))*zeodsi))*zmu0in
      zeo4=(zeo4/zmu0ic+(1.-zusn(ji))*((dp(ji,1)*zqli(ji,1))*zeodsn)   &
                       +(1.-zusi(ji))*((dp(ji,1)*zqice(ji,1))*zeodsi))*zmu0in
      zeo5=zeo*zmu0in
      za1n(ji)=exp(max(-zeo5,zargli))
      zg1=(zeo3*(zeo5-zeo1)-zeo2*zeo4)*(1./(zeo5*zeo5-zeps*zeps))
      zg2=-(zeo4*(zeo5+zeo1)+zeo2*zeo3)*(1./(zeo5*zeo5-zeps*zeps))
      za2n(ji)=zg2*(za1n(ji)-za4n(ji))-zg1*za5n(ji)*za1n(ji)
      za3n(ji)=zg1*(1.-za4n(ji)*za1n(ji))-zg2*za5n(ji)
      zal1(ji)=zb2(ji,1)
      zal2(ji)=zb1(ji,1)
      zbe1(ji)=1.-zb2(ji,1)
      zbe2(ji)=1.-zb1(ji,1)
      zga1(ji)=1.-zb4(ji,1)
      zga2(ji)=1.-zb3(ji,1)
      zde1(ji)=zb4(ji,1)
      zde2(ji)=zb3(ji,1)
  620 continue

!     VI.3   PREMIERE COUCHE, ELIMINATION (SANS PROBLEMES).

      do 630 jn=1,iaucr
      do 630 ji=iideb(jn),iiut(jn)
      zfmc(ji,1)=za3c(ji)*(zal2(ji)*zfpc(ji,1))
      zfpc(ji,2)=za1c(ji)*(zal2(ji)*zfpc(ji,1))
      zfdc(ji,2)=za2c(ji)*(zal2(ji)*zfpc(ji,1))
      zfmn(ji,1)=za3n(ji)*(zbe2(ji)*zfpc(ji,1))
      zfpn(ji,2)=za1n(ji)*(zbe2(ji)*zfpc(ji,1))
      zfdn(ji,2)=za2n(ji)*(zbe2(ji)*zfpc(ji,1))
      ztu1(ji,1)=0.
      ztu2(ji,1)=zal1(ji)*za4c(ji)
      ztu3(ji,1)=zga1(ji)*za4c(ji)
      ztu4(ji,1)=zbe1(ji)*za4n(ji)
      ztu5(ji,1)=zde1(ji)*za4n(ji)
      ztu6(ji,1)=zal1(ji)*za5c(ji)
      ztu7(ji,1)=zga1(ji)*za5c(ji)
      ztu8(ji,1)=zbe1(ji)*za5n(ji)
      ztu9(ji,1)=zde1(ji)*za5n(ji)
  630 continue

!     VI.4   BOUCLE SUR LES NIVEAUX, CALCULS PRELIMINAIRES PUIS
!            ELIMINATION.

      do 643 jk=2,nk
      do 640 jn=1,iaucr
      do 640 ji=iideb(jn),iiut(jn)
      zeo1 = zueogs(ji,jk) + zeo1sa*daer(ji,jk) + (zeoray*dp(ji,jk))
      zeo2 =                 zeo2sa*daer(ji,jk) + (zeoray*dp(ji,jk))
      zeps=sqrt(zeo1*zeo1-zeo2*zeo2)
      ztau=exp(max(-zeps,zargli))
      zrho=zeo2/(zeo1+zeps)
      za4c(ji)=ztau*(1.-(zrho*zrho))*(1./(1.-(zrho*zrho)*(ztau*ztau)))
      za5c(ji)=zrho*(1.-(ztau*ztau))*(1./(1.-(zrho*zrho)*(ztau*ztau)))
      zeo=0.5*zdeogs(ji,jk)+zeosa*daer(ji,jk)+(zeoray*dp(ji,jk))
      zmu0ic=(zeps/zeo)+dsign(max(dabs(zmu0i(ji)-(zeps/zeo)),zepres),zmu0i(ji)-(zeps/zeo))
      zeo3=((zeodsa*daer(ji,jk))*zusa(ji)+(0.5*(zeoray*dp(ji,jk))))*zmu0ic
      zeo4=((zeodsa*daer(ji,jk))*(1.-zusa(ji))+(0.5*(zeoray*dp(ji,jk))))*zmu0ic
      zeo5=zeo*zmu0ic
      za1c(ji)=exp(max(-zeo5,zargli))
      zg1=(zeo3*(zeo5-zeo1)-zeo2*zeo4)*(1./(zeo5*zeo5-zeps*zeps))
      zg2=-(zeo4*(zeo5+zeo1)+zeo2*zeo3)*(1./(zeo5*zeo5-zeps*zeps))
      za2c(ji)=zg2*(za1c(ji)-za4c(ji))-zg1*za5c(ji)*za1c(ji)
      za3c(ji)=zg1*(1.-za4c(ji)*za1c(ji))-zg2*za5c(ji)
      zeo1=zeo1+zeo1sn*(dp(ji,jk)*zqli(ji,jk))+zeo1si*(dp(ji,jk)*zqice(ji,jk))
      zeo2=zeo2+zeo2sn*(dp(ji,jk)*zqli(ji,jk))+zeo2si*(dp(ji,jk)*zqice(ji,jk))
      zeps=sqrt(zeo1*zeo1-zeo2*zeo2)
      ztau=exp(max(-zeps,zargli))
      zrho=zeo2/(zeo1+zeps)
      za4n(ji)=ztau*(1.-(zrho*zrho))*(1./(1.-(zrho*zrho)*(ztau*ztau)))
      za5n(ji)=zrho*(1.-(ztau*ztau))*(1./(1.-(zrho*zrho)*(ztau*ztau)))
      zeo=zeo+zeosn*(dp(ji,jk)*zqli(ji,jk))+zeosi*(dp(ji,jk)*zqice(ji,jk))
      zmu0in=(zeps/zeo)+dsign(max(dabs(zmu0i(ji)-(zeps/zeo)),zepres),zmu0i(ji)-(zeps/zeo))
      zeo3=(zeo3/zmu0ic+zusn(ji)*((dp(ji,jk)*zqli(ji,jk))*zeodsn)       &
                       +zusi(ji)*((dp(ji,jk)*zqice(ji,jk))*zeodsi))*zmu0in
      zeo4=(zeo4/zmu0ic+(1.-zusn(ji))*((dp(ji,jk)*zqli(ji,jk))*zeodsn)  &
                       +(1.-zusi(ji))*((dp(ji,jk)*zqice(ji,jk))*zeodsi))*zmu0in
      zeo5=zeo*zmu0in
      za1n(ji)=exp(max(-zeo5,zargli))
      zg1=(zeo3*(zeo5-zeo1)-zeo2*zeo4)*(1./(zeo5*zeo5-zeps*zeps))
      zg2=-(zeo4*(zeo5+zeo1)+zeo2*zeo3)*(1./(zeo5*zeo5-zeps*zeps))
      za2n(ji)=zg2*(za1n(ji)-za4n(ji))-zg1*za5n(ji)*za1n(ji)
      za3n(ji)=zg1*(1.-za4n(ji)*za1n(ji))-zg2*za5n(ji)
      zal1(ji)=zb2(ji,jk)
      zal2(ji)=zb1(ji,jk)
      zbe1(ji)=1.-zb2(ji,jk)
      zbe2(ji)=1.-zb1(ji,jk)
      zga1(ji)=1.-zb4(ji,jk)
      zga2(ji)=1.-zb3(ji,jk)
      zde1(ji)=zb4(ji,jk)
      zde2(ji)=zb3(ji,jk)
  640 continue
      do 641 jn=1,iaucr
      do 641 ji=iideb(jn),iiut(jn)
      ztd1=1./(1.-za5c(ji)*(zal2(ji)*ztu6(ji,jk-1)+zga2(ji)*ztu8(ji,jk-1)))
      zfmc(ji,jk)=ztd1*(za5c(ji)*(zal2(ji)*zfdc(ji,jk)+zga2(ji)   &
                 *zfdn(ji,jk))+za3c(ji)*(zal2(ji)*zfpc(ji,jk)+zga2(ji)*zfpn(ji,jk)))
      ztu1(ji,jk)=ztd1*za5c(ji)*(zal2(ji)*ztu7(ji,jk-1)+zga2(ji)*ztu9(ji,jk-1))
      ztu2(ji,jk)=(ztd1*za4c(ji))*zal1(ji)
      ztu3(ji,jk)=(ztd1*za4c(ji))*zga1(ji)
      ztd2=za5n(ji)*(zbe2(ji)*ztu6(ji,jk-1)+zde2(ji)*ztu8(ji,jk-1))
      ztd3=1./(1.-za5n(ji)*(zbe2(ji)*ztu7(ji,jk-1)+zde2(ji)*ztu9(ji,jk-1))-ztd2*ztu1(ji,jk))
      zfmn(ji,jk)=ztd3*(za5n(ji)*(zbe2(ji)*zfdc(ji,jk)+zde2(ji)    &
                 *zfdn(ji,jk))+ztd2*zfmc(ji,jk)+za3n(ji)*(zbe2(ji) &
                 *zfpc(ji,jk)+zde2(ji)*zfpn(ji,jk)))
      ztu4(ji,jk)=ztd3*(za4n(ji)*zbe1(ji)+ztd2*ztu2(ji,jk))
      ztu5(ji,jk)=ztd3*(za4n(ji)*zde1(ji)+ztd2*ztu3(ji,jk))
  641 continue
      do 642 jn=1,iaucr
      do 642 ji=iideb(jn),iiut(jn)
      zfpc(ji,jk+1)=za1c(ji)*(zal2(ji)*zfpc(ji,jk)+zga2(ji)*zfpn(ji,jk))
      zfpn(ji,jk+1)=za1n(ji)*(zbe2(ji)*zfpc(ji,jk)+zde2(ji)*zfpn(ji,jk))
      ztd4=za4c(ji)*(zal2(ji)*ztu6(ji,jk-1)+zga2(ji)*ztu8(ji,jk-1))
      ztd5=za4c(ji)*(zal2(ji)*ztu7(ji,jk-1)+zga2(ji)*ztu9(ji,jk-1))
      zfdc(ji,jk+1)=za4c(ji)*(zal2(ji)*zfdc(ji,jk)+zga2(ji)*zfdn(ji,jk)) &
                   +ztd4*zfmc(ji,jk)+ztd5*zfmn(ji,jk)+za2c(ji)*(zal2(ji) &
                   *zfpc(ji,jk)+zga2(ji)*zfpn(ji,jk))
      ztu6(ji,jk)=za5c(ji)*zal1(ji)+ztd4*ztu2(ji,jk)+ztd5*ztu4(ji,jk)
      ztu7(ji,jk)=za5c(ji)*zga1(ji)+ztd4*ztu3(ji,jk)+ztd5*ztu5(ji,jk)
      ztd6=za4n(ji)*(zbe2(ji)*ztu6(ji,jk-1)+zde2(ji)*ztu8(ji,jk-1))
      ztd7=za4n(ji)*(zbe2(ji)*ztu7(ji,jk-1)+zde2(ji)*ztu9(ji,jk-1))
      zfdn(ji,jk+1)=za4n(ji)*(zbe2(ji)*zfdc(ji,jk)+zde2(ji)*zfdn(ji,jk)) &
                   +ztd6*zfmc(ji,jk)+ztd7*zfmn(ji,jk)+za2n(ji)*(zbe2(ji) &
                   *zfpc(ji,jk)+zde2(ji)*zfpn(ji,jk))
      ztu8(ji,jk)=za5n(ji)*zbe1(ji)+ztd6*ztu2(ji,jk)+ztd7*ztu4(ji,jk)
      ztu9(ji,jk)=za5n(ji)*zde1(ji)+ztd6*ztu3(ji,jk)+ztd7*ztu5(ji,jk)
  642 continue
  643 continue

!     VI.5   TRAITEMENT EN SURFACE.

      do 650 jn=1,iaucr
      do 650 ji=iideb(jn),iiut(jn)
      zal=alb(ji)                 ! diffuse
      zalp=zalpn(ji)              ! parallel
      ztds1=1./(1.-zal*ztu6(ji,nk))
      zfmc(ji,nk+1)=ztds1*(zal*zfdc(ji,nk+1)+zalp*zfpc(ji,nk+1))
      ztus1=ztds1*zal*ztu7(ji,nk)
      ztds2=zal*ztu8(ji,nk)
      ztds3=1./(1.-zal*ztu9(ji,nk)-ztds2*ztus1)
      zfmn(ji,nk+1)=ztds3*(zal*zfdn(ji,nk+1)+ztds2*zfmc(ji,nk+1)+zalp*zfpn(ji,nk+1))
      zfmc(ji,nk+1)=zfmc(ji,nk+1)+ztus1*zfmn(ji,nk+1)
  650 continue

!     VI.6   SUBSTITUTION NIVEAU PAR NIVEAU.

      do 661 jk=nk,1,-1
      do 660 jn=1,iaucr
      do 660 ji=iideb(jn),iiut(jn)
      zfdn(ji,jk+1)=zfdn(ji,jk+1)+ztu8(ji,jk)*zfmc(ji,jk+1)+ztu9(ji,jk)*zfmn(ji,jk+1)
      zfdc(ji,jk+1)=zfdc(ji,jk+1)+ztu6(ji,jk)*zfmc(ji,jk+1)+ztu7(ji,jk)*zfmn(ji,jk+1)
      zfmn(ji,jk)=zfmn(ji,jk)+ztu4(ji,jk)*zfmc(ji,jk+1)+ztu5(ji,jk)*zfmn(ji,jk+1)
      zfmc(ji,jk)=zfmc(ji,jk)+ztu2(ji,jk)*zfmc(ji,jk+1)+ztu3(ji,jk)*zfmn(ji,jk+1)+ztu1(ji,jk)*zfmn(ji,jk)
  660 continue
  661 continue

!     VI.7   CALCUL DES FLUX, DES DELTA-FLUX ET DES TAUX NETS
!            D'ECHAUFFEMENT RADIATIF.

      do 670 jn=1,iaucr
      do 670 ji=iideb(jn),iiut(jn)
      ztrans(ji)=zfpc(ji,nk+1)+zfdc(ji,nk+1)-zfmc(ji,nk+1)
      ztrans(ji)=ztrans(ji)+zfpn(ji,nk+1)+zfdn(ji,nk+1)-zfmn(ji,nk+1)
      zfnet(ji)=ztrans(ji)
      rg(ji)=zfnet(ji)
  670 continue
      do 673 jk=nk,1,-1
      do 671 jn=1,iaucr
      do 671 ji=iideb(jn),iiut(jn)
      ztranb(ji)=ztrans(ji)
  671 continue
      do 672 jn=1,iaucr
      do 672 ji=iideb(jn),iiut(jn)
      ztrans(ji)=zfpc(ji,jk)+zfdc(ji,jk)-zfmc(ji,jk)
      ztrans(ji)=ztrans(ji)+zfpn(ji,jk)+zfdn(ji,jk)-zfmn(ji,jk)
      dtfr(ji,jk)=dtfr(ji,jk)+(ztrans(ji)-ztranb(ji))*gg/(dp(ji,jk)*cph(ji,jk))
      zfnet(ji)=ztrans(ji)
  672 continue
  673 continue
      do 674 jn=1,iaucr
      do 674 ji=iideb(jn),iiut(jn)
      rvis(ji)=zfnet(ji)
  674 continue

      return
      end subroutine radial
!###############################################################################################################
      subroutine u_ghost (bufsend, ip_to, bufrecv, ip_from, nbuf)

      use mod_moloch, only: nbuffer
      implicit none

#ifdef mpi
      include 'mpif.h'
      integer :: status(mpi_status_size)
#else
      integer :: status
#endif

      integer :: nbuf, ip_to, ip_from, ierr, comm, tag=1, isend, irecv
      real    :: bufsend(nbuf), bufrecv(nbuf)

#ifdef mpi
      comm  = mpi_comm_world

      call mpi_isend (bufsend, nbuf, mpi_real, ip_to,   tag, comm, isend, ierr)
      call mpi_irecv (bufrecv, nbuf, mpi_real, ip_from, tag, comm, irecv, ierr)
      call mpi_wait  (isend, status, ierr)
      call mpi_wait  (irecv, status, ierr)
!      call mpi_send (bufsend, nbuf, mpi_real, ip_to  , tag, comm,         ierr)
!      call mpi_recv (bufrecv, nbuf, mpi_real, ip_from, tag, comm, status, ierr)

#endif

      return
      end subroutine u_ghost
!###############################################################################################################
    subroutine wrshf (jstep)

    use mod_moloch
    implicit none
    real(4) ztg, zzeta, ztbarv, zss, zesk, zee, zt0t, slp(nlon,nlat), zalf, ze1, ze2, zqsat, zriav, zwink
    integer :: iunit=26, iunit_work=29, jlon, jlat, jstep, j, jklev, &
 nlevinp, nlevout, jlonm1, jlatp1, jklev1, jklev2, jkpbl, jkup
    character(len=15) file_out
    integer, dimension(50) :: grib2_descript = 0
    real, dimension(200)   :: pdr_loc = 0.
    integer, parameter :: ivalmiss = -9999
    integer nmsg2
    real, dimension(nlev)  :: finp, xinp
    real, dimension(100)   :: fout, xout

    nmsg2 = 24
    if (jstep.eq.1) qprectot = 0.
    if (jstep.eq.1) qprecsolid = 0.

!  computation of SLP (Pascal)

    zss = 6.0e-03     ! average stability
    do jlat = 1, nlat
    do jlon = 1, nlon
    zzeta  = zeta(jlon,jlat,1) + phig(jlon,jlat)/g
    ztg    = tvirt(jlon,jlat,1) + zss*zzeta
    ztbarv = .5*(ztg + tvirt(jlon,jlat,1))
    slp(jlon,jlat) = p(jlon,jlat,1)*exp(g*zzeta/(ztbarv*rd))
    enddo
    enddo

    do j = 1, 8
    call filt2d (slp, 1.)
    enddo

!  completing external boundaries

!      if (ip_w.eq.ip_null) then
      do jlat = 1, nlat
      jklev=n_std_lev_sl(2,jlat)
      t_std_lev(1,jlat,1:jklev) = t_std_lev(2,jlat,1:jklev)
      u_std_lev(1,jlat,1:jklev) = u_std_lev(2,jlat,1:jklev)
      v_std_lev(1,jlat,1:jklev) = v_std_lev(2,jlat,1:jklev)
      q_std_lev(1,jlat,1:jklev) = q_std_lev(2,jlat,1:jklev)
      tg(1,jlat,1:nlevg) = tg(2,jlat,1:nlevg)
      qg(1,jlat,1:nlevg) = qg(2,jlat,1:nlevg)
      enddo
!      endif
!      if (ip_e.eq.ip_null) then
      do jlat = 1, nlat
      jklev=n_std_lev_sl(nlonm1,jlat)
      t_std_lev(nlon,jlat,1:jklev) = t_std_lev(nlonm1,jlat,1:jklev)
      u_std_lev(nlon,jlat,1:jklev) = u_std_lev(nlonm1,jlat,1:jklev)
      v_std_lev(nlon,jlat,1:jklev) = v_std_lev(nlonm1,jlat,1:jklev)
      q_std_lev(nlon,jlat,1:jklev) = q_std_lev(nlonm1,jlat,1:jklev)
      tg(nlon,jlat,1:nlevg) = tg(nlonm1,jlat,1:nlevg)
      qg(nlon,jlat,1:nlevg) = qg(nlonm1,jlat,1:nlevg)
      enddo
!      endif
!      if (ip_s.eq.ip_null) then
      do jlon = 1, nlon
      jklev=n_std_lev_sl(jlon,2)
      t_std_lev(jlon,1,1:jklev) = t_std_lev(jlon,2,1:jklev)
      u_std_lev(jlon,1,1:jklev) = u_std_lev(jlon,2,1:jklev)
      v_std_lev(jlon,1,1:jklev) = v_std_lev(jlon,2,1:jklev)
      q_std_lev(jlon,1,1:jklev) = q_std_lev(jlon,2,1:jklev)
      tg(jlon,1,1:nlevg) = tg(jlon,2,1:nlevg)
      qg(jlon,1,1:nlevg) = qg(jlon,2,1:nlevg)
      enddo
!      endif
!      if (ip_n.eq.ip_null) then
      do jlon = 1, nlon
      jklev=n_std_lev_sl(jlon,nlatm1)
      t_std_lev(jlon,nlat,1:jklev) = t_std_lev(jlon,nlatm1,1:jklev)
      u_std_lev(jlon,nlat,1:jklev) = u_std_lev(jlon,nlatm1,1:jklev)
      v_std_lev(jlon,nlat,1:jklev) = v_std_lev(jlon,nlatm1,1:jklev)
      q_std_lev(jlon,nlat,1:jklev) = q_std_lev(jlon,nlatm1,1:jklev)
      tg(jlon,nlat,1:nlevg) = tg(jlon,nlatm1,1:nlevg)
      qg(jlon,nlat,1:nlevg) = qg(jlon,nlatm1,1:nlevg)
      enddo
!      endif

! Lapse rate

      do jlat = 1, nlat
      do jlon = 1, nlon
      jklev2 = 5
       do jklev = nlev/2, 4, -1
       if (zeta(jlon,jlat,jklev) <= 1500.) exit
       enddo
      jklev2 = max(jklev, 5)
      lapse_rate(jlon,jlat) = (t(jlon,jlat,jklev2)-t(jlon,jlat,3))/(zeta(jlon,jlat,jklev2)-zeta(jlon,jlat,3))
      enddo
      enddo

!-----------------------------------------------------------------------------
!  Interpolation of variables at constant geometric heigh (m) over the surface
!-----------------------------------------------------------------------------

    zalf = 0.7
    ze1 = 0.5
    ze2 = 0.5

    nlevinp = nlev

    do jlat = 1, nlat
    do jlon = 1, nlon
      jlonm1=max(jlon-1,1)
      jlatp1=min(jlat+1,nlat)
      nlevout = n_std_lev_atm - n_std_lev_sl(jlon,jlat)

      if (nlevout >= 1 ) then

        xinp(1:nlev) = zeta(jlon,jlat,1:nlev)

        do jklev = 1, nlevout
        jklev1 = jklev+n_std_lev_sl(jlon,jlat)
        xout(jklev) = std_lev_atm(jklev1)
        enddo

        finp(1:nlev) = t(jlon,jlat,1:nlev)
        call interp (zalf, ze1, ze2, nlev, xinp, finp, xout(1:nlevout), fout(1:nlevout), nlevout)
        do jklev = 1, nlevout
        jklev1 = jklev+n_std_lev_sl(jlon,jlat)
        t_std_lev(jlon,jlat,jklev1)=fout(jklev)
        enddo

        finp(1:nlev) = 0.5*(u(jlon,jlat,1:nlev)+u(jlonm1,jlat,1:nlev)) ! u computed on T-points
        call interp (zalf, ze1, ze2, nlev, xinp, finp, xout(1:nlevout), fout(1:nlevout), nlevout)
        do jklev = 1, nlevout
        jklev1 = jklev+n_std_lev_sl(jlon,jlat)
        u_std_lev(jlon,jlat,jklev1)=fout(jklev)
        enddo

        finp(1:nlev) = 0.5*(v(jlon,jlat,1:nlev)+v(jlon,jlatp1,1:nlev)) ! v computed on T-points
        call interp (zalf, ze1, ze2, nlev, xinp, finp, xout(1:nlevout), fout(1:nlevout), nlevout)
        do jklev = 1, nlevout
        jklev1 = jklev+n_std_lev_sl(jlon,jlat)
        v_std_lev(jlon,jlat,jklev1)=fout(jklev)
        enddo

        finp(1:nlev) = q(jlon,jlat,1:nlev)
        call interp (zalf, ze1, ze2, nlev, xinp, finp, xout(1:nlevout), fout(1:nlevout), nlevout)
        do jklev = 1, nlevout
        jklev1 = jklev+n_std_lev_sl(jlon,jlat)
        q_std_lev(jlon,jlat,jklev1)=fout(jklev)
        enddo

! Dew point

        do jklev = 1,nlev
        zee = p(jlon,jlat,jklev)*q(jlon,jlat,jklev)/(eps+q(jlon,jlat,jklev)*(1.-eps))
        if (t(jlon,jlat,jklev) >= tzer) then
        finp(jklev) = (tzer-33.65/17.40*log(zee/611.))/(1.-1./17.40*log(zee/611.))
        else
        finp(jklev) = (tzer- 0.75/22.45*log(zee/611.))/(1.-1./22.45*log(zee/611.))
        endif
        enddo
        call interp (zalf, ze1, ze2, nlev, xinp, finp, xout(1:nlevout), fout(1:nlevout), nlevout)
        do jklev = 1, nlevout
        jklev1 = jklev+n_std_lev_sl(jlon,jlat)
        td_std_lev(jlon,jlat,jklev1)=fout(jklev)
        enddo

! Relative humidity

        do jklev = 1, nlev
        zee = p(jlon,jlat,jklev)*q(jlon,jlat,jklev)/(eps+q(jlon,jlat,jklev)*(1.-eps))
        zt0t = tzer/t(jlon,jlat,jklev)
        zesk = ezer*exp(-ccw1*log(zt0t)+ccw2*(1.-zt0t))  ! saturated pressure over water (Pa)
        finp(jklev) = zee/zesk*100.
        enddo
        call interp (zalf, ze1, ze2, nlev, xinp, finp, xout(1:nlevout), fout(1:nlevout), nlevout)
        do jklev = 1, nlevout
        jklev1 = jklev+n_std_lev_sl(jlon,jlat)
        rh_std_lev(jlon,jlat,jklev1)=fout(jklev)
        enddo

      endif

    enddo
    enddo

! Definition of relative humidity and dew point temperature in the surface layer

    do jlat = 1, nlat
    do jlon = 1, nlon
      do jklev = 1,n_std_lev_sl(jlon,jlat)
        zt0t = tzer/t_std_lev(jlon,jlat,jklev)
        zesk = ezer*exp(-ccw1*log(zt0t)+ccw2*(1.-zt0t))  ! saturated pressure over water (Pa)
        zqsat = ps(jlon,jlat)*q_std_lev(jlon,jlat,jklev)/(eps+q_std_lev(jlon,jlat,jklev)*(1.-eps))
        rh_std_lev(jlon,jlat,jklev) = min(zqsat/zesk*100., 100.)
        zee = ps(jlon,jlat)*q_std_lev(jlon,jlat,jklev)/(eps+q_std_lev(jlon,jlat,jklev)*(1.-eps))
        if (t_std_lev(jlon,jlat,jklev) >= tzer) then
          td_std_lev(jlon,jlat,jklev) = (tzer-33.65/17.40*log(zee/611.))/(1.-1./17.40*log(zee/611.))
        else
          td_std_lev(jlon,jlat,jklev) = (tzer- 0.75/22.45*log(zee/611.))/(1.-1./22.45*log(zee/611.))
        endif
      enddo
    enddo
    enddo

    if (n_std_lev_soil >= 1 ) then

      nlevinp = nlevg+1
      nlevout=n_std_lev_soil

      xinp(1) = 0.
      do jklev = 1,nlevg
        xinp(jklev+1) = lev_soil(jklev)
      enddo

      xout(1:nlevout) = std_lev_soil(1:nlevout)

      do jlat = 1, nlat
      do jlon = 1, nlon

        finp(1) = tskin(jlon,jlat)
        do jklev = 1,nlevg
          finp(jklev+1) = tg(jlon,jlat,jklev)
        enddo
        call interp (zalf, ze1, ze2, nlevinp, xinp, finp, xout(1:nlevout), fout(1:nlevout), nlevout)
        tg_std_lev(jlon,jlat,1:nlevout)=fout(1:nlevout)

        finp(1) = qg(jlon,jlat,1)
        do jklev = 1,nlevg
          finp(jklev+1) = qg(jlon,jlat,jklev)
        enddo
        call interp (zalf, ze1, ze2, nlevinp, xinp, finp, xout(1:nlevout), fout(1:nlevout), nlevout)
        qg_std_lev(jlon,jlat,1:nlevout)=fout(1:nlevout)

      enddo
      enddo

    endif

!------------------------------------------------------------------------------------------
!   output section
!------------------------------------------------------------------------------------------

  if(myid == 0) then

    ishf = ishf+1

#ifdef oper
    write (file_out,'(a,i3.3,a)') 'moloch_',ishf,'.shf'
    open (iunit, file=file_out, form='unformatted', status='unknown')
#else
    open (iunit, file='moloch.shf', form='unformatted', status='unknown', position='append')
#endif

    nfdr(19) = nmsg2
    pdr_loc(1:100) = pdr(1:100)
    write(iunit) nfdr
    write(iunit) pdr_loc(1:200)

  endif

! grib2_description - descripting record for grib2 format coding

!  1 - model index (1 Bolam, 2 Moloch, 3 Globo)
!  2 - grid template index (1 horizontal grid, 1000 vertical cross-section)
!  3 - product template index (0 instant, 8 statistical, 32 forecast satellite, 1 instant individual ensemble, &
!      11 statistical individual ensemble)
!  4 - time unit index (0 minute, 1 hour, 2 day, 3 month, 4 year, 13 second)
!  5 - statistical elaboration type (for statistical products only): 0 average, 1 accumulation, 2 maximum, 3 minimum
!  6 - production status of data (0 Operational products, 2 Research products, 3 Re-analysis products, &
!      7 Sub-seasonal to seasonal prediction S2S)
!  7 - type of data (0 Analysis, 1 Forecast, 2 Analysis and forecast, 3 Control forecast, 4 Perturbed forecast)
!  8 - indicator of time unit for the increment between successive fields used for statistical elaboration
! 10 - level (layer) type (for forecast satellite products, code of satellite platform and sensor)
! 11 - first scaled value of level (layer)
! 12 - scale of first value of level (layer)
! 13 - second scaled value of level (layer)
! 14 - scale of second value of level (layer)
! 15 - level type of second level for a layer
! 20 - product discipline
! 21 - product category
! 22 - product parameter
! 30 - flag of bit-map presence: 0 not present, 1 present (use VALMISS)
! 31 - member number of perturbed forecast
! 32 - member index of perturbed forecast

    grib2_descript( 1)    = 2 ! Moloch model
    grib2_descript( 2)    = 1 ! horizontal grid
    grib2_descript( 3)    = 0 ! instant product
    grib2_descript( 4)    = 0 ! time unit is minute
    grib2_descript( 6)    = 0 ! operational products
    grib2_descript( 7)    = 1 ! forecast
    grib2_descript(10:15) = ivalmiss
    grib2_descript(30)    = 0 ! bit-map absent

    call collect (slp, gfield) ! #1
    if(myid == 0) then
      grib2_descript(10) = 101 ! Mean sea level
      grib2_descript(20) =   0 ! Discipline: Meteorological products
      grib2_descript(21) =   3 ! Category: Mass
      grib2_descript(22) =   1 ! Parameter: Pressure reduced to MSL (Pa)
      write (iunit) grib2_descript
      call wrec2 (iunit, gnlon, gnlat, gfield)
    endif

    grib2_descript(10) =   1 ! Ground or water surface
    grib2_descript(11:12) = 0
 
    call collect (qprectot, gfield) ! #2
    if(myid == 0) then
      grib2_descript(3)  =   8 ! statistical product
      grib2_descript(5)  =   1 ! accumulation
      grib2_descript(20) =   0 ! Discipline: Meteorological products
      grib2_descript(21) =   1 ! Category: Moisture
      grib2_descript(22) =   8 ! Parameter: Total precipitation (kg m-2)
      write (iunit) grib2_descript
      call wrec2 (iunit, gnlon, gnlat, gfield)
    endif

    call collect (cloudt, gfield) ! #3
    if(myid == 0) then
      grib2_descript(3)  =   0 ! instant product
      grib2_descript(5)  =   0
      grib2_descript(20) =   0 ! Discipline: Meteorological products
      grib2_descript(21) =   6 ! Category: Cloud
      grib2_descript(22) =   1 ! Parameter: Total cloud cover (%)
      write (iunit) grib2_descript
      call wrec2 (iunit, gnlon, gnlat, gfield*1.e2)
    endif

    if (jstep.gt.1) totsky = totsky/float(ntswshf)
    call collect (totsky, gfield) ! #4
    if(myid == 0) then
      grib2_descript(3)  =   8 ! statistical product
      grib2_descript(5)  =   0 ! average
      grib2_descript(20) =   0 ! Discipline: Meteorological products
      grib2_descript(21) =   4 ! Category: Short-wave Radiation
      grib2_descript(22) =   7 ! Parameter: Downward short-wave radiation flux (W m-2)
      write (iunit) grib2_descript
      call wrec2 (iunit, gnlon, gnlat, gfield)
    endif

!    if (jstep.gt.1) soldir = soldir/float(ntswshf)
!    call collect (soldir, gfield)
!    if(myid == 0) then
!      grib2_descript(3)  =   8 ! statistical product
!      grib2_descript(5)  =   0 ! average
!      grib2_descript(20) =   0 ! Discipline: Meteorological products
!      grib2_descript(21) =   4 ! Category: Short-wave Radiation
!      grib2_descript(22) =  13 ! Parameter: Direct solar radiation flux (W m-2)
!      write (iunit) grib2_descript
!      call wrec2 (iunit, gnlon, gnlat, gfield)
!    endif

    call collect (qprecsolid, gfield) ! #5
    if(myid == 0) then
      grib2_descript(3)  =   8 ! statistical product
      grib2_descript(5)  =   1 ! accumulation
      grib2_descript(20) =   0 ! Discipline: Meteorological products
      grib2_descript(21) =   1 ! Category: Moisture
      grib2_descript(22) =  29 ! Parameter: Total snowfall (m)
      write (iunit) grib2_descript
      call wrec2 (iunit, gnlon, gnlat, gfield*1.e-3)
    endif

    call collect (snow, gfield) ! #6
    if(myid == 0) then
      grib2_descript(3)  =   0 ! instant product
      grib2_descript(5)  =   0
      grib2_descript(20) =   0 ! Discipline: Meteorological products
      grib2_descript(21) =   1 ! Category: Moisture
      grib2_descript(22) =  13 ! Parameter: Water equivalent of accumulated snow depth (kg m-2)
      write (iunit) grib2_descript
      call wrec2 (iunit, gnlon, gnlat, gfield*1.e3)
    endif

    call collect (lapse_rate, gfield) ! #7
    if(myid == 0) then
      grib2_descript(20) =   0 ! Discipline: Meteorological products
      grib2_descript(21) =   0 ! Category: Temperature
      grib2_descript(22) =   8 ! Parameter: Lapse rate (K m-1)
      write (iunit) grib2_descript
      call wrec2 (iunit, gnlon, gnlat, gfield)
    endif

    call collect (tg(1,1,1), gfield) ! #8
    if(myid == 0) then
      grib2_descript(10) = 106  ! Depth below land surface (m)
      grib2_descript(11) = nint(lev_soil(1)*1.e3) ! First scaled value of level (layer)
      grib2_descript(12) =   3  ! Scale of first value of level (layer)
      grib2_descript(20) =   2  ! Discipline: Land surface products
      grib2_descript(21) =   0  ! Category: Vegetation/Biomass
      grib2_descript(22) =   2  ! Parameter: Soil temperature (K)
      write (iunit) grib2_descript
      call wrec2 (iunit, gnlon, gnlat, gfield)
    endif

    grib2_descript(10) =   1 ! Ground or water surface
    grib2_descript(11:12) = 0

    call collect (tskin, gfield) ! #9
    if(myid == 0) then
      grib2_descript(20) =   0 ! Discipline: Meteorological products
      grib2_descript(21) =   0 ! Category: Temperature
      grib2_descript(22) =   0 ! Parameter: Temperature (K)
      write (iunit) grib2_descript
      call wrec2 (iunit, gnlon, gnlat, gfield)
    endif

    if (jstep.gt.1) then
      call collect (shf_accum/float(ntswshf), gfield)
    else
      call collect (shf_accum, gfield) ! #10
    endif
    if(myid == 0) then
      grib2_descript(3)  =   8 ! statistical product
      grib2_descript(5)  =   0 ! average
      grib2_descript(20) =   0 ! Discipline: Meteorological products
      grib2_descript(21) =   0 ! Category: Temperature
      grib2_descript(22) =  11 ! Parameter: Sensible heat net flux (W m-2)
      write (iunit) grib2_descript
      call wrec2 (iunit, gnlon, gnlat, gfield)
      grib2_descript(3)  =   0 ! instant product
      grib2_descript(5)  =   0
    endif

    if (jstep.gt.1) then
      call collect (lhf_accum/float(ntswshf), gfield)
    else
      call collect (lhf_accum, gfield) ! #11
    endif
    if(myid == 0) then
      grib2_descript(3)  =   8 ! statistical product
      grib2_descript(5)  =   0 ! average
      grib2_descript(20) =   0 ! Discipline: Meteorological products
      grib2_descript(21) =   0 ! Category: Temperature
      grib2_descript(22) =  10 ! Parameter: Latent heat net flux (W m-2)
      write (iunit) grib2_descript
      call wrec2 (iunit, gnlon, gnlat, gfield)
      grib2_descript(3)  =   0 ! instant product
      grib2_descript(5)  =   0
    endif

    if (jstep.gt.1) then
      call collect (qf_accum/float(ntswshf), gfield)
    else
      call collect (qf_accum, gfield) ! #12
    endif
    if(myid == 0) then
      grib2_descript(3)  =   8 ! statistical product
      grib2_descript(5)  =   0 ! average
      grib2_descript(20) =   0 ! Discipline: Meteorological products
      grib2_descript(21) =   1 ! Category: Moisture
      grib2_descript(22) =   6 ! Parameter: Evaporation (kg m-2)
      write (iunit) grib2_descript
      call wrec2 (iunit, gnlon, gnlat, gfield)
      grib2_descript(3)  =   0 ! instant product
      grib2_descript(5)  =   0
    endif

    grib2_descript(3)  =   0 ! instant product
    grib2_descript(5)  =   0
    grib2_descript(10) = 103 ! Specified height level above ground (m)
    grib2_descript(12) = 0

    do jklev = 1, n_std_lev_atm
      if (int(std_lev_atm(jklev)) == 2) then            ! 2m temperature only
        call collect (t_std_lev(1,1,jklev), gfield) ! #13
        if(myid == 0) then
          grib2_descript(11) =   int(std_lev_atm(jklev))
          grib2_descript(20) =   0 ! Discipline: Meteorological products
          grib2_descript(21) =   0 ! Category: Temperature
          grib2_descript(22) =   0 ! Parameter: Temperature (K)
          write (iunit) grib2_descript
          call wrec2 (iunit, gnlon, gnlat, gfield)
        endif
      endif

      if (int(std_lev_atm(jklev)) >= 10) then
        call collect (u_std_lev(1,1,jklev), gfield) ! #14-15-16-17
        if(myid == 0) then
          grib2_descript(11) =   int(std_lev_atm(jklev))
          grib2_descript(20) =   0 ! Discipline: Meteorological products
          grib2_descript(21) =   2 ! Category: Momentum
          grib2_descript(22) =   2 ! Parameter: u-component of wind (m s-1)
          write (iunit) grib2_descript
          call wrec2 (iunit, gnlon, gnlat, gfield)
        endif
      endif

      if (int(std_lev_atm(jklev)) >= 10) then
        call collect (v_std_lev(1,1,jklev), gfield) ! #18-19-20-21
        if(myid == 0) then
          grib2_descript(11) =   int(std_lev_atm(jklev))
          grib2_descript(20) =   0 ! Discipline: Meteorological products
          grib2_descript(21) =   2 ! Category: Momentum
          grib2_descript(22) =   3 ! Parameter: v-component of wind (m s-1)
          write (iunit) grib2_descript
          call wrec2 (iunit, gnlon, gnlat, gfield)
        endif
      endif

      if (int(std_lev_atm(jklev)) == 2) then
        call collect (rh_std_lev(1,1,jklev), gfield) ! #22
        if(myid == 0) then
          grib2_descript(11) =   int(std_lev_atm(jklev))
          grib2_descript(20) =   0 ! Discipline: Meteorological products
          grib2_descript(21) =   1 ! Category: Moisture
          grib2_descript(22) =   1 ! Parameter: Relative humidity (%)
          write (iunit) grib2_descript
          call wrec2 (iunit, gnlon, gnlat, gfield)
        endif
      endif

      if (int(std_lev_atm(jklev)) == 2) then
        call collect (q_std_lev(1,1,jklev), gfield) ! #23
        if(myid == 0) then
          grib2_descript(11) =   int(std_lev_atm(jklev))
          grib2_descript(20) =   0 ! Discipline: Meteorological products
          grib2_descript(21) =   1 ! Category: Moisture
          grib2_descript(22) =   0 ! Parameter: Specific humidity (kg kg-1)
          write (iunit) grib2_descript
          call wrec2 (iunit, gnlon, gnlat, gfield)
        endif
      endif

      if (int(std_lev_atm(jklev)) == 2) then
        call collect (td_std_lev(1,1,jklev), gfield) ! #24
        if(myid == 0) then
          grib2_descript(11) =   int(std_lev_atm(jklev))
          grib2_descript(20) =   0 ! Discipline: Meteorological products
          grib2_descript(21) =   0 ! Category: Temperature
          grib2_descript(22) =   6 ! Parameter: Dew point temperature (K)
          write (iunit) grib2_descript
          call wrec2 (iunit, gnlon, gnlat, gfield)
        endif
      endif

    enddo

    grib2_descript(10) = 106 ! Depth below land surface (m)
    grib2_descript(11) =   0
    grib2_descript(12) =   3

    do jklev = 1, n_std_lev_soil

      call collect (tg_std_lev(1,1,jklev), gfield)
      if(myid == 0) then
        grib2_descript(11) =   int(std_lev_soil(jklev)*1.e3)
        grib2_descript(20) =   2 ! Discipline:  Land surface products
        grib2_descript(21) =   0 ! Category: Vegetation/Biomass
        grib2_descript(22) =   2 ! Parameter: Soil temperature (K)
        write (iunit) grib2_descript
        call wrec2 (iunit, gnlon, gnlat, gfield)
      endif

      call collect (qg_std_lev(1,1,jklev), gfield)
      if(myid == 0) then
        grib2_descript(11) =   int(std_lev_soil(jklev)*1.e3)
        grib2_descript(20) =   2 ! Discipline:  Land surface products
        grib2_descript(21) =   0 ! Category: Vegetation/Biomass
        grib2_descript(22) =   9 ! Parameter: Volumetric soil moisture content (Proportion)
        write (iunit) grib2_descript
        call wrec2 (iunit, gnlon, gnlat, gfield)
      endif

    enddo

!    call collect (gust, gfield) ! #25
!    if(myid == 0) then
!      grib2_descript(10) =   1 ! Ground or water surface
!      grib2_descript(11:12) = 0
!      grib2_descript(20) =   0 ! Discipline: Meteorological products
!      grib2_descript(21) =   2 ! Category: Momentum
!      grib2_descript(22) =  22 ! Parameter: Wind speed (gust)  (M sec-1)
!      write (iunit) grib2_descript
!      call wrec2 (iunit, gnlon, gnlat, gfield)
!    endif

  if(myid == 0) then
    flush (iunit)
    close (iunit)

#ifdef oper
    call system("sync")
!!!      call system("ls -l -L "//file_out)
!!!      call system("date")
    open (iunit_work, file=file_out(1:14)//'.txt', status='unknown')
    write (iunit_work,'(2A)') file_out,' is full and closed'
    close (iunit_work)
    call system("sync")
    print *,'Output written on file ', file_out
#else
    print *,'Output written on file moloch.shf'
#endif

  endif

    qprectot   = 0. ! reset total precipitation
    qprecsolid = 0. ! reset snow fall
    totsky     = 0. ! reset cumulated flux
    soldir     = 0. ! reset cumulated flux
    shf_accum  = 0. ! reset cumulated sensible heat flux
    lhf_accum  = 0. ! reset cumulated latent heat flux
    qf_accum   = 0. ! reset cumulated evaporation flux

    return
    end subroutine wrshf
!###############################################################################################################
subroutine wrshf_radar(jstep)
! Writting of simulated radar reflectivity

use mod_moloch
implicit none

integer, parameter :: nlevz   = 40     ! no. of constant height levels for radar interpolation
real, dimension(nlon,nlat,nlevz) :: radarz
real, dimension(nlev)  :: zinp
real, dimension(nlevz) :: zout

character(len=15) :: file_out, amhf
integer, dimension(50) :: grib2_descript = 0
real, dimension(200)   :: pdr_loc = 0.
integer, parameter :: ivalmiss = -9999
integer :: iunit=41, iunit_work=29, jstep, jlon, jlat, jklev


! Vertical interpolation of radar reflectivity from zita to zeta coordinates

 do jklev = 1,nlevz
   zout(jklev) = float(jklev)/float(nlevz)*11000. ! 11 km top
 enddo

 do jlat = 1, nlat
 do jlon = 1, nlon
   do jklev = 1, nlev
     zinp(jklev) = zeta(jlon,jlat,jklev) + phig(jlon,jlat)/g ! height above sea level
   enddo
   call interp (.5,.1,.1,nlev,zinp,rradar(jlon,jlat,1:nlev),zout(1:nlevz),radarz(jlon,jlat,1:nlevz),nlevz)
   do jklev = 1, nlevz
     if (zout(jklev).lt.phig(jlon,jlat)/g) radarz(jlon,jlat,jklev) = radarmval
   enddo
 enddo
 enddo

 if (myid == 0) then

#ifdef oper
   write (amhf,'(i3.3)') imhf
   file_out='moloch_radar_'//trim(amhf)//'.shf'
   open (iunit, file=trim(file_out), form='unformatted', status='unknown')
#else
   open (iunit, file='moloch_radar.shf', form='unformatted', status='unknown', position='append')
#endif

   pdr_loc(1:100) = pdr(1:100)
   write (iunit) nfdr
   write (iunit) pdr_loc(1:200)

 endif

! grib2_description - descripting record for grib2 format coding

!  1 - model index (1 Bolam, 2 Moloch, 3 Globo)
!  2 - grid template index (1 horizontal grid, 1000 vertical cross-section)
!  3 - product template index (0 instant, 8 statistical, 32 forecast satellite, 1 instant individual ensemble, &
!      11 statistical individual ensemble)
!  4 - time unit index (0 minute, 1 hour, 2 day, 3 month, 4 year, 13 second)
!  5 - statistical elaboration type (for statistical products only): 0 average, 1 accumulation, 2 maximum, 3 minimum
!  6 - production status of data (0 Operational products, 2 Research products, 3 Re-analysis products, &
!      7 Sub-seasonal to seasonal prediction S2S)
!  7 - type of data (0 Analysis, 1 Forecast, 2 Analysis and forecast, 3 Control forecast, 4 Perturbed forecast)
!  8 - indicator of time unit for the increment between successive fields used for statistical elaboration
! 10 - level (layer) type (for forecast satellite products, code of satellite platform and sensor)
! 11 - first scaled value of level (layer)
! 12 - scale of first value of level (layer)
! 13 - second scaled value of level (layer)
! 14 - scale of second value of level (layer)
! 15 - level type of second level for a layer
! 20 - product discipline
! 21 - product category
! 22 - product parameter
! 30 - flag of bit-map presence: 0 not present, 1 present (use VALMISS)
! 31 - member number of perturbed forecast
! 32 - member index of perturbed forecast

 grib2_descript( 1)    = 2 ! Moloch model
 grib2_descript( 2)    = 1 ! horizontal grid
 grib2_descript( 3)    = 0 ! instant product
 grib2_descript( 4)    = 0 ! time unit is minute
 grib2_descript( 6)    = 0 ! operational products
 grib2_descript( 7)    = 1 ! forecast
 grib2_descript(10:15) = ivalmiss
 grib2_descript(30)    = 0 ! bit-map absent

 grib2_descript(3)  =   0 ! instant product
 grib2_descript(5)  =   0
 grib2_descript(iunit) = 103 ! Specified height level above ground (m)
 grib2_descript(12) = 0 
 grib2_descript(20) =   0 ! Discipline: Meteorological products
 grib2_descript(21) =  15 ! Category: Radar
 grib2_descript(22) =   1 ! Parameter: Base reflectivity (dB)

 do jklev = 1, nlevz

   call collect (radarz(1,1,jklev), gfield)
   if (myid == 0) then
     grib2_descript(11) = nint(zout(jklev))
     write (iunit) grib2_descript
     call wrec2 (iunit, gnlon, gnlat, gfield)
   endif

 enddo

 if (myid == 0) then
   flush(iunit)
   close(iunit)
#ifdef oper
   call system("sync")
!!!   call system("ls -l -L "//file_out)
!!!   call system("date")
   open (iunit_work, file=(trim(file_out)//'.txt'), status='unknown')
   write (iunit_work,'(2A)') trim(file_out),' is full and closed'
   close (iunit_work)
   call system("sync")
   print *,'Output written on file ', trim(file_out)
#else
   print *,'Output written on file moloch_radar.shf'
#endif
 endif

end subroutine wrshf_radar
!###############################################################################################################
    subroutine advect_waf (dt)

!  3D advection of all variables with WAF scheme

    use mod_moloch, only : nlon, nlat, nlev, nlevp1, nlonm1, nlatm1, u, v, w, tetav, pai, q, qcw, qci, qpw, tang, &
                           qpi1, qpi2, tke, ut, vt, wt, tket, ncw, nci, ip_e, ip_n, ip_s, ip_w, ip_null, nlmic2
    implicit none
    real, dimension(nlon,nlev) :: zvt
    real, dimension(nlat,nlev) :: zut
    real dt
    integer jlon, jlat, jklev

!----------------------
!  U-wind on T points
!----------------------

#ifdef mpi
      call u_ghost (u(nlonm1,:,:), ip_e, u(1,   :,:), ip_w, nlat*nlev)
      call u_ghost (u(2,     :,:), ip_w, u(nlon,:,:), ip_e, nlat*nlev)
#endif

    do jlon = 3, nlonm1
    ut(jlon,:,:) = .5625*(u(jlon,:,:)+u(jlon-1,:,:)) -.0625*(u(jlon+1,:,:)+u(jlon-2,:,:))
    enddo

#ifdef mpi
      call u_ghost (u(nlon-2,:,:), ip_e, zut, ip_w, nlat*nlev)
      if (ip_w.eq.ip_null) then
        ut(1,:,:) =  u(1,:,:)
        ut(2,:,:) = .5*(u(1,:,:)+u(2,:,:))
      else
        ut(2,:,:) = .5625*(u(2,:,:)+u(1,:,:)) - .0625*(u(3,:,:)+zut(:,:))
      endif
      if (ip_e.eq.ip_null) then
        ut(nlon,:,:) = .5*(u(nlon,:,:)+u(nlonm1,:,:))
      endif
#else
      ut(1,:,:) =  u(1,:,:)
      ut(2,:,:) = .5*(u(1,:,:)+u(2,:,:))
      ut(nlon,:,:) = .5*(u(nlon,:,:)+u(nlonm1,:,:))
#endif

!----------------------
!  V-wind on T points
!----------------------

#ifdef mpi
      call u_ghost (v(:,nlatm1,:), ip_n, v(:,1   ,:), ip_s, nlon*nlev)
      call u_ghost (v(:,2     ,:), ip_s, v(:,nlat,:), ip_n, nlon*nlev)
#endif

    do jlat = 2, nlat-2
    vt(:,jlat,:) = .5625*(v(:,jlat+1,:)+v(:,jlat,:))-.0625*(v(:,jlat+2,:)+v(:,jlat-1,:))
    enddo

#ifdef mpi
      call u_ghost (v(:,3,:), ip_s, zvt, ip_n, nlon*nlev)
      if (ip_n.eq.ip_null) then
        vt(:,nlat  ,:) = v(:,nlat,:)
        vt(:,nlatm1,:) = .5*(v(:,nlat,:)+v(:,nlatm1,:))
      else
        vt(:,nlatm1,:) = .5625*(v(:,nlat,:)+v(:,nlatm1,:))-.0625*(zvt(:,:)+v(:,nlat-2,:))
      endif
      if (ip_s.eq.ip_null) then
        vt(:,1,:) = .5*(v(:,1,:)+v(:,2,:))
      endif
#else
      vt(:,nlat  ,:) = v(:,nlat,:)
      vt(:,nlatm1,:) = .5*(v(:,nlat,:)+v(:,nlatm1,:))
      vt(:,1,:) = .5*(v(:,1,:)+v(:,2,:))
#endif

!----------------------
!  TKE and W at levels
!----------------------

    do jklev = 2, nlev-1
    do jlat = 2, nlatm1
    do jlon = 2, nlonm1
    wt(jlon,jlat,jklev)   = .5625*(w  (jlon,jlat,jklev+1)+w  (jlon,jlat,jklev))   &
                           -.0625*(w  (jlon,jlat,jklev+2)+w  (jlon,jlat,jklev-1))
    tket(jlon,jlat,jklev) = .5625*(tke(jlon,jlat,jklev+1)+tke(jlon,jlat,jklev))   &
                           -.0625*(tke(jlon,jlat,jklev+2)+tke(jlon,jlat,jklev-1))
    enddo
    enddo
    enddo
    do jlat = 2, nlatm1
    do jlon = 2, nlonm1
    wt  (jlon,jlat,1   )=.5*(w  (jlon,jlat,2     )+w  (jlon,jlat,1   ))
    wt  (jlon,jlat,nlev)=.5*(w  (jlon,jlat,nlevp1)+w  (jlon,jlat,nlev))
    tket(jlon,jlat,1   )=.5*(tke(jlon,jlat,2     )+tke(jlon,jlat,1   ))
    tket(jlon,jlat,nlev)=.5*(tke(jlon,jlat,nlevp1)+tke(jlon,jlat,nlev))
    enddo
    enddo

!----------------------
!  advections
!----------------------

    call wafone (tetav, dt)
    call wafone (pai, dt)
    call wafone (q, dt)
    call wafone (ut, dt)
    call wafone (vt, dt)
    call wafone (wt, dt)
    call wafone (qcw, dt)
    call wafone (qci, dt)
    if (nlmic2) call wafone (ncw, dt)     !2m
    if (nlmic2) call wafone (nci, dt)     !2m
    call wafone (qpw, dt)
    call wafone (qpi1, dt)
    call wafone (qpi2, dt)
    call wafone (tket, dt)

!---------------------------------------
!  curvature terms
!---------------------------------------

    do jklev = 1, nlev
    do jlat = 2, nlatm1
    do jlon = 2, nlonm1
    ut(jlon,jlat,jklev) = ut(jlon,jlat,jklev) + ut(jlon,jlat,jklev)*vt(jlon,jlat,jklev)*tang(jlat)*dt
    vt(jlon,jlat,jklev) = vt(jlon,jlat,jklev) - ut(jlon,jlat,jklev)*ut(jlon,jlat,jklev)*tang(jlat)*dt
    enddo
    enddo
    enddo

!---------------------------------------
!  back to wind points: U (fourth order)
!---------------------------------------

#ifdef mpi
      call u_ghost (ut(nlonm1,:,:), ip_e, ut(1   ,:,:), ip_w, nlat*nlev)
      call u_ghost (ut(2     ,:,:), ip_w, ut(nlon,:,:), ip_e, nlat*nlev)
#endif

    do jlat = 2, nlatm1
    do jlon = 2, nlon-2
    u(jlon,jlat,:) = .5625*(ut(jlon,jlat,:)+ut(jlon+1,jlat,:))-.0625*(ut(jlon-1,jlat,:)+ut(jlon+2,jlat,:))
    enddo
    enddo

#ifdef mpi
      call u_ghost (ut(3,:,:), ip_w, zut, ip_e, nlat*nlev)
      if (ip_e.eq.ip_null) then
        u(nlonm1,2:nlatm1,:) = .5*(ut(nlon,2:nlatm1,:)+ut(nlonm1,2:nlatm1,:))
      else
        u(nlonm1,2:nlatm1,:) = .5625*(ut(nlonm1,2:nlatm1,:)+ut (nlon,2:nlatm1,:)) &
                              -.0625*(ut(nlon-2,2:nlatm1,:)+zut(     2:nlatm1,:))
      endif
#else
      u(nlonm1,2:nlatm1,:) = .5*(ut(nlon,2:nlatm1,:)+ut(nlonm1,2:nlatm1,:))
#endif

!---------------------------------------
!  back to wind points: V (fourth order)
!---------------------------------------

#ifdef mpi
      call u_ghost (vt(2:nlonm1,nlatm1,:), ip_n, vt(2:nlonm1,1   ,:), ip_s, (nlon-2)*nlev)
      call u_ghost (vt(2:nlonm1,2     ,:), ip_s, vt(2:nlonm1,nlat,:), ip_n, (nlon-2)*nlev)
#endif

    do jlat = 3, nlatm1
    v(2:nlonm1,jlat,1:nlev) = .5625*(vt(2:nlonm1,jlat  ,1:nlev)+vt(2:nlonm1,jlat-1,1:nlev)) &
                             -.0625*(vt(2:nlonm1,jlat+1,1:nlev)+vt(2:nlonm1,jlat-2,1:nlev))
    enddo

#ifdef mpi
      call u_ghost (vt(2:nlonm1,nlat-2,:), ip_n, zvt(2:nlonm1,:), ip_s, (nlon-2)*nlev)
      if (ip_s.eq.ip_null) then
        v(2:nlonm1,2,1:nlev) = .5*(vt(2:nlonm1,1,1:nlev)+vt(2:nlonm1,2,1:nlev))
      else
        v(2:nlonm1,2,1:nlev) = .5625*(vt(2:nlonm1,2,1:nlev)+vt (2:nlonm1,1,1:nlev)) &
                              -.0625*(vt(2:nlonm1,3,1:nlev)+zvt(2:nlonm1,  1:nlev))
      endif
#else
      v(2:nlonm1,2,1:nlev) = .5*(vt(2:nlonm1,1,1:nlev)+vt(2:nlonm1,2,1:nlev))
#endif

!----------------------
!  back to half-levels
!----------------------

    do jklev = 3, nlev-1
    do jlat = 2, nlatm1
    do jlon = 2, nlonm1
    w  (jlon,jlat,jklev) = .5625*(wt  (jlon,jlat,jklev  )+wt  (jlon,jlat,jklev-1)) &
                          -.0625*(wt  (jlon,jlat,jklev+1)+wt  (jlon,jlat,jklev-2))
    tke(jlon,jlat,jklev) = .5625*(tket(jlon,jlat,jklev  )+tket(jlon,jlat,jklev-1)) &
                          -.0625*(tket(jlon,jlat,jklev+1)+tket(jlon,jlat,jklev-2))
    enddo
    enddo
    enddo
    do jlat = 2, nlatm1
    do jlon = 2, nlonm1
    w(jlon,jlat,2   )=.5*(wt(jlon,jlat,2   )+wt(jlon,jlat,1     ))
    w(jlon,jlat,nlev)=.5*(wt(jlon,jlat,nlev)+wt(jlon,jlat,nlev-1))
    tke(jlon,jlat,2)=.5*(tket(jlon,jlat,2)+tket(jlon,jlat,1))
    enddo
    enddo

    return
    end subroutine advect_waf
!###############################################################################################################
    subroutine micro2m (jstep)

! two-moment microphysics

    use mod_moloch, only: nlon, nlat, nlev, ntop, q, qcw, qci, qpw, qpi1, qpi2, t, p, ncw, nci, w, tke,  &
                          ccw1, ccw2, cci1, cci2, cpd, cpv, rd, rv, eps, cw, ci, yliv, yliw, ylwv, tzer, &
                          pzer, ezer, pi, dtstep, fmask, fmz, dz, raini, snowi, nlmic2, rradar, tvirt,   &
                          nlradar, nhist, ps, radarmval
    implicit none

! microphysical parameters

    real, parameter :: alfacw=6., alfaci=3.                    ! distribution parameters of clouds
    real, parameter :: yncw0=3.1e7, ynci0=5.7e7                ! numbers of cloud particles (continental)
    real, parameter :: yn0r=5.e6, yn0s=1.0e7, yn0h=6.2e6       ! distribution parameters of hydrometeors
    real, parameter :: acw=pi/6.*1000., bcw=3.                 ! cloud water mass parameters
    real, parameter :: aci=100., bci=2.55                      ! cloud ice mass parameters
    real, parameter :: ar =pi/6.*1000., br =3.                 ! rain mass parameters
    real, parameter :: as =30., bs =2.7                        ! snow mass parameters
    real, parameter :: ah =ar*.80, bh =3.                      ! graupel/hail mass parameters
    real, parameter :: ykc=.3e8,  ync= 2.                      ! fall speed parameters (cloud)
    real, parameter :: ykr=900.,  ynr=.84                      ! fall speed parameters (rain)
    real, parameter :: yks=132.,  yns=.74                      ! fall speed parameters (snow)
    real, parameter :: ykh=1.7e5, ynh=1.8                      ! fall speed parameters (graupel/hail)

    real, parameter :: zqcwth=2.e-6, zqcith=2.e-6              ! autoconversion thresholds
    real, parameter :: d0w=0.58e-4, d0i=0.125e-4               ! autoconversion diameters
    real, parameter :: zqpwth=1.e-4, ztcrit=269.               ! graupel/hail formation
    real, parameter :: ee_rcw=.5, ee_scw=.8, ee_sci=.4         ! collection efficiencies
    real, parameter :: ee_hcw=.5, ee_hci=.5, ee_rs =.7         ! collection efficiencies
    real, parameter :: fvent=.8, chi=2.26e-5, yka=2.43e-2
    real, parameter :: ymu=1.718e-5, schmidt=.6
    real, parameter :: zqmin=5.e-9                             ! minimum value of specific concentration
    real, parameter :: cccwc=1.26e9, kcwc=.308, yncmaxc=1.e9   ! continental cloud nucleation
    real, parameter :: cccwm=1.e8,   kcwm=.462, yncmaxm=1.e8   ! maritime cloud nucleation
    real, parameter :: refw=0.93, refi=0.21                    ! radar reflectivity

    real, dimension(nlon,nlev) :: zfsrain, zfssnow, zfshail, zflurn, zflusn, zfluha, zrodz, &
                                  zfsqcw, zfsqci, zfsncw, zfsnci, zflucw, zfluci, zflunw, zfluni

    integer jlon, jlat, jklev, jstep
    real zexpw, zexpi, z31, zgalfw1, zgalfwb, zgalfw2, zgalfi1, zgalfib, zgalfi2
    real zgar1, zgar2, zgar3, zgar4, zgas1, zgas2, zgas3, zgas4, zgah1, zgah2, zgah3, zgah4
    real zzcdw, zzcdi, zzafm, zzncw, zznci, cc_rcw, cc_scw, cc_sci, cc_hcw, cc_hci, cc_rs
    real zp, zro, zt, zqv, zqcw, zqci, zqpw, zqpi1, zqpi2, yncw, ynci
    real zt0t, zesk, zqsw, zqsi, ssw, ssi, yncwc, yncic, zqtot, zcoe, zctot, zh, zw, z1
    real zdqcw, zqdci, zdq, zdnn, zdw, zdi, zcond, zbase, zpote, zridu, zdcon, zdfre, zdsub, zeffi
    real zsubl, zarg, zgammacw, zgammaci, zbetw, zbeti, zdqtr, zdqmlt, zdqci
    real zdqtot, zlambdrl, zlambdhl, zlambdsl
    real zlambhm1, zlambsm1, zlambrm1, zrey, zdsmi, zintfr, zintfs, zintfh
    real zgacw1, zgacw3, zgaci1, zgaci3, zbeta, cccw, kcw, yncmax
    real zqcwnew, zncwnew, zqcinew, zncinew, zqpwnew, zqpi1new, zqpi2new, zdqwat, zdqice
    real zdegr1, zdegs1, zdegh1, zsshape, ztotref, zrainref, zsnowref, zhailref
    real zmelmax, zmeltsn, ynclow, zboul, zboui, zwes
    real zrh, ztw, zt0tw, zrip, zzqsw, ztwfg

! gamma functions

      zgalfw1 = eugamma (alfacw+1.        )
      zgalfwb = eugamma (alfacw+bcw+1.    )
      zgalfw2 = eugamma (alfacw-bcw+2.    )
      zgacw1  = eugamma (ync+bcw+alfacw+1.)
      zgacw3  = eugamma (ync+alfacw+1.    )

      zgalfi1 = eugamma (alfaci+1.        )
      zgalfib = eugamma (alfaci+bci+1.    )
      zgalfi2 = eugamma (alfaci-bci+2.    )
      zgaci1  = eugamma (ync+bci+alfaci+1.)
      zgaci3  = eugamma (ync+alfaci+1.    )

      zgar1 = eugamma (ynr+3.    )
      zgar2 = eugamma (br+1.     )
      zgar3 = eugamma (.5*ynr+2.5)
      zgar4 = eugamma (br+ynr+1. )

      zgas1 = eugamma (yns+3.    )
      zgas2 = eugamma (bs+1.     )
      zgas3 = eugamma (.5*yns+2.5)
      zgas4 = eugamma (bs+yns+1. )

      zgah1 = eugamma (ynh+3.    )
      zgah2 = eugamma (bh+1.     )
      zgah3 = eugamma (.5*ynh+2.5)
      zgah4 = eugamma (bh+ynh+1. )

! constants

      zexpw = 1./bcw
      zexpi = 1./bci
      zzcdw = dtstep*(alfacw+1.)*(zgalfw1/(zgalfwb*acw))**zexpw
      zzcdi = dtstep*(alfaci+1.)*(zgalfi1/(zgalfib*aci))**zexpi
      zzafm = 2.*pi*fvent*yka/yliw
      zzncw = dtstep*zgalfw2/zgalfw1/acw*(zgalfw1/(zgalfwb*acw))**(zexpw-1.)   !2m
      zznci = dtstep*zgalfi2/zgalfi1/aci*(zgalfi1/(zgalfib*aci))**(zexpi-1.)   !2m
      z31   = .31*schmidt**(1./3.)/sqrt(ymu)

! collection coefficients

      cc_rcw = .25*pi*ee_rcw*yn0r*ykr*zgar1*dtstep
      cc_rs  = .25*pi*ee_rs *yn0r*ykr*zgar1*dtstep
      cc_scw = .25*pi*ee_scw*yn0s*yks*zgas1*dtstep
      cc_sci = .25*pi*ee_sci*yn0s*yks*zgas1*dtstep
      cc_hcw = .25*pi*ee_hcw*yn0h*ykh*zgah1*dtstep
      cc_hci = .25*pi*ee_hci*yn0h*ykh*zgah1*dtstep

!-------------------------------------------------------------------------------------

      do 1000 jlat = 2, nlat-1

      zfsqcw = -.002
      zfsqci = -.002
      zfsncw = -.001
      zfsnci = -.001
      zfsrain = -.1
      zfssnow = -.01
      zfshail = -.1

      do jklev = 1, ntop
      do jlon = 2, nlon-1

      zt    = t   (jlon,jlat,jklev)
      zp    = p   (jlon,jlat,jklev)
      zqv   = q   (jlon,jlat,jklev)
      zqcw  = qcw (jlon,jlat,jklev)
      zqci  = qci (jlon,jlat,jklev)
      zqpw  = qpw (jlon,jlat,jklev)
      zqpi1 = qpi1(jlon,jlat,jklev)
      zqpi2 = qpi2(jlon,jlat,jklev)
      if (nlmic2) then
      yncw  = ncw (jlon,jlat,jklev)                       !2m
      ynci  = nci (jlon,jlat,jklev)                       !2m
      zw = .5*(w(jlon,jlat,jklev+1)+w(jlon,jlat,jklev))   !2m
      else
!      yncw  = yncw0*(1.-fmask(jlon,jlat))+.5*yncw0*fmask(jlon,jlat)
      yncw  = yncw0
      ynci  = ynci0
      endif

      zqtot = zqv +zqcw +zqci +zqpw +zqpi1 +zqpi2
      zro   = zp/((rd*(1.-zqtot)+rv*zqv)*zt)        ! air density
      zcoe  = sqrt(1./zro)                          ! coefficient to increase terminal velocity with height
      zboul = 1.
      zboui = zboul * max(.01, min((zp-100.e2)/300.e2, 1.))  ! to reduce high clouds

! cloud nucleation as a function of land/sea mask

      if (nlmic2) then
      if (fmask(jlon,jlat).gt.0.5) then
      cccw = cccwm        !2m
      kcw  = kcwm         !2m
      yncmax = yncmaxm    !2m
      else
      cccw = cccwc        !2m
      kcw  = kcwc         !2m
      yncmax = yncmaxc    !2m
      endif
      ynclow = yncmax*.05 !2m
      endif

! calculation of enthalpy zh before microphysical processes

      zctot =  (1.-zqtot)*cpd + zqtot*ci
      zh = zctot*zt + (yliv+(cpv-ci)*(zt-tzer))*zqv + (yliw+(cw-ci)*(zt-tzer))*(zqcw+zqpw)

! saturation specific humidity with respect to water and ice

      zt0t = tzer/zt
      zesk = ezer*exp(-ccw1*log(zt0t)+ccw2*(1.-zt0t)) !  saturated pressure over water (in Pa)
      zqsw = zesk*eps/(zp+zesk*(eps-1.))
      if (zt.lt.tzer) then
      zesk = ezer*exp(-cci1*log(zt0t)+cci2*(1.-zt0t)) !  saturated pressure over ice (in Pa)
      zqsi = zesk*eps/(zp+zesk*(eps-1.))
      else
      zqsi = zqsw
      endif

! 1)  nucleation of cloud water

      if (nlmic2) then
      ssw = (zqv-zqsw)/zqsw                                                                     !2m
      if (ssw.gt.0.and.yncw.lt.yncmax.and.zw.gt.0.) then                                        !2m
      yncwc = yncw + dtstep*.01*cccw*kcw*ssw*zw   ! modified for numerical reasons              !2m
      yncwc = min (yncwc, yncmax)                                                               !2m
      zdqcw = min (1.e-12*(yncwc-yncw), zqv)                                                    !2m
      yncw  = yncw + zdqcw*1.e12                                                                !2m
      zqcw  = zqcw + zdqcw                                                                      !2m
      zqv   = zqv - zdqcw                                                                       !2m
      endif                                                                                     !2m
      endif

! 2) evaporation-condensation rate of cloud water

      if (zt.gt.233.) then
      if (zqv.gt.zqsw.or.zqcw.gt.zqmin) then
      zdw = zqv-zqsw
      z1  = zqsw*ylwv**2*chi*zro/(yka*rv*zt**2)
      zcond = 2.*pi*fvent*chi*zdw/(1.+z1) ! dm(D)/dt = D*zro*zcond
      zbase = max(zqcw, 1.e-9)*zro/yncw
      zpote = zbase**zexpw
      zdcon = zcond*zzcdw*yncw*zpote
      zridu = 1./(1.+zqsw*ylwv**2/(rv*zctot*zt**2))
        if (zdw.gt.0.) then
        zdcon =  min ( zdcon, zdw*zridu)
        else
        zdcon = -min (-zdcon, zqcw, -zdw*zridu)
        endif
      zqcw = zqcw + zdcon                ! update cloud water
      zqv  = zqv  - zdcon

! update number density of cw

      if (nlmic2) then
      yncw = yncw + zro*zdcon*zzncw/(zzcdw*zbase)           !2m
      yncw = min (yncw, yncmax)                             !2m
      yncw = max (yncw, 100.)                               !2m
      endif
      endif
      endif

! 3) spontaneous freezing (homogeneous nucleation) of cloud water

      if (zt.lt.245.and.zqcw.gt.zqmin) then
      zbase = zqcw*zro/yncw
      zpote = zbase**zexpw
      zdfre = zzafm*(tzer-zt)/zro*zzcdw*yncw*zpote
      zdfre = min (zdfre, zqcw)
      zqcw  = zqcw - zdfre
      zqci  = zqci + zdfre

! update number density of cw and ci

      if (nlmic2) then
      zdnn = zro*zdfre*zzncw/(zzcdw*zbase)       !2m
      zdnn = min (zdnn, yncw)                    !2m
      yncw = max (yncw-zdnn, 100.)               !2m
      ynci = ynci + zdnn                         !2m
      ynci = min (ynci, yncmax)                  !2m
      endif
      endif

! 4) heterogeneous nucleation of cloud ice (formation of pristine crystals)

      if (zt.lt.269.) then
      ssi = (zqv-zqsi)/zqsi
        if (ssi.gt.0) then
        yncic = 1.e4*min (1.,.25*(269.-zt) ) * exp(-.639+12.96*ssi)
          if (nlmic2) then
          if (yncic.gt.ynci) then                    !2m
          zdqci = min (1.e-12*(yncic-ynci), zqv)     !2m
          ynci = ynci + zdqci*1.e12                  !2m
          zqci  = zqci + zdqci                       !2m
          zqv = zqv - zdqci                          !2m
          endif                                      !2m
          else
          if (1.e-12*yncic.gt.zqci) then
          zdqci = min (1.e-12*yncic-zqci, zqv)
          zqci  = zqci + zdqci
          zqv   = zqv  - zdqci
          endif
          endif
        endif
      endif

! 5) sublimation rate of cloud ice in both directions (part of Bergeron-Findeisen process)

      if (zt.lt.tzer) then
      if (zqv.gt.zqsi.or.zqci.gt.zqmin) then
      zdi= zqv-zqsi
      z1 = zqsi*yliv**2*chi*zro/(yka*rv*zt**2)
      zsubl = 2.*pi*fvent*zdi*chi/(1.+z1)
      zbase = zqci*zro/ynci
      zpote = zbase**zexpi
      zridu  = 1./(1.+zqsi*yliv**2/(rv*zctot*zt**2))
        if (zdi.gt.0.) then
          if(zt.le.tzer-15..and.zt.ge.tzer-45.) then
          zeffi = 1.- 1./6.*((zt-tzer+45.)/15.)**2
          elseif (zt.gt.tzer-15.) then
          zeffi = 1./3.*((zt-tzer)/15.)**2
          else
          zeffi = 1.
          endif
        zdsub = zsubl*zzcdi*ynci*zpote*zeffi
        zdsub = min(zdsub, zdi*zridu)
        else
        zdsub = zsubl*zzcdi*ynci*zpote
        zdsub = -min(-zdsub, zqci, -zdi*zridu)
        endif
      zqv  = zqv  - zdsub
      zqci = zqci + zdsub

! update number density of ci

      if (nlmic2) then
      ynci = ynci + zro*zdsub*zznci/(zzcdi*zbase)   !2m
      ynci = min (ynci, yncmax)                     !2m
      ynci = max (ynci, 100.)                       !2m
      endif
      endif
      endif

! 6) melting of cloud ice

      if (zt.gt.tzer.and.zqci.gt.zqmin) then
      yncic = max (ynci, 1.e4)
      zbase = zqci*zro/yncic
      zpote = zbase**zexpi
      zdsmi = zzafm*(zt-tzer)/zro*zzcdi*yncic*zpote
      zdsmi = min (zdsmi, zqci)
      zqcw  = zqcw + zdsmi
      zqci  = zqci - zdsmi

! update number density of ci and cw

      if (nlmic2) then
      zdnn = zro*zdsmi*zznci/(zzcdi*zbase)    !2m
      zdnn = min (zdnn, ynci)                 !2m
      ynci = max (ynci-zdnn, 100.)            !2m
      yncw = yncw + zdnn                      !2m
      yncw = min (yncw, yncmax)               !2m
      endif
      endif

!-------------------------------------------------------------------------------------------------------------

! 7) cloud water converted into cloud ice due to cloud ice-water interaction

     if (zt.lt.tzer-.5.and.zqcw.gt.zqmin.and.zqci.gt.zqmin) then
     zdfre = min(5.*zro**2*dtstep*zqcw*zqci, .5*zqcw)
     zqcw = zqcw - zdfre
     zqci = zqci + zdfre
     endif

! 8) autoconversion of cloud water into rain

      if (nlmic2.and.yncw.gt.ynclow) yncw = yncw-dtstep/120.*(yncw-ynclow)  ! self-collection  !2m

      if (zqcw.gt.zqcwth) then
      zbetw = ((zgalfwb*acw*yncw)/(zqcw*zro*zgalfw1))**zexpw
      zarg = d0w*zbetw*zboul
      if (zarg.lt.25.) then
      zgammacw = .5*(1.-tanh((zarg-8./zarg-8.8)/4.3))  ! partial gamma(10,zarg)
      zdq  = zgammacw*zqcw * min(dtstep/27., 1.)
      zqcw = zqcw - zdq
      zqpw = zqpw + zdq
      endif
      endif

! 9) autoconversion of cloud ice into snow

      if (nlmic2.and.ynci.gt.ynclow.and.zt.lt.tzer) ynci=ynci-dtstep/120.*(ynci-ynclow) ! self-collection !2m

      if (zqci.gt.zqcith.and.zt.lt.tzer) then
      zbeti = ((zgalfib*aci*ynci)/(zqci*zro*zgalfi1))**zexpi
      zarg = d0i*zbeti*zboui
      if (zarg.lt.20.) then
      zgammaci = .5*(1.-tanh((zarg-8./zarg-4.9)/3.9))
      zdq   = zgammaci*zqci * min(dtstep/27., 1.)
      zqci  = zqci  - zdq
      zqpi1 = zqpi1 + zdq
      endif
      endif

!-------------------------------------------------------------------------------------------------------------

! 10) graupel/hail from freezing rain

      if (zqpw.gt.zqpwth.and.zt.lt.ztcrit) then
       if (zt.gt.245.) then
       zdq = zqpw*min(0.5, dtstep*5.e-5*(ztcrit-zt)**2)
       zqpw  = zqpw  - zdq
       zqpi2 = zqpi2 + zdq
       else
       zdq   = zqpw
       zqpw  = 0.
       zqpi2 = zqpi2 + zdq
       endif
      endif

! 11) melting of snow depending on wet bulb temperature

      if (zt.gt.tzer.and.(zqpi1.gt.zqmin.or.zqpi2.gt.zqmin)) then
      zrh = zqv/zqsw
       if(zrh.ge.1.) then
        ztw = zt
       else
        ztwfg = zt + (0.04*(zt-276.) + 1.)*(8.4 + 5.6e-5*(800.e2 - zp))*(zrh - 1.) ! 1st guess of wet bulb t
        zt0tw = tzer/ztwfg
        zesk = ezer*exp(-ccw1*log(zt0tw)+ccw2*(1.-zt0tw)) ! saturated pressure over water (in Pa)
        zzqsw = zesk*eps/(zp+zesk*(eps-1.))
        ztw = zt + yliv/cpd*(zqv-zzqsw) ! wet bulb temperature for ice
        ztw = 0.5*ztw + 0.5*ztwfg       ! empirical approx.
       endif
       if (ztw.gt.tzer.and.zqpi1.gt.zqmin) then
       zlambdsl = log(zro*zqpi1/(yn0s*as*zgas2))/(bs+1.)
       zlambsm1 = exp(zlambdsl)
       zrey   = z31*sqrt(zcoe*yks*zro)*zgas3*exp((.5*yns+2.5)*zlambdsl)
       zintfs = yn0s*(.78*zlambsm1**2 + zrey)    ! integral of d*f(d)*n(d)
       zdq    = 2.*pi*yka/(yliw*zro)*(ztw-tzer)*zintfs*dtstep
       zdq    = min (zdq, (ztw-tzer)*zctot/yliw, zqpi1)
       zqpi1  = max(zqpi1 - zdq, 0.)
       zqpw   = zqpw  + zdq
       endif

! 12) melting of graupel/hail depending on wet bulb temperature

       if (ztw.gt.tzer.and.zqpi2.gt.zqmin) then
       zlambdhl = log(zro*zqpi2/(yn0h*ah*zgah2))/(bh+1.)
       zlambhm1 = exp(zlambdhl)
       zrey   = z31*sqrt(zcoe*ykh*zro)*zgah3*exp((.5*ynh+2.5)*zlambdhl)
       zintfh = yn0h*(.78*zlambhm1**2 + zrey)   ! integral of d*f(d)*n(d)
       zdq    = 2.*pi*yka/(yliw*zro)*(ztw-tzer)*zintfh*dtstep
       zdq    = min (zdq, (ztw-tzer)*zctot/yliw, zqpi2)
       zqpi2  = max(zqpi2 - zdq, 0.)
       zqpw   = zqpw  + zdq
       endif
      endif

!-------------------------------------------------------------------------------------------------------------

! 13) sublimation of snow

      if (zqpi1.gt.zqmin.and.zqv.lt.zqsi) then
      zlambdsl = log(zro*zqpi1/(yn0s*as*zgas2))/(bs+1.)
      zlambsm1 = exp(zlambdsl)
      zrey   = z31*sqrt(zcoe*yks*zro)*zgas3*exp((.5*yns+2.5)*zlambdsl)
      zintfs = yn0s*(.78*zlambsm1**2 + zrey)      ! integral of d*f(d)*n(d)
      zdi    = zqsi-zqv
      z1     = zqsi*yliv**2*chi*zro/(yka*rv*zt**2)
      zdq    = 2.*pi*zdi*chi*zintfs/(1.+z1)*dtstep
      zdq    = min(zdq, .5*zqpi1, zdi)
      zqv    = zqv  + zdq
      zqpi1  = zqpi1- zdq
      endif

! 14) sublimation of graupel/hail

      if (zqpi2.gt.zqmin.and.zqv.lt.zqsi) then
      zlambdhl = log(zro*zqpi2/(yn0h*ah*zgah2))/(bh+1.)
      zlambhm1 = exp(zlambdhl)
      zrey   = z31*sqrt(zcoe*ykh*zro)*zgah3*exp((.5*ynh+2.5)*zlambdhl)
      zintfh = yn0h*(.78*zlambhm1**2 + zrey)      ! integral of d*f(d)*n(d)
      zdi    = zqsi-zqv
      z1     = zqsi*yliv**2*chi*zro/(yka*rv*zt**2)
      zdq    = 2.*pi*zdi*chi*zintfh/(1.+z1)*dtstep
      zdq    = min(zdq, .5*zqpi2, zdi)
      zqv    = zqv  + zdq
      zqpi2  = zqpi2- zdq
      endif

! 15) evaporation of rain

      if (zqpw.gt.zqmin.and.zqv.lt.zqsw) then
      zlambdrl = log(zro*zqpw/(yn0r*ar*6.))/(br+1.)
      zlambrm1 = exp(zlambdrl)
      zrey   = z31*sqrt(zcoe*ykr*zro)*zgar3*exp((.5*ynr+2.5)*zlambdrl)
      zintfr = yn0r*(.78*zlambrm1**2 + zrey)           ! integral of d*f(d)*n(d)
      zdw    = zqsw-zqv
      z1     = zqsw*ylwv**2*chi*zro/(yka*rv*zt**2)
      zdq    = 2.*pi*zdw*chi*zintfr/(1.+z1)*dtstep
      zdq    = min(zdq, .5*zqpw, zdw)
      zqv    = zqv  + zdq
      zqpw   = zqpw - zdq
      endif

!-------------------------------------------------------------------------------------------------------------
!  fall of precipitation - collection processes
!-------------------------------------------------------------------------------------------------------------

! 16) rain from cloud water and rain (indep. of temp.)

      if (zqpw.gt.zqmin.and.zqcw.gt.zqmin) then
      zlambrm1 = (zro*zqpw/(yn0r*ar*6.))**((ynr+3.)/(br+1.))
      zdq  = cc_rcw*zcoe*zqcw*zlambrm1
      zdq  = min (zdq, zqcw*.5)  ! intercepted cloud water
      zqcw = zqcw - zdq
      zqpw = zqpw + zdq
      if (nlmic2) then
      zdnn = yncw*zdq/zqcw                                   !2m
      yncw = yncw - zdnn  ! cloud spectrum does not change   !2m
      yncw = max (yncw, 100.)                                !2m
      endif
      endif

! 17) rain and snow-graupel/hail interaction
!     part of supercooled rain is converted to graupel/hail below tzer-1. (by enthalpy conservation);
!     else, intercepted snow-graupel/hail is converted to rain

      if (zt.lt.tzer-1..or.zt.gt.tzer+.5) then
        if (zqpw.gt.zqmin.and.(zqpi1+zqpi2).gt.zqmin) then
        zlambrm1 = (zro*zqpw/(yn0r*ar*6.))**((ynr+3.)/(br+1.))
        zdq    = cc_rs*zcoe*(zqpi1+zqpi2)*zlambrm1
        zdq    = min (zdq, (zqpi1+zqpi2)*.5)                  ! intercepted snow and graupel/hail
        zqpi1  = zqpi1 - zdq*zqpi1/(zqpi1+zqpi2)
        zqpi2  = zqpi2 - zdq*zqpi2/(zqpi1+zqpi2)
          if (zt.lt.tzer-1.) then
          zdqice = min(.5*zqpw, (zqpw+zdq)*(tzer-zt)*ci/yliw) ! supercooled rain converted to graupel/hail
          zqpw   = zqpw  - zdqice
          zqpi2  = zqpi2 + zdq + zdqice
          else
          zqpw   = zqpw  + zdq
          endif
        endif
      endif

! 18) cloud ice intercepted by snow below freezing

      if (zqpi1.gt.zqmin.and.zqci.gt.zqmin.and.zt.lt.tzer) then
      zlambsm1 = (zro*zqpi1/(yn0s*as*zgas2))**((yns+3.)/(bs+1.))
      zdq   = cc_sci*zcoe*zqci*zlambsm1
      zdq   = min (zdq, zqci*.5)     ! intercepted cloud ice
      zqci  = zqci  - zdq
      zqpi1 = zqpi1 + zdq
      if (nlmic2) then
      zdnn  = ynci*zdq/zqci                                  !2m
      ynci  = ynci - zdnn  ! cloud spectrum does not change  !2m
      ynci  = max (ynci, 100.)                               !2m
      endif
      endif

! 19) snow and cloud water interaction: riming of snow below tzer-1.
!     snow melting due to enthalpy of cloud water above tzer

      if (zqpi1.gt.zqmin.and.zqcw.gt.zqmin) then
      zlambsm1 = (zro*zqpi1/(yn0s*as*zgas2))**((yns+3.)/(bs+1.))
      zdq  = cc_scw*zcoe*zqcw*zlambsm1
      zdq  = min (zdq, zqcw*.5)    ! intercepted cloud water
      zrip = min(zqcw*1.3e3, 1.)   ! sets redistribution of rimed cloud water into snow and graupel/hail
      zqcw = zqcw  - zdq
      if (nlmic2) then
      zdnn = yncw*zdq/zqcw                                         !2m
      yncw = yncw - zdnn       ! cloud spectrum does not change    !2m
      yncw = max (yncw, 100.)                                      !2m
      endif
        if (zt.lt.tzer-1.) then  ! riming
        zqpi1  = zqpi1 + (1.-zrip)*zdq
        zqpi2  = zqpi2 + zrip*zdq
        else
        zdqmlt = 0.
        if (zt.gt.tzer) zdqmlt = min(.5*zqpi1, (zdq+zqpi1)*(zt-tzer)*ci/yliw)
        zqpi1  = zqpi1 - zdqmlt
        zqpw   = zqpw  + zdqmlt + zdq
        endif
      endif

! 20) graupel/hail from cloud ice and graupel/hail below freezing

      if (zqpi2.gt.zqmin.and.zqci.gt.zqmin.and.zt.lt.tzer.and.zt.gt.253) then ! negleg. for t < -20 C (Houze)
      zlambhm1 = (zro*zqpi2/(yn0h*ah*zgah2))**((ynh+3.)/(bh+1.))
      zdq = cc_hci*zcoe*zqci*zlambhm1 *30./(30.+(tzer-zt)**2)  ! modulation decreasing with t (Houze)
      zdq = min(zdq, zqci*.5)   ! intercepted cloud ice
      zqci  = zqci  - zdq
      zqpi2 = zqpi2 + zdq
      if (nlmic2) then
      zdnn  = ynci*zdq/zqci                                    !2m
      ynci  = ynci - zdnn  ! cloud spectrum does not change    !2m
      ynci  = max (ynci, 100.)                                 !2m
      endif
      endif

! 21) graupel/hail from freezing cloud water
!     rain from cloud water and melting graupel/hail

      if (zqpi2.gt.zqmin.and.zqcw.gt.zqmin) then
      zlambhm1 = (zro*zqpi2/(yn0h*ah*zgah2))**((ynh+3.)/(bh+1.))
      zdq = cc_hcw*zcoe*zqcw*zlambhm1
      zdq = min(zdq, zqcw*.5)                                 ! intercepted cloud water
      zqcw = zqcw - zdq
      if (nlmic2) then
      zdnn = yncw*zdq/zqcw                                   !2m
      yncw = yncw - zdnn  ! cloud spectrum does not change   !2m
      yncw = max (yncw, 100.)                                !2m
      endif
        if (zt.lt.tzer) then
        zdqice = min(zdq, (zqpi2+zdq)*(tzer-zt)*ci/yliw)      ! cloud water converted to graupel/hail
        zqpi2  = zqpi2 + zdqice
        zqpw   = zqpw  + zdq - zdqice
        else
        zdqmlt = min(.5*zqpi2, (zqpi2+zdq)*(zt-tzer)*ci/yliw) ! hail melting due to enthalpy of cloud water
        zqpi2  = zqpi2 - zdqmlt
        zqpw   = zqpw  + zdq + zdqmlt
        endif
      endif

! computation of terminal fall speeds of cloud, rain, snow, and graupel/hail

      if (zqcw.gt.zqmin) then
      zbeta  = (yncw*acw*zgalfwb/(zro*zqcw*zgalfw1))**(zexpw*ync)
      zfsqcw(jlon,jklev) = -ykc*zgacw1/(zgalfwb*zbeta)
      if (nlmic2) zfsncw(jlon,jklev) = zfsqcw(jlon,jklev)*zgalfwb*zgacw3/(zgacw1*zgalfw1)    !2m
      endif
      if (zqci.gt.zqmin) then
      zbeta  = (ynci*aci*zgalfib/(zro*zqci*zgalfi1))**(zexpi*ync)
      zfsqci(jlon,jklev) = -ykc*zgaci1/(zgalfib*zbeta)
      if (nlmic2) zfsnci(jlon,jklev) = zfsqci(jlon,jklev)*zgalfib*zgaci3/(zgaci1*zgalfi1)    !2m
      endif
      if (zqpw.gt.zqmin) then
      zfsrain(jlon,jklev) = -zcoe*ykr*zgar4/6.   *(zro*zqpw /(yn0r*ar*6.   ))**(ynr/(br+1.))
      endif
      if (zqpi1.gt.zqmin) then
      zfssnow(jlon,jklev) = -zcoe*yks*zgas4/zgas2*(zro*zqpi1/(yn0s*as*zgas2))**(yns/(bs+1.))
        if (zt.gt.tzer) then
        zwes = min(1.,(zt-tzer)/3.5) ! larger fall speed of melting snow
        zfssnow(jlon,jklev) = (1.-zwes)*zfssnow(jlon,jklev) + zwes*zfsrain(jlon,jklev) ! speed of melting snow
        endif
      endif
      if (zqpi2.gt.zqmin) then
      zfshail(jlon,jklev) = -zcoe*ykh*zgah4/zgah2*(zro*zqpi2/(yn0h*ah*zgah2))**(ynh/(bh+1.))
      endif

! model variables after microphysical processes
! new temperature by conservation of enthalpy

      t(jlon,jlat,jklev) = (zh -(yliv-(cpv-ci)*tzer)*zqv -(yliw-(cw-ci)*tzer)*(zqcw+zqpw))/  &
                           (zctot +(cpv-ci)*zqv +(cw-ci)*(zqcw+zqpw))

      q   (jlon,jlat,jklev) = zqv
      qcw (jlon,jlat,jklev) = zqcw
      qci (jlon,jlat,jklev) = max(zqci, 0.)
      qpw (jlon,jlat,jklev) = zqpw
      qpi1(jlon,jlat,jklev) = max(zqpi1, 0.)
      qpi2(jlon,jlat,jklev) = max(zqpi2, 0.)
      if (nlmic2) ncw (jlon,jlat,jklev) = yncw                 !2m
      if (nlmic2) nci (jlon,jlat,jklev) = ynci                 !2m

      enddo  ! end loop on longitude
      enddo  ! end loop on levels

!  fall of precipitation by a backward-upstream scheme
!  accumulation of precipitation at the ground

      do jklev = 1, ntop+1
      do jlon = 2, nlon-1
      zro = p(jlon,jlat,jklev)/(rd*t(jlon,jlat,jklev))
      zrodz(jlon,jklev)  = zro*dz/fmz(jlon,jlat,jklev)
      zflucw(jlon,jklev) = zfsqcw (jlon,jklev)*zro*dtstep
      zfluci(jlon,jklev) = zfsqci (jlon,jklev)*zro*dtstep
      if (nlmic2) zflunw(jlon,jklev) = zfsncw (jlon,jklev)*zro*dtstep     !2m
      if (nlmic2) zfluni(jlon,jklev) = zfsnci (jlon,jklev)*zro*dtstep     !2m
      zflurn(jlon,jklev) = zfsrain(jlon,jklev)*zro*dtstep
      zflusn(jlon,jklev) = zfssnow(jlon,jklev)*zro*dtstep
      zfluha(jlon,jklev) = zfshail(jlon,jklev)*zro*dtstep
      enddo
      enddo

      do jklev = ntop, 1, -1
      do jlon = 2, nlon-1
      zqcwnew  = (qcw (jlon,jlat,jklev)*zrodz(jlon,jklev)-qcw (jlon,jlat,jklev+1)*zflucw(jlon,jklev+1)) &
                /(zrodz(jlon,jklev)-zflucw(jlon,jklev))
      zqcinew  = (qci (jlon,jlat,jklev)*zrodz(jlon,jklev)-qci (jlon,jlat,jklev+1)*zfluci(jlon,jklev+1)) &
                /(zrodz(jlon,jklev)-zfluci(jlon,jklev))
      if (nlmic2) then
      zncwnew  = (ncw (jlon,jlat,jklev)*zrodz(jlon,jklev)-ncw (jlon,jlat,jklev+1)*zflunw(jlon,jklev+1)) &  !2m
                /(zrodz(jlon,jklev)-zflunw(jlon,jklev))                                                    !2m
      zncinew  = (nci (jlon,jlat,jklev)*zrodz(jlon,jklev)-nci (jlon,jlat,jklev+1)*zfluni(jlon,jklev+1)) &  !2m
                /(zrodz(jlon,jklev)-zfluni(jlon,jklev))                                                    !2m
      endif
      zqpwnew  = (qpw (jlon,jlat,jklev)*zrodz(jlon,jklev)-qpw (jlon,jlat,jklev+1)*zflurn(jlon,jklev+1)) &
                /(zrodz(jlon,jklev)-zflurn(jlon,jklev))
      zqpi1new = (qpi1(jlon,jlat,jklev)*zrodz(jlon,jklev)-qpi1(jlon,jlat,jklev+1)*zflusn(jlon,jklev+1)) &
                /(zrodz(jlon,jklev)-zflusn(jlon,jklev))
      zqpi2new = (qpi2(jlon,jlat,jklev)*zrodz(jlon,jklev)-qpi2(jlon,jlat,jklev+1)*zfluha(jlon,jklev+1)) &
                /(zrodz(jlon,jklev)-zfluha(jlon,jklev))

!  heat capacity of precipitation (cloud excluded)

      zdqwat=abs(zflurn(jlon,jklev+1)*qpw (jlon,jlat,jklev+1))
      zdqice=abs(zflusn(jlon,jklev+1)*qpi1(jlon,jlat,jklev+1))+abs(zfluha(jlon,jklev+1)*qpi2(jlon,jlat,jklev+1))
      t(jlon,jlat,jklev)=( cpd*zrodz(jlon,jklev)*t(jlon,jlat,jklev)+(zdqwat*cw+zdqice*ci)*t(jlon,jlat,jklev+1))&
                        /(cpd*zrodz(jlon,jklev)+zdqwat*cw+zdqice*ci)

      qcw (jlon,jlat,jklev) = zqcwnew
      qci (jlon,jlat,jklev) = zqcinew
      if (nlmic2) ncw (jlon,jlat,jklev) = zncwnew      !2m
      if (nlmic2) nci (jlon,jlat,jklev) = zncinew      !2m
      qpw (jlon,jlat,jklev) = zqpwnew
      qpi1(jlon,jlat,jklev) = zqpi1new
      qpi2(jlon,jlat,jklev) = zqpi2new

      enddo
      enddo

!  instantaneous rain and snow (mm) reaching the ground

      raini(2:nlon-1,jlat) = -zflurn(2:nlon-1,1)*qpw (2:nlon-1,jlat,1) -zflucw(2:nlon-1,1)*qcw(2:nlon-1,jlat,1)
      snowi(2:nlon-1,jlat) = -zflusn(2:nlon-1,1)*qpi1(2:nlon-1,jlat,1) -zfluci(2:nlon-1,1)*qci(2:nlon-1,jlat,1)

!  melting of residual snow (not graupel/hail) at the lowest level, after fall
!  latent heat of fusion subtracted from the lowest level

      do jlon = 2, nlon-1
      if (t(jlon,jlat,1).gt.tzer.and.t(jlon,jlat,1).lt.tzer+10.) then
       zt = t(jlon,jlat,1)
       zt0t = tzer/zt
       zesk = ezer*exp(-ccw1*log(zt0t)+ccw2*(1.-zt0t)) !  saturated pressure over water (in Pa)
       zqsw = zesk*eps/(p(jlon,jlat,1)+zesk*(eps-1.))
       zrh = q(jlon,jlat,1)/zqsw ! relat. humidity
        if(zrh.ge.1.) then
         ztw = zt
        else
         ztwfg = zt + (0.04*(zt-276.)+1.)*(8.4 + 5.6e-5*(800.e2 - p(jlon,jlat,1)))*(zrh - 1.) ! 1st guess of wet bulb t
         zt0tw = tzer/ztwfg
         zesk = ezer*exp(-ccw1*log(zt0tw)+ccw2*(1.-zt0tw)) ! saturated pressure over water (in Pa)
         zzqsw = zesk*eps/(p(jlon,jlat,1)+zesk*(eps-1.))
         ztw = zt + yliv/cpd*(q(jlon,jlat,1)-zzqsw) ! wet bulb temperature for ice
         ztw = 0.5*ztw + 0.5*ztwfg                  ! empirical approx.
        endif
       if (ztw.gt.tzer+0.7) then  ! offset value empirical but important for snowfall
       zmelmax = (ztw-tzer)*cpd*zrodz(jlon,1)/yliw
       zmeltsn = min (snowi(jlon,jlat), zmelmax)
       raini(jlon,jlat) = raini(jlon,jlat) + zmeltsn
       t(jlon,jlat,1) = t(jlon,jlat,1) - yliw*zmeltsn/(cpd*zrodz(jlon,1))
       snowi(jlon,jlat) = max(snowi(jlon,jlat) - zmeltsn, 0.)
       endif
      endif
      enddo

!   graupel/hail added to snowi

      snowi(2:nlon-1,jlat) = snowi(2:nlon-1,jlat) -zfluha(2:nlon-1,1)*qpi2(2:nlon-1,jlat,1)

 1000 continue

!============================================================================
!  Radar reflectivity of precipitation in dbz
!============================================================================

    if (nlradar.and.mod(jstep,nhist).eq.0) then

!  for rain

    zgar1 = eugamma (2.*br+1.)
    zdegr1 = (2.*br+1.)/(br+1.)

!  for snow

    zgas1 = eugamma (2.*bs+1.)
    zdegs1 = (2.*bs+1.)/(bs+1.)

!  for hail

    zgah1 = eugamma (2.*bh+1.)
    zdegh1 = (2.*bh+1.)/(bh+1.)

!  deviation of snow particles from spherical shape

    zsshape = (as*6./(pi*800.))**2.

    do jklev = 1, ntop
    do jlat = 2, nlat-1
    do jlon = 2, nlon-1

    zro = p(jlon,jlat,jklev)/(rd*tvirt(jlon,jlat,jklev))

    zrainref=0.
    if (qpw(jlon,jlat,jklev).gt.zqmin) then
    zrainref = yn0r*refw*zgar1*exp(-zdegr1*log(yn0r*ar*zgar2/(zro*qpw(jlon,jlat,jklev))))
    endif

    zsnowref=0.
    if (qpi1(jlon,jlat,jklev).gt.zqmin) then
    zsnowref = yn0s*refi*zsshape*zgas1*exp(-zdegs1*log(yn0s*as*zgas2/(zro*qpi1(jlon,jlat,jklev))))
    endif

    zhailref=0.
    if (qpi2(jlon,jlat,jklev).gt.zqmin) then
    zhailref = yn0h*refi*zgah1*exp(-zdegh1*log(yn0h*ah*zgah2/(zro*qpi2(jlon,jlat,jklev))))
    endif

    ztotref = zrainref + zsnowref + zhailref

! convert from m**6/m**3 to mm**6/m**3

    ztotref = ztotref*1.e18

! convert from mm**6/m**3 to dbz

    rradar(jlon,jlat,jklev) = max(radarmval, 10.*alog10(ztotref+1.e-20))

    enddo
    enddo
    enddo

    endif

    return

    CONTAINS
    function eugamma (x)
    implicit none
    integer :: nx, ix, nfrx
    real :: eugamma, x, gfact, xx, gfrac
    real :: eug(0:100)
    data eug/1.00000,0.99433,0.98884,0.98355,0.97844,0.97350, &    !1.0-1.05
                     0.96874,0.96415,0.95973,0.95546,0.95135, &    !1.06-1.10
                     0.94740,0.94359,0.93993,0.93642,0.93304, &    !1.11-1.15
                     0.92980,0.92670,0.92373,0.92089,0.91817, &    !1.16-1.20
                     0.91558,0.91311,0.91075,0.90852,0.90640, &    !1.21-1.25
                     0.90440,0.90250,0.90072,0.89904,0.89747, &    !1.26-1.30
                     0.89600,0.89464,0.89338,0.89222,0.89115, &    !1.31-1.35
                     0.89018,0.88931,0.88854,0.88785,0.88726, &    !1.36-1.40
                     0.88676,0.88636,0.88604,0.88581,0.88566, &    !1.41-1.45
                     0.88560,0.88563,0.88575,0.88595,0.88623, &    !1.46-1.50
                     0.88659,0.88704,0.88757,0.88818,0.88887, &    !1.51-1.55
                     0.88964,0.89049,0.89142,0.89243,0.89352, &    !1.56-1.60
                     0.89468,0.89592,0.89724,0.89864,0.90012, &    !1.61-1.65
                     0.90167,0.90330,0.90500,0.90678,0.90864, &    !1.66-1.70
                     0.91057,0.91258,0.91467,0.91683,0.91906, &    !1.71-1.75
                     0.92137,0.92376,0.92623,0.92877,0.93138, &    !1.76-1.80
                     0.93408,0.93685,0.93969,0.94261,0.94561, &    !1.81-1.85
                     0.94869,0.95184,0.95507,0.95838,0.96177, &    !1.86-1.90
                     0.96523,0.96877,0.97240,0.97610,0.97988, &    !1.91-1.95
                     0.98374,0.98768,0.99171,0.99581,1.00000/      !1.96-2.00

    nx=floor(x)
    xx=x
    gfact=1.0
    do ix=1,nx-1
    xx=xx-1.0
    gfact=gfact*xx
    enddo
    nfrx=nint((xx-1.)*100.)
    gfrac=eug(nfrx)
    eugamma = gfact*gfrac
    end function eugamma

    end subroutine micro2m
!###############################################################################################################
    subroutine wafone (p, dt)

    use mod_moloch, only : nlon, nlat, nlev, nlevp1, nlonm1, nlatm1, u, v, s, fmyu, clv, dx, dy, dz, a, &
                           dlon, dlat, rdeno, ip_e, ip_n, ip_s, ip_w, ip_null, nltwice
    implicit none
    real p(nlon,nlat,nlev), zpbw(nlon), wfw(nlon,nlevp1), zpby(nlon,nlat), wz(nlon,0:nlat+1,nlev)
    real p0(0:nlon+1,nlat,nlev), r, b, zphi, denr, dt, zamu, zdv, zcost, zcostx, zcosty, zhxvt, zhxvtn
    integer jlon, jlat, jklev, j1, j1m1, is

!----------------------
!  Vertical advection
!----------------------

    zcost = dt/dz
    if (nltwice) zcost = .5*dt/dz
    do jlon = 2, nlonm1
    wfw (jlon,1)=0.
    wfw (jlon,nlevp1)=0.
    enddo

    do 10 jlat = 2, nlatm1
    do jklev = 2, nlev
    do jlon = 2, nlonm1
      zamu = s(jlon,jlat,jklev)*zcost
      if (zamu.ge.0.) then
      is=1
      j1=jklev-1
      j1m1=j1-1
      if(j1m1.lt.1) j1m1=1
      else
      is=-1
      j1=jklev+1
      j1m1=j1-1
      if(j1.gt.nlev) j1=nlev
      endif
    r = rdeno(p(jlon,jlat,j1), p(jlon,jlat,j1m1), p(jlon,jlat,jklev), p(jlon,jlat,jklev-1))
    b = max(0., min(2., max(r, min(2.*r,1.))))
    zphi = is+zamu*b-is*b
    wfw(jlon,jklev) = 0.5*zamu*((1.+zphi)*p(jlon,jlat,jklev-1)+(1.-zphi)*p(jlon,jlat,jklev))
    enddo
    enddo
    do jklev = 1, nlev
    do jlon = 2, nlonm1
    zdv = (s(jlon,jlat,jklev+1)-s(jlon,jlat,jklev))*zcost
    wz(jlon,jlat,jklev) = p(jlon,jlat,jklev)+wfw(jlon,jklev)-wfw(jlon,jklev+1)+p(jlon,jlat,jklev)*zdv
    enddo
    enddo

    if (nltwice) then
    do jklev = 2, nlev
    do jlon = 2, nlonm1
      zamu = s(jlon,jlat,jklev)*zcost
      if (zamu.ge.0.) then
      is=1
      j1=jklev-1
      j1m1=j1-1
      if(j1m1.lt.1) j1m1=1
      else
      is=-1
      j1=jklev+1
      j1m1=j1-1
      if(j1.gt.nlev) j1=nlev
      endif
    r = rdeno(wz(jlon,jlat,j1), wz(jlon,jlat,j1m1), wz(jlon,jlat,jklev), wz(jlon,jlat,jklev-1))
    b = max(0., min(2., max(r, min(2.*r,1.))))
    zphi = is+zamu*b-is*b
    wfw(jlon,jklev) = 0.5*zamu*((1.+zphi)*wz(jlon,jlat,jklev-1)+(1.-zphi)*wz(jlon,jlat,jklev))
    enddo
    enddo
    do jklev = 1, nlev
    do jlon = 2, nlonm1
    zdv = (s(jlon,jlat,jklev+1)-s(jlon,jlat,jklev))*zcost
    wz(jlon,jlat,jklev) = wz(jlon,jlat,jklev)+wfw(jlon,jklev)-wfw(jlon,jklev+1)+p(jlon,jlat,jklev)*zdv
    enddo
    enddo
    endif

 10 continue

!  Define 2 y-ghostlines of vertically advected variables

#ifdef mpi
      if (ip_s.eq.ip_null) then
        wz(2:nlonm1,1,:) = wz(2:nlonm1,2,:)
        wz(2:nlonm1,0,:) = wz(2:nlonm1,1,:)
      endif
      if (ip_n.eq.ip_null) then
        wz(2:nlonm1,nlat  ,:) = wz(2:nlonm1,nlatm1,:)
        wz(2:nlonm1,nlat+1,:) = wz(2:nlonm1,nlat,:)
      endif
#else
      wz(2:nlonm1,1,:) = wz(2:nlonm1,2,:)
      wz(2:nlonm1,0,:) = wz(2:nlonm1,1,:)
      wz(2:nlonm1,nlat  ,:) = wz(2:nlonm1,nlatm1,:)
      wz(2:nlonm1,nlat+1,:) = wz(2:nlonm1,nlat,:)
#endif

#ifdef mpi
      call u_ghost (wz(2:nlonm1,2     :3     ,:), ip_s, wz(2:nlonm1,nlat  :nlat+1,:), ip_n, 2*(nlon-2)*nlev)
      call u_ghost (wz(2:nlonm1,nlat-2:nlatm1,:), ip_n, wz(2:nlonm1,0     :1     ,:), ip_s, 2*(nlon-2)*nlev)
#endif

!----------------------
!  Meridional advection
!----------------------

    zcosty = dt/dy
    do 20 jklev = 1, nlev
    do 15 jlat = 2, nlat
    do jlon = 2, nlonm1
    zamu = v(jlon,jlat,jklev)*zcosty
      if (zamu.gt.0.) then
      is=1
      j1=jlat-1
      else
      is=-1
      j1=jlat+1
      endif
    r = rdeno (wz(jlon,j1,jklev), wz(jlon,j1-1,jklev), wz(jlon,jlat,jklev), wz(jlon,jlat-1,jklev))
    b = max(0., min(2., max(r, min(2.*r,1.))))
    zphi = is+zamu*b -is*b
    zpby(jlon,jlat) = .5*zamu*((1.+zphi)*wz(jlon,jlat-1,jklev)+(1.-zphi)*wz(jlon,jlat,jklev))
    enddo
 15 continue

    do jlat = 2, nlatm1
    zhxvtn = clv(jlat+1)*fmyu(jlat)
    zhxvt  = clv(jlat  )*fmyu(jlat)
    do jlon = 2, nlonm1
    zdv = (v(jlon,jlat+1,jklev)*zhxvtn -v(jlon,jlat,jklev)*zhxvt)*zcosty
    p0(jlon,jlat,jklev) = wz(jlon,jlat,jklev) +zpby(jlon,jlat)*zhxvt -zpby(jlon,jlat+1)*zhxvtn &
                        + p(jlon,jlat,jklev)*zdv
    enddo
    enddo
 20 continue

!  Define 2 x-ghostlines of vertically and meridionally advected variables

#ifdef mpi
      if (ip_w.eq.ip_null) then
        p0(0,:,:) = p0(2,:,:)
        p0(1,:,:) = p0(2,:,:)
      endif
      if (ip_e.eq.ip_null) then
        p0(nlon  ,:,:) = p0(nlonm1,:,:)
        p0(nlon+1,:,:) = p0(nlonm1,:,:)
      endif
#else
      p0(0,:,:) = p0(2,:,:)
      p0(1,:,:) = p0(2,:,:)
      p0(nlon  ,:,:) = p0(nlonm1,:,:)
      p0(nlon+1,:,:) = p0(nlonm1,:,:)
#endif

#ifdef mpi
      call u_ghost (p0(nlon-2:nlonm1,:,:), ip_e, p0(0     :1     ,:,:), ip_w, 2*nlat*nlev)
      call u_ghost (p0(2     :3     ,:,:), ip_w, p0(nlon  :nlon+1,:,:), ip_e, 2*nlat*nlev)
#endif

!------------------
!  Zonal advection
!------------------

    do 21 jklev = 1, nlev
    do 16 jlat = 2, nlatm1
    zcostx = dt*fmyu(jlat)/dx

    do jlon = 2, nlon
    zamu = u(jlon-1,jlat,jklev)*zcostx
      if (zamu.ge.0.) then
      is=1
      j1=jlon-1
      else
      is=-1
      j1=jlon+1
      endif
    r = rdeno (p0(j1,jlat,jklev), p0(j1-1,jlat,jklev), p0(jlon,jlat,jklev), p0(jlon-1,jlat,jklev))
    b = max(0., min(2., max(r, min(2.*r,1.))))
    zphi = is+zamu*b -is*b
    zpbw(jlon) = .5*zamu*((1.+zphi)*p0(jlon-1,jlat,jklev)+(1.-zphi)*p0(jlon,jlat,jklev))
    enddo

    do jlon = 2, nlonm1
    zdv = (u(jlon,jlat,jklev)-u(jlon-1,jlat,jklev))*zcostx
    p(jlon,jlat,jklev) = p0(jlon,jlat,jklev) +zpbw(jlon) -zpbw(jlon+1) +p(jlon,jlat,jklev)*zdv
    enddo
 16 continue
 21 continue

    return
    end subroutine wafone
!###############################################################################################################
    subroutine comp_esk(esat, qsat, t, p, iflag)

! Computes esat from temperature and qsat from absolute temperature and pressure
! IFLAG 1: esat and qsat with respect to water and ice, separately, depending if t>tzer or t<tzer
!       2: esat and qsat with an interpolation at t<tzer between water and ice
!       3: esat and qsat with respect to water also for t<tzer

    use mod_moloch, only : tzer, ezer, cpv, cw, rd, rv, yliv, yliw, eps, ci, ylwv, ccw1, ccw2, cci1, cci2
    implicit none
    real  esat, qsat, t, p, zt0t, zesk, zratio
    integer iflag

    zt0t = tzer/t
    if (zt0t.le.1.) then
    zesk = ezer*exp(-ccw1*alog(zt0t)+ccw2*(1.-zt0t)) ! Partial pressure over water
    else
      if (iflag.eq.1) then
      zesk = ezer*exp(-cci1*alog(zt0t)+cci2*(1.-zt0t)) ! Partial pressure over ice
      elseif (iflag.eq.2) then
      zratio = 1.04979*(0.5 + 0.5*tanh((t-tzer+9.)/6.))
      zesk = zratio*(ezer*exp(-ccw1*alog(zt0t)+ccw2*(1.-zt0t)))+      &
            (1.-zratio)*(ezer*exp(-cci1*alog(zt0t)+cci2*(1.-zt0t)))
      elseif (iflag.eq.3) then
      zesk = ezer*exp(-ccw1*alog(zt0t)+ccw2*(1.-zt0t))
      else
      print*, "Iflag out of range in subr. comp_esk", iflag
      stop
      endif
    endif

    esat=zesk
    qsat = zesk*eps/(p+zesk*(eps-1.))
    return
    end subroutine comp_esk
#ifdef rad_ecmwf
!###############################################################################################################
   subroutine aerdef (nyrc, nmonc, ndayc, nhouc, nminc)

! aerosol and ozone definition - MOLOCH version - called at initial time and at long intervals

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

use parkind1  ,only : jpim     ,jprb
use yomlun   , only : nulnam
use yomcst   , only : rd       ,rg       ,rtt      ,rsigma   ,       &
                      rcpd     ,rpi      ,rday     ,rcpd     ,       &
                      ri0      ,rsigma   ,                           &
                      rea      ,rlvtt    ,rlstt    ,rmd      ,rkbol, &
                      rnavo    ,r        ,repsm
use yoeaerd  , only : cvdaes   ,cvdael   ,cvdaeu   ,cvdaed   ,          &
                      rcaeops  ,rcaeopl  ,rcaeopu  ,rcaeopd  ,rctrbga  ,&
                      rcvobga  ,rcstbga  ,rctrpt   ,rcaeadm  ,rcaeros  ,&
                      rcaeadk
use yoelw    , only : nsil     ,nipd     ,ntra     ,nua      ,&
                      ng1      ,ng1p1    ,wg1
use yoephli  , only : lphylin
use yoerad   , only : naer     ,nmode    ,nozocl   ,&
                      nradfr   ,nradpfr  ,nradpla  ,nrint    ,&
                      novlp    ,nlw      ,nsw      ,&
                      ntsw     ,leradhs  ,lhvolca  ,lnewaer  ,&
                      lonewsw  ,niceopt  ,nliqopt  ,lrrtm    ,lsrtm  ,&
                      nradip   ,nradlp   ,ninhom   ,nmcica
use yoerdi   , only : rcardi   ,rch4     ,rn2o     ,ro3      ,&
                      rcfc11   ,rcfc12   ,repclc
use yoerdu   , only : nuaer    ,ntraer   ,nimp     ,nout     ,&
                      r10e     ,replog   ,repsc    ,repsco   ,&
                      repscq   ,repsct   ,repscw   ,diff
use yomrip   , only : nindat   ,nsssss   ,nstadd   ,nstass             ,&
                      rtimst   ,rstati   ,rtimtr   ,rhgmt    ,reqtim   ,&
                      rsovr    ,rdeaso   ,rdecli   ,rwsovr   ,rip0     ,&
                      rcodec   ,rsidec   ,rcovsr   ,rsivsr   ,rdtsa    ,&
                      rdtsa2   ,rdts62   ,rdts22   ,rtdt     ,rtmolt   ,&
                      rsivsrlu ,rcovsrlu ,rdeclu   ,rsideclu ,rcodeclu
use yomphy3  , only : rii0
use yoethf   , only : r2es     ,r3les    ,r3ies    ,r4les              ,&
                      r4ies    ,r5les    ,r5ies    ,r5alvcp  ,r5alscp  ,&
                      ralvdcp  ,ralsdcp  ,rtwat    ,rtice    ,rticecu  ,&
                      rtwat_rtice_r      ,rtwat_rticecu_r
use yoedbug  , only : nstpdbg, kstpdbg

use mod_moloch, only : nlon, nlat, nlev, nlevp1, gnlon, gnlat, gfield, dz, h, ps, p, t,  &
                       tskin, phig, aerosol, ozon, myid, alont, alatt, g, tzer, aerotot
implicit none

integer jlon, jlat, jklev, nyrc, nmonc, ndayc, nhouc, nminc
integer klon, klev, kaer, kaero, ksw, iiyr, ierr, jf

parameter (klon=nlon, klev=nlev, kaer=6, kaero=1)

!- reference arrays

real(kind=jprb) :: paer5(klon,kaer,klev), paero5(klon,klev,kaero)
real(kind=jprb) :: paph5(klon,klev+1), pap5(klon,klev)
real(kind=jprb) :: pgelam5(klon),pgemu5(klon), pslon5(klon),pclon5(klon)
real(kind=jprb) :: pozon5(klon,klev)
real(kind=jprb) :: pth5(klon,klev+1) , pt5(klon,klev), pts5(klon)

!- extra arrays

real(kind=jprb) :: plat5(klon), plon5(klon), zozon5(klon,klev)
real(kind=jprb) :: zeta(klev), zetah(klev+1)

integer(kind=jpim) :: klw, kmode
integer(kind=jpim) :: icldmix, iflag, igeo, ilay, iprint ,      &
                      kulout, ndump, nflux, nrad, nprint,       &
                      month, iday, iyr, ilwrad, iminut, inhomf, &
                      ilcloud, io3only, icewat
integer(kind=jpim) :: jl, jk, jnu, jaer, ja, jsw, jcld, jiclq

logical :: lcloud, lforcld, ltimtst, ldebug
logical :: levoigt, llgeose, llgeosw, llgms, llindsa, llmto

real(kind=jprb) :: zalb(klon), zemis, zemiw, rgammas
real(kind=jprb) :: prii05, pcco25, zepaer,  zthetoz, zangozc, zdegrad,       &
                   ztt5, zcondw, zfract, zzmu0, zpp5, zrr5, zddqq, zee, zqq5
real(kind=jprb) :: zcoolc, zcoolt, zheatc, zheatt
real(kind=jprb) :: zxkgkg, zfacti, zfactl
real(kind=jprb) :: ptarg
real(kind=jprb) :: zdpgcp
real(kind=jprb) :: zdflwc, zdflwt, zdfswc, zdfswt, zrcday, zsigma
real(kind=jprb) :: zbla, zlatp, zlonp, zphour, zpsur, zteta, ztim, ztsur
real(kind=jprb) :: zairmwg, zco2mwg, zch4mwg, zn2omwg, zno2mwg, zo3mwg &
                 , zc11mwg, zc12mwg, zc22mwg, zcl4mwg
real(kind=jprb) :: zco2, zch4, zn2o, zo3, zno2, zc11, zc12, zc22, zcl4

!  load external functions

#include "fctast.h"
#include "fcttim.h"
#include "fcttre.h"

if (myid == 0) print*, 'Definition of aerosol and ozone'

!  basic constants

rpi   = 3.1415926535898
rtt   = tzer
rea   = 149597870000.   ! Earth-Sun mean distance
rtwat = rtt
rtice = rtt-23.
rday  = 86400.
rg    = g
rmd   = 28.9644        ! molecular weight of air
rkbol = 1.380658e-23   ! Boltzmann constant
rnavo = 6.0221367e+23  ! Avogadro number
r     = rnavo*rkbol
rd    = 1000.*r/rmd
rcpd  = 3.5*rd
ri0   = 1365.          ! Solar constant
rsigma= 5.67e-08       ! Stefan's constant

rgammas= 0.08
icldmix= 3
lcloud =.true.
lforcld=.false.
ltimtst=.false.
lhvolca=.false.
lnewaer=.true.      ! If true, aerosol monthly distributions are used
nradip = 3          ! Index for diagnosis of ice cloud effective radius
nradlp = 2          ! Index for diagnosis of liq. cloud effective radius
lrrtm  =.false.
lsrtm  =.false.
ninhom = 0          ! 0 if no inhomogeneity scaling effect
lphylin=.false.
levoigt=.false.
leradhs=.true.      ! .t. if rad. is computed on a coarser sampled grid (unused)
lonewsw=.true.      ! .t. if new sw code is active
nstpdbg= 0
kstpdbg(:)=-999
kmode=0
zdegrad = rpi/180.

call su_mcica

!------------------------------------------------------------------
! Define the configuration
!------------------------------------------------------------------

lrrtm=.true. ! .t. if rrtm140mr is used for lw radiation transfer
lsrtm=.true.
ndump=0
ldebug=.false.

klw   =  16
nlw   =  klw ! number of longwave spectral intervals
ksw   =  14
nsw   =  ksw ! number of shortwave spectral intervals
ntsw  =  14  ! maximum possible number of sw spectral intervals
ninhom = 0   ! 0 if no inhomogeneity scaling effect

nradip = 3   ! index for diagnosis of ice cloud effective radius
nradlp = 2   ! index for diagnosis of liq. cloud effective radius
nliqopt= 2   ! index for liquid water cloud optical properties
niceopt= 3   ! index for ice cloud optical properties
naer   = 1   ! configuration index for aerosols (if =1, climatological database used)
nflux  = 6
nmode  = 0   ! configuration for radiation code: flux vs. radiance
nrad   = 1
nradfr = -3  ! frequency of full radiation computations (every '-nradfr' hours)
nradpfr= 36  ! print frequency for rad. statistics (in rad. t. steps)
nradpla= 15  ! print rad. statistics every 'nradpla' rows
nrint  = 4   ! interpolation distance (in points)
nuaer  = 31  ! number of absorber amounts w or w/o aerosols
ntraer = 19  ! number of transmission functions w or w/o aerosols
novlp  = 1   ! cloud overlap configuration
nprint = 1

zairmwg = 28.970_jprb
zco2mwg = 44.011_jprb
zch4mwg = 16.043_jprb
zn2omwg = 44.013_jprb
zno2mwg = 46.006_jprb
zo3mwg  = 47.9982_jprb
zc11mwg = 137.3686_jprb
zc12mwg = 120.9140_jprb
zc22mwg =  86.4690_jprb
zcl4mwg = 153.8230_jprb

zch4 = 1.72e-06_jprb*zch4mwg/zairmwg
zn2o = 310.e-09_jprb*zn2omwg/zairmwg
zno2 = 500.e-13_jprb*zno2mwg/zairmwg
zo3  =   1.e-06_jprb*zo3mwg /zairmwg
zc11 = 280.e-12_jprb*zc11mwg/zairmwg
zc12 = 484.e-12_jprb*zc12mwg/zairmwg
zc22 =   1.e-12_jprb*zc22mwg/zairmwg
zcl4 =   1.e-12_jprb*zcl4mwg/zairmwg

igeo=0
llgeose =.false.
llgeosw =.false.
llgms   =.false.
llindsa =.false.
llmto   =.false.

repsc  = 1.e-12   ! sec. epsilon for cloud cover
repsco = 1.e-12   ! sec. epsilon for ozone amount
repscq = 1.e-12   ! sec. epsilon for water vapor
repsct = 1.e-12   ! sec. epsilon for shortwave optical thickness
repscw = 1.e-12   ! sec. epsilon for cloud liquid water
replog = 1.e-12   ! sec. epsilon for abs.amount in laplace transform
nout   = 6        ! unit number for the extra prints

!------------------------------------------------------------------

call surdi
call sulwn
call suaerl
call suaerh

call surrtab
call surrtpk
call surrtrf
call surrtftr

call rrtm_kgb1
call rrtm_kgb2
call rrtm_kgb3
call rrtm_kgb4
call rrtm_kgb5
call rrtm_kgb6
call rrtm_kgb7
call rrtm_kgb8
call rrtm_kgb9
call rrtm_kgb10
call rrtm_kgb11
call rrtm_kgb12
call rrtm_kgb13
call rrtm_kgb14
call rrtm_kgb15
call rrtm_kgb16

call rrtm_init_140gp

call suswn   (ntsw, nsw)
call suclopn (ntsw, nsw, nlev)
call suaersn  (ntsw, nsw)

! Routines specific to SRTM

if (lsrtm) then
call srtm_init
call susrtaer
call susrtcop
else
call suswn   (ntsw,nsw)
call suaersn (ntsw,nsw)
endif

! nindat:  YYYYMMDD date of the simulation, for example 19870701

  iiyr = nyrc
  nindat =  iiyr*10000 + nmonc*100 + ndayc

! nsss: no. of seconds in a day (rounded to one minute)

  nsssss = nhouc*3600 + nminc*60

! Basic constants
! 0 as last arg. of sucts: no printouts; 1: prints all constants

  kulout = 6  ! logical unit for output of sucst
  call sucst (kulout, nindat, nsssss, 0)  ! defines costants needed by rtime - also rday defined above!

! rtimtr: time in sec. after 12:00 gmt of 1 jan 2000

  rtimtr = rtime (iiyr, nmonc, ndayc, nsssss, rday) ! function in /include/fcttim.h

  rhgmt = nhouc*3600. + nminc*60. ! time of the day in sec rounded to min.
  zteta = rteta(rtimtr)
  rdeaso= rrs(zteta)
  rdecli= rds(zteta)
  reqtim= ret(zteta)
  rsovr = reqtim + rhgmt
  rwsovr= rsovr*2.0*rpi/rday
  rip0  = ri0*rea*rea/(rdeaso*rdeaso)
  rii0  = rip0
  rcodec= cos(rdecli)
  rsidec= sin(rdecli)
  rcovsr= cos(rwsovr)
  rsivsr= sin(rwsovr)

  nozocl = 2        ! nozocl=2: Fortuin-Langematz O3 climatology + aerosol
  io3only= 0
  zrcday = rday * rg / rcpd
  diff   = 1.66
  r10e   = 0.4342945
  prii05 = rii0
  zepaer = 1.e-12

  iminut=int(float(nsssss)/60.)

!- Fortuin-Langematz O3 climatology

  call suecozc ( nindat , iminut )

!- ECMWF Geleyn O3 climatology

  zthetoz=rteta(rtimtr)
  zangozc=rel(zthetoz)-1.7535
  call suecozo ( zangozc )

! climatological aerosol

  if (naer.eq.1) then
  call suecaebc                  ! black carbon aerosol (fires)
  call suecaeor                  ! organic type aerosol
  call suecaesd                  ! soil-dust aerosol
  call suecaess                  ! sea-salt aerosol
  call suecaesu                  ! sulfate-type aerosol
  call suecaec (nindat, iminut)  ! climatol. distrib. of volcanic aerosol
  endif

! zeta and zetah: sigma-like coordinate for the routines that set parameters for aerosol

    do jk = 1, nlev
    zeta(jk)  = (jk-.5)/float(nlev)
    zetah(jk) = (jk-1.)/float(nlev)
    enddo
    zetah(nlevp1) = 1.

! Latitudes and longitudes in radians

    zdegrad = rpi/180.

! Definition of cloud parameters (in module yoecld, used by radlswr, called in radintec)

    call sucld (nlev, zeta)

!-------------------------------------------------------------------
! Definition of atmospheric variables

    do jlat = 1, nlat

    do jlon = 1, nlon
    plat5(jlon)  = alatt(jlon,jlat) * zdegrad
    plon5(jlon)  = alont(jlon,jlat) * zdegrad  ! long. < 0 not accepted
    if (plon5(jlon).lt.0.) plon5(jlon) = plon5(jlon) + 2.*rpi
    pgelam5(jlon)= plon5(jlon)
    pgemu5(jlon) = sin(plat5(jlon))
    pclon5(jlon) = cos(plon5(jlon))
    pslon5(jlon) = sin(plon5(jlon))
    enddo

! Pressure at full levels (filtered in the vertical)

    do jlon = 1, nlon
    pap5(jlon,nlev) = p(jlon,jlat,1   )
    pap5(jlon,1   ) = p(jlon,jlat,nlev)
    enddo
    do jklev = 2, nlev-1
    jk = nlevp1-jklev
    do jlon = 1, nlon
    pap5(jlon,jk) = .25*(p(jlon,jlat,jklev-1)+p(jlon,jlat,jklev+1))+.5*p(jlon,jlat,jklev)
    enddo
    enddo

! Pressure at half levels

    do jlon = 1, nlon
    paph5(jlon,1)      = 0.
    paph5(jlon,nlevp1) = ps(jlon,jlat)
    enddo
    do jk = 2, nlev
    do jlon = 1, nlon
    paph5(jlon,jk) = 0.5*(pap5(jlon,jk-1)+pap5(jlon,jk))
    enddo
    enddo

! Temperature at full levels

    do jklev = 1, nlev
    jk = nlevp1-jklev
    do jlon = 1, nlon
    pt5(jlon,jk) = t(jlon,jlat,jklev)
    enddo
    enddo

! Temperature at half-levels: at top defined as lev. 1
! at bottom as tskin, elsewhere is the arithmetic mean

    do jlon = 1, nlon
    pth5(jlon,1     ) = pt5(jlon,1)
    pth5(jlon,nlevp1) = tskin(jlon,jlat)
    pts5(jlon)        = tskin(jlon,jlat)
    enddo

    do jk = 2, nlev
    do jlon = 1, nlon
    pth5(jlon,jk) = 0.5*(pt5(jlon,jk-1)+pt5(jlon,jk))
    enddo
    enddo

! Call to other set-up routines for various coefficients. These are dependent on the vertical resolution

    call suaerv (nlev, zetah, cvdaes, cvdael, cvdaeu, cvdaed, rctrbga, rcvobga, rcstbga, rcaeops, rcaeopl, &
                 rcaeopu, rcaeopd, rctrpt, rcaeadk, rcaeadm, rcaeros)

! Derive the aerosols and ozone distribution from climatology

    call radaca (1 , nlon , nlon , 1 , nlev ,                                        &
                  paph5 , pgelam5 , pgemu5 , pclon5 , pslon5 , pth5 , paer5 , zozon5)

! The computation of ozone must be done separately at each point in longitude, because radocz
! assumes (incorrectly for rotated grid) that latitude is the same for all vectors in longitude

    do jlon= 1, nlon
    call radozc (1, 1, 1, 1, nlev, 1, 1, 0, paph5(jlon,:), pgemu5(jlon), zozon5(jlon,:))
    enddo
    pozon5 = zozon5

    if (naer.eq.0) paer5 = zepaer

    do jk = 1, nlev
    do jlon = 1, nlon
    ozon   (jlon,jlat,jk  ) = pozon5(jlon,  jk)     ! ozone
    aerosol(jlon,jlat,jk,1) = paer5 (jlon,1,jk)     ! land (organic + sulfate) aerosol
    aerosol(jlon,jlat,jk,2) = paer5 (jlon,2,jk)     ! sea salt aerosol
    aerosol(jlon,jlat,jk,3) = paer5 (jlon,3,jk)     ! desert dust aerosol
    aerosol(jlon,jlat,jk,4) = paer5 (jlon,4,jk)     ! urban + black carbon aerosol
    aerosol(jlon,jlat,jk,5) = paer5 (jlon,5,jk)     ! volcanic aerosol
    aerosol(jlon,jlat,jk,6) = paer5 (jlon,6,jk)     ! stratospheric background aerosol
    enddo
    enddo

! Ad hoc increase of some aerosol (urban and sulfate) over the Po Valley in the lower troposphere

    do jk = 1, nlev
    do jlon = 1, nlon
    if(alatt(jlon,jlat).gt.44.3.and.alatt(jlon,jlat).lt.46.2.and.   &
       alont(jlon,jlat).gt. 7.0.and.alont(jlon,jlat).lt.13.4.and.   &
       phig(jlon,jlat)/g.lt.500..and.pap5(jlon,jk).gt.800.e2) then
    aerosol(jlon,jlat,jk,1) = 1.10*aerosol(jlon,jlat,jk,1)
    aerosol(jlon,jlat,jk,4) = 1.20*aerosol(jlon,jlat,jk,4)
    endif
    enddo
    enddo

    enddo ! jlat

    do jf = 1, 4
    call filt3d (ozon, 1.)
    enddo
    ozon = max (ozon, 0.)
    do jf = 1, 20
    do jk = 1, 6
    call filt3d (aerosol(1:nlon,1:nlat,1:nlev,jk), 1.)
    enddo
    enddo
    do jf = 1, 90
    call filt3d (aerosol(1:nlon,1:nlat,1:nlev, 2), 1.)
    enddo
    aerosol = max (aerosol, zepaer)

    aerotot(:,:,:) = aerosol(:,:,:,1) + aerosol(:,:,:,2) + aerosol(:,:,:,3) +  &
                     aerosol(:,:,:,4) + aerosol(:,:,:,5) + aerosol(:,:,:,6)

    goto 345  ! skip plots for verification

    call collect (ozon(1:nlon,1:nlat,3), gfield)
    if (myid == 0) call plotout (gfield, gnlon, gnlat, 99)

    call collect (aerosol(1:nlon,1:nlat,nlev-5 ,1), gfield)
    if (myid == 0) call plotout (gfield, gnlon, gnlat, 99)
    call collect (aerosol(1:nlon,1:nlat,nlev-5 ,2), gfield)
    if (myid == 0) call plotout (gfield, gnlon, gnlat, 99)
    call collect (aerosol(1:nlon,1:nlat,nlev-5 ,3), gfield)
    if (myid == 0) call plotout (gfield, gnlon, gnlat, 99)
    call collect (aerosol(1:nlon,1:nlat,nlev-5 ,4), gfield)
    if (myid == 0) call plotout (gfield, gnlon, gnlat, 99)
    call collect (aerosol(1:nlon,1:nlat,5 ,5), gfield)
    if (myid == 0) call plotout (gfield, gnlon, gnlat, 99)
    call collect (aerosol(1:nlon,1:nlat,5 ,6), gfield)
    if (myid == 0) call plotout (gfield, gnlon, gnlat, 99)
    call collect (aerotot(1:nlon,1:nlat,nlev-5), gfield)
    if (myid == 0) call plotout (gfield, gnlon, gnlat, 99)

    call collect (aerosol(1:nlon,1:nlat,nlev ,1), gfield)
    if (myid == 0) call plotout (gfield, gnlon, gnlat, 99)
    call collect (aerosol(1:nlon,1:nlat,nlev ,2), gfield)
    if (myid == 0) call plotout (gfield, gnlon, gnlat, 99)
    call collect (aerosol(1:nlon,1:nlat,nlev ,3), gfield)
    if (myid == 0) call plotout (gfield, gnlon, gnlat, 99)
    call collect (aerosol(1:nlon,1:nlat,nlev ,4), gfield)
    if (myid == 0) call plotout (gfield, gnlon, gnlat, 99)
    call collect (aerotot(1:nlon,1:nlat,nlev), gfield)
    if (myid == 0) call plotout (gfield, gnlon, gnlat, 99)

#ifdef mpi
      call mpi_finalize (ierr)
#endif
    stop
345 continue

    return
    end subroutine aerdef
#endif
!###############################################################################################################
#ifdef rad_ecmwf
    subroutine radintec (jlat, jl1, nyrc, nmonc, ndayc, nhouc, nminc)

! Interface subroutine between MOLOCH and ECMWF radiation version 35R2
! (received by J.J. Morcrette, Feb. 2012)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Input variables to radlswr:
! nlonr:   no. of final (and total) points in longitude (defined in bolam:
!          normally radiation is computed here only at alternate points)
! kmode:   control param. (unused)
! kaer:    no. of aerosol types (set to 6) - it must be set to the same value here and in aerdef
! kaero:   last dim. of paero (progn. aerosol) - set to 1 - must be set to the same value here and in aerdef
! prii05:  effective solar input, taking into account orbital annual modulation
! paer5:   optical depth of the different aerosol types
! paero5:  prognostic aerosol (unused)
! palbd5:  albedo for the different spectral intervals of diffuse short wave (visible)
!          it must include the snow effect
! palbp5:  as palbd5, but for direct (parallel) short wave
! paph5:   pressure at semi-integer levels
! pap5:    pressure at integer levels
! pccnl5:  cloud condens. nuclei over land (use depending on lccnl, undefined, so false)
! pccno5:  as above but over ocean
! pgelam5: longitude in radians
! pgemu5:  sin of latitude
! pco25:   concentration of CO2 (pa/pa) - defined in bolam
! pch45:   concentration of minor gases - defined below
! pn2o5:        "
! pno25:        "
! pc115:        "
! pc125:        "
! pc225:        "
! pcl45:        "
! pclfr5:  cloud fraction
! pdp5:    pressure thickness of each integer layer
! pemis5:  surface emissivity for long waves, except in the 8-12.5 micron window
! pemiw5:  surface emissivity for long waves within the 8-12.5 micron window
! plsm5:   land (+ sea ice) to sea fraction (i.e. solid ground - 1: all land or sea ice; 0: all liquid sea)
! pmu05:   cosine of the solar zenith angle
! pozon5:  concentration of ozone (pa/pa)
! pq5:     specific humidity
! pqiwp5:  cloud ice (kg/kg)
! pqlwp5:  cloud water (kg/kg)
! pqs5:    saturation specific humidity kg/kg (not used in radlswr)
! pqrain5: rain water (not used in radlswr)
! praint5: rain rate (not used in radlswr)
! pth5:    temperature of semi-integer levels - valuest at top and bottom not obvious
! pt5:     temperature of integer levels
! pts5:    surface (skin) temperature
! pnbas5:  index of base of convective layer (not used in radlswr)
! pntop5:  index of top of convective layer (not used in radlswr)

! Output variables from radlswr:

! pemit5:  surface total longwave emissivity (not used)
! pfct5:   clear-sky lw net fluxes
! pflt5:   total lw net fluxes
! pfcs5:   clear-sky sw net fluxes
! pfls5:   total sw net fluxes
! pfrsod5: total-sky surface sw downward flux
! psudu5:  solar radiance in sun's direction
! puvdf5:  surface downward u.v. radiation (not used)
! pparf:   photosynthetically active radiation (not used)
! pparcf5: clear-sky photosynthetically active radiation (not used)
! ptincf5: top-of-atmosphere incident solar radiation (not used)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

use parkind1  ,only : jpim     ,jprb
use yomlun   , only : nulnam

use yomcst   , only : rd       ,rg       ,rtt      ,rsigma   ,       &
                      rcpd     ,rpi      ,rday     ,rcpd     ,       &
                      ri0      ,rsigma   ,                           &
                      rea      ,rlvtt    ,rlstt    ,rmd      ,rkbol, &
                      rnavo    ,r        ,repsm
use yoeaerd  , only : cvdaes   ,cvdael   ,cvdaeu   ,cvdaed   ,          &
                      rcaeops  ,rcaeopl  ,rcaeopu  ,rcaeopd  ,rctrbga  ,&
                      rcvobga  ,rcstbga  ,rctrpt   ,rcaeadm  ,rcaeros  ,&
                      rcaeadk
use yoelw    , only : nsil     ,nipd     ,ntra     ,nua      ,&
                      ng1      ,ng1p1    ,wg1
use yoephli  , only : lphylin
use yoerad   , only : naer     ,nmode    ,nozocl   ,&
                      nradfr   ,nradpfr  ,nradpla  ,nrint    ,&
                      novlp    ,nlw      ,nsw      ,&
                      ntsw     ,leradhs  ,lhvolca  ,lnewaer  ,&
                      lonewsw  ,niceopt  ,nliqopt  ,lrrtm    ,lsrtm  ,&
                      nradip   ,nradlp   ,ninhom   ,nmcica
use yoerdi   , only : rcardi   ,rch4     ,rn2o     ,ro3      ,&
                      rcfc11   ,rcfc12   ,repclc
use yoerdu   , only : nuaer    ,ntraer   ,nimp     ,nout     ,&
                      r10e     ,replog   ,repsc    ,repsco   ,&
                      repscq   ,repsct   ,repscw   ,diff
use yomrip   , only : nindat   ,nsssss   ,nstadd   ,nstass             ,&
                      rtimst   ,rstati   ,rtimtr   ,rhgmt    ,reqtim   ,&
                      rsovr    ,rdeaso   ,rdecli   ,rwsovr   ,rip0     ,&
                      rcodec   ,rsidec   ,rcovsr   ,rsivsr   ,rdtsa    ,&
                      rdtsa2   ,rdts62   ,rdts22   ,rtdt     ,rtmolt   ,&
                      rsivsrlu ,rcovsrlu ,rdeclu   ,rsideclu ,rcodeclu
use yomphy3  , only : rii0
use yoethf   , only : r2es     ,r3les    ,r3ies    ,r4les              ,&
                      r4ies    ,r5les    ,r5ies    ,r5alvcp  ,r5alscp  ,&
                      ralvdcp  ,ralsdcp  ,rtwat    ,rtice    ,rticecu  ,&
                      rtwat_rtice_r      ,rtwat_rticecu_r
use yoedbug  , only : nstpdbg, kstpdbg

use mod_moloch, only : nlon, nlat, nlonr, nradm, nlev, nlevp1, dz, h, g, cpd, cpv, ps, p, swsdtf, &
                       t, tskin, aerosol, ozon, myid, corvis, corirr, corrdt, tzer, swsddr,       &
                       q, qcw, qci, qpw, qpi1, fcloud, ccw1, ccw2, cci1, cci2, eps, ezer, emisg1, &
                       emisg2, albedo, fsnow, fmask, fice, alsn, g, co2ppm, mcica,           &
                       alont, alatt, cloudt, slopeff

implicit none

integer jlon, jlat, jklev, nyrc, nmonc, ndayc, nhouc, nminc
integer klon, klev, kaer, kaero, ksw, iiyr, ierr, jl, jl1
real zalsn, zfsnow, zdstcp, zomd

parameter (klon=nlonr, klev=nlev, kaer=6, kaero=1, zomd=cpv/cpd-1.)

!- reference arrays

real(kind=jprb) :: paer5(klon,kaer,klev), paero5(klon,klev,kaero)
real(kind=jprb) :: palbd5(klon,14)   , palbp5(klon,14)
real(kind=jprb) :: paph5(klon,klev+1), pap5(klon,klev)
real(kind=jprb) :: pclfr5(klon,klev) , pdp5(klon,klev)
real(kind=jprb) :: pemis5(klon), pemiw5(klon), plsm5(klon), pmu05(klon)
real(kind=jprb) :: pgelam5(klon),pgemu5(klon), pslon5(klon),pclon5(klon)
real(kind=jprb) :: pozon5(klon,klev)
real(kind=jprb) :: pq5(klon,klev)
real(kind=jprb) :: pqiwp5(klon,klev) , pqlwp5(klon,klev)
real(kind=jprb) :: pqs5(klon,klev)
real(kind=jprb) :: pqrain5(klon,klev), praint5(klon,klev)
real(kind=jprb) :: pccnl5(klon)      , pccno5(klon)
real(kind=jprb) :: pth5(klon,klev+1) , pt5(klon,klev), pts5(klon)
real(kind=jprb) :: pco25(klon,klev)  , pch45(klon,klev)  , pn2o5(klon,klev), &
                   pno25(klon,klev)  , pc115(klon,klev)  , pc125(klon,klev), &
                   pc225(klon,klev)  , pcl45(klon,klev)
real(kind=jprb) :: pemit5(klon)
real(kind=jprb) :: pfct5(klon,klev+1), pflt5(klon,klev+1)
real(kind=jprb) :: pfcs5(klon,klev+1), pfls5(klon,klev+1)
real(kind=jprb) :: pfrsod5(klon)       , psudu5(klon)
real(kind=jprb) :: pnbas5(klon)        , pntop5(klon)
real(kind=jprb) :: puvdf5(klon)        , ptincf5(klon)
real(kind=jprb) :: pparf5(klon)        , pparcf5(klon)

!- extra arrays

real(kind=jprb) :: plat5(klon), plon5(klon), zozon5(klon,klev)
real(kind=jprb) :: zeta(klev), zetah(klev+1)

integer(kind=jpim) :: klw, kmode
integer(kind=jpim) :: icldmix, iflag, igeo, ilay, iprint ,      &
                      kulout, ndump, nflux, nrad, nprint,       &
                      month, iday, iyr, ilwrad, iminut, inhomf, &
                      ilcloud, io3only, icewat
integer(kind=jpim) :: jk, jnu, jaer, ja, jsw, jcld, jiclq

logical :: lcloud, lforcld, ltimtst, ldebug
logical :: levoigt, llgeose, llgeosw, llgms, llindsa, llmto

real(kind=jprb) :: zalb(klon), zemis, zemiw, rgammas
real(kind=jprb) :: prii05, pcco25, zepaer,  zthetoz, zangozc, zdegrad,      &
                   ztt5, zcondw, zfract, zzmu0, zpp5, zrr5, zddqq, zee, zqq5
real(kind=jprb) :: zcoolc, zcoolt, zheatc, zheatt
real(kind=jprb) :: zxkgkg, zfacti, zfactl
real(kind=jprb) :: ptarg
real(kind=jprb) :: zdpgcp
real(kind=jprb) :: zdflwc, zdflwt, zdfswc, zdfswt, zrcday, zsigma
real(kind=jprb) :: zbla, zlatp, zlonp, zphour, zpsur, zteta, ztim, ztsur
real(kind=jprb) :: zairmwg, zco2mwg, zch4mwg, zn2omwg, zno2mwg, zo3mwg, &
                   zc11mwg, zc12mwg, zc22mwg, zcl4mwg
real(kind=jprb) :: zco2, zch4, zn2o, zo3, zno2, zc11, zc12, zc22, zcl4

! load external functions

#include "fctast.h"
#include "fcttim.h"
#include "fcttre.h"

!  basic constants

rpi    = 3.1415926535898
rtt    = tzer
rea    = 149597870000. ! Earth-Sun mean distance
rtwat  = rtt
rtice  = rtt-23.
rday   = 86400.
rg     = g
rmd    = 28.9644       ! molecular weight of air
rkbol  = 1.380658e-23  ! Boltzmann constant
rnavo  = 6.0221367e+23 ! Avogadro number
r      = rnavo*rkbol
rd     = 1000.*r/rmd
rcpd   = 3.5*rd
ri0    = 1365.         ! Solar constant
rsigma = 5.67e-08      ! Stefan's constant

rgammas= 0.08
icldmix= 3
lcloud =.true.
lforcld=.false.
ltimtst=.false.
lhvolca=.false.
lnewaer=.true.       ! If true, aerosol monthly distributions are used
nradip = 3           ! Index for diagnosis of ice cloud effective radius
nradlp = 2           ! Index for diagnosis of liq. cloud effective radius
lrrtm  =.false.
lsrtm  =.false.
ninhom = 0           ! 0 if no inhomogeneity scaling effect
lphylin=.false.
levoigt=.false.
leradhs=.false.      ! not used
lonewsw=.true.       ! .t. if new sw code is active
nstpdbg= 0
kstpdbg(:)=-999
kmode= 0
zdegrad = rpi/180.

call su_mcica

!------------------------------------------------------------------
! Define the configuration
!------------------------------------------------------------------

lrrtm=.true.  ! .t. if rrtm140mr is used for lw radiation transfer
lsrtm=.true.
ndump= 0
ldebug=.false.

nmcica = mcica ! if =2, use of McICA method for clouds - mcica defined in the main program

klw    =  16
nlw    =  klw ! number of longwave spectral intervals
ksw    =  14
nsw    =  ksw ! number of shortwave spectral intervals
ntsw   =  14  ! maximum possible number of sw spectral intervals
ninhom =  0   ! 0 if no inhomogeneity scaling effect

nradip = 3    ! index for diagnosis of ice cloud effective radius
nradlp = 2    ! index for diagnosis of liq. cloud effective radius
nliqopt= 2    ! index for liquid water cloud optical properties
niceopt= 3    ! index for ice cloud optical properties
naer   = 1    ! configuration index for aerosols
nflux  = 6
nmode  = 0    ! configuration for radiation code: flux vs. radiance
nrad   = 1
nradfr =-3    ! frequency of full radiation computations (every '-nradfr' hours)
nradpfr= 36   ! print frequency for rad. statistics (in rad. t. steps)
nradpla= 15   ! print rad. statistics every 'nradpla' rows
nrint  =  4   ! interpolation distance (in points)
nuaer  = 31   ! number of absorber amounts w or w/o aerosols
ntraer = 19   ! number of transmission functions w or w/o aerosols
novlp  =  1   ! cloud overlap configuration
nprint =  1

zairmwg =  28.970_jprb
zco2mwg =  44.011_jprb
zch4mwg =  16.043_jprb
zn2omwg =  44.013_jprb
zno2mwg =  46.006_jprb
zo3mwg  =  47.9982_jprb
zc11mwg = 137.3686_jprb
zc12mwg = 120.9140_jprb
zc22mwg =  86.4690_jprb
zcl4mwg = 153.8230_jprb

zch4 = 1.72e-06_jprb*zch4mwg/zairmwg
zn2o = 310.e-09_jprb*zn2omwg/zairmwg
zno2 = 500.e-13_jprb*zno2mwg/zairmwg
zo3  =   1.e-06_jprb*zo3mwg /zairmwg
zc11 = 280.e-12_jprb*zc11mwg/zairmwg
zc12 = 484.e-12_jprb*zc12mwg/zairmwg
zc22 =   1.e-12_jprb*zc22mwg/zairmwg
zcl4 =   1.e-12_jprb*zcl4mwg/zairmwg

igeo=0
llgeose =.false.
llgeosw =.false.
llgms   =.false.
llindsa =.false.
llmto   =.false.

repsc  = 1.e-12    ! epsilon for cloud cover
repsco = 1.e-12    ! epsilon for ozone amount
repscq = 1.e-12    ! epsilon for water vapor
repsct = 1.e-12    ! epsilon for shortwave optical thickness
repscw = 1.e-12    ! epsilon for cloud liquid water
replog = 1.e-12    ! epsilon for abs. amount in laplace transform
nout   = 6         ! unit number for the extra prints

ptincf5 = 0.

!------------------------------------------------------------------

call surdi
call sulwn
call suaerl
call suaerh

call surrtab
call surrtpk
call surrtrf
call surrtftr

call rrtm_kgb1
call rrtm_kgb2
call rrtm_kgb3
call rrtm_kgb4
call rrtm_kgb5
call rrtm_kgb6
call rrtm_kgb7
call rrtm_kgb8
call rrtm_kgb9
call rrtm_kgb10
call rrtm_kgb11
call rrtm_kgb12
call rrtm_kgb13
call rrtm_kgb14
call rrtm_kgb15
call rrtm_kgb16

call rrtm_init_140gp

call suswn   (ntsw, nsw)
call suclopn (ntsw, nsw, nlev)
call suaersn (ntsw, nsw)

!-- routines specific to SRTM

if (lsrtm) then
call srtm_init
call susrtaer
call susrtcop
else
call suswn   (ntsw,nsw)
call suaersn (ntsw,nsw)
endif

! nindat:  YYYYMMDD date of the simulation, for example 19870701

  iiyr = nyrc
  nindat =  iiyr*10000 + nmonc*100 + ndayc

! nsssss: no. of seconds of a day (rounded to minute)

  nsssss = nhouc*3600 + nminc*60

! Basic constants
! 0 as last arg. of sucts: no printouts; 1: prints all constants

  kulout = 6  ! logical unit for output of sucst
  call sucst (kulout, nindat, nsssss, 0)  ! defines costants needed by rtime - also rday defined above!

! rtimtr: time in sec. after 12:00 gmt of 1 jan 2000

  rtimtr = rtime (iiyr, nmonc, ndayc, nsssss, rday) ! function in /include/fcttim.h

  rhgmt = nhouc*3600. + nminc*60.   ! time of the day in sec rounded to min.
  zteta=rteta(rtimtr)               ! function in /include/fctast.h
  rdeaso=rrs(zteta)                 ! (also other functions below)
  rdecli=rds(zteta)
  reqtim=ret(zteta)
  rsovr =reqtim + rhgmt
  rwsovr=rsovr*2.0*rpi/rday
  rip0=ri0*rea*rea/(rdeaso*rdeaso)
  rii0 = rip0
  rcodec=cos(rdecli)
  rsidec=sin(rdecli)
  rcovsr=cos(rwsovr)
  rsivsr=sin(rwsovr)

  diff   = 1.66
  r10e   = 0.4342945
  prii05 = rii0
  zco2   = co2ppm*1.e-06*zco2mwg/zairmwg ! CO2 concentration
  pcco25 = zco2

! Definition of atmospheric variables
! Latitudes and longitudes in radians

  zdegrad = rpi/180.

  do jl = 1, nlonr
  jlon = jl1 + (jl-1)*nradm
  plat5(jl)  = alatt(jlon,jlat) * zdegrad
  plon5(jl)  = alont(jlon,jlat) * zdegrad  ! long. < 0 not accepted
  if (plon5(jl).lt.0.) plon5(jl) = plon5(jl) + 2.*rpi
  pgelam5(jl)= plon5(jl)
  pgemu5(jl) = sin(plat5(jl))
  pclon5(jl) = cos(plon5(jl))
  pslon5(jl) = sin(plon5(jl))
  pmu05(jl)  = max(0._jprb, rsidec*pgemu5(jl)-rcodec*sqrt(1.-pgemu5(jl)**2)*cos(pgelam5(jl)+rwsovr))
  enddo

! Pressure at full levels (filtered in the vertical)

  do jl = 1, nlonr
  jlon = jl1 + (jl-1)*nradm
  pap5(jl,nlev) = p(jlon,jlat,1   )
  pap5(jl,1   ) = p(jlon,jlat,nlev)
  enddo
  do jklev = 2, nlev-1
  jk = nlevp1-jklev
  do jl = 1, nlonr
  jlon = jl1 + (jl-1)*nradm
  pap5(jl,jk) = .25*(p(jlon,jlat,jklev-1)+p(jlon,jlat,jklev+1))+.5*p(jlon,jlat,jklev)
  enddo
  enddo

! Pressure at half levels

  do jl = 1, nlonr
  jlon = jl1 + (jl-1)*nradm
  paph5(jl,1)      = 0.
  paph5(jl,nlevp1) = ps(jlon,jlat)
  enddo
  do jk = 2, nlev
  do jl = 1, nlonr
  paph5(jl,jk) = 0.5*(pap5(jl,jk-1)+pap5(jl,jk))
  enddo
  enddo

! Temperature at full levels

  do jklev = 1, nlev
  jk = nlevp1-jklev
  do jl = 1, nlonr
  jlon = jl1 + (jl-1)*nradm
  pt5(jl,jk) = t(jlon,jlat,jklev)
  enddo
  enddo

! Temperature at half-levels: at top defined as lev. 1
! at bottom as t at 2m or tskin, elsewhere is the arithmetic mean

  do jl = 1, nlonr
  jlon = jl1 + (jl-1)*nradm
  pth5(jl,1     ) = pt5(jl,1)
  pth5(jl,nlevp1) = tskin(jlon,jlat)
  pts5(jl)        = tskin(jlon,jlat)
  enddo
  do jk = 2, nlev
  do jl = 1, nlonr
  pth5(jl,jk) = 0.5*(pt5(jl,jk-1)+pt5(jl,jk))
  enddo
  enddo

  do jk = 1, nlev
  do jl = 1, nlonr
  jlon = jl1 + (jl-1)*nradm
  pozon5(jl,  jk) = ozon   (jlon,jlat,jk  )
  paer5 (jl,1,jk) = aerosol(jlon,jlat,jk,1)
  paer5 (jl,2,jk) = aerosol(jlon,jlat,jk,2)
  paer5 (jl,3,jk) = aerosol(jlon,jlat,jk,3)
  paer5 (jl,4,jk) = aerosol(jlon,jlat,jk,4)
  paer5 (jl,5,jk) = aerosol(jlon,jlat,jk,5)
  paer5 (jl,6,jk) = aerosol(jlon,jlat,jk,6)
  enddo
  enddo

!--------------------------------------------------------------------------------

  do jklev = 1, nlev
  jk = nlevp1-jklev
  do jl = 1, nlonr
  jlon = jl1 + (jl-1)*nradm
  pq5(jl,jk)     = q(jlon,jlat,jklev)
  pclfr5(jl,jk)  = fcloud(jlon,jlat,jklev)
  pqiwp5(jl,jk)  = qci(jlon,jlat,jklev) + 0.2*qpi1(jlon,jlat,jklev)   ! to be verified
  pqlwp5(jl,jk)  = qcw(jlon,jlat,jklev)
  pdp5  (jl,jk)  = max(paph5(jl,jk+1)-paph5(jl,jk),10.d0) ! allows negative thickness!

!  The LW and SW radiation schemes can handle profiles of the trace gases.
!  Here fixed concentrations are used.

  pco25(jl,jk) = zco2
  pch45(jl,jk) = zch4
  pn2o5(jl,jk) = zn2o
  pno25(jl,jk) = zno2
  pc115(jl,jk) = zc11
  pc125(jl,jk) = zc12
  pc225(jl,jk) = zc22
  pcl45(jl,jk) = zcl4
  enddo
  enddo

  do jl = 1, nlonr
  jlon = jl1 + (jl-1)*nradm
  zfsnow = fsnow(jlon,jlat)

! Change of snow albedo in a range depending on skin temperature

!  if (tskin(jlon,jlat).lt.277.) then
!  zalsn = alsn-.1 + .2*(277.-tskin(jlon,jlat))/14.
!  zalsn = min (zalsn, alsn+.1)
!  else
!  zalsn = alsn-.1
!  endif

!  zalb  (jl) = albedo(jlon,jlat)*(1.-zfsnow) + zalsn*zfsnow
  zalb  (jl) = albedo(jlon,jlat)*(1.-zfsnow) + alsn(jlon,jlat)*zfsnow
  pemis5(jl) = emisg1(jlon,jlat)*(1.-zfsnow) + 0.99*zfsnow
  pemiw5(jl) = emisg2(jlon,jlat)*(1.-zfsnow) + 0.99*zfsnow
  plsm5 (jl) = 1. - fmask(jlon,jlat) + fice(jlon,jlat)*fmask(jlon,jlat)
  enddo

! Albedo: for parallel (direct) radiation, dependency on solar zenith angle is computed over land
! (Yang et al) and over sea (Taylor et al, as ECMWF)

  do jnu = 1, ksw
  do jl = 1, nlonr
  jlon = jl1 + (jl-1)*nradm
  palbd5(jl,jnu) = zalb(jl)      ! albedo for diffuse light
  if (fmask(jlon,jlat).lt.0.5) then
  palbp5(jl,jnu) = min(1., zalb(jl)*(1. + 0.26)/(1. + 2.*0.26*pmu05(jl)))
  else
  palbp5(jl,jnu) = (1.-plsm5(jl))*0.037/(1.1*pmu05(jl)**1.4 + 0.15) + zalb(jl)*plsm5(jl)
  endif
  enddo
  enddo

  call radlswr ( 1, nlonr, nlonr, nlev, kmode, kaer, kaero                     &
               , prii05                                                        &
               , paer5 , paero5, palbd5, palbp5, paph5 , pap5                  &
               , pccnl5, pccno5, pgelam5,pgemu5                                &
               , pco25 , pch45 , pn2o5 , pno25 , pc115  , pc125 , pc225, pcl45 &
               , pclfr5, pdp5  , pemis5, pemiw5, plsm5  , pmu05 , pozon5       &
               , pq5   , pqiwp5, pqlwp5, pqs5  , pqrain5, praint5              &
               , pth5  , pt5   , pts5  , pnbas5, pntop5                        &
               , pemit5, pfct5 , pflt5 , pfcs5 , pfls5  , pfrsod5              &
               , psudu5, puvdf5, pparf5, pparcf5,ptincf5 )

  do jk = 1, nlev
  jklev = nlevp1 - jk
  do jl = 1, nlonr
  jlon = jl1 + (jl-1)*nradm
  zdstcp = g/(cpd*(1.+zomd*q(jlon,jlat,jklev))*(paph5(jl,jk+1)-paph5(jl,jk)))
  corrdt(jlon,jlat,jklev) = zdstcp*(pflt5(jl,jk)-pflt5(jl,jk+1)+pfls5(jl,jk)-pfls5(jl,jk+1))
  enddo
  enddo

! Surface fluxes of visible and infrared radiation (positive downward)

  do jl = 1, nlonr
  jlon = jl1 + (jl-1)*nradm
!  corvis(jlon,jlat) = pfls5(jl,nlevp1)  ! total net sw radiation on a flat surface
  swsdtf(jlon,jlat) = pfrsod5(jl)       ! total sky sw radiation
  swsddr(jlon,jlat) = psudu5(jl)        ! parallel solar sw radiation
  corvis(jlon,jlat) = (swsdtf(jlon,jlat) - swsddr(jlon,jlat))*(1.-palbd5(jl,1)) + &
                       swsddr(jlon,jlat)*(1.-palbp5(jl,1))*slopeff(jlon,jlat) ! total net sw radiation on topog. slopes
!  corirr(jlon,jlat) = pflt5(jl,nlevp1)
  corirr(jlon,jlat) = max(pflt5(jl,nlevp1), 0.92*(1.-0.05*cloudt(jlon,jlat))*pflt5(jl,nlevp1)) ! tuning!
  enddo

  return
  end subroutine radintec
#endif
!###############################################################################################################
      subroutine calendar (nyrin, nmonin, ndayin, nhouin, nminin, iday,        &
                           ihou, imin, nyrc, nmonc, ndayc, nhouc, nminc, ndayr)

!  Defines current calendar date and time (nyrc, nmonc, ndayc, nhouc, nminc)
!  and "julian" day of the year (ndayr) by adding the forecast validity time
!  (iday, ihou, imin) to the forecast initial date and time
!  (nyrin, nmonin, ndayin, nhouin, nminin).

      implicit none
      integer imonth(12), nyrin, nmonin, ndayin, nhouin, nminin, iday, &
       ihou, imin, nyrc, nmonc, ndayc, nhouc, nminc, ndayr, j, itday,  &
       isum

!  define day of the year ( 1 < ndayr < 366 ) of initial date

      imonth( 1) = 31
      imonth( 2) = 28
      imonth( 3) = 31
      imonth( 4) = 30
      imonth( 5) = 31
      imonth( 6) = 30
      imonth( 7) = 31
      imonth( 8) = 31
      imonth( 9) = 30
      imonth(10) = 31
      imonth(11) = 30
      imonth(12) = 31

      if ( mod(nyrin,400).eq.0 .or. (mod(nyrin,4).eq.0.and.mod(nyrin,100).ne.0) ) imonth(2) = 29

      if(ndayin.gt.imonth(nmonin)) then
      write(*, '(a,i3,a,i3,a,i5)') " This day does not exist: day =", &
            ndayin, ", month =", nmonin, ", year =", nyrin
      stop " Stop."
      endif

      ndayr = ndayin
      do j = 1, nmonin-1
      ndayr = ndayr + imonth(j)
      enddo

!  update initial hour and minutes

      nminc = nminin + imin
      nhouc = nhouin + ihou
      ndayr = ndayr  + iday
      if (nminc.ge.60) then
      nhouc = nhouc + nminc/60
      nminc = mod(nminc,60)
      endif
      if (nhouc.ge.24) then
      ndayr = ndayr + nhouc/24
      nhouc = mod(nhouc, 24)
      endif

!  update ndayr and initial day, month, year

      nyrc = nyrin
 1    itday  = 365
      if ( mod(nyrc,400).eq.0 .or. (mod(nyrc,4).eq.0.and.mod(nyrc,100).ne.0) ) itday = 366
      if (ndayr.gt.itday) then
      ndayr = ndayr - itday
      nyrc  = nyrc + 1
      else
      go to 2
      endif
      go to 1

 2    imonth(2) = 28
      if ( mod(nyrc,400).eq.0 .or. (mod(nyrc,4).eq.0.and.mod(nyrc,100).ne.0) ) imonth(2) = 29
      isum = 0
      do nmonc = 1, 12
      isum = isum + imonth(nmonc)
      if (ndayr.le.isum) go to 3
      enddo

 3    ndayc = ndayr + imonth(nmonc) - isum

      return
      end subroutine calendar
!###############################################################################################################
subroutine sea_surface

!  Sea surface scheme: definition of tskin and qskin over sea surface

 use mod_moloch, only : tg, tskin, qskin, nlon, nlat, nlevg, q, hflux, qflux, frvis, frirr, fsnow,  &
                        dtstep, ccw1, ccw2, cci1, cci2, tzer, ezer, ps, eps, roscdt, nlonm1, nlatm1, &
                        snow, fmask, fice, rgmd
implicit none

 integer :: jlon, jlat
 real :: zt0t, zesk, zeskw, zeski, zqsat, zeffw, zrc, ztotfl

 do jlat = 2, nlatm1
 do jlon = 2, nlonm1

   if (fmask(jlon,jlat).ge.0.5.and.fice(jlon,jlat).lt.0.8) then  ! over sea (including thin sea ice)
!   if (mask_soil(jlon,jlat) == 0) then  ! over sea (excluding thick sea ice)

     zt0t = tzer/tg(jlon,jlat,1)
     if (tg(jlon,jlat,1).ge.tzer) then
       zesk = ezer*exp(-ccw1*log(zt0t)+ccw2*(1.-zt0t))       ! partial pressure over water
     else
       zeskw = ezer*exp(-ccw1*log(zt0t)+ccw2*(1.-zt0t))      ! partial pressure over water
       zeski = ezer*exp(-cci1*log(zt0t)+cci2*(1.-zt0t))      ! partial pressure over ice
       zesk  = zeski*fice(jlon,jlat) + zeskw*(1.-fice(jlon,jlat))
     endif
     zqsat = zesk*eps/(ps(jlon,jlat)+zesk*(eps-1.))
     zeffw = .058/(.058 + roscdt(jlon,jlat) + 1.e-9)
     qskin(jlon,jlat) = zeffw*zqsat*0.97 + (1.-zeffw)*q(jlon,jlat,1)

     ztotfl = hflux(jlon,jlat) +2.5008e6*qflux(jlon,jlat) -frvis(jlon,jlat) -frirr(jlon,jlat)

     zrc = 1./(1. + 1.5e2*rgmd(jlon,jlat))  ! reduces the net heat flux effect with rough sea...
     if (fice(jlon,jlat).le.0.5) then
       tg(jlon,jlat,1) = tg(jlon,jlat,1) - dtstep*(.25e-7*zrc*ztotfl + .8e-5*(tg(jlon,jlat,1)-tg(jlon,jlat,nlevg)))
     else             !  reduced heat capacity and heat diffusivity over sea-ice
       tg(jlon,jlat,1) = tg(jlon,jlat,1) - dtstep*( .48e-7*ztotfl + .5e-5* (tg(jlon,jlat,1)-tg(jlon,jlat,nlevg)) )
     endif
     if(fice(jlon,jlat).lt.0.5) then
       snow(jlon,jlat)  = 0.
       fsnow(jlon,jlat) = 0.
     endif
     tskin(jlon,jlat) = tg(jlon,jlat,1)

   endif

 enddo
 enddo

return
end subroutine sea_surface
!###############################################################################################################
    subroutine paidef (p, t, q, qcw, qci, pai)

    use mod_moloch, only : nlon, nlat, nlev, pzer, rd, cpd, eps, fmzh, dz, g
    implicit none
    real(4), dimension(nlon,nlat,nlev) :: p, t, q, qcw, qci, pai, tvirt
    real zeps, zb, zdelta
    integer jlon, jlat, jklev

!--------------------------------------------------------------------------------------------------------
!  Hydrostatic initialization of PAI
!--------------------------------------------------------------------------------------------------------

    pai(:,:,1) = (p(:,:,1)/pzer)**(rd/cpd)
    zeps   = 1./eps
    tvirt = t*(1. -(q+qcw+qci) +zeps*q)

    do jklev = 2, nlev
    do jlat = 1, nlat
    do jlon = 1, nlon
    zb = 2.*g/fmzh(jlon,jlat,jklev)*dz/cpd+tvirt(jlon,jlat,jklev)-tvirt(jlon,jlat,jklev-1)
    zdelta = (zb**2+4.*tvirt(jlon,jlat,jklev-1)*tvirt(jlon,jlat,jklev))**.5
    pai(jlon,jlat,jklev) = -pai(jlon,jlat,jklev-1)/(2.*tvirt(jlon,jlat,jklev-1))*(zb-zdelta)
    enddo
    enddo
    enddo

    return
    end subroutine paidef
!###############################################################################################################
 subroutine interp(alfa, ex1, ex2, npi, xi, g, x, f, nval)

!  Interpolates with splines with tension in one dimension.
!  The spline is defined imposing that the second derivative is the average
!  of second derivatives computed at the two adjacent points.
!  At interval extremes the second derivative is assumed null.
!  This subr. also extrapolates out of the interval where the input funtion g is defined

!  Input:  function g defined at coordinates xi (caution: can be changed by this subr.)
!          g(1:npi) values at irregular but strictly growing coordinates xi(1:npi)
!  Output: f(1:nval) interpolated values at arbitrary coordinates x(1:nval)

!  alfa: spline tension parameter, comprised between 0 and 1:
!  if alfa=1, pure linear interpolation; if alfa=0, pure spline

!  ex1: param. determining extrapolation for x < xi(1)
!  ex2: param. determining extrapolation for x > xi(npi)
!  if ex1=0 or ex2=0, constant value extrapolation is used at corresponding extreme
!  if ex1=1 or ex2=1, linear extrapolation is used at corresponding extreme
!  intermediate values of ex1 and ex2 give intermediate extrapolation values

  real, dimension(npi )  :: xi, g
  real, dimension(nval)  :: x,  f

  if(alfa.lt..0.or.alfa.gt.1) then
  print*, 'Caution: in interp, alfa out of interval 0-1'
  endif
  if(ex1.lt..0.or.ex1.gt.1) then
  print*, 'Caution: in interp, ex1 out of interval 0-1'
  endif
  if(ex2.lt..0.or.ex2.gt.1) then
  print*, 'Caution: in interp, ex2 out of interval 0-1'
  endif

! Fix for the case in which coordinates of the input function are not strictly increasing
! Note that this changes the input coordinates in the calling programme

  do k=2,npi
   if(xi(k).le.xi(k-1)) then
   print*, "Caution: in interp, coordinates of input function changed because not monotonic!"
   exit
   endif
  enddo

  zeps=(xi(npi)-xi(1))*1.e-6   ! small deviation used to set apart interlaced coordinates
200 do k=2,npi
     if(xi(k).le.xi(k-1)) then
     ximed=0.5*(xi(k)+xi(k-1))
     xi(k-1)=ximed-zeps
     xi(k)=ximed+zeps
     gmed=0.5*(g(k)+g(k-1))
     g(k-1)=gmed
     g(k)=gmed
     endif
    enddo

 do k=2,npi
  if(xi(k).le.xi(k-1)) then
  goto 200
  endif
 enddo

 do 100 jval =1, nval

!  2 cases of extrapolation

 if(x(jval).lt.xi(1)) then
 f(jval) = g(1) + ex1*(g(1)-g(2))/(xi(1)-xi(2)) * (x(jval)-xi(1))
 go to 100
 elseif (x(jval).ge.xi(npi)) then
 f(jval) = g(npi) + ex2*(g(npi)-g(npi-1))/(xi(npi)-xi(npi-1)) * (x(jval)-xi(npi))
 go to 100
 endif

 ir = 0

!  ir is a reference index determining the interpolation interval
!  The interpolation expression is applied also if x=xi(j)

 do j = 1, npi
 if (x(jval).ge.xi(j)) ir = ir + 1
 enddo

 if (ir.eq.1) then
 fmm = 2*g(1) - g(2)
 xmm = 2*xi(1) - xi(2)
 fpp = g(ir+2)
 xpp = xi(ir+2)
 elseif (ir.eq.(npi-1)) then
 fpp = 2*g(npi) - g(npi-1)
 xpp = 2*xi(npi) - xi(npi-1)
 fmm = g(ir-1)
 xmm = xi(ir-1)
 else
 fmm = g(ir-1)
 xmm = xi(ir-1)
 fpp = g(ir+2)
 xpp = xi(ir+2)
 endif

 fm     = g(ir)
 xm     = xi(ir)
 fp     = g(ir+1)
 xp     = xi(ir+1)
 delx   = xp - xm
 delxp  = xpp - xp
 delxm  = xm - xmm
 delx1  = x(jval) - xm
 delx2  = xp - x(jval)
 delxs  = delx**2
 delx1s = delx1**2
 delx2s = delx2**2

!  Spline contribution to interpolation

 spl = fm*(delx2/delx + delx1*delx2s/(delxs*delxm) - delx1s*     &
       delx2/((delx+delxp)*delxs)) + fp*(delx1/delx +            &
       delx1s*delx2/(delxs*delxp) - delx1*delx2s/((delx+delxm)*  &
       delxs)) - fmm * delx1*delx2s/((delx+delxm)*delx*delxm) -  &
       fpp * delx1s*delx2/((delx+delxp)*delx*delxp)

!  Linear interpolation contribution

 clin = (fm*delx2 + fp*delx1)/delx

!  Final interpolation combined using alfa

 f(jval) = alfa*clin + (1.-alfa)*spl

 100  continue

 return
 end subroutine interp
!###############################################################################################################
    subroutine plotout (a, b, n, m, title, nhf)

! Creates a file containing 2D fields (matrix a) for plotting.
! Matrix must contain the fmask to be plotted superimposed on each field.
! Title of the graph and dimension of the matrix are written first.

   real a(n,m)
   real b(n,m)
   character*30 title

   write(nhf,*) title
   write (nhf,*) n, m
   do j=1,m
   write (nhf,37) (b(i,j),i=1,n)
   enddo
   do j=1,m
   write (nhf,38) (a(i,j),i=1,n)
   enddo

37 format (20f6.3)
38 format (10e12.5)

   return
end subroutine plotout
!###############################################################################################################
      subroutine hdiffu

!  Horizontal turbulent diffusion of u and v

      use mod_moloch, only : u, v, ip_n, ip_s, ip_e, ip_w, nlon, nlat, nlev, nlonm1, nlatm1, ntop, dx, dy, fmyu, fmyv, &
                             dtstep, chm, tke, ip_null
      implicit none
      real, dimension(nlon,nlat) :: zux, zvy, zuyvx
      real zam, zdu, zdv, zup, zum, zvp, zvm
      integer jlon, jlat, jklev

#ifdef mpi
      call u_ghost (u(nlonm1,:,:), ip_e, u(1,   :,:), ip_w, nlat*nlev)
      call u_ghost (u(:,nlatm1,:), ip_n, u(:   ,1,:), ip_s, nlon*nlev)
      call u_ghost (v(2,:     ,:), ip_w, v(nlon,:,:), ip_e, nlat*nlev)
      call u_ghost (v(:,2,     :), ip_s, v(:,nlat,:), ip_n, nlon*nlev)
#endif

      do jklev = 1, ntop

      do jlat = 2, nlatm1
      do jlon = 2, nlonm1
      zux(jlon,jlat) = (u(jlon,jlat,jklev)-u(jlon-1,jlat,jklev))*fmyu(jlat)/dx
      zvy(jlon,jlat) = (v(jlon,jlat+1,jklev)-v(jlon,jlat,jklev))/dy
      zuyvx(jlon,jlat) = (u(jlon,jlat,jklev)-u(jlon,jlat-1,jklev))/dy +         &
                         (v(jlon+1,jlat,jklev)-v(jlon,jlat,jklev))*fmyv(jlat)/dx

      zam = chm(jlon,jlat,jklev)
      zam = min (zam, .05*dx**2/dtstep)  !!  Flux limitation

      zux(jlon,jlat) = -2.*zam*zux(jlon,jlat) +(tke(jlon,jlat,jklev)+tke(jlon,jlat,jklev+1))/3. ! flussi di momento
      zvy(jlon,jlat) = -2.*zam*zvy(jlon,jlat) +(tke(jlon,jlat,jklev)+tke(jlon,jlat,jklev+1))/3. ! flussi di momento
      zuyvx(jlon,jlat) = -zam*zuyvx(jlon,jlat)
      enddo
      enddo

#ifdef mpi
      call u_ghost (zux  (2,     :), ip_w, zux  (nlon,:), ip_e, nlat)
      call u_ghost (zvy  (:,nlatm1), ip_n, zvy  (:,   1), ip_s, nlon)
      call u_ghost (zuyvx(nlonm1,:), ip_e, zuyvx(1,   :), ip_w, nlat)
      call u_ghost (zuyvx(:,     2), ip_s, zuyvx(:,nlat), ip_n, nlon)
#endif

      if (ip_w.eq.ip_null) zuyvx(1,:   ) = zuyvx(2,:     )
      if (ip_n.eq.ip_null) zuyvx(:,nlat) = zuyvx(:,nlatm1)
      if (ip_s.eq.ip_null) zvy  (:,1   ) = zvy  (:,2     )
      if (ip_e.eq.ip_null) zux  (nlon,:) = zux  (nlonm1,:)

      do jlat = 2, nlatm1
      do jlon = 2, nlonm1
      zdu = (zux(jlon+1,jlat)-zux(jlon,jlat))*fmyu(jlat)/dx + (zuyvx(jlon,jlat+1)-zuyvx(jlon,jlat))/dy
      u(jlon,jlat,jklev) = u(jlon,jlat,jklev) - dtstep*zdu
      zdv = (zuyvx(jlon,jlat)-zuyvx(jlon-1,jlat))*fmyv(jlat)/dx + (zvy(jlon,jlat)-zvy(jlon,jlat-1))/dy
      v(jlon,jlat,jklev) = v(jlon,jlat,jklev) - dtstep*zdv
      enddo
      enddo

      enddo

      return
      end subroutine hdiffu
!###############################################################################################################
      subroutine hdifft

!  Horizontal turbulent diffusion of scalars on T points

      use mod_moloch, only : ip_n, ip_s, ip_e, ip_w, nlon, nlat, nlev, nlonm1, nlatm1, ntop, dx, dy, fmyu, &
                             dtstep, tetav, q, fmz, dz, hx, hy, gzita, clv, rd, p, tvirt, chm, ip_null
      implicit none
      real, dimension(nlon,nlat,nlev) :: zftx, zfty, zfqx, zfqy
      real ztx, zty, zqx, zqy, zro, zrom1, rdz2, zgz, zgzp, zgzm
      real zah, zdt, zdq, dz2, zita, zrfmzu, zrfmzum, zrfmzvp, zrfmzv, zup, zum, zvp, zvm
      integer jlon, jlat, jklev, jklevm1

#ifdef mpi
      call u_ghost (tetav(2,:     ,:), ip_w, tetav(nlon,:,:), ip_e, nlat*nlev)
      call u_ghost (tetav(:,nlatm1,:), ip_n, tetav(:   ,1,:), ip_s, nlon*nlev)
      call u_ghost (q    (2,:     ,:), ip_w, q    (nlon,:,:), ip_e, nlat*nlev)
      call u_ghost (q    (:,nlatm1,:), ip_n, q    (:   ,1,:), ip_s, nlon*nlev)
#endif

!  Computation of horizontal fluxes

      do jklev = 1, ntop+1
      jklevm1 = max (1, jklev-1)
      rdz2 = 1./(dz*(jklev+1-jklevm1))
      zita = dz/2.+(jklev-1)*dz
      zgz  = gzita(zita)
      do jlat = 2, nlatm1
      do jlon = 2, nlonm1

      zah   = chm(jlon,jlat,jklev)/.75 ! Prandtl
      zah   = min (zah, .05*dx**2/dtstep)  !!  Flux limitation

      zro   = p(jlon,jlat,jklev)/(rd*tvirt(jlon,jlat,jklev))

      ztx = (tetav(jlon+1,jlat,jklev)-tetav(jlon,jlat,jklev))*fmyu(jlat)/dx
      ztx = ztx - zgz*hx(jlon,jlat)*.25*(fmz(jlon,jlat,jklev)+fmz(jlon+1,jlat,jklev))* &
             (tetav(jlon+1,jlat,jklev+1)+tetav(jlon,jlat,jklev+1)-tetav(jlon+1,jlat,jklevm1)-tetav(jlon,jlat,jklevm1))*rdz2

      zty = (tetav(jlon,jlat,jklev)-tetav(jlon,jlat-1,jklev))/dy
      zty = zty - zgz*hy(jlon,jlat)*.25*(fmz(jlon,jlat,jklev)+fmz(jlon,jlat-1,jklev))* &
             (tetav(jlon,jlat,jklev+1)+tetav(jlon,jlat-1,jklev+1)-tetav(jlon,jlat,jklevm1)-tetav(jlon,jlat-1,jklevm1))*rdz2

      zqx = (q(jlon+1,jlat,jklev)-q(jlon,jlat,jklev))*fmyu(jlat)/dx
      zqx = zqx - zgz*hx(jlon,jlat)*.25*(fmz(jlon,jlat,jklev)+fmz(jlon+1,jlat,jklev))* &
             (q(jlon+1,jlat,jklev+1)+q(jlon,jlat,jklev+1)-q(jlon+1,jlat,jklevm1)-q(jlon,jlat,jklevm1))*rdz2

      zqy = (q(jlon,jlat,jklev)-q(jlon,jlat-1,jklev))/dy
      zqy = zqy - zgz*hy(jlon,jlat)*.25*(fmz(jlon,jlat,jklev)+fmz(jlon,jlat-1,jklev))* &
             (q(jlon,jlat,jklev+1)+q(jlon,jlat-1,jklev+1)-q(jlon,jlat,jklevm1)-q(jlon,jlat-1,jklevm1))*rdz2

      zftx(jlon,jlat,jklev) = -zro*zah*ztx
      zfty(jlon,jlat,jklev) = -zro*zah*zty
      zfqx(jlon,jlat,jklev) = -zro*zah*zqx
      zfqy(jlon,jlat,jklev) = -zro*zah*zqy

      enddo
      enddo
      enddo

!  Divergence of horizontal fluxes

#ifdef mpi
      call u_ghost (zftx(nlonm1,:,:), ip_e, zftx(1,   :,:), ip_w, nlat*nlev)
      call u_ghost (zfty(:,     2,:), ip_s, zfty(:,nlat,:), ip_n, nlon*nlev)
      call u_ghost (zfqx(nlonm1,:,:), ip_e, zfqx(1,   :,:), ip_w, nlat*nlev)
      call u_ghost (zfqy(:,     2,:), ip_s, zfqy(:,nlat,:), ip_n, nlon*nlev)
#endif

      if (ip_w.eq.ip_null) then
      zftx(1,:,:) = zftx(2,:,:)
      zfqx(1,:,:) = zfqx(2,:,:)
      endif
      if (ip_n.eq.ip_null) then
      zfty(:,nlat,:) = zfty(:,nlatm1,:)
      zfqy(:,nlat,:) = zfqy(:,nlatm1,:)
      endif

      do jklev = 1, ntop
      jklevm1 = max (1, jklev-1)
      rdz2 = 1./(2.*dz)
      zita = dz/2.+(jklev-1)*dz
      zgzp = gzita(zita+dz)
      if (jklev.eq.1) then
      zgzm = 0.
      else
      zgzm = gzita(zita-dz)
      endif
      do jlat = 2, nlatm1
      do jlon = 2, nlonm1

      zrom1 = rd*tvirt(jlon,jlat,jklev)/p(jlon,jlat,jklev)
      zrfmzu  = 2./(fmz(jlon,jlat,jklev)+fmz(jlon+1,jlat,jklev))
      zrfmzum = 2./(fmz(jlon,jlat,jklev)+fmz(jlon-1,jlat,jklev))
      zrfmzvp = 2./(fmz(jlon,jlat,jklev)+fmz(jlon,jlat+1,jklev))
      zrfmzv  = 2./(fmz(jlon,jlat,jklev)+fmz(jlon,jlat-1,jklev))

      zup = zftx(jlon  ,jlat,jklev)*zrfmzu
      zum = zftx(jlon-1,jlat,jklev)*zrfmzum
      zvp = zfty(jlon,jlat+1,jklev)*zrfmzvp*clv(jlat+1)
      zvm = zfty(jlon,jlat  ,jklev)*zrfmzv *clv(jlat)
      zdt = fmz(jlon,jlat,jklev)*fmyu(jlat)*((zup-zum)/dx +(zvp-zvm)/dy)

      zdt = zdt - fmz(jlon,jlat,jklev)*rdz2*.5* &
                 (zgzp*(hx(jlon,jlat)*zftx(jlon,jlat,jklev+1)+ hx(jlon-1,jlat)*zftx(jlon-1,jlat,jklev+1)) - &
                  zgzm*(hx(jlon,jlat)*zftx(jlon,jlat,jklevm1)+ hx(jlon-1,jlat)*zftx(jlon-1,jlat,jklevm1)) + &
                  zgzp*(hy(jlon,jlat)*zfty(jlon,jlat,jklev+1)+ hy(jlon,jlat+1)*zfty(jlon,jlat+1,jklev+1)) - &
                  zgzm*(hy(jlon,jlat)*zfty(jlon,jlat,jklevm1)+ hy(jlon,jlat+1)*zfty(jlon,jlat+1,jklevm1)) )

      tetav(jlon,jlat,jklev) = tetav(jlon,jlat,jklev) - dtstep*zrom1*zdt

      zup = zfqx(jlon  ,jlat,jklev)*zrfmzu
      zum = zfqx(jlon-1,jlat,jklev)*zrfmzum
      zvp = zfqy(jlon,jlat+1,jklev)*zrfmzvp*clv(jlat+1)
      zvm = zfqy(jlon,jlat  ,jklev)*zrfmzv *clv(jlat)
      zdq = fmz(jlon,jlat,jklev)*fmyu(jlat)*((zup-zum)/dx +(zvp-zvm)/dy)

      zdq = zdq - fmz(jlon,jlat,jklev)*rdz2*.5* &
                 (zgzp*(hx(jlon,jlat)*zfqx(jlon,jlat,jklev+1)+ hx(jlon-1,jlat)*zfqx(jlon-1,jlat,jklev+1)) - &
                  zgzm*(hx(jlon,jlat)*zfqx(jlon,jlat,jklevm1)+ hx(jlon-1,jlat)*zfqx(jlon-1,jlat,jklevm1)) + &
                  zgzp*(hy(jlon,jlat)*zfqy(jlon,jlat,jklev+1)+ hy(jlon,jlat+1)*zfqy(jlon,jlat+1,jklev+1)) - &
                  zgzm*(hy(jlon,jlat)*zfqy(jlon,jlat,jklevm1)+ hy(jlon,jlat+1)*zfqy(jlon,jlat+1,jklevm1)) )

      q(jlon,jlat,jklev) = q(jlon,jlat,jklev) - dtstep*zrom1*zdq

      enddo
      enddo
      enddo

      q = max (q, 1.e-10)

      return
      end subroutine hdifft
!###############################################################################################################
subroutine surfradpar(fmask, fsea_ice, qg_surf, veg_map, soil_qmax, soil_qmin, &
 soil_albedo_dry, soil_albedo_wet, soil_emiss1_dry, soil_emiss1_wet, soil_emiss2_dry, soil_emiss2_wet, &
 veg_frac, veg_albedo, veg_emiss1, veg_emiss2, &
 nlon, nlat, nlevel_soil, nvt, iflag, albedo, emismap1, emismap2)

! Defines radiative parameters of land surface as a function of albedo, emissivity of soil (dry and wet) and vegetation,
! and of fraction of vegetation cover;
! and as a function ao sea ice fraction.

! Input variables:

!      fmask: fraction of land/sea (proportion, 1-sea, 0-land);
!      fsea_ice: fraction of sea ice cover (proportion);
!      qg_surf: specific volumetric soil water content at top soil level (m3/m3);
!      veg_map: vegetation (landuse) types (classification GLC2000, 22 types) defined on the grid points, 3D array,
!               1-st and 2-nd array index are grid point index, 3-rd index is vegetation types number plus 1,
!               if 3-rd index has value 1, then the variable means dominate vegetation type,
!               values 2-23 of 3-rd index means fraction (of area, proportion) of the each of 22 vegetation types;
!      soil_qmax: maximum specific volumetric soil water content (m^3/m^3) at all soil levels, 3D array;
!      soil_qmin: minimum specific volumetric soil water content (m^3/m^3) at all soil levels, 3D array;
!      soil_albedo_dry: dry soil albedo (proportion), 2D array;
!      soil_albedo_wet: wet soil albedo (proportion), 2D array;
!      soil_emiss_1_dry: dry soil emissivity in broadband window (proportion), 2D array;
!      soil_emiss_1_wet: wet soil emissivity in broadband window (proportion), 2D array;
!      soil_emiss_2_dry: dry soil emissivity in 8-12 micron window (proportion), 2D array;
!      soil_emiss_2_wet: wet soil emissivity in 8-12 micron window (proportion), 2D array;
!      veg_frac: vegetation cover fraction, 2D array;
!      veg_albedo: vegetation albedo (proportion), 2D array;
!      veg_emiss_1: vegetaion emissivity in broadband window (proportion), 2D array;
!      veg_emiss_2: vegetation emissivity in 8-12 micron window (proportion), 2D array;
!      nlon: grid point number along axis of rotated longitude;
!      nlat: grid point number along axis of rotated latidude;
!      nlevel_soil: number of soil levels;
!      nvt: total number vegetaion types;
!      iflag: flag for definition of radiative parameters over sea surface; if not 1, 
!             then the parameters not will be defined.

! Output variables:

!      albedo: albedo on the model grid (proportion);
!      emismap1: emissivity in broadband window (proportion);
!      emismap2: emissivity in 8-12 micron window (proportion).

implicit none

! Input

integer :: nlon, nlat, nlevel_soil, nvt, iflag
real, dimension(nlon,nlat) :: fmask, fsea_ice, qg_surf, soil_albedo_dry, soil_albedo_wet, &
 soil_emiss1_dry, soil_emiss1_wet, soil_emiss2_dry, soil_emiss2_wet, &
 veg_frac, veg_albedo, veg_emiss1, veg_emiss2
real, dimension(nlon,nlat,nlevel_soil) :: soil_qmax, soil_qmin
real, dimension(nlon,nlat,nvt+1) :: veg_map

! Output :

real, dimension(nlon,nlat) :: albedo, emismap1, emismap2

! Internal:

integer :: jlon, jlat, nv1, &
 iurbveg = 20  ! vegetation type index for urban area
real :: zvegfr, zsoildf, zsoilwf, &
 albedo_urban=0.6 ! albedo for urban area

! print *
! print *,"Renewal of radiation parameters of the surface"

 do jlat=1,nlat
 do jlon=1,nlon

   nv1 = int(veg_map(jlon,jlat,1)) ! Domain vegetation type at grid point

   zvegfr = veg_frac(jlon,jlat)

   if (fmask(jlon,jlat) < 0.5) then ! land

     zsoilwf = (qg_surf(jlon,jlat)-soil_qmin(jlon,jlat,1))/(soil_qmax(jlon,jlat,1)-soil_qmin(jlon,jlat,1))
     zsoildf = 1.-zsoilwf

! Tuning: increase of vegetation fraction only to determine emissivity
! and albedo (due to previous modification of vegetation fraction)

!!!     if (nv1 < nvt-1) zvegfr = min(1., zvegfr+0.2) ! for Version 19 and earlier

     albedo(jlon,jlat)   = zvegfr*veg_albedo(jlon,jlat)+ &
 (1.-zvegfr)*(zsoildf*soil_albedo_dry(jlon,jlat)+zsoilwf*soil_albedo_wet(jlon,jlat))

     emismap1(jlon,jlat) = zvegfr*veg_emiss1(jlon,jlat)+ &
 (1.-zvegfr)*(zsoildf*soil_emiss1_dry(jlon,jlat)+zsoilwf*soil_emiss1_wet(jlon,jlat))

     emismap2(jlon,jlat) = zvegfr*veg_emiss2(jlon,jlat)+ &
 (1.-zvegfr)*(zsoildf*soil_emiss2_dry(jlon,jlat)+zsoilwf*soil_emiss2_wet(jlon,jlat))

! Urban effect for radiative properties of land surface: albedo is decreased
! proportionally to urban vegetation that defines fraction urban areas

     albedo(jlon,jlat) = albedo(jlon,jlat)*(1.-veg_map(jlon,jlat,iurbveg+1))+  &
                         albedo(jlon,jlat)*veg_map(jlon,jlat,iurbveg+1)*albedo_urban

   else ! sea

! Definition of albedo and emissivity values as a function of land-sea fraction
! (fmask) and ice fraction (fsea_ice).
! 0.07 is the value of water body albedo; 0.40 of sea ice albedo.
! 0.970 is the value of water body emissivity; 0.980 of sea ice emissivity.

     if (iflag.eq.1) then
!!!       albedo(jlon,jlat) = (1.-fsea_ice(jlon,jlat))*0.07 + fsea_ice(jlon,jlat)*0.55
!!! test      emismap1(jlon,jlat) = (1.-fsea_ice(jlon,jlat))*0.990 + fsea_ice(jlon,jlat)*0.997
!!! test      emismap2(jlon,jlat) = (1.-fsea_ice(jlon,jlat))*0.990 + fsea_ice(jlon,jlat)*0.997
       albedo(jlon,jlat) = (1.-fsea_ice(jlon,jlat))*0.07 + fsea_ice(jlon,jlat)*0.40
       emismap1(jlon,jlat) = (1.-fsea_ice(jlon,jlat))*0.970 + fsea_ice(jlon,jlat)*0.980
       emismap2(jlon,jlat) = (1.-fsea_ice(jlon,jlat))*0.970 + fsea_ice(jlon,jlat)*0.980
     endif

   endif

 enddo
 enddo

return
end subroutine surfradpar
!###############################################################################################################
subroutine wrback

!  scrive file per back-trajectories (vento a 10 m e sul secondo model level, pbl)

use mod_moloch
implicit none

integer :: iunit=42, jlon, jlat

    if(myid == 0) then
      open (iunit, file='traject.mhf', form='unformatted', status='unknown', position='append')
      write(iunit) gnlon-2, gnlat-2
      write(iunit) nfdr(5:12)                                           ! data iniziale e validita'
      write(iunit) b0, h, x0, y0, dlon, dlat, pdr(5)+dlon, pdr(4)+dlat  ! coordinate verticali e griglia orizzontale
    endif

      call collect (u_std_lev(1,1,2), gfield)
      if (myid == 0) then
      do jlat = 2, gnlat-1
      write(iunit) (gfield(jlon,jlat), jlon = 2, gnlon-1)
      enddo
      endif
      call collect (v_std_lev(1,1,2), gfield)
      if (myid == 0) then
      do jlat = 2, gnlat-1
      write(iunit) (gfield(jlon,jlat), jlon = 2, gnlon-1)
      enddo
      endif
      call collect (u(1,1,2), gfield)
      if (myid == 0) then
      do jlat = 2, gnlat-1
      write(iunit) (gfield(jlon,jlat), jlon = 2, gnlon-1)
      enddo
      endif
      call collect (v(1,1,2), gfield)
      if (myid == 0) then
      do jlat = 2, gnlat-1
      write(iunit) (gfield(jlon,jlat), jlon = 2, gnlon-1)
      enddo
      endif
      call collect (pbl, gfield)
      if (myid == 0) then
      do jlat = 2, gnlat-1
      write(iunit) (gfield(jlon,jlat), jlon = 2, gnlon-1)
      enddo
      endif

    if (myid == 0) then
      print *,'traject. file written'
      close (iunit)
    endif

return
end subroutine wrback
!#####################################################################################################
     subroutine pbl_height

     use mod_moloch, only : nlon, nlat, nlev, pbl, gust, rich, zeta, u_std_lev, v_std_lev, u, v

! Definition of PBL height and gust based on PBL
 
     pbl = 1.

     do jlat = 1, nlat
     do jlon = 1, nlon

     jkpbl = 0
     do jklev = 1, nlev/4
     zriav = .5*(rich(jlon,jlat,jklev+1)+rich(jlon,jlat,jklev))
     if (zriav.gt..25) exit
     jkpbl = jklev
     enddo
     if (jkpbl.eq.0.and.rich(jlon,jlat,1).lt..25) jkpbl = 1
     if (jkpbl.gt.0) pbl(jlon,jlat) = zeta(jlon,jlat,jkpbl) + 10.
     pbl(jlon,jlat) = min (pbl(jlon,jlat), 2500.)

     gust(jlon,jlat) = sqrt(u_std_lev(jlon,jlat,1)**2+v_std_lev(jlon,jlat,1)**2)
     do jklev = 1, nlev/4
     jkup = jklev
     if (pbl(jlon,jlat).lt.zeta(jlon,jlat,jklev)) exit
     enddo
     if (jkup.gt.1) then
     do jklev = 1, jkup-1
     zwink = sqrt(u(jlon,jlat,jklev)**2+v(jlon,jlat,jklev)**2)
     gust(jlon,jlat) = max (gust(jlon,jlat), zwink)
     enddo
     endif

     enddo
     enddo

     return
     end
!##################################################################################################################
      subroutine convection_kf

! Kain-Fritsch (2004) convection (deep and shallow) scheme, revised

! July 2020: Moloch version

! Aug. 2018: max. timec decreased to avoid spurious expicit convection
! June 2012:
! changes mainly affecting shallow convection:
! dpmin lowered by 20% (also affects deep conv.); other corrections by M. Fantini (see "Mau")
! change in chmax (lowered from 1000 to 200 m to start shallow convection) - see "Andr"
! change in TIMEC (increased) for shallow convection - see "Andr"
! (the last two changes make shallow convection more effective on tropical oceans, with a shallow pbl,
! but at the same time limit the warming at about 850 hPa due to shallow conv. at mid-latitudes in summer,
! that seems excessive - perhaps shallow conv. algorithm should be revised, introducing iteration?

!  Oct. 2008 - parameterization of convection using conservation of liquid water static energy

      use mod_moloch, only : nlon, nlat, nlev, hxt, ps, u, v, t, tvirt, q, p, zeta, w, dz, &
                             r=>rd, eps, g, cp=>cpd, cv=>cvd, dx0=>dx, dy, t00=>tzer, cpv, &
                             snocon, raicon, dtdt, dqdt, dqcwdt, dqcidt, dqpwdt, dqpi1dt,  &
                             fmz, ntsrad, dtstep, nlonm1, nlatm1, nprocsy, myid, ntop
      implicit integer (i-n), real(4) (a-h, o-z)

! constants for saturation vapour pressure (Tetens expression)

      parameter (aliq=611., bliq=17.4,  cliq=bliq*273.15, dliq=33.65)

      parameter (ttfrz=268.15, tbfrz=248.15, rate=0.01)

!  constants used in the computation of gaussian integrals

      parameter (za1=0.4361836, za2=-0.1201676, za3=0.9372980, sigma=1./6., fe=1./0.202765151)
      parameter (zt1=0.500498, zc1=zt1*(za1+za2*zt1+za3*zt1**2), e45=.011108997)   ! e45=exp(-4.5)

      real, dimension(nlev) :: p0, dp, z0, t0, tv0, q0, u0, v0, w0, thta0, rh0, ems, emsd,                    &
                               umf, uer, udr, der, ddr, uer2, udr2, der2, ddr2, cldhgt, dilfrc,               &
                               qliq, qice, qlqout, qicout, pptliq, pptice, detlq, detic,                      &
                               tu, tvu, qu, td, qdu, qdd, thtau, thtad, tvqu,                                 &
                               fxm, thfxin, thfxout, qfxin, qfxout, qlfxin, qlfxout, qifxin, qifxout, qrfxin, &
                               qrfxout, qsfxin, qsfxout, tg, tvg, qg, qlg, qig, qrg, qsg, thtag

      jstart = 2
      jend   = nlatm1

!*******************************************************************
!                  environmental properties                        *
!*******************************************************************

      gdry   = -g/cp
      ep     = 1./eps - 1.
      dts    = dtstep*ntsrad
      do 999 jlat = jstart, jend
      dx     = dx0*hxt(jlat)
      dxsq   = dx*dy
      do 999 jlon = 2, nlonm1

      ishall = 0
      dpmin  = 4.5e3                  ! thickness of layer whose properties characterize the parcel
      p300   = ps(jlon,jlat)-30000.   ! pressure at 300 mb above surface
      kznew  = ntop

!  input:  temperature (t0, kelvin) ; specific humidity (q0, kg/kg) ;
!          horizontal wind speed (u0 and v0, m/s) ; pressure (p0, pascal) ;
!          height (z0, m);  vertical motion (w0, m/s)

      ml = 0
      l5 = 1
      do 15 k = 1, kznew
      p0(k)  = p(jlon,jlat,k)
      dp(k)  = p0(k)*g*dz/(r*tvirt(jlon,jlat,k)*fmz(jlon,jlat,k))  ! dp is the pressure interval at levels (hydrostatic)
      t0(k)  = t(jlon,jlat,k)
      q0(k)  = max(q(jlon,jlat,k), 1.e-10)
      u0(k)  = .5*(u(jlon,jlat,k)+u(jlon-1,jlat,k))
      v0(k)  = .5*(v(jlon,jlat,k)+v(jlon,jlat+1,k))
      zes    = aliq*exp((bliq*t0(k)-cliq)/(t0(k)-dliq))
      zqes   = 0.622*zes/(p0(k)-zes)
      q0(k)  = min (zqes, q0(k))         ! if q0 is above saturation value, reduce it to saturation level
      rh0(k) = q0(k)/zqes
      tv0(k) = t0(k)*(1.+ep*q0(k))
      w0(k)  = .5*(w(jlon,jlat,k)+w(jlon,jlat,k+1))
      z0(k)  = zeta(jlon,jlat,k)
      thta0(k) = t0(k)*(1.e5/p0(k))**(0.2854*(1.-0.28*q0(k)))  ! theta environment
      ems (k)  = dp(k)*dxsq/g                   ! ems is mass in the box: rho*volume
      emsd(k)  = 1./ems(k)

      if (p0(k).ge.500e2) l5   = k
      if (p0(k).ge.p300 ) llfc = k     ! llfc is the last level below p300
      if (t0(k).gt.t00  ) ml   = k     ! ml is the highest lev. with t above zero - melting level
  15  continue

      cldhgt = 0.    ! cloud thickness is initalized to zero at all levels

!*******************************************************************
!                 mixture                                          *
!*******************************************************************

      kmix = 1
  25  lc = kmix

      if (lc.gt.llfc) then
      chmax = 0.
      do nk = 1, llfc
        if (cldhgt(nk).gt.chmax) then
        chmax = cldhgt(nk)
        nchm = nk
        endif
      enddo

      if (chmax.gt.300.) then  ! Andr mar. 2017
!      goto 999     ! if uncommented, shallow convection is disallowed
      ishall = 1

      lc = nchm    ! pick the tallest cloud as candidate for shallow convection
      else
      goto 999     ! if shallow convection is not possible, go to the next grid point
      endif
      endif

!  assume that in order to support a deep updraft you need a layer of unstable air 50 to 100 mb deep.
!  to approximate this, isolate a group of adjacent individual model layers, with the base at level lc
!  such that the combined depth of these layers is at least dpmin

      dpthmx = 0.
      do nk = lc, kznew
      dpthmx = dpthmx + dp(nk)
      if (dpthmx.gt.dpmin) goto 64
      enddo
      goto 999
 64   kpbl = nk

!  go ahead and determine what level to start with for the next mixture in case the current mixture,
!  with base at level lc, is shallow or not buoyant.
!  instead of checking mixtures using every single layer, move up in increments of at least 15 mb

      do nk = lc+1, kznew
      if (p0(lc)-p0(nk).ge.15.e2) then
      kmix = nk
      go to 66
      endif
      enddo
      goto 999
 66   continue

!  find the thermodynamic characteristics of the layer by mass-weighting the characteristics
!  of the individual model layers

      tmix = 0.
      qmix = 0.
      zmix = 0.
      pmix = 0.
      do nk = lc, kpbl
      tmix = tmix + dp(nk)*t0(nk)
      qmix = qmix + dp(nk)*q0(nk)
      zmix = zmix + dp(nk)*z0(nk)
      pmix = pmix + dp(nk)*p0(nk)
      enddo
      tmix = tmix/dpthmx
      qmix = qmix/dpthmx
      zmix = zmix/dpthmx
      pmix = pmix/dpthmx
      emix = qmix*pmix/(0.622+qmix)    ! partial press. of mixture

!*******************************************************************
!                 lifted condensation level                        *
!*******************************************************************

!  find the temperature of the mixture at its lcl

      ztlog = log (emix/aliq)
      tdpt  = (cliq-dliq*ztlog)/(bliq-ztlog)   ! dew point
      tlcl  = tdpt -(.212+1.571e-3*(tdpt-t00)-4.36e-4*(tmix-t00))*(tmix-tdpt)
      tlcl  = min (tlcl, tmix)           ! case of over-saturated mixture
      tvlcl = tlcl*(1.+0.608*qmix)
      zlcl  = zmix + (tlcl-tmix)/gdry    ! note that the lapse rate is exactly dry adiabatic

      do nk = lc, kznew
      if (zlcl.le.z0(nk)) goto 35
      enddo
!      goto 999                  !  <<<<<======== sugg.: goto 25 if ishall=0 - Mau

      if(ishall.eq.0) then    ! Mau sugg. just above
      goto 25
      else
      goto 999
      endif

 35   klcl = nk           ! klcl is the lev. just above lcl

!  estimate environmental temperature and mixing ratio at the lcl

      dlp  = (zlcl-z0(klcl-1))/(z0(klcl)-z0(klcl-1))    ! calculate dlp using z instead of log(p)
      tenv = t0(klcl-1) + (t0(klcl)-t0(klcl-1))*dlp
      qenv = q0(klcl-1) + (q0(klcl)-q0(klcl-1))*dlp
      tven = tenv*(1.+0.608*qenv)

!*******************************************************************
!                 trigger                                          *
!*******************************************************************

!  check to see if cloud is buoyant using Fritsch-Chappell trigger function described in kain and fritsch (1992).
!  w0 is an approximate value for the running-mean grid-scale vertical velocity, which gives smoother fields
!  of convective initiation than the instantaneous value. Formula relating temperature perturbation to
!  vertical velocity has been used with the most success at grid lengths near 25 km.  For different grid-lengths,
!  adjust vertical velocity assuming linear dependence of w on grid length.

      if (zlcl.lt.2.e3) then
      wkl0 = 0.02*zlcl/2.e3     ! wkl0 is a threshold defined linearly up to 2 cm/sec
      else
      wkl0 = 0.02
      endif

      wkl = (w0(klcl-1) + (w0(klcl)-w0(klcl-1))*dlp)*dx/25.e3 - wkl0
      if (wkl.lt.0.) then
      dtlcl = 0.
      rad   = 1000.
      else
      dtlcl = 3.13*wkl**0.33
      rad   = 1000.+1000*wkl/0.1
      rad   = min (rad, 2000.)
      endif

!  if no convection, goes up and uses the new kmix defined 15 mb above the previous one

      if (tlcl+dtlcl.le.tenv.and.ishall.eq.0) goto 25               ! <<<<=======  Mau

!*******************************************************************
!                 compute updraft properties                       *
!*******************************************************************

!  convective triggering criteria has been satisfied. Compute liquid water static energy (hw)

      hw = (cp+cpv*qmix)*tmix + (1.+qmix)*g*zmix

!  calculation of initial vertical velocity of the parcel. JSK 11/26/97

      if (dtlcl.gt.0.) then
      gdt  = 2.*g*dtlcl*500./tven
      wlcl = 1.+0.5*sqrt(gdt)
      wlcl = min (wlcl,3.)
      else
      wlcl = 1.
      endif
      wtw  = wlcl**2
      au0  = 0.01*dxsq
      plcl = p0(klcl-1) + (p0(klcl)-p0(klcl-1))*dlp

!  ttemp is used during calculation of the linear glaciation process.
!  It is initially set to the temperature at which freezing is specified to begin.
!  Within the glaciation interval, it is set equal to the updraft temp. at the previous model level

      ttemp = ttfrz          ! ttfrz defined -5 deg C

!  estimate initial updraft mass flux umf(klcl-1)

      vmflcl       = wlcl *au0 * plcl/(r*tvlcl)
      umf(klcl-1)  = vmflcl  ! defined at the lower interface of layer klcl

      uer(klcl-1)  = 0.
      tvu(klcl-1)  = tvlcl
      qu(klcl-1)   = qmix
      qliq(klcl-1) = 0.
      qice(klcl-1) = 0.
      let = klcl             ! level of neutral buoyancy (equilibrium temp.)
      rei = 0.               ! rate of environmental inflow
      ee1 = 1.               ! env. entrainment fraction
      ud1 = 0.               ! updraft detrained fraction
      abe = 0.               ! available buoyant energy
      trppt = 0.             ! total rate of precipitation production

!  enter the loop for updraft calculations. calculate updraft temp, mixing ratio, vertical mass flux,
!  lateral detrainment of mass and moisture, precipitation rates at each model level

      do 60 nk = klcl, kznew

      qu(nk)     = qu(nk-1)
      qliq(nk)   = qliq(nk-1)
      qice(nk)   = qice(nk-1)

!---------------------------
!  parcel rise
!---------------------------

!  p0,z0: new pressure and height (input)
!  hw: updraft liq. water static energy (input)
!  tu: updraft temperature (output)
!  qu: updraft humidity (input-output)
!  qliq, qice: updraft old condensate (input-output)
!  qlnew: fresh condensate (output), all liquid

      call hliq (hw, z0(nk), p0(nk), tu(nk), qu(nk), qliq(nk), qice(nk), qlnew, g, cp, cpv, 1.)

!---------------------------
!  glaciation
!---------------------------

!  check to see if updraft temp is above the temperature at which glaciation is assumed to initiate;
!  if it is, calculate the fraction of remaining liquid water to freeze.
!  ttfrz is the temp at which freezing begins, ttfrz=-5 deg c
!  tbfrz is the temp below which all liquid water is frozen at each level, tbfrz=-25 deg c

      qinew = 0.
      if (tu(nk).le.ttfrz) then

!  determine the effects of liquid water freezing when temperature is below ttfrz;
!  comput. of frc1: fraction of liquid that freezes in the current layer

      if (tu(nk).gt.tbfrz) then
!      frc1 = (ttemp-tu(nk))/(ttemp-tbfrz)
      frc1 = min(1.,max((ttemp-tu(nk))/(ttemp-tbfrz), 0.)) ! Mau, Dec. 2011
      else
      frc1 = 1.          ! for colder temperatures, freeze all liquid water.
      endif
      ttemp = tu(nk)     ! this ttemp will be used to compute frc1 at the level above

      qinew = qlnew*frc1
      qlnew = qlnew-qinew
      qice(nk) = qice(nk)+qliq(nk)*frc1+qinew
      qliq(nk) = qliq(nk)-qliq(nk)*frc1+qlnew

!  freezing warms the air and it becomes unsaturated. assume that water and ice evaporate to maintain saturation

      call hliq (hw, z0(nk), p0(nk), tu(nk), qu(nk), qliq(nk), qice(nk), zzz, g, cp, cpv, 1.)
      qliq(nk) = qliq(nk)-qlnew      ! fresh and old condensate are separated again
      qice(nk) = qice(nk)-qinew
      endif

!------------------------------------------------
!  vertical velocity and precipitation production
!------------------------------------------------

!  calculate updraft vertical velocity and precipitation fallout
!  Bechtold et al 2001 QJ eq(12). Factor 1/1.5 takes into account non-hydrostatic effect

      tvu(nk) = tu(nk)*(1.+0.608*qu(nk))    ! virtual temp. of the updraft, here without water loading
      if (nk.eq.klcl) then
      dzz = z0(nk)-zlcl
      boterm  = 2.*((tvlcl+tvu(nk))/(tven+tv0(nk))-1.)*dzz*g/1.5
      else
      dzz = z0(nk)-z0(nk-1)
      boterm  = 2.*((tvu(nk-1)+tvu(nk))/(tv0(nk-1)+tv0(nk))-1.)*dzz*g/1.5  !  (t_parcel-t_env)/t_env
      endif

!  estimate the vertical velocity so that an average vertical velocity can
!  be calculated to estimate the time required for ascent between model levels

      enterm = 2.*rei*wtw/umf(nk-1)
      zqest  = .5*(qliq(nk)+qice(nk)+qlnew+qinew)
      wtw1   = max ( wtw+boterm-enterm-2.*g*dzz*zqest/1.5, 0. )
      zwavg  = .5*(sqrt(wtw)+sqrt(wtw1))

!  this precipitation fallout scheme is based on the scheme used by Ogura and Cho (1973).
!  Liquid water fallout from a parcel is calculated using dq=-rate*q*dt, but to simulate a quasi-continuous
!  process, and to eliminate a dependency on vertical resolution, this is expressed as q=q*exp(-rate*dz)
!  only 60% of the fresh condensate is allowed to participate in the conversion process

      zqconv = qliq(nk) + qice(nk) + .6*(qlnew+qinew)
!      zdq    = zqconv*(1.-exp(-rate*dzz/zwavg))  ! total loss due to precipitation

      rate1 = 0.008 + 0.0018*(273.15-tu(nk))                       ! Andr mar. 2017
      rate1 = min(rate1, 0.03)                                     ! Andr mar. 2017
      rate1 = max(rate1, 0.008)                                    ! Andr mar. 2017
      zdq  = zqconv*(1.-exp(-rate1*dzz/zwavg))  ! total loss due to precipitation

!  estimate the mean load of condensate on the updraft in the layer, calculate vertical velocity

      pptdrg = qliq(nk)+qice(nk) +.5*(qlnew+qinew-zdq)   ! average water loading
      wtw = wtw + boterm - enterm -2.*g*dzz*pptdrg/1.5

!  ratio3 is the fraction of liquid water in fresh condensate
!  ratio4 is the fraction of liquid water in the total amount of condensate involved in the precipitation process

      ratio3 = qlnew/(qlnew+qinew+1.e-10)
      ratio4 = (qliq(nk)+0.6*qlnew)/(zqconv+1.e-10)

!  determine the new liquid water and ice concentrations and precipitation loss

      qlqout(nk) = ratio4*zdq
      qicout(nk) = (1.-ratio4)*zdq
      qliq(nk)   = ratio4*(zqconv-zdq) + ratio3*0.4*(qlnew+qinew)
      qice(nk)   = (1.-ratio4)*(zqconv-zdq) + (1.-ratio3)*0.4*(qlnew+qinew)

!  if vert vel is negative, exit the updraft loop and, if cloud is tall enough, finalize updraft calculations

      if (wtw.lt.1.e-3) goto 65   ! exit updraft

!  update the abe

      tvqu(nk) = tu(nk)*(1.+0.608*qu(nk)-qliq(nk)-qice(nk)) ! virtual temp. with water loading
      if (nk.eq.klcl) then
      dilbe = ((tvlcl+tvqu(nk))/(tven+tv0(nk))-1.)*dzz*g
      else
      dilbe = ((tvqu(nk-1)+tvqu(nk))/(tv0(nk-1)+tv0(nk))-1.)*dzz*g
      endif
      abe = abe + max(dilbe,0.)

!-------------------------
!  entrainment/detrainment
!-------------------------

      rei = vmflcl*dp(nk)*0.03/rad         ! rei is the rate of environmental inflow

!  updraft liquid water static energy

      rll = 2.5008e6-2369.*(tu(nk)-273.16)
      rls = 2.8345e6-260. *(tu(nk)-273.16)
      zqtot = qu(nk) + qliq(nk) + qice(nk)
      hwu = (cp+cpv*zqtot)*tu(nk) -rll*qliq(nk) -rls*qice(nk) +(1.+zqtot)*g*z0(nk)

!  environmental liquid water static energy to be mixed with parcel's at constant pressure

      hw0 = (cp+q0(nk)*cpv)*t0(nk) +(1.+q0(nk))*g*z0(nk)

!  if cloud parcels are virtually colder than the environment, no entrainment is allowed at this level

      if (tvqu(nk).le.tv0(nk)) then
      ee1 = 0.5
      ud1 = 1.0
      uer(nk) = 0.5*rei
      udr(nk) = 1.5*rei
      else
      let = nk               ! level of equilibrium temperature (neutral buoyancy)

!  determine the fractional entrain. and detrain. rates  ee2, ud2  at the current level (.5<=ee2<=1, 0<=ud2<=1.5)

      qtmp   = .95*q0(nk) + .05*qu(nk)      ! mix 95% envir. air and 5% updraft air
      tmpliq = .05*qliq(nk)
      tmpice = .05*qice(nk)
      hwmix  = .95*hw0 + .05*hwu
      call hliq (hwmix, z0(nk), p0(nk), ttmp, qtmp, tmpliq, tmpice, qlnew, g, cp, cpv, 1.)
      tu95 = ttmp*(1.+0.608*qtmp-tmpliq-tmpice)
        if (tu95.gt.tv0(nk)) then             ! if still buoyant, entrainm.=max, detrainm.=min
        ee2 = 1.
        ud2 = 0.
        goto 50
        endif

      qtmp   = .1*q0(nk) + .9*qu(nk)        ! mix 10% envir. air and 90% updraft air
      tmpliq = .9*qliq(nk)
      tmpice = .9*qice(nk)
      hwmix  = .1*hw0 + .9*hwu
      call hliq (hwmix, z0(nk), p0(nk), ttmp, qtmp, tmpliq, tmpice, qlnew, g, cp, cpv, 1.)
      tu10 = ttmp*(1.+0.608*qtmp-tmpliq-tmpice)
      tvdiff = tu10-tvqu(nk)
        if (abs(tvdiff).lt.1.e-3) then
        ee2 = 1.
        ud2 = 0.
        goto 50
        endif

!  determine the critical mixed fraction of updraft and environmental air (eqfrc)
!  newton step: tvqu + (d tv)/(d f)*eqfrc = tv0, eqfrc fraction of env. air in neutrally buoyant mixture

      eqfrc = (tv0(nk)-tvqu(nk))*.1/tvdiff
      eqfrc = max (0. ,eqfrc)
      eqfrc = min (1. ,eqfrc)

      if (eqfrc.eq.1.) then
      ee2 = 1.
      ud2 = 0.
      elseif (eqfrc.eq.0.) then
      ee2 = 0.5
      ud2 = 1.5
      else

!  integral over gaussian distribution - the numerical approximation to the integrals are taken from
!  "Handbook of mathematical functions with formulas, graphs and math. tables" ed. by Abramowitz and Stegun,
!  Nat. Bureau of Standards, Applied mathematics series.  June, 1964., May, 1968.

      zy  = 6.*eqfrc-3.
      ey  = exp(-.5*zy**2)
        if (zy.ge.0.) then
        zt2 = 1./(1.+0.33267*zy)
        zc2 = zt2*(za1+za2*zt2+za3*zt2**2)
        ee2 = sigma*(.5*(2.506628-e45*zc1-ey*zc2) + sigma*(e45-ey)) - e45*    .5*eqfrc**2
        ud2 = sigma*(.5*(        -e45*zc1+ey*zc2) + sigma*(e45-ey)) - e45*(.5+.5*eqfrc**2-eqfrc)
        else
        zt2 = 1./(1.-0.33267*zy)
        zc2 = zt2*(za1+za2*zt2+za3*zt2**2)
        ee2 = sigma*(.5*(        -e45*zc1+ey*zc2) + sigma*(e45-ey)) - e45*    .5*eqfrc**2
        ud2 = sigma*(.5*(2.506628-e45*zc1-ey*zc2) + sigma*(e45-ey)) - e45*(.5+.5*eqfrc**2-eqfrc)
        endif
      ee2 = max (ee2*fe, .5)
      ud2 = 1.5*ud2*fe
      endif

!  net entrainment and detrainment rates are given by the average fractional values in the layer

 50   uer(nk) = .5*rei*(ee1+ee2)
      udr(nk) = .5*rei*(ud1+ud2)
      ee1 = ee2
      ud1 = ud2
      endif               !  end detrainment/entrainment section

!  if the calculated updraft detrainment rate is greater than the total
!  updraft mass flux, all cloud mass detrains, exit updraft calculations

      if (umf(nk-1)-udr(nk).lt.10.) then
      abe = abe - max(dilbe,0.)
      let = nk-1
      goto 65         ! exit updraft
      endif

      umf(nk)    = umf(nk-1)-udr(nk)+uer(nk)     ! umf(nk) is defined at the interface above
      dilfrc(nk) = umf(nk)/(umf(nk)-uer(nk))     ! dilfrc >= 1

!  detlq and detic are the rates of detrainment of liquid water and ice in the detraining updraft mass

      detlq(nk) = qliq(nk)*udr(nk)       ! total detrained water
      detic(nk) = qice(nk)*udr(nk)       ! total detrained ice
      qdu(nk)   = qu(nk)                 ! humidity which is detrained
      qu(nk)    = (qu(nk) + q0(nk)*(dilfrc(nk)-1.))/dilfrc(nk)
      qliq(nk)  = qliq(nk)/dilfrc(nk)
      qice(nk)  = qice(nk)/dilfrc(nk)
      hw = (hwu + hw0*(dilfrc(nk)-1.))/dilfrc(nk)

!  pptliq is the rate of generation (fallout) of liquid precip at a given model lvl, pptice the same for ice,
!  trppt is the total rate of production of precipitation

      pptliq(nk) = qlqout(nk)*umf(nk-1)
      pptice(nk) = qicout(nk)*umf(nk-1)
      trppt = trppt + pptliq(nk) + pptice(nk)
      if (nk.le.kpbl) uer(nk) = uer(nk) + vmflcl*dp(nk)/dpthmx
 60   continue                            ! fine loop updraft
 65   ltop = nk-1

      if (ishall.eq.0) then

!  do not allow any cloud from this layer if:
!  1.) cloud top is at model level just above lcl, or
!  2.) cloud top is within updraft source layer, or
!  3.) cloud-top detrainment layer begins within updraft source layer.

      if (ltop.le.klcl .or. ltop.le.kpbl .or. let+1.le.kpbl) goto 25  ! check next level

!  if cloud top height is less than the specified minimum for deep convection, save value
!  to consider this level as source for shallow convection, go back up to check next level
!  try specifying minimum cloud depth as a function of tlcl

      cldhgt(lc) = z0(ltop)-zlcl       ! define cloud depth
      if (tlcl.gt.293.) then
!      chmin = 4.e3
      chmin = 2.8e3 ! Andr oct. 2017
      elseif (tlcl.le.293. .and. tlcl.ge.273.) then
!      chmin = 2.e3 + 100.*(tlcl-273.)
      chmin = 2.e3 + 40.*(tlcl-273.) ! Andr oct. 2017
      elseif (tlcl.lt.273.) then
      chmin = 2.e3
      endif

      if (cldhgt(lc).lt.chmin .or. abe.lt.1.) goto 25  ! check next level
      else     ! ishall=1
      let = max (kpbl, klcl)
      endif    ! condition on ishall

!  if the let and ltop are the same, detrain all of the updraft mass flux at this level

      if (let.eq.ltop) then
      udr(ltop)   = udr(ltop)+umf(ltop)-uer(ltop)
      uer(ltop)   = 0.
      umf(ltop)   = 0.
      detlq(ltop) = qliq(ltop)*udr(ltop)*dilfrc(ltop)
      detic(ltop) = qice(ltop)*udr(ltop)*dilfrc(ltop)
      else      !  begin total detrainment at the level above the let

!  adjust mass flux profiles, detrainment rates, and precipitation fall
!  rates to reflect the linear decrease in mass flux between the let and ltop

      dptt=0.
      do  nk = let+1, ltop
      dptt = dptt + dp(nk)
      enddo
      dumfdp = umf(let)/dptt

!  entrainment is allowed at every level except for ltop, so disallow
!  entrainment at ltop and adjust entrainment rates between let and ltop;
!  maintains the ratio between mass flux at the interface above and entrainment computed before

      do nk = let+1, ltop
        if (nk.eq.ltop) then
        udr(nk)   = umf(nk-1)
        uer(nk)   = 0.
        umf(nk)   = 0.
        detlq(nk) = udr(nk)*qliq(nk)*dilfrc(nk)
        detic(nk) = udr(nk)*qice(nk)*dilfrc(nk)
        else
        umf(nk)   = umf(nk-1)-dp(nk)*dumfdp
        uer(nk)   = umf(nk)*(1.-1./dilfrc(nk))
        udr(nk)   = umf(nk-1)-umf(nk)+uer(nk)
        detlq(nk) = udr(nk)*qliq(nk)*dilfrc(nk)
        detic(nk) = udr(nk)*qice(nk)*dilfrc(nk)
        endif

        if (nk.ge.let+2) then
        trppt = trppt - pptliq(nk) - pptice(nk)
        pptliq(nk) = umf(nk-1)*qlqout(nk)
        pptice(nk) = umf(nk-1)*qicout(nk)
        trppt = trppt + pptliq(nk) + pptice(nk)
        endif
      enddo
      endif

!  extend the updraft mass flux profile down to the source layer for the updraft air

      kkk = min(kpbl,klcl-1)
      uer(1:lc-1) = 0.
      umf(1:lc-1) = 0.
      uer(lc:kkk) = vmflcl*dp(lc:kkk)/dpthmx
      umf(lc) = uer(lc)
      do nk = lc+1, kkk
      umf(nk) = umf(nk-1) + uer(nk)
      enddo
      uer(kpbl+1:klcl-1) = 0.      ! no more entrainment above the top of mixture and below klcl
      umf(kpbl+1:klcl-1) = vmflcl  ! mas flux is constant

      tu(1:lc-1)    = 0.
      qu(1:lc-1)    = 0.
      tu(lc:klcl-1) = tmix + (z0(lc:klcl-1)-zmix)*gdry   ! dry adiabatic lapse rate above the source layer
      qu(lc:klcl-1) = qmix
      udr   (1:klcl-1) = 0.
      qdu   (1:klcl-1) = 0.
      qliq  (1:klcl-1) = 0.
      qice  (1:klcl-1) = 0.
      pptliq(1:klcl-1) = 0.
      pptice(1:klcl-1) = 0.
      detlq (1:klcl-1) = 0.
      detic (1:klcl-1) = 0.

!  updraft potential temperature

      do nk = 1, ltop
      thtau(nk) = tu(nk)*(1.e5/p0(nk))**(0.2854*(1.-0.28*qu(nk)))
      enddo

!*****************************************************************
!                  compute downdraft properties                  *
!*****************************************************************

      tder  = 0.   ! total downdraft evaporation rate
      der   = 0.   ! downdraft entrainment rate
      ddr   = 0.   ! downdraft detrainment rate
      qdd   = 0.   ! downdraft detrained humidityi
      thtad = 0.   ! downdraft detrained pot. temperature

      if (ishall.eq.1)  goto 141      ! 141 exit for no downdraft

!  if lfs is not at least 50 mb above cloud base (implying that the
!  level of equilibrium temp, let, is just above cloud base) do not allow a downdraft

      kstart = kpbl+1
      lfs    = let-1
      do k = kstart+1, let-2
        if (p0(kstart)-p0(k).gt.200.e2) then
        lfs = k         ! lfs: 200 mb above p0(kstart), or let-1
        exit
        endif
      enddo
      if (p0(kstart)-p0(lfs).lt.50.e2) goto 141   ! exit downdraft

      rhbar = rh0(lfs)*dp(lfs)
      dptt = dp(lfs)
      do k = lfs-1, kstart, -1
      dptt = dptt +dp(k)
      rhbar = rhbar + rh0(k)*dp(k)
      enddo
      rhbar = min (rhbar/dptt,.99)  ! Marino
!      dmf = min(3.*(1.-rhbar),1.5) * umf(klcl)  !  downdraft mass flux (Kain, J.App.Met., 2004, eq.11)
      dmf = 2.5*(1.-0.5*rhbar-0.5*rhbar*rhbar) * umf(klcl) ! Andr
      hwd = 0.
      qd  = 0.
      do k = lfs, kstart, -1
      der(k) = dmf*dp(k)/dptt
      hwd = hwd + ((cp+q0(k)*cpv)*t0(k)+(1.+q0(k))*g*z0(k))*der(k)
      qd  = qd  + q0(k)*der(k)
      enddo
      hwd = hwd/dmf
      qd  = qd/dmf

!  calculate total liquid and iced precipitation and melting effect

      qld = 0.
      qid = 0.
      do k = klcl, ltop
      qld = qld + pptliq(k)
      qid = qid + pptice(k)
      enddo
      qld = qld/dmf
      qid = qid/dmf
      if (lc.lt.ml) then  ! only place where ml is used (level where temp=0 deg C)
      qld = qld+qid
      qid = 0.
      endif

!  define (approx.) liquid water static energy at top of the downdraft detrainment layer

      rll = 2.5008e6-2369.*(t0(kstart)-273.16)
      rls = 2.8345e6-260. *(t0(kstart)-273.16)
      hwd = hwd -rll*qld -rls*qid +cpv*t0(kstart)*(qld+qid)

!  downdraft detrainment

      dpdd = 0.
      do nd = kstart-1, 1, -1
      dpdd = dpdd + dp(nd)
      qdd(nd) = qd
      qld0 = qld
      qid0 = qid

!  descent, hence qdd > qd. The missing part should come from evap. of precip.
!  td is temp. reached due to evap. to saturation,
!  but here evap. is required not up to saturation but only to pre-determined rh

      rh = 1.-0.16/1000.*(z0(kstart)-z0(nd))   ! specify relative humidity decrease of 16%/km in downdraft
      call hliq (hwd, z0(nd), p0(nd), td(nd), qdd(nd), qld0, qid0, qlnew, g, cp, cpv, rh)

!  if the downdraft becomes virtually warmer than the environment, it does not sink any more;
!  defines ldb - downdraft base, see also Bechtold

      ztvd = td(nd)*(1.+0.608*qdd(nd))
      if (ztvd.gt.tv0(nd).or.nd.eq.1) then
      ldb = nd
      exit
      endif
      enddo             ! end downdraft detrainment

!  calculate the total downdraft evaporation rate (tder)

      do k = kstart-1, ldb, -1
      ddr(k)   = dmf*dp(k)/dpdd         ! detrainment linearly distributed (in pressure)
      tder     = tder + (qdd(k)-qd)*ddr(k)
      thtad(k) = td(k)*(1.e5/p0(k))**(0.2854*(1.-0.28*qdd(k)))
      enddo

      ddinc = min(1., trppt/(tder+1.))  ! reduce downd. mass flux if evaporation is larger than total precip.
      tder  = tder*ddinc
      trppt = trppt-tder                ! residual precipitation (excluding evaporation)
      do k = ldb, lfs
      der(k) = der(k)*ddinc
      ddr(k) = ddr(k)*ddinc
      enddo

 141  continue  ! gets here when downdraft is excluded (shallow cloud, p0(kstart)-p0(lfs)<50mb)

!  compute convective time scale (timec). The mean wind at the lcl and midtroposphere is used.

      wspd_klcl = sqrt (u0(klcl)**2 + v0(klcl)**2)
      wspd_l5   = sqrt (u0(l5  )**2 + v0(l5  )**2)
      vconv    = .5*(wspd_klcl+wspd_l5) + 0.001
      timec    = dx/vconv
      timec    = max (900. ,timec)
      timec    = min (1500.,timec)
      timec    = max (dts  ,timec)

      if (ishall.eq.1) timec = 5400. ! care:  redef. near the end - Andr oct. 2017

!  set limits on the updraft and downdraft mass fluxes so that the inflow
!  into convective drafts from a given layer is no more than is available in that layer initially

      aincmax = 1000.
      do nk = lc, ltop
      if (uer(nk)+der(nk).gt.1.e-3) then
      aincmax = min (aincmax, 0.99*ems(nk)/((uer(nk)+der(nk))*timec))
      endif
      enddo

      if (ishall.eq.1) then
      tkemax = 5.
      ainc = 0.05*tkemax*dpthmx*dxsq/(vmflcl*g*timec)
      ainc = min (ainc, aincmax)
      else
      ainc = min (1., aincmax)
      endif

      udr2(1:ltop) = udr(1:ltop)*ainc
      uer2(1:ltop) = uer(1:ltop)*ainc
      der2(1:ltop) = der(1:ltop)*ainc
      ddr2(1:ltop) = ddr(1:ltop)*ainc
      ncount = 0
      abegold = abe
      aincold = 0.
 175  ncount = ncount+1

!*****************************************************************
!           compute properties for compensational subsidence     *
!*****************************************************************

! determine mass flux fxm at bottom of each layer to satisfy mass continuity

      dtt      = timec
      thtag(1) = thta0(1)
      qg(1)    = q0(1)
      fxm(1)   = 0.
      do nk = 2, ltop
      fxm(nk) = fxm(nk-1) - uer2(nk-1)-der2(nk-1)+udr2(nk-1)+ddr2(nk-1)
      dtt = min (dtt, .75/(abs(fxm(nk))*emsd(nk)+1.e-12))
      thtag(nk) = thta0(nk)
      qg(nk)    = q0(nk)
      enddo

      nstep = nint (timec/dtt+1.)
      dtime = timec/float(nstep)

!  do an upstream/forward-in-time advection of theta and q

      do 495 ntc = 1, nstep

!  assign theta and q values at the top and bottom of each layer based on sign of fxm

      thfxin (1:ltop) = 0.
      thfxout(1:ltop) = 0.
      qfxin  (1:ltop) = 0.
      qfxout (1:ltop) = 0.
      do nk = 2, ltop  ! loop on the interfaces
      if (fxm(nk).ge.0.) then
      thfxin(nk)    = fxm(nk)*thtag(nk-1)
      qfxin(nk)     = fxm(nk)*qg(nk-1)
      thfxout(nk-1) = thfxout(nk-1) + thfxin(nk)
      qfxout(nk-1)  = qfxout(nk-1)  + qfxin(nk)
      else
      thfxout(nk)   = -fxm(nk)*thtag(nk)
      qfxout(nk)    = -fxm(nk)*qg(nk)
      thfxin(nk-1)  = thfxin(nk-1)+thfxout(nk)
      qfxin(nk-1)   = qfxin(nk-1)+qfxout(nk)
      endif
      enddo

!  update the theta and q values at each level

      do nk = 1, ltop
      thtag(nk) = thtag(nk) + ( thfxin(nk) + udr2(nk)*thtau(nk) + ddr2(nk)*thtad(nk) - thfxout(nk)- &
                  (uer2(nk)+der2(nk))*thta0(nk) )*dtime*emsd(nk)
      qg(nk)    = qg(nk) + ( qfxin(nk) + udr2(nk)*qdu(nk) + ddr2(nk)*qdd(nk) - qfxout(nk)- &
                  (uer2(nk)+der2(nk))*q0(nk) )*dtime*emsd(nk)
      enddo
  495 continue

!  check to see if mixing ratio dips below zero anywhere;  if so, borrow
!  moisture from adjacent layers to bring it back up above zero

      if (qg(1).lt.0.) qg(1) = 1.e-9
      do nk = 2, ltop
        if (qg(nk).lt.0.) then
        nk1 = nk+1
        if (nk.eq.ltop) nk1 = klcl   ! tele-transport
        tma = qg(nk1)*ems(nk1)
        tmb = qg(nk-1)*ems(nk-1)
        tmm = (qg(nk)-1.e-9)*ems(nk)
        bcoeff = -tmm/(tma**2/tmb+tmb)
        acoeff = bcoeff*tma/tmb
        tmb = tmb*(1.-bcoeff)
        tma = tma*(1.-acoeff)
        qg(nk)  = 1.e-9
        qg(nk1) = tma*emsd(nk1)
        qg(nk-1)= tmb*emsd(nk-1)
        endif
      enddo

      do nk = 1, ltop
      tg(nk) = thtag(nk)*(1.e5/p0(nk))**(-0.2854*(1.-0.28*qg(nk)))  ! convert theta to t
      enddo

      if (ishall.eq.1 .or. ncount.gt.7) goto 265     ! exit iteration

!----------------------------------------------------------
!  compute new cloud and change in available buoyant energy
!----------------------------------------------------------

!  find the thermodynamic characteristics of the layer by
!  mass-weighting the characteristics of the individual model layers

      tmix = 0.
      qmix = 0.
      do nk = lc, kpbl
      tmix = tmix + dp(nk)*tg(nk)
      qmix = qmix + dp(nk)*qg(nk)
      enddo
      tmix = tmix/dpthmx
      qmix = qmix/dpthmx

!  remove supersaturation, if necessary

      es  = aliq*exp((tmix*bliq-cliq)/(tmix-dliq))
      qss = 0.622*es/(pmix-es)
      if (qmix.gt.qss) then
      rll = 2.5008e6-2369.*(tmix-273.16)
      cpm = cp*(1.+0.887*qmix)
      dssdt = qss*(cliq-bliq*dliq)/(tmix-dliq)**2
      dq = (qmix-qss)/(1.+rll*dssdt/cpm)
      tmix = tmix + rll/cp*dq
      qmix = qmix-dq
      tlcl = tmix
      else
      qmix = max (qmix, 1.e-12)
      emix = qmix*pmix/(0.622+qmix)
      tlog = log(emix/aliq)
      tdpt = (cliq-dliq*tlog)/(bliq-tlog)
      tlcl = tdpt - (.212+1.571e-3*(tdpt-t00)-4.36e-4*(tmix-t00))*(tmix-tdpt)
      tlcl = min (tlcl,tmix)
      endif
      tvlcl = tlcl*(1.+0.608*qmix)
      zlcl  = zmix + (tlcl-tmix)/gdry

      do nk = lc, kznew    ! note: recomputation of klcl
      klcl = nk
      if (zlcl.le.z0(nk)) goto 735
      enddo
      print*, 'Bad adjust - lcl cannot be found after conv. adjustment'
      goto 999
 735  dlp = (zlcl-z0(klcl-1))/(z0(klcl)-z0(klcl-1))

!  estimate environmental temperature and mixing ratio at the lcl

      tenv = tg(klcl-1) + (tg(klcl)-tg(klcl-1))*dlp
      qenv = qg(klcl-1) + (qg(klcl)-qg(klcl-1))*dlp
      tven = tenv*(1.+0.608*qenv)
      hw   = (cp+cpv*qmix)*tmix + (1.+qmix)*g*zmix
      qgu  = qmix
      qlu  = 0.
      qiu  = 0.

!  compute adjusted abe (abeg) -  undiluted and 'heavy' updraft

      abeg = 0.
      do k = klcl, ltop
      call hliq (hw, z0(k), p0(k), tgu, qgu, qlu, qiu, qlnew, g, cp, cpv, 1.)
      qlu = qlu + qlnew
      tvqu(k) = tgu*(1.+0.608*qgu-qlu)
      tvg(k)  = tg(k)*(1.+0.608*qg(k))
        if (k.eq.klcl) then
        dilbe = ((tvlcl+tvqu(k))/(tven+tvg(k))-1.)*(z0(k)-zlcl)*g
        else
        dilbe = ((tvqu(k-1)+tvqu(k))/(tvg(k-1)+tvg(k))-1.)*(z0(k)-z0(k-1))*g
        endif
      abeg = abeg + max(dilbe,0.)
      enddo
      if (ncount.eq.1.and.abeg.gt.abe) abeg=abe*.9

!  assume at least 90% of cape (abe) is removed by convection during the period timec

      if (abeg.gt.abe) then
!      print*,'abeg.gt.abe',ncount,abeg,abe
      goto 999
      endif
      if (abeg.le..1*abe) goto 265  ! CAPE reduced to less than 10% of its initial value: exit iteration

      dabeda = (abeg-abegold)/(ainc-aincold)
      if (dabeda.ge.0.) goto 265 ! if it goes away form min. of abeg before reaching 10% of abe, exit

!  if more than 10% of the original cape remains, and cape is still decreasing,
!  increase the convective mass flux by the factor ainc computed by a Newton step

      abegold = abeg
      aincold = ainc
      ainc = ainc - (abeg-.02*abe)/dabeda

!  exit iteration with the old ainc value if mass is not enough to increase flux at lcl

      if (ainc.gt.aincmax) then
      ainc = aincold
      goto 265
      endif

 255  udr2(1:ltop) = udr(1:ltop)*ainc
      uer2(1:ltop) = uer(1:ltop)*ainc
      der2(1:ltop) = der(1:ltop)*ainc
      ddr2(1:ltop) = ddr(1:ltop)*ainc
      goto 175     ! go back up for another iteration

 265  continue

!*******************************************************************
!                       final section                              *
!*******************************************************************

!  compute hydrometeor tendencies as is done for t, q

      qlg = 0.
      qig = 0.
      qrg = 0.
      qsg = 0.

      do 290 ntc = 1, nstep

!  assign hydrometeors concentrations at the top and bottom of each layer based on the sign of fxm

      qlfxin (1:ltop) = 0.
      qlfxout(1:ltop) = 0.
      qifxin (1:ltop) = 0.
      qifxout(1:ltop) = 0.
      qrfxin (1:ltop) = 0.
      qrfxout(1:ltop) = 0.
      qsfxin (1:ltop) = 0.
      qsfxout(1:ltop) = 0.
      do nk = 2, ltop
        if (fxm(nk).ge.0.) then
        qlfxin(nk) = fxm(nk)*qlg(nk-1)
        qifxin(nk) = fxm(nk)*qig(nk-1)
        qrfxin(nk) = fxm(nk)*qrg(nk-1)
        qsfxin(nk) = fxm(nk)*qsg(nk-1)
        qlfxout(nk-1) = qlfxout(nk-1)+qlfxin(nk)
        qifxout(nk-1) = qifxout(nk-1)+qifxin(nk)
        qrfxout(nk-1) = qrfxout(nk-1)+qrfxin(nk)
        qsfxout(nk-1) = qsfxout(nk-1)+qsfxin(nk)
        else
        qlfxout(nk) = -fxm(nk)*qlg(nk)
        qifxout(nk) = -fxm(nk)*qig(nk)
        qrfxout(nk) = -fxm(nk)*qrg(nk)
        qsfxout(nk) = -fxm(nk)*qsg(nk)
        qlfxin(nk-1) = qlfxin(nk-1)+qlfxout(nk)
        qifxin(nk-1) = qifxin(nk-1)+qifxout(nk)
        qrfxin(nk-1) = qrfxin(nk-1)+qrfxout(nk)
        qsfxin(nk-1) = qsfxin(nk-1)+qsfxout(nk)
        endif
      enddo

!  update the hydrometeor concentration values at each level

      frc2 = trppt/(trppt+tder+1.e-15)  ! fraction of total condensate generated that goes into precipitation
      do nk = 1, ltop
      rainfb  = pptliq(nk)*ainc*frc2*float(ishall)
      snowfb  = pptice(nk)*ainc*frc2*float(ishall)
      qlg(nk) = qlg(nk) + (qlfxin(nk) +detlq(nk)*ainc -qlfxout(nk))*dtime*emsd(nk)
      qig(nk) = qig(nk) + (qifxin(nk) +detic(nk)*ainc -qifxout(nk))*dtime*emsd(nk)
      qrg(nk) = qrg(nk) + (qrfxin(nk) -qrfxout(nk) +rainfb)*dtime*emsd(nk)
      qsg(nk) = qsg(nk) + (qsfxin(nk) -qsfxout(nk) +snowfb)*dtime*emsd(nk)
      enddo

 290  continue

!  evaluate moisture budget

      qinit = 0.
      qfinl = 0.
      do nk = 1, ltop
      qinit = qinit + q0(nk)*ems(nk)
      qfinl = qfinl + qg(nk)*ems(nk)
      qfinl = qfinl + (qlg(nk)+qig(nk)+qrg(nk)+qsg(nk))*ems(nk)
      enddo
      qfinl = qfinl + trppt*ainc*timec*float(1-ishall)
      err2 = (qfinl-qinit)*100./qinit

      if (abs(err2).gt.0.05) then
      print*, '!! Moisture budget error in newconv !!'
!      stop 'newconv'
      goto 999
      endif

!  feedback to resolvable scale tendencies

      if (ishall.eq.1) timec = 5400.

      do k = 1, ltop
      dtdt  (jlon,jlat,k) = (tg(k)-t0(k))/timec
      dqdt  (jlon,jlat,k) = (qg(k)-q0(k))/timec
      dqcwdt(jlon,jlat,k) = max(qlg(k)/timec,0.)
      dqcidt(jlon,jlat,k) = max(qig(k)/timec,0.)
      dqpwdt(jlon,jlat,k) = max(qrg(k)/timec,0.)
      dqpi1dt(jlon,jlat,k)= max(qsg(k)/timec,0.)
        if(t0(k).gt.t00+1..and.dqcidt(jlon,jlat,k).gt.1.e-8) then
        dqcwdt(jlon,jlat,k) = dqcwdt(jlon,jlat,k) + dqcidt(jlon,jlat,k)
        dqcidt(jlon,jlat,k) = 0.
        endif
        if(t0(k).gt.t00+2..and.dqpi1dt(jlon,jlat,k).gt.1.e-8) then
        dqpwdt(jlon,jlat,k) = dqpwdt(jlon,jlat,k) + dqpi1dt(jlon,jlat,k)
        dqpi1dt(jlon,jlat,k) = 0.
        endif
      enddo

      precon = max(trppt*ainc/dxsq*dts*float(1-ishall),0.)   ! precipitation rate (mm/sec) x timestep

      if (t0(1).le.(t00+1.4).and.t0(2).le.(t00+0.2).and.t0(3).le.t00.and.t0(4).le.t00) then
      zds = min(t0(1)-t00-0.4, 1.)
      zds = max(zds, 0.)
      raicon(jlon,jlat) = zds*precon
      snocon(jlon,jlat) = (1.-zds)*precon
      else
      raicon(jlon,jlat) = precon
      snocon(jlon,jlat) = 0.
      endif

  999 continue   ! end loop on grid points

      return
      end subroutine convection_kf
!##################################################################################################################
      subroutine hliq (hw, z, p, t, q, ql, qi, qlnew, g, cp, cpv, rh)

! Computes thermodynamic properties of water, used by newconv

      parameter (aliq=611., bliq=17.4, cliq=bliq*273.15, dliq=33.65)

      ql0   = ql
      qi0   = qi
      qtot  = q+ql0+qi0
      cptot = cp+cpv*qtot
      hw1   = hw - (1.+qtot)*g*z

      t  = hw1/cptot
      if (t.lt.180.) then
      qs = 0.
      else
      es = aliq*exp((bliq*t-cliq)/(t-dliq))
      qs = rh*0.622*es/(p-es)
      endif

      if (qtot.le.qs) then                     ! unsaturated
      q  = qtot
      ql = 0.
      qi = 0.
      qlnew = 0.
      else                                     ! saturated
      rll = 2.5008e6+2369.*273.16
      rls = 2.8345e6+260. *273.16
      t1  = (hw1 +rll*ql0 +rls*qi0)/(cptot+2369.*ql0+260.*qi0)
      es  = aliq*exp((bliq*t1-cliq)/(t1-dliq))
      qs1 = rh*0.622*es/(p-es)

        if (qtot.gt.ql0+qi0+qs1) then          ! case in which new condensate is generated
          if (rh.lt.0.99) then                 ! return without condensation
          t = t1
          return
          endif
        rll = 2.5008e6-2369.*(t1-273.16)
        rls = 2.8345e6-260. *(t1-273.16)
        res1 = hw1 - cptot*t1 + rls*qi0 + rll*(qtot-qi0-qs1)
        t = t1+1.  ! first guess temperature
 10     es  = aliq*exp((bliq*t-cliq)/(t-dliq))
        q = 0.622*es/(p-es)
        rll = 2.5008e6-2369.*(t-273.16)
        rls = 2.8345e6-260. *(t-273.16)
        res = hw1 - cptot*t + rls*qi0 + rll*(qtot-qi0-q)
        deriv = (res-res1)/(t-t1)
        t1  = t
        res1 = res
        t = t - res/deriv   ! newton step
          if(abs(t-t1).gt.1.e-2) then
          goto 10
          else
          es = aliq*exp((bliq*t-cliq)/(t-dliq))
          q = 0.622*es/(p-es)
          qlnew = qtot-ql0-qi0-q
          endif
        else                          ! case in which it evaporates to saturation
        fract = ql0/(ql0+qi0+1.e-15)
        qlnew = 0.
        rll = 2.5008e6-2369.*(t1-273.16)
        rls = 2.8345e6-260. *(t1-273.16)
        res1 = hw1 - cptot*t1 + rll*(qtot-qs1)*fract + rls*(qtot-qs1)*(1.-fract)
        t = t1+1.  ! first guess temperature
 11     es  = aliq*exp((bliq*t-cliq)/(t-dliq))
        q = rh*0.622*es/(p-es)
        rll = 2.5008e6-2369.*(t-273.16)
        rls = 2.8345e6-260. *(t-273.16)
        res = hw1 - cptot*t + rll*(qtot-q)*fract + rls*(qtot-q)*(1.-fract)
        deriv = (res-res1)/(t-t1)
        t1  = t
        res1 = res
        t = t - res/deriv   ! Newton step
          if(abs(t-t1).gt.1.e-2) then
          goto 11
          else
          es = aliq*exp((bliq*t-cliq)/(t-dliq))
          q = rh*0.622*es/(p-es)
          ql = (qtot-q)*fract
          qi = (qtot-q)*(1.-fract)
          endif
        endif
      endif

      return
      end subroutine hliq
!###############################################################################################################
     subroutine divdamp (ddamp)

     use mod_moloch, only: nlon, nlat, nlev, nlonm1, nlatm1, ip_e, ip_n, ip_s, ip_w, u, v, div2, dx, dy, dt, hxt
     implicit none
     integer jlon, jlat, jklev
     real ddamp, ddamp1, zprof, zdtrdx, zdtrdy, p2(nlon,nlat)

#ifdef mpi
       call u_ghost (div2(2,:     ,:), ip_w, div2(nlon,:,:), ip_e, nlat*nlev)
       call u_ghost (div2(nlonm1,:,:), ip_e, div2(1   ,:,:), ip_w, nlat*nlev)  ! manda questo, li', ricevi quello, da la'
       call u_ghost (div2(:,nlatm1,:), ip_n, div2(:   ,1,:), ip_s, nlon*nlev)
       call u_ghost (div2(:,2     ,:), ip_s, div2(:,nlat,:), ip_n, nlon*nlev)
#endif

     do 50 jklev = 1, nlev
     zprof  = ddamp + (1.-ddamp)/(nlev-jklev+3.)

     ddamp1 = zprof*.125*dy**2/dt
     zdtrdy = ddamp1/dy
     do jlat = 2, nlatm1
     zdtrdx = ddamp1/(dx*hxt(jlat))
     do jlon = 2, nlonm1
     u(jlon,jlat,jklev) = u(jlon,jlat,jklev) + zdtrdx*(div2(jlon+1,jlat,jklev)-div2(jlon,jlat,jklev))
     v(jlon,jlat,jklev) = v(jlon,jlat,jklev) + zdtrdy*(div2(jlon,jlat,jklev)-div2(jlon,jlat-1,jklev))
     enddo
     enddo

!  Horizontal diffusion of divergence

     do jlat = 2, nlatm1
     do jlon = 2, nlonm1
     p2(jlon,jlat) = .125 *(div2(jlon,jlat-1,jklev) + div2(jlon-1,jlat,jklev) + &
                            div2(jlon+1,jlat,jklev) + div2(jlon,jlat+1,jklev)) -.5*div2(jlon,jlat,jklev)
     enddo
     enddo

     do jlat = 2, nlatm1
     do jlon = 2, nlonm1
     div2(jlon,jlat,jklev) = div2(jlon,jlat,jklev) + zprof*p2(jlon,jlat)
     enddo
     enddo

50   continue

     return
     end
!###############################################################################################################
subroutine nudging_tuv(kstep, time_data)
!
! Large scale nudging for T, U, V variables
!
! Data fields for nudging are the same used for boundary contidion definition.
! Method of nudging: Fourier transformation until wave number km (x axis) and
! lm (y axis), smooth transition is applied over the truncating point.
! Period of relaxation of nudging is defined by nudging_time variable (about 4 hours).

use mod_moloch, only : t, tb1, tb2, u, ub1, ub2, v, vb1, vb2, zeta, phig, &
 dtstep, nudging_time, ntnudg, nlon, nlat, nlev
implicit none

integer :: kstep, jlon, jlat, jklev
real :: time_data, cnudg
real, save :: cnudg0
real, dimension(nlon,nlat), save :: znudg_t, znudg_u, znudg_v

 if (kstep == 1) then
   call lpfilts (znudg_t, 0)
   cnudg0 = ntnudg*dtstep/nudging_time    ! nudging coefficient "tau", nudging_time = 4 hours 
 endif

 do jklev = 3, nlev
   znudg_t(:,:) = t(:,:,jklev) - (tb1(:,:,jklev) + time_data*(tb2(:,:,jklev)-tb1(:,:,jklev)))
   call lpfilts (znudg_t, 1)
   znudg_u(:,:) = u(:,:,jklev) - (ub1(:,:,jklev) + time_data*(ub2(:,:,jklev)-ub1(:,:,jklev)))
   call lpfilts (znudg_u, 1)
   znudg_v(:,:) = v(:,:,jklev) - (vb1(:,:,jklev) + time_data*(vb2(:,:,jklev)-vb1(:,:,jklev)))
   call lpfilts (znudg_v, 1)
   do jlat = 1, nlat
   do jlon = 1, nlon
     cnudg = cnudg0* min ((zeta(jlon,jlat,jklev)/1500.)**2, 1.)   ! nudging coefficient "f1"
     cnudg = cnudg * max (1.-abs(phig(jlon,jlat))/25000., 0.)     ! nudging coefficient "f2"
     t(jlon,jlat,jklev) = t(jlon,jlat,jklev) - cnudg*znudg_t(jlon,jlat)
     u(jlon,jlat,jklev) = u(jlon,jlat,jklev) - cnudg*znudg_u(jlon,jlat)
     v(jlon,jlat,jklev) = v(jlon,jlat,jklev) - cnudg*znudg_v(jlon,jlat)
   enddo
   enddo
 enddo

end subroutine nudging_tuv
!###############################################################################################################
    subroutine lpfilts (f, iflag)

! Low-pass filter in two dimensions (valid for zero boundary) with MPI (July 2021)
! f(nlon,nlat): input and output
! km,lm: number of spectral components in x and y directions
! bvx,bvy: basis vectors
! px,py: convolution function. Pure spectral truncation when px=py=1
!        (default is set to transform of Green function of diffusion eq.)
! sx,sy: spectral sine coefficients
! First call with iflag=0 to initialize basis vectors

use mod_moloch, only : nlon, nlat, gnlon, gnlat, km, lm, myid, comm_row, comm_col, nprocsy
#ifdef mpi
  use mod_moloch, only : comm_row, comm_col, row_color, col_color
  include 'mpif.h'
#endif

    real f(nlon,nlat), sx(nlat,2*km), sxg(nlat,2*km), sy(nlon,2*lm), syg(nlon,2*lm)
    real bvx(nlon,2*km), bvy(nlat,2*lm), px(2*km), py(2*lm)
    save bvx, bvy

    if (iflag.eq.0) then ! initialization
    do k = 1, 2*km
    px(k) = exp(-(float(k)/float(km))**2)
    enddo
    do l = 1, 2*lm
    py(l) = exp(-(float(l)/float(lm))**2)
    enddo
    infx = 1 + (myid/nprocsy)*(nlon-2)
    infy = 1 + (myid-(myid/nprocsy)*nprocsy)*(nlat-2)
    dx = 3.1415927/float(gnlon-1)
    dy = 3.1415927/float(gnlat-1)
    do k = 1, 2*km
    do jlon = 1, nlon
    bvx(jlon,k) = sqrt(2./float(gnlon-1)*px(k))*sin(k*(infx+jlon-2)*dx)
    enddo
    enddo
    do l = 1, 2*lm
    do jlat = 1, nlat
    bvy(jlat,l) = sqrt(2./float(gnlat-1)*py(l))*sin(l*(infy+jlat-2)*dy)
    enddo
    enddo
    if (myid == 0) print*, ' ***** Nudging has been initialized with km,lm = ', km, lm
    return
    endif ! end initialization

    do k = 1, 2*km
    do jlat = 1, nlat  ! x-transform
    sx(jlat,k) = 0.
      do jlon = 2, nlon-1
      sx(jlat,k) = sx(jlat,k) + f(jlon,jlat)*bvx(jlon,k)
      enddo
    enddo
    enddo

    if (nlon.eq.gnlon) then
    sxg = sx
    else
#ifdef mpi
      call mpi_allreduce (sx, sxg, 2*km*nlat, mpi_real, mpi_sum, comm_row, ierr)
#endif
    endif

    f = 0.
    do k = 1, 2*km
    do jlat = 1, nlat
    do jlon = 1, nlon
    f(jlon,jlat) = f(jlon,jlat) + sxg(jlat,k)*bvx(jlon,k)
    enddo
    enddo
    enddo

    do l = 1, 2*lm
    do jlon = 1, nlon  ! y-transform
    sy(jlon,l) = 0.
      do jlat = 2, nlat-1
      sy(jlon,l) = sy(jlon,l) + f(jlon,jlat)*bvy(jlat,l)
      enddo
    enddo
    enddo

    if (nlat.eq.gnlat) then
    syg = sy
    else
#ifdef mpi
      call mpi_allreduce (sy, syg, 2*lm*nlon, mpi_real, mpi_sum, comm_col, ierr)
#endif
    endif

    f = 0.
    do l = 1, 2*lm
    do jlat = 1, nlat
    do jlon = 1, nlon
    f(jlon,jlat) = f(jlon,jlat) + syg(jlon,l)*bvy(jlat,l)
    enddo
    enddo
    enddo

    return
    end
!###############################################################################################################
subroutine snowfall_level

! Definition of atmospheric level index (from bottom to top) when
! ice precipitation predominate liquid precipitation
! to output for thje aim of definition of snowfall altitude definition
! by prostprocessing elaboration.
!
! Approximating hypothesis: hydrometeors terminal velocity, that must be 
! considereted in precipitation flow determination, is proportional
! of specific content of precipitation (q), because of it is possible
! to calculate the precipitation flow at each atmospheric levels for
! convection precipitation.
!
! level_snowfall(nlon,nlat): output
! qcw(nlon,nlat,nlev), qci(nlon,nlat,nlev), qpw(nlon,nlat,nlev), qpi1(nlon,nlat,nlev) and qpi2(nlon,nlat,nlev): input

use mod_moloch, only : level_snowfall, qpw, qpi1, qpi2, nlon, nlat, nlev, myid

integer :: jlon, jlat, jklev
real :: zqprec, zqprec_ice, zfrac

 level_snowfall(:,:) = 0

 do jlat = 1,nlat
 do jlon = 1,nlon

  do jklev = 1,nlev

    zqprec = qpw(jlon,jlat,jklev) + qpi1(jlon,jlat,jklev) + qpi2(jlon,jlat,jklev)
    if (zqprec > 1.e-8) then
      zqprec_ice = qpi1(jlon,jlat,jklev) + qpi2(jlon,jlat,jklev)
      zfrac = zqprec_ice / zqprec
      if (zfrac > 0.5) then
        level_snowfall(jlon,jlat) = jklev
        exit
      endif
    endif

  enddo

 enddo
 enddo

return
end
!###############################################################################################################
    subroutine raglio (x0,y0, x00,y00, dlon,dlat, xc,yc, nxg,nyg, nx,ny, i1,j1)

! X00, Y00 primo punto in coordinate ruotate
! X0, Y0 centro delle coordinate ruotate in gradi
! xc, yc coordinate geografiche da trasformare in coordinate ruotate

    real*8 zpi, zfac, zx0, zy0, zx, zy, zzlat, zz

    if (abs(x0)>0.01.or.abs(y0)>0.01) then ! Case of rotated grid

    zpi  = dabs(dacos(-1.d0))
    zfac = zpi/180.d0
      zx0  = dble(x0)*zfac
      zy0  = dble(y0)*zfac
      zx   = dble(xc)*zfac
      zy   = dble(yc)*zfac
      if (zx-zx0.gt. zpi) zx = zx - 2.d0*zpi
      if (zx-zx0.lt.-zpi) zx = zx + 2.d0*zpi

      zzlat = dasin( -dcos(zy)*dsin(zy0)*dcos(zx-zx0) + dsin(zy)*dcos(zy0) )
      zz = (dsin(zy)-dcos(zy0)*dsin(zzlat))/(dsin(zy0)*dcos(zzlat))
      if (zz < -1.d0.and.zz > -1.00001d0) zz = -1.d0
      if (zz >  1.d0.and.zz <  1.00001d0) zz =  1.d0
        if (zx < zx0) then
        zrlon = -dacos(zz)/zfac
        else
        zrlon =  dacos(zz)/zfac
        endif
      zrlat = zzlat/zfac

      if (zrlon >  180.) zrlon = zrlon - 360.
      if (zrlon < -180.) zrlon = zrlon + 360.

      else ! Case of non rotated grid

      zrlon = xc
      zrlat = yc

      endif

      ic = (zrlon-x00)/dlon
      jc = (zrlat-y00)/dlat
      if (ic+nx/2.gt.nyg) ic = nxg-nx/2-1
      if (jc+ny/2.gt.nyg) jc = nyg-ny/2-1
      i1 = max(1, ic-nx/2)
      j1 = max(1, jc-ny/2)

      return
      end
!###############################################################################################################
subroutine wrspray (i1, j1, inst)

use mod_moloch
implicit none

integer :: iunit=43, iunit_work=29, inst, i1, j1, i2, j2, jlon, jlat, jklev
real(4) zwork(nlon,nlat), zlon00, zlat00
character(len=30) :: file_out="spray_000.mhf"

 write (file_out(7:9),'(i3.3)') inst

 if (myid == 0) then

      open (iunit, file=trim(file_out), form='unformatted', status='unknown')

      i2 = i1 + nxspray-1
      j2 = j1 + nyspray-1
      zlon00 = pdr(5) + (i1-1)*dlon
      zlat00 = pdr(4) + (j1-1)*dlat +dlat/2.  ! Per spray solo output su punti T
      write(iunit) nxspray, nyspray, nzspray+1
      write(iunit) nfdr(5:12)                                         ! data iniziale e validita'
      write(iunit) a0, b0, h, x0, y0, dlon, dlat, dz, zlon00, zlat00  ! coordinate verticali e griglia orizzontale

 endif


      call collect (phig, gfield)  ! surface geopotential
      if (myid == 0) then
      do jlat = j1, j2
      write(iunit) (gfield(jlon,jlat), jlon = i1, i2)
      enddo
      endif

      call collect (fmask, gfield) ! land sea mask (sea=1)
      if (myid == 0) then
      do jlat = j1, j2
      write(iunit) (gfield(jlon,jlat), jlon = i1, i2)
      enddo
      endif

      call collect (alont, gfield) ! longitudes (degrees)
      if (myid == 0) then
      do jlat = j1, j2
      write(iunit) (gfield(jlon,jlat), jlon = i1, i2)
      enddo
      endif

      call collect (alatt, gfield)  ! latitudes (degrees)
      if (myid == 0) then
      do jlat = j1, j2
      write(iunit) (gfield(jlon,jlat), jlon = i1, i2)
      enddo
      endif

      call collect (ustar, gfield)  ! U*
      if (myid == 0) then
      do jlat = j1, j2
      write(iunit) (gfield(jlon,jlat), jlon = i1, i2)
      enddo
      endif

      call collect (tstar, gfield)  ! T*
      if (myid == 0) then
      do jlat = j1, j2
      write(iunit) (gfield(jlon,jlat), jlon = i1, i2)
      enddo
      endif

      call collect (qstar, gfield)  ! Q*
      if (myid == 0) then
      do jlat = j1, j2
      write(iunit) (gfield(jlon,jlat), jlon = i1, i2)
      enddo
      endif

      call collect (rgmd, gfield)  ! momentum roughness
      if (myid == 0) then
      do jlat = j1, j2
      write(iunit) (gfield(jlon,jlat), jlon = i1, i2)
      enddo
      endif

      call collect (tg(1,1,1), gfield)  ! sst and soil temperature
      if (myid == 0) then
      do jlat = j1, j2
      write(iunit) (gfield(jlon,jlat), jlon = i1, i2)
      enddo
      endif

      call collect (tskin, gfield)
      if (myid == 0) then
      do jlat = j1, j2
      write(iunit) (gfield(jlon,jlat), jlon = i1, i2)
      enddo
      endif
      do jklev = 1, nzspray
      call collect (t(1,1,jklev), gfield)
      if (myid == 0) then
      do jlat = j1, j2
      write(iunit) (gfield(jlon,jlat), jlon = i1, i2)
      enddo
      endif
      enddo

      call collect (ps, gfield)  ! pressure
      if (myid == 0) then
      do jlat = j1, j2
      write(iunit) (gfield(jlon,jlat), jlon = i1, i2)
      enddo
      endif
      do jklev = 1, nzspray
      call collect (p(1,1,jklev), gfield)
      if (myid == 0) then
      do jlat = j1, j2
      write(iunit) (gfield(jlon,jlat), jlon = i1, i2)
      enddo
      endif
      enddo

      call collect (qskin, gfield)  ! qskin
      if (myid == 0) then
      do jlat = j1, j2
      write(iunit) (gfield(jlon,jlat), jlon = i1, i2)
      enddo
      endif
      do jklev = 1, nzspray
      call collect (q(1,1,jklev), gfield)  ! specific humidity
      if (myid == 0) then
      do jlat = j1, j2
      write(iunit) (gfield(jlon,jlat), jlon = i1, i2)
      enddo
      endif
      enddo

      if (myid == 0) then  ! wind
      gfield = 0.
      do jlat = j1, j2
      write(iunit) (gfield(jlon,jlat), jlon = i1, i2)
      enddo
      endif
      do jklev = 1, nzspray
      call collect (u(1,1,jklev), gfield)
      if (myid == 0) then
      do jlat = j1, j2
      write(iunit) (.5*(gfield(jlon,jlat)+gfield(jlon-1,jlat)), jlon = i1, i2)
      enddo
      endif
      enddo

      if (myid == 0) then
      gfield = 0.
      do jlat = j1, j2
      write(iunit) (gfield(jlon,jlat), jlon = i1, i2)
      enddo
      endif
      do jklev = 1, nzspray
      call collect (v(1,1,jklev), gfield)
      if (myid == 0) then
      do jlat = j1, j2
      write(iunit) (.5*(gfield(jlon,jlat)+gfield(jlon,jlat+1)), jlon = i1, i2)
      enddo
      endif
      enddo

      call collect (w(1,1,1), gfield)
      if (myid == 0) then
      do jlat = j1, j2
      write(iunit) (gfield(jlon,jlat), jlon = i1, i2)
      enddo
      endif
      do jklev = 1, nzspray
      zwork(:,:) = .5*(w(:,:,jklev)+w(:,:,jklev+1))
      call collect (zwork, gfield)  ! vertical velocity
      if (myid == 0) then
      do jlat = j1, j2
      write(iunit) (gfield(jlon,jlat), jlon = i1, i2)
      enddo
      endif
      enddo

      call collect (tke(1,1,1), gfield)
      if (myid == 0) then
      do jlat = j1, j2
      write(iunit) (gfield(jlon,jlat), jlon = i1, i2)
      enddo
      endif
      do jklev = 1, nzspray
      zwork(:,:) = .5*(tke(:,:,jklev)+tke(:,:,jklev+1))
      call collect (zwork, gfield)  ! turbulent kinetic energy
      if (myid == 0) then
      do jlat = j1, j2
      write(iunit) (gfield(jlon,jlat), jlon = i1, i2)
      enddo
      endif
      enddo

      call collect (cvm(1,1,1), gfield)
      if (myid == 0) then
      do jlat = j1, j2
      write(iunit) (gfield(jlon,jlat), jlon = i1, i2)
      enddo
      endif
      do jklev = 1, nzspray
      zwork(:,:) = .5*(cvm(:,:,jklev)+cvm(:,:,jklev+1))
      call collect (zwork, gfield)  ! coefficient of vertical diffusion
      if (myid == 0) then
      do jlat = j1, j2
      write(iunit) (gfield(jlon,jlat), jlon = i1, i2)
      enddo
      endif
      enddo

      if (myid == 0) then
      gfield = 0.
      do jlat = j1, j2
      write(iunit) (gfield(jlon,jlat), jlon = i1, i2)
      enddo
      endif
      do jklev = 1, nzspray
      call collect (chm(1,1,jklev), gfield)  ! coefficient of horizontal diffusion
      if (myid == 0) then
      do jlat = j1, j2
      write(iunit) (gfield(jlon,jlat), jlon = i1, i2)
      enddo
      endif
      enddo

 if (myid == 0) then
     close (iunit)
     open (iunit_work, file=trim(file_out)//'.txt', status='unknown')
     write (iunit_work,'(2a)') trim(file_out),' is full and closed'
     close (iunit_work)
     print *
     print *,'  Output written on file ', trim(file_out)
 endif

return
end subroutine wrspray
!##################################################################################################################
subroutine rdrf(istart)
! Reads restart file

use mod_moloch
implicit none

integer :: iunit=17, istart, ierr_open, jlon, jlat, jklev, comm, error, jpr, tag=0, itype
real :: zdtstep_old, ztime
integer, dimension(gnlon,gnlat)      :: igfield
character(len=30) :: file_in="moloch.rf"
integer, dimension(1) :: i1

#ifdef mpi
      include 'mpif.h'
      integer :: status(mpi_status_size)

      comm = mpi_comm_world
#endif

 if(myid == 0) then

   open (iunit, file=file_in, form='unformatted', status='old', iostat=ierr_open)

   if (ierr_open /= 0) then
     close (iunit)
     print *,'Restart, error in ',trim(file_in),' opening, stop'
#ifdef mpi
        call mpi_abort(comm, error)
#endif
     stop
   endif

   read (iunit) nstep0
   read (iunit) nfdr0
   read (iunit) pdr0
   read (iunit) nfdr
   read (iunit) pdr
   read (iunit) reqtim, rdecli

   if (istart /= 0) close (iunit)
    
 endif ! myid

#ifdef mpi
     call mpi_bcast(nstep0,        1, mpi_integer,  0, comm, error)   
     call mpi_bcast(nfdr0(1:50),  50, mpi_integer,  0, comm, error)   
     call mpi_bcast(pdr0(1:100), 100, mpi_real,     0, comm, error)   
     call mpi_bcast(nfdr(1:50),   50, mpi_integer,  0, comm, error)   
     call mpi_bcast(pdr(1:100),  100, mpi_real,     0, comm, error)   
     call mpi_bcast(reqtim,        1, mpi_double,   0, comm, error)   
     call mpi_bcast(rdecli,        1, mpi_double,   0, comm, error)   
#endif

 if (istart /= 0) return

! Basic atmospheric variables

 do jklev=1,nlev
   if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, p(1,1,jklev))
 enddo

 do jklev=1,nlev
   if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, u(1,1,jklev))
 enddo

 do jklev=1,nlev
   if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, v(1,1,jklev))
 enddo

 do jklev=1,nlevp1
   if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, w(1,1,jklev))
 enddo

 do jklev=1,nlev
   if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, t(1,1,jklev))
 enddo

 do jklev=1,nlev
   if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, q(1,1,jklev))
 enddo

 do jklev=1,nlev
   if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, qcw(1,1,jklev))
 enddo

 do jklev=1,nlev
   if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, qci(1,1,jklev))
 enddo

 do jklev=1,nlev
   if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, qpw(1,1,jklev))
 enddo

 do jklev=1,nlev
   if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, qpi1(1,1,jklev))
 enddo

 do jklev=1,nlev
   if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, qpi2(1,1,jklev))
 enddo

! Secondary atmospheric variables

 do jklev=1,nlev
   if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, pai(1,1,jklev))
 enddo

 do jklev=1,nlevp1
   if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, s(1,1,jklev))
 enddo

 do jklev=1,nlevp1
   if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, wwkw(1,1,jklev))
 enddo

 do jklev=1,nlevp1
   if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, deltaw(1,1,jklev))
 enddo

 do jklev=1,nlev
   if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, div2(1,1,jklev))
 enddo

 do jklev=1,nlev
   if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, ncw(1,1,jklev))
 enddo

 do jklev=1,nlev
   if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, nci(1,1,jklev))
 enddo

 do jklev=1,nlev
   if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, fcloud(1,1,jklev))
 enddo

 do jklev=1,nlevp1
   if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, tke(1,1,jklev))
 enddo

 do jklev=1,nlevp1
   if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, mlz(1,1,jklev))
 enddo

 do jklev=1,nlevp1
   if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, rich(1,1,jklev))
 enddo

 do jklev=1,nlev
   if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, dtdt(1,1,jklev))
 enddo

 do jklev=1,nlev
   if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, geldt(1,1,jklev))
 enddo

 do jklev=1,nlev
   if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, corrdt(1,1,jklev))
 enddo

 do jklev=1,nlev
   if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, rradar(1,1,jklev))
 enddo

! Surface prognostic and accumulation variables 

 if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, lai)

 if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, fveg)

 if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, rgm)

 if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, rgq)

 if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, fice)

 if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, iceth)

 if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, albedo)

 if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, emisg1)

 if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, emisg2)

 if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, cloudt)

 if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, prectot)

 if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, precconv)

 if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, precsolid)

 if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, tskin)

 if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, tg_surf)

 do jklev=1,nlevg
   if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, tg(1,1,jklev))
 enddo

 if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, qskin)

 if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, qg_surf)

 do jklev=1,nlevg
   if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, qg(1,1,jklev))
 enddo

 if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, fice_soil_surf)

 do jklev=1,nlevg
   if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, fice_soil(1,1,jklev))
 enddo

 if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, snow)

 do jklev=1,nlev_snow
   if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, snow_lev(1,1,jklev))
 enddo

 do jklev=1,nlev_snow
   if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, snow_t(1,1,jklev))
 enddo

 do jklev=1,nlev_snow
   if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, snow_fice(1,1,jklev))
 enddo

 do jklev=1,nlev_snow
   if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, snow_age(1,1,jklev))
 enddo

 do jklev=1,nlev_snow
   if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, snow_melt_age(1,1,jklev))
 enddo

 do jklev=1,nlev_snow
   if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, snow_dens(1,1,jklev))
 enddo

 if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, alsn)

 if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, fsnow)

 if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, roscdm)

 if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, roscdt)

 if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, cswfl)

 if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, clwfl)

 if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, cshflux)

 if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, clhflux)

 if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, t2min)

 if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, t2max)

 if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, ws10max)

 if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, runoff)

 if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, runoff_tot)

 if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, cwvflux)

 if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, cfl_heat_soil_bottom)

 if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, cfl_water_soil_bottom)

 if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, frvis)

 if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, frirr)

 if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, gelvis)

 if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, gelirr)

 if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, corvis)

 if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, corirr)

 if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, raini)

 if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, snowi)

 if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, hflux)

 if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, qflux)

 if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, totsky)

 if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, soldir)

 if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, shf_accum)

 if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, lhf_accum)

 if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, qf_accum)

 if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, slopeff)

!  Physical parameters of radiation scheme

 do jklev=1,nlev
   if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, ozon(1,1,jklev))
 enddo

 do itype=1,ntype_aerosol
 do jklev=1,nlev
   if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, aerosol(1,1,jklev,itype))
 enddo
 enddo

 do jklev=1,nlev
   if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, aerotot(1,1,jklev))
 enddo

!  Physical parameters of soil scheme

 if(myid == 0) call rrec2_int (iunit, gnlon, gnlat, igfield)
 call disperse_int (igfield, mask_soil)

 if(myid == 0) call rrec2_int (iunit, gnlon, gnlat, igfield)
 call disperse_int (igfield, ind_lev_soil_h_bottom)

 if(myid == 0) call rrec2_int (iunit, gnlon, gnlat, igfield)
 call disperse_int (igfield, ind_lev_soil_w_bottom)

 if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, kturb_surf_m)

 if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, kturb_surf_h)

 if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, kturb_surf_q)

 do jklev = 1, 20
   if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, kturb_surf_h_mem(1,1,jklev))
 enddo

 do jklev = 1, 20
   if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, kturb_surf_q_mem(1,1,jklev))
 enddo

 if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, flh_specif)

 if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, fl_wv)

 if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, fl_heat_soil_bottom)

 if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, fl_water_soil_bottom)

 if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, fl_runoff)

 if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, fl_runoff_tot)

 do jklev = 1, nlevg
   if(myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, psi_soil(1,1,jklev))
 enddo

 if(myid == 0) then
   close (iunit)
   print *
   print*,trim(file_in),' read, nstep0= ',nstep0
 endif

! ---> pochva
    mask_soil(1,   1:nlat) = 0
    mask_soil(nlon,1:nlat) = 0
    mask_soil(1:nlon,   1) = 0
    mask_soil(1:nlon,nlat) = 0
! <---

return
end subroutine rdrf
!###############################################################################################################
subroutine wrrf(kstep)
! Writes restart file

use mod_moloch
implicit none

integer :: iunit=27, kstep, jlon, jlat, jklev, i1, i2, i3, &
 ilon1=1, ilon2=gnlon, jlat1=1, jlat2=gnlat, flag_lon=0, itype
integer, dimension(gnlon,gnlat)      :: igfield
character(len=30) :: file_out

 if(myid.eq.0) then

   irf=irf+1
   i1=irf/10
   i2=i1*10
   i3=irf-i2
   write (file_out,'(A,I2.2,A)') 'moloch_',i3,'.rf'

   open (iunit, file=trim(file_out), form='unformatted', status='unknown')

   write (iunit) kstep+1
   write (iunit) nfdr0
   write (iunit) pdr0
   write (iunit) nfdr
   write (iunit) pdr
   write (iunit) reqtim, rdecli

 endif

! Basic atmospheric variables

 do jklev=1,nlev
   call collect (p(1,1,jklev), gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)
 enddo

 do jklev=1,nlev
   call collect (u(1,1,jklev), gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)
 enddo

 do jklev=1,nlev
   call collect (v(1,1,jklev), gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)
 enddo

 do jklev=1,nlevp1
   call collect (w(1,1,jklev), gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)
 enddo

 do jklev=1,nlev
   call collect (t(1,1,jklev), gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)
 enddo

 do jklev=1,nlev
   call collect (q(1,1,jklev), gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)
 enddo

 do jklev=1,nlev
   call collect (qcw(1,1,jklev), gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)
 enddo

 do jklev=1,nlev
   call collect (qci(1,1,jklev), gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)
 enddo

 do jklev=1,nlev
   call collect (qpw(1,1,jklev), gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)
 enddo

 do jklev=1,nlev
   call collect (qpi1(1,1,jklev), gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)
 enddo

 do jklev=1,nlev
   call collect (qpi2(1,1,jklev), gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)
 enddo

! Secondary atmospheric variables

 do jklev=1,nlev
   call collect (pai(1,1,jklev), gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)
 enddo

 do jklev=1,nlevp1
   call collect (s(1,1,jklev), gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)
 enddo

 do jklev=1,nlevp1
   call collect (wwkw(1,1,jklev), gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)
 enddo

 do jklev=1,nlevp1
   call collect (deltaw(1,1,jklev), gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)
 enddo

 do jklev=1,nlev
   call collect (div2(1,1,jklev), gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)
 enddo

 do jklev=1,nlev
   call collect (ncw(1,1,jklev), gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)
 enddo

 do jklev=1,nlev
   call collect (nci(1,1,jklev), gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)
 enddo

 do jklev=1,nlev
   call collect (fcloud(1,1,jklev), gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)
 enddo

 do jklev=1,nlevp1
   call collect (tke(1,1,jklev), gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)
 enddo

 do jklev=1,nlevp1
   call collect (mlz(1,1,jklev), gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)
 enddo

 do jklev=1,nlevp1
   call collect (rich(1,1,jklev), gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)
 enddo

 do jklev=1,nlev
   call collect (dtdt(1,1,jklev), gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)
 enddo

 do jklev=1,nlev
   call collect (geldt(1,1,jklev), gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)
 enddo

 do jklev=1,nlev
   call collect (corrdt(1,1,jklev), gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)
 enddo

 do jklev=1,nlev
   call collect (rradar(1,1,jklev), gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)
 enddo

! Surface prognostic and accumulation variables

 call collect (lai, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

 call collect (fveg, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

 call collect (rgm, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

 call collect (rgq, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

 call collect (fice, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

 call collect (iceth, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

 call collect (albedo, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

 call collect (emisg1, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

 call collect (emisg2, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

 call collect (cloudt, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

 call collect (prectot, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

 call collect (precconv, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

 call collect (precsolid, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

 call collect (tskin, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

 call collect (tg_surf, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

 do jklev=1,nlevg
   call collect (tg(1,1,jklev), gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)
 enddo

 call collect (qskin, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

 call collect (qg_surf, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

 do jklev=1,nlevg
   call collect (qg(1,1,jklev), gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)
 enddo

 call collect (fice_soil_surf, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

 do jklev=1,nlevg
   call collect (fice_soil(1,1,jklev), gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)
 enddo

 call collect (snow, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

 do jklev=1,nlev_snow
   call collect (snow_lev(1,1,jklev), gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)
 enddo

 do jklev=1,nlev_snow
   call collect (snow_t(1,1,jklev), gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)
 enddo

 do jklev=1,nlev_snow
   call collect (snow_fice(1,1,jklev), gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)
 enddo

 do jklev=1,nlev_snow
   call collect (snow_age(1,1,jklev), gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)
 enddo

 do jklev=1,nlev_snow
   call collect (snow_melt_age(1,1,jklev), gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)
 enddo

 do jklev=1,nlev_snow
   call collect (snow_dens(1,1,jklev), gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)
 enddo

 call collect (alsn, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

 call collect (fsnow, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

 call collect (roscdm, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

 call collect (roscdt, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

 call collect (cswfl, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

 call collect (clwfl, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

 call collect (cshflux, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

 call collect (clhflux, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

 call collect (t2min, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

 call collect (t2max, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

 call collect (ws10max, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

 call collect (runoff, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

 call collect (runoff_tot, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

 call collect (cwvflux, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

 call collect (cfl_heat_soil_bottom, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

 call collect (cfl_water_soil_bottom, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

 call collect (frvis, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

 call collect (frirr, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

 call collect (gelvis, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

 call collect (gelirr, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

 call collect (corvis, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

 call collect (corirr, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

 call collect (raini, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

 call collect (snowi, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

 call collect (hflux, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

 call collect (qflux, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

 call collect (totsky, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

 call collect (soldir, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

 call collect (shf_accum, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

 call collect (lhf_accum, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

 call collect (qf_accum, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

 call collect (slopeff, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

 do jklev=1,nlev
   call collect (ozon(1,1,jklev), gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)
 enddo

 do itype=1,ntype_aerosol
 do jklev=1,nlev
   call collect (aerosol(1,1,jklev,itype), gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)
 enddo
 enddo

 do jklev=1,nlev
   call collect (aerotot(1,1,jklev), gfield)
   if (myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)
 enddo

!  Physical parameters of soil scheme

 call collect_int (mask_soil, igfield)
 if(myid == 0) call wrec2_int (iunit, gnlon, gnlat, igfield)

 call collect_int (ind_lev_soil_h_bottom, igfield)
 if(myid == 0) call wrec2_int (iunit, gnlon, gnlat, igfield)

 call collect_int (ind_lev_soil_w_bottom, igfield)
 if(myid == 0) call wrec2_int (iunit, gnlon, gnlat, igfield)

 call collect (kturb_surf_m, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

 call collect (kturb_surf_h, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

 call collect (kturb_surf_q, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

 do jklev = 1, 20
   call collect (kturb_surf_h_mem(1,1,jklev), gfield)
   if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)
 enddo

 do jklev = 1, 20
   call collect (kturb_surf_q_mem(1,1,jklev), gfield)
   if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)
 enddo

 call collect (flh_specif, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

 call collect (fl_wv, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

 call collect (fl_heat_soil_bottom, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

 call collect (fl_water_soil_bottom, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

 call collect (fl_runoff, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

 call collect (fl_runoff_tot, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)

 do jklev = 1, nlevg
   call collect (psi_soil(1,1,jklev), gfield)
   if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield)
 enddo

 if(myid == 0) then
   close (iunit)
   print *
   print *,trim(file_out),' written at step ',kstep
 endif

return
end subroutine wrrf
!###############################################################################################################
SUBROUTINE RDRF_POCHVA(NSTEP0, NPROCX, NPROCY, PROC_IND, FIELD_GLOB_R, NX_GLOB, NY_GLOB, DTSTEP_EXT)
! Write reastart file with pochva variables

USE MODULE_POCHVA, ONLY : DTSTEP, NPOINT_X, NPOINT_Y, NLEV_SOIL_MAX, NLEV_SNOW, &
 NLEV_SOIL_HEAT,  NLEV_SOIL_WATER, LEV_SOIL, LEV_SOIL_H, D_LEV_SOIL, D_LEV_SOIL_H, SUM_D_LEV_SOIL_H, &
 TS, QS, QS_REL, PSIS, KAPPA_H, FRAC_SICE, DFRAC_SICE_DT, QSURF, TSURF, TFOREST, QFOREST, QS_REL_VEG, &
 D_LEV_SNOW, TSNOW, LEV_SNOW, FICE_SNOW, RHOSNOW, SNOW_AGE, SNOW_MELT_AGE, DLEV_SNOW_BOTTOM, FLAG_SNOW_THICK, SNOW_DIRT, &
 QS_OLD, QSI, DQS, RUN_OFF, FLUX_ENTROPY_ATM, FLUX_ENTROPY_ATM_SPEC, FLUX_ENTROPY_ATM_LAT, &
 FLUX_ENTROPY_ATM_RAD, FLUX_ENTROPY_ATM_SOIL, FLUX_ENTROPY_ATM_SNOW, FLUX_ENTROPY_ATM_FOREST, &
 FLUX_ENTROPY_SNOW_SOIL, FLUX_ENTROPY_FOREST_SOIL, &
 FLUX_WV, FLUX_WV_SOIL, FLUX_WV_SNOW, FLUX_WV_VEGLOW_DRY, FLUX_WV_VEGLOW_WET, &
 FLUX_WV_VEGLOW_WET, FLUX_WV_FOREST_DRY, FLUX_WV_FOREST_WET, &
 FLUX_PREC_LIQ, FLUX_PREC_ICE, FLUX_W_LIQ_ATM_SOIL, FLUX_W_LIQ_ATM_SNOW, FLUX_W_ICE_ATM_SNOW, &
 FLUX_W_LIQ_ATM_VEGLOW, FLUX_W_LIQ_ATM_FOREST, FLUX_WATER_SNOW_SOIL, &
 FLUX_SOIL_ENTROPY, FLUX_SOIL_WATER, FLUX_SOIL_WATER_VEG, FLUX_SOIL_WATER_BOTTOM_DIAG, &
 LAMBDAG, FLUX_SNOW_ENTROPY, FLUX_SNOW_WATER, &
 QSMAX, QSMIN, CG, RHOG, PSIS_SAT, KAPPA_SAT, PAR_B, PAR_C, PAR_FB, QS_REL_VEGWILT, QS_REL_VEGREF, &
 LAI_MAX, LAI_VEGLOW, LAI_FOREST, FVEGLOW, FFOREST, FVEGLOWWET, FFORESTWET, ROOT_DEPTH_VEGLOW, &
 QW_VEGLOW, QW_FOREST, ROOT_DEPTH_FOREST, QW_VEGLOW_MAX, QW_FOREST_MAX, DZSUM_VEGLOW, DZSUM_FOREST, &
 FSNCOVER, DZ_DZSUM_VEGLOW, DZ_DZSUM_FOREST, MASK_SOIL, MASK_SOIL_AXED,&
 D_LEV_SNOW_BASE, RHOSNOW_FRESH, RHOSNOW_FIRN, RHOSNOW_OLD, VAL_MISSING, &
 PAR_FA, T0, T_FICE_CRIT

IMPLICIT NONE

INTEGER :: IUNIT=17, NSTEP0, NPROCX, NPROCY, PROC_IND, IERR_OPEN, NX_GLOB, NY_GLOB, K, J, I
REAL :: DTSTEP_EXT, ZDTSTEP, ZTIME
CHARACTER(LEN=30) :: FILE_IN="moloch_pochva.rf", FILE_OUT

REAL, DIMENSION(NX_GLOB, NY_GLOB), INTENT(IN) ::  FIELD_GLOB_R
INTEGER, DIMENSION(NX_GLOB, NY_GLOB) ::  FIELD_GLOB_I

! For testing output "test" --->
INTEGER :: N_OUTPOINT=0
INTEGER, DIMENSION(10) :: &
 I_OUTPOINT=    (/ -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,  -1/), &
 J_OUTPOINT=    (/ -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,  -1/), &
 IPROC_OUTPOINT=(/ -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,  -1/)
! <---

#ifdef mpi
  include 'mpif.h'
  integer :: comm, error, status(mpi_status_size)

  comm = mpi_comm_world
#endif
 
 IF (PROC_IND == 0) THEN

   OPEN (IUNIT, FILE=FILE_IN, FORM='UNFORMATTED', STATUS='OLD', IOSTAT=IERR_OPEN)

   IF (IERR_OPEN /= 0) THEN
     CLOSE (IUNIT)
     PRINT *,'Restart, error in ',TRIM(FILE_IN),' opening, stop'
#ifdef mpi
        call mpi_abort(comm, error)
#endif
     STOP
   ENDIF

   READ (IUNIT) NSTEP0
   READ (IUNIT) ZDTSTEP
   READ (IUNIT) NPOINT_X, NPOINT_Y, NLEV_SOIL_MAX

   DTSTEP=DTSTEP_EXT

 ENDIF

#ifdef mpi
   CALL MPI_BCAST(NSTEP0,        1, MPI_INTEGER, 0, COMM, ERROR)   
   CALL MPI_BCAST(DTSTEP,        1, MPI_REAL,    0, COMM, ERROR)   
   CALL MPI_BCAST(NPOINT_X,      1, MPI_INTEGER, 0, COMM, ERROR)   
   CALL MPI_BCAST(NPOINT_Y,      1, MPI_INTEGER, 0, COMM, ERROR)   
   CALL MPI_BCAST(NLEV_SOIL_MAX, 1, MPI_INTEGER, 0, COMM, ERROR)   
#endif

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

 IF (PROC_IND == 0) THEN

   READ (IUNIT) LEV_SOIL(0:NLEV_SOIL_MAX)
   READ (IUNIT) LEV_SOIL_H(1:NLEV_SOIL_MAX)
   READ (IUNIT) D_LEV_SOIL(1:NLEV_SOIL_MAX)
   READ (IUNIT) D_LEV_SOIL_H(0:NLEV_SOIL_MAX-1)
   READ (IUNIT) SUM_D_LEV_SOIL_H
   READ (IUNIT) LAI_MAX

 ENDIF

#ifdef mpi
   CALL MPI_BCAST(LEV_SOIL(0:NLEV_SOIL_MAX), NLEV_SOIL_MAX+1, MPI_REAL,    0, COMM, ERROR)   
   CALL MPI_BCAST(LEV_SOIL_H(1:NLEV_SOIL_MAX), NLEV_SOIL_MAX, MPI_REAL,    0, COMM, ERROR)   
   CALL MPI_BCAST(D_LEV_SOIL(1:NLEV_SOIL_MAX), NLEV_SOIL_MAX, MPI_REAL,    0, COMM, ERROR)   
   CALL MPI_BCAST(D_LEV_SOIL_H(0:NLEV_SOIL_MAX-1), NLEV_SOIL_MAX, MPI_REAL,    0, COMM, ERROR)   

   CALL MPI_BCAST(SUM_D_LEV_SOIL_H, 1, MPI_REAL,    0, COMM, ERROR)   
   CALL MPI_BCAST(LAI_MAX, 1, MPI_REAL,    0, COMM, ERROR)   
#endif

! Allocations --->

 IF (ALLOCATED(NLEV_SOIL_HEAT)) DEALLOCATE(NLEV_SOIL_HEAT)
 ALLOCATE(NLEV_SOIL_HEAT(NPOINT_X,NPOINT_Y))
 NLEV_SOIL_HEAT(:,:)=0
 IF (ALLOCATED(NLEV_SOIL_WATER)) DEALLOCATE(NLEV_SOIL_WATER)
 ALLOCATE(NLEV_SOIL_WATER(NPOINT_X,NPOINT_Y))
 NLEV_SOIL_WATER(:,:)=0

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
 TSNOW(:,:,:)=VAL_MISSING
 IF (ALLOCATED(LEV_SNOW)) DEALLOCATE(LEV_SNOW)
 ALLOCATE(LEV_SNOW(NPOINT_X,NPOINT_Y,0:NLEV_SNOW))
 LEV_SNOW(:,:,:)=VAL_MISSING
 IF (ALLOCATED(FICE_SNOW)) DEALLOCATE(FICE_SNOW)
 ALLOCATE(FICE_SNOW(NPOINT_X,NPOINT_Y,0:NLEV_SNOW))
 FICE_SNOW(:,:,:)=VAL_MISSING
 IF (ALLOCATED(RHOSNOW)) DEALLOCATE(RHOSNOW)
 ALLOCATE(RHOSNOW(NPOINT_X,NPOINT_Y,0:NLEV_SNOW))
 RHOSNOW(:,:,:)=VAL_MISSING
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

! <---

 IF (PROC_IND == 0) CALL RREC2_INT (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_I)
 CALL DISPERSE_INT (FIELD_GLOB_I, NLEV_SOIL_HEAT)

 IF (PROC_IND == 0) CALL RREC2_INT (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_I)
 CALL DISPERSE_INT (FIELD_GLOB_I, NLEV_SOIL_WATER)

 DO K=0,NLEV_SOIL_MAX
   IF (PROC_IND == 0) CALL RREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R)
   CALL DISPERSE (FIELD_GLOB_R, QSMAX(:,:,K))
 ENDDO

 DO K=0,NLEV_SOIL_MAX
   IF (PROC_IND == 0) CALL RREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R)
   CALL DISPERSE (FIELD_GLOB_R, QSMIN(:,:,K))
 ENDDO

 DO K=0,NLEV_SOIL_MAX
   IF (PROC_IND == 0) CALL RREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R)
   CALL DISPERSE (FIELD_GLOB_R, CG(:,:,K))
 ENDDO

 DO K=0,NLEV_SOIL_MAX
   IF (PROC_IND == 0) CALL RREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R)
   CALL DISPERSE (FIELD_GLOB_R, RHOG(:,:,K))
 ENDDO

 DO K=0,NLEV_SOIL_MAX
   IF (PROC_IND == 0) CALL RREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R)
   CALL DISPERSE (FIELD_GLOB_R, PSIS_SAT(:,:,K))
 ENDDO

 DO K=0,NLEV_SOIL_MAX
   IF (PROC_IND == 0) CALL RREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R)
   CALL DISPERSE (FIELD_GLOB_R, KAPPA_SAT(:,:,K))
 ENDDO

 DO K=0,NLEV_SOIL_MAX
   IF (PROC_IND == 0) CALL RREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R)
   CALL DISPERSE (FIELD_GLOB_R, PAR_B(:,:,K))
 ENDDO

 DO K=0,NLEV_SOIL_MAX
   IF (PROC_IND == 0) CALL RREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R)
   CALL DISPERSE (FIELD_GLOB_R, PAR_C(:,:,K))
 ENDDO

 DO K=0,NLEV_SOIL_MAX
   IF (PROC_IND == 0) CALL RREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R)
   CALL DISPERSE (FIELD_GLOB_R, PAR_FB(:,:,K))
 ENDDO

 DO K=0,NLEV_SOIL_MAX
   IF (PROC_IND == 0) CALL RREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R)
   CALL DISPERSE (FIELD_GLOB_R, QS_REL_VEGWILT(:,:,K))
 ENDDO

 DO K=0,NLEV_SOIL_MAX
   IF (PROC_IND == 0) CALL RREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R)
   CALL DISPERSE (FIELD_GLOB_R, QS_REL_VEGREF(:,:,K))
 ENDDO

 IF (PROC_IND == 0) CALL RREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R)
 CALL DISPERSE (FIELD_GLOB_R, LAI_VEGLOW)

 IF (PROC_IND == 0) CALL RREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R)
 CALL DISPERSE (FIELD_GLOB_R, LAI_FOREST)

 IF (PROC_IND == 0) CALL RREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R)
 CALL DISPERSE (FIELD_GLOB_R, FVEGLOW)

 IF (PROC_IND == 0) CALL RREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R)
 CALL DISPERSE (FIELD_GLOB_R, FFOREST)

 IF (PROC_IND == 0) CALL RREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R)
 CALL DISPERSE (FIELD_GLOB_R, FVEGLOWWET)

 IF (PROC_IND == 0) CALL RREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R)
 CALL DISPERSE (FIELD_GLOB_R, FFORESTWET)

 IF (PROC_IND == 0) CALL RREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R)
 CALL DISPERSE (FIELD_GLOB_R, ROOT_DEPTH_VEGLOW)

 IF (PROC_IND == 0) CALL RREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R)
 CALL DISPERSE (FIELD_GLOB_R, ROOT_DEPTH_FOREST)

 IF (PROC_IND == 0) CALL RREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R)
 CALL DISPERSE (FIELD_GLOB_R, QW_VEGLOW)

 IF (PROC_IND == 0) CALL RREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R)
 CALL DISPERSE (FIELD_GLOB_R, QW_FOREST)

 IF (PROC_IND == 0) CALL RREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R)
 CALL DISPERSE (FIELD_GLOB_R, QW_VEGLOW_MAX)

 IF (PROC_IND == 0) CALL RREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R)
 CALL DISPERSE (FIELD_GLOB_R, QW_FOREST_MAX)

 IF (PROC_IND == 0) CALL RREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R)
 CALL DISPERSE (FIELD_GLOB_R, DZSUM_VEGLOW)

 IF (PROC_IND == 0) CALL RREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R)
 CALL DISPERSE (FIELD_GLOB_R, DZSUM_FOREST)

 DO K=0,NLEV_SOIL_MAX-1
   IF (PROC_IND == 0) CALL RREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R)
   CALL DISPERSE (FIELD_GLOB_R, DZ_DZSUM_VEGLOW(:,:,K))
 ENDDO

 IF (PROC_IND == 0) CALL RREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R)
 CALL DISPERSE (FIELD_GLOB_R, DZSUM_FOREST)

 DO K=0,NLEV_SOIL_MAX-1
   IF (PROC_IND == 0) CALL RREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R)
   CALL DISPERSE (FIELD_GLOB_R, DZ_DZSUM_FOREST(:,:,K))
 ENDDO

 IF (PROC_IND == 0) CALL RREC2_INT (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_I)
 CALL DISPERSE_INT (FIELD_GLOB_I, MASK_SOIL)

 IF (PROC_IND == 0) CALL RREC2_INT (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_I)
 CALL DISPERSE_INT (FIELD_GLOB_I, MASK_SOIL_AXED)

 IF (PROC_IND == 0) CALL RREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R)
 CALL DISPERSE (FIELD_GLOB_R, TSURF)

 DO K=0,NLEV_SOIL_MAX
   IF (PROC_IND == 0) CALL RREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R)
   CALL DISPERSE (FIELD_GLOB_R, TS(:,:,K))
 ENDDO

 DO K=0,NLEV_SOIL_MAX
   IF (PROC_IND == 0) CALL RREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R)
   CALL DISPERSE (FIELD_GLOB_R, QS(:,:,K))
 ENDDO

 DO K=0,NLEV_SOIL_MAX
   IF (PROC_IND == 0) CALL RREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R)
   CALL DISPERSE (FIELD_GLOB_R, QS_REL(:,:,K))
 ENDDO

 DO K=0,NLEV_SOIL_MAX
   IF (PROC_IND == 0) CALL RREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R)
   CALL DISPERSE (FIELD_GLOB_R, QSI(:,:,K))
 ENDDO

 DO K=0,NLEV_SOIL_MAX
   IF (PROC_IND == 0) CALL RREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R)
   CALL DISPERSE (FIELD_GLOB_R, FRAC_SICE(:,:,K))
 ENDDO

 DO K=0,NLEV_SOIL_MAX
   IF (PROC_IND == 0) CALL RREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R)
   CALL DISPERSE (FIELD_GLOB_R, DFRAC_SICE_DT(:,:,K))
 ENDDO

 IF (PROC_IND == 0) CALL RREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R)
 CALL DISPERSE (FIELD_GLOB_R, FSNCOVER)

 DO K=0,NLEV_SNOW
   IF (PROC_IND == 0) CALL RREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R)
   CALL DISPERSE (FIELD_GLOB_R, LEV_SNOW(:,:,K))
 ENDDO

 DO K=0,NLEV_SNOW
   IF (PROC_IND == 0) CALL RREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R)
   CALL DISPERSE (FIELD_GLOB_R, TSNOW(:,:,K))
 ENDDO

 DO K=0,NLEV_SNOW
   IF (PROC_IND == 0) CALL RREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R)
   CALL DISPERSE (FIELD_GLOB_R, FICE_SNOW(:,:,K))
 ENDDO

 DO K=0,NLEV_SNOW
   IF (PROC_IND == 0) CALL RREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R)
   CALL DISPERSE (FIELD_GLOB_R, RHOSNOW(:,:,K))
 ENDDO

 DO K=0,NLEV_SNOW
   DO J=1,NPOINT_Y
     DO I=1,NPOINT_X
       IF (ABS(LEV_SNOW(I,J,K)-VAL_MISSING) < 10.) LEV_SNOW(I,J,K)=VAL_MISSING
       IF (ABS(TSNOW(I,J,K)-VAL_MISSING) < 10.) TSNOW(I,J,K)=VAL_MISSING
       IF (ABS(FICE_SNOW(I,J,K)-VAL_MISSING) < 10.) FICE_SNOW(I,J,K)=VAL_MISSING
       IF (ABS(RHOSNOW(I,J,K)-VAL_MISSING) < 10.) THEN
         RHOSNOW(I,J,K)=VAL_MISSING
       ELSE
         RHOSNOW(I,J,K)=MIN(MAX(RHOSNOW(I,J,K), RHOSNOW_FRESH), RHOSNOW_FIRN)
       ENDIF
     ENDDO
   ENDDO
 ENDDO

 DO K=0,NLEV_SNOW
   IF (PROC_IND == 0) CALL RREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R)
   CALL DISPERSE (FIELD_GLOB_R, SNOW_AGE(:,:,K))
   SNOW_AGE(:,:,K)=MAX(SNOW_AGE(:,:,K), 0.)
 ENDDO

 DO K=0,NLEV_SNOW
   IF (PROC_IND == 0) CALL RREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R)
   CALL DISPERSE (FIELD_GLOB_R, SNOW_MELT_AGE(:,:,K))
   SNOW_MELT_AGE(:,:,K)=MAX(SNOW_MELT_AGE(:,:,K), 0.)
 ENDDO

 IF (PROC_IND == 0) CALL RREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R)
 CALL DISPERSE (FIELD_GLOB_R, D_LEV_SNOW)

 IF (PROC_IND == 0) CALL RREC2_INT (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_I)
 CALL DISPERSE_INT (FIELD_GLOB_I, FLAG_SNOW_THICK)

 IF (PROC_IND == 0) CALL RREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R)
 CALL DISPERSE (FIELD_GLOB_R, SNOW_DIRT)

 IF (PROC_IND == 0) CALL RREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R)
 CALL DISPERSE (FIELD_GLOB_R, DLEV_SNOW_BOTTOM)

 IF (PROC_IND == 0) CALL RREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R)
 CALL DISPERSE (FIELD_GLOB_R, TFOREST)

 IF (PROC_IND == 0) THEN
   CLOSE (IUNIT)
   PRINT *
   PRINT *,TRIM(FILE_IN),' read, nstep0 ',NSTEP0
 ENDIF

 PAR_FA=-4./ALOG((T0+T_FICE_CRIT)/T0)

 MASK_SOIL(1,         1:NPOINT_Y) = 0
 MASK_SOIL(NPOINT_X,  1:NPOINT_Y) = 0
 MASK_SOIL(1:NPOINT_X,   1)       = 0
 MASK_SOIL(1:NPOINT_X,NPOINT_Y)   = 0

! Testing output "test" --->
 DO K=1,N_OUTPOINT
   IF (PROC_IND == IPROC_OUTPOINT(K)) THEN
     I=I_OUTPOINT(K)
     J=J_OUTPOINT(K)
     WRITE (FILE_OUT,'(A,I2.2,A)') "restart_pochva_",K,".dat2"
     OPEN (IUNIT, FILE=FILE_OUT, STATUS="UNKNOWN")
       WRITE (IUNIT,*) DTSTEP, NPOINT_X, NPOINT_Y, NLEV_SOIL_MAX
       WRITE (IUNIT,*) LEV_SOIL(0:NLEV_SOIL_MAX), LEV_SOIL_H(1:NLEV_SOIL_MAX), D_LEV_SOIL(1:NLEV_SOIL_MAX), &
  D_LEV_SOIL_H(0:NLEV_SOIL_MAX-1), SUM_D_LEV_SOIL_H, LAI_MAX
       WRITE (IUNIT,*) NLEV_SOIL_HEAT(I,J)
       WRITE (IUNIT,*) NLEV_SOIL_WATER(I,J)
       WRITE (IUNIT,*) QSMAX(I,J,0:NLEV_SOIL_MAX)
       WRITE (IUNIT,*) QSMIN(I,J,0:NLEV_SOIL_MAX)
       WRITE (IUNIT,*) CG(I,J,0:NLEV_SOIL_MAX)
       WRITE (IUNIT,*) RHOG(I,J,0:NLEV_SOIL_MAX)
       WRITE (IUNIT,*) PSIS_SAT(I,J,0:NLEV_SOIL_MAX)
       WRITE (IUNIT,*) KAPPA_SAT(I,J,0:NLEV_SOIL_MAX)
       WRITE (IUNIT,*) PAR_B(I,J,0:NLEV_SOIL_MAX)
       WRITE (IUNIT,*) PAR_C(I,J,0:NLEV_SOIL_MAX)
       WRITE (IUNIT,*) PAR_FB(I,J,0:NLEV_SOIL_MAX)
       WRITE (IUNIT,*) QS_REL_VEGWILT(I,J,0:NLEV_SOIL_MAX)
       WRITE (IUNIT,*) QS_REL_VEGREF(I,J,0:NLEV_SOIL_MAX)
       WRITE (IUNIT,*) LAI_VEGLOW(I,J)
       WRITE (IUNIT,*) LAI_FOREST(I,J)
       WRITE (IUNIT,*) FVEGLOW(I,J)
       WRITE (IUNIT,*) FFOREST(I,J)
       WRITE (IUNIT,*) FVEGLOWWET(I,J)
       WRITE (IUNIT,*) FFORESTWET(I,J)
       WRITE (IUNIT,*) ROOT_DEPTH_VEGLOW(I,J)
       WRITE (IUNIT,*) ROOT_DEPTH_FOREST(I,J)
       WRITE (IUNIT,*) QW_VEGLOW(I,J)
       WRITE (IUNIT,*) QW_FOREST(I,J)
       WRITE (IUNIT,*) QW_VEGLOW_MAX(I,J)
       WRITE (IUNIT,*) QW_FOREST_MAX(I,J)
       WRITE (IUNIT,*) DZSUM_VEGLOW(I,J)
       WRITE (IUNIT,*) DZSUM_FOREST(I,J)
       WRITE (IUNIT,*) DZ_DZSUM_VEGLOW(I,J,0:NLEV_SOIL_MAX-1)
       WRITE (IUNIT,*) DZSUM_FOREST(I,J)
       WRITE (IUNIT,*) DZ_DZSUM_FOREST(I,J,0:NLEV_SOIL_MAX-1)
       WRITE (IUNIT,*) MASK_SOIL(I,J)
       WRITE (IUNIT,*) MASK_SOIL_AXED(I,J)
       WRITE (IUNIT,*) TSURF(I,J)
       WRITE (IUNIT,*) TS(I,J,0:NLEV_SOIL_MAX)
       WRITE (IUNIT,*) QS(I,J,0:NLEV_SOIL_MAX)
       WRITE (IUNIT,*) QS_REL(I,J,0:NLEV_SOIL_MAX)
       WRITE (IUNIT,*) QSI(I,J,0:NLEV_SOIL_MAX)
       WRITE (IUNIT,*) FRAC_SICE(I,J,0:NLEV_SOIL_MAX)
       WRITE (IUNIT,*) DFRAC_SICE_DT(I,J,0:NLEV_SOIL_MAX)
       WRITE (IUNIT,*) FSNCOVER(I,J)
       WRITE (IUNIT,*) LEV_SNOW(I,J,0:NLEV_SNOW)
       WRITE (IUNIT,*) TSNOW(I,J,0:NLEV_SNOW)
       WRITE (IUNIT,*) FICE_SNOW(I,J,0:NLEV_SNOW)
       WRITE (IUNIT,*) RHOSNOW(I,J,0:NLEV_SNOW)
       WRITE (IUNIT,*) SNOW_AGE(I,J,0:NLEV_SNOW)
       WRITE (IUNIT,*) SNOW_MELT_AGE(I,J,0:NLEV_SNOW)
       WRITE (IUNIT,*) D_LEV_SNOW(I,J)
       WRITE (IUNIT,*) FLAG_SNOW_THICK(I,J)
       WRITE (IUNIT,*) SNOW_DIRT(I,J)
       WRITE (IUNIT,*) DLEV_SNOW_BOTTOM(I,J)
       WRITE (IUNIT,*) TFOREST(I,J)
     CLOSE (IUNIT) 
   ENDIF
 ENDDO
! <---

RETURN
END SUBROUTINE RDRF_POCHVA
!###############################################################################################################
SUBROUTINE WRRF_POCHVA(KSTEP, PROC_IND, FIELD_GLOB_R, NX_GLOB, NY_GLOB, IRF)
! Write reastart file with pochva variables

USE MODULE_POCHVA, ONLY : DTSTEP, NPOINT_X, NPOINT_Y, NLEV_SOIL_MAX, NLEV_SNOW, &
 NLEV_SOIL_HEAT, NLEV_SOIL_WATER, LEV_SOIL, LEV_SOIL_H, D_LEV_SOIL, D_LEV_SOIL_H, SUM_D_LEV_SOIL_H, &
 TS, QS, QS_REL, PSIS, KAPPA_H, FRAC_SICE, DFRAC_SICE_DT, QSURF, TSURF, TFOREST, QFOREST, QS_REL_VEG, &
 D_LEV_SNOW, TSNOW, LEV_SNOW, FICE_SNOW, RHOSNOW, SNOW_AGE, SNOW_MELT_AGE, DLEV_SNOW_BOTTOM, FLAG_SNOW_THICK, SNOW_DIRT, QSI, & 
 QSMAX, QSMIN, CG, RHOG, PSIS_SAT, KAPPA_SAT, PAR_B, PAR_C, PAR_FB, QS_REL_VEGWILT, QS_REL_VEGREF, &
 LAI_MAX, LAI_VEGLOW, LAI_FOREST, FVEGLOW, FFOREST, FVEGLOWWET, FFORESTWET, ROOT_DEPTH_VEGLOW, &
 QW_VEGLOW, QW_FOREST, ROOT_DEPTH_FOREST, QW_VEGLOW_MAX, QW_FOREST_MAX, DZSUM_VEGLOW, DZSUM_FOREST, &
 FSNCOVER, DZ_DZSUM_VEGLOW, DZ_DZSUM_FOREST, MASK_SOIL, MASK_SOIL_AXED

IMPLICIT NONE

INTEGER :: IUNIT=27, KSTEP, PROC_IND, NX_GLOB, NY_GLOB, K, IRF, I1, I2, I3, &
 ILON1, ILON2, JLAT1, JLAT2, FLAG_LON
CHARACTER(LEN=30) :: FILE_OUT

REAL, DIMENSION(NX_GLOB, NY_GLOB), INTENT(IN) ::  FIELD_GLOB_R
INTEGER, DIMENSION(NX_GLOB, NY_GLOB) ::  FIELD_GLOB_I

 ILON1=1
 ILON2=NX_GLOB
 JLAT1=1
 JLAT2=NY_GLOB
 FLAG_LON=0
 
 IF (PROC_IND == 0) THEN

   I1=IRF/10
   I2=I1*10
   I3=IRF-I2
   WRITE (FILE_OUT,'(A,I2.2,A)') 'moloch_pochva_',I3,'.rf'
   OPEN (IUNIT, FILE=TRIM(FILE_OUT), FORM='UNFORMATTED', STATUS='UNKNOWN')

   WRITE (IUNIT) KSTEP+1
   WRITE (IUNIT) DTSTEP
   WRITE (IUNIT) NPOINT_X, NPOINT_Y, NLEV_SOIL_MAX
   WRITE (IUNIT) LEV_SOIL(0:NLEV_SOIL_MAX)
   WRITE (IUNIT) LEV_SOIL_H(1:NLEV_SOIL_MAX)
   WRITE (IUNIT) D_LEV_SOIL(1:NLEV_SOIL_MAX)
   WRITE (IUNIT) D_LEV_SOIL_H(0:NLEV_SOIL_MAX-1)
   WRITE (IUNIT) SUM_D_LEV_SOIL_H
   WRITE (IUNIT) LAI_MAX

 ENDIF

 CALL COLLECT_INT (NLEV_SOIL_HEAT, FIELD_GLOB_I)
 IF (PROC_IND == 0) CALL WREC2_INT (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_I, ILON1, ILON2, JLAT1, JLAT2, FLAG_LON)

 CALL COLLECT_INT (NLEV_SOIL_WATER, FIELD_GLOB_I)
 IF (PROC_IND == 0) CALL WREC2_INT (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_I, ILON1, ILON2, JLAT1, JLAT2, FLAG_LON)

 DO K=0,NLEV_SOIL_MAX
   CALL COLLECT (QSMAX(:,:,K), FIELD_GLOB_R)
   IF (PROC_IND == 0) CALL WREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R, ILON1, ILON2, JLAT1, JLAT2, FLAG_LON)
 ENDDO

 DO K=0,NLEV_SOIL_MAX
   CALL COLLECT (QSMIN(:,:,K), FIELD_GLOB_R)
   IF (PROC_IND == 0) CALL WREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R, ILON1, ILON2, JLAT1, JLAT2, FLAG_LON)
 ENDDO

 DO K=0,NLEV_SOIL_MAX
   CALL COLLECT (CG(:,:,K), FIELD_GLOB_R)
   IF (PROC_IND == 0) CALL WREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R, ILON1, ILON2, JLAT1, JLAT2, FLAG_LON)
 ENDDO

 DO K=0,NLEV_SOIL_MAX
   CALL COLLECT (RHOG(:,:,K), FIELD_GLOB_R)
   IF (PROC_IND == 0) CALL WREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R, ILON1, ILON2, JLAT1, JLAT2, FLAG_LON)
 ENDDO

 DO K=0,NLEV_SOIL_MAX
   CALL COLLECT (PSIS_SAT(:,:,K), FIELD_GLOB_R)
   IF (PROC_IND == 0) CALL WREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R, ILON1, ILON2, JLAT1, JLAT2, FLAG_LON)
 ENDDO

 DO K=0,NLEV_SOIL_MAX
   CALL COLLECT (KAPPA_SAT(:,:,K), FIELD_GLOB_R)
   IF (PROC_IND == 0) CALL WREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R, ILON1, ILON2, JLAT1, JLAT2, FLAG_LON)
 ENDDO

 DO K=0,NLEV_SOIL_MAX
   CALL COLLECT (PAR_B(:,:,K), FIELD_GLOB_R)
   IF (PROC_IND == 0) CALL WREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R, ILON1, ILON2, JLAT1, JLAT2, FLAG_LON)
 ENDDO

 DO K=0,NLEV_SOIL_MAX
   CALL COLLECT (PAR_C(:,:,K), FIELD_GLOB_R)
   IF (PROC_IND == 0) CALL WREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R, ILON1, ILON2, JLAT1, JLAT2, FLAG_LON)
 ENDDO

 DO K=0,NLEV_SOIL_MAX
   CALL COLLECT (PAR_FB(:,:,K), FIELD_GLOB_R)
   IF (PROC_IND == 0) CALL WREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R, ILON1, ILON2, JLAT1, JLAT2, FLAG_LON)
 ENDDO

 DO K=0,NLEV_SOIL_MAX
   CALL COLLECT (QS_REL_VEGWILT(:,:,K), FIELD_GLOB_R)
   IF (PROC_IND == 0) CALL WREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R, ILON1, ILON2, JLAT1, JLAT2, FLAG_LON)
 ENDDO

 DO K=0,NLEV_SOIL_MAX
   CALL COLLECT (QS_REL_VEGREF(:,:,K), FIELD_GLOB_R)
   IF (PROC_IND == 0) CALL WREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R, ILON1, ILON2, JLAT1, JLAT2, FLAG_LON)
 ENDDO

 CALL COLLECT (LAI_VEGLOW, FIELD_GLOB_R)
 IF (PROC_IND == 0) CALL WREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R, ILON1, ILON2, JLAT1, JLAT2, FLAG_LON)

 CALL COLLECT (LAI_FOREST, FIELD_GLOB_R)
 IF (PROC_IND == 0) CALL WREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R, ILON1, ILON2, JLAT1, JLAT2, FLAG_LON)

 CALL COLLECT (FVEGLOW, FIELD_GLOB_R)
 IF (PROC_IND == 0) CALL WREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R, ILON1, ILON2, JLAT1, JLAT2, FLAG_LON)

 CALL COLLECT (FFOREST, FIELD_GLOB_R)
 IF (PROC_IND == 0) CALL WREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R, ILON1, ILON2, JLAT1, JLAT2, FLAG_LON)

 CALL COLLECT (FVEGLOWWET, FIELD_GLOB_R)
 IF (PROC_IND == 0) CALL WREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R, ILON1, ILON2, JLAT1, JLAT2, FLAG_LON)

 CALL COLLECT (FFORESTWET, FIELD_GLOB_R)
 IF (PROC_IND == 0) CALL WREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R, ILON1, ILON2, JLAT1, JLAT2, FLAG_LON)

 CALL COLLECT (ROOT_DEPTH_VEGLOW, FIELD_GLOB_R)
 IF (PROC_IND == 0) CALL WREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R, ILON1, ILON2, JLAT1, JLAT2, FLAG_LON)

 CALL COLLECT (ROOT_DEPTH_FOREST, FIELD_GLOB_R)
 IF (PROC_IND == 0) CALL WREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R, ILON1, ILON2, JLAT1, JLAT2, FLAG_LON)

 CALL COLLECT (QW_VEGLOW, FIELD_GLOB_R)
 IF (PROC_IND == 0) CALL WREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R, ILON1, ILON2, JLAT1, JLAT2, FLAG_LON)

 CALL COLLECT (QW_FOREST, FIELD_GLOB_R)
 IF (PROC_IND == 0) CALL WREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R, ILON1, ILON2, JLAT1, JLAT2, FLAG_LON)

 CALL COLLECT (QW_VEGLOW_MAX, FIELD_GLOB_R)
 IF (PROC_IND == 0) CALL WREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R, ILON1, ILON2, JLAT1, JLAT2, FLAG_LON)

 CALL COLLECT (QW_FOREST_MAX, FIELD_GLOB_R)
 IF (PROC_IND == 0) CALL WREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R, ILON1, ILON2, JLAT1, JLAT2, FLAG_LON)

 CALL COLLECT (DZSUM_VEGLOW, FIELD_GLOB_R)
 IF (PROC_IND == 0) CALL WREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R, ILON1, ILON2, JLAT1, JLAT2, FLAG_LON)

 CALL COLLECT (DZSUM_FOREST, FIELD_GLOB_R)
 IF (PROC_IND == 0) CALL WREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R, ILON1, ILON2, JLAT1, JLAT2, FLAG_LON)

 DO K=0,NLEV_SOIL_MAX-1
   CALL COLLECT (DZ_DZSUM_VEGLOW(:,:,K), FIELD_GLOB_R)
   IF (PROC_IND == 0) CALL WREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R, ILON1, ILON2, JLAT1, JLAT2, FLAG_LON)
 ENDDO

 CALL COLLECT (DZSUM_FOREST, FIELD_GLOB_R)
 IF (PROC_IND == 0) CALL WREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R, ILON1, ILON2, JLAT1, JLAT2, FLAG_LON)

 DO K=0,NLEV_SOIL_MAX-1
   CALL COLLECT (DZ_DZSUM_FOREST(:,:,K), FIELD_GLOB_R)
   IF (PROC_IND == 0) CALL WREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R, ILON1, ILON2, JLAT1, JLAT2, FLAG_LON)
 ENDDO

 CALL COLLECT_INT (MASK_SOIL, FIELD_GLOB_I)
 IF (PROC_IND == 0) CALL WREC2_INT (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_I, ILON1, ILON2, JLAT1, JLAT2, FLAG_LON)

 CALL COLLECT_INT (MASK_SOIL_AXED, FIELD_GLOB_I)
 IF (PROC_IND == 0) CALL WREC2_INT (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_I, ILON1, ILON2, JLAT1, JLAT2, FLAG_LON)

 CALL COLLECT (TSURF, FIELD_GLOB_R)
 IF (PROC_IND == 0) CALL WREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R, ILON1, ILON2, JLAT1, JLAT2, FLAG_LON)

 DO K=0,NLEV_SOIL_MAX
   CALL COLLECT (TS(:,:,K), FIELD_GLOB_R)
   IF (PROC_IND == 0) CALL WREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R, ILON1, ILON2, JLAT1, JLAT2, FLAG_LON)
 ENDDO

 DO K=0,NLEV_SOIL_MAX
   CALL COLLECT (QS(:,:,K), FIELD_GLOB_R)
   IF (PROC_IND == 0) CALL WREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R, ILON1, ILON2, JLAT1, JLAT2, FLAG_LON)
 ENDDO

 DO K=0,NLEV_SOIL_MAX
   CALL COLLECT (QS_REL(:,:,K), FIELD_GLOB_R)
   IF (PROC_IND == 0) CALL WREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R, ILON1, ILON2, JLAT1, JLAT2, FLAG_LON)
 ENDDO

 DO K=0,NLEV_SOIL_MAX
   CALL COLLECT (QSI(:,:,K), FIELD_GLOB_R)
   IF (PROC_IND == 0) CALL WREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R, ILON1, ILON2, JLAT1, JLAT2, FLAG_LON)
 ENDDO

 DO K=0,NLEV_SOIL_MAX
   CALL COLLECT (FRAC_SICE(:,:,K), FIELD_GLOB_R)
   IF (PROC_IND == 0) CALL WREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R, ILON1, ILON2, JLAT1, JLAT2, FLAG_LON)
 ENDDO

 DO K=0,NLEV_SOIL_MAX
   CALL COLLECT (DFRAC_SICE_DT(:,:,K), FIELD_GLOB_R)
   IF (PROC_IND == 0) CALL WREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R, ILON1, ILON2, JLAT1, JLAT2, FLAG_LON)
 ENDDO

 CALL COLLECT (FSNCOVER, FIELD_GLOB_R)
 IF (PROC_IND == 0) CALL WREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R, ILON1, ILON2, JLAT1, JLAT2, FLAG_LON)

 DO K=0,NLEV_SNOW
   CALL COLLECT (LEV_SNOW(:,:,K), FIELD_GLOB_R)
   IF (PROC_IND == 0) CALL WREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R, ILON1, ILON2, JLAT1, JLAT2, FLAG_LON)
 ENDDO

 DO K=0,NLEV_SNOW
   CALL COLLECT (TSNOW(:,:,K), FIELD_GLOB_R)
   IF (PROC_IND == 0) CALL WREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R, ILON1, ILON2, JLAT1, JLAT2, FLAG_LON)
 ENDDO

 DO K=0,NLEV_SNOW
   CALL COLLECT (FICE_SNOW(:,:,K), FIELD_GLOB_R)
   IF (PROC_IND == 0) CALL WREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R, ILON1, ILON2, JLAT1, JLAT2, FLAG_LON)
 ENDDO

 DO K=0,NLEV_SNOW
   CALL COLLECT (RHOSNOW(:,:,K), FIELD_GLOB_R)
   IF (PROC_IND == 0) CALL WREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R, ILON1, ILON2, JLAT1, JLAT2, FLAG_LON)
 ENDDO

 DO K=0,NLEV_SNOW
   CALL COLLECT (SNOW_AGE(:,:,K), FIELD_GLOB_R)
   IF (PROC_IND == 0) CALL WREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R, ILON1, ILON2, JLAT1, JLAT2, FLAG_LON)
 ENDDO

 DO K=0,NLEV_SNOW
   CALL COLLECT (SNOW_MELT_AGE(:,:,K), FIELD_GLOB_R)
   IF (PROC_IND == 0) CALL WREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R, ILON1, ILON2, JLAT1, JLAT2, FLAG_LON)
 ENDDO

 CALL COLLECT (D_LEV_SNOW, FIELD_GLOB_R)
 IF (PROC_IND == 0) CALL WREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R, ILON1, ILON2, JLAT1, JLAT2, FLAG_LON)

 CALL COLLECT_INT (FLAG_SNOW_THICK, FIELD_GLOB_I)
 IF (PROC_IND == 0) CALL WREC2_INT (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_I, ILON1, ILON2, JLAT1, JLAT2, FLAG_LON)

 CALL COLLECT (SNOW_DIRT, FIELD_GLOB_R)
 IF (PROC_IND == 0) CALL WREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R, ILON1, ILON2, JLAT1, JLAT2, FLAG_LON)

 CALL COLLECT (DLEV_SNOW_BOTTOM, FIELD_GLOB_R)
 IF (PROC_IND == 0) CALL WREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R, ILON1, ILON2, JLAT1, JLAT2, FLAG_LON)

 CALL COLLECT (TFOREST, FIELD_GLOB_R)
 IF (PROC_IND == 0) CALL WREC2 (IUNIT, NX_GLOB, NY_GLOB, FIELD_GLOB_R, ILON1, ILON2, JLAT1, JLAT2, FLAG_LON)

 IF (PROC_IND == 0) THEN
   CLOSE (IUNIT)
   PRINT *
   PRINT *,TRIM(FILE_OUT),' written at step ',KSTEP
 ENDIF

RETURN
END SUBROUTINE WRRF_POCHVA
!###############################################################################################################
