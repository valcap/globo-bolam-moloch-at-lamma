! Last update 29/07/2024

! Version 24.1.1

! Sett. 2024: 
! Cambiata la scrittura di model_param_constant.bin - piu' grandezze.
! Nuovo formato mhf (scrittura record non 1d ma 2d)

! Ott. 2023: Option of compilation without ECMWF radiation scheme pachage,
! astronomical variables definition and aerosols and ozone content definitions
! are extrated from that scheme.

! May 2003: possibility to have full mhf and cut mht (in case of Globo 
! option) at the same time.
! Version for new storage system (as tintin cluster).

! August 2022: modification of output: it is possibile to output both mhf
! complete and mhf for cut domain in the same time or one of them,
! bolam.inp is modified also.

! July 2022: Modification in pochva scheme: definition of qsurf over bare soil
! to solution of problem of too humid surface air over very dry desert soil
! (top soil relative humidity <= 0.4).

! May 2022: Correction of vegetation fraction in subroutine surfradpar

! Dec. 2021: Input/output units have been reordered, for each type of data
! a specific unit number has been assigned.

! Dec. 2021: Correction of error in pochva scheme, the definition of qv_surf
! (air humidity) at the surface (mainly at bare soil surface) as the function
! of turbulent coefficient in the surface layer and top soil water content,
! error was due to turbulent coefficient that was not normalized to
! 1 m over surface.

! Nov. 2021: Substitution !$ symbols for mpi library implementation by ifdef mpi command

! Gen. 2021:
! parametrizzazione orogrdag revista, parametrizzazione blockdrag inserita,
! subroutine cloudfr revista,
! in subroutine radintec per old_rad bias (aggiunta a corirr) modificato, 
! in subroutine convection_kf timec modificato (timec = 4300.)
! Feb. 2020: shf di Globo in modo analogo a shf di Bolam
! Apr. 2019: Pochva per Globo, nuovo shf di Globo (analogo a shf di Bolam)
! Gen. 2019:
! Turbul. orogr. drag inserito in vdiff (coeff. ctofd)
! Inserita riduz. radiaz. all'inizio del run
!---------------------------------------------------------------
! Ago. 2018: 
! Nuovo schema del suolo (pochva), i nuovi dataset fisiografici,
! parametri del suolo sono diversi su vari livelli,
! manto nevoso muti-stato, il numero dei livelli di neve e' 11,
! la coordinate verticale nello manto nevoso e' kg/m/m,
! quando su un livello la neve non c'e', tutti i paramtri di neve sono -9999.;
! nuovo formato mhf: mfs_atm, mhf_soil e il file con i campi statici (model_param_constant.bin).
!
! Ago. 2018: modifiche nel timec di convection_kf per eliminare spot di convez. esplicita
! Apr. 2018: introdotta capacita' termica nel calcolo di tskin.
! Nov. 2017: rivista traspiraz. vegetaz. (messa in funz. della rad. netta) e stominr
! (problema eccessiva RH a 2m al tramonto).
!---------------------------------------------------------------
! Directives (ifdef):
! "oper" for the operational version;
! "globo" for the global domain;
! "oldrad" for using the 2004 ECMWF radiation.
! 
!##################################################################################################################

    module mod_model

    include 'dimensions.inc'

! pochva --->
    integer, parameter :: nlev_snow=11
! <---

    logical :: nlbfix, nlcadj, nlana, nlclimate, nlrst
    integer :: nstep, nstep0, ntsbou, nbc, nhist, nhist_zoom, nhist_soil_full, ndrunt, &
               ntsrc, ntswshf, nradm, ntop, nadj, nbl, nstep_sl_filter, nrst
    real(4) :: dtstep, anu2, anu2v, ddamp, hrun, hbound, hist, hist_zoom, hist_soil_full, &
               hdiag, mradconv, mswshf, mslfilter, hrst, &
               zoom_lon_ini, zoom_lon_fin, zoom_lat_ini, zoom_lat_fin
    namelist /model/ dtstep, anu2, anu2v, ddamp, &
                     hrun, hdiag, hbound, &
                     hist, mswshf, hrst, &
                     mradconv, nradm, &
                     ntop, nadj, nbl, & 
                     nlcadj, nlana, nlbfix, nlclimate, &
                     nlrst, &
                     mslfilter, &
                     hist_soil_full, hist_zoom, &
                     zoom_lon_ini, zoom_lon_fin, zoom_lat_ini, zoom_lat_fin

    integer, parameter :: nlon=(gnlon-2)/nprocsx+2, nlat=(gnlat-2)/nprocsy+2
!    integer, parameter :: nst=15, nvt=14
    integer, parameter :: nst=31, nvt=22
    integer, parameter :: nlonm1=nlon-1, nlatm1=nlat-1, nlevp1=nlev+1, nbuffer=2*max(nlon,nlat)*nlev
    integer            :: nyrin, nmonin, ndayin, nhouin, nminin
    integer            :: i_zoom_lon_ini, i_zoom_lon_fin, j_zoom_lat_ini, j_zoom_lat_fin, flag_joint=0
    real(4)            :: alon0, alat0, dlon, dlat, x0, y0, dx, dy, alfa, pzer, d1, d2, d3, d4, dt, alp0, &
                          psmin, co2ppm, rdecli0, reqtim0
    real(4), parameter :: yliv=2834170.5, yliw=333560.5, ylwv=yliv-yliw, tzer=273.15, ezer= 611.,               &
                          rd=287.05, rv=461.51, eps=rd/rv, cpd=1004.6, cvd=cpd-rd, cpv=1869.46, cw=4186.8,      &
                          ci=cw/2., ep=1./eps-1., gamma=cpd/cvd, rdrcp=rd/cpd, a=6371.e3, g=9.807,              &
                          omega=7.292e-5, pi=3.14159265, tkemin=1.e-9
!!!    real(4), parameter :: alsn = .71      ! mean snow albedo (exept glaciers)

    integer, dimension(50)        :: nfdr0, nfdr, nfdrb ! nfdrb for Globo zoom output
    real(4), dimension(200)       :: pdr0, pdr, pdrb ! pdrb for Globo zoom output
    real(4), dimension(nlev)      :: dsig, sigint, dsigalf, sigalf, huc
    real(4), dimension(nlon)      :: bndrel
    real(4), dimension(nlat)      :: hxt, hst, hxv, tangu, tangv
    real(4), dimension(nlevp1)    :: sig
    real(4), dimension(nlevg)     :: d, lev_soil
    real(4), dimension(nlon,nlat) :: ps, pst, fmask, phig, orogvar, tskin, qskin, cloudt, albedo, rgm, rgq, rgmd, &
                                     snow, fsnow, raicon, snocon, rainls, snowls, alsn, tgclim, qgclim, &
                                     water_table_depth, tg_bottom, qg_rel_bottom, qg_rel_surf_approx, &
                                     fveg, lai, proot, veg_rough, &
                                     roscd, qsat, fvegs, fwetl, cw1, cw2, cw3, qveg, deriv, hflux, qflux,  &
                                     frirr, frvis, emisg1, emisg2, fice, alont, alatt, fcorio, t2, &
                                     psb1, psb2, ustar, tstar, rought, psih, bvf, &
                                     u10, v10, q2, qprec, qprecc, qsnfall, dfrirr, dfrvis, solar,                 &
                                    snowfor, iceth, lapse_rate, hsnc, fsnowmax,               &
                                    shf_accum, lhf_accum, qf_accum, albedo_soil_dry, albedo_soil_wet, albedo_veg, &
                                    emiss1_soil_dry, emiss1_soil_wet, emiss1_veg,                                 &
                                    emiss2_soil_dry, emiss2_soil_wet, emiss2_veg,                                 &
                                    fice_rd, iceth_rd
! cumulated variables
    real(4), dimension(nlon,nlat) :: totpre, conpre, snfall, &
                                     cswfl, clwfl, chflux, cqflux, cwvflux, &
                                     cfl_heat_soil_bottom, cfl_water_soil_bottom, &
                                     runoff, runoff_tot, t2min, t2max
    real(4), dimension(nlon,nlat) :: totpre_zoom, conpre_zoom, snfall_zoom, &
                                     cswfl_zoom, clwfl_zoom, chflux_zoom, cqflux_zoom, cwvflux_zoom, &
                                     cfl_heat_soil_bottom_zoom, cfl_water_soil_bottom_zoom, &
                                     runoff_zoom, runoff_tot_zoom, t2min_zoom, t2max_zoom
    real(4), dimension(gnlon,gnlat)      :: gfield, gfield1, gfield2, gfield3
    real(4), dimension(nlon,nlat,nlevg)  :: tg, qg, qgice
    real(4), dimension(nlon,nlat,nlev)   :: tvirt, u, v, t, q, qcw, qci, qrain, qsnow, qhail, phi, omeg, dtdt,    &
                                            dqdt, dqcwdt, dqcidt, dqrndt, dqsndt, div1, div2, fcloud, ut, vt,     &
                                            ub1, vb1, tb1, qb1, qcwb1, qcib1, ub2, vb2, tb2, qb2, qcwb2, qcib2,   &
                                            trd_a, trd_c, tket
    real(4), dimension(nlon,nlat,nlevp1) :: sigdot, phih, tke, lml, rich
    real(4), dimension(nlon,nlat,nst+1)  :: suolo
    real(4), dimension(nlon,nlat,nvt+1)  :: vegeta

! mpi veriables

    integer :: myid=0, iprocs=1, infx=1, infy=1, supx, supy, &
               ip_e=-1, ip_n=-1, ip_s=-1, ip_w=-1, ip_null=-1

! pochva --->

    integer, dimension(nlon,nlat) :: mask_soil, ind_lev_soil_h_bottom, ind_lev_soil_w_bottom
    real(4), dimension(nlon,nlat) :: kturb_surf_m, kturb_surf_h, kturb_surf_q,  flh_specif=0., flh_latent=0., fl_wv=0., &
                                     h_lev_bottom, p_lev_bottom, fl_rain=0., fl_snow=0., fl_rad_tot=0., &
                                     fl_heat_soil_bottom=0., fl_water_soil_bottom=0., fl_runoff=0., fl_runoff_tot=0., &
                                     tg_surf, qg_surf, fice_soil_surf, snow_dirt
    real, dimension(nlon,nlat,20) :: kturb_surf_h_mem, kturb_surf_q_mem
    real(4), dimension(nlon,nlat,nlev_snow) :: snow_lev, snow_t, snow_fice, snow_age, snow_melt_age, snow_dens
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

! For "interpolation" of wind, temp. and humidity at prescribed levels of geometric height (m)
! above and/or below the earth surface

    integer, parameter :: n_std_lev_atm = 5, n_std_lev_soil = 0
    integer, dimension(nlon,nlat) :: n_std_lev_sl
    real(4), dimension(n_std_lev_atm) :: std_lev_atm = (/2., 10., 50., 80., 100./)
    real(4), dimension(n_std_lev_soil) :: std_lev_soil
    real(4), dimension(nlon,nlat,n_std_lev_atm) :: t_std_lev, u_std_lev, v_std_lev, q_std_lev, rh_std_lev, td_std_lev
    real(4), dimension(nlon,nlat,n_std_lev_soil) :: tg_std_lev, qg_std_lev

! costants used to compute partial pressures at saturation

    real(4), parameter :: ccw1=(cpv-cw)/rv, ccw2=ylwv/tzer/rv-ccw1, cci1=(cpv-ci)/rv, cci2=yliv/tzer/rv-cci1

! radintec: mcica=2: McICA clouds; mcica=0: no McICA clouds

    integer, parameter :: nlonr = nlon/2-1, mcica = 0, ntaer = 6
    real(4), dimension(nlon,nlat,nlev)   :: ozon, aerotot
    real(4), dimension(nlon,nlat,nlev,ntaer) :: aerosol

! output of radiation schemes

    real(4), dimension(nlon,nlat)        :: corvis, corirr, gelvis, gelirr
    real(4), dimension(nlon,nlat,nlev)   :: corrdt, geldt

    integer :: imhf=0, imhf_zoom=0, ishf=0, irf=0

! Globo variables and parameters

    integer ip_oppo
    real(4) cpole
    integer, dimension(nlat)      :: nsweep
    real(4), dimension(nlon,nfou) :: snt, cst
    real(4), dimension(nlon,nlat) :: olr, olrtot
    real(4), dimension(nlat)      :: filtt1, filtu1, filtv1

#ifdef globo
    real(4), dimension(nlon,nlat,nlev) :: uf, vf, tf, qf, qcwf, qcif
#endif

    CONTAINS

    function denomin (t1, t2)   ! computes denominators in WAF scheme
    real(4) denomin, t1, t2, zzden
    zzden = t1-t2
    if (abs(zzden).lt.1.e-15) then
    denomin = sign(1.e15,zzden)
    else
    denomin = 1./zzden
    endif
    end function denomin

    end module mod_model
!##################################################################################################################
    program bolam

    use mod_model
    use pochva_scheme, only : pochva

    implicit none

#ifdef mpi
      include 'mpif.h'
      integer      :: status(mpi_status_size)
#endif

    integer      :: ierr, comm
    integer      :: iunit_inp=10, j, n, jadj, jstep, interv, jsday,                &
                    stime1, stime2, countrate, jlon, jlat, jklev, jklev_ref=3, nyrc, nmonc, ndayc, nhouc, nminc,  &
                    nyrc1, nmonc1, ndayc1, nhouc1, nminc1, iday, ihou, imin, ndayr, ndayr1,          &
                    jlo1, jlo2, jla1, jbl, imhf_copy
    real(4)      :: zubar, zvbar, zcz, zdtrdy, zpsav, zcoef, zgamma, zexp, zzpp, zpsb, zub, zvb, zdt, &
                    zalf, zalf2, zbet, zbet2, dtstep2, zdz, z, zhea, ztime, z1, z2, ztb, zqb,         &
                    zqcwb, zqcib, anu41, zesk, zqsat, zt0t, zqgleq, ddamp1, zrid, zalsn, zmin, zmax
    real(4)      :: hlow, hmid, hhig, slow, smid, shig, snwe
    real(kind=8) :: zfac, zx0, zy0, zlatt, zlatv, zlont, zzlatt, zargt, zaarg, zxt, zdpi
    real(4), dimension(gnlat)     :: zfiltgt, zfiltgu, zfiltgv
    real(4), dimension(nlat)      :: zdtrdx
    real(4), dimension(nlon,nlat) :: zu1, zv1
    character(len=15) :: filesoil
    character(len=40) :: command
    character(len=30) :: str, name_file_mhf
    real         :: pd1, pd2, pd4, pd5, pd38, pd39
    integer      :: gnlonr, gnlatr
    real         :: rrf
    real(4), dimension(nlon,nlat) :: totpre_copy, conpre_copy, snfall_copy, &
                                     cswfl_copy, clwfl_copy, chflux_copy, cqflux_copy, cwvflux_copy, &
                                     cfl_heat_soil_bottom_copy, cfl_water_soil_bottom_copy, &
                                     runoff_copy, runoff_tot_copy, t2min_copy, t2max_copy

! pochva --->
    real         :: zdlev, dlev_base=10., dlev_min=5.
    integer :: n_outpoint=0, i, k, flag_call_pochva=1
    integer, dimension(10) :: &
 i_outpoint=    (/ 14,   54,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1/), &
 j_outpoint=    (/ 49,   42,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1/), &
 iproc_outpoint=(/ 51,   23,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1/)
    character (len=30) :: file_outpoint
! <---

! initialize the scalar case

    nstep0  = 1
    supx    = nlon
    supy    = nlat

! initialize mpi environment

#ifdef mpi
      comm = mpi_comm_world
      call mpi_init(ierr)
      if (ierr .ne. 0) then
        print *,'error starting mpi program.  terminating.'
        call mpi_abort(comm, ierr)
      endif

      call mpi_comm_size(comm, iprocs, ierr)  ! find out number of processes
#endif

    if(iprocs.ne.nprocsx*nprocsy) then
      print *,iprocs,nprocsx,nprocsy,'Error in the definition of number of processes: terminating.'
#ifdef mpi
      call mpi_abort(comm, ierr)
#endif
      stop
    endif

#ifdef mpi
      if (gnlat-2.ne.(nlat-2)*nprocsy.or.gnlon-2.ne.(nlon-2)*nprocsx) then
        print *,'errors in domain decomposition:  terminating.'
        call mpi_abort(comm, ierr)
      endif
#endif

    if (mod(nlon,2).ne.0.or.mod(nlat,2).ne.0) then
      print *,'nlon or nlat not even numbers: terminating.'   ! for ECMWF radiation computation
#ifdef mpi
        call mpi_abort(comm, ierr)
#endif
    endif

#ifdef mpi

      call mpi_comm_rank(comm, myid, ierr) ! find out my process id
      if (myid.eq.0) then
        print *

#ifdef globo
        print*,'#### GLOBO parallel version ####'
#else
        print*,'#### BOLAM parallel version ####'
#endif

      write(*,6000) nprocsx, nprocsy, gnlon, gnlat
 6000 format(/,' Number of processors =',i2,' x',i2,/,' Global x dim. =',i5,/,' Global y dim. =',i5)
      endif

#endif

! set exterior limits of the local domain

#ifdef mpi

      infx = 1+(myid/nprocsy)*(nlon-2)
      supx = infx+nlon-1
      infy = 1+(myid-(myid/nprocsy)*nprocsy)*(nlat-2)
      supy = infy+nlat-1

      write (*,6001) myid, nlon, nlat, nlev, infx, supx, infy, supy
   6001 format(' myid =',i4,'  -  nlon,nlat,nlev =',3i5,' -  position in global domain',4i4)

! define id's of neighbour processes

      ip_e = myid+nprocsy
      ip_w = myid-nprocsy

#ifdef globo
        if (ip_e.gt.iprocs-1) ip_e = ip_e - iprocs
        if (ip_w.lt.0       ) ip_w = ip_w + iprocs
#else
        if (ip_e.gt.iprocs-1) ip_e = mpi_proc_null
        if (ip_w.lt.0       ) ip_w = mpi_proc_null
#endif

      ip_s = myid-1
      ip_n = myid+1
      if (mod(myid,nprocsy).eq.0        ) ip_s = mpi_proc_null
      if (mod(myid,nprocsy).eq.nprocsy-1) ip_n = mpi_proc_null
      ip_null = mpi_proc_null
      ip_oppo = myid + iprocs/2   ! for Globo only
      if (ip_oppo.gt.iprocs-1) ip_oppo = ip_oppo - iprocs

!      call system('hostname')

#endif

!------------------------------
! read parameters from namelist
!------------------------------

    open (iunit_inp, file='bolam.inp', status='old')
    read (iunit_inp, model)
    close (iunit_inp)

#ifndef rad_ecmwf
    nradm = min(nradm, 1) ! No ECMWF radiation scheme
#endif

    if (myid.eq.0) then
    print*
    print *,'Parameters of bolam (bolam.inp):'
    print model
    endif
    nstep  = hrun*3600./dtstep+.5
    ntsbou = hbound*3600./dtstep+.5
    nbc    = (nstep-1)/ntsbou+2
    nhist  = hist*3600./dtstep+.5
    if (hist < 0.) nhist = -1
    nhist_zoom  = int(hist_zoom*3600./dtstep+.5)
    if (hist_zoom < 0.) nhist_zoom = -1
    ntswshf= mswshf*60./dtstep+.5
    if (mswshf < 0.) ntswshf = -1
    nhist_soil_full  = int(hist_soil_full*3600./dtstep+.5)
    if (hist_soil_full < 0.) nhist_soil_full = -1
    ndrunt = max( int(hdiag*3600./dtstep+.5), 1)
    ntsrc  = mradconv*60./dtstep+.5
    nstep_sl_filter = int(mslfilter*60./dtstep)
    nrst = hrst*3600./dtstep+.5
    flag_joint = 0

!-------------------------------------
! read initial condition
!-------------------------------------

    if (nlrst) then
      call rdrf (1)
      interv = nstep0/ntsbou + 1
    else
      interv = 1
    endif

#ifdef globo
    interv = 1
#endif

    call rdmhf_atm (interv, 1)
    call rdmhf_soil (interv, 1)

    call rd_param_const

    if (nlrst) call rdrf (1) 

!*************************************************
!  year/month/day/hour/min. of initial condition
!*************************************************

    nyrin  = nfdr0(5)
    nmonin = nfdr0(6)
    ndayin = nfdr0(7)
    nhouin = nfdr0(8)
    nminin = nfdr0(9)
    iday  =  nfdr0(10)
    ihou  =  nfdr0(11)
    imin  =  nfdr0(12)

    

    call calendar (nyrin, nmonin, ndayin, nhouin, nminin, iday, ihou, imin, &
                   nyrc, nmonc, ndayc, nhouc, nminc, ndayr)

!  update initial condition date (if initial condition is a restart)

    nminin = nminc
    nhouin = nhouc
    ndayin = ndayc
    nmonin = nmonc
    nyrin  = nyrc

    if (myid.eq.0) then
    print*, 'Initial date (year, month, day, hour):'
    write(*,'(12x,I8.2,3I6.2)') nyrin, nmonin, ndayin, nhouin
    endif

! CO2 concentration defined as a function of the year, assuming a linear trend

    co2ppm = 400. + 2.*(nyrin - 2015)
    co2ppm = max (co2ppm, 280.)
    if(myid.eq.0) then
    write(*, '(a, f7.2)') ' CO2 concentration in ppm =', co2ppm
    print*
    endif

!******************************************************
! read soil variables tg and qg (except the deepest layer)
! and snow cover from file if it exists;
! otherwise use present analysis soil variables
!******************************************************

    jsday = 86400./dtstep+.5

!******************************************************
!  horizontal grid parameters
!******************************************************

    dlat   = pdr0(1)
    dlon   = pdr0(2)
    alat0  = pdr0(4) + (infy-1)*dlat  ! latit. of the 1st v-point of the local domain
    alon0  = pdr0(5) + (infx-1)*dlon  ! longit. of the 1st t-point of the local domain
    y0     = pdr0(38)
    x0     = pdr0(39)
    dx     = a*dlon*pi/180.
    dy     = a*dlat*pi/180.
    dt     = dtstep/float(nadj)

#ifdef globo
    i_zoom_lon_ini = int((zoom_lon_ini- 0.)/dlon)+1
    i_zoom_lon_fin = int((zoom_lon_fin- 0.)/dlon)+1
    j_zoom_lat_ini = int((zoom_lat_ini+90.)/dlat)+1
    j_zoom_lat_fin = int((zoom_lat_fin+90.)/dlat)+1
    if (i_zoom_lon_ini > i_zoom_lon_fin) flag_joint=1
#endif

    zfac = dabs(dacos(-1.d0))/180.d0
    zx0  = dble(x0)*zfac
    zy0  = dble(y0)*zfac
    do jlat = 1, nlat
    zlatv = (dble(alat0)+                (jlat-1)*dble(dlat))*zfac
    zlatt = (dble(alat0)+dble(dlat)/2.d0+(jlat-1)*dble(dlat))*zfac
    hxt(jlat)    = dcos(zlatt)
    hst(jlat)    = dsin(zlatt)
    hxv(jlat)    = dcos(zlatv)
    tangu(jlat)  = hst(jlat)/hxt(jlat)/a
    tangv(jlat)  = dsin(zlatv)/hxv(jlat)/a
    zdtrdx(jlat) = dt/hxt(jlat)/dx

#ifdef globo
      do jlon = 1, nlon
      fcorio(jlon,jlat) = 2.*omega*hst(jlat)
      alont(jlon,jlat) = alon0+float(jlon-1)*dlon
      enddo
      alatt(:,jlat) = alat0+dlat*0.5+float(jlat-1)*dlat
#else
      do jlon = 1, nlon
      zlont  = (dble(alon0)+(jlon-1)*dble(dlon))*zfac
      zzlatt = 1.d0/zfac*dasin( dcos(zy0)*dsin(zlatt) + dsin(zy0)*dcos(zlatt)*dcos(zlont))
      zargt  = -dsin(zlatt)*dsin(zy0)+dcos(zy0)*dcos(zlatt)*dcos(zlont)
      zaarg  = zargt/dcos(zfac*zzlatt)
      alatt(jlon,jlat) = zzlatt
      if (zaarg.lt.-1..and.zaarg.gt.-1.00001) zaarg = -1.d0
      if (zaarg.gt. 1..and.zaarg.lt. 1.00001) zaarg =  1.d0
        if (zlont.lt.0.d0) then
        alont(jlon,jlat) = 1.d0/zfac*(zx0-dacos(zaarg))
        else
        alont(jlon,jlat) = 1.d0/zfac*(zx0+dacos(zaarg))
        endif
      fcorio(jlon,jlat) = 2.*omega*sin(pi*alatt(jlon,jlat)/180.)
      enddo
#endif

    enddo

! reduction of snow (in the Alpine region only) - may be necessary from time to time...

!   do jlat=2,nlat-1
!   do jlon=2,nlon-1
!   if(alatt(jlon,jlat).gt.43.8.and.alatt(jlon,jlat).lt.48.2.and.alont(jlon,jlat).gt.5..and.alont(jlon,jlat).lt.16.5) then
!   snow(jlon,jlat) = max(0., 0.55*snow(jlon,jlat) - 0.005)
!   endif
!   enddo
!   enddo

!******************************************************
!  definition of vertical levels
!******************************************************

    pzer   = pdr0(36)
    alfa   = pdr0(37)

    if (myid.eq.0) print*, 'Pzer =', pzer, 'alfa =', alfa

    sig(1) = 0.
    do jklev = 2, nlevp1
    sig(jklev) = pdr0(38+jklev)
    enddo
    do jklev =1, nlev
    dsig(jklev)    = sig(jklev+1)-sig(jklev)
    sigint(jklev)  = .5*(sig(jklev+1)+sig(jklev))
    sigalf(jklev)  = sigint(jklev)**3 * (1.+alfa*(1.-sigint(jklev))*(2.-sigint(jklev)))
    dsigalf(jklev) = dsig(jklev)*sigint(jklev)**2 * (3.+alfa*(6.-12.*sigint(jklev)+5.*sigint(jklev)**2))
    enddo

!  definition of cloud water/ice for the initial condition

    call defclt (qcw, qci, t, nlon, nlat, nlev, tzer, 1)

!*****************************************************************************************
!  initialization of boundary conditions
!  definition of boundary relaxation coefficients  - relax requires nbl to be a power of 2
!*****************************************************************************************

#ifndef globo
    call relax (nbl, .01, 1., bndrel)
    do jlon = nbl, 2, -1
    bndrel(jlon) = bndrel(jlon-1)
    enddo
    bndrel(1) = 1.

    psb1  = ps
    ub1   = u
    vb1   = v
    tb1   = t
    qb1   = q
    qcwb1 = qcw
    qcib1 = qci

    if (.not.nlbfix) then
      if (nlrst) then
        interv = nstep0/ntsbou + 2
      else
        interv = 2
      endif
      call rdmhf_atm (interv, 0)
      call defclt (qcwb2, qcib2, tb2, nlon, nlat, nlev, tzer, 2)
    endif
#endif

!  define descriptor records

    nfdr  = nfdr0
    pdr   = pdr0

!  nfdr

    nfdr(1)  = 1
    nfdr(5)  = nyrin
    nfdr(6)  = nmonin
    nfdr(7)  = ndayin
    nfdr(8)  = nhouin
    nfdr(9)  = nminin
    nfdr(13) = npolcap
    nfdr(14) = nbl
!   nfdr(15) contains no. of soil levels nlevg
    nfdr(16) = nstep
    nfdr(17) = nhist
!!!    if (nhist < 0.and.nhist_zoom > 0) nfdr(17) = nhist_zoom
    nfdr(18) = ntswshf

!  pdr

    pdr(3)  = dtstep

    nfdrb(:)  = nfdr(:)
    nfdrb(17) = nhist_zoom
    pdrb(:)   = pdr(:)
#ifdef globo
    nfdrb(1) = 8
    if (flag_joint == 1) then
    nfdrb(2) = i_zoom_lon_fin+gnlon-1-i_zoom_lon_ini
    else
    nfdrb(2) = i_zoom_lon_fin-i_zoom_lon_ini+1
    endif
    nfdrb(3) = j_zoom_lat_fin-j_zoom_lat_ini+1
    pdrb = 0.
    do jklev = 1, nlev
    pdrb(39+jklev) = sig(jklev+1)
    enddo
    pdrb(1) = dlat
    pdrb(2) = dlon
    pdrb(3) = dtstep
    pdrb(4) = float(j_zoom_lat_ini-1)*dlat-90.-dlat/2.
    if (flag_joint == 1) then
    pdrb(5) = float(i_zoom_lon_ini-1)*dlon-360.
    else
    pdrb(5) = float(i_zoom_lon_ini-1)*dlon
    endif
    pdrb(6:5+nlevg) = pdr0(6:5+nlevg)
    pdrb(36) = pzer
    pdrb(37) = alfa
#endif

!***********************************************************************
!  initializations
!***********************************************************************

    alp0   = log(1.e5)
    zdtrdy = dt/dy
    zcz    = dt*rd/cpd

!  depth soil layers

    lev_soil(1:nlevg) = pdr0(6:5+nlevg)

    d(1) = lev_soil(1)*2.
    do jklev = 2, nlevg
      d(jklev)=2.*(lev_soil(jklev)-lev_soil(jklev-1))-d(jklev-1)
    enddo

    if (.not.nlrst) then

    nstep0 = 1

!  physical parameters of soil scheme

    do jlat = 1, nlat
    do jlon = 1, nlon

    alsn(jlon,jlat) = .71      ! mean snow albedo (exept glaciers) will be changed by soil processes scheme

    fice(jlon,jlat) = max(fice(jlon,jlat), 0.)
    fice(jlon,jlat) = min(fice(jlon,jlat), 1.)
    if (fmask(jlon,jlat).lt.0.5) fice(jlon,jlat) = 0.

    qg(jlon,jlat,1:nlevg) = max( min( qg(jlon,jlat,1:nlevg), qsoil_max(jlon,jlat,1:nlevg)), qsoil_min(jlon,jlat,1:nlevg))

! pochva --->

    qg_surf(jlon,jlat) = qg(jlon,jlat,1)

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

!    snow(jlon,jlat) = snow(jlon,jlat)*1.e3 ! m of equiv. water ---> kg/m/m
!
!    if (fmask(jlon,jlat) < 0.5.and.int(suolo(jlon,jlat,1)) == nst-1) then
!! glacier
!      snow(jlon,jlat) = max (50., snow(jlon,jlat)) ! snow (kg/m/m)
!    endif
!
!    if (fmask(jlon,jlat) >= 0.5.and.fice(jlon,jlat) >= 0.8) then
!! sea ice
!      snow(jlon,jlat) = max (50., snow(jlon,jlat)) ! snow (kg/m/m)
!    endif
!
!    snow_lev(jlon,jlat,:) = valmiss
!    snow_t(jlon,jlat,:) = valmiss
!    snow_fice(jlon,jlat,:) = valmiss
!    snow_age(jlon,jlat,:) = valmiss
!    snow_melt_age(jlon,jlat,:) = valmiss
!    snow_dens(jlon,jlat,:) = valmiss
!
!    if (snow(jlon,jlat) > 0.) then
!      tskin(jlon,jlat) = min(tskin(jlon,jlat), tzer-0.1)
!      zdlev=dlev_base
!      do while (.true.)
!        if (zdlev*float(nlev_snow-1) >= snow(jlon,jlat)) exit
!        zdlev=zdlev+dlev_base
!      enddo
!      snow_lev(jlon,jlat,1)=0.
!      if (snow(jlon,jlat) >= zdlev+dlev_min) THEN
!        n=int(snow(jlon,jlat)/zdlev)+1
!        if (snow(jlon,jlat)-zdlev*float(n-1) >= dlev_min) n=n+1
!        do jklev=2,n
!          snow_lev(jlon,jlat,jklev)=snow(jlon,jlat)-float(n-jklev)*zdlev
!        enddo
!        n=n-1
!      else
!        n=1
!        snow_lev(jlon,jlat,2)=snow(jlon,jlat)
!      endif
!      snow_t(jlon,jlat,1) = tskin(jlon,jlat)
!      if (n < 2) then
!        snow_t(jlon,jlat,2) = min(tg(jlon,jlat,1), tzer-0.1)
!      else
!        snow_t(jlon,jlat,n+1) = min(tg(jlon,jlat,1), tzer-0.1)
!        z1=1./float(n)
!        do jklev=2,n
!          snow_t(jlon,jlat,jklev) = snow_t(jlon,jlat,1)*(1.-z1*float(jklev-1)) + &
! snow_t(jlon,jlat,n+1)*z1*float(jklev-1)
!        enddo
!      endif
!      snow_fice(jlon,jlat,1:n+1) = 1.
!      snow_age(jlon,jlat,1:n+1) = 10.
!      snow_melt_age(jlon,jlat,1:n+1) = 0.
!    endif
!
!    if (fmask(jlon,jlat) < 0.5.or.fice(jlon,jlat) >= 0.8) then
!! glacier
!      if (int(suolo(jlon,jlat,1)) == nst-1) then
!        do jklev = 1, nlevg
!          tg(jlon,jlat,jklev) = min(tg(jlon,jlat,jklev), tzer-0.1)
!        enddo
!        snow_age(jlon,jlat,1:n+1) = 180. ! snow age (days)
!        snow_melt_age(jlon,jlat,1:n+1) = 30. ! snow age (days)
!!print *,' Glacier ',myid,jlon,jlat
!      endif
!! sea ice
!      if (fmask(jlon,jlat) > 0.5) then
!        psi_soil(jlon,jlat,:) = 2.2 ! thermal conductivity of ice
!        snow_age(jlon,jlat,1:n+1) = 180. ! snow age (days)
!        do jklev = 1, nlevg
!          if (lev_soil(jklev) > iceth(jlon,jlat)) exit
!        enddo
!        ind_lev_soil_h_bottom(jlon,jlat) = max(jklev-1, 3)
!        ind_lev_soil_w_bottom(jlon,jlat) = 0
!        do jklev = 1, ind_lev_soil_h_bottom(jlon,jlat)
!          tg(jlon,jlat,jklev) = min(tg(jlon,jlat,jklev), tzer-0.1)
!        enddo
!        albedo_soil_dry(jlon,jlat) = 0.4
!        albedo_soil_wet(jlon,jlat) = 0.4
!        albedo(jlon,jlat) = 0.4
!      endif
!    endif
!
!    snow(jlon,jlat) = snow(jlon,jlat)*1.e-3 ! kg/m/m ---> m of equiv. water

! <---

    enddo
    enddo

! ---> pochva
    mask_soil(1,   1:nlat) = 0
    mask_soil(nlon,1:nlat) = 0
    mask_soil(1:nlon,   1) = 0
    mask_soil(1:nlon,nlat) = 0
! <---

    t2        = tskin
    geldt     = 0.
    corrdt    = 0.
    qrain     = 0.
    qsnow     = 0.
    qhail     = 0.
    tke       = tkemin ! turbulent kinetic energy
    lml       = 1.e-6  ! mixing length
    trd_a     = 0.
    trd_c     = 0.
    totpre    = 0.     ! total rain+snow (kg/m2) accumulated between two checkpoints
    conpre    = 0.     ! convective rain+snow (kg/m2) accumulated between two checkpoints
    snfall    = 0.     ! snow fall (kg/m2 of equivalent water) accumulated between two checkpoints
    raicon    = 0.     ! rain (convective) in ntsrc timesteps (kg/m2)
    snocon    = 0.     ! snow (convective) in ntsrc timesteps (kg/m2)
    rainls    = 0.     ! rain (large scale) in one timestep (kg/m2)
    snowls    = 0.
    qprec     = 0.
    qprecc    = 0.
    qsnfall   = 0.
    hflux     = 0.
    qflux     = 0.
    frvis     = 0.
    frirr     = 0.
    gelvis    = 0.
    gelirr    = 0.
    corvis    = 0.
    corirr    = 0.
    chflux    = 0.
    cqflux    = 0.
    cwvflux   = 0.
    cswfl     = 0.
    solar     = 0.
    shf_accum = 0.
    lhf_accum = 0.
    qf_accum  = 0.
    clwfl     = 0.
    fl_heat_soil_bottom  =0.
    fl_water_soil_bottom =0.
    cfl_heat_soil_bottom  =0.
    cfl_water_soil_bottom =0.
    runoff      =0.
    runoff_tot  =0.
    fl_runoff   =0.
    fl_runoff_tot=0.
    fcloud    = 1.e-10 ! box fraction covered by cloud (do not set to 0.)
    roscd     = 1.e-2
    rich      = 10.

    snow = max (snow, 0.)   ! snow height in mm of equivalent water
    do jlat=1,nlat
    do jlon=1,nlon
      if (snow(jlon,jlat) > 0.) then
        fsnow(jlon,jlat) = min( ((snow(jlon,jlat)/5.e-3)**0.5), 1.0)
      else
        fsnow(jlon,jlat) = 0.
      endif
    enddo
    enddo
!    fsnow(:,:) = min((snow(:,:)/hsnc(:,:))**0.67, fsnowmax(:,:)) ! snow fraction (re-defined in subr. soil)
    lai  = max (lai, 1.e-8)
    t2min(:,:) = t(:,:,nlev)
    t2max(:,:) = t(:,:,nlev)
    t2min(2:nlonm1,2:nlatm1) = 999.
    t2max(2:nlonm1,2:nlatm1) = 0.
    sigdot(:,:,1     ) = 0.
    sigdot(:,:,nlevp1) = 0.
    phih  (:,:,nlevp1) = phig(:,:)
    phih  (:,:,     1) = 1.e6
    rgmd(:,:) = rgm(:,:)

    totpre_zoom = 0.
    conpre_zoom = 0.
    snfall_zoom = 0.
    cswfl_zoom  = 0.
    clwfl_zoom  = 0.
    chflux_zoom = 0.
    cqflux_zoom = 0.
    cwvflux_zoom = 0.
    cfl_heat_soil_bottom_zoom = 0.
    cfl_water_soil_bottom_zoom = 0.
    runoff_zoom = 0.
    runoff_tot_zoom = 0.
    t2min_zoom = t2min
    t2max_zoom = t2max

    else ! nlrst

! Reading of restart file

      call rdrf (0)
      call rdrf_pochva(nstep0, nprocsx, nprocsy, myid, gfield, gnlon, gnlat, dtstep)
      flag_pochva_initial = 0

    endif ! .not.nlrst

!  definition of critical relative humidity for cloud formation
!  depending on level (do not set huc = 1. to avoid overflow)

    hlow = 0.88
    hmid = 0.87
    hhig = 0.92
    slow = 0.900
    smid = 0.550
    shig = 0.250
    zalf =(hlow-hmid)/(1.-slow)
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

    if (.not.nlrst) then

    tvirt = t*(1. +ep*q -qcw -qci -qrain -qsnow)

!  Computation of geopotential

    call phicomp

!  definition of albedo, emisg1, emisg2

 if (myid == 0) then
   print *
   print *,"Renewal of radiation parameters of the surface"
 endif
!!!    call surfradpar
    call surfradpar(fmask, fice, qg_surf, vegeta, qsoil_max, qsoil_min, &
 albedo_soil_dry, albedo_soil_wet, emiss1_soil_dry, emiss1_soil_wet, emiss2_soil_dry, emiss2_soil_wet, &
 fveg, albedo_veg, emiss1_veg, emiss2_veg, &
 nlon, nlat, nlevg, nvt, 1, albedo, emisg1, emisg2)

    endif

! <augh,max: tue sep 30 09:15:42 cest 2008>
#ifdef chem
    call init_chem (myid)
#endif
! </augh,max>

!******************************************************
!  Globo variables and constants
!******************************************************

    cpole = sin(dlat*pi/360.)/(a*float(gnlon-2)*(1.-cos(dlat*pi/360.)))

!  Number of sweeps in zonal advection and polar filter

#ifdef globo
    nsweep(:) = max( 1, nint(0.7071/(hxt(:)+1.e-5)) )
    if (mod(myid,nprocsy).eq.0        ) nsweep(1   ) = 1
    if (mod(myid,nprocsy).eq.nprocsy-1) nsweep(nlat) = 1
#else
    nsweep(:) = max( 1, nint(1./(hxt(:)+1.e-5)) )
#endif

#ifdef globo
    zdpi = dabs(dacos(-1.d0))
    do n = 1, nfou
    do jlon = 1, nlon
    zxt = (infx+jlon-2) * 360.d0/(gnlon-2) * zdpi/180.d0
    snt(jlon,n) = dsin(zxt*n)
    cst(jlon,n) = dcos(zxt*n)
    enddo
    enddo

    zfiltgt = 0
    zfiltgu = 0
    zfiltgv = 0
    do jlat = 2, gnlat/2
    zfiltgt(jlat) = 4.*cos(45.*pi/180.)**2/(pi*cos(-pi/2.+(jlat-1)*dlat*pi/180.))**2 ! see 2011 paper on Weather and For.
    zfiltgu(jlat) = 4./(dlon*pi/180.)**2/(3.14*(jlat-1)+1.)**2          ! approx.
    zfiltgv(jlat) = 4./(dlon*pi/180.)**2/(3.14*(jlat-1)+1.-3.14/2.)**2  ! approx.

    zfiltgt(jlat) = max (zfiltgt(jlat)-4./pi**2, 0.)  ! no filter below 45 degrees
    zfiltgu(jlat) = max (zfiltgu(jlat)-4./pi**2, 0.)
    zfiltgv(jlat) = max (zfiltgv(jlat)-4./pi**2, 0.)
    enddo
    do jlat = 1, gnlat/2-1
    zfiltgt(gnlat-jlat  ) = zfiltgt(jlat+1)
    zfiltgu(gnlat-jlat  ) = zfiltgu(jlat+1)
    zfiltgv(gnlat-jlat+1) = zfiltgv(jlat+1)
    enddo

    do jlat = 1, nlat
    filtt1(jlat) = zfiltgt(infy+jlat-1)
    filtu1(jlat) = zfiltgu(infy+jlat-1)
    filtv1(jlat) = zfiltgv(infy+jlat-1)
    enddo
    if (mod(myid,nprocsy).eq.0        ) filtt1(1             :npolcap+1) = 0
    if (mod(myid,nprocsy).eq.0        ) filtu1(1             :npolcap+1) = 0
    if (mod(myid,nprocsy).eq.0        ) filtv1(1             :npolcap+1) = 0
    if (mod(myid,nprocsy).eq.nprocsy-1) filtt1(nlat-npolcap  :nlat     ) = 0
    if (mod(myid,nprocsy).eq.nprocsy-1) filtu1(nlat-npolcap  :nlat     ) = 0
    if (mod(myid,nprocsy).eq.nprocsy-1) filtv1(nlat-npolcap+1:nlat     ) = 0

!  inizialization at poles for possible errors in initial fields

    call filt2uv (1.)
    call polaveruv
    call polavert (t)
    call polavert (q)
    q = max (q, 1.e-12)
#endif

    if (.not.nlrst) then

    call cloudfr
    call ccloud

#ifdef globo
    if (nhist_zoom > 0) then
      call wr_param_const (nfdrb, pdrb, i_zoom_lon_ini, i_zoom_lon_fin, j_zoom_lat_ini, j_zoom_lat_fin, flag_joint)
    endif
#endif

!  save initial condition

    if (nlana.and.myid == 0) then
      imhf=imhf+1
      imhf_zoom=imhf_zoom+1
    endif

#ifdef globo
    call runout_g (0, nyrc, nmonc, ndayc, nhouc, nminc)  ! run time diagnostics
    if (nlana) then
      if (nhist > 0) then
        name_file_mhf="globo"
        call wrmhf_atm (name_file_mhf, nfdr, pdr, 1, gnlon, 1, gnlat, 0)
        call wrmhf_soil(name_file_mhf, nfdr, pdr, 1, gnlon, 1, gnlat, 0, 0)
      endif
      if (nhist_zoom > 0) then
        name_file_mhf="globo_zoom"
        call wrmhf_atm (name_file_mhf, nfdrb, pdrb, i_zoom_lon_ini, i_zoom_lon_fin, j_zoom_lat_ini, j_zoom_lat_fin, flag_joint)
        call wrmhf_soil(name_file_mhf, nfdrb, pdrb, i_zoom_lon_ini, i_zoom_lon_fin, j_zoom_lat_ini, j_zoom_lat_fin, flag_joint, 0)
      endif
    endif
#else
    call runout_b (0)  ! run time diagnostics
    if (nlana) then
      name_file_mhf="bolam"
      call wrmhf_atm (name_file_mhf, nfdr, pdr, 1, gnlon, 1, gnlat, 0)
      call wrmhf_soil(name_file_mhf, nfdr, pdr, 1, gnlon, 1, gnlat, 0, 0)
    endif
#endif

    else ! nlrst

    if (myid == 0) print *,'--- Restart ---'
#ifdef globo
    call runout_g (nstep0, nyrc, nmonc, ndayc, nhouc, nminc)  ! run time diagnostics
#else
    call runout_b (nstep0)  ! run time diagnostics
#endif
    if (myid == 0) print *,'---------------'
    imhf = (nstep0-1)/nhist
    imhf_zoom = (nstep0-1)/nhist_zoom
    if (nlana) then
      imhf = imhf+1
      imhf_zoom = imhf_zoom+1
    endif
    ishf = (nstep0-1)/ntswshf+1
    irf = nstep0/nrst

    endif ! .not.nlrst

    call system_clock (stime1, countrate)  ! .....elapsed time

!***********************************************************************
!  time integration begins....
!***********************************************************************

    do 1000 jstep = nstep0, nstep

#ifdef globo
    uf = u
    vf = v
    tf = t
    qf = q
    qcwf = qcw
    qcif = qci
#endif

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

    tvirt = t*(1. +ep*q -qcw -qci -qrain -qsnow)

!***********************************************************************
    do 200 jadj = 1, nadj
!***********************************************************************

!  adjustment step -- only terms describing gravity waves are considered
!  forward integration of surface pressure
!  divergence of horizontal motion, vertical velocity and tendency of ps is made in sub. diverg

#ifdef mpi
      call u_ghost (v(:,2,:), ip_s, v(:,nlat,:), ip_n, nlon*nlev)
      if (nprocsx.eq.1) then

#ifdef globo
          u(1,:,:) = u(nlonm1,:,:)
          u(nlon,:,:) = u(2,:,:)
#endif

      else
        call u_ghost (u(nlonm1,:,:), ip_e, u(1,:,:), ip_w, nlat*nlev)
      endif
#else
#ifdef globo
        u(1,:,:) = u(nlonm1,:,:)
        u(nlon,:,:) = u(2,:,:)
#endif
#endif

    call diverg (jadj)

!  new surface pressure

    ps = ps + dt*pst

!  computation of new virtual temperature (omega-alpha term)

    tvirt = tvirt * (1.+zcz*omeg)

!  ghostline update to integrate momentum equations

#ifdef mpi

      call u_ghost (ps(:,2     ), ip_s, ps(:,nlat), ip_n, nlon)
      call u_ghost (ps(:,nlatm1), ip_n, ps(:,1   ), ip_s, nlon)
      call u_ghost (tvirt(:,nlatm1,:), ip_n, tvirt(:,1,:), ip_s, nlon*nlev)
      call u_ghost (div1 (:,nlatm1,:), ip_n, div1 (:,1,:), ip_s, nlon*nlev)
      if (nprocsx.eq.1) then

#ifdef globo
        ps(1   ,:) = ps(nlonm1,:)
        ps(nlon,:) = ps(2     ,:)
        tvirt(1   ,:,:) = tvirt(nlonm1,:,:)
        tvirt(nlon,:,:) = tvirt(2     ,:,:)
        div1 (1   ,:,:) = div1 (nlonm1,:,:)
        div1 (nlon,:,:) = div1 (2     ,:,:)
#endif

      else
        call u_ghost (ps(2     ,:), ip_w, ps(nlon,:), ip_e, nlat)
        call u_ghost (ps(nlonm1,:), ip_e, ps(1   ,:), ip_w, nlat)
        call u_ghost (tvirt(2,:,:), ip_w, tvirt(nlon,:,:), ip_e, nlat*nlev)
        call u_ghost (div1 (2,:,:), ip_w, div1 (nlon,:,:), ip_e, nlat*nlev)
      endif
#else
#ifdef globo
      ps(1   ,:) = ps(nlonm1,:)
      ps(nlon,:) = ps(2     ,:)
      tvirt(1   ,:,:) = tvirt(nlonm1,:,:)
      tvirt(nlon,:,:) = tvirt(2     ,:,:)
      div1 (1   ,:,:) = div1 (nlonm1,:,:)
      div1 (nlon,:,:) = div1 (2     ,:,:)
#endif
#endif

!  Computation of geopotential

    call phicomp

!  backward integration of momentum equations
!  divergence damping
!  the Coriolis terms must be included here

    do jklev = 1, nlev
    ddamp1 = (ddamp + (1.-ddamp)/float(jklev+2)-(1.-ddamp)/float(nlev+2))*.125*dy**2/dt ! sponge at the top of the atmosphere
!    ddamp1 = ddamp*.125*dy**2/dt ! no sponge layer
    do jlat = 2, nlatm1
    do jlon = 2, nlonm1
    zvbar = .25*(v(jlon,jlat,jklev)+v(jlon+1,jlat,jklev)+v(jlon,jlat+1,jklev)+v(jlon+1,jlat+1,jklev))
    zpsav = .5*(ps(jlon+1,jlat)+ps(jlon,jlat))
    zcoef = .5*rd*(ps(jlon+1,jlat)-ps(jlon,jlat))*sigalf(jklev)/(pzer*sigint(jklev)-(pzer-zpsav)*sigalf(jklev))
    u(jlon,jlat,jklev) = u(jlon,jlat,jklev) - ( (tvirt(jlon,jlat,jklev)+tvirt(jlon+1,jlat,jklev))*zcoef &
                       + phi(jlon+1,jlat,jklev)-phi(jlon,jlat,jklev) )*zdtrdx(jlat)                     &
                       + .5*dt*(fcorio(jlon,jlat)+fcorio(jlon+1,jlat))*zvbar                            &
                       + ddamp1*(div1(jlon+1,jlat,jklev)-div1(jlon,jlat,jklev))*zdtrdx(jlat)
    enddo
    enddo
    enddo

    do jklev = 1, nlev
    ddamp1 = (ddamp + (1.-ddamp)/float(jklev+2)-(1.-ddamp)/float(nlev+2))*.125*dy**2/dt ! sponge at the top of the atmosphere
!    ddamp1 = ddamp*.125*dy**2/dt ! no sponge layer
    do jlat = 2, nlatm1
    do jlon = 2, nlonm1
    zubar = .25*(u(jlon,jlat,jklev)+u(jlon-1,jlat,jklev)+u(jlon,jlat-1,jklev)+u(jlon-1,jlat-1,jklev))
    zpsav = .5*(ps(jlon,jlat)+ps(jlon,jlat-1))
    zcoef = .5*rd*(ps(jlon,jlat)-ps(jlon,jlat-1))*sigalf(jklev)/(pzer*sigint(jklev)-(pzer-zpsav)*sigalf(jklev))
    v(jlon,jlat,jklev) = v(jlon,jlat,jklev) - ( (tvirt(jlon,jlat,jklev)+tvirt(jlon,jlat-1,jklev))*zcoef &
                       + phi(jlon,jlat,jklev)-phi(jlon,jlat-1,jklev) )*zdtrdy                           &
                       - .5*dt*(fcorio(jlon,jlat)+fcorio(jlon,jlat-1))*zubar                            &
                       + ddamp1*(div1(jlon,jlat,jklev)-div1(jlon,jlat-1,jklev))*zdtrdy
    enddo
    enddo
    enddo

!***********************************************************************
 200 continue    ! end of gravity mode computation
!***********************************************************************

    call wafone (tvirt, dtstep)
    call wafone (q,     dtstep)
    call wafone (qcw,   dtstep)
    call wafone (qci,   dtstep)

#ifndef globo
    call wafone (qrain, dtstep)
    call wafone (qsnow, dtstep)
    call wafone (qhail, dtstep)

!---------------------------------
! TKE from half to integer levels
!---------------------------------

    do jklev = 2, nlev-1
    tket(:,:,jklev) = .5625*(tke(:,:,jklev+1)+tke(:,:,jklev))-.0625*(tke(:,:,jklev+2)+tke(:,:,jklev-1))
    enddo
    tket(:,:,1   ) = .5*(tke(:,:,2     )+tke(:,:,1   ))
    tket(:,:,nlev) = .5*(tke(:,:,nlevp1)+tke(:,:,nlev))
    call wafone (tket, dtstep)
    do jklev = 3, nlev-1
    tke(:,:,jklev) = .5625*(tket(:,:,jklev)+tket(:,:,jklev-1))-.0625*(tket(:,:,jklev+1)+tket(:,:,jklev-2))
    enddo
    tke(:,:,nlev) = .5*(tket(:,:,nlev)+tket(:,:,nlev-1))
    tke(:,:,   2) = .5*(tket(:,:,2   )+tket(:,:,1     ))
#endif

    call wafuv (dtstep)

!***********************************************************************
!  'Physics'
!***********************************************************************

    do jklev = 1, ntop  ! max. relative humidity above and at ntop is set to .2
    do jlat = 1, nlat
    do jlon = 1, nlon
    zzpp = pzer*sigint(jklev) - (pzer-ps(jlon,jlat))*sigalf(jklev)
    zt0t = tzer/t(jlon,jlat,jklev)
    if (t(jlon,jlat,jklev).ge.tzer) then
    zesk = ezer*exp(-ccw1*log(zt0t)+ccw2*(1.-zt0t))  !  partial pressure over water
    else
    zesk = ezer*exp(-cci1*log(zt0t)+cci2*(1.-zt0t))  !  partial pressure over ice
    endif
    zqsat = zesk*eps/(zzpp+zesk*(eps-1.))
    q(jlon,jlat,jklev) = min (q(jlon,jlat,jklev), .2*zqsat)
    q(jlon,jlat,jklev) = max (q(jlon,jlat,jklev), 1.e-12)
    enddo
    enddo
    enddo

    qcw  (:,:,1:ntop)=0.
    qci  (:,:,1:ntop)=0.
    qrain(:,:,1:ntop)=0.
    qsnow(:,:,1:ntop)=0.
    qhail(:,:,1:ntop)=0.

    tke = max(tke, tkemin)

!  definition of temperature and geopot. at half-levels

    t = tvirt/(1. +ep*q -qcw -qci -qrain -qsnow)

    do jklev = nlev, 2, -1
    do jlat = 1, nlat
    do jlon = 1, nlon
    phih(jlon,jlat,jklev)=.5*(phi(jlon,jlat,jklev)+phi(jlon,jlat,jklev-1))
    enddo
    enddo
    enddo

!---------------------------------
!  Update albedo and emissivity
!---------------------------------

!!!    if (mod(jstep,6*ntsrc).eq.0) call surfradpar
    if (mod(jstep,6*ntsrc) == 0) then
      call surfradpar (fmask, fice, qg_surf, vegeta, qsoil_max, qsoil_min, &
 albedo_soil_dry, albedo_soil_wet, emiss1_soil_dry, emiss1_soil_wet, emiss2_soil_dry, emiss2_soil_wet, &
 fveg, albedo_veg, emiss1_veg, emiss2_veg, &
 nlon, nlat, nlevg, nvt, 1, albedo, emisg1, emisg2)
    endif
 

!---------------------------------
!  Update veget. fraction and LAI
!---------------------------------

#ifdef globo
    if (mod(jstep,8*jsday).eq.0) then
      call collect (alont, gfield)
      call collect (alatt, gfield1)
      if (myid == 0) then
        print *
        call veget_param_time_series(gnlon-2, gnlat, &
 gfield(1:gnlon-2,:), gfield1(1:gnlon-2,:), x0, y0, dlon, dlat, pdr0(5), pdr0(4)+dlat*0.5, 1, ndayr, &
 gfield2(1:gnlon-2,:), lai_max, gfield3(1:gnlon-2,:), valmiss)
        gfield2(gnlon-1:gnlon,:) = gfield2(1:2,:)
        gfield3(gnlon-1:gnlon,:) = gfield3(1:2,:)
        print *
      endif
      call disperse (gfield2, lai)
      call disperse (gfield3, fveg)
    endif
#endif
 

!---------------------------------
!  large scale precipitation
!---------------------------------

    rainls = 0.
    snowls = 0.
    call micro

    totpre  = totpre + rainls + snowls
    snfall  = snfall + snowls
    totpre_zoom  = totpre_zoom + rainls + snowls
    snfall_zoom  = snfall_zoom + snowls
    qprec   = qprec  + rainls + snowls  ! for shf files (different accum. period)
    qsnfall = qsnfall + snowls          ! for shf files (different accum. period)

!---------------------------------
!  radiation and convection
!---------------------------------

! Tapering of surface radiation at the beginning of the forecast, due to the
! initial cloud deficit

    rrf = 0.33*exp(-(iday*1440. + ihou*60. + imin)/110.) ! Radiation reduction factor

    if (jstep == nstep0.or.mod(jstep-1,ntsrc) == 0) then

    if (nradm > 0) then ! radiation

! definition of aerosol, ozone, solar time reqtim and solar declin. rdecli,
! used in radint for Geleyn radiation - the last two quantites vary slowly so can be kept
! constant for a few days, as in the case radintec is not called)

    if (jstep == nstep0.or.mod(jstep-1, ntsrc*400).eq.0) then
#ifdef rad_ecmwf
      call aerdef (nyrc, nmonc, ndayc, nhouc, nminc, infx, infy)
#else
      call radiat_init(rdecli0, reqtim0, ozon, aerosol, aerotot, &
   nyrc, nmonc, ndayc, nhouc, nminc, &
   nlon, nlat, nlev, ntaer, sig, sigint, sigalf, pzer, dlat, dlon, &
   alatt, alont, ps, tskin, t, phig(:,:)/g, infy, infx, myid)
#endif
    endif

! time lapse of ntsrc*dtstep sec. to compute tendencies of surface radiation fluxes (useful only for solar rad.)

    call calendar (nyrc , nmonc , ndayc , nhouc , nminc , 0, 0, nint(ntsrc*dtstep/60.),  &
                   nyrc1, nmonc1, ndayc1, nhouc1, nminc1, ndayr1)

    if (myid == 0.and.jstep == nstep0) then
    if(nradm == 1) then
    print*, 'Geleyn radiation scheme'
    elseif(nradm == 2) then

#ifdef oldrad
    print*, 'Old (2004) ECMWF radiation'
#else
    write(*,'(a,i2)') ' ECMWF radiation (2012) with McICA =', mcica
#endif

    else
    print*, 'No radiation scheme applied'
    endif
    print*
    endif

    call cloudfr                                          ! cloud fraction
    call ccloud                                           ! total cloud cover

    if (nradm.ne.0) call radint (ndayr1, nhouc1, nminc1, infx, rrf) ! Geleyn radiation

#ifdef rad_ecmwf
! ------------ ECMWF radiation scheme --------------

#ifdef globo
    if (nradm == 2 .and. (jstep == nstep0.or.mod(jstep-1,ntsrc*3) == 0)) then       ! correction with ECMWF radiation
#else
    if (nradm == 2 .and. (jstep == nstep0.or.mod(jstep-1,ntsrc*5) == 0)) then       ! correction with ECMWF radiation
#endif

    corrdt = 0.
    corvis = 0.
    corirr = 0.

    do 30 jlat = 2, nlat-2, 2
    call radintec (jlat, 2, nyrc1, nmonc1, ndayc1, nhouc1, nminc1, infx, infy, rrf)  ! ECMWF radiation

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
      if (nprocsx.eq.1) then

#ifdef globo
        corrdt(nlon,jlat,:) = corrdt(2,jlat,:)
        corvis(nlon,jlat) = corvis(2,jlat)
        corirr(nlon,jlat) = corirr(2,jlat)
#endif

      else
        call u_ghost (corrdt(2,jlat,:), ip_w, corrdt(nlon,jlat,:), ip_e, nlev)
        call u_ghost (corvis(2,jlat), ip_w, corvis(nlon,jlat), ip_e, 1)
        call u_ghost (corirr(2,jlat), ip_w, corirr(nlon,jlat), ip_e, 1)
      endif
#else
#ifdef globo
      corrdt(nlon,jlat,:) = corrdt(2,jlat,:)
      corvis(nlon,jlat) = corvis(2,jlat)
      corirr(nlon,jlat) = corirr(2,jlat)
#endif
#endif

    do jklev = 1, nlev
    do jlon = 3, nlonm1, 2
    corrdt(jlon,jlat,jklev) = 0.5*(corrdt(jlon-1,jlat,jklev)+corrdt(jlon+1,jlat,jklev))
    enddo
    enddo
    do jlon = 3, nlonm1, 2
    corvis(jlon,jlat) = 0.5*(corvis(jlon-1,jlat)+corvis(jlon+1,jlat))
    corirr(jlon,jlat) = 0.5*(corirr(jlon-1,jlat)+corirr(jlon+1,jlat))
    enddo
 30 continue

!  latitudinal interpolation

#ifdef mpi
      call u_ghost (corrdt(2:nlonm1,2,:), ip_s, corrdt(2:nlonm1,nlat,:), ip_n, (nlon-2)*nlev)
      call u_ghost (corvis(2:nlonm1,2), ip_s, corvis(2:nlonm1,nlat), ip_n, nlon-2)
      call u_ghost (corirr(2:nlonm1,2), ip_s, corirr(2:nlonm1,nlat), ip_n, nlon-2)
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
    enddo
    enddo

#ifdef globo
#ifdef mpi
        if (ip_s.eq.ip_null) then
          do jlon = 1, nlon
            corvis(jlon,1) = corvis(2,2)
            corirr(jlon,1) = corirr(2,2)
            corrdt(jlon,1,:) = corrdt(2,2,:)
          enddo
        endif
#endif
#endif

    call filt2t (corrdt, 1.)   !  Smoothing of RADINTEC correction
    call filt2t1(corvis, 1.)   !  Smoothing of RADINTEC correction
    call filt2t1(corirr, 1.)   !  Smoothing of RADINTEC correction

    endif ! condition on ECMWF radiation
! ------------ End of ECMWF radiation scheme --------------
#endif

    endif ! end of radiation

    raicon = 0.   ! output quantities are set to zero
    snocon = 0.
    dtdt   = 0.
    dqdt   = 0.
    dqcwdt = 0.
    dqcidt = 0.
    dqrndt = 0.
    dqsndt = 0.

    if (nlcadj) then ! convection

    call convection_kf (ntop)
    call filt2t (dtdt  , 1.)   !  Smoothing of convection
    call filt2t (dqdt  , 1.)   !  Smoothing of convection
    call filt2t (dqcwdt, 1.)   !  Smoothing of convection
    call filt2t (dqcidt, 1.)   !  Smoothing of convection
    call filt2t1(raicon, 1.)   !  Smoothing of convection
    call filt2t1(snocon, 1.)   !  Smoothing of convection

! Reduction of convective precip. near boundaries

#ifndef globo

#ifdef mpi
      if (ip_s.eq.ip_null) then  ! south boundary
#endif
    do jbl = 2, nbl
      zrid = 1.-sqrt(bndrel(jbl))
      do jlon = 2, nlonm1
        raicon(jlon,jbl) = raicon(jlon,jbl)*zrid
        snocon(jlon,jbl) = snocon(jlon,jbl)*zrid
      enddo
    enddo
#ifdef mpi
      endif
#endif

#ifdef mpi
      if (ip_n.eq.ip_null) then  ! north boundary
#endif
    do jbl = 1, nbl-1
      zrid = 1.-sqrt(bndrel(jbl+1))
      do jlon = 2, nlonm1
        raicon(jlon,nlat-jbl) = raicon(jlon,nlat-jbl)*zrid
        snocon(jlon,nlat-jbl) = snocon(jlon,nlat-jbl)*zrid
      enddo
    enddo
#ifdef mpi
      endif
#endif

#ifdef mpi
      if (ip_w.eq.ip_null) then   ! west boundary
#endif
    do jlat = 2, nlatm1
      do jbl = 2, nbl
        zrid = 1.-sqrt(bndrel(jbl))
        raicon(jbl,jlat) = raicon(jbl,jlat)*zrid
        snocon(jbl,jlat) = snocon(jbl,jlat)*zrid
      enddo
    enddo
#ifdef mpi
      endif
#endif

#ifdef mpi
      if (ip_e.eq.ip_null) then   ! east boundary
#endif
    do jlat = 2, nlatm1
      do jbl = 1, nbl-1
        zrid = 1.-sqrt(bndrel(jbl+1))
        raicon(nlon-jbl,jlat) = raicon(nlon-jbl,jlat)*zrid
        snocon(nlon-jbl,jlat) = snocon(nlon-jbl,jlat)*zrid
      enddo
    enddo
#ifdef mpi
      endif
#endif

#endif

    endif ! end of convection

    totpre  = totpre  + raicon + snocon
    conpre  = conpre  + raicon + snocon
    snfall  = snfall  + snocon
    totpre_zoom  = totpre_zoom  + raicon + snocon
    conpre_zoom  = conpre_zoom  + raicon + snocon
    snfall_zoom  = snfall_zoom  + snocon
    qprec   = qprec   + raicon + snocon
    qprecc  = qprecc  + raicon + snocon
    qsnfall = qsnfall + snocon

    dtdt  = dtdt + geldt + corrdt

    if (jstep == 1) frvis = max(0., gelvis + corvis)
    dfrvis = (max(0., gelvis + corvis) - frvis)/float(ntsrc)
    if (jstep == 1) frirr = gelirr + corirr
    dfrirr = (gelirr + corirr - frirr)/float(ntsrc)

    endif   ! condition on ntsrc - end of radiation and convection

!  update for tendencies due to radiation and convection

    t = t + dtdt*dtstep
    q     = max (1.e-12, q + dqdt  *dtstep)
    qcw   = max (0., qcw   + dqcwdt*dtstep)
    qci   = max (0., qci   + dqcidt*dtstep)
    qrain = max (0., qrain + dqrndt*dtstep)
    qsnow = max (0., qsnow + dqsndt*dtstep)
    frvis = max (0., frvis + dfrvis)
    frirr = frirr + dfrirr

!---------------------------------
!  soil and vertical diffusion
!---------------------------------

    call tofd

#ifdef mpi
      call u_ghost (v(:,2,:), ip_s, v(:,nlat,:), ip_n, nlon*nlev)
      call u_ghost (u(nlonm1,:,:), ip_e, u(1,:,:), ip_w, nlat*nlev)
#endif

!!!    call vdiff_old (dtstep)
    call surface_layer (dtstep, jstep)

    if (flag_surf_flux == 0 ) then
! Specific and latent heat fluxes and water vapour fluxes defined by pochva
! scheme
      do jlat=1,nlat
      do jlon=1,nlon
        if (mask_soil(jlon,jlat) == 1) then
          hflux(jlon,jlat)=-flh_specif(jlon,jlat)
          qflux(jlon,jlat)=-fl_wv(jlon,jlat)
        endif
      enddo
      enddo
    endif

    call vdiff (dtstep)

!    flag_call_pochva=1
!#ifdef globo
!    if (mod(jstep-1, 2) == 0) then
!      flag_call_pochva=1
!    else
!      flag_call_pochva=0
!    endif
!#endif

 if (flag_call_pochva == 1) then

    call sea_surface

!do jlat=2,nlat-1
!do jlon=2,nlon-1
!  if (abs(tskin(jlon,jlat))>1000.) print *,"*4 oshibka * tskin ",tskin(jlon,jlat),jlon,jlat,myid
!  if (abs(qskin(jlon,jlat))>10.) print *,"*4 oshibka * qskin ",qskin(jlon,jlat),jlon,jlat,myid
!  if (abs(snow(jlon,jlat))>1000.) print *,"*4 oshibka * snow ",snow(jlon,jlat),jlon,jlat,myid
!!  if (abs(fsnow(jlon,jlat))>10.) print *,"*4 oshibka * fsnow ",fsnow(jlon,jlat),jlon,jlat,myid
!  if (isnan(tskin(jlon,jlat))) print *,"*4 oshibka * tskin ",jstep,tskin(jlon,jlat),jlon,jlat,myid
!  if (isnan(qskin(jlon,jlat))) print *,"*4 oshibka * qskin ",jstep,qskin(jlon,jlat),jlon,jlat,myid
!  if (isnan(snow(jlon,jlat))) print *,"*4 oshibka * snow ",jstep,snow(jlon,jlat),jlon,jlat,myid
!!  if (isnan(fsnow(jlon,jlat))) print *,"*4 oshibka * fsnow ",jstep,fsnow(jlon,jlat),jlon,jlat,myid
!  if (isnan(u(jlon,jlat,nlev))) print *,"*4 oshibka * u ",jstep,u(jlon,jlat,nlev),jlon,jlat,myid
!  if (isnan(v(jlon,jlat,nlev))) print *,"*4 oshibka * v ",jstep,v(jlon,jlat,nlev),jlon,jlat,myid
!  if (isnan(t(jlon,jlat,nlev))) print *,"*4 oshibka * t ",jstep,t(jlon,jlat,nlev),jlon,jlat,myid
!  if (isnan(q(jlon,jlat,nlev))) print *,"*4 oshibka * q ",jstep,q(jlon,jlat,nlev),jlon,jlat,myid
!  if (isnan(tke(jlon,jlat,nlev))) print *,"*4 oshibka * tke ",jstep,tke(jlon,jlat,nlev),jlon,jlat,myid
!  if (isnan(phi(jlon,jlat,nlev))) print *,"*4 oshibka * phi ",jstep,phi(jlon,jlat,nlev),jlon,jlat,myid
!  if (isnan(ps(jlon,jlat))) print *,"*4 oshibka * ps ",jstep,ps(jlon,jlat),jlon,jlat,myid
!  if (isnan(rainls(jlon,jlat))) print *,"*4 oshibka * rainls ",jstep,rainls(jlon,jlat),jlon,jlat,myid
!  if (isnan(raicon(jlon,jlat))) print *,"*4 oshibka * raicon ",jstep,raicon(jlon,jlat),jlon,jlat,myid
!  if (isnan(snowls(jlon,jlat))) print *,"*4 oshibka * snowls ",jstep,snowls(jlon,jlat),jlon,jlat,myid
!  if (isnan(snocon(jlon,jlat))) print *,"*4 oshibka snocon*  ",jstep,snocon(jlon,jlat),jlon,jlat,myid
!  if (isnan(frvis(jlon,jlat))) print *,"*4 oshibka * frvis ",jstep,frvis(jlon,jlat),jlon,jlat,myid
!  if (isnan(frirr(jlon,jlat))) print *,"*4 oshibka * frirr ",jstep,frirr(jlon,jlat),jlon,jlat,myid
!  if (isnan(kturb_surf_h(jlon,jlat))) print *,"*4 oshibka * kturb_surf_h ",jstep,kturb_surf_h(jlon,jlat),jlon,jlat,myid
!  if (isnan(kturb_surf_q(jlon,jlat))) print *,"*4 oshibka * kturb_surf_q ",jstep,kturb_surf_q(jlon,jlat),jlon,jlat,myid
!enddo
!enddo

! pochva --->

    h_lev_bottom(:,:) = (phi(:,:,nlev)-phig(:,:))/g
    p_lev_bottom(:,:) = pzer*sigint(nlev) - (pzer-ps(:,:))*sigalf(nlev)
    fl_rain(:,:) = (rainls+raicon/float(ntsrc))/dtstep
    fl_snow(:,:) = (snowls+snocon/float(ntsrc))/dtstep
! City Heat Island, urban vegetation index is 20 
    fl_rad_tot(:,:) = frvis(:,:)+frirr(:,:)*(1.-0.1*vegeta(:,:,21))
    snow(:,:) = snow(:,:)*1.e3 ! m of equiv. water ---> kg/m/m

    if (flag_surf_flux == 1 ) then
!  For option of using heat fluxes as input fileds for Pochva scheme
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
                 mask_soil, h_lev_bottom, ps, p_lev_bottom, t(:,:,nlev), q(:,:,nlev), &
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

! Snow albedo on glacier and on sea ice
    do jlat=2,nlatm1
    do jlon=2,nlonm1
      if (int(psi_soil(jlon,jlat,1)*10.) == 22) then ! glacier
        if (fmask(jlon,jlat) < 0.5 ) then ! glacier on land
          alsn(jlon,jlat)=0.83
        else ! sea ice
          alsn(jlon,jlat)=0.70
        endif
      endif
    enddo
    enddo
! <---


! <---

!do jlat=2,nlat-1
!do jlon=2,nlon-1
!  if (abs(tskin(jlon,jlat))>373..or.abs(tskin(jlon,jlat))<173.) print *,"*5 oshibka * tskin ",jstep,tskin(jlon,jlat),jlon,jlat,myid
!  if (abs(qskin(jlon,jlat))>10.) print *,"*5 oshibka * qskin ",jstep,qskin(jlon,jlat),jlon,jlat,myid
!  if (abs(snow(jlon,jlat))>1000.) print *,"*5 oshibka * snow ",jstep,snow(jlon,jlat),jlon,jlat,myid
!  if (abs(fsnow(jlon,jlat))>1.) print *,"*5 oshibka * fsnow ",jstep,fsnow(jlon,jlat),jlon,jlat,myid
!  if (alsn(jlon,jlat)<0.3.or.alsn(jlon,jlat)>1.) print *,"*5 oshibka * albedo ",jstep,fsnow(jlon,jlat),jlon,jlat,myid
!  if (isnan(tskin(jlon,jlat))) print *,"*5* oshibka tskin ",jstep,tskin(jlon,jlat),jlon,jlat,myid
!  if (isnan(qskin(jlon,jlat))) print *,"*5* oshibka qskin ",jstep,qskin(jlon,jlat),jlon,jlat,myid
!  if (isnan(snow(jlon,jlat))) print *,"*5 oshibka * snow ",jstep,snow(jlon,jlat),jlon,jlat,myid
!  if (isnan(fsnow(jlon,jlat))) print *,"*5 oshibka * fsnow ",jstep,fsnow(jlon,jlat),jlon,jlat,myid
!  if (isnan(v(jlon,jlat,nlev))) print *,"*5 oshibka * v ",jstep,v(jlon,jlat,nlev),jlon,jlat,myid
!  if (isnan(t(jlon,jlat,nlev))) print *,"*5 oshibka * t ",jstep,t(jlon,jlat,nlev),jlon,jlat,myid
!  if (isnan(q(jlon,jlat,nlev))) print *,"*5 oshibka * q ",jstep,q(jlon,jlat,nlev),jlon,jlat,myid
!  if (isnan(tke(jlon,jlat,nlev))) print *,"*5 oshibka * tke ",jstep,tke(jlon,jlat,nlev),jlon,jlat,myid
!  if (isnan(u(jlon,jlat,nlev))) then
!#ifdef mpi
!      call mpi_abort(comm, ierr)
!#endif
!   stop
!  endif
!enddo
!enddo

do k=1,n_outpoint
  if (myid == iproc_outpoint(k)) then
    i=i_outpoint(k)
    j=j_outpoint(k)
    write (file_outpoint,'(a,i2.2,a)') "surf_model_",K,".dat"
    open (31, file=file_outpoint, form="formatted", position="append")
    write (31,'(100e16.8)') float(jstep)*dtstep/3600., &
 ps(i,j), &                !  2
 p_lev_bottom(i,j), &      !  3
 t(i,j,nlev)-273.15, &     !  4
 q(i,j,nlev), &            !  5
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
 emisg2(i,j),&             ! 18
 fice(i,j),&               ! 19
 float(mask_soil(i,j)),&   ! 20
 sqrt(u(i,j,nlev)**2+v(i,j,nlev)**2)  ! 21
    close (31)
    write (file_outpoint,'(a,i2.2,a)') "tsoil_model_",K,".dat"
    open (31, file=file_outpoint, form="formatted", position="append")
    write (31,'(100e16.8)') float(jstep)*dtstep/3600., &
 tskin(i,j)-273.15, &      !  2
 tg_surf(i,j)-273.15, &    !  3
 tg(i,j,1:nlevg)-273.15    !  4:12
    close (31)
    write (file_outpoint,'(a,i2.2,a)') "qsoil_model_",K,".dat"
    open (31, file=file_outpoint, form="formatted", position="append")
    write (31,'(100e16.8)') float(jstep)*dtstep/3600., &
 qskin(i,j), &        !  2
 qg_surf(i,j), &           !  3
 qg(i,j,1:nlevg)           !  4:12
    close (31)
    write (file_outpoint,'(a,i2.2,a)') "qsoil_rel_model_",K,".dat"
    open (31, file=file_outpoint, form="formatted", position="append")
    write (31,'(100e16.8)') float(jstep)*dtstep/3600., &
 qskin(i,j), &        !  2 
 (qg_surf(i,j)-qsoil_min(i,j,1))/(qsoil_max(i,j,1)-qsoil_min(i,j,1)), & !  3 
 (qg(i,j,1:nlevg)-qsoil_min(i,j,1:nlevg))/(qsoil_max(i,j,1:nlevg)-qsoil_min(i,j,1:nlevg)) !  4:12 
    close (31)
    write (file_outpoint,'(a,i2.2,a)') "snow_lev_model_",K,".dat"
    open (31, file=file_outpoint, form="formatted", position="append")
    write (31,'(100e16.8)') float(jstep)*dtstep/3600., &
 snow(i,j), &             !  2
 fsnow(i,j), &            !  3
 alsn(i,j), &             !  4
 snow_lev(i,j,1:nlev_snow) !  5:17
    close (31)
    write (file_outpoint,'(a,i2.2,a)') "fice_soil_model_",K,".dat"
    open (31, file=file_outpoint, form="formatted", position="append")
    write (31,'(100e16.8)') float(jstep)*dtstep/3600., &
 fice_soil_surf(i,j), &            !  2
 fice_soil(i,j,1:nlevg)            !  3:11
    close (31)
    write (file_outpoint,'(a,i2.2,a)') "snow_t_model_",K,".dat"
    open (31, file=file_outpoint, form="formatted", position="append")
    write (31,'(100e16.8)') float(jstep)*dtstep/3600., &
 snow_t(i,j,1:nlev_snow)-273.15    !  2:14
    close (31)
    write (file_outpoint,'(a,i2.2,a)') "snow_dens_model_",K,".dat"
    open (31, file=file_outpoint, form="formatted", position="append")
    write (31,'(100e16.8)') float(jstep)*dtstep/3600., &
 snow_dens(i,j,1:nlev_snow)        !  2:14
    close (31)
  endif
enddo

 endif ! flag_call_pochva == 1

!  cumulated fluxes and 2m min/max temperature
!  incoming solar radiation is computed using albedo with snow (changes with tskin as in radint and radintec)

    do jlat = 2, nlatm1
    do jlon = 2, nlonm1
!    if (tskin(jlon,jlat).lt.277.) then
!    zalsn = alsn-.1 + .2*(277.-tskin(jlon,jlat))/14.
!    zalsn = min (zalsn, alsn+.1)
!    else
!    zalsn = alsn-.1
!    endif
!    zalsn = albedo(jlon,jlat)*(1.-fsnow(jlon,jlat)) + zalsn*fsnow(jlon,jlat)
    zalsn = albedo(jlon,jlat)*(1.-fsnow(jlon,jlat)) + alsn(jlon,jlat)*fsnow(jlon,jlat)
    solar (jlon,jlat) = solar (jlon,jlat) + frvis(jlon,jlat)/(1.-zalsn)
    cswfl (jlon,jlat) = cswfl (jlon,jlat) + frvis(jlon,jlat)*dtstep
    clwfl (jlon,jlat) = clwfl (jlon,jlat) + frirr(jlon,jlat)*dtstep
    cswfl_zoom (jlon,jlat) = cswfl_zoom (jlon,jlat) + frvis(jlon,jlat)*dtstep
    clwfl_zoom (jlon,jlat) = clwfl_zoom (jlon,jlat) + frirr(jlon,jlat)*dtstep
    chflux(jlon,jlat) = chflux(jlon,jlat) - hflux(jlon,jlat)*dtstep
    shf_accum(jlon,jlat) = shf_accum(jlon,jlat) - hflux(jlon,jlat)
    zhea = ylwv + (cpv-cw)*(tskin(jlon,jlat)-tzer)
    cqflux(jlon,jlat) = cqflux(jlon,jlat) - qflux(jlon,jlat)*dtstep*zhea
    cwvflux(jlon,jlat) = cwvflux(jlon,jlat) - qflux(jlon,jlat)*dtstep
    chflux_zoom(jlon,jlat) = chflux_zoom(jlon,jlat) - hflux(jlon,jlat)*dtstep
    cqflux_zoom(jlon,jlat) = cqflux_zoom(jlon,jlat) - qflux(jlon,jlat)*dtstep*zhea
    cwvflux_zoom(jlon,jlat) = cwvflux_zoom(jlon,jlat) - qflux(jlon,jlat)*dtstep
    cfl_heat_soil_bottom(jlon,jlat) = cfl_heat_soil_bottom(jlon,jlat) + fl_heat_soil_bottom(jlon,jlat)*dtstep
    cfl_water_soil_bottom(jlon,jlat) = cfl_water_soil_bottom(jlon,jlat) + fl_water_soil_bottom(jlon,jlat)*dtstep
    runoff(jlon,jlat) = runoff(jlon,jlat) + fl_runoff(jlon,jlat)*dtstep
    runoff_tot(jlon,jlat) = runoff_tot(jlon,jlat) + fl_runoff_tot(jlon,jlat)*dtstep
    cfl_heat_soil_bottom_zoom(jlon,jlat) = cfl_heat_soil_bottom_zoom(jlon,jlat) + fl_heat_soil_bottom(jlon,jlat)*dtstep
    cfl_water_soil_bottom_zoom(jlon,jlat) = cfl_water_soil_bottom_zoom(jlon,jlat) + fl_water_soil_bottom(jlon,jlat)*dtstep
    runoff_zoom(jlon,jlat) = runoff_zoom(jlon,jlat) + fl_runoff(jlon,jlat)*dtstep
    runoff_tot_zoom(jlon,jlat) = runoff_tot_zoom(jlon,jlat) + fl_runoff_tot(jlon,jlat)*dtstep
    lhf_accum(jlon,jlat) = lhf_accum(jlon,jlat) - qflux(jlon,jlat)*zhea
    qf_accum(jlon,jlat) = qf_accum(jlon,jlat) - qflux(jlon,jlat)*dtstep
    t2min(jlon,jlat) = min (t2min(jlon,jlat), t2(jlon,jlat))
    t2max(jlon,jlat) = max (t2max(jlon,jlat), t2(jlon,jlat))
    t2min_zoom(jlon,jlat) = min (t2min_zoom(jlon,jlat), t2(jlon,jlat))
    t2max_zoom(jlon,jlat) = max (t2max_zoom(jlon,jlat), t2(jlon,jlat))
    enddo
    enddo

!---------------------------------
!  orographic gravity wave drag and blocking drag
!---------------------------------

    call blockdrag (jstep)
    call orogdrag (jstep)

!---------------------------------
!  filter
!---------------------------------

#ifdef globo

    if (jstep == nstep0.or.mod(jstep,5*ntsrc) == 0) then  ! Polar filter of full fields

    call polaveruv
    call polavert (t)
    call polavert (tke(:,:,2:nlevp1))
    tke = max(tke, tkemin)
    call polavert (q)
    q = max (q, 1.e-12)
    call polavert (qcw)
    qcw = max (qcw, 0.)
    call polavert (qci)
    qci = max (qci, 0.)
    call polavert (qrain)
    qsnow = max (qrain, 0.)
    call polavert (qsnow)
    qsnow = max (qsnow, 0.)
    call polavert (qhail)
    qsnow = max (qhail, 0.)

    else     ! polar filter of tendencies

    u = u-uf
    v = v-vf
    call polaveruv
    u = u+uf
    v = v+vf
    t = t-tf
    call polavert (t)
    t = t+tf
    q = q-qf
    call polavert (q)
    q = max (q+qf, 1.e-12)
    qcw = qcw - qcwf
    call polavert (qcw)
    qcw = max (qcw+qcwf, 0.)
    qci = qci - qcif
    call polavert (qci)
    qci = max (qci+qcif, 0.)

    endif

#endif

    if (mod(jstep-1,4) == 0) then
    call filt2uv(anu2v)
    call filt2t (t, anu2)
    endif

    if (mod(jstep-1,8).eq.0) then
    call filt2t (tke(:,:,2:nlevp1), anu2)
    call filt2t (q    , anu2)
    call filt2t (qcw  , anu2)
    call filt2t (qci  , anu2)
    call filt2t (qrain, anu2)
    call filt2t (qsnow, anu2)
    call filt2t (qhail, anu2)
    endif

! <augh,max: tue sep 30 09:15:42 cest 2008>
#ifdef chem
    call bol2chem (jstep, myid)
#endif
! </augh,max>

#ifndef globo

!***********************************************************************
! lateral boundary conditions
!***********************************************************************

      if(nlbfix) then

      ztime = 0.

      else ! .not.nlbfix

!  defines the boundary conditions at time jstep*dtstep by linear interpolation in time.

      ztime = float(mod(jstep,ntsbou))/float(ntsbou)
      if (jstep >= (nbc-1)*ntsbou) ztime = 1.

      if (mod(jstep,ntsbou).eq.0.and.jstep.lt.(nbc-1)*ntsbou) then    !  read new boundary file
      psb1  = psb2
      ub1   = ub2
      vb1   = vb2
      tb1   = tb2
      qb1   = qb2
      qcwb1 = qcwb2
      qcib1 = qcib2

      interv = jstep/ntsbou + 1

      call rdmhf_atm (interv+1, 0)
      call defclt (qcwb2, qcib2, tb2, nlon, nlat, nlev, tzer, 2)

      call rdmhf_soil (interv, 0)

     if (nlclimate) then
!!      tg(:,:,nlevg) = tgclim(:,:)  ! uncomment here for long runs (with variable soil T)
!!      qg(:,:,nlevg) = qgclim(:,:)  ! uncomment here for long runs (with variable soil T)
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
              snow(jlon,jlat) = max (50., snow(jlon,jlat)) ! snow (kg/m/m)
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

!  linear interpolation in time between b1 and b2

#ifdef mpi
        if (ip_w.eq.ip_null) then   ! west boundary
#endif

        do jlat = 1, nlat
          jlo2 = nbl
          if (ip_s.eq.ip_null .and. jlat.lt.nbl     ) jlo2 = min (jlat,        nbl)
          if (ip_n.eq.ip_null .and. jlat.gt.nlat-nbl) jlo2 = min (nlat-jlat+1, nbl)
        do jlon = 1, jlo2

! error trapping in case of invalid value of ps read in boundary files within the BOLAM frame

         if(psb2(jlon,jlat).lt.1.e4.and.myid.eq.0.and.(.not.nlbfix)) then
          print*
          print*, "jlon, jlat, psb2:", jlon, jlat, psb2(jlon,jlat)
          print*, "Invalid value of ps read in lateral b.c.: stop!"
          print*, "If ps=-9999. forecast b.c. files are defined on frames"
          print*, "which do not completely cover the NBL lateral frames."
          print*
#ifdef mpi
            call mpi_abort(comm, ierr)
#endif
          stop
         endif
        zpsb = psb1(jlon,jlat) + ztime*(psb2(jlon,jlat)-psb1(jlon,jlat))
        z1 = bndrel(jlon)
        z2 = 1.-z1
        ps(jlon,jlat) = z1*zpsb + z2*ps(jlon,jlat)
        enddo
        enddo
      do jklev = 1, nlev
      do jlat = 1, nlat
          jlo2 = nbl
          if (ip_s.eq.ip_null .and. jlat.lt.nbl     ) jlo2 = min (jlat,        nbl)
          if (ip_n.eq.ip_null .and. jlat.gt.nlat-nbl) jlo2 = min (nlat-jlat+1, nbl)
      do jlon = 1, jlo2
      zub   = ub1  (jlon,jlat,jklev) + ztime*(ub2  (jlon,jlat,jklev)-ub1  (jlon,jlat,jklev))
      zvb   = vb1  (jlon,jlat,jklev) + ztime*(vb2  (jlon,jlat,jklev)-vb1  (jlon,jlat,jklev))
      ztb   = tb1  (jlon,jlat,jklev) + ztime*(tb2  (jlon,jlat,jklev)-tb1  (jlon,jlat,jklev))
      zqb   = qb1  (jlon,jlat,jklev) + ztime*(qb2  (jlon,jlat,jklev)-qb1  (jlon,jlat,jklev))
      zqcwb = qcwb1(jlon,jlat,jklev) + ztime*(qcwb2(jlon,jlat,jklev)-qcwb1(jlon,jlat,jklev))
      zqcib = qcib1(jlon,jlat,jklev) + ztime*(qcib2(jlon,jlat,jklev)-qcib1(jlon,jlat,jklev))
      z1 = bndrel(jlon)
      z2 = 1.-z1
      u  (jlon,jlat,jklev) = z1*zub   + z2*u  (jlon,jlat,jklev)
      v  (jlon,jlat,jklev) = z1*zvb   + z2*v  (jlon,jlat,jklev)
      t  (jlon,jlat,jklev) = z1*ztb   + z2*t  (jlon,jlat,jklev)
      q  (jlon,jlat,jklev) = z1*zqb   + z2*q  (jlon,jlat,jklev)
      qcw(jlon,jlat,jklev) = z1*zqcwb + z2*qcw(jlon,jlat,jklev)
      qci(jlon,jlat,jklev) = z1*zqcib + z2*qci(jlon,jlat,jklev)
      enddo
      enddo
      enddo

#ifdef mpi
        endif
#endif

#ifdef mpi
        if (ip_e.eq.ip_null) then   ! east boundary
#endif

        do jlat = 1, nlat
          jlo2 = 1
          if (ip_s.eq.ip_null .and. jlat.lt.nbl) jlo2 = max(nbl-jlat+1   ,1)
          if (ip_n.eq.ip_null .and. jlat.gt.nbl) jlo2 = max(nbl+jlat-nlat,1)
        do jlon = jlo2, nbl
        jlo1 = nlon-nbl+jlon
        zpsb = psb1(jlo1,jlat) + ztime*(psb2(jlo1,jlat)-psb1(jlo1,jlat))
        z1 = bndrel(nbl-jlon+1)
        z2 = 1.-z1
        ps(jlo1,jlat) = z1*zpsb + z2*ps(jlo1,jlat)
        enddo
        enddo
      do jklev = 1, nlev
      do jlat = 1, nlat
          jlo2 = 1
          if (ip_s.eq.ip_null .and. jlat.lt.nbl) jlo2 = max(nbl-jlat+1   ,1)
          if (ip_n.eq.ip_null .and. jlat.gt.nbl) jlo2 = max(nbl+jlat-nlat,1)
      do jlon = jlo2, nbl
      jlo1 = nlon-nbl+jlon
      zub   = ub1  (jlo1,jlat,jklev) + ztime*(ub2  (jlo1,jlat,jklev)-ub1  (jlo1,jlat,jklev))
      zvb   = vb1  (jlo1,jlat,jklev) + ztime*(vb2  (jlo1,jlat,jklev)-vb1  (jlo1,jlat,jklev))
      ztb   = tb1  (jlo1,jlat,jklev) + ztime*(tb2  (jlo1,jlat,jklev)-tb1  (jlo1,jlat,jklev))
      zqb   = qb1  (jlo1,jlat,jklev) + ztime*(qb2  (jlo1,jlat,jklev)-qb1  (jlo1,jlat,jklev))
      zqcwb = qcwb1(jlo1,jlat,jklev) + ztime*(qcwb2(jlo1,jlat,jklev)-qcwb1(jlo1,jlat,jklev))
      zqcib = qcib1(jlo1,jlat,jklev) + ztime*(qcib2(jlo1,jlat,jklev)-qcib1(jlo1,jlat,jklev))
      z1 = bndrel(nbl-jlon+1)
      z2 = 1.-z1
      u  (jlo1,jlat,jklev) = z1*zub   + z2*u  (jlo1,jlat,jklev)
      v  (jlo1,jlat,jklev) = z1*zvb   + z2*v  (jlo1,jlat,jklev)
      t  (jlo1,jlat,jklev) = z1*ztb   + z2*t  (jlo1,jlat,jklev)
      q  (jlo1,jlat,jklev) = z1*zqb   + z2*q  (jlo1,jlat,jklev)
      qcw(jlo1,jlat,jklev) = z1*zqcwb + z2*qcw(jlo1,jlat,jklev)
      qci(jlo1,jlat,jklev) = z1*zqcib + z2*qci(jlo1,jlat,jklev)
      enddo
      enddo
      enddo

#ifdef mpi
        endif
#endif

#ifdef mpi
        if (ip_s.eq.ip_null) then  ! south boundary
#endif

        do jlat = 1, nbl
          jlo1 = 1
          if (ip_w.eq.ip_null) jlo1 = jlat+1
          jlo2 = nlon
          if (ip_e.eq.ip_null) jlo2 = nlon-jlat
        z1 = bndrel(jlat)
        z2 = 1.-z1
        do jlon = jlo1, jlo2
        zpsb = psb1(jlon,jlat) + ztime*(psb2(jlon,jlat)-psb1(jlon,jlat))
        ps(jlon,jlat) = z1*zpsb + z2*ps(jlon,jlat)
        enddo
        enddo
      do jklev = 1, nlev
      do jlat = 1, nbl
          jlo1 = 1
          if (ip_w.eq.ip_null) jlo1 = jlat+1
          jlo2 = nlon
          if (ip_e.eq.ip_null) jlo2 = nlon-jlat
      z1 = bndrel(jlat)
      z2 = 1.-z1
      do jlon = jlo1, jlo2
      zub   = ub1  (jlon,jlat,jklev) + ztime*(ub2  (jlon,jlat,jklev)-ub1  (jlon,jlat,jklev))
      zvb   = vb1  (jlon,jlat,jklev) + ztime*(vb2  (jlon,jlat,jklev)-vb1  (jlon,jlat,jklev))
      ztb   = tb1  (jlon,jlat,jklev) + ztime*(tb2  (jlon,jlat,jklev)-tb1  (jlon,jlat,jklev))
      zqb   = qb1  (jlon,jlat,jklev) + ztime*(qb2  (jlon,jlat,jklev)-qb1  (jlon,jlat,jklev))
      zqcwb = qcwb1(jlon,jlat,jklev) + ztime*(qcwb2(jlon,jlat,jklev)-qcwb1(jlon,jlat,jklev))
      zqcib = qcib1(jlon,jlat,jklev) + ztime*(qcib2(jlon,jlat,jklev)-qcib1(jlon,jlat,jklev))
      u  (jlon,jlat,jklev) = z1*zub   + z2*u  (jlon,jlat,jklev)
      v  (jlon,jlat,jklev) = z1*zvb   + z2*v  (jlon,jlat,jklev)
      t  (jlon,jlat,jklev) = z1*ztb   + z2*t  (jlon,jlat,jklev)
      q  (jlon,jlat,jklev) = z1*zqb   + z2*q  (jlon,jlat,jklev)
      qcw(jlon,jlat,jklev) = z1*zqcwb + z2*qcw(jlon,jlat,jklev)
      qci(jlon,jlat,jklev) = z1*zqcib + z2*qci(jlon,jlat,jklev)
      enddo
      enddo
      enddo

#ifdef mpi
        endif
#endif

#ifdef mpi
        if (ip_n.eq.ip_null) then  ! north boundary
#endif
        do jlat = 1, nbl
        jla1 = nlat-nbl+jlat
          jlo1 = 1
          if (ip_w.eq.ip_null) jlo1 = nbl-jlat+2
          jlo2 = nlon
          if (ip_e.eq.ip_null) jlo2 = nlonm1-nbl+jlat
        z1 = bndrel(nbl-jlat+1)
        z2 = 1.-z1
        do jlon = jlo1, jlo2
        zpsb = psb1(jlon,jla1) + ztime*(psb2(jlon,jla1)-psb1(jlon,jla1))
        ps(jlon,jla1) = z1*zpsb + z2*ps(jlon,jla1)
        enddo
        enddo
      do jklev = 1, nlev
      do jlat = 1, nbl
      jla1 = nlat-nbl+jlat
          jlo1 = 1
          if (ip_w.eq.ip_null) jlo1 = nbl-jlat+2
          jlo2 = nlon
          if (ip_e.eq.ip_null) jlo2 = nlonm1-nbl+jlat
      z1 = bndrel(nbl-jlat+1)
      z2 = 1.-z1
      do jlon = jlo1, jlo2
      zub   = ub1  (jlon,jla1,jklev) + ztime*(ub2  (jlon,jla1,jklev)-ub1  (jlon,jla1,jklev))
      zvb   = vb1  (jlon,jla1,jklev) + ztime*(vb2  (jlon,jla1,jklev)-vb1  (jlon,jla1,jklev))
      ztb   = tb1  (jlon,jla1,jklev) + ztime*(tb2  (jlon,jla1,jklev)-tb1  (jlon,jla1,jklev))
      zqb   = qb1  (jlon,jla1,jklev) + ztime*(qb2  (jlon,jla1,jklev)-qb1  (jlon,jla1,jklev))
      zqcwb = qcwb1(jlon,jla1,jklev) + ztime*(qcwb2(jlon,jla1,jklev)-qcwb1(jlon,jla1,jklev))
      zqcib = qcib1(jlon,jla1,jklev) + ztime*(qcib2(jlon,jla1,jklev)-qcib1(jlon,jla1,jklev))
      u  (jlon,jla1,jklev) = z1*zub   + z2*u  (jlon,jla1,jklev)
      v  (jlon,jla1,jklev) = z1*zvb   + z2*v  (jlon,jla1,jklev)
      t  (jlon,jla1,jklev) = z1*ztb   + z2*t  (jlon,jla1,jklev)
      q  (jlon,jla1,jklev) = z1*zqb   + z2*q  (jlon,jla1,jklev)
      qcw(jlon,jla1,jklev) = z1*zqcwb + z2*qcw(jlon,jla1,jklev)
      qci(jlon,jla1,jklev) = z1*zqcib + z2*qci(jlon,jla1,jklev)
      enddo
      enddo
      enddo

#ifdef mpi
        endif
#endif

#endif  ! ifndef globo

!***********************************************************************
!  check point
!***********************************************************************

    nfdr(10)  = iday
    nfdr(11)  = ihou
    nfdr(12)  = imin
    nfdrb(10) = iday
    nfdrb(11) = ihou
    nfdrb(12) = imin

#ifdef globo
!!    if (jstep.eq.jsday*3.5) then ! for a restart..
!!      name_file_mhf="globo"
!!      call wrmhf_atm (name_file_mhf, nfdr, pdr, 1, gnlon, 1, gnlat, 0)
!!      call wrmhf_soil(name_file_mhf, nfdr, pdr, 1, gnlon, 1, gnlat, 0, 0)
!!    endif
#endif

  if (nhist > 0.and.mod(jstep,nhist).eq.0) then

    call cloudfr
    call ccloud

    if (myid == 0) imhf=imhf+1

#ifdef globo
      name_file_mhf="globo"
#else
      name_file_mhf="bolam"
#endif
    call wrmhf_atm (name_file_mhf, nfdr, pdr, 1, gnlon, 1, gnlat, 0)
    call wrmhf_soil(name_file_mhf, nfdr, pdr, 1, gnlon, 1, gnlat, 0, 0)

! reset total precipitation and cumulated fluxes

    totpre = 0.
    conpre = 0.
    snfall = 0.
    cswfl  = 0.
    clwfl  = 0.
    chflux = 0.
    cqflux = 0.
    cwvflux = 0.
    cfl_heat_soil_bottom = 0.
    cfl_water_soil_bottom = 0.
    runoff = 0.
    runoff_tot = 0.
    t2min(2:nlonm1,2:nlatm1) = 999.
    t2max(2:nlonm1,2:nlatm1) = 0.

  endif

#ifdef globo
  if (nhist_zoom > 0.and.mod(jstep,nhist_zoom) == 0.and.jstep <= 7*jsday) then

    call cloudfr
    call ccloud

    if (myid == 0) imhf_zoom=imhf_zoom+1

    imhf_copy   = imhf
    totpre_copy = totpre
    conpre_copy = conpre
    snfall_copy = snfall
    cswfl_copy  = cswfl
    clwfl_copy  = clwfl
    chflux_copy = chflux
    cqflux_copy = cqflux
    cwvflux_copy = cwvflux
    cfl_heat_soil_bottom_copy = cfl_heat_soil_bottom
    cfl_water_soil_bottom_copy = cfl_water_soil_bottom
    runoff_copy = runoff
    runoff_tot_copy = runoff_tot
    t2min_copy = t2min
    t2max_copy = t2max
 
    imhf   = imhf_zoom
    totpre = totpre_zoom
    conpre = conpre_zoom
    snfall = snfall_zoom
    cswfl  = cswfl_zoom
    clwfl  = clwfl_zoom
    chflux = chflux_zoom
    cqflux = cqflux_zoom
    cwvflux = cwvflux_zoom
    cfl_heat_soil_bottom = cfl_heat_soil_bottom_zoom
    cfl_water_soil_bottom = cfl_water_soil_bottom_zoom
    runoff = runoff_zoom
    runoff_tot = runoff_tot_zoom
    t2min = t2min_zoom
    t2max = t2max_zoom

    name_file_mhf="globo_zoom"
    call wrmhf_atm (name_file_mhf, nfdrb, pdrb, i_zoom_lon_ini, i_zoom_lon_fin, j_zoom_lat_ini, j_zoom_lat_fin, flag_joint)
    call wrmhf_soil(name_file_mhf, nfdrb, pdrb, i_zoom_lon_ini, i_zoom_lon_fin, j_zoom_lat_ini, j_zoom_lat_fin, flag_joint, 0)
 
    imhf   = imhf_copy
    totpre = totpre_copy
    conpre = conpre_copy
    snfall = snfall_copy
    cswfl  = cswfl_copy
    clwfl  = clwfl_copy
    chflux = chflux_copy
    cqflux = cqflux_copy
    cwvflux = cwvflux_copy
    cfl_heat_soil_bottom = cfl_heat_soil_bottom_copy
    cfl_water_soil_bottom = cfl_water_soil_bottom_copy
    runoff = runoff_copy
    runoff_tot = runoff_tot_copy
    t2min = t2min_copy
    t2max = t2max_copy

! reset total precipitation and cumulated fluxes

    totpre_zoom = 0.
    conpre_zoom = 0.
    snfall_zoom = 0.
    cswfl_zoom  = 0.
    clwfl_zoom  = 0.
    chflux_zoom = 0.
    cqflux_zoom = 0.
    cwvflux_zoom = 0.
    cfl_heat_soil_bottom_zoom = 0.
    cfl_water_soil_bottom_zoom = 0.
    runoff_zoom = 0.
    runoff_tot_zoom = 0.
    t2min_zoom(2:nlonm1,2:nlatm1) = 999.
    t2max_zoom(2:nlonm1,2:nlatm1) = 0.

  endif
#endif

#ifdef globo
    if (nhist_soil_full > 0.and.mod(jstep,nhist_soil_full) == 0.and.jstep <= 5*jsday) then
      name_file_mhf="globo"
      call wrmhf_soil(name_file_mhf, nfdr, pdr, 1, gnlon, 1, gnlat, 0, 1)
    endif
#endif

#ifdef globo
    if (mod(jstep,ndrunt).eq.0) call runout_g (jstep, nyrc, nmonc, ndayc, nhouc, nminc)
#else
    if (mod(jstep,ndrunt).eq.0) call runout_b (jstep)  ! run time diagnostics
#endif

    if (ntswshf > 0) then
    if (mod(jstep,ntswshf) == 0.or.jstep == 1) then

    call cloudfr
    call ccloud

#ifdef globo
    call wrshf_g(jstep) !  write SHF
#else
    call wrshf_b(jstep) !  write SHF
#endif

    endif
    endif

! Writting of restart file
    if (mod(jstep,nrst) == 0.and.jstep > nstep0) then
      call wrrf(jstep)
      call wrrf_pochva(jstep, myid, gfield, gnlon, gnlat, irf)
    endif

!***********************************************************************
 1000 continue    ! End of timestep
!***********************************************************************

!    call collect (fmask, gfield1) ! examples for 2-D plotting
!    str = 'u(nlev)'
!    call collect (u(1:nlon,1:nlat,nlev), gfield)
!    if (myid.eq.0) call plotout (gfield, gfield1, gnlon, gnlat, str, 99)
!    str = 'gelvis'
!    call collect (gelvis, gfield)
!    if (myid.eq.0) call plotout (gfield, gfield1, gnlon, gnlat, str, 99)

    call system_clock (stime2, countrate)
    if(myid.eq.0) print*,'System time (sec):', float(stime2-stime1)/countrate

! <augh,max: tue sep 30 09:15:42 cest 2008>
#ifdef chem
    call last_chem (myid)
#endif
! </augh,max>

! clean up mpi environment

#ifdef mpi
      call mpi_finalize (ierr)
#endif

    end program bolam
!##################################################################################################################
    subroutine phicomp

! Computes geopotential

    use mod_model
    implicit none
    integer jlon, jlat, jklev
    real(4) zc1, zsigmed, zdsgalf, zsgalf

!  phi is geop. at integer levels

    zsigmed =.5*(1.+sigint(nlev))
    zsgalf  = zsigmed**3 * (1.+alfa*(1.-zsigmed)*(2.-zsigmed))
    zdsgalf = zsigmed**2 * (3.+alfa*(6.-12.*zsigmed+5.*zsigmed**2))
    do jlat = 1, nlat
    do jlon = 1, nlon
    zc1 = ( pzer-(pzer-ps(jlon,jlat))*zdsgalf )/( pzer*zsigmed-(pzer-ps(jlon,jlat))*zsgalf )
    phi(jlon,jlat,nlev) = phig(jlon,jlat) + rd*(1.-sigint(nlev))*zc1*tvirt(jlon,jlat,nlev)
    enddo
    enddo

    do jklev = nlev-1, 1, -1
    zsigmed = .5*(sigint(jklev+1)+sigint(jklev))
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
!##################################################################################################################
    subroutine rdmhf_atm (nun, init_flag)

! Reads the Model History File MHF with prognostic atmospheric variables

    use mod_model
    implicit none
    integer :: nun, init_flag, iunit=12, iunit_work=19, &
 jlon, jlat, jklev, isleep, no_input, ierr_open, comm, error
    character(len=30) :: filerd, file_work="input_request.txt",anun

#ifdef mpi
      include 'mpif.h'

      comm = mpi_comm_world
#endif

    no_input = 0
    isleep = 0

    if (myid.eq.0) then

      print *

      if (.not.nlclimate) then
        write (anun,'(i2.2)') nun
      else
        write (anun,'(i4.4)') nun
      endif

      filerd="input_atm_"//trim(anun)//".mhf"

      do while (.true.)
        open (iunit, file=trim(filerd), form='unformatted', status='old', iostat=ierr_open)
        if (ierr_open == 0) then
#ifdef oper
          call system("sync")
          call system("ls -l -L "//filerd)
          call system("date")
          call system("sleep 1")
#endif
          exit
        else
          print *,"Input file ",trim(filerd)," not ready: wait - sleep 10 s..."
          isleep = isleep + 1

#ifndef oper
          if (isleep.eq.2) then
            write(*, '(a, i3, a)') " Input file no.", nun, " not found: program stopped!"
            print*
            no_input = 1
            exit
          endif
#endif
          call system ("sleep 10")
        endif
      enddo

      if (nlclimate) then
        open (iunit_work, file=trim(file_work))
        write (iunit_work,*) nun+1
        print *,' ----- write ',nun+1,' in ',trim(file_work),'----'
        close (iunit_work)
      endif

    endif

#ifdef mpi
      call mpi_bcast (no_input, 1, mpi_integer, 0, comm, error)
#endif

    if (no_input.eq.1) then
#ifdef mpi
        call mpi_finalize (error)
#endif
      stop
    endif

!  read descriptor records and broadcast to all processes

    if (myid.eq.0) then
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
      if (init_flag == 1) call mpi_bcast (pdr0  , 200, mpi_real   , 0, comm, error)
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

! ps  surface pressure

    if (myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
    if (init_flag == 1) call disperse (gfield, ps  )
    if (init_flag /= 1) call disperse (gfield, psb2)

!  u, v component of the wind

    do jklev = 1, nlev
    if (myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
    if (init_flag == 1) call disperse (gfield, u  (1,1,jklev))
    if (init_flag /= 1) call disperse (gfield, ub2(1,1,jklev))
    enddo

    do jklev = 1, nlev
    if (myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
    if (init_flag == 1) call disperse (gfield, v  (1,1,jklev))
    if (init_flag /= 1) call disperse (gfield, vb2(1,1,jklev))
    enddo

!  temperature

    do jklev = 1, nlev
    if (myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
    if (init_flag == 1) call disperse (gfield, t  (1,1,jklev))
    if (init_flag /= 1) call disperse (gfield, tb2(1,1,jklev))
    enddo

!  specific humidity q

    do jklev = 1, nlev
    if (myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
    if (init_flag == 1) call disperse (gfield, q  (1,1,jklev))
    if (init_flag /= 1) call disperse (gfield, qb2(1,1,jklev))
    enddo

! cloud water and ice (summed) - provisionally attributed to cloud water

    do jklev = 1, nlev
    if (myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
    if (init_flag == 1) call disperse (gfield, qcw  (1,1,jklev))
    if (init_flag /= 1) call disperse (gfield, qcwb2(1,1,jklev))
    enddo

    if (myid == 0) then
      write(*,'(2a)') " read file ",trim(filerd)
      close (iunit)
      write (*,*)
    endif

    return
    end subroutine rdmhf_atm
!##################################################################################################################
    subroutine rdmhf_soil (nun, init_flag)

! Reads the Model History File MHF

    use mod_model
    implicit none
    integer :: nun, init_flag, iunit=13, jlon, jlat, jklev, isleep, no_input, ierr_open, ierr_read, ird, comm, error
    character(len=30) :: filerd, anun
    integer, dimension(50) :: nfdr_local
    real, dimension(200) :: pdr_local
    real, dimension(nlon,nlat) :: field2d_add

#ifdef mpi
      include 'mpif.h'

      comm = mpi_comm_world
#endif

    no_input = 0
    isleep = 0

    if (myid.eq.0) then

      write (*,*)

      if (.not.nlclimate) then
        write (anun,'(i2.2)') nun
      else
        write (anun,'(i4.4)') nun
      endif

      filerd="input_soil_"//trim(anun)//".mhf"

      do while (.true.)
        open (iunit, file=trim(filerd), form='unformatted', status='old', iostat=ierr_open)

        if (ierr_open == 0) then
#ifdef oper
          call system("sync")
          call system("ls -l -L "//filerd)
          call system("date")
          call system("sleep 1")
#endif
          exit
        else
          print *,"Input file ",trim(filerd)," not ready: wait - sleep 10 s..."
          isleep = isleep + 1

#ifndef oper
          if (isleep.eq.2) then
            write(*, '(a, i3, a)') " Input file no.", nun, " not found: program stopped!"
            print*
            no_input = 1
            exit
          endif
#endif

          call system ("sleep 10")
        endif
      enddo

    endif

#ifdef mpi
      call mpi_bcast (no_input, 1, mpi_integer, 0, comm, error)
#endif

    if (no_input.eq.1) then
#ifdef mpi
        call mpi_finalize (error)
#endif
      stop
    endif

!  read descriptor records and broadcast to all processes

    if (myid.eq.0) then

      ierr_read=0

      if (init_flag == 1) then
        read(iunit) nfdr_local
        read(iunit) pdr_local
        if (any(nfdr_local(:) /= nfdr0(:)).or.any(pdr_local(:) /= pdr0(:))) then
          ierr_read = 1
          print *, ' Header parameter (nfdr, pdr) of ',trim(filerd), &
 ' input file not coincide with defined parameters.   STOP the program'
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
    if (init_flag == 1.or.nlclimate) call disperse (gfield, iceth_rd)
    if (init_flag == 1) iceth(:,:)=iceth_rd(:,:)

    if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
    if (init_flag == 1.or.nlclimate) call disperse (gfield, fice_rd)
    if (init_flag == 1) fice(:,:)=fice_rd(:,:)

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
    if (init_flag == 1) call disperse (gfield, totpre)
    if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
    if (init_flag == 1) call disperse (gfield, conpre)
    if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
    if (init_flag == 1) call disperse (gfield, snfall)

! Prognostic surface and soil/sea fields

    if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
    if (init_flag == 1) call disperse (gfield, tskin)
    if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
    if (init_flag == 1) call disperse (gfield, tg_surf)

    if (init_flag == 1) then
      do jklev = 1, nlevg
        if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
        call disperse (gfield, tg(1,1,jklev))
      enddo
    else
      do jklev = 1, nlevg-1
        if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
      enddo
      if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
      call disperse (gfield, tgclim)
    endif

    if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
    if (init_flag == 1)  then
      call disperse (gfield, qskin)
    endif
    if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
    if (init_flag == 1)  call disperse (gfield, qg_surf)

    if (init_flag == 1) then
      do jklev = 1, nlevg
        if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
        call disperse (gfield, qg(1,1,jklev))
      enddo
    else
      do jklev = 1, nlevg-1
        if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
      enddo
      if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
      call disperse (gfield, qgclim)
    endif

    if (init_flag == 1) then
      if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
      call disperse (gfield, fice_soil_surf)
      do jklev = 1, nlevg
        if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
        call disperse (gfield, fice_soil(1,1,jklev))
      enddo
      if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
      call disperse (gfield, snow)
      do jklev = 1, nlev_snow
        if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
        call disperse (gfield, snow_lev(1,1,jklev))
      enddo
      do jklev = 1, nlev_snow
        if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
        call disperse (gfield, snow_t(1,1,jklev))
      enddo
      do jklev = 1, nlev_snow
        if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
        call disperse (gfield, snow_fice(1,1,jklev))
      enddo
      do jklev = 1, nlev_snow
        if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
        call disperse (gfield, snow_age(1,1,jklev))
      enddo
      do jklev = 1, nlev_snow
        if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
        call disperse (gfield, snow_melt_age(1,1,jklev))
      enddo
      do jklev = 1, nlev_snow
        if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
        call disperse (gfield, snow_dens(1,1,jklev))
      enddo
      if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
      call disperse (gfield, alsn)
      if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
      call disperse (gfield, cswfl)
      if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
      call disperse (gfield, clwfl)
      if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
      call disperse (gfield, chflux)
      if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
      call disperse (gfield, cqflux)
      if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
      call disperse (gfield, t2min)
      if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
      call disperse (gfield, t2max)
      if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
      call disperse (gfield, runoff)
      if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
      call disperse (gfield, runoff_tot)
   
! Reading of additional 2D fields
   
      do ird = 1, 3
        if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
      enddo

    endif ! init_flag == 1

    if (myid == 0) then
      close (iunit)
      write(*,'(2a)') " read file ",trim(filerd)
      write(*,*)
    endif

    if (nun >= 3) flag_pochva_par_change=1

    return
    end subroutine rdmhf_soil
!##################################################################################################################
subroutine rd_param_const

! Reads from additional input file
! all constant (in time) model physiographical parameters

use mod_model
implicit none

character(len=30) :: filerd="model_param_constant.bin"
integer :: iunit=11, ierr_open, nlon_local, nlat_local, nlevg_local, nst_local, nvt_local, ierr=0, ird, jklev, &
 comm, error
real :: dlon_local, dlat_local, x0_local, y0_local, alon0_local, alat0_local
integer, dimension(gnlon,gnlat) :: igfield
real, dimension(1) :: z1

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

   read (iunit) nlon_local, nlat_local, nlevg_local, dlon_local, dlat_local, &
 x0_local, y0_local, alon0_local, alat0_local, nst_local, nvt_local

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
   if (nst_local /= nst) ierr=ierr+1
   if (nvt_local /= nvt) ierr=ierr+1

   if (ierr /= 0) then
     print *,"Error in header parameters in input file, ",trim(filerd),", not coincident with defined parameters"
     print *,"Model nlon, nlat, nlevg, dlon, dlat, x0, y0, alon0, alat0, nst, nvt :", &
 gnlon, gnlat, nlevg, pdr0(2), pdr0(1), pdr0(39), pdr0(38), pdr0(5), pdr0(4), nst, nvt
     print *,"Read nlon, nlat, nlevg, dlon, dlat, x0, y0, alon0, alat0, nst, nvt :", &
 nlon_local, nlat_local, nlevg_local, dlon_local, dlat_local, x0_local, y0_local, &
 alon0_local, alat0_local, nst_local, nvt_local
     print *,"   STOP"
#ifdef mpi
       call mpi_abort(comm, error)
#endif
     stop
   endif

 endif ! myid == 0

!  land sea mask and orography

 if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, fmask)
 if (myid == 0) then
   call rrec2 (iunit, gnlon, gnlat, gfield)
   gfield(:,:) = gfield(:,:)*g
 endif
 call disperse (gfield, phig)

! Standard deviation computed from orography hight variance

 if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
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
 if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, veg_rough)

!  radiation parameters of vegetation surface

 if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, albedo_veg)
 if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, emiss1_veg)
 if (myid == 0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, emiss2_veg)

! soil top and bottom parameters

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
   write(*,*) " read file ",trim(filerd)
   write(*,*)
 endif

return
end subroutine rd_param_const
!##################################################################################################################
    subroutine rrec2 (kunit, nlon, nlat, vect)

! Used by rdmhf, rdrf

    implicit none

    integer :: kunit, nlon, nlat
    real(4), dimension(nlon,nlat) :: vect

    read(kunit) vect(1:nlon, 1:nlat)

    return
    end subroutine rrec2
!##################################################################################################################
    subroutine rrec2_int (kunit, nlon, nlat, ivect)

! Used by rdmhf, rdrf

    implicit none

    integer :: kunit, nlon, nlat
    integer, dimension(nlon,nlat) :: ivect

    read(kunit) ivect(1:nlon, 1:nlat)

    return
    end subroutine rrec2_int
!##################################################################################################################
    subroutine rrec2_old (kunit, nlon, nlat, vect)

! Used by rdmhf, rdrf

    real(4) vect(nlon,nlat)

    do jlat = 1, nlat
    read(kunit) (vect(jlon,jlat), jlon = 1, nlon)
    enddo

    return
    end subroutine rrec2_old
!##################################################################################################################
    subroutine rrec2_int_old (kunit, nlon, nlat, ivect)

! Used by rdmhf, rdrf

    integer ivect(nlon,nlat)

    do jlat = 1, nlat
    read(kunit) (ivect(jlon,jlat), jlon = 1, nlon)
    enddo

    return
    end subroutine rrec2_int_old
!##################################################################################################################
subroutine wrmhf_atm(name, nfdr, pdr, ilon1, ilon2, jlat1, jlat2, flag_lon)

! Writes the Model History File (MHF) with prognostic atmospheric variables

use mod_model, only : imhf, gnlon, gnlat, nlon, nlat, nlev, myid, nlclimate, gfield, &
                      ps, u, v, t, q, qcw, qci

implicit none

character(len=30) :: name
integer :: iunit=22, iunit_work=29, ilon1, ilon2, jlat1, jlat2, flag_lon, jlon, jlat, jklev
integer, dimension(50)  :: nfdr
real(4), dimension(200) :: pdr
real, dimension(nlon,nlat) :: zwork, zfrunoff
character(len=50) :: file_out, amhf

 if(myid == 0) then

#ifdef oper
 if (.not.nlclimate) then
   write (amhf,'(i3.3)') imhf
 else
   write (amhf,'(i4.4)') imhf
 endif
 file_out=trim(name)//'_atm_'//trim(amhf)//'.mhf'
 open (iunit, file=trim(file_out), form='unformatted', status='unknown')
#else
! call system ("touch mhf-open")
 file_out=trim(name)//"_atm.mhf"
 open (iunit, file=trim(file_out), form='unformatted', status='unknown', position='append')
#endif

!  write descriptor records

 write(iunit) nfdr
 write(iunit) pdr

 endif

!  surface pressure

 call collect (ps, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

!  u, v component of the wind

 do jklev=1,nlev
 call collect (u(1,1,jklev), gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
 enddo
 do jklev=1,nlev
 call collect (v(1,1,jklev), gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
 enddo

!  temperature

 do jklev=1,nlev
 call collect (t(1,1,jklev), gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
 enddo

!  specific humidity q

 do jklev=1,nlev
 call collect (q(1,1,jklev), gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
 enddo

! cloud water and ice (summed)

 do jklev=1,nlev
 zwork(:,:) = qcw(:,:,jklev) + qci(:,:,jklev)
 call collect (zwork, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
 enddo

 if (myid == 0) then

   flush (iunit)
   close (iunit)

   print *

#ifdef oper
   call system("sync")
!!!   call system("ls -l -L "//file_out)
!!!   call system("date")
   open  (iunit_work, file=trim(file_out)//'.txt', status='unknown')
   write (iunit_work,'(2a)') trim(file_out),' is full and close'
   close (iunit_work)
   call system("sync")
#endif
   print *,'Output written on unit ',trim(file_out)
!      call system ("rm -f mhf-open")

 endif

return
end subroutine wrmhf_atm
!##################################################################################################################

subroutine wrmhf_soil(name, nfdr, pdr, ilon1, ilon2, jlat1, jlat2, flag_lon, iflag_soil_full)

! Writes the Model History File MHF

use mod_model, only : imhf, gnlon, gnlat, nlon, nlat, nlevg, nlev_snow, myid, nlclimate, gfield, &
                      lai, fveg, rgm, rgq, iceth, fice, albedo, emisg1, emisg2,&
                      cloudt, totpre, conpre, snfall, &
                      tskin, tg_surf, tg, qskin, qg_surf, qg, fice_soil_surf, fice_soil, &
                      snow, snow_lev, snow_t, snow_fice, snow_age, snow_melt_age, snow_dens, alsn, &
                      cswfl, clwfl, chflux, cqflux, t2min, t2max, runoff, runoff_tot, &
                      cwvflux, fl_heat_soil_bottom, fl_water_soil_bottom

implicit none

character(len=30) :: name
integer :: iunit=23, iunit_work=29, ilon1, ilon2, jlat1, jlat2, flag_lon, iflag_soil_full, jlon, jlat, jklev
integer, dimension(50)  :: nfdr
real(4), dimension(200) :: pdr
real, dimension(nlon,nlat) :: zwork, zfrunoff
character(len=50) :: file_out, amhf, aini, aterm

 if(myid == 0) then

#ifdef oper
   if (.not.nlclimate) then
     write (amhf,'(i3.3)') imhf
   else
     write (amhf,'(i4.4)') imhf
   endif
   file_out=trim(name)//'_soil_'//trim(amhf)//'.mhf'
   if (iflag_soil_full == 1) then
     write (aini,'(i4.4,3i2.2)') nfdr(5:8)
     write (aterm,'(i3.3,2i2.2)') nfdr(10:12)
     file_out=trim(name)//'_soil_'//trim(aini)//'_'//trim(aterm)//'.mhf'
   endif
   open (iunit, file=trim(file_out), form='unformatted', status='unknown')
#else
!   call system ("touch mhf-open")
   file_out=trim(name)//"_soil.mhf"
   open (iunit, file=trim(file_out), form='unformatted', status='unknown', position='append')
#endif

!  write descriptor records

   write(iunit) nfdr
   write(iunit) pdr

 endif

!  Physiographical parameters changing in time

 call collect (lai, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (fveg, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (rgm, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (rgq, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (iceth, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (fice, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (albedo, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (emisg1, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (emisg2, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

! Prognostic cloud and precipitation variables at the surface

 call collect (cloudt, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (totpre, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (conpre, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (snfall, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

! Prognostic surface and soil/sea fields

 call collect (tskin, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (tg_surf, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 do jklev = 1, nlevg
   call collect (tg(1,1,jklev), gfield)
   if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
 enddo

 call collect (qskin, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (qg_surf, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 do jklev = 1, nlevg
   call collect (qg(1,1,jklev), gfield)
   if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
 enddo

 call collect (fice_soil_surf, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 do jklev = 1, nlevg
   call collect (fice_soil(1,1,jklev), gfield)
   if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
 enddo

 call collect (snow, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 do jklev = 1, nlev_snow
   call collect (snow_lev(1,1,jklev), gfield)
   if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
 enddo

 do jklev = 1, nlev_snow
   call collect (snow_t(1,1,jklev), gfield)
   if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
 enddo

 do jklev = 1, nlev_snow
   call collect (snow_fice(1,1,jklev), gfield)
   if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
 enddo

 do jklev = 1, nlev_snow
   call collect (snow_age(1,1,jklev), gfield)
   if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
 enddo

 do jklev = 1, nlev_snow
   call collect (snow_melt_age(1,1,jklev), gfield)
   if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
 enddo

 do jklev = 1, nlev_snow
   call collect (snow_dens(1,1,jklev), gfield)
   if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
 enddo

 call collect (alsn, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (cswfl, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (clwfl, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (chflux, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (cqflux, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (t2min, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (t2max, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (runoff, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (runoff_tot, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

! Writting of additional 2D fields

!    zwork = 0.
!    call collect (zwork, gfield)
!    if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
!    if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
!    if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (cwvflux, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (fl_heat_soil_bottom, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (fl_water_soil_bottom, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 if (myid == 0) then

   flush (iunit)
   close (iunit)

   print *

#ifdef oper
   call system("sync")
!!!   call system("ls -l -L "//file_out)
!!!   call system("date")
   open  (iunit_work, file=trim(file_out)//'.txt', status='unknown')
   write (iunit_work,'(2a)') trim(file_out),' is full and close'
   close (iunit_work)
   call system("sync")
#endif
   print *,'Output written on unit ',trim(file_out)
!      call system ("rm -f mhf-open")

 endif

return
end subroutine wrmhf_soil
!######################################################################################################

subroutine wr_param_const(nfdr, pdr, ilon1, ilon2, jlat1, jlat2, flag_lon)

! Writes file with all constant (in time) model physiographical parameters.
! It is used only by Globo to write zoomed domain.

use mod_model, only : gnlon, gnlat, nlon, nlat, nlevg, nst, nvt, myid, gfield, &
                      fmask, phig, g, orogvar, suolo, vegeta, &
                      qsoil_max, qsoil_min, c_soil, rho_soil, psi_soil, &
                      kw_soil, par_b_soil, par_c_soil, qsrel_wilt, qsrel_ref, &
                      albedo_soil_dry, albedo_soil_wet, &
                      emiss1_soil_dry, emiss1_soil_wet, emiss2_soil_dry, emiss2_soil_wet, &
                      proot, veg_rough, albedo_veg, emiss1_veg, emiss2_veg, &
                      water_table_depth, tg_bottom, qg_rel_bottom, qg_rel_surf_approx, &
                      ind_lev_soil_h_bottom, ind_lev_soil_w_bottom, snow_dirt, lai_max 

implicit none

integer :: iunit=21, ilon1, ilon2, jlat1, jlat2, flag_lon, jlon, jlat, jklev, iwr
integer, dimension(50)  :: nfdr
real(4), dimension(200) :: pdr
integer, dimension(gnlon,gnlat) :: igfield
character(len=50) :: file_out="model_param_constant_zoom.bin"

 if(myid == 0) then

   open (iunit, file=trim(file_out), form='unformatted', status='unknown')

! T-grid parameters: gnlon, gnlat, nlevg, dlon, dlat, x0, y0, alon0, alat0, nst, nvt

   write(iunit) nfdr(2), nfdr(3), nfdr(15), pdr(2), pdr(1), &
 pdr(39), pdr(38), pdr(5), pdr(4)+pdr(1)*0.5, nst, nvt

 endif

!  land sea mask and orography

 call collect (fmask, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (phig(:,:)/g, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

! Standard deviation computed from orography hight variance

 call collect (orogvar, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

!  soil types

 do iwr = 1, nst+1
   call collect (suolo(1,1,iwr), gfield)
   if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
 enddo

!  vegetation types

 do iwr = 1, nvt+1
   call collect (vegeta(1,1,iwr), gfield)
   if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
 enddo

!  physical soil parameters at soil levels

 do jklev = 1, nlevg
   call collect (qsoil_max(1,1,jklev), gfield)
   if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
 enddo
 do jklev = 1, nlevg
   call collect (qsoil_min(1,1,jklev), gfield)
   if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
 enddo
 do jklev = 1, nlevg
   call collect (c_soil(1,1,jklev), gfield)
   if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
 enddo
 do jklev = 1, nlevg
   call collect (rho_soil(1,1,jklev), gfield)
   if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
 enddo
 do jklev = 1, nlevg
   call collect (psi_soil(1,1,jklev), gfield)
   if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
 enddo
 do jklev = 1, nlevg
   call collect (kw_soil(1,1,jklev), gfield)
   if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
 enddo
 do jklev = 1, nlevg
   call collect (par_b_soil(1,1,jklev), gfield)
   if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
 enddo
 do jklev = 1, nlevg
   call collect (par_c_soil(1,1,jklev), gfield)
   if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
 enddo
 do jklev = 1, nlevg
   call collect (qsrel_wilt(1,1,jklev), gfield)
   if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
 enddo
 do jklev = 1, nlevg
   call collect (qsrel_ref(1,1,jklev), gfield)
   if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
 enddo

!  radiation parameters of soil surface

 call collect (albedo_soil_dry, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
 call collect (albedo_soil_wet, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
 call collect (emiss1_soil_dry, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
 call collect (emiss1_soil_wet, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
 call collect (emiss2_soil_dry, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
 call collect (emiss2_soil_wet, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

!  physical parameters of vegetation

 call collect (proot, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
 call collect (veg_rough, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

!  radiation parameters of vegetation surface

 call collect (albedo_veg, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
 call collect (emiss1_veg, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
 call collect (emiss2_veg, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

! soil top and bottom parameters

 call collect (water_table_depth, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
 call collect (tg_bottom, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
 call collect (qg_rel_bottom, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
 call collect (qg_rel_surf_approx, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

!  soil level index of bottom contidion for heat transport and water transport schemes

 call collect_int (ind_lev_soil_h_bottom, igfield)
 if(myid == 0) call wrec2_int (iunit, gnlon, gnlat, igfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
 call collect_int (ind_lev_soil_w_bottom, igfield)
 if(myid == 0) call wrec2_int (iunit, gnlon, gnlat, igfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

! snow_dirt: snow "dirtibility", that is weight (proportion 0-1) of dirt
! growth of snow surface due to vegetation waste, aerosol deposition, etc., 
! it is used in snow albedo definition.

 call collect (snow_dirt, gfield)
 if(myid == 0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

!  maximum values of LAI in used dataset (for a normalization)

 if(myid == 0) write(iunit) lai_max

 if (myid == 0) then
   flush (iunit)
   close (iunit)
   print *
   print *,'Output written on unit ',trim(file_out)
   print *
 endif

return
end subroutine wr_param_const
!######################################################################################################
    subroutine wrec2 (kunit, nlon, nlat, vect, ilon1, ilon2, jlat1, jlat2, flag_lon)

! Used by wrmhf_atm, wrmhf_soil, wrshf_b, wrshf_g

    implicit none

    integer :: kunit, nlon, nlat, ilon1, ilon2, jlat, jlat1, jlat2, flag_lon, jlonf1, jlonf2
    real(4), dimension(nlon,nlat) :: vect, vect2

#ifdef globo
    vect(1   ,jlat1:jlat2) = vect(nlon-1,jlat1:jlat2)
    vect(nlon,jlat1:jlat2) = vect(2     ,jlat1:jlat2)
#endif

    if (flag_lon == 0) then

      write(kunit) vect(ilon1:ilon2, jlat1:jlat2)

    else
       
      do jlat = jlat1, jlat2
        jlonf1 = nlon-ilon1
        jlonf2 = nlon-ilon1+ilon2-1
        vect2(1:jlonf1, jlat) = vect(ilon1:nlon-1, jlat)
        vect2(jlonf1+1:jlonf2, jlat) = vect(2:ilon2, jlat)
      enddo

      write(kunit) vect2(1:jlonf2, jlat1:jlat2)

    endif

    return
    end subroutine wrec2
!##################################################################################################################
    subroutine wrec2_int (kunit, nlon, nlat, ivect, ilon1, ilon2, jlat1, jlat2, flag_lon)

! Used by wrrf

    implicit none

    integer :: kunit, nlon, nlat, ilon1, ilon2, jlat, jlat1, jlat2, flag_lon, jlonf1, jlonf2
    integer, dimension(nlon,nlat) :: ivect, ivect2

#ifdef globo
    ivect(1   ,jlat1:jlat2) = ivect(nlon-1,jlat1:jlat2)
    ivect(nlon,jlat1:jlat2) = ivect(2     ,jlat1:jlat2)
#endif

    if (flag_lon == 0) then

      write(kunit) ivect(ilon1:ilon2, jlat1:jlat2)

    else
       
      do jlat = jlat1, jlat2
        jlonf1 = nlon-ilon1
        jlonf2 = nlon-ilon1+ilon2-1
        ivect2(1:jlonf1, jlat) = ivect(ilon1:nlon-1, jlat)
        ivect2(jlonf1+1:jlonf2, jlat) = ivect(2:ilon2, jlat)
      enddo

      write(kunit) ivect2(1:jlonf2, jlat1:jlat2)

    endif

    return
    end subroutine wrec2_int
!##################################################################################################################
    subroutine wrec2_old (kunit, nlon, nlat, vect, ilon1, ilon2, jlat1, jlat2, flag_lon)

! Used by wrmhf_atm, wrmhf_soil, wrshf_b, wrshf_g

    implicit none

    integer :: kunit, nlon, nlat, ilon1, ilon2, jlat1, jlat2, flag_lon, jlat, jlon
    real(4), dimension(nlon,nlat) :: vect

#ifdef globo
    vect(1   ,jlat1:jlat2) = vect(nlon-1,jlat1:jlat2)
    vect(nlon,jlat1:jlat2) = vect(2     ,jlat1:jlat2)
#endif

    if (flag_lon == 0) then

    do jlat = jlat1, jlat2
    write(kunit) (vect(jlon,jlat), jlon = ilon1, ilon2)
    enddo

    else

    do jlat = jlat1, jlat2
    write(kunit) (vect(jlon,jlat), jlon = ilon1,nlon-1), (vect(jlon,jlat), jlon = 2,ilon2)
    enddo

    endif

    return
    end subroutine wrec2_old
!##################################################################################################################
    subroutine wrec2_int_old (kunit, nlon, nlat, ivect, ilon1, ilon2, jlat1, jlat2, flag_lon)

! Used by wrrf

    implicit none

    integer :: kunit, nlon, nlat, ilon1, ilon2, jlat1, jlat2, flag_lon, jlat, jlon
    integer, dimension(nlon,nlat) :: ivect

#ifdef globo
    ivect(1   ,jlat1:jlat2) = ivect(nlon-1,jlat1:jlat2)
    ivect(nlon,jlat1:jlat2) = ivect(2     ,jlat1:jlat2)
#endif

    if (flag_lon == 0) then

    do jlat = jlat1, jlat2
    write(kunit) (ivect(jlon,jlat), jlon = ilon1, ilon2)
    enddo

    else

    do jlat = jlat1, jlat2
    write(kunit) (ivect(jlon,jlat), jlon = ilon1,nlon-1), (ivect(jlon,jlat), jlon = 2,ilon2)
    enddo

    endif

    return
    end subroutine wrec2_int_old
!##################################################################################################################
    subroutine collect (lfield, gfield)

! Routine for parallel computing: collects data from sub-domains

!  root process receives all sub-domains to reconstruct the global domain
!  other processes send their sub-domains to root

    use mod_model, only : nprocsx, nprocsy, nlon, nlat, gnlon, gnlat, myid
    implicit none

    integer                         :: count, sx, ex, sy, ey, jpr, comm, error, tag=0
    real(4), dimension(nlon,nlat)   :: lfield
    real(4), dimension(gnlon,gnlat) :: gfield

#ifdef mpi

      include 'mpif.h'
      integer                       :: status(mpi_status_size)

      comm = mpi_comm_world

      if (myid.eq.0) then
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
!##################################################################################################################
    subroutine collect_int (ilfield, igfield)

! Routine for parallel computing: collects data from sub-domains (for integer 2d field)

!  root process receives all sub-domains to reconstruct the global domain
!  other processes send their sub-domains to root

    use mod_model, only : nprocsx, nprocsy, nlon, nlat, gnlon, gnlat, myid
    implicit none
    integer                         :: count, sx, ex, sy, ey, jpr, comm, error, tag=0
    integer, dimension(nlon,nlat)   :: ilfield
    integer, dimension(gnlon,gnlat) :: igfield

#ifdef mpi

      include 'mpif.h'
      integer                       :: status(mpi_status_size)

      comm = mpi_comm_world

      if (myid.eq.0) then
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
          CALL MPI_RECV (igfield(sx:ex,sy:ey), count, MPI_INTEGER, jpr, jpr, comm, status, error)
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
        CALL MPI_SEND (ilfield(sx:ex,sy:ey), count, MPI_INTEGER, 0, myid, comm, error)
      endif

      call mpi_barrier(comm, error)

#else

      igfield(1:nlon,1:nlat) = ilfield(1:nlon,1:nlat)

#endif

    return
    end subroutine collect_int
!##################################################################################################################
      subroutine disperse (gfield, lfield)

! Routine for parallel computing: distributes data to sub-domains

!  root process reads the global domain and distributes local domains to other processes
!  other processes receive their sub-domains from root

    use mod_model, only : myid, nprocsx, nprocsy, nlon, nlat, gnlon, gnlat
    implicit none

    integer                          :: sx, ex, sy, ey, jpr, comm, error, tag=0, count
    real(4), dimension(nlon,nlat)    :: lfield
    real(4), dimension(gnlon,gnlat)  :: gfield

#ifdef mpi

      include 'mpif.h'
      integer                        :: status(mpi_status_size)

      count = nlon*nlat
      comm = mpi_comm_world
      call mpi_barrier(comm, error)
      if (myid.eq.0) then
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
!##################################################################################################################
     subroutine disperse_int (igfield, ilfield)

! Routine for parallel computing: distributes data to sub-domains (for integer 2d fields)

!  root process reads the global domain and distributes local domains to other processes
!  other processes receive their sub-domains from root

    use mod_model, only : myid, nprocsx, nprocsy, nlon, nlat, gnlon, gnlat
    implicit none

    integer                          :: sx, ex, sy, ey, jpr, comm, error, tag=0, count
    integer, dimension(nlon,nlat)    :: ilfield
    integer, dimension(gnlon,gnlat)  :: igfield

#ifdef mpi

      include 'mpif.h'
      integer                        :: status(mpi_status_size)

      count = nlon*nlat
      comm = mpi_comm_world
      if (myid.eq.0) then
        ilfield(1:nlon,1:nlat) = igfield(1:nlon,1:nlat)
        do jpr = 1, nprocsx*nprocsy-1
          sx = 1+(jpr/nprocsy)*(nlon-2)
          ex = sx+nlon-1
          sy = 1+(jpr-(jpr/nprocsy)*nprocsy)*(nlat-2)
          ey = sy+nlat-1
          call mpi_send (igfield(sx:ex,sy:ey), count, mpi_integer, jpr, tag, comm, error)
        enddo
      else
        call mpi_recv (ilfield, count, mpi_integer, 0, tag, comm, status, error)
      endif

      call mpi_barrier(comm, error)

#else

      ilfield(1:nlon,1:nlat) = igfield(1:nlon,1:nlat)

#endif

    return
    end subroutine disperse_int
!##################################################################################################################
    subroutine diverg (kadj)

!  Definition of horizontal divergence
!  Definition of tendency of ps (pst)
!  Definition of vertical velocity (sigdot)
!  Computation of omega/p (omeg)

    use mod_model, only : nlon, nlat, nlev, nlevp1, nlonm1, nlatm1, ip_null, ip_e, ip_n, ip_s, ip_w, nadj,  &
                          ps, pst, u, v, omeg, sigdot, sig, sigint, hxt, hxv, dsig, dx, dy, pzer,           &
                          div1, div2, alfa, sigalf, gnlon, myid, pi, a, dlon, dlat, nprocsx, nprocsy, cpole
    implicit none

    integer :: jlon, jlat, jklev, kadj, jstart, jend, jpr, ierr, comm, tag1=1, tag2=2
    real(4) zpbs, zpbn, zrdx, zhxvt, zhxvtn, zdpdsig, zzpp, zdiv(nlon,nlev), zom1(nlon,nlev)
    real(4) zsgalf, zdsgalf, zdvint(nlon,nlevp1), zpbx(nlon), zpby(nlon), zpbyn(nlon), zzz(nlev)

#ifdef mpi
      include 'mpif.h'
      integer :: status(mpi_status_size)

      comm = mpi_comm_world
#endif

    pst = 0.
    if (kadj.eq.1) sigdot = 0.

#ifdef globo

!---------------
!  div1 at poles
!---------------

#ifdef mpi
      if (ip_s.eq.ip_null) then
#endif
    div1(1,1,:) = 0.
    do jlon = 2, nlonm1
      div1(1,1,:) = div1(1,1,:)+v(jlon,2,:)*cpole
    enddo
#ifdef mpi
      if (myid.ne.0) then
        call mpi_send (div1(1,1,1:nlev), nlev, mpi_real, 0, tag1, comm, ierr)
        call mpi_recv (div1(1,1,1:nlev), nlev, mpi_real, 0, tag2, comm, status, ierr)
      else
        do jpr = nprocsy, nprocsy*nprocsx-1, nprocsy
          call mpi_recv (zzz, nlev, mpi_real, jpr, tag1, comm, status, ierr)
          div1(1,1,:) = div1(1,1,:) + zzz(:)
        enddo
        do jpr = nprocsy, nprocsy*nprocsx-1, nprocsy
          call mpi_send (div1(1,1,1:nlev), nlev, mpi_real, jpr, tag2, comm, ierr)
        enddo
      endif
#endif
    do jlon = 2, nlon
    div1(jlon,1,:) = div1(1,1,:)
    enddo
#ifdef mpi
      endif
#endif

#ifdef mpi
      if (ip_n.eq.ip_null) then
#endif
    div1(1,nlat,:) = 0.
    do jlon = 2, nlonm1
      div1(1,nlat,:) = div1(1,nlat,:)-v(jlon,nlat,:)*cpole
    enddo
#ifdef mpi
      if (myid.ne.nprocsy-1) then
        call mpi_send (div1(1,nlat,1:nlev), nlev, mpi_real, nprocsy-1, tag1, comm, ierr)
        call mpi_recv (div1(1,nlat,1:nlev), nlev, mpi_real, nprocsy-1, tag2, comm, status, ierr)
      else
        do jpr = 2*nprocsy-1, nprocsy*nprocsx-1, nprocsy
          call mpi_recv (zzz, nlev, mpi_real, jpr, tag1, comm, status, ierr)
          div1(1,nlat,:) = div1(1,nlat,:) + zzz(:)
        enddo
        do jpr = 2*nprocsy-1, nprocsy*nprocsx-1, nprocsy
          call mpi_send (div1(1,nlat,1:nlev), nlev, mpi_real, jpr, tag2, comm, ierr)
        enddo
      endif
#endif
    do jlon = 2 ,nlon
      div1(jlon,nlat,:) = div1(1,nlat,:)
    enddo
#ifdef mpi
      endif
#endif

#endif ! ifdef globo

!-----------------------------
!  div1 in the interior points
!-----------------------------

    do 15 jlat = 2, nlatm1
    zrdx   = 1.         /(hxt(jlat)*dx)
    zhxvt  = hxv(jlat  )/(hxt(jlat)*dy)
    zhxvtn = hxv(jlat+1)/(hxt(jlat)*dy)
    do jklev = 1, nlev
    do jlon = 2, nlonm1
    div1(jlon,jlat,jklev) = (u(jlon,jlat,jklev)-u(jlon-1,jlat,jklev))*zrdx  &
                          + zhxvtn*v(jlon,jlat+1,jklev)-zhxvt*v(jlon,jlat,jklev)
    enddo
    enddo
 15 continue

    call wafps   ! computes div2, divergence of mass fluxes

#ifdef globo
    call filt2t1  (div1, .1)  ! only a small value is good...
    call filt2t1  (div2, .1)
    call polavert (div1)
    call polavert (div2)
#else
    if (ip_w.eq.ip_null) div1(1,2:nlatm1,:)    = div1(2,2:nlatm1,:)
    if (ip_w.eq.ip_null) div2(1,2:nlatm1,:)    = div2(2,2:nlatm1,:)
    if (ip_e.eq.ip_null) div1(nlon,2:nlatm1,:) = div1(nlonm1,2:nlatm1,:)
    if (ip_e.eq.ip_null) div2(nlon,2:nlatm1,:) = div2(nlonm1,2:nlatm1,:)
    if (ip_s.eq.ip_null) div1(:,1,:)           = div1(:,2,:)
    if (ip_s.eq.ip_null) div2(:,1,:)           = div2(:,2,:)
    if (ip_n.eq.ip_null) div1(:,nlat,:)        = div1(:,nlatm1,:)
    if (ip_n.eq.ip_null) div2(:,nlat,:)        = div2(:,nlatm1,:)

!    call filt2t (div1, 0.8)
!    call filt2t (div2, 0.8)

! The following filter seems to allow for higher (conv.) precip. maxima for bolam...

    call filt (div1, nlev, 0.77)
    call filt (div2, nlev, 0.77)
#endif

!  definition of tendency of ps and of (-) integral from top to a given half-level of horizontal div. (zdvint)

    jstart = 2
    jend = nlatm1

#ifdef globo
    if (ip_s.eq.ip_null) jstart = 1
    if (ip_n.eq.ip_null) jend = nlat
#endif

    do 20 jlat = jstart, jend

    do jklev = 1, nlev
    zdsgalf = sigint(jklev)**2 * (3.+alfa*(6.-12.*sigint(jklev)+5.*sigint(jklev)**2))
    do jlon = 2, nlonm1
    zdiv(jlon,jklev) = pzer*(1.-zdsgalf)*div1(jlon,jlat,jklev) + zdsgalf*div2(jlon,jlat,jklev)
    zom1(jlon,jklev) = div2(jlon,jlat,jklev)-ps(jlon,jlat)*div1(jlon,jlat,jklev)  ! for comput. of omega
    pst(jlon,jlat) = pst(jlon,jlat) - zdiv(jlon,jklev)*dsig(jklev)
    zdvint(jlon,jklev+1) = pst(jlon,jlat)
    enddo
    enddo

!  definition of vertical velocity sigdot on half-levels

    do jklev = 2, nlev
    zsgalf  = sig(jklev)**3 * (1.+alfa*(1.-sig(jklev))*(2.-sig(jklev)))
    zdsgalf = sig(jklev)**2 * (3.+alfa*(6.-12.*sig(jklev)+5.*sig(jklev)**2))
    do jlon = 2, nlonm1
    zdpdsig = pzer - (pzer-ps(jlon,jlat))*zdsgalf
    sigdot(jlon,jlat,jklev) = sigdot(jlon,jlat,jklev) + (zdvint(jlon,jklev)-zsgalf*pst(jlon,jlat)) &
                              /(float(nadj)*zdpdsig)
    enddo
    enddo

!  omega/p (omeg) at t-points and integer levels

    do jklev = 1, nlev
    do jlon = 2, nlonm1
    zzpp = pzer*sigint(jklev)-(pzer-ps(jlon,jlat))*sigalf(jklev)
    omeg(jlon,jlat,jklev) = ( sigalf(jklev)*zom1(jlon,jklev) + zdvint(jlon,jklev+1) &
                            + zdiv(jlon,jklev)*(sig(jklev+1)-sigint(jklev)) )/zzpp
    enddo
    enddo

 20 continue

    return
    end subroutine diverg
!##################################################################################################################
    subroutine wafuv (dt)

! WAF advection of horizontal velocity components

    use mod_model, only : nlon, gnlon, nlat, nlev, nlevp1, nlonm1, nlatm1, u, v,                       &
                          sigdot, hxt, hst, hxv, dsig, snt, cst, pi, dx, dy, a, dlon, dlat, ut,        &
                          vt, nsweep, nprocsy, nprocsx, myid, ip_e, ip_n, ip_s, ip_w, ip_null, ip_oppo, denomin
    implicit none

    real(4), dimension(nlon)               :: zpbw1, zpbw2
    real(4), dimension(nlon,nlev)          :: zvt
    real(4), dimension(nlat)               :: zut
    real(4), dimension(nlat,nlev)          :: zutb
    real(4), dimension(nlon,nlevp1)        :: wfw1, wfw2
    real(4), dimension(nlon,nlat)          :: zpby1, zpby2
    real(4), dimension(nlon,0:nlat+1,nlev) :: wz1, wz2
    real(4), dimension(0:nlon+1,nlat,nlev) :: p01, p02
    integer jlon, jlat, jklev, jlon1, j1, j1m1, is, jsweep, tag3, jstart, jend, iprocs, jpr
    real(4) dt, zc, zcost, zamu, denr1, denr2, r1, r2, b1, b2, zphi1, zphi2
    real(4) zvi(2,nlev), zvo(2,nlev), zdiv, zcostx, zcosty, zcostg, zhxvt, zhxvtn
    integer                                :: ierr, comm, tag1=1, tag2=2

#ifdef mpi
      include 'mpif.h'
      integer                              :: status(mpi_status_size)

      comm = mpi_comm_world
#endif

    jstart = 2
    jend = nlatm1

#ifdef globo
    if (ip_s.eq.ip_null) jstart = 1
    if (ip_n.eq.ip_null) jend = nlat
#endif

    iprocs = nprocsx*nprocsy

!----------------------
!  u-wind on t-points
!----------------------

#ifdef globo

#ifdef mpi
      if (nprocsx.eq.1) then
#endif
    u(1,:,:) = u(nlonm1,:,:)
    u(nlon,:,:) = u(2,:,:)
#ifdef mpi
      else
        call u_ghost (u(nlonm1,:,:), ip_e, u(1   ,:,:), ip_w, nlat*nlev)
        call u_ghost (u(2     ,:,:), ip_w, u(nlon,:,:), ip_e, nlat*nlev)
      endif
#endif

    do 2 jklev = 1, nlev
    do jlat = 2, nlatm1
    do jlon = 3, nlonm1
    ut(jlon,jlat,jklev) =.5625*(u(jlon  ,jlat,jklev)+u(jlon-1,jlat,jklev)) &
                        -.0625*(u(jlon+1,jlat,jklev)+u(jlon-2,jlat,jklev))
    enddo
    enddo

#ifdef mpi
      if (nprocsx.eq.1) then
#endif
    zut(:) = u(nlon-2,:,jklev)
#ifdef mpi
      else
        call u_ghost (u(nlon-2,:,jklev), ip_e, zut, ip_w, nlat)
      endif
#endif
    do jlat = 2, nlatm1
    ut(2,jlat,jklev) = .5625*(u(2,jlat,jklev)+u(1,jlat,jklev)) - .0625*(u(3,jlat,jklev)+zut(jlat))
    enddo
 2  continue

#else ! no globo

#ifdef mpi
      call u_ghost (u(2     ,:,:), ip_w, u(nlon,:,:), ip_e, nlat*nlev)
      call u_ghost (u(nlonm1,:,:), ip_e, u(1   ,:,:), ip_w, nlat*nlev)
#endif

    do jlon = 3, nlonm1
    ut(jlon,:,:) = .5625*(u(jlon,:,:)+u(jlon-1,:,:)) -.0625*(u(jlon+1,:,:)+u(jlon-2,:,:))
    enddo

#ifdef mpi
      call u_ghost (u(nlon-2,:,:), ip_e, zutb, ip_w, nlat*nlev)
      if (ip_w.eq.ip_null) then
#endif
    ut(1,:,:) =  u(1,:,:)
    ut(2,:,:) = .5*(u(1,:,:)+u(2,:,:))
#ifdef mpi
      else
        ut(2,:,:) = .5625*(u(2,:,:)+u(1,:,:)) - .0625*(u(3,:,:)+zutb(:,:))
      endif
      if (ip_e.eq.ip_null) then
#endif
    ut(nlon,:,:) = .5*(u(nlon,:,:)+u(nlonm1,:,:))
#ifdef mpi
      endif
#endif

#endif

!----------------------
!  v-wind on t-points
!----------------------

#ifdef globo

#ifdef mpi
      call u_ghost (v(2:nlonm1,nlatm1,:), ip_n, v  (2:nlonm1,1   ,:), ip_s, (nlon-2)*nlev)
      call u_ghost (v(2:nlonm1,2     ,:), ip_s, v  (2:nlonm1,nlat,:), ip_n, (nlon-2)*nlev)
      call u_ghost (v(2:nlonm1,3     ,:), ip_s, zvt(2:nlonm1,     :), ip_n, (nlon-2)*nlev)

      if (ip_s.eq.ip_null) then   ! row of processors including south pole

        if (nprocsx.eq.1) then
#endif

      do jlon = 2, nlonm1
      jlon1 = jlon + (nlon-2)/2
      if (jlon1.gt.nlon) jlon1 = jlon1-nlon+2
      v(jlon,1,:) = -v(jlon1,2,:)
      enddo
#ifdef mpi
        else
          call u_ghost (-v(2:nlonm1,2,:), ip_oppo, v(2:nlonm1,1,:), ip_oppo, (nlon-2)*nlev)
        endif
#endif

! Find components sin (zvi(1)) and cos (zvi(2)) of the 1st wavenumber of v at latitude
! closest to the South Pole (v(2))

    zvi = 0.
    do jklev = 1, nlev
    do jlon = 2, nlonm1
    zvi(1,jklev) = zvi(1,jklev) + v(jlon,2,jklev)*snt(jlon,1)
    zvi(2,jklev) = zvi(2,jklev) + v(jlon,2,jklev)*cst(jlon,1)
    enddo
    enddo
    zvi = zvi*dlon/180.

! Sum of projections sin and cos on all processes at the South Pole

#ifdef mpi
      if (myid.gt.0) then
        call mpi_send (zvi, 2*nlev, mpi_real, 0, tag1, comm, ierr)
        call mpi_recv (zvi, 2*nlev, mpi_real, 0, tag2, comm, status, ierr)
      else
        do jpr = nprocsy, iprocs-1, nprocsy
          call mpi_recv (zvo, 2*nlev, mpi_real, jpr, tag1, comm, status, ierr)
          zvi = zvi + zvo
        enddo
        do jpr = nprocsy, iprocs-1, nprocsy
          call mpi_send (zvi, 2*nlev, mpi_real, jpr, tag2, comm, ierr)
        enddo
      endif
#endif

! Analitical reconstruction of wind at the South Pole

    do jklev = 1, nlev
    do jlon = 2, nlonm1
    vt(jlon,1,jklev) = zvi(1,jklev)*snt(jlon,1) + zvi(2,jklev)*cst(jlon,1)
    ut(jlon,1,jklev) = zvi(1,jklev)*cst(jlon,1) - zvi(2,jklev)*snt(jlon,1)
    enddo
    enddo

#ifdef mpi
      endif ! row of processors including south pole
#endif

#ifdef mpi
      if (ip_n.eq.ip_null) then  !  row of processors including north pole
        if (nprocsx.eq.1) then
#endif
      do jlon = 2, nlonm1
      jlon1 = jlon + (nlon-2)/2
      if (jlon1.gt.nlon) jlon1 = jlon1-nlon+2
      zvt(jlon,:) = -v(jlon1,nlat,:)
      enddo
#ifdef mpi
        else
          call u_ghost (-v(2:nlonm1,nlat,:), ip_oppo, zvt(2:nlonm1,:), ip_oppo, (nlon-2)*nlev)
        endif
#endif
    zvi = 0.
    do jklev = 1, nlev
    do jlon = 2, nlonm1
    zvi(1,jklev) = zvi(1,jklev) + v(jlon,nlat,jklev)*snt(jlon,1)
    zvi(2,jklev) = zvi(2,jklev) + v(jlon,nlat,jklev)*cst(jlon,1)
    enddo
    enddo
    zvi = zvi*dlon/180.
#ifdef mpi
      if (myid.gt.nprocsy-1) then
        call mpi_send (zvi, 2*nlev, mpi_real, nprocsy-1, tag1, comm, ierr)
        call mpi_recv (zvi, 2*nlev, mpi_real, nprocsy-1, tag2, comm, status, ierr)
      else
        do jpr = 2*nprocsy-1, iprocs-1, nprocsy
          call mpi_recv (zvo, 2*nlev, mpi_real, jpr, tag1, comm, status, ierr)
          zvi = zvi + zvo
        enddo
        do jpr = 2*nprocsy-1, iprocs-1, nprocsy
          call mpi_send (zvi, 2*nlev, mpi_real, jpr, tag2, comm, ierr)
        enddo
      endif
#endif
    do jklev = 1, nlev
    do jlon = 2, nlonm1
    vt(jlon,nlat,jklev) = zvi(1,jklev)*snt(jlon,1) + zvi(2,jklev)*cst(jlon,1)
    ut(jlon,nlat,jklev) =-zvi(1,jklev)*cst(jlon,1) + zvi(2,jklev)*snt(jlon,1)
    enddo
    enddo
#ifdef mpi
      endif
#endif

    do jlat = 2, nlat-2
    vt(2:nlonm1,jlat,:) = .5625*(v(2:nlonm1,jlat+1,:)+v(2:nlonm1,jlat  ,:)) &
                         -.0625*(v(2:nlonm1,jlat+2,:)+v(2:nlonm1,jlat-1,:))
    enddo
    vt(2:nlonm1,nlatm1,:) = .5625*(v(2:nlonm1,nlat,:)+v(2:nlonm1,nlatm1,:)) &
                           -.0625*(zvt(2:nlonm1,:)+v(2:nlonm1,nlat-2,:))

#else ! no globo

#ifdef mpi
      call u_ghost (v(:,nlatm1,:), ip_n, v(:,1   ,:), ip_s, nlon*nlev)
      call u_ghost (v(:,2     ,:), ip_s, v(:,nlat,:), ip_n, nlon*nlev)
      call u_ghost (v(:,3     ,:), ip_s, zvt        , ip_n, nlon*nlev)
#endif

    do jlat = 2, nlat-2
    vt(:,jlat,:) = .5625*(v(:,jlat+1,:)+v(:,jlat,:))-.0625*(v(:,jlat+2,:)+v(:,jlat-1,:))
    enddo

#ifdef mpi
      if (ip_n.eq.ip_null) then
#endif
    vt(:,nlat  ,:) = v(:,nlat,:)
    vt(:,nlatm1,:) = .5*(v(:,nlat,:)+v(:,nlatm1,:))
#ifdef mpi
      else
        vt(:,nlatm1,:) = .5625*(v(:,nlat,:)+v(:,nlatm1,:))-.0625*(zvt(:,:)+v(:,nlat-2,:))
      endif
      if (ip_s.eq.ip_null) then
#endif
    vt(:,1,:) = .5*(v(:,1,:)+v(:,2,:))
#ifdef mpi
      endif
#endif

#endif ! globo

!---------------------
!  Vertical advections
!---------------------

    wfw1(:,1)=0.
    wfw2(:,1)=0.
    wfw1(:,nlevp1)=0.
    wfw2(:,nlevp1)=0.

    do 10 jlat = jstart, jend

    do jklev = 2, nlev
    zcost = .5*dt/dsig(jklev-1)    ! split of vert. adv.
    do jlon = 2, nlonm1
    zamu = sigdot(jlon,jlat,jklev)*zcost
      if (zamu.ge.0.) then
      is=1
      j1=jklev-1
      j1m1=j1-1
      if (j1m1.lt.1) j1m1=1
      else
      is=-1
      j1=jklev+1
      j1m1=j1-1
      if (j1.gt.nlev) j1=nlev
      endif
!    denr1 = 1./(ut(jlon,jlat,jklev)-ut(jlon,jlat,jklev-1))
!    denr2 = 1./(vt(jlon,jlat,jklev)-vt(jlon,jlat,jklev-1))
    denr1 = denomin(ut(jlon,jlat,jklev), ut(jlon,jlat,jklev-1))
    denr2 = denomin(vt(jlon,jlat,jklev), vt(jlon,jlat,jklev-1))
    r1 = (ut(jlon,jlat,j1)-ut(jlon,jlat,j1m1))*denr1
    r2 = (vt(jlon,jlat,j1)-vt(jlon,jlat,j1m1))*denr2
    b1 = max(0., min(2., max(r1, min(2.*r1,1.))))
    b2 = max(0., min(2., max(r2, min(2.*r2,1.))))
    zphi1 = is+zamu*b1 -is*b1
    zphi2 = is+zamu*b2 -is*b2
    wfw1(jlon,jklev)=.5*sigdot(jlon,jlat,jklev)*((1.+zphi1)*ut(jlon,jlat,jklev-1)+(1.-zphi1)*ut(jlon,jlat,jklev))
    wfw2(jlon,jklev)=.5*sigdot(jlon,jlat,jklev)*((1.+zphi2)*vt(jlon,jlat,jklev-1)+(1.-zphi2)*vt(jlon,jlat,jklev))
    enddo
    enddo

    do jklev = 1, nlev
    zcost = .5*dt/dsig(jklev)   ! split of vert. adv.
    do jlon = 2, nlonm1
    zdiv = (sigdot(jlon,jlat,jklev+1)-sigdot(jlon,jlat,jklev))*zcost
    wz1(jlon,jlat,jklev)=ut(jlon,jlat,jklev)+(wfw1(jlon,jklev)-wfw1(jlon,jklev+1))*zcost+ut(jlon,jlat,jklev)*zdiv
    wz2(jlon,jlat,jklev)=vt(jlon,jlat,jklev)+(wfw2(jlon,jklev)-wfw2(jlon,jklev+1))*zcost+vt(jlon,jlat,jklev)*zdiv
    enddo
    enddo

    do jklev = 2, nlev
    zcost = .5*dt/dsig(jklev-1)
    do jlon = 2, nlonm1
    zamu = sigdot(jlon,jlat,jklev)*zcost
      if (zamu.ge.0.) then
      is=1
      j1=jklev-1
      j1m1=j1-1
      if (j1m1.lt.1) j1m1=1
      else
      is=-1
      j1=jklev+1
      j1m1=j1-1
      if (j1.gt.nlev) j1=nlev
      endif
!    denr1 = 1./(wz1(jlon,jlat,jklev)-wz1(jlon,jlat,jklev-1))
!    denr2 = 1./(wz2(jlon,jlat,jklev)-wz2(jlon,jlat,jklev-1))
    denr1 = denomin(wz1(jlon,jlat,jklev), wz1(jlon,jlat,jklev-1))
    denr2 = denomin(wz2(jlon,jlat,jklev), wz2(jlon,jlat,jklev-1))
    r1 = (wz1(jlon,jlat,j1)-wz1(jlon,jlat,j1m1))*denr1
    r2 = (wz2(jlon,jlat,j1)-wz2(jlon,jlat,j1m1))*denr2
    b1 = max(0., min(2., max(r1, min(2.*r1,1.))))
    b2 = max(0., min(2., max(r2, min(2.*r2,1.))))
    zphi1 = is+zamu*b1 -is*b1
    zphi2 = is+zamu*b2 -is*b2
    wfw1(jlon,jklev)=.5*sigdot(jlon,jlat,jklev)*((1.+zphi1)*wz1(jlon,jlat,jklev-1)+(1.-zphi1)*wz1(jlon,jlat,jklev))
    wfw2(jlon,jklev)=.5*sigdot(jlon,jlat,jklev)*((1.+zphi2)*wz2(jlon,jlat,jklev-1)+(1.-zphi2)*wz2(jlon,jlat,jklev))
    enddo
    enddo

    do jklev = 1, nlev
    zcost = .5*dt/dsig(jklev)
    do jlon = 2, nlonm1
    zdiv = (sigdot(jlon,jlat,jklev+1)-sigdot(jlon,jlat,jklev))*zcost
    wz1(jlon,jlat,jklev)=wz1(jlon,jlat,jklev)+(wfw1(jlon,jklev)-wfw1(jlon,jklev+1))*zcost+ut(jlon,jlat,jklev)*zdiv
    wz2(jlon,jlat,jklev)=wz2(jlon,jlat,jklev)+(wfw2(jlon,jklev)-wfw2(jlon,jklev+1))*zcost+vt(jlon,jlat,jklev)*zdiv
    enddo
    enddo

 10 continue

!  Define 2 y-ghostlines of vertically advected variables

#ifdef globo

!  extra ghostline at poles (wind components must change sign)

#ifdef mpi
      if (ip_s.eq.ip_null) then
        if (nprocsx.eq.1) then
#endif
    do jlon = 2, nlonm1
    jlon1 = jlon + (nlon-2)/2
    if (jlon1.gt.nlon) jlon1 = jlon1-nlon+2
    wz1(jlon,0,:) = -wz1(jlon1,2,:)
    wz2(jlon,0,:) = -wz2(jlon1,2,:)
    enddo
#ifdef mpi
        else
          call u_ghost (-wz1(2:nlonm1,2,:), ip_oppo, wz1(2:nlonm1,0,:), ip_oppo, (nlon-2)*nlev)
          call u_ghost (-wz2(2:nlonm1,2,:), ip_oppo, wz2(2:nlonm1,0,:), ip_oppo, (nlon-2)*nlev)
        endif
      endif
#endif

#ifdef mpi
      if (ip_n.eq.ip_null) then
        if (nprocsx.eq.1) then
#endif
    do jlon = 2, nlonm1
    jlon1 = jlon + (nlon-2)/2
    if (jlon1.gt.nlon) jlon1 = jlon1-nlon+2
    wz1(jlon,nlat+1,:) = -wz1(jlon1,nlatm1,:)
    wz2(jlon,nlat+1,:) = -wz2(jlon1,nlatm1,:)
    enddo
#ifdef mpi
        else
          call u_ghost (-wz1(2:nlonm1,nlatm1,:), ip_oppo, wz1(2:nlonm1,nlat+1,:), ip_oppo, (nlon-2)*nlev)
          call u_ghost (-wz2(2:nlonm1,nlatm1,:), ip_oppo, wz2(2:nlonm1,nlat+1,:), ip_oppo, (nlon-2)*nlev)
        endif
      endif
#endif

#else ! no globo

#ifdef mpi
      if (ip_s.eq.ip_null) then
#endif
    wz1(2:nlonm1,0,:) = wz1(2:nlonm1,2,:)
    wz2(2:nlonm1,0,:) = wz2(2:nlonm1,2,:)
    wz1(2:nlonm1,1,:) = wz1(2:nlonm1,2,:)
    wz2(2:nlonm1,1,:) = wz2(2:nlonm1,2,:)
#ifdef mpi
      endif
      if (ip_n.eq.ip_null) then
#endif
    wz1(2:nlonm1,nlat  ,:) = wz1(2:nlonm1,nlatm1,:)
    wz2(2:nlonm1,nlat  ,:) = wz2(2:nlonm1,nlatm1,:)
    wz1(2:nlonm1,nlat+1,:) = wz1(2:nlonm1,nlatm1,:)
    wz2(2:nlonm1,nlat+1,:) = wz2(2:nlonm1,nlatm1,:)
#ifdef mpi
      endif
#endif

#endif ! globo

#ifdef mpi
      call u_ghost (wz1(2:nlonm1,2     :3     ,:), ip_s, wz1(2:nlonm1,nlat  :nlat+1,:), ip_n, 2*(nlon-2)*nlev)
      call u_ghost (wz1(2:nlonm1,nlat-2:nlatm1,:), ip_n, wz1(2:nlonm1,0     :1     ,:), ip_s, 2*(nlon-2)*nlev)
      call u_ghost (wz2(2:nlonm1,2     :3     ,:), ip_s, wz2(2:nlonm1,nlat  :nlat+1,:), ip_n, 2*(nlon-2)*nlev)
      call u_ghost (wz2(2:nlonm1,nlat-2:nlatm1,:), ip_n, wz2(2:nlonm1,0     :1     ,:), ip_s, 2*(nlon-2)*nlev)
#endif

!-----------------------
!  Horizontal advections
!-----------------------

    zcosty = dt/dy
    do 20 jklev = 1, nlev

!  Meridional advection (except poles)

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
!    denr1 = 1./(wz1(jlon,jlat,jklev)-wz1(jlon,jlat-1,jklev))
!    denr2 = 1./(wz2(jlon,jlat,jklev)-wz2(jlon,jlat-1,jklev))
    denr1 = denomin(wz1(jlon,jlat,jklev), wz1(jlon,jlat-1,jklev))
    denr2 = denomin(wz2(jlon,jlat,jklev), wz2(jlon,jlat-1,jklev))
    r1 = (wz1(jlon,j1,jklev)-wz1(jlon,j1-1,jklev))*denr1
    r2 = (wz2(jlon,j1,jklev)-wz2(jlon,j1-1,jklev))*denr2
    b1 = max(0., min(2., max(r1, min(2.*r1,1.))))
    b2 = max(0., min(2., max(r2, min(2.*r2,1.))))
    zphi1 = is+zamu*b1 -is*b1
    zphi2 = is+zamu*b2 -is*b2
    zpby1(jlon,jlat) = .5*((1.+zphi1)*wz1(jlon,jlat-1,jklev)+(1.-zphi1)*wz1(jlon,jlat,jklev))
    zpby2(jlon,jlat) = .5*((1.+zphi2)*wz2(jlon,jlat-1,jklev)+(1.-zphi2)*wz2(jlon,jlat,jklev))
    enddo
 15 continue

    do jlat = 2, nlatm1
    zhxvtn = hxv(jlat+1)/(hxt(jlat)*dy)
    zhxvt  = hxv(jlat  )/(hxt(jlat)*dy)
    do jlon = 2, nlonm1
    p01(jlon,jlat,jklev) = wz1(jlon,jlat,jklev) - dt*(v(jlon,jlat+1,jklev)*(zpby1(jlon,jlat+1)-             &
              ut(jlon,jlat,jklev))*zhxvtn -v(jlon,jlat,jklev)*(zpby1(jlon,jlat)-ut(jlon,jlat,jklev))*zhxvt)
    p02(jlon,jlat,jklev) = wz2(jlon,jlat,jklev) - dt*(v(jlon,jlat+1,jklev)*(zpby2(jlon,jlat+1)-             &
              vt(jlon,jlat,jklev))*zhxvtn -v(jlon,jlat,jklev)*(zpby2(jlon,jlat)-vt(jlon,jlat,jklev))*zhxvt)
    enddo
    enddo

#ifdef globo
#ifdef mpi
      if (ip_s.eq.ip_null) then !  South Pole
#endif
    p02(2:nlonm1,1,jklev) = wz2(2:nlonm1,1,jklev) -.25*dt/dy*(wz2(2:nlonm1,2,jklev)**2-wz2(2:nlonm1,0,jklev)**2)
#ifdef mpi
      endif
      if (ip_n.eq.ip_null) then !  North Pole
#endif
    p02(2:nlonm1,nlat,jklev) = wz2(2:nlonm1,nlat  ,jklev) - .25*dt/dy*                        &
                              (wz2(2:nlonm1,nlat+1,jklev)**2 - wz2(2:nlonm1,nlatm1,jklev)**2)
#ifdef mpi
      endif
#endif
#endif

 20 continue

#ifndef globo
#ifdef mpi
      if (ip_w.eq.ip_null) then
#endif
    p01(0,:,:) = p01(2,:,:)
    p01(1,:,:) = p01(2,:,:)
    p02(0,:,:) = p02(2,:,:)
    p02(1,:,:) = p02(2,:,:)
#ifdef mpi
      endif
      if (ip_e.eq.ip_null) then
#endif
    p01(nlon  ,:,:) = p01(nlonm1,:,:)
    p01(nlon+1,:,:) = p01(nlonm1,:,:)
    p02(nlon  ,:,:) = p02(nlonm1,:,:)
    p02(nlon+1,:,:) = p02(nlonm1,:,:)
#ifdef mpi
      endif
#endif
#endif

!  Zonal advection and centrifugal acceleration (exept poles)

    do 16 jlat = 2, nlatm1
    zcostx = dt/( hxt(jlat)*dx*float(nsweep(jlat)) )
    zcostg = .25*dt*hst(jlat)/( float(nsweep(jlat))*a*hxt(jlat) )

    do jsweep = 1, nsweep(jlat)

#ifdef mpi
      if (nprocsx.eq.1) then
#endif

#ifdef globo
    p01(0:1,jlat,1:nlev) = p01(nlon-2:nlonm1,jlat,1:nlev)
    p02(0:1,jlat,1:nlev) = p02(nlon-2:nlonm1,jlat,1:nlev)
    p01(nlon:nlon+1,jlat,1:nlev) = p01(2:3,jlat,1:nlev)
    p02(nlon:nlon+1,jlat,1:nlev) = p02(2:3,jlat,1:nlev)
#endif

#ifdef mpi
      else
        call u_ghost (p01(nlon-2:nlonm1,jlat,:), ip_e, p01(0:1        ,jlat,:), ip_w, 2*nlev)
        call u_ghost (p01(2:3          ,jlat,:), ip_w, p01(nlon:nlon+1,jlat,:), ip_e, 2*nlev)
        call u_ghost (p02(nlon-2:nlonm1,jlat,:), ip_e, p02(0:1        ,jlat,:), ip_w, 2*nlev)
        call u_ghost (p02(2:3          ,jlat,:), ip_w, p02(nlon:nlon+1,jlat,:), ip_e, 2*nlev)
      endif
#endif

    do 21 jklev = 1, nlev
    do jlon = 2, nlon
    zamu = u(jlon-1,jlat,jklev)*zcostx
      if(zamu.ge.0.) then
      is=1
      j1=jlon-1
      else
      is=-1
      j1=jlon+1
      endif
!    denr1 = 1./(p01(jlon,jlat,jklev)-p01(jlon-1,jlat,jklev))
!    denr2 = 1./(p02(jlon,jlat,jklev)-p02(jlon-1,jlat,jklev))
    denr1 = denomin(p01(jlon,jlat,jklev), p01(jlon-1,jlat,jklev))
    denr2 = denomin(p02(jlon,jlat,jklev), p02(jlon-1,jlat,jklev))
    r1 = (p01(j1,jlat,jklev)-p01(j1-1,jlat,jklev))*denr1
    r2 = (p02(j1,jlat,jklev)-p02(j1-1,jlat,jklev))*denr2
    b1 = max(0., min(2., max(r1, min(2.*r1,1.))))
    b2 = max(0., min(2., max(r2, min(2.*r2,1.))))
    zphi1 = is+zamu*b1 -is*b1
    zphi2 = is+zamu*b2 -is*b2
    zpbw1(jlon) = .5*((1.+zphi1)*p01(jlon-1,jlat,jklev)+(1.-zphi1)*p01(jlon,jlat,jklev))
    zpbw2(jlon) = .5*((1.+zphi2)*p02(jlon-1,jlat,jklev)+(1.-zphi2)*p02(jlon,jlat,jklev))
    enddo

    do jlon = 2, nlonm1
    p01(jlon,jlat,jklev) = p01(jlon,jlat,jklev) - zcostx*( u(jlon,jlat,jklev)*(zpbw1(jlon+1) -             &
                             ut(jlon,jlat,jklev))-u(jlon-1,jlat,jklev)*(zpbw1(jlon)-ut(jlon,jlat,jklev)) ) &
                             +zcostg*(zpbw1(jlon+1)+zpbw1(jlon))*(zpbw2(jlon+1)+zpbw2(jlon))
    p02(jlon,jlat,jklev) = p02(jlon,jlat,jklev) - zcostx*( u(jlon,jlat,jklev)*(zpbw2(jlon+1) -             &
                             vt(jlon,jlat,jklev))-u(jlon-1,jlat,jklev)*(zpbw2(jlon)-vt(jlon,jlat,jklev)) ) &
                             -zcostg*(zpbw1(jlon+1)+zpbw1(jlon))**2
    enddo
 21 continue
    enddo  !! end sweeps
 16 continue

!---------------------------------------
!  back to wind points: u (fourth order)
!---------------------------------------

#ifndef globo
#ifdef mpi
      if (ip_w.eq.ip_null) then
#endif
    p01(1,:,:) = p01(2,:,:)
    p02(1,:,:) = p02(2,:,:)
#ifdef mpi
      endif
      if (ip_e.eq.ip_null) then
#endif
    p01(nlon  ,:,:) = p01(nlonm1,:,:)
    p01(nlon+1,:,:) = p01(nlonm1,:,:)
    p02(nlon  ,:,:) = p02(nlonm1,:,:)
    p02(nlon+1,:,:) = p02(nlonm1,:,:)
#ifdef mpi
      endif
      if (ip_s.eq.ip_null) p02(2:nlonm1,1   ,:) = vt(2:nlonm1,1   ,:)
      if (ip_n.eq.ip_null) p02(2:nlonm1,nlat,:) = vt(2:nlonm1,nlat,:)
#endif
#endif

#ifdef mpi
      if (nprocsx.eq.1) then
#endif

#ifdef globo
    p01(1,2:nlatm1,:) = p01(nlonm1,2:nlatm1,:)
    p01(nlon:nlon+1,2:nlatm1,:) = p01(2:3,2:nlatm1,:)
#endif

#ifdef mpi
      else
        call u_ghost (p01(nlonm1,2:nlatm1,:), ip_e, p01(1          ,2:nlatm1,:), ip_w,   (nlat-2)*nlev)
        call u_ghost (p01(2:3   ,2:nlatm1,:), ip_w, p01(nlon:nlon+1,2:nlatm1,:), ip_e, 2*(nlat-2)*nlev)
      endif
#endif

    do jklev = 1, nlev
    do jlat = 2, nlatm1
    do jlon = 2, nlonm1
    u(jlon,jlat,jklev) = .5625*(p01(jlon,  jlat,jklev)+p01(jlon+1,jlat,jklev)) &
                        -.0625*(p01(jlon-1,jlat,jklev)+p01(jlon+2,jlat,jklev))
    enddo
    enddo
    enddo

!---------------------------------------
!  back to wind points: v (fourth order)
!---------------------------------------

#ifdef mpi
      call u_ghost (p02(2:nlonm1,nlatm1,:), ip_n, p02(2:nlonm1,1   ,:), ip_s, (nlon-2)*nlev)
      call u_ghost (p02(2:nlonm1,2     ,:), ip_s, p02(2:nlonm1,nlat,:), ip_n, (nlon-2)*nlev)
      call u_ghost (p02(2:nlonm1,nlat-2,:), ip_n, zvt(2:nlonm1     ,:), ip_s, (nlon-2)*nlev)
#endif

#ifdef globo

#ifdef mpi
      if (ip_s.eq.ip_null) then
#endif
    call slofou1 (p02(1:nlon,1,1:nlev), 1., .true.)
#ifdef mpi
        if (nprocsx.eq.1) then
#endif
    do jklev = 1, nlev
    do jlon = 2, nlonm1
    jlon1 = jlon + (nlon-2)/2
    if (jlon1.gt.nlon) jlon1 = jlon1 - nlon + 2
    zvt(jlon,jklev) = -p02(jlon1,2,jklev)
    enddo
    enddo
#ifdef mpi
        else
          call u_ghost (-p02(2:nlonm1,2,:), ip_oppo, zvt(2:nlonm1,:), ip_oppo, (nlon-2)*nlev)
        endif
      endif
#endif

#ifdef mpi
      if (ip_n.eq.ip_null) then
#endif
    call slofou1 (p02(1:nlon,nlat,1:nlev), 1., .true.)
#ifdef mpi
        if (nprocsx.eq.1) then
#endif
    do jklev = 1, nlev
    do jlon = 2, nlonm1
    jlon1 = jlon + (nlon-2)/2
    if (jlon1.gt.nlon) jlon1 = jlon1-nlon+2
    v(jlon,nlat,jklev) = -p02(jlon1,nlatm1,jklev)
    enddo
    enddo
#ifdef mpi
        else
          call u_ghost (-p02(2:nlonm1,nlatm1,:), ip_oppo, v(2:nlonm1,nlat,:), ip_oppo, (nlon-2)*nlev)
        endif
#endif
    v(2:nlonm1,nlat,:) = .5625*(p02(2:nlonm1,nlat,:)+p02(2:nlonm1,nlatm1,:)) &
                        -.0625*(v  (2:nlonm1,nlat,:)+p02(2:nlonm1,nlat-2,:))
#ifdef mpi
      endif
#endif

#endif ! globo

    do jlat = 3, nlatm1
    v(2:nlonm1,jlat,:) = .5625*(p02(2:nlonm1,jlat  ,:)+p02(2:nlonm1,jlat-1,:)) &
                        -.0625*(p02(2:nlonm1,jlat+1,:)+p02(2:nlonm1,jlat-2,:))
    enddo

#ifdef globo
    v(2:nlonm1,2   ,:) = .5625*(p02(2:nlonm1,2,:)+p02(2:nlonm1,1,:)) &
                        -.0625*(p02(2:nlonm1,3,:)+zvt(2:nlonm1  ,:))
#else ! no globo
#ifdef mpi
      if (ip_s.eq.ip_null) then
#endif
    v(2:nlonm1,2,:) = .5*(p02(2:nlonm1,1,:)+p02(2:nlonm1,2,:))
#ifdef mpi
      else
      v(2:nlonm1,2,:) = .5625*(p02(2:nlonm1,2,:)+p02(2:nlonm1,1,:)) &
                       -.0625*(p02(2:nlonm1,3,:)+zvt(2:nlonm1,  :))
      endif
#endif ! globo
#endif

    return
    end subroutine wafuv
!##################################################################################################################
    subroutine wafone (var, dt)

! WAF advection of scalar variables

    use mod_model, only : nlon, gnlon, nlat, nlev, nlevp1, nlonm1, nlatm1, u, v, sigdot, hxt, hxv, dsig, pi, &
            dx, dy, a, dlat, nprocsy, nprocsx, nsweep, myid, ip_e, ip_n, ip_s, ip_w, ip_null, ip_oppo, cpole, denomin

    implicit none
    real(4), dimension(nlon,nlat,nlev)     :: var, dvarx
    real(4), dimension(nlon)               :: zpbw
    real(4), dimension(nlon,nlevp1)        :: wfw
    real(4), dimension(nlon,nlat)          :: zpby
    real(4), dimension(nlon,0:nlat+1,nlev) :: wz
    real(4), dimension(0:nlon+1,nlat,nlev) :: p0
    integer :: jlon, jlat, jklev, jlon1, j1, j1m1, is, jsweep, jstart, jend, iprocs, jpr, &
 ierr, comm, tag1=1, tag2=2
    real(4) dt, zcost, zamu, denr, r, b, zphi, zp0i, zp0o
    real(4) zdiv, zcostx, zcosty, zhxvt, zhxvtn

#ifdef mpi
      include 'mpif.h'
      integer                              :: status(mpi_status_size)

      comm = mpi_comm_world
#endif

    jstart = 2
    jend = nlatm1

#ifdef globo
    if (ip_s.eq.ip_null) jstart = 1
    if (ip_n.eq.ip_null) jend = nlat
#endif

    iprocs = nprocsx*nprocsy

!---------------------
!  Vertical advections
!---------------------

    wfw(:,1)=0.
    wfw(:,nlevp1)=0.

    do 10 jlat = jstart, jend

    do jklev = 2, nlev
    zcost = .5*dt/dsig(jklev-1)
    do jlon = 2, nlonm1
    zamu = sigdot(jlon,jlat,jklev)*zcost
      if (zamu.ge.0.) then
      is=1
      j1=jklev-1
      j1m1=j1-1
      if (j1m1.lt.1) j1m1=1
      else
      is=-1
      j1=jklev+1
      j1m1=j1-1
      if (j1.gt.nlev) j1=nlev
      endif
!    denr = 1./(var(jlon,jlat,jklev)-var(jlon,jlat,jklev-1))
    denr = denomin(var(jlon,jlat,jklev), var(jlon,jlat,jklev-1))
    r = (var(jlon,jlat,j1)-var(jlon,jlat,j1m1))*denr
    b = max(0., min(2., max(r, min(2.*r,1.))))
    zphi = is+zamu*b -is*b
    wfw(jlon,jklev) = .5*sigdot(jlon,jlat,jklev)*((1.+zphi)*var(jlon,jlat,jklev-1)+(1.-zphi)*var(jlon,jlat,jklev))
    enddo
    enddo

    do jklev = 1, nlev
    zcost = .5*dt/dsig(jklev)   ! split of vert. adv.
    do jlon = 2, nlonm1
    zdiv = (sigdot(jlon,jlat,jklev+1)-sigdot(jlon,jlat,jklev))*zcost
    wz(jlon,jlat,jklev) = var(jlon,jlat,jklev)+(wfw(jlon,jklev)-wfw(jlon,jklev+1))*zcost + var(jlon,jlat,jklev)*zdiv
    enddo
    enddo

    do jklev = 2, nlev
    zcost = .5*dt/dsig(jklev-1)
    do jlon = 2, nlonm1
    zamu = sigdot(jlon,jlat,jklev)*zcost
      if (zamu.ge.0.) then
      is=1
      j1=jklev-1
      j1m1=j1-1
      if (j1m1.lt.1) j1m1=1
      else
      is=-1
      j1=jklev+1
      j1m1=j1-1
      if (j1.gt.nlev) j1=nlev
      endif
!    denr = 1./(wz(jlon,jlat,jklev)-wz(jlon,jlat,jklev-1))
    denr = denomin (wz(jlon,jlat,jklev), wz(jlon,jlat,jklev-1))
    r = (wz(jlon,jlat,j1)-wz(jlon,jlat,j1m1))*denr
    b = max(0., min(2., max(r, min(2.*r,1.))))
    zphi = is+zamu*b -is*b
    wfw(jlon,jklev) = .5*sigdot(jlon,jlat,jklev)*((1.+zphi)*wz(jlon,jlat,jklev-1)+(1.-zphi)*wz(jlon,jlat,jklev))
    enddo
    enddo

    do jklev = 1, nlev
    zcost = .5*dt/dsig(jklev)
    do jlon = 2, nlonm1
    zdiv = (sigdot(jlon,jlat,jklev+1)-sigdot(jlon,jlat,jklev))*zcost
    wz(jlon,jlat,jklev) = wz(jlon,jlat,jklev)+(wfw(jlon,jklev)-wfw(jlon,jklev+1))*zcost + var(jlon,jlat,jklev)*zdiv
    enddo
    enddo

 10 continue

#ifndef globo
#ifdef mpi
      if (ip_s.eq.ip_null) then
#endif
    wz(2:nlonm1,0,:) = wz(2:nlonm1,2,:)
    wz(2:nlonm1,1,:) = wz(2:nlonm1,2,:)
#ifdef mpi
      endif
      if (ip_n.eq.ip_null) then
#endif
    wz(2:nlonm1,nlat  ,:) = wz(2:nlonm1,nlatm1,:)
    wz(2:nlonm1,nlat+1,:) = wz(2:nlonm1,nlatm1,:)
#ifdef mpi
      endif
#endif
#endif

#ifdef mpi
      call u_ghost (wz(2:nlonm1,2     :3     ,:), ip_s, wz(2:nlonm1,nlat  :nlat+1,:), ip_n, 2*(nlon-2)*nlev)
      call u_ghost (wz(2:nlonm1,nlat-2:nlatm1,:), ip_n, wz(2:nlonm1,0     :1     ,:), ip_s, 2*(nlon-2)*nlev)
#endif

#ifdef globo

!  extra ghostline at poles

#ifdef mpi
      if (ip_s.eq.ip_null) then
        if (nprocsx.eq.1) then
#endif
    do jlon = 2, nlonm1
    jlon1 = jlon + (nlon-2)/2
    if (jlon1.ge.nlon) jlon1 = jlon1-nlon+2
    wz(jlon,0,:) =  wz(jlon1,2,:)
    enddo
#ifdef mpi
        else
          call u_ghost ( wz(2:nlonm1,2,:), ip_oppo, wz(2:nlonm1,0,:), ip_oppo, (nlon-2)*nlev)
        endif
      endif
#endif

#ifdef mpi
      if (ip_n.eq.ip_null) then
        if (nprocsx.eq.1) then
#endif
    do jlon = 2, nlonm1
    jlon1 = jlon + (nlon-2)/2
    if (jlon1.ge.nlon) jlon1 = jlon1-nlon+2
    wz(jlon,nlat+1,:) =  wz(jlon1,nlatm1,:)
    enddo
#ifdef mpi
        else
          call u_ghost ( wz(2:nlonm1,nlatm1,:), ip_oppo, wz(2:nlonm1,nlat+1,:), ip_oppo, (nlon-2)*nlev)
        endif
      endif
#endif

#endif ! globo

!-----------------------
!  Horizontal advections
!-----------------------

    zcosty = dt/dy
    do 20 jklev = 1, nlev

!  Meridional advection (except poles)

    do 15 jlat = 2, nlat
    do jlon = 2, nlonm1
    zamu = v(jlon,jlat,jklev)*zcosty
      if(zamu.gt.0.) then
      is=1
      j1=jlat-1
      else
      is=-1
      j1=jlat+1
      endif
!    denr = 1./(wz(jlon,jlat,jklev)-wz(jlon,jlat-1,jklev))
    denr = denomin(wz(jlon,jlat,jklev), wz(jlon,jlat-1,jklev))
    r = (wz(jlon,j1,jklev)-wz(jlon,j1-1,jklev))*denr
    b = max(0., min(2., max(r, min(2.*r,1.))))
    zphi = is+zamu*b -is*b
    zpby(jlon,jlat)=.5*((1.+zphi)*wz(jlon,jlat-1,jklev)+(1.-zphi)*wz(jlon,jlat,jklev))
    enddo
 15 continue

    do jlat = 2, nlatm1
    zhxvtn = hxv(jlat+1)/(hxt(jlat)*dy)
    zhxvt  = hxv(jlat  )/(hxt(jlat)*dy)
    do jlon = 2, nlonm1
    p0(jlon,jlat,jklev) = wz(jlon,jlat,jklev) - dt*(v(jlon,jlat+1,jklev)*(zpby(jlon,jlat+1)-                &
            var(jlon,jlat,jklev))*zhxvtn - v(jlon,jlat,jklev)*(zpby(jlon,jlat)-var(jlon,jlat,jklev))*zhxvt)
    enddo
    enddo

#ifdef globo

#ifdef mpi
      if (ip_s.eq.ip_null) then !  South Pole
#endif
    zp0i = 0.
    do jlon = 2, nlonm1
    zp0i = zp0i + dt*cpole*v(jlon,2,jklev)*(var(jlon,1,jklev)-zpby(jlon,2))
    enddo
#ifdef mpi
      if (myid.ne.0) then
        call mpi_send (zp0i, 1, mpi_real, 0, tag1, comm, ierr)
        call mpi_recv (zp0i, 1, mpi_real, 0, tag2, comm, status, ierr)
      else
        do jpr = nprocsy, iprocs-1, nprocsy
          call mpi_recv (zp0o, 1, mpi_real, jpr, tag1, comm, status, ierr)
          zp0i = zp0i + zp0o
        enddo
        do jpr = nprocsy, iprocs-1, nprocsy
          call mpi_send (zp0i, 1, mpi_real, jpr, tag2, comm, ierr)
        enddo
      endif
#endif
    do jlon = 2, nlonm1
    var(jlon,1,jklev) = zp0i + wz(2,1,jklev)
    enddo
#ifdef mpi
      endif
#endif

#ifdef mpi
      if (ip_n.eq.ip_null) then !  North Pole
#endif
    zp0i = 0.
    do jlon = 2, nlonm1
    zp0i = zp0i - dt*cpole*v(jlon,nlat,jklev)*(var(jlon,nlat,jklev)-zpby(jlon,nlat))
    enddo
#ifdef mpi
        if (myid.ne.nprocsy-1) then
          call mpi_send (zp0i, 1, mpi_real, nprocsy-1, tag1, comm, ierr)
          call mpi_recv (zp0i, 1, mpi_real, nprocsy-1, tag2, comm, status, ierr)
        else
          do jpr = 2*nprocsy-1, iprocs-1, nprocsy
            call mpi_recv (zp0o, 1, mpi_real, jpr, tag1, comm, status, ierr)
            zp0i = zp0i + zp0o
          enddo
          do jpr = 2*nprocsy-1, iprocs-1, nprocsy
            call mpi_send (zp0i, 1, mpi_real, jpr, tag2, comm, ierr)
          enddo
        endif
#endif
    do jlon = 2, nlonm1
    var(jlon,nlat,jklev) = zp0i + wz(2,nlat,jklev)
    enddo
#ifdef mpi
      endif
#endif

#endif ! globo

 20 continue

!  Zonal advection (except poles)

#ifndef globo

#ifdef mpi
      if (ip_w.eq.ip_null) then
#endif
    p0(0,:,:) = p0(2,:,:)
    p0(1,:,:) = p0(2,:,:)
#ifdef mpi
      endif
      if (ip_e.eq.ip_null) then
#endif
    p0(nlon  ,:,:) = p0(nlonm1,:,:)
    p0(nlon+1,:,:) = p0(nlonm1,:,:)
#ifdef mpi
      endif
#endif

#endif ! no globo

!  Define 2 x-ghostlines of vertically and meridionally advected variables

#ifdef mpi
      if (nprocsx.eq.1) then
#endif

#ifdef globo
    p0(0:1        ,:,:) = p0(nlon-2:nlonm1,:,:)
    p0(nlon:nlon+1,:,:) = p0(2:3          ,:,:)
#endif

#ifdef mpi
      else
        call u_ghost (p0(nlon-2:nlonm1,:,:), ip_e, p0(0:1        ,:,:), ip_w, 2*nlat*nlev)
        call u_ghost (p0(2:3          ,:,:), ip_w, p0(nlon:nlon+1,:,:), ip_e, 2*nlat*nlev)
      endif
#endif

    do 21 jklev = 1, nlev
    do 16 jlat = 2, nlatm1
    zcostx = dt/(hxt(jlat)*dx)
    do jlon = 2, nlon
    zamu = u(jlon-1,jlat,jklev)*zcostx
      if (zamu.ge.0.) then
      is=1
      j1=jlon-1
      else
      is=-1
      j1=jlon+1
      endif
!    denr = 1./(p0(jlon,jlat,jklev)-p0(jlon-1,jlat,jklev))
    denr = denomin(p0(jlon,jlat,jklev), p0(jlon-1,jlat,jklev))
    r = (p0(j1,jlat,jklev)-p0(j1-1,jlat,jklev))*denr
    b = max(0., min(2., max(r, min(2.*r,1.))))
    zphi = is+zamu*b -is*b
    zpbw(jlon) = .5*((1.+zphi)*p0(jlon-1,jlat,jklev)+(1.-zphi)*p0(jlon,jlat,jklev))
    enddo

    do jlon = 2, nlonm1
    dvarx(jlon,jlat,jklev) = -zcostx*( u(jlon,jlat,jklev)*(zpbw(jlon+1)-var(jlon,jlat,jklev)) &
                                      -u(jlon-1,jlat,jklev)*(zpbw(jlon)-var(jlon,jlat,jklev)) )
    enddo
 16 continue
 21 continue

#ifdef globo
    call polavert (dvarx)
#endif

    do jlat = 2, nlatm1
    var(2:nlonm1,jlat,:) = max (0., p0(2:nlonm1,jlat,:)+dvarx(2:nlonm1,jlat,:))
    enddo

    return
    end subroutine wafone
!##################################################################################################################
    subroutine surface_layer (pstep, kstep)

! Computes parameters and variables of the turbulent exchanging in the surface layer
! Computes heat and specific humidity fluxes at the ground
! Fluxes are positive upward

    use mod_model, only: nlon, nlat, nlev, nlonm1, nlatm1, ps, u, v, t, q, tvirt, roscd,                    &
                         rd, cpd, g, ep, rdrcp, phig, tskin, qskin, sigint, pzer, alp0, dx, phi, fmask, fice, &
                         rgm, rgq, rgmd, sigalf, fsnow, t2, q2, u10, v10, ustar, tstar, rought, psih, rich, &
                         kturb_surf_m, kturb_surf_h, kturb_surf_q, hflux, qflux,                            &
                         n_std_lev_atm, n_std_lev_sl, std_lev_atm, t_std_lev, u_std_lev, v_std_lev, q_std_lev, &
                         nstep_sl_filter, myid, kturb_surf_h_mem, kturb_surf_q_mem

#ifdef chem
  use mod_bol2chem
#endif

    real(4), parameter              :: zak=.4, zgam=16., zaholt=1., zbholt=2./3., zcholt=5., zdholt=0.35
    real(4), dimension(n_std_lev_atm) :: psim_sl, psih_sl, psiq_sl
    integer :: nlev_sl, kstep
    real(4) :: ztetav_bottom, zlev_bottom
    real(4), dimension(nlon,nlat) :: kturb_surf_h_filter, kturb_surf_q_filter

!  Businger functions

    psium(zz1,zz2) = log( (1.+zz1)**2*(1.+zz1**2)/((1.+zz2)**2*(1.+zz2**2)) ) -2.*(atan(zz1)-atan(zz2))
    psiuh(zz1,zz2) = 2.*log((1.+zz1**2)/(1.+zz2**2))

!  Holtslag functions

    psism(zz1,zz2)=-zaholt*zz1-zbholt*(zz1-zcholt/zdholt)*exp(-zdholt*zz1)  &
                   +zaholt*zz2+zbholt*(zz2-zcholt/zdholt)*exp(-zdholt*zz2)
    psish(zz1,zz2)=-(1.+2./3.*zaholt*zz1)**1.5-zbholt*(zz1-zcholt/zdholt)*exp(-zdholt*zz1) &
                   +(1.+2./3.*zaholt*zz2)**1.5+zbholt*(zz2-zcholt/zdholt)*exp(-zdholt*zz2)

    do jlat = 1, nlat
    do jlon = 1, nlon
      do jklev = 1, n_std_lev_atm
        if (std_lev_atm(jklev) > (phi(jlon,jlat,nlev)-phig(jlon,jlat))/g) then
         n_std_lev_sl(jlon,jlat) = jklev-1
         exit
        endif
      enddo
    enddo
    enddo

!--------------------------------------------
!  loop in latitude to the end of the routine
!--------------------------------------------

    jstart = 2
    jend = nlatm1

#ifdef globo
    if (ip_s.eq.ip_null) jstart = 1
    if (ip_n.eq.ip_null) jend = nlat
#endif

    do jlat = jstart, jend
    jlatp1 = min (jlat+1, nlat) ! necessaria per globo

!  conversion factor from t to theta
!  definition of virtual theta
!  specific humidity at saturation
!  half-level height

!------------------------------------------------------------------------------------------
!  turbulent fluxes at the ground (positive upward) - surface layer - Monin Obukhov theory
!------------------------------------------------------------------------------------------

    do jlon = 2, nlonm1

    zzpp = pzer*sigint(nlev) - (pzer-ps(jlon,jlat))*sigalf(nlev)
    zconv = exp(rdrcp*(alp0-log(zzpp)))
    ztetav_bottom = tvirt(jlon,jlat,nlev)*zconv

    zza   = (phi(jlon,jlat,nlev)-phig(jlon,jlat))/g + rgm(jlon,jlat) !  sigma=1 is located at z=rgm
    zua   = .5*(u(jlon,jlat,nlev)+u(jlon-1,jlat,nlev))
    zva   = .5*(v(jlon,jlat,nlev)+v(jlon,jlatp1,nlev))
    zmod2 = zua**2 + zva**2 + .07
    zmod  = sqrt(zmod2)
    nlev_sl = n_std_lev_sl(jlon,jlat)

!  virtual potential temperature computed with skin temperature

    zconvg = exp(rdrcp*(alp0-log(ps(jlon,jlat))))*(1.+ep*qskin(jlon,jlat))
    ztevg = tskin(jlon,jlat)*zconvg

!  bulk Richardson number

    ztebar =.5*(ztevg+ztetav_bottom)
    beta   = g/ztebar
    zri    = zza*beta*(ztetav_bottom-ztevg)/zmod2

    if(fmask(jlon,jlat).ge..5.and.fice(jlon,jlat).lt.0.1) then

!  computation of Charnock roughness

      zchar = 5.e-4
      zcoch1 = .0185*zak**2/g*zmod2
      do jiter = 1, 5
      zcoch2 = zza/zchar
      zchar = zcoch1/log(zcoch2)**2
      enddo
      zrgm = zchar

! roughness lengths over sea interpolated logarithmically between 0.<Ri<0.25 with Large & Pond values

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

    else                     ! soil or sea ice > 0.1
      zrgm = rgm(jlon,jlat)
      zrgt = rgq(jlon,jlat)
      zrgm = zrgm*(1.-0.5*fsnow(jlon,jlat)*0.4/(0.4+zrgm))  ! reduction over snow
      zrgt = zrgt*(1.-0.5*fsnow(jlon,jlat)*0.4/(0.4+zrgt))  ! reduction over snow
      zrgq = zrgt
    endif

    rgmd(jlon,jlat) = zrgm

#ifdef chem
      bolchem_rgm(jlon,jlat)=zrgm
      bolchem_rgq(jlon,jlat)=zrgq
#endif

    zalzam = log(zza/zrgm)
    zalzat = log(zza/zrgt)
    zalzaq = log(zza/zrgq)

    if(zri.gt.0.) then         ! Holtslag functions

      zal = .2 + 4.*zri
      zpsim = psism(zal, zal*zrgm/zza)
      zpsih = psish(zal, zal*zrgt/zza)
      do jiter = 1, 10
      zalm = zal
      zal = zri*(zalzam-zpsim)**2/(zalzat-zpsih)   ! za/L from eq (6.47) of Businger
      error = abs(zal-zalm)/zal
      if (error.lt.1.e-2) go to 1
      zpsim = psism(zal, zal*zrgm/zza)
      zpsih = psish(zal, zal*zrgt/zza)
      enddo
 1    continue
      zpsiq = zpsih  ! because rgt=rgq in the stable case
      ustar(jlon,jlat) = zak*zmod/(zalzam-zpsim)
      tstar(jlon,jlat) = zak/(zalzat-zpsih)*(ztetav_bottom-ztevg)
      zqstar           = zak/(zalzaq-zpsiq)*(q(jlon,jlat,nlev)-qskin(jlon,jlat))
      zcdm = ustar(jlon,jlat)**2/zmod
      zcdt = ustar(jlon,jlat)*zak/(zalzat-zpsih)
      zcdq = zcdt
           zpsim10 = psism(zal*(10.+zrgm)/zza, zal*zrgm/zza)  ! to compute 10 m wind
           zpsih2  = psish(zal*(2. +zrgt)/zza, zal*zrgt/zza)  ! to compute 2 m temp.
           zpsiq2  = zpsih2    ! because rgt=rgq in the stable case
           do jklev = 1, nlev_sl
             psim_sl(jklev) = psism(zal*(std_lev_atm(jklev)+zrgm)/zza, zal*zrgm/zza)
             psih_sl(jklev) = psish(zal*(std_lev_atm(jklev)+zrgt)/zza, zal*zrgt/zza)
             psiq_sl(jklev) = psih_sl(jklev)         ! because rgt=rgq in the stable case
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
      tstar(jlon,jlat) = zak/(zalzat-zpsih)*(ztetav_bottom-ztevg)
      zqstar           = zak/(zalzaq-zpsiq)*(q(jlon,jlat,nlev)-qskin(jlon,jlat))
      zcdm = ustar(jlon,jlat)**2/zmod
      zcdt = ustar(jlon,jlat)*zak/(zalzat-zpsih)
      zcdq = ustar(jlon,jlat)*zak/(zalzaq-zpsiq)
           zz1 = (1.-zgam*zal*(10.+zrgm)/zza)**.25   ! for 10 m wind
           zpsim10 = psium(zz1, zx2)                 ! for 10 m wind
           zz1 = (1.-zgam*zal*(2. +zrgt)/zza)**.25   ! for 2 m t
           zpsih2  = psiuh(zz1, zy2)                 ! for 2 m t
           zz1 = (1.-zgam*zal*(2. +zrgq)/zza)**.25   ! for 2 m q
           zpsiq2  = psiuh(zz1, zz2)                 ! for 2 m q
           do jklev = 1, nlev_sl
             zz1 = (1.-zgam*zal*(std_lev_atm(jklev)+zrgm)/zza)**.25
             psim_sl(jklev) = psium(zz1, zx2)
             zz1 = (1.-zgam*zal*(std_lev_atm(jklev)+zrgt)/zza)**.25
             psih_sl(jklev) = psiuh(zz1, zy2)
             zz1 = (1.-zgam*zal*(std_lev_atm(jklev)+zrgq)/zza)**.25
             psiq_sl(jklev) = psiuh(zz1, zz2)
           enddo

    endif

!--------------------------------------------------------------------------------------------------
! Output variables
!--------------------------------------------------------------------------------------------------

    zlev_bottom = (phi(jlon,jlat,nlev)-phig(jlon,jlat))/g

    kturb_surf_m(jlon,jlat) = abs( -zcdm*zlev_bottom )
    kturb_surf_h(jlon,jlat) = abs( (-zcdt/zconvg)*zlev_bottom )
    kturb_surf_q(jlon,jlat) = abs( -zcdq*zlev_bottom)

!--------------------------------------------------------------------------------------------------

    zuv10             = ustar(jlon,jlat)/zak*(log((10.+zrgm)/zrgm)-zpsim10)
    u10   (jlon,jlat) = zuv10*zua/zmod
    v10   (jlon,jlat) = zuv10*zva/zmod
    q2    (jlon,jlat) = qskin(jlon,jlat)+zqstar/zak*(log((2.+zrgq)/zrgq)-zpsiq2)
    zconv2            = zconvg/(1.+ep*qskin(jlon,jlat))*(1.+ep*q2(jlon,jlat))
    t2    (jlon,jlat) = (ztevg+tstar(jlon,jlat)/zak*(log((2.+zrgt)/zrgt)-zpsih2))/zconv2
    rought(jlon,jlat) = zrgt   ! roughness for temperature diffusion
    psih  (jlon,jlat) = zpsih

! "interpolation" at standard levels of geometric height (m) above the surface, inside the surface layer 

    do jklev = 1, nlev_sl
      zuv = ustar(jlon,jlat)/zak*(log((std_lev_atm(jklev)+zrgm)/zrgm)-psim_sl(jklev))
      u_std_lev(jlon,jlat,jklev)  = zuv*zua/zmod
      v_std_lev(jlon,jlat,jklev)  = zuv*zva/zmod
      q_std_lev(jlon,jlat,jklev)  = qskin(jlon,jlat)+zqstar/zak*(log((std_lev_atm(jklev)+zrgq)/zrgq)-psiq_sl(jklev))
      zconv2 = zconvg/(1.+ep*qskin(jlon,jlat))*(1.+ep*q_std_lev(jlon,jlat,jklev))
      t_std_lev(jlon,jlat,jklev)  = &
 (ztevg+tstar(jlon,jlat)/zak*(log((std_lev_atm(jklev)+zrgt)/zrgt)-psih_sl(jklev)))/zconv2

!  Fix over high mountains...

     zglin = max (0., min((phig(jlon,jlat)-10000.)/15000., 1.))
     u_std_lev(jlon,jlat,jklev) = (1.-zglin)*u_std_lev(jlon,jlat,jklev) + zglin*u(jlon,jlat,nlev)
     v_std_lev(jlon,jlat,jklev) = (1.-zglin)*v_std_lev(jlon,jlat,jklev) + zglin*v(jlon,jlat,nlev)
     if (t_std_lev(jlon,jlat,jklev).lt.t(jlon,jlat,nlev)) then
     zglin = max (0., min((phig(jlon,jlat)-8000.)/20000., 1.))
     t_std_lev(jlon,jlat,jklev) = (1.-zglin)*t_std_lev(jlon,jlat,jklev) + zglin*t(jlon,jlat,nlev)
     q_std_lev(jlon,jlat,jklev) = (1.-zglin)*q_std_lev(jlon,jlat,jklev) + zglin*q(jlon,jlat,nlev)
     endif

    enddo

!  coefficients for comput. of surface fluxes of heat and q

    zros = ps(jlon,jlat)/(rd*tskin(jlon,jlat)*(1.+ep*qskin(jlon,jlat)))
    roscd(jlon,jlat)=  .25*zros*zcdq + .75*roscd(jlon,jlat)
    rich(jlon,jlat,nlev+1) =  min(zri, 500.)

!--------------------------------------------------------------------------------------------------

#ifdef chem
  invra(jlon,jlat) = zcdt
  rhos(jlon,jlat) = zros
#endif

!--------------------------------------------------------------------------------------------------

 enddo ! end of loop in longitude
 enddo ! end of loop in latitude

!--------------------------------------------------------------------------------------------------

! Time filter for turbulent coefficient (heat and water vapour) in the surface layer

 if (nstep_sl_filter > 1) then

   if (kstep > nstep_sl_filter) then
     do n = 1,nstep_sl_filter-1
       kturb_surf_h_mem(:,:,n) = kturb_surf_h_mem(:,:,n+1) 
       kturb_surf_q_mem(:,:,n) = kturb_surf_q_mem(:,:,n+1) 
     enddo
     kturb_surf_h_mem(:,:,nstep_sl_filter) = kturb_surf_h(:,:)
     kturb_surf_q_mem(:,:,nstep_sl_filter) = kturb_surf_q(:,:)
     kturb_surf_h_filter(:,:) = 0.
     kturb_surf_q_filter(:,:) = 0.
     do n = 1,nstep_sl_filter
       kturb_surf_h_filter(:,:) = kturb_surf_h_filter(:,:) + kturb_surf_h_mem(:,:,n)
       kturb_surf_q_filter(:,:) = kturb_surf_q_filter(:,:) + kturb_surf_q_mem(:,:,n)
     enddo
     kturb_surf_h_filter(:,:) = kturb_surf_h_filter(:,:) / float(nstep_sl_filter)
     kturb_surf_q_filter(:,:) = kturb_surf_q_filter(:,:) / float(nstep_sl_filter)
     kturb_surf_h(:,:) = kturb_surf_h_filter(:,:)
     kturb_surf_q(:,:) = kturb_surf_q_filter(:,:)
  !   if (kstep == nstep_sl_filter+1.or.mod( (kstep*int(dtstep)), 3600) == 0) then
  !     do jlat = 2, nlat-1
  !     do jlon = 2, nlon-1
  !       zmin = minval( kturb_surf_h_mem(jlon,jlat,1:nstep_sl_filter) )
  !       zmax = maxval( kturb_surf_h_mem(jlon,jlat,1:nstep_sl_filter) )
  !       if (zmax > 1..and.zmin < 1.e-2) then
  !         kturb_surf_h(jlon,jlat) = kturb_surf_h_filter(jlon,jlat)
  !         kturb_surf_q(jlon,jlat) = kturb_surf_q_filter(jlon,jlat)
  !       endif
  !     enddo 
  !     enddo 
  !   endif
   else
     kturb_surf_h_mem(:,:,kstep) = kturb_surf_h(:,:)
     kturb_surf_q_mem(:,:,kstep) = kturb_surf_q(:,:)
   endif

 endif

!  fluxes at new time step (qflux is redefined in subr. soil but only over land)
!  (note that to effectively suppress h and/or q fluxes to the atmosphere, it is necessary
!  to set zct(:) and/or zcq(:) to zero at the beginning of vdiff - but hflux and qflux
!  still affect heat and moisture fluxes into soil and sea)

 do jlat = 2, nlatm1
 do jlon = 2, nlonm1
   zzpp = pzer*sigint(nlev) - (pzer-ps(jlon,jlat))*sigalf(nlev)
   zconv = exp(rdrcp*(alp0-log(zzpp)))
   ztetav_bottom = tvirt(jlon,jlat,nlev)*zconv
   zconvg = exp(rdrcp*(alp0-log(ps(jlon,jlat))))*(1.+ep*qskin(jlon,jlat))
   ztevg = tskin(jlon,jlat)*zconvg
   zlev_bottom = (phi(jlon,jlat,nlev)-phig(jlon,jlat))/g
   zros = ps(jlon,jlat)/(rd*tskin(jlon,jlat)*(1.+ep*qskin(jlon,jlat)))
   hflux(jlon,jlat) = kturb_surf_h(jlon,jlat)*(zros*cpd)/zlev_bottom*(ztevg-ztetav_bottom)
   qflux(jlon,jlat) = kturb_surf_q(jlon,jlat)*zros/zlev_bottom*(qskin(jlon,jlat)-q(jlon,jlat,nlev))
 enddo
 enddo

    return
    end subroutine surface_layer
!##################################################################################################################
    subroutine tofd

! Turbulent orographic form drag - Beljaars et al, QJRMS, 2004

    use mod_model, only: nlon, nlat, ntop, dtstep, phi, phig, g, orogvar, u, v

    real, parameter :: ctofd = 5.e-9   ! turbulent orographic form drag coefficient

    do jklev =1,ntop
      do jlat = 1, nlat
        do jlon = 1, nlon

          zilev = (phi(jlon,jlat,jklev)-phig(jlon,jlat))/g +.01

          if (zilev < 3000.) then  ! Turbulent orographic form drag - Beljaars et al, QJRMS, 2004
            zvar = orogvar(jlon,jlat)**2
            ztofd = ctofd*sqrt(u(jlon,jlat,jklev)**2+v(jlon,jlat,jklev)**2)*zvar/zilev**1.2
           if (zilev > 500.) ztofd = ztofd*exp(-(zilev/1500.)**1.5)
           ztofd = 1./(1.+dtstep*ztofd)
           u(jlon,jlat,jklev) = u(jlon,jlat,jklev)*ztofd
           v(jlon,jlat,jklev) = v(jlon,jlat,jklev)*ztofd
         endif

        enddo
      enddo
    enddo

    return
    end subroutine tofd
!##################################################################################################################
    subroutine vdiff (pstep)

! Computes vertical diffusion of u, v, q, thetav and tke using TKE-l closure

    use mod_model, only: nlon, nlat, nlev, ntop, nlonm1, nlatm1, nlevp1, ps, u, v, t, q, tke, lml, tvirt,        &
                         rd, cpd, g, ep, eps, rdrcp, phig, tskin, qskin, sigint, dsig, pzer, tzer, ezer, alp0,   &
                         dx, ccw1, ccw2, phi, phih, tkemin, sigalf, dsigalf, ylwv, ustar, tstar,                &
                         kturb_surf_m, kturb_surf_h, kturb_surf_q, trd_a, trd_c, qcw, qci, rich, bvf, myid

#ifdef chem
  use mod_bol2chem
#endif

    real(4), dimension(nlon,nlev)   :: za, zam, zah, zc, zcm, zch, ze, zconv, ztetav, zqs, arraytmp
    real(4), dimension(nlon,nlevp1) :: zhlev, dlogthe, prandtl
    real(4), dimension(nlon,nlat)   :: ztevg
    real(4), dimension(nlon,nlat)   :: qcwskin
    real(4), dimension(nlon)        :: zct, zcq, psisurf
    real(4), parameter              :: zak=.4
    real(4) :: zlev_bottom

    zce = .17       ! ustar**2 = zce*tke
    zmlmax = 100.   ! stable mixing length maximum value
    zcbml = .37     ! adjustable parameter in unstable mixing length

#ifdef globo
    zcbml = 1.     ! adjustable parameter in unstable mixing length
#endif

    zmlcut = dx     ! mixing length maximum value
    zsqce=sqrt(zce)*pstep

!--------------------------------------------
!  loop in latitude to the end of the routine
!--------------------------------------------

 qcwskin = 0.
 bvf = 0.

 jstart = 2
 jend = nlatm1

#ifdef globo
 if (ip_s.eq.ip_null) jstart = 1
 if (ip_n.eq.ip_null) jend = nlat
#endif

 do jlat = jstart, jend
   jlatp1 = min (jlat+1, nlat) ! necessaria per globo

!  conversion factor from t to theta
!  definition of virtual theta
!  specific humidity at saturation
!  half-level height (zhlev(nlevp1)=0)

    do jklev = 1, nlev   ! loop on half levels
    do jlon = 2, nlonm1
     zzpp = pzer*sigint(jklev) - (pzer-ps(jlon,jlat))*sigalf(jklev)
     zconv(jlon,jklev) = exp(rdrcp*(alp0-log(zzpp)))
     zt0t = tzer/t(jlon,jlat,jklev)
     zesk = ezer*exp(-ccw1*log(zt0t)+ccw2*(1.-zt0t))     ! partial pressure over water
     zqs(jlon,jklev)     = zesk*eps/(zzpp+zesk*(eps-1.))
     ztetav(jlon,jklev)  = tvirt(jlon,jlat,jklev)*zconv(jlon,jklev)
     zhlev(jlon,jklev+1) = (phih(jlon,jlat,jklev+1)-phig(jlon,jlat))/g
    enddo
    enddo

!  computation of dry and moist (saturated) static stability as Durran & Klemp (1982)
!  (specific humidity is used in place of mixing ratio)

    do jklev = 2, nlev   ! loop on half levels
    do jlon = 2, nlonm1
     dthd = 2.*(ztetav(jlon,jklev-1)-ztetav(jlon,jklev))/(ztetav(jlon,jklev-1)+ztetav(jlon,jklev) ) ! dry
     r_up = zqs(jlon,jklev-1)
     r_do = zqs(jlon,jklev  )
     r_av = 0.5*(r_up + r_do)
     t_av = 0.5*(t(jlon,jlat,jklev) + t(jlon,jlat,jklev-1))
     theta_up = t(jlon,jlat,jklev-1)*zconv(jlon,jklev-1)
     theta_do = t(jlon,jlat,jklev  )*zconv(jlon,jklev  )
     zaa = (1. + ylwv*r_av/(rd*t_av))/(1. + eps*ylwv**2*r_av/(cpd*rd*t_av**2))
     dthm = zaa*((theta_up-theta_do)*2./(theta_up + theta_do) + ylwv/(cpd*t_av)*(r_up - r_do)) &
            - q(jlon,jlat,jklev-1) - qcw(jlon,jlat,jklev-1) - qci(jlon,jlat,jklev-1)           &
            + q(jlon,jlat,jklev  ) + qcw(jlon,jlat,jklev  ) + qci(jlon,jlat,jklev  )

! average relative humidity computed giving some more weight to the layer below

     zrhm = 0.55*q(jlon,jlat,jklev)/zqs(jlon,jklev) + 0.45*q(jlon,jlat,jklev-1)/zqs(jlon,jklev-1)
     zcof = max(-24.5 + 25.*zrhm, 0.)    ! zcof=0. for rh<=0.98, zcof=0.5 for rh=1
     zcof = min(zcof, .8)
     dlogthe(jlon,jklev) = zcof*dthm + (1.-zcof)*dthd  ! effective stability near saturation
    enddo
    enddo

  do jlon = 2, nlonm1

!  virtual potential temperature computed with skin temperature

    zconvg=exp(rdrcp*(alp0-log(ps(jlon,jlat))))*(1.+ep*qskin(jlon,jlat))
    ztevg(jlon,jlat)=tskin(jlon,jlat)*zconvg
    ztebar =.5*(ztevg(jlon,jlat)+ztetav(jlon,nlev))

!--------------------------------------------------------------------------------------------------

!  coefficients c(nlev) for momentum, heat, and moisture computed on t points

    zlev_bottom = (phi(jlon,jlat,nlev)-phig(jlon,jlat))/g
    zcdm = kturb_surf_m(jlon,jlat)/zlev_bottom
    zcdt = kturb_surf_h(jlon,jlat)*zconvg/zlev_bottom
    zcdq = kturb_surf_q(jlon,jlat)/zlev_bottom

    zros = ps(jlon,jlat)/(rd*tskin(jlon,jlat)*(1.+ep*qskin(jlon,jlat)))
    zdpn = pzer*dsig(nlev) - (pzer-ps(jlon,jlat))*dsigalf(nlev)
    zzc  = -pstep*g*zros/zdpn
    zcm(jlon,nlev) = zzc*zcdm
    zct(jlon) = zzc*zcdt
    zcq(jlon) = zzc*zcdq

    ztflux          = -ustar(jlon,jlat)*tstar(jlon,jlat)

!--------------------------------------------------------
!  mixing length:   stable   -> Blackadar modified
!                   unstable -> Bougeault 1986 (modified)
!--------------------------------------------------------

    jkbot = nlev      ! half-integer level

    do 50 jklev = nlev, ntop, -1   ! loop on half-integer levels

    if (dlogthe(jlon,jklev).gt.0..or.jklev.eq.ntop) then

!  mixing length (stable case)

    zblack = zak*zhlev(jlon,jklev)*zmlmax/(zak*zhlev(jlon,jklev)+zmlmax)
    lml(jlon,jlat,jklev) = zblack

!  mixing length (unstable layer below current stable level, if present)

    do jk = jkbot, jklev+1, -1
    zlbot = zhlev(jlon,jk)-zhlev(jlon,jkbot+1)
    zltop = zhlev(jlon,jklev)-zhlev(jlon,jk)
    lml(jlon,jlat,jk) = min (zcbml*sqrt(zlbot*sqrt(zlbot*zltop)), zmlcut)
    lml(jlon,jlat,jk) = max (lml(jlon,jlat,jk), zblack)
    enddo
    jkbot = jklev-1
    endif

 50 continue

!--------------------------------------------------------
!  tke lower boundary condition
!--------------------------------------------------------

    tke(jlon,jlat,nlevp1) = ustar(jlon,jlat)**2/zce+tkemin
    zwstar = 0.
    jkubl  = nlevp1  ! top of unstable boundary layer

    if (ztflux.gt.0.) then  ! in case of unstable surface layer

    do jklev = nlev, ntop, -1  ! loop on half-integer levels
      if (dlogthe(jlon,jklev).gt.0..or.jklev.eq.ntop) then
      jkubl = jklev
      exit
      endif
    enddo

    zwstar = (g*zhlev(jlon,jkubl)*ztflux/ztebar)**(1./3.)
    tke(jlon,jlat,nlevp1) = tke(jlon,jlat,nlevp1) + .4*zwstar**2

    endif

#ifdef chem
  w_star(jlon,jlat) = zwstar
#endif

   enddo ! end of loop in longitude

!--------------------------------------------------------
!  tke equation
!--------------------------------------------------------

    do jklev = nlev, 2, -1
    do jlon = 2, nlonm1
    zrdz  = g/(phi(jlon,jlat,jklev-1)-phi(jlon,jlat,jklev))
    zdudz = .5*(u(jlon,jlat,jklev-1)+u(jlon-1,jlat,jklev-1)-u(jlon,jlat,jklev)-u(jlon-1,jlat,jklev))*zrdz
    zdvdz = .5*(v(jlon,jlat,jklev-1)+v(jlon,jlatp1,jklev-1)-v(jlon,jlat,jklev)-v(jlon,jlatp1,jklev))*zrdz
    zshear= zdudz**2 + zdvdz**2
    zbuoy = g*dlogthe(jlon,jklev)*zrdz
    zrich = min (zbuoy/(zshear + 1.e-6), 500.)
    rich(jlon,jlat,jklev) = zrich

! Brunt-Vaisala frequency squared, averaged near surface with different weights,
! to be used in orogdrag and in blockdrag

    if(jklev==nlev)   bvf(jlon,jlat) = 0.2*zbuoy
    if(jklev==nlev-1) bvf(jlon,jlat) = bvf(jlon,jlat) + 0.5*zbuoy
    if(jklev==nlev-2) bvf(jlon,jlat) = bvf(jlon,jlat) + 0.45*zbuoy

      if (zrich > 0.) then
      prandtl(jlon,jklev) = .75 + 1.5*zrich + .5*zrich*sqrt(zrich) ! (after Anderson, Bound. L. Met., 2009)
!      zfm = 1./(1. + 12.*zrich)
!      lml(jlon,jlat,jklev) = lml(jlon,jlat,jklev)*zfm  ! Blackadar modified
      lml(jlon,jlat,jklev) = min(0.52*sqrt(zce*tke(jlon,jlat,jklev)/zbuoy), lml(jlon,jlat,jklev)) ! Deardorff (1980)
      else
      prandtl(jlon,jklev) = .75
      endif
    zp = zsqce*lml(jlon,jlat,jklev)*(zshear-zbuoy/prandtl(jlon,jklev))
    zd = zsqce*zce/(lml(jlon,jlat,jklev)+.01)
    tke(jlon,jlat,jklev) = tke(jlon,jlat,jklev) + sqrt(tke(jlon,jlat,jklev))*(-zd*tke(jlon,jlat,jklev)+zp)
    if (tke(jlon,jlat,jklev).lt.tkemin) tke(jlon,jlat,jklev) = tkemin
    if (tke(jlon,jlat,jklev).gt.20.   ) tke(jlon,jlat,jklev) = 20.
    enddo
    enddo
    rich(:,jlat,1) = rich(:,jlat,2)

!  frictional heating

    do jklev = ntop+1, nlev
    do jlon = 2, nlonm1
    zlmlm = .5*(lml(jlon,jlat,jklev)+lml(jlon,jlat,jklev+1))
    zlmlm = max(zlmlm, 1.5)
    ztkem = .5*zce*(tke(jlon,jlat,jklev)+tke(jlon,jlat,jklev+1))
    zdtfric = ztkem*sqrt(ztkem)/zlmlm
    ztetav(jlon,jklev) = ztetav(jlon,jklev) + zdtfric*zconv(jlon,jklev)*pstep/cpd
    enddo
    enddo

!  vertical diffusion of tke on half-levels k, k=ntop+1, nlev
!  coefficients a(k), c(k) [a(ntop+1)=0]

    za(:,ntop+1) = 0.
    do jklev = ntop+1, nlev  ! loop on levels
    do jlon = 2, nlonm1
    zkk   = .5*(lml(jlon,jlat,jklev+1)+lml(jlon,jlat,jklev))*          &
            sqrt(zce*.5*(tke(jlon,jlat,jklev+1)+tke(jlon,jlat,jklev)))
    zp    = pzer*sigint(jklev-1)-(pzer-ps(jlon,jlat))*sigalf(jklev-1)
    zpp   = pzer*sigint(jklev  )-(pzer-ps(jlon,jlat))*sigalf(jklev  )
    zdpk  = zpp-zp
    zdpkp = pzer*dsig(jklev) - (pzer-ps(jlon,jlat))*dsigalf(jklev)
    zdzp  = zhlev(jlon,jklev)-zhlev(jlon,jklev+1)
    zc(jlon,jklev) = -pstep*zdpkp/(zdpk*zdzp**2)*zkk
      if(jklev.lt.nlev) then
      zppp  = pzer*sigint(jklev+1)-(pzer-ps(jlon,jlat))*sigalf(jklev+1)
      za(jlon,jklev+1) = zc(jlon,jklev)*zdpk/(zppp-zpp)
      endif
    enddo
    enddo

!  tridiagonal matrix inversion for tke at new time level

    do jlon = 2, nlonm1
    zrb =1./(1.-zc(jlon,ntop+1))
    ze(jlon,ntop+1)=-zc(jlon,ntop+1)*zrb
    tke(jlon,jlat,ntop+1)= tke(jlon,jlat,ntop+1)*zrb
    enddo

    do jklev = ntop+2, nlev
    do jlon = 2, nlonm1
    zb = 1.-za(jlon,jklev)-zc(jlon,jklev)
    zrden = 1./(zb+za(jlon,jklev)*ze(jlon,jklev-1))
    ze(jlon,jklev)=-zc(jlon,jklev)*zrden
    tke(jlon,jlat,jklev)=(tke(jlon,jlat,jklev)-za(jlon,jklev)*tke(jlon,jlat,jklev-1))*zrden
    enddo
    enddo

    do jklev = nlev, ntop+1, -1
    do jlon = 2, nlonm1
    tke(jlon,jlat,jklev) = ze(jlon,jklev)*tke(jlon,jlat,jklev+1)+tke(jlon,jlat,jklev)
    if (tke(jlon,jlat,jklev).lt.tkemin) tke(jlon,jlat,jklev) = tkemin
    enddo
    enddo

!-----------------------------------------------
!  vertical diffusion of u, v, thetav, q, qcw, qci
!  coefficients c(k), k=ntop, nlev-1
!  coefficients a(k), k=ntop, nlev   [a(ntop)=0]
!-----------------------------------------------

    zam(:,ntop) = 0.
    zah(:,ntop) = 0.

    do jklev = ntop, nlev-1  ! loop on levels
    do jlon = 2, nlonm1

!  flux/gradient coefficients on half-level jklev+1/2 and t points

    zkkm = lml(jlon,jlat,jklev+1)*sqrt(zce*tke(jlon,jlat,jklev+1))
    zkkh = zkkm/prandtl(jlon,jklev+1)

    zp     = pzer*sigint(jklev  ) - (pzer-ps(jlon,jlat))*sigalf(jklev  )
    zpp    = pzer*sigint(jklev+1) - (pzer-ps(jlon,jlat))*sigalf(jklev+1)
    zdpkp  = zpp-zp
    zdpk   = pzer*dsig(jklev  ) - (pzer-ps(jlon,jlat))*dsigalf(jklev  )
    zdpkp1 = pzer*dsig(jklev+1) - (pzer-ps(jlon,jlat))*dsigalf(jklev+1)
    zdz    = (phi(jlon,jlat,jklev)-phi(jlon,jlat,jklev+1))/g
    zcm(jlon,jklev  ) = -pstep*zdpkp/(zdpk*zdz**2)*zkkm
    zam(jlon,jklev+1) =  zcm(jlon,jklev)*zdpk/zdpkp1
    zch(jlon,jklev  ) = -pstep*zdpkp/(zdpk*zdz**2)*zkkh
    zah(jlon,jklev+1) =  zch(jlon,jklev)*zdpk/zdpkp1

#ifdef chem
    kh( jlon,jlat,jklev+1 ) = zkkh
#endif

    enddo
    enddo

!  save coefficients ah and ch for vertical diffusion of chemical parameters
!  ch(nlev) is undefined (zero for zero flux)

#ifdef chem
    trd_a( :,jlat,ntop:nlev ) = zah( :,ntop:nlev )
    trd_c( :,jlat,ntop:nlev-1 ) = zch( :,ntop:nlev-1 )
    trd_c( :,jlat,nlev ) = zct( : )
#endif

!  tridiagonal matrix inversion

    psisurf = 0.
    arraytmp(1:nlon,1:nlev) = u(1:nlon,jlat,1:nlev)
    call tridiag (zam, zcm, psisurf, arraytmp, nlon, nlev, ntop)           ! u
    u(1:nlon,jlat,1:nlev) = arraytmp(1:nlon,1:nlev)

    arraytmp(1:nlon,1:nlev) = v(1:nlon,jlat,1:nlev)
    call tridiag (zam, zcm, psisurf, arraytmp, nlon, nlev, ntop)           ! v
    v(1:nlon,jlat,1:nlev) = arraytmp(1:nlon,1:nlev)

    zch(:,nlev) = zct(:)
    call tridiag (zah, zch, ztevg(1,jlat), ztetav, nlon, nlev, ntop)       ! virtual theta

    zch(:,nlev) = zcq(:)
    arraytmp(1:nlon,1:nlev) = q(1:nlon,jlat,1:nlev)
    call tridiag (zah, zch, qskin(1,jlat), arraytmp, nlon, nlev, ntop)     ! q
    q(1:nlon,jlat,1:nlev) = arraytmp(1:nlon,1:nlev)

!    qcwskin(1:nlon, jlat) = qcw(1:nlon,jlat,nlev)      ! lower b.c. of no flux
    arraytmp(1:nlon,1:nlev) = qcw(1:nlon,jlat,1:nlev)
    call tridiag (zah, zch, qcwskin(1,jlat), arraytmp, nlon, nlev, ntop)   ! qcw (lower b.c. 0.)
    qcw(1:nlon,jlat,1:nlev) = arraytmp(1:nlon,1:nlev)

!    qcwskin(1:nlon, jlat) = qci(1:nlon,jlat,nlev)      ! lower b.c. of no flux
    arraytmp(1:nlon,1:nlev) = qci(1:nlon,jlat,1:nlev)
    call tridiag (zah, zch, qcwskin(1,jlat), arraytmp, nlon, nlev, ntop)   ! qci (lower b.c. 0.)
    qci(1:nlon,jlat,1:nlev) = arraytmp(1:nlon,1:nlev)

!  conversion from new virtual theta to t

    do jklev = ntop, nlev
    do jlon = 2, nlonm1
    t(jlon,jlat,jklev) = ztetav(jlon,jklev)*t(jlon,jlat,jklev)/(tvirt(jlon,jlat,jklev)*zconv(jlon,jklev))
    enddo
    enddo

 enddo !  end of loop in latitude

    return
    end subroutine vdiff
!##################################################################################################################
      subroutine runout_g (jstep, nyrc, nmonc, ndayc, nhouc, nminc)

! For globo, computes run-time diagnostics

      use mod_model
      implicit none
      real(4) zpmin, zpmax, zumax, zzm, zdummy, zday, zmass
      integer jstep, jlon, jlat, jklev, j1, j2, ierr, nyrc, nmonc, ndayc, nhouc, nminc
#ifdef mpi
        include 'mpif.h'
#endif

      zpmin =  1.e8
      zpmax = -1.e8
      zumax = -1.e8
      j1=2
      j2=nlatm1
      if (ip_s.eq.ip_null) j1=1
      if (ip_n.eq.ip_null) j2=nlat

      do jlat = j1, j2
      do jlon = 2, nlonm1
      if (ps(jlon,jlat).lt.zpmin) zpmin=ps(jlon,jlat)
      if (ps(jlon,jlat).gt.zpmax) zpmax=ps(jlon,jlat)
      enddo
      enddo

      do jklev = 1, nlev
      do jlat = j1, j2
      do jlon = 2, nlonm1
      zzm = sqrt(u(jlon,jlat,jklev)**2+v(jlon,jlat,jklev)**2)
      if (zzm.gt.zumax) zumax=zzm
      enddo
      enddo
      enddo

!  calculate global minimum and maximum

#ifdef mpi
        call mpi_reduce(zpmin, zdummy, 1, mpi_real, mpi_min, 0, mpi_comm_world, ierr)
        zpmin = zdummy
        call mpi_reduce(zpmax, zdummy, 1, mpi_real, mpi_max, 0, mpi_comm_world, ierr)
        zpmax = zdummy
        call mpi_reduce(zumax, zdummy, 1, mpi_real, mpi_max, 0, mpi_comm_world, ierr)
        zumax = zdummy
#endif

!      call tmass (ps, nlon, nlat, hxt, dlon, dlat, myid, nprocsx, nprocsy, zmass)
      call tqmass (ps, q, nlon, nlat, nlev, hxt, dlon, dlat, myid, nprocsx, nprocsy, &
                   dsig, dsigalf, pzer, g, zmass)

#ifdef mpi
        call mpi_reduce(zmass, zdummy, 1, mpi_real, mpi_sum, 0, mpi_comm_world, ierr)
        zmass = zdummy
#endif

      if (myid.eq.0) then
        zday=jstep*dtstep/86400.
        print 1010, jstep, zday, zpmin/100., zpmax/100., zumax, zmass, nyrc, nmonc, ndayc, nhouc, nminc
      endif

! 1010 format(' jstep',i8,'    day', f8.2,'    min/max Ps=',2f8.2,'    Vmax=',f6.1,'    total mass=',f10.1,/, &
 1010 format(' jstep',i8,'    day', f8.2,'    min/max Ps=',2f8.2,'    Vmax=',f6.1,'    tot. q mass=',f10.2,/, &
             ' Current Date (Year/Month/Day/Hour/min) ', i5, 4i3)
      return
      end subroutine runout_g
!##################################################################################################################
      subroutine runout_b (jstep)

! For bolam, computes run-time diagnostics

      use mod_model
      implicit none
      real(4) zpmin, zpmax, zumax, zvmax, ztsmin, ztsmax, zdummy, zday, zhour
      integer jstep, jlon, jlat, jklev, j1, j2, ierr
#ifdef mpi
        include 'mpif.h'
#endif

      zpmin  =  1.e8
      zpmax  = -1.e8
      zumax  = -1.e8
      zvmax  = -1.e8
      ztsmin =  1.e8
      ztsmax = -1.e8
      j1=2
      j2=nlat-1
      if (mod(myid,nprocsy).eq.0        ) j1=1
      if (mod(myid,nprocsy).eq.nprocsy-1) j2=nlat

      do jlat = j1, j2
      do jlon = 2, nlonm1
      if (ps(jlon,jlat).lt.zpmin) zpmin=ps(jlon,jlat)
      if (ps(jlon,jlat).gt.zpmax) zpmax=ps(jlon,jlat)
      if (tskin(jlon,jlat).lt.ztsmin) ztsmin=tskin(jlon,jlat)
      if (tskin(jlon,jlat).gt.ztsmax) ztsmax=tskin(jlon,jlat)
      do jklev = 1, nlev
        if (abs(u(jlon,jlat,jklev)).gt.zumax) zumax=abs(u(jlon,jlat,jklev))
        if (abs(v(jlon,jlat,jklev)).gt.zvmax) zvmax=abs(v(jlon,jlat,jklev))
      enddo
      enddo
      enddo

!  calculate global minimum and maximum

#ifdef mpi
        call mpi_reduce(zpmin, zdummy, 1, mpi_real, mpi_min, 0, mpi_comm_world, ierr)
        zpmin = zdummy
        call mpi_reduce(zpmax, zdummy, 1, mpi_real, mpi_max, 0, mpi_comm_world, ierr)
        zpmax = zdummy
        call mpi_reduce(ztsmin, zdummy, 1, mpi_real, mpi_min, 0, mpi_comm_world, ierr)
        ztsmin = zdummy
        call mpi_reduce(ztsmax, zdummy, 1, mpi_real, mpi_max, 0, mpi_comm_world, ierr)
        ztsmax = zdummy
        call mpi_reduce(zumax, zdummy, 1, mpi_real, mpi_max, 0, mpi_comm_world, ierr)
        zumax = zdummy
        call mpi_reduce(zvmax, zdummy, 1, mpi_real, mpi_max, 0, mpi_comm_world, ierr)
        zvmax = zdummy
#endif

      if (myid.eq.0) then
        zday  = jstep*dtstep/86400.
        zhour = jstep*dtstep/3600.
!!!        print 1010, jstep, zhour, zpmin/100., zpmax/100.
        write (*,'(a,i8,a,f8.2,3(a,2f8.2))') ' jstep =', jstep,',   hour =', zhour, &
 ',   min/max ps =', zpmin/100., zpmax/100.,',   min/max tsurf=', ztsmin-273.15, ztsmax-273.15, &
 ',   max u,v=', zumax, zvmax
        print*
      endif

 1010 format(' jstep =', i8, ',   hour =', f8.2, ',   min/max ps =', 2f8.2)
      return
      end subroutine runout_b
!##################################################################################################################
    subroutine micro

! Microphysics: explicit cloud and precipitation processes

    use mod_model, only: nlon, nlat, nlev, ntop, q, qcw, qci, qrain, qsnow, qhail, t, ps, tke,            &
                         ccw1, ccw2, cci1, cci2, cpd, cpv, rd, rv, eps, cw, ci, yliv, yliw, ylwv, tzer,   &
                         pzer, ezer, pi, dtstep, fmask, rainls, snowls, sigint, sigalf, dsig, dsigalf, g, &
                         dqcwdt, dqcidt, dqrndt, dqsndt, &
                         nlatm1, ip_n, ip_s, ip_null, myid

#ifdef chem
      use mod_bol2chem, only: precipitation_below, precipitation_in, precip_intercepted_volume
#endif

    implicit none

! microphysical parameters

    real, parameter :: alfacw=6., alfaci=3.                    ! distribution parameters of clouds
    real, parameter :: yncw0=3.2e7, ynci0=5.7e7                ! numbers of cloud particles (continental)
    real, parameter :: yn0r=5.e6, yn0s=1.0e7, yn0h=6.2e6       ! distribution parameters of hydrometeors
    real, parameter :: acw=pi/6.*1000., bcw=3.                 ! cloud water mass parameters
    real, parameter :: aci=100., bci=2.55                      ! cloud ice mass parameters
    real, parameter :: ar =pi/6.*1000., br =3.                 ! rain mass parameters
    real, parameter :: as =30., bs =2.7                        ! snow mass parameters
    real, parameter :: ah =ar*.75, bh =3.                      ! graupel/hail mass parameters
    real, parameter :: ykc=.3e8,  ync= 2.                      ! fall speed parameters (cloud)
    real, parameter :: ykr=900.,  ynr=.84                      ! fall speed parameters (rain)
    real, parameter :: yks=132.,  yns=.74                      ! fall speed parameters (snow)
    real, parameter :: ykh=1.7e5, ynh=1.8                      ! fall speed parameters (graupel/hail)

    real, parameter :: zqcwth=5.e-6, zqcith=5.e-6              ! autoconversion thresholds
    real, parameter :: d0w=0.58e-4, d0i=0.125e-4               ! autoconversion diameters
    real, parameter :: zqpwth=1.e-4, ztcrit=269.               ! graupel/hail formation
    real, parameter :: ee_rcw=.7, ee_scw=.8, ee_sci=.5         ! collection efficiencies
    real, parameter :: ee_hcw=.7, ee_hci=.7, ee_rs =.8         ! collection efficiencies
    real, parameter :: fvent=.8, chi=2.26e-5, yka=2.43e-2
    real, parameter :: ymu=1.718e-5, schmidt=.6
    real, parameter :: zqmin=1.e-8                             ! minimum value of specific concentration

    real, dimension(nlon,nlev) :: zfsrain, zfssnow, zfshail, zflurn, zflusn, zfluha, zrodz, &
                                  zfsqcw, zfsqci, zfsncw, zfsnci, zflucw, zfluci

    integer jlon, jlat, jklev, jstart, jend
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
    real zmelmax, zmeltsn, zcapt, zboul, zboui, zwes
    real zrh, ztw, zt0tw, zrip, zzqsw, ztwfg

#ifdef chem
  ! volume intercepted by rain per unit volume [precip_intercepted_volume] = [1]
  precip_intercepted_volume = 0.0
#endif

! gamma functions

      zgalfw1 = eugamma (alfacw+1.        )
      zgalfwb = eugamma (alfacw+bcw+1.    )
      zgacw1  = eugamma (ync+bcw+alfacw+1.)

      zgalfi1 = eugamma (alfaci+1.        )
      zgalfib = eugamma (alfaci+bci+1.    )
      zgaci1  = eugamma (ync+bci+alfaci+1.)

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
      z31   = .31*schmidt**(1./3.)/sqrt(ymu)

! collection coefficients

      cc_rcw = .25*pi*ee_rcw*yn0r*ykr*zgar1*dtstep
      cc_rs  = .25*pi*ee_rs *yn0r*ykr*zgar1*dtstep
      cc_scw = .25*pi*ee_scw*yn0s*yks*zgas1*dtstep
      cc_sci = .25*pi*ee_sci*yn0s*yks*zgas1*dtstep
      cc_hcw = .25*pi*ee_hcw*yn0h*ykh*zgah1*dtstep
      cc_hci = .25*pi*ee_hci*yn0h*ykh*zgah1*dtstep

!-------------------------------------------------------------------------------------

      jstart = 2
      jend = nlatm1

#ifdef globo
      if (ip_s.eq.ip_null) jstart = 1
      if (ip_n.eq.ip_null) jend = nlat
#endif

      do 1000 jlat = jstart, jend

! Initial values of terminal fall speeds
      zfsqcw = -.002
      zfsqci = -.002
      zfsrain = -.1
      zfssnow = -.01
      zfshail = -.1

      do jklev = ntop, nlev
      do jlon = 2, nlon-1

      zt    = t    (jlon,jlat,jklev)
      zp    = pzer*sigint(jklev)-(pzer-ps(jlon,jlat))*sigalf(jklev)
      zqv   = q    (jlon,jlat,jklev)
      zqcw  = qcw  (jlon,jlat,jklev)
      zqci  = qci  (jlon,jlat,jklev)
      zqpw  = qrain(jlon,jlat,jklev)
      zqpi1 = qsnow(jlon,jlat,jklev)
      zqpi2 = qhail(jlon,jlat,jklev)
      yncw  = 1.05*yncw0*(1.-fmask(jlon,jlat))+.85*yncw0*fmask(jlon,jlat)
      ynci  = 1.05*ynci0*(1.-fmask(jlon,jlat))+.85*ynci0*fmask(jlon,jlat)
      zboul = 1.
      zboui = zboul * max(.01, min((zp-100.e2)/300.e2, 1.))   ! to reduce high clouds
      zqtot = zqv +zqcw +zqci +zqpw +zqpi1 +zqpi2
      zro   = zp/((rd*(1.-zqtot)+rv*zqv)*zt)        ! air density
      zcoe  = sqrt(1./zro)                          ! coefficient to increase terminal velocity with height

! calculation of enthalpy zh before microphysical processes

      zctot =  (1.-zqtot)*cpd + zqtot*ci
      zh = zctot*zt + (yliv+(cpv-ci)*(zt-tzer))*zqv + (yliw+(cw-ci)*(zt-tzer))*(zqcw+zqpw)

! saturation specific humidity with respect to water and ice

      zt0t = tzer/zt
      zesk = ezer*exp(-ccw1*log(zt0t)+ccw2*(1.-zt0t)) !  saturated pressure over water (in Pa)
      zqsw = max( zesk*eps/(zp+zesk*(eps-1.)), 1.e-9)
      if (zt.lt.tzer) then
      zesk = ezer*exp(-cci1*log(zt0t)+cci2*(1.-zt0t)) !  saturated pressure over ice (in Pa)
      zqsi = zesk*eps/(zp+zesk*(eps-1.))
      else
      zqsi = zqsw
      endif
      zqsi = max( zqsi, 1.e-9)

! 1) (no) nucleation of cloud water

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
      if(zdcon.lt.0..and.dqcwdt(jlon,jlat,jklev).gt.1.e-8) &
         zdcon = 0.02*zdcon     ! reduce evap. of conv. cloud water
      zqcw = zqcw + zdcon       ! update cloud water
      zqv  = zqv  - zdcon
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
      endif

! 4) heterogeneous nucleation of cloud ice (formation of pristine crystals)

      if (zt.lt.269.) then
      ssi = min( (zqv-zqsi)/zqsi, 3.25) ! ---> yncic=1.e18 as maximum (it was happen value 7, so yncic=1.e72 = Floating-point exception
        if (ssi.gt.0) then
        yncic = 1.e4*min (1.,.25*(269.-zt) ) * exp(-.639+12.96*ssi)
          if (1.e-12*yncic.gt.zqci) then
          zdqci = min (1.e-12*yncic-zqci, zqv)
          zqci  = zqci + zdqci
          zqv   = zqv  - zdqci
          endif
        endif
      endif

! 5) sublimation rate of cloud ice in both directions (part of Bergeron-Findeisen process)

      if (zt.lt.tzer) then
      if (zqv.gt.zqsi.or.zqci.gt.zqmin) then
      zqci = max(zqci, 0.)
      zdi= zqv-zqsi
      z1 = zqsi*yliv**2*chi*zro/(yka*rv*zt**2)
      zsubl = 2.*pi*fvent*zdi*chi/(1.+z1)

! AUGH FIXME: see similar solution above

      zbase = max(zqci, 1.e-9)*zro/ynci ! poisk
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
      if(zdsub.lt.0..and.dqcidt(jlon,jlat,jklev).gt.1.e-8) &
         zdsub = 0.02*zdsub     ! reduce subl. of conv. cloud ice
      zqci = zqci + zdsub       ! update cloud ice
      zqv  = zqv  - zdsub
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
      endif

!-------------------------------------------------------------------------------------------------------------

! 7) cloud water converted into cloud ice due to cloud ice-water interaction

     if (zt.lt.tzer-.5.and.zqcw.gt.zqmin.and.zqci.gt.zqmin) then
     zdfre = min(5.*zro**2*dtstep*zqcw*zqci, .5*zqcw)
     zqcw = zqcw - zdfre
     zqci = zqci + zdfre
     endif

! 8) autoconversion of cloud water into rain

      if (zqcw.gt.zqcwth) then
      zbetw = ((zgalfwb*acw*yncw)/(zqcw*zro*zgalfw1))**zexpw
      zarg = d0w*zbetw*zboul
      if (zarg.lt.25.) then
      zgammacw = .5*(1.-tanh((zarg-8./zarg-8.8)/4.3))  ! partial gamma(10,zarg)
      zdq  = zgammacw*zqcw
      zqcw = zqcw - zdq
      zqpw = zqpw + zdq
      endif
      endif

! 9) autoconversion of cloud ice into snow

      if (zqci.gt.zqcith.and.zt.lt.tzer) then
      zbeti = ((zgalfib*aci*ynci)/(zqci*zro*zgalfi1))**zexpi
      zarg = d0i*zbeti*zboui
      if (zarg.lt.25.) then
      zgammaci = .5*(1.-tanh((zarg-8./zarg-4.9)/3.9))
      zdq   = zgammaci*zqci
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
!        zzqsw = zesk*eps(zp-zesk)                ! with mixing ratio
        ztw = zt + yliv/cpd*(zqv-zzqsw) ! wet bulb temperature for ice
!        ztw = zt + yliv/cpd*(zqv/(1.-zqv)-zzqsw) ! with mixing ratio
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
      if(dqsndt(jlon,jlat,jklev).gt.1.e-8) zdq = 0.05*zdq  ! reduce subl. of conv. snow
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
      if(dqrndt(jlon,jlat,jklev).gt.1.e-8) zdq = 0.05*zdq  ! reduce evap. of conv. rain
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

#ifdef chem

! volume intercepted by rain per unit volume [precip_intercepted_volume] = [1]
! upper limited by half of the cell grid volume according to the removal of cloud water (see below)
! precip_intercepted_volume(jlon,jlat,jklev) = min( (cc_rcw/ee_rcw)*zcoe*zlambrm1, 0.5 )

  precip_intercepted_volume(jlon,jlat,jklev) = precip_intercepted_volume(jlon,jlat,jklev) + &
                                               (cc_rcw/ee_rcw)*zcoe*zlambrm1
#endif

      zqcw = zqcw - zdq
      zqpw = zqpw + zdq
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
      endif

! 19) snow and cloud water interaction: riming of snow below tzer-1.
!     snow melting due to enthalpy of cloud water above tzer

      if (zqpi1.gt.zqmin.and.zqcw.gt.zqmin) then
      zlambsm1 = (zro*zqpi1/(yn0s*as*zgas2))**((yns+3.)/(bs+1.))
      zdq   = cc_scw*zcoe*zqcw*zlambsm1
      zdq   = min (zdq, zqcw*.5)  ! intercepted cloud water

#ifdef chem

! volume intercepted by rain per unit volume [precip_intercepted_volume] = [1]
! upper limited by half of the cell grid volume according to the removal of cloud water (see below)
! precip_intercepted_volume(jlon,jlat,jklev) = min( (cc_scw/ee_scw)*zcoe*zlambsm1, 0.5 )

  precip_intercepted_volume(jlon,jlat,jklev) = precip_intercepted_volume(jlon,jlat,jklev) + &
                                               (cc_scw/ee_scw)*zcoe*zlambsm1
#endif

      zrip = min(zqcw*1.3e3, 1.)  ! determines redistribution of rimes cloud water into snow and graupel/hail
      zqcw  = zqcw  - zdq
        if (zt.lt.tzer-0.5) then  ! riming
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

      if (zqpi2.gt.zqmin.and.zqci.gt.zqmin.and.zt.lt.tzer.and.zt.gt.253) then ! negleg. for t<-20 C (Houze)
      zlambhm1 = (zro*zqpi2/(yn0h*ah*zgah2))**((ynh+3.)/(bh+1.))
      zdq = cc_hci*zcoe*zqci*zlambhm1 *30./(30.+(tzer-zt)**2)  ! modulation decreasing with t (Houze)
      zdq = min(zdq, zqci*.5)   ! intercepted cloud ice
      zqci  = zqci  - zdq
      zqpi2 = zqpi2 + zdq
      endif

! 21) graupel/hail from freezing cloud water
!     rain from cloud water and melting graupel/hail

      if (zqpi2.gt.zqmin.and.zqcw.gt.zqmin) then
      zlambhm1 = (zro*zqpi2/(yn0h*ah*zgah2))**((ynh+3.)/(bh+1.))
      zdq = cc_hcw*zcoe*zqcw*zlambhm1
      zdq = min(zdq, zqcw*.5)                                 ! intercepted cloud water

#ifdef chem
  ! volume intercepted by rain per unit volume [precip_intercepted_volume] = [1]
  precip_intercepted_volume(jlon,jlat,jklev) = precip_intercepted_volume(jlon,jlat,jklev) + &
                                               (cc_hcw/ee_hcw)*zcoe*zlambhm1
  ! upper limited by half of the cell grid volume according to the removal of cloud water (see below)
  !precip_intercepted_volume(jlon,jlat,jklev) = min( (cc_scw/ee_scw)*zcoe*zlambsm1, 0.5 )
#endif

      zqcw = zqcw - zdq
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
      endif
      if (zqci.gt.zqmin) then
      zbeta  = (ynci*aci*zgalfib/(zro*zqci*zgalfi1))**(zexpi*ync)
      zfsqci(jlon,jklev) = -ykc*zgaci1/(zgalfib*zbeta)
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

      q    (jlon,jlat,jklev) = zqv
      qcw  (jlon,jlat,jklev) = zqcw
      qci  (jlon,jlat,jklev) = max(zqci, 0.)
      qrain(jlon,jlat,jklev) = zqpw
      qsnow(jlon,jlat,jklev) = max(zqpi1, 0.)
      qhail(jlon,jlat,jklev) = max(zqpi2, 0.)

      enddo  ! end loop on longitude
      enddo  ! end loop on levels

!  fall of precipitation by a backward-upstream scheme
!  accumulation of precipitation at the ground

      zflucw = 0.
      zfluci = 0.
      zflurn = 0.
      zflusn = 0.
      zfluha = 0.

      do jklev = ntop-1, nlev
      do jlon = 2, nlon-1
      zp = pzer*sigint(jklev)-(pzer-ps(jlon,jlat))*sigalf(jklev)
      zro = zp/(rd*t(jlon,jlat,jklev))
      zrodz(jlon,jklev)  = (pzer*dsig(jklev)-(pzer-ps(jlon,jlat))*dsigalf(jklev))/g
      zflucw(jlon,jklev) = zfsqcw (jlon,jklev)*zro*dtstep
      zfluci(jlon,jklev) = zfsqci (jlon,jklev)*zro*dtstep
      zflurn(jlon,jklev) = zfsrain(jlon,jklev)*zro*dtstep
      zflusn(jlon,jklev) = zfssnow(jlon,jklev)*zro*dtstep
      zfluha(jlon,jklev) = zfshail(jlon,jklev)*zro*dtstep
      enddo
      enddo

      do jklev = ntop, nlev
      do jlon = 2, nlon-1
      zqcwnew  = (qcw  (jlon,jlat,jklev)*zrodz(jlon,jklev)-qcw  (jlon,jlat,jklev-1)*zflucw(jlon,jklev-1)) &
                /(zrodz(jlon,jklev)-zflucw(jlon,jklev))
      zqcinew  = (qci  (jlon,jlat,jklev)*zrodz(jlon,jklev)-qci  (jlon,jlat,jklev-1)*zfluci(jlon,jklev-1)) &
                /(zrodz(jlon,jklev)-zfluci(jlon,jklev))
      zqpwnew  = (qrain(jlon,jlat,jklev)*zrodz(jlon,jklev)-qrain(jlon,jlat,jklev-1)*zflurn(jlon,jklev-1)) &
                /(zrodz(jlon,jklev)-zflurn(jlon,jklev))
      zqpi1new = (qsnow(jlon,jlat,jklev)*zrodz(jlon,jklev)-qsnow(jlon,jlat,jklev-1)*zflusn(jlon,jklev-1)) &
                /(zrodz(jlon,jklev)-zflusn(jlon,jklev))
      zqpi2new = (qhail(jlon,jlat,jklev)*zrodz(jlon,jklev)-qhail(jlon,jlat,jklev-1)*zfluha(jlon,jklev-1)) &
                /(zrodz(jlon,jklev)-zfluha(jlon,jklev))

!  heat capacity of precipitation (cloud excluded)

      zdqwat=abs(zflurn(jlon,jklev-1)*qrain(jlon,jlat,jklev-1))
      zdqice=abs(zflusn(jlon,jklev-1)*qsnow(jlon,jlat,jklev-1))+abs(zfluha(jlon,jklev-1)*qhail(jlon,jlat,jklev-1))
      t(jlon,jlat,jklev)=( cpd*zrodz(jlon,jklev)*t(jlon,jlat,jklev)+(zdqwat*cw+zdqice*ci)*t(jlon,jlat,jklev-1)) &
                        /(cpd*zrodz(jlon,jklev)+zdqwat*cw+zdqice*ci)

      qcw  (jlon,jlat,jklev) = zqcwnew
      qci  (jlon,jlat,jklev) = zqcinew
      qrain(jlon,jlat,jklev) = zqpwnew
      qsnow(jlon,jlat,jklev) = zqpi1new
      qhail(jlon,jlat,jklev) = zqpi2new

#ifdef chem

!  flux is negative downward while precipitation must be positive

      precipitation_in(jlon,jlat,jklev)    = -1 * qcw(jlon,jlat,jklev)   * zflucw(jlon,jklev)/dtstep
      precipitation_below(jlon,jlat,jklev) = -1 * qrain(jlon,jlat,jklev) * zflurn(jlon,jklev)/dtstep
#endif

      enddo
      enddo

!  instantaneous rain, snow and graupel/hail (mm) reaching the ground

      rainls(:,jlat) = -zflurn(:,nlev)*qrain(:,jlat,nlev) -zflucw(:,nlev)*qcw(:,jlat,nlev)
      snowls(:,jlat) = -zflusn(:,nlev)*qsnow(:,jlat,nlev) -zfluci(:,nlev)*qci(:,jlat,nlev)  &
                       -zfluha(:,nlev)*qhail(:,jlat,nlev)

!  melting of residual snow and graupel/hail at the lowest level, after fall
!  latent heat of fusion subtracted from the lowest level

      do jlon = 2, nlon-1
      if (t(jlon,jlat,nlev).gt.tzer.and.t(jlon,jlat,nlev).lt.tzer+10.) then
       zt = t(jlon,jlat,nlev)
       zt0t = tzer/zt
       zesk = ezer*exp(-ccw1*log(zt0t)+ccw2*(1.-zt0t)) ! saturated pressure over water (in Pa)
       zqsw = zesk*eps/(ps(jlon,jlat)+zesk*(eps-1.))
       zrh = q(jlon,jlat,nlev)/zqsw ! relat. humidity
        if(zrh.ge.1.) then
         ztw = zt
        else
         ztwfg = zt + (0.04*(zt-276.)+1.)*(8.4 + 5.6e-5*(800.e2 - ps(jlon,jlat)))*(zrh - 1.) ! 1st guess of wet bulb t
         zt0tw = tzer/ztwfg
         zesk = ezer*exp(-ccw1*log(zt0tw)+ccw2*(1.-zt0tw)) ! saturated pressure over water (in Pa)
         zzqsw = zesk*eps/(ps(jlon,jlat)+zesk*(eps-1.))
         ztw = zt + yliv/cpd*(q(jlon,jlat,nlev)-zzqsw) ! wet bulb temperature for ice
         ztw = 0.5*ztw + 0.5*ztwfg                     ! empirical approx.
        endif
       if (ztw.gt.tzer+0.7) then  ! offset value empirical but important for snowfall
       zcapt = g*(pzer*dsig(nlev)-(pzer-ps(jlon,jlat))*dsigalf(nlev))*cpd
       zmelmax = (ztw-tzer)*zcapt/yliw
       zmeltsn = min (snowls(jlon,jlat), zmelmax)
       rainls(jlon,jlat) = rainls(jlon,jlat) + zmeltsn
       t(jlon,jlat,nlev) = t(jlon,jlat,nlev) - yliw*zmeltsn/zcapt
       snowls(jlon,jlat) = max(snowls(jlon,jlat) - zmeltsn, 0.)
       endif
      endif
      enddo

!-----------------------------------------------------------------------
! elimination of residuals of hydrometeors
!-----------------------------------------------------------------------

      do jklev = ntop, nlev
      do jlon = 2, nlon-1
      if(qcw  (jlon,jlat,jklev).lt.1.e-9) qcw  (jlon,jlat,jklev)=0.
      if(qci  (jlon,jlat,jklev).lt.1.e-9) qci  (jlon,jlat,jklev)=0.
      if(qrain(jlon,jlat,jklev).lt.1.e-8) qrain(jlon,jlat,jklev)=0.
      if(qsnow(jlon,jlat,jklev).lt.1.e-8) qsnow(jlon,jlat,jklev)=0.
      if(qhail(jlon,jlat,jklev).lt.1.e-8) qhail(jlon,jlat,jklev)=0.
      enddo
      enddo

 1000 continue

    return

    CONTAINS
    function eugamma (x)
    implicit none
    integer :: nx, ix, nfrx
    real :: eugamma, x, gfact, xx, gfrac
    real :: eug(0:100)
    data eug/1.00000,0.99433,0.98884,0.98355,0.97844,0.97350, &    !1.0 -1.05
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

    end subroutine micro
#ifndef rad_ecmwf
!##################################################################################################################
subroutine radiat_init(rds, ret, ozon, aerosol, aerotot, &
   nyrc, nmonc, ndayc, nhouc, nminc, &
   nlon, nlat, nlev, ntaer, sig, sigint, sigalf, pzer, dlat, dlon, &
   alat, alon, ps, tsurf, t, htop, infy, infx, proc_id)

! Initialization of astronomical parameters and
! total aerosols content for Ritter-Geleyn radiation scheme;
! ozone content and content of single type aerosol are defined also.
!
! Version for Bolam/Globo
!

use module_radiat_param

implicit none

! Input and output variables

integer, intent(in) :: nyrc, nmonc, ndayc, nhouc, nminc, nlon, nlat, nlev, ntaer, infy, infx, proc_id
real, intent(in) :: pzer, dlat, dlon
real, dimension(nlev+1), intent(in) :: sig
real, dimension(nlev), intent(in) :: sigint, sigalf
real, dimension(nlon,nlat), intent(in) :: alat, alon, ps, tsurf, htop
real, dimension(nlon,nlat,nlev), intent(in) :: t

real, intent(out) :: rds, ret
real, dimension(nlon,nlat,nlev), intent(out) :: ozon, aerotot
real, dimension(nlon,nlat,nlev,ntaer), intent(out) :: aerosol

! Internal variables

integer :: nsss, nzzaa, nzzmm, iiyr, nindat, iminut, jlon, jlat, jklev, jf, jk
real :: julian_day, rtime, rteta, rel, rem, rlls, rllls, zangozc

real, dimension(nlon,kaer,nlev) :: paer5
real, dimension(nlon,nlev) :: pap5, pozon5
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

! zeta and zetah: vertical sigma coordin. for the routines that set parameters for aerosol

   do jklev = 1,nlev
    zeta (jklev) = sigint(jklev)
    zetah(jklev) = sig(jklev)
   enddo
   zetah(nlev+1) = 1.

  loop_latitude_index: do jlat =1,nlat

! Geographical coordinate parameters

    do jlon = 1, nlon
#ifndef globo
      plat5(jlon)  = alat(jlon,jlat) * zdegrad
      plon5(jlon)  = alon(jlon,jlat) * zdegrad  ! long. < 0 not accepted
#else
      plat5(jlon)  = (-90. + real(infy+jlat-2)*dlat) * zdegrad
      plon5(jlon)  = (       real(infx+jlon-2)*dlon) * zdegrad
#endif
      if (plon5(jlon) < 0.) plon5(jlon) = plon5(jlon) + 2.*rpi
      pgelam5(jlon)= plon5(jlon)
      pgemu5(jlon) = sin(plat5(jlon))
      pclon5(jlon) = cos(plon5(jlon))
      pslon5(jlon) = sin(plon5(jlon))
    enddo

! Pressure at full levels

   do jklev = 1,nlev
     do jlon = 1,nlon
       pap5(jlon,jklev) = pzer*sigint(jklev) - (pzer-ps(jlon,jlat))*sigalf(jklev)
     enddo
   enddo

! Pressure at half levels: paph5(nlon,nlev+1)

    do jlon = 1,nlon
      paph5(jlon,1     ) = 1.e-8
      paph5(jlon,nlev+1) = ps(jlon,jlat)
    enddo
    do jklev = 2,nlev
      do jlon = 1,nlon
        paph5(jlon,jklev) = 0.5*(pap5(jlon,jklev-1)+pap5(jlon,jklev))
      enddo
    enddo

! Temperature at half-levels: at top defined as lev. 1
! at bottom as tsurf, elsewhere is the arithmetic mean

    do jlon = 1,nlon
      pth5(jlon,1) = t(jlon,jlat,1)
      pth5(jlon,nlev+1) = tsurf(jlon,jlat)
    enddo

    do jklev = 2,nlev
      do jlon = 1,nlon
        pth5(jlon,jklev) = 0.5*(t(jlon,jlat,jklev-1)+t(jlon,jlat,jklev))
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

#ifndef globo

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

#endif

  enddo loop_latitude_index

  do jf = 1, 2
    call filt2t (ozon, 1.)
  enddo
  ozon = max (ozon, 0.)
  do jf = 1, 5
    do jk = 1, 6
      call filt2t (aerosol(1:nlon,1:nlat,1:nlev,jk), 1.)
    enddo
  enddo
  do jf = 1, 60
    call filt2t (aerosol(1:nlon,1:nlat,1:nlev, 2), 1.)
  enddo
  aerosol = max (aerosol, zepaer)

#ifdef globo
  call polavert (ozon)
  ozon = max (ozon, 0.)
  do jk = 1, 6
    call polavert (aerosol(1:nlon,1:nlat,1:nlev,jk))
    aerosol(:,:,:,jk) = max (aerosol(:,:,:,jk), sngl(zepaer))
  enddo
#endif

  aerotot(:,:,:) = aerosol(:,:,:,1) + aerosol(:,:,:,2) + aerosol(:,:,:,3) +  &
                   aerosol(:,:,:,4) + aerosol(:,:,:,5) + aerosol(:,:,:,6)


return
end subroutine radiat_init
#endif
!##################################################################################################################
      subroutine radint (ndayr, nhouc, nminc, infx, rrf)

!  Interface to radial

      use mod_model, only : nlon, nlat, nlev, nlonm1, nlevp1, sig, sigint, dsig, ps, t, q, qcw, qci, emisg1,  &
                            fsnow, albedo, alsn, fmask, tskin, qskin, rd, g, pzer, dtstep, cpd, cpv, pi,      &
                            geldt, gelvis, gelirr, sigalf, dsigalf, alfa, qrain, qsnow, alont, alatt, nlatm1, &
                            myid, aerotot, fice, fcloud, cloudt, dlon, nprocsy, hst, rdecli0, reqtim0, myid
      implicit none
      common /radiazio/ stefan, ro, xts, gg, sdec, cdec, zalpn(nlon),                        &
              qfs(nlon,nlev), tfs(nlon,nlev), dp(nlon,nlev), pm(nlon,nlev), p(nlon,nlev),    &
              zqli(nlon,nlev), zqice(nlon,nlev), zneb(nlon,nlev), cph(nlon,nlevp1),          &
              daer(nlon,nlev), tsf(nlon), alb(nlon), emis(nlon), zmu0(nlon),                 &
              dtfr(nlon,nlev), rat(nlon), rg(nlon), rirs(nlon), rvis(nlon)
      real(kind=8) dtfr, rat, rg, rirs, rvis, stefan, ro, xts, gg, sdec, cdec, qfs, tfs, dp, &
                   pm, p, zqli, zqice, zneb, cph, daer, tsf, alb, zalpn, emis, zmu0, zslat, zlont
      real(4) aer(nlev), zarg, zgmt, zqhl, zsgalf, zalsn, zm, zhard, rrf
      integer jlon, jlat, jklev, ndayr, nhouc, nminc, infx, jstart, jend

!  Definition of radiation parameters

      ro = 1365.       ! ECMWF rad. value (Geleyn orig.: ro = 1373.)
      ro = ro*0.982    ! further tuning vs. ECMWF vis. radiation at surface
      stefan = 5.67e-8
      gg     = g

!  Seasonal variation of solar constant (factor - Pielke, p. 211)

      zarg = 2.*pi*float(ndayr-1)/365.
      xts  = 1.00011 + .034221*cos(zarg) + .00128*sin(zarg) + .000719*cos(2.*zarg) + .000077*sin(2.*zarg)

!  Solar declination (rdecli defined in aerdef (first step) and radintec)

      sdec = sin(rdecli0)
      cdec = cos(rdecli0)

!  Vertical profile of aerosol content (aer is a vertical integral; daer is aerosol on integer levels)
!  (case in which ECMWF aerosol is not used)

      jstart = 2
      jend = nlatm1

#ifdef globo
      if (mod(myid,nprocsy).eq.0) jstart = 1
      if (mod(myid,nprocsy).eq.nprocsy-1) jend = nlat
#endif

!  Loop over all latitudes

      do 1000 jlat = jstart, jend

      do jklev = 1, nlev
      do jlon = 2, nlonm1
      daer(jlon,jklev) = 1.08*aerotot(jlon,jlat,jklev) ! tuning vs. ECMWF vis. radiation at surf.
      enddo
      enddo

!  Temp., spec. humidity, water and ice clouds
!  Pressure intervals between half-levels and pressure at integer levels
!  Fraction of cloud cover fcloud put in zneb

      do jklev = 1, nlev
      do jlon = 2, nlonm1
      tfs  (jlon,jklev) = t  (jlon,jlat,jklev)
      qfs  (jlon,jklev) = q  (jlon,jlat,jklev)
      zqli (jlon,jklev) = qcw(jlon,jlat,jklev)
      zqice(jlon,jklev) = qci(jlon,jlat,jklev) + 0.2*qsnow(jlon,jlat,jklev)
      dp   (jlon,jklev) = pzer*dsig(jklev) - (pzer-ps(jlon,jlat))*dsigalf(jklev)
      pm   (jlon,jklev) = pzer*sigint(jklev) - (pzer-ps(jlon,jlat))*sigalf(jklev)
      zneb (jlon,jklev) = fcloud(jlon,jlat,jklev)
      enddo
      enddo

!  Specific heat of moist air at half-levels

      do jklev = 2, nlev
      do jlon = 2, nlonm1
      zqhl = .5*(q(jlon,jlat,jklev) + q(jlon,jlat,jklev-1))
      cph(jlon,jklev) = cpd*(1.-zqhl) + cpv*zqhl
      enddo
      enddo
      do jlon = 2, nlonm1
      cph(jlon,1)      = cpd
      cph(jlon,nlevp1) = cpd*(1.-qskin(jlon,jlat)) + cpv*qskin(jlon,jlat)
      enddo

!  Pressure at half-levels - not defined at the uppermost half-level

      do jklev = 2, nlevp1
      zsgalf = sig(jklev)**3 * (1.+alfa*(1.-sig(jklev))*(2.-sig(jklev)))
      do jlon = 2, nlonm1
      p(jlon,jklev-1) = pzer*sig(jklev) - (pzer-ps(jlon,jlat))*zsgalf
      enddo
      enddo

!  Surface temperature, albedo and emissivity
!  Change of snow albedo (alb) in a range depending on skin temperature (as in radintec)

      do jlon = 2, nlonm1
      tsf(jlon) = tskin(jlon,jlat)
!      if (tsf(jlon).lt.277.) then
!      zalsn = alsn-.1 + .2*(277.-tsf(jlon))/14.
!      zalsn = min (zalsn, alsn+.1)
!      else
!      zalsn = alsn-.1
!      endif
      zalsn = alsn(jlon,jlat) ! test
      alb(jlon)  = albedo(jlon,jlat)*(1.-fsnow(jlon,jlat)) + zalsn*fsnow(jlon,jlat)
      emis(jlon) = emisg1(jlon,jlat)*(1.-fsnow(jlon,jlat)) + 0.87* fsnow(jlon,jlat)
      enddo

!  True solar hour (units: hours and fractions of hours)
!  Reqtim0 (in sec) accounts for astron. time; it is computed in aerdef and radintec;
!  it is required to assure syncronization with ECMWF solar radiation
!  Albedo: for parallel (direct) radiation, dependency on solar zenith angle is computed
!  over land (Yang et al) and over sea (Taylor et al, as ECMWF)

      zgmt = float(nhouc) + float(nminc)/60. + reqtim0/3600.

#ifdef globo
      zslat = hst(jlat)
#endif

      do jlon = 2, nlonm1

#ifdef globo
      zlont = (infx+jlon-2)*dlon
#else
      zslat = sin(alatt(jlon,jlat)*pi/180.)
      zlont = alont(jlon,jlat)
#endif

      zmu0(jlon) = zslat*sdec-sqrt(1.-zslat**2)*cdec*cos(.2618993878d0*(zgmt+zlont*24./360.))
      zm = max(zmu0(jlon), 1.e-12)
      if (fmask(jlon,jlat).lt.0.5) then
      zalpn(jlon) = min(1., alb(jlon)*(1. + 0.26)/(1. + 2.*0.26*zm))
      else
      zhard = 1. - fmask(jlon,jlat) + fice(jlon,jlat)*fmask(jlon,jlat)
      zalpn(jlon) = (1.-zhard)*0.037/(1.1*zm**1.4 + 0.15) + alb(jlon)*zhard
      endif
      enddo

!  Call of radiation subrout.

      call radial

!  Output variables

      do jklev = 1, nlev
      do jlon = 2, nlonm1
      geldt(jlon,jlat,jklev) = sngl(dtfr(jlon,jklev))
      enddo
      enddo

!  Surface fluxes of visible and infrared radiation (positive downward)

      do jlon = 2, nlonm1
      gelvis(jlon,jlat) = sngl(rg(jlon))
!      gelirr(jlon,jlat) = max (sngl(rat(jlon)), sngl(rat(jlon))*0.855) ! tuning vs. ECMWF ir radiation at surface
      gelirr(jlon,jlat) = max (sngl(rat(jlon)), sngl(rat(jlon))*0.855*0.92*(1.-0.05*cloudt(jlon,jlat))) ! more tuning!

      if(rrf.gt.1.e-2) then ! Initial tapering
      gelvis(jlon,jlat) = gelvis(jlon,jlat)*(1. - rrf*cloudt(jlon,jlat)**.5)
      gelirr(jlon,jlat) = gelirr(jlon,jlat)*(1. - rrf*cloudt(jlon,jlat)**.5)
      endif

      enddo

1000  continue

      return
      end subroutine radint
!##################################################################################################################
      subroutine radial

! Ritter -Geleyn radiation scheme (modified)

      use mod_model, only: nlon, nlonm1, nlev, nlevp1, co2ppm
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
!##################################################################################################################
    subroutine filt2t (p, anu)

!  Filter of 'T points' 3D fields  -  it is assumed that ghostlines are not updated

    use mod_model, only: nlon, nlat, nlev, nlonm1, nlatm1, nprocsx, ip_e, ip_n, ip_s, ip_w, ip_oppo, ip_null
    real(4) p(nlon,nlat,nlev), p2(nlon,0:nlat+1,nlev)

    if (anu.lt.1.e-6) return

    jstart = 2
    jend = nlatm1

#ifdef globo
    if (ip_s.eq.ip_null) jstart = 1
    if (ip_n.eq.ip_null) jend = nlat
#endif

!------------------------
!  Update all ghostlines
!------------------------

#ifdef mpi
      call u_ghost (p(2:nlonm1,nlatm1,:), ip_n, p(2:nlonm1,1   ,:), ip_s, (nlon-2)*nlev)
      call u_ghost (p(2:nlonm1,2     ,:), ip_s, p(2:nlonm1,nlat,:), ip_n, (nlon-2)*nlev)
      if (nprocsx.eq.1) then
#endif

#ifdef globo
    p(1   ,:,:) = p(nlonm1,:,:)
    p(nlon,:,:) = p(2     ,:,:)
#endif

#ifdef mpi
      else
        call u_ghost (p(nlonm1,:,:), ip_e, p(1   ,:,:), ip_w, nlat*nlev)
        call u_ghost (p(2     ,:,:), ip_w, p(nlon,:,:), ip_e, nlat*nlev)
      endif
#endif

!------------------------
!  Horizontal diffusion
!------------------------

    do k = 1, nlev
    do jlat = 1, nlat
    do jlon = 2, nlonm1
    p2(jlon,jlat,k) = .25*(p(jlon-1,jlat,k)+p(jlon+1,jlat,k))+.5*p(jlon,jlat,k)
    enddo
    enddo
    enddo

#ifdef globo

#ifdef mpi
      if (ip_s.eq.ip_null) then
        if (nprocsx.eq.1) then
#endif
    do jlon = 2, nlonm1
    jlon1 = jlon + (nlon-2)/2
    if (jlon1.gt.nlon) jlon1 = jlon1-nlon+2
    p2(jlon,0,:) =  p2(jlon1,2,:)
    enddo
#ifdef mpi
        else
          call u_ghost (p2(2:nlonm1,2,:), ip_oppo, p2(2:nlonm1,0,:), ip_oppo, (nlon-2)*nlev)
        endif
      endif
      if (ip_n.eq.ip_null) then
        if (nprocsx.eq.1) then
#endif
    do jlon = 2, nlonm1
    jlon1 = jlon + (nlon-2)/2
    if (jlon1.gt.nlon) jlon1 = jlon1-nlon+2
    p2(jlon,nlat+1,:) =  p2(jlon1,nlatm1,:)
    enddo
#ifdef mpi
        else
          call u_ghost (p2(2:nlonm1,nlatm1,:), ip_oppo, p2(2:nlonm1,nlat+1,:), ip_oppo, (nlon-2)*nlev)
        endif
      endif
#endif

#endif ! globo

    do k = 1, nlev
    do jlat = jstart, jend
    do jlon = 2, nlonm1
    p(jlon,jlat,k)=(1.-anu)*p(jlon,jlat,k)+anu*(.25*(p2(jlon,jlat+1,k)+p2(jlon,jlat-1,k))+.5*p2(jlon,jlat,k))
    enddo
    enddo
    enddo

    return
    end subroutine filt2t
!##################################################################################################################
    subroutine filt2t1 (p, anu)

!  Filter of 'T-points' 2D fields  -  it is assumed that ghostlines are not updated

    use mod_model, only: nlon, nlat, nlonm1, nlatm1, nprocsx, ip_e, ip_n, ip_s, ip_w, ip_oppo, ip_null
    real(4) p(nlon,nlat), p2(nlon,0:nlat+1)

    if (anu.lt.1.e-6) return

    jstart = 2
    jend = nlatm1

#ifdef globo
    if (ip_s.eq.ip_null) jstart = 1
    if (ip_n.eq.ip_null) jend = nlat
#endif

!------------------------
!  Update all ghostlines
!------------------------

#ifdef mpi
      call u_ghost (p(2:nlonm1,nlatm1), ip_n, p(2:nlonm1,1   ), ip_s, nlon-2)
      call u_ghost (p(2:nlonm1,2     ), ip_s, p(2:nlonm1,nlat), ip_n, nlon-2)
      if (nprocsx.eq.1) then
#endif

#ifdef globo
    p(1   ,:) = p(nlonm1,:)
    p(nlon,:) = p(2     ,:)
#endif

#ifdef mpi
      else
        call u_ghost (p(nlonm1,:), ip_e, p(1   ,:), ip_w, nlat)
        call u_ghost (p(2     ,:), ip_w, p(nlon,:), ip_e, nlat)
      endif
#endif

!------------------------
!  Horizontal diffusion
!------------------------

    do jlat = 1, nlat
    do jlon = 2, nlonm1
    p2(jlon,jlat) = .25*(p(jlon-1,jlat)+p(jlon+1,jlat))+.5*p(jlon,jlat)
    enddo
    enddo

#ifdef globo

#ifdef mpi
      if (ip_s.eq.ip_null) then
        if (nprocsx.eq.1) then
#endif
    do jlon = 2, nlonm1
    jlon1 = jlon + (nlon-2)/2
    if (jlon1.gt.nlon) jlon1 = jlon1-nlon+2
    p2(jlon,0) =  p2(jlon1,2)
    enddo
#ifdef mpi
        else
          call u_ghost (p2(2:nlonm1,2), ip_oppo, p2(2:nlonm1,0), ip_oppo, nlon-2)
        endif
      endif
      if (ip_n.eq.ip_null) then
        if (nprocsx.eq.1) then
#endif
    do jlon = 2, nlonm1
    jlon1 = jlon + (nlon-2)/2
    if (jlon1.gt.nlon) jlon1 = jlon1-nlon+2
    p2(jlon,nlat+1) =  p2(jlon1,nlatm1)
    enddo
#ifdef mpi
        else
          call u_ghost (p2(2:nlonm1,nlatm1), ip_oppo, p2(2:nlonm1,nlat+1), ip_oppo, nlon-2)
        endif
      endif
#endif

#endif ! globo

    do jlat = jstart, jend
    do jlon = 2, nlonm1
    p(jlon,jlat) = (1.-anu)*p(jlon,jlat)+anu*(.25*(p2(jlon,jlat+1)+p2(jlon,jlat-1))+.5*p2(jlon,jlat))
    enddo
    enddo

    return
    end subroutine filt2t1
!##################################################################################################################
    subroutine filt2uv (anu)

    use mod_model, only: nlon, nlat, nlev, nlonm1, nlatm1, nprocsx, nprocsy, dlon, pi, snt, cst, u, v, &
                         myid, ip_e, ip_n, ip_s, ip_w, ip_oppo, ip_null
    real(4) zu2(nlon,nlat), zv2(nlon,nlat+1,nlev), zvi(2,nlev), zvo(2,nlev)
    integer :: ierr, comm, tag1=1, tag2=2

#ifdef mpi
      include 'mpif.h'
      integer :: status(mpi_status_size)

      comm  = mpi_comm_world
#endif

    if (anu.lt.1.e-6) return

    iprocs = nprocsx*nprocsy
    jendv = nlatm1

#ifdef globo
    if (ip_n.eq.ip_null) jendv = nlat
#endif

!------------------------
!  Update all ghostlines
!------------------------

#ifdef mpi
      call u_ghost (u(2:nlonm1,nlatm1,:), ip_n, u(2:nlonm1,1   ,:), ip_s, (nlon-2)*nlev)
      call u_ghost (u(2:nlonm1,2     ,:), ip_s, u(2:nlonm1,nlat,:), ip_n, (nlon-2)*nlev)
      call u_ghost (v(2:nlonm1,nlatm1,:), ip_n, v(2:nlonm1,1   ,:), ip_s, (nlon-2)*nlev)
      call u_ghost (v(2:nlonm1,2     ,:), ip_s, v(2:nlonm1,nlat,:), ip_n, (nlon-2)*nlev)

      if (nprocsx.eq.1) then
#endif

#ifdef globo
    u(1   ,:,:) = u(nlonm1,:,:)
    u(nlon,:,:) = u(2     ,:,:)
    v(1   ,:,:) = v(nlonm1,:,:)
    v(nlon,:,:) = v(2     ,:,:)
#endif

#ifdef mpi
      else
        call u_ghost (u(nlonm1,:,:), ip_e, u(1   ,:,:), ip_w, nlat*nlev)
        call u_ghost (u(2     ,:,:), ip_w, u(nlon,:,:), ip_e, nlat*nlev)
        call u_ghost (v(nlonm1,:,:), ip_e, v(1   ,:,:), ip_w, nlat*nlev)
        call u_ghost (v(2     ,:,:), ip_w, v(nlon,:,:), ip_e, nlat*nlev)
      endif
#endif

!----------------------------------
!  Second order diffusion of wind
!----------------------------------

! 2-grid-interval filter of v over the whole domain

    do k = 1, nlev
    do jlat = 1, nlat
    do jlon = 2, nlonm1
    zv2(jlon,jlat,k) = .25*(v(jlon-1,jlat,k)+v(jlon+1,jlat,k))+.5*v(jlon,jlat,k)
    enddo
    enddo
    enddo

#ifdef globo

#ifdef mpi
      if (ip_s.eq.ip_null) then
        if (nprocsx.eq.1) then
#endif
    do jlon = 2, nlonm1
    jlon1 = jlon + (nlon-2)/2
    if (jlon1.gt.nlon) jlon1 = jlon1-nlon+2
    zv2(jlon,1,:) = -zv2(jlon1,2,:)
    enddo
#ifdef mpi
        else
          call u_ghost (-zv2(2:nlonm1,2,:), ip_oppo, zv2(2:nlonm1,1,:), ip_oppo, (nlon-2)*nlev)
        endif
      endif
#endif

#ifdef mpi
      if (ip_n.eq.ip_null) then
        if (nprocsx.eq.1) then
#endif
    do jlon = 2, nlonm1
    jlon1 = jlon + (nlon-2)/2
    if (jlon1.gt.nlon) jlon1 = jlon1-nlon+2
    zv2(jlon,nlat+1,:) =  -zv2(jlon1,nlat,:)
    enddo
#ifdef mpi
        else
          call u_ghost (-zv2(2:nlonm1,nlat,:), ip_oppo, zv2(2:nlonm1,nlat+1,:), ip_oppo, (nlon-2)*nlev)
        endif
      endif
#endif

#endif ! globo

    do k = 1, nlev
    do jlat = 2, jendv
    do jlon = 2, nlonm1
    v(jlon,jlat,k) = (1.-anu)*v(jlon,jlat,k)+anu*(.25*(zv2(jlon,jlat+1,k)+zv2(jlon,jlat-1,k))+.5*zv2(jlon,jlat,k))
    enddo
    enddo
    enddo

! u at poles

#ifdef globo

#ifdef mpi
      if (ip_s.eq.ip_null) then   ! row of processors including south pole
#endif

    zvi = 0.
    do k = 1, nlev
    do jlon = 2, nlonm1
    zvi(1,k) = zvi(1,k) + v(jlon,2,k)*snt(jlon,1)
    zvi(2,k) = zvi(2,k) + v(jlon,2,k)*cst(jlon,1)
    enddo
    enddo
    zvi = zvi*dlon/180.
#ifdef mpi
        if (myid.gt.0) then
          call mpi_send (zvi, 2*nlev, mpi_real, 0, tag1, comm, ierr)
          call mpi_recv (zvi, 2*nlev, mpi_real, 0, tag2, comm, status, ierr)
        else
          do jpr = nprocsy, iprocs-1, nprocsy
            call mpi_recv (zvo, 2*nlev, mpi_real, jpr, tag1, comm, status, ierr)
            zvi = zvi + zvo
          enddo
          do jpr = nprocsy, iprocs-1, nprocsy
            call mpi_send (zvi, 2*nlev, mpi_real, jpr, tag2, comm, ierr)
          enddo
        endif
#endif

    do k = 1, nlev
    u(:,1,k) = zvi(1,k)*cst(:,1) - zvi(2,k)*snt(:,1)
    enddo

#ifdef mpi
      endif ! ip_s.eq.ip_null
#endif

#ifdef mpi
      if (ip_n.eq.ip_null) then  !  row of processors including north pole
#endif

    zvi = 0.
    do k = 1, nlev
    do jlon = 2, nlonm1
    zvi(1,k) = zvi(1,k) + v(jlon,nlat,k)*snt(jlon,1)
    zvi(2,k) = zvi(2,k) + v(jlon,nlat,k)*cst(jlon,1)
    enddo
    enddo
    zvi = zvi*dlon/180.
#ifdef mpi
        if (myid.gt.nprocsy-1) then
          call mpi_send (zvi, 2*nlev, mpi_real, nprocsy-1, tag1, comm, ierr)
          call mpi_recv (zvi, 2*nlev, mpi_real, nprocsy-1, tag2, comm, status, ierr)
        else
          do jpr = 2*nprocsy-1, iprocs-1, nprocsy
            call mpi_recv (zvo, 2*nlev, mpi_real, jpr, tag1, comm, status, ierr)
            zvi = zvi + zvo
          enddo
          do jpr = 2*nprocsy-1, iprocs-1, nprocsy
            call mpi_send (zvi, 2*nlev, mpi_real, jpr, tag2, comm, ierr)
          enddo
        endif
#endif
    do k = 1, nlev
    u(:,nlat,k) = -zvi(1,k)*cst(:,1) + zvi(2,k)*snt(:,1)
    enddo

#ifdef mpi
      endif ! ip_n.eq.ip_null
#endif

#endif ! globo

! 2-grid-interval filter of u over the whole domain

    do k = 1, nlev
    do jlat = 1, nlat
    do jlon = 2, nlonm1
    zu2(jlon,jlat) = .25*(u(jlon-1,jlat,k)+u(jlon+1,jlat,k))+.5*u(jlon,jlat,k)
    enddo
    enddo
    do jlat = 2, nlatm1
    do jlon = 2, nlonm1
    u(jlon,jlat,k) = (1.-anu)*u(jlon,jlat,k) +anu*(.25*(zu2(jlon,jlat+1)+zu2(jlon,jlat-1))+.5*zu2(jlon,jlat))
    enddo
    enddo
    enddo

    return
    end subroutine filt2uv
!##################################################################################################################
    subroutine filt (p, nt, anu2)

! Horizontal diffusion (5 points)

    use mod_model, only: nlon, nlat, nlev, nlonm1, nlatm1, ip_e, ip_n, ip_s, ip_w
    implicit none

    integer :: nt, jlon, jlat, k, ierr, comm, tag1=1, tag2=2, isend1, isend2, irecv1, irecv2
    real(4) p(nlon,nlat,nlev), p2(nlon,nlat), anu2
#ifdef mpi
      include 'mpif.h'
      integer :: status(mpi_status_size)

    comm  = mpi_comm_world
#endif

!------------------------
!  Update all ghostlines
!------------------------

#ifdef mpi
      do k = 1, nt
!  qui non uso la u_ghost perche' perde un sacco di tempo... perche'?
!        call u_ghost (p(2:nlonm1,nlatm1,k), ip_n, p(2:nlonm1,1   ,k), ip_s, nlon-2)
!        call u_ghost (p(2:nlonm1,2     ,k), ip_s, p(2:nlonm1,nlat,k), ip_n, nlon-2)
        call mpi_isend (p(2:nlonm1,nlatm1,k), nlon-2, mpi_real, ip_n, tag1, comm, isend1, ierr)
        call mpi_isend (p(2:nlonm1,2     ,k), nlon-2, mpi_real, ip_s, tag2, comm, isend2, ierr)
        call mpi_irecv (p(2:nlonm1,1     ,k), nlon-2, mpi_real, ip_s, tag1, comm, irecv1, ierr)
        call mpi_irecv (p(2:nlonm1,nlat  ,k), nlon-2, mpi_real, ip_n, tag2, comm, irecv2, ierr)
        call mpi_wait  (isend1, status, ierr)
        call mpi_wait  (isend2, status, ierr)
        call mpi_wait  (irecv1, status, ierr)
        call mpi_wait  (irecv2, status, ierr)
        call u_ghost (p(nlonm1,:,k), ip_e, p(1   ,:,k), ip_w, nlat)
        call u_ghost (p(2     ,:,k), ip_w, p(nlon,:,k), ip_e, nlat)
      enddo
#endif

!-----------------------
!  Horizontal diffusion
!-----------------------

    do 50 k = 1, nt
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
    end subroutine filt
!##################################################################################################################
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
!##################################################################################################################
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
!##################################################################################################################
    subroutine tridiag (za, zc, psisurf, psi, nlon, nlev, ntop)

! Computes tridiagonal matrix inversion

    implicit none
    real(4) za(nlon,nlev), zc(nlon,nlev), psi(nlon,nlev), psisurf(nlon)
    integer nlon, nlev, ntop, jlon, jklev
    real(4) zrb, zb, zrden, ze(nlon,nlev)

!-------------------------------------------------
!  tridiagonal matrix inversion:                 -
!  a(k)*psi(k-1)+b(k)*psi(k)+c(k)*psi(k+1)=r(k)  -
!  b(k) = 1-a(k)-c(k)                            -
!  a(ntop) = 0                                   -
!  r(k) = psim(k),   ntop <= k < nlev            -
!  r(nlev) = -c(nlev)*psisurf + psim(nlev)       -
!  where psi(k) = u,v,teta,q  at new time level  -
!  psim prognostic variables at old time level   -
!  and psisurf their values at the surface.      -
!  Solving formula: psi(k)=e(k)*psi(k+1)+f(k)    -
!  e(k) independent of psi                       -
!  f(k) function of psi saved in psi             -
!-------------------------------------------------

    do jlon = 2, nlon-1
    zrb=1./(1.-za(jlon,ntop)-zc(jlon,ntop))
    ze(jlon,ntop) = -zc(jlon,ntop)*zrb
    psi(jlon,ntop) = psi(jlon,ntop)*zrb
    enddo

    do jklev = ntop+1, nlev-1
    do jlon = 2, nlon-1
    zb = 1.-za(jlon,jklev)-zc(jlon,jklev)
    zrden = 1./(zb+za(jlon,jklev)*ze(jlon,jklev-1))
    ze(jlon,jklev) = -zc(jlon,jklev)*zrden
    psi(jlon,jklev) = (psi(jlon,jklev)-za(jlon,jklev)*psi(jlon,jklev-1))*zrden
    enddo
    enddo

    do jlon = 2, nlon-1
    zb = 1.-za(jlon,nlev)-zc(jlon,nlev)
    zrden = 1./(zb+za(jlon,nlev)*ze(jlon,nlev-1))
    psi(jlon,nlev) = ( psi(jlon,nlev)-za(jlon,nlev)*psi(jlon,nlev-1)-zc(jlon,nlev)*psisurf(jlon) )*zrden
    enddo

    do jklev = nlev-1, ntop, -1
    do jlon = 2, nlon-1
    psi(jlon,jklev) = ze(jlon,jklev)*psi(jlon,jklev+1) +psi(jlon,jklev)
    enddo
    enddo

    return
    end subroutine tridiag
!##################################################################################################################
    subroutine wafps

!  Computes divergence on t points of the 2-d mass flux with WAF scheme:
!  div2(2:nlonm1,2:nlatm1,1:nlev)

    use mod_model, only : nlon, nlat, nlev, nlonm1, nlatm1, myid, nprocsx, nprocsy, ip_oppo, cpole, &
                          ip_e, ip_n, ip_s, ip_w, ip_null, ps, u, v, hxt, hxv, dx, dy, dt, div2, denomin

    real(4) denrx(nlon,nlat), denry(nlon,nlat), waflux(nlon,nlat), waflux1(nlon), psex(0:nlon+1,0:nlat+1)
    real(4) zzz(nlev)
    integer :: ierr, comm, tag1=1, tag2=2, jpr, iproc

#ifdef mpi
    include 'mpif.h'
    integer :: status(mpi_status_size)

    iprocs = nprocsx*nprocsy
    comm  = mpi_comm_world
#endif

    rdt = 1./dt

    psex(1:nlon,1:nlat) = ps(1:nlon,1:nlat)    ! psex contains ps with two ghostlines

#ifndef globo

#ifdef mpi
      if (ip_w.eq.ip_null) then
        psex(0,1:nlat) = 2.*ps(1,1:nlat)-ps(2,1:nlat)
      endif
      if (ip_e.eq.ip_null) then
        psex(nlon+1,1:nlat) = 2.*ps(nlon,1:nlat)-ps(nlonm1,1:nlat)
      endif
      if (ip_s.eq.ip_null) then
        psex(0:nlon+1,0) = 2.*psex(0:nlon+1,1)-psex(0:nlon+1,2)
      endif
      if (ip_n.eq.ip_null) then
        psex(0:nlon+1,nlat+1) = 2.*psex(0:nlon+1,nlat)-psex(0:nlon+1,nlatm1)
      endif
#else
      psex(0,1:nlat) = 2.*ps(1,1:nlat)-ps(2,1:nlat)
      psex(nlon+1,1:nlat) = 2.*ps(nlon,1:nlat)-ps(nlonm1,1:nlat)
      psex(0:nlon+1,0) = 2.*psex(0:nlon+1,1)-psex(0:nlon+1,2)
      psex(0:nlon+1,nlat+1) = 2.*psex(0:nlon+1,nlat)-psex(0:nlon+1,nlatm1)
#endif

#endif ! globo

#ifdef mpi
      if (nprocsx.eq.1) then
#endif

#ifdef globo
    psex(0,1:nlat) = psex(nlon-2,1:nlat)
    psex(nlon+1,1:nlat) = psex(3,1:nlat)
#endif

#ifdef mpi
      else
        call u_ghost (psex(nlon-2,1:nlat), ip_e, psex(0     ,1:nlat), ip_w, nlat)
        call u_ghost (psex(3     ,1:nlat), ip_w, psex(nlon+1,1:nlat), ip_e, nlat)
      endif
      call u_ghost (psex(0:nlon+1,nlat-2), ip_n, psex(0:nlon+1,0     ), ip_s, nlon+2)
      call u_ghost (psex(0:nlon+1,3     ), ip_s, psex(0:nlon+1,nlat+1), ip_n, nlon+2)
#endif

#ifdef globo

#ifdef mpi
      if (ip_s.eq.ip_null) then
        if (nprocsx.eq.1) then
#endif
    do jlon = 0, nlon+1
    jlon1 = jlon + (nlon-2)/2
    if (jlon1.gt.nlon) jlon1 = jlon1-nlon+2
    psex(jlon,0) = psex(jlon1,2)
    enddo
#ifdef mpi
        else
          call u_ghost (psex(0:nlon+1,2), ip_oppo, psex(0:nlon+1,0), ip_oppo, nlon+2)
        endif
      endif

      if (ip_n.eq.ip_null) then
        if (nprocsx.eq.1) then
#endif
    do jlon = 0, nlon+1
    jlon1 = jlon + (nlon-2)/2
    if (jlon1.gt.nlon) jlon1=jlon1-nlon+2
    psex(jlon,nlat+1) = psex(jlon1,nlatm1)
    enddo
#ifdef mpi
        else
          call u_ghost (psex(0:nlon+1,nlatm1), ip_oppo, psex(0:nlon+1,nlat+1), ip_oppo, nlon+2)
        endif
      endif
#endif

#endif ! globo

    do jlat = 2, nlatm1
    do jlon = 2, nlon
!    denrx(jlon,jlat) = 1./(psex(jlon,jlat)-psex(jlon-1,jlat))
    denrx(jlon,jlat) = denomin(psex(jlon,jlat),psex(jlon-1,jlat))
    enddo
    enddo
    do jlat = 2, nlat
    do jlon = 2, nlonm1
!    denry(jlon,jlat) = 1./(psex(jlon,jlat)-psex(jlon,jlat-1))
    denry(jlon,jlat) = denomin(psex(jlon,jlat),psex(jlon,jlat-1))
    enddo
    enddo

    do 10 jklev = 1, nlev
    do jlat = 2, nlatm1
    zcost = dt/(dx*hxt(jlat))
    do jlon = 2, nlon
    amu = u(jlon-1,jlat,jklev)*zcost
        if(amu.ge.0.) then
        is=1
        j1=jlon-1
        else
        is=-1
        j1=jlon+1
        endif
    r = (psex(j1,jlat)-psex(j1-1,jlat))*denrx(jlon,jlat)
    b = max(0., min(2., max(r, min(2.*r,1.))))
    phi = is+amu*b-is*b
    waflux1(jlon) = .5*amu*((1.+phi)*psex(jlon-1,jlat) + (1.-phi)*psex(jlon,jlat))
    enddo
    do jlon = 2, nlonm1
    div2(jlon,jlat,jklev) = waflux1(jlon+1) - waflux1(jlon)
    enddo
    enddo
 10 continue

    zcost = dt/dy
    do 20 jklev = 1, nlev
    do jlat = 2, nlat
    do jlon = 2, nlonm1
    amu = v(jlon,jlat,jklev)*zcost
        if(amu.ge.0.) then
        is=1
        j1=jlat-1
        else
        is=-1
        j1=jlat+1
        endif
    r = (psex(jlon,j1)-psex(jlon,j1-1))*denry(jlon,jlat)
    b = max(0., min(2., max(r, min(2.*r,1.))))
    phi = is+amu*b-is*b
    waflux(jlon,jlat) = .5*amu*((1.+phi)*psex(jlon,jlat-1) + (1.-phi)*psex(jlon,jlat))
    enddo
    enddo
    do jlat = 2, nlatm1
    rhxt = 1./hxt(jlat)
    areap = hxv(jlat+1)*rhxt
    aream = hxv(jlat  )*rhxt
    do jlon = 2, nlonm1
    div2(jlon,jlat,jklev) = (div2(jlon,jlat,jklev) + waflux(jlon,jlat+1)*areap - waflux(jlon,jlat)*aream)*rdt
    enddo
    enddo

#ifdef globo

    zc1 = cpole*dy*rdt
#ifdef mpi
      if (ip_s.eq.ip_null) then
#endif
    div2(1,1,jklev) = 0.
    do jlon = 2, nlonm1
    div2(1,1,jklev) = div2(1,1,jklev) + waflux(jlon,2)*zc1
    enddo
#ifdef mpi
      endif
      if (ip_n.eq.ip_null) then
#endif
    div2(1,nlat,jklev) = 0.
    do jlon = 2, nlonm1
    div2(1,nlat,jklev) = div2(1,nlat,jklev) - waflux(jlon,nlat)*zc1
    enddo
#ifdef mpi
      endif
#endif

#endif ! globo

 20 continue

#ifdef globo

#ifdef mpi
      if (ip_s.eq.ip_null) then
        if (nprocsx.gt.1) then
          if (myid.ne.0) then
            call mpi_send (div2(1,1,1:nlev), nlev, mpi_real, 0, tag1, comm, ierr)
            call mpi_recv (div2(1,1,1:nlev), nlev, mpi_real, 0, tag2, comm, status, ierr)
          else
            do jpr = nprocsy, iprocs-1, nprocsy
              call mpi_recv (zzz, nlev, mpi_real, jpr, tag1, comm, status, ierr)
              div2(1,1,:) = div2(1,1,:) + zzz(:)
            enddo
            do jpr = nprocsy, iprocs-1, nprocsy
              call mpi_send (div2(1,1,1:nlev), nlev, mpi_real, jpr, tag2, comm, ierr)
            enddo
          endif
        endif
#endif
    do jklev = 1, nlev
    do jlon = 2, nlon
    div2(jlon,1,jklev) = div2(1,1,jklev)
    enddo
    enddo
#ifdef mpi
      endif ! ip_s.eq.ip_null
#endif

#ifdef mpi
      if (ip_n.eq.ip_null) then
        if (nprocsx.gt.1) then
          if (myid.ne.nprocsy-1) then
            call mpi_send (div2(1,nlat,1:nlev), nlev, mpi_real, nprocsy-1, tag1, comm, ierr)
            call mpi_recv (div2(1,nlat,1:nlev), nlev, mpi_real, nprocsy-1, tag2, comm, status, ierr)
          else
            do jpr = 2*nprocsy-1, iprocs-1, nprocsy
              call mpi_recv (zzz, nlev, mpi_real, jpr, tag1, comm, status, ierr)
              div2(1,nlat,:) = div2(1,nlat,:) + zzz(:)
            enddo
            do jpr = 2*nprocsy-1, iprocs-1, nprocsy
              call mpi_send (div2(1,nlat,1:nlev), nlev, mpi_real, jpr, tag2, comm, ierr)
            enddo
          endif
        endif
#endif
    do jklev = 1, nlev
    do jlon = 2, nlon
    div2(jlon,nlat,jklev) = div2(1,nlat,jklev)
    enddo
    enddo
#ifdef mpi
      endif ! ip_n.eq.ip_null
#endif

#endif ! globo

    return
    end subroutine wafps
!##################################################################################################################
      subroutine relax (is, gammin, gammax, alpha)

!  Computes optimal relaxation coefficients for bolam lateral boundary conditions
!  (Lehmann, MAP, 1993,1-14).
!  Input:  is       width of boundary relaxation zone (power of 2)
!          gammin   minimal courant number (c*dt/dx)
!          gammax   maximal courant number
!  Output: alpha() weight of externally specified values in the boundary
!          zone (corresponding to optimal relax. coefficients)

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
        write (*,*) " nbl is not a power of 2: stop in subrout. relax"
        stop
        endif

      do i=n,1,-1
      kk=p(i)/q(i-1)

      do j=i,1,-1
      xxx  = q(j)
      q(j) = p(j)-kk*q(j-1)
      p(j) = xxx
      enddo

      xxx  = q(0)
      q(0) = p(0)
      p(0) = xxx
      kdt2 = kk*sqrt(gammin*gammax)
      alpha(i) = kdt2/(1.+kdt2)
      enddo

!  remark: this alpha corresponds to the leapfrog scheme,
!  whereas kdt2 is independent of the integration scheme

      return
      end subroutine relax
!##################################################################################################################
      subroutine defclt (qcw, qci, t, nlon, nlat, nlev, tzer, nt)

! Defines cloud ice and cloud water assuming a ratio depending on T
! Limits max. values of cloud ice and cloud water (except for initial condition: nt=1)

      real qcw(nlon,nlat,nlev), qci(nlon,nlat,nlev), t(nlon,nlat,nlev)

      do jklev = 1, nlev
      do jlat = 1, nlat
      do jlon = 1, nlon
      ztemp = t(jlon,jlat,jklev)
      if(ztemp.ge.tzer) then
       zratio = 1.
      elseif (ztemp.lt.tzer-28.) then
       zratio = 0.
      else
       zratio = 1.04979*(0.5 + 0.5*tanh((ztemp-tzer+9.)/6.))
      endif
      zctot=qcw(jlon,jlat,jklev)
      qcw(jlon,jlat,jklev) = zratio*zctot
      qci(jlon,jlat,jklev) = zctot - qcw(jlon,jlat,jklev)
      if(nt.ne.1) then
       qcw(jlon,jlat,jklev) = min(qcw(jlon,jlat,jklev), 0.5e-3)  ! to prevent excessive precip. at boundaries
       qci(jlon,jlat,jklev) = min(qci(jlon,jlat,jklev), 0.08e-3) ! to prevent excessive precip. at boundaries
      endif
      enddo
      enddo
      enddo

      return
      end subroutine defclt
!##################################################################################################################
      subroutine convection_kf (ltopc)

! Kain-Fritsch (2004) convection (deep and shallow) scheme, revised

! Aug. 2018: max. timec decreased to avoid spurious explicit convection
! June 2012:
! changes mainly affecting shallow convection:
! dpmin lowered by 20% (also affects deep conv.); other corrections by M. Fantini (see "Mau")
! change in chmax (lowered from 1000 to 200 m to start shallow convection) - see "Andr"
! change in TIMEC (increased) for shallow convection - see "Andr"
! (the last two changes make shallow convection more effective on tropical oceans, with a shallow pbl,
! but at the same time limit the warming at about 850 hPa due to shallow conv. at mid-latitudes in summer,
! that seems excessive - perhaps shallow conv. algorithm should be revised, introducing iteration?

!  Oct. 2008 - parameterization of convection using conservation of liquid water static energy

      use mod_model, only : nlon, nlat, nlev, hxt, phig, ps, u, v, t, q, omeg, phi, sigint, sigalf, dsigalf, &
                            dsig, pzer, r=>rd, rv, ep, g, cp=>cpd, cv=>cvd, dx0=>dx, dy, t00=>tzer, cpv,     &
                            hst, snocon, raicon, dtdt, dqdt, dqcwdt, dqcidt, dqrndt, dqsndt, ntsrc, dtstep,  &
                            nlonm1, nlatm1, nprocsy, myid

#ifdef chem
      use mod_bol2chem, only: convflag, pconv, econv
#endif

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

#ifdef chem
      ! init pconv
      pptliq = 0.0
      pptice = 0.0
      convflag = .false.
      pconv = 0.0
      econv = 0.0
#endif

      jstart = 2
      jend   = nlatm1

#ifdef globo
      if (mod(myid,nprocsy).eq.0)         jstart = 7         ! no convection at poles
      if (mod(myid,nprocsy).eq.nprocsy-1) jend   = nlat-6    ! no convection at poles
#endif

!*******************************************************************
!                  environmental properties                        *
!*******************************************************************

      gdry   = -g/cp
      dts    = dtstep*ntsrc
      do 999 jlat = jstart, jend
      dx     = dx0*hxt(jlat)
      dxsq   = dx*dy
      do 999 jlon = 2, nlonm1

      ishall = 0
      dpmin  = 4.5e3                  ! thickness of layer whose properties characterize the parcel
      p300   = ps(jlon,jlat)-30000.   ! pressure at 300 mb above surface
      kznew  = nlev-ltopc+1

!  input:  temperature (t0, kelvin) ; specific humidity (q0, kg/kg) ;
!          horizontal wind speed (u0 and v0, m/s) ; pressure (p0, pascal) ;
!          height (z0, m);  vertical motion (w0, m/s)

      ml = 0
      l5 = 1
      do 15 k = 1, kznew
      nk = nlev-k+1
      p0(k)  = pzer*sigint(nk)- (pzer-ps(jlon,jlat))*sigalf(nk)
      dp(k)  = pzer*dsig(nk)  - (pzer-ps(jlon,jlat))*dsigalf(nk)  ! dp is the pressure interval between sigma levels
      t0(k)  = t(jlon,jlat,nk)
      q0(k)  = max(q(jlon,jlat,nk), 1.e-10)
      u0(k)  = .5*(u(jlon,jlat,nk)+u(jlon-1,jlat,nk))
      v0(k)  = .5*(v(jlon,jlat,nk)+v(jlon,jlat+1,nk))
      zes    = aliq*exp((bliq*t0(k)-cliq)/(t0(k)-dliq))
      zqes   = 0.622*zes/(p0(k)-zes)
      q0(k)  = min (zqes, q0(k))                ! if q0 is above saturation value, reduce it to saturation level
      rh0(k) = q0(k)/zqes
      tv0(k) = t0(k)*(1.+ep*q0(k))
      w0(k)  = -r/g*omeg(jlon,jlat,nk)*tv0(k)             ! omeg is omega/p
      z0(k)  = (phi(jlon,jlat,nk)-phig(jlon,jlat))/g
      thta0(k) = t0(k)*(1.e5/p0(k))**(0.2854*(1.-0.28*q0(k)))  ! theta environment
      ems (k)  = dp(k)*dxsq/g                   ! ems is mass in the box: rho*volume
      emsd(k)  = 1./ems(k)

      if (p0(k).ge.500e2) l5   = k
      if (p0(k).ge.p300 ) llfc = k     ! llfc is the last level below p300
      if (t0(k).gt.t00  ) ml   = k     ! ml is the highest lev. with t above zero - melting level
  15  continue

      cldhgt = 0.    ! cloud thickness is initialized to zero at all levels

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

#ifdef globo
      wkl = (w0(klcl-1) + (w0(klcl)-w0(klcl-1))*dlp) - wkl0*hst(jlat)**2
#else
      wkl = (w0(klcl-1) + (w0(klcl)-w0(klcl-1))*dlp)*dx/25.e3 - wkl0
#endif

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

#ifdef chem
      ! pconv is the level-by-level accumulated precipitation; vertical
      ! axis of pconv index is in bolam convention (inverted with
      ! respect to convection_kf)
      if( trppt > 0.0 )then ! update pconv only if there is some precipitation
        convflag( jlon,jlat ) = .true.
        ! TODO: verify that pconv(nlev-klcl+1) is equal to trppt
        pconv( jlon,jlat, nlev-kznew+1-1 ) = 0.0
        do nk = kznew, 1, -1 ! kznew: top of comput domain; must be defined up to the lowest model level!
          pconv( jlon,jlat, nlev-nk+1 ) = pconv( jlon,jlat, nlev-nk+1-1 ) + pptliq(nk) + pptice(nk)
        end do
      end if
#endif

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

#ifdef chem
      ! local (not cumulated) re-evaporation
      econv( jlon,jlat, nlev-k+1 ) = (qdd(k)-qd)*ddr(k)
      ! subtract cumulated level-by-leve downdraft evaporation rate from pconv at each level
      ! and set the result to 0 if negative (evaporation larger than precipitation)
      pconv( jlon,jlat, nlev-k+1 ) = pconv( jlon,jlat, nlev-k+1 ) - tder
      if( pconv( jlon,jlat, nlev-k+1 ) < 0 ) pconv( jlon,jlat, nlev-k+1 ) = 0.0
#endif

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
      timec    = max (1500.,timec)
      timec    = max (dts,timec)
      timec    = min (1800.,timec)
      if (ishall.eq.1) timec = 4300. ! care: redef. near the end - optimized with globo, aug. 2020

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
      print*, '!! Moisture budget error in convection_kf !!'
!      stop 'convection_kf'
      goto 999
      endif

!  feedback to resolvable scale tendencies

      if (ishall.eq.1) timec = 4300. ! optimized with globo, aug. 2020

      do k = 1, ltop
      nk = nlev-k+1
      dtdt  (jlon,jlat,nk) = (tg(k)-t0(k))/timec
      dqdt  (jlon,jlat,nk) = (qg(k)-q0(k))/timec
      dqcwdt(jlon,jlat,nk) = max(qlg(k)/timec,0.)
      dqcidt(jlon,jlat,nk) = max(qig(k)/timec,0.)
      dqrndt(jlon,jlat,nk) = max(qrg(k)/timec,0.)
      dqsndt(jlon,jlat,nk) = max(qsg(k)/timec,0.)
        if(t0(k).gt.t00+1..and.dqcidt(jlon,jlat,nk).gt.1.e-8) then
        dqcwdt(jlon,jlat,nk) = dqcwdt(jlon,jlat,nk) + dqcidt(jlon,jlat,nk)
        dqcidt(jlon,jlat,nk) = 0.
        endif
        if(t0(k).gt.t00+2..and.dqsndt(jlon,jlat,nk).gt.1.e-8) then
        dqrndt(jlon,jlat,nk) = dqrndt(jlon,jlat,nk) + dqsndt(jlon,jlat,nk)
        dqsndt(jlon,jlat,nk) = 0.
        endif
      enddo

      precon = max(trppt*ainc/dxsq*dts*float(1-ishall),0.)   ! precipitation rate (mm/sec) x timestep

#ifdef chem
      ! Renormalize pconv and econv using ainc factor; units mm/s
      ! Excludes the shallow cases
      pconv(jlon,jlat,:) = max( pconv(jlon,jlat,:)*ainc / dxsq * (1-ishall), 0.0 )
      econv(jlon,jlat,:) = max( econv(jlon,jlat,:)*ainc / dxsq * (1-ishall), 0.0 )
      ! FIXME: "max" is useless because pconv is already set as to be positive or zero here
#endif

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

! Computes thermodynamic properties of water, used by convection_kf

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
!##################################################################################################################
      subroutine orogdrag (jstep)

! Computes orographic wave drag

      use mod_model, only: nlon, nlat, nlev, nlevp1, nlonm1, nlatm1, phig, ps, u, v, hxt, phih, &
                           pzer, sigint, sigalf, dsig, dsigalf, rd, g, dtstep, dlon, dlat,      &
                           ip_e, ip_n, ip_s, ip_w, bvf, rich, tvirt, orogvar
      implicit none
      real, dimension(nlon,nlat,nlevp1) :: fmx, fmy
      real, dimension(nlon,nlat)        :: zorxx, zoryy, ordragx, ordragy
      real, dimension(nlev)             :: zwe
      real, dimension(nlevp1)           :: fup, fdw
      integer, dimension(nlon,nlat)     :: iormask, lcrit
      integer jstep, jlon, jlat, jklev, lc
      real z, zuav, zvav, zdphi, zn, zro, zsprod, zrich, zut, zvt, hra, znfa, zinx, ziny, zdp

! zdragc is a general amplitude coefficient;
! cext and zrichden determine the extinction coefficient zwe
! (cext is the inverse of a typical atmospheric depth of absorption).

      real, parameter :: zdragc = 0.028, cext = 1.7e-3, zrichden = 1.25
      save zorxx, zoryy, iormask

      fmx = 0.
      fmy = 0.
      ordragx = 0.
      ordragy = 0.

!----------------------------------------
!  1- Initialization (first instant only)
!----------------------------------------

      if (jstep.eq.1) then
      iormask = 0

! Second derivatives of orography
! Only crests, not valleys, are considered - a const. is subtracted to avoid "noise" over sea or hills
! A function of subgrid orography variance is added
! Threshold values for drag application (using iormask)
! znfa: normaliz. factor to keep the average drag in x nearly independent of model resolution
! 0.40+0.28 is dlon+dlat for globo 31 km res. - expon. 2./3. found empirically

      znfa = ((0.40+0.28)/(dlon+dlat))**(2./3.)

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

!-------------------------------------------------------------------
!  2- Computation of momentum flux at the surface (ordragx, ordragy)
!-------------------------------------------------------------------

! All computations on T points - tendencies averaged on u,v points

! Comput. of u e v as averages over the three lowest levels
! Brunt-Vaisala freq. (squared) bvf is computed in vdiff over the 3 lowest layers)
! zuav and zvav are averaged on T points

#ifdef mpi
        call u_ghost (u(nlonm1,:,:), ip_e, u(1,:,:), ip_w, nlat*nlev)
        call u_ghost (v(:,2,:), ip_s, v(:,nlat,:), ip_n, nlon*nlev)
#endif

      do jlat = 2, nlatm1
      do jlon = 2, nlonm1
      if (iormask(jlon,jlat).eq.1) then

      zuav = 1./6.*(u(jlon-1,jlat,nlev  )+u(jlon,jlat,nlev  )+   &
                    u(jlon-1,jlat,nlev-1)+u(jlon,jlat,nlev-1)+   &
                    u(jlon-1,jlat,nlev-2)+u(jlon,jlat,nlev-2))
      zvav = 1./6.*(v(jlon,jlat,nlev  )+v(jlon,jlat+1,nlev  )+   &
                    v(jlon,jlat,nlev-1)+v(jlon,jlat+1,nlev-1)+   &
                    v(jlon,jlat,nlev-2)+v(jlon,jlat+1,nlev-2))

      zn  = sqrt(max(bvf(jlon,jlat), 0.)) ! bvf is the square of Brunt-Vaisala frequency near surface
      if (zn.gt.1.5e-2) zn = 1.5e-2 ! to limit g.w. drag for high static stability
      zro = (pzer*sigint(nlev)-(pzer-ps(jlon,jlat))*sigalf(nlev))/(rd*tvirt(jlon,jlat,nlev)) ! density at nlev
      ordragx(jlon,jlat) = zdragc*zro*zn*zuav*zorxx(jlon,jlat)
      ordragy(jlon,jlat) = zdragc*zro*zn*zvav*zoryy(jlon,jlat)

! Computation of the critical level, defined as the level where the wind speed becomes zero
! or the wind vector becomes perpendicular to the orographic drag (that can have a different
! direction with respect to the surface (3 level-mean) wind vector).
! The critical level is defined where the scalar product between the wind and the drag
! becomes (for the first time going upward) null or positive.
! lcrit contains the index of the half-level below the level where the above condition is verified.
! In case of no critical level, lcrit = 2.

      lcrit(jlon,jlat) = 2
      do jklev = nlev-3, 1, -1
      zut = .5*(u(jlon-1,jlat,jklev) + u(jlon,jlat,  jklev))
      zvt = .5*(v(jlon,  jlat,jklev) + v(jlon,jlat+1,jklev))
      zsprod = zut*ordragx(jlon,jlat) + zvt*ordragy(jlon,jlat)
        if (zsprod.ge.0..or.rich(jlon,jlat,jklev).lt.0.) then
        lcrit(jlon,jlat) = max (lcrit(jlon,jlat), jklev+1)
        exit
        endif
      enddo

      endif ! condition on iormask
      enddo
      enddo

!----------------------------------------------------------------------------------------------------
!  3- Extintion coefficient (zwe) for upward and downward mom. flux
!----------------------------------------------------------------------------------------------------

! Gravity wave drag is applied up to the critical level
! The moist Richardson number computed in vdiff (defined on half levels) is used

      do jlat = 2, nlatm1
      do jlon = 2, nlonm1
      if (iormask(jlon,jlat).eq.1) then

! Computation of zwe on levels, depending on Richardson

      lc = lcrit(jlon,jlat)
      do jklev = nlev, lc-1, -1
      zrich = .5*(rich(jlon,jlat,jklev) + rich(jlon,jlat,jklev+1))
      zrich = min (zrich, 8.)
      zrich = max (zrich, -.5)
      zdphi = phih(jlon,jlat,jklev)-phih(jlon,jlat,jklev+1)
      zwe(jklev) = cext*zdphi/(g*(2.*zrich + zrichden))
      zwe(jklev) = min (zwe(jklev), 1.) ! zwe = 2 implies total absorption in one layer
      enddo
      if (rich(jlon,jlat,lcrit(jlon,jlat)).lt.0.) zwe(lcrit(jlon,jlat)+1) = 1.

!----------------------------------------------------------------------------
!  4- Computation of upward and downward momentum flux profile on half-levels
!----------------------------------------------------------------------------

      fup = 0.
      fdw = 0.
      fup(nlevp1) = 1.
      do jklev = nlev, lc, -1
      fup(jklev) = fup(jklev+1)*(1.-.5*zwe(jklev))/(1.+.5*zwe(jklev))
      enddo
      fdw(lc) = fup(lc)    ! reflection at critical level
      do jklev = lc+1, nlevp1
      fdw(jklev) = fdw(jklev-1)*(1.-.5*zwe(jklev-1))/(1.+.5*zwe(jklev-1))
      enddo

      do jklev = lc, nlevp1   ! x and y total momentum flux profile (on halh-levels)
      fmx(jlon,jlat,jklev) = ordragx(jlon,jlat)*(fup(jklev)-fdw(jklev))
      fmy(jlon,jlat,jklev) = ordragy(jlon,jlat)*(fup(jklev)-fdw(jklev))
      enddo

      endif ! cond. on iormask
      enddo
      enddo

!-----------------------------------------------
!  5- Divergence of momentum flux and u,v update
!-----------------------------------------------

#ifdef mpi
        call u_ghost (fmx(2,:,:), ip_w, fmx(nlon,:,:), ip_e, nlat*nlev)
        call u_ghost (fmy(:,nlatm1,:), ip_n, fmy(:,1,:), ip_s, nlon*nlev)
#endif

      do jlat = 2, nlatm1
      do jlon = 2, nlonm1
      if (iormask(jlon,jlat).eq.1) then
      lc = lcrit(jlon,jlat)
      do jklev = lc, nlev
      zdp = pzer*dsig(jklev)-(pzer-ps(jlon,jlat))*dsigalf(jklev)
      zinx = .5*dtstep*g/zdp*((fmx(jlon,jlat,jklev+1)+fmx(jlon+1,jlat,jklev+1))- &
                             (fmx(jlon,jlat,jklev  )+fmx(jlon+1,jlat,jklev  )))
      ziny = .5*dtstep*g/zdp*((fmy(jlon,jlat-1,jklev+1)+fmy(jlon,jlat,jklev+1))- &
                             (fmy(jlon,jlat-1,jklev  )+fmy(jlon,jlat,jklev  )))
      u(jlon,jlat,jklev) = u(jlon,jlat,jklev) + zinx
      v(jlon,jlat,jklev) = v(jlon,jlat,jklev) + ziny
      enddo
      endif ! condition on iormask

      enddo
      enddo

      return
      end subroutine orogdrag
!##################################################################################################################
      subroutine blockdrag (jstep)

! Computes drag due to orographic flow blocking

      use mod_model, only: nlon, nlat, nlev, nlonm1, nlatm1, dlon, dlat, phig, ps, u, v, hxt, &
                           phi, rd, g, sigint, sigalf, dtstep, bvf, tvirt, orogvar, pzer,     &
                           ip_e, ip_n, ip_s, ip_w
      implicit none
      real, dimension(nlon,nlat)    :: zorx, zory, zostdx, zostdy
      integer, dimension(nlon,nlat) :: iblmaskx, iblmasky
      integer jstep, jlon, jlat, jklev
      real zro, decscu, decscv, zilev, zuinc, zvinc, znfa, zpsx, zpsy, ztx, zty
      real znx, zny, bldragx, bldragy, zbvfx, zbvfy, zuav, zvav

! zbldragc is a general amplitude coefficient

      real, parameter :: zbldragc = 2.e-4

      save zorx, zory, zostdx, zostdy, iblmaskx, iblmasky

!----------------------------------
!  Initialization (first call only)
!----------------------------------

      if (jstep.eq.1) then
      iblmaskx = 0
      iblmasky = 0

! znfa: normaliz. factor
! 0.40+0.28 is dlon+dlat for globo 31 km res. - expon. 0.5 found empirically.
! With this factor zorx slightly decreases as model resolution increases

! Computing orography increments, orog. std at in u,v points and masks

      do jlat = 2, nlatm1
      znfa = ((0.40+0.28)/(dlon*hxt(jlat)/.707+dlat))**.5
      do jlon = 2, nlonm1
      zostdx(jlon,jlat) = .5*(orogvar(jlon,jlat)+orogvar(jlon+1,jlat))
      zorx(jlon,jlat) = znfa*(phig(jlon+1,jlat)-phig(jlon,jlat))/g
      if(abs(zorx(jlon,jlat)).gt.200.) iblmaskx(jlon,jlat) = 1
      if(zostdx(jlon,jlat).gt.250.) iblmaskx(jlon,jlat) = 1
      zostdy(jlon,jlat) = .5*(orogvar(jlon,jlat)+orogvar(jlon,jlat-1))
      zory(jlon,jlat) = znfa*(phig(jlon,jlat)-phig(jlon,jlat-1))/g
      if(abs(zory(jlon,jlat)).gt.200.) iblmasky(jlon,jlat) = 1
      if(zostdy(jlon,jlat).gt.250.) iblmasky(jlon,jlat) = 1
      enddo
      enddo

      endif  ! end of computations at the first instant only

#ifdef mpi
        call u_ghost (bvf(2,     :), ip_w, bvf(nlon,:), ip_e, nlat)
        call u_ghost (bvf(:,nlatm1), ip_n, bvf(:,1   ), ip_s, nlon)
#endif

!----------------------------------------------------------------------------------
!  Computation of momentum flux at the surface (bldragx, bldragy) on u and v points
!----------------------------------------------------------------------------------

! Comput. of near-surface wind (zuav, zvav) as averages over the three lowest levels
! Brunt-Vaisala freq. squared bvf is computed in vdiff over the lowest 3 layers

      do jlat = 2, nlatm1
      do jlon = 2, nlonm1

! U-component

      if (iblmaskx(jlon,jlat).eq.1) then
      zbvfx = .5*(bvf(jlon,jlat)+bvf(jlon+1,jlat)) ! bvf is the square of the Brunt-Vaisala freq.
      znx = sqrt(max(zbvfx, 0.))
        if(znx.gt.3.e-3.and.znx*(abs(zorx(jlon,jlat))+.35*zostdx(jlon,jlat)).gt.1.) then
        zuav = (u(jlon,jlat,nlev) + u(jlon,jlat,nlev-1) + u(jlon,jlat,nlev-2))/3.
        zuav = sign(sqrt(abs(zuav)), zuav) ! to reduce drag for high speed
        zpsx = .5*(ps(jlon,jlat) + ps(jlon+1,jlat))
        ztx = .5*(tvirt(jlon,jlat,nlev) + tvirt(jlon+1,jlat,nlev))
        zro = (pzer*sigint(nlev)-(pzer-zpsx)*sigalf(nlev))/rd/ztx ! density at nlev (u points)
        bldragx = sign(max(zbldragc*zro*znx*zuav*zorx(jlon,jlat), 0.), zuav) ! Orog. slope drag applied only at upslope flow
        bldragx = bldragx + .35*zbldragc*zro*znx*zuav*zostdx(jlon,jlat)      ! Orographic std. drag added

! Application of drag at lower levels, decaying with height

        decscu = 15.*zbvfx ! inverse of expon. decay scale
        do jklev = nlev, nlev/3, -1
        if (u(jlon,jlat,jklev)*zuav.le.0.) exit ! drag applied only until u and zuav have the same sign
        zilev = .5*(phi(jlon,jlat,jklev)+phi(jlon+1,jlat,jklev)-phig(jlon,jlat)-phig(jlon+1,jlat))/g
         if (zilev.le.zostdx(jlon,jlat)) then  ! slow linear decay at the lowest levels
         zuinc = dtstep*bldragx*(1.-0.3*zilev/(zostdx(jlon,jlat)))
         else                                   ! exponential decay above
         zuinc = dtstep*bldragx*0.7*exp(-decscu*(zilev-zostdx(jlon,jlat)))
         endif  ! cond. on zilev
        if(abs(zuinc).lt.1.e-2) exit
        zuinc = sign(min(0.7*abs(u(jlon,jlat,jklev)), abs(zuinc)), zuinc) ! limitation of increment
        u(jlon,jlat,jklev) = u(jlon,jlat,jklev) - zuinc
        enddo
        endif  ! condition on znx
      endif ! condition on iblmaskx

! V-component

      if (iblmasky(jlon,jlat).eq.1) then
      zbvfy = .5*(bvf(jlon,jlat-1)+bvf(jlon,jlat))
      zny = sqrt(max(zbvfy, 0.))
        if(zny.gt.3.e-3.and.zny*(abs(zory(jlon,jlat))+.35*zostdy(jlon,jlat)).gt.1.) then
        zvav = (v(jlon,jlat,nlev)+v(jlon,jlat,nlev-1)+v(jlon,jlat,nlev-2))/3.
        zvav = sign(sqrt(abs(zvav)), zvav) ! to reduce drag for high wind speed
        zpsy = .5*(ps(jlon,jlat-1) + ps(jlon,jlat))
        zty = .5*(tvirt(jlon,jlat-1,nlev) + tvirt(jlon,jlat,nlev))
        zro = (pzer*sigint(nlev)-(pzer-zpsy)*sigalf(nlev))/rd/zty ! density at nlev
        bldragy = sign(max(zbldragc*zro*zny*zvav*zory(jlon,jlat), 0.), zvav) ! Orog. slope drag applied only at upslope flow
        bldragy = bldragy + .35*zbldragc*zro*zny*zvav*zostdy(jlon,jlat)      ! Orographic std. drag added

! Application of drag at lower levels, decaying with height

        decscv = 15.*zbvfy ! inverse of expon. decay scale
        do jklev = nlev, nlev/3, -1
        if (v(jlon,jlat,jklev)*zvav.le.0.) exit ! drag applied only until v and zvav have the same sign
        zilev = .5*(phi(jlon,jlat,jklev)+phi(jlon,jlat-1,jklev)-phig(jlon,jlat)-phig(jlon,jlat-1))/g
         if (zilev.le.zostdy(jlon,jlat)) then  ! slow linear decay at the lowest levels
         zvinc = dtstep*bldragy*(1.-0.3*zilev/(zostdy(jlon,jlat)))
         else                                   ! exponential decay above
         zvinc = dtstep*bldragy*0.7*exp(-decscv*(zilev-zostdy(jlon,jlat)))
         endif  ! cond. on zilev
        if(abs(zvinc).lt.1.e-2) exit
        zvinc = sign(min(0.7*abs(v(jlon,jlat,jklev)), abs(zvinc)), zvinc) ! limitation of increment
        v(jlon,jlat,jklev) = v(jlon,jlat,jklev) - zvinc
        enddo
        endif ! condition on zny
      endif ! condition on iblmasky

      enddo
      enddo

      return
      end subroutine blockdrag
!##################################################################################################################
#ifdef mpi

      subroutine u_ghost (bufsend, ip_to, bufrecv, ip_from, nbuf)

      use mod_model, only: nbuffer
      implicit none
      include 'mpif.h'
      integer :: nbuf, ip_to, ip_from, ierr, comm, tag=1, status(mpi_status_size), isend, irecv
      real    :: bufsend(nbuf), bufrecv(nbuf)
      comm  = mpi_comm_world

      call mpi_isend (bufsend, nbuf, mpi_real, ip_to,   tag, comm, isend, ierr)
      call mpi_irecv (bufrecv, nbuf, mpi_real, ip_from, tag, comm, irecv, ierr)
      call mpi_wait  (isend, status, ierr)
      call mpi_wait  (irecv, status, ierr)
!      call mpi_send (bufsend, nbuf, mpi_real, ip_to  , tag, comm,         ierr)
!      call mpi_recv (bufrecv, nbuf, mpi_real, ip_from, tag, comm, status, ierr)

      return
      end subroutine u_ghost

#endif
!##################################################################################################################
    subroutine wrshf_b (jstep)

! For bolam, writes SHF files: files containing selected surface fields written at run-time

    use mod_model
    implicit none
    integer :: iunit=26, iunit_work=29, jlon, jlat, jstep, j, jklev, &
 jklev1, jklev2, nlevint, nlevout, iplev, jlonm1, jlatp1, &
 ilon1=1, ilon2=gnlon, jlat1=1, jlat2=gnlat, flag_lon=0
    integer, parameter :: nplev=3
    real zzpp, zzee, ztg, zt0t, zw, zesk, zqsat, ztbarv, zss, zalf, zalp1, zalp2, zalp3, zalp4, ze1, ze2
    real, dimension(nlevp1)    :: zlnp, zaux
    real, dimension(1)         :: zaux1
    real, dimension(100)       :: zaux2, zaux3
    real, dimension(nlon,nlat) :: ztz0, slp, q2rel, z1000, z850, z500, t850
    real, dimension(nplev) :: plev=(/1000.e2, 850.e2, 500.e2/), zalp
    real, dimension(nlon,nlat,nplev) :: zplev, tplev
    character(len=15)          :: file_out
    integer, dimension(50)     :: grib2_descript = 0
    integer, parameter         :: ivalmiss = -9999

    if (jstep.eq.1) qprec   = 0.
    if (jstep.eq.1) qprecc  = 0.
    if (jstep.eq.1) qsnfall = 0.

!  computation of slp (Pascal)

    zss = 6.5e-03     ! Average stability with climatological value of 6.5 K/km
    do jlat = 1, nlat
    do jlon = 1, nlon
    zzpp = pzer*sigint(nlev-1)-(pzer-ps(jlon,jlat))*sigalf(nlev-1)
    ztz0(jlon,jlat) = t(jlon,jlat,nlev-1) + zss*(phi(jlon,jlat,nlev-1)/g)
    ztg  = t(jlon,jlat,nlev-1) + zss*log(ps(jlon,jlat)/zzpp)*rd/g*t(jlon,jlat,nlev)*(1.+ep*q(jlon,jlat,nlev))
    ztbarv = .5*(ztz0(jlon,jlat)+ztg) * (1.+ep*q(jlon,jlat,nlev))
    slp(jlon,jlat) = ps(jlon,jlat)*exp(phig(jlon,jlat)/(ztbarv*rd))
    enddo
    enddo
    call filt2t1 (slp, 1.)
    call filt2t1 (slp, 1.)
    call filt2t1 (slp, 1.)
    call filt2t1 (slp, 1.)

!  completing external boundaries

!      if (ip_w.eq.ip_null) then
      do jlat = 1, nlat
      u10(1,jlat) = u10(2,jlat)
      v10(1,jlat) = v10(2,jlat)
      t2(1,jlat) = t2(2,jlat)
      q2(1,jlat) = q2(2,jlat)
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
      u10(nlon,jlat) = u10(nlonm1,jlat)
      v10(nlon,jlat) = v10(nlonm1,jlat)
      t2(nlon,jlat) = t2(nlonm1,jlat)
      q2(nlon,jlat) = q2(nlonm1,jlat)
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
      u10(jlon,1) = u10(jlon,2)
      v10(jlon,1) = v10(jlon,2)
      t2(jlon,1) = t2(jlon,2)
      q2(jlon,1) = q2(jlon,2)
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
      u10(jlon,nlat) = u10(jlon,nlatm1)
      v10(jlon,nlat) = v10(jlon,nlatm1)
      t2(jlon,nlat) = t2(jlon,nlatm1)
      q2(jlon,nlat) = q2(jlon,nlatm1)
      jklev=n_std_lev_sl(jlon,nlatm1)
      t_std_lev(jlon,nlat,1:jklev) = t_std_lev(jlon,nlatm1,1:jklev)
      u_std_lev(jlon,nlat,1:jklev) = u_std_lev(jlon,nlatm1,1:jklev)
      v_std_lev(jlon,nlat,1:jklev) = v_std_lev(jlon,nlatm1,1:jklev)
      q_std_lev(jlon,nlat,1:jklev) = q_std_lev(jlon,nlatm1,1:jklev)
      tg(jlon,nlat,1:nlevg) = tg(jlon,nlatm1,1:nlevg)
      qg(jlon,nlat,1:nlevg) = qg(jlon,nlatm1,1:nlevg)
      enddo
!      endif

! Definition of relative humidity at 2 m

      do jlat = 1, nlat
      do jlon = 1, nlon
      zt0t = tzer/t2(jlon,jlat)
      zesk = ezer*exp(-ccw1*log(zt0t)+ccw2*(1.-zt0t))  ! saturated pressure over water (Pa)
      zqsat = ps(jlon,jlat)*q2(jlon,jlat)/(eps+q2(jlon,jlat)*(1.-eps))
      q2rel(jlon,jlat) = min(zqsat/zesk*100., 100.)
      enddo
      enddo

! Lapse rate

      do jlat = 1, nlat
      do jlon = 1, nlon
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
        lapse_rate=(t(jlon,jlat,jklev2)-t(jlon,jlat,jklev1))*g/(phi(jlon,jlat,jklev2)-phi(jlon,jlat,jklev1))
      enddo
      enddo

!-----------------------------------------------------------
!  Interpolation of variables at constant pressure surfaces
!-----------------------------------------------------------

    zalf = .7

    do iplev = 1,nplev ! isobaric levels

      zalp(iplev) = log(plev(iplev))

      do jlat = 1, nlat
      do jlon = 1, nlon

        do jklev = 1, nlev
          zzpp = pzer*sigint(jklev)-(pzer-ps(jlon,jlat))*sigalf(jklev)
          zlnp(jklev) = log(zzpp)
        enddo
        zlnp(nlevp1) = log(slp(jlon,jlat)) ! auxiliary sea level

! In case zlnp(nlevp1) <= zlnp(nlev) (lowest model level below sea level) the auxiliary level
! at sea level pressure cannot be used - this concerns vert. interpolation of t and phi only

        if (zlnp(nlevp1).le.zlnp(nlev)) then
          nlevint = nlev
        else
          nlevint = nlevp1
        endif

        zaux1(1) = zalp(iplev)

!  Interpolation of geopotential (converted to geopotential meters)

        ze1 = 1.
        ze2 = 1.
        zaux(1:nlev) = phi(jlon,jlat,1:nlev)/9.8
        zaux(nlevp1) = 0.
        call interp (zalf, ze1, ze2, nlevint, zlnp, zaux, zaux1, zplev(jlon,jlat,iplev), 1)

!  Interpolation of t (Kelvin)

        ze1 = .4
        ze2 = .4
        zaux(1:nlev) = t(jlon,jlat,1:nlev)
        zaux(nlevp1) = ztz0(jlon,jlat)
        call interp (zalf, ze1, ze2, nlevint, zlnp, zaux, zaux1, tplev(jlon,jlat,iplev), 1)

      enddo
      enddo

      call filt2t1 (zplev(1:nlon,1:nlat,iplev), 1.)
      call filt2t1 (tplev(1:nlon,1:nlat,iplev), 1.)

    enddo ! isobaric levels

!-----------------------------------------------------------------------------
!  Interpolation of variables at constant geometric heigh (m) above the surface
!-----------------------------------------------------------------------------

    zalf = 0.7
    ze1 = 0.5
    ze2 = 0.5

    nlevint = nlev

    do jlat = 1, nlat
    do jlon = 1, nlon
      jlonm1=max(jlon-1,1)
      jlatp1=min(jlat+1,nlat)
      nlevout = n_std_lev_atm - n_std_lev_sl(jlon,jlat)

      if (nlevout >= 1 ) then

        do jklev = 1,nlev
          jklev1 = nlev-jklev+1
          zlnp(jklev1) = (phi(jlon,jlat,jklev)-phig(jlon,jlat))/g
        enddo

        do jklev = 1,nlevout
          jklev1 = jklev+n_std_lev_sl(jlon,jlat)
          zaux2(jklev) = std_lev_atm(jklev1)
        enddo

        do jklev = 1,nlev
          jklev1 = nlev-jklev+1
          zaux(jklev1) = t(jlon,jlat,jklev)
        enddo
        call interp (zalf, ze1, ze2, nlevint, zlnp, zaux, zaux2(1:nlevout), zaux3(1:nlevout), nlevout)
        do jklev = 1,nlevout
          jklev1 = jklev+n_std_lev_sl(jlon,jlat)
          t_std_lev(jlon,jlat,jklev1)=zaux3(jklev)
        enddo

        do jklev = 1,nlev
          jklev1 = nlev-jklev+1
          zaux(jklev1) = 0.5*(u(jlon,jlat,jklev)+u(jlonm1,jlat,jklev)) ! u computed on T-points
        enddo
        call interp (zalf, ze1, ze2, nlevint, zlnp, zaux, zaux2(1:nlevout), zaux3(1:nlevout), nlevout)
        do jklev = 1,nlevout
          jklev1 = jklev+n_std_lev_sl(jlon,jlat)
          u_std_lev(jlon,jlat,jklev1)=zaux3(jklev)
        enddo

        do jklev = 1,nlev
          jklev1 = nlev-jklev+1
          zaux(jklev1) = 0.5*(v(jlon,jlat,jklev)+v(jlon,jlatp1,jklev)) ! v computed on T-points
        enddo
        call interp (zalf, ze1, ze2, nlevint, zlnp, zaux, zaux2(1:nlevout), zaux3(1:nlevout), nlevout)
        do jklev = 1,nlevout
          jklev1 = jklev+n_std_lev_sl(jlon,jlat)
          v_std_lev(jlon,jlat,jklev1)=zaux3(jklev)
        enddo

        do jklev = 1,nlev
          jklev1 = nlev-jklev+1
          zaux(jklev1) = q(jlon,jlat,jklev)
        enddo
        call interp (zalf, ze1, ze2, nlevint, zlnp, zaux, zaux2(1:nlevout), zaux3(1:nlevout), nlevout)
        do jklev = 1,nlevout
          jklev1 = jklev+n_std_lev_sl(jlon,jlat)
          q_std_lev(jlon,jlat,jklev1)=zaux3(jklev)
        enddo

! Dew point

        do jklev = 1,nlev
          jklev1 = nlev-jklev+1
          zzpp = pzer*sigint(jklev) - (pzer-ps(jlon,jlat))*sigalf(jklev)
          zzee = zzpp*q(jlon,jlat,jklev)/(eps+q(jlon,jlat,jklev)*(1.-eps))
          if (t(jlon,jlat,jklev) >= tzer) then
            zaux(jklev1) = (tzer-33.65/17.40*log(zzee/611.))/(1.-1./17.40*log(zzee/611.))
          else
            zaux(jklev1) = (tzer- 0.75/22.45*log(zzee/611.))/(1.-1./22.45*log(zzee/611.))
          endif
        enddo
        call interp (zalf, ze1, ze2, nlevint, zlnp, zaux, zaux2(1:nlevout), zaux3(1:nlevout), nlevout)
        do jklev = 1,nlevout
          jklev1 = jklev+n_std_lev_sl(jlon,jlat)
          td_std_lev(jlon,jlat,jklev1)=zaux3(jklev)
        enddo

! Relative humidity

        do jklev = 1,nlev
          jklev1 = nlev-jklev+1
          zzpp = pzer*sigint(jklev) - (pzer-ps(jlon,jlat))*sigalf(jklev)
          zzee = zzpp*q(jlon,jlat,jklev)/(eps+q(jlon,jlat,jklev)*(1.-eps))
          zt0t = tzer/t(jlon,jlat,jklev)
          zesk = ezer*exp(-ccw1*log(zt0t)+ccw2*(1.-zt0t))  ! saturated pressure over water (Pa)
          zaux(jklev1) = zzee/zesk*100.
        enddo
        call interp (zalf, ze1, ze2, nlevint, zlnp, zaux, zaux2(1:nlevout), zaux3(1:nlevout), nlevout)
        do jklev = 1,nlevout
          jklev1 = jklev+n_std_lev_sl(jlon,jlat)
          rh_std_lev(jlon,jlat,jklev1)=zaux3(jklev)
        enddo

      endif

    enddo
    enddo

! Definition of relative humidity and dew point temperature in the surface layer

    do jlat = 1, nlat
    do jlon = 1, nlon
      do jklev = 1,n_std_lev_sl(jlon,jlat)
        if (t_std_lev(jlon,jlat,jklev) /= 0.) then
          zt0t = tzer/t_std_lev(jlon,jlat,jklev)
          zesk = ezer*exp(-ccw1*log(zt0t)+ccw2*(1.-zt0t))  ! saturated pressure over water (Pa)
          zqsat = ps(jlon,jlat)*q_std_lev(jlon,jlat,jklev)/(eps+q_std_lev(jlon,jlat,jklev)*(1.-eps))
          rh_std_lev(jlon,jlat,jklev) = min(zqsat/zesk*100., 100.)
          zzee = ps(jlon,jlat)*q_std_lev(jlon,jlat,jklev)/(eps+q_std_lev(jlon,jlat,jklev)*(1.-eps))
          if (t_std_lev(jlon,jlat,jklev) >= tzer) then
            td_std_lev(jlon,jlat,jklev) = (tzer-33.65/17.40*log(zzee/611.))/(1.-1./17.40*log(zzee/611.))
          else
            td_std_lev(jlon,jlat,jklev) = (tzer- 0.75/22.45*log(zzee/611.))/(1.-1./22.45*log(zzee/611.))
          endif
        else
            td_std_lev(jlon,jlat,jklev) = 0.
        endif
      enddo
    enddo
    enddo

    if (n_std_lev_soil >= 1 ) then

      nlevint = nlevg+1
      nlevout=n_std_lev_soil

      zlnp(1) = 0.
      do jklev = 1,nlevg
        zlnp(jklev+1) = lev_soil(jklev)
      enddo

      zaux2(1:nlevout) = std_lev_soil(1:nlevout)

      do jlat = 1, nlat
      do jlon = 1, nlon

        zaux(1) = tskin(jlon,jlat)
        do jklev = 1,nlevg
          zaux(jklev+1) = tg(jlon,jlat,jklev)
        enddo
        call interp (zalf, ze1, ze2, nlevint, zlnp, zaux, zaux2(1:nlevout), zaux3(1:nlevout), nlevout)
        tg_std_lev(jlon,jlat,1:nlevout)=zaux3(1:nlevout)

        zaux(1) = qg(jlon,jlat,1)
        do jklev = 1,nlevg
          zaux(jklev+1) = qg(jlon,jlat,jklev)
        enddo
        call interp (zalf, ze1, ze2, nlevint, zlnp, zaux, zaux2(1:nlevout), zaux3(1:nlevout), nlevout)
        qg_std_lev(jlon,jlat,1:nlevout)=zaux3(1:nlevout)

      enddo
      enddo

    endif

!   writing of slp, u10, v10, t2, q2rel, qprec....

    if(myid.eq.0) then

#ifdef oper
      ishf=ishf+1
      if (.not.nlclimate) then
        write (file_out,'(a,i3.3,a)') 'bolam_',ishf,'.shf'
      else
        write (file_out,'(a,i5.5,a)') 'bolam_',ishf,'.shf'
      endif
      open (iunit, file=file_out, form='unformatted', status='unknown')
#else
      open (iunit, file='bolam.shf', form='unformatted', status='unknown', position='append')
#endif

      write(iunit) nfdr
      write(iunit) pdr

    endif

! grib2_description - descriptor record for grib2 format coding

!  1 - model index (1 Bolam, 2 Moloch, 3 Globo)
!  2 - grid template index (1 horizontal grid, 1000 vertical cross-section)
!  3 - product template index (0 instant, 8 statistical, 32 forecast satellite, 1 instant individual ensemble, 11 statistical individual ensemble)
!  4 - time unit index (0 minute, 1 hour, 2 day, 3 month, 4 year, 13 second)
!  5 - statistical elaboration type (for statistical products only) : 0 average, 1 accumulation, 2 maximum, 3 minimum
!  6 - production status of data (0 Operational products, 2 Research products, 3 Re-analysis products, 7 Sub-seasonal to seasonal prediction S2S)
!  7 - type of data (0 Analysis, 1 Forecast, 2 Analysis and forecast, 3 Control forecast, 4 Perturbed forecast)
!  8 - indicator of time unit for the increment between successive fields used for statistical elaboration
! 10 - level (layer) type (for forecast satellite products code of satellite platform and sensor)
! 11 - first scaled value of level (layer)
! 12 - scale of first value of level (layer)
! 13 - second scaled value of level (layer)
! 14 - scale of second value of level (layer)
! 15 - level type of second level for a layer
! 20 - product discipline
! 21 - product category
! 22 - product parameter
! 30 - flag of bit-map presence: 0 not present, 1 present (uses VALMISS)
! 31 - member number of perturbed forecast
! 32 - member index of perturbed forecast

    grib2_descript( 1)    = 1 ! Bolam model
    grib2_descript( 2)    = 1 ! horizontal grid
    grib2_descript( 3)    = 0 ! instant product
    grib2_descript( 4)    = 1 ! time unit is hour
    grib2_descript( 6)    = 0 ! operational products
    grib2_descript( 7)    = 1 ! forecast
    grib2_descript(10:15) = ivalmiss
    grib2_descript(30)    = 0 ! bit-map absent

    call collect (slp, gfield)
    if(myid.eq.0) then
      grib2_descript(10) = 101 ! Mean sea level
      grib2_descript(20) =   0 ! Discipline: Meteorological products
      grib2_descript(21) =   3 ! Category: Mass
      grib2_descript(22) =   1 ! Parameter: Pressure reduced to MSL (Pa)
      write (iunit) grib2_descript
      call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
    endif

!    call collect (u10, gfield)
!    if(myid.eq.0) then
!      grib2_descript(10) = 103 ! Specified height level above ground (m)
!      grib2_descript(11) =  10 ! 10 m
!      grib2_descript(12) =   0
!      grib2_descript(20) =   0 ! Discipline: Meteorological products
!      grib2_descript(21) =   2 ! Category: Momentum
!      grib2_descript(22) =   2 ! Parameter: u-component of wind (m s-1)
!      write (iunit) grib2_descript
!      call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
!    endif
!
!    call collect (v10, gfield)
!    if(myid.eq.0) then
!      grib2_descript(10) = 103 ! Specified height level above ground (m)
!      grib2_descript(11) =  10 ! 10 m
!      grib2_descript(12) =   0
!      grib2_descript(20) =   0 ! Discipline: Meteorological products
!      grib2_descript(21) =   2 ! Category: Momentum
!      grib2_descript(22) =   3 ! Parameter: v-component of wind (m s-1)
!      write (iunit) grib2_descript
!      call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
!    endif
!
!    call collect (t2, gfield)
!    if(myid.eq.0) then
!      grib2_descript(10) = 103 ! Specified height level above ground (m)
!      grib2_descript(11) =   2 ! 2 m
!      grib2_descript(12) =   0
!      grib2_descript(20) =   0 ! Discipline: Meteorological products
!      grib2_descript(21) =   0 ! Category: Temperature
!      grib2_descript(22) =   0 ! Parameter: Temperature (K)
!      write (iunit) grib2_descript
!      call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
!    endif
!
!    call collect (q2rel, gfield)
!    if(myid.eq.0) then
!      grib2_descript(10) = 103 ! Specified height level above ground (m)
!      grib2_descript(11) =   2 ! 2 m
!      grib2_descript(12) =   0
!      grib2_descript(20) =   0 ! Discipline: Meteorological products
!      grib2_descript(21) =   1 ! Category: Moisture
!      grib2_descript(22) =   1 ! Parameter: Relative humidity (%)
!      write (iunit) grib2_descript
!      call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
!    endif

    grib2_descript(10) =   1 ! Ground or water surface
    grib2_descript(11:12) = 0

    call collect (qprec, gfield)
    if(myid.eq.0) then
      grib2_descript(3)  =   8 ! statistical product
      grib2_descript(5)  =   1 ! accumulation
      grib2_descript(20) =   0 ! Discipline: Meteorological products
      grib2_descript(21) =   1 ! Category: Moisture
      grib2_descript(22) =   8 ! Parameter: Total precipitation (kg m-2)
      write (iunit) grib2_descript
      call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
      grib2_descript(3) = 0 ! instant product
      grib2_descript(5) = 0
    endif

    call collect (cloudt, gfield)
    if(myid.eq.0) then
      grib2_descript(20) =   0 ! Discipline: Meteorological products
      grib2_descript(21) =   6 ! Category: Cloud
      grib2_descript(22) =   1 ! Parameter: Total cloud cover (%)
      write (iunit) grib2_descript
      call wrec2 (iunit, gnlon, gnlat, gfield*1.e2, ilon1, ilon2, jlat1, jlat2, flag_lon)
    endif

    if (jstep.gt.1) solar = solar/float(ntswshf)
    call collect (solar, gfield)
    if(myid.eq.0) then
      grib2_descript(3)  =   8 ! statistical product
      grib2_descript(5)  =   0 ! average
      grib2_descript(20) =   0 ! Discipline: Meteorological products
      grib2_descript(21) =   4 ! Category: Short-wave Radiation
      grib2_descript(22) =   7 ! Parameter: Downward short-wave radiation flux (W m-2)
      write (iunit) grib2_descript
      call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
    endif

    call collect (qsnfall, gfield)
    if(myid.eq.0) then
      grib2_descript(3)  =   8 ! statistical product
      grib2_descript(5)  =   1 ! accumulation
      grib2_descript(20) =   0 ! Discipline: Meteorological products
      grib2_descript(21) =   1 ! Category: Moisture
      grib2_descript(22) =  29 ! Parameter: Total snowfall (m)
      write (iunit) grib2_descript
      call wrec2 (iunit, gnlon, gnlat, gfield*1.e-3, ilon1, ilon2, jlat1, jlat2, flag_lon)
      grib2_descript(3) = 0 ! instant product
      grib2_descript(5) = 0
    endif

    call collect (snow, gfield)
    if(myid.eq.0) then
      grib2_descript(20) =   0 ! Discipline: Meteorological products
      grib2_descript(21) =   1 ! Category: Moisture
      grib2_descript(22) =  13 ! Parameter: Water equivalent of accumulated snow depth (kg m-2)
      write (iunit) grib2_descript
      call wrec2 (iunit, gnlon, gnlat, gfield*1.e3, ilon1, ilon2, jlat1, jlat2, flag_lon)
    endif

    call collect (lapse_rate, gfield)
    if(myid.eq.0) then
      grib2_descript(20) =   0 ! Discipline: Meteorological products
      grib2_descript(21) =   0 ! Category: Temperature
      grib2_descript(22) =   8 ! Parameter: Lapse rate (K m-1)
      write (iunit) grib2_descript
      call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
    endif

    call collect (tskin, gfield)
    if(myid.eq.0) then
      grib2_descript(20) =   0 ! Discipline: Meteorological products
      grib2_descript(21) =   0 ! Category: Temperature
      grib2_descript(22) =   0 ! Parameter: Temperature (K)
      write (iunit) grib2_descript
      call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
    endif

    if (jstep.gt.1) then
      call collect (shf_accum/float(ntswshf), gfield)
    else
      call collect (shf_accum, gfield)
    endif
    if(myid.eq.0) then
      grib2_descript(3)  =   8 ! statistical product
      grib2_descript(5)  =   0 ! average
      grib2_descript(20) =   0 ! Discipline: Meteorological products
      grib2_descript(21) =   0 ! Category: Temperature
      grib2_descript(22) =  11 ! Parameter: Sensible heat net flux (W m-2)
      write (iunit) grib2_descript
      call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
      grib2_descript(3)  =   0 ! instant product
      grib2_descript(5)  =   0
    endif

    if (jstep.gt.1) then
      call collect (lhf_accum/float(ntswshf), gfield)
    else
      call collect (lhf_accum, gfield)
    endif
    if(myid.eq.0) then
      grib2_descript(3)  =   8 ! statistical product
      grib2_descript(5)  =   0 ! average
      grib2_descript(20) =   0 ! Discipline: Meteorological products
      grib2_descript(21) =   0 ! Category: Temperature
      grib2_descript(22) =  10 ! Parameter: Latent heat net flux (W m-2)
      write (iunit) grib2_descript
      call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
      grib2_descript(3)  =   0 ! instant product
      grib2_descript(5)  =   0
    endif

    if (jstep.gt.1) then
      call collect (qf_accum, gfield)
    else
      call collect (qf_accum, gfield)
    endif
    if(myid.eq.0) then
      grib2_descript(3)  =   8 ! statistical product
      grib2_descript(5)  =   1 ! accumulation
      grib2_descript(20) =   0 ! Discipline: Meteorological products
      grib2_descript(21) =   1 ! Category: Moisture
      grib2_descript(22) =   6 ! Parameter: Evaporation (kg m-2)
      write (iunit) grib2_descript
      call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
      grib2_descript(3)  =   0 ! instant product
      grib2_descript(5)  =   0
    endif

!    grib2_descript(10) = 100 ! Isobaric surface  (Pa)
!    grib2_descript(11:12) = 0
!
!    iplev = 2
!    call collect (tplev(1:nlon,1:nlat,2), gfield)
!    if(myid.eq.0) then
!      grib2_descript(11) = nint(plev(iplev)*1.e-2) ! 850 hPa
!      grib2_descript(12) =  -2
!      grib2_descript(20) =   0 ! Discipline: Meteorological products
!      grib2_descript(21) =   0 ! Category: Temperature
!      grib2_descript(22) =   0 ! Parameter: Temperature (K)
!      write (iunit) grib2_descript
!      call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
!    endif

    iplev = 1
    call collect (zplev(1:nlon,1:nlat,iplev), gfield)
    if(myid.eq.0) then
      grib2_descript(11) = nint(plev(iplev)*1.e-2) ! 1000 hPa
      grib2_descript(12) =  -2
      grib2_descript(20) =   0 ! Discipline: Meteorological products
      grib2_descript(21) =   3 ! Category: Mass
      grib2_descript(22) =   5 ! Parameter: Geopotential height (gpm)
      write (iunit) grib2_descript
      call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
    endif

!
!    iplev = 2
!    call collect (zplev(1:nlon,1:nlat,iplev), gfield)
!    if(myid.eq.0) then
!      grib2_descript(11) = nint(plev(iplev)*1.e-2) ! 850 hPa
!      grib2_descript(12) =  -2
!      grib2_descript(20) =   0 ! Discipline: Meteorological products
!      grib2_descript(21) =   3 ! Category: Mass
!      grib2_descript(22) =   5 ! Parameter: Geopotential height (gpm)
!      write (iunit) grib2_descript
!      call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
!    endif
!
!    iplev = 3
!    call collect (zplev(1:nlon,1:nlat,iplev), gfield)
!    if(myid.eq.0) then
!      grib2_descript(11) = nint(plev(iplev)*1.e-2) ! 500 hPa
!      grib2_descript(12) =  -2
!      grib2_descript(20) =   0 ! Discipline: Meteorological products
!      grib2_descript(21) =   3 ! Category: Mass
!      grib2_descript(22) =   5 ! Parameter: Geopotential height (gpm)
!      write (iunit) grib2_descript
!      call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
!    endif

    grib2_descript(10) = 103 ! Specified height level above ground (m)
    grib2_descript(11:12) = 0

    do jklev = 1,n_std_lev_atm

      if (int(std_lev_atm(jklev)) == 2) then
        call collect (t_std_lev(1,1,jklev), gfield)
        if(myid.eq.0) then
          grib2_descript(11) =   int(std_lev_atm(jklev))
          grib2_descript(20) =   0 ! Discipline: Meteorological products
          grib2_descript(21) =   0 ! Category: Temperature
          grib2_descript(22) =   0 ! Parameter: Temperature (K)
          write (iunit) grib2_descript
          call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
        endif
      endif

      if (int(std_lev_atm(jklev)) >= 10) then
        call collect (u_std_lev(1,1,jklev), gfield)
        if(myid.eq.0) then
          grib2_descript(11) =   int(std_lev_atm(jklev))
          grib2_descript(20) =   0 ! Discipline: Meteorological products
          grib2_descript(21) =   2 ! Category: Momentum
          grib2_descript(22) =   2 ! Parameter: u-component of wind (m s-1)
          write (iunit) grib2_descript
          call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
        endif
      endif

      if (int(std_lev_atm(jklev)) >= 10) then
        call collect (v_std_lev(1,1,jklev), gfield)
        if(myid.eq.0) then
          grib2_descript(11) =   int(std_lev_atm(jklev))
          grib2_descript(20) =   0 ! Discipline: Meteorological products
          grib2_descript(21) =   2 ! Category: Momentum
          grib2_descript(22) =   3 ! Parameter: v-component of wind (m s-1)
          write (iunit) grib2_descript
          call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
        endif
      endif

      if (int(std_lev_atm(jklev)) == 2) then
        call collect (rh_std_lev(1,1,jklev), gfield)
        if(myid.eq.0) then
          grib2_descript(11) =   int(std_lev_atm(jklev))
          grib2_descript(20) =   0 ! Discipline: Meteorological products
          grib2_descript(21) =   1 ! Category: Moisture
          grib2_descript(22) =   1 ! Parameter: Relative humidity (%)
          write (iunit) grib2_descript
          call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
        endif
      endif

      if (int(std_lev_atm(jklev)) == 2) then
        call collect (q_std_lev(1,1,jklev), gfield)
        if(myid.eq.0) then
          grib2_descript(11) =   int(std_lev_atm(jklev))
          grib2_descript(20) =   0 ! Discipline: Meteorological products
          grib2_descript(21) =   1 ! Category: Moisture
          grib2_descript(22) =   0 ! Parameter: Specific humidity (kg kg-1)
          write (iunit) grib2_descript
          call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
        endif
      endif

      if (int(std_lev_atm(jklev)) == 2) then
        call collect (td_std_lev(1,1,jklev), gfield)
        if(myid.eq.0) then
          grib2_descript(11) =   int(std_lev_atm(jklev))
          grib2_descript(20) =   0 ! Discipline: Meteorological products
          grib2_descript(21) =   0 ! Category: Temperature
          grib2_descript(22) =   6 ! Parameter: Dew point temperature (K)
          write (iunit) grib2_descript
          call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
        endif
      endif

    enddo

    grib2_descript(10) = 106 ! Depth below land surface (m)
    grib2_descript(11) =   0
    grib2_descript(12) =   3

    do jklev = 1,n_std_lev_soil

      call collect (tg_std_lev(1,1,jklev), gfield)
      if(myid.eq.0) then
        grib2_descript(11) =   int(std_lev_soil(jklev)*1.e3)
        grib2_descript(20) =   2 ! Discipline:  Land surface products
        grib2_descript(21) =   0 ! Category: Vegetation/Biomass
        grib2_descript(22) =   2 ! Parameter: Soil temperature (K)
        write (iunit) grib2_descript
        call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
      endif

      call collect (qg_std_lev(1,1,jklev), gfield)
      if(myid.eq.0) then
        grib2_descript(11) =   int(std_lev_soil(jklev)*1.e3)
        grib2_descript(20) =   2 ! Discipline:  Land surface products
        grib2_descript(21) =   0 ! Category: Vegetation/Biomass
        grib2_descript(22) =   9 ! Parameter: Volumetric soil moisture content (Proportion)
        write (iunit) grib2_descript
        call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
      endif

    enddo

    if(myid.eq.0) then

      flush (iunit)
      close (iunit)

#ifdef oper
      call system("sync")
!!!      call system("ls -l -L "//file_out)
!!!      call system("date")
      open  (iunit_work, file=file_out(1:13)//'.txt', status='unknown')
      write (iunit_work,'(2a)') file_out,' is full and close'
      close (iunit_work)
      call system("sync")
      print*, file_out, ' written'
#else
      print*,'bolam.shf written'
#endif

    endif

    qprec     = 0. ! reset precipitation
    qprecc    = 0.
    qsnfall   = 0.
    solar     = 0.
    shf_accum = 0.
    lhf_accum = 0.
    qf_accum  = 0.

    return
    end subroutine wrshf_b
!##################################################################################################################
    subroutine wrshf_g (jstep)

! For globo, writes SHF files: files containing selected surface fields written at run-time

    use mod_model
    implicit none
    integer :: iunit=26, iunit_work=29, jlon, jlat, jstep, j, jklev, nlevint, &
 ilon1=1, ilon2=gnlon, jlat1=1, jlat2=gnlat, flag_lon=0
    real zzpp, ztg, zt0t, zw, zesk, zqsat, ztbarv, zss, zalf, zalp1, zalp2, zalp3, zalp4, ze1, ze2
    real, dimension(nlevp1)    :: zlnp, zaux
    real, dimension(nlon,nlat) :: ztz0, slp, q2rel, z850, z500, z50, t850, t500, u250, v250, u850, v850
    character(len=15)          :: file_out
    integer, dimension(50)     :: grib2_descript = 0
    integer, parameter         :: ivalmiss = -9999

!  Computation of slp (Pascal)

    zss = 6.0e-3     ! average stability in value k/km
    do jlat = 1, nlat
    do jlon = 1, nlon
    zzpp = pzer*sigint(nlev-1)-(pzer-ps(jlon,jlat))*sigalf(nlev-1)
    ztz0(jlon,jlat) = t(jlon,jlat,nlev-1) + zss*(phi(jlon,jlat,nlev-1)/g)
    ztg  = t(jlon,jlat,nlev-1) + zss*log(ps(jlon,jlat)/zzpp)*rd/g*t(jlon,jlat,nlev)*(1.+ep*q(jlon,jlat,nlev))
    ztbarv = .5*(ztz0(jlon,jlat)+ztg) * (1.+ep*q(jlon,jlat,nlev))
    slp(jlon,jlat) = ps(jlon,jlat)*exp(phig(jlon,jlat)/(ztbarv*rd))
    enddo
    enddo
    call filt2t1 (slp, 1.)
    call filt2t1 (slp, 1.)
    call filt2t1 (slp, 1.)
    call filt2t1 (slp, 1.)

!  Computation of relative humidity at 2m (in per cent)

    do jlat = 1, nlat
    do jlon = 1, nlon
    zt0t = tzer/t2(jlon,jlat)
    zesk = ezer*exp(-ccw1*log(zt0t)+ccw2*(1.-zt0t))  !! partial pressure over water
    zqsat = zesk*eps/(ps(jlon,jlat)+zesk*(eps-1.))
    q2rel(jlon,jlat) = q2(jlon,jlat)/zqsat*100.
    enddo
    enddo

!-----------------------------------------------------------
!  Interpolation of variables on constant pressure surfaces
!-----------------------------------------------------------

    zalf = .7
    zalp1 = log(500.e2)
    zalp2 = log(850.e2)
    zalp3 = log(250.e2)
    zalp4 = log(50.e2)

    do jlat = 1, nlat
    do jlon = 1, nlon

    do jklev = 1, nlev
    zzpp = pzer*sigint(jklev)-(pzer-ps(jlon,jlat))*sigalf(jklev)
    zlnp(jklev) = log(zzpp)
    enddo
    zlnp(nlevp1) = log(slp(jlon,jlat)) !! livello ausiliario sea level

! In case zlnp(nlevp1) <= zlnp(nlev) (lowest model level below sea level) the auxiliary level
! at sea level pressure cannot be used - this concerns vert. interpolation of t and phi only

    if (zlnp(nlevp1).le.zlnp(nlev)) then
    nlevint = nlev
    else
    nlevint = nlevp1
    endif

!  Interpolation of geopotential (converted to geopotential meters)

    ze1 = 1.
    ze2 = 1.
    zaux(1:nlev) = phi(jlon,jlat,1:nlev)/9.8
    zaux(nlevp1) = 0.
    call interp (zalf, ze1, ze2, nlevint, zlnp, zaux, zalp1, z500(jlon,jlat), 1)
    call interp (zalf, ze1, ze2, nlevint, zlnp, zaux, zalp4, z50 (jlon,jlat), 1)

!  Interpolation of t (converted to Celsius)

    ze1 = .4
    ze2 = .4
!    zaux(1:nlev) = t(jlon,jlat,1:nlev) - 273.15
!    zaux(nlevp1) = ztz0(jlon,jlat) - 273.15
    zaux(1:nlev) = t(jlon,jlat,1:nlev)
    zaux(nlevp1) = ztz0(jlon,jlat)
    call interp (zalf, ze1, ze2, nlevint, zlnp, zaux, zalp2, t850(jlon,jlat), 1)
    call interp (zalf, ze1, ze2, nlevint, zlnp, zaux, zalp1, t500(jlon,jlat), 1)

!  Interpolation of u, v

    ze1 = .4
    ze2 = .4
    zaux(1:nlev) = u(jlon,jlat,1:nlev)
    zaux(nlevp1) = 0.
    call interp (zalf, ze1, ze2, nlevint, zlnp, zaux, zalp2, u850(jlon,jlat), 1)
    zaux(1:nlev) = v(jlon,jlat,1:nlev)
    zaux(nlevp1) = 0.
    call interp (zalf, ze1, ze2, nlevint, zlnp, zaux, zalp2, v850(jlon,jlat), 1)

    zaux(1:nlev) = u(jlon,jlat,1:nlev)
    zaux(nlevp1) = 0.
    call interp (zalf, ze1, ze2, nlevint, zlnp, zaux, zalp3, u250(jlon,jlat), 1)
    zaux(1:nlev) = v(jlon,jlat,1:nlev)
    zaux(nlevp1) = 0.
    call interp (zalf, ze1, ze2, nlevint, zlnp, zaux, zalp3, v250(jlon,jlat), 1)

    enddo
    enddo

    call filt2t1 (t850, 1.)
    call filt2t1 (t500, 1.)
    call filt2t1 (z500, 1.)
    call filt2t1 (z50 , 1.)

!-----------------------------------------------------------
!  Writting
!-----------------------------------------------------------

! grib2_description - descriptor record for grib2 format coding

!  1 - model index (1 Bolam, 2 Moloch, 3 Globo)
!  2 - grid template index (1 horizontal grid, 1000 vertical cross-section)
!  3 - product template index (0 instant, 8 statistical, 32 forecast satellite, 1 instant individual ensemble, 11 statistical individual ensemble)
!  4 - time unit index (0 minute, 1 hour, 2 day, 3 month, 4 year, 13 second)
!  5 - statistical elaboration type (for statistical products only) : 0 average, 1 accumulation, 2 maximum, 3 minimum
!  6 - production status of data (0 Operational products, 2 Research products, 3 Re-analysis products, 7 Sub-seasonal to seasonal prediction S2S)
!  7 - type of data (0 Analysis, 1 Forecast, 2 Analysis and forecast, 3 Control forecast, 4 Perturbed forecast)
!  8 - indicator of time unit for the increment between successive fields used for statistical elaboration
! 10 - level (layer) type (for forecast satellite products code of satellite platform and sensor)
! 11 - first scaled value of level (layer)
! 12 - scale of first value of level (layer)
! 13 - second scaled value of level (layer)
! 14 - scale of second value of level (layer)
! 15 - level type of second level for a layer
! 20 - product discipline
! 21 - product category
! 22 - product parameter
! 30 - flag of bit-map presence: 0 not present, 1 present (uses VALMISS)
! 31 - member number of perturbed forecast
! 32 - member index of perturbed forecast

    grib2_descript( 1)    = 3 ! Globo model
    grib2_descript( 2)    = 1 ! horizontal grid
    grib2_descript( 3)    = 0 ! instant product
    grib2_descript( 4)    = 1 ! time unit is hour
    grib2_descript( 6)    = 0 ! operational products
    grib2_descript( 7)    = 1 ! forecast
    grib2_descript(10:15) = ivalmiss
    grib2_descript(30)    = 0 ! bit-map absent

#ifndef oper
    if (myid.eq.0) then
      open (iunit, file='globo.shf', form='unformatted', status='unknown', position='append')

      write (iunit) nfdr
      write (iunit) pdr

    endif
#endif

    if (jstep.eq.1) then

#ifdef oper
      if (myid.eq.0) then
        file_out="globo_const.shf"
        open (iunit, file=trim(file_out), form='unformatted', status='unknown')

        write (iunit) nfdr
        write (iunit) pdr

      endif
#endif

      call collect (phig, gfield)
      if (myid.eq.0) then
        grib2_descript(10) =   1 ! Ground or water surface
        grib2_descript(20) =   2 ! Discipline: Land surface products
        grib2_descript(21) =   0 ! Category: Vegetation
        grib2_descript(22) =   7 ! Parameter: Model terrain height (m)
        write (iunit) grib2_descript
        gfield(:,:) = gfield(:,:)/g
        call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
      endif

      call collect (fmask, gfield)
      if (myid.eq.0) then
        grib2_descript(10) =   1 ! Ground or water surface
        grib2_descript(20) =   2 ! Discipline: Land surface products
        grib2_descript(21) =   0 ! Category: Vegetation
        grib2_descript(22) =   0 ! Parameter: Land cover (0=land, 1=sea) (Proportion)
        write (iunit) grib2_descript
        gfield(:,:) = 1.-gfield(:,:)
        call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
      endif

#ifdef oper
      if(myid.eq.0) then
        close (iunit)
        open  (iunit_work, file=trim(file_out)//'.txt', status='unknown')
        write (iunit_work,'(2a)') trim(file_out),' is full and close'
        close (iunit_work)
        print*, trim(file_out), ' written'
      endif
#endif

    endif ! jstep.eq.1

#ifdef oper
    if (myid.eq.0) then
      ishf=ishf+1
      if (.not.nlclimate) then
        write (file_out,'(a,i3.3,a)') 'globo_',ishf,'.shf'
      else
        write (file_out,'(a,i5.5,a)') 'globo_',ishf,'.shf'
      endif
      open (iunit, file=trim(file_out), form='unformatted', status='unknown')

      write (iunit) nfdr
      write (iunit) pdr

    endif
#endif

    call collect (slp, gfield)
    if (myid.eq.0) then
      grib2_descript(10) =   1 ! Ground or water surface
!      grib2_descript(10) = 101 ! Mean sea level
      grib2_descript(20) =   0 ! Discipline: Meteorological products
      grib2_descript(21) =   3 ! Category: Mass
      grib2_descript(22) =   1 ! Parameter: Pressure reduced to MSL (Pa)
      write (iunit) grib2_descript
      call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
    endif

    call collect (tg(1,1,1), gfield)
    if (myid.eq.0) then
      grib2_descript(10) = 106 ! Depth below land surface (m)
      grib2_descript(11) = lev_soil(1)*1.e3
      grib2_descript(12) =   3
      grib2_descript(20) =   2 ! Discipline: Land surface products
      grib2_descript(21) =   0 ! Category: Vegetation/Biomass
      grib2_descript(22) =   2 ! Parameter: Soil temperature (K)
      write (iunit) grib2_descript
      call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
    endif

    grib2_descript(10) =   1 ! Ground or water surface
    grib2_descript(11:12) = 0

    grib2_descript(3)  =   8 ! statistical product
    grib2_descript(5)  =   1 ! accumulation

    call collect (qprec, gfield)
    if (myid.eq.0) then
      grib2_descript(20) =   0 ! Discipline: Meteorological products
      grib2_descript(21) =   1 ! Category: Moisture
      grib2_descript(22) =   8 ! Parameter: Total precipitation (kg m-2)
      write (iunit) grib2_descript
      call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
    endif

    call collect (qprecc, gfield)
    if (myid.eq.0) then
      grib2_descript(20) =   0 ! Discipline: Meteorological products
      grib2_descript(21) =   1 ! Category: Moisture
      grib2_descript(22) =  10 ! Parameter: Convective precipitation (kg m-2)
      write (iunit) grib2_descript
      call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
    endif

    call collect (qsnfall, gfield)
    if (myid.eq.0) then
      grib2_descript(20) =   0 ! Discipline: Meteorological products
      grib2_descript(21) =   1 ! Category: Moisture
      grib2_descript(22) =  29 ! Parameter: Total snowfall (m)
      write (iunit) grib2_descript
      gfield(:,:) = gfield(:,:)*1.e-3
      call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
    endif

    grib2_descript(3) = 0 ! instant product
    grib2_descript(5) = 0

    call collect (snow, gfield)
    if (myid.eq.0) then
      grib2_descript(20) =   0 ! Discipline: Meteorological products
      grib2_descript(21) =   1 ! Category: Moisture
      grib2_descript(22) =  13 ! Parameter: Water equivalent of accumulated snow depth (kg m-2)
      write (iunit) grib2_descript
      gfield(:,:) = gfield(:,:)*1.e3
      call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
    endif

    grib2_descript(10) = 103 ! Specified height level above ground (m)

    call collect (u10, gfield)
    if (myid.eq.0) then
      grib2_descript(11) =  10 ! 10 m
      grib2_descript(12) =   0
      grib2_descript(20) =   0 ! Discipline: Meteorological products
      grib2_descript(21) =   2 ! Category: Momentum
      grib2_descript(22) =   2 ! Parameter: u-component of wind (m s-1)
      write (iunit) grib2_descript
      call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
    endif

    call collect (v10, gfield)
    if (myid.eq.0) then
      grib2_descript(11) =  10 ! 10 m
      grib2_descript(12) =   0
      grib2_descript(20) =   0 ! Discipline: Meteorological products
      grib2_descript(21) =   2 ! Category: Momentum
      grib2_descript(22) =   3 ! Parameter: v-component of wind (m s-1)
      write (iunit) grib2_descript
      call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
    endif

    call collect (t2, gfield)
    if (myid.eq.0) then
      grib2_descript(11) =   2 ! 2 m
      grib2_descript(12) =   0
      grib2_descript(20) =   0 ! Discipline: Meteorological products
      grib2_descript(21) =   0 ! Category: Temperature
      grib2_descript(22) =   0 ! Parameter: Temperature (K)
      write (iunit) grib2_descript
      call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
    endif

    call collect (q2rel, gfield)
    if (myid.eq.0) then
      grib2_descript(11) =   2 ! 2 m
      grib2_descript(12) =   0
      grib2_descript(20) =   0 ! Discipline: Meteorological products
      grib2_descript(21) =   1 ! Category: Moisture
      grib2_descript(22) =   1 ! Parameter: Relative humidity (%)
      write (iunit) grib2_descript
      call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
    endif

    grib2_descript(10) =   1 ! Ground or water surface
    grib2_descript(11:12) = 0

    call ccloud
    call collect (cloudt, gfield)
    if (myid.eq.0) then
      grib2_descript(20) =   0 ! Discipline: Meteorological products
      grib2_descript(21) =   6 ! Category: Cloud
      grib2_descript(22) =   1 ! Parameter: Total cloud cover (%)
      write (iunit) grib2_descript
      gfield(:,:) = gfield(:,:)*1.e2
      call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
    endif

    call collect (fice, gfield)
    if (myid.eq.0) then
      grib2_descript(20) =   10! Discipline: Oceanographic products
      grib2_descript(21) =   2 ! Category: Ice
      grib2_descript(22) =   0 ! Parameter: Ice cover (Proportion)
      write (iunit) grib2_descript
      call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
    endif

    grib2_descript(10) = 100 ! Isobaric surface  (Pa)

    call collect (t850, gfield)
    if (myid.eq.0) then
      grib2_descript(11) = 85000 ! 850 hPa
      grib2_descript(12) =   0
      grib2_descript(20) =   0 ! Discipline: Meteorological products
      grib2_descript(21) =   0 ! Category: Temperature
      grib2_descript(22) =   0 ! Parameter: Temperature
      write (iunit) grib2_descript
      call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
    endif

    call collect (z500, gfield)
    if (myid.eq.0) then
      grib2_descript(11) = 50000 ! 500 hPa
      grib2_descript(12) =   0
      grib2_descript(20) =   0 ! Discipline: Meteorological products
      grib2_descript(21) =   3 ! Category: Mass
      grib2_descript(22) =   5 ! Parameter: Geopotential height (gpm)
      write (iunit) grib2_descript
      call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
    endif

    call collect (z50 , gfield)
    if (myid.eq.0) then
      grib2_descript(11) =  5000 ! 50 hPa
      grib2_descript(12) =   0
      grib2_descript(20) =   0 ! Discipline: Meteorological products
      grib2_descript(21) =   3 ! Category: Mass
      grib2_descript(22) =   5 ! Parameter: Geopotential height (gpm)
      write (iunit) grib2_descript
      call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
    endif

    grib2_descript(10) =   1 ! Ground or water surface
    grib2_descript(11:12) = 0

    if (jstep.gt.1) solar = solar*dtstep/float(ntswshf)
    call collect (solar, gfield)
    if (myid.eq.0) then
      grib2_descript(3)  =   8 ! statistical product
      grib2_descript(5)  =   1 ! accumulation
      grib2_descript(20) =   0 ! Discipline: Meteorological products
      grib2_descript(21) =   4 ! Category: Short-wave Radiation
      grib2_descript(22) =   7 ! Parameter: Downward short-wave radiation flux (W m-2)
      write (iunit) grib2_descript
      call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
      grib2_descript(3) = 0 ! instant product
      grib2_descript(5) = 0
    endif

    grib2_descript(10) = 100 ! Isobaric surface  (Pa)

    call collect (u250, gfield)
    if (myid.eq.0) then
      grib2_descript(11) = 25000 ! 250 hPa
      grib2_descript(12) =   0
      grib2_descript(20) =   0 ! Discipline: Meteorological products
      grib2_descript(21) =   2 ! Category: Momentum
      grib2_descript(22) =   2 ! Parameter: u-component of wind (m s-1)
      write (iunit) grib2_descript
      call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
    endif

    call collect (v250, gfield)
    if (myid.eq.0) then
      grib2_descript(11) = 25000 ! 250 hPa
      grib2_descript(12) =   0
      grib2_descript(20) =   0 ! Discipline: Meteorological products
      grib2_descript(21) =   2 ! Category: Momentum
      grib2_descript(22) =   3 ! Parameter: v-component of wind (m s-1)
      write (iunit) grib2_descript
      call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
    endif

    call collect (t500, gfield)
    if (myid.eq.0) then
      grib2_descript(11) = 50000 ! 500 hPa
      grib2_descript(12) =   0
      grib2_descript(20) =   0 ! Discipline: Meteorological products
      grib2_descript(21) =   0 ! Category: Temperature
      grib2_descript(22) =   0 ! Parameter: Temperature (K)
      write (iunit) grib2_descript
      call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
    endif

    call collect (u850, gfield)
    if (myid.eq.0) then
      grib2_descript(11) = 85000 ! 850 hPa
      grib2_descript(12) =   0
      grib2_descript(20) =   0 ! Discipline: Meteorological products
      grib2_descript(21) =   2 ! Category: Momentum
      grib2_descript(22) =   2 ! Parameter: u-component of wind (m s-1)
      write (iunit) grib2_descript
      call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
    endif

    call collect (v850, gfield)
    if (myid.eq.0) then
      grib2_descript(11) = 85000 ! 850 hPa
      grib2_descript(12) =   0
      grib2_descript(20) =   0 ! Discipline: Meteorological products
      grib2_descript(21) =   2 ! Category: Momentum
      grib2_descript(22) =   3 ! Parameter: v-component of wind (m s-1)
      write (iunit) grib2_descript
      call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
    endif

    if(myid.eq.0) then

      flush (iunit)
      close (iunit)

#ifdef oper
      call system("sync")
!!!      call system("ls -l -L "//file_out)
!!!      call system("date")
      open  (iunit_work, file=trim(file_out)//'.txt', status='unknown')
      write (iunit_work,'(2a)') trim(file_out),' is full and close'
      close (iunit_work)
      call system("sync")
      print*, trim(file_out), ' written'
#else
      print*,'globo.shf written'
#endif

    endif

! reset total precipitation and cumulated fluxes

    if (jstep.ne.1) then
      qprec  = 0.
      qprecc = 0.
      qsnfall= 0.
      solar  = 0.
      olrtot = 0.
    endif

    return
    end subroutine wrshf_g
!##################################################################################################################
    subroutine comp_esk (esat, qsat, t, p, iflag)

! Computes esat from temperature and qsat from absolute temperature and pressure.
! Iflag 1: esat and qsat with respect to water and ice, separately, depending if t>tzer or t<tzer
!       2: esat and qsat with an interpolation at t<tzer between water and ice
!       3: esat and qsat with respect to water also for t<tzer

    use mod_model, only : tzer, ezer, cpv, cw, rd, rv, yliv, yliw, eps, ci, ylwv, ccw1, ccw2, cci1, cci2
    implicit none
    real  esat, qsat, t, p, zt0t, zesk, zratio
    integer iflag

    zt0t = tzer/t
    if (zt0t.le.1.) then
     zesk = ezer*exp(-ccw1*alog(zt0t)+ccw2*(1.-zt0t)) ! partial pressure over water
     else
      if (iflag.eq.1) then
      zesk = ezer*exp(-cci1*alog(zt0t)+cci2*(1.-zt0t)) ! partial pressure over ice
      elseif (iflag.eq.2) then
      zratio = 1.04979*(0.5 + 0.5*tanh((t-tzer+9.)/6.))
      zesk = zratio*(ezer*exp(-ccw1*alog(zt0t)+ccw2*(1.-zt0t)))+      &
        (1.-zratio)*(ezer*exp(-cci1*alog(zt0t)+cci2*(1.-zt0t)))
      elseif (iflag.eq.3) then
      zesk = ezer*exp(-ccw1*alog(zt0t)+ccw2*(1.-zt0t))
      else
      print*, "iflag out of range in comp_esk", iflag
      stop
      endif
    endif

    esat=zesk
    qsat = zesk*eps/(p+zesk*(eps-1.))
    return
    end subroutine comp_esk
!##################################################################################################################
#ifdef rad_ecmwf
#ifndef oldrad
   subroutine aerdef (nyrc, nmonc, ndayc, nhouc, nminc, infx, infy)

! For new (2012) ECMWF rad.
! Aerosol and ozone definition - called at initial time and at long intervals

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

use mod_model, only : nlon, nlat, nlev, gnlon, gnlat, gfield, gfield1, &
                      dlon, dlat, sig, sigint, sigalf, ps, pzer, t,    &
                      tskin, aerosol, ozon, myid, alont, alatt, tzer,  &
                      g, phig, aerotot, rdecli0, reqtim0, fmask
implicit none

integer jlon, jlat, jklev, nyrc, nmonc, ndayc, nhouc, nminc, infx, infy
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

character (len=30) :: str

!  load external functions

#include "fctast.h"
#include "fcttim.h"
#include "fcttre.h"

if (myid.eq.0) print*, 'Definition of aerosol and ozone'

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
  rdecli0 = rdecli
  reqtim0 = reqtim
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

! zeta and zetah: vertical sigma coordin. for the routines that set parameters for aerosol

    do jklev = 1, nlev
    zeta (jklev) = sigint(jklev)
    zetah(jklev) = sig(jklev)
    enddo
    zetah(nlev+1) = 1.

! Latitudes and longitudes in radians

    zdegrad = rpi/180.

! Definition of cloud parameters (in module yoecld, used by radlswr, called in radintec)

    call sucld (nlev, zeta)

!-------------------------------------------------------------------
! Definition of atmospheric variables

    do jlat = 1, nlat

    do jlon = 1, nlon

#ifdef globo
    plat5(jlon)  = (-90. + (infy+jlat-2)*dlat) * zdegrad
    plon5(jlon)  = (       (infx+jlon-2)*dlon) * zdegrad
#else
    plat5(jlon)  = alatt(jlon,jlat) * zdegrad
    plon5(jlon)  = alont(jlon,jlat) * zdegrad  ! long. < 0 not accepted
#endif

    if (plon5(jlon).lt.0.) plon5(jlon) = plon5(jlon) + 2.*rpi
    pgelam5(jlon)= plon5(jlon)
    pgemu5(jlon) = sin(plat5(jlon))
    pclon5(jlon) = cos(plon5(jlon))
    pslon5(jlon) = sin(plon5(jlon))
    enddo

! Pressure at full levels

    do jklev = 1, nlev
    do jlon = 1, nlon
    pap5(jlon,jklev) = pzer*sigint(jklev) - (pzer-ps(jlon,jlat))*sigalf(jklev)
    enddo
    enddo

! Pressure at half levels

    do jlon = 1, nlon
    paph5(jlon,1     ) = 1.e-8
    paph5(jlon,nlev+1) = ps(jlon,jlat)
    enddo
    do jklev = 2, nlev
    do jlon = 1, nlon
    paph5(jlon,jklev) = 0.5*(pap5(jlon,jklev-1)+pap5(jlon,jklev))
    enddo
    enddo

! Temperature at full levels

    do jklev = 1, nlev
    do jlon = 1, nlon
    pt5(jlon,jklev) = t(jlon,jlat,jklev)
    enddo
    enddo

! Temperature at half-levels: at top defined as lev. 1
! at bottom as tskin, elsewhere is the arithmetic mean

    do jlon = 1, nlon
    pth5(jlon,1     ) = pt5(jlon,1)
    pth5(jlon,nlev+1) = tskin(jlon,jlat)
    pts5(jlon)        = tskin(jlon,jlat)
    enddo

    do jklev = 2, nlev
    do jlon = 1, nlon
    pth5(jlon,jklev) = 0.5*(pt5(jlon,jklev-1)+pt5(jlon,jklev))
    enddo
    enddo

! Call to other set-up routines for various coefficients. These are dependent on the vertical resolution

    call suaerv (nlev, zetah, cvdaes, cvdael, cvdaeu, cvdaed, rctrbga, rcvobga, rcstbga, rcaeops, rcaeopl, &
                 rcaeopu, rcaeopd, rctrpt, rcaeadk, rcaeadm, rcaeros)

! Derive the aerosols and ozone distribution from climatology

    call radaca ( 1 , nlon , nlon , 1 , nlev ,                                       &
                  paph5 , pgelam5 , pgemu5 , pclon5 , pslon5 , pth5 , paer5 , zozon5)

! The computation of ozone must be done separately at each point in longitude, because radocz
! assumes (incorrectly for rotated grid) that latitude is the same for all vectors in longitude

    do jlon= 1, nlon
    call radozc (1, 1, 1, 1, nlev, 1, 1, 0, paph5(jlon,:), pgemu5(jlon), zozon5(jlon,:))
    enddo
    pozon5 = zozon5

    if (naer.eq.0) paer5 = zepaer

    do jklev = 1, nlev
    do jlon = 1, nlon
    ozon   (jlon,jlat,jklev  ) = pozon5(jlon,  jklev)     ! ozone
    aerosol(jlon,jlat,jklev,1) = paer5 (jlon,1,jklev)     ! land (organic + sulfate) aerosol
    aerosol(jlon,jlat,jklev,2) = paer5 (jlon,2,jklev)     ! sea salt aerosol
    aerosol(jlon,jlat,jklev,3) = paer5 (jlon,3,jklev)     ! desert dust aerosol
    aerosol(jlon,jlat,jklev,4) = paer5 (jlon,4,jklev)     ! urban + black carbon aerosol
    aerosol(jlon,jlat,jklev,5) = paer5 (jlon,5,jklev)     ! volcanic aerosol
    aerosol(jlon,jlat,jklev,6) = paer5 (jlon,6,jklev)     ! stratospheric background aerosol
    enddo
    enddo

#ifndef globo

! Ad hoc increase of some aerosol (urban and sulfate) over the Po Valley in the lower troposphere

    do jklev = 1, nlev
    do jlon = 1, nlon
    if(alatt(jlon,jlat).gt.44.3.and.alatt(jlon,jlat).lt.46.2.and.   &
       alont(jlon,jlat).gt. 7.0.and.alont(jlon,jlat).lt.13.4.and.   &
       phig(jlon,jlat)/g.lt.500..and.pap5(jlon,jklev).gt.800.e2) then
    aerosol(jlon,jlat,jklev,1) = 1.10*aerosol(jlon,jlat,jklev,1)
    aerosol(jlon,jlat,jklev,4) = 1.20*aerosol(jlon,jlat,jklev,4)
    endif
    enddo
    enddo
#endif

    enddo ! jlat

    do jf = 1, 2
    call filt2t (ozon, 1.)
    enddo
    ozon = max (ozon, 0.)
    do jf = 1, 5
    do jk = 1, 6
    call filt2t (aerosol(1:nlon,1:nlat,1:nlev,jk), 1.)
    enddo
    enddo
    do jf = 1, 60
    call filt2t (aerosol(1:nlon,1:nlat,1:nlev, 2), 1.)
    enddo
    aerosol = max (aerosol, zepaer)

#ifdef globo
    call polavert (ozon)
    ozon = max (ozon, 0.)
    do jk = 1, 6
    call polavert (aerosol(1:nlon,1:nlat,1:nlev,jk))
    aerosol(:,:,:,jk) = max (aerosol(:,:,:,jk), sngl(zepaer))
    enddo
#endif

    aerotot(:,:,:) = aerosol(:,:,:,1) + aerosol(:,:,:,2) + aerosol(:,:,:,3) +  &
                     aerosol(:,:,:,4) + aerosol(:,:,:,5) + aerosol(:,:,:,6)

    goto 345  ! skip plots for verification

!    call collect (fmask, gfield1) ! examples for 2-D plotting
!    if (myid.eq.0) call plotout (gfield, gfield1, gnlon, gnlat, str, 99)
!    str = 'albedo'
!    call collect (albedo, gfield)
!    if (myid.eq.0) call plotout (gfield, gfield1, gnlon, gnlat, str, 99)

    call collect (fmask, gfield1)
    str = 'ozon(3)'
    call collect (ozon(1:nlon,1:nlat,3), gfield)
    if (myid.eq.0) call plotout (gfield, gfield1, gnlon, gnlat, str, 99)
    str = 'aerosol1(nlev-5)'
    call collect (aerosol(1:nlon,1:nlat,nlev-5  ,1), gfield)
    if (myid.eq.0) call plotout (gfield, gfield1, gnlon, gnlat, str, 99)
    str = 'aerosol2(nlev-5)'
    call collect (aerosol(1:nlon,1:nlat,nlev-5  ,2), gfield)
    if (myid.eq.0) call plotout (gfield, gfield1, gnlon, gnlat, str, 99)
    str = 'aerosol3(nlev-5)'
    call collect (aerosol(1:nlon,1:nlat,nlev-5  ,3), gfield)
    if (myid.eq.0) call plotout (gfield, gfield1, gnlon, gnlat, str, 99)
    str = 'aerosol4(nlev-5)'
    call collect (aerosol(1:nlon,1:nlat,nlev-5  ,4), gfield)
    if (myid.eq.0) call plotout (gfield, gfield1, gnlon, gnlat, str, 99)
    str = 'aerosol5(nlev-5)'
    call collect (aerosol(1:nlon,1:nlat,nlev-5  ,5), gfield)
    if (myid.eq.0) call plotout (gfield, gfield1, gnlon, gnlat, str, 99)
    str = 'aerosol6(nlev-5)'
    call collect (aerosol(1:nlon,1:nlat,nlev-5  ,6), gfield)
    if (myid.eq.0) call plotout (gfield, gfield1, gnlon, gnlat, str, 99)
    str = 'aerotot(nlev-5)'
    call collect (aerotot(1:nlon,1:nlat,nlev-5), gfield)
    if (myid.eq.0) call plotout (gfield, gfield1, gnlon, gnlat, str, 99)
    str = 'aerosol1(nlev)'
    call collect (aerosol(1:nlon,1:nlat,nlev  ,1), gfield)
    if (myid.eq.0) call plotout (gfield, gfield1, gnlon, gnlat, str, 99)
    str = 'aerosol2(nlev)'
    call collect (aerosol(1:nlon,1:nlat,nlev  ,2), gfield)
    if (myid.eq.0) call plotout (gfield, gfield1, gnlon, gnlat, str, 99)
    str = 'aerosol3(nlev)'
    call collect (aerosol(1:nlon,1:nlat,nlev  ,3), gfield)
    if (myid.eq.0) call plotout (gfield, gfield1, gnlon, gnlat, str, 99)
    str = 'aerosol4(nlev)'
    call collect (aerosol(1:nlon,1:nlat,nlev  ,4), gfield)
    if (myid.eq.0) call plotout (gfield, gfield1, gnlon, gnlat, str, 99)
    str = 'aerotot(nlev)'
    call collect (aerotot(1:nlon,1:nlat,nlev), gfield)
    if (myid.eq.0) call plotout (gfield, gfield1, gnlon, gnlat, str, 99)

#ifdef mpi
      call mpi_finalize (ierr)
#endif
    stop

345 continue

    return
    end subroutine aerdef
#endif
#endif
!##################################################################################################################
#ifdef rad_ecmwf
#ifdef oldrad
      subroutine aerdef (nyrc, nmonc, ndayc, nhouc, nminc, infx, infy)

! For old (2004) ECMWF rad.
! aerosol and ozone definition
! INTEGER_M and REAL_B must be uppercase

#include "tsmbkind.h"

  use yomcst,    only : rd, rg, rtt, rsigma, rcpd, rpi, rday, rea, repsm, ri0
  use yoeaerd,   only : cvdaes, cvdael, cvdaeu, cvdaed, rcaeops, rcaeopl, rcaeopu, rcaeopd, rctrbga,          &
                        rcvobga, rcstbga, rctrpt, rcaeadm, rcaeros, rcaeadk
  use yoerad,    only : naer, nmode, nozocl, nradfr, nradpfr, nradpla, nrint, nhowinh, novlp, nradf2c, nradc2f, &
                        nlw, nsw, ntsw, lerad6h, leradhs, lhvolca, lonewsw, lowasyf, lowhsss, loifuec, lrrtm,   &
                        lradlp, linhom, raovlp, rbovlp, niceopt, nliqopt, nradip, nradlp, rminice, lnewaer
  use yoerdi,    only : rcardi, rch4, rn2o, ro3, rcfc11, rcfc12, repclc
  use yoerdu,    only : nuaer, ntraer, nout, rcday, r10e, replog, repsc, repsco, repscq, repsct, repscw, diff
  use yoesw,     only : lo3only

  use mod_model, only : nlon, nlat, nlev, gnlon, gnlat, gfield, nlevp1, dlon, dlat, sig, sigint, sigalf,  &
                        ps, pzer, t, tskin, aerosol, ozon, aerotot, myid, rdecli0, reqtim0, alont, alatt

  parameter (jp_lon=nlon, jp_lev=nlev, jp_lw=6, jp_sw=6, jp_nua=24, jp_mode=1, jp_aer=6, jp_levp1=jp_lev+1)

!- reference arrays

  REAL_B :: paer5(jp_lon,jp_aer,jp_lev)
  REAL_B :: paph5(jp_lon,jp_levp1), pap5(jp_lon,jp_lev)
  REAL_B :: pth5(jp_lon,jp_levp1) , pt5(jp_lon,jp_lev), pts5(jp_lon)
  REAL_B :: pgemu5(jp_lon), pslon5(jp_lon), pclon5(jp_lon), pspgem(jp_lon)
  REAL_B :: pozon5(jp_lon,jp_lev)
  REAL_B :: plat5(jp_lon), plon5(jp_lon), zozon5(jp_lon,jp_lev)

  REAL_B :: plon6(jp_lon), zeta(jp_lev), zetah(jp_levp1)

! dummy integer variables

  INTEGER_M :: kaer, kidia, kfdia, ktdia, klw, ksw, kmode
  INTEGER_M :: iflag, ilay, kulout, ndump, nflux
  INTEGER_M :: nrad, nindat, nsss, month, iday, iyr, ilwrad, ifin, iminut, inhomf
  INTEGER_M :: io3only

! dummy real variables

  REAL_B :: prii05, zepaer, zthetoz, zangozc, zdegrad, ztt5, zcondw, zfract, zpp5, zrr5, zddqq, zee, zqq5
  REAL_B :: rtimtr, zcoolc, zcoolt, zheatc, zheatt, zteta

! load external functions

#include "fctast.h"   ! def. funzioni astronomiche
#include "fcttim.h"   ! def. funzioni temporali
#include "fcttre.h"   ! def. funz. termodinamiche

!------------------------------------------------------------------
! DEFINE THE CONFIGURATION
!------------------------------------------------------------------

    kidia = 1
    kfdia = nlon
    ktdia = 1
    kmode = jp_mode
    kaer  = jp_aer
    klw   = jp_lw

! ECMWF climatological aerosols (tanre=1, gads=2, ECMWF oper)

    naer = 2
    if (naer.eq.1) lnewaer=.false.
    if (naer.eq.2) lnewaer=.true.        ! use of new aerosol climatology
    if (myid.eq.0) write(*,'(a, i2)') ' Aerosol definition with naer = ', naer

! ksw number of SW spectral intervals: 2, 4 or 6

    ksw     = jp_sw
    nflux   = 6
    nmode   = 0
    nrad    = 1
    nradfr  = -3
    nradpfr = 36
    nradpla = 15
    nrint   = 4
    nradf2c = 1
    nradc2f = 1
    nuaer   = 31
    ntraer  = 19
    novlp   = 1
    nlw     = 6
    ntsw    = 6
    nsw     = ksw

! nindat:  YYYYMMDD date of the simulation, for example 19870701

    iiyr = nyrc
    nindat =  iiyr*10000 + nmonc*100 + ndayc

! nsss: no. sec. del giorno (arrotondati al minuto)

    nsss = nhouc*3600 + nminc*60

! basic constants
! 0 come ultimo argomento: non fa stampe; 1: stampa tutte le cost.

    kulout = 6  ! unita' logica di output per la sucst
    call sucst (kulout, nindat, nsss, 0)     ! def. costanti tra cui rday necess. alla rtime

! rtimtr: time in sec. after 11:57 gmt of 1 jan 2000

    rtimtr = rtime (iiyr, nmonc, ndayc, nsss) ! function in /include/fcttim.h

! calcolo cost. astronomiche come ECMWF

    zteta  = rteta(rtimtr) ! zteta tempo in fraz. di anno dalle 11:57 del 1 genn.
    rdecli = rds(zteta)    ! decl. terra
    reqtim = ret(zteta)    ! time equation
    rdecli0 = rdecli
    reqtim0 = reqtim

! Load Ozone climatology
! NOZOCL =-1 tests the vertical quadrature
! NOZOCL =0 whatever is read for O3 as input profile
! NOZOCL =1 old ECMWF O3 and climatol. aerosols
! NOZOCL =2 Fortuin-Langematz O3 climatology + aero
! NOZOCL =3 old ECMWF O3 and no aerosol
! NOZOCL =4 Fortuin-Langematz O3 and no aerosol

    nozocl=2

! io3only =0 oper.; =1 only O3 absorption in UV-Vis.

    io3only = 0
    lo3only =.false.

    call surdi
    call sulwn
    call suswn    (ntsw,ksw)
    call suaerl
    call suaerh
    call suaersn  (ntsw,ksw)
    call suclopn  (ntsw,ksw,nlev)

    if (nozocl.eq.-1) then
    rch4  =1.e-18
    rn2o  =1.e-18
    ro3   =1.e-18
    rcfc11=1.e-18
    rcfc12=1.e-18
    endif

    zepaer = 1.e-12
    zepaers= 1.e-12

! Fortuin-Langematz O3 climatology
! IMINUT e' il no. di minuti dall'inizio del giorno (non IMIN!)

    iminut = nsss/60
    call suecozc (nindat, iminut)

! ECMWF Geleyn O3 climatology

    zthetoz = rteta(rtimtr)
    zangozc = rel(zthetoz)-1.7535
    call suecozo (zangozc)

! nuova climat. aerosol

    if (naer.eq.2) then
    call suecaebc                  ! black carbon aerosol (fires)
    call suecaeor                  ! organic type aerosol
    call suecaesd                  ! soil-dust aerosol
    call suecaess                  ! sea-salt aerosol
    call suecaesu                  ! sulfate-type aerosol
    call suecaec (nindat, iminut)  ! new climatol. distrib. of aerosol
    endif

! definiz. delle variabili atmosferiche

    do 100 jlat = 1, nlat

! latitudes and longitudes in radians

    zdegrad = rpi/180.

    do jlon = 1, nlon

#ifdef globo
    plat5(jlon)  = (-90. + (infy+jlat-2)*dlat) * zdegrad
    plon5(jlon)  = (       (infx+jlon-2)*dlon) * zdegrad
#else
    plat5(jlon)  = alatt(jlon,jlat) * zdegrad
    plon5(jlon)  = alont(jlon,jlat) * zdegrad  ! long. < 0 not accepted
#endif

    plon6(jlon)  = plon5(jlon) + rpi
    if (plon6(jlon).gt.2.*rpi) plon6(jlon) = plon6(jlon)-2.*rpi
    pgemu5(jlon) = sin(plat5(jlon))
    pspgem(jlon) = sqrt(1.-pgemu5(jlon)**2)
    pclon5(jlon) = cos(plon5(jlon))
    pslon5(jlon) = sin(plon5(jlon))
    enddo

! Pressione sui livelli

    do jklev = 1, nlev
    do jlon = 1, nlon
    pap5(jlon,jklev) = pzer*sigint(jklev) - (pzer-ps(jlon,jlat))*sigalf(jklev)
    enddo
    enddo

! Pressione sui semi-livelli

    do jlon = 1, nlon
    paph5(jlon,1     ) = 0.
    paph5(jlon,nlevp1) = ps(jlon,jlat)
    enddo
    do jklev = 2, nlev
    do jlon = 1, nlon
    paph5(jlon,jklev) = 0.5*(pap5(jlon,jklev-1)+pap5(jlon,jklev))
    enddo
    enddo

! Temperatura sui livelli

    do jklev = 1, nlev
    do jlon = 1, nlon
    pt5(jlon,jklev) = t(jlon,jlat,jklev)
    enddo
    enddo

! temperatura sui semi-livelli: al top si def. uguale a quella del liv. 1
! al bottom e' la Tskin, altrove si fa la media aritmetica

    do jlon = 1, nlon
    pth5(jlon,1     ) = pt5(jlon,1)
    pth5(jlon,nlevp1) = tskin(jlon,jlat)
    pts5(jlon)        = tskin(jlon,jlat)
    enddo

    do jklev = 2, nlev
    do jlon = 1, nlon
    pth5(jlon,jklev) = 0.5*(pt5(jlon,jklev-1)+pt5(jlon,jklev))
    enddo
    enddo

! ZETA e ZETAH sono coordin. vert. sigma per le due routines che
! seguono e che settano parametri per l'aerosol - sono espresse in SIGMA

    do jklev = 1, nlev
    zeta (jklev) = sigint(jklev)
    zetah(jklev) = sig(jklev)
    enddo
    zetah(nlevp1) = 1.

! call to other set-up routines for various coefficients. These are dependent on the vertical resolution

    call sucld (nlev, zeta)

    call suaerv (nlev, zetah, cvdaes, cvdael, cvdaeu, cvdaed, rctrbga, rcvobga, rcstbga, rcaeops, rcaeopl, &
                 rcaeopu, rcaeopd, rctrpt, rcaeadk, rcaeadm, rcaeros)

! derive the aerosols and ozone distribution from climatology

    if (naer.eq.1) then
    call radaca ( kidia , kfdia , kfdia  , ktdia , nlev , &
                  paph5 , plon5 , pgemu5 , pclon5 , pslon5 , pth5 , paer5 , zozon5)
    endif
    if (naer.eq.2) then
    call radaca ( kidia , kfdia , kfdia  , ktdia , nlev , &
                  paph5 , plon6 , -pgemu5 , -pclon5 , pslon5 , pth5 , paer5 , zozon5)
    endif

    if (nozocl.eq.1 .or. nozocl.eq.3) pozon5 = zozon5

    if (nozocl.eq.-1) then
    paer5  = zepaer
    pozon5 = 0.
    endif

    if (nozocl.eq.2 .or. nozocl.eq.4) then
    call radozc (kidia, kfdia, kfdia, ktdia, nlev, 1, kfdia, 0, paph5, pgemu5, zozon5)
    pozon5 = zozon5
    endif

    if (nozocl.gt.2) then
    paer5 = zepaer
    else
    do jklev = 1, nlev
    do jlon = 1, nlon
    paer5(jlon,5,jklev) = zepaer  ! volcanic aerosol set to epsilon in absence of eruption
    enddo
    enddo
    endif

    paer5 = max (zepaer, paer5)     ! security check on aerosol amounts

    do jklev = 1, nlev
    do jlon = 1, nlon
    ozon   (jlon,jlat,jklev  ) = pozon5(jlon,  jklev)
    aerosol(jlon,jlat,jklev,1) = paer5 (jlon,1,jklev)
    aerosol(jlon,jlat,jklev,2) = paer5 (jlon,2,jklev)
    aerosol(jlon,jlat,jklev,3) = paer5 (jlon,3,jklev)
    aerosol(jlon,jlat,jklev,4) = paer5 (jlon,4,jklev)
    aerosol(jlon,jlat,jklev,5) = paer5 (jlon,5,jklev)
    aerosol(jlon,jlat,jklev,6) = paer5 (jlon,6,jklev)
    enddo
    enddo

#ifndef globo

! Ad hoc increase of some aerosol (urban and sulfate) over the Po Valley in the lower troposphere

    do jklev = 1, nlev
    do jlon = 1, nlon
    if(alatt(jlon,jlat).gt.44.3.and.alatt(jlon,jlat).lt.46.2.and.   &
       alont(jlon,jlat).gt. 7.0.and.alont(jlon,jlat).lt.13.4.and.   &
       phig(jlon,jlat)/g.lt.500..and.pap5(jlon,jklev).gt.800.e2) then
    aerosol(jlon,jlat,jklev,1) =  1.10*aerosol(jlon,jlat,jklev,1)
    aerosol(jlon,jlat,jklev,4) =  1.20*aerosol(jlon,jlat,jklev,4)
    endif
    enddo
    enddo
#endif

100 continue

    do jf = 1, 2
    call filt2t (ozon, 1.)
    enddo
    ozon = max (ozon, 0.)
    do jf = 1, 5
    do jk = 1, 6
    call filt2t (aerosol(1:nlon,1:nlat,1:nlev,jk), 1.)
    enddo
    enddo
    do jf = 1, 60
    call filt2t (aerosol(1:nlon,1:nlat,1:nlev,2), 1.)
    enddo
    aerosol = max (aerosol, zepaer)

#ifdef globo
    call polavert (ozon)
    ozon = max (ozon, 0.)
    do jk = 1, 6
    call polavert (aerosol(1:nlon,1:nlat,1:nlev,jk))
    aerosol(:,:,:,jk) = max (aerosol(:,:,:,jk), sngl(zepaer))
    enddo
#endif

    aerotot(:,:,:) = aerosol(:,:,:,1) + aerosol(:,:,:,2) + aerosol(:,:,:,3) +  &
                     aerosol(:,:,:,4) + aerosol(:,:,:,5) + aerosol(:,:,:,6)

    return
    end subroutine aerdef
#endif
#endif
!##################################################################################################################
#ifdef rad_ecmwf
#ifndef oldrad
   subroutine radintec (jlat, jl1, nyrc, nmonc, ndayc, nhouc, nminc, infx, infy, rrf)

! For new (2012) ECMWF rad.
! Interface subroutine between Bolam and ECMWF radiation version 35R2
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
! pfrsod5: total-sky surface sw downward flux (not used)

! psudu5:  solar radiance in sun's direction (not used)
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

use mod_model, only : nlon, nlat, nlonr, nradm, nlev, dlon, dlat, cpd, cpv, sig, sigint, sigalf, ps, &
                      pzer, t, tskin, t2, aerosol, ozon, myid, corvis, corirr, corrdt, tzer, q, qcw, &
                      qci, qrain, qsnow, qhail, fcloud, ccw1, ccw2, cci1, cci2, eps, ezer, emisg1,   &
                      emisg2, albedo, fsnow, fmask, fice, alsn, g, co2ppm, mcica, alont, alatt,      &
                      cloudt, rdecli0, reqtim0

#ifdef gas
  use mod_chem, only: top_level
  use mod_gas, only: phk_attn_fact
#endif

implicit none

integer jlon, jlat, jklev, nyrc, nmonc, ndayc, nhouc, nminc, infx, infy
integer klon, klev, kaer, kaero, ksw, iiyr, ierr, jl, jl1
real zalsn, zfsnow, zdstcp, zomd, rrf

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
  rdecli0 = rdecli
  reqtim0 = reqtim
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

#ifdef globo
  plat5(jl)  = (-90. + (infy+jlat-2)*dlat) * zdegrad
  plon5(jl)  = (       (infx+jlon-2)*dlon) * zdegrad
#else
  plat5(jl)  = alatt(jlon,jlat) * zdegrad
  plon5(jl)  = alont(jlon,jlat) * zdegrad  ! long. < 0 not accepted
#endif

  if (plon5(jl).lt.0.) plon5(jl) = plon5(jl) + 2.*rpi
  pgelam5(jl)= plon5(jl)
  pgemu5(jl) = sin(plat5(jl))
  pclon5(jl) = cos(plon5(jl))
  pslon5(jl) = sin(plon5(jl))
  pmu05(jl)  = max(0._jprb, rsidec*pgemu5(jl)-rcodec*sqrt(1.-pgemu5(jl)**2)*cos(pgelam5(jl)+rwsovr))
  enddo

! Pressure at full levels

  do jklev = 1, nlev
  do jl = 1, nlonr
  jlon = jl1 + (jl-1)*nradm
  pap5(jl,jklev) = pzer*sigint(jklev) - (pzer-ps(jlon,jlat))*sigalf(jklev)
  enddo
  enddo

! Pressure at half levels

  do jl = 1, nlonr
  jlon = jl1 + (jl-1)*nradm
  paph5(jl,1     ) = 0.
  paph5(jl,nlev+1) = ps(jlon,jlat)
  enddo
  do jklev = 2, nlev
  do jl = 1, nlonr
  paph5(jl,jklev) = 0.5*(pap5(jl,jklev-1)+pap5(jl,jklev))
  enddo
  enddo

! Temperature at full levels

  do jklev = 1, nlev
  do jl = 1, nlonr
  jlon = jl1 + (jl-1)*nradm
  pt5(jl,jklev) = t(jlon,jlat,jklev)
  enddo
  enddo

! Temperature at half-levels: at top defined as lev. 1
! at bottom as t at 2m or tskin, elsewhere is the arithmetic mean

  do jl = 1, nlonr
  jlon = jl1 + (jl-1)*nradm
  pth5(jl,1     ) = pt5(jl,1)
!  pth5(jl,nlev+1) = t2(jlon,jlat)
  pth5(jl,nlev+1) = tskin(jlon,jlat)
  pts5(jl)        = tskin(jlon,jlat)
  enddo
  do jklev = 2, nlev
  do jl = 1, nlonr
  pth5(jl,jklev) = 0.5*(pt5(jl,jklev-1)+pt5(jl,jklev))
  enddo
  enddo

  do jklev = 1, nlev
  do jl = 1, nlonr
  jlon = jl1 + (jl-1)*nradm
  pozon5(jl,  jklev) = ozon   (jlon,jlat,jklev  )
  paer5 (jl,1,jklev) = aerosol(jlon,jlat,jklev,1)
  paer5 (jl,2,jklev) = aerosol(jlon,jlat,jklev,2)
  paer5 (jl,3,jklev) = aerosol(jlon,jlat,jklev,3)
  paer5 (jl,4,jklev) = aerosol(jlon,jlat,jklev,4)
  paer5 (jl,5,jklev) = aerosol(jlon,jlat,jklev,5)
  paer5 (jl,6,jklev) = aerosol(jlon,jlat,jklev,6)
  enddo
  enddo

!--------------------------------------------------------------------------------

  do jklev = 1, nlev
  do jl = 1, nlonr
  jlon = jl1 + (jl-1)*nradm
  pq5(jl,jklev)     = q(jlon,jlat,jklev)
  pclfr5(jl,jklev)  = fcloud(jlon,jlat,jklev)
  pqiwp5(jl,jklev)  = qci(jlon,jlat,jklev)+0.2*qsnow(jlon,jlat,jklev)
  pqlwp5(jl,jklev)  = qcw(jlon,jlat,jklev)
  pdp5  (jl,jklev)  = paph5(jl,jklev+1)-paph5(jl,jklev)

!  The LW and SW radiation schemes can handle profiles of the trace gases.
!  Here fixed concentrations are used.

  pco25(jl,jklev) = zco2
  pch45(jl,jklev) = zch4
  pn2o5(jl,jklev) = zn2o
  pno25(jl,jklev) = zno2
  pc115(jl,jklev) = zc11
  pc125(jl,jklev) = zc12
  pc225(jl,jklev) = zc22
  pcl45(jl,jklev) = zcl4
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
  zalsn = alsn(jlon,jlat)

  zalb  (jl) = albedo(jlon,jlat)*(1.-zfsnow) + zalsn*zfsnow
  pemis5(jl) = emisg1(jlon,jlat)*(1.-zfsnow) + 0.87*zfsnow
  pemiw5(jl) = emisg2(jlon,jlat)*(1.-zfsnow) + 0.87*zfsnow
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

  do jklev = 1, nlev
  do jl = 1, nlonr
  jlon = jl1 + (jl-1)*nradm
  zdstcp = g/(cpd*(1.+zomd*q(jlon,jlat,jklev))*(paph5(jl,jklev+1)-paph5(jl,jklev)))
  corrdt(jlon,jlat,jklev) = zdstcp*                                                 &
              (pflt5(jl,jklev)-pflt5(jl,jklev+1)+pfls5(jl,jklev)-pfls5(jl,jklev+1))
  enddo
  enddo

! Surface fluxes of visible and infrared radiation (positive downward)

  do jl = 1, nlonr
  jlon = jl1 + (jl-1)*nradm
  corvis(jlon,jlat) = pfls5(jl,nlev+1)
!  corirr(jlon,jlat) = pflt5(jl,nlev+1)
  corirr(jlon,jlat) = max(pflt5(jl,nlev+1), 0.92*(1.-0.05*cloudt(jlon,jlat))*pflt5(jl,nlev+1)) ! tuning!

  if (rrf.gt.1.e-2) then ! Initial tapering
  corvis(jlon,jlat) = corvis(jlon,jlat)*(1. - rrf*cloudt(jlon,jlat)**.5)
  corirr(jlon,jlat) = corirr(jlon,jlat)*(1. - rrf*cloudt(jlon,jlat)**.5)
  endif

#ifdef gas
  do jklev = top_level, nlev
    if( pfcs5(jl,jklev) /= 0 )then
      phk_attn_fact(jlon,jlat,jklev)=pfls5(jl,jklev)/pfcs5(jl,jklev)
    else
      phk_attn_fact(jlon,jlat,jklev)=0
    end if
  end do
#endif

  enddo

  return
  end subroutine radintec
#endif
#endif
!##################################################################################################################
#ifdef rad_ecmwf
#ifdef oldrad
  subroutine radintec (jlat, jl1, nyrc, nmonc, ndayc, nhouc, nminc, infx, infy)

! For old (2004) ECMWF rad.
! interface subroutine between GLOBO and Morcrette ECMWF radiation v. 2004 (Morcrette, mar. 2005)
! INTEGER_M e REAL_B devono rimanere maiuscoli

#include "tsmbkind.h"

  use yomcst,    only : rd, rg, rtt, rsigma, rcpd, rpi, rday, rea, repsm, ri0
  use yoeaerd,   only : cvdaes, cvdael, cvdaeu, cvdaed, rcaeops, rcaeopl, rcaeopu, rcaeopd, rctrbga,          &
                        rcvobga, rcstbga, rctrpt, rcaeadm, rcaeros, rcaeadk
  use yoerad,    only : nmode, nozocl, nradfr, nradpfr, nradpla, nrint, nhowinh, novlp, nradf2c, nradc2f,     &
                        nlw, nsw, ntsw, lerad6h, leradhs, lhvolca, lonewsw, lowasyf, lowhsss, loifuec, lrrtm, &
                        lradlp, linhom, raovlp, rbovlp, niceopt, nliqopt, nradip, nradlp, rminice
  use yoerdi,    only : rcardi, rch4, rn2o, ro3, rcfc11, rcfc12, repclc
  use yoerdu,    only : nuaer, ntraer, nout, rcday, r10e, replog, repsc, repsco, repscq, repsct, repscw, diff
  use yoesw,     only : lo3only
  use yoethf,    only : rtwat, rtice
  use yoedbug,   only : ldebug
  use yoephli,   only : lphylin

  use mod_model, only : nlon, nlat, nlonr, nradm, nlev, nlevp1, dlon, dlat, cpd, cpv, sig, sigint, sigalf,    &
                        ps, fmask, phig, albedo, alsn, emisg1, emisg2, fsnow, tskin, qskin, corrdt, corvis, corirr, &
                        fcloud, t, q, qcw, qci, qrain, qsnow, ezer, pzer, g, eps, ccw1, ccw2, cci1, cci2,     &
                        co2ppm, fice, nprocsy, myid, aerosol, ozon, olr, ntop, alont, alatt, rdecli0, reqtim0

  parameter (jp_lon=nlonr, jp_lev=nlev, jp_lw=6, jp_sw=6, jp_nua=24, jp_mode=1, jp_aer=6, jp_levp1=jp_lev+1)

!- reference arrays

  REAL_B :: paer5(jp_lon,jp_aer,jp_lev)
  REAL_B :: palbd5(jp_lon,jp_sw)  , palbp5(jp_lon,jp_sw)
  REAL_B :: paph5(jp_lon,jp_levp1), pap5(jp_lon,jp_lev)
  REAL_B :: pclfr5(jp_lon,jp_lev) , pdp5(jp_lon,jp_lev)
  REAL_B :: pemis5(jp_lon), pemiw5(jp_lon), plsm5(jp_lon), pmu05(jp_lon)
  REAL_B :: pgemu5(jp_lon), pslon5(jp_lon), pclon5(jp_lon),pspgem(jp_lon)
  REAL_B :: pozon5(jp_lon,jp_lev) , pq5(jp_lon,jp_lev)
  REAL_B :: pqiwp5(jp_lon,jp_lev) , pqlwp5(jp_lon,jp_lev)
  REAL_B :: psqiw5(jp_lon,jp_lev) , psqlw5(jp_lon,jp_lev)
  REAL_B :: prlvri5(jp_lon,jp_lev), prlvrl5(jp_lon,jp_lev)
  REAL_B :: pqs5(jp_lon,jp_lev)   , pqrain5(jp_lon,jp_lev), praint5(jp_lon,jp_lev)
  REAL_B :: pth5(jp_lon,jp_levp1) , pt5(jp_lon,jp_lev), pts5(jp_lon)

  REAL_B :: pemit5(jp_lon)
  REAL_B :: pfct5(jp_lon,jp_levp1), pflt5(jp_lon,jp_levp1)
  REAL_B :: pfcs5(jp_lon,jp_levp1), pfls5(jp_lon,jp_levp1)

! box-type arrays (usati entro la RADLSW)

  REAL_B :: pfdct5(jp_lon,jp_levp1), pfdlt5(jp_lon,jp_levp1)
  REAL_B :: pfdcs5(jp_lon,jp_levp1), pfdls5(jp_lon,jp_levp1)
  REAL_B :: pfuct5(jp_lon,jp_levp1), pfult5(jp_lon,jp_levp1)
  REAL_B :: pfucs5(jp_lon,jp_levp1), pfuls5(jp_lon,jp_levp1)

! box-type results (usati entro la RADLSW)

  REAL_B :: aswbox(jp_lon,100), olrbox(jp_lon,100)
  REAL_B :: slwbox(jp_lon,100), sswbox(jp_lon,100)
  REAL_B :: taubox(jp_lon,100), cldbox(jp_lon,100,jp_lev)

  REAL_B :: pfrsod5(jp_lon), psudu5(jp_lon)
  REAL_B :: puvdf5(jp_lon) , pparf5(jp_lon)
  REAL_B :: pnbas5(jp_lon) , pntop5(jp_lon)

!- extra arrays

  REAL_B :: zalb(jp_lon), plat5(jp_lon), plon5(jp_lon)

! dummy integer variables

  INTEGER_M :: kaer, kidia, kfdia, ktdia, klw, ksw, kmode
  INTEGER_M :: iflag, ilay, kbox, kulout, ndump, nflux
  INTEGER_M :: nrad, nindat, nsss, month, iday, iyr, ilwrad, ifin, inhomf
  INTEGER_M :: ilcloud, io3only, icewat, nbox

! dummy real variables

  REAL_B :: prii05, pcco25, zthetoz, zangozc, zdegrad, ztt5, zcondw, zfract, zpp5, zrr5, zddqq, zee, zqq5
  REAL_B :: rtimtr, zcoolc, zcoolt, zheatc, zheatt, zteta

! load external functions

#include "fctast.h"   ! def. funzioni astronomiche
#include "fcttim.h"   ! def. funzioni temporali
#include "fcttre.h"   ! def. funz. termodinamiche

!------------------------------------------------------------------
! DEFINE THE CONFIGURATION
!------------------------------------------------------------------

  kidia = 1
  kfdia = nlonr
  ktdia = 1
  kmode = jp_mode
  kaer  = jp_aer
  klw   = jp_lw
  kbox  = 0           ! no box-type computation
  nbox  = 1

  lhvolca = .false.
  lradlp  = .false.
  lo3only = .false.
  lrrtm   = .false.
  linhom  = .false.
  lphylin = .false.
  lerad6h = .true.
  leradhs = .true.
  lonewsw = .true.

  kulout = 6         ! unita' logica di output per la sucst
  ndump  = 3         ! per evitare stampe
  ldebug = .false.   ! per evitare stampe

! ksw number of SW spectral intervals: 2, 4 or 6

  ksw = jp_sw     ! preced. era 4 - cambia poco l'efficienza, ma con 6 si abbassa la temp. media!

! lascio i param. seg. come sono (validi anche nella vers. 2004)

  nflux   = 6
  nmode   = 0
  nrad    = 1
  nradfr  = -3
  nradpfr = 36
  nradpla = 15
  nrint   = 4
  nradf2c = 1
  nradc2f = 1
  nuaer   = 31
  ntraer  = 19
  novlp   = 1
  nlw     = 6
  ntsw    = 6
  nsw     = ksw
  repsc   = 1.e-12
  repsco  = 1.e-12
  repscq  = 1.e-12
  repsct  = 1.e-12
  repscw  = 1.e-12
  replog  = 1.e-12
  nout    = 6

! call to set-up routines for various coefficients (independent from the vertical resolution)

  call surdi ! inizializza common
  call sulwn ! inizializza common
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

! Initialization routine for RRTM
! Reduce absorption coefficient data from 256 to 140 g-points

  call rrtm_init_140gp

!------------------------------------------------------------------
! DEFINE THE ATMOSPHERIC CASE
!------------------------------------------------------------------

! nindat:  YYYYMMDD date of the simulation, for example 19870701

  iiyr = nyrc
  nindat =  iiyr*10000 + nmonc*100 + ndayc

! nsss: no. sec. del giorno (ma arrotondati al minuto!)

  nsss = nhouc*3600 + nminc*60

! basic constants
! 0 come ultimo argomento: non fa stampe; 1: stampa tutte le cost.

  call sucst (kulout, nindat, nsss, 0)     ! def. costanti tra cui rday necess. alla rtime

! rtimtr: time in sec. after 11:57 gmt of 1 jan (<0 in the morning of 1 jan!?)

  rtimtr = rtime(iiyr, nmonc, ndayc, nsss) ! function in /include/fcttim.h

! call to set-up routines for various coefficients
! These are independent from the vertical resolution

! ilcloud=0 clear-sky; =1 cloudy computations

  ilcloud = 1

  if (ilcloud.eq.1) then

! Choose Optical Properties
! IOPTPROP=0        Water:Fouquart  Ice=Ebert-Curry
! IOPTPROP=1        Water:Slingo    Ice=Ebert-Curry
! IOPTPROP=2        Water:Slingo    Ice=Fu-Liou

  ioptprop = 0   ! operational

    if (ioptprop.eq.0) then
    lowasyf=.false.
    loifuec=.false.
    elseif (ioptprop.eq.1) then
    lowasyf=.true.
    loifuec=.false.
    elseif (ioptprop.eq.2) then
    lowasyf=.true.
    loifuec=.true.
    endif

! ICEWAT=0 Liquid Water; 1=Ice; 2=both

  icewat = 2

! ILWRAD =0 Morcrette, 1991 operational before 20000627'
!        =1 Mlawer et al., 1997, ECMWF-operational since 2000'
!        =2 Morcrette, 1991 original as in ERA-15'

  ilwrad = 1  ! Operativo ECMWF da giu. 2000

! Attenzione: qui nel seguito messe opzioni prese da esempi e valide solo nel caso ILWRAD=1

  lrrtm=.true.
  novlp=1 ! cloud overlap

! NLIQOPT =0  Water LW: Smith-Shi, 1992; SW: Fouquart, 1987
!         =1  Water LW: Savijarvi, 1997; SW: Slingo  , 1989
!         =2  Water LW: Lindner,Li,2000; SW: Slingo  , 1989

  nliqopt=0

! NICEOPT =0  Ice LW: Smith,Shi  , 1992; SW: Ebert-Curry, 1992
!         =1  Ice LW: Ebert,Curry, 1992; SW: Ebert-Curry, 1992
!         =2  Ice LW: Fu,Liou    , 1993; SW: Fu,Liou    , 1993
!         =3  Ice LW: Fu et al.  , 1998; SW: Fu

  niceopt=1

! IRADLP =0 effective radius - liquid as f(Pressure)
!        =1 fixed 10 microns over land, 13 over ocean
!        =2 computed from LWC Martin et al, 1994

  iradlp=2
  nradlp=iradlp
    if (iradlp.le.1) then
    lradlp=.false.
    else
    lradlp=.true.
    endif

! IRADIP =0 fixed effective radius - ice 40 microns
!        =1   f(T)   40 - 130 microns
!        =2   f(T)   30 -  60 microns Jakob-Klein
!        =3   f(T,IWC)  Sun-Rikus, 99
! Non usare opz.IRADIP=3 perchè richiede di setttare MINICE (v. esempi)

  iradip=2
  nradip=iradip

! INHOMF  =0 cloud tau is taken as is
! INHOMF  =1 Tiedtke, 1995  tau x 0.7

  inhomf=1

    if (inhomf.eq.0) then
    linhom=.false.
    elseif (inhomf.eq.1) then
    linhom=.true.

! NHOWINH =1 MT 0.7 factor
! NHOWINH =2 exp(-(sig/tau)^2)
! NHOWINH =3 Cairns, 2000

    nhowinh=1
    endif

  endif

! ICORSFPR =0 nothing done
! ICORSFPR =1 corrected for Delta P
! ICORSFPR =2 corrected for Delta P, ajusted T
! ICORSFPR =3 corrected for Delta P, ajusted T and q

  icorsfpr=0

  call surdi
  call sulwn
  call suswn    (ntsw,ksw)
  call suaerl
  call suaerh
  call suaersn  (ntsw,ksw)
  call suclopn  (ntsw,ksw,nlev)

!- basic constants
! e qui di seguito si vanno a risettare alcune delle cost. prima def. da SUCST.

  rtt    = 273.15
  rtwat  = rtt
  rtice  = rtt-23.
  rg     = g           ! globo
  rcpd   = cpd         ! globo
  ri0    = 1373.
  rsigma = 5.67e-8
  rcday  = rday*rg/rcpd
  diff   = 1.66        ! parametro in yoerdu
  r10e   = 0.4342945   ! parametro in yoerdu
  zomd   = cpv/cpd - 1.

!********************** CO2 *******************************
  pcco25 = co2ppm*1.e-06*44./29.
  if (nozocl.eq.-1) pcco25 = 1.e-18
!**********************************************************

! latitudes and longitudes in radians

  zdegrad = rpi/180.

  do jl = kidia, kfdia
  jlon = jl1+(jl-1)*nradm

#ifdef globo
  plat5(jl)  = (-90. + (infy+jlat-2)*dlat) * zdegrad
  plon5(jl)  = (       (infx+jlon-2)*dlon) * zdegrad
#else
  plat5(jl)  = alatt(jlon,jlat) * zdegrad
  plon5(jl)  = alont(jlon,jlat) * zdegrad  ! long. < 0 not accepted
#endif

  enddo

! definiz. delle variabili atmosferiche

! Pressione sui livelli

  do jklev = 1, nlev
  do jl = kidia, kfdia
  jlon = jl1+(jl-1)*nradm
  pap5(jl,jklev) = pzer*sigint(jklev) - (pzer-ps(jlon,jlat))*sigalf(jklev)
  enddo
  enddo

! Pressione sui semi-livelli

  do jl = kidia, kfdia
  jlon = jl1+(jl-1)*nradm
  paph5(jl,1     ) = 0.
  paph5(jl,nlevp1) = ps(jlon,jlat)
  enddo

  do jklev = 2, nlev
  do jl = kidia, kfdia
  paph5(jl,jklev) = 0.5*(pap5(jl,jklev-1)+pap5(jl,jklev))
  enddo
  enddo

! Temperatura sui livelli

  do jklev = 1, nlev
  do jl = kidia, kfdia
  jlon = jl1+(jl-1)*nradm
  pt5(jl,jklev) = t(jlon,jlat,jklev)
  enddo
  enddo

! temperatura sui semi-livelli: al top si def. uguale a quella del liv. 1
! al bottom e' la Tskin, altrove si fa la media aritmetica

  do jl = kidia, kfdia
  jlon = jl1+(jl-1)*nradm
  pth5(jl,1     ) = pt5(jl,1)
  pth5(jl,nlevp1) = tskin(jlon,jlat)
  pts5(jl)        = tskin(jlon,jlat)
  enddo

  do jklev = 2, nlev
  do jl = kidia, kfdia
  pth5(jl,jklev) = 0.5*(pt5(jl,jklev-1)+pt5(jl,jklev))
  enddo
  enddo

! q e qsat sui livelli

  do jklev = 1, nlev
  do jl = kidia, kfdia
  jlon = jl1+(jl-1)*nradm
  pq5(jl,jklev) = q(jlon,jlat,jklev)
  enddo
  enddo

  do jklev = 1, nlev
  do jl = kidia, kfdia
  ztempsp = pt5(jl,jklev)
  zt0t = rtt/ztempsp
    if(ztempsp.ge.rtt) then
    zesk = ezer*exp(-ccw1*alog(zt0t)+ccw2*(1.-zt0t))  ! partial pressure over water
    else
    zesk = ezer*exp(-cci1*alog(zt0t)+cci2*(1.-zt0t))  ! partial pressure over ice
    endif
  pqs5(jl,jklev) = zesk*eps/(pap5(jl,jklev)+zesk*(eps-1.))
  pq5(jl,jklev) = min (pq5(jl,jklev), 1.05*pqs5(jl,jklev))   ! max q = 1.05*sat.
  enddo
  enddo

! coordinate geografiche - lat e lon in radianti
! PNBAS55 e PNTOP5 sono base e top nubi convettive, ma non usati

  do jl = kidia, kfdia
  pgemu5(jl) = sin(plat5(jl))
  pspgem(jl) = sqrt(1.-pgemu5(jl)**2)
  pclon5(jl) = cos(plon5(jl))
  pslon5(jl) = sin(plon5(jl))
  pnbas5(jl) = 1.
  pntop5(jl) = 1.
  enddo

! PALBD5(JP_LON) SW surface albedo (diffuse)
! PALBP5(JP_LON) SW surface albedo (direct)
! PEMIS5(JP_LON) LW emissivity outside the LW window region
! PEMIW5(JP_LON) LW emissivity within the 8-12.5 micron window region
! PLSM5(JP_LON)  land/sea mask   (1.= land   0.=ocean)
! IN EMISG1 AND EMISG2, SEA IS INCLUDED, SNOW IS NOT

  do jl = kidia, kfdia
  jlon = jl1+(jl-1)*nradm
  plsm5 (jl) = 1. -(fmask(jlon,jlat)-fice(jlon,jlat))
  zalb  (jl) = albedo(jlon,jlat)*(1.-fsnow(jlon,jlat)) + alsn(jlon,jlat)*fsnow(jlon,jlat)
  pemis5(jl) = emisg1(jlon,jlat)
  pemiw5(jl) = emisg2(jlon,jlat)
  enddo

! albedo

  do jnu = 1, ksw
  do jl = kidia, kfdia
  palbd5(jl,jnu) = zalb(jl)
  palbp5(jl,jnu) = (zalb(jl)-.07*(1.-plsm5(jl)))*1.26/(1.+.52*pmu05(jl)) + (1.-plsm5(jl))*.0397/(1.1*pmu05(jl)**1.4+.15)
  palbp5(jl,jnu) = min (palbp5(jl,jnu), .999)
  enddo
  enddo

! calcolo cost. astronomiche come ECMWF

  zteta  = rteta(rtimtr) ! zteta tempo in fraz. di anno dalle 11:57 del 1 genn.
  rdeaso = rrs(zteta)    ! dist. terra-sole
  rdecli = rds(zteta)    ! decl. terra

! reqtim dovrebbe una correz. astron. (in s) per la durata del giorno
! rispetto all'andamento stagionale puramente sinusoidale
! ret: 'eq. del tempo' (funz. astron. in fctst.h)

  reqtim = ret(zteta)
  rdecli0 = rdecli
  reqtim0 = reqtim
  rsovr  = reqtim+float(nsss)          ! nsss no. di sec. dalle 00 di ogni giorno
  rwsovr = rsovr*2.*rpi/rday           ! rday=86400., rpi e' pigreco
  prii05 = ri0*rea*rea/(rdeaso*rdeaso) ! "cost. solare" con var. stag.
  rcodec = cos(rdecli)
  rsidec = sin(rdecli)
  rcovsr = cos(rwsovr)
  rsivsr = sin(rwsovr)

  do jl = kidia, kfdia
  pmu05(jl) = max(rsidec*pgemu5(jl) - rcodec*rcovsr*pspgem(jl)*pclon5(jl) &
            + rcodec*rsivsr*pspgem(jl)*pslon5(jl), 0.d0)
  enddo

  do jklev = 1, nlev
  do jl = kidia, kfdia
  jlon = jl1+(jl-1)*nradm
  pclfr5(jl,jklev)  = fcloud(jlon,jlat,jklev)
  pqiwp5(jl,jklev)  = qci(jlon,jlat,jklev)!+qsnow(jlon,jlat,jklev)
  pqlwp5(jl,jklev)  = qcw(jlon,jlat,jklev)
  pdp5  (jl,jklev)  = paph5(jl,jklev+1)-paph5(jl,jklev)
!  pqrain5(jl,jklev) = qrain(jlon,jlat,jklev) ! non usato!
  pqrain5(jl,jklev) = 0.
  praint5(jl,jklev) = 0.

  pozon5(jl,  jklev) = ozon   (jlon,jlat,jklev  )
  paer5 (jl,1,jklev) = aerosol(jlon,jlat,jklev,1)
  paer5 (jl,2,jklev) = aerosol(jlon,jlat,jklev,2)
  paer5 (jl,3,jklev) = aerosol(jlon,jlat,jklev,3)
  paer5 (jl,4,jklev) = aerosol(jlon,jlat,jklev,4)
  paer5 (jl,5,jklev) = aerosol(jlon,jlat,jklev,5)
  paer5 (jl,6,jklev) = aerosol(jlon,jlat,jklev,6)

  enddo
  enddo

  do jl = 1, kidia, kfdia
  nextra = nint(3.*((abs(plat5(jl))/(.5*3.14))**2))
  pclfr5(jl,1:ntop+nextra) = 0.
  pqiwp5(jl,1:ntop+nextra) = 0.
  pqlwp5(jl,1:ntop+nextra) = 0.
  enddo

! Le quantita' seg. dip. da rapporti di mix. ratio, settate come negli esempi, ma sembrano ininfluenti

  psqiw5  = 1.
  psqlw5  = 1.
  prlvri5 = 0.
  prlvrl5 = 0.

  call radlsw (kidia, kfdia, kfdia, ktdia, nlev, kmode, kaer, kbox, nbox, ndump, ilwrad, prii05, &
               paer5, palbd5, palbp5, paph5, pap5, pcco25, pclfr5, pdp5, pemis5, pemiw5, plsm5, pmu05, pozon5, &
               pq5, pqiwp5, pqlwp5, psqiw5, psqlw5, pqs5, pqrain5, praint5, &
               prlvri5, prlvrl5, pth5, pt5, pts5, pnbas5, pntop5, &
               pemit5, pfct5, pflt5, pfcs5, pfls5, pfrsod5, psudu5, puvdf5, pparf5, &
               pfdct5, pfuct5, pfdlt5, pfult5, pfdcs5, pfucs5, pfdls5, pfuls5, &
               aswbox, olrbox, slwbox, sswbox, taubox, cldbox)

!---------------------------------------------------------------
! OUTPUT VARIABLES
!---------------------------------------------------------------

! si parte dal secondo liv. in alto per evitare il problema al top

  do jklev = 2, nlev
  do jl = kidia, kfdia
  jlon = jl1 + (jl-1)*nradm
  zdstcp = g/(cpd*(1.+zomd*q(jlon,jlat,jklev))*(paph5(jl,jklev+1)-paph5(jl,jklev)))
  corrdt(jlon,jlat,jklev) = zdstcp*(pflt5(jl,jklev) -pflt5(jl,jklev+1) +pfls5(jl,jklev) -pfls5(jl,jklev+1))
  enddo
  enddo

! si def. (arbitr.) al liv. 1 una frazione della tendenza calcolata per il liv. 2

  do jl = kidia, kfdia
  jlon = jl1+(jl-1)*nradm
  zdstcp = g/(cpd*(1.+zomd*q(jlon,jlat,2))*(paph5(jl,3)-paph5(jl,2)))
  corrdt(jlon,jlat,1) = 0.7*zdstcp*(pflt5(jl,2) -pflt5(jl,3) +pfls5(jl,2)-pfls5(jl,3))
  enddo

! surface fluxes of visible and infrared radiation (positive downward)

  do jl = kidia, kfdia
  jlon = jl1+(jl-1)*nradm
  corvis(jlon,jlat) = pfls5(jl,nlevp1)
  corirr(jlon,jlat) = pflt5(jl,nlevp1)  + 12. ! known bias
  olr   (jlon,jlat) = pflt5(jl,1)
  enddo

  return
  end subroutine radintec
#endif
#endif
!##################################################################################################################
      subroutine cloudfr

! Defines fcloud: fraction of cloudy sky, as a function of
! explicit cloud water, cloud ice, snow and relative humidity.
! It must be stricly 0 > fcloud < 1.
! Cloud fraction is reduced depending on local (moist) Richardson number
! and pbl stability.

      use mod_model, only : nlon, nlat, nlev, nlonm1, nlatm1, nlevp1, ntop, u, v, t, &
                            q, qcw, qci, tvirt, fcloud, qsnow, rich, sigint, ps,     &
                            sigalf, rgm, huc, tskin, qskin, phig, phi, phih, fmask,  &
                            rd, cpv, cpd, g, ep, eps, rdrcp, pzer, tzer, ezer, alp0, &
                            ccw1, ccw2, ylwv
      implicit none
      real zstabg(nlon,nlat), zindn(nlon,nlat), ztevg(nlon,nlat)
      real zconv(nlon,nlev), zqs(nlon,nlev), ztetav(nlon,nlev)
      real dlogthe(nlon,nlevp1), zhlev(nlon,nlevp1)
      real zstabr, zdstabg, zqcrit, zepneb, zomd, fc1, fc2, cloudd, zwe, ritop, fcrit, cltop
      real zppp, ztmp, zesat, zqs1, zqs2, zqsa, zrh, zrh1, zrc
      real zzpp, zt0t, zesk, dthd, dthm, r_up, r_do, r_av, t_av, theta_up, theta_do, zrhm
      real zaa, zza, zcof, zua, zva, zmod2, zconvg, ztebar, beta, zri, zrdz, zdtdz
      real zdudz, zdvdz, zshear, zbuoy, zrich
      integer jlat, jlon, jklev, jk, ktop, jlatp1, ifl
      parameter (zepneb = 1.e-10, zqcrit = 2.2e-4)

      zomd = cpv/cpd-1.

! re-computation of the (moist) Richardson number (as in vdiff), needed to reduce cloud fraction
! (the Richardson no. computed in vdiff is not updated when this subroutine is called -
! this implies differences of cloud fraction with those computed in postbolam).

    do jlat = 1, nlat
    jlatp1 = min (jlat+1, nlat) ! necessary for globo

    do jklev = 1, nlev   ! loop on half levels
    do jlon = 2, nlonm1
     zzpp = pzer*sigint(jklev) - (pzer-ps(jlon,jlat))*sigalf(jklev)
     zconv(jlon,jklev) = exp(rdrcp*(alp0-log(zzpp)))
     zt0t = tzer/t(jlon,jlat,jklev)
     zesk = ezer*exp(-ccw1*log(zt0t)+ccw2*(1.-zt0t))     ! partial pressure over water
     zqs(jlon,jklev)     = zesk*eps/(zzpp+zesk*(eps-1.))
     ztetav(jlon,jklev)  = tvirt(jlon,jlat,jklev)*zconv(jlon,jklev)
     zhlev(jlon,jklev+1) = (phih(jlon,jlat,jklev+1)-phig(jlon,jlat))/g
    enddo
    enddo

!  computation of dry and moist (saturated) static stability as Durran & Klemp (1982)
!  (specific humidity is used in place of mixing ratio)

    do jklev = 2, nlev   ! loop on half levels
    do jlon = 2, nlonm1
     dthd = 2.*(ztetav(jlon,jklev-1)-ztetav(jlon,jklev))/(ztetav(jlon,jklev-1)+ztetav(jlon,jklev) ) ! dry
     r_up = zqs(jlon,jklev-1)
     r_do = zqs(jlon,jklev  )
     r_av = 0.5*(r_up + r_do)
     t_av = 0.5*(t(jlon,jlat,jklev) + t(jlon,jlat,jklev-1))
     theta_up = t(jlon,jlat,jklev-1)*zconv(jlon,jklev-1)
     theta_do = t(jlon,jlat,jklev  )*zconv(jlon,jklev  )
     zaa = (1. + ylwv*r_av/(rd*t_av))/(1. + eps*ylwv**2*r_av/(cpd*rd*t_av**2))
     dthm = zaa*((theta_up-theta_do)*2./(theta_up + theta_do) + ylwv/(cpd*t_av)*(r_up - r_do)) &
            - q(jlon,jlat,jklev-1) - qcw(jlon,jlat,jklev-1) - qci(jlon,jlat,jklev-1)           &
            + q(jlon,jlat,jklev  ) + qcw(jlon,jlat,jklev  ) + qci(jlon,jlat,jklev  )

! average relative humidity computed giving some more weight to the layer below

     zrhm = 0.55*q(jlon,jlat,jklev)/zqs(jlon,jklev) + 0.45*q(jlon,jlat,jklev-1)/zqs(jlon,jklev-1)
     zcof = max(-24.5 + 25.*zrhm, 0.)    ! zcof=0. for rh<=0.98, zcof=0.5 for rh=1
     zcof = min(zcof, .8)
     dlogthe(jlon,jklev) = zcof*dthm + (1.-zcof)*dthd  ! effective stability near saturation
    enddo
    enddo

    do jlon = 2, nlonm1
    zza   = (phi(jlon,jlat,nlev)-phig(jlon,jlat))/g + rgm(jlon,jlat) !  sigma=1 is located at z=rgm
    zua   = .5*(u(jlon,jlat,nlev)+u(jlon-1,jlat,nlev))
    zva   = .5*(v(jlon,jlat,nlev)+v(jlon,jlatp1,nlev))
    zmod2 = zua**2 + zva**2 + 0.07

!  virtual potential temperature computed with skin temperature

    zconvg = exp(rdrcp*(alp0-log(ps(jlon,jlat))))*(1.+ep*qskin(jlon,jlat))
    ztevg(jlon,jlat) = tskin(jlon,jlat)*zconvg

!  bulk Richardson number

    ztebar =.5*(ztevg(jlon,jlat)+ztetav(jlon,nlev))
    beta   = g/ztebar
    zri    = zza*beta*(ztetav(jlon,nlev)-ztevg(jlon,jlat))/zmod2
    rich(jlon,jlat,nlev+1) =  min(zri, 500.)
    enddo

    do jklev = nlev, 2, -1
    do jlon = 2, nlonm1
    zrdz  = g/(phi(jlon,jlat,jklev-1)-phi(jlon,jlat,jklev))
    zdtdz = dlogthe(jlon,jklev)*zrdz
    zdudz = .5*(u(jlon,jlat,jklev-1)+u(jlon-1,jlat,jklev-1)-u(jlon,jlat,jklev)-u(jlon-1,jlat,jklev))*zrdz
    zdvdz = .5*(v(jlon,jlat,jklev-1)+v(jlon,jlatp1,jklev-1)-v(jlon,jlat,jklev)-v(jlon,jlatp1,jklev))*zrdz
    zshear= zdudz**2 + zdvdz**2
    zbuoy = g*dlogthe(jlon,jklev)*zrdz
    zrich = min (zbuoy/(zshear + 1.e-6), 500.)
    rich(jlon,jlat,jklev) = zrich
    enddo
    enddo

    enddo    !  end of loop in latitude
    rich(:,:,1) = rich(:,:,2)

! Initialization for computing Geleyn stability in the pbl:
! weighted average of static energy between the surface and the lowest model level.
! More weight is given to surface t and q (with respect to air) over sea than over land,
! to reduce clouds more over sea than over land in unstable pbl.

      do jlat = 1, nlat
      do jlon = 2, nlonm1
      zwe = 0.1 + fmask(jlon,jlat)*0.6
      zstabg(jlon,jlat) = zwe*(cpd*(1. + zomd*qskin(jlon,jlat))*tskin(jlon,jlat) + phig(jlon,jlat)) +         &
                         (1.-zwe)*(cpd*(1. + zomd*q(jlon,jlat,nlev))*t(jlon,jlat,nlev) + phi(jlon,jlat,nlev))
      zindn(jlon,jlat) = 1.
      enddo
      enddo

! Preliminary definition of fcloud

      do jklev = ntop, nlev
      do jlat = 1, nlat
      do jlon = 2, nlonm1
      zppp = pzer*sigint(jklev)-(pzer-ps(jlon,jlat))*sigalf(jklev)
      ztmp = t(jlon,jlat,jklev)
      call comp_esk(zesat,zqs1,ztmp,zppp,2)  ! blended saturation
      call comp_esk(zesat,zqs2,ztmp,zppp,3)  ! sat. to water below 0
      zqsa = 0.70*zqs1+0.30*zqs2             ! ad hoc to limit contrib. to high clouds
      zrh1 = min(1., q(jlon,jlat,jklev)/zqsa)
      zrh = (zrh1-huc(jklev))/(1.-huc(jklev))
      zrh = min(1., zrh)
      zrh = max(0., zrh)

      fcloud(jlon,jlat,jklev) = (max(0.15*zrh, (qcw(jlon,jlat,jklev) + 0.77*qci(jlon,jlat,jklev) + &
                                 0.2*qsnow(jlon,jlat,jklev))/zqcrit))**0.7
      enddo
      enddo
      enddo

! Reduction of clouds in unstable layers

      fcrit = 0.3
      do jklev = nlev, ntop, -1
      do jlat = 1, nlat
      do jlon = 2, nlonm1

! Geleyn stability, based on (dry, virtual) static energy profile (non-local, for pbl only)

      zdstabg = zstabg(jlon,jlat) - cpd*(1. + zomd*q(jlon,jlat,jklev))*t(jlon,jlat,jklev)-phi(jlon,jlat,jklev)
      zindn(jlon,jlat) = zindn(jlon,jlat)*max(0., sign(1., zdstabg))

! Local stability derived from the Richardson no.

      zstabr = 0.5*rich(jlon,jlat,jklev+1) + 0.5*rich(jlon,jlat,jklev)

      fc1 = 1.
      fc2 = 1.
      if(zstabr.lt..25.and.zstabr.gt.0.) fc1 = 0.875 + 0.5*zstabr
      if(zstabr.gt.0..and.rich(jlon,jlat,jklev+1).lt.0.) fc1 = 0.8
      if(zstabr.le.0.) fc1 = 0.7

! Over the sea: computation of cloud depth and cloud top to distinguish between cumuli
! and stratocumuli, the latter to be more strongly reduced by Geleyn method

      if(fmask(jlon,jlat).ge.0.5) then ! sea
      cloudd = 0.
      ktop = nlev
      ifl = 0
      if(zindn(jlon,jlat).gt..9) then
       do jk = jklev-1, nlev/2, -1
       if(fcloud(jlon,jlat,jk).gt.fcrit) then
       cloudd = cloudd + 0.5*(phi(jlon,jlat,jk-1)-phi(jlon,jlat,jk+1))/g
       ktop = jk
       ifl = 1
       else
       exit
       endif
       enddo
       if(ifl.eq.1) then
        ritop = rich(jlon,jlat,ktop)  ! Ri at half lev. above cloud top
        cltop = (phi(jlon,jlat,ktop)-phig(jlon,jlat))/g  ! height of cloud top
        if(cloudd.gt.480..or.cltop.gt.1550..or.ritop.lt.0.5) zindn(jlon,jlat) = 0.55 ! cumuli, not stratocumuli
       endif
      endif
#ifdef globo
      zrc = 0.52
#else
      zrc = 0.40
#endif
      else  ! land
      zrc = 0.22
      endif !  sea or land

      fc2 = 1. - zrc*zindn(jlon,jlat)
      fcloud(jlon,jlat,jklev) = min(fc1, fc2)*fcloud(jlon,jlat,jklev)
      fcloud(jlon,jlat,jklev) = max(zepneb, min(.999995, fcloud(jlon,jlat,jklev)))
      enddo
      enddo
      enddo

      return
      end subroutine cloudfr
!##################################################################################################################
   subroutine ccloud

! Computes total cloud cover - uses subr. ccolumn

   use mod_model, only : nlon, nlat, nlev, nlonm1, nlatm1, cloudt, fcloud, ntop, myid, nprocsy
   real zprod(nlon), zcol(nlev-ntop+1)

   jstart = 2
   jend = nlatm1

#ifdef globo
   if (mod(myid,nprocsy).eq.0) jstart = 1
   if (mod(myid,nprocsy).eq.nprocsy-1) jend = nlat
#endif

   do jlat = jstart, jend
   do jlon = 2, nlonm1
   zcol(1:nlev-ntop+1) = fcloud(jlon,jlat,ntop:nlev)
   call ccolumn(zcol, nlev-ntop+1, zprod(jlon))
   cloudt(jlon,jlat) = max (0., 1.-zprod(jlon))
   cloudt(jlon,jlat) = min (cloudt(jlon,jlat), 1.)
   enddo
   enddo

   return
   end subroutine ccloud
!##################################################################################################################
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
!##################################################################################################################
    subroutine sea_surface

!  Slab sea scheme
!  Sea includes thin sea ice - thick sea ice is treated as land glacier

    use mod_model, only : tg, qgice, tskin, qskin, nlevg, t, q, hflux, qflux, frvis, frirr, fsnow, snow,  &
                          dtstep, nlev, ccw1, ccw2, cci1, cci2, tzer, ezer, ps, eps, &
                          roscd, nlonm1, nlatm1, fmask, fice, rgmd, mask_soil, kturb_surf_h, kturb_surf_q, tvirt, phi, phig, &
                          pzer, g, cpd, sigalf, sigint, rd, rdrcp, alp0, ep, myid
    implicit none

    real, parameter :: ztice=tzer-0.5
    integer jlon, jlat
    real zt0t, ztotfl, zesk, zeskw, zeski, zqsat, zrc, zeffw, &
         zhflux, zqflux, zrho_skin, zzpp, zconv, zconvg, ztetav, ztevg

    do 1000 jlat = 2, nlatm1
    do 1000 jlon = 2, nlonm1

!--------------------------------
!    if (fmask(jlon,jlat).ge.0.5.and.fice(jlon,jlat).lt.0.8) then  ! over sea (including thin sea ice)
    if (mask_soil(jlon,jlat) == 0) then  ! over sea (excluding thick sea ice)
!--------------------------------

      zt0t = tzer/tg(jlon,jlat,1)
      if (tg(jlon,jlat,1).ge.tzer) then
      zesk = ezer*exp(-ccw1*log(zt0t)+ccw2*(1.-zt0t))       ! partial pressure over water
      else
      zeskw = ezer*exp(-ccw1*log(zt0t)+ccw2*(1.-zt0t))      ! partial pressure over water
      zeski = ezer*exp(-cci1*log(zt0t)+cci2*(1.-zt0t))      ! partial pressure over ice
      zesk  = zeski*fice(jlon,jlat) + zeskw*(1.-fice(jlon,jlat))
      endif
      zqsat = zesk*eps/(ps(jlon,jlat)+zesk*(eps-1.))

#ifdef globo
      zeffw = .043/(.043 + roscd(jlon,jlat) + 1.e-9)
      qskin(jlon,jlat) = zeffw*zqsat + (1.-zeffw)*q(jlon,jlat,nlev)
!      qskin(jlon,jlat) = zqsat*.963*(1.-fice(jlon,jlat)) + zqsat*fice(jlon,jlat)
#else
      zeffw = .058/(.058 + roscd(jlon,jlat) + 1.e-9)
      qskin(jlon,jlat) = zeffw*zqsat*0.97 + (1.-zeffw)*q(jlon,jlat,nlev)
#endif

      zconvg = exp(rdrcp*(alp0-log(ps(jlon,jlat))))*(1.+ep*qskin(jlon,jlat))
      ztevg = tskin(jlon,jlat)*zconvg
      zzpp = pzer*sigint(nlev) - (pzer-ps(jlon,jlat))*sigalf(nlev)
      zconv = exp(rdrcp*(alp0-log(zzpp)))
      ztetav = tvirt(jlon,jlat,nlev)*zconv
      zrho_skin = ps(jlon,jlat)/(rd*tskin(jlon,jlat)*(1.+ep*qskin(jlon,jlat)))
      zhflux = kturb_surf_h(jlon,jlat)*cpd*zrho_skin/((phi(jlon,jlat,nlev)-phig(jlon,jlat))/g)*(ztevg-ztetav)
      zqflux = kturb_surf_q(jlon,jlat)*zrho_skin/((phi(jlon,jlat,nlev)-phig(jlon,jlat))/g)*(qskin(jlon,jlat)-q(jlon,jlat,nlev))

!      ztotfl = hflux(jlon,jlat) +2.5008e6*qflux(jlon,jlat) -frvis(jlon,jlat) -frirr(jlon,jlat)
      ztotfl = zhflux +2.5008e6*zqflux -frvis(jlon,jlat) -frirr(jlon,jlat)

      zrc = 1./(1. + 1.5e2*rgmd(jlon,jlat))
      if (fice(jlon,jlat).le.0.5) then
      tg(jlon,jlat,1) = tg(jlon,jlat,1) - dtstep*( .25e-7*zrc*ztotfl + .8e-5*(tg(jlon,jlat,1)-tg(jlon,jlat,nlevg)) )
      else             !  reduced heat capacity and heat diffusivity over sea-ice
      tg(jlon,jlat,1) = tg(jlon,jlat,1) - dtstep*( .48e-7*ztotfl + .5e-5* (tg(jlon,jlat,1)-tg(jlon,jlat,nlevg)) )
      endif
        if(fice(jlon,jlat).lt.0.5) then
        snow(jlon,jlat)  = 0.
        fsnow(jlon,jlat) = 0.
        endif
      tskin(jlon,jlat) = tg(jlon,jlat,1)

      endif

 1000 continue

      return
      end subroutine sea_surface
!##################################################################################################################
 subroutine interp (alfa, ex1, ex2, npi, xi, g, x, f, nval)

!  Interpolates with splines with tension in one dimension.
!  The spline is defined imposing that the second derivative is the average
!  of second derivatives computed at the two adjacent points.
!  At interval extremes the second derivative is assumed null.
!  This subrout. also extrapolates out of the interval where the input funtion g is defined

!  Input:  function g defined at coordinates xi (caution: can be changed by this subrout.)
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
!##################################################################################################################
      subroutine tqmass (ps, q, nlon, nlat, nlev, hxt, dlon, dlat, myid, nprocsx, nprocsy, &
                         dsig, dsigalf, pzer, g, tom)

! For globo

      real(4) ps(nlon,nlat), hxt(nlat), q(nlon,nlat,nlev), dsig(nlev), dsigalf(nlev)

      pi = abs(acos(-1.))
      tom = 0.

      do jklev = 2, nlev
      do jlat = 2, nlat-1
      do jlon = 2, nlon-1
      zdpkp = pzer*dsig(jklev) - (pzer-ps(jlon,jlat))*dsigalf(jklev)
!      tom = tom + ps(jlon,jlat)*hxt(jlat)*dlon*dlat*(pi/180.)**2
      tom = tom + hxt(jlat)*dlon*dlat*(pi/180.)**2*q(jlon,jlat,jklev)*zdpkp/g
      enddo
      enddo
      enddo

      if (mod(myid,nprocsy).eq.0) then
      do jklev = 2, nlev
      zdpkp = pzer*dsig(jklev) - (pzer-ps(1,1))*dsigalf(jklev)
      tom = tom + q(1,1,jklev)*zdpkp/g*2.*pi*(1.-cos(dlat*pi/360.))/float(nprocsx)
      enddo
      endif

      if (mod(myid,nprocsy).eq.nprocsy-1) then
      do jklev = 2, nlev
      zdpkp = pzer*dsig(jklev) - (pzer-ps(1,nlat))*dsigalf(jklev)
!      tom = tom + ps(1,nlat)*2.*pi*(1.-cos(dlat*pi/360.))/float(nprocsx)
      tom = tom + q(1,nlat,jklev)*zdpkp/g*2.*pi*(1.-cos(dlat*pi/360.))/float(nprocsx)
      enddo
      endif

      tom = tom/(4.*pi)

      return
      end subroutine tqmass
!##################################################################################################################
    subroutine polavert (p)

!  For globo, polar average of T variables - it is assumed that ghostlines are not updated

    use mod_model, only: gnlon, nlon, nlat, nlev, nlonm1, nlatm1, nprocsx, npolcap, filtt1,  &
                         hxt, myid, ip_e, ip_n, ip_s, ip_w, ip_null

    real(4) p(nlon,nlat,nlev), p1(nlon,nlev)
    real, parameter :: zsq2 = (gnlon-2)/sqrt(2.)

!------------------------
!  Polar filter
!------------------------

    if (ip_s.eq.ip_null) then
    do jlat = 1, npolcap+1
    znt = zsq2*hxt(jlat)  ! see 2011 paper on Weather and Forecasting
    call slofou1 (p(1:nlon,jlat,1:nlev), znt, .false.)
    enddo
    endif
    if (ip_n.eq.ip_null) then
    do jlat = nlat-npolcap, nlat
    znt = zsq2*hxt(jlat)
    call slofou1 (p(1:nlon,jlat,1:nlev), znt, .false.)
    enddo
    endif

!-------------------------------
!  Update meridional ghostlines
!-------------------------------

#ifdef mpi
      if (nprocsx.eq.1) then
#endif
    p(1,:,:) = p(nlonm1,:,:)
    p(nlon,:,:) = p(2,:,:)
#ifdef mpi
      else
        call u_ghost (p(nlonm1,:,:), ip_e, p(1   ,:,:), ip_w, nlat*nlev)
        call u_ghost (p(2     ,:,:), ip_w, p(nlon,:,:), ip_e, nlat*nlev)
      endif
#endif

    do jlat = 1, nlat
    nfilt = filtt1(jlat) + .99999
    anu = 0.
    if (nfilt.gt.0) anu = filtt1(jlat)/float(nfilt)
    do n = 1, nfilt
      do jklev = 1, nlev
      do jlon = 2, nlonm1
      p1(jlon,jklev) = .25*(p(jlon-1,jlat,jklev)+p(jlon+1,jlat,jklev))-.5*p(jlon,jlat,jklev)
      enddo
      enddo
#ifdef mpi
      if (nprocsx.eq.1) then
#endif
    p1(1   ,1:nlev) = p1(nlonm1,1:nlev)
    p1(nlon,1:nlev) = p1(2     ,1:nlev)
#ifdef mpi
      else
        call u_ghost (p1(nlonm1,1:nlev), ip_e, p1(1   ,1:nlev), ip_w, nlev)
        call u_ghost (p1(2     ,1:nlev), ip_w, p1(nlon,1:nlev), ip_e, nlev)
      endif
#endif
    p(1:nlon,jlat,1:nlev) = p(1:nlon,jlat,1:nlev) + anu*p1(1:nlon,1:nlev)
    enddo
    enddo

    return
    end subroutine polavert
!##################################################################################################################
    subroutine polaveruv

! For globo

    use mod_model, only: gnlon, nlon, nlat, nlev, nlonm1, nlatm1, nprocsx, npolcap, u, v, filtu1, filtv1, &
                         hxt, myid, ip_e, ip_n, ip_s, ip_w, ip_null

    real(4) p1(nlon,nlev)
    real, parameter :: zsq2 = (gnlon-2)/sqrt(2.)

!---------------
!  Polar average
!---------------

    if (ip_s.eq.ip_null) then
    do jlat = 1, npolcap+1    ! U points
    znt = zsq2*hxt(jlat) + 1.
    call slofou1 (u(1:nlon,jlat,1:nlev), znt, .true.)
    enddo
    do jlat = 2, npolcap+1    ! V points
    znt = .5*zsq2*(hxt(jlat)+hxt(jlat-1)) + 1.
    call slofou1 (v(1:nlon,jlat,1:nlev), znt, .true.)
    enddo
    endif
    if (ip_n.eq.ip_null) then
    do jlat = nlat-npolcap, nlat      ! U points
    znt = zsq2*hxt(jlat) + 1.
    call slofou1 (u(1:nlon,jlat,1:nlev), znt, .true.)
    enddo
    do jlat = nlat-npolcap+1, nlat    ! V points
    znt = .5*zsq2*(hxt(jlat)+hxt(jlat-1)) + 1.
    call slofou1 (v(1:nlon,jlat,1:nlev), znt, .true.)
    enddo
    endif

#ifdef mpi
      if (nprocsx.eq.1) then
#endif
    u(1   ,:,:) = u(nlonm1,:,:)
    u(nlon,:,:) = u(2     ,:,:)
    v(1   ,:,:) = v(nlonm1,:,:)
    v(nlon,:,:) = v(2     ,:,:)
#ifdef mpi
      else
        call u_ghost (u(nlonm1,:,:), ip_e, u(1   ,:,:), ip_w, nlat*nlev)
        call u_ghost (u(2     ,:,:), ip_w, u(nlon,:,:), ip_e, nlat*nlev)
        call u_ghost (v(nlonm1,:,:), ip_e, v(1   ,:,:), ip_w, nlat*nlev)
        call u_ghost (v(2     ,:,:), ip_w, v(nlon,:,:), ip_e, nlat*nlev)
      endif
#endif

    do jlat = 1, nlat
    nfilt = filtu1(jlat) + .99999
    anu = 0.
    if (nfilt.gt.0) anu = filtu1(jlat)/float(nfilt)
    do n = 1, nfilt
      do jklev = 1, nlev
      do jlon = 2, nlonm1
      p1(jlon,jklev) = .25*(u(jlon-1,jlat,jklev)+u(jlon+1,jlat,jklev))-.5*u(jlon,jlat,jklev)
      enddo
      enddo
#ifdef mpi
      if (nprocsx.eq.1) then
#endif
    p1(1   ,1:nlev) = p1(nlonm1,1:nlev)
    p1(nlon,1:nlev) = p1(2     ,1:nlev)
#ifdef mpi
      else
        call u_ghost (p1(nlonm1,1:nlev), ip_e, p1(1   ,1:nlev), ip_w, nlev)
        call u_ghost (p1(2     ,1:nlev), ip_w, p1(nlon,1:nlev), ip_e, nlev)
      endif
#endif
    u(1:nlon,jlat,1:nlev) = u(1:nlon,jlat,1:nlev) + anu*p1(1:nlon,1:nlev)
    enddo

    nfilt = filtv1(jlat) + .99999
    anu = 0.
    if (nfilt.gt.0) anu = filtv1(jlat)/float(nfilt)
    do n = 1, nfilt
      do jklev = 1, nlev
      do jlon = 2, nlonm1
      p1(jlon,jklev) = .25*(v(jlon-1,jlat,jklev)+v(jlon+1,jlat,jklev))-.5*v(jlon,jlat,jklev)
      enddo
      enddo
#ifdef mpi
      if (nprocsx.eq.1) then
#endif
    p1(1   ,1:nlev) = p1(nlonm1,1:nlev)
    p1(nlon,1:nlev) = p1(2     ,1:nlev)
#ifdef mpi
      else
        call u_ghost (p1(nlonm1,1:nlev), ip_e, p1(1   ,1:nlev), ip_w, nlev)
        call u_ghost (p1(2     ,1:nlev), ip_w, p1(nlon,1:nlev), ip_e, nlev)
      endif
#endif
    v(1:nlon,jlat,1:nlev) = v(1:nlon,jlat,1:nlev) + anu*p1(1:nlon,1:nlev)
    enddo
    enddo

    return
    end subroutine polaveruv
!##################################################################################################################
      subroutine slofou1 (p, nt, nlwind)

! For globo

      use mod_model, only: nlon, gnlon, nlonm1, nlev, nprocsx, nprocsy, myid, npolcap, dlon, snt, cst
      implicit none
      logical nlwind
      integer :: jlon, jklev, ntot, n, jpr, iprocs, ip_first, ierr, comm, tag1, tag2
      real(4), dimension(nlon,nlev)   :: p
      real(4), dimension(10*(npolcap+1),nlev)   :: zps, zs
      real(4), dimension(0:10*(npolcap+1),nlev) :: zpc, zc
      real(4) psmoo(gnlon), nt

#ifdef mpi
      include 'mpif.h'
      integer :: status(mpi_status_size)

      comm  = mpi_comm_world
#endif
      iprocs = nprocsx*nprocsy
      ip_first = mod(myid,nprocsy)
#ifdef mpi
        tag1 = 1
        tag2 = 2
#endif

      ntot = nint(3.*nt)
      do n = 1, ntot
      psmoo(n) = exp(-float(n)**2/nt**2)*dlon/180.
      enddo
      if (nlwind.and.ntot.eq.3) ntot = 1  ! only for wind at poles
      if (nlwind) psmoo(1) = dlon/180.

      zpc(0,:) = 0.
      do jklev = 1, nlev
      do jlon = 2, nlonm1
      zpc(0,jklev) = zpc(0,jklev) + p(jlon,jklev)
      enddo
      zpc(0,jklev) = zpc(0,jklev)/float(gnlon-2)
      enddo
      if (nlwind.and.ntot.eq.1) zpc(0,:) = 0.  ! only for wind at poles

      do 10 n = 1, ntot
      do jklev = 1, nlev
      zps(n,jklev) = 0.
      zpc(n,jklev) = 0.
      do jlon = 2 , nlonm1
      zps(n,jklev) = zps(n,jklev)+p(jlon,jklev)*snt(jlon,n)
      zpc(n,jklev) = zpc(n,jklev)+p(jlon,jklev)*cst(jlon,n)
      enddo
      zps(n,jklev) = zps(n,jklev)*psmoo(n)
      zpc(n,jklev) = zpc(n,jklev)*psmoo(n)
      enddo
 10   continue

!  global average of zpc, zps and broadcast

      if (nprocsx.gt.1) then
#ifdef mpi
        if (myid.ge.nprocsy) then
          call mpi_send (zpc(0:ntot,1:nlev), (ntot+1)*nlev, mpi_real, ip_first, tag1, comm, ierr)
          call mpi_send (zps(1:ntot,1:nlev), (ntot  )*nlev, mpi_real, ip_first, tag1, comm, ierr)
          call mpi_recv (zpc(0:ntot,1:nlev), (ntot+1)*nlev, mpi_real, ip_first, tag2, comm, status, ierr)
          call mpi_recv (zps(1:ntot,1:nlev), (ntot  )*nlev, mpi_real, ip_first, tag2, comm, status, ierr)
        else
          do jpr = ip_first+nprocsy, iprocs-1, nprocsy
            call mpi_recv (zc(0:ntot,1:nlev), (ntot+1)*nlev, mpi_real, jpr, tag1, comm, status, ierr)
            call mpi_recv (zs(1:ntot,1:nlev), (ntot  )*nlev, mpi_real, jpr, tag1, comm, status, ierr)
            zpc(0:ntot,1:nlev) = zpc(0:ntot,1:nlev) + zc(0:ntot,1:nlev)
            zps(1:ntot,1:nlev) = zps(1:ntot,1:nlev) + zs(1:ntot,1:nlev)
          enddo
          do jpr = ip_first+nprocsy, iprocs-1, nprocsy
            call mpi_send (zpc(0:ntot,1:nlev), (ntot+1)*nlev, mpi_real, jpr, tag2, comm, ierr)
            call mpi_send (zps(1:ntot,1:nlev), (ntot  )*nlev, mpi_real, jpr, tag2, comm, ierr)
          enddo
        endif
#endif
      endif

      do jklev = 1, nlev
      p(1:nlon,jklev) = zpc(0,jklev)
      do n = 1, ntot
      p(1:nlon,jklev) = p(1:nlon,jklev) + zps(n,jklev)*snt(1:nlon,n) + zpc(n,jklev)*cst(1:nlon,n)
      enddo
      enddo

      return
      end subroutine slofou1
!##################################################################################################################
SUBROUTINE VEGET_PARAM_TIME_SERIES(NLON,NLAT,ALON,ALAT,X0,Y0,DLON,DLAT,LONINI,LATINI,FLAG_GLOB,&
 JDAY,LAI_DEF,LAI_MAX,FVEG_DEF,VAL_MISSING)

! LAI - Leaf Area Index
! LAI_MAX - Maximum value of LAI (refer value)
! FVEG - Vegetation Cover Fraction
! FPAR - Active Photosynthetically Absorbed Radiotion Fraction

! Dataset CYCLOPES: http://postel.mediasfrance.org/en/DOWNLOAD/Biogeophysical-Products/ 

IMPLICIT NONE

!      NLON: grid point number along axis of rotated longitude;
!      NLAT: grid point number along axis of rotated latidude;
!      X0, Y0: Geographical coordinates of the centre of rotation (degree);
!      DLON: grid step along the domain axis of rotated longitude (degree);
!      DLAT: grid step along the domain axis of rotated latitude (degree);
!      LONINI: longitude of fist (south-west coner) domain point in rotated coordinates (degree);
!      LATINI: langitude of fist (south-west coner) domain point in rotated coordinates (degree);
!      ALON: geographical longitude of the grid points (degree, -180...180);
!      ALAT: geographical latitude of the grid points (degree, -90...90);
! FLAG_GLOB: if 0, then not global input grid, if 1, then input grid is global not rotated
!      JDAY: Julian day minus 1
!      VAL_MISSING: values used for not difened point.

INTEGER :: NLON, NLAT, FLAG_GLOB, JDAY
REAL :: X0, Y0, DLON, DLAT, LONINI, LATINI, LAI_MAX, VAL_MISSING
REAL, DIMENSION(NLON,NLAT) :: ALON, ALAT, LON_GEO, LAT_GEO, LAI_DEF, FVEG_DEF
REAL, PARAMETER :: DLON_READ=1./112., DLAT_READ=1./112.,&
                   LONINI_READ=0.+DLON_READ*0.5, LONFIN_READ=360.-DLON_READ*0.5,&
                   LATINI_READ=-60.+DLAT_READ*0.5, LATFIN_READ=80.-DLAT_READ*0.5
INTEGER, PARAMETER :: NX_READ=40320, NY_READ=(LATFIN_READ-LATINI_READ)/DLAT_READ
INTEGER*1, DIMENSION(NX_READ) :: VALUE_READ
REAL, DIMENSION(:,:,:), ALLOCATABLE :: LAI_READ, FVEG_READ, FPAR_READ
REAL, DIMENSION(:,:), ALLOCATABLE :: LON_READ_GEO, LAT_READ_GEO, LON_READ_MOD, LAT_READ_MOD 

INTEGER, DIMENSION(36) :: DAY_LIST=(/5, 15, 25, 36, 46, 56, 64, 74, 84, 95, 105, 115, 125, 135, 145, 156, 166, 176, 186, 196,&
 206, 217, 227, 237, 248, 258, 268, 278, 288, 298, 309, 319, 329, 339, 349, 359/)
CHARACTER(LEN=40), DIMENSION(2) :: FILEINPUT
INTEGER, DIMENSION(2) :: IUNIT=(/14, 15/)
INTEGER, DIMENSION(NLON,NLAT) :: NPOINT_LAI, NPOINT_FVEG
INTEGER, PARAMETER :: NX_PART=2
INTEGER, DIMENSION(NX_PART) :: IINI, IFIN
INTEGER :: DAYR, DAYL, IERR_OPEN, NLON_READ, NLAT_READ, ISTART, IFINISH, JSTART, JFINISH,&
 IFILE, IPAR, I, J, IR, JR, IX_PART, III, JJJ, JINI, JFIN, FLAG_PERIOD=0, IND, IMOD, JMOD
REAL :: LON_GEO_MIN, LON_GEO_MAX, LAT_GEO_MIN, LAT_GEO_MAX, LON_START_READ, LAT_START_READ, &
 WEIGH_DAYR, WEIGH_DAYL

 LAI_MAX=7.5

 PRINT *
 PRINT *,'Updating of vegetation LAI and fraction'
 PRINT *
 PRINT *,'Reading and processing data about vegetation parameters: LAI and FVEG using corrisponding datasets'

! Defintion days interval

 DO I=1,36
   IF (JDAY<DAY_LIST(I)) EXIT
 ENDDO

 IF (I==1) THEN
   DAYL=DAY_LIST(36)
 ELSE
   DAYL=DAY_LIST(I-1)
 ENDIF
 IF (I <= 36) THEN
   DAYR=DAY_LIST(I)
 ELSE
   DAYR=DAY_LIST(1)
 ENDIF

 IF (DAYR > DAYL) THEN
   WEIGH_DAYR=FLOAT(JDAY-DAYL)/FLOAT(DAYR-DAYL)
 ELSE
   IF (JDAY >= DAYL ) THEN
     WEIGH_DAYR=FLOAT(JDAY-DAYL)/FLOAT(DAYR+365-DAYL)
   ELSE
     WEIGH_DAYR=FLOAT(JDAY+365-DAYL)/FLOAT(DAYR+365-DAYL)
   ENDIF
 ENDIF
 WEIGH_DAYL=1.-WEIGH_DAYR

! Definition of parameters of input data grid

 LAT_GEO(:,:)=ALAT(:,:)
 LON_GEO(:,:)=ALON(:,:)
 DO J=1,NLAT
 DO I=1,NLON
   IF (ABS(ABS(LAT_GEO(I,J))-90.) < 1.e-10) LAT_GEO(I,J) = SIGN(89.999999, LAT_GEO(I,J))
   IF (ALON(I,J) < -180.) LON_GEO(I,J) = ALON(I,J)+360.
   IF (ALON(I,J) >  180.) LON_GEO(I,J) = ALON(I,J)-360.
 ENDDO
 ENDDO
! Dataset longitude: 0..360
 IF (FLAG_GLOB == 1.AND.ABS(LONINI-LONINI_READ) < 10.) THEN
   LON_GEO(:,:)=ALON(:,:)
 ENDIF

 LON_GEO_MIN=MINVAL(LON_GEO)
 LON_GEO_MAX=MAXVAL(LON_GEO)
 LAT_GEO_MIN=MINVAL(LAT_GEO)
 LAT_GEO_MAX=MAXVAL(LAT_GEO)

 ISTART=INT(LON_GEO_MIN*0.2)*5-10
 IFINISH=INT(LON_GEO_MAX*0.2)*5+15
 JSTART=MAX( INT(LAT_GEO_MIN*0.2)*5-10, -90)
 JFINISH=MIN( INT(LAT_GEO_MAX*0.2)*5+15, 90)

 IF (FLAG_GLOB == 0) THEN

! Control of read data domain along latitude axis

   IF (FLOAT(JSTART) >= LATINI_READ.AND.FLOAT(JFINISH) <= LATFIN_READ) THEN

! Domain extremes locate inside of latitude limits of dataset global domain

     JINI=NINT((FLOAT(JSTART)-LATINI_READ)/DLAT_READ)+1
     JFIN=NINT((FLOAT(JFINISH)-LATINI_READ)/DLAT_READ)+1

   ELSE

! Domain extremes located outside of latitude limits of dataset global domain

     IF (FLOAT(JSTART) < LATINI_READ.AND.FLOAT(JFINISH) <= LATFIN_READ) THEN

       JINI=1
       JFIN=INT((FLOAT(JFINISH)-LATINI_READ)/DLAT_READ)+1
       IF (FLOAT(JSTART) < LATINI_READ) THEN
         print *,' Caution: LAI, Vegetation fraction datasets does not cover the south part of the model domain'
         print *,' Model grid south latitude extreme',FLOAT(JSTART)
         print *,' Dataset grid south latitude extreme',LATINI_READ
       ENDIF

     ELSE IF (FLOAT(JFINISH) > LATFIN_READ.AND.FLOAT(JSTART) >= LATFIN_READ) THEN

       JINI = INT((FLOAT(JSTART)-LATINI_READ)/DLAT_READ)+1
       JFIN = NY_READ
       IF (FLOAT(JFINISH) > LATFIN_READ) THEN
         print *,' Caution: LAI, Vegetation fraction datasets does not cover the north part of the model domain'
         print *,' Model grid north latitude extreme',FLOAT(JFINISH)
         print *,' Dataset grid north latitude extrem',LATFIN_READ
       ENDIF

     ELSE

! Domain covers the whole latitude range

       JINI = 1
       JFIN = NY_READ
       IF (FLOAT(JSTART) < LATINI_READ.OR.FLOAT(JFINISH) > LATFIN_READ) THEN
         print *,' Caution: LAI, Vegetation fraction datasets does not cover the north and the south parts of the model domain'
         print *,' Model grid latitude extremes',FLOAT(JSTART),FLOAT(JFINISH)
         print *,' Dataset grid latitude extremes',LATINI_READ,LATFIN_READ
       ENDIF

     ENDIF

   ENDIF

 ELSE

! Domain covers the whole latitude range

   JINI = 1
   JFIN = NY_READ
   print *,' Caution: LAI, Vegetation fraction datasets does not cover the north and the south parts of the model domain'
   print *,' Model grid latitude extremes',LAT_GEO(1,1),LAT_GEO(1,NLAT)
   print *,' Dataset grid latitude extremes',LATINI_READ,LATFIN_READ

 ENDIF 

 LAT_START_READ = LATINI_READ+FLOAT(JINI-1)*DLAT_READ

! Control of read data domain along longitude axis, possible cut and paste operation for periodic grid

 IF (FLAG_GLOB == 0) THEN

   IF (ISTART >= INT(LONINI_READ).AND.IFINISH <= INT(LONFIN_READ)) THEN

! Domain extremes locate inside of longitude limits of dataset global domain

     IINI(1) = INT((FLOAT(ISTART)-LONINI_READ)/DLON_READ)+1
     IFIN(1) = INT((FLOAT(IFINISH)-LONINI_READ)/DLON_READ)+1

     IINI(2) = 0
     IFIN(2) = -1

     FLAG_PERIOD = 0

   ELSE

! Domain extremes locate outside of longitude limits of dataset global domain -> cut and paste

     IF (ISTART < INT(LONINI_READ).AND.IFINISH <= INT(LONFIN_READ)) THEN

       IINI(1) = INT((FLOAT(ISTART)+360.-LONINI_READ)/DLON_READ)+1
       IFIN(1) = NX_READ

       IINI(2) = 1
       IFIN(2) = INT((FLOAT(IFINISH)-LONINI_READ)/DLON_READ)+1

       FLAG_PERIOD = 1

     ELSEIF (IFINISH > INT(LONFIN_READ).AND.ISTART >= INT(LONINI_READ)) THEN

       IINI(1) = INT((FLOAT(ISTART)-LONINI_READ)/DLON_READ)+1
       IFIN(1) = NX_READ

       IINI(2) = 1
       IFIN(2) = INT((FLOAT(IFINISH)-360.-LONINI_READ)/DLON_READ)+1

       FLAG_PERIOD = 1

     ELSE

! Domain covers the whole longitude circle

       IINI(1) = 1
       IFIN(1) = NX_READ

       IINI(2) = 0
       IFIN(2) = -1

       FLAG_PERIOD = 2

     ENDIF

   ENDIF

 ELSE ! global output grid

! Domain covers the whole longitude circle

   IF (ABS(LONINI-LONINI_READ) < 10.) THEN

! Input (dataset) and output (model) global grid have the same longitude origin

     IINI(1) = 1
     IFIN(1) = NX_READ

     IINI(2) = 0
     IFIN(2) = -1

     FLAG_PERIOD = 2

   ELSE
! Input (dataset) and output (model) global grid have the longitude origin shited by 180 degees

     IINI(1) = INT((LONINI-LONINI_READ)/DLON_READ)+1
     IFIN(1) = NX_READ

     IINI(2) = 1
     IFIN(2) = IINI(1)-1

     FLAG_PERIOD = 3

   ENDIF

 ENDIF

 LON_START_READ = LONINI_READ+FLOAT(IINI(1)-1)*DLON_READ

 NLAT_READ = JFIN-JINI+1
 NLON_READ = IFIN(1)-IINI(1)+1
 IF (FLAG_PERIOD == 1 ) THEN
   NLON_READ = NLON_READ+IFIN(2)-IINI(2)+1
 ENDIF
 IF (FLAG_PERIOD == 2.OR.FLAG_PERIOD == 3) THEN
   NLON_READ = NX_READ+2
   IF (FLAG_PERIOD == 2) LON_START_READ = LONINI_READ-DLON_READ
   IF (FLAG_PERIOD == 3) LON_START_READ = LON_START_READ-DLON_READ
 ENDIF

 IF (LON_START_READ > 180.) LON_START_READ=LON_START_READ-360.

 ALLOCATE(LAI_READ(NLON_READ,NLAT_READ,3))
 ALLOCATE(FVEG_READ(NLON_READ,NLAT_READ,3))
!!! ALLOCATE(FPAR_READ(NLON_READ,NLAT_READ,3))
 ALLOCATE(LON_READ_GEO(NLON_READ,NLAT_READ))
 ALLOCATE(LAT_READ_GEO(NLON_READ,NLAT_READ))
 ALLOCATE(LON_READ_MOD(NLON_READ,NLAT_READ))
 ALLOCATE(LAT_READ_MOD(NLON_READ,NLAT_READ))

 LAI_READ(:,:,:)=VAL_MISSING
 FVEG_READ(:,:,:)=VAL_MISSING
!!! FPAR_READ(:,:,:)=VAL_MISSING
 LON_READ_GEO(:,:)=VAL_MISSING
 LAT_READ_GEO(:,:)=VAL_MISSING
 LON_READ_MOD(:,:)=VAL_MISSING
 LAT_READ_MOD(:,:)=VAL_MISSING

! Reading of datasets

 DO IPAR=1,2 ! 1 - LAI, 2 - FVEG, 3 - FPAR

   IF (IPAR==1) THEN
     WRITE (FILEINPUT(1),'(A,I3.3,A)') "lai_global_latlon_1km_day",DAYL,".bin"
     WRITE (FILEINPUT(2),'(A,I3.3,A)') "lai_global_latlon_1km_day",DAYR,".bin"
   ENDIF
   IF (IPAR==2) THEN
     WRITE (FILEINPUT(1),'(A,I3.3,A)') "fveg_global_latlon_1km_day",DAYL,".bin"
     WRITE (FILEINPUT(2),'(A,I3.3,A)') "fveg_global_latlon_1km_day",DAYR,".bin"
   ENDIF
!!!   IF (IPAR==3) THEN
!!!     WRITE (FILEINPUT(1),'(A,I3.3,A)') "fpar_global_latlon_1km_day",DAYL,".bin"
!!!     WRITE (FILEINPUT(2),'(A,I3.3,A)') "fpar_global_latlon_1km_day",DAYR,".bin"
!!!   ENDIF

   DO IFILE=1,2

     PRINT *
     PRINT *,'Reading the global vegetation parameter dataset from ',FILEINPUT(IFILE)

     OPEN(IUNIT(IFILE),FILE=FILEINPUT(IFILE),STATUS='OLD',FORM='UNFORMATTED',ACCESS='DIRECT',RECL=NX_READ,IOSTAT=IERR_OPEN)

     IF (IERR_OPEN == 0) THEN

       J=0
       DO JJJ=JINI,JFIN
         J=J+1
      
         READ (IUNIT(IFILE), REC=JJJ) VALUE_READ
  
         IF (IPAR == 1) THEN ! LAI
  
           I=0
           IF (FLAG_PERIOD == 2.OR.FLAG_PERIOD == 3) I=1
           DO IX_PART=1,NX_PART
             DO III=IINI(IX_PART),IFIN(IX_PART)
               I=I+1
               IF (VALUE_READ(III) >= 0) THEN
                 LAI_READ(I,J,IFILE)=FLOAT(VALUE_READ(III))
               ELSE
                 LAI_READ(I,J,IFILE)=FLOAT(VALUE_READ(III))+256.
               ENDIF
               LON_READ_GEO(I,J)=LON_START_READ+FLOAT(I-1)*DLON_READ
               IF (FLAG_GLOB == 0.AND.LON_READ_GEO(I,J)>180.) LON_READ_GEO(I,J)=LON_READ_GEO(I,J)-360. !???
             ENDDO
           ENDDO ! IX_PART
           IF (FLAG_PERIOD == 2.OR.FLAG_PERIOD == 3) THEN
             LAI_READ(1,J,IFILE)=LAI_READ(NLON_READ-1,J,IFILE)
             LAI_READ(NLON_READ,J,IFILE)=LAI_READ(2,J,IFILE)
             LON_READ_GEO(1,J)=LON_START_READ
             LON_READ_GEO(NLON_READ,J)=LON_START_READ+FLOAT(NLON_READ-1)*DLON_READ
             IF (FLAG_GLOB == 0) THEN
               IF (LON_READ_GEO(1,J)>180.) LON_READ_GEO(1,J)=LON_READ_GEO(1,J)-360. !???
               IF (LON_READ_GEO(NLON_READ,J)>180.) LON_READ_GEO(NLON_READ,J)=LON_READ_GEO(NLON_READ,J)-360. !???
             ENDIF
           ENDIF
        
           LAT_READ_GEO(:,J)=LAT_START_READ+FLOAT(J-1)*DLAT_READ
  
         ELSE IF (IPAR == 2) THEN ! FVEG
  
           I=0
           IF (FLAG_PERIOD == 2.OR.FLAG_PERIOD == 3) I=1
           DO IX_PART=1,NX_PART
             DO III=IINI(IX_PART),IFIN(IX_PART)
               I=I+1
               IF (VALUE_READ(III) >= 0) THEN
                 FVEG_READ(I,J,IFILE)=FLOAT(VALUE_READ(III))
               ELSE
                 FVEG_READ(I,J,IFILE)=FLOAT(VALUE_READ(III))+256.
               ENDIF
             ENDDO
           ENDDO ! IX_PART
           IF (FLAG_PERIOD == 2.OR.FLAG_PERIOD == 3) THEN
             FVEG_READ(1,J,IFILE)=FVEG_READ(NLON_READ-1,J,IFILE)
             FVEG_READ(NLON_READ,J,IFILE)=FVEG_READ(2,J,IFILE)
           ENDIF
        
!!!         ELSE IF (IPAR == 3) THEN ! FPAR
!!!  
!!!           I=0
!!!           IF (FLAG_PERIOD == 2.OR.FLAG_PERIOD == 3) I=1
!!!           DO IX_PART=1,NX_PART
!!!             DO III=IINI(IX_PART),IFIN(IX_PART)
!!!               I=I+1
!!!               IF (VALUE_READ(III) >= 0) THEN
!!!                 FPAR_READ(I,J,IFILE)=FLOAT(VALUE_READ(III))
!!!               ELSE
!!!                 FPAR_READ(I,J,IFILE)=FLOAT(VALUE_READ(III))+256.
!!!               ENDIF
!!!             ENDDO
!!!           ENDDO ! IX_PART
!!!           IF (FLAG_PERIOD == 2.OR.FLAG_PERIOD == 3) THEN
!!!             FPAR_READ(1,J,IFILE)=FPAR_READ(NLON_READ-1,J,IFILE)
!!!             FPAR_READ(NLON_READ,J,IFILE)=FPAR_READ(2,J,IFILE)
!!!           ENDIF
        
         ENDIF
      
       ENDDO ! JJJ
      
!       IF (IPAR == 1) THEN
!         IF (FLAG_PERIOD == 2 ) THEN
!           LON_READ_GEO(1,:) = LON_START_READ
!           LON_READ_GEO(NLON_READ,:) = LON_START_READ+FLOAT(NLON_READ-1)*DLON_READ
!         ENDIF
!       ENDIF
  
       CLOSE(IUNIT(IFILE))

     ELSE ! Not found file

       PRINT *
       PRINT *,'Not found input file ',FILEINPUT(IFILE)
       IF (IPAR == 1) PRINT *,'Vegetation LAI not updated'
       IF (IPAR == 2) PRINT *,'Vegetation fraction not updated'
       PRINT *
       RETURN

     ENDIF

   ENDDO ! IFILE

 ENDDO ! IPAR

 DO IFILE=1,2
   DO J=1,NLAT_READ
   DO I=1,NLON_READ
     IF (LAI_READ(I,J,IFILE) < 255.) THEN
       LAI_READ(I,J,IFILE)=LAI_READ(I,J,IFILE)/30.
     ELSE
       LAI_READ(I,J,IFILE)=VAL_MISSING
     ENDIF
     IF (FVEG_READ(I,J,IFILE) < 255.) THEN
       FVEG_READ(I,J,IFILE)=FVEG_READ(I,J,IFILE)/250.
     ELSE
       FVEG_READ(I,J,IFILE)=VAL_MISSING
     ENDIF
!!!     IF (FPAR_READ(I,J,IFILE) < 255.) THEN
!!!       FPAR_READ(I,J,IFILE)=FPAR_READ(I,J,IFILE)/250.
!!!     ELSE
!!!       FPAR_READ(I,J,IFILE)=VAL_MISSING
!!!     ENDIF
   ENDDO
   ENDDO
 ENDDO

! Time interpolation

 PRINT *
 PRINT *,'Time interpolation'

 IND=3
 DO J=1,NLAT_READ
 DO I=1,NLON_READ
   IF (LAI_READ(I,J,1)/=VAL_MISSING.AND.LAI_READ(I,J,2)/=VAL_MISSING) THEN
     LAI_READ(I,J,IND)=LAI_READ(I,J,1)*WEIGH_DAYL+LAI_READ(I,J,2)*WEIGH_DAYR
   ELSEIF (LAI_READ(I,J,1)/=VAL_MISSING) THEN
     LAI_READ(I,J,IND)=LAI_READ(I,J,1)
   ELSEIF (LAI_READ(I,J,2)/=VAL_MISSING) THEN
     LAI_READ(I,J,IND)=LAI_READ(I,J,2)
   ENDIF
   IF (FVEG_READ(I,J,1)/=VAL_MISSING.AND.FVEG_READ(I,J,2)/=VAL_MISSING) THEN
     FVEG_READ(I,J,IND)=FVEG_READ(I,J,1)*WEIGH_DAYL+FVEG_READ(I,J,2)*WEIGH_DAYR
   ELSEIF (FVEG_READ(I,J,1)/=VAL_MISSING) THEN
     FVEG_READ(I,J,IND)=FVEG_READ(I,J,1)
   ELSEIF (FVEG_READ(I,J,2)/=VAL_MISSING) THEN
     FVEG_READ(I,J,IND)=FVEG_READ(I,J,2)
   ENDIF
 ENDDO
 ENDDO

! Space interpolation = average value for read pixel valuees around model grid point

 PRINT *,'Space interpolation'

 CALL ANTI_ROT_GRID(X0,Y0,LON_READ_GEO,LAT_READ_GEO,LON_READ_MOD,LAT_READ_MOD,NLON_READ,NLAT_READ)

 DO J=1,NLAT_READ
 DO I=1,NLON_READ
  IF (LON_READ_GEO(I,J)==VAL_MISSING) LON_READ_MOD(I,J)=VAL_MISSING
  IF (LAT_READ_GEO(I,J)==VAL_MISSING) LAT_READ_MOD(I,J)=VAL_MISSING
 ENDDO
 ENDDO

 NPOINT_LAI(:,:)=0
 NPOINT_FVEG(:,:)=0
 DO J=1,NLAT_READ
 DO I=1,NLON_READ
   IF (LON_READ_MOD(I,J)/=VAL_MISSING.AND.LAT_READ_MOD(I,J)/=VAL_MISSING) THEN
     IMOD=NINT((LON_READ_MOD(I,J)-LONINI)/DLON)+1
     JMOD=NINT((LAT_READ_MOD(I,J)-LATINI)/DLAT)+1
     IF (IMOD>=1.AND.IMOD<=NLON.AND.JMOD>=1.AND.JMOD<=NLAT) THEN
       IF (LAI_READ(I,J,IND)/=VAL_MISSING) THEN
         IF (NPOINT_LAI(IMOD,JMOD)>0) THEN
           LAI_DEF(IMOD,JMOD)=LAI_DEF(IMOD,JMOD)+LAI_READ(I,J,IND)
         ELSE
           LAI_DEF(IMOD,JMOD)=LAI_READ(I,J,IND)
         ENDIF
         NPOINT_LAI(IMOD,JMOD)=NPOINT_LAI(IMOD,JMOD)+1
       ENDIF
       IF (FVEG_READ(I,J,IND)/=VAL_MISSING) THEN
         IF (NPOINT_FVEG(IMOD,JMOD)>0) THEN
           FVEG_DEF(IMOD,JMOD)=FVEG_DEF(IMOD,JMOD)+FVEG_READ(I,J,IND)
         ELSE
           FVEG_DEF(IMOD,JMOD)=FVEG_READ(I,J,IND)
         ENDIF
         NPOINT_FVEG(IMOD,JMOD)=NPOINT_FVEG(IMOD,JMOD)+1
       ENDIF
     ENDIF
   ENDIF
 ENDDO 
 ENDDO 

 DO J=1,NLAT
 DO I=1,NLON
   IF (NPOINT_LAI(I,J)>0) THEN
     LAI_DEF(I,J)=LAI_DEF(I,J)/FLOAT(NPOINT_LAI(I,J))
     LAI_DEF(I,J)=MAX( MIN( LAI_DEF(I,J), LAI_MAX), 0.)
   ELSE
     LAI_DEF(I,J)=0.
   ENDIF
   IF (NPOINT_FVEG(I,J)>0) THEN
     FVEG_DEF(I,J)=FVEG_DEF(I,J)/FLOAT(NPOINT_LAI(I,J))
     FVEG_DEF(I,J)=MAX( MIN( FVEG_DEF(I,J), 1.), 0.)
   ELSE
     FVEG_DEF(I,J)=0.
   ENDIF
 ENDDO
 ENDDO

END SUBROUTINE VEGET_PARAM_TIME_SERIES
!##################################################################################################################
subroutine anti_rot_grid(x0,y0,alon,alat,xr,yr,nlon,nlat)

! Calculation of the rotated coordinates (put in XR, YR) with rotation centre X0, Y0
! (given in geographical coordinates) for input grid points given in
! geographical coordinates (defined in ALON, ALAT)
! Some computations require double precision

implicit none

real :: x0, y0
integer :: nlon, nlat, jlat, jlon
real, dimension(nlon,nlat) :: alat, alon, xr, yr
real*8 :: zfac, zx0, zy0, zx, zy, zlon, zlat, zz, pi

if (abs(x0)>0.01.or.abs(y0)>0.01) then

! Case of rotated grid

  pi  = dabs(dacos(-1.d0))
  zfac = pi/180.d0
  zx0  = dble(x0)*zfac
  zy0  = dble(y0)*zfac
  do jlat=1,nlat
  do jlon=1,nlon
    if (jlon==1.and.jlat==1) call sleep (0) ! fix for a problem of ifort 14 (but needs compil. with -O2 or less)
    zx = dble(alon(jlon,jlat))*zfac
    if (zx-zx0 > pi) zx = zx - 2.d0*pi
    if (zx-zx0 < -pi) zx = zx + 2.d0*pi
    zy = dble(alat(jlon,jlat))*zfac
    zlat = dasin( -dcos(zy)*dsin(zy0)*dcos(zx-zx0)+dsin(zy)*dcos(zy0) )
    zz = dsin(zy)-dcos(zy0)*dsin(zlat)
    zz = zz/(dsin(zy0)*dcos(zlat))
    if (zz < -1.d0.and.zz > -1.00001d0) zz = -1.d0
    if (zz > 1.d0.and.zz < 1.00001d0) zz = 1.d0
    if (zx < zx0) then
      zlon = -dacos(zz)
    else
      zlon = dacos(zz)
    endif
    zx = zlon/zfac
    zy = zlat/zfac
    xr(jlon,jlat) = sngl(zx)
    yr(jlon,jlat) = sngl(zy)
  enddo
  enddo

else

! Case of non rotated grid

  xr(:,:) = alon(:,:)
  yr(:,:) = alat(:,:)

endif

return
end subroutine anti_rot_grid
!##################################################################################################################
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
!##################################################################################################################
subroutine rdrf(istart)
! Reads restart file

use mod_model
#ifdef rad_ecmwf
use yomrip, only : reqtim, rdecli
#endif
implicit none

integer :: iunit=17, istart, ierr_open, jlon, jlat, jklev, comm, error, jpr, tag=0
real :: zdtstep_old, ztime
integer, dimension(gnlon,gnlat)      :: igfield
character(len=30) :: file_in="bolam.rf"
integer, dimension(1) :: i1

#ifdef mpi
      include 'mpif.h'
      integer :: status(mpi_status_size)

      comm = mpi_comm_world
#endif

 if(myid.eq.0) then

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
   read (iunit) reqtim0, rdecli0
#ifdef rad_ecmwf
   reqtim = dble(reqtim0)
   rdecli = dble(rdecli0)
#endif

   if (istart /= 0) close (iunit)
    
 endif ! myid

#ifdef mpi
     call mpi_bcast(nstep0,        1, mpi_integer,  0, comm, error)   
     call mpi_bcast(nfdr0(1:50),  50, mpi_integer,  0, comm, error)   
     call mpi_bcast(pdr0(1:200), 200, mpi_real,     0, comm, error)   
     call mpi_bcast(nfdr(1:50),   50, mpi_integer,  0, comm, error)   
     call mpi_bcast(pdr(1:200),  200, mpi_real,     0, comm, error)   
     call mpi_bcast(reqtim0,       1, mpi_real,     0, comm, error)   
     call mpi_bcast(rdecli0,       1, mpi_real,     0, comm, error)   
#ifdef rad_ecmwf
     call mpi_bcast(reqtim,        1, mpi_double,   0, comm, error)   
     call mpi_bcast(rdecli,        1, mpi_double,   0, comm, error)   
#endif
#endif

 if (istart /= 0) return

! Basic variables

 if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, phig)

 if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, fmask)

 if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, ps)

 do jklev=1,nlev
   if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, u(1,1,jklev))
 enddo

 do jklev=1,nlev
   if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, v(1,1,jklev))
 enddo

 do jklev=1,nlev
   if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, t(1,1,jklev))
 enddo

 do jklev=1,nlev
   if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, q(1,1,jklev))
 enddo

 do jklev=1,nlev
   if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, qcw(1,1,jklev))
 enddo

 do jklev=1,nlev
   if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, qci(1,1,jklev))
 enddo

 do jklev=1,nlev
   if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, qrain(1,1,jklev))
 enddo

 do jklev=1,nlev
   if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, qsnow(1,1,jklev))
 enddo

 do jklev=1,nlev
   if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, qhail(1,1,jklev))
 enddo

 if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, tskin)

 do jklev = 1, nlevg
   if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, tg(1,1,jklev))
 enddo

 if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, qskin)

 do jklev = 1, nlevg
   if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, qg(1,1,jklev))
 enddo

 if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, cloudt)

 if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, totpre)

 if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, conpre)

 if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, snfall)

 if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, snow)

 if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, cswfl)

 if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, clwfl)

 if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, chflux)

 if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, cqflux)

 if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, t2min)

 if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, t2max)

 if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, albedo)

 if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, emisg1)

 if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, emisg2)

 if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, rgm)

 if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, rgmd)

 if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, rgq)

 if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, runoff)

 if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, runoff_tot)

 if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, fice)

 if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, iceth)

 if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, alsn)

 if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, cwvflux)

 if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, fl_heat_soil_bottom)

 if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, fl_water_soil_bottom)

!  Physical parameters of soil scheme

 if(myid.eq.0) call rrec2_int (iunit, gnlon, gnlat, igfield)
 call disperse_int (igfield, mask_soil)

 if(myid.eq.0) call rrec2_int (iunit, gnlon, gnlat, igfield)
 call disperse_int (igfield, ind_lev_soil_h_bottom)

 if(myid.eq.0) call rrec2_int (iunit, gnlon, gnlat, igfield)
 call disperse_int (igfield, ind_lev_soil_w_bottom)

 if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, kturb_surf_m)

 if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, kturb_surf_h)

 if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, kturb_surf_q)

 do jklev = 1, 20
   if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, kturb_surf_h_mem(1,1,jklev))
 enddo

 do jklev = 1, 20
   if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, kturb_surf_q_mem(1,1,jklev))
 enddo

 if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, flh_specif)

 if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, fl_wv)

 if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, fl_heat_soil_bottom)

 if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, fl_water_soil_bottom)

 if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, fl_runoff)

 if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, fl_runoff_tot)

 if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, cfl_heat_soil_bottom)

 if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, cfl_water_soil_bottom)

 do jklev = 1, nlev_snow
   if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, snow_lev(1,1,jklev))
 enddo

 do jklev = 1, nlev_snow
   if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, snow_t(1,1,jklev))
 enddo

 do jklev = 1, nlev_snow
   if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, snow_age(1,1,jklev))
 enddo

 do jklev = 1, nlev_snow
   if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, snow_melt_age(1,1,jklev))
 enddo

 do jklev = 1, nlev_snow
   if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, snow_dens(1,1,jklev))
 enddo

 do jklev = 1, nlevg
   if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, psi_soil(1,1,jklev))
 enddo

 do jklev = 1, nlevg
   if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, qgice(1,1,jklev))
 enddo

 if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, fveg)

 if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, lai)

! Other variables

 if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, t2)

 do jklev=1,nlev
   if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, geldt(1,1,jklev))
 enddo

 do jklev=1,nlev
   if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, corrdt(1,1,jklev))
 enddo

 do jklev=1,nlev
   if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, phi(1,1,jklev))
 enddo

 do jklev=1,nlevp1
   if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, phih(1,1,jklev))
 enddo

 do jklev=1,nlevp1
   if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, sigdot(1,1,jklev))
 enddo

 do jklev=1,nlevp1
   if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, tke(1,1,jklev))
 enddo

 do jklev=1,nlevp1
   if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, lml(1,1,jklev))
 enddo

 do jklev=1,nlev
   if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, trd_a(1,1,jklev))
 enddo

 do jklev=1,nlev
   if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, trd_c(1,1,jklev))
 enddo

 if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, raicon)

 if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, snocon)

 if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, rainls)

 if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, snowls)

 if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, qprec)

 if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, qprecc)

 if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, qsnfall)

 if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, hflux)

 if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, qflux)

 if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, frvis)

 if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, frirr)

 if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, gelvis)

 if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, gelirr)

 if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, corvis)

 if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, corirr)

 if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, solar)

 do jklev=1,nlev
   if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
   call disperse (gfield, fcloud(1,1,jklev))
 enddo

 if(myid.eq.0) call rrec2 (iunit, gnlon, gnlat, gfield)
 call disperse (gfield, fsnow)

 if(myid.eq.0) then
   close (iunit)
   print *
   print*,trim(file_in),' read, nstep0= ',nstep0
 endif

 do jlon=1,nlon
 do jlat=1,nlat
   if (fice(jlon,jlat).ge.0.8)  suolo(jlon,jlat,14) = 0.
 enddo
 enddo

 roscd       =1.e-2
!!! sigdot(:,:,1     )=0.
!!! sigdot(:,:,nlevp1)=0.
!!! phih  (:,:,nlevp1)=phig(:,:)
!!! phih  (:,:,     1)=1.e6
!!! rgmd(:,:) = rgm(:,:)

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

use mod_model
implicit none

integer :: iunit=27, kstep, jlon, jlat, jklev, i1, i2, i3, &
 ilon1=1, ilon2=gnlon, jlat1=1, jlat2=gnlat, flag_lon=0
integer, dimension(gnlon,gnlat)      :: igfield
character(len=30) :: file_out

 if(myid.eq.0) then

   irf=irf+1
   i1=irf/10
   i2=i1*10
   i3=irf-i2
   write (file_out,'(A,I2.2,A)') 'bolam_',i3,'.rf'

   open (iunit, file=trim(file_out), form='unformatted', status='unknown')

   write (iunit) kstep+1
   write (iunit) nfdr0
   write (iunit) pdr0
   write (iunit) nfdr
   write (iunit) pdr
   write (iunit) reqtim0, rdecli0

 endif

! Basic variables

 call collect (phig, gfield)
 if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (fmask, gfield)
 if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (ps, gfield)
 if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 do jklev=1,nlev
   call collect (u(1,1,jklev), gfield)
   if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
 enddo

 do jklev=1,nlev
   call collect (v(1,1,jklev), gfield)
   if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
 enddo

 do jklev=1,nlev
   call collect (t(1,1,jklev), gfield)
   if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
 enddo

 do jklev=1,nlev
   call collect (q(1,1,jklev), gfield)
   if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
 enddo

 do jklev=1,nlev
   call collect (qcw(1,1,jklev), gfield)
   if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
 enddo

 do jklev=1,nlev
   call collect (qci(1,1,jklev), gfield)
   if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
 enddo

 do jklev=1,nlev
   call collect (qrain(1,1,jklev), gfield)
   if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
 enddo

 do jklev=1,nlev
   call collect (qsnow(1,1,jklev), gfield)
   if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
 enddo

 do jklev=1,nlev
   call collect (qhail(1,1,jklev), gfield)
   if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
 enddo

 call collect (tskin, gfield)
 if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 do jklev = 1, nlevg
   call collect (tg(1,1,jklev), gfield)
   if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
 enddo

 call collect (qskin, gfield)
 if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 do jklev = 1, nlevg
   call collect (qg(1,1,jklev), gfield)
   if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
 enddo

 call collect (cloudt, gfield)
 if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (totpre, gfield)
 if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (conpre, gfield)
 if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (snfall, gfield)
 if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (snow, gfield)
 if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (cswfl, gfield)
 if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (clwfl, gfield)
 if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (chflux, gfield)
 if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (cqflux, gfield)
 if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (t2min, gfield)
 if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (t2max, gfield)
 if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (albedo, gfield)
 if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (emisg1, gfield)
 if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (emisg2, gfield)
 if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (rgm, gfield)
 if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (rgmd, gfield)
 if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (rgq, gfield)
 if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (runoff, gfield)
 if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (runoff_tot, gfield)
 if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (fice, gfield)
 if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (iceth, gfield)
 if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (alsn, gfield)
 if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (cwvflux, gfield)
 if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (fl_heat_soil_bottom, gfield)
 if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (fl_water_soil_bottom, gfield)
 if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

!  Physical parameters of soil scheme

 call collect_int (mask_soil, igfield)
 if(myid.eq.0) call wrec2_int (iunit, gnlon, gnlat, igfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect_int (ind_lev_soil_h_bottom, igfield)
 if(myid.eq.0) call wrec2_int (iunit, gnlon, gnlat, igfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect_int (ind_lev_soil_w_bottom, igfield)
 if(myid.eq.0) call wrec2_int (iunit, gnlon, gnlat, igfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (kturb_surf_m, gfield)
 if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (kturb_surf_h, gfield)
 if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (kturb_surf_q, gfield)
 if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 do jklev = 1, 20
   call collect (kturb_surf_h_mem(1,1,jklev), gfield)
   if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
 enddo

 do jklev = 1, 20
   call collect (kturb_surf_q_mem(1,1,jklev), gfield)
   if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
 enddo

 call collect (flh_specif, gfield)
 if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (fl_wv, gfield)
 if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (fl_heat_soil_bottom, gfield)
 if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (fl_water_soil_bottom, gfield)
 if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (fl_runoff, gfield)
 if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (fl_runoff_tot, gfield)
 if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (cfl_heat_soil_bottom, gfield)
 if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (cfl_water_soil_bottom, gfield)
 if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 do jklev = 1, nlev_snow
   call collect (snow_lev(1,1,jklev), gfield)
   if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
 enddo

 do jklev = 1, nlev_snow
   call collect (snow_t(1,1,jklev), gfield)
   if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
 enddo

 do jklev = 1, nlev_snow
   call collect (snow_age(1,1,jklev), gfield)
   if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
 enddo

 do jklev = 1, nlev_snow
   call collect (snow_melt_age(1,1,jklev), gfield)
   if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
 enddo

 do jklev = 1, nlev_snow
   call collect (snow_dens(1,1,jklev), gfield)
   if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
 enddo

 do jklev = 1, nlevg
   call collect (psi_soil(1,1,jklev), gfield)
   if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
 enddo

 do jklev = 1, nlevg
   call collect (qgice(1,1,jklev), gfield)
   if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
 enddo

 call collect (fveg, gfield)
 if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (lai, gfield)
 if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

! Other variables

 call collect (t2, gfield)
 if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 do jklev=1,nlev
   call collect (geldt(1,1,jklev), gfield)
   if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
 enddo

 do jklev=1,nlev
   call collect (corrdt(1,1,jklev), gfield)
   if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
 enddo

 do jklev=1,nlev
   call collect (phi(1,1,jklev), gfield)
   if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
 enddo

 do jklev=1,nlevp1
   call collect (phih(1,1,jklev), gfield)
   if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
 enddo

 do jklev=1,nlevp1
   call collect (sigdot(1,1,jklev), gfield)
   if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
 enddo

 do jklev=1,nlevp1
   call collect (tke(1,1,jklev), gfield)
   if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
 enddo

 do jklev=1,nlevp1
   call collect (lml(1,1,jklev), gfield)
   if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
 enddo

 do jklev=1,nlev
   call collect (trd_a(1,1,jklev), gfield)
   if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
 enddo

 do jklev=1,nlev
   call collect (trd_c(1,1,jklev), gfield)
   if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
 enddo

 call collect (raicon, gfield)
 if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (snocon, gfield)
 if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (rainls, gfield)
 if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (snowls, gfield)
 if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (qprec, gfield)
 if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (qprecc, gfield)
 if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (qsnfall, gfield)
 if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (hflux, gfield)
 if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (qflux, gfield)
 if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (frvis, gfield)
 if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (frirr, gfield)
 if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (gelvis, gfield)
 if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (gelirr, gfield)
 if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (corvis, gfield)
 if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (corirr, gfield)
 if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 call collect (solar, gfield)
 if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 do jklev=1,nlev
   call collect (fcloud(1,1,jklev), gfield)
   if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)
 enddo

 call collect (fsnow, gfield)
 if(myid.eq.0) call wrec2 (iunit, gnlon, gnlat, gfield, ilon1, ilon2, jlat1, jlat2, flag_lon)

 if(myid.eq.0) then
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
CHARACTER(LEN=30) :: FILE_IN="bolam_pochva.rf", FILE_OUT

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
   WRITE (FILE_OUT,'(A,I2.2,A)') 'bolam_pochva_',I3,'.rf'
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
