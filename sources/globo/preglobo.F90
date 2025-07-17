! Preglobo
! Introdure lbot !!!!
! v. preglobo_s2s.F90


! Input data GFS-NOAA or IFS-ECMWF in grib2 format
!
! Ultima modifica: 03/04/2025
!
! Apr. 2025: New definition of vertical coordinate (vert_coord_c, aka, bika)
!
! Sett. 2024: 
! Cambiamento in subroutine physiographic_param (def_soil.F90) e scrittura
! di model_param_constant.bin - piu' ordine, eliminazione dei passaggi 
! inutile. 
! Nuovo formato mhf (scrittura e lettura di mhf di Bolam record non 1d ma 2d)
!
! Feb. 2019:
! Nuovo schema del suolo (pochva), i nuovi dataset fisiografici (nuovo def_soil.F90),
! parametri del suolo sono diversi su vari livelli,
! manto nevoso muti-stato, il numero dei livelli di neve e' 11,
! la coordinate verticale nello manto nevoso e' kg/m/m,
! quando su un livello la neve non c'e', tutti i paramtri di neve sono -9999.;
! nuovo formato mhf: mfs_atm, mhf_soil e il file con i campi statici (model_param_constant.bin).
!##################################################################################################
! ********  list of required variables in input data *********
!
!            %%%%     gfs-noaa/ncep    %%%%
! temperature, u-component and v-component of wind; geopotential at all
! isobaric levels in the atmosphere; relative humidity at all isobaric levels
! below 100 hpa; (optionally) cloud total liquid + ice water content at the
! same levels);
! soil temperature and soil specific volumetric water content at 4 levels
! below ground surface;
! at ground or sea surface: sea-land fraction, geopotential, temperature,
! equivalent water content in the snow cover, sea ice fraction.
!
!            %%%%      ifs-ecmwf       %%%%
! case of data on isobaric atmospheric levels:
! temperature; u-component and v-component of wind; geopotential;
! specific humidity; (optionally) cloud liquid and cloud ice water content,
! at all isobaric levels, geopotential of the surface "model level",
! consistent with ps

! case of data on hybrid model atmospheric levels:
! temperature; u-component and v-component of wind; specific humidity;
! (optionally) cloud liquid and cloud ice water content), at all selected
! hybrid levels; natural logarithm of surface pressure at model levels;
! soil temperature and soil specific volumetric water content at 4 levels
! below ground surface, or at the top and bottom soil levels only;
! at ground or sea surface: sea-land fraction, geopotential, soil type;
! temperature, equivalent water content in the snow cover, sea ice fraction,
! sea ice temperature level 1.
!##################################################################################################
   module param
   integer :: nlon, nlat, nlev, nlevp1, nlevg, npolcap, ntot, ivert_coord_type=1
   integer, parameter :: nlevg0=20, nlevsnow=11

! alfa, pzer: parameters defining the hybrid vertical coordinate
! soil_lev: depth of soil levels (m)

   real :: dlon, dlat, alfa, pzer, vert_coord_c=0.
   real, dimension(nlevg0) :: soil_lev0
   real, dimension(:), allocatable :: soil_lev

! ------  input data grid parameters ------

! iana, jana, dxa, dya: dimensions and resol. of the input ("analysis") grid.
! nlev_atm_inp_max: max. possible no. of atmospheric levels (used to allocate arrays only).
! nlevp: actual no. of atmospheric levels (read from input data).
! nlev_soil_inp_max : max. possible number of soil levels in input data (used to allocate arrays only).
! nlev_soil_inp : actual number of soil levels in input data.
! x1a, y1a: coord. in deg. of the sw corner of the input data grid grid.

  integer :: iana, iana0, jana, nlev_atm_inp_max, nlev_atm_inp, nlev_atm_inpp1, nlev_soil_inp, nlev_soil_inp_max, &
             nauxp, flag_cut_paste
  real :: dxa, dya, x1a, y1a

 logical :: surf_elaborate=.true.

! nst is soil types number, nvt is vegetation type number

 integer :: nst=31, nvt=22

! flags and other output parameters of the mhf

 integer, dimension(50) :: nfdr
 real, dimension(200)   :: pdr

 integer, dimension(5) :: idate0
 integer, dimension(4) :: iperiod_inp
 integer, dimension(3) :: iperiod

! p: pressure of full levels; pl: logarithm of pressure; phl: pressure of half-levels

 real, dimension(:,:,:), allocatable :: p, pl, phl

! log. of press. on auxiliary levels and model hybrid levels

 real, dimension(:), allocatable :: plaux, s, sigint, plsig, hxt, hxv, aka, bika, akah, bikah

! lapse rate

 real :: gammac, gamma_inp, gamma0=6.5e-3
 real, dimension(:,:), allocatable :: gamma

! no. of days of each months

 integer, dimension(12) :: imon=(/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)

! physical constants

 real :: pi, eps, ep, rcp

 real, parameter :: yliv=2834170.5, yliw=333560.5, ylwv=yliv-yliw, t0=273.15, e0= 611., g0=9.807
 real, parameter :: rd=287.05, rv=461.51, cp=1004.6, cpv=1869.46, cw=4186.8, ci=cw/2., ra=6371.e+3
 real, parameter :: ccw1=(cpv-cw)/rv, ccw2=ylwv/t0/rv-ccw1, cci1=(cpv-ci)/rv, cci2=yliv/t0/rv-cci1

 namelist/param_preglobo/ nlon, nlat, nlev, npolcap, &
                          alfa, pzer, &
                          nlev_atm_inp_max, nlev_soil_inp_max, &
                          ivert_coord_type, &
                          soil_lev0

 real :: val_missing=-9999.

end module param
!##################################################################################################
module gribana

! Model used for input data: 'ifs' (IFS, ECMWF) or 'gfs' (GFS, NOAA, NCEP, USA)

 character(len=3) :: model_input

! input fields and related parameters:

 integer :: npar3d=15, npar2d=50, npar3d_soil=6, level_type=0

 real, dimension(:,:,:,:), allocatable :: field3d, field3d_soil
 real, dimension(:,:,:), allocatable   :: field2d

 real, dimension(:,:,:), allocatable :: ge, te, ue, ve, qe, qce, qcwe, qcie
 real, dimension(:,:), allocatable   :: phige, phige2d, psel, snowe, fmaske, tskine, &
                                        soile, pmsle, cctote, u10e, v10e, t2e, td2e, &
                                        ficee, ficeee, ticee, ticeee
 real, dimension(:,:,:), allocatable :: tge, qge

 real, dimension(:,:), allocatable :: alon_inp, alat_inp

 integer, dimension(:,:), allocatable :: mask_frame

! for 3d input data fields at model hybrid or sigma levels:
! ak and bk vectors are used to define pressure on ifs model levels

 real, dimension(:), allocatable :: ak, bk

! matrices with coordinates of the input grid

 real, dimension(:), allocatable :: xe, ye, lev_list, lev_list_soil

! flag indicating presence (1) or absence (0) of cloud water/ice in input data

 integer :: iflag_cloude=0

 integer :: ifl_xf

! for conversion of volumetric water soil content into relative soil content:
! minimum and maximum soil water content as a function of soil type in input data

 real, dimension(:), allocatable   :: qgmine, qgmaxe
 real, dimension(:,:), allocatable :: qgmine2d, qgmaxe2d

end module gribana
!##################################################################################################
module mod_fields

! output fields at surface and hybrid levels

 real, dimension(:,:,:), allocatable :: u, v, t, q, qc, press
 real, dimension(:,:), allocatable   :: phig, fmask, tskin, qvsurf, ps, cloud, totpre, &
                                        conpre, snfall, snow, fice, rgm, rgq, tice, iceth, fmaski, &
 tgsurf, qgsurf, water_table_depth, tg_bottom, qg_rel_bottom, qg_rel_surf_approx, &
 albedo, emismap1, emismap2, snow_albedo, &
 soil_albedo_dry, soil_albedo_wet, soil_emiss1_dry, soil_emiss1_wet, soil_emiss2_dry, soil_emiss2_wet, &
 veg_lai, veg_frac, veg_root_depth, veg_roughness, veg_albedo, veg_emiss1, veg_emiss2, snow_dirt, fice_soil_surf

 real, dimension(:,:,:), allocatable :: tg, qg, tg_first_guess, qg_rel_first_guess, fice_soil, &
 soil_map, veg_map, &
 soil_qmax, soil_qmin, soil_c, soil_rho, soil_psi, soil_k, soil_par_b, soil_par_c, soil_qrel_wilt, soil_qrel_ref, &
 snow_lev, snow_t, snow_fice, snow_age, snow_melt_age, snow_dens
 real :: veg_lai_max
 integer, dimension(:,:), allocatable :: ind_lev_soil_h_bottom, ind_lev_soil_w_bottom

! matrices with grid coordinates, for grid rotation

 real, dimension(:,:), allocatable   :: xxt, xxu, xxv, yyt, yyu, yyv, xxtg, yytg

! intermediate 3d fields (interpolated horizontally at isobaric levels).
! wind components are computed at u and v grid points, on both u and v grids

 real, dimension(:,:,:), allocatable :: gp, tp, upu, upv, vpu, vpv, upup, vpvp, qp

! Copies
 integer :: iflag_soil_snow_frc
 real, dimension(:,:), allocatable :: tskin_copy, iceth_frc, fice_frc, &
 tskin_frc, tgsurf_frc, qvsurf_frc, qgsurf_frc, fice_soil_surf_frc, snow_frc
 real, dimension(:,:,:), allocatable :: tg_copy, tg_frc, qg_frc, fice_soil_frc, &
 snow_lev_frc, snow_t_frc, snow_fice_frc, snow_age_frc, snow_melt_age_frc, snow_dens_frc

end module mod_fields
!##################################################################################################
   program preglobo

! Assimilation of global data in lat-lon coordinates

   use param
   use gribana
   use mod_fields

   implicit none

! length of file_name_work string must be 30 in according to declaration in subroutine read_grib2_data

   character (len=80) :: file_name_work

real, dimension(:,:), allocatable ::  workf, work, qgmin

! Intermediate (work) 2d fields interpolated on the model grid.
! they are derived (some only virtually) from the analyses,
! except htopi that is derived from the global orography dataset

real, dimension(:,:), allocatable   :: phigi, phigi2d, psil, snowi, htopi
real, dimension(:,:,:), allocatable :: tgi, qgi, tgii, qgii

! Sea surface temperature (sst) derived from analysis

real, dimension(:,:), allocatable :: sste, ssti

! Intermediate work fields

real, dimension(:,:,:), allocatable  :: tv, wspeed
real, dimension(:), allocatable      :: gaux
real, dimension(:,:), allocatable    :: tsair, fhard, diftop, twater
integer, dimension(3)                :: ijk

! Pressure at the surface (defined at t, u, v input grid point respectively)

real, dimension(:,:), allocatable :: pst, psu, psv

! Standard deviation of orography

 real, dimension(:,:), allocatable :: zhtopvar

! work matrix used for vertical interpolation

real, dimension(:), allocatable :: uk, vk, tk, qk, upk, vpk, qpk, tvk, pl2, &
                                    tsoile, qsoile, zsoil, zsoil_work, tsoili, qsoili

integer, dimension(:), allocatable :: iv, ivg

integer :: i, j, k, i0, j0, iii, ki, k1, k2, k3, k4, kf, kt, jklev, igrib, icentre_code, model_code, ndayr, jmon, nsmooth, &
           isubcentre_code, nyrc, nmonc, ndayc, nhouc, nminc, &
           flag_constant_fields, nlon_limit, nlsnow, nk

real :: x0a, y0a, x0=0., y0=0., zday, zcoeday, zlon_lim, zzz, zlatt, zlatv,            &
        alon0, alat0, alatcr0, alatcr25, gamma1, zlaps0, zt0t, zratio, zesk,     &
        al, alf, ex1, ex2, zphl1, zphl2, zphl3, zphl4, zpl1, zpl2, peso, peso2,        &
        delp1, delp2, delp3, tv1, tv2, tvs, psmin, zdiftp, zwf, zge, wei, tvm, zalmax, &
        psig, zqs, zqsw, zqsi, zes, zesw, zesi, qsat, wsat, w, rh, wm, ztlake,         &
        ws, topcr, tland, zfsnow, zterm, zterq, wf, ztland, ttme,  zqgmin, zqgmax,     &
        zlon, zlat, zdelt, qgav1, qgav2, qgav3, qgav4, zeskw, zdq,                 &
        zrgm, fatt, zeski, zqsat, zqtrue, zcaz, zzh1, zzh2, zs, z1, z2

   pi  = abs(acos(-1.))
   eps = rd/rv
   ep  = rv/rd-1.
   rcp = rd/cp

   open (11, file='preglobo.inp', status='old')
   read (11, param_preglobo)
   close (11)
   print *,'parameters of Preglobo (defined in preglobo.inp):'
   print param_preglobo
   print *,'  '

   nlevp1 = nlev+1
   ntot   = nlon*nlat

   dlon  = 360./float(nlon-2)
   dlat  = 180./float(nlat-1)
   alon0 = 0.
   alat0 = -90.

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

   open (11, file='input_data.grib2', status='old', err=50)
   close (11)

   file_name_work = 'input_data.grib2'
   call read_grib2_data(1, file_name_work, 1, .true., flag_cut_paste, icentre_code, isubcentre_code, model_code,&
             iana0, jana, nlev_atm_inp, nlev_atm_inp_max, nlev_soil_inp, nlev_soil_inp_max,          &
             x0a, y0a, x1a, y1a, dxa, dya, idate0, iperiod_inp, lev_list, lev_list_soil, level_type, &
             npar3d, npar2d, npar3d_soil, field3d, field2d, field3d_soil, ak, bk, val_missing)


   iana = iana0 + 1

   print *,'Parameters of the input grid are:'
   print *,'nx=', iana, 'ny=', jana,'  x1=', x1a, 'y1=', y1a
   print*, 'dlon=', dxa, 'dlat=', dya
   print *

   if (abs(x0a) > 0.001 .or. abs(y0a) > 0.001) then
   print *,"The input grid uses rotated coordinates => stop"
   stop
   endif

   if (dxa*float(iana) < 357. .or. dya*float(jana) < 179.) then
   print *,"Input data grid is not global => stop"
   stop
   endif

! Definition of input data origin (model_input)

   print*
   if (icentre_code == 98) then
   model_input='ifs'
   print *,'Input data are provided by IFS model, ECMWF'
   elseif (icentre_code == 7) then
   model_input='gfs'
   print *,'Input data are provided by GFS model, NOAA, NCEP, USA'
   else
   print *,'Emission Center of input data (',icentre_code,') not managed => stop'
   stop
   endif
   print*

!---------------------------------------------------------------------------
!  allocation of arrays
!---------------------------------------------------------------------------

 allocate(field3d(iana0,jana,nlev_atm_inp_max,npar3d))
 allocate(field3d_soil(iana0,jana,nlev_soil_inp_max,npar3d_soil))
 allocate(field2d(iana0,jana,npar2d))

 allocate(ge     (iana,jana,nlev_atm_inp_max))
 allocate(te     (iana,jana,nlev_atm_inp_max))
 allocate(ue     (iana,jana,nlev_atm_inp_max))
 allocate(ve     (iana,jana,nlev_atm_inp_max))
 allocate(qe     (iana,jana,nlev_atm_inp_max))
 allocate(qce    (iana,jana,nlev_atm_inp_max))
 allocate(qcwe   (iana,jana,nlev_atm_inp_max))
 allocate(qcie   (iana,jana,nlev_atm_inp_max))
 allocate(tge    (iana,jana,nlev_soil_inp_max))
 allocate(qge    (iana,jana,nlev_soil_inp_max))

 allocate(soile  (iana,jana))
 allocate(tskine (iana,jana))
 allocate(phige  (iana,jana))
 allocate(phige2d(iana,jana))
 allocate(fmaske (iana,jana))
 allocate(qgmine2d(iana,jana))
 allocate(qgmaxe2d(iana,jana))
 allocate(psel   (iana,jana))
 allocate(pmsle  (iana,jana))
 allocate(snowe  (iana,jana))
 allocate(ficee  (iana,jana))
 allocate(ficeee (iana,jana))
 allocate(ticee  (iana,jana))
 allocate(ticeee (iana,jana))
 allocate(xe     (iana))
 allocate(ye     (jana))
 allocate(alon_inp(iana,jana))
 allocate(alat_inp(iana,jana))
 allocate(mask_frame(iana,jana))
 allocate(cctote (iana,jana))
 allocate(u10e   (iana,jana))
 allocate(v10e   (iana,jana))
 allocate(t2e    (iana,jana))
 allocate(td2e   (iana,jana))
 allocate(sste   (iana,jana))

 allocate(p      (nlon,nlat,nlev_atm_inp_max))
 allocate(pl     (nlon,nlat,nlev_atm_inp_max))
 allocate(phl    (nlon,nlat,nlev_atm_inp_max+1))
 allocate(plaux  (2*nlev_atm_inp_max-1))
 allocate(sigint (nlev))
 allocate(aka    (nlev))
 allocate(bika   (nlev))
 allocate(s      (nlevp1))
 allocate(akah   (nlevp1))
 allocate(bikah  (nlevp1))
 allocate(plsig  (nlev))
 allocate(hxt    (nlat))
 allocate(hxv    (nlat))

 allocate(ak     (nlev_atm_inp_max+1))
 allocate(bk     (nlev_atm_inp_max+1))
 allocate(lev_list(nlev_atm_inp_max))
 allocate(lev_list_soil(nlev_soil_inp_max))

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

 allocate(u      (nlon,nlat,nlev))
 allocate(v      (nlon,nlat,nlev))
 allocate(t      (nlon,nlat,nlev))
 allocate(q      (nlon,nlat,nlev))
 allocate(qc     (nlon,nlat,nlev))
 allocate(press  (nlon,nlat,nlev))

 allocate(phig      (nlon,nlat))
 allocate(fmask     (nlon,nlat))
 allocate(fmaski    (nlon,nlat))
 allocate(tskin     (nlon,nlat))
 allocate(tskin_copy(nlon,nlat))
 allocate(tskin_frc (nlon,nlat))
 allocate(tgsurf    (nlon,nlat))
 allocate(tgsurf_frc(nlon,nlat))
 allocate(qvsurf    (nlon,nlat))
 allocate(qvsurf_frc(nlon,nlat))
 allocate(qgsurf    (nlon,nlat))
 allocate(qgsurf_frc(nlon,nlat))
 allocate(ps        (nlon,nlat))
 allocate(cloud     (nlon,nlat))
 allocate(totpre    (nlon,nlat))
 allocate(conpre    (nlon,nlat))
 allocate(snfall    (nlon,nlat))
 allocate(snow      (nlon,nlat))
 allocate(snow_frc  (nlon,nlat))
 allocate(fice      (nlon,nlat))
 allocate(fice_frc  (nlon,nlat))
 allocate(iceth     (nlon,nlat))
 allocate(iceth_frc (nlon,nlat))
 allocate(rgm       (nlon,nlat))
 allocate(rgq       (nlon,nlat))
 allocate(tice      (nlon,nlat))
 allocate(albedo    (nlon,nlat))
 allocate(emismap1  (nlon,nlat))
 allocate(emismap2  (nlon,nlat))
 allocate(snow_albedo(nlon,nlat))
 allocate(fice_soil_surf    (nlon,nlat))
 allocate(fice_soil_surf_frc(nlon,nlat))

 allocate(tg_first_guess    (nlon,nlat,nlevg))
 allocate(tg           (nlon,nlat,nlevg))
 allocate(tg_copy      (nlon,nlat,nlevg))
 allocate(tg_frc       (nlon,nlat,nlevg))
 allocate(qg_rel_first_guess(nlon,nlat,nlevg))
 allocate(qg           (nlon,nlat,nlevg))
 allocate(qg_frc       (nlon,nlat,nlevg))
 allocate(fice_soil    (nlon,nlat,nlevg))
 allocate(fice_soil_frc(nlon,nlat,nlevg))

 allocate(snow_lev         (nlon,nlat,nlevsnow))
 allocate(snow_lev_frc     (nlon,nlat,nlevsnow))
 allocate(snow_t           (nlon,nlat,nlevsnow))
 allocate(snow_t_frc       (nlon,nlat,nlevsnow))
 allocate(snow_fice        (nlon,nlat,nlevsnow))
 allocate(snow_fice_frc    (nlon,nlat,nlevsnow))
 allocate(snow_age         (nlon,nlat,nlevsnow))
 allocate(snow_age_frc     (nlon,nlat,nlevsnow))
 allocate(snow_melt_age    (nlon,nlat,nlevsnow))
 allocate(snow_melt_age_frc(nlon,nlat,nlevsnow))
 allocate(snow_dens        (nlon,nlat,nlevsnow))
 allocate(snow_dens_frc    (nlon,nlat,nlevsnow))
 allocate(snow_dirt        (nlon,nlat))

 allocate(ind_lev_soil_h_bottom(nlon,nlat))
 allocate(ind_lev_soil_w_bottom(nlon,nlat))

 allocate(water_table_depth  (nlon,nlat))
 allocate(tg_bottom          (nlon,nlat))
 allocate(qg_rel_bottom      (nlon,nlat))
 allocate(qg_rel_surf_approx (nlon,nlat))

 allocate(soil_map (nlon,nlat,nst+1))
 allocate(veg_map  (nlon,nlat,nvt+1))

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

 allocate(workf  (iana,jana))
 allocate(work   (nlon,nlat))

 allocate(phigi  (nlon,nlat))
 allocate(phigi2d(nlon,nlat))
 allocate(psil   (nlon,nlat))
 allocate(snowi  (nlon,nlat))
 allocate(htopi  (nlon,nlat))
 allocate(zhtopvar(nlon,nlat))
 allocate(tgi    (nlon,nlat,nlevg))
 allocate(qgi    (nlon,nlat,nlevg))
 allocate(tgii   (nlon,nlat,nlev_soil_inp_max))
 allocate(qgii   (nlon,nlat,nlev_soil_inp_max))
 allocate(ssti   (nlon,nlat))

 allocate(tv     (nlon,nlat,2*nlev_atm_inp_max-1))
 allocate(gaux   (2*nlev_atm_inp_max-1))
 allocate(tsair  (nlon,nlat))
 allocate(fhard  (nlon,nlat))
 allocate(diftop (nlon,nlat))
 allocate(twater (nlon,nlat))
 allocate(wspeed (nlon,nlat,nlev))
 allocate(gamma  (nlon,nlat))
 allocate(qgmin  (nlon,nlat))

 allocate(pst    (nlon,nlat))
 allocate(psu    (nlon,nlat))
 allocate(psv    (nlon,nlat))

 allocate(uk     (nlev))
 allocate(vk     (nlev))
 allocate(tk     (nlev))
 allocate(qk     (nlev))
 allocate(upk    (nlev_atm_inp_max))
 allocate(vpk    (nlev_atm_inp_max))
 allocate(qpk    (nlev_atm_inp_max))
 allocate(tvk    (2*nlev_atm_inp_max-1))
 allocate(pl2    (nlev_atm_inp_max))

 allocate(tsoile (nlev_soil_inp_max+1))
 allocate(qsoile (nlev_soil_inp_max+1))
 allocate(zsoil      (nlevg))
 allocate(zsoil_work (nlevg))
 allocate(tsoili     (nlevg))
 allocate(qsoili     (nlevg))

 allocate(iv     (nlev))
 allocate(ivg    (nlevg))

!---------------------------------------------------------------------------

   nfdr = 0
   pdr  = 0.

   snow_albedo(:,:) = 0.71

! Definition of hybrid half-integer and integer atmosferic levels

   s(1) = 0.

   vert_coord_c = 0.

! 1. Case of simple sigma standard coordinates = equal-mass layers

   if (ivert_coord_type <= 1) then
     do k = 2, nlevp1
       s(k) = float(k-1)/float(nlev)
     enddo
   endif

! 2. Case of bolam standard coordinates

   if (ivert_coord_type == 2) then
     do k = 2, nlevp1
       zs = float(k-1)/float(nlev)
       s(k) = 0.78*zs + 1.44*zs**3 - 1.22*zs**4
     enddo
   endif

! 3. Case of vert_coord_c using: vert_coord_c flattens hybrid surfaces more rapidly above mountains at stratospheric levels

   if (ivert_coord_type == 3) then
     vert_coord_c = 0.09**3*(1.+alfa*(1.-0.09)*(2.-0.09))/(1.-0.09)
     do k = 2, nlevp1
       zs = float(k-1)/float(nlev)
       s(k) = 0.78*zs + 1.44*zs**3 - 1.22*zs**4
       s(k) = s(k) * max ((1.+tanh((zs-0.5)/0.7))/(1.+tanh((1.-0.5)/0.7)), 0.02)
     enddo
   endif

! At integer levels:

   do k = 1, nlev
     sigint(k) = 0.5*(s(k+1)+s(k))
     bika(k) = max( sigint(k)**3 * (1.+alfa*(1.-sigint(k))*(2.-sigint(k))) - vert_coord_c*(1.-sigint(k)), 0.)
     aka (k)  = pzer*(sigint(k) - bika(k))
   enddo

! At semi-integer levels:

   bikah(1) = 0.
   bikah(nlevp1) = 1.
   do k = 2, nlev
     bikah(k) = 0.5*(bika(k)+bika(k-1))
   enddo

! Definition of horizontal grid variables

   do j = 2, nlat-1
   zlatt  = pi/180.*(-90.+(j-1)*dlat)
   zlatv  = zlatt-dlat/2.*pi/180.
   hxt(j) = cos(zlatt)
   hxv(j) = cos(zlatv)
   enddo
   hxt(1   ) = 0.
   hxt(nlat) = 0.
   zlatv = pi/180.*(90.-dlat/2.)
   hxv(nlat) = cos(zlatv)

! initialization of some parameters and flags of mhf

   nfdr(1) = 1
   nfdr(2) = nlon
   nfdr(3) = nlat
   nfdr(4) = nlev
   nfdr(13)= npolcap
   nfdr(15)= nlevg

   pdr(1)  = dlat
   pdr(2)  = dlon
   pdr(4)  = -90.-dlat/2.
   pdr(5)  = 0.
   pdr(6:5+nlevg) = soil_lev(1:nlevg)
   pdr(36) = pzer
   pdr(37) = alfa
   pdr(38) = 0.
   pdr(39) = 0.
   pdr(199) = vert_coord_c

   do k = 2, nlevp1
   if (38+k.gt.200) then ! 200 is the present dimension of pdr
   print*, "The number of levels of the output model exceeds the upper limit"
   stop
   endif
   pdr(38+k) = s(k)
   enddo

! definition of coordinate grid used in the model

   zlon_lim = 360.
   if (x1a < 0.) zlon_lim = 180.

   do j = 1, nlat
   do i = 1, nlon
   xxt(i,j) = float(i-1)*dlon
   xxu(i,j) = xxt(i,j)+dlon/2.
   if (xxt(i,j) >= zlon_lim) xxt(i,j) = xxt(i,j)-360.
   if (xxu(i,j) >= zlon_lim) xxu(i,j) = xxu(i,j)-360.
   xxv(i,j) = xxt(i,j)
   yyt(i,j) = -90. + float(j-1)*dlat
   yyv(i,j) = yyt(i,j)-dlat/2.
   if (yyt(i,j) >=  90.-1.e-4) yyt(i,j)= 90.-1.e-4
   if (yyv(i,j) <  -90.+1.e-4) yyv(i,j)=-90.+1.e-4
   if (yyt(i,j) <  -90.+1.e-4) yyt(i,j)=-90.+1.e-4
   yyu(i,j) = yyt(i,j)
   enddo
   enddo

   xxtg(:,:) = xxt(:,:)
   yytg(:,:) = yyt(:,:)

! definition of geographical coordinates (in deg.) of the input data grid

   do i = 1, iana
   xe(i) = x1a + float(i-1)*dxa
   alon_inp(i,:) = xe(i)
   enddo
   do j = 1, jana
   ye(j) = y1a + float(j-1)*dya
   alat_inp(:,j) = ye(j)
   enddo

!--------------------------------------------------------------------------------
! call of the procedure for reading and decoding input data files in grib2 format
!--------------------------------------------------------------------------------

   file_name_work = 'input_data.grib2'
   call read_grib2_data(1, file_name_work, 1, .false., flag_cut_paste, icentre_code, isubcentre_code, model_code,&
             iana0, jana, nlev_atm_inp, nlev_atm_inp_max, nlev_soil_inp, nlev_soil_inp_max,          &
             x0a, y0a, x1a, y1a, dxa, dya, idate0, iperiod_inp, lev_list, lev_list_soil, level_type, &
             npar3d, npar2d, npar3d_soil, field3d, field2d, field3d_soil, ak, bk, val_missing)

   iana = iana0 + 1
   nlev_atm_inpp1 = nlev_atm_inp + 1
   nauxp = 2*nlev_atm_inp - 1

! preliminary elaboration of input data

   if (model_input=='ifs') call conv_ifs_data
   if (model_input=='gfs') call conv_gfs_data
!   call plotout (fmaske, iana, jana, 99)

   print*, iperiod(1), iperiod(2), iperiod(3)

 call calendar (idate0(1), idate0(2), idate0(3),  idate0(4) ,idate0(5), iperiod(1), iperiod(2), iperiod(3), &
                nyrc, nmonc, ndayc, nhouc, nminc, ndayr)

 nfdr(5) = nyrc
 nfdr(6) = nmonc
 nfdr(7) = ndayc
 nfdr(8) = nhouc
 nfdr(9) = nminc

 print *,' Initial condition date (yyyy/mm/dd/hh/min): '
 print *, nfdr(5:9)
 print *,' Ndayr: ', ndayr

 print*
 if (level_type == 1) then
   print *,'type of atmospheric levels: isobaric'
 elseif (level_type == 2) then
   print *,'type of atmospheric levels: hybrid'
 else
   print *,'type of atmospheric levels not defined'
 endif
 print*

 do j = 1, nlat
 zlaps0 = gamma0*sqrt(hxt(j))
 do i = 1, nlon
 gamma(i,j) = zlaps0*(1.-0.2*sin(2.*yyt(i,j)*pi/180.)*cos(ndayr*2.*pi/365.))
 enddo
 enddo

! Seatemp extends t from the sea towards the land (without changing tge(1) )
! first over land t is prescribed as a combination of tge(2) and tge(3) is to reduce
! sst problems when sst is used later to define lake temperature.

   !!!!!!!ticeee=val_missing!!!!!!!!!!

!   call plotout (tskine, iana, jana, 99)
!   call plotout (tge(:,:,1), iana, jana, 99)
!   call plotout (tge(:,:,2), iana, jana, 99)
!   call plotout (tge(:,:,3), iana, jana, 99)
!   call plotout (ficeee, iana, jana, 99)
!   call plotout (fmaske, iana, jana, 99)
!   call plotout (ticeee, iana, jana, 99)

   zcaz = ticeee(1,1)
   do j = 1, jana
   do i = 1, iana

   if (fmaske(i,j) >= .5) then
   workf(i,j) = tge(i,j,1)
   else
   workf(i,j) = 0.5*tge(i,j,2) + 0.5*tge(i,j,3)
   endif

   if (ficeee(i,j).lt.0.05) then
   ticeee(i,j) = 271.4
   else
   if (int(zcaz).eq.int(val_missing)) ticeee(i,j) = min(workf(i,j),271.4)
   endif

   enddo
   enddo

   call seatemp (workf, fmaske, sste, iana, jana, 5, 1, .5)

! Ice sea fraction, ice temperature and surface level temperature

   if (model_input=='ifs') then

     call seatemp(ficeee, fmaske, ficee, iana, jana, 5, 1, .5)
     call seatemp(ticeee, fmaske, ticee, iana, jana, 5, 1, .5)

   elseif (model_input=='gfs') then

     do j = 1, jana
     do i = 1, iana
     if (fmaske(i,j) > 0.5.and.tge(i,j,1) < 271.4.and.ficeee(i,j) < 0.5) &
         ficeee(i,j) = max (ficeee(i,j), min (fmaske(i,j), (271.4-tge(i,j,1))/10.))
     if (fmaske(i,j) > 0.5.and.tge(i,j,1) > 271.4.and.ficeee(i,j) > 0.5) &
         tge(i,j,1) = min (271.4, tge(i,j,1))
     enddo
     enddo

! Definition of ticeee from tge (same as TskinE over sea, including sea ice)

     ticeee(:,:) = tge(:,:,1)

     do j = 1, jana
     do i = 1, iana
     ticeee(i,j) = max (ticeee(i,j), te(i,j,nlev_atm_inp)-10.)
     ticeee(i,j) = min (ticeee(i,j), 271.4)
     enddo
     enddo

     ficee(:,:) = ficeee(:,:)
     call seatemp (ticeee, ficee, ticee, iana, jana, 5, 1, .5)

   endif

! landtemp extends t and qg from the land towards the sea:
! these variables are redefined, so they must not be used to define quantities over the sea

   call landtemp (tskine, fmaske, iana, jana, 4, 0, .8)
   do k = 1, nlev_soil_inp
   call landtemp (tge(1,1,k), fmaske, iana, jana, 4, 0, .8)
   call landtemp (qge(1,1,k), fmaske, iana, jana, 4, 0, .8)
   enddo

! --------------------------------------------------
! horizontal interpolation of surface fields
! --------------------------------------------------

! note: phige and phige2d are the same for gfs data and different for ifs data, case of hybrid levels

   al = 0.8
   if (level_type == 2) call interp_spline_2d (phige,iana,jana,xe,ye,xxt,yyt,ntot,phigi,al) ! hybrid levels
   print *, "phige2d max = ",maxval(phige2d) 
   print *, "phige2d min = ",minval(phige2d) 
   call interp_spline_2d (phige2d,iana,jana,xe,ye,xxt,yyt,ntot,phigi2d,al)

   do k = 1, nlev_soil_inp
   call interp_spline_2d (tge(1,1,k), iana, jana, xe, ye, xxt, yyt, ntot, tgii(1,1,k), al)
! qge - relative (!) soil water content (proportion) in input data;
! qgii - relative (!) soil water content (proportion) interpolated at model grid
   call interp_spline_2d (qge(1,1,k), iana, jana, xe, ye, xxt, yyt, ntot, qgii(1,1,k), al)
   enddo

   if (minval(tskine) > 180.) then ! there are input data for surface temperature
     call interp_spline_2d (tskine,iana,jana,xe,ye,xxt,yyt,ntot,tskin,al)
   else
     tskin(:,:) = tgii(:,:,1)
   endif
   call interp_spline_2d (sste,iana,jana,xe,ye,xxt,yyt,ntot,ssti,al)
   call interp_spline_2d (ticee,iana,jana,xe,ye,xxt,yyt,ntot,tice,al)
   call interp_spline_2d (fmaske,iana,jana,xe,ye,xxt,yyt,ntot,fmaski,1.)

   al = 1.
   call interp_spline_2d (snowe, iana, jana, xe, ye, xxt, yyt, ntot, snowi, al)
   call interp_spline_2d (ficee, iana, jana, xe, ye, xxt, yyt, ntot, fice , al)

! -----------------------------------------------------------------------

! Definition of all model physiographical parameters: land-sea fraction,
! orography parameters, soil and vegetation parameters

! If model_param_constant.bin exists (flag_constant_fields=0), then reads constant model
! physiographical parameters from this file,
! else (flag_constant_fields/=0) parameters are defined using 
! subroutine physiographic_param (in def_soil.F90)

! Physiographical parameters variable in time (LAI, vegetation frac,
! soil temperature and soil water content vertical approximated profies)
! are difined by subroutine physiographic_param (in def_soil.F90)

   call read_param_const(x0, y0, alon0, alat0, htopi, zhtopvar, flag_constant_fields)

   nlon_limit=nlon-2

   print *

   call physiographic_param (nlon_limit, nlat, nst, nvt, nlevg, &
 alon0, alat0, dlon, dlat, x0, y0, xxtg(1:nlon_limit,:), yytg(1:nlon_limit,:), &
 tskin(1:nlon_limit,:), soil_lev, ndayr, flag_constant_fields, &
 htopi(1:nlon_limit,:), zhtopvar(1:nlon_limit,:), fmask(1:nlon_limit,:), &
 water_table_depth(1:nlon_limit,:), tg_bottom(1:nlon_limit,:), &
 qg_rel_bottom(1:nlon_limit,:), qg_rel_surf_approx(1:nlon_limit,:), &
 ind_lev_soil_h_bottom(1:nlon_limit,:), ind_lev_soil_w_bottom(1:nlon_limit,:), &
 tg_first_guess(1:nlon_limit,:,:), qg_rel_first_guess(1:nlon_limit,:,:), &
 soil_map(1:nlon_limit,:,:), veg_map(1:nlon_limit,:,:), &
 soil_qmax(1:nlon_limit,:,:), soil_qmin(1:nlon_limit,:,:), soil_c(1:nlon_limit,:,:), &
 soil_rho(1:nlon_limit,:,:), soil_psi(1:nlon_limit,:,:), soil_k(1:nlon_limit,:,:), &
 soil_par_b(1:nlon_limit,:,:), soil_par_c(1:nlon_limit,:,:), &
 soil_qrel_wilt(1:nlon_limit,:,:), soil_qrel_ref(1:nlon_limit,:,:), &
 soil_albedo_dry(1:nlon_limit,:), soil_albedo_wet(1:nlon_limit,:), &
 soil_emiss1_dry(1:nlon_limit,:), soil_emiss1_wet(1:nlon_limit,:), &
 soil_emiss2_dry(1:nlon_limit,:), soil_emiss2_wet(1:nlon_limit,:), &
 veg_lai(1:nlon_limit,:), veg_frac(1:nlon_limit,:), &
 veg_root_depth(1:nlon_limit,:), veg_roughness(1:nlon_limit,:), &
 veg_albedo(1:nlon_limit,:), veg_emiss1(1:nlon_limit,:), veg_emiss2(1:nlon_limit,:), veg_lai_max, &
 snow_dirt(1:nlon_limit,:))

   htopi(nlon-1:nlon,:) = htopi(1:2,:)
   zhtopvar(nlon-1:nlon,:) = zhtopvar(1:2,:)
   fmask(nlon-1:nlon,:) = fmask(1:2,:)
   tg_bottom(nlon-1:nlon,:) = tg_bottom(1:2,:)
   qg_rel_bottom(nlon-1:nlon,:) = qg_rel_bottom(1:2,:)
   qg_rel_surf_approx(nlon-1:nlon,:) = qg_rel_surf_approx(1:2,:)
   water_table_depth(nlon-1:nlon,:) = water_table_depth(1:2,:)
   ind_lev_soil_h_bottom(nlon-1:nlon,:) = ind_lev_soil_h_bottom(1:2,:)
   ind_lev_soil_w_bottom(nlon-1:nlon,:) = ind_lev_soil_w_bottom(1:2,:)
   soil_map(nlon-1:nlon,:,:) = soil_map(1:2,:,:)
   veg_map(nlon-1:nlon,:,:) = veg_map(1:2,:,:)
   soil_qmax(nlon-1:nlon,:,:) = soil_qmax(1:2,:,:)
   soil_qmin(nlon-1:nlon,:,:) = soil_qmin(1:2,:,:)
   soil_c(nlon-1:nlon,:,:) = soil_c(1:2,:,:)
   soil_rho(nlon-1:nlon,:,:) = soil_rho(1:2,:,:)
   soil_psi(nlon-1:nlon,:,:) = soil_psi(1:2,:,:)
   soil_k(nlon-1:nlon,:,:) = soil_k(1:2,:,:)
   soil_par_b(nlon-1:nlon,:,:) = soil_par_b(1:2,:,:)
   soil_par_c(nlon-1:nlon,:,:) = soil_par_c(1:2,:,:)
   soil_qrel_wilt(nlon-1:nlon,:,:) = soil_qrel_wilt(1:2,:,:)
   soil_qrel_ref(nlon-1:nlon,:,:) = soil_qrel_ref(1:2,:,:)
   soil_albedo_dry(nlon-1:nlon,:) = soil_albedo_dry(1:2,:)
   soil_albedo_wet(nlon-1:nlon,:) = soil_albedo_wet(1:2,:)
   soil_emiss1_dry(nlon-1:nlon,:) = soil_emiss1_dry(1:2,:)
   soil_emiss1_wet(nlon-1:nlon,:) = soil_emiss1_wet(1:2,:)
   soil_emiss2_dry(nlon-1:nlon,:) = soil_emiss2_dry(1:2,:)
   soil_emiss2_wet(nlon-1:nlon,:) = soil_emiss2_wet(1:2,:)
   veg_lai(nlon-1:nlon,:) = veg_lai(1:2,:)
   veg_frac(nlon-1:nlon,:) = veg_frac(1:2,:)
   veg_root_depth(nlon-1:nlon,:) = veg_root_depth(1:2,:)
   veg_roughness(nlon-1:nlon,:) = veg_roughness(1:2,:)
   veg_albedo(nlon-1:nlon,:) = veg_albedo(1:2,:)
   veg_emiss1(nlon-1:nlon,:) = veg_emiss1(1:2,:)
   veg_emiss2(nlon-1:nlon,:) = veg_emiss2(1:2,:)
   snow_dirt(nlon-1:nlon,:) = snow_dirt(1:2,:)

   tg_first_guess(nlon-1:nlon,:,:) = tg_first_guess(1:2,:,:)
   qg_rel_first_guess(nlon-1:nlon,:,:) = qg_rel_first_guess(1:2,:,:)
   veg_lai(nlon-1:nlon,:) = veg_lai(1:2,:)
   veg_frac(nlon-1:nlon,:) = veg_frac(1:2,:)
 
   call polavert (phig, nlon, nlat, npolcap)
   call polavert (zhtopvar, nlon, nlat, npolcap)
   zhtopvar(:,:) = max (zhtopvar(:,:), 0.)
   call polavert (veg_roughness, nlon, nlat, npolcap)
   veg_roughness(:,:) = max (veg_roughness(:,:), 1.e-8)

   phig(:,:) = htopi(:,:)*g0

   print *

! -----------------------------------------------------------------------

   work = fice
   call seatemp (work, fmask, fice, nlon, nlat, 5, 1, .5)
   work = tice
   call seatemp (work, fmask, tice, nlon, nlat, 5, 1, .5)
   work = tice
   call seatemp (work, fice, tice, nlon, nlat, 5, 1, .5)
   do k = 1, nlev_soil_inp
   work(:,:) = tgii(:,:,k)
   call landtemp (tgii(1,1,k), fmask, nlon, nlat, 4, 0, .8)
   do j = 1, nlat
   do i = 1, nlon
   if (fmask(i,j).gt.0.5) tgii(i,j,k) = work(i,j)
   enddo
   enddo
   call landtemp (qgii(1,1,k), fmask, nlon, nlat, 4, 0, .8)
   enddo

! --------------------------------------------------
! Vertical interpolation of soil parameters (temp. and water cont.)
! --------------------------------------------------

! LEV_LIST_SOIL: levels of input model soil
! ZSOILI and SLT: levels and depths of output model soil layers

   zsoil(1:nlevg) = sqrt(soil_lev(1:nlevg)) ! to make the depth coordinate more linear

! Vertical interpolation of soil variables from input to output model levels

  alf = 0.
  ex1 = 1.
  ex2 = 1.
  z1=sum(lev_list_soil(1:nlev_soil_inp))
  z2=sum(soil_lev(1:nlevg))

  if (abs(z1-z2)>1.e-1.or.nlev_soil_inp<nlevg) then  ! case vertical interpolation is needed

    zsoil(1:nlevg) = sqrt(soil_lev(1:nlevg))
    zsoil_work(1:nlev_soil_inp) = sqrt(lev_list_soil(1:nlev_soil_inp))

    do j = 1, nlat
    do i = 1, nlon
  
! Temperature of soil and sea water

     tsoile(1:nlev_soil_inp) = tgii(i,j,1:nlev_soil_inp)

     if (fmask(i,j) < 0.5.or.fice(i,j) >= 0.8) then

! Soil or thick sea ice

        tgi(i,j,1:nlevg) = tg_first_guess(i,j,1:nlevg)
  
        nk = ind_lev_soil_h_bottom(i,j)

        if (soil_lev(nk) >= lev_list_soil(nlev_soil_inp)+0.5) then

! The bottom temperature soil level is below the bottom soil level in input data:
! temperature at upper soil levels is defined by vertical interpolation of
! input data and temperature at lower soil level is defined by first guess temperature

          do jklev = 1,nk
            if (soil_lev(jklev) >= lev_list_soil(nlev_soil_inp)+0.5) exit
          enddo
          nk = jklev
          zsoil_work(nlev_soil_inp+1) = zsoil(nk)
          tsoile(nlev_soil_inp+1) = tg_first_guess(i,j,nk)
          call near (zsoil(1:nk), nk, zsoil_work(1:nlev_soil_inp+1), nlev_soil_inp+1, ivg(1:nk))
          call interp_spline_1d (tsoili(1:nk), zsoil(1:nk), nk, tsoile(1:nlev_soil_inp+1), zsoil_work(1:nlev_soil_inp+1), &
   nlev_soil_inp+1, ivg(1:nk), alf, ex1, ex2)
          tgi(i,j,1:nk) = tsoili(1:nk)

        else

! The bottom temperature soil level is upper then bottom soil level in input data:
! temperature at all soil levels is defined by vertical interpolation of input data

          call near (zsoil(1:nk), nk, zsoil_work(1:nlev_soil_inp), nlev_soil_inp, ivg(1:nk))
          call interp_spline_1d (tsoili(1:nk), zsoil(1:nk), nk, tsoile(1:nlev_soil_inp), zsoil_work(1:nlev_soil_inp), &
   nlev_soil_inp, ivg(1:nk), alf, ex1, ex2)
          tgi(i,j,1:nk) = tsoili(1:nk)

        endif

      else

! Sea water

        call near (zsoil(1:nlevg), nlevg, zsoil_work(1:nlev_soil_inp), nlev_soil_inp, ivg(1:nlevg))
        call interp_spline_1d (tsoili(1:nlevg), zsoil(1:nlevg), nlevg, tsoile(1:nlev_soil_inp), zsoil_work(1:nlev_soil_inp), &
   nlev_soil_inp, ivg(1:nlevg), alf, ex1, ex2)
        tgi(i,j,1:nlevg) = tsoili(1:nlevg)

      endif

! Relative soil water content

      qsoile(1:nlev_soil_inp) = qgii(i,j,1:nlev_soil_inp)

      nk = ind_lev_soil_w_bottom(i,j)

      if (nk > 0) then ! No water body, no glacier
    
        qgi(i,j,1:nlevg) = qg_rel_first_guess(i,j,1:nlevg)

        if (soil_lev(nk) >= lev_list_soil(nlev_soil_inp)+0.5) then
  
! The bottom soil water content level is below the bottom soil level in input data:
! water content at upper soil levels is defined by vertical interpolation of
! input data and water content at lower soil level is defined by first guess temperature
  
          do jklev = 1,nk
            if (soil_lev(jklev) >= lev_list_soil(nlev_soil_inp)+0.5) exit
          enddo
          nk = jklev
          zsoil_work(nlev_soil_inp+1) = zsoil(nk)
          qsoile(nlev_soil_inp+1) = qg_rel_first_guess(i,j,nk)
          call near (zsoil(1:nk), nk, zsoil_work(1:nlev_soil_inp+1), nlev_soil_inp+1, ivg(1:nk))
          call interp_spline_1d (qsoili(1:nk), zsoil(1:nk), nk, qsoile(1:nlev_soil_inp+1), zsoil_work(1:nlev_soil_inp+1), &
     nlev_soil_inp+1, ivg(1:nk), alf, ex1, ex2)
          qgi(i,j,1:nk) = qsoili(1:nk)
  
        else
  
! The bottom soil water content level is upper then bottom soil level in input data:
! water content at all soil levels is defined by vertical interpolation of input data
  
          call near (zsoil(1:nk), nk, zsoil_work(1:nlev_soil_inp), nlev_soil_inp, ivg(1:nk))
          call interp_spline_1d (qsoili(1:nk), zsoil(1:nk), nk, qsoile(1:nlev_soil_inp), zsoil_work(1:nlev_soil_inp), &
     nlev_soil_inp, ivg(1:nk), alf, ex1, ex2)
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

! reset of possible unfeasible values

   snowi(:,:)  = max (snowi(:,:), 0.)
   fice(:,:)   = max (fice(:,:), 0.)
   fice(:,:)   = min (fice(:,:), fmask(:,:))
   tice(:,:)   = min (tice(:,:), 271.4)
   qgi(:,:,:)  = min (max(qgi(:,:,:),0.), 1.) ! relative value

!i=2; j=588
!i=552; j=890
!print *
!print *,'poisk fmask',fmask(i,j)
!print *,'poisk htopi',phig(i,j)/g0
!print *,'poisk zhtopvar',zhtopvar(i,j)
!print *,'poisk tg_bottom, qg_rel_bottom, qg_rel_surf_approx, water_table_depth',&
! tg_bottom(i,j), qg_rel_bottom(i,j), qg_rel_surf_approx(i,j), water_table_depth(i,j)
!print *,'poisk ind_lev_soil_h_bottom ind_lev_soil_w_bottom',ind_lev_soil_h_bottom(i,j),ind_lev_soil_w_bottom(i,j)
!print *,'poisk soil_map',soil_map(i,j,:)
!print *,'poisk veg_map',veg_map(i,j,:)
!print *,'poisk soil_qmax',soil_qmax(i,j,:)
!print *,'poisk soil_qmin',soil_qmin(i,j,:)
!print *,'poisk soil_c',soil_c(i,j,:)
!print *,'poisk soil_rho',soil_rho(i,j,:)
!print *,'poisk soil_psi',soil_psi(i,j,:)
!print *,'poisk soil_k',soil_k(i,j,:)
!print *,'poisk soil_par_b',soil_par_b(i,j,:)
!print *,'poisk soil_qrel_wilt',soil_qrel_wilt(i,j,:)
!print *,'poisk soil_qrel_ref',soil_qrel_ref(i,j,:)
!print *,'poisk soil_albedo_dry',soil_albedo_dry(i,j)
!print *,'poisk soil_albedo_dry',soil_albedo_dry(i,j)
!print *,'poisk soil_albedo_wet',soil_albedo_wet(i,j)
!print *,'poisk soil_emiss1_dry',soil_emiss1_dry(i,j)
!print *,'poisk soil_emiss1_wet',soil_emiss1_wet(i,j)
!print *,'poisk soil_emiss2_dry',soil_emiss2_dry(i,j)
!print *,'poisk soil_emiss2_wet',soil_emiss2_wet(i,j)
!print *,'poisk veg_root_depth',veg_root_depth(i,j)
!print *,'poisk veg_roughness',veg_roughness(i,j)
!print *,'poisk veg_albedo',veg_albedo(i,j)
!print *,'poisk veg_emiss1',veg_emiss1(i,j)
!print *,'poisk veg_emiss2',veg_emiss2(i,j)
!print *,'poisk snow_dirt',snow_dirt(i,j)
!print *,'poisk veg_lai_max',veg_lai_max
!print *,'poisk tg_first_guess',tg_first_guess(i,j,:)
!print *,'poisk qg_rel_first_guess',qg_rel_first_guess(i,j,:)
!print *,'poisk veg_lai',veg_lai(i,j)
!print *,'poisk veg_frac',veg_frac(i,j)
!print *,'poisk tgi',tgi(i,j,:)
!print *,'poisk qgi',qgi(i,j,:)
!print *

! --------------------------------------------------
!  only in the case of hybrid levels: horizontal interpolation of surface pressure
! --------------------------------------------------

   if (level_type == 2) then
   al = 0.8
   call interp_spline_2d (psel,iana,jana,xe,ye,xxt,yyt,ntot,psil,al)
   endif

! ----------------------------------------------
! Definition of the fields on vertical levels
! ----------------------------------------------

! Note that the same 3d matrices p and pl are used for pressure on input model levels,
! for both cases of input pressure levels (in this case only the vertical index varies)
! and of hybrid levels (considered later)
! Note also that the no. of press. levels can be different at different instants

 if (level_type==1) then  ! isobaric levels

! definition of auxiliary isobaric levels (logaritmic interp.):
! p(1,1,k)      standard levels
! pl(1,1,k)     logp(k)
! plaux(k)      log paux (interpolated)

   i0=1
   j0=1
   do k=1,nlev_atm_inp
     p(i0,j0,k)=lev_list(k)
     pl(i0,j0,k)=alog(p(i0,j0,k))
     plaux(2*k-1)=pl(i0,j0,k)
   enddo
   do k=1,nlev_atm_inp-1
     plaux(2*k)=(pl(i0,j0,k)+pl(i0,j0,k+1))/2.
   enddo

! Recomputation of geopotential values at isobaric levels below ground surface:
! sometime these values are not consistent with the average temperature in the
! layer above, so using them for computing the temperature at auxiliary levels can
! produce unrealistic values; therefore, the geopotential is recomputed using the
! virtual temperature of the layer above.
! the correction is applied depending on the level difference.

   do j=1,jana
   do i=1,iana

! Search of the index of the first level above ground

     ki=0
     do k=1,nlev_atm_inp
       if(phige2d(i,j) < ge(i,j,k)) ki=k
     enddo

! Geopotential correction starting from above (the updated geopotential is
! computed starting from the first level located below ground ki+1).
! a compromise solution is applied: weighted average between the old and the
! new geopotential, corrected as a function of the distance between the gfs orography
! and the level considered (zwf=0.: no correct.; zwf=1.: full correct.)

     if (ki < nlev_atm_inp) then
       do k=ki+1,nlev_atm_inp
         zdiftp = max(phige2d(i,j)-ge(i,j,k), 0.)
         zwf = min(zdiftp/1500.,1.)
         zge=ge(i,j,k)
         ge(i,j,k) = ge(i,j,k-1) + alog(p(i0,j0,k-1)/p(i0,j0,k))*rd*.5* &
 ((1.+ep*qe(i,j,k))*te(i,j,k)+(1.+ep*qe(i,j,k-1))*te(i,j,k-1))
         ge(i,j,k) = zwf*ge(i,j,k)+(1.-zwf)*zge
       enddo
     endif

   enddo
   enddo

   if (model_input=='gfs') then

! horizontal filtering of the upper level (atmospheric) input fields
! phige used as workspace

     print*, "Horizontal smoothing applied to gfs fields at stratospheric levels"
     wei = 0.5
     do k = 1, nlev_atm_inp
     zzz = -8000.*alog(p(1,1,k)/1.e5)
     nsmooth = int((zzz/1000.)**2/45.)
     nsmooth = nsmooth/4

!       print*, "level", nlev_atm_inp-k+1, zzz, "nsmooth", nsmooth

     do kf = 1, nsmooth
     call filt2 (ge(1,1,k), iana, jana, 1.)
     call filt2 (te(1,1,k), iana, jana, 1.)
     call filt2 (ue(1,1,k), iana, jana, 1.)
     call filt2 (ve(1,1,k), iana, jana, 1.)
     call filt2 (qe(1,1,k), iana, jana, 1.)
     enddo
     enddo

   endif ! gfs

   gamma_inp = 0.
   iii = 0
   do j = 1, jana
   do i = 1, iana
   if (mask_frame(i,j)==1) then
   iii = iii + 1
   gamma_inp = gamma_inp + (te(i,j,nlev_atm_inp)-te(i,j,nlev_atm_inp-1)) &
                          /(ge(i,j,nlev_atm_inp-1)-ge(i,j,nlev_atm_inp))*g0
   endif
   enddo
   enddo
   gamma_inp = gamma_inp/float(iii)

!   print*, "Mean Lapse rate of input model near the surface", gamma_inp
!   print*, "Mean Lapse rate computed on", iii, "points"

 else  ! hybrid levels

! computation of p on half levels, with the expression:
! p(k+1/2) = a(k+1/2) + b(k+1/2)*ps

   do k=1,nlev_atm_inp
     k1=int(lev_list(k))
     phl(:,:,k)=ak(k1) + bk(k1)*exp(psil(:,:))
   enddo
   k1=int(lev_list(nlev_atm_inp))+1
   phl(:,:,nlev_atm_inp+1)=ak(k1) + bk(k1)*exp(psil(:,:))

! definition of p on full levels,
! with the expression: p(k)=0.5*(p(k+1/2)+p(k-1/2))

   do k=1,nlev_atm_inp
     p(:,:,k) = 0.5*(phl(:,:,k) + phl(:,:,k+1))
     pl(:,:,k) = alog(p(:,:,k))
   enddo

! gamma_inp is computed from second and fourth levels above the ground (to avoid extreme values
! near the surface)

   gamma_inp=0.
   k1=int(lev_list(nlev_atm_inp))
   k2=int(lev_list(nlev_atm_inp-1))
   k3=int(lev_list(nlev_atm_inp-2))
   k4=int(lev_list(nlev_atm_inp-3))
   iii=0
   do j=1,jana
   do i=1,iana
     if (mask_frame(i,j)==1) then
       iii=iii+1
       zphl1=ak(k1) + bk(k1)*exp(psel(i,j))
       zphl2=ak(k2) + bk(k2)*exp(psel(i,j))
       zphl3=ak(k3) + bk(k3)*exp(psel(i,j))
       zphl4=ak(k4) + bk(k4)*exp(psel(i,j))
       zpl1=alog(0.5*(zphl1 + zphl2))
       zpl2=alog(0.5*(zphl3 + zphl4))
       gamma_inp=gamma_inp+g0/rd*(alog(te(i,j,nlev_atm_inp-3))-alog(te(i,j,nlev_atm_inp-1)))/(zpl2-zpl1)
     endif
   enddo
   enddo
   gamma_inp=gamma_inp/float(iii)
!   print*, "Mean Lapse Rate of input model near the surface", gamma_inp
!   print*, "Mean Lapse rate computed on", iii, "points"

 endif

! comput. of (minus) spatially averaged istantaneous lapse rate of the input model near the surface

! gamma = max(0.35*gammac + 0.65*gamma_inp, 0.)

! print*, "Blended Lapse Rate", gamma

! ----------------------------------------------
! horizontal interpolation of the upper level fields
! ----------------------------------------------
! al: spline tension parameter

 al = 0.2
 if (level_type == 1) then ! isobaric levels
   do k=1,nlev_atm_inp
     call interp_spline_2d(ge(1,1,k),iana,jana,xe,ye,xxt,yyt,ntot,gp(1,1,k),al)
   enddo
 endif

 do k=1,nlev_atm_inp
   al = 0.2
   call interp_spline_2d(te(1,1,k),iana,jana,xe,ye,xxt,yyt,ntot,tp(1,1,k),al)
   call interp_spline_2d(ue(1,1,k),iana,jana,xe,ye,xxu,yyu,ntot,upu(1,1,k),al)
   call interp_spline_2d(ue(1,1,k),iana,jana,xe,ye,xxv,yyv,ntot,upv(1,1,k),al)
   call interp_spline_2d(ve(1,1,k),iana,jana,xe,ye,xxu,yyu,ntot,vpu(1,1,k),al)
   call interp_spline_2d(ve(1,1,k),iana,jana,xe,ye,xxv,yyv,ntot,vpv(1,1,k),al)
   al = 0.4
   call interp_spline_2d(qe(1,1,k),iana,jana,xe,ye,xxt,yyt,ntot,qp(1,1,k),al)
 enddo

! elimination of possible negative values of specific humidity

   qp =max (1.e-9, qp)

   upup = upu
   vpvp = vpv

! computation of geopotential at the ground surface

   htopi = phig/g0

!----------------------------------
! start of vertical interpolation
!----------------------------------

! definition of tv (virtual temperature) and ps (surface pressure)

 if (level_type==1) then ! isobaric levels

   i0=1
   j0=1
   do j=1,nlat
   do i=1,nlon

! Definition of virtual t at the auxiliary pressure levels

! The layer average temperature is imposed to be consistent with the hydrostatic eq.

     tv(i,j,1)=(1+ep*qp(i,j,1))*tp(i,j,1)
     gaux(1)=gp(i,j,1)
     do k=1,nlev_atm_inp-1
       k1=2*k+1
       tv(i,j,k1)=(1.+ep*qp(i,j,k+1))*tp(i,j,k+1)
       gaux(k1)=gp(i,j,k+1)
       k2=2*k

! Alternative expression for tvm: after tests on gfs data compar. 26 and 47 levels (apr. 2010),
! a definition of delp depending on p is introduced, with a transition from dp/p at low levels
! to d(ln(p)) at higher levels - this empirically reduces "noise" al low levels but maintains
! realistic corrections at higher levels

       peso=0.5*(p(i0,j0,k)+p(i0,j0,k+1))/1.1e5          ! weight for the denomin. (p/prif)
       delp1=pl(i0,j0,k)-pl(i0,j0,k+1)                   ! delta(log(p))
       delp2=2.*(p(i0,j0,k)-p(i0,j0,k+1))/(p(i0,j0,k)+p(i0,j0,k+1)) ! deltap/pmed
       delp3=peso*delp2+(1.-peso)*delp1
       tvm=-(gp(i,j,k)-gp(i,j,k+1))/delp3/rd

! Below: blending depending again on p, since empirically it seems that the correction based
! on the use of gph data is more realistic at higher levels than at lower

       tv1=1.6*tvm-0.3*(tv(i,j,k2-1)+tv(i,j,k2+1))
       tv2=0.5*(tv(i,j,k2-1)+tv(i,j,k2+1))
       peso2=sqrt(peso)
       tv(i,j,k2)=peso2*tv2+(1.-peso2)*tv1
       gaux(k2)=gaux(k1)-rd*(plaux(k2)-plaux(k1))*(tv(i,j,k2)+tv(i,j,k1))/2.
     enddo

! Computation of surface pressure using t virt. at the auxiliary levels computed above.
! the hydrostatic relation is used to extrapolate (or interpolate) tvirt using
! geopotential as independent variable (instead of log(p) to avoid a second order equation).
! The geopotential is computed first at the auxiliary intermediate levels
! (consistent with the comput. of tvirt at the same levels, using thickness between standard levels).

     ki=0
     do k=1,nauxp
       if (phig(i,j) < gaux(k)) ki=k
     enddo
     if (ki == nauxp)  ki=nauxp-1

     if (ki == 0) then
       print *,'ki=0, error in computing geopotential!'
       print *,'i = ',i,'  j = ',j
       print *,'phig(i,j)=',phig(i,j),' gaux(1)=',gaux(1)
       stop
     endif

! ts (virt.) extrapolated from air temp. (def. in tsair)
! t is kept constant in case of t inversions, when the model surface
! located below the lowest analysis level)

     tvs=tv(i,j,ki+1)+(tv(i,j,ki+1)-tv(i,j,ki))*(phig(i,j)-gaux(ki+1))/(gaux(ki+1)-gaux(ki))
     if (tv(i,j,ki+1) < tv(i,j,ki).and.(phig(i,j)-gaux(ki+1)) < 0.) tvs=tv(i,j,ki+1)
     tsair(i,j)=tvs

! Definition of surface pressure - pst contains the logarithm of ps on t points

     pst(i,j) = plaux(ki+1)-(2./rd)*(phig(i,j)-gaux(ki+1))/(tvs+tv(i,j,ki+1))
     ps(i,j)  = exp(pst(i,j))

! end of loop on the grid points

   enddo
   enddo

 else  ! from here case of hybrid levels

! Definition of geopotential at input model levels, not available in input data.
! it is computed using phigi and pressure data.
! geopotential is used to compute the surface pressure only (not for temperature definition).
! fields, which are resuls of horizontal interpolation, are used.

! Definition of virtual temperature

   do k=1,nlev_atm_inp
     tv(:,:,k)=(1.+ep*qp(:,:,k))*tp(:,:,k)
   enddo

! Definition of geopotential at internal levels (gp)

   gp(:,:,nlev_atm_inp) = phigi(:,:)-rd*tv(:,:,nlev_atm_inp)*(pl(:,:,nlev_atm_inp)-psil(:,:))

   kf=nlev_atm_inp-1
   do k=kf,1,-1
     gp(:,:,k) = gp(:,:,k+1)-rd*0.5*(tv(:,:,k)+tv(:,:,k+1))*(pl(:,:,k)-pl(:,:,k+1))
   enddo

   do j=1,nlat
   do i=1,nlon

! definition of surface pressure ps
! search of the index of the first level above ground

     ki=0
     do k=1,nlev_atm_inp
       if (phig(i,j) < gp(i,j,k)) ki=k
     enddo

     if (ki == 0) then
       print *,'ki=0: error in geopotential computation! stop.'
       print *,'i = ',i,'  j = ',j
       print *,'phig(i,j)=',phig(i,j),' gp(i,j,1)=',gp(i,j,1)
       stop
     endif

! ts is extrapolated from air temperature (in tsair).
! if the grid point is under the orography, ts is computed as a weighted average
! between the result of the extrapolation with lapse rate gamma and
! the result of the linear extrapolation from higher levels.

     if (ki == nlev_atm_inp) then

! In this case the ground surface is below the lowest level of the input model:
! in cases of inversion extrapol. is problematic, therefore only the lapse rate gamma
! is used; otherwise a weighted average of the lapse rate gamma
! and linear extrapolation from levels above is applied

       if ((tv(i,j,ki)-tv(i,j,ki-1)) > 0.) then
         tsair(i,j)=tv(i,j,ki)+0.7*(gamma(i,j)*(gp(i,j,ki)-phig(i,j))/g0) &
         +0.3*((tv(i,j,ki)-tv(i,j,ki-1))*(phig(i,j)-gp(i,j,ki))/(gp(i,j,ki)-gp(i,j,ki-1)))
       else
         tsair(i,j)=tv(i,j,ki)+gamma(i,j)*(gp(i,j,ki)-phig(i,j))/g0
       endif

! Extrapolation of surface pressure ps at t points
! pst contains the logarithm of ps on t points

       pst(i,j) = pl(i,j,ki)+(gp(i,j,ki)-phig(i,j))/(rd*0.5*(tv(i,j,ki)+tsair(i,j)))
       ps (i,j) = exp(pst(i,j))

     else

! In this case the ground surface is above the lowest level of the input model:
! interpolation is applied

       tsair(i,j)=tv(i,j,ki+1)+(tv(i,j,ki+1)-tv(i,j,ki))*(phig(i,j)-gp(i,j,ki+1))/(gp(i,j,ki+1)-gp(i,j,ki))

! Extrapolation of surface pressure ps at t points
! pst contains the logarithm of ps on t points

       pst(i,j) = pl(i,j,ki+1)+(gp(i,j,ki+1)-phig(i,j))/(rd*0.5*(tsair(i,j)+tv(i,j,ki+1)))
       ps (i,j) = exp(pst(i,j))

     endif

   enddo
   enddo

 endif ! end of distinction between level types

! Definition of pressure at integer levels

 do k = 1, nlev
 do j = 1, nlat
 do i = 1, nlon
   press(i,j,k) = aka(k) + ps(i,j)*bika(k)
 enddo
 enddo
 enddo

   psmin  = minval(ps(:,:))
   print*
   print*, "psmin =", psmin, "alfa =", alfa, "Pzer =", pzer, "Vert. coord. parameter =",vert_coord_c
!print *,minloc(ps(:,:))
!stop

! definition of pressure at the ground surface on u, v points

   do j = 1, nlat
   do i = 1, nlon-1
   psu(i,j) = (pst(i,j)+pst(i+1,j))/2.
   enddo
   psu(nlon,j) = psu(1,j)
   enddo

   do i = 1, nlon
   psv(i,1) = pst(i,2)  ! non importa se e' sbagliato, tanto e' fuori dall'universo
   do j = 2, nlat
   psv(i,j) = (pst(i,j)+pst(i,j-1))/2.
   enddo
   enddo

! ----------------------------------------------------------
! Vertical interpolation of fields at the atmospheric levels
! Start of the main loop on grid points
! ----------------------------------------------------------

   do j = 1, nlat
   do i = 1, nlon

! Computation of pressure at the hybrid levels on the t points (for t and q)

   call levels (pst(i,j), nlev, plsig, aka, bika)

! Definition of auxiliary input vectors (qpk contains sqrt(q); tvk contains virtual t)

   do k = 1, nlev_atm_inp
   upk(k) = upup(i,j,k)
   vpk(k) = vpvp(i,j,k)
   qpk(k) = sqrt(qp(i,j,k))
   enddo

   if (level_type==1) then  ! isobaric levels

     i0=1
     j0=1
     pl2(1:nlev_atm_inp)=pl(i0,j0,1:nlev_atm_inp)
     tvk(1:nauxp)=tv(i,j,1:nauxp)

! Interpolation of tv
! tv is kept constant in case of extrapolation toward the ground surface with t inversion

     call near(plsig,nlev,plaux,nauxp,iv)
     alf=.65
     ex1=0.5
     ex2=0.6
     call interp_spline_1d(tk,plsig,nlev,tvk,plaux,nauxp,iv,alf,ex1,ex2)
     do k=1,nlev
       if(plsig(k) > plaux(nauxp).and.tvk(nauxp) < tvk(nauxp-1)) tk(k)=tvk(nauxp)
     enddo

   else   ! hybrid levels

     do k=1,nlev_atm_inp
       tvk(k)=tp(i,j,k)
       pl2(k)=pl(i,j,k)
     enddo

! Interpolation of t
! t is kept constant in case of extrapolation toward the ground surface with t inversion
! (in this case t is recomputed below)

     call near(plsig,nlev,pl2,nlev_atm_inp,iv)
     alf=.65
     ex1=0.5
     ex2=0.6
     call interp_spline_1d(tk,plsig,nlev,tvk,pl2,nlev_atm_inp,iv,alf,ex1,ex2)
     do k=1,nlev
       if (plsig(k) > pl2(nlev_atm_inp).and.tvk(nlev_atm_inp) < tvk(nlev_atm_inp-1)) tk(k)=tvk(nlev_atm_inp)
     enddo

   endif

! From here: for all level types
! Interpolation of sqrt of specific humidity

   call near (plsig,nlev,pl2,nlev_atm_inp,iv)
   alf=.65
   ex1=.3
   ex2=.3
   call interp_spline_1d(qk,plsig,nlev,qpk,pl2,nlev_atm_inp,iv,alf,ex1,ex2)

! q recomputation

   q(i,j,1:nlev)=max(1.e-7,qk(1:nlev)**2)

   if (level_type == 1) then  ! isobaric levels

! Conversion of t, that in case of isob. levels contains virt. temp., into temperature

     t(i,j,1:nlev)=tk(1:nlev)/(1+ep*q(i,j,1:nlev))

   else  ! hybrid levels

! Redefinition of t using prescribed lapse rate in the case when
! a level is located under the lowest level of the input model

     ki=0
     do k=1,nlev
       if (plsig(k) <= pl(i,j,nlev_atm_inp)) ki=k
     enddo

     if (ki < nlev) then

! Temperature extrapolation: weighted average between the result of vert. interp.
! and extrapolation (using the constant lapse rate gamma);
! in case of inversion, the constant lapse rate is used only

       do kt=ki+1,nlev
       gamma1=gamma(i,j)*(0.5*rd/g0)*(pl(i,j,nlev_atm_inp)-plsig(kt))
         if ((tp(i,j,nlev_atm_inp-1)-tp(i,j,nlev_atm_inp)) > 0.) then
         tk(kt)= tp(i,j,nlev_atm_inp)*(1.-gamma1)/(1.+gamma1)
         else
         tk(kt)=.5*(tk(kt)+tp(i,j,nlev_atm_inp)*(1.-gamma1)/(1.+gamma1))
         endif
       enddo

     endif

     t(i,j,1:nlev)=tk(1:nlev)

   endif

! From here: for all level types

! Controlq is called only in the case in which cloud water/ice is absent in input;
! it checks that q does not exceed (significantly) saturation and
! increases slightly q values close to saturation

    do k = 1, nlev
    psig = exp(plsig(k))
    if (iflag_cloude == 0) call controlq (q(i,j,k), t(i,j,k), psig, eps)
    enddo

    if (level_type == 2) then  ! hybrid levels only

! Reset of specific humidity using the condition of constant relative humidity
! at downward extrapolated levels;
! but if lapse rate gamma is large, extrapolation with constant q is applied

    if (ki < nlev) then

       call qsat_tetens(tp(i,j,nlev_atm_inp), exp(pl(i,j,nlev_atm_inp)), eps, zqs, zqsw, zqsi, zes, zesw, zesi)
       if (tp(i,j,nlev_atm_inp) >= t0) then
       qsat=zqsw
       else
       qsat=zqsi
       endif
       wsat = qsat/(1.-qsat)
       w = qp(i,j,nlev_atm_inp)/(1.-qp(i,j,nlev_atm_inp))
       rh = min(w/wsat, 0.995)

       do kt=ki+1,nlev
         call qsat_tetens(t(i,j,kt), exp(plsig(kt)), eps, zqs, zqsw, zqsi, zes, zesw, zesi)
         if (t(i,j,kt) >= t0) then
           qsat=zqsw
         else
           qsat=zqsi
         endif
         wsat = qsat/(1.-qsat)
         w =  rh*wsat
         q(i,j,kt) = w/(1.+w)
         if (gamma(i,j) > 7.e-3) q(i,j,kt)=qp(i,j,nlev_atm_inp)
       enddo

     endif

   endif

! From here: for all level types

! Interpolation of wind components;
! wind is kept almost constant at a level under nlev_atm_inp

   call levels (psu(i,j), nlev, plsig, aka, bika)
   call near (plsig, nlev, pl2,nlev_atm_inp, iv)
   alf = .65
   ex1 = .2
   ex2 = .2
   call interp_spline_1d (uk, plsig, nlev, upk, pl2, nlev_atm_inp, iv, alf, ex1, ex2)

   call levels (psv(i,j), nlev, plsig, aka, bika)
   call near (plsig, nlev, pl2, nlev_atm_inp, iv)
   call interp_spline_1d (vk, plsig, nlev, vpk, pl2, nlev_atm_inp, iv, alf, ex1, ex2)

   u(i,j,1:nlev) = uk(1:nlev)
   v(i,j,1:nlev) = vk(1:nlev)

! End of vertical interpolations (loop on grid points)

   enddo
   enddo

! Comput. of max. wind speed (only for printout)

   wspeed = sqrt(u**2+v**2)
   wm  = maxval(wspeed)
   ijk = maxloc(wspeed)
   write(*,'(a,f8.2,a,3i5)') " Max. wind speed:", wm," at grid point",ijk(1),ijk(2),ijk(3)

!------------------------------------------------------------
! Definition of fields at the ground surface and in the soil
!------------------------------------------------------------

   zzh1 = 20. + 70.*(dlat/0.12)**.5
   zzh2 = 2.*zzh1

!  Start of the loop on the grid points

   do j = 1, nlat
   do i = 1, nlon

! Computation of a weighting factor depending on the difference between the
! input model topography interpolated on the model grid (phigi2d) and topography
! defined in output grid (htopi). diftop > 0 if  phigi2d > htopi.
! in case of "valley" with respect to the input model orography,
! and in winter only, the lapse rate is further reduced

   diftop(i,j) = (phigi2d(i,j)/g0)-htopi(i,j)

! Definition of twater: temperature of surface water (sea, lakes, rivers).
! the value is selected between ssti (that is correct for open sea, but not for lakes)
! and a combination of temp. at levels 2 and 3 (that is better for lakes)

   ztlake = .5*tgi(i,j,2) + .5*tgi(i,j,3) + gamma(i,j)*diftop(i,j)

   if (htopi(i,j) < zzh1) then
     twater(i,j) = ssti(i,j)
   elseif (htopi(i,j) > zzh2) then
     twater(i,j) = ztlake
   else
     twater(i,j) = ssti(i,j)+(htopi(i,j)-zzh1)*(ztlake-ssti(i,j))/(zzh2-zzh1)
   endif

! Definition of the surface temperature.

   tsair(i,j) = tsair(i,j)/(1.+ep*q(i,j,nlev))

! Definition of tg
! for definition of orographic corrections of soil temperature at various
! levels, different lapse rates are applied, with seasonal variations decreasing
! from near the surface to deeper levels.

   if (fmask(i,j) >= .5) then
   tg(i,j,:) = twater(i,j)*(1.-fice(i,j)) + tice(i,j)*fice(i,j)
   else
   tg(i,j,1) = tgi(i,j,1) + gamma(i,j)*diftop(i,j)
   tg(i,j,2) = tgi(i,j,2) + gamma(i,j)*diftop(i,j)
   tg(i,j,3) = tgi(i,j,3) + 0.5*(gamma(i,j)+6.e-3)*diftop(i,j)
   tg(i,j,4:nlevg) = tgi(i,j,4:nlevg) + 6.e-3*diftop(i,j)
   endif

! Re-definition of FICE

   if (twater(i,j)<271.4) then
   fice(i,j) = max (fice(i,j), min (1., (271.4-twater(i,j))/10. ))
   fice(i,j) = min (fice(i,j), fmask(i,j))
   endif

! Definition of iceth as a function of fice

   iceth(i,j) = -3.5 + 5.*fice(i,j)  ! sets 0.5<iceth<1.5 for 0.8<fice<1.
   if (fice(i,j).lt.0.8) then
   iceth(i,j) = 0.                   ! treshold must be consistent with the BOLAM model
   else
   iceth(i,j) = max(iceth(i,j), 0.5) ! treshold and value must be consistent with the BOLAM model
   endif
   iceth(i,j) = min(iceth(i,j), 3.0)

! Case of thick ice: all soil TG defined as in soil/glacier scheme

   if (fice(i,j).gt.0.8) then
   if (htopi(i,j).lt.170.) then 
   do k = 2, nlevg
   tg(i,j,k) = tg(i,j,1) + (zsoil(k)-zsoil(1))*(271.4-tg(i,j,1))/(zsoil(nlevg)-zsoil(1))
   enddo
   else
   do k = 2, nlevg
   tg(i,j,k) = tg(i,j,1) + (zsoil(k)-zsoil(1))*(273.1-tg(i,j,1))/(zsoil(nlevg)-zsoil(1))
   enddo
   endif
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
        soil_albedo_dry(i,j) = 0.7
        soil_albedo_wet(i,j) = 0.7
        albedo(i,j) = 0.7
        tskin(i,j) = min(tskin(i,j), 273.0)
        fice_soil_surf(i,j) = 1.
        do jklev = 1, ind_lev_soil_h_bottom(i,j)
          tg(i,j,jklev) = min(tg(i,j,jklev), 273.0)
          fice_soil(i,j,jklev) = 1.
        enddo
      else ! land glacier
        tskin(i,j) = min(tskin(i,j), 273.1)
        tg(i,j,1) = min(tg(i,j,1), 272.5)
        tg(i,j,2) = min(tg(i,j,2), 272.0)
        tg(i,j,3) = min(tg(i,j,3), 271.5)
        tg(i,j,4:nlevg) = min(tg(i,j,4:nlevg), 271.0)
      endif

    endif

! For the case of gfs (noaa-ncep) input data: decrease of temperature values at low levels 
! and in the soil (less at deepest levels) over the Po valley
! the max. correction is at the end of dec. - correction is null from apr. to sept.

   if (model_input=='gfs') then
   zlon = xxt(i,j)
   zlat = yyt(i,j)
    if (zlat>44.3 .and. zlat<46.2 .and. zlon>7. .and. zlon<13.4 .and. fmask(i,j)<.35 .and. htopi(i,j)<500.) then
    zdelt = -1.5*max (0., cos(ndayr*2.*pi/365.))
    tg(i,j,1)  = tg(i,j,1)  + zdelt
    tg(i,j,2)  = tg(i,j,2)  + zdelt
    tg(i,j,3)  = tg(i,j,3)  + zdelt*.7
    tg(i,j,4:nlevg)  = tg(i,j,4:nlevg)  + zdelt*.45
     do k = nlev/2, nlev
     if (press(i,j,k) > 945.e2 .and. press(i,j,k) <= 960.e2) t(i,j,k) = t(i,j,k) + 0.5*zdelt
     if (press(i,j,k) > 960.e2) t(i,j,k) = t(i,j,k) + zdelt
     enddo
    endif
   endif

! Snow definition (in meters of equivalent water)
! Further snow reduction over land as a function of tg1 and tsair

   snow(i,j) = max (snowi(i,j), 0.)
   if (tg(i,j,1)>277.5)  snow(i,j) = 0.
   if (snow(i,j)<0.0005) snow(i,j) = 0.
   ttme = .8*tg(i,j,1) + .2*tsair(i,j)
   if (ttme>278.) snow(i,j) = 0.
   snow(i,j) = min (snowi(i,j), 0.111)

! Defintion of soil temperature at the soil top

    if (snow(i,j) < 0.01) then
      tgsurf(i,j)=tskin(i,j)
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

! Conversion of relative into volumetric water content (m**3/m**3)

   qgmin(i,j) = soil_qmin(i,j,3)
   qg(i,j,1:nlevg) = qgi(i,j,1:nlevg)*(soil_qmax(i,j,1:nlevg)-soil_qmin(i,j,1:nlevg)) + soil_qmin(i,j,1:nlevg)

! Soil water content definition for water bodies: sea, lakes, rivers, glaciers

     if (fmask(i,j)>.5 .or. int(soil_map(i,j,1)) == nst-1) qg(i,j,:) = 1.

   qgsurf(i,j) = qg(i,j,1)

! Definition of roughness length

   rgm(i,j) = max( (1.-fmask(i,j))*veg_roughness(i,j) + 1.e-3*(fmask(i,j)-fice(i,j)) + 5.e-2*fice(i,j), 1.e-3)
   rgq(i,j) = rgm(i,j)

! Correction of wind at the lowest level as a function of roughness

   zrgm = min (rgm(i,j), 8.)
   u(i,j,nlev) = u(i,j,nlev)/(1.+.3*zrgm)
   v(i,j,nlev) = v(i,j,nlev)/(1.+.3*zrgm)

! Denition of snow variables at snow levels
! 1-st level is snow top, bottom level is soil surface

    snow_lev(i,j,:) = val_missing
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
      snow_t(i,j,1) = min(tskin(i,j), 273.)
      if (nlsnow < 2) then
        snow_t(i,j,1) = min(tg(i,j,1), 273.)
      else
        z1=1./float(nlsnow+1)
        snow_t(i,j,nlsnow+1) = min(tg(i,j,1), 273.)
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

   enddo    ! End of the loop on grid points
   enddo    ! End of the loop on grid points

   rgm = max (rgm, 1.e-3)
   rgq = max (rgq, 1.e-3)

      do k = 1, nlev
      do j = 1, nlat
      do i = 1, nlon
    if (press(i,j,k) >= 500.) then
      zt0t = t0/t(i,j,k)
      if (zt0t.le.1.) then
      zesk = e0*exp(-ccw1*alog(zt0t)+ccw2*(1.-zt0t)) !! Partial pressure over water
      else
      zesk = e0*exp(-cci1*alog(zt0t)+cci2*(1.-zt0t)) !! Partial pressure over ice
      endif

      zqsat = zesk*eps/(press(i,j,k)+zesk*(eps-1.))
      zqtrue = min (zqsat, q(i,j,k))
      zdq = max (q(i,j,k)-zqtrue, 0.)
      qc(i,j,k) = zdq
      q(i,j,k)  = zqtrue
    else
      q(i,j,k)  = 1.e-8
      qc(i,j,k) = 0.
    endif
if (i==330.and.j==199) print *,'poisk  qc ',k,press(i,j,k),t(i,j,k)-273.15,q(i,j,k),qc(i,j,k),zesk
!if (i==330.and.j==199.and.k==1) print *,'poisk 1 qc ',qc(i,j,k),q(i,j,k),zqsat,press(i,j,k),t(i,j,k),zesk
!if (k==1.and.qc(i,j,k) > 1.e-8) print *,'poisk 2 qc ',i,j,qc(i,j,k)
      enddo
      enddo
      enddo

     qvsurf (:,:) = q(:,:,nlev)

!--------------------------------------------------------------------------------------------

! Using of forecast soil/snow data in place of analysis data 

  call read_forecast_mhf_soil(iflag_soil_snow_frc)
!  call read_forecast_mhf_soil_access_direct(iflag_soil_snow_frc)
!  call read_forecast_mhf_ascii_soil(iflag_soil_snow_frc)

  write (*,*)
  if (iflag_soil_snow_frc == 1) then

    do j=1,nlat
    do i=1,nlon

      tskin_copy(i,j)=tskin(i,j)
      tg_copy(i,j,1:nlevg)=tg(i,j,1:nlevg)

! poisk
!      if (fmask(i,j) < 0.5 ) then ! land or temperary sea glacier
      if ( (fmask(i,j) < 0.5).or.(fmask(i,j) >= 0.5.and.fice(i,j) >= 0.8.and.fice_frc(i,j) >= 0.8) ) then ! land or temperary sea glacier

        tskin(i,j) = tskin_frc(i,j)
        tgsurf(i,j) = tgsurf_frc(i,j)
        tg(i,j,1:nlevg) = tg_frc(i,j,1:nlevg)
        qvsurf(i,j) = qvsurf_frc(i,j)
        qgsurf(i,j) = qgsurf_frc(i,j)
        qg(i,j,1:nlevg) = qg_frc(i,j,1:nlevg)
! poisk 05.07.2021 --->
!  if (yytg(i,j) > 10..and.yytg(i,j) < 32.) then
!    do k=1,5
!      qg(i,j,k)=min(qg(i,j,k)*1.5, soil_qmax(i,j,k)*0.5)
!    enddo
!  endif
! <---
        fice_soil_surf(i,j) = fice_soil_surf_frc(i,j)
        fice_soil(i,j,1:nlevg) = fice_soil_frc(i,j,1:nlevg)
!if (fmask(i,j)<0.01.and.maxval(fice_soil(i,j,1:6))>0.8.and.minval(fice_soil(i,j,1:6))<0.01) then
!  write (31,*) i,j,xxtg(i,j),yytg(i,j),tg(i,j,1:nlevg),fice_soil(i,j,1:nlevg)
!endif
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
          soil_albedo_dry(i,j) = 0.7
          soil_albedo_wet(i,j) = 0.7

          do jklev = 1, nlevg
            if (soil_lev(jklev) > iceth(i,j)) exit
          enddo
          ind_lev_soil_h_bottom(i,j) = max(jklev-1, 3)
          ind_lev_soil_w_bottom(i,j) = 0

          tskin(i,j) = min(tskin(i,j), 273.0)
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
            snow_t(i,j,1) = min(tskin(i,j), 273.)
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

          soil_psi(i,j,:) = 0. ! thermal conductivity of ice
          soil_albedo_dry(i,j) = 0.07
          soil_albedo_wet(i,j) = 0.07

          ind_lev_soil_h_bottom(i,j) = 0
          ind_lev_soil_w_bottom(i,j) = 0

          tskin(i,j) = tskin_copy(i,j)
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
! Write mhf
!--------------------------------------------------------------------------------------------

   if (flag_constant_fields == 1) then
     call write_param_const(x0, y0, alon0, alat0, htopi, zhtopvar)
   endif
 
   call wrmhf_atm
   call wrmhf_soil

   print*
   print *,"Globo initial condition written"
   print*

!   call plotout (rgm, nlon, nlat, 99)
!   call plotout (rgq, nlon, nlat, 99)
!   call plotout (fice, nlon, nlat, 99)
!   call plotout (fmask, nlon, nlat, 99)
!   call plotout (snow, nlon, nlat, 99)
!   call plotout (tg(:,:,1), nlon, nlat, 99)
!   call plotout (tg(:,:,nlevg), nlon, nlat, 99)
!   call plotout (qg(:,:,1), nlon, nlat, 99)
!   call plotout (qgmin, nlon, nlat, 99)
!   call plotout (qg(:,:,3), nlon, nlat, 99)
!   call plotout (qg(:,:,5), nlon, nlat, 99)
!   call plotout (qg(:,:,7), nlon, nlat, 99)

   stop

50 print *,"input file not available"

   stop
   end program preglobo
!##################################################################################################
subroutine wrec2 (kunit, nlon, nlat, vect)

implicit none

integer :: kunit, nlon, nlat
real, dimension(nlon,nlat) :: vect

 write(kunit) vect(1:nlon,1:nlat)

return
end subroutine wrec2
!##################################################################################################
subroutine wrec2_int (kunit, nlon, nlat, ivect)

implicit none

integer :: kunit, nlon, nlat
integer, dimension(nlon,nlat) :: ivect

 write(kunit) ivect(1:nlon,1:nlat)

return
end subroutine wrec2_int
!##################################################################################################
subroutine rrec2 (kunit, nlon, nlat, vect)

implicit none

integer :: kunit, nlon, nlat
real, dimension(nlon,nlat) :: vect

 read(kunit) vect(1:nlon,1:nlat)

return
end subroutine rrec2
!##################################################################################################
subroutine rrec2_int (kunit, nlon, nlat, ivect)

implicit none

integer :: kunit, nlon, nlat
integer, dimension(nlon,nlat) :: ivect

 read(kunit) ivect(1:nlon,1:nlat)

return
end subroutine rrec2_int
!##################################################################################################
subroutine wrmhf_atm

! Writes a MHF file of BOLAM with atmospheric variables

use param
use mod_fields

implicit none

real, dimension(nlon,nlat) :: field2d_add
integer :: nf, iunit = 11, i, j, k
character(len=30) :: file_output='input_atm_01.mhf'

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

return
end
!##################################################################################################
subroutine wrmhf_soil

! Writes a MHF file of BOLAM with surface, sea, soil variables

use param
use mod_fields

implicit none

real, dimension(nlon,nlat) :: field2d_add
integer :: iunit = 11, i, j, k, iwr
character(len=30) :: file_output='input_soil_01.mhf'

 field2d_add(:,:) = 0.

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

! cloud, totpre, conpre, snfall

 do iwr = 1, 4
   call wrec2 (iunit, nlon, nlat, field2d_add(1:nlon,1:nlat))
 enddo

 call wrec2 (iunit, nlon, nlat, tskin(1:nlon,1:nlat))

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

! cswfl, clwfl, chflux, cqflux, t2min, t2max, runoff, runoff_tot

 do iwr = 1, 8
   call wrec2 (iunit, nlon, nlat, field2d_add(1:nlon,1:nlat))
 enddo

! Writing of additional 2D fields

 field2d_add(:,:) = 0.
 do iwr = 1,3
   call wrec2 (iunit, nlon, nlat, field2d_add(1:nlon,1:nlat))
 enddo

 close (iunit)

 write (*,*)
 write (*,*) 'Output mhf file ',trim(file_output),' written' 
 write (*,*)
 
return
end
!##################################################################################################
   subroutine levels (psl, nlev, plsig, aka, bika)

! given the suface pressure psl=logps, computes the logp on hybrid levels

   real, dimension(nlev) :: plsig, aka, bika

   do k = 1, nlev
     plsig(k) = log ( aka(k) + exp(psl)*bika(k) )
   enddo

   return
   end
!##################################################################################################
   subroutine conv_ifs_data

! procedure to convert meteorological fields derived from ecmwf-ifs

use param, only : iana, iana0, jana, nlev_atm_inp_max, nlev_soil_inp_max, nlev_atm_inp, nlev_soil_inp, &
 val_missing, pi, idate0, iperiod, iperiod_inp, surf_elaborate
use gribana

implicit none

integer :: ind_field, i, j, k, iii, jjj, ifl, np, iste, npt, soil_flag
real :: zqcmin, qcrit, zzz, zflatt, day1, day2
real, dimension(iana,jana,nlev_atm_inp_max) :: rh
integer, dimension(iana,jana) :: imaske
real, dimension(iana,jana,nlev_soil_inp_max) :: worki, workf
real, dimension(nlev_soil_inp_max) :: zdtg
real, save :: day
integer, save :: iday, nday, nste
integer, dimension(12) :: imon=(/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
real :: depth, depthsum
real :: tgav1,tgav2,tgav3,tgav4,qgav1,qgav2,qgav3,qgav4

 mask_frame(:,:)=1

 if (iperiod_inp(1)==1) then         ! time unit is hour
   iperiod(1)=iperiod_inp(2)/24                             ! days
   iperiod(2)=iperiod_inp(2)-iperiod(1)*24                  ! hours
   iperiod(3)=0                                             ! minutes
 elseif (iperiod_inp(1)==0) then     ! time unit is minute
   iperiod(1)=iperiod_inp(2)/24/60                          ! days
   iperiod(2)=(iperiod_inp(2)-iperiod(1)*24*60)/60          ! hours
   iperiod(3)=iperiod_inp(2)-iperiod(1)*24*60-iperiod(2)*60 ! minutes
 endif

! definition of meteorological fields declared in main program
! list to be checked ! see subroutine read_ecmwf_grib_data

 ge(1:iana0,1:jana,1:nlev_atm_inp) = field3d(1:iana0,1:jana,1:nlev_atm_inp,1) ! geopotential in m**2 s**-2
 te(1:iana0,1:jana,1:nlev_atm_inp) = field3d(1:iana0,1:jana,1:nlev_atm_inp,2)
 ue(1:iana0,1:jana,1:nlev_atm_inp) = field3d(1:iana0,1:jana,1:nlev_atm_inp,3)
 ve(1:iana0,1:jana,1:nlev_atm_inp) = field3d(1:iana0,1:jana,1:nlev_atm_inp,4)
 qe(1:iana0,1:jana,1:nlev_atm_inp) = field3d(1:iana0,1:jana,1:nlev_atm_inp,5)
 rh(1:iana0,1:jana,1:nlev_atm_inp) = field3d(1:iana0,1:jana,1:nlev_atm_inp,6)

! note that in the following the last index denotes:
! 7: total cloud water content (sum of water and ice - parameter normally not available from ecmwf-ifs);
! 8: cloud liquid water; 9: cloud ice.

 zqcmin=1.e-7
 iflag_cloude=0
 qce (:,:,:) = 0.
 qcwe(:,:,:) = 0.
 qcie(:,:,:) = 0.

! check availability of total cloud water (liquid + ice) in input data

 if (any(int(field3d(:,:,1:nlev_atm_inp,7)) == int(val_missing)).or.maxval(field3d(:,:,1:nlev_atm_inp,7)) < zqcmin) then
 print*, "no total cloud water (liquid + ice) available in input"
 iflag_cloude=0
 else
 iflag_cloude=1
 print*, "total cloud water (liquid + ice) available in input"
 qce (1:iana0,1:jana,1:nlev_atm_inp) = field3d(1:iana0,1:jana,1:nlev_atm_inp,7)
 endif

! check availability of cloud liquid water in input data

 if (any(int(field3d(:,:,1:nlev_atm_inp,8)) == int(val_missing)).or.maxval(field3d(:,:,1:nlev_atm_inp,8)) < zqcmin) then
 print*, "no cloud liquid water available in input"
 iflag_cloude=max(iflag_cloude, 0)
 else
 iflag_cloude=1
 print*, "cloud liquid water available in input"
 qcwe(1:iana0,1:jana,1:nlev_atm_inp) = field3d(1:iana0,1:jana,1:nlev_atm_inp,8)
 endif

! check availability of cloud ice in input data

 if (any(int(field3d(:,:,1:nlev_atm_inp,9)) == int(val_missing)).or.maxval(field3d(:,:,1:nlev_atm_inp,9)) < zqcmin) then
 print*, "no cloud ice available in input"
 iflag_cloude=max(iflag_cloude, 0)
 else
 iflag_cloude=1
 print*, "cloud ice available in input"
 qcie(1:iana0,1:jana,1:nlev_atm_inp) = field3d(1:iana0,1:jana,1:nlev_atm_inp,9)
 endif

   if (iflag_cloude==1) qce(:,:,1:nlev_atm_inp) = max (qce(:,:,1:nlev_atm_inp), &
                        qcwe(:,:,1:nlev_atm_inp) + qcie(:,:,1:nlev_atm_inp))

   ge  (iana, 1:jana, 1:nlev_atm_inp) = ge  (1, 1:jana, 1:nlev_atm_inp)
   te  (iana, 1:jana, 1:nlev_atm_inp) = te  (1, 1:jana, 1:nlev_atm_inp)
   ue  (iana, 1:jana, 1:nlev_atm_inp) = ue  (1, 1:jana, 1:nlev_atm_inp)
   ve  (iana, 1:jana, 1:nlev_atm_inp) = ve  (1, 1:jana, 1:nlev_atm_inp)
   qe  (iana, 1:jana, 1:nlev_atm_inp) = qe  (1, 1:jana, 1:nlev_atm_inp)
   rh  (iana, 1:jana, 1:nlev_atm_inp) = rh  (1, 1:jana, 1:nlev_atm_inp)
   qce (iana, 1:jana, 1:nlev_atm_inp) = qce (1, 1:jana, 1:nlev_atm_inp)
   qcwe(iana, 1:jana, 1:nlev_atm_inp) = qcwe(1, 1:jana, 1:nlev_atm_inp)
   qcie(iana, 1:jana, 1:nlev_atm_inp) = qcie(1, 1:jana, 1:nlev_atm_inp)

! if available, the 2-d geopotential at the ground is set, otherwise the
! "3-d" geop. (model level) at the ground is used

 if (any(int(field2d(:,:,1))==int(val_missing))) then
   phige2d(1:iana0,1:jana) = field2d(1:iana0,1:jana,2)
 else
   phige2d(1:iana0,1:jana) = field2d(1:iana0,1:jana,1)
 endif

 phige(1:iana0,1:jana)   = field2d(1:iana0,1:jana,2)
 fmaske(1:iana0,1:jana)  = field2d(1:iana0,1:jana,3)
 soile(1:iana0,1:jana)   = field2d(1:iana0,1:jana,21)

 phige2d(iana,1:jana) = phige2d(1,1:jana)
 phige(iana,1:jana) = phige(1,1:jana)
 fmaske(iana,1:jana) = fmaske(1,1:jana)
 soile(iana,1:jana) = soile(1,1:jana)

! conversion of fmask: 0-sea, 1-land to 1-sea, 0-land

 fmaske(:,:)=min(max(1.-fmaske(:,:),0.),1.)

 psel(1:iana0,1:jana)    = field2d(1:iana0,1:jana,4)
 pmsle(1:iana0,1:jana)   = field2d(1:iana0,1:jana,15)
 cctote(1:iana0,1:jana)  = field2d(1:iana0,1:jana,16)
 u10e(1:iana0,1:jana)    = field2d(1:iana0,1:jana,17)
 v10e(1:iana0,1:jana)    = field2d(1:iana0,1:jana,18)
 t2e(1:iana0,1:jana)     = field2d(1:iana0,1:jana,19)
 td2e(1:iana0,1:jana)    = field2d(1:iana0,1:jana,20)

 psel  (iana,1:jana) = psel  (1,1:jana)
 pmsle (iana,1:jana) = pmsle (1,1:jana)
 cctote(iana,1:jana) = cctote(1,1:jana)
 u10e  (iana,1:jana) = u10e  (1,1:jana)
 v10e  (iana,1:jana) = v10e  (1,1:jana)
 t2e   (iana,1:jana) = t2e   (1,1:jana)
 td2e  (iana,1:jana) = td2e  (1,1:jana)

   if (level_type == 1) then
     print *, 'input data on isobaric vertical levels with the following values:'
     do k=1,nlev_atm_inp
       print '(i3,a2,f5.0,a4)',k,'  ',lev_list(k)*0.01,'e+02'
     enddo
     print *
   else
     print *,'input data on hybrid vertical levels:'
     print*, "actual levels available (left col.) and ecmwf level index (right col.)"
     do k=1,nlev_atm_inp
       print '(i3,a,i3)',k,'   ',int(lev_list(k))
     enddo
     print*
   endif

   print *,'soil levels and mid-level depth in input data (cm):'
   depthsum=0.
   do k=1,nlev_soil_inp
     depth=2.*(lev_list_soil(k) - depthsum)
     depthsum=depthsum+depth
     print '(i3,2f8.1)',k,lev_list_soil(k)*1.e2, depth*1.e2
   enddo
   print*

! definition of qe as summ of specific humidity and total cloud water content

 qe(:,:,:) = qe(:,:,:) + qce(:,:,:)
 qe(:,:,:) = max(qe(:,:,:), zqcmin)

! definition of maximum and minumum volumetric soil content
! for input parameters (on the input data grid)
! conversion of volumetric water soil content into relative soil content

! parameters of soil in input data:
! nste - number of soil texture types
! qgmine - minimum volumetric water content (m**3/m**3)
! qgmaxe - maximum (saturation, porosity) volumetric water content (m**3/m**3)

! ifs-ecmwf global model

  iday=idate0(1)*10000+idate0(2)*100+idate0(3)
  soil_flag=1
#ifdef interim
  soil_flag=0
#endif

  if (iday >= 20071201.and.soil_flag==1) then

! after 01/12/2007:
! soil texture types: 1 - coarse, 2 - medium, 3 - medium-fine,
!                     4 - fine, 5 - very fine, 6 - organic
!                     7 - unknown!

    nste=7
    allocate(qgmine(nste))
    allocate(qgmaxe(nste))

    qgmine(1)=0.000
    qgmine(2)=0.000
    qgmine(3)=0.000
    qgmine(4)=0.000
    qgmine(5)=0.000
    qgmine(6)=0.000
    qgmine(7)=0.000
    qgmaxe(1)=0.403
    qgmaxe(2)=0.439
    qgmaxe(3)=0.430
    qgmaxe(4)=0.520
    qgmaxe(5)=0.614
    qgmaxe(6)=0.766
    qgmaxe(7)=0.590

    if (any(int(soile(:,:))==val_missing)) then
    write (*,*) 'soil type not available in input data file'
    stop
    endif

    do j=1,jana
    do i=1,iana
    iste=int(soile(i,j))
      if (iste>=1.and.iste<=nste) then
      qgmine2d(i,j)=qgmine(iste)
      qgmaxe2d(i,j)=qgmaxe(iste)
      else
      qgmine2d(i,j)=0.
      qgmaxe2d(i,j)=0.
      if (iste > nste) write (*,*) "soil texture type in input is outside the allowed range at point:",i,j,iste
      endif
    enddo
    enddo

  else

! before 01/12/2007:
! unique soil texture type: 1 - loamy

    soile(:,:)=0.
    nste=1
    allocate(qgmine(nste))
    allocate(qgmaxe(nste))
    qgmine(1)=0.000
    qgmaxe(1)=0.472

    qgmine2d(:,:)=qgmine(1)
    qgmaxe2d(:,:)=qgmaxe(1)

  endif

 if (surf_elaborate) then

 tskine(1:iana0,1:jana)  = field2d(1:iana0,1:jana,5)
 snowe(1:iana0,1:jana)   = field2d(1:iana0,1:jana,14)*1.e-3 ! conversion from kg m**-2 (mm) into m of water
 ficeee(1:iana0,1:jana)  = min( max( field2d(1:iana0,1:jana,22), 0.), 1.)
 tge(1:iana0,1:jana,1:nlev_soil_inp) = field3d_soil(1:iana0,1:jana,1:nlev_soil_inp,1)
 qge(1:iana0,1:jana,1:nlev_soil_inp) = max( field3d_soil(1:iana0,1:jana,1:nlev_soil_inp,2), 0.)
 ticeee(1:iana0,1:jana) = field3d_soil(1:iana0,1:jana,1,6)

 tskine(iana,1:jana) = tskine(1,1:jana)
 snowe(iana,1:jana) = snowe(1,1:jana)
 ficeee(iana,1:jana) = ficeee(1,1:jana)
 tge(iana,1:jana,1:nlev_soil_inp) = tge(1,1:jana,1:nlev_soil_inp)
 qge(iana,1:jana,1:nlev_soil_inp) = qge(1,1:jana,1:nlev_soil_inp)
 ticeee(iana,1:jana) = ticeee(1,1:jana)

 do j=1,jana
 do i=1,iana
 if (ficeee(i,j).gt.0.9) snowe(i,j) = .025
 qge(i,j,1:nlev_soil_inp)=min( max( (qge(i,j,1:nlev_soil_inp)-qgmine2d(i,j))/(qgmaxe2d(i,j)-qgmine2d(i,j)), 0.), 1.)
 enddo
 enddo

 endif ! surf_elaborate

 return
 end
!##################################################################################################
  subroutine conv_gfs_data

! procedure to convert meteorological fields derived from input data files for gfs data

  use param, only : iana, iana0, jana, nlev_atm_inp_max, nlev_soil_inp_max, nlev_atm_inp, nlev_soil_inp, &
                    eps, val_missing, pi, idate0, iperiod, iperiod_inp, surf_elaborate
  use gribana

  implicit none

  integer :: i, j, k, khumini, nday, npt
  integer, dimension(3)  :: iii
  integer, dimension(12) :: imon=(/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
  integer, save          :: nste
  real, dimension(iana,jana,nlev_atm_inp_max) :: rh, wspeed
  real :: wm, qsat, qsatw, qsati, esat, esatw, esati, eee, rhcwi, day, zflatt, zrid, depth, depthsum
  real, dimension(nlev_atm_inp_max)  :: p
  real, dimension(nlev_soil_inp_max) :: zdtg, zdqgrel

  mask_frame(:,:)=1

! note: the definitions below are used only for diagnostics
! date, time etc.

 if (iperiod_inp(1)==1) then         ! time unit is hour
   iperiod(1)=iperiod_inp(2)/24                             ! days
   iperiod(2)=iperiod_inp(2)-iperiod(1)*24                  ! hours
   iperiod(3)=0                                             ! minutes
 elseif (iperiod_inp(1)==0) then     ! time unit is minute
   iperiod(1)=iperiod_inp(2)/24/60                          ! days
   iperiod(2)=(iperiod_inp(2)-iperiod(1)*24*60)/60          ! hours
   iperiod(3)=iperiod_inp(2)-iperiod(1)*24*60-iperiod(2)*60 ! minutes
 endif

! definition of meteorological fields declared in main program
! list to be checked ! see subroutine read_grib2_data

 ge(1:iana0,1:jana,1:nlev_atm_inp) = field3d(1:iana0,1:jana,1:nlev_atm_inp,1)
 te(1:iana0,1:jana,1:nlev_atm_inp) = field3d(1:iana0,1:jana,1:nlev_atm_inp,2)
 ue(1:iana0,1:jana,1:nlev_atm_inp) = field3d(1:iana0,1:jana,1:nlev_atm_inp,3)
 ve(1:iana0,1:jana,1:nlev_atm_inp) = field3d(1:iana0,1:jana,1:nlev_atm_inp,4)
! qe(1:iana0,1:jana,1:nlev_atm_inp) = field3d(1:iana0,1:jana,1:nlev_atm_inp,5)
 rh(1:iana0,:,1:nlev_atm_inp) = field3d(1:iana0,1:jana,1:nlev_atm_inp,6)

 do k=1,nlev_atm_inp
   if (lev_list(k)>15000.) goto 100 ! for pressure levels > 150 hpa no modification is introduced
 enddo
100 khumini=k

 if (any(int(field3d(:,:,khumini:nlev_atm_inp,7))==-9999).or.maxval(field3d(:,:,khumini:nlev_atm_inp,7))<1.e-7) then
 iflag_cloude=0
 qce(:,:,:) = 0.
 print*, "no cloud total water (liquid+ice) available in input"
 else
 iflag_cloude=1
 qce(1:iana0,1:jana,khumini:nlev_atm_inp) = field3d(1:iana0,1:jana,khumini:nlev_atm_inp,7)
 qce(:,:,1:khumini) = 0.   ! definition of missing values of cloud water+ice at high levels
 print*, "cloud total water (liquid+ice) available in input"
 endif

 ge(iana,1:jana,1:nlev_atm_inp) = ge(1,1:jana,1:nlev_atm_inp)
 te(iana,1:jana,1:nlev_atm_inp) = te(1,1:jana,1:nlev_atm_inp)
 ue(iana,1:jana,1:nlev_atm_inp) = ue(1,1:jana,1:nlev_atm_inp)
 ve(iana,1:jana,1:nlev_atm_inp) = ve(1,1:jana,1:nlev_atm_inp)
! qe(iana,1:jana,1:nlev_atm_inp) = qe(1,1:jana,1:nlev_atm_inp)
 rh(iana,1:jana,1:nlev_atm_inp) = rh(1,1:jana,1:nlev_atm_inp)
 qce(iana,1:jana,1:nlev_atm_inp) = qce(1,1:jana,1:nlev_atm_inp)

! list of isobaric level:

 p(1:nlev_atm_inp)=lev_list(1:nlev_atm_inp)

 print*, "number of isobaric levels nlev_atm_inp", nlev_atm_inp
 print *,'isobaric levels for which t is available in input data:'
 do k = 1, nlev_atm_inp
 print '(i3,a5,f5.0,a4)',k,'  p =',p(k)*0.01,'e+02'
 enddo
 print*

   print *,'soil levels and depths in input data (cm):'
   depthsum=0.
   do k=1,nlev_soil_inp
   depth=2.*(lev_list_soil(k) - depthsum)
   depthsum=depthsum+depth
   print '(i3,2f8.1)',k,lev_list_soil(k)*1.e2, depth*1.e2
   enddo
   print*

! for the gfs model, phige2d and phige are the same

  phige2d(1:iana0,1:jana) = field2d(1:iana0,1:jana,1)
  phige2d(iana   ,1:jana) = phige2d(1,1:jana)
  phige(:,:) = phige2d(:,:)

  if (any(int(field2d(:,:,3))==int(val_missing)).and.any(int(field2d(:,:,46))==int(val_missing))) then
  print*, "The land-sea mask field is missing in GFS input grib data: stop!"
  stop
  endif

! Land-sea mask LANDN (coding discipline=2, category=0, parameter=218, index=46 here), valid
! after 19.07.2017, takes precedence over the "old" LAND (discipline=2, category=0,
! parameter=0, index=3 here), which is the standard land-sea mask to be used before 19.07.2017.

  if (any(int(field2d(:,:,3))/=int(val_missing))) then
  fmaske(1:iana0,1:jana)  = field2d(1:iana0,1:jana,3)
  endif

  if (any(int(field2d(:,:,46))/=int(val_missing))) then
  fmaske(1:iana0,1:jana)  = field2d(1:iana0,1:jana,46)
  endif

  fmaske(iana,1:jana) = fmaske(1,1:jana)

! conversion of fmask: 0-sea, 1-land to 1-sea, 0-land

  fmaske(:,:) = min(max(1.-fmaske(:,:),0.),1.)

 if (surf_elaborate) then

   tge(1:iana0,1:jana,1:nlev_soil_inp) = field3d_soil(1:iana0,1:jana,1:nlev_soil_inp,1)
   qge(1:iana0,1:jana,1:nlev_soil_inp) = field3d_soil(1:iana0,1:jana,1:nlev_soil_inp,2)

   tskine(1:iana0,1:jana)  = field2d(1:iana0,1:jana,5)
   snowe(1:iana0,1:jana)   = field2d(1:iana0,1:jana,14)*1.e-3 ! conversion from kg m**-2 (mm) into m of water
   ficeee(1:iana0,1:jana)   = min( max( field2d(1:iana0,1:jana,22), 0.), 1.)

   tge(iana,1:jana,1:nlev_soil_inp) = tge(1,1:jana,1:nlev_soil_inp)
   qge(iana,1:jana,1:nlev_soil_inp) = qge(1,1:jana,1:nlev_soil_inp)
   tskine(iana,1:jana) = tskine(1,1:jana)
   snowe(iana,1:jana) = snowe(1,1:jana)
   ficeee(iana,1:jana) = ficeee(1,1:jana)

! definition of temperature and soil water content in the case of missing data (-9999)
! (gfs t and qg over sea not defined)

  do k=1,nlev_soil_inp
  do j=1,jana
  do i=1,iana
  if (int(tge(i,j,k))==int(val_missing)) tge(i,j,k)=tskine(i,j)
  if (int(qge(i,j,k))==int(val_missing)) qge(i,j,k)=0.33 ! to be used to define qg on small islands
  enddo
  enddo
  enddo

! definition of maximum and minumum volumetric soil content
! for input parameters (on the input data grid)
! conversion of volumetric water soil content into relative soil content
! parameters of soil in input data:
! nste - number of soil texture types
! qgmine - minimum volumetric water content (m**3/m**3)
! qgmaxe - maximum (saturation, porosity) volumetric water conent (m**3/m**3)
! gfs (noaa-ncep) global model
! 0.470 is the max. water soil content ((m**3/m**3), this value has been
! empirically verified by data during period 01.06.2008-30.06.2009

    nste=1
    allocate(qgmine(nste))
    allocate(qgmaxe(nste))
    qgmine(1)=0.000
    qgmaxe(1)=0.470

    qgmine2d(:,:)=qgmine(1)
    qgmaxe2d(:,:)=qgmaxe(1)

    do k=1,nlev_soil_inp
    qge(:,:,k)=(qge(:,:,k)-qgmine2d(:,:))/(qgmaxe2d(:,:)-qgmine2d(:,:))
    enddo

  endif

! conversion from geop. height (gph) to geopotential (m**2 s**-2)

!!!!! ge(:,:,:)=ge(:,:,:)*9.8

! relative humidity is correctly defined only up to 100-150 hpa only

 do k = 1, nlev_atm_inp
 if(p(k).le.100.e2) rh(:,:,k) = 2.
 if(p(k).le.70.e2)  rh(:,:,k) = 0.1
 if(p(k).le.20.e2)  rh(:,:,k) = 0.001
 if(p(k).le.10.e2)  rh(:,:,k) = 0.0001
 enddo

! conversion of relative humidity into specific humidity

 do k = 1, nlev_atm_inp
 do j = 1, jana
 do i = 1, iana
 rh(i,j,k) = max(rh(i,j,k),0.)
 rh(i,j,k) = min(rh(i,j,k),102.)
 call qsat_tetens (te(i,j,k), p(k), eps, qsat, qsatw, qsati, esat, esatw, esati)
 eee=esat*rh(i,j,k)*1.e-2
 qe(i,j,k) = eps*eee/(eps*eee-eee+p(k))
 qe(i,j,k) = max (qe(i,j,k), 1.e-7)

! reduction of cloud water+ice

 qce(i,j,k)= min (qce(i,j,k), 1.e-3)    ! to avoid excessive condensate
 qe(i,j,k) = qe(i,j,k) + qce(i,j,k)     ! cloud water+ice added to specific humidity
 enddo
 enddo
 enddo

 qge(:,:,:) = max(min(qge(:,:,:),1.),0.) ! reset to min-max of relative values

 return
 end
!##################################################################################################
   subroutine filt2 (p, im, jm, anu)

   real p(im,jm), p2(im,0:jm+1)

!  2-grid interval filter over the whole domain

   do jlat = 1, jm
   do jlon = 2, im-1
   p2(jlon,jlat) = .25*(p(jlon-1,jlat)+p(jlon+1,jlat))+.5*p(jlon,jlat)
   enddo
   p2(1 ,jlat) = p2(im-1,jlat)
   p2(im,jlat) = p2(2,jlat)
   enddo

   do jlon = 1, im
   jlon1 = jlon + im/2
   if (jlon1 > im) jlon1 = jlon1-im
   p2(jlon,jm+1) = p2(jlon1,jm-1)
   p2(jlon,0   ) = p2(jlon1,2   )
   enddo

   do jlat = 1, jm
   p(1:im,jlat) = (1.-anu)*p(1:im,jlat) + anu*    &
                  (.25*p2(1:im,jlat+1)+.25*p2(1:im,jlat-1)+.5*p2(1:im,jlat))
   enddo

   return
   end
!=======================================================================
subroutine read_param_const(x0, y0, alon0, alat0, htop, htopvar, flag_data)

! If model_param_constant.bin exists, then reads constant (in time) model
! physiographical parameters from this file

use param
use mod_fields, only : fmask, rgm, rgq, soil_map, veg_map, &
 water_table_depth, tg_bottom, qg_rel_bottom, qg_rel_surf_approx, &
 soil_qmax, soil_qmin, soil_c, soil_rho, soil_psi, soil_k, soil_par_b, soil_par_c, soil_qrel_wilt, soil_qrel_ref, &
 soil_albedo_dry, soil_albedo_wet, soil_emiss1_dry, soil_emiss1_wet, soil_emiss2_dry, soil_emiss2_wet, &
 veg_lai, veg_frac, veg_root_depth, veg_roughness, veg_albedo, veg_emiss1, veg_emiss2, &
 ind_lev_soil_h_bottom, ind_lev_soil_w_bottom, veg_lai_max, &
 snow_dirt

implicit none
integer :: flag_data, iunit=10, nlon_read, nlat_read, nlevg_read, nst_read, nvt_read, &
 ierr_open, flag, i, j, k
real :: x0, y0, x0_read, y0_read, alon0, alat0, alon0_read, alat0_read, dlon_read, dlat_read
real, dimension(nlon,nlat) :: htop, htopvar
character(len=30) :: file_data="model_param_constant.bin"

 flag_data=0

 open (iunit, file=trim(file_data), form="unformatted", status="old", iostat=ierr_open)

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

 read(iunit) veg_lai_max

 close (iunit)

 write (*,*)
 write (*,*) 'File with model constant parameter fields ',trim(file_data),' read'
 write (*,*)

return
end subroutine read_param_const
!=======================================================================
subroutine write_param_const(x0, y0, alon0, alat0, htop, htopvar)

! Writes constant (in time) model physiographical parameters in model_param_constant.bin file

use param
use mod_fields, only : phig, fmask, rgm, rgq, soil_map, veg_map, &
  water_table_depth, tg_bottom, qg_rel_bottom, qg_rel_surf_approx, &
 soil_qmax, soil_qmin, soil_c, soil_rho, soil_psi, soil_k, soil_par_b, soil_par_c, soil_qrel_wilt, soil_qrel_ref, &
 soil_albedo_dry, soil_albedo_wet, soil_emiss1_dry, soil_emiss1_wet, soil_emiss2_dry, soil_emiss2_wet, &
 veg_lai, veg_frac, veg_root_depth, veg_roughness, veg_albedo, veg_emiss1, veg_emiss2, &
 ind_lev_soil_h_bottom, ind_lev_soil_w_bottom, veg_lai_max, &
 snow_dirt

implicit none
integer :: iunit=11, i, j, k
real :: x0, y0, alon0, alat0
real, dimension(nlon,nlat) :: htop, htopvar
character(len=30) :: file_data="model_param_constant.bin"

 open (iunit, file=trim(file_data), form="unformatted", status="unknown")

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

 close (iunit)

 write (*,*)
 write (*,*) 'Output file ',trim(file_data),' written' 
 write (*,*)

return
end subroutine write_param_const
!=======================================================================
subroutine read_forecast_mhf_soil(iflag)

! Reads a MHF file of BOLAM with surface and soil variables 
! simulated by a forecast run with siutable validity date and hour

use param
use mod_fields, only : iceth_frc, fice_frc, tskin_frc, tgsurf_frc, tg_frc, &
 qvsurf_frc, qgsurf_frc, qg_frc, fice_soil_surf_frc, fice_soil_frc, &
 snow_frc, snow_lev_frc, snow_t_frc, snow_fice_frc, snow_age_frc, snow_melt_age_frc, snow_dens_frc, snow_albedo

implicit none

integer :: iflag, iunit=11, iopen_err=0, ird, ierr=0, i, j, k, ndayr, &
 year_frc, month_frc, day_frc, hour_frc, minute_frc
character(len=30) :: file_inp="globo_forecast_soil.mhf",str_date0,str_frc_term,str_date_frc
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

 call rrec2 (iunit, nlon, nlat, tskin_frc(1:nlon,1:nlat))

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

 call rrec2 (iunit, nlon, nlat, snow_albedo(1:nlon,1:nlat))

! do ird=1,11 ! cswfl, clwfl, chflux, cqflux, t2min, t2max, runoff, runoff_tot, 3 additional 2D fields 
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
use mod_fields, only : iceth_frc, fice_frc, tskin_frc, tgsurf_frc, tg_frc, &
 qvsurf_frc, qgsurf_frc, qg_frc, fice_soil_surf_frc, fice_soil_frc, &
 snow_frc, snow_lev_frc, snow_t_frc, snow_fice_frc, snow_age_frc, snow_melt_age_frc, snow_dens_frc, snow_albedo

implicit none

integer :: iflag, iunit=11, iopen_err=0, ird, ierr=0, i, j, k, ndayr, &
 year_frc, month_frc, day_frc, hour_frc, minute_frc
character(len=30) :: file_inp="globo_forecast_soil.mhf", str_date0, str_frc_term, str_date_frc
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
   read(iunit) (tskin_frc(i,j),i=1,nlon)
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

 do j = 1,nlat
   read(iunit) (snow_albedo(i,j),i=1,nlon)
 enddo

! do ird=1,11 ! cswfl, clwfl, chflux, cqflux, t2min, t2max, runoff, runoff_tot, 3 additional 2D fields 
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
use mod_fields, only : iceth, fice, tskin, tgsurf, tg, qvsurf, qgsurf, qg, fice_soil_surf, fice_soil, &
 snow, snow_lev, snow_t, snow_fice, snow_age, snow_melt_age, snow_dens, &
 iceth_frc, fice_frc, tskin_frc, tgsurf_frc, tg_frc, qvsurf_frc, qgsurf_frc, qg_frc, fice_soil_surf_frc, fice_soil_frc, &
 snow_frc, snow_lev_frc, snow_t_frc, snow_fice_frc, snow_age_frc, snow_melt_age_frc, snow_dens_frc

implicit none

integer :: iflag, iunit=11, iopen_err=0, rec_ind, rec_lon, ird, ierr=0, i, j, k, ndayr, &
 year_frc, month_frc, day_frc, hour_frc, minute_frc
character(len=30) :: file_inp="globo_forecast_soil.bin",str_date0,str_frc_term,str_date_frc
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
   tskin_frc(1:nlon,j) = array_1d(1:rec_lon)
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
use mod_fields, only : iceth_frc, fice_frc, tskin_frc, tgsurf_frc, tg_frc, &
 qvsurf_frc, qgsurf_frc, qg_frc, fice_soil_surf_frc, fice_soil_frc, &
 snow_frc, snow_lev_frc, snow_t_frc, snow_fice_frc, snow_age_frc, snow_melt_age_frc, snow_dens_frc, snow_albedo

implicit none

integer :: iflag, iunit=11, iopen_err=0, ird, ierr=0, i, j, k, ndayr, &
 year_frc, month_frc, day_frc, hour_frc, minute_frc
character(len=30) :: file_inp="globo_forecast_soil.dat",str_date0,str_frc_term,str_date_frc
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
   read(iunit,1001) (tskin_frc(i,j),i=1,nlon)
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
 do j = 1,nlat
   read(iunit,1001) (snow_albedo(i,j),i=1,nlon)
 enddo

! do ird=1,11 ! snow cswfl, clwfl, chflux, cqflux, t2min, t2max, runoff, runoff_tot, 3 additional 2D fields 
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
!    include 'read_grib2_data.F90'
