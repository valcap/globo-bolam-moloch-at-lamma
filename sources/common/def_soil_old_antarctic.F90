! Last update 29/08/2024
!
! 2024.08: puting the not constant variables in order;
! new variables in output: t_bottom, qrel_bottom, veg_qrel_surf, water_table_depth
!
! May 2022: Correction of vegetation fraction in subroutine surfradpar,
! albedo of sea ice, emissivities of sea and sea ice
!
! 02/07/2021: change of dencity and heat capacity of upper horizonts of 
! dry soil for all soil types to forcing of heat isolation effect of upper
! soil layers (subroutine definition_soil_phys_param) 
!
! 08/07/2020: introduction of the coefficient (0-1) snow surface "dirtibility":
! growth of snow surface dirty due to vegetation waste, aerosol deposition, etc.,
! it is used in snow albedo definition.
!
! 23/08/2019: correction of dry and wet soil albedo and two emissivity for desert
! soil types (North Africa and Middle East, Yermosols, Xerosols, Lithosols,
! Dunes, Rock debris); correction of water table depth for soil types in
! South Mediterrenian zone (cambisols, Luvisols, Fluvisols) and for Cropland
! vegetation types; change of definition of water table depth using it's depth
! following soil and vegetation types, appending of dependence on orography
! height; reduction of hydraulic conductivity for 10 times.
! 
! 20/02/2019: version for ISAC models version 20
! 
! Nuovo suolo "Pochva", tutti i parametri di suolo sono 3d,
! utiluzzo diretto dei tipi pedologici FAO, senza passare in "texture types",
! nuovo dataset di vegetazione GLC2000 a posto di GLC1990,
! Lai e fazione di vegetazione vengono definiti (ciclo stagionale) non con
! le formule empiriche ma utilizzando i dataset appositi (36 file per l'anno)

!=======================================================================

subroutine physiographic_param (nlon, nlat, nst, nvt, nlevel_soil, &
 alon0, alat0, dlon, dlat, x0, y0, alon, alat, tsurf, soil_lev, jday, flag_variable_const, &
 htop, htopvar, fmask, water_table, t_bottom, qrel_bottom, veg_qrel_surf, &
 index_lev_tbottom, index_lev_qbottom, tsoil, qsoil_rel, &
 soil_map, veg_map, &
 soil_qmax, soil_qmin, soil_c, soil_rho, soil_psi, soil_k, &
 soil_par_b, soil_par_c, soil_qrel_wilt, soil_qrel_ref, &
 soil_albedo_dry, soil_albedo_wet, &
 soil_emiss1_dry, soil_emiss1_wet, soil_emiss2_dry, soil_emiss2_wet, &
 veg_lai, veg_frac, veg_root_depth, veg_roughness, veg_albedo, veg_emiss1, veg_emiss2, veg_lai_max, &
 snow_dirt)

! Procedure defines all physiographic parameters on the presciried grid:
! topography height, variance of topography height,
! land-sea fraction, soil and vegetation parameters;
! soil temperature and water content on the bottom level and also at
! all soil levels in following some climatic hypothesis and surface data.

! Input variables:

!      nlon: grid point number along axis of rotated longitude;
!      nlat: grid point number along axis of rotated latidude;
!      nst: total number soil types (classification FAO);
!      nvt: total number vegetaion types (classification GLC2000);
!      nlevel_soil: total number of soil levels;
!      dlon: grid step along the domain axis of rotated longitude (degree);
!      dlat: grid step along the domain axis of rotated latitude (degree);
!      alon: geographical longitude of the grid points (degree, -180...180);
!      alat: geographical latitude of the grid points (degree, -90...90);
!      tsurf: surface temperature (K), may be replaced by air temperatureat 2 m;
!      level_soil: soil levels values (m);
!      jday: Julian day minus 1 (0-364).;
!      flag_variable_const: flag for definition of variable in time parameters only (veg. LAI, veg. fraction, tsoil and qsoil_rel profiles): 
!                 1 - definition of all output variables,
!             not 1 - definition of variable parameters only.

! Output variables:

!      fmask: land-sea franction: 1-sea, 0-land (proportion);
!      htop: orography height (m);
!      htopvar: orography variance (later, std. dev.) (m);
!      water_table: depth of water table (m);
!      t_bottom: temperature at soil bottom depth for heat precesses (K);
!      qrel_bottom: relative soil water content at soil bottom depth for water precesses (proportion);
!      veg_qrel_surf: relative soil water content at soil surface (proportion);
!      index_lev_tbottom: index of bottom soil level for heat trasport at each grid point;
!      index_lev_qbottom: index of bottom soil level for water trasport at each grid point;
!      tsoil: soil temperature at all soil levels (K);
!      qsoil_rel: relative soil water content at all soil levels (proportion); 
!      soil_map:  soil types (classification FAO, 31 types) defined on the grid points, 3D array,
!                1-st and 2-nd array index are grid point index, 3-rd index is soil type number plus 1,
!                if 3-rd index has value 1, then the variable means dominate soil type,
!                values 2-32 of 3-rd index means fraction (of area, proportion) of the each of 31 soil types;
!      veg_map: vegetation (landuse) types (classification GLC2000, 22 types) defined on the grid points, 3D array,
!               1-st and 2-nd array index are grid point index, 3-rd index is vegetation types number plus 1,
!               if 3-rd index has value 1, then the variable means dominate vegetation type,
!               values 2-23 of 3-rd index means fraction (of area, proportion) of the each of 22 vegetation types;
!      soil_qmax: maximum specific volumetric soil water content (m^3/m^3) at all soil levels, 3D array;
!      soil_qmin: minimum specific volumetric soil water content (m^3/m^3) at all soil levels, 3D array;
!      soil_c: dry soil thermal capacity (J/K/m^3) at all soil levels, 3D array;
!      soil_rho: dry soil density (kg/m^3) at all soil levels, 3D array;
!      soil_psi: idraulic potential of saturated soil (m) at all soil levels, 3D array;
!      soil_k: idraulic conductivity of saturated soil (m/s) at all soil levels, 3D array;
!      soil_par_b: parameter in soil water transport equation (dimensionless) at all soil levels, 3D array;
!      soil_par_c: parameter in soil water transport equation (dimensionless) at all soil levels, 3D array;
!      soil_qrel_wilt: relative soil water content (proportion) at which the vegetation wilt at all soil levels, 3D array;
!      soil_qrel_ref: relative soil water content (proportion) at which the vegetation stop evapotraspiration increase at all soil levels, 3D array;
!      soil_albedo_dry: dry soil albedo (proportion), 2D array;
!      soil_albedo_wet: wet soil albedo (proportion), 2D array;
!      soil_albedo_dry: dry soil albedo (proportion), 2D array;
!      soil_albedo_wet: wet soil albedo (proportion), 2D array;
!      soil_emiss_1_dry: dry soil emissivity in broadband window (proportion), 2D array;
!      soil_emiss_1_wet: wet soil emissivity in broadband window (proportion), 2D array;
!      soil_emiss_2_dry: dry soil emissivity in 8-12 micron window (proportion), 2D array;
!      soil_emiss_2_wet: wet soil emissivity in 8-12 micron window (proportion), 2D array;
!      veg_lai: Leaf Area Index (LAI) of vegetation (dimensionless), 2D array;
!      veg_frac: fraction of vegetation (proportion), 2D array;
!      veg_root_depth: root depth of vegetation (m), 2D array;
!      veg_roughness: vegetation cover roughness (m), 2D array;
!      veg_albedo: vegetation albedo (proportion), 2D array;
!      veg_emiss_1: vegetaion emissivity in broadband window (proportion), 2D array;
!      veg_emiss_2: vegetation emissivity in 8-12 micron window (proportion), 2D array.
!      veg_lai_max: LAI maximum value permited for used vegetation dataset (refer value), real variable;
!      snow_dirt: snow "dirtibility", that is coefficient (proportion 0-1) of dirt
! growth of snow surface dirty due to vegetation waste, aerosol deposition, etc.,
! it is used in snow albedo definition.

implicit none

! Input:

integer :: nlon, nlat, nst, nvt, nlevel_soil, jday, flag_variable_const
real :: alon0, alat0, dlon, dlat, x0, y0
real, dimension(nlon,nlat) :: alon, alat, tsurf
real, dimension(nlevel_soil) :: soil_lev

! Ouput:

real, dimension(nlon,nlat) :: htop, htopvar, fmask, water_table, t_bottom, qrel_bottom, veg_qrel_surf
real, dimension(nlon,nlat,nst+1) :: soil_map
real, dimension(nlon,nlat,nvt+1) :: veg_map
real, dimension(nlon,nlat,nlevel_soil) :: soil_qmax, soil_qmin, soil_c, soil_rho, &
 soil_psi, soil_k, soil_par_b, soil_par_c, soil_qrel_wilt, soil_qrel_ref
real, dimension(nlon,nlat) :: soil_albedo_dry, soil_albedo_wet, &
 soil_emiss1_dry, soil_emiss1_wet, soil_emiss2_dry, soil_emiss2_wet, &
 veg_lai, veg_frac, veg_root_depth, veg_roughness, veg_albedo, veg_emiss1, veg_emiss2, snow_dirt
real :: veg_lai_max
integer, dimension(nlon,nlat) :: index_lev_tbottom, index_lev_qbottom
real, dimension(nlon,nlat,nlevel_soil) :: tsoil, qsoil_rel

! Internal:

integer, parameter :: nvt_old=14
real, dimension(nlon,nlat) :: fmask1, fmask_ecmwf, t_climate, t_amplitude, &
 soil_water_table, veg_water_table, depth_bottom, q_bottom
real, dimension(nlon,nlat,nvt_old+1) :: veg_map_old
real, dimension(nlon,nlat,nlevel_soil) :: qrel
real :: val_missing=-9999., pi, x, xx, t1, dt_work
integer :: flag_glob=0, iopen_err, i, j, k, nbottom, npolcap
real, dimension(1) :: jday_work, a, b, c, d, dt, jday_work_sh, a_sh, b_sh, c_sh, d_sh, dt_sh

! For analytical hypothetical vertical profile of soil temperature:
real, dimension(5) :: jday0=(/  0.,  90., 181., 274., 366./)
real, dimension(5) :: a0=   (/ 0.5,  1.5,  0.5,  1.5,  0.5/)
real, dimension(5) :: b0=   (/ 1.5,  0.5,  0.5, -0.5,  1.5/)
real, dimension(5) :: c0=   (/ 1.0,  0.5,  1.0,  0.5,  1.0/)
real, dimension(5) :: d0=   (/-1.0, -0.5,  1.0,  0.5, -1.0/)
real, dimension(5) :: dt0=  (/-10.,  -3.,  10.,   3., -10./)
integer, dimension(500) :: iv
real, dimension(nlevel_soil) :: func

character (len=20) :: fl_land='land_par.bin'

integer, dimension(5) :: idate0=(/2018,09,20, 0, 0/)
integer, dimension(3) :: iperiod=(/0, 0, 0/)
integer, dimension(12) :: imon=(/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)


!---------------------------------------------------------------------------

 pi = abs(acos(-1.))

 a0(:)=a0(:)*pi
 b0(:)=b0(:)*pi

 flag_glob=0
 if (abs(x0) < 0.1.and.abs(y0) < 0.1) then
   if (abs(alon0) < dlon.and.alat0 < -90.+dlat.and.&
alon0+float(nlon-1)*dlon > 360.-dlon*2..and.alat0+float(nlat-1)*dlat > 90.-dlat) then
     flag_glob=1
     npolcap=8
     if (dlat < 0.20) npolcap=10 ! 22 km
     if (dlat < 0.18) npolcap=12 ! 19-20 km
     if (dlat < 0.15) npolcap=14 ! 17-18 km
     if (dlat < 0.12) npolcap=16 ! 14-15 km
   endif
 endif

 j=jday+1
 do i=1,12
   j=j-imon(i)
   if (j <= 0) exit
 enddo
 idate0(2)=i

 j=jday+1
 do i=1,idate0(2)-1
   j=j-imon(i)
 enddo
 idate0(3)=j

! Reading of geographical files, in case they are pre-computed - if they do not exist, the same
! quantities are then computed below

 if (flag_variable_const == 1) then

   open (10,file=fl_land,form='unformatted',status='old',iostat=iopen_err)

   if (iopen_err == 0) then ! File with phisiographic parameters exists

     read (10) htop
     read (10) htopvar
     read (10) fmask
     read (10) t_climate
     read (10) t_amplitude
     read (10) soil_map
     read (10) veg_map
     read (10) veg_map_old
     close (10)

   else ! File with phisiographic parameters do not exist, and this parameters must be definded


! Reading of global orography dataset and definition of orography
! and its standard deviation on the model grid

     call orogdef (alon,alat,htop,htopvar,alon0,alat0,dlon,dlat,nlon,nlat,flag_glob)
     if (flag_glob == 1) then
       call polavert (htop, nlon, nlat, npolcap)
       call polavert (htopvar, nlon, nlat, npolcap)
       htopvar(:,:) = max (htopvar(:,:), 0.)
     endif

! Reading of global soil (FAO), definition of soil types (FAO classification) destribution on the model grid

     call dataset_soil_fao(nlon, nlat, dlon, dlat, alon, alat, flag_glob, soil_map)

! Reading of global landuse - vegetation datasets (GLC1990) and
! definition of vegetation types (GLC1990) destribution on the model grid,
! and definition of land-sea fraction

     call dataset_vegetation(nlon, nlat, dlon, dlat, alon, alat, flag_glob, veg_map_old, fmask1)

! Reading of global landuse - vegetation datasets (GLC2000) and
! definition of vegetation types (GLC2000) destribution on the model grid,
! and definition of lan-sea fraction

     call dataset_vegetation_cyclopes(nlon, nlat, x0, y0, dlon, dlat, alon0, alat0, alon, alat, flag_glob, &
 nvt, val_missing, veg_map, fmask)

! Dataset GLC2000 cover the area from -60 to 80 degree of latitude only,
! and land-sea fraction are defined on this area only,
! land-sea fraction in polare zones are defined with using land-sea fraction accoding
! GLC1990 landuse dataset

!     if (any(veg_map(:,:,1) == val_missing)) fmask(:,:)=fmask1(:,:)

     fmask(:,:)=fmask1(:,:) ! fmask1 is defined using a dataset with higher resoltution

! For the antarctic area: vegetation dataset (GLC1990) do not include perennial sea ice close to the antarctic coast;
! this omission is corrected with land-sea fraction of ifs_ecmwf (fmask_ecmwf).

! Reading of global land-sea fraction data of ifs-ecmwf with resolution 0.25 degree and 
! definition of land-sea fraction according to this data on the model grid

     call dataset_ecmwf_lsm(nlon, nlat, alon, alat, flag_glob, fmask_ecmwf)
 
! Reading of global climate temperature and amplitude of climate temperature,
! following 36 year (1979-2014) re-analysis data of IFS-ECMWF with resolution 0.75 degree and 
! interpolation of this variables on the model grid

     call dataset_t_climate(nlon, nlat, alon, alat, flag_glob, t_climate, t_amplitude)

! Correction of soil, vegetation and land-sea fraction maps defined in the model grid,
! to coordinate among themselves 

     print *
     print *,"Coordination among various defined physiographic datasets, corrections of them"

     call soil_veg_lsf_correct(nlon, nlat, nst, nvt, 2, alat, fmask_ecmwf, &
 soil_map, veg_map, fmask, t_climate, t_amplitude)

! Writing of computed geographical quantities - needed for the operational suite only!

     open (10,file=fl_land,form='unformatted',status='unknown')
     write (10) htop
     write (10) htopvar
     write (10) fmask
     write (10) t_climate
     write (10) t_amplitude
     write (10) soil_map
     write (10) veg_map
     write (10) veg_map_old
     close (10)

   endif ! File with phisiographic parameters exists or do not exist

! call outgraph(80102,2,nlon,nlat,x0,y0,alon0,alat0,dlon,dlat,&
! 2, 0,  0,  1, 1, 0,idate0(1:5),iperiod(1:3),fmask(:,:),1.,0.)
! call outgraph(80102,2,nlon,nlat,x0,y0,alon0,alat0,dlon,dlat,&
! 2, 0,193,  1, 1, 0,idate0(1:5),iperiod(1:3),soil_map(:,:,1),1.,0.)
! call outgraph(80102,2,nlon,nlat,x0,y0,alon0,alat0,dlon,dlat,&
! 2, 0,195,  1, 1, 0,idate0(1:5),iperiod(1:3),veg_map(:,:,1),1.,0.)
! call outgraph(80102,2,nlon,nlat,x0,y0,alon0,alat0,dlon,dlat,&
! 2, 0,194,  1, 1, 0,idate0(1:5),iperiod(1:3),veg_map_old(:,:,1),1.,0.)

! call outgraph(80102,2,nlon,nlat,x0,y0,alon0,alat0,dlon,dlat,&
! 2, 0,  0,  1, 1, 0,idate0(1:5),iperiod(1:3),fmask(:,:),1.,0.)
! call outgraph(80102,2,nlon,nlat,x0,y0,alon0,alat0,dlon,dlat,&
! 2, 0,  7,  1, 1, 0,idate0(1:5),iperiod(1:3),htop(:,:),1.,0.)
! call outgraph(80102,2,nlon,nlat,x0,y0,alon0,alat0,dlon,dlat,&
! 2, 0,  7,  1, 1, 0,idate0(1:5),iperiod(1:3),htopvar(:,:),1.,0.)
! call outgraph(80102,2,nlon,nlat,x0,y0,alon0,alat0,dlon,dlat,&
! 0, 0,  0,  1, 1, 0,idate0(1:5),iperiod(1:3),t_climate(:,:),1.,-273.15)
! call outgraph(80102,2,nlon,nlat,x0,y0,alon0,alat0,dlon,dlat,&
! 0, 0,  0,  1, 1, 0,idate0(1:5),iperiod(1:3),t_amplitude(:,:),1.,0.)

! Definition of physical parameters of soil using soil types (FAO classification) defined on the model grid

 call definition_soil_phys_param(soil_map, nlevel_soil, soil_lev, nlon, nlat, &
 soil_qmax, soil_qmin, soil_c, soil_rho, soil_psi, soil_k, soil_par_b, soil_par_c, soil_qrel_wilt, soil_qrel_ref, &
 soil_water_table, soil_albedo_dry, soil_albedo_wet, soil_emiss1_dry, soil_emiss1_wet, soil_emiss2_dry, soil_emiss2_wet)

! call outgraph(80102,2,nlon,nlat,x0,y0,alon0,alat0,dlon,dlat,&
! 2, 0,  0,  1, 1, 0,idate0(1:5),iperiod(1:3),fmask(:,:),1.,0.)
!do k=1,nlevel_soil
!   call outgraph(80102,2,nlon,nlat,x0,y0,alon0,alat0,dlon,dlat,&
! 2, 0,193,  1, 1, 0,idate0(1:5),iperiod(1:3),soil_psi(:,:,k),1.,0.)
! print *,k,maxval(soil_psi(:,:,k)),minval(soil_psi(:,:,k))
!enddo

! call outgraph(80102,2,nlon,nlat,x0,y0,alon0,alat0,dlon,dlat,&
! 2, 0,  0,  1, 1, 0,idate0(1:5),iperiod(1:3),fmask(:,:),1.,0.)
! call outgraph(80102,2,nlon,nlat,x0,y0,alon0,alat0,dlon,dlat,&
! 2, 0,193,  1, 1, 0,idate0(1:5),iperiod(1:3),soil_emiss1_dry,1.,0.)
! call outgraph(80102,2,nlon,nlat,x0,y0,alon0,alat0,dlon,dlat,&
! 2, 0,193,  1, 1, 0,idate0(1:5),iperiod(1:3),soil_emiss1_wet,1.,0.)
! call outgraph(80102,2,nlon,nlat,x0,y0,alon0,alat0,dlon,dlat,&
! 2, 0,193,  1, 1, 0,idate0(1:5),iperiod(1:3),soil_emiss2_dry,1.,0.)
! call outgraph(80102,2,nlon,nlat,x0,y0,alon0,alat0,dlon,dlat,&
! 2, 0,193,  1, 1, 0,idate0(1:5),iperiod(1:3),soil_emiss2_wet,1.,0.)
! call outgraph(80102,2,nlon,nlat,x0,y0,alon0,alat0,dlon,dlat,&
! 2, 0,193,  1, 1, 0,idate0(1:5),iperiod(1:3),soil_albedo_dry,1.,0.)
! call outgraph(80102,2,nlon,nlat,x0,y0,alon0,alat0,dlon,dlat,&
! 2, 0,193,  1, 1, 0,idate0(1:5),iperiod(1:3),soil_albedo_wet,1.,0.)
! call outgraph(80102,2,nlon,nlat,x0,y0,alon0,alat0,dlon,dlat,&
! 2, 0,193,  1, 1, 0,idate0(1:5),iperiod(1:3),soil_water_table,1.,0.)
!print *,maxval(soil_emiss1_dry),minval(soil_emiss1_dry)
!print *,maxval(soil_emiss1_wet),minval(soil_emiss1_wet)
!print *,maxval(soil_emiss2_dry),minval(soil_emiss2_dry)
!print *,maxval(soil_emiss2_wet),minval(soil_emiss2_wet)
!print *,maxval(soil_albedo_dry),minval(soil_albedo_dry)
!print *,maxval(soil_water_table),minval(soil_water_table)

! Definition of physical parameters of vegetation using vegetation types GLC1990 defined on the model map

! call definition_veg1_phys_param(veg_map_old, alat, nlon, nlat, idate0(2), idate0(3), &
! veg_lai, veg_lai_max, veg_frac, veg_root_depth, veg_roughness, veg_albedo, veg_emiss1, veg_emiss2)

! Definition of constant physical parameters of vegetation using vegetation
! types GLC2000 defined on the model map

 snow_dirt(:,:)=1. ! Snow may be dirty (snow albedo is low) due to vegetation waste, aerosol deposition, etc.

 call definition_veg2_const_phys_param(veg_map, t_climate, t_amplitude, nlon, nlat, &
 veg_water_table, veg_qrel_surf, veg_root_depth, veg_roughness, &
 veg_albedo, veg_emiss1, veg_emiss2, snow_dirt)

 if (flag_glob == 1) then
   call polavert (veg_roughness, nlon, nlat, npolcap)
   veg_roughness(:,:) = max (veg_roughness(:,:), 1.e-8)
 endif

! Some correction of snow surface "dirtibility"

 do j=1,nlat
 do i=1,nlon
   if (soil_psi(i,j,1) > 0.) snow_dirt(i,j)=0. ! Glacier
   if (fmask(i,j) >= 0.5) snow_dirt(i,j)=0.2   ! Sea ice
!! Antarctic coasts: snow "dirtibility" 0.3 for not glacier surface
!   if (soil_psi(i,j,1) < 0.) snow_dirt(i,j)=0.3
 enddo
 enddo

! Definition of soil bottom condition

 call definition_soil_bottom(soil_map, veg_map, t_amplitude, t_climate, htop, veg_water_table, soil_water_table, water_table, &
 veg_qrel_surf, soil_qmax, soil_qmin, soil_c, soil_rho, soil_psi, soil_par_b, &
 alat, soil_lev, nlon, nlat, nst, nvt, nlevel_soil, &
 index_lev_tbottom, index_lev_qbottom, depth_bottom, t_bottom, q_bottom, qrel_bottom)

 endif ! flag_variable_const

! Definition of physical parameters, that have season cycle,
! of vegetation, using vegetation types GLC2000 defined on the model map

 call veget_param_time_series(nlon, nlat, alon, alat, x0, y0, dlon, dlat, alon0, alat0, flag_glob, &
 jday, veg_lai, veg_lai_max, veg_frac, val_missing)

! Correction over glacier and water body

 do j=1,nlat
 do i=1,nlon
   if (int(soil_map(i,j,1)) >= nst-1) then
     veg_lai(i,j)=0.
     veg_frac(i,j)=0.
   endif
 enddo
 enddo
!
! call outgraph(80102,2,nlon,nlat,x0,y0,alon0,alat0,dlon,dlat,&
! 2, 0,  0,  1, 1, 0,idate0(1:5),iperiod(1:3),fmask(:,:),1.,0.)
! call outgraph(80102,2,nlon,nlat,x0,y0,alon0,alat0,dlon,dlat,&
! 2, 0,193,  1, 1, 0,idate0(1:5),iperiod(1:3),veg_lai,1.,0.)
! call outgraph(80102,2,nlon,nlat,x0,y0,alon0,alat0,dlon,dlat,&
! 2, 0,193,  1, 1, 0,idate0(1:5),iperiod(1:3),veg_frac,1.,0.)

! call outgraph(80102,2,nlon,nlat,x0,y0,alon0,alat0,dlon,dlat,&
! 2, 0,193,  1, 1, 0,idate0(1:5),iperiod(1:3),veg_root_depth,1.,0.)
! call outgraph(80102,2,nlon,nlat,x0,y0,alon0,alat0,dlon,dlat,&
! 2, 0,193,  1, 1, 0,idate0(1:5),iperiod(1:3),veg_roughness,1.,0.)
! call outgraph(80102,2,nlon,nlat,x0,y0,alon0,alat0,dlon,dlat,&
! 2, 0,193,  1, 1, 0,idate0(1:5),iperiod(1:3),veg_emiss1,1.,0.)
! call outgraph(80102,2,nlon,nlat,x0,y0,alon0,alat0,dlon,dlat,&
! 2, 0,193,  1, 1, 0,idate0(1:5),iperiod(1:3),veg_emiss2,1.,0.)
! call outgraph(80102,2,nlon,nlat,x0,y0,alon0,alat0,dlon,dlat,&
! 2, 0,193,  1, 1, 0,idate0(1:5),iperiod(1:3),veg_albedo,1.,0.)
! call outgraph(80102,2,nlon,nlat,x0,y0,alon0,alat0,dlon,dlat,&
! 2, 0,193,  1, 1, 0,idate0(1:5),iperiod(1:3),veg_water_table,1.,0.)
!print *,maxval(veg_lai),minval(veg_lai)
!print *,veg_lai_max
!print *,maxval(veg_frac),minval(veg_frac)
!print *,maxval(veg_root_depth),minval(veg_root_depth)
!print *,maxval(veg_roughness),minval(veg_roughness)
!print *,maxval(veg_emiss1),minval(veg_emiss1)
!print *,maxval(veg_emiss2),minval(veg_emiss2)
!print *,maxval(veg_albedo),minval(veg_albedo)

!! All surface radiation parameters:
! call outgraph(80102,2,nlon,nlat,x0,y0,alon0,alat0,dlon,dlat,&
! 2, 0,  0,  1, 1, 0,idate0(1:5),iperiod(1:3),fmask(:,:),1.,0.)
! call outgraph(80102,2,nlon,nlat,x0,y0,alon0,alat0,dlon,dlat,&
! 2, 0,193,  1, 1, 0,idate0(1:5),iperiod(1:3),soil_emiss1_dry,1.,0.)
! call outgraph(80102,2,nlon,nlat,x0,y0,alon0,alat0,dlon,dlat,&
! 2, 0,193,  1, 1, 0,idate0(1:5),iperiod(1:3),soil_emiss1_wet,1.,0.)
! call outgraph(80102,2,nlon,nlat,x0,y0,alon0,alat0,dlon,dlat,&
! 2, 0,193,  1, 1, 0,idate0(1:5),iperiod(1:3),soil_emiss2_dry,1.,0.)
! call outgraph(80102,2,nlon,nlat,x0,y0,alon0,alat0,dlon,dlat,&
! 2, 0,193,  1, 1, 0,idate0(1:5),iperiod(1:3),soil_emiss2_wet,1.,0.)
! call outgraph(80102,2,nlon,nlat,x0,y0,alon0,alat0,dlon,dlat,&
! 2, 0,193,  1, 1, 0,idate0(1:5),iperiod(1:3),veg_emiss1,1.,0.)
! call outgraph(80102,2,nlon,nlat,x0,y0,alon0,alat0,dlon,dlat,&
! 2, 0,193,  1, 1, 0,idate0(1:5),iperiod(1:3),veg_emiss2,1.,0.)
! call outgraph(80102,2,nlon,nlat,x0,y0,alon0,alat0,dlon,dlat,&
! 2, 0,193,  1, 1, 0,idate0(1:5),iperiod(1:3),soil_albedo_dry,1.,0.)
! call outgraph(80102,2,nlon,nlat,x0,y0,alon0,alat0,dlon,dlat,&
! 2, 0,193,  1, 1, 0,idate0(1:5),iperiod(1:3),soil_albedo_wet,1.,0.)
! call outgraph(80102,2,nlon,nlat,x0,y0,alon0,alat0,dlon,dlat,&
! 2, 0,193,  1, 1, 0,idate0(1:5),iperiod(1:3),veg_albedo,1.,0.)
! call outgraph(80102,2,nlon,nlat,x0,y0,alon0,alat0,dlon,dlat,&
! 2, 0,193,  1, 1, 0,idate0(1:5),iperiod(1:3),snow_dirt,1.,0.)
!print *,maxval(soil_emiss1_dry),minval(soil_emiss1_dry)
!print *,maxval(soil_emiss1_wet),minval(soil_emiss1_wet)
!print *,maxval(soil_emiss2_dry),minval(soil_emiss2_dry)
!print *,maxval(soil_emiss2_wet),minval(soil_emiss2_wet)
!print *,maxval(veg_emiss1),minval(veg_emiss1)
!print *,maxval(veg_emiss2),minval(veg_emiss2)
!print *,maxval(soil_albedo_dry),minval(soil_albedo_dry)
!print *,maxval(veg_albedo),minval(veg_albedo)

!! All constant variables 2d:
! call outgraph(80102,2,nlon,nlat,x0,y0,alon0,alat0,dlon,dlat,&
! 2, 0,  0,  1, 1, 0,idate0(1:5),iperiod(1:3),fmask(:,:),1.,0.)
! call outgraph(80102,2,nlon,nlat,x0,y0,alon0,alat0,dlon,dlat,&
! 2, 0,  7,  1, 1, 0,idate0(1:5),iperiod(1:3),htop(:,:),1.,0.)
! call outgraph(80102,2,nlon,nlat,x0,y0,alon0,alat0,dlon,dlat,&
! 2, 0,  7,  1, 1, 0,idate0(1:5),iperiod(1:3),htopvar(:,:),1.,0.)
! call outgraph(80102,2,nlon,nlat,x0,y0,alon0,alat0,dlon,dlat,&
! 2, 0,193,  1, 1, 0,idate0(1:5),iperiod(1:3),t_climate,1.,-273.15)
! call outgraph(80102,2,nlon,nlat,x0,y0,alon0,alat0,dlon,dlat,&
! 2, 0,193,  1, 1, 0,idate0(1:5),iperiod(1:3),t_amplitude,1.,0.)
! call outgraph(80102,2,nlon,nlat,x0,y0,alon0,alat0,dlon,dlat,&
! 2, 0,193,  1, 1, 0,idate0(1:5),iperiod(1:3),soil_water_table,1.,0.)
! call outgraph(80102,2,nlon,nlat,x0,y0,alon0,alat0,dlon,dlat,&
! 2, 0,193,  1, 1, 0,idate0(1:5),iperiod(1:3),veg_water_table,1.,0.)
! call outgraph(80102,2,nlon,nlat,x0,y0,alon0,alat0,dlon,dlat,&
! 2, 0,193,  1, 1, 0,idate0(1:5),iperiod(1:3),water_table,1.,0.)
!! call outgraph(80102,2,nlon,nlat,x0,y0,alon0,alat0,dlon,dlat,&
!! 2, 0,193,  1, 1, 0,idate0(1:5),iperiod(1:3),veg_qrel_surf,1.,0.)
! call outgraph(80102,2,nlon,nlat,x0,y0,alon0,alat0,dlon,dlat,&
! 2, 0,193,  1, 1, 0,idate0(1:5),iperiod(1:3),depth_bottom,1.,0.)
!! call outgraph(80102,2,nlon,nlat,x0,y0,alon0,alat0,dlon,dlat,&
!! 2, 0,193,  1, 1, 0,idate0(1:5),iperiod(1:3),t_bottom,1.,-273.15)
!! call outgraph(80102,2,nlon,nlat,x0,y0,alon0,alat0,dlon,dlat,&
!! 2, 0,193,  1, 1, 0,idate0(1:5),iperiod(1:3),q_bottom,1.,0.)
! call outgraph(80102,2,nlon,nlat,x0,y0,alon0,alat0,dlon,dlat,&
! 2, 0,193,  1, 1, 0,idate0(1:5),iperiod(1:3),qrel_bottom,1.,0.)
! call outgraph(80102,2,nlon,nlat,x0,y0,alon0,alat0,dlon,dlat,&
! 2, 0,193,  1, 1, 0,idate0(1:5),iperiod(1:3),veg_root_depth,1.,0.)
! call outgraph(80102,2,nlon,nlat,x0,y0,alon0,alat0,dlon,dlat,&
! 2, 0,193,  1, 1, 0,idate0(1:5),iperiod(1:3),veg_roughness,1.,0.)
!print *,maxval(htop),minval(htop)
!print *,maxval(htopvar),minval(htopvar)
!print *,maxval(t_climate),minval(t_climate)
!print *,maxval(t_amplitude),minval(t_amplitude)
!print *,maxval(soil_water_table),minval(soil_water_table)
!print *,maxval(veg_water_table),minval(veg_water_table)
!print *,maxval(veg_qrel_surf),minval(veg_qrel_surf)
!print *,maxval(depth_bottom),minval(depth_bottom)
!!print *,maxval(t_bottom),minval(t_bottom)
!!print *,maxval(q_bottom),minval(q_bottom)
!print *,maxval(qrel_bottom),minval(qrel_bottom)
!print *,maxval(veg_root_depth),minval(veg_root_depth)
!print *,maxval(veg_roughness),minval(veg_roughness)

! Variables NOT constant in time: temperature soil and water content soil
! vertical profiles, using the following variables:
! jday, soil_lev, fmask, alat, soil_map, soil_qmax, soil_qmin,
! t_bottom, tsurf, index_lev_tbottom,
! qrel_bottom, veg_qrel_surf, index_lev_qbottom

! Definition of temperature at all soil levels using temperature at surface,
! bottom level temperature and analytical hypotetical vertical profile

! Interpolation of profile parameters for prescribed Julian day using values
! defined for various season (a0, b0, c0, d0, td0)

! North hemisphere
  
 jday_work(1)=float(jday)
  
 call near(jday_work(1:1), 1, jday0, 5, iv(1:1))
  
 call interp_spline_1d(a(1:1), jday_work(1:1), 1, a0, jday0, 5, iv(1:1), 1., 0., 0.)
 call interp_spline_1d(b(1:1), jday_work(1:1), 1, b0, jday0, 5, iv(1:1), 1., 0., 0.)
 call interp_spline_1d(c(1:1), jday_work(1:1), 1, c0, jday0, 5, iv(1:1), 1., 0., 0.)
 call interp_spline_1d(d(1:1), jday_work(1:1), 1, d0, jday0, 5, iv(1:1), 1., 0., 0.)
 call interp_spline_1d(dt(1:1), jday_work(1:1), 1, dt0, jday0, 5, iv(1:1), 1., 0., 0.)
  
! South hemisphere
  
 jday_work_sh(1)=jday_work(1)+183.
 if (jday_work_sh(1) >= 365.) jday_work_sh(1)=jday_work_sh(1)-365.
  
 call near(jday_work_sh(1:1), 1, jday0, 5, iv(1:1))
  
 call interp_spline_1d(a_sh(1:1), jday_work(1:1), 1, a0, jday0, 5, iv(1:1), 1., 0., 0.)
 call interp_spline_1d(b_sh(1:1), jday_work(1:1), 1, b0, jday0, 5, iv(1:1), 1., 0., 0.)
 call interp_spline_1d(c_sh(1:1), jday_work(1:1), 1, c0, jday0, 5, iv(1:1), 1., 0., 0.)
 call interp_spline_1d(d_sh(1:1), jday_work(1:1), 1, d0, jday0, 5, iv(1:1), 1., 0., 0.)
 call interp_spline_1d(dt_sh(1:1), jday_work(1:1), 1, dt0, jday0, 5, iv(1:1), 1., 0., 0.)
  
! Definition of soil temperature profile at each grid point
  
 do j=1,nlat
 do i=1,nlon
  
   if (fmask(i,j) < 0.5 ) then
  
     tsoil(i,j,:)=t_bottom(i,j)
  
     nbottom=index_lev_tbottom(i,j)
  
     do k=1,nbottom
  
! Scale of depth is the variable of analytic function
  
       x=(soil_lev(k)-soil_lev(1))/(soil_lev(nbottom)-soil_lev(1))
  
! Analytical function of temperature profile form as a function of depth scale
  
       if (alat (i,j) >= 0.) then
         xx=x*a(1)+b(1)
         func(k)=cos(xx)*c(1)+d(1)
       else
         xx=x*a_sh(1)+b_sh(1)
         func(k)=cos(xx)*c_sh(1)+d_sh(1)
       endif
  
     enddo
  
     if (alat(i,j) >= 0.) then
       dt_work=dt(1)
     else
       dt_work=dt_sh(1)
     endif
  
     t1=tsurf(i,j) ! Temperature at soil top level
  
     if (dt_work >= 0.) then
       if (t1 < t_bottom(i,j)) t1=t_bottom(i,j)
     else
       if (t1 > t_bottom(i,j)) t1=t_bottom(i,j)
     endif
  
     do k=1,nbottom-1
       tsoil(i,j,k)=func(k)*(t_bottom(i,j)-t1)/(func(nbottom)-func(1))+t_bottom(i,j)
     enddo

! Glacier
     if (soil_map(i,j,1) == nst-1) then
       do k=1,nlevel_soil
         tsoil(i,j,k)=min(tsoil(i,j,k), 272.)
       enddo
     endif

   else
  
     tsoil(i,j,:)=t_bottom(i,j)
     tsoil(i,j,1)=tsurf(i,j)
  
   endif
  
 enddo
 enddo
  
! Definition of relative and specific volumentric soil water content at all grid point at all soil levels
  
 do j=1,nlat
 do i=1,nlon

   if (fmask(i,j) < 0.5.and.int(soil_map(i,j,1)) < nst-1 ) then

! No Water body, no Glacier

     qrel(i,j,:)=qrel_bottom(i,j)
     qrel(i,j,1)=veg_qrel_surf(i,j)

     nbottom=index_lev_qbottom(i,j)

     do k=2,nbottom-1
       x=(soil_lev(k)-soil_lev(1))/(soil_lev(nbottom)-soil_lev(1))
       qrel(i,j,k)=qrel(i,j,nbottom)*x+qrel(i,j,1)*(1.-x)
       qrel(i,j,k)=min( max( qrel(i,j,k), 0.), 1.)
     enddo

     qsoil_rel(i,j,1:nbottom) = qrel(i,j,1:nbottom)

     qsoil_rel(i,j,nbottom+1:nlevel_soil)=qsoil_rel(i,j,nbottom)

   else

! Water body or Glacier

     qsoil_rel(i,j,:)=0.

   endif

 enddo
 enddo

return
end subroutine physiographic_param

!=======================================================================
! Group  of subroutines defining orography, soil and vegetation properties,
! to be included in preprocessing and nesting of BOLAM and MOLOCH.
!=======================================================================
subroutine orogdef(alon,alat,orog,orog_variance,alon0,alat0,dlon,dlat,nlon,nlat,flag_glob)

! Reads orography dataset from global lat-lon regular database with resolution
! 1/120 degree, derived from www.ngdc.noaa.gov/mgg/topo/globe.html.
! From this site 16 original binary files were downloaded. the 16 files were
! merged into one file written in the same format: 1 data value uses 16 bit.
! The merged files is orogr_global_latlon_1km.bin. It contains
! 21600 lines and 43200 columns, starting from the point -180 LON AND 90 LAT,
! corresponding to the North-West extreme of the NW-most
! database pixel.

! For correct reading of input database files, for ifort compiler it is necessary to specify
! "-assume byterecl" as compilation option

! Model grid orography is defined below using one of two methods, as a function of
! model horiz. resolution dlat:
! If dlat > 4.*dlat_read (orog. data resolution, 1/120 deg., equivalent to 0.033 deg.),
! then a 2-D average weighted with a Gaussian function is computed with parameters:
! radius of gaussian width (zradinfl) and a factor determining the area of gaussian
! filter application (zdistfr).
! If dlat < 4.*dlat_read, a bilinear interpolation is applied.

! In both cases the procedure computes the model grid orography std. deviation
! with different algorithms but always applying gaussian-weighted average.
! An "envelope orography" enhancement is applied using orog. std. deviation.

! Input:
! alon, alat:   model grid coordinates in geographical coordinates
! dx, dy:       model grid resolution (degrees)
! nlon, nlat:   number of model grid points
! flag_glob:    if 0, then not global input grid, if 1, then input grid is global not rotated

! Output:
! orog: orography on the model grid
! orog_variance: orography variance (later, std. dev.) on the model grid

implicit none

integer :: nlon, nlat, flag_glob
real :: dlon, dlat, alon0, alat0
real, dimension(nlon,nlat) :: orog, orog_variance, alon, alat

! Input dataset parameters:

real, parameter :: dlon_read=1./120., dlat_read=1./120.,                            &
                   lonini_read=-180.+dlon_read*0.5, lonfin_read=180.-dlon_read*0.5, &
                   latini_read=90.-dlon_read*0.5, latfin_read=-90.+dlon_read*0.5

! dlon_read, dlat_read: resolutions of orography databases
! lonini_read, lonfin_read: longitude of the west and east extremes of orography database
! latini_read, latfin_read: latitude of the north and south extremes of orography database

integer, parameter :: nx_read=43200, ny_read=(latini_read-latfin_read)/dlat_read


real, dimension(nlon,nlat) :: lat_geo, lon_geo, variance2, grad, wrk
integer*2, dimension(:,:), allocatable :: orog_read
integer*2, dimension(nx_read) :: orog_read0
real, dimension(:,:), allocatable :: zfunc
character(len=27) :: afile='orogr_global_latlon_1km.bin'

real :: zpi, val_missing=-9999, lon_geo_min, lon_geo_max, lat_geo_min, lat_geo_max, &
        zradinfl, zdistfr, lat_start_read, lon_start_read, zdistfac, zdistinfl,     &
        zdistmax, zparam, zlatfac, zlatnb, zlatsb, zlonwb, zloneb, zfuncsum,        &
        zlatbase, zlonbase, zdist, wei, dlat1, dlat2, dlon1, dlon2, zfy1, zfy2,     &
        zlatpoint, zlonpoint, orogaver, zlat, zlon, zehf

integer :: istart, ifinish, jstart, jfinish, nlon_read, nlat_read, jini, jfin, ix_part, &
           flag_period=0, jjj, iii, j, i, jlat, jlon, jnb, jsb, iwb, ieb, jj, ii,       &
           jl, il, nnit, nsmooth
integer, parameter :: nx_part=2
integer, dimension(nx_part) :: iini, ifin

INTEGER, DIMENSION(5) :: IDATE0=(/2018,9,20,0,0/)
INTEGER, DIMENSION(3) :: IPERIOD=(/0,0,0/)

 zpi=abs(acos(-1.))

 lat_geo(:,:) = alat(:,:)
 lon_geo(:,:) = alon(:,:)
 do jlat=1,nlat
 do jlon=1,nlon
   if (abs(abs(lat_geo(jlon,jlat))-90.) < 1.e-10) lat_geo(jlon,jlat) = sign(89.999999, lat_geo(jlon,jlat))
! Dataset longitude: -180..180
   if (alon(jlon,jlat) < -180.) lon_geo(jlon,jlat) = alon(jlon,jlat)+360.
   if (alon(jlon,jlat) >  180.) lon_geo(jlon,jlat) = alon(jlon,jlat)-360.
 enddo
 enddo
 if (flag_glob == 1.and.abs(alon0-lonini_read) > 10.) then
! Reset (below) dataset longitude: 0..360
   lon_geo(:,:) = alon(:,:)
 endif

! Original dataset definition and reading:

 lon_geo_min = minval(lon_geo)
 lon_geo_max = maxval(lon_geo)
 lat_geo_min = minval(lat_geo)
 lat_geo_max = maxval(lat_geo)

 istart = int(lon_geo_min*0.2)*5-10
 ifinish = int(lon_geo_max*0.2)*5+15
 jstart = max( int(lat_geo_min*0.2)*5-10, -90 )
 jfinish = min( int(lat_geo_max*0.2)*5+15, 90 )

! Gaussian filter parameters:
! Gaussian half width

 zradinfl = 1.5*dlat/(1.+5.*dlat)

! Factor determining radius of gaussian filter application (around a grid point).
! This radius is equal to ZRADINFL*ZDISTFR:

 zdistfr = 3.
 if (dlat > 0.2) zdistfr = 2.4 ! to speed up execution

 open(16,file=afile,status='old',form='unformatted',access='direct',recl=nx_read*2)

! Reading of the database

 print *,'Reading the global topography dataset...'

! Control of read data domain along latitude axis

 if (flag_glob == 0) then

   if (float(jstart) >= latfin_read.and.float(jfinish) <= latini_read) then

! Domain extremes located inside of latitude limits of dataset global domain

     jini = int((latini_read-float(jstart))/dlat_read)+1
     jfin = int((latini_read-float(jfinish))/dlat_read)+1

   else

! Domain extremes located outside of latitude limits of dataset global domain

     if (float(jstart) < latfin_read.and.float(jfinish) <= latini_read) then

       jini = ny_read
       jfin = int((latini_read-float(jfinish))/dlat_read)+1
       if (abs(jstart) /= 90) then
         print *,'Caution: the orography dataset does not cover the south part of model domain'
         print *,'Model grid south latitude extreme',float(jstart)
         print *,'Dataset grid south latitude extreme',latfin_read
       endif

     elseif (float(jfinish) > latini_read.and.float(jstart) >= latfin_read) then

       jini = int((latini_read-float(jstart))/dlat_read)+1
       jfin = 1
       if (abs(jfinish) /= 90) then
         print *,'Caution: the orography dataset does not cover the north part of model domain'
         print *,'Model grid north latitude extreme',float(jfinish)
         print *,'Dataset grid north latitude extreme',latini_read
       endif

     else

! Domain covers the whole latitude range

       jini=ny_read
       jfin=1
       if (abs(jstart) /= 90.or.abs(jfinish) /= 90) then
         print *,'Caution: the orography dataset does not cover the north and south parts of model domain'
         print *,'Model grid latitude extremes',float(jstart),float(jfinish)
         print *,'Dataset grid latitude extremes',latfin_read,latini_read
       endif

     endif

   endif

 else

! Domain covers the whole latitude range

   jini=ny_read
   jfin=1

 endif

 lat_start_read = latini_read-float(jini-1)*dlat_read

! Control of read data domain along longitude axis, possible cut and paste operation for periodic grid

 if (flag_glob == 0) then

   if (istart >= int(lonini_read).and.ifinish <= int(lonfin_read)) then

! Domain extremes locate inside of longitude limits of dataset global domain

     iini(1) = int((float(istart)-lonini_read)/dlon_read)+1
     ifin(1) = int((float(ifinish)-lonini_read)/dlon_read)+1

     iini(2) = 0
     ifin(2) = -1

     flag_period = 0

   else

! Limited domain locates outside of longitude limits of dataset global domain -> cut and paste

     if (istart < int(lonini_read).and.ifinish <= int(lonfin_read)) then

       iini(1) = int((float(istart)+360.-lonini_read)/dlon_read)+1
       ifin(1) = nx_read

       iini(2) = 1
       ifin(2) = int((float(ifinish)-lonini_read)/dlon_read)+1

       flag_period = 1

     elseif (ifinish > int(lonfin_read).and.istart >= int(lonini_read)) then

       iini(1) = int((float(istart)-lonini_read)/dlon_read)+1
       ifin(1) = nx_read

       iini(2) = 1
       ifin(2) = int((float(ifinish)-360.-lonini_read)/dlon_read)+1

       flag_period = 1

     else

! Domain covers the whole longitude circle

       iini(1) = 1
       ifin(1) = nx_read

       iini(2) = 0
       ifin(2) = -1

       flag_period = 2

     endif

   endif

 else ! global output grid

! Domain covers the whole longitude circle

   if (abs(alon0-lonini_read) < 10.) then

! Input (dataset) and output (model) global grid have the same longitude origin

     iini(1) = 1
     ifin(1) = nx_read

     iini(2) = 0
     ifin(2) = -1

     flag_period = 2

   else

! Input (dataset) and output (model) global grid have the longitude origin shited by 180 degees

     iini(1) = int((alon0-lonini_read)/dlon_read)+1
     ifin(1) = nx_read

     iini(2) = 1
     ifin(2) = iini(1)-1

     flag_period = 3

   endif

 endif

 lon_start_read = lonini_read+float(iini(1)-1)*dlon_read

 nlat_read = jini-jfin+1
 nlon_read = ifin(1)-iini(1)+1
 if (flag_period == 1 ) then
   nlon_read = nlon_read+ifin(2)-iini(2)+1
 endif
 if (flag_period == 2.or.flag_period == 3 ) then
   nlon_read = nx_read+2
   if (flag_period == 2 ) lon_start_read = lonini_read-dlon_read
   if (flag_period == 3 ) lon_start_read = lon_start_read-dlon_read
 endif

 allocate(orog_read(nlon_read,nlat_read))
 orog_read(:,:) = int(val_missing)

 j=0
 do jjj=jini, jfin, -1
   j = j+1

   read (16, rec=jjj) orog_read0

   i=0
   if (flag_period == 2.or.flag_period == 3 ) i=1
   do ix_part=1,nx_part
     do iii=iini(ix_part), ifin(ix_part)
       i = i+1
       orog_read(i,j) = orog_read0(iii)
       if (orog_read(i,j) == -500) orog_read(i,j) = 0 ! Sea
     enddo
   enddo ! ix_part
   if (flag_period == 2.or.flag_period == 3 ) then
     orog_read(1,j) = orog_read(nlon_read-1,j)
     orog_read(nlon_read,j) = orog_read(2,j)
   endif

 enddo ! jjj

 close (16)

! Definition of model rotated grid orography ZOROG from original
! database grid (lat-lon regular) IOROGBASE

 zdistfac = 6371.e3*zpi/180.

 if (flag_period == 1) then
   do jlat=1,nlat
   do jlon=1,nlon
     if (lon_geo(jlon,jlat) < 0.) lon_geo(jlon,jlat) = lon_geo(jlon,jlat)+360.
   enddo
   enddo
 endif

 if (dlat > dlat_read*4.) then

! Average weight (with gaussian function) calculation:

   print*, "Model orography computed with gaussian filter"
   print*, "Gaussian half-width (deg.)", zradinfl

   zdistinfl = zradinfl*zdistfac
   zdistmax = zradinfl*zdistfr
   zparam = 1./zdistinfl

   do jlat=1,nlat
   do jlon=1,nlon

     zlatfac = cos(max( min( lat_geo(jlon,jlat), 89.), -89.)*zpi/180.)

     zlatsb = lat_geo(jlon,jlat)-zdistmax
     zlatnb = lat_geo(jlon,jlat)+zdistmax
     zlonwb = lon_geo(jlon,jlat)-zdistmax/zlatfac
     zloneb = lon_geo(jlon,jlat)+zdistmax/zlatfac

!     jsb = int((zlatsb-lat_start_read)/dlat_read)+1   ! old
!     jnb = int((zlatnb-lat_start_read)/dlat_read)+1
!     iwb = int((zlonwb-lon_start_read)/dlon_read)+1
!     ieb = int((zloneb-lon_start_read)/dlon_read)+1

     jsb = nint((zlatsb-lat_start_read)/dlat_read)+1  ! Nov. 2014
     jnb = nint((zlatnb-lat_start_read)/dlat_read)+1  !   "
     iwb = nint((zlonwb-lon_start_read)/dlon_read)+1  !   "
     ieb = nint((zloneb-lon_start_read)/dlon_read)+1  !   "

     jsb = max(min(jsb, nlat_read), 1)
     jnb = max(min(jnb, nlat_read), 1)
     iwb = max(min(iwb, nlon_read), 1)
     ieb = max(min(ieb, nlon_read), 1)

     allocate(zfunc((ieb-iwb+10),(jnb-jsb+10)))

     orog(jlon,jlat) = 0.
     zfuncsum = 0.

     jl = 0
     do jj=jsb,jnb
       jl = jl+1
       il = 0
       do ii=iwb,ieb
         il = il+1
         zlatbase = lat_start_read+dlat_read*float(jj-1)
         zlonbase = lon_start_read+dlon_read*float(ii-1)
         zdist = sqrt(((lon_geo(jlon,jlat)-zlonbase)*zlatfac*zdistfac)**2+((lat_geo(jlon,jlat)- &
         zlatbase)*zdistfac)**2)
         zfunc(il,jl) = exp(-(zdist*zparam)**2)
         orog(jlon,jlat) = orog(jlon,jlat)+float(orog_read(ii,jj))*zfunc(il,jl)
         zfuncsum = zfuncsum+zfunc(il,jl)
       enddo
     enddo

     orog(jlon,jlat) = orog(jlon,jlat)/zfuncsum

! Orography variance calculation:

     orog_variance(jlon,jlat) = 0.
     jl = 0
     do jj=jsb,jnb
       jl = jl+1
       il = 0
       do ii=iwb,ieb
         il = il+1
         orog_variance(jlon,jlat) = orog_variance(jlon,jlat)+(float(orog_read(ii,jj))-orog(jlon,jlat))**2* &
                      zfunc(il,jl)
       enddo
     enddo

! Standard deviation computed from variance

     orog_variance(jlon,jlat) = sqrt(orog_variance(jlon,jlat)/zfuncsum)

     deallocate(zfunc)

   enddo
   enddo

   do jlat=2,nlat-1
   do jlon=2,nlon-1
     grad(jlon,jlat) = sqrt((orog(jlon+1,jlat)-orog(jlon-1,jlat))**2+(orog(jlon,jlat+1)-  &
                       orog(jlon,jlat-1))**2)
   enddo
   enddo
   do jlon=2,nlon-1
     grad(jlon,1) = grad(jlon,2)
     grad(jlon,nlat) = grad(jlon,nlat-1)
   enddo
   do jlat=1,nlat
     grad(1,jlat) = grad(2,jlat)
     grad(nlon,jlat) = grad(nlon-1,jlat)
   enddo

! Correction of orogr. std. deviation with orog. gradient to reduce std. dev. due to cliffs

   do jlat=1,nlat
   do jlon=1,nlon
     variance2(jlon,jlat) = orog_variance(jlon,jlat)/(0.00075*grad(jlon,jlat)+1.)
   enddo
   enddo

   nnit = max(15 - nint(33.*dlat), 6) ! funct. of resol.
   print*, "No. of filter iterations for orography and orogr. std. dev.", nnit
   call hd4(orog_variance,nlon,nlat,dlon,dlat,nnit)
   call hd4(variance2,nlon,nlat,dlon,dlat,nnit)

! "Envelope orography" definition using modified std. deviation
!  and further weak modif. as a funct. of resolution

   do jlat=1,nlat
   do jlon=1,nlon
     orog(jlon,jlat) = orog(jlon,jlat)+0.7*variance2(jlon,jlat)
     orog(jlon,jlat) = orog(jlon,jlat)*(1. + 0.08*dlat)
     orog_variance(jlon,jlat) = orog_variance(jlon,jlat)*.7 ! (only for roughness def.)
   enddo
   enddo

! Smoothing of orography field (fourth order and weak second order)

   call hd4(orog,nlon,nlat,dlon,dlat,nnit)
   wei = 0.15
   nsmooth = 1
   call smooth_soil(orog,wrk,nlon,nlat,wei,nsmooth)

! Further smoothing of std. dev. to be used for roughness comput.

   nsmooth = nnit/5
   print*, "No. of 2nd filter iterations for orography std. dev.", NSMOOTH
   wei = 0.5
   call smooth_soil(orog_variance,wrk,nlon,nlat,wei,nsmooth)

 else ! Case dlat <= dlat_read*4.

! Case of bi-linear interpolation:

   print*, "Model orography computed by interpolation"

   do jlat=1,nlat
   do jlon=1,nlon

!     jj = int((lat_geo(jlon,jlat)-lat_start_read)/dlat_read)+1  ! old
!     ii = int((lon_geo(jlon,jlat)-lon_start_read)/dlon_read)+1

     jj = nint((lat_geo(jlon,jlat)-lat_start_read)/dlat_read)  ! Nov. 2014
     ii = nint((lon_geo(jlon,jlat)-lon_start_read)/dlon_read)

     ii = max(min(ii, nlon_read), 1)  ! (to avoid problem at the South Pole)
     jj = max(min(jj, nlat_read), 1)  !         "                  "

     dlat1 = lat_geo(jlon,jlat)-lat_start_read+dlat_read*float(jj-1)
     dlat2 = lat_start_read+dlat_read*float(jj)-lat_geo(jlon,jlat)
     dlon1 = lon_geo(jlon,jlat)-lon_start_read+dlon_read*float(ii-1)
     dlon2 = lon_start_read+dlon_read*float(ii)-lon_geo(jlon,jlat)

     zfy1 = (float(orog_read(ii,jj  ))*dlon2+float(orog_read(ii+1,jj  ))*dlon1)/(dlon1+dlon2)
     zfy2 = (float(orog_read(ii,jj+1))*dlon2+float(orog_read(ii+1,jj+1))*dlon1)/(dlon1+dlon2)
     orog(jlon,jlat) = (zfy1*dlat2+zfy2*dlat1)/(dlat1+dlat2)
   enddo
   enddo

! Orography variance calculation:

   zdistinfl = 2.*dlat
   zdistmax = zdistinfl*zdistfr
   zparam = 1./zdistinfl

   do jlat=1,nlat
   do jlon=1,nlon

     zlatpoint = alat0+dlat*float(jlat-1)
     zlonpoint = alon0+dlon*float(jlon-1)

     zlatfac = cos(max( min( zlatpoint, 89.), -89.)*zpi/180.)

     zlatsb = alat0+dlat*float(jlat-1)-zdistmax
     zlatnb = alat0+dlat*float(jlat-1)+zdistmax
     zlonwb = alon0+dlon*float(jlon-1)-zdistmax/zlatfac
     zloneb = alon0+dlon*float(jlon-1)+zdistmax/zlatfac

     jsb = max(int((zlatsb-alat0)/dlat)+1, 1)
     jnb = min(int((zlatnb-alat0)/dlat)+1, nlat)
     iwb = max(int((zlonwb-alon0)/dlon)+1, 1)
     ieb = min(int((zloneb-alon0)/dlon)+1, nlon)

     allocate(zfunc((ieb-iwb+10),(jnb-jsb+10)))

     orogaver = 0.
     zfuncsum = 0.
     jl = 0
       do jj=jsb,jnb
       jl = jl+1
       il = 0
       do ii=iwb,ieb
         il = il+1
         zlat = alat0+dlat*float(jj-1)
         zlon = alon0+dlon*float(ii-1)
         zdist = sqrt((zlonpoint-zlon)**2+(zlatpoint-zlat)**2)
         zfunc(il,jl) = exp(-(zdist*zparam)**2)
         orogaver = orogaver+orog(ii,jj)*zfunc(il,jl)
         zfuncsum = zfuncsum+zfunc(il,jl)
       enddo
     enddo

     orogaver = orogaver/zfuncsum

     orog_variance(jlon,jlat) = 0.
     jl = 0
     do jj=jsb,jnb
       jl = jl+1
       il = 0
       do ii=iwb,ieb
         il = il+1
         orog_variance(jlon,jlat) = orog_variance(jlon,jlat)+(orog(ii,jj)-orogaver)**2*zfunc(il,jl)
       enddo
     enddo
     orog_variance(jlon,jlat) = sqrt(orog_variance(jlon,jlat)/zfuncsum)

     deallocate(zfunc)

   enddo
   enddo

! Orography (small) enhancement

   zehf = 1. + 0.6*dlat
   print*, "Enhancement of filtered orography", zehf
   orog(:,:) = orog(:,:)*zehf

! Smoothing of interpolated orography

   nnit = max(32 - nint(500.*dlat), 6)      ! funct. of resol.
   print*, "Iterations for filtering orography", nnit
   call hd4(orog,nlon,nlat,dlon,dlat,nnit)

   wei = 0.5
   nsmooth = max(nnit/9, 1)
   print*, "NSMOOTH for orography", nsmooth
   call smooth_soil(orog,wrk,nlon,nlat,wei,nsmooth)

! Rescaling and smoothing of orog. std. deviation, used to define roughness

   orog_variance(:,:) = orog_variance(:,:)*0.65

   wei = 0.5
   nsmooth = max(nnit/7, 1)
   print*, "NSMOOTH for orog. st. dev.", nsmooth
   call smooth_soil(orog_variance,wrk,nlon,nlat,wei,nsmooth)

 endif ! Case dlat

 orog_variance(:,:) = max(orog_variance(:,:), 0.)

return
end
!=======================================================================
subroutine dataset_soil_fao(nlon, nlat, dlon, dlat, alon, alat, flag_glob, soilmap)

! Procedure defines soil types (classification FAO) over the model grid

! Internal basic variables:

!      nst: total number soil types;

! Input variables:

!      nlon: grid point number along axis of rotated longitude;
!      nlat: grid point number along axis of rotated latidude;
!      dlon: grid step along the domain axis of rotated longitude (degree);
!      dlat: grid step along the domain axis of rotated latitude (degree);
!      alon: geographical longitude of the grid points (degree, -180...180);
!      alat: geographical latitude of the grid points (degree, -90...90);
!      flag_glob:    if 0, then not global input grid, if 1, then input grid is global not rotated

! Output variables:

!      soilmap:  soil types (classification FAO, 31 types) defined on the grid points, 3D array,
!                1-st and 2-nd array index are grid point index, 3-rd index is soil type number plus 1,
!                if 3-rd index has value 1, then the variable means dominate soil type,
!                values 2-32 of 3-rd index means fraction (of area, proportion) of the each of 31 soil types;

! Soil types are defined starting from the global soil dataset in file:
!             soil_fao_global_latlon_8km.bin
! and legend (descripting) file:
!             worldexp.dat

! The global soil dataset (with resolution 1/12 degree) was obtained from the global
! soil map of FAO and corresponding global dataset at www.fao.org,
! dataset grid: lon -180 --> 180, lat 90 --> -90

implicit none

integer, parameter :: nst=31 ! NST must be 31!

integer :: nlon, nlat, flag_glob
real :: dlat, dlon
real, dimension(nlon,nlat) ::  alat, alon
real, dimension(nlon,nlat,nst+1) :: soilmap

! Work variables:

real, parameter :: dlon_read_s=1./12., dlat_read_s=1./12., &
                   lonini_read_s=-180.+dlon_read_s*0.5, lonfin_read_s=180.-dlon_read_s*0.5, &
                   latini_read_s=90.-dlon_read_s*0.5, latfin_read_s=-90.+dlon_read_s*0.5

! dlon_read_s, dlat_read_s, dlon_read_v, dlat_read_v: resolutions of soil and landcover databases
! lonini_read_s, lonfin_read_s: longitude of the west and east extremes of soil database
! latini_read_s, latfin_read_s: latitude of the north and south extremes of soil database

integer, parameter :: nx_read_s=4320, ny_read_s=(latini_read_s-latfin_read_s)/dlat_read_s

integer, parameter :: nsnum=7000, nsoil=8 ! 7000 is the number of soil units

real, dimension(nlon,nlat) :: lat_geo, lon_geo
integer*2, dimension(:,:), allocatable :: soil_read
integer*2 :: soil_read0(nx_read_s)
character (len=50) :: file_input='soil_fao_global_latlon_8km.bin'
real, dimension(nsoil,nsnum) :: fsoil
real, dimension(nst) :: ffsoil
integer, dimension(nsoil,nsnum) :: soil_type
integer, dimension(nsnum) :: isnum
character(len=2), dimension(nsoil,nsnum) :: asoil
real :: pi, val_missing=-9999., lat_geo_min, lat_geo_max, lon_geo_min, lon_geo_max, &
        lat_start_read_s, lon_start_read_s, zfmax
integer :: nlon_read_s, nlat_read_s, istart, ifinish, jstart, jfinish, jini, jfin, j, i, jjj, iii, &
           flag_period_s=0, ix_part, isoil, is, jsoil, js, jsmax, insoil, jlon, jlat
real, dimension(4) :: zdist
integer, parameter :: nx_part = 2
integer, dimension(nx_part) :: iini, ifin
integer, dimension(1) :: iind

type work_format
  character(len=2) :: asoil
  real :: fsoil
  real, dimension(10) :: zzz
end type work_format

type soil_unit_descript
  integer :: iunit
  character(len=13) :: descript1
  character(len=7) :: descript2
  type (work_format), dimension(8) :: read
end type soil_unit_descript

type (soil_unit_descript), dimension(nsnum) :: soil_unit

 pi = abs(acos(-1.))

 lat_geo(:,:) = alat(:,:)
 lon_geo(:,:) = alon(:,:)
 do jlat=1,nlat
 do jlon=1,nlon
   if (abs(abs(lat_geo(jlon,jlat))-90.) < 1.e-10) lat_geo(jlon,jlat) = sign(89.999999, lat_geo(jlon,jlat))
! Dataset longitude: -180..180
   if (alon(jlon,jlat) < -180.) lon_geo(jlon,jlat) = alon(jlon,jlat)+360.
   if (alon(jlon,jlat) >  180.) lon_geo(jlon,jlat) = alon(jlon,jlat)-360.
 enddo
 enddo
 if (flag_glob == 1.and.abs(alon(1,1)-lonini_read_s) > 10.) then
! Reset (below) dataset longitude: 0..360
   lon_geo(:,:) = alon(:,:)
 endif

! Coordinates of the centre of the extremes for reading

 lon_geo_min = minval(lon_geo)
 lon_geo_max = maxval(lon_geo)
 lat_geo_min = minval(lat_geo)
 lat_geo_max = maxval(lat_geo)

 istart = int(lon_geo_min*0.2)*5-10
 ifinish = int(lon_geo_max*0.2)*5+15
 jstart = max( int(lat_geo_min*0.2)*5-10, -90 )
 jfinish = min( int(lat_geo_max*0.2)*5+15, 90 )

! Reading of soil original dataset

 print *
 print *,'Reading the global soil dataset from ',file_input

 open(16,file=file_input,status='old',form='unformatted',access='direct',recl=nx_read_s*2)

! Control of read data domain along latitude axis

 if (flag_glob == 0) then

   if (float(jstart) >= latfin_read_s.and.float(jfinish) <= latini_read_s) then

! Domain extremes locate inside of latitude limits of dataset global domain

     jini = int((latini_read_s-float(jstart))/dlat_read_s)+1
     jfin = int((latini_read_s-float(jfinish))/dlat_read_s)+1

   else

! Domain extremes locate outside of latitude limits of dataset global domain

     if (float(jstart) < latfin_read_s.and.float(jfinish) <= latini_read_s) then

       jini = ny_read_s
       jfin = int((latini_read_s-float(jfinish))/dlat_read_s)+1
       if (abs(jstart) /= 90) then
         print *,' Caution: the soil dataset does not cover the south part of the model domain'
         print *,' Model grid south latitude extreme',float(jstart)
         print *,' Dataset grid south latitude extreme',latfin_read_s
       endif

     elseif (float(jfinish) > latini_read_s.and.float(jstart) >= latfin_read_s) then

       jini = int((latini_read_s-float(jstart))/dlat_read_s)+1
       jfin = 1
       if (abs(jfinish) /= 90) then
         print *,' Caution: the soil dataset does not cover the north part of the model domain'
         print *,' Model grid north latitude extreme',float(jfinish)
         print *,' Dataset grid north latitude extreme',latini_read_s
       endif

     else

! Domain covers the whole latitude range

       jini = ny_read_s
       jfin = 1
       if (abs(jstart) /= 90.or.abs(jfinish) /= 90) then
         print *,' Caution: the soil dataset does not cover the north and south parts of the model domain'
         print *,' Model grid latitude extremes',float(jstart),float(jfinish)
         print *,' Dataset grid latitude extremes',latfin_read_s,latini_read_s
       endif

     endif

   endif

 else

! Domain covers the whole latitude range

   jini=ny_read_s
   jfin=1

   print *,' Caution: the soil dataset does not cover the north and south parts of the model domain'
   print *,' Model grid latitude extremes',lat_geo(1,1),lat_geo(1,nlat)
   print *,' Dataset grid latitude extremes',latfin_read_s,latini_read_s

 endif

 lat_start_read_s = latini_read_s-float(jini-1)*dlat_read_s

! Control of read data domain along longitude axis, possible cut and paste operation for periodic grid

 if (flag_glob == 0) then

   if (istart >= int(lonini_read_s).and.ifinish <= int(lonfin_read_s)) then

! Domain extremes located inside of longitude limits of dataset global domain

     iini(1) = int((float(istart)-lonini_read_s)/dlon_read_s)+1
     ifin(1) = int((float(ifinish)-lonini_read_s)/dlon_read_s)+1

     iini(2)=0
     ifin(2)=-1

     flag_period_s = 0

   else

! Domain extremes located outside of longitude limits of dataset global domain -> cut and paste

     if (istart < int(lonini_read_s).and.ifinish <= int(lonfin_read_s)) then

       iini(1) = int((float(istart)+360.-lonini_read_s)/dlon_read_s)+1
       ifin(1) = nx_read_s

       iini(2) = 1
       ifin(2) = int((float(ifinish)-lonini_read_s)/dlon_read_s)+1

       flag_period_s = 1

     elseif (ifinish > int(lonfin_read_s).and.istart >= int(lonini_read_s)) then

       iini(1) = int((float(istart)-lonini_read_s)/dlon_read_s)+1
       ifin(1) = nx_read_s

       iini(2) = 1
       ifin(2) = int((float(ifinish)-360.-lonini_read_s)/dlon_read_s)+1

       flag_period_s = 1

     else

! Domain covers the whole longitude circle

       iini(1) = 1
       ifin(1) = nx_read_s

       iini(2) = 0
       ifin(2) = -1

       flag_period_s = 2

     endif

   endif

 else ! global output grid

! Domain covers the whole longitude circle

   if (abs(alon(1,1)-lonini_read_s) < 10.) then

! Input (dataset) and output (model) global grid have the same longitude origin

     iini(1) = 1
     ifin(1) = nx_read_s

     iini(2) = 0
     ifin(2) = -1

     flag_period_s = 2

   else

! Input (dataset) and output (model) global grid have the longitude origin shited by 180 degees

     iini(1) = int((alon(1,1)-lonini_read_s)/dlon_read_s)+1
     ifin(1) = nx_read_s

     iini(2) = 1
     ifin(2) = iini(1)-1

     flag_period_s = 3

   endif

 endif

 lon_start_read_s = lonini_read_s+float(iini(1)-1)*dlon_read_s

 nlat_read_s = jini-jfin+1
 nlon_read_s = ifin(1)-iini(1)+1
 if (flag_period_s == 1 ) then
   nlon_read_s = nlon_read_s+ifin(2)-iini(2)+1
 endif
 if (flag_period_s == 2.or.flag_period_s == 3) then
   nlon_read_s = nx_read_s+2
   if (flag_period_s == 2) lon_start_read_s = lonini_read_s-dlon_read_s
   if (flag_period_s == 3) lon_start_read_s = lon_start_read_s-dlon_read_s
 endif

 allocate(soil_read(nlon_read_s,nlat_read_s))
 soil_read(:,:) = int(val_missing)

 j = 0
 do jjj=jini, jfin, -1
   j = j+1

   read (16, rec=jjj) soil_read0

   i = 0
   if (flag_period_s == 2.or.flag_period_s == 3 ) i = 1
   do ix_part=1,nx_part
     do iii=iini(ix_part), ifin(ix_part)
       i = i+1
       soil_read(i,j) = soil_read0(iii)
     enddo
   enddo ! ix_part
   if (flag_period_s == 2.or.flag_period_s == 3 ) then
     soil_read(1,j) = soil_read(nlon_read_s-1,j)
     soil_read(nlon_read_s,j) = soil_read(2,j)
   endif

 enddo ! jjj

 close (16)

 print *,'Soil dataset processing'

! "Interpolation" of soil types from the original soil map grid (ISOLBAS)
! to the model grid (SOILMAP)

! The soil type is defined by a code number from 1 to 6998 (prescribed by FAO)
! in the 2-D grid vector 'ISOLBAS'
! The interpretation of the code value is obtained from the auxiliary file
! 'worldexp.dat' that defines the 6988 codes: each code is associated with a min. of 1
! and a max. of 8 soil types chosen among 31 types of the FAO classification.
! The usual FAO types are 26. Additional types are:
! 'GL': glaciers; 'DS': dunes or drifting sand; 'RK': rock debris or desert detritus;
! 'ST': salt flats; 'WR' OR 'ND': water bodies (sea, lakes, rivers).

! Total 31 FAO soil types:

!  1 Acrisols
!  2 Cambisols
!  3 Chernozems
!  4 Podzoluvisols
!  5 Rendzinans
!  6 Ferrasols
!  7 Gleysols
!  8 Phaeozems 
!  9 Lithosols
! 10 Fluvisols
! 11 Kastanozems
! 12 Luvisols
! 13 Greyzems
! 14 Nitosols
! 15 Histosols
! 16 Podzols
! 17 Arenosols
! 18 Regosols
! 19 Solonetz
! 20 Andosols
! 21 Rankers
! 22 Vertisols
! 23 Planosols
! 24 Xerosols
! 25 Yermosols
! 26 Solonchaks
! 27 Dunes or shifting sands
! 28 Salt flats
! 29 Rock debris or desert detritus
! 30 Glaciers
! 31 Water body

! Reading soil unit description from worldexp.dat
! (the format is complicate but follows the original setup, except that
! commas have been added for fortran reading)

 isnum(1:nsnum) = 0

 open (18,file='worldexp.dat',status='old')

 do i=1,nsnum
   read (18,'(i4,1x,a,1x,a,8(1x,a,1x,f8.0,10(1x,f7.0)))',end=211) soil_unit(i)
   isnum(i) = soil_unit(i)%iunit
   do j=1,nsoil
     asoil(j,isnum(i)) = soil_unit(i)%read(j)%asoil
     fsoil(j,isnum(i)) = soil_unit(i)%read(j)%fsoil
   enddo
 enddo ! i=1,nsnum

 211  close (18)

! Recoding soil name: from FAO types to our codes (types from 1 to 15)

 do iii=1,nsnum
   if (isnum(iii) /= 0) then
     i = isnum(iii)
     do j=1,nsoil
       if (asoil(j,i) /= '  ') then
         if (asoil(j,i)(1:1) == 'A') soil_type(j,i) = 1
         if (asoil(j,i)(1:1) == 'B') soil_type(j,i) = 2
         if (asoil(j,i)(1:1) == 'C') soil_type(j,i) = 3
         if (asoil(j,i)(1:1) == 'D') soil_type(j,i) = 4
         if (asoil(j,i)(1:1) == 'E') soil_type(j,i) = 5
         if (asoil(j,i)(1:1) == 'F') soil_type(j,i) = 6
         if (asoil(j,i)(1:1) == 'G') soil_type(j,i) = 7
         if (asoil(j,i)(1:1) == 'H') soil_type(j,i) = 8
         if (asoil(j,i)(1:1) == 'I') soil_type(j,i) = 9
         if (asoil(j,i)(1:1) == 'J') soil_type(j,i) = 10
         if (asoil(j,i)(1:1) == 'K') soil_type(j,i) = 11
         if (asoil(j,i)(1:1) == 'L') soil_type(j,i) = 12
         if (asoil(j,i)(1:1) == 'M') soil_type(j,i) = 13
         if (asoil(j,i)(1:1) == 'N') soil_type(j,i) = 14
         if (asoil(j,i)(1:1) == 'O') soil_type(j,i) = 15
         if (asoil(j,i)(1:1) == 'P') soil_type(j,i) = 16
         if (asoil(j,i)(1:1) == 'Q') soil_type(j,i) = 17
         if (asoil(j,i)(1:1) == 'R') soil_type(j,i) = 18
         if (asoil(j,i)(1:1) == 'S') soil_type(j,i) = 19
         if (asoil(j,i)(1:1) == 'T') soil_type(j,i) = 20
         if (asoil(j,i)(1:1) == 'U') soil_type(j,i) = 21
         if (asoil(j,i)(1:1) == 'V') soil_type(j,i) = 22
         if (asoil(j,i)(1:1) == 'W') soil_type(j,i) = 23
         if (asoil(j,i)(1:1) == 'X') soil_type(j,i) = 24
         if (asoil(j,i)(1:1) == 'Y') soil_type(j,i) = 25
         if (asoil(j,i)(1:1) == 'Z') soil_type(j,i) = 26
         if (asoil(j,i)(1:2) == 'DS') soil_type(j,i) = 27
         if (asoil(j,i)(1:2) == 'ST') soil_type(j,i) = 28
         if (asoil(j,i)(1:2) == 'RK') soil_type(j,i) = 29
         if (asoil(j,i)(1:2) == 'GL') soil_type(j,i) = 30
         if (asoil(j,i)(1:2) == 'WR') soil_type(j,i) = nst ! nst must be the last!
         if (asoil(j,i)(1:2) == 'ND') soil_type(j,i) = nst
       else
         soil_type(j,i) = 0
       endif
     enddo
   endif
 enddo     ! i=1,nsnum

 if (flag_period_s == 1) then
   do jlat=1,nlat
   do jlon=1,nlon
     if (lon_geo(jlon,jlat) < 0.) lon_geo(jlon,jlat) = lon_geo(jlon,jlat)+360.
   enddo
   enddo
 endif

 do j=1,nlat
 do i=1,nlon

! Search of the nearest database pixel (by pixel centre coordinate):

!   isoil = int((lon_geo(i,j)-lon_start_read_s)/dlon_read_s)+1   ! old
!   jsoil = int((lat_geo(i,j)-lat_start_read_s)/dlat_read_s)+1

!   isoil = nint((lon_geo(i,j)-lon_start_read_s)/dlon_read_s)+1   ! Oxana Nov. 2014
!   jsoil = nint((lat_geo(i,j)-lat_start_read_s)/dlat_read_s)+1

   isoil = nint((lon_geo(i,j)-lon_start_read_s)/dlon_read_s + 0.25)+1 ! Nov. 2014 (Europe only?)
   jsoil = nint((lat_geo(i,j)-lat_start_read_s)/dlat_read_s + 0.25)   ! Nov. 2014 (Europe only?)

   zdist(1) = (lon_start_read_s+dlon_read_s*float(isoil-1)-lon_geo(i,j))**2 +  &
              (lat_start_read_s+dlat_read_s*float(jsoil-1)-lat_geo(i,j))**2
   zdist(2) = (lon_start_read_s+dlon_read_s*float(isoil)  -lon_geo(i,j))**2 +  &
              (lat_start_read_s+dlat_read_s*float(jsoil-1)-lat_geo(i,j))**2
   zdist(3) = (lon_start_read_s+dlon_read_s*float(isoil)  -lon_geo(i,j))**2 +  &
              (lat_start_read_s+dlat_read_s*float(jsoil)  -lat_geo(i,j))**2
   zdist(4) = (lon_start_read_s+dlon_read_s*float(isoil-1)-lon_geo(i,j))**2 +  &
              (lat_start_read_s+dlat_read_s*float(jsoil)  -lat_geo(i,j))**2

   iind = minloc(zdist)

   if (iind(1) == 2.or.iind(1) == 3) isoil = isoil+1
   if (iind(1) == 3.or.iind(1) == 4) jsoil = jsoil+1

   isoil = max(min(isoil, nlon_read_s), 1)
   jsoil = max(min(jsoil, nlat_read_s), 1)

! Definition of 3-d vector SOILMAP, determining soil types (dominant type
! and fractions of all types in order 1-15)
! A geographical 'interpolation' is made from the soil_read grid (FAO cells):
! model grid points in the same FAO cell assume the same soil types.

! Code 0: missing value code

   is = soil_read(isoil,jsoil)
   if (is==0) is = 6997! water body
   do js=1,nst
     ffsoil(js) = 0.
   enddo
   do js=1,nst
     do jjj=1,nsoil
       insoil = soil_type(jjj,is)
       if (insoil == nst+1) insoil = nst
       if (insoil == js) ffsoil(js) = ffsoil(js)+fsoil(jjj,is)*0.01
     enddo
   enddo
   do js=1,nst
     soilmap(i,j,js+1) = ffsoil(js)
   enddo
   zfmax = 0.
   do js=1,nst
     if (ffsoil(js) > zfmax) then
       jsmax = js
       zfmax = ffsoil(js)
     endif
   enddo
   soilmap(i,j,1) = float(jsmax)

 enddo
 enddo

! End of soil type definition

!------------------------------------------------------------

! Antarctica in soilmap:

 do jlat=1,nlat
 do jlon=1,nlon
   if (lat_geo(jlon,jlat) < -60.and.nint(soilmap(jlon,jlat,1)) /= nst) then
     soilmap(jlon,jlat,1) = nst-1 ! glacier
     soilmap(jlon,jlat,2:nst+1) = 0.
     soilmap(jlon,jlat,nst) = 1.
   endif
 enddo
 enddo

return
end subroutine dataset_soil_fao
!=======================================================================
subroutine dataset_vegetation(nlon, nlat, dlon, dlat, alon, alat, flag_glob, vegmap, fmask)

! Procedure defines vegetation (landuse) types (classification GLC1990) over the model grid, 
! and land-sea fraction over the model grid using landuse dataset

! Internal basic variables:

!      nvt: total number vegetaion types;

! Input variables:

!      nlon: grid point number along axis of rotated longitude;
!      nlat: grid point number along axis of rotated latidude;
!      dlon: grid step along the domain axis of rotated longitude (degree);
!      dlat: grid step along the domain axis of rotated latitude (degree);
!      alon: geographical longitude of the grid points (degree, -180...180);
!      alat: geographical latitude of the grid points (degree, -90...90),
! flag_glob: if 0, then not global input grid, if 1, then input grid is global not rotated

! Output variables:

!      vegmap: vegetation (landuse) types (classification GLC1990, 14 types) defined on the grid points, 3D array,
!               1-st and 2-nd array index are grid point index, 3-rd index is vegetation types number plus 1,
!               if 3-rd index has value 1, then the variable means dominate vegetation type,
!               values 2-15 of 3-rd index means fraction (of area, proportion) of the each of 14 vegetation types
!       fmask: land-sea franction: 1-sea, 0-land (proportion).

! Vegetation types are defined starting from the global dataset in file:
!             veget_global_latlon_1km.bin 

! The global vegetation cover dataset (with resolution 1/120 degree) was obtained
! (with no modifications) from http://glcf.umd.edu/ following the links:
! Data & Products ->
!   AVHRR: Land Cover Classification ->
!     Download via web page with links to FTP Server ->
!       Global Coverage Data Sets: 1 Kilometer pixel resolution: Goodes projection, binary (BSQ)
! or directly from:
! ftp://ftp.glcf.umd.edu/glcf/Global_Land_Cover/Global/1deg/gl-latlong-1deg-landcover.bsq.gz
! Dataset grid: lon -180 --> 180, lat 90 --> -90

implicit none

integer, parameter :: nvt=14 ! NVT must be 14!

integer :: nlon, nlat, flag_glob
real :: dlat, dlon
real, dimension(nlon,nlat) ::  alat, alon, fmask
real, dimension(nlon,nlat,nvt+1) :: vegmap

! Work variables:

real, parameter :: dlon_read_v=1./120., dlat_read_v=1./120., &
                   lonini_read_v=-180.+dlon_read_v*0.5, lonfin_read_v=180.-dlon_read_v*0.5, &
                   latini_read_v=90.-dlat_read_v*0.5, latfin_read_v=-90.+dlat_read_v*0.5

! lonini_read_v, lonfin_read_v: longitude of the west and east extremes of landcover database
! latini_read_v, latfin_read_v: latitude of the north and south extremes of land database

integer, parameter :: nx_read_v=43200, ny_read_v=(latini_read_v-latfin_read_v)/dlat_read_v

real, dimension(nlon,nlat) ::  lat_geo, lon_geo, work
integer*2, dimension (:,:), allocatable :: veg_read
integer*1 :: veg_read0(nx_read_v)
character (len=50) :: file_input='veget_global_latlon_1km.bin'
integer nitype(nvt)

real :: pi, radinfl, val_missing=-9999., lat_geo_min, lat_geo_max, lon_geo_min, lon_geo_max, &
        lat_start_read_v, lon_start_read_v, &
        zzfac, zcos, zlonwb, zloneb, zlatsb, zlatnb, ziwb, zieb, zjsb, zjnb, &
        znpoint, wei
integer :: nlon_read_v, nlat_read_v, &
           istart, ifinish, jstart, jfinish, jini, jfin, j, jj, jjj, i, ii, iii, flag_period_v=0, &
           ix_part, jlon, jlat, itsqr, iwb, ieb, jsb, jnb, itype, nmax, nmax1, maxtype, maxtype1, npoint, &
           iveg, jveg, iv, nsmooth, if1, iter
real, dimension(4) :: zdist
integer, parameter :: nx_part = 2
integer, dimension(nx_part) :: iini, ifin
integer, dimension(1) :: iind

 pi = abs(acos(-1.))

 lat_geo(:,:) = alat(:,:)
 lon_geo(:,:) = alon(:,:)
 do jlat=1,nlat
 do jlon=1,nlon
   if (abs(abs(lat_geo(jlon,jlat))-90.) < 1.e-10) lat_geo(jlon,jlat) = sign(89.999999, lat_geo(jlon,jlat))
! Dataset longitude: -180..180
   if (alon(jlon,jlat) < -180.) lon_geo(jlon,jlat) = alon(jlon,jlat)+360.
   if (alon(jlon,jlat) >  180.) lon_geo(jlon,jlat) = alon(jlon,jlat)-360.
 enddo
 enddo
 if (flag_glob == 1.and.abs(alon(1,1)-lonini_read_v) > 10.) then
! Reset (below) dataset longitude: 0..360
   lon_geo(:,:) = alon(:,:)
 endif

! For vegetation dataset: "interpolation" to the model grid with resolution,
! coarser then original datset
! radinfl - radius around the grid point within which the veg. dataset pixels are considered:
! about 0.3*dlat for resolution coarser then about 5 km
! about 0.4*dlat for resolution close to vegetation dataset resolution (0.9 km)

 radinfl = dlat*0.3*(1.+0.5*(0.5-0.5*tanh((dlat-0.02)/0.02)))

! Coordinates of the centre of the extremes for reading

 lon_geo_min = minval(lon_geo)
 lon_geo_max = maxval(lon_geo)
 lat_geo_min = minval(lat_geo)
 lat_geo_max = maxval(lat_geo)

 istart = int(lon_geo_min*0.2)*5-10
 ifinish = int(lon_geo_max*0.2)*5+15
 jstart = max( int(lat_geo_min*0.2)*5-10, -90 )
 jfinish = min( int(lat_geo_max*0.2)*5+15, 90 )

! Reading of vegetation original map

 print *
 print *,'Reading the global vegetation dataset from ',file_input

 open(17,file=file_input,status='old',form='unformatted',access='direct',recl=nx_read_v)

! Control of read data domain along latitude axis

 if (flag_glob == 0) then

   if (float(jstart) >= latfin_read_v.and.float(jfinish) <= latini_read_v) then

! Domain extremes located inside of latitude limits of dataset global domain

     jini = int((latini_read_v-float(jstart))/dlat_read_v)+1
     jfin = int((latini_read_v-float(jfinish))/dlat_read_v)+1

 else

! Domain extremes located outside of latitude limits of dataset global domain

     if (float(jstart) < latfin_read_v.and.float(jfinish) <= latini_read_v) then

       jini = ny_read_v
       jfin = int((latini_read_v-float(jfinish))/dlat_read_v)+1
       if (abs(jstart) /= 90) then
         print *,' Caution: vegetation (landuse) dataset does not cover the south part of the model domain'
         print *,' Model grid south latitude extreme',float(jstart)
         print *,' Dataset grid south latitude extreme',latfin_read_v
       endif

     elseif (float(jfinish) > latini_read_v.and.float(jstart) >= latfin_read_v) then

       jini = int((latini_read_v-float(jstart))/dlat_read_v)+1
       jfin = 1
       if (abs(jfinish) /= 90) then
         print *,' Caution: vegetation (landuse) dataset does not cover the north part of the model domain'
         print *,' Model grid north latitude extreme',float(jfinish)
         print *,' Dataset grid north latitude extreme',latini_read_v
       endif

   else

! Domain covers the whole latitude range

       jini = ny_read_v
       jfin = 1
       if (abs(jstart) /= 90.or.abs(jfinish) /= 90) then
         print *,' Caution: vegetation (landuse) dataset does not cover the north and south parts of the model domain'
         print *,' Model grid latitude extremes',float(jstart),float(jfinish)
         print *,' Dataset grid latitude extremes',latfin_read_v,latini_read_v
       endif

     endif

   endif

 else

! Domain covers the whole latitude range

   jini = ny_read_v
   jfin = 1
   print *,' Caution: vegetation (landuse) dataset does not cover the north and south parts of the model domain'
   print *,' Model grid latitude extremes',lat_geo(1,1),lat_geo(1,nlat)
   print *,' Dataset grid latitude extremes',latfin_read_v,latini_read_v

 endif

 lat_start_read_v = latini_read_v-float(jini-1)*dlat_read_v

! Control of read data domain along longitude axis, possible cut and paste operation for periodic grid

 if (flag_glob == 0) then

   if (istart >= int(lonini_read_v).and.ifinish <= int(lonfin_read_v)) then

! Domain extremes located inside of longitude limits of dataset global domain

     iini(1) = int((float(istart)-lonini_read_v)/dlon_read_v)+1
     ifin(1) = int((float(ifinish)-lonini_read_v)/dlon_read_v)+1

     iini(2) = 0
     ifin(2) = -1

     flag_period_v = 0

   else

! Domain extremes located outside of longitude limits of dataset global domain -> cut and paste

     if (istart < int(lonini_read_v).and.ifinish <= int(lonfin_read_v)) then

       iini(1) = int((float(istart)+360.-lonini_read_v)/dlon_read_v)+1
       ifin(1) = nx_read_v

       iini(2) = 1
       ifin(2) = int((float(ifinish)-lonini_read_v)/dlon_read_v)+1

       flag_period_v = 1

     elseif (ifinish > int(lonfin_read_v).and.istart >= int(lonini_read_v)) then

       iini(1) = int((float(istart)-lonini_read_v)/dlon_read_v)+1
       ifin(1) = nx_read_v

       iini(2) = 1
       ifin(2) = int((float(ifinish)-360.-lonini_read_v)/dlon_read_v)+1

       flag_period_v = 1

     else

! Domain covers the whole longitude circle

       iini(1) = 1
       ifin(1) = nx_read_v

       iini(2) = 0
       ifin(2) = -1

       flag_period_v = 2

     endif

   endif

 else ! global output grid

! Domain covers the whole longitude circle

   if (abs(alon(1,1)-lonini_read_v) < 10.) then

! Input (dataset) and output (model) global grid have the same longitude origin

     iini(1) = 1
     ifin(1) = nx_read_v

     iini(2) = 0
     ifin(2) = -1

     flag_period_v = 2

   else

! Input (dataset) and output (model) global grid have the longitude origin shited by 180 degees

     iini(1) = int((alon(1,1)-lonini_read_v)/dlon_read_v)+1
     ifin(1) = nx_read_v

     iini(2) = 1
     ifin(2) = iini(1)-1

     flag_period_v = 3

   endif

 endif

 lon_start_read_v = lonini_read_v+float(iini(1)-1)*dlon_read_v

 nlat_read_v = jini-jfin+1
 nlon_read_v = ifin(1)-iini(1)+1
 if (flag_period_v == 1 ) then
   nlon_read_v = nlon_read_v+ifin(2)-iini(2)+1
 endif
 if (flag_period_v == 2.or.flag_period_v == 3 ) then
   nlon_read_v = nx_read_v+2
   if (flag_period_v == 2) lon_start_read_v = lonini_read_v-dlon_read_v
   if (flag_period_v == 3) lon_start_read_v = lon_start_read_v-dlon_read_v
 endif

 allocate(veg_read(nlon_read_v,nlat_read_v))
 veg_read(:,:) = int(val_missing)

 j=0
 do jjj=jini, jfin, -1
   j = j+1

   read (17, rec=jjj) veg_read0

   i = 0
   if (flag_period_v == 2.or.flag_period_v == 3 ) i = 1
   do ix_part=1,nx_part
     do iii=iini(ix_part), ifin(ix_part)
       i = i+1
       veg_read(i,j) = veg_read0(iii)
     enddo
   enddo ! ix_part
   if (flag_period_v == 2.or.flag_period_v == 3 ) then
     veg_read(1,j) = veg_read(nlon_read_v-1,j)
     veg_read(nlon_read_v,j) = veg_read(2,j)
   endif

 enddo ! jjj

 close (17)

 print *,'Vegetation data processing'

!------------------------------------------------------------
! Beginning vegetation type definition

! "Interpolation" of vegetation types from the original vegetation map grid
! IVEGBAS: 2-D with vegetation types (landcover types) code from 0 to 13:
! (Note: code number 0 in the original dataset has been substituted by 14 (=nsv) )
! See also subroutine deflai for a more precise definition of vegetation types below

! 1  evergreen needleleaf forest
! 2  evergreen broadleaf forest
! 3  deciduous needleleaf forest
! 4  deciduous broadleaf forest
! 5  mixed cover
! 6  woodland
! 7  wooded grassland
! 8  closed shrubland
! 9  open shrubland
! 10 grassland
! 11 cropland
! 12 bare ground ---> 13
! 13 urban and built-up ---> 12
!  0  water body (sea, lakes, rivers) ---> 14

! Code -9999: missing value code

 do j=1,nlat_read_v
 do i=1,nlon_read_v
   if (veg_read(i,j) == val_missing) then
     print*, 'Missing value code found in veget. dataset'
     print*, 'Defined as bare ground'
     veg_read(i,j) = nvt-1
   endif
   if (veg_read(i,j) == 0) veg_read(i,j) = nvt
 enddo
 enddo

! 'Interpolation' from the original 1/120 deg grid of veg_read to the model grid (VEGMAP):
! considering the higher resolution of IVEGBAS, in VEGMAP the most frequent veget. type
! (first value) in the model cell and fractions of all vegetation types are stored.

! FMASK (0-1, where 1 is all water) is defined here using the index water body
! and the fraction of water body compared to all types found in each model cell

 lon_geo(:,:) = alon(:,:)
 do jlat=1,nlat
 do jlon=1,nlon
! Dataset longitude: -180..180
   if (alon(jlon,jlat) < -180.) lon_geo(jlon,jlat) = alon(jlon,jlat)+360.
   if (alon(jlon,jlat) >  180.) lon_geo(jlon,jlat) = alon(jlon,jlat)-360.
 enddo
 enddo
 if (flag_glob == 1.and.abs(alon(1,1)-lonini_read_v) > 10.) then
! Reset (below) dataset longitude: 0..360
   lon_geo(:,:) = alon(:,:)
 endif

 lon_geo(:,:) = alon(:,:)
 do jlat=1,nlat
 do jlon=1,nlon
! Dataset longitude: -180..180
   if (alon(jlon,jlat) < -180.) lon_geo(jlon,jlat) = alon(jlon,jlat)+360.
   if (alon(jlon,jlat) >  180.) lon_geo(jlon,jlat) = alon(jlon,jlat)-360.
 enddo
 enddo
 if (flag_period_v == 1) then
   do jlat=1,nlat
   do jlon=1,nlon
     if (lon_geo(jlon,jlat) < 0.) lon_geo(jlon,jlat) = lon_geo(jlon,jlat)+360.
   enddo
   enddo
 endif
 if (flag_period_v == 3) then
! Reseted dataset longitude: 0..360
   lon_geo(:,:) = alon(:,:)
 endif

 zzfac = pi/180.

! Definition of 'interpolation' mode flag

   if (dlat > dlat_read_v*1.2) then

   print*, "Vegetation and land-sea fraction defined counting pixels of types in a surrounding area"

   itsqr = 1
   do j=1,nlat
   do i=1,nlon

     zcos = cos( (max( min(lat_geo(i,j), 89.), -89.))*zzfac )
     zlonwb = lon_geo(i,j)-radinfl/zcos
     zloneb = lon_geo(i,j)+radinfl/zcos
     zlatsb = lat_geo(i,j)-radinfl
     zlatnb = lat_geo(i,j)+radinfl

     ziwb = (zlonwb-lon_start_read_v)/dlon_read_v
!     iwb = max(min(nint(ziwb)+1, nlon_read_v), 1)     ! old
     iwb = max(min(nint(ziwb+1.5), nlon_read_v), 1)    ! Nov. 2014

     zieb = (zloneb-lon_start_read_v)/dlon_read_v
!     ieb = max(min(nint(zieb)+1, nlon_read_v), 1)      ! old
     ieb = max(min(nint(zieb+1.5), nlon_read_v), 1)     ! Nov. 2014

     zjsb = (zlatsb-lat_start_read_v)/dlat_read_v
     jsb = max(min(nint(zjsb)+1, nlat_read_v), 1)

     zjnb = (zlatnb-lat_start_read_v)/dlat_read_v
     jnb = max(min(nint(zjnb)+1, nlat_read_v), 1)

     do itype=1,nvt
       nitype(itype) = 0
     end do

     do jj=jsb,jnb
     do ii=iwb,ieb
       itype = veg_read(ii,jj)
       nitype(itype) = nitype(itype)+1
     enddo
     enddo

     nmax = 0
     if (itsqr == 1) then
       do itype=1,nvt
         if (nitype(itype) > nmax) then
           nmax = nitype(itype)
           maxtype = itype
          endif
       enddo
       itsqr = itsqr+1
     else
       do itype=nvt,1,-1
         if (nitype(itype) > nmax) then
           nmax = nitype(itype)
           maxtype = itype
         endif
       enddo
       itsqr = itsqr-1
     endif

     nmax1 = 0
     do itype=1,nvt
       if (itype /= maxtype) then
         if (nitype(itype) > nmax1) then
           nmax1 = nitype(itype)
           maxtype1 = itype
         endif
       endif
     enddo

     if (maxtype == nvt) then
       if (nmax1 > 0) maxtype = maxtype1
     endif

     vegmap(i,j,1) = float(maxtype)

! Fraction of vetegation type exept for water body

     npoint = 0
     do itype=1,nvt-1
       npoint = npoint+nitype(itype)
     enddo

     if (npoint > 0) then ! not water body

       znpoint = 1./float(npoint)
       do itype=1,nvt-1
         vegmap(i,j,itype+1) = float(nitype(itype))*znpoint
       enddo
       vegmap(i,j,nvt+1) = 0.

     else ! water body

       vegmap(i,j,:) = 0.
       vegmap(i,j,1) = float(nvt)
       vegmap(i,j,nvt+1) = 1.

     endif

! Definition of fmask

     npoint = 0
     do itype=1,nvt
       npoint = npoint+nitype(itype)
     enddo

      if (npoint > 0) then
        fmask(i,j) = float(nitype(nvt))/float(npoint)
      else
        fmask(i,j) = 1.
     endif

   enddo
   enddo

 else     ! dlat <= dlat_read_v*1.2

   print*, "Vegetation and land-sea fraction computed by direct grid interpolation"

   do j=1,nlat
   do i=1,nlon

! Search of the nearest database pixel (by pixel centre coordinate):

!     iveg = nint((lon_geo(i,j)-lon_start_read_v)/dlon_read_v)+1      ! old
     iveg = nint((lon_geo(i,j)-lon_start_read_v)/dlon_read_v+0.5)+1   ! Nov. 2014
     jveg = nint((lat_geo(i,j)-lat_start_read_v)/dlat_read_v)+1

     zdist(1) = (lon_start_read_v+dlon_read_v*float(iveg-1)-lon_geo(i,j))**2 +  &
                (lat_start_read_v+dlat_read_v*float(jveg-1)-lat_geo(i,j))**2
     zdist(2) = (lon_start_read_v+dlon_read_v*float(iveg)  -lon_geo(i,j))**2 +  &
                (lat_start_read_v+dlat_read_v*float(jveg-1)-lat_geo(i,j))**2
     zdist(3) = (lon_start_read_v+dlon_read_v*float(iveg)  -lon_geo(i,j))**2 +  &
                (lat_start_read_v+dlat_read_v*float(jveg)  -lat_geo(i,j))**2
     zdist(4) = (lon_start_read_v+dlon_read_v*float(iveg-1)-lon_geo(i,j))**2 +  &
                (lat_start_read_v+dlat_read_v*float(jveg)  -lat_geo(i,j))**2

     iind = minloc(zdist)

     if (iind(1) == 2.or.iind(1) == 3) iveg = iveg+1
     if (iind(1) == 3.or.iind(1) == 4) jveg = jveg+1

     iveg = max(min(iveg, nlon_read_v), 1)
     jveg = max(min(jveg, nlat_read_v), 1)

     iv = veg_read(iveg,jveg)
     vegmap(i,j,1) = float(iv)

     do itype=1,nvt
       vegmap(i,j,itype+1) = 0.
     enddo
     vegmap(i,j,iv+1) = 1.

! Definition of fmask

     if (iv == nvt) then
       fmask(i,j) = 1.
     else
       fmask(i,j) = 0.
     endif

   enddo
   enddo

   wei = 0.5
   nsmooth = 1
   call smooth_soil(fmask,work,nlon,nlat,wei,nsmooth)

 endif    ! dlat > dlat_read_v*1.2

 fmask(:,:) = fmask(:,:)**1.7   ! ad hoc correction to slightly enlarge land areas

! Reduction of isolated (very small) lakes or thin rivers

 do iter=1,2

   work(1:nlon,1:nlat) = fmask(1:nlon,1:nlat)

   do i=2,nlon-1
   do j=2,nlat-1
     if1 = 0
     do ii=-1,1
     do jj=-1,1
       if (work(i+ii,j+jj) > 0.5) if1 = if1+1 ! >, not >= !
     enddo
     enddo
     if(fmask(i,j) >= 0.5.and.if1 <= 2) fmask(i,j) = .45 ! >=, not > !
   enddo
   enddo

 enddo

return
end subroutine dataset_vegetation
!=======================================================================
subroutine dataset_ecmwf_lsm(nlon, nlat, alon, alat, flag_glob, lsm_ecmwf)

! Procedure defined land-sea fraction over the model grid using IFS-ECMWF data about land-sea mask 
! read from input file ecmwf_glob_0_25_lsm.bin (resolution 0.25 degree)

! Input variables:

!      nlon: grid point number along axis of rotated longitude;
!      nlat: grid point number along axis of rotated latidude;
!      alon: geographical longitude of the grid points (degree, -180...180);
!      alat: geographical latitude of the grid points (degree, -90...90) 
! flag_glob: if 0, then not global input grid, if 1, then input grid is global not rotated

! Output variables:

!      lsm_ecmwf: land-sea franction: 1-sea, 0-land (proportion) according to IFS-ECMWF data woth 0.25 degree resolution

! Land-sea mask data of IFS-ECMWF
! defined on a lat-lon regular grid, with resolution of 0.25 deg.
! lon -180 --> 180, lat -90 --> 90,
! lsm=10: land, lsm=0: sea

implicit none

integer :: nlon, nlat, flag_glob
real, dimension(nlon,nlat) ::  alat, alon, lat_geo, lon_geo
real, dimension (:,:), allocatable :: lsm_read
real, dimension (:), allocatable :: lon_read, lat_read
integer*1, dimension(:,:), allocatable :: lsm_read0
character (len=50) :: file_input='ecmwf_glob_0_25_lsm.bin'

real, dimension(nlon,nlat) :: lsm_ecmwf

! Work variables:

integer :: nx_read_lsm, ny_read_lsm
real :: dlon_read_lsm, dlat_read_lsm, &
        lonini_read_lsm, latini_read_lsm, lonfin_read_lsm, latfin_read_lsm

! dlon_read_lsm, dlat_read_lsm: resolutions IFS-ECMWF LSM databases
! lonini_read_lsm, lonfin_read_lsm: longitude of the west and east extremes of IFS-ECMWF LSM database
! latini_read_lsm, latfin_read_lsm: latitude of the north and south extremes of IFS-ECMWF LSM databases

real :: val_missing=-9999., lat_geo_min, lat_geo_max, lon_geo_min, lon_geo_max, &
        lat_start_read_lsm, lon_start_read_lsm
integer :: nlon_read_lsm, nlat_read_lsm, &
           istart, ifinish, jstart, jfinish, jini, jfin, j, jjj, i, iii, &
           ix_part, flag_period_lsm, jlon, jlat
integer, parameter :: nx_part = 2
integer, dimension(nx_part) :: iini, ifin

 lat_geo(:,:) = alat(:,:)
 lon_geo(:,:) = alon(:,:)
 do jlat=1,nlat
 do jlon=1,nlon
   if (abs(abs(lat_geo(jlon,jlat))-90.) < 1.e-10) lat_geo(jlon,jlat) = sign(89.999999, lat_geo(jlon,jlat))
   if (alon(jlon,jlat) < -180.) lon_geo(jlon,jlat) = alon(jlon,jlat)+360.
   if (alon(jlon,jlat) >  180.) lon_geo(jlon,jlat) = alon(jlon,jlat)-360.
 enddo
 enddo

! Coordinates of the centre of the extremes for reading

 lon_geo_min = minval(lon_geo)
 lon_geo_max = maxval(lon_geo)
 lat_geo_min = minval(lat_geo)
 lat_geo_max = maxval(lat_geo)

 istart = int(lon_geo_min*0.2)*5-10
 ifinish = int(lon_geo_max*0.2)*5+15
 jstart = max( int(lat_geo_min*0.2)*5-10, -90 )
 jfinish = min( int(lat_geo_max*0.2)*5+15, 90 )

 print *
 print *,'Reading IFS-ECMWF land-sea mask database from ',file_input

 open (18, file=file_input, form='unformatted')

 read (18) ny_read_lsm, nx_read_lsm
 read (18) latini_read_lsm, lonini_read_lsm, dlat_read_lsm, dlon_read_lsm
 latfin_read_lsm=latini_read_lsm+float(ny_read_lsm-1)*dlat_read_lsm
 lonfin_read_lsm=lonini_read_lsm+float(nx_read_lsm-1)*dlon_read_lsm

 allocate (lsm_read0(nx_read_lsm,ny_read_lsm))

 read (18) lsm_read0

 close (18)

 if (flag_glob == 0) then

! Control of read data domain along latitude axis

   if (float(jstart) >= latini_read_lsm.and.float(jfinish) <= latfin_read_lsm) then

! Domain extremes locate inside of latitude limits of dataset global domain

     jini = int((float(jstart)-latini_read_lsm)/dlat_read_lsm)+1
     jfin = int((float(jfinish)-latini_read_lsm)/dlat_read_lsm)+1

   else

! Domain extremes located outside of latitude limits of dataset global domain

     if (float(jstart) < latini_read_lsm.and.float(jfinish) <= latfin_read_lsm) then

       jini = 1
       jfin = int((float(jfinish)-latini_read_lsm)/dlat_read_lsm)+1
       if (abs(jstart) /= 90) then
         print *,' Caution: IFS-ECMWF land-sea mask dataset does not cover the south part of the model domain'
         print *,' Model grid south latitude extreme',float(jstart)
         print *,' Dataset grid south latitude extreme',latini_read_lsm
       endif

     elseif (float(jfinish) > latfin_read_lsm.and.float(jstart) >= latfin_read_lsm) then

       jini = int((float(jstart)-latini_read_lsm)/dlat_read_lsm)+1
       jfin = ny_read_lsm
       if (abs(jfinish) /= 90) then
         print *,' Caution: IFS-ECMWF land-sea mask dataset does not cover the north part of the model domain'
         print *,' Model grid north latitude extreme',float(jfinish)
         print *,' Dataset grid north latitude extrem',latfin_read_lsm
       endif

     else

! Domain covers the whole latitude range

       jini = 1
       jfin = ny_read_lsm
       if (abs(jstart) /= 90.or.abs(jfinish) /= 90) then
         print *,' Caution: IFS-ECMWF land-sea mask dataset does not cover the north and south parts of the model domain'
         print *,' Model grid latitude extremes',float(jstart),float(jfinish)
         print *,' Dataset grid latitude extremes',latini_read_lsm,latfin_read_lsm
       endif

     endif

   endif

 else

! Domain covers the whole latitude range

   jini = 1
   jfin = ny_read_lsm

 endif

 lat_start_read_lsm = latini_read_lsm+float(jini-1)*dlat_read_lsm

! Control of read data domain along longitude axis, possible cut and paste operation for periodic grid

 if (flag_glob == 0) then

   if (istart >= int(lonini_read_lsm).and.ifinish <= int(lonfin_read_lsm)) then

! Domain extremes locate inside of longitude limits of dataset global domain

     iini(1) = int((float(istart)-lonini_read_lsm)/dlon_read_lsm)+1
     ifin(1) = int((float(ifinish)-lonini_read_lsm)/dlon_read_lsm)+1

     iini(2) = 0
     ifin(2) = -1

     flag_period_lsm = 0

   else

! Domain extremes locate outside of longitude limits of dataset global domain -> cut and paste

     if (istart < int(lonini_read_lsm).and.ifinish <= int(lonfin_read_lsm)) then

       iini(1) = int((float(istart)+360.-lonini_read_lsm)/dlon_read_lsm)+1
       ifin(1) = nx_read_lsm

       iini(2) = 1
       ifin(2) = int((float(ifinish)-lonini_read_lsm)/dlon_read_lsm)+1

       flag_period_lsm = 1

     elseif (ifinish > int(lonfin_read_lsm).and.istart >= int(lonini_read_lsm)) then

       iini(1) = int((float(istart)-lonini_read_lsm)/dlon_read_lsm)+1
       ifin(1) = nx_read_lsm

       iini(2) = 1
       ifin(2) = int((float(ifinish)-360.-lonini_read_lsm)/dlon_read_lsm)+1

       flag_period_lsm = 1

     else

! Domain covers the whole longitude circle

       iini(1) = 1
       ifin(1) = nx_read_lsm

       iini(2) = 0
       ifin(2) = -1

       flag_period_lsm = 2

     endif

   endif

 else ! global output grid

! Domain covers the whole longitude circle

   if (abs(alon(1,1)-lonini_read_lsm) < 10.) then

! Input (dataset) and output (model) global grid have the same longitude origin

     iini(1) = 1
     ifin(1) = nx_read_lsm

     iini(2) = 0
     ifin(2) = -1

     flag_period_lsm = 2

   else

! Input (dataset) and output (model) global grid have the longitude origin shited by 180 degees

     iini(1) = int((alon(1,1)-lonini_read_lsm)/dlon_read_lsm)+1
     ifin(1) = nx_read_lsm

     iini(2) = 1
     ifin(2) = iini(1)-1

     flag_period_lsm = 3

   endif

 endif

 lon_start_read_lsm = lonini_read_lsm+float(iini(1)-1)*dlon_read_lsm

 nlat_read_lsm = jfin-jini+1
 nlon_read_lsm = ifin(1)-iini(1)+1
 if (flag_period_lsm == 1 ) then
   nlon_read_lsm = nlon_read_lsm+ifin(2)-iini(2)+1
 endif
 if (flag_period_lsm == 2.or.flag_period_lsm == 3 ) then
   nlon_read_lsm = nx_read_lsm+2
   if (flag_period_lsm == 2) lon_start_read_lsm = lonini_read_lsm-dlon_read_lsm
   if (flag_period_lsm == 3) lon_start_read_lsm = lon_start_read_lsm-dlon_read_lsm
 endif

 allocate(lsm_read(nlon_read_lsm,nlat_read_lsm))
 allocate(lon_read(nlon_read_lsm))
 allocate(lat_read(nlat_read_lsm))
 lsm_read(:,:) = val_missing
 lon_read(:) = val_missing
 lat_read(:) = val_missing

 j = 0
 do jjj=jini, jfin
   j = j+1

   i = 0
   if (flag_period_lsm == 2.or.flag_period_lsm == 3) i=1
   do ix_part=1,nx_part
     do iii=iini(ix_part), ifin(ix_part)
       i = i+1
       lsm_read(i,j) = float(1-lsm_read0(iii,jjj)/10)
       lon_read(i) = lon_start_read_lsm+float(i-1)*dlon_read_lsm
     enddo
   enddo ! ix_part
   if (flag_period_lsm == 2 ) then
     lsm_read(1,j) = lsm_read(nlon_read_lsm-1,j)
     lsm_read(nlon_read_lsm,j) = lsm_read(2,j)
   endif

   lat_read(j) = lat_start_read_lsm+float(j-1)*dlat_read_lsm
 enddo ! jjj
 if (flag_period_lsm == 2.or.flag_period_lsm == 3) then
   lon_read(1) = lon_start_read_lsm
   lon_read(nlon_read_lsm) = lon_start_read_lsm+float(nlon_read_lsm-1)*dlon_read_lsm
 endif

! Interpolation of the read LSM field to the model grid using model rotated coordinate system

 lon_geo(:,:) = alon(:,:)
 do jlat=1,nlat
 do jlon=1,nlon
   if (abs(abs(lat_geo(jlon,jlat))-90.) < 1.e-10) lon_geo(jlon,jlat) = 0.
   if (alon(jlon,jlat) < -180.) lon_geo(jlon,jlat) = alon(jlon,jlat)+360.
   if (alon(jlon,jlat) >  180.) lon_geo(jlon,jlat) = alon(jlon,jlat)-360.
 enddo
 enddo
 if (flag_glob == 1.and.flag_period_lsm == 3) then
! Reseted dataset longitude: 0..360
   lon_geo(:,:) = alon(:,:)
 endif
 if (flag_period_lsm == 1) then
   do jlat=1,nlat
   do jlon=1,nlon
     if (lon_geo(jlon,jlat) < 0.) lon_geo(jlon,jlat) = lon_geo(jlon,jlat)+360.
   enddo
   enddo
 endif

 call interp_spline_2d(lsm_read, nlon_read_lsm, nlat_read_lsm, lon_read, lat_read, &
                       lon_geo, lat_geo, nlon*nlat, lsm_ecmwf, 1.)

return
end subroutine dataset_ecmwf_lsm
!=======================================================================
SUBROUTINE DATASET_VEGETATION_CYCLOPES(NLON,NLAT,X0,Y0,DLON,DLAT,LONINI,LATINI,ALON,ALAT,FLAG_GLOB,&
 NVT,VAL_MISSING,VEGET,LSF)

! Procedure defines vegetation (landuse) types (classification GLC2000) over the model grid,
! and land-sea fraction over the model grid using landuse dataset
! between -60 and 80 degree latitude only !!!

! Input variables:

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
!      NVT: total number vegetaion types;
!      VAL_MISSING: values used for not difened point.

! Output variables:

!      VEGET: vegetation (landuse) types (classification GLC2000, 22 types) defined on the grid points, 3D array,
!               1-st and 2-nd array index are grid point index, 3-rd index is vegetation types number plus 1,
!               if 3-rd index has value 1, then the variable means dominate vegetation type,
!               values 2-23 of 3-rd index means fraction (of area, proportion) of the each of 22 vegetation types
!      LSF: land-sea franction: 1-sea, 0-land (proportion).

! Vegetation types are defined starting from the global dataset in file:
!             vegtype_GLC2000_global_latlon_1km.bin

! The global vegetation cover dataset (with resolution 1/112 degree) is
! vegetation types using in CYCLOPES dataset,
! this dataset is following GLC2000 classification (http://bioval.jrc.ec.europa.eu/products/glc2000/legend.php) 

! GLC2000 classification:

!  1 - Tree cover, broadleaved evergreen, closed to open (>15%)
!  2 - Tree cover, broadleaved deciduous, closed (>40%)
!  3 - Tree cover, broadleaved deciduous, open (15-40%)
!  4 - Tree cover, needleleaved evergreen,  closed to open (>15%)
!  5 - Tree cover, needleleaved decidous, closed to open (>15%)
!  6 - Tree cover, mixed leaftype,  closed to open (>15%)
!  7 - Tree cover, closed to open (>15%), regularly flooded, fresh or brackish water: Swamp Forests
!  8 - Tree cover, closed to open (>15%),  regularly flooded, saline water: Mangrove Forests
!  9 - Mosaic of tree cover and other natural vegetation (Crop component possible)
! 10 - Tree Cover, burnt (mainly boreal forests)
! 11 - Shrubcover, closed to open (>15%) , evergreen(broadleaved or needleleaved)
! 12 - Shrubcover, closed to open (>15%), deciduous (broadleaved)
! 13 - Herbaceous cover, closed to open (>15%)
! 14 - Sparse Herbaceous or sparse Shrub cover
! 15 - Regularly flooded ( > 2 month) Shrub or Herbaceous cover, closed to open 
! 16 - Cropland (upland crops or inundated/ flooded crops as e.g. rice)
! 17 - Mosaic of Cropland / Tree cover/ Other Natural Vegetation
! 18 - Mosaic of Cropland / Shrub or Herbaceous cover
! 19 - Bare Areas
! 20 - Urban Areas
! 21 - Snow or Ice (natural or artificial)
! 22 - Water Bodies (natural or artificial)

IMPLICIT NONE

INTEGER :: NLON, NLAT, FLAG_GLOB, NVT
REAL :: X0, Y0, DLON, DLAT, LONINI, LATINI, VAL_MISSING
REAL, DIMENSION(NLON,NLAT,NVT+1) :: VEGET
REAL, DIMENSION(NLON,NLAT) :: LSF, ALON, ALAT, LON_GEO, LAT_GEO
REAL, PARAMETER :: DLON_READ=1./112., DLAT_READ=1./112.,&
                   LONINI_READ=0.+DLON_READ*0.5, LONFIN_READ=360.-DLON_READ*0.5,&
                   LATINI_READ=-60.+DLAT_READ*0.5, LATFIN_READ=80.-DLAT_READ*0.5
INTEGER, PARAMETER :: NX_READ=40320, NY_READ=(LATFIN_READ-LATINI_READ)/DLAT_READ
INTEGER*1, DIMENSION(NX_READ) :: VALUE_READ
INTEGER, DIMENSION(:,:), ALLOCATABLE :: VEGET_READ
REAL, DIMENSION(:,:), ALLOCATABLE :: LON_READ_GEO, LAT_READ_GEO, LON_READ_MOD, LAT_READ_MOD 

CHARACTER(LEN=50) :: FILEINPUT="vegtype_GLC2000_global_latlon_1km.bin"

INTEGER, DIMENSION(NLON,NLAT) :: NPOINT
INTEGER, DIMENSION(NLON,NLAT,NVT) :: NPOINT_VEGTYPE
INTEGER, DIMENSION(NVT) :: NTYPE
INTEGER :: IUNIT=10, NLON_READ, NLAT_READ, IVAL_MISSING, ISTART, IFINISH, JSTART, JFINISH,&
 I, J, IR, JR, IX_PART, III, JJJ, JINI, JFIN, FLAG_PERIOD=0, &
 ITSQR, IMOD, JMOD, ITYPE, IWB, IEB, JSB, JNB, NMAX, MAXTYPE, NMAX2, MAXTYPE2, NP
REAL :: PI, LON_GEO_MIN, LON_GEO_MAX, LAT_GEO_MIN, LAT_GEO_MAX, LON_START_READ, LAT_START_READ, &
 RADINFL, ZFAC, ZCOS, ZLONWB, ZLONEB, ZLATSB, ZLATNB, ZIWB, ZIEB, ZJSB, ZJNB, ZNP
INTEGER, PARAMETER :: NX_PART=2
INTEGER, DIMENSION(NX_PART) :: IINI, IFIN
INTEGER, DIMENSION(:,:), ALLOCATABLE :: IWORK
INTEGER, DIMENSION(1) :: IARRAY

 PI=ABS(ACOS(-1.))
 ZFAC=PI/180.
 IVAL_MISSING=INT(VAL_MISSING)

 PRINT *
 PRINT *,"Reading the global vegetation dataset from ",FILEINPUT

 LON_GEO_MIN=MINVAL(ALON)
 LON_GEO_MAX=MAXVAL(ALON)

 LAT_GEO(:,:)=ALAT(:,:)
 LON_GEO(:,:)=ALON(:,:)

 IF (LON_GEO_MIN >= 0..AND.LON_GEO_MAX > 0.) THEN
   LON_GEO(:,:)=ALON(:,:)
 ELSE
   DO J=1,NLAT
   DO I=1,NLON
     IF (ALON(I,J) < -180.) LON_GEO(I,J) = ALON(I,J)+360.
     IF (ALON(I,J) >  180.) LON_GEO(I,J) = ALON(I,J)-360.
   ENDDO
   ENDDO
 ENDIF

 DO J=1,NLAT
 DO I=1,NLON
   IF (ABS(ABS(LAT_GEO(I,J))-90.) < 1.e-10) LAT_GEO(I,J) = SIGN(89.999999, LAT_GEO(I,J))
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
         print *,' Caution: GLC2000 landuse dataset does not cover the south part of the model domain'
         print *,' Model grid south latitude extreme',FLOAT(JSTART)
         print *,' Dataset grid south latitude extreme',LATINI_READ
       ENDIF

     ELSE IF (FLOAT(JFINISH) > LATFIN_READ.AND.FLOAT(JSTART) >= LATFIN_READ) THEN

       JINI = INT((FLOAT(JSTART)-LATINI_READ)/DLAT_READ)+1
       JFIN = NY_READ
       IF (FLOAT(JFINISH) > LATFIN_READ) THEN
         print *,' Caution: GLC2000 landuse dataset does not cover the north part of the model domain'
         print *,' Model grid north latitude extreme',FLOAT(JFINISH)
         print *,' Dataset grid north latitude extrem',LATFIN_READ
       ENDIF

     ELSE

! Domain covers the whole latitude range

       JINI = 1
       JFIN = NY_READ
       IF (FLOAT(JSTART) < LATINI_READ.OR.FLOAT(JFINISH) > LATFIN_READ) THEN
         print *,' Caution: GLC2000 landuse dataset does not cover the north and the south parts of the model domain'
         print *,' Model grid south latitude extremes',FLOAT(JSTART),FLOAT(JFINISH)
         print *,' Dataset grid south latitude extremes',LATINI_READ,LATFIN_READ
       ENDIF

     ENDIF

   ENDIF

 ELSE

! Domain covers the whole latitude range

   JINI = 1
   JFIN = NY_READ
   print *,' Caution: GLC2000 landuse dataset does not cover the north and the south parts of the model domain'
   print *,' Model grid south latitude extremes',LAT_GEO(1,1),LAT_GEO(1,NLAT)
   print *,' Dataset grid south latitude extremes',LATINI_READ,LATFIN_READ

 ENDIF

 LAT_START_READ = LATINI_READ+FLOAT(JINI-1)*DLAT_READ

 IF (FLAG_GLOB == 0) THEN

! Control of read data domain along longitude axis, possible cut and paste operation for periodic grid

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
 IF (FLAG_PERIOD == 2.OR.FLAG_PERIOD == 3 ) THEN
   NLON_READ = NX_READ+2
   IF (FLAG_PERIOD == 2) LON_START_READ = LONINI_READ-DLON_READ
   IF (FLAG_PERIOD == 3) LON_START_READ = LON_START_READ-DLON_READ
 ENDIF

 IF (LON_START_READ > 180.) LON_START_READ=LON_START_READ-360.

 ALLOCATE(VEGET_READ(NLON_READ,NLAT_READ))
 ALLOCATE(LON_READ_GEO(NLON_READ,NLAT_READ))
 ALLOCATE(LAT_READ_GEO(NLON_READ,NLAT_READ))
 ALLOCATE(LON_READ_MOD(NLON_READ,NLAT_READ))
 ALLOCATE(LAT_READ_MOD(NLON_READ,NLAT_READ))

 VEGET(:,:,:)=VAL_MISSING
 VEGET_READ(:,:)=IVAL_MISSING
 LON_READ_GEO(:,:)=VAL_MISSING
 LAT_READ_GEO(:,:)=VAL_MISSING
 LON_READ_MOD(:,:)=VAL_MISSING
 LAT_READ_MOD(:,:)=VAL_MISSING
 LSF(:,:)=VAL_MISSING

! Reading of dataset

 OPEN(IUNIT,FILE=FILEINPUT,STATUS='OLD',FORM='UNFORMATTED',ACCESS='DIRECT',RECL=NX_READ)

 J=0
 DO JJJ=JINI,JFIN
   J=J+1

   READ (IUNIT, REC=JJJ) VALUE_READ

   I=0
   IF (FLAG_PERIOD == 2.OR.FLAG_PERIOD == 3) I=1
   DO IX_PART=1,NX_PART
     DO III=IINI(IX_PART),IFIN(IX_PART)
       I=I+1
       VEGET_READ(I,J)=INT(FLOAT(VALUE_READ(III)))
       LON_READ_GEO(I,J)=LON_START_READ+FLOAT(I-1)*DLON_READ
       IF (VEGET_READ(I,J) < 0) VEGET_READ(I,J)=VEGET_READ(I,J)+256
       IF (VEGET_READ(I,J) >= 255.OR.VEGET_READ(I,J) < 0) VEGET_READ(I,J)=IVAL_MISSING
     ENDDO
   ENDDO ! IX_PART
   IF (FLAG_PERIOD == 2.OR.FLAG_PERIOD == 3) THEN
     VEGET_READ(1,J)=VEGET_READ(NLON_READ-1,J)
     VEGET_READ(NLON_READ,J)=VEGET_READ(2,J)
     LON_READ_GEO(1,J)=LON_START_READ
     LON_READ_GEO(NLON_READ,J)=LON_START_READ+FLOAT(NLON_READ-1)*DLON_READ
   ENDIF

   LAT_READ_GEO(:,J)=LAT_START_READ+FLOAT(J-1)*DLAT_READ

 ENDDO ! JJJ

 CLOSE(IUNIT)

! Space interpolation = definition of fraction of eache vegetation types around
! model grid point and definition domain type

 PRINT *,"Vegetation data processing"

 CALL ANTI_ROT_GRID(X0,Y0,LON_READ_GEO,LAT_READ_GEO,LON_READ_MOD,LAT_READ_MOD,NLON_READ,NLAT_READ)

 DO J=1,NLAT_READ
 DO I=1,NLON_READ
  IF (LON_READ_GEO(I,J)==VAL_MISSING) LON_READ_MOD(I,J)=VAL_MISSING
  IF (LAT_READ_GEO(I,J)==VAL_MISSING) LAT_READ_MOD(I,J)=VAL_MISSING
 ENDDO
 ENDDO

! Caution: change of order types in classification !!!

 ALLOCATE(IWORK(NLON_READ,NLAT_READ))
 IWORK(1:NLON_READ,1:NLAT_READ)=VEGET_READ(1:NLON_READ,1:NLAT_READ)

 DO J=1,NLAT_READ
 DO I=1,NLON_READ
! Significance of indexes 50, 23, 24 and negative value is found by empirical mode (there are not in dataset legend) 
   IF (IWORK(I,J)==50.OR.IWORK(I,J)==23.OR.IWORK(I,J)==20.OR.IWORK(I,J)<=0) VEGET_READ(I,J)=NVT ! Water Body
   IF (IWORK(I,J)==24) VEGET_READ(I,J)=19 ! Bare Areas
   IF (IWORK(I,J)==21) VEGET_READ(I,J)=NVT-1 ! Snow or Ice
   IF (IWORK(I,J)==22) VEGET_READ(I,J)=20 ! Urban Areas
 ENDDO 
 ENDDO 

 DEALLOCATE(IWORK)

! Definition of 'interpolation' mode

 IF (DLAT > DLAT_READ*1.2) THEN ! Resoltuion of dataset is higher then model grid resolution

   PRINT *, "Vegetation and land-sea fraction defined counting pixels of types in a surrounding area"

   RADINFL = DLAT*0.3*(1.+0.5*(0.5-0.5*TANH((DLAT-0.02)/0.02)))

   ITSQR=1
   DO J=1,NLAT
   DO I=1,NLON

     ZCOS=COS( (MAX( MIN(LAT_GEO(I,J), 89.), -89.))*ZFAC)     
     ZLONWB=LON_GEO(I,J)-RADINFL/ZCOS
     ZLONEB=LON_GEO(I,J)+RADINFL/ZCOS
     ZLATSB=LAT_GEO(I,J)-RADINFL
     ZLATNB=LAT_GEO(I,J)+RADINFL

     ZIWB=(ZLONWB-LON_START_READ)/DLON_READ
     IWB=MAX(MIN(NINT(ZIWB)+1, NLON_READ), 1)     ! old
!     IWB=MAX(MIN(NINT(ZIWB+1.5), NLON_READ), 1)   ! Nov. 2014

     ZIEB=(ZLONEB-LON_START_READ)/DLON_READ
     IEB=MAX(MIN(NINT(ZIEB)+1, NLON_READ), 1)     ! old
!     IEB=MAX(MIN(NINT(ZIEB+1.5), NLON_READ), 1)   ! Nov. 2014

     ZJSB=(ZLATSB-LAT_START_READ)/DLAT_READ
     JSB=MAX(MIN(NINT(ZJSB)+1, NLAT_READ), 1)

     ZJNB=(ZLATNB-LAT_START_READ)/DLAT_READ
     JNB=MAX(MIN(NINT(ZJNB)+1, NLAT_READ), 1)

     NTYPE(1:NVT)=0

     DO JJJ=JSB,JNB
     DO III=IWB,IEB
       ITYPE=INT(VEGET_READ(III,JJJ))
       IF (ITYPE /= IVAL_MISSING) NTYPE(ITYPE)=NTYPE(ITYPE)+1
     ENDDO
     ENDDO

! The most frequent veg. type in a grid box

     NMAX=0
     IF (ITSQR == 1) THEN ! "Odd" grid box
       DO ITYPE=1,NVT
         IF (NTYPE(ITYPE) > NMAX) THEN
           NMAX=NTYPE(ITYPE)
           MAXTYPE=ITYPE
         ENDIF
       ENDDO
       ITSQR=ITSQR+1
     ELSE ! "Even" grid box
       DO ITYPE=NVT,1,-1
         IF (NTYPE(ITYPE) > NMAX) THEN
           NMAX=NTYPE(ITYPE)
           MAXTYPE=ITYPE
         ENDIF
       ENDDO
       ITSQR=ITSQR-1
     ENDIF

! Veg. type, the second for frequency in a box grid

     NMAX2=0
     DO ITYPE=1,NVT
       IF (ITYPE /= MAXTYPE) THEN
         IF (NTYPE(ITYPE) > NMAX2) THEN
           NMAX2=NTYPE(ITYPE)
           MAXTYPE2=ITYPE
         ENDIF
       ENDIF
     ENDDO

     IF (MAXTYPE == NVT) THEN ! The most frequent veg. type in a grid box is glacier or water body
       IF (NMAX2 > 0) MAXTYPE=MAXTYPE2 ! If not water body veg. type esists, then this type is considerated as the most frequent veg. type
     ENDIF

     VEGET(I,J,1)=FLOAT(MAXTYPE)

! Fraction of vetegation type exept for glacier and water body

     IF (MAXTYPE < NVT-1) THEN

       NP=0
       DO ITYPE=1,NVT-2
         NP=NP+NTYPE(ITYPE)
       ENDDO

       IF (NP > 0) THEN

         ZNP=1./FLOAT(NP)
         DO ITYPE=1,NVT-2
           VEGET(I,J,ITYPE+1)=FLOAT(NTYPE(ITYPE))*ZNP
         ENDDO
         DO ITYPE=NVT-1,NVT
           VEGET(I,J,ITYPE+1)=0.
         ENDDO

       ELSE

         VEGET(I,J,2:NVT+1)=0.
         VEGET(I,J,MAXTYPE+1)=1.

       ENDIF

     ELSE ! Glacier or water body

       VEGET(I,J,2:NVT+1)=0.
       IF (NTYPE(NVT-1) > NTYPE(NVT)) THEN ! Glacier
         ITYPE=NVT-1
       ELSE ! Water body
         ITYPE=NVT
       ENDIF
       VEGET(I,J,ITYPE+1)=1.
       VEGET(I,J,1)=FLOAT(ITYPE)

     ENDIF

! Land-sea fraction

     NP=0
     DO ITYPE=1,NVT
       NP=NP+NTYPE(ITYPE)
     ENDDO

     IF (NP > 0) THEN
       ITYPE=NVT ! Water body veg. type
       LSF(I,J)=FLOAT(NTYPE(ITYPE))/FLOAT(NP) ! Fraction of water body (sea)
     ELSE
       LSF(I,J)=1.
     ENDIF

   ENDDO
   ENDDO

 ELSE ! Resolution of dataset is close to model grid resoltion

   NPOINT(:,:)=0
   NPOINT_VEGTYPE(:,:,:)=0
  
   DO J=1,NLAT_READ
   DO I=1,NLON_READ
     IF (LON_READ_MOD(I,J)/=VAL_MISSING.AND.LAT_READ_MOD(I,J)/=VAL_MISSING) THEN
       IMOD=NINT((LON_READ_MOD(I,J)-LONINI)/DLON)+1
       JMOD=NINT((LAT_READ_MOD(I,J)-LATINI)/DLAT)+1
       IF (IMOD>=1.AND.IMOD<=NLON.AND.JMOD>=1.AND.JMOD<=NLAT) THEN
         IF (VEGET_READ(I,J)/=IVAL_MISSING) THEN
           ITYPE=INT(VEGET_READ(I,J))
           NPOINT_VEGTYPE(IMOD,JMOD,ITYPE)=NPOINT_VEGTYPE(IMOD,JMOD,ITYPE)+1
           IF (VEGET_READ(I,J)<NVT-1) NPOINT(IMOD,JMOD)=NPOINT(IMOD,JMOD)+1 ! No glacier or water body
         ENDIF
       ENDIF
     ENDIF
   ENDDO 
   ENDDO 
  
   DO J=1,NLAT
   DO I=1,NLON
     IF (NPOINT(I,J)>0) THEN
       ZNP=1./FLOAT(NPOINT(I,J))
       DO ITYPE=1,NVT-2
         VEGET(I,J,ITYPE+1)=FLOAT(NPOINT_VEGTYPE(I,J,ITYPE))*ZNP
       ENDDO
       DO ITYPE=NVT-1,NVT
         VEGET(I,J,ITYPE+1)=0.
       ENDDO
       LSF(I,J)=0. ! Fraction of water body (sea) is equal 0
     ELSE
       VEGET(I,J,2:NVT+1)=0.
       IF (NPOINT_VEGTYPE(I,J,NVT-1) > NPOINT_VEGTYPE(I,J,NVT)) THEN ! Glacier
         ITYPE=NVT-1 
         LSF(I,J)=1. ! Fraction of water body (sea) is equal 0
       ELSE ! Water body
         ITYPE=NVT 
         LSF(I,J)=1. ! Fraction of water body (sea) is equal 1
       ENDIF 
       VEGET(I,J,ITYPE+1)=1.
     ENDIF
     IARRAY=MAXLOC(VEGET(I,J,2:NVT+1))
     ITYPE=IARRAY(1)
     VEGET(I,J,1)=FLOAT(ITYPE)
   ENDDO
   ENDDO

 ENDIF ! Correlation between dataset resolution and model grid resolution

RETURN

END SUBROUTINE DATASET_VEGETATION_CYCLOPES
!=======================================================================
subroutine dataset_t_climate(nlon, nlat, alon, alat, flag_glob, t_climate, t_amplitude)

! Procedure defines climate temperature and amplitude of climate temperature 
! over the model grid using global using 36 year (1979-2014) re-analysis data
! of IFS-ECMWF with resolution 0.75 degree (ERA-Interim dataset) 

! Input variables:

!      nlon: grid point number along axis of rotated longitude;
!      nlat: grid point number along axis of rotated latidude;
!      alon: geographical longitude of the grid points (degree, -180...180);
!      alat: geographical latitude of the grid points (degree, -90...90) 
! flag_orog: if 0, then not global input grid, if 1, then input grid is global not rotated

! Output variables:

!      t_climate: climate temperature (K) average in 36 year (1979-2014)
!      t_amplitude: amplitude of temperature (K) (maximum in 36 years minus minimum in 36 years)

! Data of IFS-ECMWF
! defined on a lat-lon regular grid, with resolution of 0.75 deg.
! lon 0 --> 359.25, lat -90 --> 90

implicit none

integer :: nlon, nlat, flag_glob
real, dimension(nlon,nlat) ::  alat, alon, lat_geo, lon_geo
real, dimension (:,:), allocatable :: tclim_read0, tamplid_read0, tclim_read, tamplid_read
real, dimension (:), allocatable :: lon_read, lat_read
character (len=50) :: file_input='ecmwf_ifs_1979_2014_t2_clim.bin'

real, dimension(nlon,nlat) :: t_climate, t_amplitude

! Work variables:

integer :: nx_read, ny_read
real :: dlon_read, dlat_read, lonini_read, latini_read, lonfin_read, latfin_read

! dlon_read, dlat_read: resolutions IFS-ECMWF data
! lonini_read, lonfin_read: longitude of the west and east extremes of IFS-ECMWF data
! latini_read, latfin_read: latitude of the north and south extremes of IFS-ECMWF data

real :: val_missing=-9999., lat_geo_min, lat_geo_max, lon_geo_min, lon_geo_max, &
        lat_start_read, lon_start_read
integer :: nlon_read, nlat_read, &
           istart, ifinish, jstart, jfinish, jini, jfin, j, jjj, i, iii, &
           ix_part, flag_period, jlon, jlat
integer, parameter :: nx_part = 2
integer, dimension(nx_part) :: iini, ifin

 lat_geo(:,:) = alat(:,:)
 lon_geo(:,:) = alon(:,:)
 do jlat=1,nlat
 do jlon=1,nlon
   if (abs(abs(lat_geo(jlon,jlat))-90.) < 1.e-10) lat_geo(jlon,jlat) = sign(89.999999, lat_geo(jlon,jlat))
   if (alon(jlon,jlat) < -180.) lon_geo(jlon,jlat) = alon(jlon,jlat)+360.
   if (alon(jlon,jlat) >  180.) lon_geo(jlon,jlat) = alon(jlon,jlat)-360.
 enddo
 enddo

! Coordinates of the centre of the extremes for reading

 lon_geo_min = minval(lon_geo)
 lon_geo_max = maxval(lon_geo)
 lat_geo_min = minval(lat_geo)
 lat_geo_max = maxval(lat_geo)

 istart = int(lon_geo_min*0.2)*5-10
 ifinish = int(lon_geo_max*0.2)*5+15
 jstart = max( int(lat_geo_min*0.2)*5-10, -90 )
 jfinish = min( int(lat_geo_max*0.2)*5+15, 90 )

 print *
 print *,'Reading ECMWF REA-Interim climate temperature database from ',file_input

 open (18, file=file_input, form='unformatted')

 read (18) ny_read, nx_read
 read (18) latini_read, lonini_read, dlat_read, dlon_read
 latfin_read=latini_read+float(ny_read-1)*dlat_read
 lonfin_read=lonini_read+float(nx_read-1)*dlon_read

 allocate (tclim_read0(nx_read,ny_read))
 allocate (tamplid_read0(nx_read,ny_read))

 read (18) tclim_read0
 read (18) tamplid_read0

 close (18)

! Control of read data domain along latitude axis

 if (flag_glob == 0) then

   if (float(jstart) >= latini_read.and.float(jfinish) <= latfin_read) then

! Domain extremes locate inside of latitude limits of dataset global domain

     jini = int((float(jstart)-latini_read)/dlat_read)+1
     jfin = int((float(jfinish)-latini_read)/dlat_read)+1

   else

! Domain extremes located outside of latitude limits of dataset global domain

     if (float(jstart) < latini_read.and.float(jfinish) <= latfin_read) then

       jini = 1
       jfin = int((float(jfinish)-latini_read)/dlat_read)+1
       if (abs(jstart) /= 90) then
         print *,' Caution: ECMWF ERA-Interim climate temperature dataset does not cover the south part of the model domain'
         print *,' Model grid south latitude extreme',float(jstart)
         print *,' Dataset grid south latitude extreme',latini_read
       endif

     elseif (float(jfinish) > latfin_read.and.float(jstart) >= latfin_read) then

       jini = int((float(jstart)-latini_read)/dlat_read)+1
       jfin = ny_read
       if (abs(jfinish) /= 90) then
         print *,' Caution: ECMWF ERA-Interim climate temperature dataset does not cover the north part of the model domain'
         print *,' Model grid north latitude extreme',float(jfinish)
         print *,' Dataset grid north latitude extrem',latfin_read
       endif

     else

! Domain covers the whole latitude range

       jini = 1
       jfin = ny_read
       if (abs(jstart) /= 90.or.abs(jfinish) /= 90) then
         print *,' Caution: ECMWF ERA-Interim climate temperature dataset does not cover the north and south parts of the model domain'
         print *,' Model grid latitude extremes',float(jstart),float(jfinish)
         print *,' Dataset grid latitude extremes',latini_read,latfin_read
       endif

     endif

   endif

 else

! Domain covers the whole latitude range

   jini = 1
   jfin = ny_read

 endif

 lat_start_read = latini_read+float(jini-1)*dlat_read

! Control of read data domain along longitude axis, possible cut and paste operation for periodic grid

 if (flag_glob == 0) then

   if (istart >= int(lonini_read).and.ifinish <= int(lonfin_read)) then

! Domain extremes locate inside of longitude limits of dataset global domain

     iini(1) = int((float(istart)-lonini_read)/dlon_read)+1
     ifin(1) = int((float(ifinish)-lonini_read)/dlon_read)+1

     iini(2) = 0
     ifin(2) = -1

     flag_period = 0

   else

! Domain extremes locate outside of longitude limits of dataset global domain -> cut and paste

     if (istart < int(lonini_read).and.ifinish <= int(lonfin_read)) then

       iini(1) = int((float(istart)+360.-lonini_read)/dlon_read)+1
       ifin(1) = nx_read

       iini(2) = 1
       ifin(2) = int((float(ifinish)-lonini_read)/dlon_read)+1
       if (ifin(2) <= 0) ifin(2)=int((float(ifinish)+360.-lonini_read)/dlon_read)+1

       flag_period = 1

     elseif (ifinish > int(lonfin_read).and.istart >= int(lonini_read)) then

       iini(1) = int((float(istart)-lonini_read)/dlon_read)+1
       if (iini(1) >= nx_read) iini(1)=int((float(istart)-360.-lonini_read)/dlon_read)+1
       ifin(1) = nx_read

       iini(2) = 1
       ifin(2) = int((float(ifinish)-360.-lonini_read)/dlon_read)+1

       flag_period = 1

     else

! Domain covers the whole longitude circle

       iini(1) = 1
       ifin(1) = nx_read

       iini(2) = 0
       ifin(2) = -1

       flag_period = 2

     endif

   endif

 else ! global output grid

! Domain covers the whole longitude circle

   if (abs(alon(1,1)-lonini_read) < 10.) then

! Input (dataset) and output (model) global grid have the same longitude origin

     iini(1) = 1
     ifin(1) = nx_read

     IINI(2) = 0
     IFIN(2) = -1

     flag_period = 2

   else

! Input (dataset) and output (model) global grid have the longitude origin shited by 180 degees

     iini(1) = int((alon(1,1)-lonini_read)/dlon_read)+1
     ifin(1) = nx_read

     iini(2) = 1
     ifin(2) = iini(1)-1

     flag_period = 3

   endif

 endif

 lon_start_read = lonini_read+float(iini(1)-1)*dlon_read

 nlat_read = jfin-jini+1
 nlon_read = ifin(1)-iini(1)+1
 if (flag_period == 1 ) then
   nlon_read = nlon_read+ifin(2)-iini(2)+1
 endif
 if (flag_period == 2.or.flag_period == 3) then
   nlon_read = nx_read+2
   if (flag_period == 2) lon_start_read = lonini_read-dlon_read
   if (flag_period == 3) lon_start_read = lon_start_read-dlon_read
 endif

 if (lon_start_read > 180.) lon_start_read=lon_start_read-360.

 allocate(tclim_read(nlon_read,nlat_read))
 allocate(tamplid_read(nlon_read,nlat_read))
 allocate(lon_read(nlon_read))
 allocate(lat_read(nlat_read))
 tclim_read(:,:) = val_missing
 tamplid_read(:,:) = val_missing
 lon_read(:) = val_missing
 lat_read(:) = val_missing

 j = 0
 do jjj=jini, jfin
   j = j+1

   i = 0
   if (flag_period == 2.or.flag_period == 3) i=1
   do ix_part=1,nx_part
     do iii=iini(ix_part), ifin(ix_part)
       i = i+1
       tclim_read(i,j) = tclim_read0(iii,jjj)
       tamplid_read(i,j) = tamplid_read0(iii,jjj)
       lon_read(i) = lon_start_read+float(i-1)*dlon_read
     enddo
   enddo ! ix_part
   if (flag_period == 2.or.flag_period == 3) then
     tclim_read(1,j) = tclim_read(nlon_read-1,j)
     tclim_read(nlon_read,j) = tclim_read(2,j)
     tamplid_read(1,j) = tamplid_read(nlon_read-1,j)
     tamplid_read(nlon_read,j) = tamplid_read(2,j)
   endif

   lat_read(j) = lat_start_read+float(j-1)*dlat_read
 enddo ! jjj
 if (flag_period == 2.or.flag_period == 3) then
   lon_read(1) = lon_start_read
   lon_read(nlon_read) = lon_start_read+float(nlon_read-1)*dlon_read
 endif

! Interpolation of the read data to the model grid using model rotated coordinate system

 lon_geo(:,:) = alon(:,:)
 do jlat=1,nlat
 do jlon=1,nlon
   if (alon(jlon,jlat) < -180.) lon_geo(jlon,jlat) = alon(jlon,jlat)+360.
   if (alon(jlon,jlat) >  180.) lon_geo(jlon,jlat) = alon(jlon,jlat)-360.
 enddo
 enddo
 if (flag_glob == 1.and.flag_period == 2) then
! Dataset longitude: 0..360
   lon_geo(:,:) = alon(:,:)
 endif

 call interp_spline_2d(tclim_read, nlon_read, nlat_read, lon_read, lat_read, &
                       lon_geo, lat_geo, nlon*nlat, t_climate, 1.)

 call interp_spline_2d(tamplid_read, nlon_read, nlat_read, lon_read, lat_read, &
                       lon_geo, lat_geo, nlon*nlat, t_amplitude, 1.)

return
end subroutine dataset_t_climate
!=======================================================================
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
INTEGER, DIMENSION(2) :: IUNIT=(/11, 12/)
INTEGER, DIMENSION(NLON,NLAT) :: NPOINT_LAI, NPOINT_FVEG
INTEGER, PARAMETER :: NX_PART=2
INTEGER, DIMENSION(NX_PART) :: IINI, IFIN
INTEGER :: DAYR, DAYL, IERR_OPEN, NLON_READ, NLAT_READ, ISTART, IFINISH, JSTART, JFINISH,&
 IFILE, IPAR, I, J, IR, JR, IX_PART, III, JJJ, JINI, JFIN, FLAG_PERIOD=0, IND, IMOD, JMOD
REAL :: LON_GEO_MIN, LON_GEO_MAX, LAT_GEO_MIN, LAT_GEO_MAX, LON_START_READ, LAT_START_READ, &
 WEIGH_DAYR, WEIGH_DAYL

 LAI_MAX=7.5

 PRINT *
 PRINT *,'Reading and processing data about vegetation parameters: LAI and FVEG'

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

 LON_GEO_MIN=MINVAL(ALON)
 LON_GEO_MAX=MAXVAL(ALON)

 LAT_GEO(:,:)=ALAT(:,:)
 LON_GEO(:,:)=ALON(:,:)

 IF (LON_GEO_MIN >= 0..AND.LON_GEO_MAX > 0.) THEN
   LON_GEO(:,:)=ALON(:,:)
 ELSE
   DO J=1,NLAT
   DO I=1,NLON
     IF (ALON(I,J) < -180.) LON_GEO(I,J) = ALON(I,J)+360.
     IF (ALON(I,J) >  180.) LON_GEO(I,J) = ALON(I,J)-360.
   ENDDO
   ENDDO
 ENDIF

 DO J=1,NLAT
 DO I=1,NLON
   IF (ABS(ABS(LAT_GEO(I,J))-90.) < 1.e-10) LAT_GEO(I,J) = SIGN(89.999999, LAT_GEO(I,J))
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
         print *,' Model grid south latitude extremes',FLOAT(JSTART),FLOAT(JFINISH)
         print *,' Dataset grid south latitude extremes',LATINI_READ,LATFIN_READ
       ENDIF

     ENDIF

   ENDIF

 ELSE

! Domain covers the whole latitude range

   JINI = 1
   JFIN = NY_READ
   print *,' Caution: LAI, Vegetation fraction datasets does not cover the north and the south parts of the model domain'
   print *,' Model grid south latitude extremes',LAT_GEO(1,1),LAT_GEO(1,NLAT)
   print *,' Dataset grid south latitude extremes',LATINI_READ,LATFIN_READ

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

 IF (LON_START_READ > 180..AND.LON_GEO_MIN < 180.) LON_START_READ=LON_START_READ-360.

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
               IF (FLAG_GLOB == 0.AND.LON_READ_GEO(I,J)>180..AND.LON_GEO_MIN < 180.) LON_READ_GEO(I,J)=LON_READ_GEO(I,J)-360. !???
             ENDDO
           ENDDO ! IX_PART
           IF (FLAG_PERIOD == 2.OR.FLAG_PERIOD == 3) THEN
             LAI_READ(1,J,IFILE)=LAI_READ(NLON_READ-1,J,IFILE)
             LAI_READ(NLON_READ,J,IFILE)=LAI_READ(2,J,IFILE)
             LON_READ_GEO(1,J)=LON_START_READ
             LON_READ_GEO(NLON_READ,J)=LON_START_READ+FLOAT(NLON_READ-1)*DLON_READ
             IF (FLAG_GLOB == 0) THEN
               IF (LON_READ_GEO(1,J)>180..AND.LON_GEO_MIN < 180.) LON_READ_GEO(1,J)=LON_READ_GEO(1,J)-360. !???
               IF (LON_READ_GEO(NLON_READ,J)>180..AND.LON_GEO_MIN < 180.) LON_READ_GEO(NLON_READ,J)=LON_READ_GEO(NLON_READ,J)-360. !???
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

       PRINT *,'Not found input file ',FILEINPUT(IFILE),', stop'
       STOP

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
!=======================================================================
subroutine soil_veg_lsf_correct(nlon, nlat, nst, nvt, veget_classif, lat_geo, lsf_ecmwf, &
 soilmap, vegmap, lsf, t_climate, t_amplitude)

! Procedure makes some correction of soil, vegetation and land-sea fraction
! maps defined in the model grid, to coordinate among themselves

! Input variables:

!      nlon: grid point number along axis of rotated longitude;
!      nlat: grid point number along axis of rotated latidude;
!      nst: total number of soil types;
!      nvt: total number of vegetation types;
!      veget_classif: index of vegetetion classification 1 - GLC1990, 2 - GLC2000
!      lat_geo: geographical latitude of the grid points (degree, -90...90);
!      lsf_ecwmf: land-sea fraction according to IFS-ECMWF data, (proportion, 1.0 - sea, 0.0 -land)
!      t_climate: climate temperature (K) average in 36 year (1979-2014)
!      t_amplitude: amplitude of temperature (K) (maximum in 36 years minus minimum in 36 years)

! Input and output variables: 

!      soilmap:  soil types (classification FAO) defined on the grid points, 3D array,
!                1-st and 2-nd array index are grid point index, 3-rd index is soil type number plus 1,
!                if 3-rd index has value 1, then the variable means dominate soil type,
!                values 2...nst+1 of 3-rd index means fraction (of area, proportion) of the each of nst soil types;

!      vegmap:  vegetation types (classification GLC1990 or GLC2000) defined on the grid points, 3D array,
!                1-st and 2-nd array index are grid point index, 3-rd index is soil type number plus 1,
!                if 3-rd index has value 1, then the variable means dominate soil type,
!                values 2...nvt+1 of 3-rd index means fraction (of area, proportion) of the each of nvt vegetation types;

!     lsf: land-sea fraction (proportion, 1.0 - sea, 0.0 -land)

implicit none

integer :: nlon, nlat, nst, nvt, veget_classif
real, dimension(nlon,nlat) :: lat_geo, lsf_ecmwf, lsf, t_climate, t_amplitude
real, dimension(nlon,nlat,nst+1) :: soilmap
real, dimension(nlon,nlat,nvt+1) :: vegmap

! Work variables:

integer, dimension(nlon,nlat) :: imap1
integer, dimension(100) :: nist
integer :: jlon,jlat, i, j, iini, ifin, jini, jfin, itype, itypemax, ntypemax, &
 icondition, itype_correct, icorrect, idefault, iter, &
 iurbveg, iurbsoil, iantarctic_soil, iantarctic_veg, &
 iglacier_soil, iglacier_veg
real :: lat_geo_min, zfmax

! --------------------------------------------------------------------------

! Correction of vegetation type map, defined on the grid, "expanding" vegetation over coastline to sea
! for the purpose to have assured definition of vegetation along costline also 

 itype_correct = nvt ! water body
 if (veget_classif == 1) idefault = 5 ! default vegetation type - Mixed cover
 if (veget_classif == 2) idefault = 9 ! default vegetation type - Mosaic of tree cover and other natural vegetation (Crop component possible) 
 imap1(1:nlon,1:nlat) = nint(vegmap(1:nlon,1:nlat,1))

 do jlat=1,nlat
 do jlon=1,nlon

   if (imap1(jlon,jlat) == nvt.and.lsf(jlon,jlat) <= 0.5) then ! vegetation type is not defined in coastline zone 

     do iter=1,nlon/5

       jini = max(jlat-iter,1)
       jfin = min(jlat+iter,nlat)
       iini = max(jlon-iter,1)
       ifin = min(jlon+iter,nlon)
       nist(1:nvt) = 0
       do j=jini,jfin
       do i=iini,ifin
         if (i == jlon.and.j == jlat) then
           continue
         else
           nist(imap1(i,j)) = nist(imap1(i,j))+1
         endif
       enddo
       enddo
       itypemax = 0
       ntypemax = 0
       do itype=1,nvt
         if (itype == itype_correct) cycle
         if (nist(itype) > itypemax) then
           itypemax = nist(itype)
           ntypemax = itype
         endif
       enddo

       if (ntypemax > 0) exit

     enddo

     if (ntypemax == 0) ntypemax = idefault

     vegmap(jlon,jlat,1:nvt+1) = 0.
     vegmap(jlon,jlat,1) = float(ntypemax)
     vegmap(jlon,jlat,ntypemax+1) = 1.

   endif

 enddo
 enddo

! --------------------------------------------------------------------------

! Correction of soil type map, defined on the grid using vegetation map,
! becuase of vegetation database has a very high resolution
! in comparison with the soil database

 idefault = 16 ! default soil type - Podzols

! Case of GCL1990 vegetation dataset

 if (veget_classif == 1) then

! 1) "Expanding" soil over coastline to sea
! for the purpose to have assured definition of vegetation along costline also 
! 2) Definition of a soil type different from "water body" in the ponts when
! soil map has "water body", but vegetation map has not "water body"
! 3) Definition of a soil type different from "glaciers" in the ponts when
! soil map has "glaciers", but vegetation map has not "bare soil"

   do icorrect=1,3

     imap1(1:nlon,1:nlat) = nint(soilmap(1:nlon,1:nlat,1))

     do jlat=1,nlat
     do jlon=1,nlon
  
! Condition of a correction

       icondition=0

       if (icorrect == 1) then
         itype_correct=nst ! "water body" is soil type that must be corrected
         if (imap1(jlon,jlat) == nst.and.lsf(jlon,jlat) <= 0.5) icondition=1
       endif
       if (icorrect == 2) then
         itype_correct=nst ! "water body" soil type that must be corrected
         if (imap1(jlon,jlat) == nst.and.nint(vegmap(jlon,jlat,1)) /= nvt) icondition=1
       endif
       if (icorrect == 3) then
         itype_correct=nst-1 ! "glacier" soil type that must be corrected
         if (imap1(jlon,jlat) == nst-1.and.nint(vegmap(jlon,jlat,1)) /= 12) icondition=1
       endif

       if (icondition == 1) then
  
         do iter=1,nlon/5
    
           jini = max(jlat-iter,1)
           jfin = min(jlat+iter,nlat)
           iini = max(jlon-iter,1)
           ifin = min(jlon+iter,nlon)
           nist(1:nst) = 0
           do j=jini,jfin
           do i=iini,ifin
             if (i == jlon.and.j == jlat) then
               continue
             else
               nist(imap1(i,j)) = nist(imap1(i,j))+1
             endif
           enddo
           enddo
           itypemax = 0
           ntypemax = 0
           do itype=1,nst
             if (itype == itype_correct) cycle
             if (nist(itype) > itypemax) then
               itypemax = nist(itype)
               ntypemax = itype
             endif
           enddo
    
           if (ntypemax > 0) exit
    
         enddo
  
         if (ntypemax == 0) ntypemax = idefault
  
         soilmap(jlon,jlat,1:nst+1) = 0.
         soilmap(jlon,jlat,1) = float(ntypemax)
         soilmap(jlon,jlat,ntypemax+1) = 1.
  
       endif
  
     enddo
     enddo

   enddo ! Various corrections

!! 4) Definition of water body points in soil map,
!! if this points are not "water body" in vegetation map
!
!   do jlat=1,nlat
!   do jlon=1,nlon
!
!     if (nint(vegmap(jlon,jlat,1)) == nvt.and.nint(soilmap(jlon,jlat,1)) /= nst) then ! vegetation type is not defined, but soil type is defined
!       soilmap(jlon,jlat,1:nst+1) = 0.
!       soilmap(jlon,jlat,1) = nst
!       soilmap(jlon,jlat,nst+1) = 1.
!     endif
!
!   enddo
!   enddo

 endif ! veget. classification 1 (GLC1990)

! Case of GCL2000 vegetation dataset

 if (veget_classif == 2) then 

! 1) "Expanding" soil over coastline to sea
! for the purpose to have assured definition of vegetation along costline also 
! 2) Definition of a soil type different from "water body" in the ponts when
! soil map has "water body", but vegetation map has not "water body"
! 3) Definition of a soil type different from "glaciers" in the ponts when
! soil map has "glaciers", but vegetation map has not "Snow or Ice"

   do icorrect=1,3

     imap1(1:nlon,1:nlat) = nint(soilmap(1:nlon,1:nlat,1))

     do jlat=1,nlat
     do jlon=1,nlon
  
! Condition of a correction

       icondition=0

       if (icorrect == 1) then
         itype_correct=nst ! "water body" is soil type that must be corrected
         if (imap1(jlon,jlat) == nst.and.lsf(jlon,jlat) <= 0.5) icondition=1
       endif
       if (icorrect == 2) then
         itype_correct=nst ! "water body" soil type that must be corrected
         if (imap1(jlon,jlat) == nst.and.nint(vegmap(jlon,jlat,1)) /= nvt) icondition=1
       endif
       if (icorrect == 3) then
         itype_correct=nst-1 ! "glacier" soil type that must be corrected
         if (imap1(jlon,jlat) == nst-1.and.nint(vegmap(jlon,jlat,1)) /= nvt-1) icondition=1
       endif

       if (icondition == 1) then
  
         do iter=1,nlon/5
    
           jini = max(jlat-iter,1)
           jfin = min(jlat+iter,nlat)
           iini = max(jlon-iter,1)
           ifin = min(jlon+iter,nlon)
           nist(1:nst) = 0
           do j=jini,jfin
           do i=iini,ifin
             if (i == jlon.and.j == jlat) then
               continue
             else
               nist(imap1(i,j)) = nist(imap1(i,j))+1
             endif
           enddo
           enddo
           itypemax = 0
           ntypemax = 0
           do itype=1,nst
             if (itype == itype_correct) cycle
             if (nist(itype) > itypemax) then
               itypemax = nist(itype)
               ntypemax = itype
             endif
           enddo
    
           if (ntypemax > 0) exit
    
         enddo
  
         if (ntypemax == 0) ntypemax = idefault
  
         soilmap(jlon,jlat,1:nst+1) = 0.
         soilmap(jlon,jlat,1) = float(ntypemax)
         soilmap(jlon,jlat,ntypemax+1) = 1.
  
       endif
  
     enddo
     enddo

   enddo ! Various corrections

! 4) Definition of glacier points in soil map,
! if this points are not "Snow or Ice" in vegetation map;
!! 5) Definition of water body points in soil map,
!! if this points are not "water body" in vegetation map;

   do jlat=1,nlat
   do jlon=1,nlon

     if (nint(vegmap(jlon,jlat,1)) == nvt-1.and.nint(soilmap(jlon,jlat,1)) /= nst-1) then ! 2-nd correction
       soilmap(jlon,jlat,1:nst+1) = 0.
       soilmap(jlon,jlat,1) = nst-1
       soilmap(jlon,jlat,nst) = 1.
     endif

!     if (nint(vegmap(jlon,jlat,1)) == nvt.and.nint(soilmap(jlon,jlat,1)) /= nst) then ! 1-st correction
!       soilmap(jlon,jlat,1:nst+1) = 0.
!       soilmap(jlon,jlat,1) = nst
!       soilmap(jlon,jlat,nst+1) = 1.
!     endif

   enddo
   enddo

 endif ! veget. classification 2 (GLC2000)

! --------------------------------------------------------------------------

! Urban soil type introduced as Lithosols (mountain sols) on the base of vegetation urban type

 if (veget_classif == 1) iurbveg = 13
 if (veget_classif == 2) iurbveg = 20
 iurbsoil = 8

 do jlat=1,nlat
 do jlon=1,nlon
   if (vegmap(jlon,jlat,iurbveg+1) > 0.01) then
     do itype=2,nst+1
       soilmap(jlon,jlat,itype) = soilmap(jlon,jlat,itype)*(1.-vegmap(jlon,jlat,iurbveg+1))
     enddo
     soilmap(jlon,jlat,iurbsoil+1) = soilmap(jlon,jlat,iurbsoil+1)+vegmap(jlon,jlat,iurbveg+1)
     zfmax = 0.
     do itype=1,nst
       if (soilmap(jlon,jlat,itype+1) > zfmax) then
         itypemax = itype
         zfmax = soilmap(jlon,jlat,itype+1)
       endif
     enddo
     soilmap(jlon,jlat,1) = float(itypemax)
   endif
 enddo
 enddo

! --------------------------------------------------------------------------

! Correction of the coastline of Antarctica with including of perennial sea ice
! with using land-sea fraction from IFS-ECMWF data interpolated to the model grid

 lat_geo_min=minval(lat_geo(:,:))

! New Antarctic --->
! if (lat_geo_min < -60.) then ! Antarctic area
!
!   if (veget_classif == 1) iantarctic_veg=12 ! "bare soil"
!   if (veget_classif == 2) iantarctic_veg=nvt-1 ! "snow or ice"
!
!   do jlat=1,nlat
!   do jlon=1,nlon
!     if (lat_geo(jlon,jlat) < -60..and.lsf(jlon,jlat) > lsf_ecmwf(jlon,jlat)) then
!       lsf(jlon,jlat) = lsf_ecmwf(jlon,jlat)
!     endif
!! Antarctic is covered by glacier with the exception of some coast zone
!     if (lat_geo(jlon,jlat) < -60..and.lsf(jlon,jlat) < 0.5) then
!       iantarctic_soil=nst-1 ! Glacier
!       if (t_climate(jlon,jlat) > 248..and.t_amplitude(jlon,jlat) < 55.) then
!!       if (t_climate(jlon,jlat) > 248.15) then
!! Antarctic "warm" coast zone: average year surface temperature is over then -25 degree,
!! and extreme year surface temperature amplitude is less then 60 degree
!         iantarctic_soil=18 ! Regosols: crust of weathering, rock
!       endif
!       soilmap(jlon,jlat,:) = 0.
!       soilmap(jlon,jlat,1) = float(iantarctic_soil)
!       soilmap(jlon,jlat,iantarctic_soil+1) = 1.
!       vegmap(jlon,jlat,:) = 0.
!       vegmap(jlon,jlat,1) = float(iantarctic_veg)
!       vegmap(jlon,jlat,iantarctic_veg+1) = 1.
!      endif
!   enddo
!   enddo
!
! endif ! Antarctic area
! <---

! Old Antarctic --->
 if (lat_geo_min < -60.) then ! Antarctic area

   iglacier_soil=nst-1 ! "glacier"
   if (veget_classif == 1) iglacier_veg=12 ! "bare soil"
   if (veget_classif == 2) iglacier_veg=nvt-1 ! "snow or ice"

   do jlat=1,nlat
   do jlon=1,nlon
     if (lat_geo(jlon,jlat) < -60..and.lsf(jlon,jlat) > lsf_ecmwf(jlon,jlat)) then
       lsf(jlon,jlat) = lsf_ecmwf(jlon,jlat)
     endif
     if (lat_geo(jlon,jlat) < -60..and.lsf(jlon,jlat) < 0.5) then
       soilmap(jlon,jlat,:) = 0.
       soilmap(jlon,jlat,1) = float(iglacier_soil)
       soilmap(jlon,jlat,iglacier_soil+1) = 1.
       vegmap(jlon,jlat,:) = 0.
       vegmap(jlon,jlat,1) = float(iglacier_veg)
       vegmap(jlon,jlat,iglacier_veg+1) = 1.
      endif
   enddo
   enddo

 endif ! Antarctic area
! <---


! --------------------------------------------------------------------------

return
end subroutine soil_veg_lsf_correct
!=======================================================================
subroutine definition_soil_phys_param(soil_fao,nlevel_soil,level_soil,nlon,nlat, &
 soil_qmax, soil_qmin, soil_c, soil_rho, soil_psi, soil_k, soil_par_b, soil_par_c, soil_qrel_wilt, soil_qrel_ref, &
 soil_water_table, soil_albedo_dry, soil_albedo_wet, soil_emiss_1_dry, soil_emiss_1_wet, soil_emiss_2_dry, soil_emiss_2_wet)

! Procedure defines physical parameters of soil using soil types (classification FAO)
! defined on the grid points

! Internal basic variables:

!      nst: total number soil soil types;

! Input variables:

!      soil_fao: soil types (classification FAO, 31 types) defined on the grid points, 3D array,
!                1-st and 2-nd array index are grid point index, 3-rd index is soil type number plus 1,
!                if 3-rd index has value 1, then the variable means dominate soil type,
!                values 2-32 of 3-rd index means fraction (of area, proportion) of the each of 31 soil types;

!      nlevel_soil: number of soil levels;
!      level_soil: soil levels values (m);
!      nlon: grid point number along axis of rotated longitude;
!      nlat: grid point number along axis of rotated latidude.

! Output variables:

!      soil_qmax: maximum specific volumetric soil water content (m^3/m^3) at all soil levels, 3D array;
!      soil_qmin: minimum specific volumetric soil water content (m^3/m^3) at all soil levels, 3D array;
!      soil_c: dry soil thermal capacity (J/K/m^3) at all soil levels, 3D array;
!      soil_rho: dry soil density (kg/m^3) at all soil levels, 3D array;
!      soil_psi: idraulic potential of saturated soil (m) at all soil levels, 3D array;
!      soil_k: idraulic conductivity of saturated soil (m/s) at all soil levels, 3D array;
!      soil_par_b: parameter in soil water transport equation (dimensionless) at all soil levels, 3D array;
!      soil_par_c: parameter in soil water transport equation (dimensionless) at all soil levels, 3D array;
!      soil_qrel_wilt: relative soil water content (proportion) at which the vegetation wilt at all soil levels, 3D array;
!      soil_qrel_ref: relative soil water content (proportion) at which the vegetation stop evapotraspiration increase at all soil levels, 3D array;
!      soil_water_table: approximate depth of water table top (m), 2D array;
!      soil_albedo_dry: dry soil albedo (proportion), 2D array;
!      soil_albedo_wet: wet soil albedo (proportion), 2D array;
!      soil_emiss_1_dry: dry soil emissivity in broadband window (proportion), 2D array;
!      soil_emiss_1_wet: wet soil emissivity in broadband window (proportion), 2D array;
!      soil_emiss_2_dry: dry soil emissivity in 8-12 micron window (proportion), 2D array;
!      soil_emiss_2_wet: wet soil emissivity in 8-12 micron window (proportion), 2D array.

implicit none

integer, parameter :: nst=31

! Input:

integer :: nlon, nlat, nlevel_soil
real, dimension(nlon,nlat,nst+1) :: soil_fao
real, dimension(nlevel_soil) :: level_soil

! Ouput:

real, dimension(nlon,nlat,nlevel_soil) :: soil_qmax, soil_qmin, soil_c, soil_rho, &
 soil_psi, soil_k, soil_par_b, soil_par_c, soil_qrel_wilt, soil_qrel_ref
real, dimension(nlon,nlat) :: soil_water_table, soil_albedo_dry, soil_albedo_wet, &
 soil_emiss_1_dry, soil_emiss_1_wet, soil_emiss_2_dry, soil_emiss_2_wet

! Work:

integer, parameter :: soil_horizon_number=5, ival_missing=-9999
real, dimension(soil_horizon_number,nst) :: horizon_bottom, &
 qgmax, qgmin, rogdr, cgdry, psig0, ykghy, bghyd, cghyd, qg_rel_wilt, qg_rel_refw
real, dimension(nst) :: water_table_depth, soildalbed, soilwalbed, &
 soildemis1, soilwemis1, soildemis2, soilwemis2
real, dimension(nlevel_soil,nst) :: qgmax_lev, qgmin_lev, rogdr_lev, cgdry_lev, &
 psig0_lev, ykghy_lev, bghyd_lev, cghyd_lev, qg_rel_wilt_lev, qg_rel_refw_lev
real, dimension(nlon,nlat) :: zzz1, zzz2
integer :: itype, nsmooth, ilev, ihorizon, i, ip, jp, k, n11
real :: zfrac

! Definition of phisical parameters for various soil type at all soil horizons,
! various soil type may have variuos horizon number (maximum 5) with
! various thikness

! Depth (m) of bottom of soil horisons
!                                Horizon
!                   1       2       3       4       5          Soil type
data &
 horizon_bottom/   0.10,   0.40,   0.90,   1.50,   100., &  !   1 Acrisols
                   0.05,   0.30,   0.80,   1.50,   100., &  !   2 Cambisols
                   0.05,   0.80,   1.50,   2.00,   100., &  !   3 Chernozems
                   0.10,   0.30,   0.70,   1.40,   100., &  !   4 Podzoluvisols
                   0.10,   0.50,   100., -9999., -9999., &  !   5 Rendzinans
                   0.70,   2.00,    80.,   100., -9999., &  !   6 Ferrasols
                   0.20,   0.50,   1.50,   100., -9999., &  !   7 Gleysols
                   0.50,   1.00,   2.00,   100., -9999., &  !   8 Phaeozems
                   0.10,   100., -9999., -9999., -9999., &  !   9 Lithosols
                   0.20,   1.50,   100., -9999., -9999., &  !  10 Fluvisols
                   0.20,   0.60,   1.60,   100., -9999., &  !  11 Kastanozems
                   0.05,   0.20,   0.80,   1.50,   100., &  !  12 Luvisols
                   0.20,   0.80,   1.50,   100., -9999., &  !  13 Greyzems
                   0.40,   1.30,   2.00,   100., -9999., &  !  14 Nitosols
                   0.80,   1.20,   2.00,   100., -9999., &  !  15 Histosols
                   0.10,   0.40,   0.80,   1.40,   100., &  !  16 Podzols
                   0.80,   2.00,   100., -9999., -9999., &  !  17 Arenosols
                   0.10,   1.00,   100., -9999., -9999., &  !  18 Regosols
                   0.20,   0.50,   1.00,   100., -9999., &  !  19 Solonetz
                   0.30,   1.50,   100., -9999., -9999., &  !  20 Andosols
                   0.20,   0.50,   100., -9999., -9999., &  !  21 Rankers
                   0.20,   0.50,   1.00,   100., -9999., &  !  22 Vertisols
                   0.10,   2.00,   100., -9999., -9999., &  !  23 Planosols
                   0.50,   0.80,   1.50,   100., -9999., &  !  24 Xerosols
                   0.20,   1.00,   100., -9999., -9999., &  !  25 Yermosols
                   0.30,   1.00,   100., -9999., -9999., &  !  26 Solonchaks
                   2.00,   100., -9999., -9999., -9999., &  !  27 Dunes or shifting sands
                   1.00,   100., -9999., -9999., -9999., &  !  28 Salt flats
                   1.00,   100., -9999., -9999., -9999., &  !  29 Rock debris or desert detritus
                   100., -9999., -9999., -9999., -9999., &  !  30 Glaciers
                   100., -9999., -9999., -9999., -9999./    !  31 Water body

! Maximum volumetric water content of the soil (m**3/m**3) (qgmax*row is soil porosity)
!                                Horizon
!                   1       2       3       4       5          Soil type
data qgmax/       0.600,  0.480,  0.470,  0.490,  0.420, &  !   1 Acrisols
                  0.550,  0.400,  0.450,  0.490,  0.400, &  !   2 Cambisols
                  0.550,  0.500,  0.450,  0.400,  0.400, &  !   3 Chernozems
                  0.600,  0.400,  0.450,  0.420,  0.400, &  !   4 Podzoluvisols
                  0.450,  0.400,  0.400, -9999., -9999., &  !   5 Rendzinans
                  0.500,  0.500,  0.450,  0.420, -9999., &  !   6 Ferrasols
                  0.450,  0.400,  0.480,  0.420, -9999., &  !   7 Gleysols
                  0.450,  0.480,  0.490,  0.400, -9999., &  !   8 Phaeozems
                  0.450,  0.420, -9999., -9999., -9999., &  !   9 Lithosols
                  0.500,  0.450,  0.420, -9999., -9999., &  !  10 Fluvisols
                  0.420,  0.450,  0.400,  0.400, -9999., &  !  11 Kastanozems
                  0.550,  0.400,  0.480,  0.450,  0.400, &  !  12 Luvisols
                  0.450,  0.480,  0.490,  0.400, -9999., &  !  13 Greyzems
                  0.480,  0.490,  0.490,  0.420, -9999., &  !  14 Nitosols
                  0.900,  0.490,  0.490,  0.420, -9999., &  !  15 Histosols
                  0.600,  0.400,  0.480,  0.420,  0.400, &  !  16 Podzols
                  0.400,  0.430,  0.400, -9999., -9999., &  !  17 Arenosols
                  0.430,  0.460,  0.400, -9999., -9999., &  !  18 Regosols
                  0.410,  0.470,  0.490,  0.420, -9999., &  !  19 Solonetz
                  0.500,  0.700,  0.400, -9999., -9999., &  !  20 Andosols
                  0.450,  0.400,  0.400, -9999., -9999., &  !  21 Rankers
                  0.350,  0.400,  0.450,  0.420, -9999., &  !  22 Vertisols
                  0.450,  0.470,  0.420, -9999., -9999., &  !  23 Planosols
                  0.450,  0.480,  0.460,  0.400, -9999., &  !  24 Xerosols
                  0.450,  0.420,  0.400, -9999., -9999., &  !  25 Yermosols
                  0.450,  0.480,  0.420, -9999., -9999., &  !  26 Solonchaks
                  0.400,  0.400, -9999., -9999., -9999., &  !  27 Dunes or shifting sands
                  0.420,  0.400, -9999., -9999., -9999., &  !  28 Salt flats
                  0.420,  0.400, -9999., -9999., -9999., &  !  29 Rock debris or desert detritus
                  1.000, -9999., -9999., -9999., -9999., &  !  30 Glaciers
                  1.000, -9999., -9999., -9999., -9999./    !  31 Water body

! Minimum volumetric water content of the soil (m**3/m**3)
!                                Horizon
!                   1       2       3       4       5          Soil type
data qgmin/       0.150,  0.145,  0.140,  0.142,  0.025, &  !   1 Acrisols
                  0.100,  0.030,  0.060,  0.120,  0.020, &  !   2 Cambisols
                  0.100,  0.100,  0.080,  0.100,  0.080, &  !   3 Chernozems
                  0.150,  0.050,  0.130,  0.110,  0.050, &  !   4 Podzoluvisols
                  0.050,  0.020,  0.050, -9999., -9999., &  !   5 Rendzinans
                  0.100,  0.150,  0.140,  0.100, -9999., &  !   6 Ferrasols
                  0.040,  0.020,  0.140,  0.100, -9999., &  !   7 Gleysols
                  0.060,  0.100,  0.110,  0.050, -9999., &  !   8 Phaeozems
                  0.030,  0.020, -9999., -9999., -9999., &  !   9 Lithosols
                  0.100,  0.140,  0.025, -9999., -9999., &  !  10 Fluvisols
                  0.080,  0.130,  0.110,  0.050, -9999., &  !  11 Kastanozems
                  0.100,  0.050,  0.130,  0.120,  0.050, &  !  12 Luvisols
                  0.060,  0.100,  0.110,  0.050, -9999., &  !  13 Greyzems
                  0.135,  0.145,  0.140,  0.025, -9999., &  !  14 Nitosols
                  0.100,  0.140,  0.140,  0.025, -9999., &  !  15 Histosols
                  0.150,  0.050,  0.140,  0.110,  0.050, &  !  16 Podzols
                  0.020,  0.120,  0.020, -9999., -9999., &  !  17 Arenosols
                  0.100,  0.120,  0.020, -9999., -9999., &  !  18 Regosols
                  0.050,  0.130,  0.140,  0.050, -9999., &  !  19 Solonetz
                  0.020,  0.010,  0.020, -9999., -9999., &  !  20 Andosols
                  0.050,  0.020,  0.050, -9999., -9999., &  !  21 Rankers
                  0.100,  0.130,  0.140,  0.100, -9999., &  !  22 Vertisols
                  0.080,  0.140,  0.100, -9999., -9999., &  !  23 Planosols
                  0.120,  0.140,  0.130,  0.100, -9999., &  !  24 Xerosols
                  0.020,  0.050,  0.050, -9999., -9999., &  !  25 Yermosols
                  0.120,  0.140,  0.050, -9999., -9999., &  !  26 Solonchaks
                  0.020,  0.100, -9999., -9999., -9999., &  !  27 Dunes or shifting sands
                  0.050,  0.100, -9999., -9999., -9999., &  !  28 Salt flats
                  0.030,  0.100, -9999., -9999., -9999., &  !  29 Rock debris or desert detritus
                  0.000, -9999., -9999., -9999., -9999., &  !  30 Glaciers
                  0.000, -9999., -9999., -9999., -9999./    !  31 Water body

! Dry soil density (kg/m**3)
!                                Horizon
!                   1       2       3       4       5          Soil type
data rogdr/        900.,  1100.,  1400.,  1600.,  2300., &  !   1 Acrisols
                   500.,   900.,  1300.,  1600.,  2300., &  !   2 Cambisols
                   500.,  1000.,  1300.,  1500.,  2300., &  !   3 Chernozems
                   500.,  1200.,  1500.,  1500.,  2300., &  !   4 Podzoluvisols
                   500.,  1500.,  2300., -9999., -9999., &  !   5 Rendzinans
                  1100.,  1300.,  1600.,  2300., -9999., &  !   6 Ferrasols
                  1100.,  1400.,  1600.,  2300., -9999., &  !   7 Gleysols
                  1000.,  1400.,  1600.,  2300., -9999., &  !   8 Phaeozems
                  1200.,  2300., -9999., -9999., -9999., &  !   9 Lithosols
                  1000.,  1400.,  2300., -9999., -9999., &  !  10 Fluvisols
                  1100.,  1300.,  1500.,  2300., -9999., &  !  11 Kastanozems
                   500.,  1100.,  1400.,  1500.,  2300., &  !  12 Luvisols
                  1000.,  1400.,  1600.,  2300., -9999., &  !  13 Greyzems
                  1100.,  1400.,  1600.,  2300., -9999., &  !  14 Nitosols
                   300.,  1400.,  1600.,  2300., -9999., &  !  15 Histosols
                   500.,  1300.,  1600.,  1500.,  2300., &  !  16 Podzols
                  1100.,  1500.,  2300., -9999., -9999., &  !  17 Arenosols
                  1200.,  1600.,  2300., -9999., -9999., &  !  18 Regosols
                  1200.,  1500.,  1600.,  2300., -9999., &  !  19 Solonetz
                   900.,  1000.,  2300., -9999., -9999., &  !  20 Andosols
                  1100.,  1500.,  2300., -9999., -9999., &  !  21 Rankers
                  1100.,  1300.,  1500.,  2300., -9999., &  !  22 Vertisols
                   600.,  1500.,  2300., -9999., -9999., &  !  23 Planosols
                  1100.,  1400.,  1500.,  2300., -9999., &  !  24 Xerosols
                  1100.,  1400.,  2300., -9999., -9999., &  !  25 Yermosols
                  1100.,  1500.,  2300., -9999., -9999., &  !  26 Solonchaks
                  1200.,  2300., -9999., -9999., -9999., &  !  27 Dunes or shifting sands
                  1400.,  2300., -9999., -9999., -9999., &  !  28 Salt flats
                  1400.,  2300., -9999., -9999., -9999., &  !  29 Rock debris or desert detritus
                   915., -9999., -9999., -9999., -9999., &  !  30 Glaciers
                  1000., -9999., -9999., -9999., -9999./    !  31 Water body

! Dry soil heat capacity (J/kg/K)
!                                Horizon
!                   1       2       3       4       5          Soil type
data cgdry/        900.,   900.,   900.,   800.,   800., &  !   1 Acrisols
                   900.,   850.,   900.,   900.,   800., &  !   2 Cambisols
                   900.,   900.,   850.,   850.,   800., &  !   3 Chernozems
                   900.,   850.,   800.,   800.,   800., &  !   4 Podzoluvisols
                   900.,   850.,   800., -9999., -9999., &  !   5 Rendzinans
                   850.,   900.,   900.,   800., -9999., &  !   6 Ferrasols
                   850.,   750.,   700.,   800., -9999., &  !   7 Gleysols
                   850.,   850.,   700.,   800., -9999., &  !   8 Phaeozems
                   850.,   800., -9999., -9999., -9999., &  !   9 Lithosols
                   850.,   800.,   800., -9999., -9999., &  !  10 Fluvisols
                   850.,   850.,   800.,   800., -9999., &  !  11 Kastanozems
                   900.,   850.,   800.,   800.,   800., &  !  12 Luvisols
                   850.,   850.,   700.,   800., -9999., &  !  13 Greyzems
                   900.,   800.,   800.,   800., -9999., &  !  14 Nitosols
                  1000.,   900.,   900.,   900., -9999., &  !  15 Histosols
                   900.,   850.,   800.,   800.,   800., &  !  16 Podzols
                   850.,   800.,   800., -9999., -9999., &  !  17 Arenosols
                   850.,   800.,   800., -9999., -9999., &  !  18 Regosols
                   850.,   850.,   850.,   800., -9999., &  !  19 Solonetz
                   800.,   800.,   800., -9999., -9999., &  !  20 Andosols
                   850.,   800.,   800., -9999., -9999., &  !  21 Rankers
                   850.,   800.,   800.,   800., -9999., &  !  22 Vertisols
                   900.,   850.,   800., -9999., -9999., &  !  23 Planosols
                   850.,   850.,   800.,   800., -9999., &  !  24 Xerosols
                   800.,   800.,   800., -9999., -9999., &  !  25 Yermosols
                   800.,   800.,   800., -9999., -9999., &  !  26 Solonchaks
                   800.,   800., -9999., -9999., -9999., &  !  27 Dunes or shifting sands
                   800.,   800., -9999., -9999., -9999., &  !  28 Salt flats
                   850.,   800., -9999., -9999., -9999., &  !  29 Rock debris or desert detritus
                 2093.4, -9999., -9999., -9999., -9999., &  !  30 Glaciers
                 4186.8, -9999., -9999., -9999., -9999./    !  31 Water body

! Saturated moisture potential of the soil (at QG=QGMAX) (m)
!                                Horizon
!                   1       2       3       4       5          Soil type
data psig0/      -0.350, -0.400, -0.450, -0.500, -0.300, &  !   1 Acrisols
                 -0.300, -0.400, -0.400, -0.450, -0.300, &  !   2 Cambisols
                 -0.300, -0.250, -0.200, -0.150, -0.200, &  !   3 Chernozems
                 -0.330, -0.130, -0.500, -0.350, -0.100, &  !   4 Podzoluvisols
                 -0.200, -0.150, -0.100, -9999., -9999., &  !   5 Rendzinans
                 -0.300, -0.500, -0.500, -0.300, -9999., &  !   6 Ferrasols
                 -0.200, -0.150, -0.400, -0.300, -9999., &  !   7 Gleysols
                 -0.500, -0.300, -0.150, -0.100, -9999., &  !   8 Phaeozems
                 -0.150, -0.100, -9999., -9999., -9999., &  !   9 Lithosols
                 -0.500, -0.500, -0.300, -9999., -9999., &  !  10 Fluvisols
                 -0.200, -0.150, -0.150, -0.200, -9999., &  !  11 Kastanozems
                 -0.300, -0.300, -0.650, -0.450, -0.300, &  !  12 Luvisols
                 -0.500, -0.300, -0.150, -0.100, -9999., &  !  13 Greyzems
                 -0.300, -0.600, -0.600, -0.300, -9999., &  !  14 Nitosols
                 -0.100, -0.600, -0.600, -0.300, -9999., &  !  15 Histosols
                 -0.350, -0.130, -0.600, -0.350, -0.100, &  !  16 Podzols
                 -0.100, -0.150, -0.300, -9999., -9999., &  !  17 Arenosols
                 -0.150, -0.200, -0.300, -9999., -9999., &  !  18 Regosols
                 -0.250, -0.350, -0.450, -0.300, -9999., &  !  19 Solonetz
                 -0.050, -0.020, -0.050, -9999., -9999., &  !  20 Andosols
                 -0.200, -0.150, -0.100, -9999., -9999., &  !  21 Rankers
                 -0.400, -0.500, -0.500, -0.300, -9999., &  !  22 Vertisols
                 -0.200, -0.600, -0.300, -9999., -9999., &  !  23 Planosols
                 -0.300, -0.400, -0.350, -0.300, -9999., &  !  24 Xerosols
                 -0.100, -0.080, -0.100, -9999., -9999., &  !  25 Yermosols
                 -0.300, -0.350, -0.300, -9999., -9999., &  !  26 Solonchaks
                 -0.100, -0.300, -9999., -9999., -9999., &  !  27 Dunes or shifting sands
                 -0.200, -0.300, -9999., -9999., -9999., &  !  28 Salt flats
                 -0.100, -0.300, -9999., -9999., -9999., &  !  29 Rock debris or desert detritus
                   2.20, -9999., -9999., -9999., -9999., &  !  30 Glaciers
                   0.57, -9999., -9999., -9999., -9999./    !  31 Water body

! Saturated hydraulic conductivity (m/s)
!                                Horizon
!                   1       2        3       4        5          Soil type
data ykghy/      2.0e-7, 2.5e-7,  3.0e-7, 4.5e-7, 0.5e-7, &  !   1 Acrisols
                 5.0e-7, 7.5e-7,  8.5e-7,10.0e-7, 1.0e-7, &  !   2 Cambisols
                 5.0e-7, 9.0e-7, 10.0e-7, 2.0e-7, 0.5e-7, &  !   3 Chernozems
                 3.0e-7, 8.0e-7, 12.5e-7, 1.5e-7, 0.5e-7, &  !   4 Podzoluvisols
                 4.5e-7, 5.0e-7,  1.0e-7, -9999., -9999., &  !   5 Rendzinans
                 1.5e-7, 6.5e-7,  1.5e-7, 1.0e-7, -9999., &  !   6 Ferrasols
                 4.5e-7,11.0e-7,  3.5e-7, 1.0e-7, -9999., &  !   7 Gleysols
                 7.0e-7,10.0e-7,  5.0e-7, 1.0e-7, -9999., &  !   8 Phaeozems
                12.0e-7, 1.0e-7,  -9999., -9999., -9999., &  !   9 Lithosols
                 1.2e-7, 1.5e-7,  1.0e-7, -9999., -9999., &  !  10 Fluvisols
                 8.0e-7,11.0e-7,  4.0e-7, 0.5e-7, -9999., &  !  11 Kastanozems
                 5.0e-7, 6.0e-7,  8.0e-7,10.0e-7, 1.0e-7, &  !  12 Luvisols
                 7.0e-7,10.0e-7,  5.0e-7, 1.0e-7, -9999., &  !  13 Greyzems
                 1.0e-7, 2.0e-7,  2.5e-7, 0.5e-7, -9999., &  !  14 Nitosols
                 6.0e-7, 2.0e-7,  2.5e-7, 0.5e-7, -9999., &  !  15 Histosols
                 3.0e-7, 8.0e-7, 10.0e-7, 1.5e-7, 0.5e-7, &  !  16 Podzols
                12.0e-7, 3.5e-7,  1.0e-7, -9999., -9999., &  !  17 Arenosols
                 6.0e-7, 11.e-7,  1.0e-7, -9999., -9999., &  !  18 Regosols
                 2.0e-7, 4.5e-7,  6.0e-7, 1.0e-7, -9999., &  !  19 Solonetz
                12.0e-7,20.0e-7,  2.0e-7, -9999., -9999., &  !  20 Andosols
                 8.0e-7,10.0e-7,  1.0e-7, -9999., -9999., &  !  21 Rankers
                 3.0e-7, 6.0e-7, 10.0e-7, 1.0e-7, -9999., &  !  22 Vertisols
                 2.0e-7, 3.0e-7,  1.0e-7, -9999., -9999., &  !  23 Planosols
                 8.0e-7,12.5e-7,  1.5e-7, 0.5e-7, -9999., &  !  24 Xerosols
                 5.5e-7,11.0e-7,  0.5e-7, -9999., -9999., &  !  25 Yermosols
                 3.0e-7, 4.0e-7,  1.0e-7, -9999., -9999., &  !  26 Solonchaks
                15.0e-7, 1.0e-7,  -9999., -9999., -9999., &  !  27 Dunes or shifting sands
                 1.0e-7, 1.0e-7,  -9999., -9999., -9999., &  !  28 Salt flats
                12.0e-7, 1.0e-7,  -9999., -9999., -9999., &  !  29 Rock debris or desert detritus
                     0., -9999.,  -9999., -9999., -9999., &  !  30 Glaciers
                     0., -9999.,  -9999., -9999., -9999./    !  31 Water body

! Exponent B (dimensionless) for moisture potential equation:
! PSI=PSIG*(QGMAX/QG)**B
!                                Horizon
!                   1       2       3       4       5          Soil type
data bghyd/        8.00,   9.50,   8.00,  11.00,   3.00, &  !   1 Acrisols
                   5.00,   4.50,   7.50,   9.00,   3.00, &  !   2 Cambisols
                   5.00,   5.00,   4.50,   4.00,   3.00, &  !   3 Chernozems
                   5.00,   4.50,   8.50,   7.50,   3.00, &  !   4 Podzoluvisols
                   5.00,   4.10,   3.00, -9999., -9999., &  !   5 Rendzinans
                   5.00,   9.00,  10.00,   5.00, -9999., &  !   6 Ferrasols
                   4.90,   4.20,  11.00,   5.00, -9999., &  !   7 Gleysols
                   5.30,   4.90,   4.10,   3.00, -9999., &  !   8 Phaeozems
                   4.50,   3.00, -9999., -9999., -9999., &  !   9 Lithosols
                   5.00,  10.50,   3.00, -9999., -9999., &  !  10 Fluvisols
                   4.70,   4.50,   4.00,   3.00, -9999., &  !  11 Kastanozems
                   5.00,   6.00,   8.00,   4.00,   3.00, &  !  12 Luvisols
                   5.30,   4.90,   4.10,   3.00, -9999., &  !  13 Greyzems
                   8.50,  11.50,  11.00,   5.00, -9999., &  !  14 Nitosols
                   4.00,  11.50,  11.00,   5.00, -9999., &  !  15 Histosols
                   5.00,   4.50,  10.50,   7.50,    3.0, &  !  16 Podzols
                   4.20,  10.00,   3.00, -9999., -9999., &  !  17 Arenosols
                   6.00,   7.00,   3.00, -9999., -9999., &  !  18 Regosols
                   7.00,  10.00,  11.00,     3., -9999., &  !  19 Solonetz
                   4.10,   4.00,   3.00, -9999., -9999., &  !  20 Andosols
                   4.50,   4.00,   3.00, -9999., -9999., &  !  21 Rankers
                   7.50,   8.00,   8.00,   3.00, -9999., &  !  22 Vertisols
                   5.50,  11.00,   3.00, -9999., -9999., &  !  23 Planosols
                   8.00,   9.00,   8.50,   3.00, -9999., &  !  24 Xerosols
                   4.00,   4.50,   3.00, -9999., -9999., &  !  25 Yermosols
                   8.00,  10.00,   3.00, -9999., -9999., &  !  26 Solonchaks
                   4.00,   3.00, -9999., -9999., -9999., &  !  27 Dunes or shifting sands
                   6.00,   3.00, -9999., -9999., -9999., &  !  28 Salt flats
                   4.50,   3.00, -9999., -9999., -9999., &  !  29 Rock debris or desert detritus
                   0.00, -9999., -9999., -9999., -9999., &  !  30 Glaciers
                   0.00, -9999., -9999., -9999., -9999./    !  31 Water body

! Exponent C (dimensionless) for ... equation:
! 
!                                Horizon
!                   1       2       3       4       5          Soil type
data cghyd/        3.00,   3.00,   3.00,   3.00,   3.00, &  !   1 Acrisols
                   3.00,   3.00,   3.00,   3.00,   3.00, &  !   2 Cambisols
                   3.00,   3.00,   3.00,   3.00,   3.00, &  !   3 Chernozems
                   3.00,   3.00,   3.00,   3.00,   3.00, &  !   4 Podzoluvisols
                   3.00,   3.00,   3.00, -9999., -9999., &  !   5 Rendzinans
                   3.00,   3.00,   3.00,   3.00, -9999., &  !   6 Ferrasols
                   3.00,   3.00,   3.00,   3.00, -9999., &  !   7 Gleysols
                   3.00,   3.00,   3.00,   3.00, -9999., &  !   8 Phaeozems
                   3.00,   3.00, -9999., -9999., -9999., &  !   9 Lithosols
                   3.00,   3.00,   3.00, -9999., -9999., &  !  10 Fluvisols
                   3.00,   3.00,   3.00,   3.00, -9999., &  !  11 Kastanozems
                   3.00,   3.00,   3.00,   3.00,   3.00, &  !  12 Luvisols
                   3.00,   3.00,   3.00,   3.00, -9999., &  !  13 Greyzems
                   3.00,   3.00,   3.00,   3.00, -9999., &  !  14 Nitosols
                   3.00,   3.00,   3.00,   3.00, -9999., &  !  15 Histosols
                   3.00,   3.00,   3.00,   3.00,   3.00, &  !  16 Podzols
                   3.00,   3.00,   3.00, -9999., -9999., &  !  17 Arenosols
                   3.00,   3.00,   3.00, -9999., -9999., &  !  18 Regosols
                   3.00,   3.00,   3.00,   3.00, -9999., &  !  19 Solonetz
                   3.00,   3.00,   3.00, -9999., -9999., &  !  20 Andosols
                   3.00,   3.00,   3.00, -9999., -9999., &  !  21 Rankers
                   3.00,   3.00,   3.00,   3.00, -9999., &  !  22 Vertisols
                   3.00,   3.00,   3.00, -9999., -9999., &  !  23 Planosols
                   3.00,   3.00,   3.00,   3.00, -9999., &  !  24 Xerosols
                   3.00,   3.00,   3.00, -9999., -9999., &  !  25 Yermosols
                   3.00,   3.00,   3.00, -9999., -9999., &  !  26 Solonchaks
                   3.00,   3.00, -9999., -9999., -9999., &  !  27 Dunes or shifting sands
                   3.00,   3.00, -9999., -9999., -9999., &  !  28 Salt flats
                   3.00,   3.00, -9999., -9999., -9999., &  !  29 Rock debris or desert detritus
                   0.00, -9999., -9999., -9999., -9999., &  !  30 Glaciers
                   0.00, -9999., -9999., -9999., -9999./    !  31 Water body

! Relative (0-1) water content of the soil at wilting point (proportion)
!                                Horizon
!                   1       2       3       4       5          Soil type
data &
 qg_rel_wilt/      0.45,   0.42,   0.35,   0.40,   0.10, &  !   1 Acrisols
                   0.40,   0.30,   0.35,   0.40,   0.10, &  !   2 Cambisols
                   0.40,   0.20,   0.15,   0.15,   0.10, &  !   3 Chernozems
                   0.30,   0.17,   0.30,   0.25,   0.10, &  !   4 Podzoluvisols
                   0.20,   0.10,   0.10, -9999., -9999., &  !   5 Rendzinans
                   0.25,   0.40,   0.50,   0.20, -9999., &  !   6 Ferrasols
                   0.18,   0.15,   0.40,   0.20, -9999., &  !   7 Gleysols
                   0.25,   0.20,   0.10,   0.10, -9999., &  !   8 Phaeozems
                   0.15,   0.10, -9999., -9999., -9999., &  !   9 Lithosols
                   0.22,   0.40,   0.10, -9999., -9999., &  !  10 Fluvisols
                   0.17,   0.30,   0.25,   0.10, -9999., &  !  11 Kastanozems
                   0.40,   0.20,   0.35,   0.30,   0.10, &  !  12 Luvisols
                   0.25,   0.20,   0.10,   0.10, -9999., &  !  13 Greyzems
                   0.35,   0.40,   0.40,   0.20, -9999., &  !  14 Nitosols
                   0.50,   0.40,   0.40,   0.20, -9999., &  !  15 Histosols
                   0.30,   0.17,   0.40,   0.25,   0.10, &  !  16 Podzols
                   0.15,   0.30,   0.10, -9999., -9999., &  !  17 Arenosols
                   0.20,   0.30,   0.10, -9999., -9999., &  !  18 Regosols
                   0.25,   0.30,   0.40,   0.10, -9999., &  !  19 Solonetz
                   0.10,   0.10,   0.10, -9999., -9999., &  !  20 Andosols
                   0.20,   0.15,   0.10, -9999., -9999., &  !  21 Rankers
                   0.35,   0.40,   0.40,   0.10, -9999., &  !  22 Vertisols
                   0.20,   0.40,   0.10, -9999., -9999., &  !  23 Planosols
                   0.30,   0.40,   0.35,   0.10, -9999., &  !  24 Xerosols
                   0.10,   0.15,   0.10, -9999., -9999., &  !  25 Yermosols
                   0.30,   0.40,   0.10, -9999., -9999., &  !  26 Solonchaks
                   0.10,   0.10, -9999., -9999., -9999., &  !  27 Dunes or shifting sands
                   0.60,   0.10, -9999., -9999., -9999., &  !  28 Salt flats
                   0.15,   0.10, -9999., -9999., -9999., &  !  29 Rock debris or desert detritus
                   0.00, -9999., -9999., -9999., -9999., &  !  30 Glaciers
                   0.00, -9999., -9999., -9999., -9999./    !  31 Water body

! Relative (0-1) water content of soil at reference point for evapotranspiration (proportion)
!                                Horizon
!                   1       2       3       4       5          Soil type
data &
 qg_rel_refw/      0.70,   0.70,   0.75,   0.68,   0.65, &  !   1 Acrisols
                   0.70,   0.70,   0.75,   0.68,   0.65, &  !   2 Cambisols
                   0.70,   0.80,   0.65,   0.68,   0.65, &  !   3 Chernozems
                   0.68,   0.67,   0.65,   0.68,   0.65, &  !   4 Podzoluvisols
                   0.70,   0.67,   0.65, -9999., -9999., &  !   5 Rendzinans
                   0.80,   0.70,   0.70,   0.65, -9999., &  !   6 Ferrasols
                   0.70,   0.68,   0.67,   0.65, -9999., &  !   7 Gleysols
                   0.70,   0.68,   0.67,   0.65, -9999., &  !   8 Phaeozems
                   0.67,   0.65, -9999., -9999., -9999., &  !   9 Lithosols
                   0.80,   0.68,   0.65, -9999., -9999., &  !  10 Fluvisols
                   0.70,   0.65,   0.68,   0.65, -9999., &  !  11 Kastanozems
                   0.70,   0.70,   0.70,   0.68,   0.65, &  !  12 Luvisols
                   0.70,   0.68,   0.67,   0.65, -9999., &  !  13 Greyzems
                   0.75,   0.68,   0.67,   0.65, -9999., &  !  14 Nitosols
                   0.65,   0.68,   0.67,   0.65, -9999., &  !  15 Histosols
                   0.68,   0.67,   0.65,   0.68,   0.65, &  !  16 Podzols
                   0.67,   0.70,   0.65, -9999., -9999., &  !  17 Arenosols
                   0.67,   0.66,   0.65, -9999., -9999., &  !  18 Regosols
                   0.68,   0.67,   0.66,   0.65, -9999., &  !  19 Solonetz
                   0.66,   0.65,   0.65, -9999., -9999., &  !  20 Andosols
                   0.68,   0.67,   0.65, -9999., -9999., &  !  21 Rankers
                   0.72,   0.70,   0.68,   0.65, -9999., &  !  22 Vertisols
                   0.68,   0.70,   0.65, -9999., -9999., &  !  23 Planosols
                   0.75,   0.68,   0.67,   0.65, -9999., &  !  24 Xerosols
                   0.67,   0.66,   0.65, -9999., -9999., &  !  25 Yermosols
                   0.68,   0.67,   0.65, -9999., -9999., &  !  26 Solonchaks
                   0.67,   0.65, -9999., -9999., -9999., &  !  27 Dunes or shifting sands
                   0.65,   0.65, -9999., -9999., -9999., &  !  28 Salt flats
                   0.67,   0.65, -9999., -9999., -9999., &  !  29 Rock debris or desert detritus
                   0.00, -9999., -9999., -9999., -9999., &  !  30 Glaciers
                   0.00, -9999., -9999., -9999., -9999./    !  31 Water body

! Dry soil (qg=qgmin) albedo
!                                Soil type
!                 Acrisols    Cambisols   Chernozems   Podzoluvisols Rendzinans
data soildalbed/    0.20,        0.20,       0.10,         0.20,       0.10, & 
!                Ferrasols    Gleysols    Phaeozems    Lithosols     Fluvisols
!!!                    0.15,        0.20,       0.10,         0.22,       0.15, & 
                    0.15,        0.20,       0.10,         0.15,       0.15, & 
!                Kastanozems  Luvisols    Greyzems     Nitosols      Histosols 
                    0.15,        0.20,       0.10,         0.18,       0.15, & 
!                Podzols      Arenosols   Regosols     Solonetz      Andosols
                    0.20,        0.25,       0.25,         0.17,       0.05, & 
!                Rankers      Vertisols   Planosols    Xerosols      Yermosols
!!!                    0.20,        0.17,       0.20,         0.24,       0.22, & 
                    0.20,        0.17,       0.20,         0.17,       0.16, & 
!                Solonchaks   Dunes       Salt flats   Rock debris   Glaciers
!!!                    0.22,        0.30,       0.35,         0.25,       0.30, &
                    0.20,        0.20,       0.25,         0.20,       0.30, &
!                Water body
                    0.07/

! Wet soil (qg=qgmax) albedo
!                                Soil type
!                 Acrisols    Cambisols   Chernozems   Podzoluvisols Rendzinans
data soilwalbed/    0.10,        0.10,       0.05,         0.10,       0.05, & 
!                Ferrasols    Gleysols    Phaeozems    Lithosols     Fluvisols
!!!                    0.08,        0.13,       0.05,         0.14,       0.08, & 
                    0.08,        0.13,       0.05,         0.08,       0.08, & 
!                Kastanozems  Luvisols    Greyzems     Nitosols      Histosols 
                    0.08,        0.10,       0.05,         0.10,       0.05, & 
!                Podzols      Arenosols   Regosols     Solonetz      Andosols
                    0.10,        0.15,       0.15,         0.14,       0.03, & 
!                Rankers      Vertisols   Planosols    Xerosols      Yermosols
!!!                    0.10,        0.13,       0.10,         0.14,       0.12, & 
                    0.10,        0.13,       0.10,         0.08,       0.07, & 
!                Solonchaks   Dunes       Salt flats   Rock debris   Glaciers
!!!                    0.18,        0.18,       0.20,         0.18,       0.30, &
                    0.18,        0.13,       0.20,         0.18,       0.30, &
!                Water body
                    0.07/

! Dry soil (qg=qgmin) emissivity in broadband window
!                                Soil type
!                 Acrisols    Cambisols   Chernozems   Podzoluvisols Rendzinans
data soildemis1/   0.960,       0.960,      0.975,        0.960,      0.970, & 
!                Ferrasols    Gleysols    Phaeozems    Lithosols     Fluvisols
!!!                   0.965,       0.960,      0.975,        0.940,      0.965, & 
                   0.965,       0.960,      0.975,        0.960,      0.965, & 
!                Kastanozems  Luvisols    Greyzems     Nitosols      Histosols 
                   0.970,       0.960,      0.975,        0.970,      0.970, & 
!                Podzols      Arenosols   Regosols     Solonetz      Andosols
                   0.960,       0.940,      0.940,        0.970,      0.980, & 
!                Rankers      Vertisols   Planosols    Xerosols      Yermosols
!!!                   0.960,       0.960,      0.960,        0.950,      0.960, & 
                   0.960,       0.960,      0.960,        0.960,      0.965, & 
!                Solonchaks   Dunes       Salt flats   Rock debris   Glaciers
!!!                   0.950,       0.930,      0.940,        0.950,      0.970, &
                   0.950,       0.945,      0.945,        0.960,      0.970, &
!                Water body
                   0.970/

! Wet soil (qg=qgmax) emissivity in broadband window
!                                Soil type
!                 Acrisols    Cambisols   Chernozems   Podzoluvisols Rendzinans
data soilwemis1/   0.985,       0.985,      0.995,        0.985,      0.990, & 
!                Ferrasols    Gleysols    Phaeozems    Lithosols     Fluvisols
!!!                   0.985,       0.985,      0.995,        0.980,      0.985, & 
                   0.985,       0.985,      0.995,        0.985,      0.985, & 
!                Kastanozems  Luvisols    Greyzems     Nitosols      Histosols 
                   0.990,       0.985,      0.995,        0.990,      0.990, & 
!                Podzols      Arenosols   Regosols     Solonetz      Andosols
                   0.985,       0.970,      0.970,        0.990,      0.995, & 
!                Rankers      Vertisols   Planosols    Xerosols      Yermosols
!!!                   0.985,       0.980,      0.985,        0.975,      0.980, & 
                   0.985,       0.980,      0.985,        0.980,      0.985, & 
!                Solonchaks   Dunes       Salt flats   Rock debris   Glaciers
!!!                   0.975,       0.960,      0.980,        0.960,      0.970, &
                   0.980,       0.975,      0.975,        0.970,      0.970, &
!                Water body
                   0.970/

! Dry soil (qg=qgmin) emissivity in 8-12 micron window
!                                Soil type
!                 Acrisols    Cambisols   Chernozems   Podzoluvisols Rendzinans
data soildemis2/   0.950,       0.950,      0.960,        0.950,      0.955, & 
!                Ferrasols    Gleysols    Phaeozems    Lithosols     Fluvisols
!!!                   0.950,       0.950,      0.960,        0.935,      0.950, & 
                   0.950,       0.950,      0.960,        0.945,      0.950, & 
!                Kastanozems  Luvisols    Greyzems     Nitosols      Histosols 
                   0.955,       0.950,      0.960,        0.955,      0.960, & 
!                Podzols      Arenosols   Regosols     Solonetz      Andosols
                   0.950,       0.920,      0.940,        0.960,      0.970, & 
!                Rankers      Vertisols   Planosols    Xerosols      Yermosols
!!!                   0.950,       0.940,      0.950,        0.930,      0.940, & 
                   0.950,       0.940,      0.950,        0.945,      0.950, & 
!                Solonchaks   Dunes       Salt flats   Rock debris   Glaciers
!!!                   0.940,       0.890,      0.920,        0.920,      0.980, &
                   0.940,       0.930,      0.925,        0.930,      0.980, &
!                Water body
                   0.970/

! Wet soil (qg=qgmax) emissivity in 8-12 micron window
!                                Soil type
!                 Acrisols    Cambisols   Chernozems   Podzoluvisols Rendzinans
data soilwemis2/   0.970,       0.970,      0.985,        0.955,      0.990, & 
!                Ferrasols    Gleysols    Phaeozems    Lithosols     Fluvisols
!!!                   0.970,       0.955,      0.985,        0.940,      0.970, & 
                   0.970,       0.955,      0.985,        0.965,      0.970, & 
!                Kastanozems  Luvisols    Greyzems     Nitosols      Histosols 
                   0.980,       0.970,      0.985,        0.975,      0.985, & 
!                Podzols      Arenosols   Regosols     Solonetz      Andosols
                   0.970,       0.940,      0.960,        0.980,      0.990, & 
!                Rankers      Vertisols   Planosols    Xerosols      Yermosols
!!!                   0.970,       0.950,      0.970,        0.950,      0.960, & 
                   0.970,       0.950,      0.970,        0.960,      0.970, & 
!                Solonchaks   Dunes       Salt flats   Rock debris   Glaciers
!!!                   0.960,       0.920,      0.950,        0.930,      0.980, &
                   0.960,       0.950,      0.955,        0.960,      0.980, &
!                Water body
                   0.970/

! Approximate depth (m) of water table top
!                                Soil type
!                 Acrisols    Cambisols   Chernozems   Podzoluvisols Rendzinans
data &
!!! water_table_depth/  10.,          4.,        20.,           2.,        50., & 
 water_table_depth/  10.,         20.,        25.,           2.,        50., & 
!                Ferrasols    Gleysols    Phaeozems    Lithosols     Fluvisols
!!!                     20.,          2.,        20.,         100.,         2., & 
                     20.,          5.,        20.,         100.,         5., & 
!                Kastanozems  Luvisols    Greyzems     Nitosols      Histosols 
!!!                     20.,          4.,        10.,          15.,        0.5, & 
                     20.,         20.,        10.,          15.,        0.5, & 
!                Podzols      Arenosols   Regosols     Solonetz      Andosols
                      3.,         50.,       100.,          50.,       100., & 
!                Rankers      Vertisols   Planosols    Xerosols      Yermosols
                     10.,         10.,        1.5,          50.,        50., & 
!                Solonchaks   Dunes       Salt flats   Rock debris   Glaciers
                     50.,        100.,        50.,          50.,         0., &
!                Water body
                      0./

 print *
 print *,"Definition of physical parameters of soil (FAO soil classification), parameters change depending on soil level"

 nsmooth = 0
! nsmooth = 1
 if (nsmooth > 0) print*, "nsmooth in presoil", nsmooth

! Definition of physical parameters of soil

! Definition of soil parameters at soil levels for various soil types using data
! about soil horizons location

! Correction of hydraulic conductivity following the results of experiment:
! real atmospheric condition in arid zone (Greece), no (!) presipitation
! during 3 summer month, it is necessary reduce hydraulic conductivity in 10 (!)
! times to obtain the decreasing of top soil water content

 do itype=1,nst
   do ihorizon=1,soil_horizon_number
     if (int(ykghy(ihorizon,itype)) /= -9999) ykghy(ihorizon,itype)=ykghy(ihorizon,itype)*0.1
   enddo
 enddo

 do itype=1,nst
   do ilev=1,nlevel_soil
     do ihorizon=1,soil_horizon_number
       if (level_soil(ilev)-horizon_bottom(ihorizon,itype) <= 0.) exit
     enddo
     qgmax_lev(ilev,itype) = qgmax(ihorizon,itype)
     qgmin_lev(ilev,itype) = qgmin(ihorizon,itype)
     rogdr_lev(ilev,itype) = rogdr(ihorizon,itype)
     cgdry_lev(ilev,itype) = cgdry(ihorizon,itype)
     psig0_lev(ilev,itype) = psig0(ihorizon,itype)
     ykghy_lev(ilev,itype) = ykghy(ihorizon,itype)
     bghyd_lev(ilev,itype) = bghyd(ihorizon,itype)
     cghyd_lev(ilev,itype) = cghyd(ihorizon,itype)
     qg_rel_wilt_lev(ilev,itype) = qg_rel_wilt(ihorizon,itype)
     qg_rel_refw_lev(ilev,itype) = qg_rel_refw(ihorizon,itype)
   enddo
 enddo

 soil_qmax(:,:,:) = 0.
 soil_qmin(:,:,:) = 0.
 soil_c(:,:,:) = 0.
 soil_rho(:,:,:) = 0.
 soil_psi(:,:,:) = 0.
 soil_k(:,:,:) = 0.
 soil_par_b(:,:,:) = 0.
 soil_par_c(:,:,:) = 0.
 soil_qrel_wilt(:,:,:) = 0.
 soil_qrel_ref(:,:,:) = 0.
 soil_water_table(:,:) = 0.
 soil_albedo_dry(:,:) = 0.
 soil_albedo_wet(:,:) = 0.
 soil_emiss_1_dry(:,:) = 0.
 soil_emiss_1_wet(:,:) = 0.
 soil_emiss_2_dry(:,:) = 0.
 soil_emiss_2_wet(:,:) = 0.

 do jp=1,nlat
 do ip=1,nlon

   do itype=1,nst
     n11 = itype
     zfrac = soil_fao(ip,jp,n11+1)

! Weighted (with soil fraction) average of the various soil parameters.
! The following fills glacier and water points by soil type "Podzols", for subsequent spatial smoothing only

      if (n11 >= nst-1.and.nsmooth /= 0) n11 = 16

      do ilev=1,nlevel_soil
        soil_qmax(ip,jp,ilev) = soil_qmax(ip,jp,ilev) + qgmax_lev(ilev,n11)*zfrac
        soil_qmin(ip,jp,ilev) = soil_qmin(ip,jp,ilev) + qgmin_lev(ilev,n11)*zfrac
        soil_c(ip,jp,ilev) = soil_c(ip,jp,ilev) + cgdry_lev(ilev,n11)*zfrac
        soil_rho(ip,jp,ilev) = soil_rho(ip,jp,ilev) + rogdr_lev(ilev,n11)*zfrac
        soil_psi(ip,jp,ilev) = soil_psi(ip,jp,ilev) + psig0_lev(ilev,n11)*zfrac
        soil_k(ip,jp,ilev) = soil_k(ip,jp,ilev) + ykghy_lev(ilev,n11)*zfrac
        soil_par_b(ip,jp,ilev) = soil_par_b(ip,jp,ilev) + bghyd_lev(ilev,n11)*zfrac
        soil_par_c(ip,jp,ilev) = soil_par_c(ip,jp,ilev) + cghyd_lev(ilev,n11)*zfrac
        soil_qrel_wilt(ip,jp,ilev) = soil_qrel_wilt(ip,jp,ilev) + qg_rel_wilt_lev(ilev,n11)*zfrac
        soil_qrel_ref(ip,jp,ilev) = soil_qrel_ref(ip,jp,ilev) + qg_rel_refw_lev(ilev,n11)*zfrac
      enddo

      soil_water_table(ip,jp) = soil_water_table(ip,jp) + water_table_depth(n11)*zfrac
      soil_albedo_dry(ip,jp) = soil_albedo_dry(ip,jp) + soildalbed(n11)*zfrac
      soil_albedo_wet(ip,jp) = soil_albedo_wet(ip,jp) + soilwalbed(n11)*zfrac
      soil_emiss_1_dry(ip,jp) = soil_emiss_1_dry(ip,jp) + soildemis1(n11)*zfrac
      soil_emiss_1_wet(ip,jp) = soil_emiss_1_wet(ip,jp) + soilwemis1(n11)*zfrac
      soil_emiss_2_dry(ip,jp) = soil_emiss_2_dry(ip,jp) + soildemis2(n11)*zfrac
      soil_emiss_2_wet(ip,jp) = soil_emiss_2_wet(ip,jp) + soilwemis2(n11)*zfrac

   enddo

 enddo
 enddo

! Horizontal smoothing of soil parameters at all soil levels

 if (nsmooth.ne.0) then

   do ilev=1,nlevel_soil
     zzz1(:,:) = soil_qmax(:,:,ilev)
     call smooth_soil(zzz1,zzz2,nlon,nlat,0.2,nsmooth)
     soil_qmax(:,:,ilev) = zzz2(:,:)
   enddo

   do ilev=1,nlevel_soil
     zzz1(:,:) = soil_qmin(:,:,ilev)
     call smooth_soil(zzz1,zzz2,nlon,nlat,0.2,nsmooth)
     soil_qmin(:,:,ilev) = zzz2(:,:)
   enddo

   do ilev=1,nlevel_soil
     zzz1(:,:) = soil_c(:,:,ilev)
     call smooth_soil(zzz1,zzz2,nlon,nlat,0.2,nsmooth)
     soil_c(:,:,ilev) = zzz2(:,:)
   enddo

   do ilev=1,nlevel_soil
     zzz1(:,:) = soil_rho(:,:,ilev)
     call smooth_soil(zzz1,zzz2,nlon,nlat,0.2,nsmooth)
     soil_rho(:,:,ilev) = zzz2(:,:)
   enddo

   do ilev=1,nlevel_soil
     zzz1(:,:) = soil_psi(:,:,ilev)
     call smooth_soil(zzz1,zzz2,nlon,nlat,0.2,nsmooth)
     soil_psi(:,:,ilev) = zzz2(:,:)
   enddo

   do ilev=1,nlevel_soil
     zzz1(:,:) = soil_k(:,:,ilev)
     call smooth_soil(zzz1,zzz2,nlon,nlat,0.2,nsmooth)
     soil_k(:,:,ilev) = zzz2(:,:)
   enddo

   do ilev=1,nlevel_soil
     zzz1(:,:) = soil_par_b(:,:,ilev)
     call smooth_soil(zzz1,zzz2,nlon,nlat,0.2,nsmooth)
     soil_par_b(:,:,ilev) = zzz2(:,:)
   enddo

   do ilev=1,nlevel_soil
     zzz1(:,:) = soil_par_c(:,:,ilev)
     call smooth_soil(zzz1,zzz2,nlon,nlat,0.2,nsmooth)
     soil_par_c(:,:,ilev) = zzz2(:,:)
   enddo

   do ilev=1,nlevel_soil
     zzz1(:,:) = soil_qrel_wilt(:,:,ilev)
     call smooth_soil(zzz1,zzz2,nlon,nlat,0.2,nsmooth)
     soil_qrel_wilt(:,:,ilev) = zzz2(:,:)
   enddo

   do ilev=1,nlevel_soil
     zzz1(:,:) = soil_qrel_ref(:,:,ilev)
     call smooth_soil(zzz1,zzz2,nlon,nlat,0.2,nsmooth)
     soil_qrel_ref(:,:,ilev) = zzz2(:,:)
   enddo

   zzz1(:,:) = soil_water_table(:,:)
   call smooth_soil(zzz1,zzz2,nlon,nlat,0.2,nsmooth)
   soil_water_table(:,:) = zzz2(:,:)

   zzz1(:,:) = soil_albedo_dry(:,:)
   call smooth_soil(zzz1,zzz2,nlon,nlat,0.2,nsmooth)
   soil_albedo_dry(:,:) = zzz2(:,:)

   zzz1(:,:) = soil_albedo_wet(:,:)
   call smooth_soil(zzz1,zzz2,nlon,nlat,0.2,nsmooth)
   soil_albedo_wet(:,:) = zzz2(:,:)

   zzz1(:,:) = soil_emiss_1_dry(:,:)
   call smooth_soil(zzz1,zzz2,nlon,nlat,0.2,nsmooth)
   soil_emiss_1_dry(:,:) = zzz2(:,:)

   zzz1(:,:) = soil_emiss_1_wet(:,:)
   call smooth_soil(zzz1,zzz2,nlon,nlat,0.2,nsmooth)
   soil_emiss_1_wet(:,:) = zzz2(:,:)

   zzz1(:,:) = soil_emiss_2_dry(:,:)
   call smooth_soil(zzz1,zzz2,nlon,nlat,0.2,nsmooth)
   soil_emiss_2_dry(:,:) = zzz2(:,:)

   zzz1(:,:) = soil_emiss_2_wet(:,:)
   call smooth_soil(zzz1,zzz2,nlon,nlat,0.2,nsmooth)
   soil_emiss_2_wet(:,:) = zzz2(:,:)

 endif

! Water body and glacier parameters

 do jp=1,nlat
 do ip=1,nlon

   n11 = soil_fao(ip,jp,1)
   if (n11 >= nst-1) then
     soil_qmax(ip,jp,:) = qgmax(1,n11)
     soil_qmin(ip,jp,:) = qgmin(1,n11)
     soil_c(ip,jp,:) = cgdry(1,n11)
     soil_rho(ip,jp,:) = rogdr(1,n11)
     soil_psi(ip,jp,:) = psig0(1,n11)
     soil_k(ip,jp,:) = ykghy(1,n11)
     soil_par_b(ip,jp,:) = bghyd(1,n11)
     soil_par_c(ip,jp,:) = cghyd(1,n11)
     soil_qrel_wilt(ip,jp,:) = qg_rel_wilt(1,n11)
     soil_qrel_ref(ip,jp,:) = qg_rel_refw(1,n11)
     soil_water_table(ip,jp) = water_table_depth(n11)
     soil_albedo_dry(ip,jp) = soildalbed(n11)
     soil_albedo_wet(ip,jp) = soilwalbed(n11)
     soil_emiss_1_dry(ip,jp) = soildemis1(n11)
     soil_emiss_1_wet(ip,jp) = soilwemis1(n11)
     soil_emiss_2_dry(ip,jp) = soildemis2(n11)
     soil_emiss_2_wet(ip,jp) = soilwemis2(n11)
   else
     soil_qmin(ip,jp,:) = min(soil_qmin(ip,jp,:), soil_qmax(ip,jp,:)-1.e-8)
     soil_qrel_wilt(ip,jp,:) = min(soil_qrel_wilt(ip,jp,:), soil_qrel_ref(ip,jp,:)-1.e-8)
   endif
 enddo
 enddo

return
end subroutine definition_soil_phys_param
!=======================================================================
subroutine definition_veg1_phys_param(vegeta,alat,nlon,nlat,imonthc,idayc, &
 veg_lai, veg_lai_max, veg_frac, veg_root_depth, veg_roughness, veg_albedo, veg_emiss_1, veg_emiss_2)

! Procedure defines physical parameters vegetation using! vegetation type (classification GLC1990)
! defined on the grid points

! Internal basic variables:

!      nvt: total number of vegetation (landuse) types;

! Input variables:

!      vegeta: vegetation (landuse) types (classification GLC1990, 14 types) defined on the grid points, 3D array,
!               1-st and 2-nd array index are grid point index, 3-rd index is vegetation types number plus 1,
!               if 3-rd index has value 1, then the variable means dominate vegetation type,
!               values 2-15 of 3-rd index means fraction (of area, proportion) of the each of 14 vegetation types;    

!      alat: geographical latitude of the grid points (degree);
!      nlon: grid point number along axis of rotated longitude;
!      nlat: grid point number along axis of rotated latidude;
!      imonthc: month number of the date (1-12);
!      idayc: day number of the date.

! Output variables:

!      veg_lai: Leaf Area Index (LAI) of vegetation (dimensionless), 2D array;
!      veg_lai_max: LAI maximum value permited for used vegetation dataset (refer value), real variable;
!      veg_frac: fraction of vegetation (proportion), 2D array;
!      veg_root_depth: root depth of vegetation (m), 2D array;
!      veg_roughness: vegetation cover roughness (m), 2D array;
!      veg_albedo: vegetation albedo (proportion), 2D array;
!      veg_emiss_1: vegetaion emissivity in broadband window (proportion), 2D array;
!      veg_emiss_2: vegetation emissivity in 8-12 micron window (proportion), 2D array.

implicit none

integer, parameter :: nvt=14

! Input:

integer :: nlon, nlat, imonthc, idayc
real, dimension(nlon,nlat,nvt+1) :: vegeta
real, dimension(nlon,nlat) :: alat

! Ouput:

real, dimension(nlon,nlat) :: veg_lai, veg_frac, veg_root_depth, veg_roughness, veg_albedo, veg_emiss_1, veg_emiss_2
real :: veg_lai_max

! Work:

integer, parameter :: ival_missing=-9999
real, dimension(nvt) :: vegrtzn, vegrogh, vegalb, vegemis1, vegemis2, &
 veglain, vegfrcn, veglais, vegfrcs, veglaic, vegfrcc
real, dimension(nlon,nlat) :: zzz1, zzz2
integer, dimension(12) :: imon
integer :: ndayn, itype, ndays, nsmooth, i, ip, jp, nv11
real :: zndayc, zdayref, zalatref, zfrac, z1

! Depth of plant's root zone (m)

data vegrtzn/ &
! 1-Evergreen needleleaf forest 2-Evergreen broadleaf forest 3-Deciduous needleleaf forest
       1.00,                        1.25,                      1.00, &
! 4-Deciduous broadleaf forest  5-Mixed cover                6-Woodland
       1.25,                        0.80,                      0.90, &
! 7-Wooded grassland            8-Closed shrubland           9-Open shrubland
       0.70,                        0.50,                      0.45, &
!10-Grassland                  11-Cropland                  12-Bare soil
       0.35,                        0.40,                      0.00, &
!13-Urban and built-up         14-Water body
       0.60,                        0.00/

! Vegetation cover roughness

data vegrogh/ &
! 1-Evergreen needleleaf forest 2-Evergreen broadleaf forest 3-Deciduous needleleaf forest
       0.80,                        1.00,                      0.80, &
! 4-Deciduous broadleaf forest  5-Mixed cover                6-Woodland
       0.90,                        0.50,                      0.60, &
! 7-Wooded grassland            8-Closed shrubland           9-Open shrubland
       0.40,                        0.18,                      0.14, &
!10-Grassland                  11-Cropland                  12-Bare soil
       0.07,                        0.12,                      0.05, &
!13-Urban and built-up         14-Water body
       1.20,                        0.001/

! Vegetation albedo

data vegalb/ &
! 1-Evergreen needleleaf forest 2-Evergreen broadleaf forest 3-Deciduous needleleaf forest
       0.14,                        0.13,                      0.15, &
! 4-Deciduous broadleaf forest  5-Mixed cover                6-Woodland
       0.15,                        0.13,                      0.15, &
! 7-Wooded grassland            8-Closed shrubland           9-Open shrubland
       0.18,                        0.19,                      0.20, &
!10-Grassland                  11-Cropland                  12-Bare soil
       0.19,                        0.19,                      0.20, &
!13-Urban and built-up         14-Water body
       0.15,                        0.07/

! Vegetation emissivity in broadband window

data vegemis1/ &
! 1-Evergreen needleleaf forest 2-Evergreen broadleaf forest 3-Deciduous needleleaf forest
       0.997,                       0.995,                     0.991, &
! 4-Deciduous broadleaf forest  5-Mixed cover                6-Woodland
       0.992,                       0.990,                     0.987, &
! 7-Wooded grassland            8-Closed shrubland           9-Open shrubland
       0.984,                       0.986,                     0.982, &
!10-Grassland                  11-Cropland                  12-Bare soil
       0.985,                      0.983,                      0.963, &
!13-Urban and built-up         14-Water body
       0.986,                      0.990/

! Vegetation emissivity in 8-12 micron window

data vegemis2/ &
! 1-Evergreen needleleaf forest 2-Evergreen broadleaf forest 3-Deciduous needleleaf forest
       0.993,                       0.985,                     0.983, &
! 4-Deciduous broadleaf forest  5-Mixed cover                6-Woodland
       0.988,                       0.982,                     0.975, &
! 7-Wooded grassland            8-Closed shrubland           9-Open shrubland
       0.975,                       0.975,                     0.970, &
!10-Grassland                  11-Cropland                  12-Bare soil
       0.978,                       0.983,                     0.958, &
!13-Urban and built-up         14-Water body
       0.975,                       0.990/

data imon /31,28,31,30,31,30,31,31,30,31,30,31/

 print *
 print *,"Definition of physical parameters of vegetation (GLC1990 classification)"

! DEFLAY defines vegetation parameters (LAI and veget. fraction) that are
! function of the day of the year

! Date for the Northern Hemisphere season

 ndayn = 0
 do i=1,imonthc-1
   ndayn = ndayn+imon(i)
 enddo
 ndayn = ndayn+idayc

! Date for the Southern Hemisphere season

 ndays = 0
 do i=1,imonthc+6-1
   if (i.le.12) then
     ndays = ndays+imon(i)
   else
     ndays = ndays+imon(i-12)
   endif
 enddo
 ndays = ndays+idayc
 if (ndays.gt.365) ndays = ndays-365

 zndayc = ndayn
 call deflai(veglain,vegfrcn,zndayc)
 zndayc = ndays
 call deflai(veglais,vegfrcs,zndayc)

 zdayref = 183. ! year half, 2 July
 zalatref = 10.

 nsmooth = 0
! nsmooth = 1
 if (nsmooth > 0) print*, "nsmooth in presoil", nsmooth


! Definition of vegetation parameters for various vegeation types

 veg_lai(:,:) = 0.
 veg_frac(:,:) = 0.
 veg_root_depth(:,:) = 0.
 veg_roughness(:,:) = 0.
 veg_albedo(:,:) = 0.
 veg_emiss_1(:,:) = 0.
 veg_emiss_2(:,:) = 0.

 do jp=1,nlat
 do ip=1,nlon

   do itype=1,nvt
     nv11 = itype
     zfrac = vegeta(ip,jp,nv11+1)

! Filling water points by the type "Mixed cover", for spatial smoothing only

     if (nv11.eq.nvt.and.nsmooth.ne.0) nv11 = 5

! Definition of changing (in time) parameter of vegetation
! as a function of date and latitude (hemisphere).
! Season variability is decreased in tropical zone between
! -zalatref and zalatref latitudes. At the equator the date is
! constant and equal to zdayref.

     if (alat(ip,jp).ge.zalatref.or.alat(ip,jp).le.-zalatref) then
       if (alat(ip,jp).ge.0.) then
         veg_lai(ip,jp) = veg_lai(ip,jp)+veglain(nv11)*zfrac
         veg_frac(ip,jp) = veg_frac(ip,jp)+vegfrcn(nv11)*zfrac
       else
         veg_lai(ip,jp) = veg_lai(ip,jp)+veglais(nv11)*zfrac
         veg_frac(ip,jp) = veg_frac(ip,jp)+vegfrcs(nv11)*zfrac
       endif
     else
       if (alat(ip,jp).ge.0.) then
         zndayc = ndayn
         zndayc = zdayref+(zndayc-zdayref)*alat(ip,jp)/zalatref
       else
         zndayc = ndays
         zndayc = zdayref-(zndayc-zdayref)*alat(ip,jp)/zalatref
       endif
       call deflai(veglaic,vegfrcc,zndayc)
       veg_lai(ip,jp) = veg_lai(ip,jp)+veglaic(nv11)*zfrac
       veg_frac(ip,jp) = veg_frac(ip,jp)+vegfrcc(nv11)*zfrac
     endif

     veg_root_depth(ip,jp) = veg_root_depth(ip,jp)+vegrtzn(nv11)*zfrac
     veg_roughness(ip,jp) = veg_roughness(ip,jp)+vegrogh(nv11)*zfrac
     veg_albedo(ip,jp) = veg_albedo(ip,jp)+vegalb(nv11)*zfrac
     veg_emiss_1(ip,jp) = veg_emiss_1(ip,jp)+vegemis1(nv11)*zfrac
     veg_emiss_2(ip,jp) = veg_emiss_2(ip,jp)+vegemis2(nv11)*zfrac

   enddo

 enddo
 enddo

! Horizontal smoothing of vegetation parameters

 if (nsmooth.ne.0) then

   zzz1(:,:) = veg_lai(:,:)
   call smooth_soil(zzz1,zzz2,nlon,nlat,0.2,nsmooth)
   veg_lai(:,:) = zzz2(:,:)

   zzz1(:,:) = veg_frac(:,:)
   call smooth_soil(zzz1,zzz2,nlon,nlat,0.2,nsmooth)
   veg_frac(:,:) = zzz2(:,:)

   zzz1(:,:) = veg_root_depth(:,:)
   call smooth_soil(zzz1,zzz2,nlon,nlat,0.2,nsmooth)
   veg_root_depth(:,:) = zzz2(:,:)

   zzz1(:,:) = veg_roughness(:,:)
   call smooth_soil(zzz1,zzz2,nlon,nlat,0.2,nsmooth)
   veg_roughness(:,:) = zzz2(:,:)

   zzz1(:,:) = veg_emiss_1(:,:)
   call smooth_soil(zzz1,zzz2,nlon,nlat,0.2,nsmooth)
   veg_emiss_1(:,:) = zzz2(:,:)

   zzz1(:,:) = veg_emiss_2(:,:)
   call smooth_soil(zzz1,zzz2,nlon,nlat,0.2,nsmooth)
   veg_emiss_2(:,:) = zzz2(:,:)

 endif

! Water body parameters

 do jp=1,nlat
 do ip=1,nlon

   nv11 = vegeta(ip,jp,1)
   if (nv11 == nvt) then
     veg_lai(ip,jp) = veglain(nvt)
     veg_frac(ip,jp) = vegfrcn(nvt)
     veg_root_depth(ip,jp) = vegrtzn(nvt)
     veg_roughness(ip,jp) = vegrogh(nvt)
     veg_albedo(ip,jp) = vegalb(nvt)
     veg_emiss_1(ip,jp) = vegemis1(nvt)
     veg_emiss_2(ip,jp) = vegemis2(nvt)
   endif
 enddo
 enddo

! LAI maximum value permited for used vegetation dataset (refer value)

 z1=maxval(veg_lai(:,:))
 veg_lai_max = max(z1, 15.)

return
end subroutine definition_veg1_phys_param
!=======================================================================
subroutine definition_veg2_const_phys_param(vegeta,t_climate,t_amplitude,nlon,nlat, &
 veg_water_table, veg_qrel_surf, veg_root_depth, veg_roughness, veg_albedo, veg_emiss_1, veg_emiss_2, snow_dirt)

! Procedure defines physical parameters vegetation using vegetation type (classification GLC2000)
! defined on the grid points

! Internal basic variables:

!      nvt: total number of vegetation (landuse) types;

! Input variables:

!      vegeta: vegetation (landuse) types (classification GLC2000, 22 types) defined on the grid points, 3D array,
!               1-st and 2-nd array index are grid point index, 3-rd index is vegetation types number plus 1,
!               if 3-rd index has value 1, then the variable means dominate vegetation type,
!               values 2-23 of 3-rd index means fraction (of area, proportion) of the each of 22 vegetation types;    

!      t_climate: climate temperature (K);
!      t_amplitude: climate amplitude of temperature (K);
!      nlon: grid point number along axis of rotated longitude;
!      nlat: grid point number along axis of rotated latidude;

! Output variables:

!      veg_water_table: depth of ground water table (m) following to vegetation type;
!      veg_qrel_surf: relative soil water content (proportion) in top soil layer following to vegetation type;
!      veg_root_depth: root depth of vegetation (m), 2D array;
!      veg_roughness: vegetation cover roughness (m), 2D array;
!      veg_albedo: vegetation albedo (proportion), 2D array;
!      veg_emiss_1: vegetaion emissivity in broadband window (proportion), 2D array;
!      veg_emiss_2: vegetation emissivity in 8-12 micron window (proportion), 2D array.

implicit none

integer, parameter :: nvt=22

! Input:

integer :: nlon, nlat
real, dimension(nlon,nlat,nvt+1) :: vegeta
real, dimension(nlon,nlat) :: t_climate, t_amplitude

! Ouput:

real, dimension(nlon,nlat) :: veg_water_table, veg_qrel_surf, &
 veg_root_depth, veg_roughness, veg_albedo, veg_emiss_1, veg_emiss_2, snow_dirt

! Work:

integer, parameter :: ival_missing=-9999
real, dimension(nvt) :: vegwt1, vegwt2, qrelsrf1, qrelsrf2, vegrtzn, vegrogh, vegalb, vegemis1, vegemis2, snow_dirt_loc
real, dimension(nlon,nlat) :: zzz1, zzz2
integer :: itype, nsmooth, i, ip, jp, nv11
real :: zfrac, z1

! Ground water table depth (m) in "tropical" zones (for some type only: 12, 13, 14, 19)

data vegwt1/ &
! 1-Tree broadl.evergreen >15%  2-Tree broadl.deciduous >40% 3-Tree broadl.deciduous 15-40%
         2.,                          4.,                       10., &
! 4-Tree needl. evergreen >15%  5-Tree needl.decidous 15-40% 6-Tree mixed leaftype > 15%
         3.,                          3.,                        3., &
! 7-Tree >15% Swamp Forests     8-Tree >15% Mangrove Forests 9-Mosaic tree and other
         1.,                          0.,                        4., &
!10-Tree burnt forests         11-Shrub >15% evergreen      12-Shrub >15% deciduous broadl.
         4.,                          6.,                       15., &
!13-Herbaceous >15%            14-Sparse herbaceous/shurb   15-Flooded herbaceous/shurb > 15%
        20.,                         50.,                        2., &
!16-Cropland                   17-Mosaic crop/tree/other    18-Mosaic cropland/shub/herbaceous
!!!         8.,                          6.,                        8., &
        12.,                         10.,                       10., &
!19-Bare Areas                 20-Urban Areas               21-Snow or Ice
        50.,                         10.,                        0., &
!22-Water Bodies
         0./

! Ground water table depth (m) in "boreal" zones (for some type only: 12, 13, 14, 19)

data vegwt2/ &
! 1-Tree broadl.evergreen >15%  2-Tree broadl.deciduous >40% 3-Tree broadl.deciduous 15-40%
         2.,                          4.,                       10., &
! 4-Tree needl. evergreen >15%  5-Tree needl.decidous 15-40% 6-Tree mixed leaftype > 15%
         3.,                          3.,                        3., &
! 7-Tree >15% Swamp Forests     8-Tree >15% Mangrove Forests 9-Mosaic tree and other
         1.,                          0.,                        4., &
!10-Tree burnt forests         11-Shrub >15% evergreen      12-Shrub >15% deciduous broadl.
         4.,                          6.,                        2., &
!13-Herbaceous >15%            14-Sparse herbaceous/shurb   15-Flooded herbaceous/shurb > 15%
         2.,                          1.,                        2., &
!16-Cropland                   17-Mosaic crop/tree/other    18-Mosaic cropland/shub/herbaceous
!!!         8.,                          6.,                        8., &
        10.,                         10.,                       10., &
!19-Bare Areas                 20-Urban Areas               21-Snow or Ice
         1.,                         10.,                        0., &
!22-Water Bodies
         0./

! Relative soil water content (proportion) in top soil layer in "tropical" zones (for some type only: 12, 13, 14, 19)

data qrelsrf1/ &
! 1-Tree broadl.evergreen >15%  2-Tree broadl.deciduous >40% 3-Tree broadl.deciduous 15-40%
        0.7,                         0.4,                       0.3, &
! 4-Tree needl. evergreen >15%  5-Tree needl.decidous 15-40% 6-Tree mixed leaftype > 15%
        0.6,                         0.6,                       0.6, &
! 7-Tree >15% Swamp Forests     8-Tree >15% Mangrove Forests 9-Mosaic tree and other
        0.8,                         1.0,                       0.5, &
!10-Tree burnt forests         11-Shrub >15% evergreen      12-Shrub >15% deciduous broadl.
        0.5,                         0.3,                       0.3, &
!13-Herbaceous >15%            14-Sparse herbaceous/shurb   15-Flooded herbaceous/shurb > 15%
        0.3,                         0.2,                       0.8, &
!16-Cropland                   17-Mosaic crop/tree/other    18-Mosaic cropland/shub/herbaceous
        0.3,                         0.4,                       0.3, &
!19-Bare Areas                 20-Urban Areas               21-Snow or Ice
        0.1,                         0.5,                       1.0, &
!22-Water Bodies
        1.0/

! Relative soil water content (proportion) in top soil layer in "boreal" zones (for some type only: 12, 13, 14, 19)

data qrelsrf2/ &
! 1-Tree broadl.evergreen >15%  2-Tree broadl.deciduous >40% 3-Tree broadl.deciduous 15-40%
        0.7,                         0.4,                       0.3, &
! 4-Tree needl. evergreen >15%  5-Tree needl.decidous 15-40% 6-Tree mixed leaftype > 15%
        0.6,                         0.6,                       0.6, &
! 7-Tree >15% Swamp Forests     8-Tree >15% Mangrove Forests 9-Mosaic tree and other
        0.8,                         1.0,                       0.5, &
!10-Tree burnt forests         11-Shrub >15% evergreen      12-Shrub >15% deciduous broadl.
        0.5,                         0.3,                       0.6, &
!13-Herbaceous >15%            14-Sparse herbaceous/shurb   15-Flooded herbaceous/shurb > 15%
        0.6,                         0.9,                       0.8, &
!16-Cropland                   17-Mosaic crop/tree/other    18-Mosaic cropland/shub/herbaceous
        0.3,                         0.4,                       0.3, &
!19-Bare Areas                 20-Urban Areas               21-Snow or Ice
        0.9,                         0.5,                       1.0, &
!22-Water Bodies
        1.0/

! Depth of plant's root zone (m)

data vegrtzn/ &
! 1-Tree broadl.evergreen >15%  2-Tree broadl.deciduous >40% 3-Tree broadl.deciduous 15-40%
       1.50,                        1.50,                      1.50, &
! 4-Tree needl. evergreen >15%  5-Tree needl.decidous 15-40% 6-Tree mixed leaftype > 15%
       1.50,                        1.50,                      1.50, &
! 7-Tree >15% Swamp Forests     8-Tree >15% Mangrove Forests 9-Mosaic tree and other
       1.50,                        1.50,                      1.00, &
!10-Tree burnt forests         11-Shrub >15% evergreen      12-Shrub >15% deciduous broadl.
       1.50,                        0.50,                      0.50, &
!13-Herbaceous >15%            14-Sparse herbaceous/shurb   15-Flooded herbaceous/shurb > 15%
       0.15,                        0.30,                      0.15, &
!16-Cropland                   17-Mosaic crop/tree/other    18-Mosaic cropland/shub/herbaceous
       0.20,                        0.30,                      0.20, &
!19-Bare Areas                 20-Urban Areas               21-Snow or Ice
       0.00,                        0.30,                      0.00, &
!22-Water Bodies
       0.00/

data vegrogh/ &
! 1-Tree broadl.evergreen >15%  2-Tree broadl.deciduous >40% 3-Tree broadl.deciduous 15-40%
       1.50,                        1.50,                      1.70, &
! 4-Tree needl. evergreen >15%  5-Tree needl.decidous 15-40% 6-Tree mixed leaftype > 15%
       1.70,                        1.70,                      1.50, &
! 7-Tree >15% Swamp Forests     8-Tree >15% Mangrove Forests 9-Mosaic tree and other
       1.50,                        1.50,                      1.50, &
!10-Tree burnt forests         11-Shrub >15% evergreen      12-Shrub >15% deciduous broadl.
       1.60,                        0.20,                      0.20, &
!13-Herbaceous >15%            14-Sparse herbaceous/shurb   15-Flooded herbaceous/shurb > 15%
       0.05,                        0.10,                      0.05, &
!16-Cropland                   17-Mosaic crop/tree/other    18-Mosaic cropland/shub/herbaceous
       0.10,                        0.20,                      0.15, &
!19-Bare Areas                 20-Urban Areas               21-Snow or Ice
       0.05,                        0.70,                      0.05, &
!22-Water Bodies
       0.001/

! Vegetation albedo

data vegalb/ &
! 1-Tree broadl.evergreen >15%  2-Tree broadl.deciduous >40% 3-Tree broadl.deciduous 15-40%
       0.13,                        0.15,                      0.15, &
! 4-Tree needl. evergreen >15%  5-Tree needl.decidous 15-40% 6-Tree mixed leaftype > 15%
       0.14,                        0.15,                      0.14, &
! 7-Tree >15% Swamp Forests     8-Tree >15% Mangrove Forests 9-Mosaic tree and other
       0.13,                        0.14,                      0.16, &
!10-Tree burnt forests         11-Shrub >15% evergreen      12-Shrub >15% deciduous broadl.
       0.08,                        0.18,                      0.19, &
!13-Herbaceous >15%            14-Sparse herbaceous/shurb   15-Flooded herbaceous/shurb > 15%
       0.19,                        0.19,                      0.18, &
!16-Cropland                   17-Mosaic crop/tree/other    18-Mosaic cropland/shub/herbaceous
       0.20,                        0.19,                      0.19, &
!19-Bare Areas                 20-Urban Areas               21-Snow or Ice
       0.15,                        0.15,                      0.55, &
!22-Water Bodies
       0.07/

! Vegetation emissivity in broadband window

data vegemis1/ &
! 1-Tree broadl.evergreen >15%  2-Tree broadl.deciduous >40% 3-Tree broadl.deciduous 15-40%
      0.995,                       0.992,                     0.992, &
! 4-Tree needl. evergreen >15%  5-Tree needl.decidous 15-40% 6-Tree mixed leaftype > 15%
      0.997,                       0.992,                     0.993, &
! 7-Tree >15% Swamp Forests     8-Tree >15% Mangrove Forests 9-Mosaic tree and other
      0.994,                       0.994,                     0.990, &
!10-Tree burnt forests         11-Shrub >15% evergreen      12-Shrub >15% deciduous broadl.
      0.960,                       0.982,                     0.986, &
!13-Herbaceous >15%            14-Sparse herbaceous/shurb   15-Flooded herbaceous/shurb > 15%
      0.985,                       0.985,                     0.988, &
!16-Cropland                   17-Mosaic crop/tree/other    18-Mosaic cropland/shub/herbaceous
      0.985,                       0.985,                     0.985, &
!19-Bare Areas                 20-Urban Areas               21-Snow or Ice
      0.970,                       0.970,                     0.980, &
!22-Water Bodies
      0.970/

! Vegetation emissivity in 8-12 micron window

data vegemis2/ &
! 1-Tree broadl.evergreen >15%  2-Tree broadl.deciduous >40% 3-Tree broadl.deciduous 15-40%
      0.985,                       0.986,                     0.986, &
! 4-Tree needl. evergreen >15%  5-Tree needl.decidous 15-40% 6-Tree mixed leaftype > 15%
      0.992,                       0.985,                     0.987, &
! 7-Tree >15% Swamp Forests     8-Tree >15% Mangrove Forests 9-Mosaic tree and other
      0.984,                       0.984,                     0.985, &
!10-Tree burnt forests         11-Shrub >15% evergreen      12-Shrub >15% deciduous broadl.
      0.950,                       0.973,                     0.977, &
!13-Herbaceous >15%            14-Sparse herbaceous/shurb   15-Flooded herbaceous/shurb > 15%
      0.978,                       0.978,                     0.979, &
!16-Cropland                   17-Mosaic crop/tree/other    18-Mosaic cropland/shub/herbaceous
      0.975,                       0.980,                     0.978, &
!19-Bare Areas                 20-Urban Areas               21-Snow or Ice
      0.960,                       0.960,                     0.980, &
!22-Water Bodies
      0.970/

! Weight (proportion 0-1) of dirt growth of snow surface due to vegetation
! waste, aerosol deposition, etc., used in snow albedo definition

data snow_dirt_loc/ &
! 1-Tree broadl.evergreen >15%  2-Tree broadl.deciduous >40% 3-Tree broadl.deciduous 15-40%
         1.,                          1.,                        1., &
! 4-Tree needl. evergreen >15%  5-Tree needl.decidous 15-40% 6-Tree mixed leaftype > 15%
         1.,                          1.,                        1., &
! 7-Tree >15% Swamp Forests     8-Tree >15% Mangrove Forests 9-Mosaic tree and other
         1.,                          1.,                       0.8, &
!10-Tree burnt forests         11-Shrub >15% evergreen      12-Shrub >15% deciduous broadl.
         1.,                         0.7,                       0.7, &
!13-Herbaceous >15%            14-Sparse herbaceous/shurb   15-Flooded herbaceous/shurb > 15%
        0.5,                         0.6,                       0.5, &
!16-Cropland                   17-Mosaic crop/tree/other    18-Mosaic cropland/shub/herbaceous
        0.6,                         0.8,                       0.7, &
!19-Bare Areas                 20-Urban Areas               21-Snow or Ice
        0.0,                          1.,                       0.0, &
!22-Water Bodies
        0.0/

 print *
 print *,"Definition of physical parameters of vegetation (GLC2000 classification)"

 nsmooth = 0
! nsmooth = 1
 if (nsmooth > 0) print*, "nsmooth in presoil", nsmooth

! Definition of vegetation parameters for various vegeation types

 veg_water_table(:,:) = 0.
 veg_qrel_surf(:,:) = 0.
 veg_roughness(:,:) = 0.
 veg_albedo(:,:) = 0.
 veg_emiss_1(:,:) = 0.
 veg_emiss_2(:,:) = 0.
 snow_dirt(:,:) = 0.

 do jp=1,nlat
 do ip=1,nlon

   do itype=1,nvt
     nv11 = itype
     zfrac = vegeta(ip,jp,nv11+1)

! Filling water and glacier points by the type "Mosaic cover", for spatial smoothing only

     if (nv11 >= nvt-1.and.nsmooth.ne.0) nv11 = 9
  
     if (t_climate(ip,jp) < 277.) then ! very cold zones = "boreal" zones
       veg_water_table(ip,jp) = veg_water_table(ip,jp)+vegwt2(nv11)*zfrac 
       veg_qrel_surf(ip,jp) = veg_qrel_surf(ip,jp)+qrelsrf2(nv11)*zfrac
     else
       if (t_climate(ip,jp) < 285.) then ! cold zones
         if (t_amplitude(ip,jp) < 50.) then ! cold and humid zones = "boreal" zones
           veg_water_table(ip,jp) = veg_water_table(ip,jp)+vegwt2(nv11)*zfrac 
           veg_qrel_surf(ip,jp) = veg_qrel_surf(ip,jp)+qrelsrf2(nv11)*zfrac
         else ! cold and dry zones = not tropical deserts
           veg_water_table(ip,jp) = veg_water_table(ip,jp)+vegwt1(nv11)*zfrac
           veg_qrel_surf(ip,jp) = veg_qrel_surf(ip,jp)+qrelsrf1(nv11)*zfrac
         endif
       else ! warm zones = "tropical" zones 
         veg_water_table(ip,jp) = veg_water_table(ip,jp)+vegwt1(nv11)*zfrac
         veg_qrel_surf(ip,jp) = veg_qrel_surf(ip,jp)+qrelsrf1(nv11)*zfrac
       endif
     endif

     veg_root_depth(ip,jp) = veg_root_depth(ip,jp)+vegrtzn(nv11)*zfrac
     veg_roughness(ip,jp) = veg_roughness(ip,jp)+vegrogh(nv11)*zfrac
     veg_albedo(ip,jp) = veg_albedo(ip,jp)+vegalb(nv11)*zfrac
     veg_emiss_1(ip,jp) = veg_emiss_1(ip,jp)+vegemis1(nv11)*zfrac
     veg_emiss_2(ip,jp) = veg_emiss_2(ip,jp)+vegemis2(nv11)*zfrac
     snow_dirt(ip,jp) = snow_dirt(ip,jp)+snow_dirt_loc(nv11)*zfrac

   enddo

 enddo
 enddo

! Horizontal smoothing of vegetation parameters

 if (nsmooth.ne.0) then

   zzz1(:,:) = veg_water_table(:,:)
   call smooth_soil(zzz1,zzz2,nlon,nlat,0.2,nsmooth)
   veg_water_table(:,:) = zzz2(:,:)

   zzz1(:,:) = veg_qrel_surf(:,:)
   call smooth_soil(zzz1,zzz2,nlon,nlat,0.2,nsmooth)
   veg_qrel_surf(:,:) = zzz2(:,:)

   zzz1(:,:) = veg_root_depth(:,:)
   call smooth_soil(zzz1,zzz2,nlon,nlat,0.2,nsmooth)
   veg_root_depth(:,:) = zzz2(:,:)

   zzz1(:,:) = veg_roughness(:,:)
   call smooth_soil(zzz1,zzz2,nlon,nlat,0.2,nsmooth)
   veg_roughness(:,:) = zzz2(:,:)

   zzz1(:,:) = veg_emiss_1(:,:)
   call smooth_soil(zzz1,zzz2,nlon,nlat,0.2,nsmooth)
   veg_emiss_1(:,:) = zzz2(:,:)

   zzz1(:,:) = veg_emiss_2(:,:)
   call smooth_soil(zzz1,zzz2,nlon,nlat,0.2,nsmooth)
   veg_emiss_2(:,:) = zzz2(:,:)

   zzz1(:,:) = snow_dirt(:,:)
   call smooth_soil(zzz1,zzz2,nlon,nlat,0.2,nsmooth)
   snow_dirt(:,:) = zzz2(:,:)

 endif

! Water body and Snow or Ice parameters

 do jp=1,nlat
 do ip=1,nlon

   nv11 = nint(vegeta(ip,jp,1))
   if (nv11 >= nvt-1) then
     veg_water_table(ip,jp) = vegwt1(nv11)
     veg_qrel_surf(ip,jp) = qrelsrf1(nv11)
     veg_root_depth(ip,jp) = vegrtzn(nv11)
     veg_roughness(ip,jp) = vegrogh(nv11)
     veg_albedo(ip,jp) = vegalb(nv11)
     veg_emiss_1(ip,jp) = vegemis1(nv11)
     veg_emiss_2(ip,jp) = vegemis2(nv11)
     snow_dirt(ip,jp) = snow_dirt_loc(nv11)
   endif
   snow_dirt(ip,jp)=min( max( snow_dirt(ip,jp), 0.), 1.)
 enddo
 enddo

return
end subroutine definition_veg2_const_phys_param
!=======================================================================
subroutine definition_soil_bottom(soil_type, veg_type, t_amplitude, t_climate, horog, &
 veg_water_table, soil_water_table, water_table, &
 qrel_surf, soil_qmax, soil_qmin, soil_c, soil_rho, soil_psi, soil_par_b, &
 lat_geo, level_soil, nlon, nlat, nst, nvt, nlevel_soil_tot, &
 nlevel_soil, nlevel_soil_water, depth, t_bottom, q_bottom, qrel_bottom)

! Procedure defines bottom condition of soil:
! index of bottom soil level at each grid point and bottom relative soil water content;
! at the bottom soil level climate temperature and bottom relative soil water content
! will be used as the bottom condition.

! Input variables:

!      soil_type:  soil types (classification FAO, 31 types) defined on the grid points, 3D array,
!                1-st and 2-nd array index are grid point index, 3-rd index is soil type number plus 1,
!                if 3-rd index has value 1, then the variable means dominate soil type,
!                values 2-32 of 3-rd index means fraction (of area, proportion) of the each of 31 soil types;

!      veg_type:   vegetation (landuse) types (classification GLC2000, 22 types) defined on the grid points, 3D array,
!               1-st and 2-nd array index are grid point index, 3-rd index is vegetation types number plus 1,
!               if 3-rd index has value 1, then the variable means dominate vegetation type,
!               values 2-23 of 3-rd index means fraction (of area, proportion) of the each of 14 vegetation types

!      t_amplitude: climate amplitude of temperature (K);
!      t_climate: climate temperature (K);
!      horog: orography height (m);
!      veg_water_table: depth of ground water table (m) following to vegetation type;
!      soil_water_table: depth of ground water table (m) following to soil type;
!      water_table: depth of ground water table (m);
!      qrel_surf: relative soil water content (proportion) in top soil layer;
!      soil_qmax: maximum specific volumetric soil water content (m^3/m^3) at all soil levels;
!      soil_qmin: minimum specific volumetric soil water content (m^3/m^3) at all soil levels;
!      soil_c: dry soil thermal capacity (J/K/m^3) at all soil levels;
!      soil_rho: dry soil density (kg/m^3) at all soil levels;
!      soil_psi: idraulic potential of saturated soil (m) at all soil levels;
!      soil_par_b: parameter in soil water transport equation (dimensionless) at all soil levels;

!      lat_geo:   geographical latitude of the grid points (degree, -90...90)

!      level_soil:   soil levels values (m) at total (maximum possible) soil levels;

!      nlon:  grid point number along axis of rotated longitude;
!      nlat:  grid point number along axis of rotated latidude;
!      nst:   total number of soil types;
!      nvt:   total number of vegetation types;
!      nlevel_soil_tot: number (maximum possible) of soil levels.

! Output variables:

!      nlevel_soil: index of bottom (climate) soil level for heat transport at each grid point;
!      nlevel_soil_water: index of bottom soil level for water transport at each grid point;
!      depth: depth (m) of bottom (climate) soil level at each grid point;
!      t_bottom: soil temperature (K) at bottom soil level;
!      q_bottom: soil water content (m^3/m^3) at bottom soil level;
!      qrel_bottom: relative soil water content (proportion) at bottom soil level..

implicit none

! Input:

integer :: nlon, nlat, nst, nvt, nlevel_soil_tot
real, dimension(nlon,nlat,nst+1) :: soil_type
real, dimension(nlon,nlat,nvt+1) :: veg_type
real, dimension(nlon,nlat) :: t_amplitude, t_climate, horog, veg_water_table, soil_water_table, qrel_surf, lat_geo
real, dimension(nlon,nlat,nlevel_soil_tot) :: soil_qmax, soil_qmin, soil_c, soil_rho, soil_psi, soil_par_b
real, dimension(nlevel_soil_tot) :: level_soil

! Output:

integer, dimension(nlon,nlat) :: nlevel_soil, nlevel_soil_water
real, dimension(nlon,nlat) :: depth, water_table, t_bottom, q_bottom, qrel_bottom

! Internal:

real, parameter :: water_table_veg_weight=0.7
integer, parameter :: nveg_water_table=3
integer, dimension(nveg_water_table) :: vegtype_water_table=(/16, 17, 18/)
real, dimension(nveg_water_table) :: weight_work
integer :: i, ii, j, lev, lev_ref, iflag, itype
real :: depth_ref=1.5, period=31536000., pi, amplitude_ref=1., depth_prec_season=4., depth_min=2., &
 qrel_ref, q, psi, lambda, kappa, cw=4186.8, rhow=1000., frac_bare, zdepth, z1, weight, weight_orog

 print *
 print *,"Definition of bottom condition for soil variables"

 pi=abs(acos(-1.))

! Reference soil level

 do lev=1,nlevel_soil_tot
   if (level_soil(lev) >= depth_ref) exit
 enddo
 lev_ref=max( min(lev, nlevel_soil_tot), 1)

 do j=1,nlat
 do i=1,nlon

   if (soil_psi(i,j,1) < 0.) then ! Not water body or glacier

! Reference relative soil water content

     qrel_ref=max(qrel_surf(i,j), 0.5)

!  Reference soil water content

     q=qrel_ref*(soil_qmax(i,j,lev_ref)-soil_qmin(i,j,lev_ref))+soil_qmin(i,j,lev_ref)

! Reference soil hydraulic potential

     psi=soil_psi(i,j,lev_ref)*((soil_qmax(i,j,lev_ref)/q)**soil_par_b(i,j,lev_ref))

! Thermal conductivity (dynamic) is the function of the soil hydraulic potential
! Soil thermal conductivity at reference level is used for the whole soil column

     lambda=min( 0.5+1400.*exp(-6.-alog10(abs(psi))), 2.5)

! Kinematic thermal conductivity of soil

     kappa=lambda/(soil_c(i,j,lev_ref)*soil_rho(i,j,lev_ref)+q*cw*rhow)

! Depth of attenuation of season temperature oscilation to refenence value

     depth(i,j)=(-alog((amplitude_ref/t_amplitude(i,j)**2)))/(2.*sqrt(pi/(kappa*period)))

     depth(i,j)=max(depth(i,j), depth_min)

! Control of presence of vegeation with is typical for expressed season precipitation regime

     iflag=0

     itype=2 ! Tree cover, broadleaved deciduous, closed (>40%)
     if (veg_type(i,j,itype+1) > 0.3) iflag=1

     itype=3 ! Tree cover, broadleaved deciduous, open (15-40%)
     if (veg_type(i,j,itype+1) > 0.3) iflag=1

     itype=11 ! Shrubcover, closed to open (>15%), evergreen(broadleaved or needleleaved)
     if (veg_type(i,j,itype+1) > 0.3) iflag=1

     if (iflag == 1) depth(i,j)=max(depth(i,j), depth_prec_season)

   else ! Water body or glacier

     depth(i,j)=100.

     iflag=0
     itype=nst-1 ! glacier
     if (soil_type(i,j,itype+1) > 0.5) iflag=1
     if (iflag == 1) depth(i,j)=10.

   endif

! Index of soil level with be used as bottom condition level for heat transport simulation

   do lev=1,nlevel_soil_tot
     if (level_soil(lev) >= depth(i,j)) exit
   enddo
   nlevel_soil(i,j)=min(lev, nlevel_soil_tot) 

! Index of soil level with be used as bottom condition level for water transport simulation

   if (soil_psi(i,j,1) < 0.) then ! Not water body or glacier
!     itype=19 ! Bare Areas 
!     frac_bare=veg_type(i,j,itype+1)
!     water_table(i,j)=veg_water_table(i,j)*(1.-frac_bare)+soil_water_table(i,j)*frac_bare 
!     water_table(i,j)=min(veg_water_table(i,j), soil_water_table(i,j)) 
     weight=water_table_veg_weight
     weight_work(1:nveg_water_table)=water_table_veg_weight
     do ii=1,nveg_water_table
       itype=vegtype_water_table(ii)+1
       weight_work(ii)=veg_type(i,j,itype)*(0.2-water_table_veg_weight)+water_table_veg_weight
     enddo
     weight=minval(weight_work(1:nveg_water_table))
     water_table(i,j)=veg_water_table(i,j)*weight+soil_water_table(i,j)*(1.-weight) 
     if (horog(i,j) <= 500.) then
       weight_orog=0.
     elseif (horog(i,j) >=1500.) then
       weight_orog=1.
     else
       weight_orog=1./1000.*horog(i,j)-500./1000.
     endif
     water_table(i,j)=water_table(i,j)*(1.-weight_orog)+50.*weight_orog
   else ! Water body or glacier
     water_table(i,j)=0.
   endif

   do lev=1,nlevel_soil_tot
     if (level_soil(lev) >= water_table(i,j)) exit
   enddo
   nlevel_soil_water(i,j)=min(lev, nlevel_soil_tot) 
   nlevel_soil_water(i,j)=min(nlevel_soil_water(i,j), nlevel_soil(i,j)) 

! Bottom condition of temperature and water content

   if (soil_psi(i,j,1) < 0.) then ! Not water body or glacier

     t_bottom(i,j)=t_climate(i,j)

     zdepth=level_soil(nlevel_soil(i,j))

     if (water_table(i,j) > zdepth) then
       z1=zdepth/water_table(i,j)
       qrel_bottom(i,j)=min(1.*z1+qrel_surf(i,j)*(1.-z1), 1.)
     else
       qrel_bottom(i,j)=1.
     endif

     lev=nlevel_soil_water(i,j)
     q_bottom(i,j)=qrel_bottom(i,j)*(soil_qmax(i,j,lev)-soil_qmin(i,j,lev))+soil_qmin(i,j,lev)

   else ! Water body or glacier

     t_bottom(i,j)=t_climate(i,j)

     iflag=0
     itype=nst-1 ! glacier
     if (soil_type(i,j,itype+1) > 0.5) iflag=1
     if (iflag == 1) t_bottom(i,j)=min(t_climate(i,j), 263.15)

     qrel_bottom(i,j)=1.
     nlevel_soil_water(i,j)=0

   endif

 enddo
 enddo

return
end subroutine definition_soil_bottom
!=======================================================================
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

 print *
 print *,"Definition of radiation parameters of the surface"

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
!=======================================================================
      subroutine presoil(suolo,vegeta,soilvegpar,alat,dlat,nlon,nlat,imonthc,idayc)

! Defines physical parameters of soil and vegetation on the model grid.
! They are stored in the 3-D matrix soilvegpar(nlon,nlat,50).

! NST: number of soil types, NVT is vegetation type number

      integer, parameter :: nst=15, nvt=14

      real, dimension(nlon,nlat,nst+1) :: suolo
      real, dimension(nlon,nlat,nvt+1) :: vegeta
      real, dimension(nlon,nlat,50) :: soilvegpar
      real, dimension(nlon,nlat) :: alat, zzz1, zzz2
      integer, dimension(12) :: imon

      real, dimension(nst) :: qgmax, qgmin, rogdr, cgdry, psig0, bghyd, ykghy,                         &
                              qglwt1, qglwt2, qglwt3, qglwt4, qglwt5, qglwt6,  qglwt7, qglwt8, qglwt9, &
                              qgwilt, qgrefw, soildalbed, soilwalbed,                                  &
                              soildemis1, soilwemis1, soildemis2, soilwemis2, zparqwf
      real, dimension(nvt) :: vegrtzn, vegrogh, vegalb, vegemis1, vegemis2, &
             veglain, vegfrcn, veglais, vegfrcs, veglaic, vegfrcc

! Texture soil types required by the soil process parametrization:

!  1 sand,        2 loamy sand,       3 sandy loam,       4 silty loam,
!  5 loam,        6 sandy clay loam,  7 silty clay loam,  8 clay loam,
!  9 sandy clay, 10 silty clay,      11 clay,            12 peat,
! 13 stone,      14 glaciers,        15 water body.

! Each of the following data vectors contains one parameter type for each of the above
! soil types (including water body, used by Dmitriy Pressman).
! Values are derived mostly from Pielke's (2001) book - if not available were 'invented'

! Saturated moisture potential of the soil (at QG=QGMAX) (m)
!*** PSIG0(13) (stone) defined by Oksana Drofa

      data psig0/-0.121,-0.090,-0.218,-0.786,-0.478,-0.299,-0.356,    &
          -0.630,-0.153,-0.490,-0.405,-0.356,-0.090,-5.171,-1193.463/

! Exponent B (dimensionless) for moisture potential equation:
! PSI=PSIG*(QGMAX/QG)**B
!*** BGHYD(13) defined by Oksana Drofa

      data bghyd/4.05,4.38,4.90,5.30,5.39,7.12,7.75, &
       8.52,10.40,10.40,11.40,7.75,4.00,0.00,0.00/

! Dry soil density (kg/m**3)
!*** ROGDR(1), ROGDR(10), ROGDR(11) are defined to obtain values Rog*Cg
!*** from Pielke; ROGDR(13) defined by Oksana Drofa

      data rogdr/1.55e3,1.57e3,1.51e3,1.37e3,1.464e3,1.55e3,1.40e3, &
           1.40e3,1.53e3,1.57e3,1.60e3,0.40e3,1.70e3,0.90e3,1.00e3/

! Dry soil heat capacity (J/kg/K)
!*** CGDRY(1), CGDRY(10), CGDRY(11) are defined to have values Rog*Cg
!*** from Pielke; and CGDRY(13) defined by Oksana Drofa;

      data cgdry/0.948e3,0.898e3,0.887e3,0.927e3,0.826e3,  &
                 0.761e3,0.951e3,0.878e3,0.771e3,0.732e3,  &
                 0.681e3,2.100e3,0.950e3,2.0935e3,4.186e3/

! Saturated hydraulic conductivity (m/s)

! Original (Pielke) values

!      data ykghy/1.760e-4,1.563e-4,0.341e-4,0.072e-4,0.070e-4,  &
!                 0.063e-4,0.017e-4,0.025e-4,0.022e-4,0.010e-4,  &
!                 0.013e-4,0.080e-4,2.300e-4,0.000e0,0.000e0/

! Values redefined after rescaling with min. and max. values used at ECMWF

!      data ykghy/6.97e-6,6.38e-6,2.15e-6,7.74e-7,7.61e-7,  &
!                 7.13e-7,3.38e-7,4.16e-7,3.88e-7,2.60e-7,  &
!                 2.95e-7,8.26e-7,8.50e-6,0.000e0,0.000e0/

! Intermediate values (sept. 2010)

!      data ykghy/6.0e-5,2.5e-5,7.0e-6,1.7e-6,1.70e-6,  &
!                 1.6e-6,7.0e-7,7.5e-7,6.8e-7,4.1e-7,   &
!                 4.8e-7,1.80e-6,5.0e-7,0.0e0,0.0e0/

! Revised values (feb. 2014)

      data ykghy/2.7e-5, 1.6e-5, 7.0e-6, 1.4e-6, 1.4e-6, &
                 1.3e-6, 7.0e-7, 7.5e-7, 4.9e-7, 4.7e-7, &
                 4.8e-7, 1.4e-6, 5.0e-7, 0.0e0,  0.0e0/

! Maximum volumetric water content of the soil (m**3/m**3) (qgmax*row is soil porosity)
!*** QGMAX(13) defined by Oksana Drofa

      data qgmax/0.395,0.410,0.435,0.485,0.451,0.420,0.477,0.476,  &
                 0.426,0.492,0.482,0.600,0.300,1.000,1.000/

! Minimum volumetric water content of the soil (m**3/m**3)
!*** QGMIN(1), QGMIN(10), QGMIN(11), QGMIN(13) defined by Oksana Drofa

      data qgmin/0.015,0.023,0.050,0.090,0.060,0.085,0.113,0.137,  &
                 0.134,0.140,0.145,0.188,0.020,0.000,0.000/

! Volumetric water content of the soil at wilting point (m**3/m**3)
!*** QGWILT(1),QGWILT(10),QGWILT(11),QGWILT(13) defined by Oksana Drofa

      data qgwilt/0.0677,0.0750,0.1142,0.1794,0.1547,0.1749,0.2181,  &
           0.2498,0.2193,0.2832,0.2864,0.3947,0.0300,0.0000,0.0000/

! Volumetric water content of soil at reference point for evapotranspiration (m**3/m**3)
!*** These QGREFW values are given by Dmitriy Y. Pressman
!*** QGREFW(1),QGREFW(10),QGREFW(11),QGREFW(13) defined by Oksana Drofa

      data qgrefw/0.270,0.283,0.312,0.412,0.329,0.315,0.387,0.382,  &
                  0.338,0.380,0.370,0.450,0.200,0.000,0.000/

! QGLWi are maximum volumetric liquid (fluid) water contents of the soil
! (m**3/m**3) at different temperatures below 0 C (defined below)
! ZPARQWF is a parameter used for the QGLWi definition

      data zparqwf/1.00,0.95,0.80,0.90,0.90,0.80,0.80,0.60,  &
                   0.50,0.40,0.30,0.70,1.20,1.00,1.00/

! Dry soil (qg=qgmin) emissivity in broadband window

      data soildemis1/0.930,0.950,0.960,0.965,0.975,0.965,0.970,0.960,  &
                     0.950,0.955,0.940,0.965,0.950,0.980,0.970/

! Dry soil (qg=qgmin) emissivity in 8-12 micron window

      data soildemis2/0.860,0.880,0.950,0.950,0.960,0.955,0.965,0.945,  &
                     0.940,0.945,0.932,0.960,0.945,0.980,0.970/

! Wet soil (qg=qgmax) emissivity in broadband window

      data soilwemis1/0.960,0.972,0.985,0.985,0.995,0.985,0.990,0.985, &
                     0.975,0.978,0.968,0.988,0.965,0.980,0.970/

! Wet soil (qg=qgmax) emissivity in 8-12 micron window

      data soilwemis2/0.890,0.930,0.955,0.970,0.985,0.975,0.980,0.965, &
                    0.960,0.965,0.950,0.982,0.958,0.980,0.970/

! Dry soil (qg=qgmin) albedo

      data soildalbed/0.45,0.35,0.25,0.15,0.10,0.18,0.17,0.22,  &
                     0.25,0.28,0.35,0.15,0.18,0.55,0.07/

! Wet soil (qg=qgmax) albedo

      data soilwalbed/0.25,0.20,0.18,0.08,0.05,0.10,0.13,0.15,  &
                     0.20,0.15,0.20,0.05,0.08,0.55,0.07/

! Vegetation types (landcover types) that are requed by the soil model:
! (see subroutine deflai for a more precise definition of vegetation types below)

!  1 evergreen needleleaf forest,  2 evergreen broadleaf forest,
!  3 deciduous needleleaf forest,  4 deciduous broadleaf forest,
!  5 mixed cover,                  6 woodland,
!  7 wooded grassland,             8 closed shrubland,
!  9 open shrubland,              10 grassland,
! 11 cropland,                    12 bare ground,
! 13 urban and built-up,          14 water body (sea, lakes, rivers).

! Depth of plant's root zone (m)

      data vegrtzn/1.00,1.25,1.00,1.25,0.80,0.90,0.70,0.50,  &
                   0.45,0.35,0.40,0.00,0.60,0.00/

! Vegetation cover roughness

      data vegrogh/0.80,1.00,0.80,0.90,0.50,0.60,0.40,  &
                   0.18,0.14,0.07,0.12,0.05,1.20,0.001/

! Vegetation emissivity in broadband window

      data vegemis1/0.997,0.995,0.991,0.992,0.990,0.987,0.984,  &
                    0.986,0.982,0.985,0.983,0.963,0.986,0.990/

! Vegetation emissivity in 8-12 micron window

      data vegemis2/0.993,0.985,0.983,0.988,0.982,0.975,0.975,  &
                    0.975,0.970,0.978,0.983,0.958,0.975,0.990/

! Vegetation albedo

      data vegalb/0.14,0.13,0.15,0.15,0.13,0.15,0.18,  &
                  0.19,0.20,0.19,0.19,0.20,0.15,0.07/

      data imon /31,28,31,30,31,30,31,31,30,31,30,31/

! DEFLAY defines vegetation parameters (LAI and veget. fraction) that are
! function of the day of the year

! Date for the Northern Hemisphere season

      ndayn = 0
      do i=1,imonthc-1
       ndayn = ndayn+imon(i)
      enddo
      ndayn = ndayn+idayc

! Date for the Southern Hemisphere season

      ndays = 0
      do i=1,imonthc+6-1
       if (i.le.12) then
         ndays = ndays+imon(i)
       else
         ndays = ndays+imon(i-12)
       endif
      enddo
      ndays = ndays+idayc
      if (ndays.gt.365) ndays = ndays-365

      zndayc = ndayn
      call deflai(veglain,vegfrcn,zndayc)
      zndayc = ndays
      call deflai(veglais,vegfrcs,zndayc)

      zdayref = 183. ! year half, 2 July
      zalatref = 10.

! Definition of QGLWTi (partly fitting empirical values using analytical functions)
! Maximum volumetric liquid (fluid) water content of the soil
! (m**3/m**3) at different temperatures below 0 C

      do k=1,13
      qglwt1(k) = qgmax(k)*0.60 ! T  0
      qglwt2(k) = qgmax(k)*0.50 ! T -1
      qglwt3(k) = qgmax(k)*0.40 ! T -2
      qglwt4(k) = qgmax(k)*0.30 ! T -3
      qglwt5(k) = qgmax(k)*0.25 ! T -4
      qglwt6(k) = qgmax(k)*0.20 ! T -5
      qglwt7(k) = qgmax(k)*0.10 ! T -10
      qglwt8(k) = qgmax(k)*0.05 ! T -15
      qglwt9(k) = qgmax(k)*0.00 ! T -20
      enddo
      do k=14,nst
      qglwt9(k) = 0.
      enddo

!      nsmooth = 0
      nsmooth = 1
      print*, "nsmooth in presoil", nsmooth

! Definition of physical parameters of soil

      do 25 jp=1,nlat
      do 25 ip=1,nlon
      do k=1,50
      soilvegpar(ip,jp,k) = 0.
      enddo
      do k=1,nst
      n11 = k
      zfrac = suolo(ip,jp,n11+1)

! Weighted (with soil fraction) average of the various soil parameters.
! The following fills glacier and water points by sandy loam, for subsequent spatial smoothing only

      if (n11 >= nst-1.and.nsmooth /= 0) n11 = 3

      soilvegpar(ip,jp, 1) = soilvegpar(ip,jp, 1)+psig0(n11)*zfrac
      soilvegpar(ip,jp, 2) = soilvegpar(ip,jp, 2)+bghyd(n11)*zfrac
      soilvegpar(ip,jp, 3) = soilvegpar(ip,jp, 3)+rogdr(n11)*zfrac
      soilvegpar(ip,jp, 4) = soilvegpar(ip,jp, 4)+cgdry(n11)*zfrac
      soilvegpar(ip,jp, 5) = soilvegpar(ip,jp, 5)+ykghy(n11)*zfrac
      soilvegpar(ip,jp, 6) = soilvegpar(ip,jp, 6)+qgmax(n11)*zfrac
      soilvegpar(ip,jp, 7) = soilvegpar(ip,jp, 7)+qgmin(n11)*zfrac
      soilvegpar(ip,jp, 8) = soilvegpar(ip,jp, 8)+qgwilt(n11)*zfrac
      soilvegpar(ip,jp, 9) = soilvegpar(ip,jp, 9)+qgrefw(n11)*zfrac

! soilvegpar(:,:,10) - soilvegpar(:,:,18) define the volumetric content (m**3/m**3)
! of liquid water in soil at temperature values:
! 0, -1, -2, -3, -4, -5, -10, -15, -20 deg. Celsius

      soilvegpar(ip,jp,10) = soilvegpar(ip,jp,10)+qglwt1(n11)*zfrac
      soilvegpar(ip,jp,11) = soilvegpar(ip,jp,11)+qglwt2(n11)*zfrac
      soilvegpar(ip,jp,12) = soilvegpar(ip,jp,12)+qglwt3(n11)*zfrac
      soilvegpar(ip,jp,13) = soilvegpar(ip,jp,13)+qglwt4(n11)*zfrac
      soilvegpar(ip,jp,14) = soilvegpar(ip,jp,14)+qglwt5(n11)*zfrac
      soilvegpar(ip,jp,15) = soilvegpar(ip,jp,15)+qglwt6(n11)*zfrac
      soilvegpar(ip,jp,16) = soilvegpar(ip,jp,16)+qglwt7(n11)*zfrac
      soilvegpar(ip,jp,17) = soilvegpar(ip,jp,17)+qglwt8(n11)*zfrac
      soilvegpar(ip,jp,18) = soilvegpar(ip,jp,18)+qglwt9(n11)*zfrac

      soilvegpar(ip,jp,19) = soilvegpar(ip,jp,19)+soildemis1(n11)*zfrac
      soilvegpar(ip,jp,20) = soilvegpar(ip,jp,20)+soildemis2(n11)*zfrac
      soilvegpar(ip,jp,21) = soilvegpar(ip,jp,21)+soilwemis1(n11)*zfrac
      soilvegpar(ip,jp,22) = soilvegpar(ip,jp,22)+soilwemis2(n11)*zfrac
      soilvegpar(ip,jp,23) = soilvegpar(ip,jp,23)+soildalbed(n11)*zfrac
      soilvegpar(ip,jp,24) = soilvegpar(ip,jp,24)+soilwalbed(n11)*zfrac
      enddo ! k=1,nst
   25 continue

      if(nsmooth.ne.0) then
      do k=1,24
      do jp=1,nlat
      do ip=1,nlon
      zzz1(ip,jp) = soilvegpar(ip,jp,k)
      enddo
      enddo

! Smoothing of adjacent parameter values

      call smooth_soil(zzz1,zzz2,nlon,nlat,0.2,nsmooth)
      do jp=1,nlat
      do ip=1,nlon
      soilvegpar(ip,jp,k) = zzz2(ip,jp)
      enddo
      enddo
      enddo
      endif

! Definition of physical parameters of vegetation (similar procedure as for soil)

      do jp=1,nlat
      do ip=1,nlon
      do k=1,nvt
      nv11 = k
      zfrac = vegeta(ip,jp,nv11+1)

! Filling water points by mixed cover, for spatial smoothing only

      if (nv11.eq.nvt.and.nsmooth.ne.0) nv11 = 10

! Definition of changing (in time) parameter of vegetation
! as a function of date and latitude (hemisphere).
! Season variability is decreased in tropical zone between
! -zalatref and zalatref latitudes. At the equator the date is
! constant and equal to zdayref.

      if (alat(ip,jp).ge.zalatref.or.alat(ip,jp).le.-zalatref) then
        if (alat(ip,jp).ge.0.) then
          soilvegpar(ip,jp,31) = soilvegpar(ip,jp,31)+veglain(nv11)*zfrac
          soilvegpar(ip,jp,32) = soilvegpar(ip,jp,32)+vegfrcn(nv11)*zfrac
        else
          soilvegpar(ip,jp,31) = soilvegpar(ip,jp,31)+veglais(nv11)*zfrac
          soilvegpar(ip,jp,32) = soilvegpar(ip,jp,32)+vegfrcs(nv11)*zfrac
        endif
      else
        if (alat(ip,jp).ge.0.) then
          zndayc = ndayn
          zndayc = zdayref+(zndayc-zdayref)*alat(ip,jp)/zalatref
        else
          zndayc = ndays
          zndayc = zdayref-(zndayc-zdayref)*alat(ip,jp)/zalatref
        endif
        call deflai(veglaic,vegfrcc,zndayc)
        soilvegpar(ip,jp,31) = soilvegpar(ip,jp,31)+veglaic(nv11)*zfrac
        soilvegpar(ip,jp,32) = soilvegpar(ip,jp,32)+vegfrcc(nv11)*zfrac
      endif

      soilvegpar(ip,jp,33) = soilvegpar(ip,jp,33)+vegrtzn(nv11)*zfrac
      soilvegpar(ip,jp,34) = soilvegpar(ip,jp,34)+vegrogh(nv11)*zfrac
      soilvegpar(ip,jp,35) = soilvegpar(ip,jp,35)+vegemis1(nv11)*zfrac
      soilvegpar(ip,jp,36) = soilvegpar(ip,jp,36)+vegemis2(nv11)*zfrac
      soilvegpar(ip,jp,37) = soilvegpar(ip,jp,37)+vegalb(nv11)*zfrac
      enddo
      enddo
      enddo

      if(nsmooth.ne.0) then
      do k=31,37
      do jp=1,nlat
      do ip=1,nlon
      zzz1(ip,jp) = soilvegpar(ip,jp,k)
      enddo
      enddo
      call smooth_soil(zzz1,zzz2,nlon,nlat,0.2,nsmooth)
      do jp=1,nlat
      do ip=1,nlon
      soilvegpar(ip,jp,k) = zzz2(ip,jp)
      enddo
      enddo
      enddo
      endif

! Water body and glacier parameters

      do jp=1,nlat
      do ip=1,nlon
      n11 = suolo(ip,jp,1)
      if (n11.ge.nst-1) then
        soilvegpar(ip,jp, 1) = psig0(n11)
        soilvegpar(ip,jp, 2) = bghyd(n11)
        soilvegpar(ip,jp, 3) = rogdr(n11)
        soilvegpar(ip,jp, 4) = cgdry(n11)
        soilvegpar(ip,jp, 5) = ykghy(n11)
        soilvegpar(ip,jp, 6) = qgmax(n11)
        soilvegpar(ip,jp, 7) = qgmin(n11)
        soilvegpar(ip,jp, 8) = qgwilt(n11)
        soilvegpar(ip,jp, 9) = qgrefw(n11)
        soilvegpar(ip,jp,10) = qglwt1(n11)
        soilvegpar(ip,jp,11) = qglwt2(n11)
        soilvegpar(ip,jp,12) = qglwt3(n11)
        soilvegpar(ip,jp,13) = qglwt4(n11)
        soilvegpar(ip,jp,14) = qglwt5(n11)
        soilvegpar(ip,jp,15) = qglwt6(n11)
        soilvegpar(ip,jp,16) = qglwt7(n11)
        soilvegpar(ip,jp,17) = qglwt8(n11)
        soilvegpar(ip,jp,18) = qglwt9(n11)
        soilvegpar(ip,jp,19) = soildemis1(n11)
        soilvegpar(ip,jp,20) = soildemis2(n11)
        soilvegpar(ip,jp,21) = soilwemis1(n11)
        soilvegpar(ip,jp,22) = soilwemis2(n11)
        soilvegpar(ip,jp,23) = soildalbed(n11)
        soilvegpar(ip,jp,24) = soilwalbed(n11)
        soilvegpar(ip,jp,31) = veglain(nvt)
        soilvegpar(ip,jp,32) = vegfrcn(nvt)
        soilvegpar(ip,jp,33) = vegrtzn(nvt)
        soilvegpar(ip,jp,34) = vegrogh(nvt)
        soilvegpar(ip,jp,35) = soildemis1(nvt)
        soilvegpar(ip,jp,36) = soildemis2(nvt)
        soilvegpar(ip,jp,37) = soildalbed(nvt)
      endif
      soilvegpar(ip,jp, 7) = min(soilvegpar(ip,jp, 7), soilvegpar(ip,jp, 6)-1.e-8)
      soilvegpar(ip,jp, 8) = min(soilvegpar(ip,jp, 8), soilvegpar(ip,jp, 9)-1.e-8)
      enddo
      enddo

      return
      end
!=======================================================================
       subroutine surfradpar_old(fmask,soil_fao,vegeta,soilvegpar,qg1,fice,alat,alat0,dlon,dlat, &
                             nlon,nlat,nst,nvt,iflag,albedo,emismap1,emismap2)

! Defines radiative parameters of land surface as a function of soil and vegetation type
! spatial distribution, soil wetness, and soil and vegetation radiative parameters
! Output: ALBEDO, EMISMAP1, EMISMAP2

! Input:

       integer :: nlon, nlat, nst, nvt, iflag
       real, dimension(nlon,nlat) :: fmask, qg1, fice, alat
       real, dimension(nlon,nlat,nst+1) :: soil_fao
       real, dimension(nlon,nlat,nvt+1) :: vegeta
       real, dimension(nlon,nlat,50) :: soilvegpar
       real :: alat0, dlon, dlat

! Output :

       real, dimension(nlon,nlat) :: albedo, emismap1, emismap2

! Work:

       real, dimension(nlon,nlat) :: zzz

! Emissivity over desert is low - to better define it, it is necessary to define here
! provisionally a 'desert';
! 50000 m is the mimimum spatial distance for desert recognition as land
! surface without vegetation cover (VEGETA(I,J,1)=12)

      pi = abs(acos(-1.))
      ddx = dlon*6371.e+3*pi/180.
      ddy = dlat*6371.e+3*pi/180.

      iurbveg = 13 ! urban soil index

      do jlat=1,nlat
      do jlon=1,nlon

      nv1 = vegeta(jlon,jlat,1)
      ns1 = soil_fao(jlon,jlat,1)

      zvegdol   = soilvegpar(jlon,jlat,32)  ! vegetation fraction
      zemveg1   = soilvegpar(jlon,jlat,35)  ! veget. broadband window emissiv.
      zemveg2   = soilvegpar(jlon,jlat,36)  ! veget. 8-12 micron window emissiv.
      zemsoild1 = soilvegpar(jlon,jlat,19)  ! dry soil broadband window emissiv.
      zemsoild2 = soilvegpar(jlon,jlat,20)  ! dry soil 8-12 micron window emissiv.
      zemsoilw1 = soilvegpar(jlon,jlat,21)  ! wet soil broadband window emissiv.
      zemsoilw2 = soilvegpar(jlon,jlat,22)  ! wet soil 8-12 micron window emissiv.
      zvegalb   = soilvegpar(jlon,jlat,37)  ! vegetation albedo
      zsoildalb = soilvegpar(jlon,jlat,23)  ! dry soil albedo
      zsoilwalb = soilvegpar(jlon,jlat,24)  ! wet soil albedo
      zqsmin    = soilvegpar(jlon,jlat, 7)  ! minimum water content of soil
      zqsmax    = soilvegpar(jlon,jlat, 6)  ! maximum water content of soil

      if (fmask(jlon,jlat).le.0.5) then
        zqs = qg1(jlon,jlat)
      else
        zqs = zqsmax
      endif

      zsoildf = (zqsmax-zqs)/(zqsmax-zqsmin)
      zsoilwf = 1.-zsoildf

! Tuning: increase of vegetation fraction only to determine emissivity
! and albedo (due to previous modification of vegetation fraction)

      if (ns1.lt.nst-1) zvegdol = min(1.,zvegdol+0.2)

      zemismap1 = zvegdol*zemveg1+(1.-zvegdol)*(zsoildf*zemsoild1+zsoilwf*zemsoilw1)
      zemismap2 = zvegdol*zemveg2+(1.-zvegdol)*(zsoildf*zemsoild2+zsoilwf*zemsoilw2)
      zalbedmap = zvegdol*zvegalb+(1.-zvegdol)*(zsoildf*zsoildalb+zsoilwf*zsoilwalb)

! Desert

      if (alat(jlon,jlat).le.50..and.alat(jlon,jlat).ge.-50..and.ns1.lt.13.and.nv1.eq.12) then

       zemisdes1 = zsoildf*0.93+zsoilwf*0.96 ! desert value
       zemisdes2 = zsoildf*0.85+zsoilwf*0.88 ! desert value
       zalbeddes = zsoildf*0.38+zsoilwf*0.23 ! desert value

       zcoslat = cos( (alat0+(jlat-1)*dlat) *pi/180.)
       nni = nint(50000./(ddx*zcoslat))
       nni = max(nni,1)
       nnj = nint(50000./ddy)
       nnj = max(nnj,1)

! Identification of desert

         i = jlon
         i1 = max(1,i-nni+1)      ! boundary
         i1 = min(i1,nlon-nni+1)  ! boundary
         do jj=1,nnj
            ip = 0

          do ii=i1,i1+nni-1
          if (int(vegeta(ii,jlat,1)).eq.12) ip = ip+1
          enddo

         if (ip.eq.nni) then
         zemismap1 = zemisdes1
         zemismap2 = zemisdes2
         zalbedmap = zalbeddes
         goto 561
         endif
         i1 = i1+1
         i1 = min(i1,nlon-nni+1)  ! boundary
         enddo

         i = jlat
         i1 = max(1,i-nnj+1)      ! boundary
         i1 = min(i1,nlat-nnj+1)  ! boundary

         do jj=1,nni
         ip = 0
          do ii=i1,i1+nnj-1
          if (int(vegeta(jlon,ii,1)).eq.12) ip = ip+1
          enddo

         if (ip.eq.nnj) then
         zemismap1 = zemisdes1
         zemismap2 = zemisdes2
         zalbedmap = zalbeddes
         goto 561
         endif
         i1 = i1+1
         i1 = min(i1,nlat-nnj+1)  ! boundary
         enddo
       endif
 561   continue

! Urban effect for radiative properties of land surface: albedo is decreased
! proportionally to urban vegetation that defines fraction urban areas

       zalbedmap = zalbedmap*(1.-vegeta(jlon,jlat,iurbveg+1))+  &
                   zalbedmap*vegeta(jlon,jlat,iurbveg+1)*0.6

! Definition of albedo and emissivity values as a function of land-sea fraction
! (FMASK) and ice fraction (FICE).
! 0.07 is the value of water body albedo; 0.55 of sea ice albedo.
! 0.990 is the value of water body emissivity; 0.997 of sea ice emissivity.

      if (fmask(jlon,jlat).le.0.5) then ! land
        albedo(jlon,jlat) = zalbedmap
        emismap1(jlon,jlat) = zemismap1
        emismap2(jlon,jlat) = zemismap2
      else ! sea
        if (iflag.eq.1) then
        albedo(jlon,jlat) = (1.-fice(jlon,jlat))*0.07 + fice(jlon,jlat)*0.55
        emismap1(jlon,jlat) = (1.-fice(jlon,jlat))*0.990 + fice(jlon,jlat)*0.997
        emismap2(jlon,jlat) = (1.-fice(jlon,jlat))*0.990 + fice(jlon,jlat)*0.997
        endif
      endif

      enddo
      enddo

      return
      end
!=======================================================================
      subroutine deflai(ylai,vegdol,zndayc)

! Defines vegetation parameters (lai and veg. fraction) that are
! function of day of the year zndayc - data obtained by Guido Fioravanti
! and subsequently revised

      parameter (nvt=14)

      dimension ylai(nvt),vegdol(nvt)

      real datlai(13,nvt),datlai1(12,nvt)
      real datfrv(13,nvt),datfrv1(12,nvt)
      real xxday(13),aavlai(12,nvt),bbvlai(12,nvt),ccvlai(12,nvt), ddvlai(12,nvt),f2vlai(nvt),alvlai(nvt),bevlai(nvt)
      real zindm(12),zlai(13),zfrv(13)

      data xxday/16.,45.,75.,105.,136.,166.,197.,228.,258.,289.,319.,350.,381./

      data datlai1/8.8, 9.2, 9.8,10.1,10.4,10.8    & ! ev.green needl.
                 ,10.5,10.2,10.1, 9.8, 9.2, 8.8,   &
                  12.0,13.0,13.0,14.0,14.0,14.0    & ! ev.green broadl.
                 ,14.0,13.0,13.0,13.0,12.0,12.0,   &
                   0.2, 0.2, 0.2, 0.5, 2.0,10.0    & ! decid. needl.
                 ,10.4,10.0, 8.0, 5.0, 1.0, 0.2,   &
                   0.0, 0.0, 0.2, 2.1, 4.5, 6.8    & ! decid. broadl.
                 , 7.2, 6.5, 5.0, 2.2, 0.5, 0.0,   &
                   1.0, 1.0, 1.5, 6.1, 7.4, 8.8    & ! mixed
                 , 8.8, 8.4, 7.6, 6.0, 5.0, 1.0,   &
                   4.5, 4.5, 4.7, 5.0, 5.5, 5.5    & ! woodland
                 , 5.4, 5.3, 5.1, 5.0, 4.7, 4.5,   &
                   3.0, 3.0, 3.5, 4.0, 5.0, 5.0    & ! wood. grass.
                 , 4.5, 4.0, 3.8, 3.5, 3.3, 3.0,   &
                   0.0, 0.0, 0.2, 0.6, 0.9, 1.7    & ! clos. shrub.
                 , 2.5, 2.5, 1.7, 0.9, 0.3, 0.0,   &
                   0.0, 0.0, 0.3, 0.2, 0.2, 0.3    & ! open. shrub.
                 , 0.4, 0.8, 1.2, 0.8, 0.3, 0.0,   &
                   0.0, 0.0, 0.5, 1.1, 1.7, 3.7    & ! grassland
                 , 4.8, 4.2, 2.0, 1.2, 0.5, 0.0,   &
                   0.0, 0.0, 1.0, 1.1, 1.8, 3.7    & ! cropland
                 , 4.8, 4.2, 2.0, 1.2, 0.5, 0.0,   &
                   0.0, 0.0, 0.0, 0.0, 0.0, 0.0    & ! bare soil
                 , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   &
                   0.3, 0.3, 0.6, 1.8, 2.5, 4.1    & ! urban
                 , 5.0, 4.6, 2.8, 1.9, 1.0, 0.3,   &
                   0.0, 0.0, 0.0, 0.0, 0.0, 0.0    & ! water body
                 , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      data datfrv1/0.8, 0.8, 0.8, 0.8, 0.8, 0.8    & ! ev.green needl.
                 , 0.8, 0.8, 0.8, 0.8, 0.8, 0.8,   &
                   1.0, 1.0, 1.0, 1.0, 1.0, 1.0    & ! ev.green broadl.
                 , 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,   &
                   0.2, 0.2, 0.2, 0.5, 0.7, 0.8    & ! decid. needl.
                 , 0.8, 0.8, 0.8, 0.6, 0.3, 0.2,   &
                   0.1, 0.1, 0.1, 0.4, 0.7, 0.8    & ! decid. broadl.
                 , 0.9, 0.9, 0.8, 0.6, 0.2, 0.1,   &
                   0.1, 0.1, 0.2, 0.4, 0.7, 0.8    & ! mixed
                 , 0.8, 0.8, 0.8, 0.7, 0.3, 0.1,   &
                   0.7, 0.7, 0.7, 0.8, 0.9, 0.9    & ! woodland
                 , 0.8, 0.7, 0.7, 0.7, 0.7, 0.7,   &
                   0.6, 0.6, 0.7, 0.9, 0.9, 0.9    & ! wood. grass.
                 , 0.8, 0.7, 0.7, 0.7, 0.6, 0.6,   &
                   0.2, 0.2, 0.2, 0.5, 0.6, 0.7    & ! clos. shrub.
                 , 0.7, 0.7, 0.6, 0.5, 0.4, 0.2,   &
                   0.1, 0.1, 0.1, 0.4, 0.5, 0.5    & ! open. shrub.
                 , 0.5, 0.5, 0.5, 0.3, 0.2, 0.1,   &
                   0.4, 0.4, 0.4, 0.6, 0.7, 0.8    & ! grassland
                 , 0.8, 0.8, 0.7, 0.6, 0.5, 0.4,   &
                   0.0, 0.1, 0.5, 0.8, 0.9, 0.9    & ! cropland
                 , 0.8, 0.5, 0.3, 0.1, 0.1, 0.0,   &
                   0.0, 0.0, 0.0, 0.0, 0.0, 0.0    & ! bare soil
                 , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   &
                   0.1, 0.1, 0.1, 0.2, 0.2, 0.2    & ! urban
                 , 0.2, 0.2, 0.2, 0.2, 0.1, 0.1,   &
                   0.0, 0.0, 0.0, 0.0, 0.0, 0.0    & ! water body
                 , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      do ivtype=1,nvt
      do iii=1,12
      datlai(iii,ivtype) = datlai1(iii,ivtype)
      datfrv(iii,ivtype) = datfrv1(iii,ivtype)
      enddo
      datlai(13,ivtype) = datlai(1,ivtype)
      datfrv(13,ivtype) = datfrv(1,ivtype)
      end do

! Spline building for vegetation lai which depends on day of year

      do ivtype=1,nvt
      znday = zndayc
      znday = min(znday,365.)

      do i=1,13
      zlai(i) = datlai(i,ivtype)
      zfrv(i) = datfrv(i,ivtype)
      end do

      call yinterp (0.7,0.7,0.7,13,xxday,zlai,znday,zzlai)
      call yinterp (0.7,0.7,0.7,13,xxday,zfrv,znday,zzfrv)

      zzlai = zzlai*0.8 ! arbitrary reduction of lai
      ylai(ivtype) = max(zzlai,0.)
      vegdol(ivtype) = max(zzfrv,0.)

      end do      ! ivtype=1,nvt

      return
      end
!=======================================================================
      subroutine yinterp(alf,ex1,ex2,npi,xi,g,x,f)

!  1-D interpolation with splines with tension of irregularly distributed data
!  Also extrap. out of the interval.
!  Returns one value for each call (scalar version)

!  Input:
!  alf: spline tension param. (0 <= 1) - 0: pure spline interp.; 1: pure linear interp.
!  ex1: param. to extrap. values for x < xi(1);
!  if =0, constant values extrap. equal to g(1); if =1, linear extrap.
!  intermediate values are possible.
!  ex2: as ex1, but to extrap. values for x > xi(npi).
!  npi: no. of points where values of g are known.
!  xi(npi): vector with values of the indep. variable where values of g are known.
!  They MUST be in GROWING order.
!  g(npi): vector with known values to be interpolated
!  x: values of the indep. variable for which the the interpolated value f(x) is computed

!  Output:
!  f: scalar containing the output value f(x) to be computed

      dimension g(npi), xi(npi)

      if(alf.lt..0.or.alf.gt.1) then
      print*, 'Caution: alf out of interval 0-1'
      endif
      if(ex1.lt..0.or.ex1.gt.1) then
      print*, 'Caution: ex1 out of interval 0-1'
      endif
      if(ex2.lt..0.or.ex2.gt.1) then
      print*, 'Caution.: ex2 out of interval 0-1'
      endif

!  2 cases of extrapolation

      if(x.lt.xi(1)) then
      f = g(1) + ex1 * (g(1)-g(2))/(xi(1)-xi(2)) * (x-xi(1))
      return
      elseif (x.ge.xi(npi)) then
      f = g(npi) + ex2 * (g(npi)-g(npi-1))/(xi(npi)-xi(npi-1)) * (x-xi(npi))
      return
      endif

      ir = 0
      do 20  j = 1, npi
      if (x.ge.xi(j)) ir = ir + 1
 20   continue

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

      fm = g(ir)
      xm = xi(ir)
      fp = g(ir+1)
      xp = xi(ir+1)
      delx = xp - xm
      delxp = xpp - xp
      delxm = xm - xmm
      delx1 = x - xm
      delx2 = xp - x
      delxs = delx**2
      delx1s = delx1**2
      delx2s = delx2**2

! Spline contribution

      spl = fm*(delx2/delx + delx1*delx2s/(delxs*delxm) - delx1s*  &
       delx2/((delx+delxp)*delxs)) + fp*(delx1/delx +              &
       delx1s*delx2/(delxs*delxp) - delx1*delx2s/((delx+delxm)*    &
       delxs)) - fmm * delx1*delx2s/((delx+delxm)*delx*delxm) -    &
       fpp * delx1s*delx2/((delx+delxp)*delx*delxp)

! Linear interpolation contribution

      clin = (fm*delx2 + fp*delx1)/delx

! Final value

      f = alf*clin + (1.-alf)*spl

      return
      end
!=======================================================================
      subroutine hd4(f,im,jm,dx,dy,nstep)

! Fourth order smoothing on cartesian coordinates

      real f(im,jm),f2(im,jm)

      zrdy2 = 1.
      zdt = .5/(4.+4./zrdy2)**2

      do jstep=1,nstep
      do j=2,jm-1
      do i=2,im-1
      f2(i,j) = (f(i-1,j)+f(i+1,j)-2.*f(i,j))+(f(i,j+1)+f(i,j-1)-2.*f(i,j))*zrdy2
      enddo
      enddo

      do j=3,jm-2
      do i=3,im-2
      f(i,j) = f(i,j)-zdt*( (f2(i-1,j)+f2(i+1,j)-2.*f2(i,j))+  &
               (f2(i,j+1)+f2(i,j-1)-2.*f2(i,j))*zrdy2 )
      enddo
      enddo
      enddo

      return
      end
!=======================================================================
      subroutine smooth_soil(vt,vout,ni,nj,w,nsmooth)

!  Performs smoothing of matrix VT(NI,NJ) using a 5-point filter (laplacian).
!  It is assumed that the grid distances in x and y are similar
!  The filtering degree is prop. to param. W (between 0 and 1)
!  (for stronger filtering it is advised to increase the number of filtering iterations
!  NSMOOTH, keeping W not larger than about 0.5).
!  Output in VOUT(NI,NJ) and also in VT (caution: VT is redefined!)

      dimension vt(ni,nj), vout(ni,nj)

      do ns=1,nsmooth

      do j=2,nj-1
      do i=2,ni-1
      vout(i,j) = (1.-w)*vt(i,j)+.25*w*(vt(i+1,j)+vt(i-1,j)+vt(i,j+1)+vt(i,j-1))
      enddo
      enddo

!  Border points (1-D filter)

      do i=2, ni-1
      vout(i,1) = (1.-.7*w)*vt(i,1)+.35*w*(vt(i+1,1)+vt(i-1,1))
      vout(i,nj) = (1.-.7*w)*vt(i,nj)+.35*w*(vt(i+1,nj)+vt(i-1,nj))
      enddo
      do j=2, nj-1
      vout(1,j) = (1.-.7*w)*vt(1,j)+.35*w*(vt(1,j+1)+vt(1,j-1))
      vout(ni,j) = (1.-.7*w)*vt(ni,j)+.35*w*(vt(ni,j+1)+vt(ni,j-1))
      enddo

!  Corner points (unchanged)

      vout(1,1) = vt(1,1)
      vout(1,nj) = vt(1,nj)
      vout(ni,nj) = vt(ni,nj)
      vout(ni,1) = vt(ni,1)

      do j=1,nj
      do i=1,ni
      vt(i,j) = vout(i,j)
      enddo
      enddo

      enddo
      return
      end
!=======================================================================
