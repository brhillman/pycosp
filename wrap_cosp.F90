! (c) British Crown Copyright 2008, the Met Office.
! All rights reserved.
! $Revision: 80 $, $Date: 2013-09-12 09:14:54 -0700 (Thu, 12 Sep 2013) $
! $URL: https://cfmip-obs-sim.googlecode.com/svn/devel/trunk/cosp_test.F90 $
! 
! Redistribution and use in source and binary forms, with or without modification, are permitted 
! provided that the following conditions are met:
! 
!     * Redistributions of source code must retain the above copyright notice, this list 
!       of conditions and the following disclaimer.
!     * Redistributions in binary form must reproduce the above copyright notice, this list
!       of conditions and the following disclaimer in the documentation and/or other materials 
!       provided with the distribution.
!     * Neither the name of the Met Office nor the names of its contributors may be used 
!       to endorse or promote products derived from this software without specific prior written 
!       permission.
! 
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR 
! IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND 
! FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR 
! CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER 
! IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT 
! OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

!
! History:
! Feb 2008 - A. Bodas-Salcedo - Initial version
! Dec 2010 - A. Bodas-Salcedo - Added capability for processing multiple files
! Nov 2015 - B. Hillman - created wrapper from COSP_TEST.F90 to isolate cosp
!            derived types from calling function so that I could wrap for python
!

#include "cosp_defs.h"

! ****************************************************************************
! Inputs
!       npoints : integer
!           Number of gridboxes (hidden in python wrapper)
!       nlevels: integer
!           Number of vertical levels (hidden in python wrapper)
!       ncolumns : integer
!           Number of subcolumns to generate from each gridbox/crm column
!           If points represent large-scale gridboxes
!           then ncolumns should be greater than 1 (somewhere around 50-100?)
!       lat : real
!           Latitude coordinate for each gridbox (npoints)
!       lon : real
!           Longitude coordinate for each gridbox (npoints)
!       skt : real
!           Surface skin temperature in kelvin (npoints)
!       landmask : integer(npoints)
!           Land-sea mask; 0 -> ocean, 1 -> land
!       sunlit : integer(npoints)
!           Sunlit flag; 0 -> nighttime, 1 -> sunlit
!       p : real(npoints, nlevels)
!           Pressure at model level midpoints (Pa?)
!       ph : real(npoints, nlevels)
!           Pressure at model level interfaces (Pa? bottom or top?)
!       zlev : real(npoints, nlevels)
!           Height at model level midpoints (m)
!       zlev_half : real(npoints, nlevels)
!           Height at model level interfaces (m; bottom or top?)
!       t : real(npoints, nlevels)
!           Temperature at model levels (K)
!       sh : real(npoints, nlevels)
!           Specific humidity
! ****************************************************************************
subroutine run_cosp(npoints,nlevels,ncolumns, &
                    lat, lon, skt, landmask, sunlit, &
                    p, ph, zlev, zlev_half, &
                    t, sh, rh, tca, cca, &
                    mr_lsliq, mr_lsice, mr_ccliq, mr_ccice, &
                    mr_lsrain, mr_lssnow, mr_lsgrpl, &
                    mr_ccrain, mr_ccsnow, &
                    dtau_s, dtau_c, dem_s, dem_c, &
                    clmisr, clisccp, cfad_ze, ze_col, ze_non, &
                    verbose, use_radar_lut, use_vgrid)

    use mod_cosp_constants
    use mod_cosp_types
    use mod_cosp
    use mod_cosp_io, only: construct_cosp_config

    IMPLICIT NONE
    !f2py integer, intent(in), hide :: npoints, nlevels
    !f2py integer, intent(in) :: ncolumns
    !f2py integer, intent(in), optional :: verbose=-1
    !f2py integer, intent(in), optional :: use_radar_lut=0
    !f2py logical, intent(in), optional :: use_vgrid=1

    ! TODO: add all of the namelist variables as optional arguments, and
    ! get rid of namelist altogether?
    ! f2py integer, intent(in), optional :: npoints_it=500, nlr=40
    ! f2py logical, intent(in), optional :: use_vgrid=.true., csat_vgrid=.true.
    ! f2py real, optional :: radar_freq=94.0, k2=-1
    ! f2py integer, optional :: surface_radar=0, use_mie_tables=0, use_gas_abs=1, do_ray=0, melt_lay=0
    ! f2py logical, optional :: use_reff=.true., use_precipitation_fluxes=.false.
    ! f2py integer, optional :: nprmts_max_hydro=12, naero=1, nprmts_max_aero=1, lidar_ice_type=0, overlap=3, isccp_topheight=1, isccp_topheight_direction=2
    ! f2py, integer, optional :: platform=1, satellite=15, instrument=0, nchannels=8
    ! f2py, integer, optional :: channels=(/1, 3, 5, 6, 8, 10, 11, 13/)
    ! f2py, real, optional :: surfem=(/0, 0, 0, 0, 0, 0, 0, 0/)
    ! f2py, real, optional :: zenang=50.0, co2=5.241e-4, ch4=9.139e-7, n2o=4.665e-7, co=2.098e-7, emsfc_lw=0.99

    ! Arguments
    integer, intent(in) :: npoints, nlevels, ncolumns
    real, intent(in), dimension(npoints) :: lat, lon, skt, landmask, sunlit
    real, intent(in), dimension(npoints,nlevels) :: p, ph, &
        zlev, zlev_half, t, sh, rh
    real, intent(in), dimension(npoints,nlevels) :: tca, cca, &
        mr_lsliq, mr_lsice, mr_ccliq, mr_ccice, &
        mr_lsrain, mr_lssnow, mr_lsgrpl, mr_ccrain, mr_ccsnow, &
        dtau_s, dtau_c, dem_s, dem_c
    real, intent(out) :: clmisr(npoints,7,16)
    real, intent(out) :: clisccp(npoints,7,7)
    real, intent(out) :: cfad_ze(npoints,15,40)
    real, intent(out) :: ze_col(npoints,nlevels)
    real, intent(out) :: ze_non(npoints,nlevels)
    integer, intent(in) :: verbose
    integer, intent(in) :: use_radar_lut
    logical, intent(in) :: use_vgrid
                    
    ! Local variables
    type(cosp_config) :: cfg   ! Configuration options
    type(cosp_gridbox) :: gbx ! Gridbox information. Input for COSP
    type(cosp_subgrid) :: sgx     ! Subgrid outputs
    type(cosp_sgradar) :: sgradar ! Output from radar simulator
    type(cosp_sglidar) :: sglidar ! Output from lidar simulator
    type(cosp_isccp)   :: isccp   ! Output from ISCCP simulator
    type(cosp_modis)   :: modis   ! Output from MODIS simulator
    type(cosp_misr)    :: misr    ! Output from MISR simulator
    type(cosp_rttov)   :: rttov   ! Output from RTTOV 
    type(cosp_vgrid)   :: vgrid   ! Information on vertical grid of stats
    type(cosp_radarstats) :: stradar ! Summary statistics from radar simulator
    type(cosp_lidarstats) :: stlidar ! Summary statistics from lidar simulator
    double precision :: time, time_bnds(2)

    ! Variables from namelist
    integer :: NPOINTS_IT=500
                ! Max number of gridpoints to be processed in one iteration
    !logical :: USE_VGRID=.true. 
                ! Use fixed vertical grid for outputs? 
                ! (if .true. then you need to define number of levels with Nlr)
    integer :: NLR=40
                ! Number of levels in statistical outputs 
                ! (only used if USE_VGRID=.true.)
    logical :: CSAT_VGRID=.true. 
                ! CloudSat vertical grid? 
                ! (if .true. then CloudSat standard grid is used for the outputs.
                ! USE_VGRID needs also be .true.)
    !------------------------------------------------------------------------------
    !--------------- Inputs related to radar simulations
    !------------------------------------------------------------------------------
    real :: RADAR_FREQ=94.0     ! CloudSat radar frequency (GHz)
    integer :: SURFACE_RADAR=0  ! surface=1, spaceborne=0
    integer :: use_mie_tables=0 ! use a precomputed lookup table? yes=1,no=0
    integer :: use_gas_abs=1    ! include gaseous absorption? yes=1,no=0
    integer :: do_ray=0         ! calculate/output Rayleigh refl=1, not=0
    integer :: melt_lay=0       ! melting layer model off=0, on=1
    real :: k2=-1               ! |K|^2, -1=use frequency dependent default
    logical :: use_reff=.true.  ! True if you want effective radius to be used 
                                ! by radar simulator (always used by lidar)
    logical :: use_precipitation_fluxes=.false.  
                ! True if precipitation fluxes are input to the algorithm 
    !------------------------------------------------------------------------------
    !---------------- Inputs related to lidar simulations
    !------------------------------------------------------------------------------
    integer :: Nprmts_max_hydro=12  ! Max number of parameters 
                                    ! for hydrometeor size distributions
    integer :: Naero=1              ! Number of aerosol species (Not used)
    integer :: Nprmts_max_aero=1    ! Max number of parameters for aerosol 
                                    ! size distributions (Not used)
    integer :: lidar_ice_type=0     ! Ice particle shape in lidar calculations 
                                    ! (0=ice-spheres ; 1=ice-non-spherical)
    integer :: OVERLAP=3            !  overlap type: 1=max, 2=rand, 3=max/rand
    !------------------------------------------------------------------------------
    !---------------- Inputs related to ISCCP simulator
    !------------------------------------------------------------------------------
    integer :: ISCCP_TOPHEIGHT=1  
                !  1 = adjust top height using both a computed
                !  infrared brightness temperature and the visible
                !  optical depth to adjust cloud top pressure. Note
                !  that this calculation is most appropriate to compare
                !  to ISCCP data during sunlit hours.
                !  2 = do not adjust top height, that is cloud top
                !  pressure is the actual cloud top pressure
                !  in the model
                !  3 = adjust top height using only the computed
                !  infrared brightness temperature. Note that this
                !  calculation is most appropriate to compare to ISCCP
                !  IR only algortihm (i.e. you can compare to nighttime
                !  ISCCP data with this option)
    integer :: ISCCP_TOPHEIGHT_DIRECTION=2   
                ! direction for finding atmosphere pressure level
                ! with interpolated temperature equal to the radiance
                ! determined cloud-top temperature
                ! 1 = find the *lowest* altitude (highest pressure) level
                ! with interpolated temperature equal to the radiance
                ! determined cloud-top temperature
                ! 2 = find the *highest* altitude (lowest pressure) level
                ! with interpolated temperature equal to the radiance 
                ! determined cloud-top temperature. This is the 
                ! default value since V4.0 of the ISCCP simulator.
                ! ONLY APPLICABLE IF top_height EQUALS 1 or 3
    !------------------------------------------------------------------------------
    !-------------- RTTOV inputs
    !------------------------------------------------------------------------------
    integer :: Platform=1    ! satellite platform
    integer :: Satellite=15  ! satellite
    integer :: Instrument=0  ! instrument
    integer :: Nchannels=8   ! Number of channels to be computed
    integer, dimension(8) :: Channels=(/1,3,5,6,8,10,11,13/) 
                ! Channel numbers (please be sure that you supply Nchannels)
    real, dimension(8) :: Surfem=(/0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/)
                ! Surface emissivity (please be sure that you supply Nchannels)
    real :: ZenAng=50.0     ! Satellite Zenith Angle
    real :: CO2=5.241e-04   ! Mixing ratios of trace gases
    real :: CH4=9.139e-07
    real :: N2O=4.665e-07
    real :: CO=2.098e-07

    ! set surface emissivity
    real :: emsfc_lw = 0.99

    integer :: i, j
    logical :: load_radar_lut = .false.
    logical :: save_radar_lut = .false.

    !---------------- End of declaration of variables --------------

    ! use radar look-up tables?
    if (use_radar_lut > 0) then
        load_radar_lut = .true.
        save_radar_lut = .true.
    end if
   
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Read COSP namelists; BRH edit: set manually instead
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    call construct_cosp_config(cfg)
  
    time = 8*1.D0/8.D0
    time_bnds = (/time-0.5*3.D0/24.D0,time+0.5*3.D0/24.D0/)

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Allocate memory for gridbox type
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if (verbose > 0) print *, 'Allocating memory for gridbox type...'
    call construct_cosp_gridbox(time,time_bnds,radar_freq,surface_radar,use_mie_tables,use_gas_abs, &
                                do_ray,melt_lay,k2, &
                                Npoints,Nlevels,Ncolumns,N_HYDRO,Nprmts_max_hydro,Naero,Nprmts_max_aero,Npoints_it, &
                                lidar_ice_type,isccp_topheight,isccp_topheight_direction,overlap,emsfc_lw, &
                                use_precipitation_fluxes,use_reff, &
                                Platform,Satellite,Instrument,Nchannels,ZenAng, &
                                channels(1:Nchannels),surfem(1:Nchannels),co2,ch4,n2o,co,gbx,load_lut=load_radar_lut)
                                
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Define new vertical grid
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if (verbose > 0) print *, 'Defining new vertical grid...'
    call construct_cosp_vgrid(gbx,Nlr,use_vgrid,csat_vgrid,vgrid)

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Allocate memory for other types
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if (verbose > 0) print *, 'Allocating memory for other types...'
    call construct_cosp_subgrid(Npoints, Ncolumns, Nlevels, sgx)
    call construct_cosp_sgradar(cfg,Npoints,Ncolumns,Nlevels,N_HYDRO,sgradar)
    call construct_cosp_radarstats(cfg,Npoints,Ncolumns,vgrid%Nlvgrid,N_HYDRO,stradar)
    call construct_cosp_sglidar(cfg,Npoints,Ncolumns,Nlevels,N_HYDRO,PARASOL_NREFL,sglidar)
    call construct_cosp_lidarstats(cfg,Npoints,Ncolumns,vgrid%Nlvgrid,N_HYDRO,PARASOL_NREFL,stlidar)
    call construct_cosp_isccp(cfg,Npoints,Ncolumns,Nlevels,isccp)
    call construct_cosp_modis(cfg,Npoints,modis)
    call construct_cosp_misr(cfg,Npoints,misr)
    call construct_cosp_rttov(cfg,Npoints,Nchannels,rttov)

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Populate input structure
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if (verbose > 0) print *, 'Populating input structure...'
    gbx%longitude(:) = lon(:)
    gbx%latitude(:) = lat(:)
    gbx%p(:,:) = p(:,:)
    gbx%ph(:,:) = ph(:,:)
    gbx%zlev(:,:) = zlev(:,:)
    gbx%zlev_half(:,:) = zlev_half(:,:)
    gbx%T(:,:) = T(:,:)
    gbx%q(:,:) = rh(:,:)
    gbx%sh(:,:) = sh(:,:)
    gbx%psfc(:) = ph(:,1)
    gbx%skt(:)  = skt(:)
    gbx%land(:) = landmask(:)
    gbx%mr_ozone(:,:) = 0
    gbx%u_wind(:) = 0
    gbx%v_wind(:) = 0
    gbx%sunlit(:) = sunlit(:)

    gbx%cca = cca(:,:)
    gbx%tca = tca(:,:)

    gbx%mr_hydro(:,:,I_LSCLIQ) = mr_lsliq(:,:)
    gbx%mr_hydro(:,:,I_LSCICE) = mr_lsice(:,:)
    gbx%mr_hydro(:,:,I_CVCLIQ) = mr_ccliq(:,:)
    gbx%mr_hydro(:,:,I_CVCICE) = mr_ccice(:,:)

    ! use mixing ratios instead of fluxes
    gbx%mr_hydro(:,:,I_LSRAIN) = mr_lsrain(:,:)
    gbx%mr_hydro(:,:,I_LSSNOW) = mr_lssnow(:,:)
    gbx%mr_hydro(:,:,I_LSGRPL) = mr_lsgrpl(:,:)
    gbx%mr_hydro(:,:,I_CVRAIN) = mr_ccrain(:,:)
    gbx%mr_hydro(:,:,I_CVSNOW) = mr_ccsnow(:,:)
    
    ! use defaults?
    !gbx%Reff = Reff
    gbx%Reff(:,:,:) = 0.0

    ! optical properties
    gbx%dtau_s(:,:) = dtau_s(:,:)
    gbx%dtau_c(:,:) = dtau_c(:,:)
    gbx%dem_s(:,:) = dem_s(:,:)
    gbx%dem_c(:,:) = dem_c(:,:)

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Call simulator
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if (verbose > 0) print *, 'Calling simulator...'
    call cosp( &
        overlap,Ncolumns,cfg,vgrid,gbx,sgx, &
        sgradar,sglidar,isccp,misr,modis,stradar,stlidar &
    )

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Copy outputs
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if (verbose > 0) print *, 'Copy outputs...'
    clmisr(:,:,:) = misr%fq_misr(:,:,:)
    clisccp(:,:,:) = isccp%fq_isccp(:,:,:)
    cfad_ze(:,:,:) = stradar%cfad_ze(:,:,:)

    ! note: these aren't quite right for ncolumns /= 1 I don't think...
    ze_col(:,:) = sum(sgradar%ze_tot(:,:,:), dim=2) / ncolumns
    ze_non(:,:) = sum(sgradar%ze_non(:,:,:), dim=2) / ncolumns

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Deallocate memory in derived types
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if (verbose > 0) print *, 'Deallocating memory...'
    call free_cosp_gridbox(gbx, save_lut=save_radar_lut)
    call free_cosp_subgrid(sgx)
    call free_cosp_sgradar(sgradar)
    call free_cosp_radarstats(stradar)
    call free_cosp_sglidar(sglidar)
    call free_cosp_lidarstats(stlidar)
    call free_cosp_isccp(isccp)
    call free_cosp_misr(misr)
    call free_cosp_modis(modis)
    call free_cosp_rttov(rttov)
    call free_cosp_vgrid(vgrid)
   
end subroutine


!***************************************************************************
!
! subroutines to get histogram bin centers and edges
! for COSP summary statistics
!
!***************************************************************************
subroutine get_tau_bins(tau, tau_bnds)
    use mod_cosp_constants
    implicit none
    real, intent(out) :: tau(7), tau_bnds(2, 7)
    tau = ISCCP_TAU
    tau_bnds = ISCCP_TAU_BNDS
end subroutine


subroutine get_plev_bins(plev, plev_bnds)
    use mod_cosp_constants
    implicit none
    real, intent(out) :: plev(7), plev_bnds(2, 7)
    plev = ISCCP_PC
    plev_bnds = ISCCP_PC_BNDS
end subroutine


subroutine get_cth_bins(cth, cth_bnds)
    use mod_cosp_constants
    implicit none
    real, intent(out) :: cth(16), cth_bnds(2, 16)
    cth = MISR_CTH
    cth_bnds = MISR_CTH_BNDS 
end subroutine


subroutine get_dbze_bins(dbze, dbze_bnds)
    use mod_cosp_constants
    implicit none
    real, intent(out) :: dbze(15), dbze_bnds(2, 15)
    integer :: i
    dbze_bnds(1, :) = (/ (CFAD_ZE_MIN + (i-1) * CFAD_ZE_WIDTH, i = 1, DBZE_BINS) /)
    dbze_bnds(2, :) = (/ (CFAD_ZE_MIN + i * CFAD_ZE_WIDTH, i = 1, DBZE_BINS) /)
    dbze = sum(dbze_bnds, dim=1)/2
end subroutine


subroutine get_alt_bins(nlevels, zstep, alt, alt_bnds)
    use mod_cosp_constants
    implicit none

    !f2py integer, intent(in), optional :: nlevels=40
    !f2py real, intent(in), optional :: zstep=480
    integer, intent(in) :: nlevels
    real, intent(in) :: zstep
    real, intent(out) :: alt(40), alt_bnds(2, 40)
    integer :: i
    do i = 1,40
        alt_bnds(1, i) = (i - 1) * zstep
        alt_bnds(2, i) = i * zstep
    end do
    alt = (alt_bnds(1, :) + alt_bnds(2, :)) / 2.0
end subroutine
