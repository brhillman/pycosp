! (c) British Crown Copyright 2008, the Met Office.
! All rights reserved.
! $Revision: 80 $, $Date: 2013-09-12 09:14:54 -0700 (Thu, 12 Sep 2013) $
! $URL: https://cfmip-obs-sim.googlecode.com/svn/devel/trunk/cosp_io.F90 $
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
! Jul 2008 - A. Bodas-Salcedo - Initial version
! Oct 2008 - S. Bony - In nc_write_cosp_1d and nc_write_cosp_2d :
!                      the label of layered cloud fractions was wrong -> corrected
!                      (before: low was actually mid, mid was high, high was total,
!                      total was low)
! Sep 2009 - A. Bodas-Salcedo - CMIP5 variable names implemented
!
! Jan 2013 - G. Cesana - Add new phase variables for outputs

MODULE MOD_COSP_IO
  USE MOD_COSP_CONSTANTS
  USE MOD_COSP_TYPES
  use MOD_COSP_Modis_Simulator
  IMPLICIT NONE

contains

! BRH added routine; do not read cosp_output_nl
subroutine construct_cosp_config(cfg)
    use mod_cosp_types
    implicit none

    type(cosp_config), intent(out) :: cfg
    
    ! set flags
    cfg%Lradar_sim = .true.
    cfg%Lisccp_sim = .true.
    cfg%Lmisr_sim = .true.

    ! not running these simulators
    cfg%Llidar_sim = .false.
    cfg%Lmodis_sim = .false.

    ! cloudsat flags
    cfg%LcfadDbze94 = .true.
    cfg%Lclcalipso2 = .false.
    cfg%Lcltlidarradar = .false.
    cfg%Ldbze94 = .false.

    ! lidar flags
    cfg%Latb532 = .false.
    cfg%LcfadLidarsr532 = .false.
    cfg%Lclcalipso2      = .false. ! Needs radar & lidar
    cfg%Lclcalipso       = .false.
    cfg%Lclhcalipso      = .false.
    cfg%Lcllcalipso      = .false.
    cfg%Lclmcalipso      = .false.
    cfg%Lcltcalipso      = .false.
    cfg%Lcltlidarradar   = .false.
    cfg%LparasolRefl    = .false.
    cfg%LlidarBetaMol532 = .false.
    cfg%Lcltlidarradar = .false. ! Needs radar & lidar
    cfg%Lclcalipsoliq    = .false.
    cfg%Lclcalipsoice    = .false.
    cfg%Lclcalipsoun    = .false.
    cfg%Lclcalipsotmp    = .false.
    cfg%Lclcalipsotmpun    = .false.
    cfg%Lclcalipsotmpliq    = .false.
    cfg%Lclcalipsotmpice    = .false.
    cfg%Lclhcalipsoliq      = .false.
    cfg%Lcllcalipsoliq      = .false.
    cfg%Lclmcalipsoliq     = .false.
    cfg%Lcltcalipsoliq      = .false.
    cfg%Lclhcalipsoice      = .false.
    cfg%Lcllcalipsoice      = .false.
    cfg%Lclmcalipsoice      = .false.
    cfg%Lcltcalipsoice      = .false.
    cfg%Lclhcalipsoun      = .false.
    cfg%Lcllcalipsoun      = .false.
    cfg%Lclmcalipsoun      = .false.
    cfg%Lcltcalipsoun      = .false.

    ! isccp flags
    cfg%Lclisccp        = .true.
    cfg%Lalbisccp       = .false.
    cfg%Lboxptopisccp   = .false.
    cfg%Lboxtauisccp    = .false.
    cfg%Lpctisccp       = .false.
    cfg%Ltauisccp       = .false.
    cfg%Lcltisccp       = .false.
    cfg%Lmeantbisccp    = .false.
    cfg%Lmeantbclrisccp = .false.
    cfg%Lfracout = .false.

    ! misr flags
    cfg%LclMISR = .true.

    ! rttov flags
    cfg%Ltbrttov = .false.

    ! modis flags
    cfg%Lcltmodis=.false.
    cfg%Lclwmodis=.false.
    cfg%Lclimodis=.false.
    cfg%Lclhmodis=.false.
    cfg%Lclmmodis=.false.
    cfg%Lcllmodis=.false.
    cfg%Ltautmodis=.false.
    cfg%Ltauwmodis=.false.
    cfg%Ltauimodis=.false.
    cfg%Ltautlogmodis=.false.
    cfg%Ltauwlogmodis=.false.
    cfg%Ltauilogmodis=.false.
    cfg%Lreffclwmodis=.false.
    cfg%Lreffclimodis=.false.
    cfg%Lpctmodis=.false.
    cfg%Llwpmodis=.false.
    cfg%Liwpmodis=.false.
    cfg%Lclmodis=.false.

    ! stats flags
    cfg%Lstats = .true.

    ! Flag to control output to file; not using cmor to write outputs
    cfg%Lwrite_output = .false.

    ! initialize output list
    ! not needed; not using cmor to write outputs

end subroutine


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------- SUBROUTINE READ_COSP_OUTPUT_NL -------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 SUBROUTINE READ_COSP_OUTPUT_NL(cosp_nl,cfg)
  character(len=*),intent(in) :: cosp_nl
  type(cosp_config),intent(out) :: cfg
  ! Local variables
  integer :: i
  logical :: Lradar_sim,Llidar_sim,Lisccp_sim,Lmodis_sim,Lmisr_sim,Lrttov_sim, &
             Lalbisccp,Latb532,Lboxptopisccp,Lboxtauisccp,LcfadDbze94, &
             LcfadLidarsr532,Lclcalipso2,Lclcalipso,Lclhcalipso,Lclisccp,Lcllcalipso, &
             Lclmcalipso,Lcltcalipso,Lcltlidarradar,Lpctisccp,Ldbze94,Ltauisccp,Lcltisccp, &
             Lclcalipsoliq,Lclcalipsoice,Lclcalipsoun, &
             Lclcalipsotmp,Lclcalipsotmpliq,Lclcalipsotmpice,Lclcalipsotmpun, &
             Lcltcalipsoliq,Lcltcalipsoice,Lcltcalipsoun, &
             Lclhcalipsoliq,Lclhcalipsoice,Lclhcalipsoun, &
             Lclmcalipsoliq,Lclmcalipsoice,Lclmcalipsoun, &
             Lcllcalipsoliq,Lcllcalipsoice,Lcllcalipsoun, &
             Ltoffset,LparasolRefl,LclMISR,Lmeantbisccp,Lmeantbclrisccp, &
             Lfracout,LlidarBetaMol532,Ltbrttov, &
             Lcltmodis,Lclwmodis,Lclimodis,Lclhmodis,Lclmmodis,Lcllmodis,Ltautmodis,Ltauwmodis,Ltauimodis,Ltautlogmodis, &
             Ltauwlogmodis,Ltauilogmodis,Lreffclwmodis,Lreffclimodis,Lpctmodis,Llwpmodis, &
             Liwpmodis,Lclmodis

  namelist/COSP_OUTPUT/Lradar_sim,Llidar_sim,Lisccp_sim,Lmodis_sim,Lmisr_sim,Lrttov_sim, &
             Lalbisccp,Latb532,Lboxptopisccp,Lboxtauisccp,LcfadDbze94, &
             LcfadLidarsr532,Lclcalipso2,Lclcalipso,Lclhcalipso,Lclisccp, &
             Lcllcalipso,Lclmcalipso,Lcltcalipso,Lcltlidarradar,Lpctisccp,Ldbze94,Ltauisccp, &
             Lclcalipsoliq,Lclcalipsoice,Lclcalipsoun, &
             Lclcalipsotmp,Lclcalipsotmpliq,Lclcalipsotmpice,Lclcalipsotmpun, &
             Lcltcalipsoliq,Lcltcalipsoice,Lcltcalipsoun, &
             Lclhcalipsoliq,Lclhcalipsoice,Lclhcalipsoun, &
             Lclmcalipsoliq,Lclmcalipsoice,Lclmcalipsoun, &
             Lcllcalipsoliq,Lcllcalipsoice,Lcllcalipsoun, &
             Lcltisccp,Ltoffset,LparasolRefl,LclMISR,Lmeantbisccp,Lmeantbclrisccp, &
             Lfracout,LlidarBetaMol532,Ltbrttov, &
             Lcltmodis,Lclwmodis,Lclimodis,Lclhmodis,Lclmmodis,Lcllmodis,Ltautmodis,Ltauwmodis,Ltauimodis,Ltautlogmodis, &
             Ltauwlogmodis,Ltauilogmodis,Lreffclwmodis,Lreffclimodis,Lpctmodis,Llwpmodis, &
             Liwpmodis,Lclmodis

  do i=1,N_OUT_LIST
    cfg%out_list(i)=''
  enddo
  open(10,file=cosp_nl,status='old')
  read(10,nml=cosp_output)
  close(10)
  
  ! Deal with dependencies
  if (.not.Lradar_sim) then
    LcfadDbze94   = .false.
    Lclcalipso2    = .false.
    Lcltlidarradar = .false. ! Needs radar & lidar
    Ldbze94        = .false.
    Lclcalipso2    = .false. ! Needs radar & lidar
  endif
  if (.not.Llidar_sim) then
    Latb532 = .false.
    LcfadLidarsr532 = .false.
    Lclcalipso2      = .false. ! Needs radar & lidar
    Lclcalipso       = .false.
    Lclhcalipso      = .false.
    Lcllcalipso      = .false.
    Lclmcalipso      = .false.
    Lcltcalipso      = .false.
    Lcltlidarradar   = .false.
    LparasolRefl    = .false.
    LlidarBetaMol532 = .false.
    Lcltlidarradar = .false. ! Needs radar & lidar
    Lclcalipsoliq    = .false.
    Lclcalipsoice    = .false.
    Lclcalipsoun    = .false.
    Lclcalipsotmp    = .false.
    Lclcalipsotmpun    = .false.
    Lclcalipsotmpliq    = .false.
    Lclcalipsotmpice    = .false.
    Lclhcalipsoliq      = .false.
    Lcllcalipsoliq      = .false.
    Lclmcalipsoliq     = .false.
    Lcltcalipsoliq      = .false.
    Lclhcalipsoice      = .false.
    Lcllcalipsoice      = .false.
    Lclmcalipsoice      = .false.
    Lcltcalipsoice      = .false.
    Lclhcalipsoun      = .false.
    Lcllcalipsoun      = .false.
    Lclmcalipsoun      = .false.
    Lcltcalipsoun      = .false.
  endif
  if (.not.Lisccp_sim) then
    Lalbisccp       = .false.
    Lboxptopisccp   = .false.
    Lboxtauisccp    = .false.
    Lclisccp        = .false.
    Lpctisccp       = .false.
    Ltauisccp       = .false.
    Lcltisccp       = .false.
    Lmeantbisccp    = .false.
    Lmeantbclrisccp = .false.
  endif
  if (.not.Lmisr_sim) then
    LclMISR = .false.
  endif
  if (.not.Lrttov_sim) then
    Ltbrttov = .false.
  endif
  if ((.not.Lradar_sim).and.(.not.Llidar_sim).and. &
      (.not.Lisccp_sim).and.(.not.Lmisr_sim)) then
    Lfracout = .false.
  endif
  if (.not.Lmodis_sim) then
    Lcltmodis=.false.
    Lclwmodis=.false.
    Lclimodis=.false.
    Lclhmodis=.false.
    Lclmmodis=.false.
    Lcllmodis=.false.
    Ltautmodis=.false.
    Ltauwmodis=.false.
    Ltauimodis=.false.
    Ltautlogmodis=.false.
    Ltauwlogmodis=.false.
    Ltauilogmodis=.false.
    Lreffclwmodis=.false.
    Lreffclimodis=.false.
    Lpctmodis=.false.
    Llwpmodis=.false.
    Liwpmodis=.false.
    Lclmodis=.false.
  endif
  if (Lmodis_sim) Lisccp_sim = .true.
  
  cfg%Lstats = .false.
  if ((Lradar_sim).or.(Llidar_sim).or.(Lisccp_sim)) cfg%Lstats = .true.
  
  ! Copy instrument flags to cfg structure
  cfg%Lradar_sim = Lradar_sim
  cfg%Llidar_sim = Llidar_sim
  cfg%Lisccp_sim = Lisccp_sim
  cfg%Lmodis_sim = Lmodis_sim
  cfg%Lmisr_sim  = Lmisr_sim
  cfg%Lrttov_sim = Lrttov_sim
  
  ! Flag to control output to file
  cfg%Lwrite_output = .false.
  if (cfg%Lstats.or.cfg%Lmisr_sim.or.cfg%Lrttov_sim) then
    cfg%Lwrite_output = .true.
  endif
  
  ! Output diagnostics
  i = 1
  if (Lalbisccp)        cfg%out_list(i) = 'albisccp'
  i = i+1
  if (Latb532)          cfg%out_list(i) = 'atb532'
  i = i+1
  if (Lboxptopisccp)    cfg%out_list(i) = 'boxptopisccp'
  i = i+1
  if (Lboxtauisccp)     cfg%out_list(i) = 'boxtauisccp'
  i = i+1
  if (LcfadDbze94)      cfg%out_list(i) = 'cfadDbze94'
  i = i+1
  if (LcfadLidarsr532)  cfg%out_list(i) = 'cfadLidarsr532'
  i = i+1
  if (Lclcalipso2)      cfg%out_list(i) = 'clcalipso2'
  i = i+1
  if (Lclcalipso)       cfg%out_list(i) = 'clcalipso'
  i = i+1
  if (Lclhcalipso)      cfg%out_list(i) = 'clhcalipso'
  i = i+1
  if (Lclisccp)         cfg%out_list(i) = 'clisccp'
  i = i+1
  if (Lcllcalipso)      cfg%out_list(i) = 'cllcalipso'
  i = i+1
  if (Lclmcalipso)      cfg%out_list(i) = 'clmcalipso'
  i = i+1
  if (Lcltcalipso)      cfg%out_list(i) = 'cltcalipso'
  i = i+1

  if (Lcllcalipsoice)      cfg%out_list(i) = 'cllcalipsoice'
  i = i+1
  if (Lclmcalipsoice)      cfg%out_list(i) = 'clmcalipsoice'
  i = i+1
  if (Lclhcalipsoice)      cfg%out_list(i) = 'clhcalipsoice'
  i = i+1
  if (Lcltcalipsoice)      cfg%out_list(i) = 'cltcalipsoice'
  i = i+1
  if (Lcllcalipsoliq)      cfg%out_list(i) = 'cllcalipsoliq'
  i = i+1
  if (Lclmcalipsoliq)      cfg%out_list(i) = 'clmcalipsoliq'
  i = i+1
  if (Lclhcalipsoliq)      cfg%out_list(i) = 'clhcalipsoliq'
  i = i+1
  if (Lcltcalipsoliq)      cfg%out_list(i) = 'cltcalipsoliq'
  i = i+1
  if (Lcllcalipsoun)      cfg%out_list(i) = 'cllcalipsoun'
  i = i+1
  if (Lclmcalipsoun)      cfg%out_list(i) = 'clmcalipsoun'
  i = i+1
  if (Lclhcalipsoun)      cfg%out_list(i) = 'clhcalipsoun'
  i = i+1
  if (Lcltcalipsoun)      cfg%out_list(i) = 'cltcalipsoun'
  i = i+1

  if (Lclcalipsoice)       cfg%out_list(i) = 'clcalipsoice'
  i = i+1
  if (Lclcalipsoliq)       cfg%out_list(i) = 'clcalipsoliq'
  i = i+1
  if (Lclcalipsoun)       cfg%out_list(i) = 'clcalipsoun'
  i = i+1

  if (Lclcalipsotmp)       cfg%out_list(i) = 'clcalipsotmp'
  i = i+1
  if (Lclcalipsotmpice)       cfg%out_list(i) = 'clcalipsotmpice'
  i = i+1
  if (Lclcalipsotmpliq)       cfg%out_list(i) = 'clcalipsotmpliq'
  i = i+1
  if (Lclcalipsotmpun)       cfg%out_list(i) = 'clcalipsotmpun'
  i = i+1
  if (Lcltlidarradar)   cfg%out_list(i) = 'cltlidarradar'
  i = i+1
  if (Lpctisccp)        cfg%out_list(i) = 'pctisccp'
  i = i+1
  if (Ldbze94)          cfg%out_list(i) = 'dbze94'
  i = i+1
  if (Ltauisccp)        cfg%out_list(i) = 'tauisccp'
  i = i+1
  if (Lcltisccp)        cfg%out_list(i) = 'cltisccp'
  i = i+1
  if (Ltoffset)         cfg%out_list(i) = 'toffset'
  i = i+1
  if (LparasolRefl)     cfg%out_list(i) = 'parasolRefl'
  i = i+1
  if (LclMISR)          cfg%out_list(i) = 'clMISR'
  i = i+1
  if (Lmeantbisccp)     cfg%out_list(i) = 'meantbisccp'
  i = i+1
  if (Lmeantbclrisccp)  cfg%out_list(i) = 'meantbclrisccp'
  i = i+1
  if (Lfracout)         cfg%out_list(i) = 'fracout'
  i = i+1
  if (LlidarBetaMol532) cfg%out_list(i) = 'lidarBetaMol532'
  i = i+1
  if (Ltbrttov)         cfg%out_list(i) = 'tbrttov'
  i = i+1
  if (Lcltmodis)        cfg%out_list(i) = 'cltmodis'
  i = i+1
  if (Lclwmodis)        cfg%out_list(i) = 'clwmodis'
  i = i+1
  if (Lclimodis)        cfg%out_list(i) = 'climodis'
  i = i+1
  if (Lclhmodis)        cfg%out_list(i) = 'clhmodis'
  i = i+1
  if (Lclmmodis)        cfg%out_list(i) = 'clmmodis'
  i = i+1
  if (Lcllmodis)        cfg%out_list(i) = 'cllmodis'
  i = i+1
  if (Ltautmodis)       cfg%out_list(i) = 'tautmodis'
  i = i+1
  if (Ltauwmodis)       cfg%out_list(i) = 'tauwmodis'
  i = i+1
  if (Ltauimodis)       cfg%out_list(i) = 'tauimodis'
  i = i+1
  if (Ltautlogmodis)    cfg%out_list(i) = 'tautlogmodis'
  i = i+1
  if (Ltauwlogmodis)    cfg%out_list(i) = 'tauwlogmodis'
  i = i+1
  if (Ltauilogmodis)    cfg%out_list(i) = 'tauilogmodis'
  i = i+1
  if (Lreffclwmodis)    cfg%out_list(i) = 'reffclwmodis'
  i = i+1
  if (Lreffclimodis)    cfg%out_list(i) = 'reffclimodis'
  i = i+1
  if (Lpctmodis)        cfg%out_list(i) = 'pctmodis'
  i = i+1
  if (Llwpmodis)        cfg%out_list(i) = 'lwpmodis'
  i = i+1
  if (Liwpmodis)        cfg%out_list(i) = 'iwpmodis'
  i = i+1
  if (Lclmodis)         cfg%out_list(i) = 'clmodis'
    
  if (i /= N_OUT_LIST) then
     print *, 'COSP_IO: wrong number of output diagnostics'
     print *, i,N_OUT_LIST
     stop
  endif

  ! Copy diagnostic flags to cfg structure
  ! ISCCP simulator  
  cfg%Lalbisccp = Lalbisccp
  cfg%Latb532 = Latb532
  cfg%Lboxptopisccp = Lboxptopisccp
  cfg%Lboxtauisccp = Lboxtauisccp
  cfg%Lmeantbisccp = Lmeantbisccp
  cfg%Lmeantbclrisccp = Lmeantbclrisccp
  cfg%Lclisccp = Lclisccp
  cfg%Lpctisccp = Lpctisccp
  cfg%Ltauisccp = Ltauisccp
  cfg%Lcltisccp = Lcltisccp
  ! CloudSat simulator  
  cfg%Ldbze94 = Ldbze94
  cfg%LcfadDbze94 = LcfadDbze94
  ! CALIPSO/PARASOL simulator  
  cfg%LcfadLidarsr532 = LcfadLidarsr532
  cfg%Lclcalipso2 = Lclcalipso2
  cfg%Lclcalipso = Lclcalipso
  cfg%Lclhcalipso = Lclhcalipso
  cfg%Lcllcalipso = Lcllcalipso
  cfg%Lclmcalipso = Lclmcalipso
  cfg%Lcltcalipso = Lcltcalipso
  cfg%Lclhcalipsoice = Lclhcalipsoice
  cfg%Lcllcalipsoice = Lcllcalipsoice
  cfg%Lclmcalipsoice = Lclmcalipsoice
  cfg%Lcltcalipsoice = Lcltcalipsoice
  cfg%Lclhcalipsoliq = Lclhcalipsoliq
  cfg%Lcllcalipsoliq = Lcllcalipsoliq
  cfg%Lclmcalipsoliq = Lclmcalipsoliq
  cfg%Lcltcalipsoliq = Lcltcalipsoliq
  cfg%Lclhcalipsoun = Lclhcalipsoun
  cfg%Lcllcalipsoun = Lcllcalipsoun
  cfg%Lclmcalipsoun = Lclmcalipsoun
  cfg%Lcltcalipsoun = Lcltcalipsoun
  cfg%Lclcalipsoice = Lclcalipsoice
  cfg%Lclcalipsoliq = Lclcalipsoliq
  cfg%Lclcalipsoun = Lclcalipsoun
  cfg%Lclcalipsotmp = Lclcalipsotmp
  cfg%Lclcalipsotmpice = Lclcalipsotmpice
  cfg%Lclcalipsotmpliq = Lclcalipsotmpliq
  cfg%Lclcalipsotmpun = Lclcalipsotmpun
  cfg%Lcltlidarradar = Lcltlidarradar
  cfg%LparasolRefl = LparasolRefl
  ! MISR simulator  
  cfg%LclMISR = LclMISR
  ! Other
  cfg%Ltoffset = Ltoffset
  cfg%Lfracout = Lfracout
  cfg%LlidarBetaMol532 = LlidarBetaMol532
  ! RTTOV
  cfg%Ltbrttov = Ltbrttov
  ! MODIS simulator  
  cfg%Lcltmodis=Lcltmodis
  cfg%Lclwmodis=Lclwmodis
  cfg%Lclimodis=Lclimodis
  cfg%Lclhmodis=Lclhmodis
  cfg%Lclmmodis=Lclmmodis
  cfg%Lcllmodis=Lcllmodis
  cfg%Ltautmodis=Ltautmodis
  cfg%Ltauwmodis=Ltauwmodis
  cfg%Ltauimodis=Ltauimodis
  cfg%Ltautlogmodis=Ltautlogmodis
  cfg%Ltauwlogmodis=Ltauwlogmodis
  cfg%Ltauilogmodis=Ltauilogmodis
  cfg%Lreffclwmodis=Lreffclwmodis
  cfg%Lreffclimodis=Lreffclimodis
  cfg%Lpctmodis=Lpctmodis
  cfg%Llwpmodis=Llwpmodis
  cfg%Liwpmodis=Liwpmodis
  cfg%Lclmodis=Lclmodis
 END SUBROUTINE READ_COSP_OUTPUT_NL

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------- SUBROUTINE ERROR_CONTROL ------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_ERROR(routine_name,message,errcode) 
    character(len = *), intent(in) :: routine_name
    character(len = *), intent(in) :: message
    integer,optional :: errcode
    
    write(6, *) " ********** Failure in ", trim(routine_name)
    write(6, *) " ********** ", trim(message)
    if (present(errcode)) write(6, *) " ********** errcode: ", errcode
    flush(6)
    stop
  END SUBROUTINE COSP_ERROR

END MODULE MOD_COSP_IO
