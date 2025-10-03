!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
module galwemrad_forcingMod
!BOP
! !MODULE: galwemrad_forcingMod

! !DESCRIPTION:
!  This module contains variables and data structures that are used
!  for the implementation of GALWEM radiation data used as supplemental
!  forcing within LIS. GALWEM radiation replaces the AGRMET radiation in polar
!  stereographic projection at hourly, which ended on 2025/05/16.
!  The standard GALWEM radiation data is at 0.25 degree downloaded from USAF.
!  The native 17km resolution GALWEM radiation data is archived by LIS-USAF team.
!  Current implementation directly uses downward shortwave and longwave fluxes.
!  Alternateve derivation of the flux fields using the cloud information 
!  similar to AGRMET approach maybe developed in future.
!  

! REVISION HISTORY:
! 11 Mar 2022; Yeosang Yoon; Initial Specification
! 08 Sep 2022; Yeosang Yoon, Add codes to read GALWEM 25 DEG dataset
! 11 Jan 2024; Eric Kemp, added third entries for galtime and metdata
!              for temporary storage.
! 13 Jun 2025: Hiroko Beaudoing, adipted GALWEM radiation implenemtation for
!                                radiation supplemental forcing option.
! !USES:
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_galwemrad      !defines the native resolution of the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: galwemrad_struc

  type, public ::  galwemrad_type_dec
     real                              :: ts
     integer                           :: nc, nr, vector_len   
     real*8                            :: fcsttime1, fcsttime2, fcsttime3
     character(len=LIS_CONST_PATH_LEN) :: odir      !GALWEM radiation forcing Directory
     character*20                      :: runmode
     integer                           :: resol     !GALWEM radiation resolution (17km or 25deg)

     integer, allocatable   :: gindex(:,:)
     integer                :: mi

     integer, allocatable   :: n111(:)
     integer, allocatable   :: n121(:)
     integer, allocatable   :: n211(:)
     integer, allocatable   :: n221(:)
     real, allocatable      :: w111(:),w121(:)
     real, allocatable      :: w211(:),w221(:)
     
     integer, allocatable   :: n112(:,:)
     integer, allocatable   :: n122(:,:)
     integer, allocatable   :: n212(:,:)
     integer, allocatable   :: n222(:,:)
     real, allocatable      :: w112(:,:),w122(:,:)
     real, allocatable      :: w212(:,:),w222(:,:)

     integer, allocatable   :: n113(:)
     
     integer                :: fcst_hour
     integer                :: init_yr, init_mo, init_da, init_hr
     real, allocatable      :: metdata1(:,:) 
     real, allocatable      :: metdata2(:,:)
     real, allocatable      :: metdata3(:,:)
     integer                :: nmodels   

  end type galwemrad_type_dec

  type(galwemrad_type_dec), allocatable :: galwemrad_struc(:)
!EOP
contains

!BOP
!
! !ROUTINE: init_galwemrad
! \label{init_galwemrad}
! 
! !INTERFACE:
  subroutine init_galwemrad(findex)
! !USES: 
    use LIS_coreMod,    only : LIS_rc, LIS_domain
    use LIS_timeMgrMod, only : LIS_update_timestep
    use LIS_logMod,     only : LIS_logunit, LIS_endrun

    implicit none
! !USES: 
    integer, intent(in)  :: findex
! 
! !DESCRIPTION: 
!  Defines the native resolution of the input forcing for GALWEM
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LIS interpolation
!  schemes (see Section~\ref{interp}).
!
!EOP

    integer :: n
    real    :: gridDesci(LIS_rc%nnest,50)
    
    write(LIS_logunit,*) "[INFO] Initializing the GALWEM radiation inputs "

    ! Forecast mode -- NOT Available at this time for this forcing reader:
    if( LIS_rc%forecastMode.eq.1 ) then
       write(LIS_logunit,*) '[ERR] Currently the GALWEM radiation forcing reader'
       write(LIS_logunit,*) '[ERR] is not set up to run in forecast mode.'
       write(LIS_logunit,*) '[ERR] LIS forecast run-time ending.'
       call LIS_endrun()
    endif

    allocate(galwemrad_struc(LIS_rc%nnest))
    call readcrd_galwemrad()

    do n=1, LIS_rc%nnest
       galwemrad_struc(n)%ts = 3600  !check
       call LIS_update_timestep(LIS_rc, n, galwemrad_struc(n)%ts)
    enddo

    do n=1, LIS_rc%nnest
       if(galwemrad_struc(n)%resol == 17) then ! galwemrad-17km
          galwemrad_struc(:)%nc = 1536  
          galwemrad_struc(:)%nr = 1152
       elseif(galwemrad_struc(n)%resol == 25) then ! galwemrad-25deg
          galwemrad_struc(:)%nc = 1440
          galwemrad_struc(:)%nr = 721
       else
          write(LIS_logunit,*) '[ERR] Currently the GALWEM radiation forcing reader'
          write(LIS_logunit,*) '[ERR] supports 17 km and 25 deg datasets.'
          write(LIS_logunit,*) '[ERR] LIS radiation run-time ending.'
          call LIS_endrun()
       endif
    enddo

    ! 2 - key met field
    LIS_rc%met_nf(findex) = 2  

    do n=1,LIS_rc%nnest
     
       ! Check if starting hour of LIS run matches 25/05/16/12Z:   ??check HKB??
!       if((LIS_rc%shr .ne.  0) .and. (LIS_rc%shr .ne.  6) .and. &
!          (LIS_rc%shr .ne. 12) .and. (LIS_rc%shr .ne. 18)) then
!          write(LIS_logunit,*) "[ERR] GALWEM radiation type begins"
!          write(LIS_logunit,*) "[ERR] at 25/05/16/12Z for a radiation window, so the "
!          write(LIS_logunit,*) "[ERR] 'Starting hour:' should be set to 25/5/16/12 in"
!          write(LIS_logunit,*) "[ERR]  your lis.config file.."
!          call LIS_endrun()
!       endif
       
       allocate(galwemrad_struc(n)%metdata1(LIS_rc%met_nf(findex),LIS_rc%ngrid(n)))
       allocate(galwemrad_struc(n)%metdata2(LIS_rc%met_nf(findex),LIS_rc%ngrid(n)))
       allocate(galwemrad_struc(n)%metdata3(LIS_rc%met_nf(findex),LIS_rc%ngrid(n)))

       ! Initialize the radiation initial date-time and grib record:
       galwemrad_struc(n)%init_yr = LIS_rc%syr
       galwemrad_struc(n)%init_mo = LIS_rc%smo
       galwemrad_struc(n)%init_da = LIS_rc%sda
       galwemrad_struc(n)%init_hr = 6*(LIS_rc%shr/6)  !multiple of 6 hours
       galwemrad_struc(n)%fcst_hour = LIS_rc%shr-galwemrad_struc(n)%init_hr

       galwemrad_struc(n)%metdata1 = 0
       galwemrad_struc(n)%metdata2 = 0
       galwemrad_struc(n)%metdata3 = 0
       gridDesci = 0

       if(galwemrad_struc(n)%resol == 17) then   !galwemrad-17km 
          gridDesci(n,1) = 0
          gridDesci(n,2) = galwemrad_struc(n)%nc !gnc
          gridDesci(n,3) = galwemrad_struc(n)%nr !gnr
          gridDesci(n,4) = -89.921875         !lat(1,1)
          gridDesci(n,5) = -179.882813        !lon(1,1)
          gridDesci(n,6) = 128
          gridDesci(n,7) = 89.921875          !lat(gnc,gnr)
          gridDesci(n,8) = 179.882813         !lon(gnc,gnr)
          gridDesci(n,9) = 0.234375           !dx
          gridDesci(n,10) = 0.15625           !dy
          gridDesci(n,20) = 0
       endif

       if(galwemrad_struc(n)%resol == 25) then   !galwemrad-25deg
          gridDesci(n,1) = 0
          gridDesci(n,2) = galwemrad_struc(n)%nc !gnc
          gridDesci(n,3) = galwemrad_struc(n)%nr !gnr
          gridDesci(n,4) = -90.0              !lat(1,1)
          gridDesci(n,5) = -180.0             !lon(1,1)
          gridDesci(n,6) = 128
          gridDesci(n,7) = 90.0               !lat(gnc,gnr)
          gridDesci(n,8) = 179.75             !lon(gnc,gnr)
          gridDesci(n,9) = 0.25               !dx
          gridDesci(n,10) = 0.25              !dy
          gridDesci(n,20) = 0
       endif

       galwemrad_struc(n)%mi = galwemrad_struc(n)%nc*galwemrad_struc(n)%nr
       galwemrad_struc(n)%fcsttime1 = 3000.0
       galwemrad_struc(n)%fcsttime2 = 0.0
       galwemrad_struc(n)%fcsttime3 = 0.0
    enddo

    do n=1,LIS_rc%nnest
       !Setting up weights for Interpolation
       if(trim(LIS_rc%met_interp(findex)).eq."bilinear") then
          allocate(galwemrad_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(galwemrad_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(galwemrad_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(galwemrad_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(galwemrad_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(galwemrad_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(galwemrad_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(galwemrad_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

          call bilinear_interp_input(n,gridDesci(n,:),&
               galwemrad_struc(n)%n111,galwemrad_struc(n)%n121,&
               galwemrad_struc(n)%n211,galwemrad_struc(n)%n221,&
               galwemrad_struc(n)%w111,galwemrad_struc(n)%w121,&
               galwemrad_struc(n)%w211,galwemrad_struc(n)%w221)

       elseif(trim(LIS_rc%met_interp(findex)).eq."budget-bilinear") then
          allocate(galwemrad_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(galwemrad_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(galwemrad_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(galwemrad_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(galwemrad_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(galwemrad_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(galwemrad_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(galwemrad_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

          call bilinear_interp_input(n,gridDesci(n,:),&
               galwemrad_struc(n)%n111,galwemrad_struc(n)%n121,&
               galwemrad_struc(n)%n211,galwemrad_struc(n)%n221,&
               galwemrad_struc(n)%w111,galwemrad_struc(n)%w121,&
               galwemrad_struc(n)%w211,galwemrad_struc(n)%w221)

          allocate(galwemrad_struc(n)%n112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(galwemrad_struc(n)%n122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(galwemrad_struc(n)%n212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(galwemrad_struc(n)%n222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(galwemrad_struc(n)%w112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(galwemrad_struc(n)%w122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(galwemrad_struc(n)%w212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(galwemrad_struc(n)%w222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))

          call conserv_interp_input(n,gridDesci(n,:),&
               galwemrad_struc(n)%n112,galwemrad_struc(n)%n122,&
               galwemrad_struc(n)%n212,galwemrad_struc(n)%n222,&
               galwemrad_struc(n)%w112,galwemrad_struc(n)%w122,&
               galwemrad_struc(n)%w212,galwemrad_struc(n)%w222)
       elseif(trim(LIS_rc%met_interp(findex)).eq."neighbor") then
          allocate(galwemrad_struc(n)%n113(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

          call neighbor_interp_input(n,gridDesci(n,:),&
               galwemrad_struc(n)%n113)
       endif
    enddo

  end subroutine init_galwemrad
end module galwemrad_forcingMod
