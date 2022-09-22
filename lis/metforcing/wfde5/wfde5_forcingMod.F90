!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
module wfde5_forcingMod
!BOP
! !MODULE: wfde5_forcingMod
!
! !DESCRIPTION:
!  This module contains variables and data structures that are used
!  for the implementation of the WFDE5 forcing data.
!  The data is global 0.5 degree dataset in latlon
!  projection, and at 1 hourly intervals. The derived
!  data type {\tt wfde5\_struc}
!  includes the variables that specify the runtime options, and the
!  weights and neighbor information to be used for spatial interpolation.
!  They are described below:
!  \begin{description}
!  \item[ncold]
!    Number of columns (along the east west dimension) for the input data
!  \item[nrold]
!    Number of rows (along the north south dimension) for the input data
!  \item[nmif]
!    Number of forcing variables in the WFDE5 data
!  \item[wfde5time1]
!    The nearest, previous 1 hour instance of the incoming
!    data (as a real time).
!  \item[wfde5time2]
!    The nearest, next 1 hour instance of the incoming
!    data (as a real time).
!  \item[wfde5dir]
!    Directory containing the input data
!  \item[mi]
!    Number of points in the input grid
!  \item[n111,n121,n211,n221]
!    Arrays containing the neighbor information of the input grid
!    for each grid point in LIS, for bilinear interpolation.
!  \item[w111,w121,w211,w221]
!    Arrays containing the weights of the input grid
!    for each grid point in LIS, for bilinear interpolation.
!  \item[n122,n122,n212,n222]
!    Arrays containing the neighbor information of the input grid
!    for each grid point in LIS, for conservative interpolation.
!  \item[w112,w122,w212,w222]
!    Arrays containing the weights of the input grid
!    for each grid point in LIS, for conservative interpolation.
!  \item[n113]
!    Arrays containing the neighbor information of the input grid
!    for each grid point in LIS, for n. neighbor interpolation.
!  \item[findtime1, findtime2]
!   boolean flags to indicate which time is to be read for
!   temporal interpolation.
!  \end{description}
!
! !USES:
  use ESMF
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_wfde5      !defines the native resolution of
                             !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: wfde5_struc

!EOP
  type, public ::  wfde5_type_dec

     real         :: ts
     integer      :: ncold, nrold
     character(len=LIS_CONST_PATH_LEN) :: wfde5dir   !WFDE5 Forcing Directory
     character(len=LIS_CONST_PATH_LEN) :: presrc     !precip forcing source
     character(len=LIS_CONST_PATH_LEN) :: wfde5alt_file   ! altitude file
     character*50 :: met_interp
     real*8       :: wfde5time1,wfde5time2
     logical      :: reset_flag
     integer      :: mo1,mo2

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
     integer                :: findtime1, findtime2
     logical                :: startFlag, dayFlag

     real, allocatable      :: tair1(:,:)
     real, allocatable      :: qair1(:,:)
     real, allocatable      :: wind1(:,:)
     real, allocatable      :: ps1(:,:)
     real, allocatable      :: rainf1(:,:)
     real, allocatable      :: snowf1(:,:)
     real, allocatable      :: swd1(:,:)
     real, allocatable      :: lwd1(:,:)

     real, allocatable      :: tair2(:,:)
     real, allocatable      :: qair2(:,:)
     real, allocatable      :: wind2(:,:)
     real, allocatable      :: ps2(:,:)
     real, allocatable      :: rainf2(:,:)
     real, allocatable      :: snowf2(:,:)
     real, allocatable      :: swd2(:,:)
     real, allocatable      :: lwd2(:,:)

     integer            :: nvars
     integer            :: uselml

     real*8             :: ringtime
     type(ESMF_Time)    :: startTime
     
     integer            :: nIter, st_iterid,en_iterid

     real, allocatable :: metdata1(:,:,:) 
     real, allocatable :: metdata2(:,:,:) 

  end type wfde5_type_dec

  type(wfde5_type_dec), allocatable :: wfde5_struc(:)

contains

!BOP
!
! !ROUTINE: init_wfde5
! \label{init_wfde5}
!
! !REVISION HISTORY:
! 23 Mar 2022: Hiroko Beaudoing, initial code 
!
! !INTERFACE:
  subroutine init_wfde5(findex)

! !USES:
    use LIS_coreMod,               only : LIS_rc, LIS_isatAfinerResolution
    use LIS_timeMgrMod
    use LIS_logMod
    use LIS_spatialDownscalingMod, only : LIS_init_pcpclimo_native
    use LIS_forecastMod
#if(defined USE_NETCDF3 || defined USE_NETCDF4)      
  use netcdf
#endif

    implicit none
! !AGRUMENTS:
    integer, intent(in) :: findex
!
! !DESCRIPTION:
!  Defines the native resolution of the input forcing for WFDE5
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LIS interpolation
!  schemes (see Section~\ref{interp}).
!
!  The routines invoked are:
!  \begin{description}
!   \item[readcrd\_wfde5](\ref{readcrd_wfde5}) \newline
!     reads the runtime options specified for WFDE5 data
!   \item[bilinear\_interp\_input](\ref{bilinear_interp_input}) \newline
!    computes the neighbor, weights for bilinear interpolation
!   \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    computes the neighbor, weights for conservative interpolation
!  \end{description}
!EOP
    real :: gridDesci(LIS_rc%nnest,50)
    integer :: updoy, yr1,mo1,da1,hr1,mn1,ss1
    real :: upgmt
    integer :: n
    integer :: ftn


    allocate(wfde5_struc(LIS_rc%nnest))

    do n=1,LIS_rc%nnest
       
       wfde5_struc(n)%ncold = 720
       wfde5_struc(n)%nrold = 360
       wfde5_struc(n)%mo1 = -1
       wfde5_struc(n)%mo2 = -1

       ! 745 should be ndays*24hr?  
       allocate(wfde5_struc(n)%tair1(LIS_rc%lnc(n)*LIS_rc%lnr(n),745))
       allocate(wfde5_struc(n)%qair1(LIS_rc%lnc(n)*LIS_rc%lnr(n),745))
       allocate(wfde5_struc(n)%wind1(LIS_rc%lnc(n)*LIS_rc%lnr(n),745))
       allocate(wfde5_struc(n)%ps1(LIS_rc%lnc(n)*LIS_rc%lnr(n),745))
       allocate(wfde5_struc(n)%rainf1(LIS_rc%lnc(n)*LIS_rc%lnr(n),745))
       allocate(wfde5_struc(n)%swd1(LIS_rc%lnc(n)*LIS_rc%lnr(n),745))
       allocate(wfde5_struc(n)%lwd1(LIS_rc%lnc(n)*LIS_rc%lnr(n),745))

       allocate(wfde5_struc(n)%tair2(LIS_rc%lnc(n)*LIS_rc%lnr(n),745))
       allocate(wfde5_struc(n)%qair2(LIS_rc%lnc(n)*LIS_rc%lnr(n),745))
       allocate(wfde5_struc(n)%wind2(LIS_rc%lnc(n)*LIS_rc%lnr(n),745))
       allocate(wfde5_struc(n)%ps2(LIS_rc%lnc(n)*LIS_rc%lnr(n),745))
       allocate(wfde5_struc(n)%rainf2(LIS_rc%lnc(n)*LIS_rc%lnr(n),745))
       allocate(wfde5_struc(n)%swd2(LIS_rc%lnc(n)*LIS_rc%lnr(n),745))
       allocate(wfde5_struc(n)%lwd2(LIS_rc%lnc(n)*LIS_rc%lnr(n),745))

    enddo

    call readcrd_wfde5()
    LIS_rc%met_nf(findex) = 8

    wfde5_struc%reset_flag = .false.

    do n=1, LIS_rc%nnest
       wfde5_struc(n)%ts = 3600  !check
       call LIS_update_timestep(LIS_rc, n, wfde5_struc(n)%ts)
    enddo

    gridDesci = 0

    do n=1,LIS_rc%nnest
       gridDesci(n,1) = 0
       gridDesci(n,2) = wfde5_struc(n)%ncold
       gridDesci(n,3) = wfde5_struc(n)%nrold
       gridDesci(n,4) = -89.75
       gridDesci(n,5) = -179.75
       gridDesci(n,6) = 128
       gridDesci(n,7) = 89.75
       gridDesci(n,8) = 179.75
       gridDesci(n,9) = 0.5
       gridDesci(n,10) = 0.5
       gridDesci(n,20) = 0

       wfde5_struc(n)%mi = wfde5_struc(n)%ncold*wfde5_struc(n)%nrold

       ! Check resolution and set up weights for Interpolation
       if ( LIS_isatAfinerResolution(n,gridDesci(n,9)) ) then
         wfde5_struc(n)%met_interp = LIS_rc%met_interp(findex)

         write(LIS_logunit,*) 'MSG: The WFDE5 forcing resolution is ' // &
                              ' coaser than the running domain.'
         write(LIS_logunit,*) '     Interpolating with the ' // &
                               trim(wfde5_struc(n)%met_interp) // ' method.'

         if(trim(wfde5_struc(n)%met_interp).eq."bilinear") then
          allocate(wfde5_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(wfde5_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(wfde5_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(wfde5_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(wfde5_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(wfde5_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(wfde5_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(wfde5_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          call bilinear_interp_input(n, gridDesci(n,:),&
               wfde5_struc(n)%n111,wfde5_struc(n)%n121,&
               wfde5_struc(n)%n211,wfde5_struc(n)%n221,&
               wfde5_struc(n)%w111,wfde5_struc(n)%w121,&
               wfde5_struc(n)%w211,wfde5_struc(n)%w221)

         elseif(trim(wfde5_struc(n)%met_interp).eq."budget-bilinear") then
          allocate(wfde5_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(wfde5_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(wfde5_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(wfde5_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(wfde5_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(wfde5_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(wfde5_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          allocate(wfde5_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          call bilinear_interp_input(n, gridDesci(n,:),&
               wfde5_struc(n)%n111,wfde5_struc(n)%n121,&
               wfde5_struc(n)%n211,wfde5_struc(n)%n221,&
               wfde5_struc(n)%w111,wfde5_struc(n)%w121,&
               wfde5_struc(n)%w211,wfde5_struc(n)%w221)

          allocate(wfde5_struc(n)%n112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(wfde5_struc(n)%n122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(wfde5_struc(n)%n212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(wfde5_struc(n)%n222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(wfde5_struc(n)%w112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(wfde5_struc(n)%w122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(wfde5_struc(n)%w212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          allocate(wfde5_struc(n)%w222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
          call conserv_interp_input(n, gridDesci(n,:),&
               wfde5_struc(n)%n112,wfde5_struc(n)%n122,&
               wfde5_struc(n)%n212,wfde5_struc(n)%n222,&
               wfde5_struc(n)%w112,wfde5_struc(n)%w122,&
               wfde5_struc(n)%w212,wfde5_struc(n)%w222)

         elseif(trim(wfde5_struc(n)%met_interp).eq."neighbor") then
          allocate(wfde5_struc(n)%n113(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          call neighbor_interp_input(n, gridDesci(n,:),&
               wfde5_struc(n)%n113)

         endif
       else
         wfde5_struc(n)%met_interp = LIS_rc%met_upscale(findex)

         write(LIS_logunit,*) 'MSG: The WFDE5 forcing resolution is finer ' // &
                              'than the running domain.'
         write(LIS_logunit,*) '     Upscaling with the ' // &
                               trim(wfde5_struc(n)%met_interp) // ' method.'

         select case( wfde5_struc(n)%met_interp )
         case( "average" )
            allocate(wfde5_struc(n)%n111(wfde5_struc(n)%mi))

            call upscaleByAveraging_input(gridDesci,                   &
                                          LIS_rc%gridDesc(n,:),        &
                                          wfde5_struc(n)%mi,           &
                                          LIS_rc%lnc(n)*LIS_rc%lnr(n), &
                                          wfde5_struc(n)%n111)
         case default
            write(LIS_logunit,*) '[ERR] Interpolation option '// &
                                trim(LIS_rc%met_interp(findex))//&
                                 ' for WFDE5 forcing is not supported'
            call LIS_endrun()
         end select
       endif

       call LIS_registerAlarm("WFDE5 forcing alarm",&
            86400.0,86400.0)
       wfde5_struc(n)%startFlag = .true.
       wfde5_struc(n)%dayFlag = .true.

       wfde5_struc(n)%nvars = 8

       ! Forecast mode:
       if(LIS_rc%forecastMode.eq.1) then 
          
          if(mod(LIS_rc%nensem(n),&
               LIS_forecast_struc(1)%niterations).ne.0) then 
             write(LIS_logunit,*) '[ERR] The number of ensembles must be a multiple'
             write(LIS_logunit,*) '[ERR] of the number of iterations '
             write(LIS_logunit,*) '[ERR] nensem = ',LIS_rc%nensem(n)
             write(LIS_logunit,*) '[ERR] niter = ',LIS_forecast_struc(1)%niterations
             call LIS_endrun()
          endif

          wfde5_struc(n)%st_iterid = LIS_forecast_struc(1)%st_iterId
          wfde5_struc(n)%en_iterId = LIS_forecast_struc(1)%niterations
          wfde5_struc(n)%nIter = LIS_forecast_struc(1)%niterations
          
          allocate(wfde5_struc(n)%metdata1(LIS_forecast_struc(1)%niterations,&
               LIS_rc%met_nf(findex),&
               LIS_rc%ngrid(n)))
          allocate(wfde5_struc(n)%metdata2(LIS_forecast_struc(1)%niterations,&
               LIS_rc%met_nf(findex),&
               LIS_rc%ngrid(n)))
          
       ! Regular retrospective or non-forecast mode:
       else

          wfde5_struc(n)%st_iterid = 1
          wfde5_struc(n)%en_iterId = 1
          wfde5_struc(n)%nIter = 1
          
          allocate(wfde5_struc(n)%metdata1(1,LIS_rc%met_nf(findex),&
               LIS_rc%ngrid(n)))
          allocate(wfde5_struc(n)%metdata2(1,LIS_rc%met_nf(findex),&
               LIS_rc%ngrid(n)))
          
       endif

       wfde5_struc(n)%metdata1 = 0
       wfde5_struc(n)%metdata2 = 0


       ! Set up precipitation climate downscaling:
       if(LIS_rc%pcp_downscale(findex).ne.0) then
          call LIS_init_pcpclimo_native(n,findex,&
               wfde5_struc(n)%ncold,&
               wfde5_struc(n)%nrold)
       endif

       if ( LIS_rc%met_ecor(findex) == "lapse-rate" .or. &
            LIS_rc%met_ecor(findex) == "lapse-rate and slope-aspect" ) then

          call read_wfde5_elev(n,findex)
       endif

    enddo   ! End nest loop
    
    
  end subroutine init_wfde5
end module wfde5_forcingMod

