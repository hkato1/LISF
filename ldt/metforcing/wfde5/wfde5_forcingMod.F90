!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
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
!    Number of forcing variables in the ECMWF data
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
!    for each grid point in LDT, for bilinear interpolation.
!  \item[w111,w121,w211,w221]
!    Arrays containing the weights of the input grid
!    for each grid point in LDT, for bilinear interpolation.
!  \item[n122,n122,n212,n222]
!    Arrays containing the neighbor information of the input grid
!    for each grid point in LDT, for conservative interpolation.
!  \item[w112,w122,w212,w222]
!    Arrays containing the weights of the input grid
!    for each grid point in LDT, for conservative interpolation.
!  \item[n113]
!    Arrays containing the neighbor information of the input grid
!    for each grid point in LDT, for n. neighbor interpolation.
!  \item[findtime1, findtime2]
!   boolean flags to indicate which time is to be read for
!   temporal interpolation.
!  \end{description}
!
! !USES:
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
     character*100 :: wfde5dir   !WFDE5 Forcing Directory
     character*40 :: presrc     !precip forcing source
     character*100 :: wfde5alt_file
     real*8       :: wfde5time1,wfde5time2
     logical      :: reset_flag
     integer      :: mo1,mo2

     real :: gridDesc(50)
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

     integer            :: nvars
     integer            :: uselml

     real*8             :: ringtime
     
     integer            :: nIter, st_iterid,en_iterid

     real, allocatable      :: tair1(:,:)
     real, allocatable      :: qair1(:,:)
     real, allocatable      :: wind1(:,:)
     real, allocatable      :: ps1(:,:)
     real, allocatable      :: rainf1(:,:)
     real, allocatable      :: snowf1(:,:)
     real, allocatable      :: dirswd1(:,:)
     real, allocatable      :: difswd1(:,:)
     real, allocatable      :: swd1(:,:)
     real, allocatable      :: lwd1(:,:)

     real, allocatable      :: tair2(:,:)
     real, allocatable      :: qair2(:,:)
     real, allocatable      :: wind2(:,:)
     real, allocatable      :: ps2(:,:)
     real, allocatable      :: rainf2(:,:)
     real, allocatable      :: snowf2(:,:)
     real, allocatable      :: dirswd2(:,:)
     real, allocatable      :: difswd2(:,:)
     real, allocatable      :: swd2(:,:)
     real, allocatable      :: lwd2(:,:)

     real, allocatable :: metdata1(:,:) 
     real, allocatable :: metdata2(:,:) 

  end type wfde5_type_dec

  type(wfde5_type_dec), allocatable :: wfde5_struc(:)

contains

!BOP
!
! !ROUTINE: init_wfde5
! \label{init_wfde5}
!
! !REVISION HISTORY:
! 29 Mar 2022: Hiroko Beaudoing, initial code 
!
! !INTERFACE:
  subroutine init_wfde5(findex)

! !USES:
    use LDT_coreMod
    use LDT_timeMgrMod
    use LDT_logMod

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
!  schemes used by NCEP and followed in the LDT interpolation
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

    integer :: updoy, yr1,mo1,da1,hr1,mn1,ss1
    real :: upgmt
    integer :: n
    integer :: ftn
    integer :: G2Pid


    allocate(wfde5_struc(LDT_rc%nnest))

    do n=1,LDT_rc%nnest
       wfde5_struc(n)%ncold = 720
       wfde5_struc(n)%nrold = 360
       wfde5_struc(n)%mo1 = -1
       wfde5_struc(n)%mo2 = -1

       LDT_rc%met_nc(findex) = wfde5_struc(n)%ncold
       LDT_rc%met_nr(findex) = wfde5_struc(n)%nrold
    
       allocate(wfde5_struc(n)%tair1(LDT_rc%lnc(n)*LDT_rc%lnr(n),745))
       allocate(wfde5_struc(n)%qair1(LDT_rc%lnc(n)*LDT_rc%lnr(n),745))
       allocate(wfde5_struc(n)%wind1(LDT_rc%lnc(n)*LDT_rc%lnr(n),745))
       allocate(wfde5_struc(n)%ps1(LDT_rc%lnc(n)*LDT_rc%lnr(n),745))
       allocate(wfde5_struc(n)%rainf1(LDT_rc%lnc(n)*LDT_rc%lnr(n),745))
       allocate(wfde5_struc(n)%swd1(LDT_rc%lnc(n)*LDT_rc%lnr(n),745))
       allocate(wfde5_struc(n)%lwd1(LDT_rc%lnc(n)*LDT_rc%lnr(n),745))

       allocate(wfde5_struc(n)%tair2(LDT_rc%lnc(n)*LDT_rc%lnr(n),745))
       allocate(wfde5_struc(n)%qair2(LDT_rc%lnc(n)*LDT_rc%lnr(n),745))
       allocate(wfde5_struc(n)%wind2(LDT_rc%lnc(n)*LDT_rc%lnr(n),745))
       allocate(wfde5_struc(n)%ps2(LDT_rc%lnc(n)*LDT_rc%lnr(n),745))
       allocate(wfde5_struc(n)%rainf2(LDT_rc%lnc(n)*LDT_rc%lnr(n),745))
       allocate(wfde5_struc(n)%swd2(LDT_rc%lnc(n)*LDT_rc%lnr(n),745))
       allocate(wfde5_struc(n)%lwd2(LDT_rc%lnc(n)*LDT_rc%lnr(n),745))

    enddo

    call readcrd_wfde5()
    LDT_rc%met_nf(findex) = 8
    LDT_rc%met_ts(findex) = 3600
    LDT_rc%met_zterp(findex) = .true. 


    wfde5_struc%reset_flag = .false.

    do n=1, LDT_rc%nnest
       wfde5_struc(n)%ts = 3600  !check
       call LDT_update_timestep(LDT_rc, n, wfde5_struc(n)%ts)
    enddo

    do n=1,LDT_rc%nnest
       wfde5_struc(n)%gridDesc = 0
       wfde5_struc(n)%gridDesc(1) = 0
       wfde5_struc(n)%gridDesc(2) = wfde5_struc(n)%ncold
       wfde5_struc(n)%gridDesc(3) = wfde5_struc(n)%nrold
       wfde5_struc(n)%gridDesc(4) = -89.75
       wfde5_struc(n)%gridDesc(5) = -179.75
       wfde5_struc(n)%gridDesc(6) = 128
       wfde5_struc(n)%gridDesc(7) = 89.75
       wfde5_struc(n)%gridDesc(8) = 179.75
       wfde5_struc(n)%gridDesc(9) = 0.5
       wfde5_struc(n)%gridDesc(10) = 0.5
       wfde5_struc(n)%gridDesc(20) = 0

       LDT_rc%met_gridDesc(findex,1:20) = wfde5_struc(n)%gridDesc(1:20)

       wfde5_struc(n)%mi = wfde5_struc(n)%ncold*wfde5_struc(n)%nrold

       ! Setting up weights for Interpolation
       if(trim(LDT_rc%met_gridtransform(findex)).eq."bilinear") then
          allocate(wfde5_struc(n)%n111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(wfde5_struc(n)%n121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(wfde5_struc(n)%n211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(wfde5_struc(n)%n221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(wfde5_struc(n)%w111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(wfde5_struc(n)%w121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(wfde5_struc(n)%w211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(wfde5_struc(n)%w221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          call bilinear_interp_input(n, wfde5_struc(n)%gridDesc(:),&
               wfde5_struc(n)%n111,wfde5_struc(n)%n121,&
               wfde5_struc(n)%n211,wfde5_struc(n)%n221,&
               wfde5_struc(n)%w111,wfde5_struc(n)%w121,&
               wfde5_struc(n)%w211,wfde5_struc(n)%w221)

       elseif(trim(LDT_rc%met_gridtransform(findex)).eq."budget-bilinear") then
          allocate(wfde5_struc(n)%n111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(wfde5_struc(n)%n121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(wfde5_struc(n)%n211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(wfde5_struc(n)%n221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(wfde5_struc(n)%w111(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(wfde5_struc(n)%w121(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(wfde5_struc(n)%w211(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          allocate(wfde5_struc(n)%w221(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          call bilinear_interp_input(n, wfde5_struc(n)%gridDesc(:),&
               wfde5_struc(n)%n111,wfde5_struc(n)%n121,&
               wfde5_struc(n)%n211,wfde5_struc(n)%n221,&
               wfde5_struc(n)%w111,wfde5_struc(n)%w121,&
               wfde5_struc(n)%w211,wfde5_struc(n)%w221)

          allocate(wfde5_struc(n)%n112(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(wfde5_struc(n)%n122(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(wfde5_struc(n)%n212(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(wfde5_struc(n)%n222(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(wfde5_struc(n)%w112(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(wfde5_struc(n)%w122(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(wfde5_struc(n)%w212(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          allocate(wfde5_struc(n)%w222(LDT_rc%lnc(n)*LDT_rc%lnr(n),25))
          call conserv_interp_input(n, wfde5_struc(n)%gridDesc(:),&
               wfde5_struc(n)%n112,wfde5_struc(n)%n122,&
               wfde5_struc(n)%n212,wfde5_struc(n)%n222,&
               wfde5_struc(n)%w112,wfde5_struc(n)%w122,&
               wfde5_struc(n)%w212,wfde5_struc(n)%w222)

       elseif(trim(LDT_rc%met_gridtransform(findex)).eq."neighbor") then
          allocate(wfde5_struc(n)%n113(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
          call neighbor_interp_input(n, wfde5_struc(n)%gridDesc(:),&
               wfde5_struc(n)%n113)

       else
          write(LDT_logunit,*) '[ERR] Interpolation option '// &
               trim(LDT_rc%met_gridtransform(findex))//&
               ' for WFDE5 forcing is not supported'
          call LDT_endrun()
       endif

       call LDT_registerAlarm("WFDE5 forcing alarm",&
            86400.0,86400.0)
       wfde5_struc(n)%startFlag = .true.
       wfde5_struc(n)%dayFlag = .true.

       wfde5_struc(n)%nvars = 8

       allocate(wfde5_struc(n)%metdata1(LDT_rc%met_nf(findex),&
            LDT_rc%ngrid(n)))
       allocate(wfde5_struc(n)%metdata2(LDT_rc%met_nf(findex),&
            LDT_rc%ngrid(n)))

       wfde5_struc(n)%metdata1 = 0
       wfde5_struc(n)%metdata2 = 0

    enddo   ! End nest loop
    
    
  end subroutine init_wfde5
end module wfde5_forcingMod

