!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT----------------------
module gpcp_forcingMod
!BOP
! !MODULE: gpcp_forcingMod
! 
! !DESCRIPTION: 
!  This module contains variables and data structures that are used
!  for the implementation of the precipitation data from the
!  NASA/GSFC's Global Precipitation Climatology Project (GPCP)
!  Daily Precipitation Analysis Climate Data Record (CDR) version 13rA1. 
!  The GPCP daily product blends data from multiple satellites
!  including SSM/I, SSM/IS, TOVS, AIRS, GOES, MeteoSat, GMS, MTSat, and 
!  GPCP monthly analysis.
!  The 3-hourly product is obtained by disaggregating 
!  GPCP1DD using the GDAS forcing is used in this implementation. 
! 
!   upto 2000/1/24          :   T126 (384x190)  grid
!   2001/1/25  - 2002/10/29 :   T170 (512x256)  grid
!   2002/10/30 - 2005/5/31  :   T254 (768x384)  grid
!   2005/6/1  - 2010/7/28   :   T382 (1152x576) grid
!   2010/7/29 - 2015/01/14  :   T574 (1760x880) grid
!   2015/1/15  onwards      :   T1534 (3072x1536) grid
!
!  The implementation in LIS has the derived data type {\tt gpcp\_struc}
!  that includes the variables to specify the runtime options, and 
!  the weights and neighbor information for spatial interpolation.
! 
!  They are desribed below: 
! \begin{description}
!  \item[ncold]
!    Number of columns (along the east west dimension) for the input data
!  \item[nrold]
!    Number of rows (along the north south dimension) for the input data
!  \item[gpcpdir]
!    Directory containing the input data
!  \item[gpcptime]
!    The nearest, hourly instance of the incoming 
!    data (as a real time).
!  \item[griduptime1]
!    The time to switch the input resolution to T170
!  \item[mi]
!    Number of points in the input grid
!  \item[n11,n121,n211,n221]
!    Arrays containing the neighbor information of the input grid 
!    for each grid point in LIS, for bilinear interpolation. 
!  \item[w111,w121,w211,w221]
!    Arrays containing the weights of the input grid 
!    for each grid point in LIS, for bilinear interpolation.
!  \item[n112,n122,n212,n222]
!    Arrays containing the neighbor information of the input grid 
!    for each grid point in LIS, for conservative interpolation. 
!  \item[w112,w122,w212,w222]
!    Arrays containing the weights of the input grid 
!    for each grid point in LIS, for conservative interpolation.
!  \end{description}
!
! !USES: 
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none
  
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: init_gpcp      !defines the native resolution of 
                           !the input data
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: gpcp_struc

!EOP
 
  type, public ::  gpcp_type_dec
     real                   :: ts
     integer                :: ncold
     integer                :: nrold  
     character(len=LIS_CONST_PATH_LEN)           :: gpcpdir  
     character*50           :: met_interp
     real*8                 :: gpcptime
     real*8                 :: griduptime1
     real*8                 :: griduptime2
     real*8                 :: griduptime3
     real*8                 :: griduptime4
     real*8                 :: griduptime5
     real*8                 :: griduptime6
     logical                :: gridchange1
     logical                :: gridchange2
     logical                :: gridchange3
     logical                :: gridchange4
     logical                :: gridchange5
     logical                :: gridchange6
     integer                :: mi

   ! Bilinear interp weights and indices
     integer, allocatable   :: n111(:)
     integer, allocatable   :: n121(:)
     integer, allocatable   :: n211(:)
     integer, allocatable   :: n221(:)
     real, allocatable      :: w111(:),w121(:)
     real, allocatable      :: w211(:),w221(:)

   ! Budget-Bilinear (conserve) interp weights and indices
     integer, allocatable   :: n112(:,:)
     integer, allocatable   :: n122(:,:)
     integer, allocatable   :: n212(:,:)
     integer, allocatable   :: n222(:,:)
     real,    allocatable   :: w112(:,:),w122(:,:)
     real,    allocatable   :: w212(:,:),w222(:,:)

     real, allocatable :: metdata1(:,:)
     real, allocatable :: metdata2(:,:)

  end type gpcp_type_dec

  type(gpcp_type_dec), allocatable :: gpcp_struc(:)

contains
  
!BOP
!
! !ROUTINE: init_gpcp
! \label{init_gpcp}
!
! !REVISION HISTORY: 
! 11Dec2003: Sujay Kumar; Initial Specification
! 02Dec2014: KR Arsenault: Added new grid change update (~2012)
! 03Jun2016: H Beaudoing: Adopted to the 3-hourly GPCP
! 
! !INTERFACE:
  subroutine init_gpcp(findex)
! !USES: 
    use LIS_coreMod,    only : LIS_rc, LIS_domain
    use LIS_timeMgrMod, only : LIS_date2time, LIS_update_timestep
    use LIS_logMod,     only : LIS_logunit, LIS_endrun

    implicit none
    integer, intent(in) :: findex

! 
! !DESCRIPTION: 
!  Defines the native resolution of the input forcing for GPCP
!  data. The grid description arrays are based on the decoding
!  schemes used by NCEP and followed in the LIS interpolation
!  schemes (see Section~\ref{interp}).
!
!  The routines invoked are: 
!  \begin{description}
!   \item[readcrd\_gpcp](\ref{readcrd_gpcp}) \newline
!     reads the runtime options specified for GPCP data
!   \item[LIS\_date2time](\ref{LIS_date2time}) \newline
!     converts date to the real time format
!   \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    computes the neighbor, weights for conservative interpolation
!   \item[upscaleByAveraging\_input](\ref{upscaleByAveraging_input}) \newline
!    computes the neighbors for upscaling by averaging
!  \end{description}
!
!EOP

    real    :: gridDesci(50)
    integer :: updoy, yr1,mo1,da1,hr1,mn1,ss1
    real    :: upgmt
    integer :: n 

    allocate(gpcp_struc(LIS_rc%nnest))

!- Read in LIS Config file entries:
    call readcrd_gpcp()

    do n=1, LIS_rc%nnest
       gpcp_struc(n)%ts = 3*60*60   ! 3-hr
       call LIS_update_timestep(LIS_rc, n, gpcp_struc(n)%ts)
    enddo

    LIS_rc%met_nf(findex) = 1

 !- Set interp arrays for reinterpolation later:
    do n=1,LIS_rc%nnest

       allocate(gpcp_struc(n)%metdata1(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))
       allocate(gpcp_struc(n)%metdata2(LIS_rc%met_nf(findex),&
            LIS_rc%ngrid(n)))

       gpcp_struc(n)%metdata1 = 0
       gpcp_struc(n)%metdata2 = 0

!------------------------------------------------------------------------ 
!- Initialize first available grid domain/res for GPCP/GDAS:
!  upto 2000/1/24          :   T126 (384x190)  grid
!- GPCP/GDAS Grid description:
       gpcp_struc(n)%ncold = 384
       gpcp_struc(n)%nrold = 190

       gridDesci = 0
       gridDesci(1) = 4
       gridDesci(2) = 384
       gridDesci(3) = 190
       gridDesci(4) = -89.277
       gridDesci(5) = 0
       gridDesci(6) = 128
       gridDesci(7) = 89.277
       gridDesci(8) = -0.9375
       gridDesci(9) = 0.9375
       gridDesci(10) = 95
       gridDesci(11) = 64
       gridDesci(20) = 0
       gpcp_struc(n)%mi = gridDesci(2)*gridDesci(3)

     ! This grid is good for some time in the 1990's.
     ! Look up the exact dates.
       yr1 = 1991     !grid update time
       mo1 = 01
       da1 = 01
       hr1 = 12
       mn1 = 0; ss1 = 0
       call LIS_date2time(gpcp_struc(n)%griduptime1,updoy,upgmt,yr1,mo1,&
            da1,hr1,mn1,ss1 )

! ==  2001/1/25  - 2002/10/29 :   T170 (512x256)  grid
       yr1 = 2000     !grid update time
       mo1 = 01
       da1 = 25
       hr1 = 0
       mn1 = 0; ss1 = 0
       call LIS_date2time(gpcp_struc(n)%griduptime2,updoy,upgmt,yr1,mo1,&
            da1,hr1,mn1,ss1 )

! ==  2002/10/30 - 2005/5/31  :   T254 (768x384)  grid
       yr1 = 2002     !grid update time
       mo1 = 10
       da1 = 30
       hr1 = 0
       mn1 = 0; ss1 = 0
       call LIS_date2time(gpcp_struc(n)%griduptime3,updoy,upgmt,yr1,mo1,&
            da1,hr1,mn1,ss1 )
       
! ==  2005/6/1  - 2010/7/28   :   T382 (1152x576) grid
       yr1 = 2005     !grid update time
       mo1 = 06
       da1 = 01
       hr1 = 0
       mn1 = 0; ss1 = 0
       call LIS_date2time(gpcp_struc(n)%griduptime4,updoy,upgmt,yr1,mo1,&
            da1,hr1,mn1,ss1 )

! ==  2010/7/29 - 2015/01/14  :   T574 (1760x880) grid
       yr1 = 2010     !grid update time
       mo1 = 07
       da1 = 29
       hr1 = 0      
       mn1 = 0; ss1 = 0
       call LIS_date2time(gpcp_struc(n)%griduptime5,updoy,upgmt,yr1,mo1,&
            da1,hr1,mn1,ss1 )

! ==  2015/1/15  onwards      :   T1534 (3072x1536) grid
       yr1 = 2015     !grid update time
       mo1 = 01
       da1 = 15
       hr1 = 0      
       mn1 = 0; ss1 = 0
       call LIS_date2time(gpcp_struc(n)%griduptime6,updoy,upgmt,yr1,mo1,&
            da1,hr1,mn1,ss1 )

       gpcp_struc(n)%gridchange1 = .true.
       gpcp_struc(n)%gridchange2 = .true.
       gpcp_struc(n)%gridchange3 = .true.
       gpcp_struc(n)%gridchange4 = .true.
       gpcp_struc(n)%gridchange5 = .true.
       gpcp_struc(n)%gridchange6 = .true.

       ! Setting up weights for Interpolation
       call gpcp_reset_interp_input(n, findex, gridDesci)
    enddo
  end subroutine init_gpcp

end module gpcp_forcingMod
