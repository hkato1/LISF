!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: get_gpcp
! \label{get_gpcp}
!
! !REVISION HISTORY:
! 03 Jun 2016: H Beaudoing: Adopted CMAP routine for GPCP 3-hourly data
!
! !INTERFACE:
subroutine get_gpcp(n,findex)
! !USES:
  use LIS_coreMod,     only : LIS_rc, LIS_domain
  use LIS_timeMgrMod,  only : LIS_tick, LIS_get_nstep
  use LIS_logMod,      only : LIS_logunit, LIS_endrun
  use gpcp_forcingMod, only : gpcp_struc
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  integer, intent(in) :: findex

!  
! !DESCRIPTION:
!  Opens, reads, and interpolates 3-hrly, GPCP forcing. 
!  At the beginning of a simulation, the code 
!  reads the most recent past data (nearest 3 hour interval), and
!  the nearest future data. These two datasets are used to 
!  temporally interpolate the data to the current model timestep. 
!  The GDAS grid resolution changes are adjusted to daily.
!
!   upto 2000/1/24          :   T126 (384x190)  grid
!   2001/1/25  - 2002/10/29 :   T170 (512x256)  grid
!   2002/10/30 - 2005/5/31  :   T254 (768x384)  grid
!   2005/6/1  - 2010/7/28   :   T382 (1152x576) grid
!   2010/7/29 - 2015/01/14  :   T574 (1760x880) grid
!   2015/1/15  onwards      :   T1534 (3072x1536) grid
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!  \item[LIS\_tick](\ref{LIS_tick}) \newline
!    determines the GPCP data times
!  \item[gpcpfile](\ref{gpcpfile}) \newline
!    Puts together appropriate file name for 6 hour intervals
!  \item[read\_gpcp](\ref{read_gpcp}) \newline
!    Interpolates GPCP data to LIS grid
!   \item[gpcp\_reset\_interp\_input](\ref{gpcp_reset_interp_input} \newline
!    resets the neighbours and weights arrays upon a grid change
!  \end{description}
!EOP
   
!==== Local Variables=======================
  integer :: ferror_gpcp   ! Error flags for precip data sources
  integer :: doy1, yr1, mo1, da1, hr1, mn1, ss1
  integer :: doy5, yr5, mo5, da5, hr5, mn5, ss5
  real*8  :: ctime,ftime_gpcp    ! Current LDAS time and end boundary times for precip data sources 
  integer :: order
  real    :: gmt1,gmt5,ts1,ts5   ! GMT times for current LDAS time and end boundary times for precip data sources
  real    :: gridDesci(50)
  character(len=LIS_CONST_PATH_LEN) :: filename ! Filename variables for precip data sources
!=== End Variable Definition =======================

!------------------------------------------------------------------------
! Determine required observed precip data times 
! (current, accumulation end time)
! Model current time
!------------------------------------------------------------------------
  yr1 = LIS_rc%yr  !current time
  mo1 = LIS_rc%mo
  da1 = LIS_rc%da
  hr1 = LIS_rc%hr
  mn1 = LIS_rc%mn
  ss1 = 0
  ts1 = 0
  call LIS_tick( ctime, doy1, gmt1, yr1, mo1, da1, hr1, mn1, ss1, ts1 )   
!------------------------------------------------------------------------ 
! GPCP product end time
!------------------------------------------------------------------------
  yr5 = LIS_rc%yr  !end accumulation time data
  mo5 = LIS_rc%mo
  da5 = LIS_rc%da
  hr5 = 3*(LIS_rc%hr/3)
  mn5 = 0
  ss5 = 0
  ts5 = 3*60*60
  call LIS_tick( ftime_gpcp, doy5, gmt5, yr5, mo5, da5, hr5, mn5, ss5, ts5 )

!------------------------------------------------------------------------
! Reinitialize the weights and neighbors for each required grid change
!------------------------------------------------------------------------

! 1979-2000 T126 Grid:
  if ( ctime > gpcp_struc(n)%griduptime1 .and. &
       ctime < gpcp_struc(n)%griduptime2 .and. &
       gpcp_struc(n)%gridchange1 ) then

     write(LIS_logunit,*) "** "
     write(LIS_logunit,*) "MSG: get_gpcp -- changing gpcp grid to 1979 -- 2000"
     write(LIS_logunit,*) "** "

     gpcp_struc(n)%ncold = 384
     gpcp_struc(n)%nrold = 190
     gpcp_struc(n)%mi = gpcp_struc(n)%ncold * gpcp_struc(n)%nrold
     gridDesci = 0
     gridDesci(1) = 4
     gridDesci(2) = 384
     gridDesci(3) = 190
     gridDesci(4) = -89.277
     gridDesci(5) = 0.0
     gridDesci(6) = 128
     gridDesci(7) = 89.277
     gridDesci(8) = -0.9375
     gridDesci(9) = 0.9375
     gridDesci(10)= 95
     gridDesci(11) = 64
     gridDesci(20)= 0

     call gpcp_reset_interp_input(n, findex, gridDesci)
     gpcp_struc(n)%gridchange1 = .false.

! 2000-2002 T170 Grid:
  elseif ( ctime >= gpcp_struc(n)%griduptime2 .and. &
           ctime < gpcp_struc(n)%griduptime3 .and. &
           gpcp_struc(n)%gridchange2 ) then

     write(LIS_logunit,*) "** "
     write(LIS_logunit,*) "MSG: get_gpcp -- changing gpcp grid to 2000 -- 2002"
     write(LIS_logunit,*) "** "

     gpcp_struc(n)%ncold = 512
     gpcp_struc(n)%nrold = 256
     gpcp_struc(n)%mi = gpcp_struc(n)%ncold * gpcp_struc(n)%nrold
     gridDesci = 0
     gridDesci(1) = 4
     gridDesci(2) = 512
     gridDesci(3) = 256
     gridDesci(4) = -89.463
     gridDesci(5) = 0
     gridDesci(6) = 128
     gridDesci(7) = 89.463
     gridDesci(8) = -0.703125
     gridDesci(9) = 0.703125
     gridDesci(10) = 128
     gridDesci(11) = 64
     gridDesci(20) = 0.0

     call gpcp_reset_interp_input(n, findex, gridDesci)
     gpcp_struc(n)%gridchange2 = .false.

! 2002-2005 T254 Grid:
  elseif ( ctime >= gpcp_struc(n)%griduptime3 .and. &
           ctime < gpcp_struc(n)%griduptime4 .and. &
           gpcp_struc(n)%gridchange3 ) then

     write(LIS_logunit,*) "** "
     write(LIS_logunit,*) "MSG: get_gpcp -- changing gpcp grid to 2002 -- 2005"
     write(LIS_logunit,*) "** "

     gpcp_struc(n)%ncold = 768
     gpcp_struc(n)%nrold = 384
     gpcp_struc(n)%mi = gpcp_struc(n)%ncold * gpcp_struc(n)%nrold
     gridDesci = 0
     gridDesci(1) = 4
     gridDesci(2) = 768
     gridDesci(3) = 384
     gridDesci(4) = -89.642
     gridDesci(5) = 0
     gridDesci(6) = 128
     gridDesci(7) = 89.642
     gridDesci(8) = -0.46875
     gridDesci(9) = 0.46875
     gridDesci(10) = 192
     gridDesci(11) = 64
     gridDesci(20) = 0.0

     call gpcp_reset_interp_input(n, findex, gridDesci)
     gpcp_struc(n)%gridchange3 = .false.

! 2005-2012 T382 Grid:
  elseif ( ctime >= gpcp_struc(n)%griduptime4 .and. &
           ctime < gpcp_struc(n)%griduptime5 .and. &
           gpcp_struc(n)%gridchange4 ) then

     write(LIS_logunit,*) "** "
     write(LIS_logunit,*) "MSG: get_gpcp -- changing gpcp grid to 2005 -- 2012 "
     write(LIS_logunit,*) "** "

     gpcp_struc(n)%ncold = 1152
     gpcp_struc(n)%nrold = 576
     gpcp_struc(n)%mi = gpcp_struc(n)%ncold * gpcp_struc(n)%nrold
     gridDesci = 0
     gridDesci(1) = 4
     gridDesci(2) = 1152
     gridDesci(3) = 576
     gridDesci(4) = -89.761
     gridDesci(5) = 0
     gridDesci(6) = 128
     gridDesci(7) = 89.761
     gridDesci(8) = -0.3125
     gridDesci(9) = 0.3125
     gridDesci(10) = 288
     gridDesci(11) = 64
     gridDesci(20) = 0.0

     call gpcp_reset_interp_input(n, findex, gridDesci)
     gpcp_struc(n)%gridchange4 = .false.

! 2010/7 -- T574 Grid:
  elseif ( ctime >= gpcp_struc(n)%griduptime5 .and. &
           ctime < gpcp_struc(n)%griduptime6 .and. &
           gpcp_struc(n)%gridchange5 ) then

     write(LIS_logunit,*) "** "
     write(LIS_logunit,*) "MSG: get_gpcp -- changing gpcp grid to 2010 -- 2015"
     write(LIS_logunit,*) "** "

     gpcp_struc(n)%ncold = 1760
     gpcp_struc(n)%nrold = 880
     gpcp_struc(n)%mi = gpcp_struc(n)%ncold * gpcp_struc(n)%nrold
     gridDesci = 0
     gridDesci(1) = 4
     gridDesci(2) = 1760
     gridDesci(3) = 880
     gridDesci(4) = -89.844
     gridDesci(5) = 0
     gridDesci(6) = 128
     gridDesci(7) = 89.844
     gridDesci(8) = -0.204545454545455
     gridDesci(9) = 0.204545454545455
     gridDesci(10) = 440
     gridDesci(11) = 64
     gridDesci(20) = 0.0

     call gpcp_reset_interp_input(n, findex, gridDesci)
     gpcp_struc(n)%gridchange5 = .false.

  elseif ( ctime >= gpcp_struc(n)%griduptime6 .and. &
           gpcp_struc(n)%gridchange6 ) then

     write(LIS_logunit,*) "** "
     write(LIS_logunit,*) "MSG: get_gpcp -- changing gpcp grid to 2015 -- "
     write(LIS_logunit,*) "** "

     gpcp_struc(n)%ncold = 3072
     gpcp_struc(n)%nrold = 1536
     gpcp_struc(n)%mi = gpcp_struc(n)%ncold * gpcp_struc(n)%nrold
     gridDesci = 0
     gridDesci(1) = 4
     gridDesci(2) = 3072
     gridDesci(3) = 1536
     gridDesci(4) = -89.910
     gridDesci(5) = 0
     gridDesci(6) = 128
     gridDesci(7) = 89.910
     gridDesci(8) = -0.1171875
     gridDesci(9) = 0.1171875
     gridDesci(10) = 768
     gridDesci(11) = 64
     gridDesci(20) = 0.0

     call gpcp_reset_interp_input(n, findex, gridDesci)
     gpcp_struc(n)%gridchange6 = .false.

  endif

!------------------------------------------------------------------------
! Ensure that data is found during first time step
!------------------------------------------------------------------------
  if( LIS_get_nstep(LIS_rc,n) == 1 .or. &
       LIS_rc%rstflag(n) == 1) then 
     LIS_rc%rstflag(n) = 0
  endif

!------------------------------------------------------------------------
! Check for and get GPCP Precipitation data
!------------------------------------------------------------------------
   filename=""
   ferror_gpcp = 0
   order = 2

 ! LIS timestep < GPCP time interval (3hr)
   if( LIS_rc%ts < gpcp_struc(n)%ts ) then
     if( ctime > gpcp_struc(n)%gpcptime ) then
      ! Put together GPCP filename:
        call gpcpfile( filename, gpcp_struc(n)%gpcpdir, yr5, mo5, da5, hr5 )
        write(LIS_logunit,*) "Getting new GPCP precip data: ",trim(filename)
        call read_gpcp( n, filename, findex, order, ferror_gpcp, hr5 )
        gpcp_struc(n)%gpcptime = ftime_gpcp
     endif  !need new time2

 ! LIS timestep >= GPCP time interval (3hr)
   elseif( LIS_rc%ts >= gpcp_struc(n)%ts ) then
    ! Put together GPCP filename:
      call gpcpfile( filename, gpcp_struc(n)%gpcpdir, yr1, mo1, da1, hr1 )
      write(LIS_logunit,*) "Getting new GPCP precip data: ",trim(filename)
      call read_gpcp( n, filename, findex, order, ferror_gpcp, hr1 )
   endif


end subroutine get_gpcp

!BOP
! !ROUTINE: gpcpfile
! \label{gpcpfile}
!
! !INTERFACE:
subroutine gpcpfile( filename, gpcpdir, yr, mo, da, hr)

  implicit none
! !ARGUMENTS: 
  character(len=*)   :: filename
  character(len=*)   :: gpcpdir
  character(len=100) :: temp
  integer            :: yr, mo, da, hr
! !DESCRIPTION:
!   This subroutine puts together GPCP file name for 
!   3 hour file intervals.
! 
!  The arguments are:
!  \begin{description}
!  \item[gpcpdir]
!    Name of the GPCP directory
!  \item[yr]
!    year 
!  \item[mo]
!   month
!  \item[da]
!   day of month
!  \item[hr]
!   hour of day
!   \item[filename]
!   file name of the timestamped GPCP file
!  \end{description}
!
!EOP

  integer :: i, c
  integer :: uyr, umo, uda, uhr, umn, uss, ts1
!  character*1 :: fbase(80), fdir(8), ftime(10), fsubs(13), fsubs2(4)
  character*1 :: fbase(80), fdir(8), ftime(10), fsubs(18), fsubs2(4)

!=== End Variable Definition ===============

!------------------------------------------------------------------------
! Make variables for the time used to create the file
! We don't want these variables being passed out
!------------------------------------------------------------------------
  uyr = yr
  umo = mo
  uda = da
  uhr = 3*(hr/3)  !hour needs to be a multiple of 3 hours
  umn = 0
  uss = 0
  ts1 = -24*60*60 !one day interval to roll back date.

  filename = ''

  write(UNIT=temp, fmt='(a40)') gpcpdir
  read(UNIT=temp, fmt='(80a1)') (fbase(i), i=1,80)

  write(UNIT=temp, fmt='(a1, i4, i2, a1)') '/', uyr, umo, '/'
  read(UNIT=temp, fmt='(8a1)') fdir
  do i = 1, 8
     if ( fdir(i) == ' ' ) fdir(i) = '0'
  end do

! disaggregated 1DD gpcp
!  write(UNIT=temp, fmt='(a13)') 'gpcp_1d_v1.2_'
!  read (UNIT=temp, fmt='(13a1)') (fsubs(i), i=1,13)
! disaggregated 1DD gpcp
  write(UNIT=temp, fmt='(a18)') 'gpcp_daily_v13rA1_'
  read (UNIT=temp, fmt='(18a1)') (fsubs(i), i=1,18)

  write(UNIT=temp, fmt='(i4, i2, i2, i2)') uyr, umo, uda, uhr
  read(UNIT=temp, fmt='(10a1)') ftime
  do i = 1, 10
     if ( ftime(i) == ' ' ) ftime(i) = '0'
  end do

  write(UNIT=temp, fmt='(a4)') '.grb'
  read (UNIT=temp, fmt='(4a1)') (fsubs2(i), i=1,4)
  c = 0
  do i = 1, 80
     if ( (fbase(i) == ' ') .and. (c == 0) ) c = i-1
  end do

  write(UNIT=temp, fmt='(80a1)') (fbase(i), i=1,c), (fdir(i), i=1,8),  &
!                       (fsubs(i), i=1,13),(ftime(i), i=1,10), &
                       (fsubs(i), i=1,18),(ftime(i), i=1,10), &
                       (fsubs2(i), i=1,4)

  read(UNIT=temp, fmt='(a80)') filename

  return
end subroutine gpcpfile

