!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: get_wfde5
! \label{get_wfde5}
!
!
! !REVISION HISTORY:
! 28 mar 2022: Hiroko Beaudoing, initial code
!
! !INTERFACE:
subroutine get_wfde5(n, findex)
! !USES:
  use LIS_coreMod
  use LIS_timeMgrMod
  use LIS_logMod
  use LIS_metforcingMod
  use wfde5_forcingMod
  use LIS_constantsMod, only: LIS_CONST_PATH_LEN

  implicit none

! !ARGUMENTS:
  integer, intent(in) :: n
  integer, intent(in) :: findex
!
! !DESCRIPTION:
!  Opens, reads, and interpolates 1-hourly WFDE5 forcing.
!
!  The WFDE5 forcing data are organized into monthly files, where each
!  file contains 24 one-hourly records over days of the month per forcing field.
!
!  In general, metforcing readers read the forcing data before the current
!  time, referred to as bookend1, and after the current time, referred to as
!  bookend2.  Then the readers temporally interpolate between bookend1 and
!  bookend2.  
!
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[findex]
!    forcing dataset index
!  \end{description}
!
!  The routines invoked are:
!  \begin{description}
!  \item[LIS\_tick](\ref{LIS_tick}) \newline
!    call to advance or retract time
!  \item[wfde5files](\ref{wfde5files}) \newline
!    Puts together appropriate file name for 1 hour intervals
!  \item[read\_wfde5](\ref{read_wfde5}) \newline
!    call to read the WFDE5 data and perform spatial interpolation
!  \end{description}
!EOP
  integer           :: order
  integer           :: ferror
  character(len=LIS_CONST_PATH_LEN),dimension(8) :: fname
  integer           :: c, r,kk,f,try
  integer           :: yr1, mo1, da1, hr1, mn1, ss1, doy1
  integer           :: yr2, mo2, da2, hr2, mn2, ss2, doy2
  real*8            :: time1, time2, timenow
  real*8            :: dtime1, dtime2
  real              :: gmt1, gmt2
  real              :: ts1, ts2

  integer           :: hr_int1, hr_int2
  integer           :: movetime  ! Flag to move bookend2 files to bookend1

! _________________________________________________________

  if( LIS_rc%nts(n).gt.3600 ) then   ! > 1-hr timestep
     write(LIS_logunit,*) '[ERR] When running LIS with WFDE5, the clock '
     write(LIS_logunit,*) '[ERR] should run with a timestep less than or '
     write(LIS_logunit,*) '[ERR] equal to one hour.'
     call LIS_endrun()
  endif

  wfde5_struc(n)%findtime1 = 0
  wfde5_struc(n)%findtime2 = 0
  movetime = 0

!=== Determine Required WFDE5 Data Times (The previous hour and the future hour)
  yr1=LIS_rc%yr
  mo1=LIS_rc%mo
  da1=LIS_rc%da
  hr1=LIS_rc%hr
  mn1=LIS_rc%mn
  ss1=0
  ts1=0
  call LIS_tick(timenow,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
 
  if(LIS_rc%ts.gt.3600) then 
     write(LIS_logunit,*) '[ERR] The model timestep is > forcing data timestep'
     write(LIS_logunit,*) '[ERR] LIS does not support this mode currently'
     write(LIS_logunit,*) '[ERR] Program stopping ...'
     call LIS_endrun()
  endif

  if(mod(nint(LIS_rc%ts),3600).eq.0) then 
     if(timenow.ge.wfde5_struc(n)%wfde5time2) then 
        yr1 = LIS_rc%yr
        mo1=LIS_rc%mo
        da1=LIS_rc%da
        hr1=LIS_rc%hr
        mn1=0
        ss1=0
        ts1=-60*60
        call LIS_tick(time1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
        
        yr2=LIS_rc%yr    !next hour
        mo2=LIS_rc%mo
        da2=LIS_rc%da
        hr2=LIS_rc%hr
        mn2=0
        ss2=0
        ts2=0
        call LIS_tick(time2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)
        movetime = 1
        wfde5_struc(n)%findtime2 = 1
     endif
  else
     if(timenow.ge.wfde5_struc(n)%wfde5time2) then 
        yr1 = LIS_rc%yr
        mo1=LIS_rc%mo
        da1=LIS_rc%da
        hr1=LIS_rc%hr
        mn1=0
        ss1=0
        ts1=0
        call LIS_tick(time1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)

        yr2=LIS_rc%yr    !next hour
        mo2=LIS_rc%mo
        da2=LIS_rc%da
        hr2=LIS_rc%hr
        mn2=0
        ss2=0
        ts2=60*60
        call LIS_tick(time2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)

        movetime = 1
        wfde5_struc(n)%findtime2 = 1
     endif
  endif

  if(LIS_rc%tscount(n).eq.1 .or.LIS_rc%rstflag(n).eq.1  ) then  
     wfde5_struc(n)%findtime1=1
     wfde5_struc(n)%findtime2=1
     movetime=0
     LIS_rc%rstflag(n) = 0
  endif
  
  if(movetime.eq.1) then
     wfde5_struc(n)%wfde5time1=wfde5_struc(n)%wfde5time2
     do f=1,LIS_rc%met_nf(findex)
        do c=1,LIS_rc%ngrid(n)
           wfde5_struc(n)%metdata1(:,f,c)=wfde5_struc(n)%metdata2(:,f,c)
        enddo
     enddo
  endif    !end of movetime=1
  
  if(wfde5_struc(n)%findtime1.eq.1) then
!=== the following looks back 10 days, at the same hour to fill data gaps.
     ferror=0
     try=0  
     ts1=-60*60*24
     do 
        if ( ferror /= 0 ) exit
        try=try+1
        do kk= wfde5_struc(n)%st_iterid, wfde5_struc(n)%en_iterid
           order = 1
           call wfde5files(n,kk,findex,wfde5_struc(n)%wfde5dir, yr1, mo1, da1, &
                wfde5_struc(n)%presrc, fname)
           call read_wfde5(n, kk,order, yr1,mo1, da1, hr1, &
                findex, fname, ferror)
        enddo

        if(ferror.ge.1) wfde5_struc(n)%wfde5time1=time1
        call LIS_tick(dtime1,doy1,gmt1,yr1,mo1,da1,hr1,mn1,ss1,ts1)
        if(try.gt.11)then
           write(LIS_logunit,*)'[ERR] WFDE5 data gap exceeds 10 days on file 1'
           call LISrun()
        endif
     enddo
!=== end of data search
  endif   

  if(wfde5_struc(n)%findtime2.eq.1) then
!=== the following looks back 10 days, at the same hour to fill data gaps.
     ferror=0
     try=0  
     ts2=-60*60*24
     do 
        if ( ferror /= 0 ) exit
        try=try+1

     !- Obtaining WFDE5 File:
        do kk= wfde5_struc(n)%st_iterid, wfde5_struc(n)%en_iterid
          order = 2
          call wfde5files(n,kk,findex,wfde5_struc(n)%wfde5dir, yr2, mo2, da2, &
               wfde5_struc(n)%presrc,fname)
          call read_wfde5(n, kk,order, yr2,mo2, da2, hr2, &
               findex, fname, ferror)
        end do

        if(ferror.ge.1) then
           wfde5_struc(n)%wfde5time2=time2
        endif
        call LIS_tick(dtime2,doy2,gmt2,yr2,mo2,da2,hr2,mn2,ss2,ts2)
        if(try.gt.11)then
           write(LIS_logunit,*)'[ERR] WFDE5 data gap exceeds 10 days on file 2'
           call LIS_endrun()
        endif
     enddo
  endif 
end subroutine get_wfde5


!BOP
! !ROUTINE: wfde5files
! \label{wfde5files}
!
! !INTERFACE:
subroutine wfde5files(n, kk, findex, wfde5dir, yr, mo, da, presrc, fname)

! !USES:
  use LIS_coreMod
  use LIS_logMod
  use LIS_forecastMod
  use LIS_timeMgrMod

  implicit none
! !ARGUMENTS:
  integer                       :: n 
  integer                       :: kk
  integer                       :: findex
  character(len=*), intent(in)  :: wfde5dir
  integer, intent(in)           :: yr,mo,da
  character(len=*), intent(in)  :: presrc
  character(len=*), dimension(8), intent(out) :: fname

! !DESCRIPTION:
!   This subroutine puts together WFDE5 file names for
!   monthly netcdf files
!
!  The arguments are:
!  \begin{description}
!  \item[wfde5dir]
!    Name of the WFDE5 directory
!  \item[yr]
!    year
!  \item[mo]
!   month
!  \item[da]
!   day of month
!  \item[fname]
!   name of the timestamped WFDE5 file
!  \end{description}
!
!EOP

  character*4  :: cyear
  character*2  :: cmonth
  character*8  :: cdate
  integer      :: hr, mn, ss
  real*8       :: time
  integer      :: doy
  real         :: gmt

  hr = 0 
  mn = 0 
  ss = 0 

  if(LIS_rc%forecastMode.eq.0) then !hindcast run
     write(unit=cyear, fmt='(i4.4)') yr
     write(unit=cmonth,fmt='(i2.2)') mo
     
     fname(1) = trim(wfde5dir)//'/'//trim(cyear)//'/PSurf_WFDE5_CRU_'//trim(cyear)//trim(cmonth)//'_v2.0.nc'
     fname(2) = trim(wfde5dir)//'/'//trim(cyear)//'/Tair_WFDE5_CRU_'//trim(cyear)//trim(cmonth)//'_v2.0.nc'
     fname(3) = trim(wfde5dir)//'/'//trim(cyear)//'/Qair_WFDE5_CRU_'//trim(cyear)//trim(cmonth)//'_v2.0.nc'
     fname(4) = trim(wfde5dir)//'/'//trim(cyear)//'/Wind_WFDE5_CRU_'//trim(cyear)//trim(cmonth)//'_v2.0.nc'
     if ( presrc.eq."GPCC" ) then
       fname(5) = trim(wfde5dir)//'/'//trim(cyear)//'/Rainf_WFDE5_CRU+GPCC_'//trim(cyear)//trim(cmonth)//'_v2.0.nc'
       fname(6) = trim(wfde5dir)//'/'//trim(cyear)//'/Snowf_WFDE5_CRU+GPCC_'//trim(cyear)//trim(cmonth)//'_v2.0.nc'
     else
       fname(5) = trim(wfde5dir)//'/'//trim(cyear)//'/Rainf_WFDE5_CRU_'//trim(cyear)//trim(cmonth)//'_v2.0.nc'
       fname(6) = trim(wfde5dir)//'/'//trim(cyear)//'/Snowf_WFDE5_CRU_'//trim(cyear)//trim(cmonth)//'_v2.0.nc'
     endif
     fname(7) = trim(wfde5dir)//'/'//trim(cyear)//'/SWdown_WFDE5_CRU_'//trim(cyear)//trim(cmonth)//'_v2.0.nc'
     fname(8) = trim(wfde5dir)//'/'//trim(cyear)//'/LWdown_WFDE5_CRU_'//trim(cyear)//trim(cmonth)//'_v2.0.nc'
    

  else !forecast mode
     !sample yr, mo, da

     call LIS_sample_forecastDate(n, kk, findex, yr,mo,da)

     write(unit=cyear, fmt='(i4.4)') yr
     write(unit=cmonth,fmt='(i2.2)') mo
     
     fname(1) = trim(wfde5dir)//'/'//trim(cyear)//'/PSurf_WFDE5_CRU_'//trim(cyear)//trim(cmonth)//'_v2.0.nc'
     fname(2) = trim(wfde5dir)//'/'//trim(cyear)//'/Tair_WFDE5_CRU_'//trim(cyear)//trim(cmonth)//'_v2.0.nc'
     fname(3) = trim(wfde5dir)//'/'//trim(cyear)//'/Qair_WFDE5_CRU_'//trim(cyear)//trim(cmonth)//'_v2.0.nc'
     fname(4) = trim(wfde5dir)//'/'//trim(cyear)//'/Wind_WFDE5_CRU_'//trim(cyear)//trim(cmonth)//'_v2.0.nc'
     if ( presrc.eq."GPCC" ) then
       fname(5) = trim(wfde5dir)//'/'//trim(cyear)//'/Rainf_WFDE5_CRU+GPCC_'//trim(cyear)//trim(cmonth)//'_v2.0.nc'
       fname(6) = trim(wfde5dir)//'/'//trim(cyear)//'/Snowf_WFDE5_CRU+GPCC_'//trim(cyear)//trim(cmonth)//'_v2.0.nc'
     else
       fname(5) = trim(wfde5dir)//'/'//trim(cyear)//'/Rainf_WFDE5_CRU_'//trim(cyear)//trim(cmonth)//'_v2.0.nc'
       fname(6) = trim(wfde5dir)//'/'//trim(cyear)//'/Snowf_WFDE5_CRU_'//trim(cyear)//trim(cmonth)//'_v2.0.nc'
     endif
     fname(7) = trim(wfde5dir)//'/'//trim(cyear)//'/SWdown_WFDE5_CRU_'//trim(cyear)//trim(cmonth)//'_v2.0.nc'
     fname(8) = trim(wfde5dir)//'/'//trim(cyear)//'/LWdown_WFDE5_CRU_'//trim(cyear)//trim(cmonth)//'_v2.0.nc'
     
  endif
end subroutine wfde5files

