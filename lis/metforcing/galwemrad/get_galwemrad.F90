!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: get_galwemrad
! \label{get_galwemrad}
!
!
! !REVISION HISTORY:
! 15 Mar 2022: Yeosang Yoon, initial code
! 08 Sep 2022; Yeosang Yoon, Add codes to read GALWEM 25 DEG dataset
! 25 May 2023: Eric Kemp, extended GALWEM 0.25 deg to 240 hours.
! 01 Jun 2023: Eric Kemp, fixes to time code to support runs other than
!              00Z, to correct time intervals between GALWEM-GD files,
!              and to stop searching for (nonexistent) GALWEM-GD file
!              when LIS reaches it's end time.
! 11 Jan 2024: Eric Kemp. Reworked code to be fault tolerant, to roll
!              back to an earlier GALWEM run if necessary, and to better
!              ensure LIS outputs history files before an early
!              termination.
! 13 Jun 2025: Hiroko Beaudoing, adipted GALWEM forecast implenemtation for 
!                                radiation supplemental forcing option.
! !INTERFACE:
subroutine get_galwemrad(n, findex)
! !USES:
  use galwemrad_forcingMod
  use LIS_constantsMod, only: LIS_CONST_PATH_LEN
  use LIS_coreMod, only: LIS_masterproc, LIS_rc, LIS_endofrun
  use LIS_logMod, only: LIS_logunit, LIS_abort, LIS_alert, LIS_endrun
  use LIS_metforcingMod
  use LIS_mpiMod
  use LIS_timeMgrMod, only: LIS_julhr_date, LIS_tick

  implicit none

! !ARGUMENTS:
  integer, intent(in) :: n
  integer, intent(in) :: findex
!
! !DESCRIPTION:
!  Opens, reads, and interpolates GALWEM forecast forcing.

!EOP
  integer           :: order,ferror
  integer           :: yr1, mo1, da1, hr1, mn1, ss1, doy1
  integer           :: yr2, mo2, da2, hr2, mn2, ss2, doy2
  real*8            :: time1, time2
  real              :: gmt1, gmt2
  real              :: ts1, ts2
  integer           :: fc_hr

  integer           :: hr_int1, hr_int2
  integer           :: valid_hour
  integer           :: fcsthr_intv
  integer           :: fcst_hour
  integer           :: openfile
  integer :: filecount
  integer :: first_fcsthr, next_fcsthr
  integer :: ierr
  logical :: lrc
  character(len=LIS_CONST_PATH_LEN) :: message(20)

  ! GALWEM radiation fields are at hourly at 0.25 degree

  if (LIS_rc%ts .gt. 3600) then
     write(LIS_logunit,*) '[ERR] The model timestep is > 1 hour ...'
     write(LIS_logunit,*) '[ERR] LIS does not support this mode currently.'
     message = ''
     message(1) = '[ERR] Program: LIS'
     message(2) = ' Routine: get_galwemrad.'
     message(3) = ' Model timestep > 1 hour'
     message(4) = ' LIS does not support this mode.'

#if (defined SPMD)
     call MPI_Barrier(LIS_MPI_COMM, ierr)
#endif
     if (LIS_masterproc) then
        call LIS_alert( 'LIS.get_galwemrad          ', 1, message)
        call LIS_abort(message)
     else
        call sleep(10) ! Make sure LIS_masterproc finishes LIS_abort
        call LIS_endrun()
     end if
  endif

  openfile=0

  if(LIS_rc%tscount(n).eq.1 .or.LIS_rc%rstflag(n).eq.1) then  !beginning of run
     LIS_rc%rstflag(n) = 0
  endif

  !EMK...Perform audit of available GALWEM-GD GRIB files for current fcst
  !run.  First two files will be saved in galwemrad_struc(n) as metdata1
  !and metdata2; also saves fcst_hour.
  if (LIS_rc%tscount(n) .eq. 1 .or. LIS_rc%rstflag(n) .eq. 1) then
     call run_audit24(n, findex, &
          LIS_rc%syr, LIS_rc%smo, LIS_rc%sda, galwemrad_struc(n)%init_hr, &
          galwemrad_struc(n)%fcst_hour, filecount)
     if (filecount < 2) then !Roll back to prior GALWEM run
        yr1 = LIS_rc%syr
        mo1 = LIS_rc%smo
        da1 = LIS_rc%sda
        hr1 = LIS_rc%shr
        mn1 = 0
        ss1 = 0
        ts1 = -43200
        call LIS_tick(time1, doy1, gmt1, &
             yr1, mo1, da1, hr1, mn1, ss1, ts1)
        call run_audit24(n, findex, &
             yr1, mo1, da1, hr1, &
             12, filecount)
        if (filecount < 2) then
           write(LIS_logunit,*) &
                "[ERR] No GALWEM files found, LIS will be halted!"
           message = ''
           message(1) = '[ERR] Program: LIS'
           message(2) = ' Routine: get_galwemrad.'
           message(3) = ' No GALWEM files found.'
           message(4) = ' LIS will be halted.'

#if (defined SPMD)
           call MPI_Barrier(LIS_MPI_COMM, ierr)
#endif
           if (LIS_masterproc) then
              call LIS_alert( 'LIS.get_galwemrad          ', 1, message)
              call LIS_abort(message)
           else
              call sleep(10) ! Make sure LIS_masterproc finishes LIS_abort
              call LIS_endrun()
           end if

        else
           galwemrad_struc(n)%init_yr = yr1
           galwemrad_struc(n)%init_mo = mo1
           galwemrad_struc(n)%init_da = da1
           galwemrad_struc(n)%init_hr = hr1
        end if
     end if
  end if

  ! Get the bookends at the start of the run
  if (LIS_rc%tscount(n) .eq. 1 .or. LIS_rc%rstflag(n) .eq. 1) then
     openfile = 1 ! Data already stored in galwemrad_struc(n) by run_audit
  end if

  ! If this is not the first LIS timestep, see if we need a new bookend
  if (openfile == 0 .and. LIS_rc%mn==0 .and. &
       LIS_rc%tscount(n) .ne. 1) then

     ! See if we need to update the bookends
     yr1 = LIS_rc%yr
     mo1 = LIS_rc%mo
     da1 = LIS_rc%da
     hr1 = LIS_rc%hr
     mn1 = 0
     ss1 = 0
     ts1 = 0
     call LIS_tick(time1, doy1, gmt1, yr1, mo1, da1, hr1, mn1, ss1, ts1)
     galwemrad_struc(n)%init_yr = yr1
     galwemrad_struc(n)%init_mo = mo1
     galwemrad_struc(n)%init_da = da1
     galwemrad_struc(n)%init_hr = 6*(hr1/6)

     if (time1 >= galwemrad_struc(n)%fcsttime2) then

        ! Update the bookends.  Note: find_next_fcsthr will save
        ! data from next file in galwemrad_struc(n)%metdata3
        call find_next_fcsthr24(n, findex, &
             galwemrad_struc(n)%init_yr, &
             galwemrad_struc(n)%init_mo, &
             galwemrad_struc(n)%init_da, &
             galwemrad_struc(n)%init_hr, &
             galwemrad_struc(n)%fcst_hour, next_fcsthr, ierr)

        if (ierr .ne. 0) then
           ! We couldn't find a newer GALWEM file.  Set flag to
           ! terminate after this timestep, and return.
           write(LIS_logunit,*) &
                '[WARN] Cannot find next GALWEM RAD file!'
           write(LIS_logunit,*) &
                '[WARN] on ',galwemrad_struc(n)%init_mo,galwemrad_struc(n)%init_da,galwemrad_struc(n)%init_hr
           flush(LIS_logunit)
           !HKB lrc = LIS_endofrun(.true.) ! Force LIS_endofrun to return true
           message = ''
           message(1) = '[WARN] Program: LIS'
           message(2) = ' Routine: get_galwemrad.'
           message(3) = ' Cannot find next GALWEM RAD file.'
           if (LIS_masterproc) then
              call LIS_alert( 'LIS.get_galwemrad          ', 1, message)
           end if
           galwemrad_struc(n)%metdata1 = LIS_rc%udef
           galwemrad_struc(n)%metdata2 = LIS_rc%udef
           galwemrad_struc(n)%fcst_hour = next_fcsthr
           galwemrad_struc(n)%fcsttime1 = galwemrad_struc(n)%fcsttime2
           galwemrad_struc(n)%fcsttime2 = galwemrad_struc(n)%fcsttime3
           return
        else
           ! We have a new second bookend.  Copy the saved older bookend
           ! to the first end.
           galwemrad_struc(n)%metdata1 = galwemrad_struc(n)%metdata2
           galwemrad_struc(n)%fcsttime1 = galwemrad_struc(n)%fcsttime2

           galwemrad_struc(n)%fcst_hour = next_fcsthr
           galwemrad_struc(n)%fcsttime2 = galwemrad_struc(n)%fcsttime3
           galwemrad_struc(n)%metdata2 = galwemrad_struc(n)%metdata3

        end if
     end if
  end if
  openfile = 0

end subroutine get_galwemrad

!BOP
!
! !ROUTINE: getgalwemradfilename
! \label{getgalwemradfilename}
!
! !INTERFACE:
subroutine getgalwemradfilename(n,rootdir,yr,mo,da,hr,fc_hr,filename)

  use galwemrad_forcingMod
  use LIS_logMod, only: LIS_logunit

  implicit none

  ! !ARGUMENTS:
  integer,          intent(in)  :: n
  character(len=*), intent(in)  :: rootdir
  integer,          intent(in)  :: yr,mo,da,hr
  integer,          intent(in)  :: fc_hr
  character(len=*), intent(out) :: filename

!EOP
  character(8) :: ftime
  character(2) :: chr
  character(3) :: fchr

  character(len=100) :: fname

  write (UNIT=chr, FMT='(i2.2)') hr
  write (UNIT=fchr, FMT='(i3.3)') fc_hr

  ! GALWEM-17km
  if(galwemrad_struc(n)%resol == 17) then
     fname = 'PS.557WW_SC.U_DI.C_GP.GALWEM-GD_GR.C17KM_AR.GLOBAL_DD.'
  endif

  ! GALWEM-25deg
  if(galwemrad_struc(n)%resol == 25) then
     fname = 'PS.557WW_SC.U_DI.C_DC.GRID_GP.GALWEM-GD_SP.COMPLEX_GR.C0P25DEG_AR.GLOBAL_PA.NASA-GODDARD_DD.'
  endif

  write (UNIT=ftime, FMT='(i4, i2.2, i2.2)') yr, mo, da

  filename = trim(rootdir) // '/' // ftime // '/' // &
             trim(fname) // ftime // '_CY.' // chr // '_FH.' // fchr // '_DF.GR2'
end subroutine getgalwemradfilename

! Run audit on GALWEM GRIB files with specified start date and time.
subroutine run_audit24(n, findex, &
     yr, mo, da, hr, first_fcsthr, filecount)

  ! Imports
  use galwemrad_forcingMod, only: galwemrad_struc
  use LIS_constantsMod, only: LIS_CONST_PATH_LEN
  use LIS_logMod, only: LIS_logunit
  use LIS_timeMgrMod, only: LIS_tick

  ! Defaults
  implicit none

  ! Arguments
  integer, intent(in) :: n
  integer, intent(in) :: findex
  integer, intent(in) :: yr
  integer, intent(in) :: mo
  integer, intent(in) :: da
  integer, intent(in) :: hr
  integer, intent(in) :: first_fcsthr
  integer, intent(out) :: filecount

  ! Locals
  integer, parameter :: MAX_FCSTHRS = 5
  character(len=LIS_CONST_PATH_LEN) :: fname
  integer :: fcsthr, delta, ferror, order
  integer :: yr1, mo1, da1, hr1, mn1, ss1, doy1
  real :: ts1, gmt1
  real*8 :: time1
  logical :: found_inq
  integer :: i

  ! Loop through list of forecast hours, and see how many files exist
  delta = 1
  filecount = 0
  do i = 0, MAX_FCSTHRS
     fcsthr = first_fcsthr + (i * delta)
     yr1 = yr
     mo1 = mo
     da1 = da
     hr1 = hr
     call getgalwemradfilename(n, galwemrad_struc(n)%odir, yr1, mo1, da1, hr1, &
          fcsthr, fname)

     ferror = 0
     inquire(file=trim(fname), exist=found_inq)
     if (found_inq) then
        if (i == 0) then
           order = 1
        else if (filecount < 2) then
           order = 2
        else
           order = 3
        end if
        call read_galwemrad(n, findex, order, fname, ferror)
        if (ferror == 0) then
           mn1 = 0
           ss1 = 0
           ts1 = fcsthr * 3600
           call LIS_tick(time1, doy1, gmt1, &
                yr1, mo1, da1, hr1, mn1, ss1, ts1)
           if (order == 1) then
              galwemrad_struc(n)%fcsttime1 = time1
              write(LIS_logunit,*) &
                   '[INFO] Read ', trim(fname)
           else if (order == 2) then
              galwemrad_struc(n)%fcsttime2 = time1
              galwemrad_struc(n)%fcst_hour = fcsthr
              write(LIS_logunit,*) &
                   '[INFO] Read ', trim(fname)
           else
              galwemrad_struc(n)%fcsttime3 = time1
           end if
        end if
     end if
     if (found_inq .and. ferror .eq. 0) then
        filecount = filecount + 1
     else
        if (i == 0) then
           write(LIS_logunit,*) &
                '[WARN] Cannot find GALWEM RAD file for LIS start date'
           write(LIS_logunit,*) &
                '[WARN] Will need to use earlier GALWEM run for forcing'
           filecount = 0
           exit
        end if
     end if

     if (filecount >= 2) exit
  end do

  !if (filecount > 0) then
  !   write(LIS_logunit,*)'[INFO] Found ', filecount, ' files'
  !end if

end subroutine run_audit24

! Find next GALWEM grib file given current selected time and
! forecast hour
!hour 00 look for CY.00_FH.000
!hour 01 look for CY.00_FH.001
!hour 02 look for CY.00_FH.002
!hour 03 look for CY.00_FH.003
!hour 04 look for CY.00_FH.004
!hour 05 look for CY.00_FH.005
!hour 06 look for CY.06_FH.000
!hour 07 look for CY.06_FH.001
!hour 08 look for CY.06_FH.002
!hour 09 look for CY.06_FH.003
!hour 10 look for CY.06_FH.004
!hour 11 look for CY.06_FH.005
!hour 12 look for CY.12_FH.000
!hour 13 look for CY.12_FH.001
!hour 14 look for CY.12_FH.002
!hour 15 look for CY.12_FH.003
!hour 16 look for CY.12_FH.004
!hour 17 look for CY.12_FH.005
!hour 18 look for CY.18_FH.000
!hour 19 look for CY.18_FH.001
!hour 20 look for CY.18_FH.002
!hour 21 look for CY.18_FH.003
!hour 22 look for CY.18_FH.004
!hour 23 look for CY.18_FH.005
subroutine find_next_fcsthr24(n, findex, &
     yr, mo, da, hr, cur_fcsthr, &
     next_fcsthr, ierr)

  ! Imports
  use galwemrad_forcingMod, only: galwemrad_struc
  use LIS_constantsMod, only: LIS_CONST_PATH_LEN
  use LIS_logMod, only: LIS_logunit
  use LIS_timeMgrMod, only: LIS_tick

  ! Defaults
  implicit none

  ! Arguments
  integer, intent(in) :: n
  integer, intent(in) :: findex
  integer, intent(in) :: yr
  integer, intent(in) :: mo
  integer, intent(in) :: da
  integer, intent(in) :: hr
  integer, intent(in) :: cur_fcsthr
  integer, intent(out) :: next_fcsthr
  integer, intent(out) :: ierr

  ! Locals
  integer, parameter :: MAX_FCSTHRS = 1
  character(len=LIS_CONST_PATH_LEN) :: fname
  integer :: fcsthr, delta, order, ferror, inithr
  integer :: yr1, mo1, da1, hr1, mn1, ss1, doy1
  integer :: yr2, mo2, da2, hr2, mn2, ss2, doy2
  real :: gmt1, ts1, gmt2, ts2
  real*8 :: time1, time2
  logical :: found_inq
  integer :: i

  ! Loop through list of forecast hours, and see how many files exist
  ierr = 1
  delta = 1
  do i = 1, MAX_FCSTHRS
     yr1 = yr
     mo1 = mo
     da1 = da
     hr1 = hr
     mn1 = 0
     ss1 = 0
     inithr = 6*(hr/6)
     if (cur_fcsthr .ge. 5 ) then
      select case (inithr)
        case (6)
         fcsthr = cur_fcsthr - inithr + i
        case (12)
         fcsthr = cur_fcsthr - inithr + 6 + i
        case (18)
         fcsthr = cur_fcsthr - inithr + 12 + i
        case default
         fcsthr = 0
      end select
     else
      fcsthr = cur_fcsthr + i
     endif
     call getgalwemradfilename(n, galwemrad_struc(n)%odir, yr1, mo1, da1, inithr, &
          fcsthr, fname)
     inquire(file=trim(fname), exist=found_inq)
     ferror = 0
     if (found_inq) then
        order = 3
        call read_galwemrad(n, findex, order, fname, ferror)
        if (ferror == 0) then
           ts1 = fcsthr * 3600
           call LIS_tick(time1, doy1, gmt1, &
                yr1, mo1, da1, hr1, mn1, ss1, ts1)
           galwemrad_struc(n)%fcsttime3 = time1
        end if
     end if
     if (found_inq .and. ferror .eq. 0) then
        write(LIS_logunit,*)'[INFO] Read ', trim(fname)
        next_fcsthr = fcsthr
        ierr = 0
        exit
     end if
  end do

end subroutine find_next_fcsthr24
