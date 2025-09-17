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
!BOP
!
! !ROUTINE: read_galwemrad
! \label{read_galwemrad}
!
! !REVISION HISTORY:
! 15 Mar 2022: Yeosang Yoon, initial code
! 11 Jan 2024; Eric Kemp, revisions for fault tolerance.
! 13 Jun 2025: Hiroko Beaudoing, adipted GALWEM forecast implenemtation for
!                                radiation supplemental forcing option.
!
! !INTERFACE:
subroutine read_galwemrad(n, findex, order, gribfile, rc)

  ! !USES:
  use LIS_constantsMod,  only : LIS_CONST_PATH_LEN
  use LIS_coreMod,       only : LIS_rc
  use LIS_logMod
  use galwemrad_forcingMod, only : galwemrad_struc

#if (defined USE_GRIBAPI)
  use grib_api
#endif

  implicit none

!ARGUMENTS:
  integer,           intent(in)    :: n
  integer,           intent(in)    :: findex  ! Forcing index
  integer,           intent(in)    :: order
  character(len=*),  intent(in)    :: gribfile

!DESCRIPTION:

!EOP
  integer         :: ftn, igrib, ierr
  integer         :: center
  character*100   :: gtype
  integer         :: file_julhr
  integer         :: yr1, mo1, da1, hr1
  integer         :: iginfo      ( 40 )
  real            :: gridres_dlat, gridres_dlon
  integer         :: ifguess, jfguess
  integer         :: dataDate, dataTime
  integer         :: fc_hr

  real            :: swdown(LIS_rc%lnc(n),LIS_rc%lnr(n))    !Downward shortwave flux at the ground [W/m^2] 
  real            :: lwdown(LIS_rc%lnc(n),LIS_rc%lnr(n))    !Downward longwave radiation at the ground [W/m^2] 
  integer, intent(out) :: rc
  logical :: found_inq

   ! Initialize return code to "no error".  We will change it below if
   ! necessary.
   rc = 0

   ! Open first guess grib data using library utility.  Just read
   ! the first file only, as all data will be of the same type.
   ! (GALWEM) because the search script ensures that it is.
   !

   ! EMK...Before using ECCODES/GRIB_API, see if the GRIB file exists
   ! using a simple inquire statement.  This avoids ECCODES/GRIB_API
   ! writing error messages to stdout/stderr, which may lead to runtime
   ! problems.
   inquire(file=trim(gribfile),exist=found_inq)
   if (.not. found_inq) then
      write(LIS_logunit,*)'[ERR] Cannot find file '//trim(gribfile)
      rc = 1
      return
   end if

#if (defined USE_GRIBAPI)

   call grib_open_file(ftn,trim(gribfile),'r',ierr)
   if ( ierr .ne. 0 ) then
      write(LIS_logunit,*) '[ERR] Failed to open - '//trim(gribfile)
      rc = 1
      return
   end if

   ! Read in the first grib record, unpack the header and extract
   ! section 1 and section 2 information.
   call grib_new_from_file(ftn,igrib,ierr)
   if ( ierr .ne. 0 ) then
      write(LIS_logunit,*) '[ERR] failed to read - '//trim(gribfile)
      call grib_close_file(ftn)
      rc = 1
      return
   endif

   call grib_get(igrib,'centre',center,ierr)
   if ( ierr .ne. 0 ) then
      write(LIS_logunit,*) '[ERR] in grib_get: ' // &
           'centre in read_galwemrad'
      call grib_release(igrib,ierr)
      call grib_close_file(ftn)
      rc = 1
      return
   endif

   call grib_get(igrib,'gridType',gtype,ierr)
   if ( ierr .ne. 0 ) then
      write(LIS_logunit,*) '[ERR] in grid_get: ' // &
           'gridtype in read_galwemrad'
      call grib_release(igrib,ierr)
      call grib_close_file(ftn)
      rc = 1
      return
   endif

   if(trim(gtype).ne."regular_ll") then
      write(LIS_logunit,*)'[ERR] GALWEM file not on lat/lon grid!'
      call grib_release(igrib,ierr)
      call grib_close_file(ftn)
      rc = 1
      return
   endif

   call grib_get(igrib,'Ni',iginfo(1),ierr)
   if ( ierr .ne. 0 ) then
      write(LIS_logunit,*) '[ERR] in grid_get: Ni in read_galwemrad'
      call grib_release(igrib,ierr)
      call grib_close_file(ftn)
      rc = 1
      return
   endif

   call grib_get(igrib,'Nj',iginfo(2),ierr)
   if ( ierr .ne. 0 ) then
      write(LIS_logunit,*) '[ERR] in grid_get: Nj in read_galwemrad'
      call grib_release(igrib,ierr)
      call grib_close_file(ftn)
      rc = 1
      return
   endif

   call grib_get(igrib,'jDirectionIncrementInDegrees',gridres_dlat,ierr)
   if ( ierr .ne. 0 ) then
      write(LIS_logunit,*) '[ERR] in grid_get: ' // &
           'jDirectionIncrementInDegrees ' // &
           'in read_galwemrad'
      call grib_release(igrib,ierr)
      call grib_close_file(ftn)
      rc = 1
      return
   endif

   ! EMK...Added dlon
   call grib_get(igrib, 'iDirectionIncrementInDegrees', gridres_dlon, ierr)
   if ( ierr .ne. 0 ) then
      write(LIS_logunit,*) '[ERR] in grid_get: ' // &
           'iDirectionIncrementInDegrees ' // &
           'in read_galwemrad'
      call grib_release(igrib, ierr)
      call grib_close_file(ftn)
      rc = 1
      return
   endif

   call grib_get(igrib,'dataDate',dataDate,ierr)
   if ( ierr .ne. 0 ) then
      write(LIS_logunit,*) '[ERR] in grid_get: ' // &
           'dataDate in read_galwemrad'
      call grib_release(igrib,ierr)
      call grib_close_file(ftn)
      rc = 1
      return
   endif

   call grib_get(igrib,'dataTime',dataTime,ierr)
   if ( ierr .ne. 0 ) then
      write(LIS_logunit,*) '[ERR] in grid_get: ' // &
           'dataTime in read_galwemrad'
      call grib_release(igrib,ierr)
      call grib_close_file(ftn)
      rc = 1
      return
   endif

   ! Here we tentatively have a file we can use. Close it for now, and
   ! prepare to pull the appropriate variables.
   call grib_release(igrib,ierr)
   call grib_close_file(ftn)

   ifguess = iginfo(1)
   jfguess = iginfo(2)

   call fldbld_read_galwemrad(n, findex, order, gribfile, ifguess, jfguess, &
            swdown, lwdown, rc)
   if (rc .ne. 0) return

   call assign_processed_galwemradforc(n,order,1,swdown)
   call assign_processed_galwemradforc(n,order,2,lwdown)

#endif

end subroutine read_galwemrad

subroutine fldbld_read_galwemrad(n, findex, order, gribfile, ifguess, jfguess,&
                              swdown, lwdown, rc)
 
  ! !USES:
  use LIS_constantsMod, only: LIS_CONST_PATH_LEN
  use LIS_coreMod, only : LIS_rc
  use LIS_logMod,  only : LIS_logunit, LIS_abort, LIS_alert, LIS_verify

#if (defined USE_GRIBAPI)
  use grib_api
#endif

  implicit none
! !ARGUMENTS:
  integer,        intent(in)    :: n
  integer, intent(in)           :: findex  ! Forcing index
  integer,        intent(in)    :: order
  character(len=*),  intent(in) :: gribfile
  integer,        intent(in)    :: ifguess
  integer,        intent(in)    :: jfguess
  real,           intent(out)   :: swdown(LIS_rc%lnc(n),LIS_rc%lnr(n))    !Downward shortwave flux at the ground [W/m^2]
  real,           intent(out)   :: lwdown(LIS_rc%lnc(n),LIS_rc%lnr(n))    !Downward longwave radiation at the ground [W/m^2]
  integer,        intent(out)   :: rc
!
! !DESCRIPTION:
!
!     To read UK Unified Model (GALWEM) data in GRIB-2 format.
!
!EOP
  character*9                   :: cstat
  character(len=LIS_CONST_PATH_LEN) :: message     ( 20 )
  character(len=7)              :: grib_msg
  character(len=7)              :: check_galwemrad_message
  integer                       :: count_swdown, count_lwdown
  integer                       :: ierr
  integer                       :: istat1
  integer                       :: igrib
  integer                       :: ftn
  integer                       :: kk, nvars
  integer                       :: prod_def_tmpl_num
  integer                       :: surface_val_2
  integer                       :: param_disc_val, param_cat_val, &
                                   param_num_val, surface_val, level_val
  real, allocatable             :: dum1d   ( : )

  real, allocatable             :: fg_swdown  ( : , : )
  real, allocatable             :: fg_lwdown  ( : , : )

  logical                       :: found_inq

  ! Executable code begins here ...

  rc = 0 ! Initialize as "no error"

  ! EMK...Before using ECCODES/GRIB_API, see if the GRIB file exists
  ! using a simple inquire statement.  This avoids ECCODES/GRIB_API
  ! writing error messages to stdout/stderr, which may lead to runtime
  ! problems.
  inquire(file=trim(gribfile),exist=found_inq)
  if (.not. found_inq) then
     write(LIS_logunit,*) '[WARN] Cannot find file '//trim(gribfile)
     rc = 1
     return
  end if

#if (defined USE_GRIBAPI)

  ! If a problem occurs here, we can just return immediately since no
  ! memory has been allocated yet.
  call grib_open_file(ftn,trim(gribfile),'r',ierr)
  if ( ierr .ne. 0 ) then
     write(LIS_logunit,*) '[WARN] Failed to open - '//trim(gribfile)
     rc = 1
     return
  end if

  allocate ( fg_swdown  (ifguess, jfguess) )
  allocate ( fg_lwdown  (ifguess, jfguess) )

  allocate ( dum1d   (ifguess*jfguess) )

  ! From this point, we must deallocate memory before returning.
  ! Unfortunately this means using a GOTO statement if a problem is
  ! encountered, but such is life.
  count_swdown  = 0
  count_lwdown = 0

  call grib_count_in_file(ftn,nvars,ierr)
  if ( ierr .ne. 0 ) then
     write(LIS_logunit,*) '[WARN] in grib_count_in_file in ' // &
          'fldbld_read_galwemrad'
     goto 100
  end if

  ! Tentatively loop through every field in GRIB file looking for the variables
  ! we want. The code below will exit the loop early if a problem is found *or*
  ! once all the required variables are found and read in.
  do kk=1,nvars

     call grib_new_from_file(ftn,igrib,ierr)
     if ( ierr .ne. 0 ) then
        write(LIS_logunit,*) '[WARN] failed to read - '//trim(gribfile)
        goto 100
     end if

     call grib_get(igrib,'discipline',param_disc_val,ierr)
     if ( ierr .ne. 0 ) then
        write(LIS_logunit,*) '[WARN] in grib_get: parameterNumber in ' // &
             'fldbld_read_galwemrad'
        call grib_release(igrib,ierr)
        goto 100
     end if

     ! EMK...Now read GRIB2 Product Definition Template Number to
     ! ensure we have an instantaneous variable at a horizontal level.
     call grib_get(igrib,'productDefinitionTemplateNumber', &
          prod_def_tmpl_num,ierr)
     if ( ierr .ne. 0 ) then
        write(LIS_logunit,*) &
             '[WARN] in grib_get: productDefinitionTemplateNumber in ' // &
             'fldbld_read_galwemrad'
        call grib_release(igrib,ierr)
        goto 100
     end if

     call grib_get(igrib,'parameterCategory',param_cat_val,ierr)
     if ( ierr .ne. 0 ) then
        write(LIS_logunit,*) &
             '[WARN] in grib_get: parameterCategory in ' // &
             'fldbld_read_galwemrad'
        call grib_release(igrib,ierr)
        goto 100
     end if

     call grib_get(igrib,'parameterNumber',param_num_val,ierr)
     if ( ierr .ne. 0 ) then
        write(LIS_logunit,*) &
             '[WARN] in grib_get: parameterNumber in ' // &
             'fldbld_read_galwemrad'
        call grib_release(igrib,ierr)
        goto 100
     end if

     call grib_get(igrib,'typeOfFirstFixedSurface',surface_val,ierr)
     if ( ierr .ne. 0 ) then
        write(LIS_logunit,*) &
             '[WARN] in grib_get: level in ' // &
             'fldbld_read_galwemrad'
        call grib_release(igrib,ierr)
        goto 100
     end if

     call grib_get(igrib,'level',level_val,ierr)
     if ( ierr .ne. 0 ) then
        write(LIS_logunit,*) &
             '[WARN] in grib_get: level in ' // &
             'fldbld_read_galwemrad'
        call grib_release(igrib,ierr)
        goto 100
     end if

     ! EMK...Now read GRIB2 type of second level to ensure we are not
     ! working with a layer-averaged field.
     call grib_get(igrib,'typeOfSecondFixedSurface',surface_val_2,ierr)
     if ( ierr .ne. 0 ) then
        write(LIS_logunit,*) &
             '[WARN] in grib_get: level in ' // &
             'fldbld_read_galwemrad'
        call grib_release(igrib,ierr)
        goto 100
     end if

     ! We have enough information to determine what GRIB parameter this
     ! is.
     grib_msg = check_galwemrad_message(param_disc_val, prod_def_tmpl_num, &
          param_cat_val, param_num_val, surface_val, level_val, surface_val_2)

     ! Skip this field if GRIB parameter is not required.
     if (grib_msg == 'none') then
        call grib_release(igrib,ierr)
        if (ierr .ne. 0) then
           write(LIS_logunit,*)'[WARN], in grib_release: in ' //&
                'fldbld_read_galwemrad'
           goto 100
        end if
        cycle ! Not a message we are interested in.
     end if

     call grib_get(igrib,'values',dum1d,ierr)
     if ( ierr .ne. 0 ) then
        write(LIS_logunit,*) &
             '[WARN] in grib_get: values in ' // &
             'fldbld_read_galwemrad'
        call grib_release(igrib,ierr)
        goto 100
     end if

     select case (grib_msg)
     case ('swdown')  ! downward shortwave
        fg_swdown = reshape(dum1d, (/ifguess,jfguess/))
        count_swdown = count_swdown + 1
     case ('lwdown')  ! downward longwave radiation
        fg_lwdown = reshape(dum1d, (/ifguess,jfguess/))
        count_lwdown = count_lwdown + 1

     case default ! Internal error, we shouldn't be here
        write(LIS_logunit,*)'[WARN] Unknown grib_message ',grib_msg
        write(LIS_logunit,*)'Aborting...'
        flush(LIS_logunit)
        write(cstat,'(i9)',iostat=istat1) ierr
        message(1) = 'Program: LIS'
        message(2) = '  Subroutine:  fldbld_read_galwemrad.'
        message(3) = '  Error reading first guess file:'
        message(4) = '  ' // trim(gribfile)
        if( istat1 .eq. 0 )then
           message(5) = '  Status = ' // trim(cstat)
        endif
        call LIS_abort( message)
     end select

     ! Finished with this field
     call grib_release(igrib,ierr)
     if (ierr .ne. 0) then
        write(LIS_logunit,*)'[WARN], in grib_release: in ' //&
             'fldbld_read_galwemrad'
        goto 100
     end if

     ! Jump out of loop early if we have everything
     if ( (count_swdown .eq. 1 ) .and. (count_lwdown  .eq. 1) ) then 
        exit
     end if

  enddo ! Loop through all GRIB file fields

  ! Make sure we have all required fields.
     if ( (count_swdown .ne. 1 ) .or. (count_lwdown  .ne. 1) ) then
     write(LIS_logunit,*)'[WARN] Missing data from GALWEM GRIB file!'
     goto 100
  end if

  ! Interpolate the fields to the LIS grid
  ! swdown
  call interp_galwemrad(n, findex, ifguess, jfguess, .false., fg_swdown, swdown)
  ! lwdown
  call interp_galwemrad(n, findex, ifguess, jfguess, .false., fg_lwdown, lwdown)

  ! At this point, we have everything.  Close the file and return.
  call grib_close_file(ftn)
  rc = 0
  return

  ! Jump down here to clean up memory before returning after finding a
  ! problem.
  100 continue
  call grib_close_file(ftn)

  deallocate ( dum1d      )
  deallocate ( fg_swdown  )
  deallocate ( fg_lwdown  )
  rc = 1
#endif

end subroutine fldbld_read_galwemrad

function check_galwemrad_message(param_disc_val, prod_def_tmpl_num, &
     param_cat_val, param_num_val, surface_val, level_val, surface_val_2)
! !USES:
! none

   implicit none
! !ARGUMENTS:
   integer, intent(in) :: param_disc_val, prod_def_tmpl_num, &
        param_cat_val, &
        param_num_val, surface_val, level_val, surface_val_2
   character(len=7)    :: check_galwemrad_message
!EOP


   ! EMK...Only use instantaneous variables
   if (prod_def_tmpl_num .ne. 0) then
      check_galwemrad_message = 'none'
      return
   end if
   ! EMK...Only use single level fields, not layers
   if (surface_val_2 .ne. 255) then
      check_galwemrad_message = 'none'
      return
   end if

   if     ( param_disc_val == 0 .and. &
            param_cat_val  == 0 .and. &
            param_num_val  == 0 .and. &
            surface_val    == 103 .and. &
            level_val == 2) then
      check_galwemrad_message = 't2'      ! Temperature interpolated to 2 metres [K] 
   elseif ( param_disc_val == 0 .and. &
            param_cat_val  == 1 .and. &
            param_num_val  == 0 .and. &
            surface_val    == 103 .and. &
            level_val == 2) then
      check_galwemrad_message = 'q2'      ! Instantaneous specific humidity interpolated to 2 metres [kg/kg]
   elseif ( param_disc_val == 0 .and. &
            param_cat_val  == 4 .and. &
            param_num_val  == 7 .and. &
            surface_val    == 1 ) then
      check_galwemrad_message = 'swdown'  ! Downward shortwave flux at the ground [W/m^2]
   elseif ( param_disc_val == 0 .and. &
            param_cat_val  == 5 .and. &
            param_num_val  == 3 .and. &
            surface_val    == 1 ) then
      check_galwemrad_message = 'lwdown'  ! Downward longwave radiation at the ground [W/m^2] 
   elseif ( param_disc_val == 0 .and. &
            param_cat_val  == 2 .and. &
            param_num_val  == 2 .and. &
            surface_val    == 103 .and. &
            level_val == 10) then
      check_galwemrad_message = 'u10'     ! Instantaneous zonal wind interpolated to 10 metres [m/s] 
   elseif ( param_disc_val == 0 .and. &
            param_cat_val  == 2 .and. &
            param_num_val  == 3 .and. &
            surface_val    == 103 .and. &
            level_val == 10) then
      check_galwemrad_message = 'v10'     ! Instantaneous meridional wind interpolated to 10 metres[m/s]
   elseif ( param_disc_val == 0 .and. &
            param_cat_val  == 3 .and. &
            param_num_val  == 0 .and. &
            surface_val    == 1 ) then
      check_galwemrad_message = 'ps'      ! Instantaneous Surface Pressure [Pa]
   elseif ( param_disc_val == 0 .and. &
            param_cat_val  == 1 .and. &
            param_num_val  == 52 .and. &
            surface_val    == 1 ) then
      check_galwemrad_message = 'prectot' ! Total precipitation [kg/m2/s]
   else
      check_galwemrad_message = 'none'
   endif
end function check_galwemrad_message

subroutine interp_galwemrad(n, findex, ifguess, jfguess, pcp_flag, input, output)

! !USES:
  use LIS_coreMod,       only : LIS_rc, LIS_domain
  use LIS_logMod,        only : LIS_logunit, LIS_endrun
  use galwemrad_forcingMod, only : galwemrad_struc

  implicit none
! !ARGUMENTS:
  integer, intent(in)  :: n
  integer, intent(in)   :: findex
  integer, intent(in)  :: ifguess
  integer, intent(in)  :: jfguess
  logical, intent(in)  :: pcp_flag
  real,    intent(in)  :: input  ( ifguess,jfguess )
  real,    intent(out) :: output ( LIS_rc%lnc(n),LIS_rc%lnr(n) )
!
! !DESCRIPTION:
!
! This routine interpolates the GALWEM data to the LIS grid.

!EOP

  integer   :: mi, mo
  integer   :: k
  integer   :: i,j
  integer   :: iret
  integer   :: midway
  character(len=50) :: method
  real, allocatable, dimension(:,:)    :: var
  logical*1, allocatable, dimension(:) :: lb
  logical*1, allocatable, dimension(:) :: lo

  allocate(var(ifguess,jfguess))
  allocate(lb(ifguess*jfguess))
  allocate(lo(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

  mi = ifguess * jfguess
  mo = LIS_rc%lnc(n)*LIS_rc%lnr(n)
  lb = .true.

  ! translate from (0,360) to (-180,180).
  midway = ifguess/2
  do j = 1, jfguess
     do i = 1, ifguess
        if ( (i+midway) < ifguess ) then
           var(i,j) = input((i+midway),j)
        else
           var(i,j) = input((i-midway+1),j)
        endif
     enddo
  enddo

  ! Interpolate to LIS grid
  select case( LIS_rc%met_interp(findex) )

    case( "bilinear" )
     call bilinear_interp(LIS_rc%gridDesc(n,:),lb,                     &
          var,lo,output,mi,mo,                                         &
          LIS_domain(n)%lat,LIS_domain(n)%lon,                         &
          galwemrad_struc(n)%w111,galwemrad_struc(n)%w121,                   &
          galwemrad_struc(n)%w211,galwemrad_struc(n)%w221,                   &
          galwemrad_struc(n)%n111,galwemrad_struc(n)%n121,                   &
          galwemrad_struc(n)%n211,galwemrad_struc(n)%n221,LIS_rc%udef,iret)

    case( "budget-bilinear" )
     if (pcp_flag) then
        call conserv_interp(LIS_rc%gridDesc(n,:),lb,                   &
             var,lo,output,mi,mo,                                      &
             LIS_domain(n)%lat, LIS_domain(n)%lon,                     &
             galwemrad_struc(n)%w112,galwemrad_struc(n)%w122,                &
             galwemrad_struc(n)%w212,galwemrad_struc(n)%w222,                &
             galwemrad_struc(n)%n112,galwemrad_struc(n)%n122,                &
             galwemrad_struc(n)%n212,galwemrad_struc(n)%n222,LIS_rc%udef,iret)
     else
        call bilinear_interp(LIS_rc%gridDesc(n,:),lb,                  &
             var,lo,output,mi,mo,                                      &
             LIS_domain(n)%lat, LIS_domain(n)%lon,                     &
             galwemrad_struc(n)%w111,galwemrad_struc(n)%w121,                &
             galwemrad_struc(n)%w211,galwemrad_struc(n)%w221,                &
             galwemrad_struc(n)%n111,galwemrad_struc(n)%n121,                &
             galwemrad_struc(n)%n211,galwemrad_struc(n)%n221,LIS_rc%udef,iret)
     endif

  case( "neighbor" )
     call neighbor_interp(LIS_rc%gridDesc(n,:),lb,                     &
          var,lo,output,mi,mo,                                         &
          LIS_domain(n)%lat, LIS_domain(n)%lon,                        &
          galwemrad_struc(n)%n113,LIS_rc%udef,iret)

  case DEFAULT
     write(LIS_logunit,*) 'ERR: Unexpected interpolation method'
     write(LIS_logunit,*) '     in interp_galwemrad_first_guess'
     write(LIS_logunit,*) '     ', trim(LIS_rc%met_interp(findex))
     call LIS_endrun()  
  end select

  deallocate(var)
  deallocate(lb)
  deallocate(lo)

end subroutine interp_galwemrad

!BOP
!
! !ROUTINE: assign_processed_galwemradforc
! \label{assign_processed_galwemradforc}
!
! !INTERFACE:
subroutine assign_processed_galwemradforc(n,order,var_index,galwemradforc)
! !USES:
  use LIS_coreMod
  use galwemrad_forcingMod, only : galwemrad_struc
!
! !DESCRIPTION:
!  This routine assigns the interpolated GALWEM forcing data
!  to the module data structures to be used later for
!  time interpolation
!
!EOP
  implicit none

  integer :: n
  integer :: order
  integer :: var_index
  real    :: galwemradforc(LIS_rc%lnc(n),LIS_rc%lnr(n))
  integer :: c,r

  do r=1,LIS_rc%lnr(n)
     do c=1,LIS_rc%lnc(n)
        if(LIS_domain(n)%gindex(c,r).ne.-1) then
           if(order.eq.1) then
              galwemrad_struc(n)%metdata1(var_index,&
                   LIS_domain(n)%gindex(c,r)) = &
                   galwemradforc(c,r)
           elseif(order.eq.2) then
              galwemrad_struc(n)%metdata2(var_index,&
                   LIS_domain(n)%gindex(c,r)) = &
                   galwemradforc(c,r)
           elseif(order.eq.3) then
              galwemrad_struc(n)%metdata3(var_index,&
                   LIS_domain(n)%gindex(c,r)) = &
                   galwemradforc(c,r)
           endif
        endif
     enddo
  enddo
end subroutine assign_processed_galwemradforc

