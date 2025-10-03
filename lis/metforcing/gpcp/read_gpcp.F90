!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
!
! !ROUTINE: read_gpcp
! \label{read_gpcp}
!
! !REVISION HISTORY:
!  06 Jun 2016: Hiroko Beaudoing; Adopted CMAP reader to GPCP
!  
! !INTERFACE:
subroutine read_gpcp( n, fname, findex, order, ferror_gpcp, filehr)
! !USES:
  use LIS_coreMod,       only : LIS_rc, LIS_domain
  use LIS_logMod,        only : LIS_logunit, LIS_warning, LIS_verify
  use gpcp_forcingMod,   only : gpcp_struc
  use LIS_constantsMod,  only : LIS_CONST_PATH_LEN
#if (defined USE_GRIBAPI)
  use grib_api
#endif

  implicit none
! !ARGUMENTS:
  integer, intent(in) :: n
  character(len=LIS_CONST_PATH_LEN)   :: fname          
  integer, intent(in) :: findex
  integer, intent(in) :: order
  integer             :: ferror_gpcp
  integer             :: filehr
!
! !DESCRIPTION:
!  For the given time, reads parameters from
!  GPCP data and interpolates to the LIS domain.
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[fname]
!    name of the 3 hour GPCP file
!  \item[ferror\_gpcp]
!    flag to indicate success of the call (=0 indicates success)
!  \item[filehr]
!    current file hour
!  \end{description}
! 
!  The routines invoked are: 
!  \begin{description}
!  \item[interp\_gpcp](\ref{interp_gpcp}) \newline
!    spatially interpolates the GPCP data
!  \end{description}
!EOP

  integer                :: ftn
  real                   :: precip_regrid(LIS_rc%lnc(n),LIS_rc%lnr(n))
  integer                :: ngpcp
  real, allocatable      :: gpcpin(:)
  logical*1,allocatable  :: lb(:)
  integer                :: index1
  integer                :: grid_size
  integer                :: i,j,t
  logical                :: file_exists
  integer                :: igrib
  integer                :: iv, iv_total
  integer                :: pds5_val, pds7_val
  integer                :: pds5(1), pds7(1)
  logical                :: var_status(1)
  logical                :: var_found
  integer                :: kk,var_index
  real                   :: missingValue
  integer                :: nvars
  integer                :: iret,rc,status
  
!=== End Variable Definition =======================

#if(defined USE_GRIBAPI) 
  pds5 = (/ 059 /) !parameter
  pds7 = (/ 000 /) !htlev2

  iv_total = 1
  ferror_gpcp = 1 
  inquire (file=trim(fname), exist=file_exists)
! File exists:
  if (file_exists) then   
     
   ! Open the GPCP grib file:
     call grib_open_file(ftn,trim(fname),'r',iret)
     if(iret.ne.0) then 
        write(LIS_logunit,*) &
             'Could not open file: ',trim(fname)
        ferror_gpcp = 0
        return
     endif
     call grib_count_in_file(ftn,nvars,iret)
     call LIS_verify(iret, 'error in grib_count_in_file in read_gpcp')
     if(iret.ne.0) then
        ferror_gpcp = 0
        return
     endif
     
     ngpcp = gpcp_struc(n)%ncold*gpcp_struc(n)%nrold
     allocate(gpcpin(ngpcp))
     allocate(lb(ngpcp)) 

   ! Search for appropiate GPCP PPT variable:
     do kk=1,nvars
        call grib_new_from_file(ftn,igrib,iret)
        call LIS_warning(iret,'error in grib_new_from_file in read_gpcp')
        
        if(iret.ne.0) then 
           write(LIS_logunit,*) &
                'Could not retrieve entries in file: ',trim(fname)
           ferror_gpcp = 0
           deallocate(lb)
           deallocate(gpcpin)
           return           
        endif

        ! Trap the old "Could not find correct forcing parameter in file"
        ! error from LIS 6.  This error occurred right before a GPCP
        ! grid change.  LIS would try to read ahead, but the new data
        ! would be on the new grid, so LIS would misread it resulting in
        ! the above error message.  The LIS would not use GPCP forcing
        ! and defaluts back to the base forcing field.
        ! Trap this by checking the number of values in one of the
        ! GRIB fields.
        call grib_get_size(igrib,'values',grid_size)
        if ( grid_size /= ngpcp ) then
           write(LIS_logunit,*) &
           'ERR: Number of values does not match due to gridchange', trim(fname)
           ferror_gpcp = 0
           deallocate(lb)
           deallocate(gpcpin)
           return
        endif

        call grib_get(igrib,'indicatorOfParameter',pds5_val,rc)
        call LIS_verify(rc, 'error in grib_get: indicatorOfParameter in read_gpcp')
        
        call grib_get(igrib,'level',pds7_val,rc)
        call LIS_verify(rc, 'error in grib_get: level in read_gpcp')
        
        var_found = .false. 
        
        do iv=1,iv_total
           if((pds5_val.eq.pds5(iv)).and.&
                (pds7_val.eq.pds7(iv))) then 
              var_found = .true.
              var_index = iv
              var_status(iv) = .true. 
              exit
           endif
        enddo
        gpcpin = -9999.0
        call grib_get(igrib,'values',gpcpin,rc)
        call LIS_warning(rc, 'error in grib_get:values in read_gpcp')
        
        if(rc.ne.0) then 
           write(LIS_logunit,*) &
                'Could not retrieve entries in file: ',trim(fname)
           ferror_gpcp = 0
           deallocate(lb)
           deallocate(gpcpin)
           return           
        endif

        call grib_get(igrib,'missingValue',missingValue,rc)
        call LIS_verify(rc, 'error in grib_get:missingValue in read_gpcp')

        call grib_release(igrib,rc)
        call LIS_verify(rc, 'error in grib_release in read_gpcp')
        
        if(var_found) then 
           lb = .false.
           do t=1,ngpcp
              if(gpcpin(t).ne.missingValue) lb(t) = .true. 
           enddo

         ! Spatially interpolate the GPCP field to target LIS grid:
           call interp_gpcp(n,findex,ngpcp,gpcpin,lb,LIS_rc%gridDesc, &
                            LIS_rc%lnc(n),LIS_rc%lnr(n),precip_regrid)
           !write(LIS_logunit,*) 'writing fort.888'
           !write(888) precip_regrid

           do j = 1,LIS_rc%lnr(n)
              do i = 1,LIS_rc%lnc(n)
                    index1 = LIS_domain(n)%gindex(i,j)
                    if(LIS_domain(n)%gindex(i,j) .ne. -1) then
                       if(order.eq.2) then 
                          gpcp_struc(n)%metdata2(1,index1) = &
                               precip_regrid(i,j) !*3600.0
                       endif
                    endif
              enddo
           enddo
        endif
     end do
  
     call grib_close_file(ftn)

     deallocate(lb)
     deallocate(gpcpin)
     
     do kk=1,iv_total
        if(.not.var_status(kk)) then 
           write(LIS_logunit,*) &
                'Could not retrieve entries in file: ',trim(fname)
           ferror_gpcp = 0
           if(order.eq.2) then
              gpcp_struc(n)%metdata2(1,:) = LIS_rc%udef
           endif
           
           return
        endif
     enddo
  else
     write(LIS_logunit,*) &
          'Could not find file: ',trim(fname)
     ferror_gpcp = 0
     if(order.eq.2) then
          gpcp_struc(n)%metdata2(1,:) = LIS_rc%udef
     endif
  endif
#endif     
     
end subroutine read_gpcp





