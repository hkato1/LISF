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
! !DESCRIPTION: 
! 
! !REVISION HISTORY: 
!  7 Jan 2016: Sujay Kumar, Initial implementation
! 17 Mar 2016: Augusto Getirana, changes in input file name generation and surface runoff and baseflow variables - this will reduce the number of times input files are read
! 17 Oct 2022: Hiroko Beaudoing, added handling of accumulated runoff inputs
! 
! !USES: 
subroutine readLISrunoffdata(n,surface_runoff, baseflow)
  use LIS_coreMod
  use LIS_logMod
  use LISrunoffdataMod
  use LIS_fileIOMod

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none

  integer,          intent(in) :: n
  real                         :: surface_runoff(LIS_rc%lnc(n),LIS_rc%lnr(n))
  real                         :: baseflow(LIS_rc%lnc(n),LIS_rc%lnr(n))
  real                         :: total_evapotranspiration(LIS_rc%lnc(n),LIS_rc%lnr(n))

  !caveats 
  ! 1) assumes the LIS outputs are the same output interval as that of
  ! the HYMAP model timestep. No temporal aggregation is done. 
  !
  ! 3) Assumes that the units of Qs and Qsb are in kg/m2s
  !
  ! 4) Assumes that the LIS outputs are in the same grid/map projection/
  ! resolution.
  ! 
  ! 5) LIS outputs are in NETCDF format. 
  !
  !Added total evapotranspiration (Evap)
  integer                       :: c,r,t
  integer, allocatable      :: nqs(:,:)

  !ag - 17Mar2016
  !real,   allocatable       :: qs(:,:),qs_t(:)
  real                  :: qs2d(LISrunoffdata_struc(n)%nc,LISrunoffdata_struc(n)%nr)
  real                  :: qsb2d(LISrunoffdata_struc(n)%nc,LISrunoffdata_struc(n)%nr)
  real,   allocatable   :: qs_t(:),qsb_t(:),evap_t(:)
  logical*1             :: lb(LISrunoffdata_struc(n)%nc*LISrunoffdata_struc(n)%nr)
  real                  :: var_input(LISrunoffdata_struc(n)%nc*LISrunoffdata_struc(n)%nr)
  logical*1             :: lo(LIS_rc%lnc(n)*LIS_rc%lnr(n))
  real                  :: var_out(LIS_rc%lnc(n)*LIS_rc%lnr(n))
  integer               :: ios, nid,qsid,qsbid,evapid
  character*100         :: filename
  logical               :: file_exists
  logical               :: check_Flag
  real                  :: dt    ! time step length (sec) of output interval
  !create LIS filename

  call LIS_create_output_filename(n, &
       filename, "netcdf", check_flag, &
       model_name='SURFACEMODEL', &
       odir=LISrunoffdata_struc(n)%odir,&
       writeint=LISrunoffdata_struc(n)%outInterval)
  
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  
  if(trim(LISrunoffdata_struc(n)%previous_filename)/=trim(filename))then
  
    LISrunoffdata_struc(n)%previous_filename=filename
    
    inquire(file=filename, exist=file_exists)
    if(file_exists) then 
      write(LIS_logunit,*) 'Reading '//trim(filename)
      !ag - 17Mar2016
      !allocate(qs(LIS_rc%lnc(n),LIS_rc%lnr(n)))
      !allocate(qsb(LIS_rc%lnc(n),LIS_rc%lnr(n)))
     
      ios = nf90_open(path=filename,&
           mode=NF90_NOWRITE,ncid=nid)
      call LIS_verify(ios,'Error in readLISrunoffdata')
     
!HKB Qs can be average or accumlation
      LISrunoffdata_struc(n)%accum = .false.
      ios = nf90_inq_varid(nid,'Qs_tavg',qsid)
      call LIS_warning(ios,'failed to inquire Qs_tavg field in readLISrunoffdata')
      if ( ios /= 0 ) then
       ios = nf90_inq_varid(nid,'Qs_acc',qsid)
       call LIS_verify(ios,'failed to inquire Qs_acc field in readLISrunoffdata')
       write(LIS_logunit,*) 'Qs in the units of accumulation'
       LISrunoffdata_struc(n)%accum = .true.
      endif
     
      ios = nf90_inq_varid(nid,'Qsb_tavg',qsbid)
      call LIS_warning(ios,'failed to inquire Qsb_tavg field in readLISrunoffdata')
      if ( ios /= 0 ) then
       ios = nf90_inq_varid(nid,'Qsb_acc',qsbid)
       call LIS_verify(ios,'failed to inquire Qsb_acc field in readLISrunoffdata')
       write(LIS_logunit,*) 'Qsb in the units of accumulation'
       LISrunoffdata_struc(n)%accum = .true.
      endif
    
      ios = nf90_inq_varid(nid,'Evap_tavg',evapid)
      call LIS_verify(ios,'failed to inquire Evap_tavg field in readLISrunoffdata')

      if(LISrunoffdata_struc(n)%domainCheck) then 
         if(LIS_rc%wopt.eq."2d gridspace") then 
            ios = nf90_get_var(nid,qsid,LISrunoffdata_struc(n)%qs, &
                 start=(/LIS_ews_halo_ind(n,LIS_localPet+1),LIS_nss_halo_ind(n,LIS_localPet+1)/),&
                 count=(/LIS_rc%lnc(n),LIS_rc%lnr(n)/))
            call LIS_verify(ios, 'failed to read Qs field in readLISrunoffdata')
            
            ios = nf90_get_var(nid,qsbid,LISrunoffdata_struc(n)%qsb,&
                 start=(/LIS_ews_halo_ind(n,LIS_localPet+1),LIS_nss_halo_ind(n,LIS_localPet+1)/),&
                 count=(/LIS_rc%lnc(n),LIS_rc%lnr(n)/))
            call LIS_verify(ios, 'failed to read Qsb field in readLISrunoffdata')
            
            if(LISrunoffdata_struc(n)%accum) then
             ! convert units from kg/m2 to kg/m2/s
             dt = LISrunoffdata_struc(n)%outInterval
              where(LISrunoffdata_struc(n)%qs/=LIS_rc%udef)
               LISrunoffdata_struc(n)%qs = LISrunoffdata_struc(n)%qs/dt
              end where
              where(LISrunoffdata_struc(n)%qsb/=LIS_rc%udef)
               LISrunoffdata_struc(n)%qsb = LISrunoffdata_struc(n)%qsb/dt
              end where
            endif

            ios = nf90_get_var(nid,evapid,LISrunoffdata_struc(n)%evap,&
                 start=(/LIS_ews_halo_ind(n,LIS_localPet+1),LIS_nss_halo_ind(n,LIS_localPet+1)/),&
                 count=(/LIS_rc%lnc(n),LIS_rc%lnr(n)/))
            call LIS_verify(ios, 'failed to read Evap_tavg field in readLISrunoffdata')
         else
            write(LIS_logunit,*) "Stand-alone HYMAP is only supported for '2d gridspace' outputs currently"
            call LIS_endrun()
         endif
         
         call LIS_verify(nf90_close(nid))

      else
         if(LIS_rc%wopt.eq."2d gridspace") then
            
            ios = nf90_get_var(nid,qsid,qs2d)
            call LIS_verify(ios, 'failed to read Qs field in readLISrunoffdata')
            
            ios = nf90_get_var(nid,qsbid,qsb2d)
            call LIS_verify(ios, 'failed to read Qsb field in readLISrunoffdata')


            if(LIS_isAtAfinerResolution(n,LISrunoffdata_struc(n)%datares)) then

               lb = .true. 
               do r=1,LISrunoffdata_struc(n)%nr
                  do c=1,LISrunoffdata_struc(n)%nc
                     var_input(c+(r-1)*LISrunoffdata_struc(n)%nc) = qs2d(c,r)
                     if(qs2d(c,r).lt.0) then
                        lb(c+(r-1)*LISrunoffdata_struc(n)%nc) = .false.
                     endif
                  enddo
               enddo
                     
               call neighbor_interp(LIS_rc%gridDesc,lb,var_input,  &
                    lo,var_out,LISrunoffdata_struc(n)%nc*LISrunoffdata_struc(n)%nr,&
                    LIS_rc%lnc(n)*LIS_rc%lnr(n),             &
                    LIS_domain(n)%lat, LIS_domain(n)%lon,  &
                    LISrunoffdata_struc(n)%n11,                         &
                    LIS_rc%udef,ios)
               
               do r=1,LIS_rc%lnr(n)
                  do c=1,LIS_rc%lnc(n)
                     LISrunoffdata_struc(n)%qs(c,r) = var_out(c+(r-1)*LIS_rc%lnc(n))
                     if(LISrunoffdata_struc(n)%accum) then
                       ! convert units from kg/m2 to kg/m2/s
                       dt = LISrunoffdata_struc(n)%outInterval
                       if (LISrunoffdata_struc(n)%qs(c,r)/=LIS_rc%udef) then
                         LISrunoffdata_struc(n)%qs(c,r) = LISrunoffdata_struc(n)%qs(c,r)/dt
                       endif
                     endif
                  enddo
               enddo

               lb = .true. 
               do r=1,LISrunoffdata_struc(n)%nr
                  do c=1,LISrunoffdata_struc(n)%nc
                     var_input(c+(r-1)*LISrunoffdata_struc(n)%nc) = qsb2d(c,r)
                     if(qsb2d(c,r).lt.0) then
                        lb(c+(r-1)*LISrunoffdata_struc(n)%nc) = .false.
                     endif
                  enddo
               enddo
                     
               call neighbor_interp(LIS_rc%gridDesc,lb,var_input,  &
                    lo,var_out,LISrunoffdata_struc(n)%nc*LISrunoffdata_struc(n)%nr,&
                    LIS_rc%lnc(n)*LIS_rc%lnr(n),             &
                    LIS_domain(n)%lat, LIS_domain(n)%lon,  &
                    LISrunoffdata_struc(n)%n11,                         &
                    LIS_rc%udef,ios)
               
               do r=1,LIS_rc%lnr(n)
                  do c=1,LIS_rc%lnc(n)
                     LISrunoffdata_struc(n)%qsb(c,r) = var_out(c+(r-1)*LIS_rc%lnc(n))
                     if(LISrunoffdata_struc(n)%accum) then
                       ! convert units from kg/m2 to kg/m2/s
                       dt = LISrunoffdata_struc(n)%outInterval
                       if (LISrunoffdata_struc(n)%qsb(c,r)/=LIS_rc%udef) then
                         LISrunoffdata_struc(n)%qsb(c,r) = LISrunoffdata_struc(n)%qsb(c,r)/dt
                       endif
                     endif
                  enddo
               enddo
               
            else

               lb = .true. 
               do r=1,LISrunoffdata_struc(n)%nr
                  do c=1,LISrunoffdata_struc(n)%nc
                     var_input(c+(r-1)*LISrunoffdata_struc(n)%nc) = qs2d(c,r)
                     if(qs2d(c,r).lt.0) then
                        lb(c+(r-1)*LISrunoffdata_struc(n)%nc) = .false.
                     endif
                  enddo
               enddo
               
               call upscaleByAveraging(&
                    LISrunoffdata_struc(n)%nc*LISrunoffdata_struc(n)%nr, &
                    LIS_rc%lnc(n)*LIS_rc%lnr(n), &
                    LIS_rc%udef, &
                    LISrunoffdata_struc(n)%n11, lb, &
                    var_input, lo, var_out)

                do r=1,LIS_rc%lnr(n)
                   do c=1,LIS_rc%lnc(n)
                      LISrunoffdata_struc(n)%qs(c,r) = var_out(c+(r-1)*LIS_rc%lnc(n))
                      if(LISrunoffdata_struc(n)%accum) then
                        ! convert units from kg/m2 to kg/m2/s
                        dt = LISrunoffdata_struc(n)%outInterval
                        if (LISrunoffdata_struc(n)%qs(c,r)/=LIS_rc%udef) then
                          LISrunoffdata_struc(n)%qs(c,r) = LISrunoffdata_struc(n)%qs(c,r)/dt
                        endif
                      endif
                   enddo
                enddo

                lb = .true.
                do r=1,LISrunoffdata_struc(n)%nr
                   do c=1,LISrunoffdata_struc(n)%nc
                      var_input(c+(r-1)*LISrunoffdata_struc(n)%nc) = qsb2d(c,r)
                     if(qsb2d(c,r).lt.0) then
                        lb(c+(r-1)*LISrunoffdata_struc(n)%nc) = .false.
                     endif
                  enddo
               enddo
               
               call upscaleByAveraging(&
                    LISrunoffdata_struc(n)%nc*LISrunoffdata_struc(n)%nr, &
                    LIS_rc%lnc(n)*LIS_rc%lnr(n), &
                    LIS_rc%udef, &
                    LISrunoffdata_struc(n)%n11, lb, &
                    var_input, lo, var_out)

                do r=1,LIS_rc%lnr(n)
                   do c=1,LIS_rc%lnc(n)
                      LISrunoffdata_struc(n)%qsb(c,r) = var_out(c+(r-1)*LIS_rc%lnc(n))
                      if(LISrunoffdata_struc(n)%accum) then
                        ! convert units from kg/m2 to kg/m2/s
                        dt = LISrunoffdata_struc(n)%outInterval
                        if (LISrunoffdata_struc(n)%qsb(c,r)/=LIS_rc%udef) then
                          LISrunoffdata_struc(n)%qsb(c,r) = LISrunoffdata_struc(n)%qsb(c,r)/dt
                        endif
                      endif
                   enddo
                enddo
               
             endif
            
         else
            write(LIS_logunit,*) "Stand-alone HYMAP is only supported for '2d gridspace' outputs currently"
            call LIS_endrun()
         endif
      
         call LIS_verify(nf90_close(nid))
      endif
   else
      write(LIS_logunit,*) 'Failed to find '//trim(filename)
      call LIS_endrun()
   endif
endif

#endif

  !ag - 17Mar2016
  where(LISrunoffdata_struc(n)%qs/=LIS_rc%udef)
    surface_runoff = LISrunoffdata_struc(n)%qs
    baseflow = LISrunoffdata_struc(n)%qsb
    total_evapotranspiration = LISrunoffdata_struc(n)%evap
  else where
    surface_runoff = 0.0
    baseflow = 0.0
    total_evapotranspiration = 0.0
  end where

end subroutine readLISrunoffdata
