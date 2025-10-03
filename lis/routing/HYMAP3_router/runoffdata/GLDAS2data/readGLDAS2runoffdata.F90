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
!BOP
! 
! !DESCRIPTION:
! This routine is for reading in GLDAS 2.0 data in native NetCDF format.
! 
! !REVISION HISTORY: 
! 30 Jan 2016: Hiroko Beaudoing, Initial implementation
! 31 Aug 2016: Augusto Getirana, Fix file name format
! 25 Sep 2025: Hiroko Beaudoing, updated for HYMAP3 
! 
!caveats
! 1) assumes the LIS outputs are the same output interval as that of
! the HYMAP model timestep. No temporal aggregation is done.
! 2) Assumes that the units of Qs and Qsb are accumulation in kg/m2/3hr 
! 3) Assumes that the LIS outputs are in the same grid/map projection/
! resolution.
! 4) LIS outputs are in NETCDF format.
!
! !USES: 
!
!EOP 
subroutine readGLDAS2runoffdata(n,surface_runoff, baseflow)

  use LIS_coreMod
  use LIS_timeMgrMod
  use LIS_logMod
  use GLDAS2runoffdataMod
  use LIS_fileIOMod
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none

  integer,          intent(in) :: n
  real                         :: surface_runoff(LIS_rc%lnc(n),LIS_rc%lnr(n))
  real                         :: baseflow(LIS_rc%lnc(n),LIS_rc%lnr(n))
  integer                       :: nc,nr
  integer                       :: c,r
  !ag - 17Mar2016
  real                  :: qs2d(GLDAS2runoffdata_struc(n)%nc,GLDAS2runoffdata_struc(n)%nr)
  real                  :: qsb2d(GLDAS2runoffdata_struc(n)%nc,GLDAS2runoffdata_struc(n)%nr)
  logical*1             :: lb(GLDAS2runoffdata_struc(n)%nc*GLDAS2runoffdata_struc(n)%nr)
  real                  :: var_input(GLDAS2runoffdata_struc(n)%nc*GLDAS2runoffdata_struc(n)%nr)
  logical*1             :: lo(LIS_rc%lnc(n)*LIS_rc%lnr(n))
  real                  :: var_out(LIS_rc%lnc(n)*LIS_rc%lnr(n))
  integer                       :: ios, nid,qsid,qsbid
  character*100                 :: filename
  logical                       :: file_exists
  logical                       :: check_Flag
  real                          :: undef
  integer                       :: doy, yr, mo, da, hr, mn, ss, ts
  real*8                        :: time
  real                          :: gmt

  external :: neighbor_interp
  external :: upscaleByAveraging
  
  yr =LIS_rc%yr    !Next Hour
  mo =LIS_rc%mo
  da =LIS_rc%da
  hr=LIS_rc%hr-imod(LIS_rc%hr,int(GLDAS2runoffdata_struc(n)%outInterval/3600.))
  !print*,int(GLDAS2runoffdata_struc(n)%outInterval/3600.),LIS_rc%hr,hr
  mn =0
  ss =0
  ts =GLDAS2runoffdata_struc(n)%outInterval

  call LIS_tick(time,doy,gmt,yr,mo,da,hr,mn,ss,real(ts))

  call create_GLDAS2_filename(GLDAS2runoffdata_struc(n)%odir,&
       GLDAS2runoffdata_struc(n)%model_name,&
       GLDAS2runoffdata_struc(n)%datares,&
       yr, mo, da, doy, hr, filename)

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  if(trim(GLDAS2runoffdata_struc(n)%previous_filename)/=trim(filename))then
  
    GLDAS2runoffdata_struc(n)%previous_filename=filename

    inquire(file=filename, exist=file_exists)
    if(file_exists) then 
      write(LIS_logunit,*) '[INFO] Reading '//trim(filename)
        
      ios = nf90_open(path=filename,&
           mode=NF90_NOWRITE,ncid=nid)
      call LIS_verify(ios,'Error in readGLDAS2runoffdata')

      ios = nf90_inq_varid(nid,'Qs_acc',qsid)
      call LIS_verify(ios,'failed to inquire Qs_acc field in readGLDAS2runoffdata')

      ios = nf90_inq_varid(nid,'Qsb_acc',qsbid)
      call LIS_verify(ios,'failed to inquire Qsb_acc field in readGLDAS2runoffdata')

      ios = nf90_get_att(nid,qsid,'missing_value',undef)
      call LIS_verify(ios,'failed to read missing_value in readGLDAS2runoffdata')

      if(GLDAS2runoffdata_struc(n)%domainCheck) then
         if(LIS_rc%wopt.eq."2d gridspace") then
            ios = nf90_get_var(nid,qsid,GLDAS2runoffdata_struc(n)%qs, &
                 start=(/LIS_ews_halo_ind(n,LIS_localPet+1),LIS_nss_halo_ind(n,LIS_localPet+1)/),&
                 count=(/LIS_rc%lnc(n),LIS_rc%lnr(n)/))
            call LIS_verify(ios, 'failed to read Qs_acc field in readGLDAS2runoffdata')

            ios = nf90_get_var(nid,qsbid,GLDAS2runoffdata_struc(n)%qsb,&
                 start=(/LIS_ews_halo_ind(n,LIS_localPet+1),LIS_nss_halo_ind(n,LIS_localPet+1)/),&
                 count=(/LIS_rc%lnc(n),LIS_rc%lnr(n)/))
            call LIS_verify(ios, 'failed to read Qsb_acc field in readGLDAS2runoffdata')
         else
            write(LIS_logunit,*) "Stand-alone HYMAP3 is only supported for '2d gridspace' outputs currently"
            call LIS_endrun()
         endif

         call LIS_verify(nf90_close(nid))

      else
         if(LIS_rc%wopt.eq."2d gridspace") then

            ios = nf90_get_var(nid,qsid,qs2d)
            call LIS_verify(ios, 'failed to read Qs_acc field in readGLDAS2runoffdata')

            ios = nf90_get_var(nid,qsbid,qsb2d)
            call LIS_verify(ios, 'failed to read Qsb_acc field in readGLDAS2runoffdata')


            if(LIS_isAtAfinerResolution(n,GLDAS2runoffdata_struc(n)%datares)) then
               lb = .true.
               do r=1,GLDAS2runoffdata_struc(n)%nr
                  do c=1,GLDAS2runoffdata_struc(n)%nc
                     var_input(c+(r-1)*GLDAS2runoffdata_struc(n)%nc) = qs2d(c,r)
                     if(qs2d(c,r).lt.0) then
                        lb(c+(r-1)*GLDAS2runoffdata_struc(n)%nc) = .false.
                     endif
                  enddo
               enddo

               call neighbor_interp(LIS_rc%gridDesc,lb,var_input,  &
                    lo,var_out, &
                    GLDAS2runoffdata_struc(n)%nc*GLDAS2runoffdata_struc(n)%nr,&
                    LIS_rc%lnc(n)*LIS_rc%lnr(n),             &
                    LIS_domain(n)%lat, LIS_domain(n)%lon,  &
                    GLDAS2runoffdata_struc(n)%n11,     &
                    LIS_rc%udef,ios)

               do r=1,LIS_rc%lnr(n)
                  do c=1,LIS_rc%lnc(n)
                     GLDAS2runoffdata_struc(n)%qs(c,r) = var_out(c+(r-1)*LIS_rc%lnc(n))
                  enddo
               enddo

               lb = .true.
               do r=1,GLDAS2runoffdata_struc(n)%nr
                  do c=1,GLDAS2runoffdata_struc(n)%nc
                     var_input(c+(r-1)*GLDAS2runoffdata_struc(n)%nc) = qsb2d(c,r)
                     if(qsb2d(c,r).lt.0) then
                        lb(c+(r-1)*GLDAS2runoffdata_struc(n)%nc) = .false.
                     endif
                  enddo
               enddo

               call neighbor_interp(LIS_rc%gridDesc,lb,var_input,  &
                    lo,var_out, &
                    GLDAS2runoffdata_struc(n)%nc*GLDAS2runoffdata_struc(n)%nr,&
                    LIS_rc%lnc(n)*LIS_rc%lnr(n),             &
                    LIS_domain(n)%lat, LIS_domain(n)%lon,  &
                    GLDAS2runoffdata_struc(n)%n11,     &
                    LIS_rc%udef,ios)

               do r=1,LIS_rc%lnr(n)
                  do c=1,LIS_rc%lnc(n)
                     GLDAS2runoffdata_struc(n)%qsb(c,r) = var_out(c+(r-1)*LIS_rc%lnc(n))
                  enddo
               enddo

            else

               lb = .true.
               do r=1,GLDAS2runoffdata_struc(n)%nr
                  do c=1,GLDAS2runoffdata_struc(n)%nc
                     var_input(c+(r-1)*GLDAS2runoffdata_struc(n)%nc) = qs2d(c,r)
                     if(qs2d(c,r).lt.0) then
                        lb(c+(r-1)*GLDAS2runoffdata_struc(n)%nc) = .false.
                     endif
                  enddo
               enddo

               call upscaleByAveraging(&
                    GLDAS2runoffdata_struc(n)%nc*GLDAS2runoffdata_struc(n)%nr, &
                    LIS_rc%lnc(n)*LIS_rc%lnr(n), &
                    LIS_rc%udef, &
                    GLDAS2runoffdata_struc(n)%n11, lb, &
                    var_input, lo, var_out)

                do r=1,LIS_rc%lnr(n)
                   do c=1,LIS_rc%lnc(n)
                      GLDAS2runoffdata_struc(n)%qs(c,r) = var_out(c+(r-1)*LIS_rc%lnc(n))
                   enddo
                enddo

                lb = .true.
                do r=1,GLDAS2runoffdata_struc(n)%nr
                   do c=1,GLDAS2runoffdata_struc(n)%nc
                      var_input(c+(r-1)*GLDAS2runoffdata_struc(n)%nc) = qsb2d(c,r)
                     if(qsb2d(c,r).lt.0) then
                        lb(c+(r-1)*GLDAS2runoffdata_struc(n)%nc) = .false.
                     endif
                  enddo
               enddo

               call upscaleByAveraging(&
                    GLDAS2runoffdata_struc(n)%nc*GLDAS2runoffdata_struc(n)%nr, &
                    LIS_rc%lnc(n)*LIS_rc%lnr(n), &
                    LIS_rc%udef, &
                    GLDAS2runoffdata_struc(n)%n11, lb, &
                    var_input, lo, var_out)

                do r=1,LIS_rc%lnr(n)
                   do c=1,LIS_rc%lnc(n)
                      GLDAS2runoffdata_struc(n)%qsb(c,r) = var_out(c+(r-1)*LIS_rc%lnc(n))
                   enddo
                enddo

             endif

         else
            write(LIS_logunit,*) "Stand-alone HYMAP3 is only supported for '2d gridspace' outputs currently"
            call LIS_endrun()
         endif

         call LIS_verify(nf90_close(nid))
      endif   ! domainCheck

    else
       write(LIS_logunit,*) 'Failed to find '//trim(filename)
       call LIS_endrun()
    endif   ! file_exists
  endif
#endif

! convert units from kg m-2 3hr-1 to kg m-2 sec-1.
where(GLDAS2runoffdata_struc(n)%qs/=LIS_rc%udef)
  surface_runoff = GLDAS2runoffdata_struc(n)%qs/GLDAS2runoffdata_struc(n)%outInterval
  baseflow = GLDAS2runoffdata_struc(n)%qsb/GLDAS2runoffdata_struc(n)%outInterval
else where
  surface_runoff = 0.0
  baseflow = 0.0
end where

end subroutine readGLDAS2runoffdata

!BOP
! 
! !ROUTINE: create_GLDAS2_filename
! \label{create_GLDAS2_filename}
!
! !INTERFACE: 
subroutine create_GLDAS2_filename(odir,model_name, datares,&
     yr,mo,da, doy,hr,filename)

  use LIS_logMod

! 
! !USES:   
  implicit none
!
! !ARGUMENTS: 
  character(len=*)             :: odir
  character(len=*)             :: model_name
  real                         :: datares
  integer                      :: yr
  integer                      :: mo
  integer                      :: da
  integer                      :: doy
  integer                      :: hr
  character(len=*)             :: filename
!
! !DESCRIPTION:
! 
! This routine creates a timestamped filename for the GLDAS2 data
! based on the given date (year, model name, month)
!
!  The arguments are: 
!  \begin{description}
!   \item[odir]            GLDAS2 base directory
!   \item[model\_name]     name of the model used in the GLDAS run
!   \item[yr]              year of data
!   \item[mo]              month of data
!   \item[filename]        Name of the GLDAS2 file
!  \end{description}
! 
!EOP
  
  character*4             :: fyr
  character*3             :: fdoy
  character*2             :: fmo, fhr, fda
  character*100           :: list_name

  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fdoy, fmt='(i3.3)') doy
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da
  write(unit=fhr, fmt='(i2.2)') hr
  
  if(datares .eq. 0.25) then   
     !ag - 31Aug2016
     list_name = 'ls '//trim(odir)//'/'//trim(fyr)//'/'//trim(fyr)//trim(fmo)//&
          '/GLDAS_'//trim(model_name)//'025_3H.A'//&
          trim(fyr)//trim(fmo)//trim(fda)//'.'//trim(fhr)//&
          '*.020.nc4 > GLDAS2_file'
     filename=trim(odir)//'/'//trim(fyr)//'/'//trim(fyr)//trim(fmo)//&
          '/GLDAS_'//trim(model_name)//'025_3H.A'//&
          trim(fyr)//trim(fmo)//trim(fda)//'.'//trim(fhr)//&
          '00.020.nc4'
  elseif(datares.eq. 1.0) then 
     list_name = 'ls '//trim(odir)//'/'//trim(fyr)//'/'//trim(fdoy)//&
          '/GLDAS_'//trim(model_name)//'10_3H.A'//&
          trim(fyr)//trim(fmo)//trim(fda)//'.'//trim(fhr)//&
          '*.020.nc4 > GLDAS2_file'
     filename=trim(odir)//'/'//trim(fyr)//trim(fmo)//&
          '/GLDAS_'//trim(model_name)//'025_3H.A'//&
          trim(fyr)//trim(fmo)//trim(fda)//'.'//trim(fhr)//&
          '00.020.nc4'

  endif
  
!  call system(trim(list_name))
!  ftn = LIS_getNextUnitNumber()
!  open(ftn,file='GLDAS2_file',status='old',iostat=ierr)
!  do while(ierr.eq.0) 
!     read(ftn,'(a)',iostat=ierr) filename
!     if(ierr.ne.0) then 
!        exit
!     endif
!  enddo
!  call LIS_releaseUnitNumber(ftn)

end subroutine create_GLDAS2_filename

!BOP
!
! !ROUTINE: interp_GLDAS2runoffdata
!  \label{interp_GLDAS2runoffdata}
!
! !INTERFACE:
  subroutine interp_GLDAS2runoffdata(n, nc,nr,var_input,var_output)
!
! !USES:
    use LIS_coreMod
    use GLDAS2runoffdataMod
      
    implicit none
!
! !INPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!   This subroutine spatially interpolates the GLDAS2 variable to the
!   model (LIS) output grid and resolution, using a bilinear interpolation
!   approach.
!
!   The arguments are:
!   \begin{description}
!    \item[nc]      number of columns in the input (GLDAS2) grid
!    \item[nr]      number of rows in the input (GLDAS2) grid
!    \item[var_input] input variable to be interpolated
!    \item[lb]        input bitmap (true//false)
!    \item[var_output] resulting interpolated field
!   \end{description}
!
! !FILES USED:
!
! !REVISION HISTORY:
!
!EOP
!BOP
!
! !ARGUMENTS:
    integer            :: n
    integer            :: nc
    integer            :: nr
    real               :: var_input(nc*nr)
    logical*1          :: lb(nc*nr)
    real               :: var_output(LIS_rc%lnc(n), LIS_rc%lnr(n))
    !EOP
    integer            :: ios
    integer            :: c,r
    logical*1          :: lo(LIS_rc%lnc(n)*LIS_rc%lnr(n))
    real               :: go(LIS_rc%lnc(n)*LIS_rc%lnr(n))

    external :: neighbor_interp
    external :: upscaleByAveraging

    var_output = LIS_rc%udef
    lb = .false.
    do r = 1,nr
       do c = 1,nc
          if (var_input(c+(r-1)*nc).ne.LIS_rc%udef) then
             lb(c+(r-1)*nc) = .true.
          endif
       enddo
    enddo

    if(LIS_isAtAfinerResolution(n,GLDAS2runoffdata_struc(n)%datares)) then
       call neighbor_interp(LIS_rc%gridDesc,lb,var_input,  &
            lo,go,nc*nr,LIS_rc%lnc(n)*LIS_rc%lnr(n),             &
            LIS_domain(n)%lat, LIS_domain(n)%lon,  & 
            GLDAS2runoffdata_struc(n)%n11,                         & 
            LIS_rc%udef,ios)
    else
       call upscaleByAveraging(&
            nc*nr, &
            LIS_rc%lnc(n)*LIS_rc%lnr(n), &
            LIS_rc%udef, &
            GLDAS2runoffdata_struc(n)%n11, lb, &
            var_input, lo, go)
    endif
    do r = 1,LIS_rc%lnr(n)
       do c = 1,LIS_rc%lnc(n)
          var_output(c,r) = go(c+(r-1)*LIS_rc%lnc(n))
       enddo
    enddo

  end subroutine interp_GLDAS2runoffdata


