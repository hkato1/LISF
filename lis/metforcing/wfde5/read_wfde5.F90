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
! !ROUTINE: read_wfde5
! \label{read_wfde5}
! 
! !REVISION HISTORY:
! 28 mar 2022: Hiroko Beaudoing, initial code 
!
! !INTERFACE:
subroutine read_wfde5(n, kk,order, year, month, day, hour, findex,          &
     fname, ferror)
! !USES:
  use LIS_coreMod,       only : LIS_rc, LIS_domain, LIS_masterproc
  use LIS_logMod
  use LIS_FORC_AttributesMod
  use LIS_metforcingMod, only : LIS_forc
  use wfde5_forcingMod,  only : wfde5_struc
  use LIS_constantsMod,  only : LIS_CONST_PATH_LEN
#if (defined USE_NETCDF3 || defined USE_NETCDF4) 
  use netcdf
#endif

  implicit none
! !ARGUMENTS: 
  integer, intent(in)          :: n
  integer, intent(in)          :: kk
  integer, intent(in)          :: order
  integer, intent(in)          :: year
  integer, intent(in)          :: month
  integer, intent(in)          :: day
  integer, intent(in)          :: hour
  integer, intent(in)          :: findex
  character(len=*), dimension(8), intent(in) :: fname
  integer, intent(out)         :: ferror          

!
! !DESCRIPTION:
!  For the given time, reads parameters from
!  WFDE5 data, transforms into 9 LIS forcing 
!  parameters and interpolates to the LIS domain. \newline
!
!  PSurf: Atmospheric pressure (Pa)        --fname(1)
!  Tair: Atmospheric temperature (K) (2m)  --fname(2)
!  Qair: Atmospheric humidity (kg/kg)      --fname(3)
!  Wind: Wind speed (m/s) (10m)            --fname(4)
!  Rainf: Rain (kg/m2/s)                   --fname(5)
!  Snowf: Snow (kg/m2/s)                   --fname(6)
!  LWdown: Long-wave radiation (W/m2)      --fname(7)
!  SWdown: Short-wave radiation (W/m2)     --fname(8)
!
!  The arguments are: 
!  \begin{description}
!  \item[order]
!    flag indicating which data to be read (order=1, read the previous 
!    1 hourly instance, order=2, read the next 1 hourly instance)
!  \item[n]
!    index of the nest
!  \item[fname]
!    names of the WFDE5 reanalysis file for each variable
!  \item[tscount]
!    time step count
!  \item[ferror]
!    return error code (0 indicates success)
!  \end{description}
! 
!  The routines invoked are: 
!  \begin{description}
!  \item[bilinear\_interp](\ref{bilinear_interp}) \newline
!    spatially interpolate the forcing data using bilinear interpolation
!  \item[conserv\_interp](\ref{conserv_interp}) \newline
!    spatially interpolate the forcing data using conservative interpolation
!  \end{description}
!EOP
  
  integer   :: ftn
  integer   :: tmpId, qId, windId, lwdId, psId, rainfId, snowfID, swdId
  integer   :: tindex,tindex2
  logical   :: file_exists
  integer   :: c,r,t,k,l,iret,status
  integer   :: cc, rr
  integer   :: mo,rec_size
  logical   :: read_lnd
  logical   :: read_flag
  real      :: undef  ! input data FillValue
  character(len=LIS_CONST_PATH_LEN) :: infile
  
  real, allocatable      :: ps(:,:,:)
  real, allocatable      :: tair(:,:,:)
  real, allocatable      :: qair(:,:,:)
  real, allocatable      :: wind(:,:,:)
  real, allocatable      :: rainf(:,:,:)
  real, allocatable      :: snowf(:,:,:)
  real, allocatable      :: swd(:,:,:)
  real, allocatable      :: lwd(:,:,:)
  real, allocatable      :: temp(:,:,:)

  integer :: days(12)
  data days /31,28,31,30,31,30,31,31,30,31,30,31/
! __________________________________________________________________________

  read_flag = .false.
  ! Check if it is the switch of a month
  if(order.eq.1) then
     if(wfde5_struc(n)%mo1.ne.month) then
        wfde5_struc(n)%mo1 = month
        read_flag = .true.
     endif
  else
     if(wfde5_struc(n)%mo2.ne.month) then
        wfde5_struc(n)%mo2 = month
        read_flag = .true.
     endif
  endif
  ferror = 1

  if(read_flag) then 
#if (defined USE_NETCDF4) 

     if((mod(year,4) .eq. 0 .and. mod(year, 100).ne.0) &!leap year
          .or.(mod(year,400) .eq.0)) then 
        days(2) = 29
     else 
        days(2) = 28
     endif
     
     mo = LIS_rc%lnc(n)*LIS_rc%lnr(n)

! Read single layer file fields:
     inquire(file=fname(8),exist=file_exists) 
     if(file_exists) then 
        if(order.eq.1) then
           wfde5_struc(n)%ps1     = LIS_rc%udef
           wfde5_struc(n)%tair1   = LIS_rc%udef
           wfde5_struc(n)%qair1   = LIS_rc%udef
           wfde5_struc(n)%wind1   = LIS_rc%udef
           wfde5_struc(n)%rainf1  = LIS_rc%udef
           wfde5_struc(n)%swd1    = LIS_rc%udef
           wfde5_struc(n)%lwd1    = LIS_rc%udef
        else
           wfde5_struc(n)%ps2     = LIS_rc%udef
           wfde5_struc(n)%tair2   = LIS_rc%udef
           wfde5_struc(n)%qair2   = LIS_rc%udef
           wfde5_struc(n)%wind2   = LIS_rc%udef
           wfde5_struc(n)%rainf2  = LIS_rc%udef
           wfde5_struc(n)%swd2    = LIS_rc%udef
           wfde5_struc(n)%lwd2    = LIS_rc%udef
        endif

        rec_size = days(month)*24 

        infile = fname(1)
        write(LIS_logunit,*)'[INFO] Reading WFDE5 file (bookend,', order,' -',trim(infile), ')', rec_size

        call LIS_verify(nf90_open(path=trim(infile), mode=NF90_NOWRITE, &
             ncid=ftn), 'nf90_open failed in read_wfde5')

        allocate(ps(wfde5_struc(n)%ncold,wfde5_struc(n)%nrold,rec_size))
        
        call LIS_verify(nf90_inq_varid(ftn,'PSurf',psId), &
             'nf90_inq_varid failed for psurf in read_wfde5')
        call LIS_verify(nf90_get_var(ftn,psId, ps),&
             'nf90_get_var failed for ps in read_wfde5') 
        call LIS_verify(nf90_get_att(ftn, psId,'_FillValue', undef), &
             'nf90_get_att failed for ps:_FillValue in read_wfde5')

        if(order.eq.1) then 
           do l=1,rec_size
              call interp_wfde5_var(n,findex,month,ps(:,:,l),   .false.,&
                   undef, wfde5_struc(n)%ps1(:,l))              
           enddo
        else
           do l=1,rec_size
              call interp_wfde5_var(n,findex,month,ps(:,:,l),    .false.,&
                   undef, wfde5_struc(n)%ps2(:,l))
           enddo
        endif
        deallocate(ps)
        call LIS_verify(nf90_close(ftn), &
             'failed to close file in read_wfde5')

        infile = fname(2)
        write(LIS_logunit,*)'[INFO] Reading WFDE5 file (bookend,', order,' -',trim(infile), ')'

        call LIS_verify(nf90_open(path=trim(infile), mode=NF90_NOWRITE, &
             ncid=ftn), 'nf90_open failed in read_wfde5')
        allocate(tair(wfde5_struc(n)%ncold,wfde5_struc(n)%nrold,rec_size))
        call LIS_verify(nf90_inq_varid(ftn,'Tair',tmpId), &
             'nf90_inq_varid failed for Tair in read_wfde5')
        call LIS_verify(nf90_get_var(ftn,tmpId, tair),&
             'nf90_get_var failed for tair in read_wfde5') 
        call LIS_verify(nf90_get_att(ftn, tmpId,'_FillValue', undef), &
             'nf90_get_att failed for tair:_FillValue in read_wfde5')

        if(order.eq.1) then 
           do l=1,rec_size
              call interp_wfde5_var(n,findex,month,tair(:,:,l), .false.,&
                   undef, wfde5_struc(n)%tair1(:,l))
           enddo
        else
           do l=1,rec_size
              call interp_wfde5_var(n,findex,month,tair(:,:,l),  .false.,&
                   undef, wfde5_struc(n)%tair2(:,l))
           enddo
        endif
        deallocate(tair)
        call LIS_verify(nf90_close(ftn), &
             'failed to close file in read_wfde5')

        infile = fname(3)
        write(LIS_logunit,*)'[INFO] Reading WFDE5 file (bookend,', order,' -',trim(infile), ')'

        call LIS_verify(nf90_open(path=trim(infile), mode=NF90_NOWRITE, &
             ncid=ftn), 'nf90_open failed in read_wfde5')
        allocate(qair(wfde5_struc(n)%ncold,wfde5_struc(n)%nrold,rec_size))
        call LIS_verify(nf90_inq_varid(ftn,'Qair',qId), &
             'nf90_inq_varid failed for Qair in read_wfde5')
        call LIS_verify(nf90_get_var(ftn,qId, qair),&
             'nf90_get_var failed for qair in read_wfde5') 
        call LIS_verify(nf90_get_att(ftn, qId,'_FillValue', undef), &
             'nf90_get_att failed for qair:_FillValue in read_wfde5')

        if(order.eq.1) then 
           do l=1,rec_size
              call interp_wfde5_var(n,findex,month,qair(:,:,l), .false.,&
                   undef, wfde5_struc(n)%qair1(:,l))
           enddo
        else
           do l=1,rec_size
              call interp_wfde5_var(n,findex,month,qair(:,:,l),  .false.,&
                   undef, wfde5_struc(n)%qair2(:,l))
           enddo
        endif
        deallocate(qair)
        call LIS_verify(nf90_close(ftn), &
             'failed to close file in read_wfde5')

        infile = fname(4)
        write(LIS_logunit,*)'[INFO] Reading WFDE5 file (bookend,', order,' -',trim(infile), ')'

        call LIS_verify(nf90_open(path=trim(infile), mode=NF90_NOWRITE, &
             ncid=ftn), 'nf90_open failed in read_wfde5')
        allocate(wind(wfde5_struc(n)%ncold,wfde5_struc(n)%nrold,rec_size))
        call LIS_verify(nf90_inq_varid(ftn,'Wind',windId), &
             'nf90_inq_varid failed for Wind in read_wfde5')
        call LIS_verify(nf90_get_var(ftn,windId, wind),&
             'nf90_get_var failed for wind in read_wfde5') 
        call LIS_verify(nf90_get_att(ftn, windId,'_FillValue', undef), &
             'nf90_get_att failed for wind:_FillValue in read_wfde5')
        if(order.eq.1) then 
           do l=1,rec_size
              call interp_wfde5_var(n,findex,month,wind(:,:,l), .false.,&
                   undef, wfde5_struc(n)%wind1(:,l))
           enddo
        else
           do l=1,rec_size
              call interp_wfde5_var(n,findex,month,wind(:,:,l),  .false.,&
                   undef, wfde5_struc(n)%wind2(:,l))
           enddo
        endif
        deallocate(wind)
        call LIS_verify(nf90_close(ftn), &
             'failed to close file in read_wfde5')
        
        ! special rec_size in 1979/01 for Rainf/Snowf/SWdown/LWdown
        if ( year.eq.1979 .and. month.eq.1 ) then
         rec_size = 737
         write(LIS_logunit,*) 'special record size for 1979-01'
        endif

        infile = fname(5)
        write(LIS_logunit,*)'[INFO] Reading WFDE5 file (bookend,', order,' -',trim(infile), ')'

        call LIS_verify(nf90_open(path=trim(infile), mode=NF90_NOWRITE, &
             ncid=ftn), 'nf90_open failed in read_wfde5')
        allocate(rainf(wfde5_struc(n)%ncold,wfde5_struc(n)%nrold,rec_size))
        allocate(temp(wfde5_struc(n)%ncold,wfde5_struc(n)%nrold,rec_size))
        call LIS_verify(nf90_inq_varid(ftn,'Rainf',rainfId), &
             'nf90_inq_varid failed for Rainf in read_wfde5')
        call LIS_verify(nf90_get_var(ftn,rainfId, temp),&
             'nf90_get_var failed for rainf in read_wfde5') 
        call LIS_verify(nf90_get_att(ftn, rainfId,'_FillValue', undef), &
             'nf90_get_att failed for rainf:_FillValue in read_wfde5')
        call LIS_verify(nf90_close(ftn), &
             'failed to close file in read_wfde5')

        infile = fname(6)
        write(LIS_logunit,*)'[INFO] Reading WFDE5 file (bookend,', order,' -',trim(infile), ')'

        call LIS_verify(nf90_open(path=trim(infile), mode=NF90_NOWRITE, &
             ncid=ftn), 'nf90_open failed in read_wfde5')
        allocate(snowf(wfde5_struc(n)%ncold,wfde5_struc(n)%nrold,rec_size))
        call LIS_verify(nf90_inq_varid(ftn,'Snowf',snowfId), &
             'nf90_inq_varid failed for Snowf in read_wfde5')
        call LIS_verify(nf90_get_var(ftn,snowfId, snowf),&
             'nf90_get_var failed for snowf in read_wfde5') 
        call LIS_verify(nf90_get_att(ftn, snowfId,'_FillValue', undef), &
             'nf90_get_att failed for snowf:_FillValue in read_wfde5')
        call LIS_verify(nf90_close(ftn), &
             'failed to close file in read_wfde5')

        rainf = 0.0
        do r=1,wfde5_struc(n)%nrold
         do c=1,wfde5_struc(n)%ncold
           do l=1,rec_size
              if(temp(c,r,l).ne.undef.and.&
                   snowf(c,r,l).ne.undef) then 
                 rainf(c,r,l) = temp(c,r,l) + snowf(c,r,l)
              else
                 rainf(c,r,l) = undef
              endif
           enddo
         enddo
        enddo
        deallocate(snowf)
        deallocate(temp)
        if(order.eq.1) then 
           do l=1,rec_size
              call interp_wfde5_var(n,findex,month,rainf(:,:,l),.true.,&
                   undef, wfde5_struc(n)%rainf1(:,l))
           enddo
        else
           do l=1,rec_size
              call interp_wfde5_var(n,findex,month,rainf(:,:,l), .true.,&
                   undef, wfde5_struc(n)%rainf2(:,l))
           enddo
        endif
        deallocate(rainf)

        infile = fname(7)
        write(LIS_logunit,*)'[INFO] Reading WFDE5 file (bookend,', order,' -',trim(infile), ')'

        call LIS_verify(nf90_open(path=trim(infile), mode=NF90_NOWRITE, &
             ncid=ftn), 'nf90_open failed in read_wfde5')
        allocate(swd(wfde5_struc(n)%ncold,wfde5_struc(n)%nrold,rec_size))
        call LIS_verify(nf90_inq_varid(ftn,'SWdown',swdId), &
             'nf90_inq_varid failed for SWdown in read_wfde5')
        call LIS_verify(nf90_get_var(ftn,swdId, swd),&
             'nf90_get_var failed for SWdown in read_wfde5') 
        call LIS_verify(nf90_get_att(ftn, swdId,'_FillValue', undef), &
             'nf90_get_att failed for swd:_FillValue in read_wfde5')

        if(order.eq.1) then 
           do l=1,rec_size
              call interp_wfde5_var(n,findex,month,swd(:,:,l),  .false.,&
                   undef, wfde5_struc(n)%swd1(:,l))
           enddo
        else
           do l=1,rec_size
              call interp_wfde5_var(n,findex,month,swd(:,:,l),   .false.,&
                   undef, wfde5_struc(n)%swd2(:,l))
           enddo
        endif
        deallocate(swd)
        call LIS_verify(nf90_close(ftn), &
             'failed to close file in read_wfde5')

        infile = fname(8)
        write(LIS_logunit,*)'[INFO] Reading WFDE5 file (bookend,', order,' -',trim(infile), ')'

        call LIS_verify(nf90_open(path=trim(infile), mode=NF90_NOWRITE, &
             ncid=ftn), 'nf90_open failed in read_wfde5')
        allocate(lwd(wfde5_struc(n)%ncold,wfde5_struc(n)%nrold,rec_size))
        call LIS_verify(nf90_inq_varid(ftn,'LWdown',lwdId), &
             'nf90_inq_varid failed for LWdown in read_wfde5')
        call LIS_verify(nf90_get_var(ftn,lwdId, lwd),&
             'nf90_get_var failed for lwd in read_wfde5')
        call LIS_verify(nf90_get_att(ftn, lwdId,'_FillValue', undef), &
             'nf90_get_att failed for lwd:_FillValue in read_wfde5')

        if(order.eq.1) then 
           do l=1,rec_size
              call interp_wfde5_var(n,findex,month,lwd(:,:,l),  .false.,&
                   undef, wfde5_struc(n)%lwd1(:,l))
           enddo
        else
           do l=1,rec_size
              call interp_wfde5_var(n,findex,month,lwd(:,:,l),   .false.,&
                   undef, wfde5_struc(n)%lwd2(:,l))
           enddo
        endif
        deallocate(lwd)
        call LIS_verify(nf90_close(ftn), &
             'failed to close file in read_wfde5')
        
                
     else
        write(LIS_logunit,*) '[ERR] ',trim(fname(1))//' does not exist'
        call LIS_endrun()
        
     endif
  endif

  tindex = (day - 1)*24 + hour + 1
  tindex2 = tindex - 7
  
  if(order.eq.1) then 
    if ( year.eq.1979 .and. month.eq.1 ) then
     call assign_processed_wfde5forc(n,kk,order,1,wfde5_struc(n)%tair1(:,tindex))
     call assign_processed_wfde5forc(n,kk,order,2,wfde5_struc(n)%qair1(:,tindex))
     call assign_processed_wfde5forc(n,kk,order,3,wfde5_struc(n)%swd1(:,tindex2))
     call assign_processed_wfde5forc(n,kk,order,4,wfde5_struc(n)%lwd1(:,tindex2))
     call assign_processed_wfde5forc(n,kk,order,5,wfde5_struc(n)%wind1(:,tindex))
     call assign_processed_wfde5forc(n,kk,order,6,wfde5_struc(n)%ps1(:,tindex))
     call assign_processed_wfde5forc(n,kk,order,7,wfde5_struc(n)%rainf1(:,tindex2))
    else
     call assign_processed_wfde5forc(n,kk,order,1,wfde5_struc(n)%tair1(:,tindex))
     call assign_processed_wfde5forc(n,kk,order,2,wfde5_struc(n)%qair1(:,tindex))
     call assign_processed_wfde5forc(n,kk,order,3,wfde5_struc(n)%swd1(:,tindex))
     call assign_processed_wfde5forc(n,kk,order,4,wfde5_struc(n)%lwd1(:,tindex))
     call assign_processed_wfde5forc(n,kk,order,5,wfde5_struc(n)%wind1(:,tindex))
     call assign_processed_wfde5forc(n,kk,order,6,wfde5_struc(n)%ps1(:,tindex))
     call assign_processed_wfde5forc(n,kk,order,7,wfde5_struc(n)%rainf1(:,tindex))
    endif
  else
    if ( year.eq.1979 .and. month.eq.1 ) then
     call assign_processed_wfde5forc(n,kk,order,1,wfde5_struc(n)%tair2(:,tindex))
     call assign_processed_wfde5forc(n,kk,order,2,wfde5_struc(n)%qair2(:,tindex))
     call assign_processed_wfde5forc(n,kk,order,3,wfde5_struc(n)%swd2(:,tindex2))
     call assign_processed_wfde5forc(n,kk,order,4,wfde5_struc(n)%lwd2(:,tindex2))
     call assign_processed_wfde5forc(n,kk,order,5,wfde5_struc(n)%wind2(:,tindex))
     call assign_processed_wfde5forc(n,kk,order,6,wfde5_struc(n)%ps2(:,tindex))
     call assign_processed_wfde5forc(n,kk,order,7,wfde5_struc(n)%rainf2(:,tindex2))
    else
     call assign_processed_wfde5forc(n,kk,order,1,wfde5_struc(n)%tair2(:,tindex))
     call assign_processed_wfde5forc(n,kk,order,2,wfde5_struc(n)%qair2(:,tindex))
     call assign_processed_wfde5forc(n,kk,order,3,wfde5_struc(n)%swd2(:,tindex))
     call assign_processed_wfde5forc(n,kk,order,4,wfde5_struc(n)%lwd2(:,tindex))
     call assign_processed_wfde5forc(n,kk,order,5,wfde5_struc(n)%wind2(:,tindex))
     call assign_processed_wfde5forc(n,kk,order,6,wfde5_struc(n)%ps2(:,tindex))
     call assign_processed_wfde5forc(n,kk,order,7,wfde5_struc(n)%rainf2(:,tindex))
    endif
  endif

#endif

end subroutine read_wfde5


!BOP
! 
! !ROUTINE: interp_wfde5_var
! \label{interp_wfde5_var}
! 
! !INTERFACE: 
subroutine interp_wfde5_var(n,findex, month, input_var, &
     pcp_flag, undef, output_var)

! !USES: 
  use LIS_coreMod
  use LIS_logMod
  use LIS_spatialDownscalingMod
  use wfde5_forcingMod, only : wfde5_struc
#if(defined USE_NETCDF3 || defined USE_NETCDF4)      
  use netcdf
#endif
  implicit none

! !ARGUMENTS:
  integer, intent(in)    :: n
  integer, intent(in)    :: findex
  integer, intent(in)    :: month
  real,    intent(in)    :: input_var(wfde5_struc(n)%ncold,wfde5_struc(n)%nrold)
  real,    intent(out)   :: output_var(LIS_rc%lnc(n)*LIS_rc%lnr(n))
  logical, intent(in)    :: pcp_flag
  real, intent(in)       :: undef

!
! !DESCRIPTION: 
!  This subroutine spatially interpolates a WFDE5 field
!  to the LIS running domain
! 
!EOP

  integer   :: t,c,r,k,iret
  integer   :: pcp1Id, pcp2Id, pcp3Id, pcp4Id,pcp5Id, pcp6Id
  real      :: f (wfde5_struc(n)%ncold*wfde5_struc(n)%nrold)
  logical*1 :: lb(wfde5_struc(n)%ncold*wfde5_struc(n)%nrold)
  logical*1 :: lo(LIS_rc%lnc(n)*LIS_rc%lnr(n))
  integer   :: input_size
! _____________________________________________________________

  input_size = wfde5_struc(n)%ncold*wfde5_struc(n)%nrold
  output_var = LIS_rc%udef

!-----------------------------------------------------------------------    
! Apply corrections
!-----------------------------------------------------------------------  
     
  lb = .false.
  k = 0
  do r=1,wfde5_struc(n)%nrold
     do c=1,wfde5_struc(n)%ncold           
        k= k + 1
        if(input_var(c,r).ne.undef) then 
           f(k) = input_var(c,r)
           lb(k) = .true.
        else
           f(k) = LIS_rc%udef
           lb(k) = .false.
        endif
     enddo
  enddo
  
!-----------------------------------------------------------------------    
! Apply downscaling
!-----------------------------------------------------------------------    
  if(pcp_flag.and.LIS_rc%pcp_downscale(findex).ne.0) then 
     call LIS_generatePcpClimoRatioField(n,findex,"WFDE5",&
          month, & 
          input_size, &
          f, &
          lb)     
  endif
          
!-----------------------------------------------------------------------
! Interpolate to LIS grid
!-----------------------------------------------------------------------
  if(pcp_flag.and.&
       trim(wfde5_struc(n)%met_interp).eq."budget-bilinear") then 
     call conserv_interp(LIS_rc%gridDesc(n,:),lb,f,lo,&
          output_var(:), &
          wfde5_struc(n)%mi,LIS_rc%lnc(n)*LIS_rc%lnr(n),& 
          LIS_domain(n)%lat, LIS_domain(n)%lon,&
          wfde5_struc(n)%w112,wfde5_struc(n)%w122,&
          wfde5_struc(n)%w212,wfde5_struc(n)%w222,&
          wfde5_struc(n)%n112,wfde5_struc(n)%n122,&
          wfde5_struc(n)%n212,wfde5_struc(n)%n222,&
          LIS_rc%udef, iret)
     
  elseif(trim(wfde5_struc(n)%met_interp).eq."bilinear".or.&
       trim(wfde5_struc(n)%met_interp).eq."budget-bilinear") then 

     call bilinear_interp(LIS_rc%gridDesc(n,:),lb,f,lo,&
          output_var(:), &
          wfde5_struc(n)%mi,LIS_rc%lnc(n)*LIS_rc%lnr(n), & 
          LIS_domain(n)%lat, LIS_domain(n)%lon,&
          wfde5_struc(n)%w111,wfde5_struc(n)%w121,&
          wfde5_struc(n)%w211,wfde5_struc(n)%w221,&
          wfde5_struc(n)%n111,wfde5_struc(n)%n121,&
          wfde5_struc(n)%n211,wfde5_struc(n)%n221,&
          LIS_rc%udef, iret)
     
  elseif(trim(wfde5_struc(n)%met_interp).eq."neighbor") then 
     call neighbor_interp(LIS_rc%gridDesc(n,:),lb,f,lo,&
          output_var(:),wfde5_struc(n)%mi,&
          LIS_rc%lnc(n)*LIS_rc%lnr(n),&
          LIS_domain(n)%lat, LIS_domain(n)%lon,&
          wfde5_struc(n)%n113,LIS_rc%udef,iret)

  elseif(trim(wfde5_struc(n)%met_interp).eq."average") then 
     call upscaleByAveraging(wfde5_struc(n)%mi, &
          LIS_rc%lnc(n)*LIS_rc%lnr(n), LIS_rc%udef, &
          wfde5_struc(n)%n111, lb, f, lo, output_var(:))

  else
     write(LIS_logunit,*) '[ERR] Spatial interpolation option '//&
          trim(wfde5_struc(n)%met_interp)//' not supported for WFDE5'
     call LIS_endrun()
  endif
  
  if( pcp_flag.and.LIS_rc%pcp_downscale(findex).ne.0 ) then 
     
     call LIS_pcpClimoDownscaling(n, findex, month,&
          LIS_rc%lnc(n)*LIS_rc%lnr(n), output_var(:), lo)
     
  endif

end subroutine interp_wfde5_var

!BOP
! 
! !ROUTINE: assign_processed_wfde5forc
! \label{assign_processed_wfde5forc}
! 
! !INTERFACE: 
subroutine assign_processed_wfde5forc(n,kk,order,var_index,wfde5forc)
! !USES: 
  use LIS_coreMod
  use wfde5_forcingMod, only : wfde5_struc
!
! !DESCRIPTION: 
!  This routine assigns the interpolated WFDE5 forcing data
!  to the module data structures to be used later for 
!  time interpolation 
!
!EOP
  implicit none

  integer :: n
  integer :: kk
  integer :: order
  integer :: var_index
  real    :: wfde5forc(LIS_rc%lnc(n)*LIS_rc%lnr(n))
  

  integer :: c,r

  do r=1,LIS_rc%lnr(n)
     do c=1,LIS_rc%lnc(n)
        if(LIS_domain(n)%gindex(c,r).ne.-1) then 
           if(order.eq.1) then 
              wfde5_struc(n)%metdata1(kk,var_index,&
                   LIS_domain(n)%gindex(c,r)) = &
                   wfde5forc(c+(r-1)*LIS_rc%lnc(n))
           elseif(order.eq.2) then 
              wfde5_struc(n)%metdata2(kk,var_index,&
                   LIS_domain(n)%gindex(c,r)) = &
                   wfde5forc(c+(r-1)*LIS_rc%lnc(n))
           endif
        endif
     enddo
  enddo
end subroutine assign_processed_wfde5forc



