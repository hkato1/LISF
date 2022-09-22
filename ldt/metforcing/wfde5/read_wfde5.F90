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
!BOP
!
! !ROUTINE: read_wfde5
! \label{read_wfde5}
! 
! !REVISION HISTORY:
! 29 mar 2022: Hiroko Beaudoing, initial code 
!
! !INTERFACE:
subroutine read_wfde5(n, kk,order, year, month, day, hour, findex,          &
     fname, ferror)
! !USES:
  use LDT_coreMod,       only : LDT_rc, LDT_domain, LDT_masterproc
  use LDT_logMod
  use LDT_FORC_AttributesMod
  use LDT_metforcingMod, only : LDT_forc
  use wfde5_forcingMod, only : wfde5_struc
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
  character(len=*),dimension(8), intent(in) :: fname
  integer, intent(out)         :: ferror          

!
! !DESCRIPTION:
!  For the given time, reads parameters from
!  WFDE5 data, transforms into 9 LDT forcing 
!  parameters and interpolates to the LDT domain. \newline
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
!  \item[name]
!    name of the 1 hour WFDE5 reanalysis file
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
  integer   :: tmpId, qId, windId, lwdId, psId, rainfId, snowfId, swdId
  integer   :: tindex
  logical   :: file_exists
  integer   :: c,r,t,k,l,iret
  integer   :: mo,rec_size
  logical   :: read_lnd
  logical   :: read_flag
  real      :: undef  ! input data FillValue
  character(len=100) :: infile
  
  real, allocatable      :: ps(:,:,:)
  real, allocatable      :: tair(:,:,:)
  real, allocatable      :: qair(:,:,:)
  real, allocatable      :: wind(:,:,:)
  real, allocatable      :: rainf(:,:,:)
  real, allocatable      :: snowf(:,:,:)
  real, allocatable      :: swd(:,:,:)
  real, allocatable      :: lwd(:,:,:)

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
     
     mo = LDT_rc%lnc(n)*LDT_rc%lnr(n)

! Read single layer file fields if lwd file exists:
     inquire(file=fname(8),exist=file_exists) 
     if(file_exists) then 
        if(order.eq.1) then
           wfde5_struc(n)%ps1     = LDT_rc%udef
           wfde5_struc(n)%tair1   = LDT_rc%udef
           wfde5_struc(n)%qair1   = LDT_rc%udef
           wfde5_struc(n)%wind1   = LDT_rc%udef
           wfde5_struc(n)%rainf1  = LDT_rc%udef
           wfde5_struc(n)%swd1    = LDT_rc%udef
           wfde5_struc(n)%lwd1    = LDT_rc%udef
        else
           wfde5_struc(n)%ps2     = LDT_rc%udef
           wfde5_struc(n)%tair2   = LDT_rc%udef
           wfde5_struc(n)%qair2   = LDT_rc%udef
           wfde5_struc(n)%wind2   = LDT_rc%udef
           wfde5_struc(n)%rainf2  = LDT_rc%udef
           wfde5_struc(n)%swd2    = LDT_rc%udef
           wfde5_struc(n)%lwd2    = LDT_rc%udef
        endif


        rec_size = days(month)*24 + 1

        infile = fname(1)
        write(LDT_logunit,*)'[INFO] Reading WFDE5 file (bookend,', order,' -',trim(infile), ')'
        call LDT_verify(nf90_open(path=trim(infile), mode=NF90_NOWRITE, &
             ncid=ftn), 'nf90_open failed in read_wfde5')

        allocate(ps(wfde5_struc(n)%ncold,wfde5_struc(n)%nrold,rec_size))
        
        call LDT_verify(nf90_inq_varid(ftn,'PSurf',psId), &
             'nf90_inq_varid failed for psurf in read_wfde5')
        call LDT_verify(nf90_get_var(ftn,psId, ps),&
             'nf90_get_var failed for ps in read_wfde5') 
        call LDT_verify(nf90_get_att(ftn, psId,'_FillValue', undef), &
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
        call LDT_verify(nf90_close(ftn), 'failed to close file in read_wfde5')

        infile = fname(2)
        write(LDT_logunit,*)'[INFO] Reading WFDE5 file (bookend,', order,' -',trim(infile), ')'
        call LDT_verify(nf90_open(path=trim(infile), mode=NF90_NOWRITE, &
             ncid=ftn), 'nf90_open failed in read_wfde5')
        allocate(tair(wfde5_struc(n)%ncold,wfde5_struc(n)%nrold,rec_size))
        call LDT_verify(nf90_inq_varid(ftn,'Tair',tmpId), &
             'nf90_inq_varid failed for Tair in read_wfde5')
        call LDT_verify(nf90_get_var(ftn,tmpId, tair),&
             'nf90_get_var failed for tair in read_wfde5') 
        call LDT_verify(nf90_get_att(ftn,tmpId,'_FillValue', undef), &
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
        call LDT_verify(nf90_close(ftn), 'failed to close file in read_wfde5')

        infile = fname(3)
        write(LDT_logunit,*)'[INFO] Reading WFDE5 file (bookend,', order,' -',trim(infile), ')'
        call LDT_verify(nf90_open(path=trim(infile), mode=NF90_NOWRITE, &
             ncid=ftn), 'nf90_open failed in read_wfde5')
        allocate(qair(wfde5_struc(n)%ncold,wfde5_struc(n)%nrold,rec_size))
        call LDT_verify(nf90_inq_varid(ftn,'Qair',qId), &
             'nf90_inq_varid failed for Qair in read_wfde5')
        call LDT_verify(nf90_get_var(ftn,qId, qair),&
             'nf90_get_var failed for qair in read_wfde5') 
        call LDT_verify(nf90_get_att(ftn,qId,'_FillValue', undef), &
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
        call LDT_verify(nf90_close(ftn), 'failed to close file in read_wfde5')

        infile = fname(4)
        write(LDT_logunit,*)'[INFO] Reading WFDE5 file (bookend,', order,' -',trim(infile), ')'
        call LDT_verify(nf90_open(path=trim(infile), mode=NF90_NOWRITE, &
             ncid=ftn), 'nf90_open failed in read_wfde5')
        allocate(wind(wfde5_struc(n)%ncold,wfde5_struc(n)%nrold,rec_size))
        call LDT_verify(nf90_inq_varid(ftn,'Wind',windId), &
             'nf90_inq_varid failed for Wind in read_wfde5')
        call LDT_verify(nf90_get_var(ftn,windId, wind),&
             'nf90_get_var failed for wind in read_wfde5') 
        call LDT_verify(nf90_get_att(ftn,windId,'_FillValue', undef), &
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
        call LDT_verify(nf90_close(ftn), 'failed to close file in read_wfde5')
        
        infile = fname(5)
        write(LDT_logunit,*)'[INFO] Reading WFDE5 file (bookend,', order,' -',trim(infile), ')'
        call LDT_verify(nf90_open(path=trim(infile), mode=NF90_NOWRITE, &
             ncid=ftn), 'nf90_open failed in read_wfde5')
        allocate(rainf(wfde5_struc(n)%ncold,wfde5_struc(n)%nrold,rec_size))
        call LDT_verify(nf90_inq_varid(ftn,'Rainf',rainfId), &
             'nf90_inq_varid failed for Rainf in read_wfde5')
        call LDT_verify(nf90_get_var(ftn,rainfId, rainf),&
             'nf90_get_var failed for rainf in read_wfde5') 
        call LDT_verify(nf90_get_att(ftn,rainfId,'_FillValue', undef), &
             'nf90_get_att failed for rainf:_FillValue in read_wfde5')
        call LDT_verify(nf90_close(ftn), 'failed to close file in read_wfde5')

        infile = fname(6)
        write(LDT_logunit,*)'[INFO] Reading WFDE5 file (bookend,', order,' -',trim(infile), ')'
        call LDT_verify(nf90_open(path=trim(infile), mode=NF90_NOWRITE, &
             ncid=ftn), 'nf90_open failed in read_wfde5')
        allocate(snowf(wfde5_struc(n)%ncold,wfde5_struc(n)%nrold,rec_size))
        call LDT_verify(nf90_inq_varid(ftn,'Snowf',snowfId), &
             'nf90_inq_varid failed for Snowf in read_wfde5')
        call LDT_verify(nf90_get_var(ftn,snowfId, snowf),&
             'nf90_get_var failed for snowf in read_wfde5') 
        call LDT_verify(nf90_get_att(ftn,snowfId,'_FillValue', undef), &
             'nf90_get_att failed for snowf:_FillValue in read_wfde5')
        call LDT_verify(nf90_close(ftn), 'failed to close file in read_wfde5')

        do r=1,wfde5_struc(n)%nrold
         do c=1,wfde5_struc(n)%ncold
           do l=1,rec_size
              if(rainf(c,r,l).ne.undef.and.&
                   snowf(c,r,l).ne.undef) then 
                 rainf(c,r,l) = rainf(c,r,l) + snowf(c,r,l)
              endif
           enddo
         enddo
        enddo
        deallocate(snowf)
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
        write(LDT_logunit,*)'[INFO] Reading WFDE5 file (bookend,', order,' -',trim(infile), ')'
        call LDT_verify(nf90_open(path=trim(infile), mode=NF90_NOWRITE, &
             ncid=ftn), 'nf90_open failed in read_wfde5')
        allocate(swd(wfde5_struc(n)%ncold,wfde5_struc(n)%nrold,rec_size))
        call LDT_verify(nf90_inq_varid(ftn,'SWdown',swdId), &
             'nf90_inq_varid failed for SWdown in read_wfde5')
        call LDT_verify(nf90_get_var(ftn,swdId, swd),&
             'nf90_get_var failed for swd in read_wfde5') 
        call LDT_verify(nf90_get_att(ftn,swdId,'_FillValue', undef), &
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
        call LDT_verify(nf90_close(ftn), 'failed to close file in read_wfde5')

        infile = fname(8)
        write(LDT_logunit,*)'[INFO] Reading WFDE5 file (bookend,', order,' -',trim(infile), ')'
        call LDT_verify(nf90_open(path=trim(infile), mode=NF90_NOWRITE, &
             ncid=ftn), 'nf90_open failed in read_wfde5')
        allocate(lwd(wfde5_struc(n)%ncold,wfde5_struc(n)%nrold,rec_size))
        call LDT_verify(nf90_inq_varid(ftn,'LWdown',lwdId), &
             'nf90_inq_varid failed for LWdown in read_wfde5')
        call LDT_verify(nf90_get_var(ftn,lwdId, lwd),&
             'nf90_get_var failed for lwd in read_wfde5')
        call LDT_verify(nf90_get_att(ftn,lwdId,'_FillValue', undef), &
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
        
        call LDT_verify(nf90_close(ftn), &
             'failed to close file in read_wfde5')
                
     else
        write(LDT_logunit,*) '[ERR] ',trim(fname(8))//' does not exist'
        call LDT_endrun()
        
     endif
  endif

  tindex = (day - 1)*24 + hour + 1
  
  if(order.eq.1) then 
     call assign_processed_wfde5forc(n,kk,order,1,wfde5_struc(n)%tair1(:,tindex))
     call assign_processed_wfde5forc(n,kk,order,2,wfde5_struc(n)%qair1(:,tindex))
     call assign_processed_wfde5forc(n,kk,order,3,wfde5_struc(n)%swd1(:,tindex))
     call assign_processed_wfde5forc(n,kk,order,4,wfde5_struc(n)%lwd1(:,tindex))
     call assign_processed_wfde5forc(n,kk,order,5,wfde5_struc(n)%wind1(:,tindex))
     call assign_processed_wfde5forc(n,kk,order,6,wfde5_struc(n)%ps1(:,tindex))
     call assign_processed_wfde5forc(n,kk,order,7,wfde5_struc(n)%rainf1(:,tindex))
  else
     call assign_processed_wfde5forc(n,kk,order,1,wfde5_struc(n)%tair2(:,tindex))
     call assign_processed_wfde5forc(n,kk,order,2,wfde5_struc(n)%qair2(:,tindex))
     call assign_processed_wfde5forc(n,kk,order,3,wfde5_struc(n)%swd2(:,tindex))
     call assign_processed_wfde5forc(n,kk,order,4,wfde5_struc(n)%lwd2(:,tindex))
     call assign_processed_wfde5forc(n,kk,order,5,wfde5_struc(n)%wind2(:,tindex))
     call assign_processed_wfde5forc(n,kk,order,6,wfde5_struc(n)%ps2(:,tindex))
     call assign_processed_wfde5forc(n,kk,order,7,wfde5_struc(n)%rainf2(:,tindex))
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
  use LDT_coreMod
  use LDT_logMod
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
  real,    intent(out)   :: output_var(LDT_rc%lnc(n)*LDT_rc%lnr(n))
  logical, intent(in)    :: pcp_flag
  real,    intent(in)    :: undef

!
! !DESCRIPTION: 
!  This subroutine spatially interpolates a WFDE5 field
!  to the LDT running domain
! 
!EOP

  integer   :: t,c,r,k,iret
  integer   :: doy
  integer   :: ftn
  integer   :: pcp1Id, pcp2Id, pcp3Id, pcp4Id,pcp5Id, pcp6Id
  real      :: f (wfde5_struc(n)%ncold*wfde5_struc(n)%nrold)
  logical*1 :: lb(wfde5_struc(n)%ncold*wfde5_struc(n)%nrold)
  logical*1 :: lo(LDT_rc%lnc(n)*LDT_rc%lnr(n))
  integer   :: input_size
! _____________________________________________________________

  input_size = wfde5_struc(n)%ncold*wfde5_struc(n)%nrold
  output_var = LDT_rc%udef

!-----------------------------------------------------------------------    
! Apply corrections
!-----------------------------------------------------------------------  
     
  lb = .false.
  do r=1,wfde5_struc(n)%nrold
     do c=1,wfde5_struc(n)%ncold           
        k= k + 1
        if(input_var(c,r).ne.undef) then 
           f(k) = input_var(c,r)
           lb(k) = .true.
        else
           f(k) = LDT_rc%udef
           lb(k) = .false.
        endif
     enddo
  enddo
  
!-----------------------------------------------------------------------    
! Apply downscaling
!-----------------------------------------------------------------------    
     

  if(pcp_flag.and.&
       trim(LDT_rc%met_gridtransform(findex)).eq."budget-bilinear") then 
     
     call conserv_interp(LDT_rc%gridDesc(n,:),lb,f,lo,&
          output_var(:), &
          wfde5_struc(n)%mi,LDT_rc%lnc(n)*LDT_rc%lnr(n),& 
          LDT_domain(n)%lat, LDT_domain(n)%lon,&
          wfde5_struc(n)%w112,wfde5_struc(n)%w122,&
          wfde5_struc(n)%w212,wfde5_struc(n)%w222,&
          wfde5_struc(n)%n112,wfde5_struc(n)%n122,&
          wfde5_struc(n)%n212,wfde5_struc(n)%n222,&
          LDT_rc%udef, iret)
     
  elseif(trim(LDT_rc%met_gridtransform(findex)).eq."bilinear".or.&
       trim(LDT_rc%met_gridtransform(findex)).eq."budget-bilinear") then 

     call bilinear_interp(LDT_rc%gridDesc(n,:),lb,f,lo,&
          output_var(:), &
          wfde5_struc(n)%mi,LDT_rc%lnc(n)*LDT_rc%lnr(n), & 
          LDT_domain(n)%lat, LDT_domain(n)%lon,&
          wfde5_struc(n)%w111,wfde5_struc(n)%w121,&
          wfde5_struc(n)%w211,wfde5_struc(n)%w221,&
          wfde5_struc(n)%n111,wfde5_struc(n)%n121,&
          wfde5_struc(n)%n211,wfde5_struc(n)%n221,&
          LDT_rc%udef, iret)
     
  elseif(trim(LDT_rc%met_gridtransform(findex)).eq."neighbor") then 
     call neighbor_interp(LDT_rc%gridDesc(n,:),lb,f,lo,&
          output_var(:),wfde5_struc(n)%mi,&
          LDT_rc%lnc(n)*LDT_rc%lnr(n),&
          LDT_domain(n)%lat, LDT_domain(n)%lon,&
          wfde5_struc(n)%n113,LDT_rc%udef,iret)
  else
     write(LDT_logunit,*) '[ERR] Spatial interpolation option '//&
          trim(LDT_rc%met_gridtransform(findex))//&
          ' not supported for WFDE5'
     call LDT_endrun()
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
  use LDT_coreMod
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
  real    :: wfde5forc(LDT_rc%lnc(n)*LDT_rc%lnc(n))
  

  integer :: c,r

  do r=1,LDT_rc%lnr(n)
     do c=1,LDT_rc%lnc(n)
        if(LDT_domain(n)%gindex(c,r).ne.-1) then 
           if(order.eq.1) then 
              wfde5_struc(n)%metdata1(var_index,&
                   LDT_domain(n)%gindex(c,r)) = &
                   wfde5forc(c+(r-1)*LDT_rc%lnc(n))
           elseif(order.eq.2) then 
              wfde5_struc(n)%metdata2(var_index,&
                   LDT_domain(n)%gindex(c,r)) = &
                   wfde5forc(c+(r-1)*LDT_rc%lnc(n))
           endif
        endif
     enddo
  enddo
end subroutine assign_processed_wfde5forc






