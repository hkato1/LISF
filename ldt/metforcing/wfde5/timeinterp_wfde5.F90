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
! !ROUTINE: timeinterp_wfde5
! \label{timeinterp_wfde5}
!
! !REVISION HISTORY:
! 29 mar 2022: Hiroko Beaudoing, initial code 
!
! !INTERFACE:

subroutine timeinterp_wfde5(n,findex)

! !USES:
  use ESMF
  use LDT_coreMod,          only : LDT_rc, LDT_domain, LDT_localPet
  use LDT_metforcingMod,    only : LDT_forc, LDT_FORC_Base_State
  use LDT_FORC_AttributesMod
  use LDT_constantsMod,     only : LDT_CONST_SOLAR
  use LDT_timeMgrMod,       only : LDT_tick
  use LDT_logMod,           only : LDT_logunit, LDT_verify, LDT_endrun
  use wfde5_forcingMod,     only : wfde5_struc

  implicit none

! !ARGUMENTS:
  integer, intent(in):: n
  integer, intent(in):: findex
!
! !DESCRIPTION:
!  Temporally interpolates the forcing data to the current model 
!  timestep. Downward shortwave radiation is interpolated using a
!  zenith-angled based approach. Precipitation 
!  is not temporally interpolated, and the previous 1 hourly value
!  is used. All other variables are linearly interpolated between 
!  the 1 hourly blocks. 
! 
!  The routines invoked are: 
!  \begin{description}
!   \item[LDT\_time2date](\ref{LDT_time2date}) \newline
!    converts the time to a date format
!   \item[zterp](\ref{zterp}) \newline
!    zenith-angle based interpolation
!  \end{description}
!EOP
  integer :: t,zdoy,k,kk
  integer :: index1
  integer :: bdoy,byr,bmo
  integer :: bda,bhr,bmn
  integer :: bss
  real*8  :: btime
  real    :: wt1,wt2,czb,cze,czm,gmt1,gmt2
  real    :: zw1,zw2,bts
  integer            :: status
  integer            :: mfactor,m
  type(ESMF_Field)   :: tmpField,q2Field,uField,vField,swdField,lwdField
  type(ESMF_Field)   :: psurfField,pcpField,snowfField
  real,pointer       :: tmp(:),q2(:),uwind(:),vwind(:)
  real,pointer       :: swd(:),lwd(:),psurf(:),pcp(:)
  real,pointer       :: snowf(:)

  btime=wfde5_struc(n)%wfde5time1
  byr = LDT_rc%yr
  bmo = LDT_rc%mo
  bda = LDT_rc%da
  bhr = LDT_rc%hr
  bmn = 30
  bss = 0
  if (LDT_rc%mn.lt.30) then
     bts = -(60*60)
  else
     bts = 0
  endif
  call LDT_tick(btime,bdoy,gmt1,byr,bmo,bda,bhr,bmn,bss,bts)
  btime=wfde5_struc(n)%wfde5time2
  byr=LDT_rc%yr    !next hour
  bmo=LDT_rc%mo
  bda=LDT_rc%da
  bhr=LDT_rc%hr
  bmn=30
  bss=0
  if (LDT_rc%mn.lt.30) then
     bts=0
  else
     bts=60*60
  endif
  call LDT_tick(btime,bdoy,gmt2,byr,bmo,bda,bhr,bmn,bss,bts)

!-----------------------------------------------------------------------
!  Interpolate Data in Time
!-----------------------------------------------------------------------

  wt1=(wfde5_struc(n)%wfde5time2-LDT_rc%time)/ & 
      (wfde5_struc(n)%wfde5time2-wfde5_struc(n)%wfde5time1)
  wt2=1.0-wt1

  call ESMF_StateGet(LDT_FORC_Base_State(n,findex),LDT_FORC_Tair%varname(1),tmpField,&
       rc=status)
  call LDT_verify(status, 'Error: Enable Tair in the forcing variables list')

  call ESMF_StateGet(LDT_FORC_Base_State(n,findex),LDT_FORC_Qair%varname(1),q2Field,&
       rc=status)
  call LDT_verify(status, 'Error: Enable Qair in the forcing variables list')

  call ESMF_StateGet(LDT_FORC_Base_State(n,findex),LDT_FORC_SWdown%varname(1),swdField,&
       rc=status)
  call LDT_verify(status, 'Error: Enable SWdown in the forcing variables list')

  call ESMF_StateGet(LDT_FORC_Base_State(n,findex),LDT_FORC_LWdown%varname(1),lwdField,&
       rc=status)
  call LDT_verify(status, 'Error: Enable LWdown in the forcing variables list')

  call ESMF_StateGet(LDT_FORC_Base_State(n,findex),LDT_FORC_Wind_E%varname(1),uField,&
       rc=status)
  call LDT_verify(status, 'Error: Enable Wind_E in the forcing variables list')

  call ESMF_StateGet(LDT_FORC_Base_State(n,findex),LDT_FORC_Wind_N%varname(1),vField,&
       rc=status)
  call LDT_verify(status, 'Error: Enable Wind_N in the forcing variables list')

  call ESMF_StateGet(LDT_FORC_Base_State(n,findex),LDT_FORC_Psurf%varname(1),psurfField,&
       rc=status)
  call LDT_verify(status, 'Error: Enable Psurf in the forcing variables list')

  call ESMF_StateGet(LDT_FORC_Base_State(n,findex),LDT_FORC_Rainf%varname(1),pcpField,&
       rc=status)
  call LDT_verify(status, 'Error: Enable Rainf in the forcing variables list')

  call ESMF_StateGet(LDT_FORC_Base_State(n,findex),&
       LDT_FORC_Snowf%varname(1),snowfField,&
       rc=status)
  call LDT_verify(status, &
       'Error: Enable Snowf in the forcing variables list')

  call ESMF_FieldGet(tmpField,localDE=0,farrayPtr=tmp,rc=status)
  call LDT_verify(status)

  call ESMF_FieldGet(q2Field,localDE=0,farrayPtr=q2,rc=status)
  call LDT_verify(status)
  
  call ESMF_FieldGet(swdField,localDE=0,farrayPtr=swd,rc=status)
  call LDT_verify(status)

  call ESMF_FieldGet(lwdField,localDE=0,farrayPtr=lwd,rc=status)
  call LDT_verify(status)

  call ESMF_FieldGet(uField,localDE=0,farrayPtr=uwind,rc=status)
  call LDT_verify(status)

  call ESMF_FieldGet(vField,localDE=0,farrayPtr=vwind,rc=status)
  call LDT_verify(status)

  call ESMF_FieldGet(psurfField,localDE=0,farrayPtr=psurf,rc=status)
  call LDT_verify(status)

  call ESMF_FieldGet(pcpField,localDE=0,farrayPtr=pcp,rc=status)
  call LDT_verify(status)

  call ESMF_FieldGet(snowfField,localDE=0,farrayPtr=snowf,rc=status)
  call LDT_verify(status)

  do t=1,LDT_rc%ntiles(n)
     index1 = LDT_domain(n)%tile(t)%index
     zdoy=LDT_rc%doy
     call zterp(1,LDT_domain(n)%grid(index1)%lat,&
          LDT_domain(n)%grid(index1)%lon,&
          gmt1,gmt2,LDT_rc%gmt,zdoy,zw1,zw2,czb,cze,czm,LDT_rc)
     
     if (wfde5_struc(n)%metdata1(3,index1).ne.LDT_rc%udef.and.&
          wfde5_struc(n)%metdata2(3,index1).ne.LDT_rc%udef) then 
        swd(t) = wfde5_struc(n)%metdata1(3,index1)*zw1+&
             wfde5_struc(n)%metdata2(3,index1)*zw2
        if (swd(t).gt.LDT_CONST_SOLAR) then
           write(unit=LDT_logunit,fmt=*) &
                   '[WARN] sw radiation too high in WFDE5!!'
           write(unit=LDT_logunit,fmt=*)'[WARN] it is',swd(t)
           write(unit=LDT_logunit,fmt=*)'[WARN] wfde5data1=',&
                wfde5_struc(n)%metdata1(3,index1)
           write(unit=LDT_logunit,fmt=*)'[WARN] wfde5data2=',&
                wfde5_struc(n)%metdata2(3,index1)
           write(unit=LDT_logunit,fmt=*)'[WARN] zw1=',zw1,'zw2=',zw2
           swd(t) = LDT_CONST_SOLAR
           write(unit=LDT_logunit,fmt=*)'[WARN] forcing set to ',swd(t) 
        endif
     endif
     
     if ((swd(t).ne.LDT_rc%udef).and.(swd(t).lt.0)) then
        if (swd(t).gt.-0.00001) then 
           swd(t) = 0.0 
        else
           write(LDT_logunit,*) &
                '[ERR] timeinterp_wfde5 -- Stopping because ', & 
                'forcing not udef but lt0,'
           write(LDT_logunit,*)'[ERR] timeinterp_wfde5 -- ', & 
                t,swd(t),wfde5_struc(n)%metdata2(3,index1), & 
                ' (',LDT_localPet,')'
           call LDT_endrun
        endif
     endif
     
     if (swd(t).gt.LDT_CONST_SOLAR) then
        swd(t)=wfde5_struc(n)%metdata2(3,index1)
     endif
  enddo
 
!-----------------------------------------------------------------------
! precip variable - constant rate over the WFDE5 hour
!-----------------------------------------------------------------------  
  do t=1,LDT_rc%ntiles(n)
     index1 = LDT_domain(n)%tile(t)%index
     pcp(t) = wfde5_struc(n)%metdata1(7,index1)
     
     if ( pcp(t) < 0 ) then
        pcp(t) = 0
     endif
     snowf(t) = wfde5_struc(n)%metdata1(8,index1)
  enddo

!-----------------------------------------------------------------------
! Linearly interpolate everything else 
!-----------------------------------------------------------------------

  do t=1,LDT_rc%ntiles(n)
     index1 = LDT_domain(n)%tile(t)%index
     
     tmp(t) = wfde5_struc(n)%metdata1(1,index1)*wt1 + &
          wfde5_struc(n)%metdata2(1,index1)*wt2
     q2(t)  = wfde5_struc(n)%metdata1(2,index1)*wt1 + &
          wfde5_struc(n)%metdata2(2,index1)*wt2
     lwd(t)  = wfde5_struc(n)%metdata1(4,index1)*wt1 + &
          wfde5_struc(n)%metdata2(4,index1)*wt2
     
     uwind(t) = wfde5_struc(n)%metdata1(5,index1)*wt1+&
          wfde5_struc(n)%metdata2(5,index1)*wt2
     vwind(t) = 0.0
     psurf(t) = wfde5_struc(n)%metdata1(6,index1)*wt1 + &
          wfde5_struc(n)%metdata2(6,index1)*wt2
     
  enddo

end subroutine timeinterp_wfde5
  
