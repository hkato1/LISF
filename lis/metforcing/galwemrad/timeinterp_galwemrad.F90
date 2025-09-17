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
! !ROUTINE: timeinterp_galwemrad
! \label{timeinterp_galwemrad}
!
! !REVISION HISTORY:
!
! 05 Apr 2022: Yeosang Yoon, Initial specification
! 13 Jun 2025: Hiroko Beaudoing, adipted GALWEM forecast implenemtation for
!                                radiation supplemental forcing option.
!
! !INTERFACE:
subroutine timeinterp_galwemrad(n,findex)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_constantsMod
  use LIS_metforcingMod
  use LIS_FORC_AttributesMod
  use LIS_timeMgrMod
  use LIS_logMod
  use galwemrad_forcingMod
  use LIS_forecastMod

  implicit none
! !ARGUMENTS:
  integer, intent(in) :: n
  integer, intent(in) :: findex
!
! !DESCRIPTION:
!  Temporally interpolates the forcing data to the current model
!  timestep. Downward shortwave radiation is interpolated using a
!  zenith-angled based approach. Precipitation and longwave radiation
!  are not temporally interpolated, and the previous value is used.
!  All other variables are linearly interpolated between the blocks.
!
!  The routines invoked are:
!  \begin{description}
!   \item[LIS\_time2date](\ref{LIS_time2date}) \newline
!    converts the time to a date format
!   \item[LIS\_tick](\ref{LIS_tick}) \newline
!    advances or retracts time by the specified amount
!   \item[zterp](\ref{zterp}) \newline
!    zenith-angle based interpolation
!  \end{description}
!EOP
  integer :: zdoy
  real    :: zw1, zw2
  real    :: czm, cze, czb
  real    :: wt1, wt2,swt1,swt2
  real    :: gmt1, gmt2, tempbts
  integer :: t,index1
  integer :: bdoy,byr,bmo,bda,bhr,bmn
  real*8  :: btime,newtime1,newtime2
  real    :: tempgmt1,tempgmt2
  integer :: tempbdoy,tempbyr,tempbmo,tempbda,tempbhr,tempbmn
  integer :: tempbss
  integer            :: status
  type(ESMF_Field)   :: swdField,lwdField
  real,pointer       :: swd(:),lwd(:)
! ________________________________________

  btime=galwemrad_struc(n)%fcsttime1
  call LIS_time2date(btime,bdoy,gmt1,byr,bmo,bda,bhr,bmn)

  tempbdoy=bdoy
  tempgmt1=gmt1
  tempbyr=byr
  tempbmo=bmo
  tempbda=bda
  tempbhr=bhr
  if (tempbhr.eq.24) tempbhr=0
  tempbmn=bmn
  tempbss=0
  tempbts=0
  call LIS_tick(newtime1,tempbdoy,tempgmt1,&
       tempbyr,tempbmo,tempbda,tempbhr,tempbmn, &
       tempbss,tempbts)

  btime=galwemrad_struc(n)%fcsttime2
  call LIS_time2date(btime,bdoy,gmt2,byr,bmo,bda,bhr,bmn)
  tempbdoy=bdoy
  tempgmt2=gmt2
  tempbyr=byr
  tempbmo=bmo
  tempbda=bda
  tempbhr=bhr
  if (tempbhr.eq.24) tempbhr=0
  tempbmn=bmn
  tempbss=0
  tempbts=0
  call LIS_tick(newtime2,tempbdoy,tempgmt2,&
       tempbyr,tempbmo,tempbda,tempbhr,tempbmn,&
       tempbss,tempbts)

!  Interpolate Data in Time
  wt1=(galwemrad_struc(n)%fcsttime2-LIS_rc%time)/ &
       (galwemrad_struc(n)%fcsttime2-galwemrad_struc(n)%fcsttime1)
  wt2=1.0-wt1
  swt1=(newtime2-LIS_rc%time)/(newtime2-newtime1)
  swt2=1.0-swt1


  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_SWdown%varname(1),swdField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable SWdown in the forcing variables list')

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),LIS_FORC_LWdown%varname(1),lwdField,&
       rc=status)
  call LIS_verify(status, 'Error: Enable LWdown in the forcing variables list')

  call ESMF_FieldGet(swdField,localDE=0,farrayPtr=swd,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(lwdField,localDE=0,farrayPtr=lwd,rc=status)
  call LIS_verify(status)

  ! Downward shortwave radiation (average):
  do t=1,LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index
     zdoy=LIS_rc%doy

     ! compute and apply zenith angle weights
     call zterp(1,LIS_domain(n)%grid(index1)%lat,&
          LIS_domain(n)%grid(index1)%lon,&
          gmt1,gmt2,LIS_rc%gmt,zdoy,zw1,zw2,czb,cze,czm,LIS_rc)

     if (galwemrad_struc(n)%metdata1(1,index1).ne.LIS_rc%udef.and.&
          galwemrad_struc(n)%metdata2(1,index1).ne.LIS_rc%udef) then
        swd(t) = galwemrad_struc(n)%metdata1(1,index1)*zw1+&
             galwemrad_struc(n)%metdata2(1,index1)*zw2

        ! In cases of small cos(zenith) angles, use linear weighting to avoid overly large weights
        if((swd(t).gt.galwemrad_struc(n)%metdata1(1,index1).and. &
             swd(t).gt.galwemrad_struc(n)%metdata2(1,index1)).and. &
             (czb.lt.0.1.or.cze.lt.0.1))then
           swd(t) = galwemrad_struc(n)%metdata1(1,index1)*swt1+ &
                galwemrad_struc(n)%metdata2(1,index1)*swt2
        endif

        if (swd(t).gt.LIS_CONST_SOLAR) then
           write(unit=LIS_logunit,fmt=*)'[WARN] sw radiation too high!!'
           write(unit=LIS_logunit,fmt=*)'[WARN] it is',swd(t)
           write(unit=LIS_logunit,fmt=*)'[WARN] galwemraddata1=',&
                galwemrad_struc(n)%metdata1(1,index1)
           write(unit=LIS_logunit,fmt=*)'[WARN] galwemraddata2=',&
                galwemrad_struc(n)%metdata2(1,index1)
           write(unit=LIS_logunit,fmt=*)'[WARN] zw1=',zw1,'zw2=',zw2
           swd(t) = LIS_CONST_SOLAR
           write(unit=LIS_logunit,fmt=*)'[WARN] forcing set to ',swd(t)
        endif
     endif

     if ((swd(t).ne.LIS_rc%udef).and.(swd(t).lt.0)) then
        if (swd(t).gt.-0.00001) then
           swd(t) = 0.0
        else
           write(LIS_logunit,*) &
                '[ERR] timeinterp_galwemrad -- Stopping because ', &
                'forcing not udef but lt 0,'
           write(LIS_logunit,*)'[ERR] timeinterp_galwemrad -- ', &
                t,swd(t),galwemrad_struc(n)%metdata2(1,index1), &
                ' (',LIS_localPet,')'
           call LIS_endrun
        endif
     endif
  enddo

!-----------------------------------------------------------------------
! Linearly interpolate everything else
!-----------------------------------------------------------------------

  do t=1,LIS_rc%ntiles(n)
     index1 = LIS_domain(n)%tile(t)%index

     ! Downward longwave field
     if((galwemrad_struc(n)%metdata1(2,index1).ne.LIS_rc%udef).and.&
          (galwemrad_struc(n)%metdata2(2,index1).ne.LIS_rc%udef)) then
        lwd(t)  = galwemrad_struc(n)%metdata1(2,index1)*wt1 + &
                  galwemrad_struc(n)%metdata2(2,index1)*wt2
     endif
  enddo

end subroutine timeinterp_galwemrad
