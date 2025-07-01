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
! !ROUTINE: timeinterp_agrradps
! \label{timeinterp_agrradps}
! 
! !REVISION HISTORY: 
!  20 Jan 2006:   Yudong Tian: Initial Implementation
!  10 Oct 2006:   Sujay Kumar: Switched to using ESMF_State for storing
!                              forcing data. 
! !INTERFACE:
subroutine timeinterp_agrradps(n,findex)
! !USES:
    use ESMF
    use LIS_coreMod,         only : LIS_rc, LIS_domain, LIS_localPet
    use LIS_FORC_AttributesMod 
    use LIS_metforcingMod,   only : LIS_FORC_Base_State, LIS_forc
    use LIS_constantsMod
    use LIS_logMod
    use agrradps_forcingMod, only : agrradps_struc

    implicit none
! !ARGUMENTS: 
    integer, intent(in) :: n
    integer, intent(in) :: findex

! !DESCRIPTION:
!  Temporally interpolates the forcing data to the current model 
!  timestep. 
!  All variables except precipitation is linearly interpolated. 
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!   index of the nest
!  \end{description}
!
!EOP
    integer :: t,index
    integer :: f
    integer :: tid, tid1, tid2
    real    :: wt1,wt2
    
    integer          :: status
    type(ESMF_Field) :: swdField
    type(ESMF_Field) :: lwdField
    real, pointer    :: swd(:),lwd(:)
!EOP
  
  wt1 = (agrradps_struc(n)%agrtime2-LIS_rc%time) / & 
        (agrradps_struc(n)%agrtime2-agrradps_struc(n)%agrtime1)
  wt2 = 1.0 - wt1
  

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),    &
                     trim(LIS_FORC_SWdown%varname(1)), &
                     swdField,                         &
                     rc=status)
  call LIS_verify(status)

  call ESMF_StateGet(LIS_FORC_Base_State(n,findex),    &
                     trim(LIS_FORC_LWdown%varname(1)), &
                     lwdField,                         &
                     rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(swdField,localDE=0,farrayPtr=swd,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(lwdField,localDE=0,farrayPtr=lwd,rc=status)
  call LIS_verify(status)

  do t = 1,LIS_rc%ntiles(n)
     index = LIS_domain(n)%tile(t)%index
     if((agrradps_struc(n)%metdata1(1,index).ne.LIS_rc%udef) .and. &
          (agrradps_struc(n)%metdata2(1,index).ne.LIS_rc%udef)) then  
        swd(t) = &
             wt1 * agrradps_struc(n)%metdata1(1,index) +  & 
             wt2 *agrradps_struc(n)%metdata2(1,index)

           if (swd(t).gt.LIS_CONST_SOLAR) then
              write(unit=LIS_logunit,fmt=*) &
                   '[WARN] sw radiation too high in AGRMET at tile',t
              write(unit=LIS_logunit,fmt=*)'[WARN] it is',swd(t)
              write(unit=LIS_logunit,fmt=*)'[WARN] data1=',&
                   agrradps_struc(n)%metdata1(1,index) 
              write(unit=LIS_logunit,fmt=*)'[WARN] data2=',&
                   agrradps_struc(n)%metdata2(1,index)
              write(unit=LIS_logunit,fmt=*)'[WARN] wt1=',wt1,'wt2=',wt2
              swd(t) = LIS_CONST_SOLAR
              write(unit=LIS_logunit,fmt=*)'[WARN] set to ',swd(t)
           endif
           if ((swd(t).ne.LIS_rc%udef).and.(swd(t).lt.0)) then
              if (swd(t).gt.-0.00001) then
                 swd(t) = 0.0
              else
                 write(LIS_logunit,*) &
                   '[ERR] timeinterp_agrradps -- Stopping because ', &
                   'SWdown forcing not udef but lt 0,'
                 write(LIS_logunit,*)'[ERR] timeinterp_agrradsp -- ', &
                   t,swd(t),agrradps_struc(n)%metdata2(1,index), &
                   ' (',LIS_localPet,')'
                 call LIS_endrun
              endif
           endif
        
     else
        swd(t) = LIS_rc%udef
     endif
  enddo
  
  do t = 1,LIS_rc%ntiles(n)
     index = LIS_domain(n)%tile(t)%index
     if((agrradps_struc(n)%metdata1(2,index).ne.LIS_rc%udef) .and. &
          (agrradps_struc(n)%metdata2(2,index).ne.LIS_rc%udef)) then  
        lwd(t) = &
             wt1 * agrradps_struc(n)%metdata1(2,index) +  & 
             wt2 * agrradps_struc(n)%metdata2(2,index)
        ! sanity check
        if (lwd(t).gt.750.0 .or. ((lwd(t).ne.LIS_rc%udef).and.(lwd(t).lt.0))) then
           write(LIS_logunit,*) &
            '[ERR] timeinterp_agrradps -- Stopping because ', &
            'LWdown is out of expected range,'
           write(LIS_logunit,*)'[ERR] timeinterp_agrradsp -- ', &
             t,lwd(t),agrradps_struc(n)%metdata1(2,index), &
             ' (',LIS_localPet,')',agrradps_struc(n)%metdata2(2,index)
           call LIS_endrun
        endif
     else
        lwd(t) = LIS_rc%udef
     endif
  enddo

end subroutine timeinterp_agrradps
