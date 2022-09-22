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
!
! !ROUTINE: readcrd_wfde5
! \label{readcrd_wfde5}
!
! !REVISION HISTORY:
! 25 Mar 2022: Hiroko Beaudoing, initial code 
!
! !INTERFACE:    
subroutine readcrd_wfde5()
! !USES:
  use ESMF
  use LIS_coreMod,      only : LIS_rc, LIS_config
  use LIS_logMod
  use LIS_timeMgrMod,   only : LIS_calendar
  use wfde5_forcingMod, only : wfde5_struc
!
! !DESCRIPTION:
!
!  This routine reads the options specific to WFDE5 forcing
!  from the LIS configuration file. 
!  
!EOP
  implicit none

  integer :: n,t,rc
  integer :: status
  type(ESMF_Time)         :: LISstartTime

  call ESMF_ConfigFindLabel(LIS_config,"WFDE5 forcing directory:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,wfde5_struc(n)%wfde5dir,&
          rc=rc)
     call LIS_verify(rc,&
          'WFDE5 forcing directory: not defined')
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"WFDE5 forcing precip source:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,wfde5_struc(n)%presrc,&
          rc=rc)
     call LIS_verify(rc,&
          'WFDE5 forcing precip source: not defined')
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"WFDE5 surface altitude file:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,wfde5_struc(n)%wfde5alt_file,&
          rc=rc)
     call LIS_verify(rc,&
          'WFDE5 surface altitude file: not defined')
  enddo

  do n=1,LIS_rc%nnest
     write(LIS_logunit,*) '[INFO] Using WFDE5 forcing'
     write(LIS_logunit,*) '[INFO] WFDE5 forcing directory: ',&
          trim(wfde5_struc(n)%wfde5DIR)

     wfde5_struc(n)%wfde5time1 = 3000.0
     wfde5_struc(n)%wfde5time2 = 0.0

     ! check if the start time needs to be adjusted
     call ESMF_TimeSet(wfde5_struc(n)%startTime, yy=1979, &
           mm = 1, dd = 1, h=7, m = 0, s=0, calendar=LIS_calendar, rc=status)
     call LIS_verify(status, &
                 'readcard_wfde5: ESMF_TimeSet wfde5_struc%startTime')

     call ESMF_TimeSet(LISstartTime, yy=LIS_rc%syr, &
           mm = LIS_rc%smo, dd = LIS_rc%sda, h=LIS_rc%shr, &
           m = LIS_rc%smn, s=LIS_rc%sss, calendar=LIS_calendar, &
           rc=status)
     call LIS_verify(status, 'readcrd_wfde5: ESMF_TimeSet LISstartTime')

     if (LISstartTime .lt. wfde5_struc(n)%startTime) then
       write(LIS_logunit,*) &
             'NOTE: WFDE5 data not available till January 1, 1979 0700 GMT'
       write(LIS_logunit,*) 'Please modify lis.config'
       write(LIS_logunit,*) 'Program stopping ... '
       call LIS_endrun()
     endif

  enddo
end subroutine readcrd_wfde5
