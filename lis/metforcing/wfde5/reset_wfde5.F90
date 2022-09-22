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
! !MODULE: reset_wfde5
! \label{reset_wfde5}
! 
! !REVISION HISTORY: 
! 25 Mar 2022: Hiroko Beaudoing, initial code 
! 
! !INTERFACE:
subroutine reset_wfde5
! !USES:
  use LIS_coreMod,  only : LIS_rc
  use LIS_timeMgrMod, only : LIS_date2time
  use wfde5_forcingMod
!
! !DESCRIPTION:
!  Routine to cleanup allocated structures for wfde5 forcing. 
!
!EOP  
  implicit none
  integer :: n 

  do n=1,LIS_rc%nnest
     wfde5_struc(n)%startFlag = .true. 
     wfde5_struc(n)%dayFlag = .true. 
     wfde5_struc(n)%wfde5time1 = 3000.0
     wfde5_struc(n)%wfde5time2 = 0.0
     wfde5_struc(n)%ringtime = 0.0
     wfde5_struc(n)%reset_flag = .true.
  enddo
end subroutine reset_wfde5
