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
! !MODULE: reset_galwemrad
!  \label{reset_galwemrad}
!
! !REVISION HISTORY:
! 05 Apr 2022; Yeosang Yoon, Initial Code
!
! !INTERFACE:
subroutine reset_galwemrad()
! !USES:
  use LIS_coreMod,       only : LIS_rc
  use galwemrad_forcingMod, only : galwemrad_struc
!
! !DESCRIPTION:
!  Routine to reset GALWEM forcing related memory allocations.
!
!EOP
  implicit none

  integer   :: n
  integer   :: findex

  do n=1,LIS_rc%nnest
     galwemrad_struc(n)%fcsttime1 = 3000.0
     galwemrad_struc(n)%fcsttime2 = 0.0
  enddo

end subroutine reset_galwemrad
