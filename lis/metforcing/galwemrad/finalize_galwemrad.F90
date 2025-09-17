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
! !MODULE: finalize_galwemrad
! \label{finalize_galwemrad}
!
! !REVISION HISTORY:
! 04 Apr 2022; Yeosang Yoon, Initial Code
!
! !INTERFACE:
subroutine finalize_galwemrad
! !USES:
  use galwemrad_forcingMod, only : galwemrad_struc
!
! !DESCRIPTION:
!  Routine to cleanup allocated structures for GALWEM forcing.
!
!EOP
  implicit none

  deallocate(galwemrad_struc)

end subroutine finalize_galwemrad
