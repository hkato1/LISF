!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: finalize_gpcp
! \label{finalize_gpcp}
! 
! !REVISION HISTORY: 
!  06 Jun 2016: Hiroko Beaudoing; Adopted CMAP routine to GPCP.
! 
! !INTERFACE:
subroutine finalize_gpcp(findex)

! !USES:
  use LIS_coreMod,     only : LIS_rc
  use gpcp_forcingMod, only : gpcp_struc
! !DESCRIPTION: 
!  Routine to cleanup GPCP forcing related memory allocations.   
!
!EOP
  implicit none
  
  integer, intent(in) :: findex

  integer :: n
  integer :: rc

  do n=1,LIS_rc%nnest
     deallocate(gpcp_struc(n)%n111, stat=rc)
     deallocate(gpcp_struc(n)%n121, stat=rc)
     deallocate(gpcp_struc(n)%n211, stat=rc)
     deallocate(gpcp_struc(n)%n221, stat=rc)
     deallocate(gpcp_struc(n)%w111, stat=rc)
     deallocate(gpcp_struc(n)%w121, stat=rc)
     deallocate(gpcp_struc(n)%w211, stat=rc)
     deallocate(gpcp_struc(n)%w221, stat=rc)

     deallocate(gpcp_struc(n)%n112, stat=rc)
     deallocate(gpcp_struc(n)%n122, stat=rc)
     deallocate(gpcp_struc(n)%n212, stat=rc)
     deallocate(gpcp_struc(n)%n222, stat=rc)
     deallocate(gpcp_struc(n)%w112, stat=rc)
     deallocate(gpcp_struc(n)%w122, stat=rc)
     deallocate(gpcp_struc(n)%w212, stat=rc)
     deallocate(gpcp_struc(n)%w222, stat=rc)
  enddo

  deallocate(gpcp_struc)

end subroutine finalize_gpcp
