!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: vic412_getsws_hymap2
!  \label{vic412_getsws_hymap2}
!
! !REVISION HISTORY:
! 12 Sep 2019: Augusto Getirana; implementation of two-way coupling
!
! !INTERFACE:
subroutine vic412_getsws_hymap2(n)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc, LIS_masterproc
  use LIS_routingMod, only : LIS_runoff_state
  use LIS_logMod
  use LIS_historyMod
  use vic412_lsmMod, only : vic412_struc

  implicit none
! !ARGUMENTS: 
  integer,  intent(in)   :: n 
!
! !DESCRIPTION:
!   This routine defines the surface water storage variables in NoahMP
!   to be updated based on feedback from HYMAP2
!  
!EOP
  integer                :: t
  integer                :: c,r
  integer                :: status
  integer                :: enable2waycpl
  
  enable2waycpl = 0 

  if(enable2waycpl==1) then 
     write(LIS_logunit,*) '[ERR] The vic412_getsws_hymap2 is not implemented'
     call LIS_endrun()
  endif

end subroutine vic412_getsws_hymap2
