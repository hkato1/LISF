!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: readcrd_gpcp
! \label{readcrd_gpcp}
!
! !REVISION HISTORY:
! 11 Dec 2003; Sujay Kumar, Initial Code
! 3 Jun 2016; Hiroko Beaudoing, adopted for GPCP 3-hourly data
!
! !INTERFACE:    
subroutine readcrd_gpcp()
! !USES:
  use ESMF 
  use gpcp_forcingMod, only : gpcp_struc
  use LIS_coreMod, only : LIS_config,LIS_rc
  use LIS_logMod, only : LIS_logunit

!
! !DESCRIPTION:
!
!  This routine reads the options specific to GPCP forcing from 
!  the LIS configuration file. 
!  
!EOP
  implicit none
  
  integer :: n,rc

  call ESMF_ConfigFindLabel(LIS_config, "GPCP forcing directory:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,gpcp_struc(n)%gpcpdir,rc=rc)
  enddo

  do n=1,LIS_rc%nnest
     write(LIS_logunit,*)"Using GPCP forcing"
     write(LIS_logunit,*)" GPCP forcing directory: ",trim(gpcp_struc(n)%gpcpdir)
!------------------------------------------------------------------------
! Setting global observed precip times to zero to ensure 
! data is read in during first time step
!------------------------------------------------------------------------
     gpcp_struc(n)%gpcptime = 0.0
  enddo

  gpcp_struc(:)%ncold = 384
  gpcp_struc(:)%nrold = 190

end subroutine readcrd_gpcp
