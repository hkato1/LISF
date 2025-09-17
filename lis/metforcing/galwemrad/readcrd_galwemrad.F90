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
! !ROUTINE: readcrd_galwemrad
! \label{readcrd_galwemrad}
!
! !REVISION HISTORY:
! 11 Mar 2022; Yeosang Yoon, Initial Code
! 08 Sep 2022; Yeosang Yoon, Add codes to read GALWEM 25 DEG dataset
! 13 Jun 2024; Hiroko Beaudoing, adopted GALWEM routines for radiation 
!                                supplemental forcing option
!
! !INTERFACE:    
subroutine readcrd_galwemrad()
! !USES:
  use LIS_logMod
  use LIS_coreMod
  use galwemrad_forcingMod, only : galwemrad_struc
  use ESMF
!
! !DESCRIPTION:
!
!  This routine reads the options specific to GALWEM radiation forcing at 0.25 deg
!  from the LIS configuration file. 
!  
!EOP

  implicit none

  integer :: n,rc

  call ESMF_ConfigFindLabel(LIS_config,"GALWEM radiation forcing directory:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,galwemrad_struc(n)%odir,rc=rc)
     call LIS_verify(rc,'GALWEM radiation forcing directory: not defined')
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"GALWEM radiation run mode:",rc=rc)
  call LIS_verify(rc, 'GALWEM radiation run mode: not defined ')
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,galwemrad_struc(n)%runmode,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"GALWEM radiation resolution:",rc=rc)
  call LIS_verify(rc, 'GALWEM radiation resolution: not defined ')
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,galwemrad_struc(n)%resol,rc=rc)
  enddo

  do n=1,LIS_rc%nnest
     write(LIS_logunit,*) '[INFO] Using GALWEM radiation forcing'
     write(LIS_logunit,*) '[INFO] GALWEM radiation forcing directory: ', trim(galwemrad_struc(n)%odir)
     write(LIS_logunit,*) '[INFO] GALWEM radiation run mode: ',galwemrad_struc(n)%runmode
     write(LIS_logunit,*) '[INFO] GALWEM radiation resolution: ',galwemrad_struc(n)%resol
  enddo
end subroutine readcrd_galwemrad
