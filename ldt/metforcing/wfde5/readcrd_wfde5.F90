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
! 29 Mar 2022: Hiroko Beaudoing, initial code 
!
! !INTERFACE:    
subroutine readcrd_wfde5()
! !USES:
  use ESMF
  use LDT_coreMod, only : LDT_rc, LDT_config
  use LDT_logMod
  use wfde5_forcingMod, only : wfde5_struc
!
! !DESCRIPTION:
!
!  This routine reads the options specific to WFDE5 forcing
!  from the LDT configuration file. 
!  
!EOP
  implicit none

  integer :: n,t,rc

  call ESMF_ConfigFindLabel(LDT_config,"WFDE5 forcing directory:",rc=rc)
  do n=1,LDT_rc%nnest
     call ESMF_ConfigGetAttribute(LDT_config,wfde5_struc(n)%wfde5dir,&
          rc=rc)
     call LDT_verify(rc,&
          'WFDE5 forcing directory: not defined')
  enddo

  call ESMF_ConfigFindLabel(LDT_config,"WFDE5 forcing precip source:",rc=rc)
  do n=1,LDT_rc%nnest
     call ESMF_ConfigGetAttribute(LDT_config,wfde5_struc(n)%presrc,&
          rc=rc)
     call LDT_verify(rc,&
          'WFDE5 forcing precip source: not defined')
  enddo

  call ESMF_ConfigFindLabel(LDT_config,"WFDE5 surface altitude file:",rc=rc)
  do n=1,LDT_rc%nnest
     call ESMF_ConfigGetAttribute(LDT_config,wfde5_struc(n)%wfde5alt_file,&
          rc=rc)
     call LDT_verify(rc,&
          'WFDE5 surface altitude file: not defined')
  enddo

  do n=1,LDT_rc%nnest
     write(LDT_logunit,*) '[INFO] Using WFDE5 forcing'
     write(LDT_logunit,*) '[INFO] WFDE5 forcing directory: ',&
          wfde5_struc(n)%wfde5dir

     wfde5_struc(n)%wfde5time1 = 3000.0
     wfde5_struc(n)%wfde5time2 = 0.0

  enddo
end subroutine readcrd_wfde5
