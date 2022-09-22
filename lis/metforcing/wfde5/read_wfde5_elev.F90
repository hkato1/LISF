!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
! !ROUTINE: read_wfde5_elev
!  \label{read_wfde5_elev}
!
! !REVISION HISTORY:
!
!  28 mar 2022; Hiroko Beaudoing; Initial Specificaton
!
! !INTERFACE:
subroutine read_wfde5_elev(n,findex)
! !USES:
  use LIS_coreMod
  use LIS_metforcingMod
  use LIS_logMod
  use wfde5_forcingMod
  use LIS_fileIOMod,     only : LIS_read_param
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  integer, intent(in) :: findex
! !DESCRIPTION:
!
!  Opens, reads, and interpolates WFDE5 model elevation to the LIS
!  grid. The data will be used to perform any topographical 
!  adjustments to the forcing.
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!   index of the nest
!  \end{description}
! 
!  The routines invoked are: 
!   \begin{description}
!    \item[ij\_to\_latlon](\ref{ij_to_latlon}) \newline
!     computes the lat lon values in LIS grid projection
!   \end{description}
!EOP
   integer :: i, err
   logical :: file_exists
   integer :: nid,elevId
   integer :: c,r
   real    :: elev(LIS_rc%gnc(n),LIS_rc%gnr(n))
   real    :: elev_subset(LIS_rc%lnc(n),LIS_rc%lnr(n))
   real    :: undef
  
#if (defined USE_NETCDF3 || defined USE_NETCDF4)

   inquire(file=wfde5_struc(n)%wfde5alt_file, exist=file_exists)
   if(file_exists) then 

      write(LIS_logunit,*) " Reading WFDE5 elevation data ... " 
      call LIS_read_param(n,"ELEV_WFDE5",elev)

      elev_subset(:,:) = elev(LIS_ews_halo_ind(n,LIS_localPet+1):&         
          LIS_ewe_halo_ind(n,LIS_localPet+1), &
          LIS_nss_halo_ind(n,LIS_localPet+1): &
          LIS_nse_halo_ind(n,LIS_localPet+1))

      do r=1,LIS_rc%lnr(n)
         do c=1,LIS_rc%lnc(n)
            if(LIS_domain(n)%gindex(c,r).ne.-1) then 
             if ( elev_subset(c,r) .ne. undef ) then
               LIS_forc(n,findex)%modelelev(LIS_domain(n)%gindex(c,r)) =&
                    elev_subset(c,r)
             endif
            endif
         enddo
      enddo
   endif

   write(LIS_logunit,*) & 
        "Finish reading original WFDE5 elevation data"
#endif

end subroutine read_wfde5_elev
 
