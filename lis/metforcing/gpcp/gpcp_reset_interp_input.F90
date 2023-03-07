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
! !ROUTINE: gpcp_reset_interp_input
!  \label{gpcp_reset_interp_input}
!
! !REVISION HISTORY:
!  01 Feb 2016: James Geiger: Initial specification
!  03 Jun 2016: Hiroko Beaudoing: adopted for GPCP 3-hourly data
!
! !INTERFACE:
subroutine gpcp_reset_interp_input(n, findex, gridDesci)
! !USES:

   use LIS_coreMod,        only : LIS_rc, LIS_isatAfinerResolution
   use LIS_logMod,         only : LIS_logunit, LIS_endrun
   use gpcp_forcingMod,    only : gpcp_struc

   implicit none
! !ARGUMENTS: 
   integer, intent(in) :: n
   integer, intent(in) :: findex
   real, intent(in)    :: gridDesci(50)

! !DESCRIPTION:
! Resets the neighbours and weights arrays used for spatially
! interpolating the GPCP forcing data to the LIS running domain.
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[findex]
!    forcing index
!  \item[gridDesci]
!    array of magic numbers describing the GPCP forcing domain
!  \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!   \item[bilinear\_interp\_input](\ref{bilinear_interp_input}) \newline
!    computes the neighbor, weights for bilinear interpolation
!   \item[conserv\_interp\_input](\ref{conserv_interp_input}) \newline
!    computes the neighbor, weights for conservative interpolation
!   \item[upscaleByAveraging\_input](\ref{upscaleByAveraging_input}) \newline
!    computes the neighbors for upscaling by averaging
!   \item[LIS\_isatAfinerResolution](\ref{LIS_isatAfinerResolution}) \newline
!    determines whether LIS' running domain is at a finer resolution
!    than the given resolution
!  \end{description}
!EOP

   integer :: rc

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

   if ( LIS_isatAfinerResolution(n,gridDesci(9)) ) then

      gpcp_struc(n)%met_interp = LIS_rc%met_interp(findex)

      write(LIS_logunit,*) 'MSG: The GPCP forcing resolution is coarser ' // &
                           'than the running domain.'
      write(LIS_logunit,*) '     Interpolating with the ' // &
                           trim(gpcp_struc(n)%met_interp) // ' method.'

      if ( gpcp_struc(n)%met_interp == "bilinear" ) then

         allocate(gpcp_struc(n)%n111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(gpcp_struc(n)%n121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(gpcp_struc(n)%n211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(gpcp_struc(n)%n221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(gpcp_struc(n)%w111(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(gpcp_struc(n)%w121(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(gpcp_struc(n)%w211(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
         allocate(gpcp_struc(n)%w221(LIS_rc%lnc(n)*LIS_rc%lnr(n)))

         call bilinear_interp_input(n,gridDesci,        &
                 gpcp_struc(n)%n111,gpcp_struc(n)%n121, &
                 gpcp_struc(n)%n211,gpcp_struc(n)%n221, &
                 gpcp_struc(n)%w111,gpcp_struc(n)%w121, &
                 gpcp_struc(n)%w211,gpcp_struc(n)%w221)

      elseif ( gpcp_struc(n)%met_interp == "budget-bilinear" ) then

         allocate(gpcp_struc(n)%n112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
         allocate(gpcp_struc(n)%n122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
         allocate(gpcp_struc(n)%n212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
         allocate(gpcp_struc(n)%n222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
         allocate(gpcp_struc(n)%w112(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
         allocate(gpcp_struc(n)%w122(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
         allocate(gpcp_struc(n)%w212(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))
         allocate(gpcp_struc(n)%w222(LIS_rc%lnc(n)*LIS_rc%lnr(n),25))

         call conserv_interp_input(n, gridDesci,                           &
                 gpcp_struc(n)%n112,gpcp_struc(n)%n122,gpcp_struc(n)%n212, &
                 gpcp_struc(n)%n222,gpcp_struc(n)%w112,gpcp_struc(n)%w122, &
                 gpcp_struc(n)%w212,gpcp_struc(n)%w222)
      endif
   else
      gpcp_struc(n)%met_interp = LIS_rc%met_upscale(findex)

      write(LIS_logunit,*) 'MSG: The GPCP forcing resolution is finer ' // &
                           'than the running domain.'
      write(LIS_logunit,*) '     Upscaling with the ' // &
                           trim(gpcp_struc(n)%met_interp) // ' method.'
   
   
      select case( gpcp_struc(n)%met_interp )
      case( "average" )
         if ( allocated(gpcp_struc(n)%n111) ) then
            deallocate(gpcp_struc(n)%n111)
         endif
         allocate(gpcp_struc(n)%n111(gpcp_struc(n)%mi))
   
         call upscaleByAveraging_input(gridDesci,                   &
                                       LIS_rc%gridDesc(n,:),        &
                                       gpcp_struc(n)%mi,            &
                                       LIS_rc%lnc(n)*LIS_rc%lnr(n), &
                                       gpcp_struc(n)%n111)
      case default
         write(LIS_logunit,*) 'The specified spatial interpolation option '
         write(LIS_logunit,*) 'is not supported for GPCP.'
         write(LIS_logunit,*) 'LIS is stopping.'
         call LIS_endrun()
      end select
   endif
end subroutine gpcp_reset_interp_input
