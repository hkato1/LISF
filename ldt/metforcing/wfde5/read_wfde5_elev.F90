!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
! !ROUTINE: read_wfde5_elev
!  \label{read_wfde5_elev}
!
! !REVISION HISTORY:
!
!  29 mar 2022: Hiroko Beaudoing; Added WFDE5 surface altitude reader
!
! !INTERFACE:
subroutine read_wfde5_elev( n, findex, wfde5elev, elevdiff )

! !USES:
  use LDT_coreMod,       only : LDT_rc, LDT_domain
  use LDT_metforcingMod, only : LDT_forc
  use wfde5_forcingMod, only : wfde5_struc
  use LDT_logMod,        only : LDT_logunit, LDT_verify, &
                                LDT_endrun
  use LDT_fileIOMod,     only : LDT_transform_paramgrid
#if (defined USE_NETCDF3 || defined USE_NETCDF4) 
  use netcdf
#endif

  implicit none

! !ARGUMENTS: 
  integer, intent(in) :: n 
  integer, intent(in) :: findex
!- Terrain height will be set to run domain:
  real, intent(inout) :: wfde5elev(LDT_rc%lnc(n),LDT_rc%lnr(n),1)
  real, intent(inout) :: elevdiff(LDT_rc%met_nc(findex), LDT_rc%met_nr(findex))

! !DESCRIPTION:
!
!  Opens, reads, and interpolates ERA5 model elevation to the LDT
!  grid. The data will be used to perform any topographical 
!  adjustments to the forcing.
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[findex]
!   index of the forcing dataset selected
!  \end{description}
! 
!  The routines invoked are: 
!   \begin{description}
!   \end{description}
!EOP
   logical :: file_exists
   integer :: ftn_const
   integer :: i,c,r,k,iret

   integer :: varId
   real    :: read_alt(wfde5_struc(n)%ncold, wfde5_struc(n)%nrold)
   real    :: undef   ! input FillValue

   ! Grid transform fields:
   integer   :: inpts, outpts
   real      :: elev1d(wfde5_struc(n)%ncold*wfde5_struc(n)%nrold)
   logical*1 :: lb(wfde5_struc(n)%ncold*wfde5_struc(n)%nrold)
   real      :: elev_regrid(LDT_rc%lnc(n)*LDT_rc%lnr(n))
   logical*1 :: lb_regrid(LDT_rc%lnc(n)*LDT_rc%lnr(n))

! _____________________________________________________________________________

#if (defined USE_NETCDF3) 
  write(LDT_logunit,*) "ERR: WFDE5 surface altitude reader requires NetCDF4"
  call LDT_endrun()
#endif

   wfde5elev = LDT_rc%udef
   elevdiff = LDT_rc%udef

   ! Check if WFDE5 grid is selected but downscaling options, 
   !  like bilinear, is incorrectly selected.
   if( wfde5_struc(n)%gridDesc(9)  == LDT_rc%gridDesc(n,9) .and. &
       wfde5_struc(n)%gridDesc(10) == LDT_rc%gridDesc(n,10).and. &
       LDT_rc%gridDesc(n,1) == 0 .and. &
       LDT_rc%met_gridtransform_parms(findex) .ne. "neighbor" ) then
      write(LDT_logunit,*) "[WARN] The WFDE5 grid was selected for the"
      write(LDT_logunit,*) "  LDT run domain; however, 'bilinear', 'budget-bilinear',"
      write(LDT_logunit,*) "  or some other unknown option was selected to spatially"
      write(LDT_logunit,*) "  downscale the grid, which will cause errors during runtime."
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun()
   endif

   inquire(file = trim(wfde5_struc(n)%wfde5alt_file), exist=file_exists)
   if(.not. file_exists) then
      write(LDT_logunit,*) "[ERR] The WFDE5 surface altitude file ",&
            trim(wfde5_struc(n)%wfde5alt_file)," is not found."
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun
   endif

! -------------------------------------------------------------------
! Open and Read-in Forcing Terrain Hght File - Bring to LIS run domain
! -------------------------------------------------------------------
   write(LDT_logunit,*) "[INFO] Reading the WFDE5 surface altitude file: ", &
        wfde5_struc(n)%wfde5alt_file

#if (defined USE_NETCDF4) 

   call LDT_verify(nf90_open(path=trim(wfde5_struc(n)%wfde5alt_file), &
            mode=NF90_NOWRITE, ncid=ftn_const), &
           'nf90_open failed in read_wfde5_elev')

   call LDT_verify(nf90_inq_varid(ftn_const,'ASurf',varId), &
           'nf90_inq_varid failed for elev in read_wfde5_elev')

   ! Reading in PHIS field:
   call LDT_verify(nf90_get_var(ftn_const, varId, read_alt), &
           'nf90_get_var failed for elev in read_wfde5_elev')
   call LDT_verify(nf90_get_att(ftn_const, varId,'_FillValue', undef), &
           'nf90_get_att failed for elev:_FillValue in read_wfde5_elev')


   ! Initialize arrays for grid transformation:
   inpts  = LDT_rc%met_nc(findex)*LDT_rc%met_nr(findex)
   outpts = LDT_rc%lnr(n)*LDT_rc%lnc(n)
   lb     = .true.
   lb_regrid = .true.
   elev1d = -9999.0

   ! Convert 2D to 1D array for the interplation call:
   do r = 1, LDT_rc%met_nr(findex)
      do c = 1, LDT_rc%met_nc(findex)
         k= c+(r-1)*LDT_rc%met_nc(findex)
         elev1d(k) = read_alt(c,r)
         if ( elev1d(k) == undef ) then
           elev1d(k)  = LDT_rc%udef
           lb(k) = .false.
         endif
      enddo
   enddo

   ! Interp elevation field to output field:
   call LDT_transform_paramgrid(n, LDT_rc%met_gridtransform_parms(findex), &
            LDT_rc%met_gridDesc(findex,:), inpts, 1, elev1d, lb, &
            outpts, elev_regrid, lb_regrid )

   i = 0
   do r = 1, LDT_rc%lnr(n)
      do c = 1, LDT_rc%lnc(n)
         i = i + 1
         wfde5elev(c,r,1) = elev_regrid(i) 
      end do
   end do


#endif
end subroutine read_wfde5_elev
