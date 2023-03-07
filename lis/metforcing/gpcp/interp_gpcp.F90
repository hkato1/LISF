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
! !ROUTINE: interp_gpcp
! \label{interp_gpcp}
!
! !REVISION HISTORY:
! 3 Jun 2016; Hiroko Beaudoing, adopted for GPCP 3-hourly data
!
! !INTERFACE: 
subroutine interp_gpcp(n,findex,ngpcp,f,lb,lis_gds,nc,nr, &
                       varfield)
! !USES:  
  use LIS_coreMod,     only : LIS_rc, LIS_domain
  use gpcp_forcingMod, only : gpcp_struc
  use LIS_logMod,      only : LIS_logunit, LIS_endrun

  implicit none
! !ARGUMENTS:
  integer, intent(in)    :: n 
  integer, intent(in)    :: findex
  integer                :: nc
  integer                :: nr
  integer                :: ngpcp
  real                   :: lis_gds(50)
  real                   :: f(ngpcp)
  logical*1              :: lb(ngpcp)
  real, dimension(nc,nr) :: varfield
!
! !DESCRIPTION:
!   This subroutine interpolates a given GPCP field 
!   to the LIS grid. 
!  The arguments are: 
!  \begin{description}
! \item[n]
!  index of the nest
! \item[ngpcp]
!  number of elements in the input grid
! \item[f]
!  input data array to be interpolated
! \item[lb]
!  input bitmap
! \item[lis\_gds]
!  array description of the LIS grid
! \item[nc]
!  number of columns (in the east-west dimension) in the LIS grid
! \item[nr]
!  number of rows (in the north-south dimension) in the LIS grid
! \item[varfield]
!  output interpolated field
!  \end{description} 
! 
!
!  The routines invoked are: 
!  \begin{description}
!  \item[bilinear\_interp](\ref{bilinear_interp}) \newline
!    spatially interpolate the forcing data using bilinear interpolation
!  \item[conserv\_interp](\ref{conserv_interp}) \newline
!    spatially interpolate the forcing data using conservative interpolation
!   \item[upscaleByAveraging](\ref{upscaleByAveraging}) \newline
!     upscales finer resolution forcing data to coarser resolution running
!     domain by averaging
! \end{description}
!EOP
  integer :: iret
  integer :: mo
  integer :: count1,i,j

  real, dimension(nc*nr) :: lis1d

  logical*1 :: lo(nc*nr)

!=== End variable declarations

  mo = nc*nr

!--------------------------------------------------------------------  
! Initialize output bitmap. 
!--------------------------------------------------------------------  
  lo = .true.
  gpcp_struc(n)%mi = ngpcp

  if ( gpcp_struc(n)%met_interp == "bilinear" ) then
     call bilinear_interp(lis_gds,lb,f,lo,lis1d,gpcp_struc(n)%mi,mo,&
          LIS_domain(n)%lat, LIS_domain(n)%lon,&
          gpcp_struc(n)%w111,gpcp_struc(n)%w121,&
          gpcp_struc(n)%w211,gpcp_struc(n)%w221,&
          gpcp_struc(n)%n111,gpcp_struc(n)%n121,&
          gpcp_struc(n)%n211,gpcp_struc(n)%n221,LIS_rc%udef, iret)
  elseif ( gpcp_struc(n)%met_interp == "budget-bilinear" ) then
     call conserv_interp(lis_gds,lb,f,lo,lis1d,gpcp_struc(n)%mi,mo,&
          LIS_domain(n)%lat, LIS_domain(n)%lon,gpcp_struc(n)%w112,&
          gpcp_struc(n)%w122,gpcp_struc(n)%w212,gpcp_struc(n)%w222,&
          gpcp_struc(n)%n112,gpcp_struc(n)%n122,gpcp_struc(n)%n212,&
          gpcp_struc(n)%n222,LIS_rc%udef,iret)
  elseif ( gpcp_struc(n)%met_interp == "average" ) then
    call upscaleByAveraging(gpcp_struc(n)%mi, mo, LIS_rc%udef, &
         gpcp_struc(n)%n111, lb, f, lo, lis1d)
  else
     write(LIS_logunit,*) 'The specified spatial interpolation option '
     write(LIS_logunit,*) 'is not supported for GPCP.'
     write(LIS_logunit,*) 'LIS is stopping.'
     call LIS_endrun()
  endif

!--------------------------------------------------------------------    
! Assign the interpolated precip field to the output variable
!--------------------------------------------------------------------    
  count1 = 0
  do j = 1, nr
     do i = 1, nc
        varfield(i,j) = lis1d(i+count1)
     enddo
     count1 = count1 + nc
  enddo

end subroutine interp_gpcp
