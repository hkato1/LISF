!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
module LIS_LMSFMod
!BOP
!
! !MODULE: LIS_LMSFMod
!
! !DESCRIPTION:
!  The code in this file implements routines to read scaling factors.
! 
!  \subsubsection{Overview}
!  This module provides routines for reading scaling factors. 
!
! !REVISION HISTORY:
!
!  30 April 2015: Augusto Getirana; Initial implementation
!
  use LIS_fileIOMod

  implicit none

  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  public :: LIS_LMSF_init          ! initializes data structures and memory
  public :: read_scaling_factor
  public :: read_first_scaling_factor
  public :: read_next_scaling_factor
  public :: LIS_LMSF_apply        
  public :: LIS_LMSF_finalize     
!------------------------------------------------------------------------------
! !PUBLIC TYPES:
!------------------------------------------------------------------------------
  public :: LIS_LMSF  ! derived datatype that stores the scaling factor data
!EOP
  type, public :: lmsf_type_dec
     real, allocatable :: tmp(:)
     real, allocatable :: q2(:)
     real, allocatable :: swd(:)
     real, allocatable :: lwd(:)
     real, allocatable :: uwind(:)
     real, allocatable :: vwind(:)
     real, allocatable :: psurf(:)
     real, allocatable :: pcp(:)
     real, allocatable :: cpcp(:)
     
  end type lmsf_type_dec
  
  type(lmsf_type_dec), allocatable :: LIS_LMSF(:)
contains

!BOP
! 
! !ROUTINE: LIS_LMSF_init
! \label{LIS_LMSF_init}
! 
! !INTERFACE:
  subroutine LIS_LMSF_init()
! !USES:
    use ESMF
    use LIS_coreMod,  only : LIS_rc, LIS_config
    use LIS_logMod
    use LIS_fileIOMod, only : LIS_readDomainConfigSpecs
  use LIS_FORC_AttributesMod 
    use LIS_timeMgrMod
! 
! !DESCRIPTION:
!
! Reads the configurable options related to scaling factors 
! and allocates memory for data structures for reading datasets
!
!EOP
    implicit none
    integer :: n, i
    integer :: rc

    allocate(LIS_LMSF(LIS_rc%nnest))
    do n=1,LIS_rc%nnest
      if(LIS_rc%scalingfactorfile(n).ne."none") then
       if (LIS_rc%scalingfactorType(n).eq."aclimo") then
        call read_scaling_factor(n)
       else if(LIS_rc%scalingfactorType(n).eq."monthly") then
        call read_first_scaling_factor(n)
       else
        write(LIS_logunit,*) 'Scaling factor Type not supported. '
        write(LIS_logunit,*) 'program stopping ...'
        call LIS_endrun
       endif
      write(LIS_logunit,*) 'Finished reading scaling factor maps'
      endif
    enddo
  end subroutine LIS_LMSF_init
!BOP
!
! !ROUTINE: read_scaling_factor
!  \label{read_scaling_factor}
!
! !REVISION HISTORY:
!  30 April 2015: Augusto Getirana; Initial Specification
!
! !INTERFACE:
subroutine read_scaling_factor(n)
! !USES:
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
  use LIS_coreMod
  use LIS_logMod,         only : LIS_logunit, &
       LIS_endrun, LIS_verify
  use LIS_FORC_AttributesMod  
  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n

! !DESCRIPTION:
!  This subroutine reads the scaling_factor data
!
!  The arguments are:
!  \begin{description}
!   \item[n]
!    index of n
!   \item[locallc]
!    landlc for the region of interest
!   \end{description}
!
!EOP      

  integer          :: ios,nid,scfId, lcid,ncId, nrId
  integer          :: nc,nr
  integer          :: r,c,i,j
  real, allocatable    :: gscale(:,:),lscale(:,:)
  logical          :: file_exists
  integer          :: sftairid,sfqairid,sfswdownid,sflwdownid,sfwindid,&
                      sfpsurfid,sfrainfid

#if (defined USE_NETCDF3 || defined USE_NETCDF4)

  inquire(file=LIS_rc%scalingfactorfile(n), exist=file_exists)
  if(file_exists) then 

     write(LIS_logunit,*)'Reading scaling factor maps...'
     write(LIS_logunit,*)trim(LIS_rc%scalingfactorfile(n))
     ios = nf90_open(path=LIS_rc%scalingfactorfile(n),&
          mode=NF90_NOWRITE,ncid=nid)
     call LIS_verify(ios,'Error in nf90_open in read_scaling_factor')
     
     ios = nf90_inq_dimid(nid,"east_west",ncId)
     call LIS_verify(ios,'Error in nf90_inq_dimid in read_scaling_factor')

     ios = nf90_inq_dimid(nid,"north_south",nrId)
     call LIS_verify(ios,'Error in nf90_inq_dimid in read_scaling_factor')

     ios = nf90_inquire_dimension(nid,ncId, len=nc)
     call LIS_verify(ios,'Error in nf90_inquire_dimension in read_scaling_factor')

     ios = nf90_inquire_dimension(nid,nrId, len=nr)
     call LIS_verify(ios,'Error in nf90_inquire_dimension in read_scaling_factor')
     
     allocate(LIS_LMSF(n)%tmp(LIS_rc%ntiles(n)))
     allocate(LIS_LMSF(n)%q2(LIS_rc%ntiles(n)))
     allocate(LIS_LMSF(n)%swd(LIS_rc%ntiles(n)))
     allocate(LIS_LMSF(n)%lwd(LIS_rc%ntiles(n)))
     allocate(LIS_LMSF(n)%uwind(LIS_rc%ntiles(n)))
     allocate(LIS_LMSF(n)%vwind(LIS_rc%ntiles(n)))
     allocate(LIS_LMSF(n)%psurf(LIS_rc%ntiles(n)))
     allocate(LIS_LMSF(n)%pcp(LIS_rc%ntiles(n)))
     !ag - LIS_FORC_CRainf%selectOpt is not yet defined 
     !if(LIS_FORC_CRainf%selectOpt.eq.1)then
       allocate(LIS_LMSF(n)%cpcp(LIS_rc%ntiles(n)))
     !end if

     allocate(gscale(LIS_rc%gnc(n),LIS_rc%gnr(n)))
     allocate(lscale(LIS_rc%lnc(n),LIS_rc%lnr(n)))


     ios = nf90_inq_varid(nid,'SF_Tair',sftairid)
     call LIS_verify(ios,'SF_Tair field not found in the LIS LSMF file')
     ios = nf90_get_var(nid,sftairid,gscale)
     call LIS_verify(ios,'Error in nf90_get_var in read_scaling_factor')

     lscale(:,:) = gscale(&
       LIS_ews_halo_ind(n,LIS_localPet+1):&         
       LIS_ewe_halo_ind(n,LIS_localPet+1), &
       LIS_nss_halo_ind(n,LIS_localPet+1): &
       LIS_nse_halo_ind(n,LIS_localPet+1))
       
     !Transferring current data to 1-D array 
     do i=1,LIS_rc%ntiles(n)
       r = LIS_domain(n)%tile(i)%row   
       c = LIS_domain(n)%tile(i)%col
       LIS_LMSF(n)%tmp(i) = lscale(c,r)
     enddo	
  
     ios = nf90_inq_varid(nid,'SF_Qair',sfqairid)
     call LIS_verify(ios,'SF_Qair field not found in the LIS LSMF file')
     ios = nf90_get_var(nid,sfqairid,gscale)
     call LIS_verify(ios,'Error in nf90_get_att in read_scaling_factor')

     lscale(:,:) = gscale(&
       LIS_ews_halo_ind(n,LIS_localPet+1):&         
       LIS_ewe_halo_ind(n,LIS_localPet+1), &
       LIS_nss_halo_ind(n,LIS_localPet+1): &
       LIS_nse_halo_ind(n,LIS_localPet+1))
       
     !Transferring current data to 1-D array 
     do i=1,LIS_rc%ntiles(n)
       r = LIS_domain(n)%tile(i)%row   
       c = LIS_domain(n)%tile(i)%col
       LIS_LMSF(n)%q2(i) = lscale(c,r)
     enddo	

     ios = nf90_inq_varid(nid,'SF_SWdown',sfswdownid)
     call LIS_verify(ios,'SF_SWdown field not found in the LIS LSMF file')
     ios = nf90_get_var(nid,sfswdownid,gscale)
     call LIS_verify(ios,'Error in nf90_get_att in read_scaling_factor')

     lscale(:,:) = gscale(&
       LIS_ews_halo_ind(n,LIS_localPet+1):&         
       LIS_ewe_halo_ind(n,LIS_localPet+1), &
       LIS_nss_halo_ind(n,LIS_localPet+1): &
       LIS_nse_halo_ind(n,LIS_localPet+1))
       
     !Transferring current data to 1-D array 
     do i=1,LIS_rc%ntiles(n)
       r = LIS_domain(n)%tile(i)%row   
       c = LIS_domain(n)%tile(i)%col
       LIS_LMSF(n)%swd(i) = lscale(c,r)
     enddo	

     ios = nf90_inq_varid(nid,'SF_LWdown',sflwdownid)
     call LIS_verify(ios,'SF_LWdown field not found in the LIS LSMF file')
     ios = nf90_get_var(nid,sflwdownid,gscale)
     call LIS_verify(ios,'Error in nf90_get_att in read_scaling_factor')

     lscale(:,:) = gscale(&
       LIS_ews_halo_ind(n,LIS_localPet+1):&         
       LIS_ewe_halo_ind(n,LIS_localPet+1), &
       LIS_nss_halo_ind(n,LIS_localPet+1): &
       LIS_nse_halo_ind(n,LIS_localPet+1))
       
     !Transferring current data to 1-D array 
     do i=1,LIS_rc%ntiles(n)
       r = LIS_domain(n)%tile(i)%row   
       c = LIS_domain(n)%tile(i)%col
       LIS_LMSF(n)%lwd(i) = lscale(c,r)
     enddo	

     ios = nf90_inq_varid(nid,'SF_Wind',sfwindid)
     call LIS_verify(ios,'SF_Wind field not found in the LIS LSMF file')
     ios = nf90_get_var(nid,sfwindid,gscale)
     call LIS_verify(ios,'Error in nf90_get_att in read_scaling_factor')

     lscale(:,:) = gscale(&
       LIS_ews_halo_ind(n,LIS_localPet+1):&         
       LIS_ewe_halo_ind(n,LIS_localPet+1), &
       LIS_nss_halo_ind(n,LIS_localPet+1): &
       LIS_nse_halo_ind(n,LIS_localPet+1))
       
     !Transferring current data to 1-D array 
     do i=1,LIS_rc%ntiles(n)
       r = LIS_domain(n)%tile(i)%row   
       c = LIS_domain(n)%tile(i)%col
       LIS_LMSF(n)%uwind(i) = lscale(c,r)
       LIS_LMSF(n)%vwind(i) = lscale(c,r)
     enddo	

     ios = nf90_inq_varid(nid,'SF_PSurf',sfpsurfid)
     call LIS_verify(ios,'SF_PSurf field not found in the LIS LSMF file')
     ios = nf90_get_var(nid,sfpsurfid,gscale)
     call LIS_verify(ios,'Error in nf90_get_att in read_scaling_factor')

     lscale(:,:) = gscale(&
       LIS_ews_halo_ind(n,LIS_localPet+1):&         
       LIS_ewe_halo_ind(n,LIS_localPet+1), &
       LIS_nss_halo_ind(n,LIS_localPet+1): &
       LIS_nse_halo_ind(n,LIS_localPet+1))
       
     !Transferring current data to 1-D array 
     do i=1,LIS_rc%ntiles(n)
       r = LIS_domain(n)%tile(i)%row   
       c = LIS_domain(n)%tile(i)%col
       LIS_LMSF(n)%psurf(i) = lscale(c,r)
     enddo	

     ios = nf90_inq_varid(nid,'SF_Rainf',sfrainfid)
     call LIS_verify(ios,'SF_Rainf field not found in the LIS LSMF file')
     ios = nf90_get_var(nid,sfrainfid,gscale)
     call LIS_verify(ios,'Error in nf90_get_att in read_scaling_factor')
     
     lscale(:,:) = gscale(&
       LIS_ews_halo_ind(n,LIS_localPet+1):&         
       LIS_ewe_halo_ind(n,LIS_localPet+1), &
       LIS_nss_halo_ind(n,LIS_localPet+1): &
       LIS_nse_halo_ind(n,LIS_localPet+1))
   
     !Transferring current data to 1-D array 
     do i=1,LIS_rc%ntiles(n)
       r = LIS_domain(n)%tile(i)%row   
       c = LIS_domain(n)%tile(i)%col
       LIS_LMSF(n)%pcp(i) = lscale(c,r)
       LIS_LMSF(n)%cpcp(i) = lscale(c,r)
     enddo	

      ios = nf90_close(nid)
     call LIS_verify(ios,'Error in nf90_close in read_scaling_factor')

     deallocate(gscale)
     deallocate(lscale)
  else
     write(LIS_logunit,*) 'scaling factor map:',LIS_rc%scalingfactorfile(n), ' does not exist'
     write(LIS_logunit,*) 'program stopping ...'
     call LIS_endrun
  endif
#endif

end subroutine read_scaling_factor
!BOP
!
! !ROUTINE: read_first_scaling_factor
!  \label{read_first_scaling_factor}
!
! !REVISION HISTORY:
!  30 Mar 2016: Hiroko Beaudoing; Initial Specification
!
! !INTERFACE:
subroutine read_first_scaling_factor(n)
! !USES:
  use LIS_coreMod
  use LIS_logMod
  use LIS_FORC_AttributesMod  
  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n

! !DESCRIPTION:
!  This subroutine reads the first instance of dynamic scaling_factor 
!  data in binary format, per forcing field
!
!  The arguments are:
!  \begin{description}
!   \item[n]
!    index of n
!   \item[locallc]
!    landlc for the region of interest
!   \end{description}
!
!EOP      

  integer          :: nc,nr
  integer          :: r,c,i,j
  real, allocatable    :: gscale(:,:),lscale(:,:)
  logical          :: file_exists
  integer          :: yr, mo, ftn
  character(len=100) :: fname1,fname2,fname3,fname4,fname5,fname6,fname7
  character(len=4)   :: fyr
  character(len=2)   :: fmo

  yr = LIS_rc%yr
  mo = LIS_rc%mo

  write(unit=fyr,fmt='(i4.4)') yr
  write(unit=fmo,fmt='(i2.2)') mo

  if ( LIS_rc%scalingfactorType(n) .eq. "monthly" ) then
  fname1=trim(LIS_rc%scalingfactorfile(n))//'_'//trim(fyr)//trim(fmo)//"_tair.1gd4r"
  fname2=trim(LIS_rc%scalingfactorfile(n))//'_'//trim(fyr)//trim(fmo)//"_qair.1gd4r"
  fname3=trim(LIS_rc%scalingfactorfile(n))//'_'//trim(fyr)//trim(fmo)//"_sw.1gd4r"
  fname4=trim(LIS_rc%scalingfactorfile(n))//'_'//trim(fyr)//trim(fmo)//"_lw.1gd4r"
  fname5=trim(LIS_rc%scalingfactorfile(n))//'_'//trim(fyr)//trim(fmo)//"_wind.1gd4r"
  fname6=trim(LIS_rc%scalingfactorfile(n))//'_'//trim(fyr)//trim(fmo)//"_pres.1gd4r"
  fname7=trim(LIS_rc%scalingfactorfile(n))//'_'//trim(fyr)//trim(fmo)//"_rain.1gd4r"
  else
   write(LIS_logunit,*) 'scaling factor file type reading not yet supported'
   write(LIS_logunit,*) 'program stopping ...'
   call LIS_endrun
  endif

  allocate(LIS_LMSF(n)%tmp(LIS_rc%ntiles(n)))
  allocate(LIS_LMSF(n)%q2(LIS_rc%ntiles(n)))
  allocate(LIS_LMSF(n)%swd(LIS_rc%ntiles(n)))
  allocate(LIS_LMSF(n)%lwd(LIS_rc%ntiles(n)))
  allocate(LIS_LMSF(n)%uwind(LIS_rc%ntiles(n)))
  allocate(LIS_LMSF(n)%vwind(LIS_rc%ntiles(n)))
  allocate(LIS_LMSF(n)%psurf(LIS_rc%ntiles(n)))
  allocate(LIS_LMSF(n)%pcp(LIS_rc%ntiles(n)))
  !ag - LIS_FORC_CRainf%selectOpt is not yet defined 
  !if(LIS_FORC_CRainf%selectOpt.eq.1)then
    allocate(LIS_LMSF(n)%cpcp(LIS_rc%ntiles(n)))
  !end if

  allocate(gscale(LIS_rc%gnc(n),LIS_rc%gnr(n)))
  allocate(lscale(LIS_rc%lnc(n),LIS_rc%lnr(n)))

  inquire(file=trim(fname1), exist=file_exists)
   if(file_exists) then 

     write(LIS_logunit,*)'Reading initial SF_tair map...'
     write(LIS_logunit,*)trim(fname1)
     ftn = LIS_getNextUnitNumber()
     open(ftn,file=fname1,form='unformatted',access='direct', &
         recl=LIS_rc%gnc(n)*LIS_rc%gnr(n)*4)
       read (ftn,rec=1) gscale
     call LIS_releaseUnitNumber(ftn)

     lscale(:,:) = gscale(&
       LIS_ews_halo_ind(n,LIS_localPet+1):&         
       LIS_ewe_halo_ind(n,LIS_localPet+1), &
       LIS_nss_halo_ind(n,LIS_localPet+1): &
       LIS_nse_halo_ind(n,LIS_localPet+1))
       
     !Transferring current data to 1-D array 
     do i=1,LIS_rc%ntiles(n)
       r = LIS_domain(n)%tile(i)%row   
       c = LIS_domain(n)%tile(i)%col
       LIS_LMSF(n)%tmp(i) = lscale(c,r)
     enddo	

   else

    LIS_LMSF(n)%tmp = 1.0

   endif ! file_exist tair

  inquire(file=trim(fname2), exist=file_exists)
   if(file_exists) then 

     write(LIS_logunit,*)'Reading initial SF_qair map...'
     write(LIS_logunit,*)trim(fname2)
     ftn = LIS_getNextUnitNumber()
     open(ftn,file=fname2,form='unformatted',access='direct', &
         recl=LIS_rc%gnc(n)*LIS_rc%gnr(n)*4)
       read (ftn,rec=1) gscale
     call LIS_releaseUnitNumber(ftn)

     lscale(:,:) = gscale(&
       LIS_ews_halo_ind(n,LIS_localPet+1):&         
       LIS_ewe_halo_ind(n,LIS_localPet+1), &
       LIS_nss_halo_ind(n,LIS_localPet+1): &
       LIS_nse_halo_ind(n,LIS_localPet+1))
       
     !Transferring current data to 1-D array 
     do i=1,LIS_rc%ntiles(n)
       r = LIS_domain(n)%tile(i)%row   
       c = LIS_domain(n)%tile(i)%col
       LIS_LMSF(n)%q2(i) = lscale(c,r)
     enddo	

   else

    LIS_LMSF(n)%q2 = 1.0

   endif ! file_exist qair

  inquire(file=trim(fname3), exist=file_exists)
   if(file_exists) then 

     write(LIS_logunit,*)'Reading initial SF_swdown map...'
     write(LIS_logunit,*)trim(fname3)
     ftn = LIS_getNextUnitNumber()
     open(ftn,file=fname3,form='unformatted',access='direct', &
         recl=LIS_rc%gnc(n)*LIS_rc%gnr(n)*4)
       read (ftn,rec=1) gscale
     call LIS_releaseUnitNumber(ftn)

     lscale(:,:) = gscale(&
       LIS_ews_halo_ind(n,LIS_localPet+1):&         
       LIS_ewe_halo_ind(n,LIS_localPet+1), &
       LIS_nss_halo_ind(n,LIS_localPet+1): &
       LIS_nse_halo_ind(n,LIS_localPet+1))
       
     !Transferring current data to 1-D array 
     do i=1,LIS_rc%ntiles(n)
       r = LIS_domain(n)%tile(i)%row   
       c = LIS_domain(n)%tile(i)%col
       LIS_LMSF(n)%swd(i) = lscale(c,r)
     enddo	

   else

    LIS_LMSF(n)%swd = 1.0

   endif ! file_exist swdown

  inquire(file=trim(fname4), exist=file_exists)
   if(file_exists) then 

     write(LIS_logunit,*)'Reading initial SF_lwdown map...'
     write(LIS_logunit,*)trim(fname4)
     ftn = LIS_getNextUnitNumber()
     open(ftn,file=fname4,form='unformatted',access='direct', &
         recl=LIS_rc%gnc(n)*LIS_rc%gnr(n)*4)
       read (ftn,rec=1) gscale
     call LIS_releaseUnitNumber(ftn)

     lscale(:,:) = gscale(&
       LIS_ews_halo_ind(n,LIS_localPet+1):&         
       LIS_ewe_halo_ind(n,LIS_localPet+1), &
       LIS_nss_halo_ind(n,LIS_localPet+1): &
       LIS_nse_halo_ind(n,LIS_localPet+1))
       
     !Transferring current data to 1-D array 
     do i=1,LIS_rc%ntiles(n)
       r = LIS_domain(n)%tile(i)%row   
       c = LIS_domain(n)%tile(i)%col
       LIS_LMSF(n)%lwd(i) = lscale(c,r)
     enddo	

   else

    LIS_LMSF(n)%lwd = 1.0

   endif ! file_exist lwdown

  inquire(file=trim(fname5), exist=file_exists)
   if(file_exists) then 

     write(LIS_logunit,*)'Reading initial SF_wind map...'
     write(LIS_logunit,*)trim(fname5)
     ftn = LIS_getNextUnitNumber()
     open(ftn,file=fname5,form='unformatted',access='direct', &
         recl=LIS_rc%gnc(n)*LIS_rc%gnr(n)*4)
       read (ftn,rec=1) gscale
     call LIS_releaseUnitNumber(ftn)

     lscale(:,:) = gscale(&
       LIS_ews_halo_ind(n,LIS_localPet+1):&         
       LIS_ewe_halo_ind(n,LIS_localPet+1), &
       LIS_nss_halo_ind(n,LIS_localPet+1): &
       LIS_nse_halo_ind(n,LIS_localPet+1))
       
     !Transferring current data to 1-D array 
     do i=1,LIS_rc%ntiles(n)
       r = LIS_domain(n)%tile(i)%row   
       c = LIS_domain(n)%tile(i)%col
       LIS_LMSF(n)%uwind(i) = lscale(c,r)
       LIS_LMSF(n)%vwind(i) = lscale(c,r)
     enddo	

   else

    LIS_LMSF(n)%uwind = 1.0
    LIS_LMSF(n)%vwind = 1.0

   endif ! file_exist wind

  inquire(file=trim(fname6), exist=file_exists)
   if(file_exists) then 

     write(LIS_logunit,*)'Reading initial SF_psurf map...'
     write(LIS_logunit,*)trim(fname6)
     ftn = LIS_getNextUnitNumber()
     open(ftn,file=fname6,form='unformatted',access='direct', &
         recl=LIS_rc%gnc(n)*LIS_rc%gnr(n)*4)
       read (ftn,rec=1) gscale
     call LIS_releaseUnitNumber(ftn)

     lscale(:,:) = gscale(&
       LIS_ews_halo_ind(n,LIS_localPet+1):&         
       LIS_ewe_halo_ind(n,LIS_localPet+1), &
       LIS_nss_halo_ind(n,LIS_localPet+1): &
       LIS_nse_halo_ind(n,LIS_localPet+1))
       
     !Transferring current data to 1-D array 
     do i=1,LIS_rc%ntiles(n)
       r = LIS_domain(n)%tile(i)%row   
       c = LIS_domain(n)%tile(i)%col
       LIS_LMSF(n)%psurf(i) = lscale(c,r)
     enddo	

   else

    LIS_LMSF(n)%psurf = 1.0

   endif ! file_exist psurf

  inquire(file=trim(fname7), exist=file_exists)
   if(file_exists) then 

     write(LIS_logunit,*)'Reading initial SF_rain map...'
     write(LIS_logunit,*)trim(fname7)
     ftn = LIS_getNextUnitNumber()
     open(ftn,file=fname7,form='unformatted',access='direct', &
         recl=LIS_rc%gnc(n)*LIS_rc%gnr(n)*4)
       read (ftn,rec=1) gscale
     call LIS_releaseUnitNumber(ftn)

     lscale(:,:) = gscale(&
       LIS_ews_halo_ind(n,LIS_localPet+1):&         
       LIS_ewe_halo_ind(n,LIS_localPet+1), &
       LIS_nss_halo_ind(n,LIS_localPet+1): &
       LIS_nse_halo_ind(n,LIS_localPet+1))
       
     !Transferring current data to 1-D array 
     do i=1,LIS_rc%ntiles(n)
       r = LIS_domain(n)%tile(i)%row   
       c = LIS_domain(n)%tile(i)%col
       LIS_LMSF(n)%pcp(i) = lscale(c,r)
       LIS_LMSF(n)%cpcp(i) = lscale(c,r)
     enddo	

   else

    LIS_LMSF(n)%pcp = 1.0
    LIS_LMSF(n)%cpcp = 1.0

   endif ! file_exist rain

  deallocate(gscale)
  deallocate(lscale)

end subroutine read_first_scaling_factor
!BOP
!
! !ROUTINE: read_next_scaling_factor
!  \label{read_next_scaling_factor}
!
! !REVISION HISTORY:
!  30 Mar 2016: Hiroko Beaudoing; Initial Specification
!
! !INTERFACE:
subroutine read_next_scaling_factor(n)
! !USES:
  use LIS_coreMod
  use LIS_logMod
  use LIS_FORC_AttributesMod  
  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n

! !DESCRIPTION:
!  This subroutine reads the dynamic scaling_factor if it's update time
!  Data is in binary format, per forcing field
!
!  The arguments are:
!  \begin{description}
!   \item[n]
!    index of n
!   \item[locallc]
!    landlc for the region of interest
!   \end{description}
!
!EOP      

  integer          :: nc,nr
  integer          :: r,c,i,j
  real, allocatable    :: gscale(:,:),lscale(:,:)
  logical          :: file_exists
  integer          :: yr, mo, ftn
  character(len=100) :: fname1,fname2,fname3,fname4,fname5,fname6,fname7
  character(len=4)   :: fyr
  character(len=2)   :: fmo

  yr = LIS_rc%yr
  mo = LIS_rc%mo

  write(unit=fyr,fmt='(i4.4)') yr
  write(unit=fmo,fmt='(i2.2)') mo

  if ( LIS_rc%scalingfactorType(n) .eq. "monthly" ) then
  fname1=trim(LIS_rc%scalingfactorfile(n))//'_'//trim(fyr)//trim(fmo)//"_tair.1gd4r"
  fname2=trim(LIS_rc%scalingfactorfile(n))//'_'//trim(fyr)//trim(fmo)//"_qair.1gd4r"
  fname3=trim(LIS_rc%scalingfactorfile(n))//'_'//trim(fyr)//trim(fmo)//"_sw.1gd4r"
  fname4=trim(LIS_rc%scalingfactorfile(n))//'_'//trim(fyr)//trim(fmo)//"_lw.1gd4r"
  fname5=trim(LIS_rc%scalingfactorfile(n))//'_'//trim(fyr)//trim(fmo)//"_wind.1gd4r"
  fname6=trim(LIS_rc%scalingfactorfile(n))//'_'//trim(fyr)//trim(fmo)//"_pres.1gd4r"
  fname7=trim(LIS_rc%scalingfactorfile(n))//'_'//trim(fyr)//trim(fmo)//"_rain.1gd4r"
  else
   write(LIS_logunit,*) 'scaling factor file type reading not yet supported'
   write(LIS_logunit,*) 'program stopping ...'
   call LIS_endrun
  endif

  allocate(gscale(LIS_rc%gnc(n),LIS_rc%gnr(n)))
  allocate(lscale(LIS_rc%lnc(n),LIS_rc%lnr(n)))

  inquire(file=trim(fname1), exist=file_exists)
   if(file_exists) then 

     write(LIS_logunit,*)'Reading next SF_tair map...'
     write(LIS_logunit,*)trim(fname1)
     ftn = LIS_getNextUnitNumber()
     open(ftn,file=fname1,form='unformatted',access='direct', &
         recl=LIS_rc%gnc(n)*LIS_rc%gnr(n)*4)
       read (ftn,rec=1) gscale
     call LIS_releaseUnitNumber(ftn)

     lscale(:,:) = gscale(&
       LIS_ews_halo_ind(n,LIS_localPet+1):&         
       LIS_ewe_halo_ind(n,LIS_localPet+1), &
       LIS_nss_halo_ind(n,LIS_localPet+1): &
       LIS_nse_halo_ind(n,LIS_localPet+1))
       
     !Transferring current data to 1-D array 
     do i=1,LIS_rc%ntiles(n)
       r = LIS_domain(n)%tile(i)%row   
       c = LIS_domain(n)%tile(i)%col
       LIS_LMSF(n)%tmp(i) = lscale(c,r)
     enddo	

   else

    LIS_LMSF(n)%tmp = 1.0

   endif ! file_exist tair

  inquire(file=trim(fname2), exist=file_exists)
   if(file_exists) then 

     write(LIS_logunit,*)'Reading next SF_qair map...'
     write(LIS_logunit,*)trim(fname2)
     ftn = LIS_getNextUnitNumber()
     open(ftn,file=fname2,form='unformatted',access='direct', &
         recl=LIS_rc%gnc(n)*LIS_rc%gnr(n)*4)
       read (ftn,rec=1) gscale
     call LIS_releaseUnitNumber(ftn)

     lscale(:,:) = gscale(&
       LIS_ews_halo_ind(n,LIS_localPet+1):&         
       LIS_ewe_halo_ind(n,LIS_localPet+1), &
       LIS_nss_halo_ind(n,LIS_localPet+1): &
       LIS_nse_halo_ind(n,LIS_localPet+1))
       
     !Transferring current data to 1-D array 
     do i=1,LIS_rc%ntiles(n)
       r = LIS_domain(n)%tile(i)%row   
       c = LIS_domain(n)%tile(i)%col
       LIS_LMSF(n)%q2(i) = lscale(c,r)
     enddo	

   else

    LIS_LMSF(n)%q2 = 1.0

   endif ! file_exist qair

  inquire(file=trim(fname3), exist=file_exists)
   if(file_exists) then 

     write(LIS_logunit,*)'Reading next SF_swdown map...'
     write(LIS_logunit,*)trim(fname3)
     ftn = LIS_getNextUnitNumber()
     open(ftn,file=fname3,form='unformatted',access='direct', &
         recl=LIS_rc%gnc(n)*LIS_rc%gnr(n)*4)
       read (ftn,rec=1) gscale
     call LIS_releaseUnitNumber(ftn)

     lscale(:,:) = gscale(&
       LIS_ews_halo_ind(n,LIS_localPet+1):&         
       LIS_ewe_halo_ind(n,LIS_localPet+1), &
       LIS_nss_halo_ind(n,LIS_localPet+1): &
       LIS_nse_halo_ind(n,LIS_localPet+1))
       
     !Transferring current data to 1-D array 
     do i=1,LIS_rc%ntiles(n)
       r = LIS_domain(n)%tile(i)%row   
       c = LIS_domain(n)%tile(i)%col
       LIS_LMSF(n)%swd(i) = lscale(c,r)
     enddo	

   else

    LIS_LMSF(n)%swd = 1.0

   endif ! file_exist swdown

  inquire(file=trim(fname4), exist=file_exists)
   if(file_exists) then 

     write(LIS_logunit,*)'Reading next SF_lwdown map...'
     write(LIS_logunit,*)trim(fname4)
     ftn = LIS_getNextUnitNumber()
     open(ftn,file=fname4,form='unformatted',access='direct', &
         recl=LIS_rc%gnc(n)*LIS_rc%gnr(n)*4)
       read (ftn,rec=1) gscale
     call LIS_releaseUnitNumber(ftn)

     lscale(:,:) = gscale(&
       LIS_ews_halo_ind(n,LIS_localPet+1):&         
       LIS_ewe_halo_ind(n,LIS_localPet+1), &
       LIS_nss_halo_ind(n,LIS_localPet+1): &
       LIS_nse_halo_ind(n,LIS_localPet+1))
       
     !Transferring current data to 1-D array 
     do i=1,LIS_rc%ntiles(n)
       r = LIS_domain(n)%tile(i)%row   
       c = LIS_domain(n)%tile(i)%col
       LIS_LMSF(n)%lwd(i) = lscale(c,r)
     enddo	

   else

    LIS_LMSF(n)%lwd = 1.0

   endif ! file_exist lwdown

  inquire(file=trim(fname5), exist=file_exists)
   if(file_exists) then 

     write(LIS_logunit,*)'Reading next SF_wind map...'
     write(LIS_logunit,*)trim(fname5)
     ftn = LIS_getNextUnitNumber()
     open(ftn,file=fname5,form='unformatted',access='direct', &
         recl=LIS_rc%gnc(n)*LIS_rc%gnr(n)*4)
       read (ftn,rec=1) gscale
     call LIS_releaseUnitNumber(ftn)

     lscale(:,:) = gscale(&
       LIS_ews_halo_ind(n,LIS_localPet+1):&         
       LIS_ewe_halo_ind(n,LIS_localPet+1), &
       LIS_nss_halo_ind(n,LIS_localPet+1): &
       LIS_nse_halo_ind(n,LIS_localPet+1))
       
     !Transferring current data to 1-D array 
     do i=1,LIS_rc%ntiles(n)
       r = LIS_domain(n)%tile(i)%row   
       c = LIS_domain(n)%tile(i)%col
       LIS_LMSF(n)%uwind(i) = lscale(c,r)
       LIS_LMSF(n)%vwind(i) = lscale(c,r)
     enddo	

   else

    LIS_LMSF(n)%uwind = 1.0
    LIS_LMSF(n)%vwind = 1.0

   endif ! file_exist wind

  inquire(file=trim(fname6), exist=file_exists)
   if(file_exists) then 

     write(LIS_logunit,*)'Reading next SF_psurf map...'
     write(LIS_logunit,*)trim(fname6)
     ftn = LIS_getNextUnitNumber()
     open(ftn,file=fname6,form='unformatted',access='direct', &
         recl=LIS_rc%gnc(n)*LIS_rc%gnr(n)*4)
       read (ftn,rec=1) gscale
     call LIS_releaseUnitNumber(ftn)

     lscale(:,:) = gscale(&
       LIS_ews_halo_ind(n,LIS_localPet+1):&         
       LIS_ewe_halo_ind(n,LIS_localPet+1), &
       LIS_nss_halo_ind(n,LIS_localPet+1): &
       LIS_nse_halo_ind(n,LIS_localPet+1))
       
     !Transferring current data to 1-D array 
     do i=1,LIS_rc%ntiles(n)
       r = LIS_domain(n)%tile(i)%row   
       c = LIS_domain(n)%tile(i)%col
       LIS_LMSF(n)%psurf(i) = lscale(c,r)
     enddo	

   else

    LIS_LMSF(n)%psurf = 1.0

   endif ! file_exist psurf

  inquire(file=trim(fname7), exist=file_exists)
   if(file_exists) then 

     write(LIS_logunit,*)'Reading next SF_rain map...'
     write(LIS_logunit,*)trim(fname7)
     ftn = LIS_getNextUnitNumber()
     open(ftn,file=fname7,form='unformatted',access='direct', &
         recl=LIS_rc%gnc(n)*LIS_rc%gnr(n)*4)
       read (ftn,rec=1) gscale
     call LIS_releaseUnitNumber(ftn)

     lscale(:,:) = gscale(&
       LIS_ews_halo_ind(n,LIS_localPet+1):&         
       LIS_ewe_halo_ind(n,LIS_localPet+1), &
       LIS_nss_halo_ind(n,LIS_localPet+1): &
       LIS_nse_halo_ind(n,LIS_localPet+1))
       
     !Transferring current data to 1-D array 
     do i=1,LIS_rc%ntiles(n)
       r = LIS_domain(n)%tile(i)%row   
       c = LIS_domain(n)%tile(i)%col
       LIS_LMSF(n)%pcp(i) = lscale(c,r)
       LIS_LMSF(n)%cpcp(i) = lscale(c,r)
     enddo	

   else

    LIS_LMSF(n)%pcp = 1.0
    LIS_LMSF(n)%cpcp = 1.0

   endif ! file_exist rain

  deallocate(gscale)
  deallocate(lscale)

end subroutine read_next_scaling_factor
!BOP
!
! !ROUTINE: LIS_LMSF_apply
!  \label{LIS_LMSF_apply}
!
! !REVISION HISTORY:
!  30 April 2015: Augusto Getirana; Initial Specification
!
! !INTERFACE:
subroutine LIS_LMSF_apply(n)
! !USES:
  use ESMF
  use LIS_coreMod,         only : LIS_rc,LIS_surface,LIS_domain
  use LIS_constantsMod,    only : LIS_CONST_SOLAR
  use LIS_FORC_AttributesMod
  use LIS_metforcingMod,  only : LIS_forc, LIS_FORC_State
  use LIS_timeMgrMod,      only : LIS_time2date, LIS_isAlarmRinging
  use LIS_logMod,          only : LIS_logunit, LIS_verify, LIS_endrun
  
  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n

! !DESCRIPTION:
!  This subroutine applies the the scaling factor (bias correction) to forcings
!
!  The arguments are:
!  \begin{description}
!   \item[n]
!    index of n
!   \item[locallc]
!    landlc for the region of interest
!   \end{description}
!
!EOP      

!  integer          :: ios1
!  integer          :: ios,nid,scfId, lcid,ncId, nrId
  integer          :: nc,nr
  integer          :: t,tid
  integer          :: findex
  integer            :: status
  type(ESMF_Field)   :: tmpField,q2Field,uField,vField,swdField,lwdField
  type(ESMF_Field)   :: psurfField,pcpField,cpcpField
  real,pointer       :: tmp(:),q2(:),uwind(:),vwind(:)
  real,pointer       :: swd(:),lwd(:),psurf(:),pcp(:),cpcp(:)

  if(LIS_rc%scalingfactorfile(n).ne."none") then
!hb - 30 Mar 2016
!only pass alarm name (if interval "monthly" is added,
!then middle of month is used as in GFRAC)
    if(LIS_rc%scalingfactorType(n).eq."monthly") then
          if(LIS_rc%da.eq.1 .and. LIS_rc%hr.eq.0 .and. LIS_rc%mn.eq.0) then
           call read_next_scaling_factor(n)
          endif
    endif
    call ESMF_StateGet(LIS_FORC_State(n),LIS_FORC_Tair%varname(1),tmpField,&
         rc=status)
    call LIS_verify(status, 'LIS_LMSF_apply error: Enable Tair in the forcing variables list')
 
    call ESMF_StateGet(LIS_FORC_State(n),LIS_FORC_Qair%varname(1),q2Field,&
         rc=status)
    call LIS_verify(status, 'LIS_LMSF_apply error: Enable Qair in the forcing variables list')

    call ESMF_StateGet(LIS_FORC_State(n),LIS_FORC_SWdown%varname(1),swdField,&
         rc=status)
    call LIS_verify(status, 'LIS_LMSF_apply error: Enable SWdown in the forcing variables list')
  
    call ESMF_StateGet(LIS_FORC_State(n),LIS_FORC_LWdown%varname(1),lwdField,&
         rc=status)
    call LIS_verify(status, 'LIS_LMSF_apply error: Enable LWdown in the forcing variables list')

    call ESMF_StateGet(LIS_FORC_State(n),LIS_FORC_Wind_E%varname(1),uField,&
         rc=status)
    call LIS_verify(status, 'LIS_LMSF_apply error: Enable Wind_E in the forcing variables list')
  
    call ESMF_StateGet(LIS_FORC_State(n),LIS_FORC_Wind_N%varname(1),vField,&
         rc=status)
    call LIS_verify(status, 'LIS_LMSF_apply error: Enable Wind_N in the forcing variables list')
  
    call ESMF_StateGet(LIS_FORC_State(n),LIS_FORC_Psurf%varname(1),psurfField,&
         rc=status)
    call LIS_verify(status, 'LIS_LMSF_apply error: Enable Psurf in the forcing variables list')
  
    call ESMF_StateGet(LIS_FORC_State(n),LIS_FORC_Rainf%varname(1),pcpField,&
         rc=status)
    call LIS_verify(status, 'LIS_LMSF_apply error: Enable Rainf in the forcing variables list')
  
    if (LIS_FORC_CRainf%selectOpt.eq.1) then
       call ESMF_StateGet(LIS_FORC_State(n),LIS_FORC_CRainf%varname(1),cpcpField,&
            rc=status)
       call LIS_verify(status, 'LIS_LMSF_apply error: Enable CRainf in the forcing variables list')
    endif

    call ESMF_FieldGet(swdField,localDE=0,farrayPtr=swd,rc=status)
    call LIS_verify(status)

    call ESMF_FieldGet(lwdField,localDE=0,farrayPtr=lwd,rc=status)
    call LIS_verify(status)
         
    call ESMF_FieldGet(pcpField,localDE=0,farrayPtr=pcp,rc=status)
    call LIS_verify(status)
  
    if (LIS_FORC_CRainf%selectOpt.eq.1) then
       call ESMF_FieldGet(cpcpField,localDE=0,farrayPtr=cpcp,rc=status)
       call LIS_verify(status)
    endif

    call ESMF_FieldGet(tmpField,localDE=0,farrayPtr=tmp,rc=status)
    call LIS_verify(status)
  
    call ESMF_FieldGet(q2Field,localDE=0,farrayPtr=q2,rc=status)
    call LIS_verify(status)

    call ESMF_FieldGet(uField,localDE=0,farrayPtr=uwind,rc=status)
    call LIS_verify(status)
  
    call ESMF_FieldGet(vField,localDE=0,farrayPtr=vwind,rc=status)
    call LIS_verify(status)

    call ESMF_FieldGet(psurfField,localDE=0,farrayPtr=psurf,rc=status)
    call LIS_verify(status)

!print*,pcp(1),LIS_rc%ntiles(n)
!nc=76089
nc=10924
    do t=1,LIS_rc%ntiles(n)
!print*,t
!if(t==nc)print*,tmp(t),LIS_LMSF(n)%tmp(t)
      tmp(t)=tmp(t)*LIS_LMSF(n)%tmp(t)
!if(t==nc)print*,q2(t),LIS_LMSF(n)%q2(t)
      q2(t)=q2(t)*LIS_LMSF(n)%q2(t)
!if(t==nc)print*,swd(t),LIS_LMSF(n)%swd(t)
      swd(t)=swd(t)*LIS_LMSF(n)%swd(t)
!if(t==nc)print*,lwd(t),LIS_LMSF(n)%lwd(t)
      lwd(t)=lwd(t)*LIS_LMSF(n)%lwd(t)
!if(t==nc)print*,uwind(t),LIS_LMSF(n)%uwind(t)
      uwind(t)=uwind(t)*LIS_LMSF(n)%uwind(t)
!if(t==nc)print*,vwind(t),LIS_LMSF(n)%vwind(t)
      vwind(t)=vwind(t)*LIS_LMSF(n)%vwind(t)
!if(t==nc)print*,psurf(t),LIS_LMSF(n)%psurf(t)
      psurf(t)=psurf(t)*LIS_LMSF(n)%psurf(t)
!if(t==nc)print*,pcp(t),LIS_LMSF(n)%pcp(t)
      pcp(t)=pcp(t)*LIS_LMSF(n)%pcp(t)
!if(t==nc)print*,cpcp(t),LIS_LMSF(n)%cpcp(t)
      if (LIS_FORC_CRainf%selectOpt.eq.1) then
        cpcp(t)=cpcp(t)*LIS_LMSF(n)%cpcp(t)
      endif
    end do  
endif

end subroutine LIS_LMSF_apply


!BOP
! !ROUTINE: LIS_LMSF_finalize
! \label{LIS_LMSF_finalize}
! 
! !INTERFACE: LMSF_finalize
  subroutine LIS_LMSF_finalize
!
! !USES:
    use LIS_coreMod,  only : LIS_rc
    implicit none
! 
! !DESCRIPTION:
! 
! Deallocates objects created in this module
! 
!EOP
    integer :: n

    do n=1,LIS_rc%nnest
      if(LIS_rc%scalingfactorfile(n).ne."none") then
        deallocate(LIS_LMSF(n)%tmp)
        deallocate(LIS_LMSF(n)%q2)
        deallocate(LIS_LMSF(n)%swd)
        deallocate(LIS_LMSF(n)%lwd)
        deallocate(LIS_LMSF(n)%uwind)
        deallocate(LIS_LMSF(n)%vwind)
        deallocate(LIS_LMSF(n)%psurf)
        deallocate(LIS_LMSF(n)%pcp)
        deallocate(LIS_LMSF(n)%cpcp)
      endif
    enddo
    
  end subroutine LIS_LMSF_finalize
 
end module LIS_LMSFMod
