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
module GLDAS2runoffdataMod
!BOP
! 
! !MODULE: GLDAS2runoffdataMod
! 
! !DESCRIPTION: 
!  This module contains the data structures and routines to handle 
!  runoff data from GLDAS 2.0. This implementation handles both
!  1 degree and 0.25 degree products from different LSMs.
!
! !REVISION HISTORY: 
! 8 Jan 2016: Sujay Kumar, initial implementation
! 25 Sep 2025: Hiroko Beaudoing, updated for HYMAP3
! 

! !USES: 
  use ESMF
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  
  implicit none
  
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: GLDAS2runoffdata_init
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  
  public :: GLDAS2runoffdata_struc
  
  type, public :: GLDAS2runoffdatadec
     
     real                    :: outInterval 
     character(len=LIS_CONST_PATH_LEN)            :: odir 
     character(len=LIS_CONST_PATH_LEN)            :: domfile
     
     character(len=LIS_CONST_PATH_LEN)       :: previous_filename
     logical             :: domaincheck
     real, allocatable   :: qs(:,:),qsb(:,:)

     character*20            :: model_name
     real                    :: datares
     integer                 :: nc, nr
     integer, allocatable    :: n11(:)
  end type GLDAS2runoffdatadec

  type(GLDAS2runoffdatadec), allocatable :: GLDAS2runoffdata_struc(:)

contains
 
!BOP
!
! !ROUTINE: GLDAS2runoffdata_init
! \label{GLDAS2runoffdata_init}
! 
  subroutine GLDAS2runoffdata_init

  use LIS_coreMod
  use LIS_logMod
  use LIS_timeMgrMod
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

    integer              :: n 
    integer              :: status
    integer              :: ftn
    character*100        :: lis_map_proj
    character*10         :: time
    integer              :: nc,nr
    integer              :: ncId, nrId
    real                 :: lat1,lat2
    real                 :: lon1,lon2
    real                 :: dx,dy
    real                 :: gridDesc(50)

    external :: neighbor_interp_input
    external :: upscaleByAveraging_input


    allocate(GLDAS2runoffdata_struc(LIS_rc%nnest))
       
    call ESMF_ConfigFindLabel(LIS_config,&
         "GLDAS2 runoff data output directory:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,GLDAS2runoffdata_struc(n)%odir,rc=status)
       call LIS_verify(status,&
            "GLDAS2 runoff data output directory: not defined")
       
    enddo

    call ESMF_ConfigFindLabel(LIS_config,&
         "GLDAS2 runoff data model name:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,&
            GLDAS2runoffdata_struc(n)%model_name,rc=status)
       call LIS_verify(status,&
            "GLDAS2 runoff data model name: not defined")
       
    enddo

    call ESMF_ConfigFindLabel(LIS_config,&
         "GLDAS2 runoff data spatial resolution (degree):",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,&
            GLDAS2runoffdata_struc(n)%datares,rc=status)
       call LIS_verify(status,&
            "GLDAS2 runoff data spatial resolution (degree): not defined")
       
    enddo

    call ESMF_ConfigFindLabel(LIS_config,&
         "GLDAS2 runoff data output interval:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,time,rc=status)
       call LIS_verify(status,&
            "GLDAS2 runoff data output interval: not defined")

       call LIS_parseTimeString(time,GLDAS2runoffdata_struc(n)%outInterval)
    
    enddo

    call ESMF_ConfigFindLabel(LIS_config,&
         "GLDAS2 runoff data domain file:",rc=status)
    do n=1, LIS_rc%nnest
       call ESMF_ConfigGetAttribute(LIS_config,&
            GLDAS2runoffdata_struc(n)%domfile,rc=status)
       call LIS_verify(status,&
            "GLDAS2 runoff data domain file: not defined")
       
    enddo

#if (defined USE_NETCDF3 || defined USE_NETCDF4)

    do n=1,LIS_rc%nnest
       call LIS_verify(nf90_open(path=GLDAS2runoffdata_struc(n)%domfile,&
            mode=NF90_NOWRITE,ncid=ftn),&
            'Error opening file '//trim(GLDAS2runoffdata_struc(n)%domfile))

       call LIS_verify(nf90_get_att(ftn,NF90_GLOBAL,'MAP_PROJECTION',&
            lis_map_proj),&
            'Error in nf90_get_att: MAP_PROJECTION')

       if(lis_map_proj.eq."EQUIDISTANT CYLINDRICAL") then

          call LIS_verify(nf90_inq_dimid(ftn,'east_west',ncId),&
               'Error in nf90_inq_dimid: east_west')

          call LIS_verify(nf90_inq_dimid(ftn,'north_south',nrId),&
               'Error in nf90_inq_dimid: north_south')

          call LIS_verify(nf90_inquire_dimension(ftn,ncId, &
               len=GLDAS2runoffdata_struc(n)%nc),&
               'Error in nf90_inquire_dimension: ncId')

          call LIS_verify(nf90_inquire_dimension(ftn,nrId, &
               len=GLDAS2runoffdata_struc(n)%nr),&
               'Error in nf90_inquire_dimension: nrId')

          GLDAS2runoffdata_struc(n)%domainCheck = .true.

          if(GLDAS2runoffdata_struc(n)%nc.ne.LIS_rc%lnc(n).or.&
               GLDAS2runoffdata_struc(n)%nr.ne.LIS_rc%lnr(n)) then

             GLDAS2runoffdata_struc(n)%domainCheck = .false.
             call LIS_verify(nf90_get_att(ftn,NF90_GLOBAL,'SOUTH_WEST_CORNER_LAT',&
                  lat1),&
                  'Error in nf90_get_att: SOUTH_WEST_CORNER_LAT')

             call LIS_verify(nf90_get_att(ftn,NF90_GLOBAL,'SOUTH_WEST_CORNER_LON',&
                  lon1),&
                  'Error in nf90_get_att: SOUTH_WEST_CORNER_LON')

             call LIS_verify(nf90_get_att(ftn,NF90_GLOBAL,'DX',&
                  dx),&
                  'Error in nf90_get_att: DX')

             call LIS_verify(nf90_get_att(ftn,NF90_GLOBAL,'DY',&
                  dy),&
                  'Error in nf90_get_att: DY')

             gridDesc = 0.0

             lat2 = (GLDAS2runoffdata_struc(n)%nr-1)*dx + lat1
             lon2 = (GLDAS2runoffdata_struc(n)%nc-1)*dy + lon1
          
             gridDesc(1) = 0  
             gridDesc(2) = GLDAS2runoffdata_struc(n)%nc
             gridDesc(3) = GLDAS2runoffdata_struc(n)%nr
             gridDesc(4) = lat1
             gridDesc(5) = lon1
             gridDesc(7) = lat2
             gridDesc(8) = lon2
             gridDesc(6) = 128
             gridDesc(9) = dx
             gridDesc(10) = dy
             gridDesc(20) = 64
          
             GLDAS2runoffdata_struc(n)%datares = min(dx,dy)
          
       
             if(LIS_isAtAfinerResolution(n, &
                  GLDAS2runoffdata_struc(n)%datares)) then
          
                  allocate(GLDAS2runoffdata_struc(n)%n11(LIS_rc%lnc(n)*LIS_rc%lnr(n)))
          
                  call neighbor_interp_input(n,gridDesc, &
                       GLDAS2runoffdata_struc(n)%n11)


             else
                 
                  nc = GLDAS2runoffdata_struc(n)%nc
                  nr = GLDAS2runoffdata_struc(n)%nr

                  allocate(GLDAS2runoffdata_struc(n)%n11(nc*nr))
                  
                  call upscaleByAveraging_input(gridDesc,&
                       LIS_rc%gridDesc(n,:),&
                       nc*nr,&
                       LIS_rc%lnc(n)*LIS_rc%lnr(n),&
                       GLDAS2runoffdata_struc(n)%n11)
             endif
          endif
       else
          write(LIS_logunit,*) '[ERR] currently only LIS data in lat/lon projection'
          write(LIS_logunit,*) '[ERR] is supported'
          call LIS_endrun()

       endif
    enddo
#endif

    !ag - 17Mar2016
    do n=1, LIS_rc%nnest
      GLDAS2runoffdata_struc(n)%previous_filename='none'
      allocate(GLDAS2runoffdata_struc(n)%qs(LIS_rc%lnc(n),LIS_rc%lnr(n)))
      allocate(GLDAS2runoffdata_struc(n)%qsb(LIS_rc%lnc(n),LIS_rc%lnr(n)))
      GLDAS2runoffdata_struc(n)%qs=LIS_rc%udef
      GLDAS2runoffdata_struc(n)%qsb=LIS_rc%udef
    enddo
    
  end subroutine GLDAS2runoffdata_init
end module GLDAS2runoffdataMod
