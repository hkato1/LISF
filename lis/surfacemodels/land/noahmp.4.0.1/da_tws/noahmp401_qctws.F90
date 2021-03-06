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
! !ROUTINE: noahmp401_qctws
! \label{noahmp401_qctws}
!
! !REVISION HISTORY:
! 14 Mar 2017: Sujay Kumar; Initial Specification
! 29 May 2020: Bailing Li; Created for Noah-MP4.0.1 
!
! !INTERFACE:
subroutine noahmp401_qctws(n, LSM_State)

! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use LIS_logMod,  only  : LIS_verify
  use noahmp401_lsmMod
  !use module_sf_noahmplsm_401
  use NOAHMP_TABLES_401, ONLY : SMCMAX_TABLE,SMCWLT_TABLE

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: LSM_State
!
! !DESCRIPTION:
!
!  Returns the soilmoisture related state prognostic variables for
!  data assimilation
! 
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest \newline
!  \item[LSM\_State] ESMF State container for LSM state variables \newline
!  \end{description}
!EOP
  integer                :: t
  integer                :: status
  real, pointer          :: soilm1(:)
  real, pointer          :: soilm2(:)
  real, pointer          :: soilm3(:)
  real, pointer          :: soilm4(:)
  real                   :: smmax
  real                   :: smmin
  
  type(ESMF_Field)       :: sm1Field
  type(ESMF_Field)       :: sm2Field
  type(ESMF_Field)       :: sm3Field
  type(ESMF_Field)       :: sm4Field

!Wanshu
  type(ESMF_Field)       :: gwField
  real, pointer          :: gws(:)
  real                   :: gwsmax, gwsmin
  real                   :: MIN_THRESHOLD,MAX_threshold,sm_threshold
  integer                :: SOILTYP
!------- 

  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 1",sm1Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet for Soil Moisture Layer 1 failed in noahmp401_qctws")

  call ESMF_FieldGet(sm1Field,localDE=0,farrayPtr=soilm1,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet for Soil Moisture Layer 1 failed in noahmp401_qctws")

  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 2",sm2Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet for Soil Moisture Layer 2 failed in noahmp401_qctws")

  call ESMF_FieldGet(sm1Field,localDE=0,farrayPtr=soilm2,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet for Soil Moisture Layer 2 failed in noahmp401_qctws")

  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 3",sm3Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet for Soil Moisture Layer 3 failed in noahmp401_qctws")

  call ESMF_FieldGet(sm1Field,localDE=0,farrayPtr=soilm3,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet for Soil Moisture Layer 3 failed in noahmp401_qctws")
  
  call ESMF_StateGet(LSM_State,"Soil Moisture Layer 4",sm4Field,rc=status)
  call LIS_verify(status,&
       "ESMF_StateGet for Soil Moisture Layer 4 failed in noahmp401_qctws")

  call ESMF_FieldGet(sm1Field,localDE=0,farrayPtr=soilm4,rc=status)
  call LIS_verify(status,&
       "ESMF_FieldGet for Soil Moisture Layer 4 failed in noahmp401_qctws")

  !Wanshu
  call ESMF_StateGet(LSM_State,"Groundwater Storage",gwField,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for gw in noahmp401_qctws')
  
  call ESMF_FieldGet(gwField,localDE=0,farrayPtr=gws,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for gw in noahmp401_qctws')

  call ESMF_AttributeGet(gwField,"Max Value",gwsmax,rc=status)
  call LIS_verify(status,&
       "ESMF_AttributeGet: Max Value failed in noahmp401_qctws")

  call ESMF_AttributeGet(gwField,"Min Value",gwsmin,rc=status)
  call LIS_verify(status,&
       "ESMF_AttributeGet: Min Value failed in noahmp401_qctws")
  !-------

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
  !Bailing Li: max min soil moisture should be retrieved based on soil type
     SOILTYP = NOAHMP401_struc(n)%noahmp401(t)%soiltype
     MAX_THRESHOLD = SMCMAX_TABLE(SOILTYP) 
     MIN_THRESHOLD = SMCWLT_TABLE(SOILTYP) 
     sm_threshold = MAX_THRESHOLD - 0.02


     if(soilm1(t).gt.sm_threshold) soilm1(t) = sm_threshold 
     if(soilm1(t).lt.MIN_THRESHOLD) soilm1(t) = MIN_THRESHOLD

     if(soilm2(t).gt.sm_threshold) soilm2(t) = sm_threshold
     if(soilm2(t).lt.MIN_THRESHOLD) soilm2(t) = MIN_THRESHOLD

     if(soilm3(t).gt.sm_threshold) soilm3(t) = sm_threshold
     if(soilm3(t).lt.MIN_THRESHOLD) soilm3(t) = MIN_THRESHOLD

     if(soilm4(t).gt.sm_threshold) soilm4(t) = sm_threshold
     if(soilm4(t).lt.MIN_THRESHOLD) soilm4(t) = MIN_THRESHOLD
     !Wanshu
     if(gws(t).gt.gwsmax) gws(t) = gwsmax
     if(gws(t).lt.gwsmin) gws(t) = gwsmin
     !------
  enddo

end subroutine noahmp401_qctws
