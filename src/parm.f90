! ====================================================================
!     Subroutine module   
!     
!     Purpose: store the variables that will be used in the model 
!              frequently.
! ==================================================================== 
   
    MODULE parm
    IMPLICIT NONE

    INTEGER, PUBLIC, PARAMETER :: KR = SELECTED_REAL_KIND(r=20, p=20) !Accuracy of Real Number.
    INTEGER, PUBLIC, PARAMETER :: KI = SELECTED_INT_KIND(4)           !Accuracy of Int Number.
    INTEGER, PUBLIC, PARAMETER :: MMPLD=50          ! Max print time and boundary condition time.
    INTEGER, PUBLIC, PARAMETER :: NMATD=8           ! Max numbers of soil materials (Vertical)
    INTEGER, PUBLIC, PARAMETER :: NREGD=8           ! Max numbers of sub-regions.
    INTEGER, PUBLIC, PARAMETER :: NPARD=10          ! Max numbers of hydraulic parameters.
    INTEGER, PUBLIC, PARAMETER :: NOBSD=10          ! Max numbers of observation nodes.
    INTEGER, PUBLIC, PARAMETER :: NLAYERD=130       ! Max numbers of soil layers.
    INTEGER, PUBLIC, PARAMETER :: numc=1            ! Max numbers of plant species.

    LOGICAL (KIND=4) :: lCheck         ! If check the input data and write the check.out file.
    LOGICAL (KIND=4) :: lWat           ! True, transient water flow; False, steady water flow.
    LOGICAL (KIND=4) :: lChem          ! If calculate the solute transport.
    LOGICAL (KIND=4) :: lTemp          ! If calculate the temperature.
    LOGICAL (KIND=4) :: lCrop          ! If calculate the crop root uptake.
    LOGICAL (KIND=4) :: IfPM           ! True, calculate the ET0 with P-M equation; False, direct input.

    INTEGER (KIND=KI) :: Terr          ! Error check index.
    INTEGER (KIND=KI) :: Nlayer        ! Numbers of real soil layers.
    INTEGER (KIND=KI) :: Wlayer
    INTEGER (KIND=KI) :: NObs          ! Numbers of the observation points.
    INTEGER (KIND=KI) :: Nmat          ! Numbers of soil materials.
    INTEGER (KIND=KI) :: NReg          ! Numbers of the sub-region, for mass balance calculate.
    INTEGER (KIND=KI) :: MPL           ! Numbers of print times.
    INTEGER (KIND=KI) :: Bup           ! Upper boundary condition. 1, atmospheric boundary, open 01.wea;
                                       ! 0, given flux (irrigation), open met.in;
                                       ! <0, precipitation and irrigation, open both files.
    INTEGER (KIND=KI) :: Bdn           ! Lower boundary condition. 0, 0 flux boundary; -1, given water content,
                                       ! open met.in; -2, given groundwater depth, open met in.
    INTEGER (KIND=KI) :: Npar          ! Numbers of soil hydraulic parameters.
    INTEGER (KIND=KI) :: Drng          ! Which Drainage Function 1/2/3/4/5/6/7.
    INTEGER (kind=KI) :: Dfit          ! The Empirical Formula of Diffusion -3/-2/-1/1/2/3.
    INTEGER (KIND=KI) :: ddn           ! Diffision process related.
    INTEGER (KIND=KI) :: Nup           ! Numbers of flux upper boundary node.
    INTEGER (KIND=KI) :: Ndn           ! Numbers of flux lower boundary node.
    INTEGER (KIND=KI) :: MaxAL         ! Numbers of atmospheric data-records.
    INTEGER (KIND=KI) :: Plevel        ! Print time counter.
    INTEGER (KIND=KI) :: Tlevel        ! Calculate time counter.
    INTEGER (KIND=KI) :: Alevel        ! Atmosphere time counter.
    INTEGER (KIND=KI), ALLOCATABLE :: Obs(:)    ! The observation points.
    INTEGER (KIND=KI), ALLOCATABLE :: MATuz(:)  ! The material serial number of each layer.
    INTEGER (KIND=KI), ALLOCATABLE :: REGuz(:)  ! The sub-region serial numer of each layer.

    CHARACTER (LEN=8) :: date          ! The date of the calculation.
    CHARACTER (LEN=5) :: LUnit         ! Length Unit. 
    CHARACTER (LEN=5) :: TUnit         ! Time Unit.
    CHARACTER (LEN=5) :: MUnit         ! Mass Unit.
        
    REAL (KIND=KR) :: Tol=1E-10_KR     ! The tolerance while calculation.
    REAL (KIND=KR) :: CosAlf           ! The Cosine value of the angle.
    REAL (KIND=KR) :: dt               ! Time step.
    REAL (KIND=KR) :: dtOld            ! The initial input time step.
    REAL (KIND=KR) :: t                ! Time.
    REAL (KIND=KR) :: tinit            ! The initial time.
    REAL (KIND=KR) :: tEnd             ! The endding time.
    REAL (KIND=KR) :: tAtm             ! The time for the atmospheric condition.
    REAL (KIND=KR) :: qair             ! Precipitation.
    REAL (KIND=KR) :: qirr             ! Irrigation.
    REAL (KIND=KR) :: TotalPerco       ! Percolation flux.
    REAL (KIND=KR) :: qairt            ! Sum of Precipitation.
    REAL (KIND=KR) :: qirrt            ! Sum of Irrigation.
    REAL (KIND=KR) :: qbtmt            ! Sum of Bottom Flux.
    REAL (KIND=KR) :: CumE             ! Sum of Evaporation.
    REAL (KIND=KR) :: CumT             ! Sum of Transpiration.
    REAL (KIND=KR) :: sinkt            ! Sum of source/sink.
    REAL (KIND=KR) :: voluz            ! The initial water volumn.
    REAL (KIND=KR) :: voluz1           ! The new water volumn.
    REAL (KIND=KR) :: Epa	           ! Actual Evaporation.
    REAL (KIND=KR) :: Tra	           ! Actual Transpiration.
    REAL (KIND=KR) :: ptab
    REAL (KIND=KR) :: xConv            ! Length conversion coefficient.
    REAL (KIND=KR) :: tConv            ! Time conversion coefficient.
    REAL (KIND=KR) :: mConv            ! Mass conversion coefficient.
    REAL (KIND=KR) :: rTop             ! If lWat=false, upper flux
    REAL (KIND=KR) :: rBot             ! If lWat=false, lower flux
    REAL (KIND=KR) :: rET            ! If lWat=false, potential evapotranspiration
    
    REAL (KIND=KR), ALLOCATABLE :: TPrint(:)        ! The print time.
    REAL (KIND=KR), ALLOCATABLE :: ths(:)           ! ths for each materials.
    REAL (KIND=KR), ALLOCATABLE :: thF(:)           ! Field capacity.
    REAL (KIND=KR), ALLOCATABLE :: thw(:)           ! Wilting point.
    REAL (KIND=KR), ALLOCATABLE :: up(:,:)          ! Met.in file, Upper input.
    REAL (KIND=KR), ALLOCATABLE :: dn(:,:)          ! Met.in file, Lower input.
    REAL (KIND=KR), ALLOCATABLE :: sp(:,:)          ! The Evaporation cumulative distribution function.
    REAL (KIND=KR), ALLOCATABLE :: Par(:,:)         ! The Hydraulic parameters.
    REAL (KINd=KR), ALLOCATABLE :: evatra(:,:),precip(:,:) ! The input value.
    ! Nlayer dimension array
    REAL (KIND=KR), ALLOCATABLE :: dz(:)            ! Thickness of Each Layer.
    REAL (KIND=KR), ALLOCATABLE :: zx(:)            ! The z Coordinates of Every Node Descend from Up to Down.
    REAL (KIND=KR), ALLOCATABLE :: th(:)            ! Soil Water Content.
    REAL (KIND=KR), ALLOCATABLE :: SInk1d(:)        ! The Actural Evapotranspiration for Each Layer.
    REAL (KIND=KR), ALLOCATABLE :: Epi(:)           ! Evaporation for Each Layer.
    REAL (KIND=KR), ALLOCATABLE :: Tri(:)           ! Transpiration for Each Layer.
    REAL (KIND=KR), ALLOCATABLE :: Slope(:)         ! Empirical Formula of Diffusion.
    REAL (KIND=KR), ALLOCATABLE :: Intercept(:)     ! Empirical Formula of Diffusion.
    REAL (KIND=KR), ALLOCATABLE :: par_n(:)         ! Empirical Formula of n.

!   Solute transport
    LOGICAL (KIND=4) :: MIM             ! If mobile-immobile.
    LOGICAL (KIND=4) :: DiCr            ! If Dissolution and Crystallization process.
    
    REAL (KIND=KR) :: ct
    REAL (KIND=KR) :: mvoluz            ! The mobile water volumn.
    REAL (KIND=KR) :: imvoluz           ! The immobile water volumn.
    REAL (KIND=KR) :: gdep              ! Groundwater level (L).
    REAL (KIND=KR) :: dth               ! Given th at lower boundary.
    
    REAL (KIND=KR), ALLOCATABLE :: qflux(:,:) ! The water flux of each layer.
    REAL (KIND=KR), ALLOCATABLE :: thini(:)   ! The initial soil water content.
    REAL (KIND=KR), ALLOCATABLE :: tht(:)     ! 
    REAL (KIND=KR), ALLOCATABLE :: thm(:)     ! Mobile soil water content.
    REAL (KIND=KR), ALLOCATABLE :: thim(:)    ! Immobile Soil Water Content.
    REAL (KIND=KR), ALLOCATABLE :: thmini(:)  ! The initial mobile soil water content.
    REAL (KIND=KR), ALLOCATABLE :: thimini(:) ! The initial immobile soil water content.
    

    INTEGER (KIND=KI) :: CBup      ! Solute up boundary condition.
    INTEGER (KIND=KI) :: CBdn      ! Solute lower boundary condition.
    INTEGER (KIND=KI) :: CNup      ! Numbers of Solute Upper Boundary Node.
    INTEGER (KIND=KI) :: CNdn      ! Numbers of Solute Lower Boundary Node.
    INTEGER (KIND=KI) :: Cind      ! The factor of solute calculation time step.
    INTEGER (KIND=KI) :: index
    
    REAL (KIND=KR) :: w            ! The weight in solute convection equation.
    REAL (KIND=KR) :: Concup       ! The concentration of the upper boundary.
    REAL (KIND=KR) :: Concdn       ! The Concentration of the lower boundary.
    REAL (KIND=KR) :: Svoluz       ! The solute volumn.
    REAL (KIND=KR) :: Smvoluz      ! The mobile region solute volumn.
    REAL (KIND=KR) :: Simvoluz     ! The immobile region solute volumn.
    REAL (KIND=KR) :: Sup          ! The solution comes from upper boundary.
    REAL (KIND=KR) :: Sdn          ! The solution comes from lower boundary.
    REAL (KIND=KR) :: SSinkt       ! The solution source/sink term.
    REAL (KIND=KR) :: Csat         ! Solubility saturation

    REAL (KIND=KR), ALLOCATABLE :: rou(:)       ! density, M3/L.
    REAL (KIND=KR), ALLOCATABLE :: kx(:)        ! Adsorption Coefficient, L3/M-1.
    REAL (KIND=KR), ALLOCATABLE :: miuw(:)      ! Water First Order Reaction Rate Coefficient, T-1.
    REAL (KIND=KR), ALLOCATABLE :: mius(:)      ! Solid First Order Reaction Rate Coefficient, T-1.
    REAL (KIND=KR), ALLOCATABLE :: gamaw(:)     ! Water Zero Order Reaction Rate Coefficient, ML-3T-1.
    REAL (KIND=KR), ALLOCATABLE :: gamas(:)     ! Solid Zero Order Reaction Rate Coefficient, T-1.
    REAL (KIND=KR), ALLOCATABLE :: tao(:)       ! Dissolution Rate of Crystalline Salt, T-1.
    REAL (KIND=KR), ALLOCATABLE :: Cup(:,:)     ! Solute Upper Boundary.
    REAL (KIND=KR), ALLOCATABLE :: Cdn(:,:)     ! Solute Lower Boundary.
    REAL (KIND=KR), ALLOCATABLE :: ChPar(:,:)   ! The Hydraulic parameters.
    REAL (KIND=KR), ALLOCATABLE :: Conc(:)      ! Solute Concentration. MIM True, Mobile Solute Concentration.
    REAL (KIND=KR), ALLOCATABLE :: Concini(:)   ! The Initial Solute Concentration.
    REAL (KIND=KR), ALLOCATABLE :: Conc1(:)
    REAL (KIND=KR), ALLOCATABLE :: Conc2(:)
    REAL (KIND=KR), ALLOCATABLE :: Conc3(:)
    REAL (KIND=KR), ALLOCATABLE :: Conc4(:)
    REAL (KIND=KR), ALLOCATABLE :: Conim(:)     ! Immobile Solute Concentration.
    REAL (KIND=KR), ALLOCATABLE :: Conimini(:)  ! The Initial Immobile Solute Concentration.
    REAL (KIND=KR), ALLOCATABLE :: Ssalt(:)     ! The Initial Mobile Soil Water Content.
    REAL (KIND=KR), ALLOCATABLE :: Ssalim(:)    ! The Initial Immobile Soil Water Content.
    REAL (KIND=KR), ALLOCATABLE :: Ssink1d(:)
    REAL (KIND=KR), ALLOCATABLE :: Ssink1dm(:) 
    REAL (KIND=KR), ALLOCATABLE :: Ssink1dim(:)
    
    END MODULE parm