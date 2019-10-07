! ===================================================================!
!   UBMOD_S -- One dimensional mass balance model for unsaturated    !
!              water flow and solute transport. Version 1.00.        !
!                                                                    !
!   Designed by Wei Mao, Yan Zhu and Jinzhong Yang.                  !
!                                                                    !
!   Cited: Mao W, Yang J, Zhu Y, et al. An efficient soil water      !
!          balance model based on hybrid numerical and statistical   !
!          methods[J]. Journal of hydrology, 2018, 559: 721-735,     !
!          doi: 10.1016/j.jhydrol.2018.02.074.                       !
!                                                                    !
!   Feel free to contact us if you have any question.                !
!       Email: weimao@whu.edu.cn, zyan0701@163.com                   !
!                                                                    !
!                            Last modified: Sep, 2019, by Wei Mao.   !
! ===================================================================!
! ====================================================================
!     Input files:
!     1. Essential files.
!         SELECTOR.IN     The basic input information.
!         uz.IN           The discrete information.
!     2. Optional files.
!         cropdat.dat     The simple crop model.
!         01.wea          Meteorological data.
! ====================================================================
! ====================================================================
!     Output files:
!         
! ====================================================================
!   storage Routing Method~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   *************************end of documentation.	
    PROGRAM Main_Program
    USE parm
    IMPLICIT NONE
    REAL (KIND=KR) :: t1, t2, Zero
    INTEGER (KIND=4) :: lengthpath, ierr
    CHARACTER (100) :: filename, datapath

    Zero = 0.0_KR
    Terr = 0_KI
    
    filename = 'LEVEL_01.DIR'
    OPEN(10,file=filename, status='old',err=901)
    READ(10,'(a)',err=904) datapath
    CLOSE(10)
    lengthpath = Len_Trim(datapath)
!---Open the interface files. The relpath is used here.
! ====================================================================
!   input files
    OPEN(33,file=datapath(1:lengthpath)//'/Selector.in', status='old',err=902)
    OPEN(32,file=datapath(1:lengthpath)//'/Profile.in' , status='old',err=902)
!   output files
    OPEN(80,file=datapath(1:lengthpath)//'/Obsnode.out' ,status='unknown',err=903) ! Observation.
    OPEN(81,file=datapath(1:lengthpath)//'/TPrint.out'  ,status='unknown',err=903) ! Profile results.
    OPEN(89,file=datapath(1:lengthpath)//'/A_Level.out' ,status='unknown',err=903) ! Accumulated results.
    OPEN(90,file=datapath(1:lengthpath)//'/Runtime.out' ,status='unknown',err=903) ! Run time.
    OPEN(91,file=datapath(1:lengthpath)//'/Balance.out' ,status='unknown',err=903) ! Mass Balance.
    OPEN(99,file=datapath(1:lengthpath)//'/Error.out'   ,status='unknown',err=903) ! Error message.
! ====================================================================
      
!-----Begin of the program.
! ====================================================================
!     subroutine about input information.
!     call for basic information. 
    CALL Selector_In
    IF (Terr.ne.0) GOTO (905, 906, 907, 908, 920, 921, 922) Terr
!     call for node information in 1D.
    CALL Profile_In
    IF (Terr.ne.0) GOTO 909
! ====================================================================

!-----preparation of the calculation
! ====================================================================
!   CPU time.
    CALL CPU_time (t1)
!   The initial water amount in model.
    CALL Balance_Initial
    IF (Terr.ne.0) GOTO 916
!     Diffusion model.
    IF (lWat) CALL Diffusion_Model
!     Call for reference Evaportranspiration and division of E&T.
    CALL Boundary_Cond(datapath)
    IF (Terr.ne.0) GOTO (911,912,913,918,910,905) Terr
! ====================================================================

!-----Begin time loop.
100 CONTINUE

! ====================================================================
!   The steady water flow.
    IF (.NOT.lWat) THEN
        CALL Steady_Flow
    ELSE

! ====================================================================
!     Four main processes.
!     Firstly, divide infiltration and surface runoff.
!     There is no surface hydrological process by now.
          !if () then
          !    call SCS
          !elseif () then
          !    call Green-Ampt
          !else
          !    qair = Q_infiltration
          !endif
        !   Set upper boundary condition.
        CALL Set_Input

! ====================================================================
!     Secondly, advective movement driven by gravitational potential.
!       "Tipping-bucket" method.
        CALL Water_Redis
        IF (Terr.ne.0) GOTO (930, 933) Terr

! ====================================================================
!     Thirdly, source/sink term.
!     open the Files that stored E&T and the rain, and the writen ETa.
        IF(bup >= Zero) CALL Water_SetET
        IF (Terr.ne.0) GOTO (931) Terr

! ====================================================================
!     Last, Diffusive soil water movement driven by matric potential.
        CALL Water_Diff
        IF (Terr.ne.0) GOTO (932) Terr
    ENDIF
    
    IF (lChem .and. Cind > 1) THEN
        th = tht
        t = t - dt
        ct = t + dt
        Cind = Cind
        dt = dt/Cind
        index = Cind
        GOTO 100
    ENDIF

!---Solute transport module.
! ====================================================================
    IF (lChem) THEN
!   1st, the solute convection
        CALL Solute_Conv

!   2nd, the solute source/sink term
        CALL Solute_Sour

!   3nd, the exchange between mobile and immobile regions.
        IF (MIM) THEN
            CALL Solute_Exch
        ELSE
            Conc3 = Conc2
        ENDIF
    
!   4th, the diffusive term driven by concentration difference.
        CALL Solute_DiffD
    ENDIF

! ====================================================================
!     Output control.
!     Output the hydraulic head and soil moisture in 1D model.
    CALL Obs_Out
    IF (Terr.ne.0) GOTO (915) Terr
!   Call for the water balance in 1D and 3D model.
    CALL BalanceT
    IF (Terr.ne.0) GOTO 916
!   call for new time and time step.
    WRITE(*,*)"t=",sngl(t)
    
!   P-Level information
    IF (abs(TPrint(Plevel)-t) < Tol) THEN
        CALL thOut
        IF (Terr.ne.0) GOTO (917) Terr
        Plevel = Plevel + 1
    ENDIF
    
    IF(abs(t-Tend) <= Tol) THEN
        GOTO 200
    ENDIF
    
    IF (dt < dtOld) THEN
        dt = dtOld
    ENDIF
    dtOld = dt
    CALL Tcontrol
    TLevel = TLevel + 1
    t = t+dt
    
    !IF (bup == 1) THEN
    !    WRITE(150,'(A18,F10.3,2F10.6)')date,t,Epa,Tra
    !    CALL DateAdd(date,1,date)
    !    Epa=Zero
    !    Tra=Zero
    !ENDIF

    GOTO 100
    
200 CALL CPU_time (t2)
    WRITE(90,*,err=914)'Real time [sec]',t2-t1

    CLOSE(80) 
    CLOSE(81) 
    CLOSE(89)
    CLOSE(90)
    CLOSE(91)
    CLOSE(99)
    CLOSE(110) 
    CLOSE(130)
    CLOSE(150) 

    STOP

! Input or output error.
901 ierr=1
    GOTO 999
902 ierr=2
    GOTO 999
903 ierr=3
    GOTO 999
904 ierr=4
    GOTO 999
905 ierr=5
    GOTO 999
906 ierr=6
    GOTO 999
907 ierr=7
    GOTO 999
908 ierr=8
    GOTO 999
909 ierr=9
    GOTO 999
910 ierr=10
    GOTO 999
911 ierr=11
    GOTO 999
912 ierr=12
    GOTO 999
913 ierr=13
    GOTO 999
914 ierr=14
    GOTO 999
915 ierr=15
    GOTO 999
916 ierr=16
    GOTO 999
917 ierr=17
    GOTO 999
918 ierr=18
    GOTO 999
919 ierr=19
    GOTO 999
! Data check error
920 ierr=20
    GOTO 999
921 ierr=21
    GOTO 999
922 ierr=22
    GOTO 999
923 ierr=23
    GOTO 999
    
930 ierr=30
    GOTO 999
931 ierr=31
    GOTO 999
932 ierr=32
    GOTO 999
933 ierr=33
    GOTO 999

999 CALL Error_Out(ierr)
    PAUSE
    STOP
    
    END PROGRAM Main_Program
    
    SUBROUTINE Error_Out(ierr)
    IMPLICIT NONE
    INTEGER (KIND=4) :: ierr
    CHARACTER (LEN=100), DIMENSION(40) :: cErr

    cErr( 1)='Open file error in file LEVEL_01.DIR !'
    cErr( 2)='Open file error for input files !'
    cErr( 3)='Error when opening the output file !'
    cErr( 4)='Error when reading from an input file LEVEL_01.DIR !'
    cErr( 5)='Error when reading from an input file Selector.in Basic Information !'
    cErr( 6)='Error when reading from an input file Selector.in Material Information !'
    cErr( 7)='Error when reading from an input file Selector.in Time Information !'
    cErr( 8)='Error when reading from an input file Selector.in Solute Information !'
    cErr( 9)='Error when reading from an input file Profile.in !'
    cErr(10)='Error when reading from an input file Met.in !'
    cErr(11)='Error when opening or reading from an input file 01.wea !'
    cErr(12)='Error when opening or reading from an input file cropdat.dat !'
    cErr(13)='Error when opening or reading from an input file crp/et0/eti !'
    cErr(14)='Error when writing to the output file Runtime.out !'
    cErr(15)='Error when writing to the output file A_Level.out !'
    cErr(16)='Error when writing to the output file Balance.out !'
    cErr(17)='Error when writing to the output file TPrint.out !'
    cErr(18)='Error when writing to the output file eta.dat !'
    ! 
    cErr(20)='Input option is wrong !'
    cErr(21)='Dimension of the Array is exceeded !'
    cErr(22)='The soil water hydraulic parameters is wrong !'
    
    cErr(30)='Mass balance error in Water_Redis module !'
    cErr(31)='Mass balance error in Water_SetET module !'
    cErr(32)='Mass balance error in Water_Diff module !'
    cErr(32)='Too much infiltration water in Water_Redis module !'
    
    WRITE(*,*) cErr(ierr)
    WRITE(99,*) cErr(ierr)
    
    END SUBROUTINE Error_Out