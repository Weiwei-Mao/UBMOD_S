! ====================================================================
!   Subroutine Solute_Conv
!
!   Purpose: Solute convection.
! ====================================================================                        
! =========================Incoming variables=========================
!	Conc		The initial solute concentration.
! =========================Outcoming variables========================
!	Conc1       The solute concentration after convection process.
! ====================================================================
    
    SUBROUTINE Solute_Conv
    USE parm
    IMPLICIT NONE
    INTEGER (KIND=4) :: i
    INTEGER (KIND=KI) :: m
    REAL (KIND=KR) :: Concen
    REAL (KIND=KR) :: Conca, Concb
    REAL (KIND=KR) :: Ssum0,Ssum1,delta
    REAL (KIND=KR), DIMENSION(Nlayer) :: thp,thq,thtim,thtm
    REAL (KIND=KR), DIMENSION(Nlayer,3) :: A
    REAL (KIND=KR), DIMENSION(Nlayer) :: B
    
    !thim = 1.628D-1
    !thm = 2.072D-1
    !th = 3.7D-1
    !qflux = 3.4392D0/480/0.56
    A = 0.0_KR
    B = 0.0_KR

    Concini = Conc
    
    IF (MIM) THEN
        DO i = 1,Nlayer
            m = MATuz(i)
            thtim(i) = par(2,m)*par(6,m)
            thtm(i) = tht(i) - thtim(i)
            IF (thtm(i) < 0) THEN
                thtm(i) = 0
                thtim(i) = tht(i)
            ENDIF
        ENDDO
    ENDIF

!   thp, last time step; thq, this time step.
    IF (MIM) THEN
        thp = thtm
        thq = thm
    ELSE
        thp = tht
        thq = th
    ENDIF
    
    IF (MIM) THEN
        Ssum0 = sum(Conc(1:Nlayer)*thtm(1:Nlayer)*dz(1:Nlayer))+ &
              & sum(Conim(1:Nlayer)*thtim(1:Nlayer)*dz(1:Nlayer))
    ELSE
        Ssum0 = sum(Conc(1:Nlayer)*thp(1:Nlayer)*dz(1:Nlayer))
    ENDIF
    
    DO i = 1, Nlayer
        
        IF (MIM) THEN
            IF (thp(i) == 0 .or. thq(i) == 0) THEN
                thp(i) = tht(i)
                thq(i) = th(i)
                Concen = Conc(i)
                Conc(i) = (Conim(i)*thtim(i)+Conc(i)*thtm(i))/tht(i)
            ENDIF
        ENDIF
        
!       the qflux considers dt, the equations download do not need to divide dt.
        IF (i == 1) THEN
            IF (qflux(1,2) >= 0) THEN
                A(1,2) = thq(i)+w*qflux(i,2)/dz(i)
                B(1) = thp(i)*Conc(i)+qirr*dt*Concup/dz(i)-(1-w)*qflux(i,2)*Conc(i)/dz(i)
            ELSEIF (qflux(1,2) < 0) THEN
                A(1,2) = thq(i)
                A(1,3) = w*qflux(i,2)/dz(i)
                B(1) = thp(i)*Conc(i)+qirr*dt*Concup/dz(i)-(1-w)*qflux(i,2)*Conc(i+1)/dz(i)
            ENDIF
        ENDIF
        IF (i >= 2 .and. i <= Nlayer-1) THEN
            IF (qflux(i,1) >= 0 .and. qflux(i,2) >= 0) THEN
                A(i,2) = thq(i)+w*qflux(i,2)/dz(i)
                A(i,1) = -w*qflux(i,1)/dz(i)
                B(i) = thp(i)*Conc(i)+(1-w)*qflux(i,1)*Conc(i-1)/dz(i)-(1-w)*qflux(i,2)*Conc(i)/dz(i)
            ELSEIF (qflux(i,1) >= 0 .and. qflux(i,2) < 0) THEN
                A(i,1) = -w*qflux(i,1)/dz(i)
                A(i,2) = thq(i)
                A(i,3) = w*qflux(i,2)/dz(i)
                B(i) = thp(i)*Conc(i)+(1-w)*qflux(i,1)*Conc(i-1)/dz(i)-(1-w)*qflux(i,2)*Conc(i+1)/dz(i)
            ELSEIF (qflux(i,1) < 0 .and. qflux(i,2) >= 0) THEN
                A(i,2) = thq(i)-w*qflux(i,1)/dz(i)+w*qflux(i,2)/dz(i)
                B(i) = thp(i)*Conc(i)+(1-w)*qflux(i,1)*Conc(i)/dz(i)-(1-w)*qflux(i,2)*Conc(i)/dz(i)
            ELSEIF (qflux(i,1) < 0 .and. qflux(i,2) < 0) THEN
                A(i,2) = thq(i)-w*qflux(i,1)/dz(i)
                A(i,3) = w*qflux(i,2)/dz(i)
                B(i) = thp(i)*Conc(i)+(1-w)*qflux(i,1)*Conc(i)/dz(i)-(1-w)*qflux(i,2)*Conc(i+1)/dz(i)
            ENDIF
        ENDIF
        IF (i == Nlayer) THEN
            IF (qflux(i,1) >= 0 .and. qflux(i,2) >= 0) THEN
                A(i,2) = thq(i)+w*qflux(i,2)/dz(i)
                A(i,1) = -w*qflux(i,1)/dz(i)
                B(i) = thp(i)*Conc(i)+(1-w)*qflux(i,1)*Conc(i-1)/dz(i)-(1-w)*qflux(i,2)*Conc(i)/dz(i)
            ELSEIF (qflux(i,1) >= 0 .and. qflux(i,2) < 0) THEN
                A(i,1) = -w*qflux(i,1)/dz(i)
                A(i,2) = thq(i)
                B(i) = thp(i)*Conc(i)+(1-w)*qflux(i,1)*Conc(i-1)/dz(i)-qflux(i,2)*Concdn/dz(i)
            ELSEIF (qflux(i,1) < 0 .and. qflux(i,2) >= 0) THEN
                A(i,2) = thq(i)-w*qflux(i,1)/dz(i)+w*qflux(i,2)/dz(i)
                B(i) = thp(i)*Conc(i)+(1-w)*qflux(i,1)*Conc(i)/dz(i)-(1-w)*qflux(i,2)*Conc(i)/dz(i)
            ELSEIF (qflux(i,1) < 0 .and. qflux(i,2) < 0) THEN
                A(i,2) = thq(i)-w*qflux(i,1)/dz(i)
                B(i) = thp(i)*Conc(i)+(1-w)*qflux(i,1)*Conc(i)/dz(i)-qflux(i,2)*Concdn/dz(i)
            ENDIF
        ENDIF
        
!!       Change Conc(i) back
!        IF (MIM) THEN
!            IF (thp(i) == 0 .or. thq(i) == 0) THEN
!                Conc(i) = Concen
!            ENDIF
!        ENDIF
        
    ENDDO
    
    CALL Chase(A,B,Conc1,Nlayer)
    
    IF (MIM) THEN
        DO i = 1,Nlayer
            IF (thtm(i) == 0 .and. thm(i) == 0) THEN
                Conim(i) = Conc1(i)
                Conc1(i)  = 0
            ELSEIF (thtm(i) == 0 .and. thm(i) /= 0) THEN
                Conim(i) = Conc1(i)
            ELSEIF (thtm(i) /= 0 .and. thm(i) == 0) THEN
                Conim(i) = Conc1(i)
                Conc1(i)  = 0
            ENDIF
        ENDDO
        Ssum1 = sum(Conc1(1:Nlayer)*thm(1:Nlayer)*dz(1:Nlayer))+ &
              & sum(Conim(1:Nlayer)*thim(1:Nlayer)*dz(1:Nlayer))
    ELSE
        Ssum1 = sum(Conc1(1:Nlayer)*thq(1:Nlayer)*dz(1:Nlayer))
    ENDIF
    
    !IF (t == 115) THEN
    !    WRITE(99,"(130F20.14)")Conim
    !    WRITE(99,"(130F20.14)")Conc1
    !ENDIF

!   Mass balance
    IF (qflux(Nlayer,2) >= 0) THEN
        delta = Ssum0+qirr*dt*Concup-qflux(Nlayer,2)*(w*Conc1(Nlayer)+(1-w)*Conc(Nlayer))-Ssum1
    ELSE
        delta = Ssum0+qirr*dt*Concup-qflux(Nlayer,2)*Concdn-Ssum1
    ENDIF

    IF (abs(delta) > Tol) THEN
        WRITE(*,*) 'Mass Balance Error, SUBTOUTINE Solute_Conv'
        WRITE(99,*) 'Mass Balance Error, SUBTOUTINE Solute_Conv'
        WRITE(99,"(300F10.4)") A
        PAUSE
!        STOP
    ENDIF

    END SUBROUTINE Solute_Conv


! ====================================================================
!   Subroutine Solute_Sour
!
!   Purpose: Solute source/sink term.
! ====================================================================                        
! =========================Incoming variables=========================
!	Conc1		The solute concentration after convection process.
! =========================Outcoming variables========================
!	Conc2		The solute concentration after source/sink term.
! ====================================================================
    SUBROUTINE Solute_Sour
    USE parm
    REAL (KIND=KR) :: step
    REAL (KIND=KR) :: con,con1,con2,ds
    REAL (KIND=KR) :: Solid,Solid1,Solid2
    REAL (KIND=KR) :: partao,partao1,partao2
    REAL (KIND=KR), dimension(Nlayer) :: R,Rm,Rim
    REAL (KIND=KR), dimension(Nlayer) :: F,Fm,Fim
    REAL (KIND=KR), dimension(Nlayer) :: G,Gm,Gim
    REAL (KIND=KR), dimension(Nlayer) :: Conim_new
    
    step = dt
    
!   Mobile-immobile model.
    IF (MIM) THEN
        DO i = 1,Nlayer
            m = MATuz(i)
            
!           The dissolution volumn of the crytallization salt.
            IF (DiCr) THEN
                IF (thm(i) /= 0) THEN
                    Con1 = Csat-(Csat-Conc1(i))/exp(tao(m)*step/thm(i))
                    Solid1 = Ssalt(i)-dz(i)*thm(i)*(Con1-Conc1(i))
                ELSE
                    Con1 = 0
                    Solid1 = Ssalt(i)
                ENDIF
                Con2 = Csat-(Csat-Conim(i))/exp(tao(m)*step/thim(i))
                Solid2 = Ssalim(i)-dz(i)*thim(i)*(Con2-Conim(i))
                IF (Solid1 < 0) THEN
                    Con1 = Ssalt(i)/dz(i)/thm(i)+Conc1(i)
                    partao1 = log((Csat-Conc1(i))/(Csat-Con1))*thm(i)/step
                ENDIF
                IF (Solid2 < 0) THEN
                    Con2 = Ssalim(i)/dz(i)/thim(i)+Conim(i)
                    partao2 = log((Csat-Conim(i))/(Csat-Con2))*thim(i)/step
                ENDIF
            ELSE
                partao1 = tao(m)
                partao2 = tao(m)
            ENDIF

!           The source/sink of immobile region.
            Rim(i) = 1+rou(m)*kx(m)/thim(i)
            Fim(i) = miuw(m)*thim(i)+mius(m)*rou(m)*kx(m)-tao(m)
            Gim(i) = gamaw(m)*thim(i)+gamas(m)*rou(m)+tao(m)*Csat
            
            Ssink1dim(i) = (-Fim(i)*Conim(i)+Gim(i))*step*dz(i)/Rim(i)
            
            Conim_new(i) = Conim(i)-Ssink1dim(i)/dz(i)/thim(i)

!           The source/sink of mobile region.
            IF (thm(i) == 0) THEN
                Ssink1dm(i) = 0
                Conc2(i) = 0
            ELSE
                Rm(i)  = 1+rou(m)*kx(m)/thm(i)
                Fm(i)  = miuw(m)*thm(i)+mius(m)*rou(m)*kx(m)-tao(m)
                Gm(i)  = gamaw(m)*thm(i)+gamas(m)*rou(m)+tao(m)*Csat
                
                Ssink1dm(i) = (-Fm(i)*Conc1(i)+Gm(i))*step*dz(i)/Rm(i)
                
                Conc2(i) = Conc1(i)-Ssink1dm(i)/dz(i)/thm(i)
            ENDIF

!           If the Concentration is bigger than saturation.
            IF (DiCr) THEN
                IF (Conim(i) > Csat) THEN
                    Ssalim(i) = Ssalim(i)+dz(i)*thim(i)*(Conim(i)-Csat)
                    Conim(i) = Csat
                ENDIF
                IF (Conc2(i) > Csat) THEN
                    Ssalt(i) = Ssalt(i)+dz(i)*thm(i)*(Conc2(i)-Csat)
                    Conc2(i) = Csat
                ENDIF
            ENDIF
            
        ENDDO
!   None mobile-immobile region.
    ELSE
        DO i = 1,Nlayer
            m = MATuz(i)
            
!           The dissolution volumn of the crytallization salt.
            IF (DiCr) THEN
                Con = Csat-(Csat-Conc1(i))/exp(tao(m)*step/th(i))
                ds = dz(i)*th(i)*(Con-Conc1(i))
                IF (Ssalt(i) < ds) THEN
                    Con = Ssalt(i)/dz(i)/th(i)+Conc1(i)
                    partao = log((Csat-Conc1(i))/(Csat-Con))*th(i)/step
                    ds = Ssalt(i)
                ENDIF
                Ssalt(i) = Ssalt(i)-ds
            ELSE
                ds = 0
            ENDIF

!           The dissolution/adsorption, zero/first order reaction term.
            R(i) = 1+rou(m)*kx(m)/th(i)
            F(i) = miuw(m)*th(i)+mius(m)*rou(m)*kx(m)
            G(i) = gamaw(m)*th(i)+gamas(m)*rou(m)
            
            Ssink1d(i) = (-F(i)*Conc1(i)+G(i))*step*dz(i)/R(i)

            Conc2(i) = Conc1(i)-Ssink1d(i)/dz(i)/th(i)+ds/dz(i)/th(i)

!           If the Concentration is bigger than saturation.
            IF (DiCr) THEN
                IF (Conc2(i) > Csat) THEN
                    Ssalt(i) = Ssalt(i)+dz(i)*th(i)*(Conc2(i)-Csat)
                    Conc2(i) = Csat
                ENDIF
            ENDIF
        ENDDO
    ENDIF

    END SUBROUTINE Solute_Sour

! ====================================================================
!   Subroutine Solute_Exch
!
!   Purpose: Solute exchange between immobile and mobile region.
! ====================================================================                        
! =========================Incoming variables=========================
!	Conc2		The solute concentration after source/sink process.
! =========================Outcoming variables========================
!	Conc3		The solute concentration after exchange.
! ====================================================================
    SUBROUTINE Solute_Exch
    USE parm
    REAL (KIND=KR) :: step
    REAL (KIND=KR) :: Ssum0,Ssum1,delta
    REAL (KIND=KR), dimension(Nlayer) :: Conim_new,Conc3_new
    
    Ssum0 = sum(Conim(1:Nlayer)*thim(1:Nlayer)*dz(1:Nlayer))+ &
          & sum(Conc2(1:Nlayer)*thm(1:Nlayer)*dz(1:Nlayer))
    step = dt/10
    
    DO j = 1,10
    DO i = 1,Nlayer
        IF (thm(i) == 0) THEN
            Conc3(i) = 0
            Conim_new(i) = Conim(i)
        ELSE
            m=MATuz(i)
            Conim_new(i) = ChPar(3,m)*step*(Conc2(i)-Conim(i))/thim(i)+Conim(i)
            Conc3_new(i) = Conc2(i)-ChPar(3,m)*step*(Conc2(i)-Conim(i))/thm(i)
        ENDIF
    ENDDO
    Conc2(i) = Conc3_new(i)
    Conim(i) = Conim_new(i)
    ENDDO
    
    Conim = Conim_new
    Conc3 = Conc3_new

    Ssum1 = sum(Conim_new(1:Nlayer)*thim(1:Nlayer)*dz(1:Nlayer))+ &
          & sum(Conc3(1:Nlayer)*thm(1:Nlayer)*dz(1:Nlayer))
    delta = Ssum0 - Ssum1
    IF (abs(delta) > Tol) THEN
        WRITE(*,*) 'Mass Balance Error, SUBTOUTINE Solute_Exch'
        WRITE(99,*) 'Mass Balance Error, SUBTOUTINE Solute_Exch'
        PAUSE
        STOP
    ENDIF

    END SUBROUTINE Solute_Exch   

! ====================================================================
!   Subroutine Solute_DiffD
!
!   Purpose: Solute diffusion. Use the total salinity.
! ====================================================================                        
! =========================Incoming variables=========================
!	Conc3		The solute concentration after exchange.
! =========================Outcoming variables========================
!	Conc4		The solute concentration after diffusion.
! ====================================================================
    SUBROUTINE Solute_DiffD
    USE parm
    INTEGER (KIND=4) :: m,m1,m2, dnn
    REAL (KIND=KR) :: step,deltat
    REAL (KIND=KR) :: delta
    REAL (KIND=KR), dimension(Nlayer,3) :: A
    REAL (KIND=KR), dimension(Nlayer) :: B,S,D,v
    REAL (KIND=KR), dimension(Nlayer) :: thq,Conc_ini
    REAL (KIND=KR), dimension(Nlayer,2) :: DDD
    
    A = 0.0_KR
    B = 0.0_KR
    
    step = dt
    deltat = dt
    
    dnn = 1
    
    IF (MIM) THEN
        thq = thm
    ELSE
        thq = th
    ENDIF
    
    DO i = 1,Nlayer
        IF (i == 1) THEN
            m = MATuz(i)
            m1 = MATuz(i+1)
            DDD(i,1) = 0
            DDD(i,2) = sqrt(ChPar(1,m)*ChPar(1,m1))*abs(qflux(i,2))/step+ &
                     & sqrt(ChPar(2,m)*ChPar(2,m1)*thq(1)**(10./3)*thq(2)**(10./3))/par(2,m)/par(2,m1)
        ELSEIF (i > 1 .and. i < Nlayer) THEN
            m  = MATuz(i)
            m1 = MATuz(i-1)
            m2 = MATuz(i+1)
            DDD(i,1) = sqrt(ChPar(1,m)*ChPar(1,m1))*abs(qflux(i,1))/step+ &
                     & sqrt(ChPar(2,m)*ChPar(2,m1)*thq(i)**(10./3)*thq(i-1)**(10./3))/par(2,m)/par(2,m1)
            DDD(i,2) = sqrt(ChPar(1,m)*ChPar(1,m2))*abs(qflux(i,2))/step+ &
                     & sqrt(ChPar(2,m)*ChPar(2,m2)*thq(i)**(10./3)*thq(i+1)**(10./3))/par(2,m)/par(2,m2)
        ELSEIF (i == Nlayer) THEN
            m = MATuz(i)
            m1 = MATuz(i-1)
            DDD(Nlayer,1) = sqrt(ChPar(1,m)*ChPar(1,m1))*abs(qflux(i,1))/step+ &
                          & sqrt(ChPar(2,m)*ChPar(2,m1)*thq(i)**(10./3)*thq(i-1)**(10./3))/par(2,m)/par(2,m1)
            DDD(Nlayer,2) = 0
        ENDIF
    ENDDO

!   If Mobile Water Content is Zero.
    IF (MIM) THEN
        DO i = 1,Nlayer
            IF (thm(i) == 0) THEN
                IF (i == 1) THEN
                    DDD(i,2)   = 0
                    DDD(i+1,1) = 0
                ELSEIF (i == Nlayer) THEN
                    DDD(i,1)   = 0
                    DDD(i-1,2) = 0
                ELSE
                    DDD(i,1) = 0
                    DDD(i,2) = 0
                    DDD(i-1,2) = 0
                    DDD(i+1,1) = 0
                ENDIF
            ENDIF
        ENDDO
    ENDIF
    
!   Initial Condition.
    Conc_ini = Conc3
    
    DO j = 1,dnn
        DO i = 1,Nlayer
            IF (i == 1) THEN
!               Upper Boundary Condition
                IF (MIM) THEN
                    A(1,2) = dz(1)*thm(1)/deltat+DDD(i,2)/(dz(1)+dz(2))
                    A(1,3) = -DDD(i,2)/(dz(1)+dz(2))
                    B(1) = dz(1)*Conc_ini(1)*thm(1)/deltat+DDD(i,2)*(Conc_ini(2)-Conc_ini(1))/(dz(1)+dz(2))
                ELSE
                    A(1,2) = dz(1)*th(1)/deltat+DDD(i,2)/(dz(1)+dz(2))
                    A(1,3) = -DDD(i,2)/(dz(1)+dz(2))
                    B(1) = dz(1)*Conc_ini(1)*th(1)/deltat+DDD(i,2)*(Conc_ini(2)-Conc_ini(1))/(dz(1)+dz(2))
                ENDIF
      
            ELSEIF (i>=2 .and. i<=Nlayer-1) THEN
      
                IF (MIM) THEN
                    A(i,1) = -DDD(i,1)/(dz(i-1)+dz(i))
                    A(i,2) = dz(i)*thm(i)/deltat+DDD(i,1)/(dz(i-1)+dz(i))+DDD(i,2)/(dz(i)+dz(i+1))
                    A(i,3) = -DDD(i,2)/(dz(i)+dz(i+1))
      
                    B(i) = dz(i)*thm(i)*Conc_ini(i)/deltat+DDD(i,1)*(Conc_ini(i-1)-Conc_ini(i))/(dz(i-1)+dz(i)) &
                       & + DDD(i,2)*(Conc_ini(i+1)-Conc_ini(i))/(dz(i)+dz(i+1))
                ELSE
                    A(i,1) = -DDD(i,1)/(dz(i-1)+dz(i))
                    A(i,2) = dz(i)*th(i)/deltat+DDD(i,1)/(dz(i-1)+dz(i))+DDD(i,2)/(dz(i)+dz(i+1))
                    A(i,3) = -DDD(i,2)/(dz(i)+dz(i+1))
      
                    B(i) = dz(i)*th(i)*Conc_ini(i)/deltat+DDD(i,1)*(Conc_ini(i-1)-Conc_ini(i))/(dz(i-1)+dz(i)) &
                       & + DDD(i,2)*(Conc_ini(i+1)-Conc_ini(i))/(dz(i)+dz(i+1))
                ENDIF
      
            ELSEIF (i == Nlayer) THEN
!               Lower Boundary Condition
      
                IF (MIM) THEN
                    A(Nlayer,1) = -DDD(i,1)/(dz(Nlayer)+dz(Nlayer-1))
                    A(Nlayer,2) = dz(Nlayer)*thm(Nlayer)/deltat+DDD(i,1)/(dz(Nlayer)+dz(Nlayer-1))
                    B(Nlayer) = dz(Nlayer)*Conc_ini(Nlayer)*thm(Nlayer)/deltat+DDD(i,1)*(Conc_ini(Nlayer-1)-Conc_ini(Nlayer))/(dz(Nlayer-1)+dz(Nlayer))
                ELSE
                    A(Nlayer,1) = -DDD(i,1)/(dz(Nlayer)+dz(Nlayer-1))
                    A(Nlayer,2) = dz(Nlayer)*th(Nlayer)/deltat+DDD(i,1)/(dz(Nlayer)+dz(Nlayer-1))
                    B(Nlayer) = dz(Nlayer)*Conc_ini(Nlayer)*th(Nlayer)/deltat+DDD(i,1)*(Conc_ini(Nlayer-1)-Conc_ini(Nlayer))/(dz(Nlayer-1)+dz(Nlayer))
                ENDIF
            ENDIF        
        ENDDO
      
        !   If Mobile Water Content is Zero.
        IF (MIM) THEN
            DO i = 1,Nlayer
                IF (thm(i) == 0) THEN
                    A(i,i) = 1
                ENDIF
            ENDDO
        ENDIF
            
        CALL Chase(A,B,Conc4,Nlayer)
      
        Conc_ini = Conc4
    ENDDO

    delta = sum(Conc3(1:Nlayer)*thq(1:Nlayer))-sum(Conc4(1:Nlayer)*thq(1:Nlayer))

    IF (abs(delta) > Tol) THEN
        WRITE(*,*) 'Mass Balance Error, SUBTOUTINE Solute_diff'
        WRITE(99,*) 'Mass Balance Error, SUBTOUTINE Solute_diff'
        PAUSE
        STOP
    ENDIF

    Conc = Conc4

    END SUBROUTINE Solute_DiffD

