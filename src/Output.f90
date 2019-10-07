! ====================================================================
!   Subroutine BalanceT
!
!   Purpose: The mass balance results, including water and solute.
! ====================================================================                        
! =========================Incoming variables=========================
!   voluz       The water volumn at the beginning.
!   voluz1      The water volumn at t time.
!   Svoluz      The solute volumn at the beginning.
!   Smvoluz     The mobile region solute volumn at the beginning.
!   Simvoluz    The immobile region solute volumn at the beginning.
!   dvol        The difference between voluz and voluz1.
!   qairt       Sum of precipitaiton.
!   qirrt       Sum of irrigation.
!   qbtmt       Sum of bottom flux.
!   sinkt       Actual evapotranspiration.
!   CumE        Sum of actual evaporation.
!   CumT        Sum of actual transpiration.
!   wbalt       Absolute water balance error.
!   wbalr       Relative water balance error.
!   sink1d(:)   Actual evaporation and transpiration at time t for
!               each layer.
! =========================Outcoming variables========================
!   None
! =========================Related files==============================
!   89(balance1d.dat)    The whole mass balance results.
! ====================================================================
    SUBROUTINE BalanceT
    USE parm
    IMPLICIT NONE

    REAL (KIND=KR) :: svoluz1,smvoluz1,simvoluz1,mvoluz1,imvoluz1
    REAL (KIND=KR) :: dvol,Sdvol,qiirt
    REAL (KIND=KR) :: wbalt,wbalr,Swbalt,Swbalr

    IF (MIM) THEN
!	    Calculate the total water storage[m]
        voluz1   = sum(th(1:Nlayer)*dz(1:Nlayer))
        mvoluz1  = sum(thm(1:Nlayer)*dz(1:Nlayer))
        imvoluz1 = sum(thim(1:Nlayer)*dz(1:Nlayer))
        IF (abs(voluz1-mvoluz1-imvoluz1) > Tol) THEN
            WRITE(*,*) 'Mass Balance Error, SUBTOUTINE BalanceT'
            WRITE(99,*) 'Mass Balance Error, SUBTOUTINE BalanceT'
            PAUSE
            STOP
        ENDIF
!	    Calculate the variation of the storage between t=(0,t)
        dvol  = voluz1-voluz
        qairt = qairt+qair*dt
        qirrt = qirrt+qirr*dt
        qbtmt = qbtmt+qflux(Nlayer,2)
        sinkt = sinkt+sum(sink1d)
        CumE  = CumE+Epa
        CumT  = CumT+Tra
        wbalt = dvol-qairt-qirrt+qbtmt+sinkt
        wbalr = abs(wbalt*1D2)/max(abs(dvol),qair+qirr+qflux(Nlayer,2)+sinkt)
        IF (wbalr > 99.9999) wbalr = 0
        
        IF (lchem) THEN
            Svoluz1   = sum(Conc(1:Nlayer)*thm(1:Nlayer)*dz(1:Nlayer)+Conim(1:Nlayer)*thim(1:Nlayer)*dz(1:Nlayer))
            Sdvol     = Svoluz1-Svoluz
            Smvoluz1  = sum(Conc(1:Nlayer)*thm(1:Nlayer)*dz(1:Nlayer))
            Simvoluz1 = sum(Conim(1:Nlayer)*thim(1:Nlayer)*dz(1:Nlayer))
            Sup     = Sup+qirr*Concup
            Ssinkt  = Ssinkt+sum(Ssink1d)
            IF (qflux(Nlayer,2) >= 0) THEN ! different direction, different concentration
                Sdn = Sdn+qflux(Nlayer,2)*Conc1(Nlayer)
            ELSE
                Sdn = Sdn+qflux(Nlayer,2)*Concdn
            ENDIF
            Swbalt  = Sdvol-Sup+Sdn+Ssinkt
            Swbalr  = abs(Swbalt*1D2)/max(abs(Sdvol),Sup+Sdn+Ssinkt)
            IF (Swbalr > 99.9999) Swbalr = 0
            WRITE(89,"(f7.3,18E20.6E3)") t,voluz1,dvol,Svoluz1,Smvoluz1,Simvoluz1,Sdvol,qairt,qirrt,Sup, &
                                     &  qbtmt,Sdn,CumE,CumT,Ssinkt,wbalt,wbalr,Swbalt,Swbalr
        ELSE
	        WRITE(89,"(f7.3,11E20.6E3)") t,voluz1,mvoluz1,imvoluz1,dvol,qairt,qiirt,qbtmt,CumE,CumT,wbalt,wbalr
        ENDIF
        
    ELSE
!	    Calculate the total water storage[m]
        voluz1  = sum(th(1:Nlayer)*dz(1:Nlayer))
!	    Calculate the variation of the storage between t=(0,t)
        dvol  = voluz1-voluz
        qairt = qairt+qair*dt
        qirrt = qirrt+qirr*dt
        qbtmt = qbtmt+qflux(Nlayer,2)
        sinkt = sinkt+sum(sink1d)
        CumE  = CumE+Epa
        CumT  = CumT+Tra
        wbalt = dvol-qairt-qirrt+qbtmt+sinkt
        wbalr = abs(wbalt*100.0_KR)/max(abs(dvol),qair+qirr+qflux(Nlayer,2)+sinkt)  !reference swms
        IF (abs(wbalr-100.0_KR)<Tol) wbalr = 0.0_KR

        IF (lchem) THEN
            Svoluz1 = sum(Conc(1:Nlayer)*th(1:Nlayer)*dz(1:Nlayer))
            Sdvol   = Svoluz1-Svoluz
            Sup     = Sup+qirr*dt*Concup
            Ssinkt  = Ssinkt+sum(Ssink1d)
            IF (qflux(Nlayer,2) >= 0) THEN
                Sdn = Sdn+qflux(Nlayer,2)*(w*Conc1(Nlayer)+(1-w)*Concini(Nlayer))
            ELSE
                Sdn = Sdn+qflux(Nlayer,2)*Concdn
            ENDIF
            Swbalt  = Sdvol-Sup+Sdn+Ssinkt
            Swbalr  = abs(Swbalt*100.0_KR)/max(abs(Sdvol),Sup+Sdn+Ssinkt)
            IF ((Swbalr-100.0_KR)<Tol) Swbalr = 0.0_KR
            WRITE(89,"(f7.3,16E20.6E3)") t,voluz1,dvol,Svoluz1,Sdvol,qairt,qirrt,Sup, &
                                     &  qbtmt,Sdn,CumE,CumT,Ssinkt,wbalt,wbalr,Swbalt,Swbalr
        ELSE
            IF (bdn == -2) THEN
                WRITE(89,"(f7.3,10E20.6E3)") t,voluz1,dvol,qairt,qirrt,qbtmt,CumE,CumT,gdep,wbalt,wbalr
            ELSE
                WRITE(89,"(f7.3,9E20.6E3)") t,voluz1,dvol,qairt,qirrt,qbtmt,CumE,CumT,wbalt,wbalr
            ENDIF
        ENDIF

    ENDIF

    RETURN
    
901 Terr=1
    RETURN
    END SUBROUTINE BalanceT

! ====================================================================
!   Subroutine Result_Out
!
!   Purpose: Print the water content and solute concentration of each
!            layer.
! ====================================================================                        
! =========================Incoming variables=========================
!   None
! =========================Outcoming variables========================
!   None
! =========================Related files==============================
!   80(thObs.dat)       All of the results.
!   81(profile.dat)     The results at print time.    
! ====================================================================
    SUBROUTINE Obs_Out
    USE parm
    IMPLICIT NONE
    INTEGER (KIND=4) :: i

    IF (MIM) THEN
        IF (lchem) THEN
            WRITE(80,"(1X,1000F10.4)") t,(th(Obs(i)),Conc(Obs(i)),Conim(Obs(i)),(thm(Obs(i))* &
                & Conc(Obs(i))+thim(Obs(i))*Conim(Obs(i)))/th(Obs(i)),i=1,NObs)
        ELSE
            WRITE(80,"(1X,1000F10.4)") t,(th(Obs(i)),thim(Obs(i)),i=1,NObs)
        ENDIF
    ELSE
        IF (lchem) THEN
            WRITE(80,"(1X,1000F10.4)") t,(th(Obs(i)),Conc(Obs(i)),i=1,NObs)
        ELSE
            WRITE(80,"(1X,1000F10.4)") t,(th(Obs(i)),i=1,NObs)
        ENDIF
    ENDIF

    !WRITE(81,*)'Zone T="', t  
    !DO j=1,Nlayer
    !    WRITE(81,*)th(j),zx(j)+dz(j)/2,dz(j)
    !ENDDO
    RETURN

901 Terr=1
    RETURN

    END SUBROUTINE Obs_Out
    
    SUBROUTINE thOut
    USE parm
    IMPLICIT NONE
    INTEGER (KIND=4) :: i,j
    
    DO i = 1,MPL
        IF ((t-TPrint(i))<Tol) THEN
            WRITE(81,"(A10,F10.6)")'Zone T=', t
            IF (MIM .and. DiCr) THEN
                WRITE(81,"()")''
            ELSEIF (MIM .and. (.NOT.DiCr)) THEN
                WRITE(81,"()")
            ELSEIF ((.NOT.MIM) .and. DiCr) THEN
                WRITE(81,"()")
            ELSEIF ((.NOT.MIM) .and. (.NOT.DiCr)) THEN
                WRITE(81,"(A80)")'       Theta                Conc               Depth                 dz         '
            ENDIF
            DO j=1,Nlayer
                IF (MIM) THEN
                    IF (DiCr) THEN
                        WRITE(81,"(8F20.14)")thim(j),Conim(j),Ssalim(j),thm(j),Conc(j),Ssalt(j),zx(j)+dz(j)/2,dz(j)
                    ELSE
                        WRITE(81,"(6F20.14)")thim(j),Conim(j),thm(j),Conc(j),zx(j)+dz(j)/2,dz(j)
                    ENDIF
                ELSE
                    IF (DiCr) THEN
                        WRITE(81,"(5F20.14)")th(j),Conc(j),Ssalt(j),zx(j)+dz(j)/2,dz(j)
                    ELSE
                        WRITE(81,"(4F20.14)")th(j),Conc(j),zx(j)+dz(j)/2,dz(j)
                    ENDIF
                ENDIF
            ENDDO
            EXIT
        ENDIF
    ENDDO
    RETURN

901 Terr=1
    RETURN
    END SUBROUTINE
