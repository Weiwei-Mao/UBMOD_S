! ====================================================================
!   SUBROUTINE BalanceI   
!     
!   PROPOSE: the initial statistics.
! ====================================================================
! =========================Incoming variables=========================
!   par(:,:)    The hydraulic parameters.
!   th(:,:)     The soil water content.
!   dz(:)       The thickness of each layer.
!   matuz(:)    The material serial number.
! =========================Outcoming variables========================
!   voluz       The initial water volumn [m3].
! =========================related files==============================
!   balance1d.dat
! =========================related functions==========================
!   None.
! ====================================================================
    SUBROUTINE Balance_Initial
    USE parm
    IMPLICIT NONE
    
    voluz  = 0.0_KR
    qairt  = 0.0_KR
    qirrt  = 0.0_KR
    qbtmt  = 0.0_KR
    CumE   = 0.0_KR
    CumT   = 0.0_KR
    sinkt  = 0.0_KR
!   Solute
    Svoluz   = 0.0_KR
	Sup      = 0.0_KR
	Sdn      = 0.0_KR
	Ssinkt   = 0.0_KR
    Smvoluz  = 0.0_KR
    Simvoluz = 0.0_KR

    IF (MIM) THEN
        voluz   = sum(thini(1:Nlayer)*dz(1:Nlayer))
        mvoluz  = sum(thmini(1:Nlayer)*dz(1:Nlayer))
        imvoluz = sum(thimini(1:Nlayer)*dz(1:Nlayer))
        
        IF (lchem) THEN
            Smvoluz  = sum(thmini(1:Nlayer)*dz(1:Nlayer)*Concini(1:Nlayer))
            Simvoluz = sum(thimini(1:Nlayer)*dz(1:Nlayer)*Conimini(1:Nlayer))
            Svoluz   = Smvoluz+Simvoluz
            WRITE(89,'(A150)')'Variables=t,mvoluz,mDvoluz,imvoluz,imDvoluz,Smvoluz,SmDvoluz, &
                             & Simvoluz,SimDvoluz,qair,Sup,qbtm,Sdn,CumE,CumT,SSink,Error, &
                             & Error%,SError,SError%'
        ELSE
            WRITE(89,'(A80)')'Variables=t,voluz,mvoluz,imvoluz,Dvoluz,qair,qirr,qbtm,CumE, &
                             & CumT,Error,Error%'
        ENDIF
    ELSE
        voluz = sum(th(1:Nlayer)*dz(1:Nlayer))
        IF (lchem) THEN
            Svoluz = sum(thini(1:Nlayer)*dz(1:Nlayer)*Concini(1:Nlayer))
            IF (bdn == -1) THEN
                Svoluz = sum(thini(1:Nlayer-1)*dz(1:Nlayer-1)*Concini(1:Nlayer-1))
            ENDIF
            WRITE(89,'(A100)')'Variables=t,voluz,Dvoluz,Svoluz,SDvoluz,qair,qirr,Sup,qbtm,Sdn, &
                             & CumE,CumT,SSink,Error,Error%,SError,SError%'
        ELSE
            WRITE(89,'(A60)')'Variables=t,voluz,Dvoluz,qair,qirr,qbtm,CumE,CumT,Error,Error%' 
        ENDIF
    ENDIF
    RETURN
    
901 Terr=1
    RETURN
    END SUBROUTINE Balance_Initial

! ====================================================================
!   SUBROUTINE Diffusion_Model   
!     
!   PROPOSE: Diffusion model, by MaoWei
! ====================================================================
! =========================Incoming variables=========================
!   Ks          Saturated soil hydrualic conductivity.
! =========================Outcoming variables========================
!   Slope       Slope.
!   Intercept   Intercept.
! =========================related files==============================
!   None.
! =========================related functions==========================
!   None.
! ====================================================================
    SUBROUTINE Diffusion_Model
    USE parm
    IMPLICIT NONE
    
    INTEGER (KIND=KI) :: m
    
    DO m = 1,NMat
        ! Curve fitting.
        IF (Dfit == 1) THEN ! quadratic curve, Dfit>0, theta
            Intercept(m) = -2.9602_KR*log10(par(5,m))**2-0.6179_KR*log10(par(5,m))-4.2946_KR
            Slope(m) = 9.4519_KR*log10(par(5,m))**2+6.6871_KR*log10(par(5,m))+11.115_KR
        ELSEIF (Dfit == 2) THEN ! first segmentation fitting
            IF (log10(par(5,m))>=-0.99866) THEN
                Intercept(m) = -0.6245_KR*log10(par(5,m))**2-0.4761_KR*log10(par(5,m))-5.2872_KR
            ELSE
                Intercept(m) = 9.8771_KR*log10(par(5,m))+4.4293_KR
            ENDIF
            IF (log10(par(5,m))>=-0.99866) THEN
                Slope(m) = 1.6683_KR*log10(par(5,m))**2+6.0643_KR*log10(par(5,m))+14.49_KR
            ELSE
                Slope(m) = -27.225_KR*log10(par(5,m))-17.452_KR
            ENDIF                
        ELSEIF (Dfit == 3) THEN ! second segmentation fitting
            IF (log10(par(5,m))>=-0.9896858) THEN
                Intercept(m) = 0.8771_KR*log10(par(5,m))**2+1.0098_KR*log10(par(5,m))-5.0988_KR
            ELSE
                Intercept(m) = 9.8771_KR*log10(par(5,m))+4.4293_KR
            ENDIF
            IF (log10(par(5,m))>=-0.9990038) THEN
                Slope(m) = -3.9735_KR*log10(par(5,m))**2+0.0609_KR*log10(par(5,m))+13.489_KR
            ELSE
                Slope(m) = -27.225_KR*log10(par(5,m))-17.452_KR
            ENDIF
        ELSEIF (Dfit == -1) THEN ! quadratic curve, Dfit<0, S
            Slope(m) = 2.8657_KR*log10(par(5,m))**2+2.5014_KR*log10(par(5,m))+4.2306_KR
            Intercept(m) = -2.4228_KR*log10(par(5,m))**2-0.4907_KR*log10(par(5,m))-3.5813_KR
        ELSEIF (Dfit == -2) THEN ! segmentation fitting, without amendment, left 6, right 6
            IF (log10(par(5,m))>=-1.004802) THEN
                Intercept(m) = -0.9483_KR*log10(par(5,m))**2-0.4217_KR*log10(par(5,m))-4.1855_KR
            ELSE
                Intercept(m) = 7.8942_KR*log10(par(5,m))+3.2129_KR
            ENDIF
            IF (log10(par(5,m))>=-1.016789) THEN
                Slope(m) = 1.6523_KR*log10(par(5,m))**2+2.5041_KR*log10(par(5,m))+4.961_KR
            ELSE
                Slope(m) = -7.2474_KR*log10(par(5,m))-3.511_KR
            ENDIF
        ELSEIF (Dfit == -3) THEN
            IF (log10(par(5,m))>=-0.8989257) THEN
                Intercept(m) = 0.12_KR*log10(par(5,m))**2+0.511_KR*log10(par(5,m))-4.1381_KR
            ELSE
                Intercept(m) = 5.2687_KR*log10(par(5,m))+0.2358_KR
            ENDIF
            IF (log10(par(5,m))>=-0.9171914) THEN
                Slope(m) = -0.4014_KR*log10(par(5,m))**2+0.5861_KR*log10(par(5,m))+4.518_KR
            ELSE
                Slope(m) = -4.8045_KR*log10(par(5,m))-0.7639_KR
            ENDIF
        ENDIF       

        ! par_n(m) = 15.812*thf(m)**2-10.416*thf(m)+2.8891
        ! par_n(m) = 0.1671*log10(par(5,m))**2+0.7103*log10(par(5,m))+1.8925
        par_n(m) = 0.9509_KR*exp(0.72223_KR*log10(par(5,m)))+0.9_KR

    ENDDO
    
    END SUBROUTINE Diffusion_Model

! ====================================================================
!   SUBROUTINE Boundary  
!     
!   PROPOSE: prepare boundary condition
! ====================================================================
! =========================Incoming variables=========================
!   Bup          Up Boundary Condition -1/0/1/2.
!   Bdn          Down Boundary Condition 0/2.
! =========================Outcoming variables========================
!   None.
! =========================related files==============================
!   None.
! =========================related functions==========================
!   None.
! ==================================================================== 
    SUBROUTINE Boundary_Cond(datapath)
    USE parm
    IMPLICIT NONE
    CHARACTER (LEN=100) :: datapath
    INTEGER (KIND=4) :: lenpath, i, j
    REAL (KIND=KR) :: Nouse, sumres
    REAL (KIND=KR), DIMENSION(Nlayer) :: res

    lenpath = Len_Trim(datapath)
    IF (lWat) THEN
        IF (bup > 0) THEN
            CALL Etp(1,datapath)
            IF (Terr.ne.0) RETURN
            OPEN(110,file=datapath(1:lenpath)//'/'//'01.et0',status='old',err=901)
            READ(110,*,err=901)
            OPEN(130,file=datapath(1:lenpath)//'/'//'01.eti',status='old',err=901)
            OPEN(150,file=datapath(1:lenpath)//'/'//'eta.dat',status='unknown',err=902)
            WRITE(150,*,err=902)"Variables=DoY,   Ea,   Ta,   Date"
            READ(130,*,err=901)
                    
            IF (.NOT. ALLOCATED(precip)) ALLOCATE(precip(2,MaxAL))
            IF (.NOT. ALLOCATED(Evatra)) ALLOCATE(Evatra(2*Nlayer,MaxAL))
            
            DO i = 1,MaxAL
                READ(110,*,err=901) Nouse,precip(1,i),Nouse,Nouse,Nouse,precip(2,i)
                READ(130,*,err=901) Nouse,Nouse,(Evatra(j,i),j=1,2*Nlayer)
            ENDDO
            precip(2,:) = precip(2,:)/1000_KI
            Evatra = Evatra/1000_KI
            CLOSE(110)
            CLOSE(130)
        
        ELSEIF (bup == 0) THEN
            OPEN(120,file=datapath(1:lenpath)//'/'//'met.in',status='old',err=903)
            READ(120,*,err=903)
            READ(120,*,err=903)Nup
            READ(120,*,err=903)
            IF (.NOT. ALLOCATED(up)) ALLOCATE(up(2,Nup))
            READ(120,*,err=903)(up(1,i),i=1,Nup)
            READ(120,*,err=903)
            READ(120,*,err=903)(up(2,i),i=1,Nup)
            DO i = 1,Nup
                up(1,i) = up(1,i)/tConv
                up(2,i) = up(2,i)/xConv
            ENDDO
        
        ELSEIF (bup < 0) THEN
            CALL Etp(1,datapath)
            IF (Terr.ne.0) RETURN
            OPEN(110,file=datapath(1:lenpath)//'/'//'01.et0',status='old',err=901)
            READ(110,*,err=901)
            OPEN(130,file=datapath(1:lenpath)//'/'//'01.eti',status='old',err=901)
            OPEN(150,file=datapath(1:lenpath)//'/'//'eta.dat',status='unknown',err=902)
            WRITE(150,*,err=902)"Variables=DoY,   Ea,   Ta,   Date"
            READ(130,*,err=901)
            IF (.NOT. ALLOCATED(precip)) ALLOCATE(precip(2,MaxAL))
            IF (.NOT. ALLOCATED(Evatra)) ALLOCATE(Evatra(2*Nlayer,MaxAL))
            DO i = 1,MaxAL
                READ(110,*,err=901) Nouse,precip(1,i),Nouse,Nouse,Nouse,precip(2,i)
                READ(130,*,err=901) Nouse,Nouse,(Evatra(j,i),j=1,2*Nlayer)
            ENDDO
            precip(2,:) = precip(2,:)/1000_KI
            Evatra = Evatra/1000_KI
            CLOSE(110)
            CLOSE(130)
        
            OPEN(120,file=datapath(1:lenpath)//'/'//'met.in',status='old',err=903)
            READ(120,*,err=903)
            READ(120,*,err=903)Nup
            READ(120,*,err=903)
            IF (.NOT. ALLOCATED(up)) ALLOCATE(up(2,Nup))
            READ(120,*,err=903)(up(1,i),i=1,Nup)
            READ(120,*,err=903)
            READ(120,*,err=903)(up(2,i),i=1,Nup)
            DO i = 1,Nup
                up(1,i) = up(1,i)/tConv
                up(2,i) = up(2,i)/xConv
            ENDDO
        ENDIF
        
!       Lower boundary condition
!       Given lower soil water content.
        IF (bdn == -1) THEN
            IF (bup > 0) OPEN(120,file=datapath(1:lenpath)//'/'//'met.in',status='old',err=903)
            READ(120,*,err=903)
            READ(120,*,err=903)Ndn
            READ(120,*,err=903)
            IF (.NOT. ALLOCATED(dn)) ALLOCATE(up(2,Ndn))
            READ(120,*,err=903)(dn(1,i),i=1,Ndn)   !
            READ(120,*,err=903)
            READ(120,*,err=903)(dn(2,i),i=1,Ndn)
            DO i = 1,Ndn
                dn(1,i) = dn(1,i)/tConv
            ENDDO
!       Given groundwater table depth.
        ELSEIF (bdn == -2) THEN
            IF (bup > 0) OPEN(120,file=datapath(1:lenpath)//'/'//'met.in',status='old',err=903)
            READ(120,*,err=903)
            READ(120,*,err=903) gdep
            gdep = gdep/xConv
        ENDIF
    ELSE
        IF (rET>Tol) THEN
            CALL Steady_ET(1,datapath,res)
        ELSE
            res=0.0_KR
        ENDIF
        qirr = rTop
        sumres = sum(res)
        DO i = 1,Nlayer
            IF (rTop-rBot-sumres<Tol) THEN
                qflux(i,1) = rTop*dt
                qflux(i,2) = rBot*dt
            ELSEIF (rTop+rBot-sumres<Tol) THEN
                qflux(i,1) = rTop*dt
                qflux(i,2) = rBot*dt
            ELSEIF (rTop+sumres-rBot<Tol) THEN
                qflux(i,1) = rTop*dt
                qflux(i,2) = rBot*dt
            ELSE
                Terr = 6
            ENDIF
            Sink1d(i) = res(i)*dt
        ENDDO
    ENDIF
!   Solute boundary condition
    IF (lchem) THEN
        IF (bup>0 .or. (.NOT.lWat)) OPEN(120,file=datapath(1:lenpath)//'/'//'met.in',status='old',err=903)
        READ(120,*,err=903)
        READ(120,*,err=903)CNup
        READ(120,*,err=903)
        IF (.NOT. ALLOCATED(Cup)) ALLOCATE(Cup(2,CNup))
        READ(120,*,err=903)(Cup(1,i),i=1,CNup)
        READ(120,*,err=903)
        READ(120,*,err=903)(Cup(2,i),i=1,CNup)
        DO i = 1,CNup
            Cup(1,i) = Cup(1,i)/tConv
            Cup(2,i) = Cup(2,i)*xConv**3/mConv
        ENDDO
        READ(120,*,err=903)
        READ(120,*,err=903)CNdn
        READ(120,*,err=903)
        IF (.NOT. ALLOCATED(Cdn)) ALLOCATE(Cdn(2,CNdn))
        READ(120,*,err=903)(Cdn(1,i),i=1,CNup)
        READ(120,*,err=903)
        READ(120,*,err=903)(Cdn(2,i),i=1,CNup)
        DO i = 1,CNdn
            Cdn(1,i) = Cdn(1,i)/tConv
            Cdn(2,i) = Cdn(2,i)*xConv**3/mConv
        ENDDO
    ENDIF
    CLOSE(120)

    RETURN
    
901 Terr=3
    RETURN
902 Terr=4
    RETURN
903 Terr=5
    RETURN
    END SUBROUTINE Boundary_Cond
    
! ====================================================================
!   Subroutine Set_Input
!     
!   Purpose: Set the boundary input values
! ====================================================================
    SUBROUTINE Set_Input
    USE parm
    IMPLICIT NONE
    INTEGER (KIND=KI) :: k, kk, i, m
    REAL (KIND=KR) :: rc, gdepold
    REAL (KIND=KR), DIMENSION(Nlayer) :: tho
    
    IF (bup > 0) THEN
        ! Read precipitation, and set irrigation to zero.
        CALL FindY_Step(MaxAL,t,qair,precip,k)
        qirr = 0.0_KR
        Epi(1:Nlayer) = Evatra(1:NLayer,k)
        Tri(1:Nlayer) = Evatra(NLayer+1:2*Nlayer,k)
        tAtm = k
        
        IF (lchem) THEN
            Concup = 0.0_KR
            CALL FindY_Step(CNdn,t,Concdn,Cdn,kk)
        ENDIF

    ELSEIF (bup == 0) THEN
        ! Read irrigation, and set precipitation to zero.
        qair = 0.0_KR
        CALL FindY_Step(Nup,t,qirr,up,k)
        IF (lchem) THEN
            CALL FindY_Step(CNup,t,Concup,Cup,kk)
            CALL FindY_Step(CNdn,t,Concdn,Cdn,kk)
        ENDIF
        
    ELSEIF (bup < 0) THEN
        CALL FindY_Step(MaxAL,t,qair,precip,k)
        CALL FindY_Step(Nup,t,qirr,up,kk)
        Epi(1:Nlayer) = Evatra(1:NLayer,k)
        Tri(1:Nlayer) = Evatra(NLayer+1:2*Nlayer,k)
        tAtm = k
        
        IF (lchem) THEN
            CALL FindY_Step(CNup,t,Concup,Cup,kk)
            CALL FindY_Step(CNdn,t,Concdn,Cdn,kk)
        ENDIF
    ENDIF
    
    IF (bdn == -1) THEN
        CALL FindY_Step(Ndn,t,dth,dn,kk)
    ELSEIF (bdn == -2) THEN
        IF (Tlevel == 1) THEN
            ! Keep the Wlayer steady, change the Nlayer to adjust the groundwater.
            Wlayer = Nlayer
        ELSE
            rc = qflux(Nlayer,2)
            IF (rc >= 0) THEN ! Unsaturated zone to saturated zone, groundwater table upward.
                tho = th
                DO i = Nlayer,1,-1
                    m = MATuz(i)
                    tho(i) = th(i) + rc/dz(i)
                    IF (tho(i) > par(2,m)) THEN
                        rc = (tho(i)-par(2,m))*dz(i)
                        tho(i) = par(2,m)
                    ELSE
                        gdep = gdep-sum(dz(i:Nlayer))+(dz(i)-rc/(par(2,m)-th(i)))
                        EXIT
                    ENDIF
                ENDDO
            ELSE
                gdepold = gdep
                DO i = Nlayer+1,Wlayer
                    m = MATuz(i)
                    tho(i) = par(2,m) + rc/dz(i)
                    IF (tho(i) < thF(m)) THEN
                        rc = rc + (par(2,m)-thF(m))*dz(i)
                        tho(i) = thF(m)
                    ELSE
                        gdep = gdep+sum(dz(Nlayer+1:i))-(dz(i)+rc/(par(2,m)-thF(m)))
                        EXIT
                    ENDIF
                ENDDO
            ENDIF
            ! Adjust dz and th
            IF (gdep > zx(Nlayer) .and. gdep < zx(Nlayer+1)) THEN
                m = MATuz(Nlayer)
                th(Nlayer) = (th(Nlayer)*dz(Nlayer)+thF(m)*(gdep-gdepold)+par(2,m)*(zx(Nlayer+1)-gdep))/(zx(Nlayer+1)-zx(Nlayer))
                dz(Nlayer) = zx(Nlayer+1)-zx(Nlayer)
            ELSE
                m = MATuz(Nlayer)
                th(Nlayer) = (th(Nlayer)*dz(Nlayer)+par(2,m)*(zx(Nlayer+1)-zx(Nlayer)-dz(Nlayer)))/(zx(Nlayer+1)-zx(Nlayer))
                dz(Nlayer) = zx(Nlayer+1)-zx(Nlayer)
                DO i = Nlayer+1,Wlayer
                    m = MATuz(i)
                    th(i) = thF(m)
                ENDDO
            ENDIF
        ENDIF
        
        IF (gdep > maxval(zx)) THEN
            WRITE(*,*) "The groundwater is too deep!"
            PAUSE
            STOP
        ENDIF
        
        DO i = 1,Wlayer
            IF (zx(i) < gdep .and. zx(i+1) >= gdep) THEN
                Nlayer = i
            ENDIF
        ENDDO
        dz(Nlayer) = gdep - zx(Nlayer)
        dz(Nlayer+1) = dz(Nlayer-1)
    ENDIF
   
    END SUBROUTINE Set_Input

! ====================================================================
!   Subroutine FindY_Step 
!     
!   Purpose: interpolation
! ====================================================================
! =========================Incoming variables=========================
!   n           The total numbers of the node.
!   x           The input value.
!   xy(2,:)     The first row is the value corresponding to x,
!               the second row is the value corresponding to y.
! =========================Outcoming variables========================
!   y           The output value.
! ====================================================================
    SUBROUTINE FindY_Step(n,x,y,xy,flag)
    USE parm

    INTEGER (kind=KI) :: n,flag
    REAL (kind=KR) :: x
    REAL (kind=KR) :: y
    REAL (kind=KR), DIMENSION(2,n) :: xy
    
    IF (x <= xy(1,1)) THEN
        Y=xy(2,1)
        flag = 1
        RETURN
    ENDIF
    IF (x >= xy(1,n)) THEN
        y=xy(2,n)
        flag = n
        RETURN
    ENDIF
    DO i=1,n-1
        IF (x > xy(1,i) .and. x <= (xy(1,i+1)+Tol)) THEN
            y=xy(2,i+1)
            flag = i+1
            RETURN
        ENDIF   
    ENDDO

    END SUBROUTINE FindY_Step
