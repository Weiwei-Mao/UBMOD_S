! ====================================================================
!     subroutine ETp
!     
!     purpose: 1. calculate the ET0 and Ep and Tp in the ??.et0 file.
!              2. calculate the Crop growth in the simulation period
!                 in the ??.crp file.
!              3. calculate the Epi[mm] and Tpi[mm] in every 1d layers
!                 in the ??.eti file.
! ====================================================================
! =========================Incoming variables=========================
!       num            The number of the 1D soil columns[-]
!       sp(:,:)        The discreted function for soil evaporation[-]
!       dx(:,:)        The z coordinates of  the discreted layers[m]
! =========================Outcoming variables========================
!       ??.ETi         file stored Epi[mm] and Tpi[mm] in every layer
! =========================related files==============================
!       ??.et0         ET0 file with Ep[mm] and Tp[mm]
!       ??.crpcrop     file with root depth[m]
!       cropdat.dat    crop datebase with the discreted function of 
!                      crop transpiration	
!       ??.ETi         file stored Epi[mm] and Tpi[mm] in every layer 
! =========================related functions==========================
!       cdatebase      The datebase to find discreted function for 
!                      soil evaporation and crop transpiration
!       crop           Find root depth[m]
!       Fdense         Find the contribution coefficients of the 
!                      discreted layers
!       Ep_Tp          Call for Ep[mm] and Tp[mm]
! ====================================================================
    SUBROUTINE ETp(num,datapath)
    USE parm
    IMPLICIT NONE

    CHARACTER (LEN=6) :: cha1, cha2, cha3
    CHARACTER (LEN=8) :: date1
    CHARACTER (LEN=8), DIMENSION(5) :: da
    CHARACTER (LEN=100) :: datapath
    INTEGER (KIND=4) :: i,j,k,ii,jj,num,lenpath
    REAL (KIND=KR), DIMENSION(Nlayer) :: dz1,sed,dense
    REAL (KIND=KR), DIMENSION(21) :: crpdat 
    REAL (KIND=KR), DIMENSION(Nlayer) :: work
    REAL (KIND=KR) :: ET0, EP, TP, Nouse, rd, dtp

    lenpath = Len_Trim(datapath)
    cha1="??.et0"
    cha2="??.crp"
    cha3="??.eti"

    date1 = date

    DO i=1,Nlayer
        dz1(i)=zx(i+1)  
    ENDDO

    IF (lCrop) THEN
        CALL Crop(datapath,num,numc,date1,MaxAL)
        IF (Terr.ne.0) RETURN
    ENDIF
    CALL ep_tp(datapath,num,date1,MaxAL,IfPM)
    IF (Terr.ne.0) RETURN

    WRITE(cha1(1:2),"(i2.2)")num
    WRITE(cha3(1:2),"(i2.2)")num
    OPEN(4,file=datapath(1:lenpath)//'/'//cha1,status="unknown",err=901)
    OPEN(6,file=datapath(1:lenpath)//'/'//cha3,status="unknown",err=901)
    WRITE(6,*,err=901)"date ordinal Es(1:Nlayer)/mm    TS(1:Nlayer)/mm"
    READ(4,*,err=901)cha1
    IF (lCrop) THEN
        WRITE(cha2(1:2),"(i2.2)")num
        OPEN(5,file=datapath(1:lenpath)//'/'//cha2,status="unknown",err=901)
        READ(5,*,err=901)cha1
    ENDIF

    DO j=1,MaxAL
        READ(4,*,err=901)date1,Nouse,ET0,Ep,Tp

        dtp=max(0.05_KR*ET0, 0.05_KR)
        dtp=dtp/4.0_KR
        DO ii=1,NMat
            DO jj=1,4
                work(jj)=sp(jj,ii)
            ENDDO
        ENDDO
        CALL Fdense(Nlayer,dz1,dtp,work,sed,4)
        
        IF (lCrop) THEN
            READ(5,*,err=901)Nouse,Nouse,rd
            
            CALL cdatebase(datapath,numc,da,crpdat)
            IF (Terr.ne.0) RETURN
            DO jj=1,10
                work(jj)=crpdat(jj+10)
            ENDDO
            ptab = crpdat(21)     
            CALL Fdense(Nlayer,dz1,0.1_KR*rd,work,dense,10)
        ENDIF
        DO k=1,Nlayer
            sed(k)=sed(k)*ep
            IF (lCrop) THEN
                dense(k)=dense(k)*tp
            ELSE
                dense(k)=0.0_KR
            ENDIF
        ENDDO
        WRITE(6,"(A9,I6,260f6.3)",err=901)date1,i,sed(1:Nlayer),dense(1:Nlayer)
    ENDDO
    CLOSE(4)
    IF (lCrop) CLOSE(5)
    CLOSE(6)
    RETURN
    
901 Terr=3
    RETURN
    END SUBROUTINE ETp

    
! ====================================================================
!     Subroutine Steady ET  
!     
!     Purpose: Calculate the ET under steady condition.
! ====================================================================
    SUBROUTINE Steady_ET(num,datapath,res)
    USE parm
    IMPLICIT NONE
    CHARACTER (LEN=100) :: datapath
    CHARACTER (LEN=8), DIMENSION(5) :: da
    REAL (KIND=KR) :: dtp, ckc, rd, clai, claif, ETc, Ep, Tp
    REAL (KIND=KR) :: f,ckc1,ckc2,ckc3,rd1,rd2,clai1,clai2,clai3,clai4
    REAL (KIND=KR), DIMENSION(21) :: crpdat
    REAL (KIND=KR), DIMENSION(Nlayer) :: work,dz1,sed,dense,res
    INTEGER (KIND=4) :: ii, jj, num, jd1, Jdate
    INTEGER (KIND=4), DIMENSION(5) :: jd
    
    IF (lCrop) THEN
        CALL cdatebase(datapath,num,da,crpdat)
        IF (Terr.NE.0) RETURN
        ! Calculate the potential Ep and Tp.
        f=crpdat(1)
        ckc1=crpdat(2)
        ckc2=crpdat(3)
        ckc3=crpdat(4)
        rd1=crpdat(5)
        rd2=crpdat(6)
        clai1=crpdat(7)
        clai2=crpdat(8)
        clai3=crpdat(9)
        clai4=crpdat(10)
        
        jd1=Jdate(date)
        DO ii=1,5
            jd(ii)=Jdate(da(ii))
        ENDDO
        IF (jd1<jd(1).or.jd1>jd(5))THEN
            ckc=0.9 
            rd=0
            clai=0
        ELSEIF (jd1.ge.jd(1).and.jd1<jd(2)) THEN
            ckc=ckc1
            rd=rd1+(rd2-rd1)*(jd1-jd(1))/(jd(3)-jd(1))
            clai=CLAI1+(CLAI2-CLAI1)*(jd1-jd(1))/(jd(2)-jd(1))
        ELSEIF (jd1.ge.jd(2).and.jd1<jd(3)) THEN
            ckc=ckc1+(ckc2-ckc1)*(jd1-jd(2))/(jd(3)-jd(2))
            rd=rd1+(rd2-rd1)*(jd1-jd(1))/(jd(3)-jd(1))
            clai=CLAI2+(CLAI3-CLAI2)*(jd1-jd(2))/(jd(3)-jd(2))
        ELSEIF (jd1.ge.jd(3).and.jd1<jd(4)) THEN
            ckc=ckc2
            rd=rd2
            clai=clai3
        ELSEIF (jd1.ge.jd(4).and.jd1.le.jd(5)) THEN
	        ckc=ckc2-(ckc2-ckc3)*(jd1-jd(4))/(jd(5)-jd(4))
            rd=rd2
            clai=CLAI3-(CLAI3-CLAI4)*(jd1-jd(4))/(jd(5)-jd(4))
        ENDIF
        Claif=exp(-f*clai)
        ETc=ckc*rET
        Ep = Claif*ETc
        Tp = ETc-Ep
    ELSE
        Ep = rET
        Tp = 0.0_KR
    ENDIF
        
    ! Redistribution to each layer.
    dtp = max(0.05_KR*rET, 0.05_KR)
    dtp = dtp/4.0_KR
    DO ii=1,NMat
        DO jj=1,4
            work(jj)=sp(jj,ii)
        ENDDO
    ENDDO
    DO ii=1,Nlayer
        dz1(ii)=zx(ii+1)  
    ENDDO
    CALL Fdense(Nlayer,dz1,dtp,work,sed,4)
    
    IF (lCrop) THEN
        DO jj=1,10
            work(jj)=crpdat(jj+10)
        ENDDO
        ptab = crpdat(21)
        CALL Fdense(Nlayer,dz1,0.1_KR*rd,work,dense,10)
    ENDIF
    DO ii=1,Nlayer
        sed(ii)=sed(ii)*ep
        IF (lCrop) THEN
            dense(ii)=dense(ii)*tp
        ELSE
            dense(ii)=0.0_KR
        ENDIF
    ENDDO
    res = sed+dense

    RETURN
    END SUBROUTINE Steady_ET


! ====================================================================
!   Subroutine ep_tp 
!     
!   Purpose: calculatin the potential Ep and Tp
! ====================================================================
! =========================Incoming variables=========================
!   num            The number of the crop[-]
!   dateini        The beginning day of simulation[yyyymmdd]
!   interval       The duration or length of the simulaiton day[-]
! =========================Outcoming variables========================
!   ??.ET0         A file stored ET0[m],Ep[m],Tp[m],PE[m] 
! =========================related files==============================
!   ??.crp         The file stored information about Kc and LAIF is 
!                  needed.
!   ??.ET0         A file stored ET0[m],Ep[m],Tp[m],PE[m] 
! =========================related functions==========================
!   None.
! ====================================================================
    SUBROUTINE ep_tp(datapath,num,dateini,interval,IfPM)
    USE parm, ONLY : KR, KI, Terr, lCrop
    IMPLICIT NONE
    
    LOGICAL (KIND=4) :: IfPM
    CHARACTER (LEN=6) :: cha1, cha2, cha3
    CHARACTER (LEN=8) :: dateini
    CHARACTER (LEN=100) :: datapath
    INTEGER (KIND=KI) ::interval
    !INTEGER (kind=KI) :: yearini, monthini, dayini
    INTEGER (KIND=4) :: num, ierr, i, year, month, day, julie, DM, lenpath
    REAL (KIND=KR) :: tmean,tmax,tmin,nact,rhmean,uh,ea,deta,dr,si,tansi,&
    tmean1,ws,ra,np,rns,eamax,eamin,ed,rnl,rn,g,lamta,u2,gama,ato_p,et0,dt,PE
    REAL (KIND=KR) :: fai,hhh,hu
    REAL (KIND=KR) :: aa, ckc, claif, etc, ep, tp

    lenpath = Len_Trim(datapath)

    IF (lCrop) THEN
        cha1="??.wea"
        cha2="??.et0"
        cha3="??.crp"
        WRITE(cha1(1:2),"(i2.2)")num
        WRITE(cha2(1:2),"(i2.2)")num
        WRITE(cha3(1:2),"(i2.2)")num
        OPEN(1,file=datapath(1:lenpath)//'/'//cha1,status="old",err=901)
        OPEN(3,file=datapath(1:lenpath)//'/'//cha3,status="old",err=902)
        OPEN(2,file=datapath(1:lenpath)//'/'//cha2,status="unknown",err=902)
        WRITE(2,*)"year month day ordinal et0 ep tp	pe"
        read(3,*,err=902) !/lai and kc/
        
        IF (IfPM) THEN
            read(1,*,err=901)
            read(1,*,err=901)fai,hhh,hu
            read(1,*,err=901)
        ELSE
            read(1,*,err=901)
        ENDIF
	    	
        DO i=1,interval
            IF (IfPM) THEN
                READ(1,*,err=901)year,month,day,tmean,tmax,tmin,nact,rhmean,uh,ato_p,PE
                CALL Refet(year,month,day,Tmean,tmax,tmin,nact,RHmean,uh,hu,ato_p,hhh,fai,ET0,tmean1)
            ELSEIF (.not.IfPM) THEN
                READ(1,*,err=901)year,month,day,et0,PE
            ENDIF
      
            READ(3,*,err=902)aa,aa,ckc,aa,aa,CLAif
            ETc=ckc*ET0 
            
!            CLAIF = 0.25    !2017-06-02
            
            Ep=cLAIf*ETc    !CLAI exp[-f*LAI]
            Tp=ETc-Ep 
            WRITE(2,"(i4,i2.2,i2.2,i6,4f8.3)",err=902)year,month,day,i,et0,ep,tp,PE
        ENDDO
        
        CLOSE(1)
        CLOSE(2)
        CLOSE(3)
    ELSE
        cha1="??.wea"
        cha2="??.et0"
        WRITE(cha1(1:2),"(i2.2)")num
        WRITE(cha2(1:2),"(i2.2)")num
        OPEN(1,file=datapath(1:lenpath)//'/'//cha1,status="old",err=901)
        OPEN(2,file=datapath(1:lenpath)//'/'//cha2,status="unknown",err=902)
        WRITE(2,*)"year month day ordinal et0 ep tp	pe"
        
        IF (IfPM) THEN
            read(1,*,err=901)
            read(1,*,err=901)fai,hhh,hu
            read(1,*,err=901)
        ELSE
            read(1,*,err=901)
        ENDIF
        
        DO i=1,interval
            IF (IfPM) THEN
                READ(1,*,err=901)year,month,day,tmean,tmax,tmin,nact,rhmean,uh,ato_p,PE
                CALL Refet(year,month,day,Tmean,tmax,tmin,nact,RHmean,uh,hu,ato_p,hhh,fai,ET0,tmean1)
            ELSEIF (.not.IfPM) THEN
                READ(1,*,err=901)year,month,day,et0,PE
            ENDIF

            Ep=ET0
            Tp=0.0_KR 
            WRITE(2,"(i4,i2.2,i2.2,i6,4f8.3)",err=902)year,month,day,i,et0,ep,tp,PE
        ENDDO
        
        CLOSE(1)
        CLOSE(2)
    ENDIF
    RETURN
    
901 Terr=1
    RETURN
902 Terr=3
    RETURN

    END SUBROUTINE ep_tp

! ====================================================================
!   Subroutine CROP  
!     
!   Purpose: Calculation Kc, RootDepth, LAI, CLAI for each day, by an 
!            interpolation method.
! ====================================================================
! =========================Incoming variables=========================
!   num        The number of the 1D soil columns[-]
!   numc       The number of the kind of crop[-]
!   dateini    The beginning day of simulation[yyyymmdd]
!   interval   The duration or length of the simulaiton day[-]
! =========================Outcoming variables========================
!   ??.crp     A file stored Kc[-],RootDepth[m],LAI[-],CLAI[-]. 
! =========================related files==============================
!   ??.crp     cropdat.dat    crop datebase.
! =========================related functions==========================
!   None.
! ====================================================================	
    SUBROUTINE Crop(datapath,num,numc,dateini,interval)
    USE parm, ONLY : KR, KI, Terr
    IMPLICIT NONE

    CHARACTER (len=8) :: Dateini, Date, Date1, Db
    CHARACTER (len=8), DIMENSION(5) :: da
    CHARACTER (len=6) :: nameF
    CHARACTER (len=10) :: cha1
    CHARACTER (LEN=100) :: datapath
    INTEGER (kind=KI) interval	! the simulation period
    INTEGER (KIND=4), DIMENSION(5) :: jd
    INTEGER (KIND=4) :: numc, nyear, nyear1, ndyear, JDATE, jd1, jmax, nyrr, dd, ifrun
    REAL (KIND=KR), DIMENSION(21) :: crpdat
    REAL (KIND=KR) :: f, ckc1, ckc2, ckc3, rd1, rd2, clai1, clai2, clai3, clai4
    REAL (KIND=KR) :: ckc, rd, clai, claif
    INTEGER (KIND=4) :: lenpath, num, i, j, k
    
    lenpath = Len_Trim(datapath)   
    namef="??.crp"
    WRITE(namef(1:2),"(i2.2)")num
    OPEN(11,file=datapath(1:lenpath)//'/'//nameF,status="unknown",err=901)
    WRITE(11,*)"date    ordinal Kc  rootdepth   LAI LAIF    rp(1:Nlayer)"
    CALL cdatebase(datapath,numc,da,crpdat)
    IF (Terr.ne.0) RETURN

    f=crpdat(1)
    ckc1=crpdat(2)
    ckc2=crpdat(3)
    ckc3=crpdat(4)
    rd1=crpdat(5)
    rd2=crpdat(6)
    clai1=crpdat(7)
    clai2=crpdat(8)
    clai3=crpdat(9)
    clai4=crpdat(10)

    ! Reset the year of the dataset to the calculation year.
    READ(da(1)(1:4),"(i4)")nyear
    READ(Dateini(1:4),"(i4)")nyear1
    ndyear=nyear1-nyear
    DO i=1,5
        READ(da(i)(1:4),"(i4)")nyear
        WRITE(da(i)(1:4),"(i4)")nyear+ndyear
    ENDDO

    DO i=1,5
        jd(i)=Jdate(da(i))
    ENDDO
    
    date=dateini
    jd1=Jdate(date)

    DO i=1,interval
        IF (jd1<jd(1).or.jd1>jd(5))THEN
            ckc=0.9 
            rd=0
            clai=0
        ELSEIF (jd1.ge.jd(1).and.jd1<jd(2)) THEN
            ckc=ckc1
            rd=rd1+(rd2-rd1)*(jd1-jd(1))/(jd(3)-jd(1))
            clai=CLAI1+(CLAI2-CLAI1)*(jd1-jd(1))/(jd(2)-jd(1))
        ELSEIF (jd1.ge.jd(2).and.jd1<jd(3)) THEN
            ckc=ckc1+(ckc2-ckc1)*(jd1-jd(2))/(jd(3)-jd(2))
            rd=rd1+(rd2-rd1)*(jd1-jd(1))/(jd(3)-jd(1))
            clai=CLAI2+(CLAI3-CLAI2)*(jd1-jd(2))/(jd(3)-jd(2))
        ELSEIF (jd1.ge.jd(3).and.jd1<jd(4)) THEN
            ckc=ckc2
            rd=rd2
            clai=clai3
        ELSEIF (jd1.ge.jd(4).and.jd1.le.jd(5)) THEN
	        ckc=ckc2-(ckc2-ckc3)*(jd1-jd(4))/(jd(5)-jd(4))
            rd=rd2
            clai=CLAI3-(CLAI3-CLAI4)*(jd1-jd(4))/(jd(5)-jd(4))
        ENDIF
        Claif=exp(-f*Clai)

        WRITE(11,"(A9,I6,14f6.3)")date,i,ckc,rd,clai,claif !
        CALL dateadd(date,1,date1)
        date=date1
        jd1=jd1+1
        
        ! If the calculation date is out of the database range !
        jmax=jd(5) ! the maximum day number
        IF (jd1>jmax) THEN
            READ(da(1)(1:4),"(i4)")nyear
            READ(da(5)(1:4),"(i4)")nyear1
            IF (nyear == nyear1) THEN
                IF (ifrun(nyear) == 0) THEN
                    nyrr = 365
                ELSE
                    Db = da(1)
                    WRITE(Db(5:8),"(i4)")0229
                    dd = Jdate(Db)
                    IF (jd(1)>dd) THEN
                        nyrr = 365
                    ELSE
                        nyrr = 366
                    ENDIF
                ENDIF
            ELSE
                IF ((ifrun(nyear)==0) .and. (ifrun(nyear+1)==0)) THEN
                    nyrr = 365
                ELSEIF (ifrun(nyear)==1) THEN
                    Db = da(1)
                    WRITE(Db(5:8),"(i4)")0229
                    dd = Jdate(Db)
                    IF (jd(1)>dd) THEN
                        nyrr = 365
                    ELSE
                        nyrr = 366
                    ENDIF
                ELSE
                    Db = da(5)
                    WRITE(Db(5:8),"(i4)")0229
                    dd = Jdate(Db)
                    IF (dd>jd(5)) THEN
                        nyrr = 365
                    ELSE
                        nyrr = 366
                    ENDIF
                ENDIF
            ENDIF
            jd1=jd1-Nyrr
        ENDIF                  

    ENDDO
    CLOSE(11)
    RETURN

901 Terr=3
    RETURN
    
    END SUBROUTINE Crop

! ====================================================================
!   Subroutine Fdense  
!     
!   Purpose: Calculation the root density function and evaporation
!            cumulative distribution function.
! ====================================================================
! =========================Incoming variables=========================
!   n1         the number for interpolation,to T is 10, to E is 4.[-]
!   n          the layer number.
!   dz1(:)     the depth of the lower boundary of the layer in its
!              controled area.[m] 
!   re         the actual length of root or depth of evaporation.[m]
!   work(:)    the discreted function for interpolation.[-]
! =========================Outcoming variables========================
!   dense(:)   the interpolation results.[-]
! =========================related files==============================
!   None.
! =========================related functions==========================
!   None.
! ====================================================================
    SUBROUTINE Fdense(mlayer,dz1,re,work,dense,N1)
    USE parm
    IMPLICIT NONE

    REAL (KIND=KR), DIMENSION(Nlayer) :: dz1
    REAL (KIND=KR), DIMENSION(Nlayer) :: dense
    REAL (KIND=KR), DIMENSION(10) :: work
    REAL (KIND=KR) :: r2, re, r1
    INTEGER (KIND=KI) :: mlayer
    INTEGER (KIND=4) ::  i, j, L, N1
    
    DO i=1,mlayer
        dense(i)=0.0_KR
    ENDDO
    
    r2=0.0_KR

    DO j=1,mlayer
        L=int(dz1(j)/re)
        
        IF (l-0.0_KR<Tol) THEN
            r1=dz1(j)/re*work(L+1)
        ELSEIF (L>0.0_KR.and.L<N1) THEN
            r1=work(L)+(dz1(j)-re*L)/re*(work(L+1)-work(L))		
        ELSEIF (L.GE.N1) THEN
            r1=1.0_KR
        ENDIF

        dense(j)=r1
        dense(j)=r1-r2
        r2=r1
    ENDDO
    
    IF (re-0.0_KR<Tol) THEN
        dense = 0.0_KR
    ENDIF
    
    END SUBROUTINE Fdense

! ====================================================================
!   Subroutine cdatebase  
!     
!   Purpose: read the information in crop datebase for crop NumC
! ====================================================================
! =========================Incoming variables=========================
!   Numc            number of the crop[-]
! =========================Outcoming variables========================
!   da(:)           date of the 5 stage of crop[yyyymmdd]
!   crpdat(:)       other numerical information such as Kc,LAI,f,
!                   root depth[m].
! =========================related files==============================
!   cropdat.dat     crop datebase with the discreted function of 
!                   transpiration
! =========================related functions==========================
!   None.
! ====================================================================
    SUBROUTINE cdatebase(datapath,Numc,da,crpdat)
    USE parm, ONLY : KR, KI, Terr
    IMPLICIT NONE

    REAL (kind=KR), DIMENSION(21) :: crpdat
    CHARACTER (len=8), DIMENSION(5) :: da
    CHARACTER (len=10) :: nameC	
    CHARACTER (LEN=100) :: datapath
    INTEGER (KIND=4) :: lenpath, i, numc

    lenpath = Len_Trim(datapath)
    OPEN(1,file=datapath(1:lenpath)//'/cropdat.dat',status="old",err=901)
    
    DO i=1,numc-1
        READ(1,*,err=901)
    ENDDO
    READ(1,*,err=901)nameC,da(1:5),crpdat(1:21)

    CLOSE(1)
    RETURN
    
901 Terr=2
    RETURN
    END SUBROUTINE cdatebase

! ====================================================================
!     subroutine refet  
!     
!     purpose: this is a subrutine to cacalate the reference et.the 
!              equation comes from the PM equation(1992) recommended 
!              by FAO(1998)
! ====================================================================
! =========================Incoming variables=========================
!   year        year[yyyy]
!   month       month[mm]
!   day         day[dd]
!   Tmean       the average temperature[C]
!   Tmax        the maximum temperature[C]
!   Tmin        the minimum temperature[C]
!   Nact        the number of hours of actual sunshine[h] 
!   RHmean      the average humidty[%]
!   Uh          the speed of wind[m/s]
!   ato_p       the atmospheric pressure[Kpa]
!   hhh         the sea level elevation[m]
!   fai         the latitude[radian]
!   Hu          the Wind speed measuring height[m].
! =========================Outcoming variables========================
!   Et0         the reference evapotranspiration[m/d]
! =========================related files==============================
!   ??.wea      weather information for ET.
! =========================related functions==========================
!   None.
! ====================================================================
    SUBROUTINE Refet(year,month,day,Tmean,tmax,tmin,nact,RHmean,uh,hu,ato_p,hhh,fai,ET0,tmean1)
    USE parm, ONLY : KR, KI
    IMPLICIT NONE

    INTEGER (KIND=4), DIMENSION(12) :: Dm ! Calculation ordinal number
    INTEGER (KIND=4) :: i,year,month,day,julie
    REAL (KIND=KR) :: tmean,tmax,tmin,nact,rhmean,uh,ea,deta,dr,si,tansi, &
         tmean1,ws,ra,np,rns,eamax,eamin,ed,rnl,rn,g,lamta,u2,gama,ato_p,et0
    REAL (KIND=KR) :: fai, hhh, hu, amta
    
    DATA DM/0,31,59,90,120,151,181,212,243,273,304,334/

    julie=day+DM(MONTH)

    ea=0.611*EXP(17.27*tmean/(tmean+237.3))
    deta=4098*ea/(tmean+237.3)**2
    dr=1+0.033*COS(0.0172*julie)
    si=0.409*SIN(0.0172*julie-1.39)
    tansi=tan(si)
    ws=acos(-TAN(fai)*tansi)
    ra=37.6*dr*(ws*SIN(fai)*SIN(si)+COS(fai)*COS(si)*SIN(ws))
    Np=7.64*ws
    rns=0.77*(0.25+0.5*nact/Np)*ra
    eamax=0.611*EXP(17.27*tmax/(tmax+237.3))
    eamin=0.611*EXP(17.27*tmin/(tmin+237.3))
    ed=(eamax+eamin)/200*RHmean
    Rnl=2.45E-9*(0.9*nact/Np+0.1)*(0.34-0.14*ed**0.5)*((tmax+273)**4+(tmin+273)**4)
    rn=rns-rnl

    IF (tmean1.eq.0) THEN
        g=0.
    ELSE
        g=0.38*(tmean-tmean1)
    ENDIF
    tmean1=tmean
    amta=2.501-(2.361e-3)*tmean
    u2=4.87*uh/Log(67.8*hu-5.42)
    gama=0.00163*ato_p/lamta
    et0=(0.408*deta*(rn-G)+gama*900/(273+tmean)*u2*(ea-ed))/(deta+gama*(1+0.34*u2))
    et0=max(et0,0.0)    ! modified by MW 2016-09-26 17:36

    END SUBROUTINE Refet
     
