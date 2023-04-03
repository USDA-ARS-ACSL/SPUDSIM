      

      subroutine SoluteUptake()
      include 'public.ins'
      IncLude 'Puplant.ins'
      Include 'Puweath.ins'
      real MMUpN, bi(3),ci(3)
      Common  / SUP /  TotWSink,TotSSink,WincrSink,TotSSINK2
      Common  / Biom /  BM
      Common  / slpt /  Iflag
      Common / Mass_root_M/Isink,Rroot,ConstI(2), constK(2),
     !    Cmin0(2),DlngR(NMatD),
     !	Disp(NumElD,2),RootExchange(2),Csink_M(NumElD),Cr_M(NumElD,2)
      Common  / SoluteM /  ThOld(NumNPD), Dmol(NumSD)
      Common  / PotNitr /  PotNitrogen_t
      Character * 40  UptakeN(6)
      data UptakeN / 
     !'Passive Mass Uptake (Mihaelson-Menton kinetics'
     !,'Passive Mass Uptake (Mihaelson-Menton kinetics)&convetcive flux'
     !,'Active  Mass Uptake  Cr = Croot',
     !' ',
     !' ',
     !' ' / 
c       real qc(NumElD),DispM(2,NumElD)
C *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  * 
	If(lInput.eq.1) then
        Open(94,file = 'TotSink.res')
		TotWSink = 0.
		TotSSINK = 0.
		TotSSINK2 = 0.0
		Iflag = 0
		MMUpN = 0.0
   
C 
C  Reading of the input files and initial calculations 
C
		im = 430
		il = 0
		Open(40,file = SoluteNitrogenFile, status = 'old',ERR = 10)   
		im = im + 1
		il = il + 1
		Read(40, * ,ERR = 10)
		im = im + 1
		il = il + 1
		Read(40, * ,ERR = 10)
		im = im + 1
		il = il + 1
		Read(40, * ,ERR = 10) Isink,Rroot
		im = im + 1
		il = il + 1
		Read(40, * ,ERR = 10)
		im = im + 1
		il = il + 1
		Do M = 1,2
		 Read(40, * ,ERR = 10) ConstI(m),constK(m), Cmin0(M)
		 il = il + 1
		enddo
		Read(40, * ,ERR = 10)
		im = im + 1
		il = il + 1
		Do M = 1,NMat
		 Read(40, * ,ERR = 10) DlngR(M)
		 il = il + 1
		enddo
		close(40)
	  Open(101, file = 'tempC.txt')
		write(101, * )      UptakeN(Isink)
		write(101,103)'    e   Csink       CsinkN      Croot       Croot  
     $  CMean' 
		write(101,103)'           passive     Total    old     young    conc'
103		format(a) 
      Endif
      if(ITIME.EQ.1) iflag =  1
      SIncrSink =  0.
      WIncrsink = 0.
      SincrSink2 = 0.0


     
c new mass uptake	
C
C calculate sum of N in all forms in soil using triangluar elements
C
C units of conc ug/liter 
	 Do n = 1,NumEl
           cSink(n,1) = 0.
           NUS = 4
           DmE=0
           if(KX(n,3).eq.KX(n,4)) NUS = 3
*          Loop on subelements
           AE = 0.0
           M = MatNumE(n)
		 do k = 1,NUS - 2
             i = KX(n,1)
             j = KX(n,k  +  1)
             l = KX(n,k  +  2)        
		     Ci(1) = x(l) - x(j)
             Ci(2) = x(i) - x(l)
             Ci(3) = x(j) - x(i)
             Bi(1) = y(j) - y(l)
             Bi(2) = y(l) - y(i)
             Bi(3) = y(i) - y(j)
             AE = AE  +  (Ci(3) * Bi(2) - Ci(2) * Bi(3)) / 2.
             CTot=Conc(i,1)+ Conc(j,1)+Conc(l,1)
             CSink(n,1) = CSink(n,1)  +  (Conc(i,1)  +  
     &             Conc(j,1)  +  Conc(l,1))/3.
    
            
		   DmE = DmE+  Dmol(1) * (ThNew(i) * Tau(hNew(i)) + 
     !	   ThNew(j) * Tau(hNew(j)) + ThNew(l) * Tau(hNew(l))) / 3.0     		       			
 		Enddo
      	Csink_M(n) = CSink(n,1)   ! mean of concentration for element              
		CSink(n,1) = CSink(n,1) * Sink(n)       
          SIncrSink2 = SIncrSink2 + cSink(n,1) * step * 14. / 62. * AE
		Disp(n,1) = DmE/(NUS-2.0) + DlngR(M) * VUP(n,1)
		Disp(n,2) = DmE/(NUS-2.0) + DlngR(M) * VUP(n,2)		
	 Enddo	  
C calculate nutrigen root uptake    
	 Call massRootflux
c -  -  -  -  -  -  -  -   -  -  -  -  -  -  -  -  -  -  -     
       do n = 1,NumEl
         
         SIncrSink = SIncrSink + CSink(n,1) * Step* Area(n) *  14. / 62.
        
         WIncrSink = WIncrSink + Sink(n) * Area(n) * Step
       enddo
       TotWSink = TotWSink + WIncrSink   ! water sink
       TotSSINK = TotSSINK + SIncrSink   ! solute sink based on csink after dispersive calcs
       TotSSINK2 = TOTSSink2 + SIncrSink2 ! solute sink based on the active uptake before calling sub
       If(ITIME.eq.24.and.iflag.eq.1) then
         Write(94, * ) Jday,TotWSink,TotSSink,TotSSInk2
         iflag = 0
       Endif
       ! for debugging
       !SIncrSink=SIncrSink2   ! set to old method for test
      return
10     Call errmes(im,il)
      end
     
      subroutine massRootflux
      include 'public.ins'
      IncLude 'Puplant.ins'
      Include 'Puweath.ins'
      Common / Mass_root_M / Isink,Rroot,ConstI(2),constK(2),Cmin0(2),
     1    DlngR(NMatD),
     1	Disp(NumElD,2),RootExchange(2),CSink_M(NumElD),Cr_M(NumElD,2) 
	Common  / PotNitr /  PotNitrogen_t
	real * 8 qsinkC(NumElD,2),partC(NumElD,2), partCr(NumElD,2)
      integer e,Iroot(NumElD,2)
	Real *8 determ,b,a,c,bet,bet2,bet2a,alf,alf2,part1,part2,cr,betR
	  
	 pi = 3.14
	 pi2 = 6.28
	 pi4 = 12.56
	 Rroot2 = Rroot * Rroot
	 alpha2=0.05
C FUP is root density, VUP is water uptake per unit root length
	 do j = 1,2
	  do e = 1, NumEl
         Iroot(e,j) = 0
	   qsinkC(e,j) = 0.0 
	   if(FUP(e,j).gt.0) THEN
	    alf = VUP(e,j) / (pi2 * alpha2*Disp(e,j))  ! /2 pi D? this is a
	    if (Disp(e,j).le.0.0) alf=0
	    SC = 0.017  -  (PSIS(e) * 0.5)  
	    betR = SC / Rroot 
	    bet = 1.0 / sqrt(pi * (FUP(e,1) + FUP(e,2)))
	    bet=bet/Rroot 
	    bet = amin1(betR,bet)
	    bet2 = bet * bet
          Iroot(e,j) = 1  
           if(alf.le.0.0) then
	      part1 = pi4 * Disp(e,j)/((bet2 * dlog(bet2) / 
     !		  (bet2 - 1.0d0) - 1.0d0)) 
		  partC(e,j) = part1 * Csink_M(e)
		  partCr(e,j) = part1
c		  qc(e) = pi4 * D(e) * (C(e) - CR(e) * (bet2 * ln(bet2) / (bet2 - 1.0) - 1)
	     else    if (alf.eq.2) then
	      part1 = (bet2 - 1.0d0)
		  part2 = dlog(bet2)
		  partC(e,j) = VUP(e,j) * part1 * Csink_M(e) / (part1 + part2) ! was + part1 part2
		  partCr(e,j) = VUP(e,j) * part2 / (part1 + part2)
c		  qc(e) = qw(e) * (part1 * C(e) - part2 * CR(e)) / 
c     !		    (part1 + part2)
	     else
		  alf2 = alf / 2.0d0
		  bet2a = bet2 * bet ** ( - alf)
		  part1 = (bet2 - 1.0d0) * (1.0d0 - alf2)
		  part2 = (bet2a - 1.0d0)
		  partC(e,j) = VUP(e,j) * part1 * Csink_M(e) / (part1 + part2)  !ug cm-3 cm root day-1 
		  partCr(e,j) =  VUP(e,j) * part2 / (part1 + part2)             !ug cm-3 cm root day-1 
c		  qc(e) = qw(e) * ((part1 * C(e) - part2 * CR(e)) / 
c     !	        (part1 + part2)  
           endif       	 
           if (Isink.eq.1.or.Isink.eq.2) then !Isink = 1 passive mass uptake (Mihaelson - Menton kinetics)
            a = partCr(e,j)                  
            if (Isink.eq.2) a = a +  VUP(e,j)!Isink = 2 passive mass uptake (Mihaelson - Menton kinetics)  + convetcive flux
	      c = - (partC(e,j) - a * Cmin0(j))
            b = ConstI(j) + c + a * constK(j)
	      c = c * constK(j)
	      determ = b * b - 4.0d0 * a * c
	      if (determ.ge.0.0) then
	       determ = dsqrt(determ)
	       Cr =  ( - b + determ) / (2.0d0 * a)
c	       qsinkC(e,j) =  partC(e,j) - partCr(e,j) * Cr
	       if(Cr.le.0.0d0) then
	        a = a
	       endif
	      else
	       Cr = 0.0
            endif
	      Cr = dmax1(Cr,0.0d0) + Cmin0(j)
	      qsinkC(e,j) = partC(e,j) - partCr(e,j) * Cr
	      qsinkC(e,j) = amax1(qsinkC(e,j),0.0)
            Cr_M(e,j) = Cr + Cmin0(j) 
	    else   if (Isink.eq.3) then ! active mass uptake Cr = Croot
		  AC = AC + partC(e,j) * FUP(e,j) * Area(e)
		  ACR = ACR + partCr(e,j) * Fup(e,j) * Area(e)
          else   if (Isink.eq.4) then ! active mass uptake rootflux = Convective + Despersive  + Mihael - Menson flux
	      AC = AC + partC(e,j) -  (VUP(e,j) 
     !                  + RootExchange(j)) * Cmin0(j)
	      partCr(e,j) = partCr(e,j) + (VUP(e,j) + RootExchange(j)) 
	    else   if (Isink.eq.5) then ! active mass uptake rootflux = Despersive  + Mihael - Menson flux
	      AC = AC + partC(e,j) -  RootExchange(j) * Cmin0(j)
	      partCr(e,j) = partCr(e,j) + RootExchange(j)
	    else   if (Isink.eq.6) then ! active mass uptake rootflux = Despersive  + Convective)
	      AC = AC + partC(e,j) * FUP(e,j)
	      ACR = ACR + partCr(e,j) * Fup(e,j)
	    endif
         endif 
	  enddo
	 enddo 	
	 if (isink.eq.3) then
	  Cr1 =  -(PotNitrogen_t - AC) / ACR   ! PotNitrogen_t should be mg day
	  Cr1 = amax1(Cr1,0.0)
10	  Cr = Cr1        
	  do j = 1,2
	   do e = 1, NumEl
          if (Iroot(e,j).eq.1) THEN
	      qsinkC(e,j) =  partC(e,j) - partCr(e,j) * Cr		   
	      AC = AC + partC(e,j) * FUP(e,j) * Area(e)
	      ACR = ACR + partCr(e,j) * Fup(e,j) * Area(e)
	    endif
	   enddo
	  enddo 
	  Cr1 = (PotNitrogen_t - AC) / ACR
	  Cr1 = amax1(Cr1,0.0)
	  if (abs(Cr - Cr1).GT.0.01 * Cr)then
	   AC = 0.0
	   ACR = 0.0
	   go to 10	
	  endif 
	  do e = 1, NumEl
         do j = 1,2
	    qsinkC(e,j) = partC(e,j) - partCr(e,j) * Cr	
	    qsinkC(e,j) = amax1(qsinkC(e,j),0.0)
	    Cr_M(e,j) = Cr
	   enddo
	  enddo
	 endif   
	 do e = 1, NumEl
	  sumSinkN = 0.0
	  do j = 1,2
	   sumSinkN = sumSinkN + qsinkC(e,j) * Fup(e,j)
	  enddo	  
	  If(CSink(e,1).gt.0 )then
	    write(101,102)time, e, CSink(e,1),sumSinkN,
     !	 (Cr_M(e,j), j = 1,2) , Csink_m(e) 
          continue
	  endif
	  CSink(e,1) = amax1(sumSinkN,0.0)
	 enddo
102    format(g14.8,1x,I5,5(2x,e10.3))
	return
	end
      Subroutine Matric(N,A,B,C,F)
	 real A( * ),B( * ),F( * ),C( * ),ALF(100)
	 N1 = N - 1
       R = 0.0
	 DO I = 1,N1
	  ALF(I) = B(I) / A(I)
	  R = R + ALF(I) * C(I)
	 ENDDO
	 CR = F(N) / R
	 DO I = 1,N1
        A(I) = F(I) / (ALF(I) + A(I))
	 ENDDO
	 RETURN
	END 
	


