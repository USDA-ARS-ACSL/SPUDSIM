cdt 7/9/07 added code to calculate loss of N in flux out of the domain 
Cdt N in 2DSOIL is nitrate (NO3) and the units are ug per cm3 or milligram/liter
       Subroutine Nitrogen_Mass_Balance()
        include 'public.ins'
        include 'nitvar.ins'
        include 'puplant.ins'
        Dimension Bi(3),Ci(3)
        Character*10 Date
        Real AE,Profile_N,Manure_N,Litter_N,
     &    Min_N,Ammon_N,Org_N,Denitr,
     &    NO3_Mean,NhumusMean,NLitterMean,NManureMean,NAmmoniaMean,All_N
        Integer ModNum
C Variables to hold mulch C and N totals
        Real humusC, litterC, manureC,
     &        CHumusMean, CLitterMean, CManureMean
        
        common /N_BAL/ModNum,CFlux,CFluxPrevious

      If (lInput.eq.1) then
        open(91,file=MassBalanceFileOut,status='unknown',recl=270)
        write(91,'(14A12)') 'time,','Date,','Min_N,','Humus_N,',
     !     'Humus_C,','Manure_N,','Manure_C,','Litter_N,',
     !     'Litter_C,','Ammon_N,','All_N,','water,','Denitr,',
     !     'CFlux'
 
        NumMod=NumMod+1
	ModNum=NumMod
	tNext(ModNum)=time
	CFlux=0
	CFluxPrevious=0
 
       Endif
          do i=1,NumBp
	     n=KXB(i)
C Cflux is loss of N in mg
	     if ((CodeW(n).eq.(-7)).or.(CodeW(n).eq.(2))) then
               
cccz	        Cflux=Cflux+Q(n)*conc(n,1)*step*14/62
cccz using a starndard conversion factor "NO3mass_2_Nmass"
cccz why "mg"??????????????????????????????????????????????????????
              Cflux=Cflux+Q(n)*conc(n,1)*step*NO3mass_2_Nmass 
           endif
C for case of downward drainage and a constrant BC 
C (does not account for upward flow yet)
           if ((CodeW(n).eq.(1))) then
             if (Q(n).lt.0) then
cccz                 Cflux=Cflux+Q(n)*conc(n,1)*step*14/62 ! only consider chemical out- flux>0
cccz using a starndard conversion factor "NO3mass_2_Nmass"
cccz why "mg"??????????????????????????????????????????????????????
                 Cflux=Cflux+Q(n)*conc(n,1)*step*NO3mass_2_Nmass  ! only consider chemical out- flux>0
             endif
           endif
	  EndDo           
           
        t=time
        if (Abs(time-tNext(ModNum)).lt.0.001*Step.or.lInput.ne.0) then
           tNext(ModNum)=tNext(ModNum)+1.0
           Profile_N=0.0
           Min_N=0.0
           Org_N=0.0
           Litter_N=0.0
           Manure_N=0.0
           Ammon_N=0.0
           W_Sum=0.0
           Denitr=0.0
           humusC=0.0
           litterC=0.0
           manureC=0.0
           
           
	   Sum=0.
           W_Sum=0.
	   Do n=1,NumEl
             NUS=4
             if(KX(n,3).eq.KX(n,4)) NUS=3
             Sum1=0.
             NDenitrifyMean=0.
             NO3_Mean=0.0
             NhumusMean=0.0
             NLitterMean=0.0
             NManureMean=0.0
             NAmmoniaMean=0.0
             CHumusMean=0.0
             CLitterMean=0.0
             CManureMean=0.0
  
*         Loop on subelements
             do k=1,NUS-2
               i=KX(n,1)
               j=KX(n,k+1)
               l=KX(n,k+2)
               Ci(1)=x(l)-x(j)
               Ci(2)=x(i)-x(l)
               Ci(3)=x(j)-x(i)
               Bi(1)=y(j)-y(l)
               Bi(2)=y(l)-y(i)
               Bi(3)=y(i)-y(j)
               AE=(Ci(3)*Bi(2)-Ci(2)*Bi(3))/2.
cccz Conc is NO3 as ug/cm3
               NO3_Mean=NO3_Mean
     &             +AE*(Conc(i,1)*ThNew(i)
     &                 +Conc(j,1)*ThNew(j)
     &                 +Conc(l,1)*ThNew(l))
     &                 /3.*NO3mass_2_Nmass  ! convert to N

               Sum1=Sum1+AE*(ThNew(i)+ThNew(j)+ThNew(l))/3.

cccz need to present all the N mass based using soil bulk density   
c  organic N is in ug/g of soil, use bd to get ug/cm3 of soil in order to sum over volume
               NhumusMean=NhumusMean+AE*(Nh(i)*Blkdn(MatNumN(i))
     &             +Nh(j)*Blkdn(MatNumN(j))
     &             +Nh(l)*Blkdn(MatNumN(l)))/3.
               NLitterMean=NLitterMean+AE*(NL(i)*
     &                 Blkdn(MatNumN(i))
     &             +NL(j)*Blkdn(MatNumN(j))
     &             +NL(l)*Blkdn(MatNumN(l)))/3.
               NManureMean=NManureMean+AE*(Nm(i)*Blkdn(MatNumN(i))
     &             +Nm(j)*Blkdn(MatNumN(j))
     &             +Nm(l)*Blkdn(MatNumN(l)))/3.
     
               ChumusMean=ChumusMean+AE*(Ch(i)*Blkdn(MatNumN(i))
     &             +Ch(j)*Blkdn(MatNumN(j))
     &             +Ch(l)*Blkdn(MatNumN(l)))/3.
               CLitterMean=CLitterMean+AE*(CL(i)*
     &                 Blkdn(MatNumN(i))
     &             +CL(j)*Blkdn(MatNumN(j))
     &             +CL(l)*Blkdn(MatNumN(l)))/3.
               CManureMean=CManureMean+AE*(Cm(i)*Blkdn(MatNumN(i))
     &             +Cm(j)*Blkdn(MatNumN(j))
     &             +Cm(l)*Blkdn(MatNumN(l)))/3.

               NAmmoniaMean=NAmmoniaMean+AE*(NNH4(i)*
     &               Blkdn(MatNumN(i))
     &             +NNH4(j)*Blkdn(MatNumN(j))
     &             +NNH4(l)*Blkdn(MatNumN(l)))/3.
               NDenitrifyMean=NDenitrifyMean+AE*(Denit(i)*
     &                    Blkdn(MatNumN(i))
     &             +denit(j)*Blkdn(MatNumN(j))
     &             +denit(l)*Blkdn(MatNumN(l)))/3.             
               Enddo
             Min_N=Min_N+NO3_Mean
             Org_N=Org_N+NhumusMean
             Litter_N=Litter_N+NLitterMean
             Manure_N=Manure_N+NManureMean
             humusC=humusC + CHumusMean
             litterC=litterC + CLitterMean
             manureC=manureC + CManureMean
             Ammon_N=Ammon_N+NAmmoniaMean
             W_Sum=W_Sum+Sum1
             Denitr=Denitr+NDenitrifyMean
             Enddo

             
             
C calculate sum of N in all forms in soil
C factor is appropriate for mg/cm3 Total N is mg per slab (grid width x 1cm)
C Mineral ad organic N is ug/cm3. Do a summation over the simulation domain - total ug in a plant slab
C Fact now works while there is a plant 

          fact=1/(0.01*RowSp/100.*EOMult) !m2 of slab
          fact=fact*10000/1000/1000/1000    !m2 of slab ->ha, ug-->Kg
          All_N=(Min_N+Org_N+Litter_N+Manure_N
     !           +Ammon_N)*fact
         iday=int(t)
	   call caldat(iday,mm,id,iyyy) 
          write (date,'(i2.2,A1,i2.2,A1,i4.4)') mm,'/',id,'/',iyyy  
          write(91,10) 
     !       time,date,Min_N*fact,
     !       Org_N*fact,humusC*fact,Manure_N*fact,
     !       manureC*fact, Litter_N*fact,      
     !       litterC*fact,Ammon_N*fact,
     !       All_N,W_Sum,Denitr*fact,CFLux*fact
	   
          CFluxPrevious=CFlux
        endif
10    Format (1F12.4,',', A12,',', 11(F12.4, ','),F12.4)
      return
      end


