*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
*
*     Chemical application module - example 1 from section 14.1  * 
*
*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      Subroutine Mngm()
      Include 'public.ins'
      Common /Mngmnt/ tAppl, cAppl,NumAp,nAppl(NumBPD),ModNum
      If (lInput.eq.1) then
        Open(40,file='Mngm.dat')
        Read(40,*)
        Read(40,*)
        Read(40,*,ERR=20) tAppl, cAppl,NumAp
        Read(40,*)
        Read(40,*,ERR=20) (nAppl(i),i=1,NumAp)
        Close(40)
        NumMod=NumMod+1
        ModNum=NumMod
        tNext(ModNum)=tAppl
      Endif
C
      If(Abs(Time-tNext(ModNum)).lt.0.1*Step) then
        Do i=1,NumAp
          Conc(nAppl(i),1)=cAppl/thNew(nAppl(i))
        Enddo
        tNext(ModNum)=1.E+32
      Endif
C
      Return
20    Stop 'Mngm data error'
      End



