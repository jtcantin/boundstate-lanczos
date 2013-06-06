c-----------------------------------------------------------------------
      subroutine calcnb(Jbig,jmax,ik,ip,numbas)
      implicit double precision(a-h,o-z)
c ... the loop to calculate the number of basis functions
      
      numbas=0
      do jsmall=0,jmax
c ...   judge the starting and end point for ksmall
        ksmmax=kstcal(jsmall,ik)
        do ksmall=-ksmmax,ksmmax,2
          Kbigst=0
          if(mod(Jbig+ip,2).eq.1.and.ksmall.eq.0)Kbigst=1
          Kbiged=min(jsmall,Jbig)
c         write(6,*)Kbigst,Kbiged
          do Kbig=Kbigst,Kbiged
c ...       loop 1: the dumb if statement below is to avoid the situation
c ...       of K=0 and k<0, which should have been considered before coding
c ...       but not.
            if(.not.(Kbig.eq.0.and.ksmall.lt.0)) then
              numbas=numbas+1
              write(6,*)jsmall,ksmall,Kbig,numbas
            endif
          enddo
        enddo
      enddo

      return
      end

c-----------------------------------------------------------------------
      integer function kstcal(jsmall,ik)
      implicit double precision(a-h,o-z)
c ... given jsmall and ik, calculate the starting value of ksmall
      if((mod(jsmall,2).eq.0.and.ik.eq.0).or.
     +   (mod(jsmall,2).eq.1.and.ik.eq.1)) then
        kstcal=jsmall
      else
        kstcal=jsmall-1
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine prepare(numbas,icode,maxj,jmax,Jbig,ik,ip,Kbgind,
     +                   coriol,rotor,Ah2o,Bh2o,Ch2o,h2_2mu,ijcori,
     +                   ncorio,ijrot,nrotor,maxfac,fact,Tcos,Tsin,
     +                   nlgrid,ncgrid,wgtgl,glgrid,gcgrid,nchi,nthe,
     +                   nrad,vpes,gchi,gthe,grad,nrpont,rsmall,rlarge,
     +                   rgrid,nrgrid,vmat,rkin,nKbig,ib000,ib100,
     +                   ib101,ib110,ib111,jb111,potfil)
      implicit double precision(a-h,o-z)

      dimension icode(numbas),Kbgind(Jbig+2),coriol(numbas*4),
     +          rotor(numbas*4),ijcori(numbas*4,2),ijrot(numbas*4,2),
     +          fact(0:maxfac),Tcos(nlgrid*ncgrid,numbas),
     +          Tsin(nlgrid*ncgrid,numbas),wgtgl(nlgrid),glgrid(nlgrid),
     +          gcgrid(ncgrid),vpes(-4:nchi+5,-4:nthe+5,nrad),
     +          gchi(-4:nchi+5),gthe(-4:nthe+5),grad(nrad),
     +          rgrid(nrpont),vmat(nrpont*nlgrid*ncgrid),
     +          rkin(nrpont*(nrpont+1)/2)
      parameter(eps=1.0d-14,one=1.d0,zero=0.d0)
      parameter(pi=3.14159265358979323846d+00)
      logical safe
      character potfil*30

c ... truncate the potfil to the correct length
      do ich=2,29
        if(potfil(ich-1:ich+1).eq.'pot')lenfil=ich+1
      enddo
c     write(6,'(a)')potfil(1:lenfil)
c ... safety check
      if(jmax.gt.maxj.or.Jbig.ge.maxj) stop'too big jmax or Jbig'

      call calfac(fact,maxfac)

      write(6,9000)

c ... the loop to generate code for potential matrix elements
      safe=.true.

      ib=0
      do jsmall=0,jmax
c ...   judge the starting and end point for ksmall
        ksmmax=kstcal(jsmall,ik)
        do ksmall=-ksmmax,ksmmax,2
          Kbigst=0
          if(mod(Jbig+ip,2).eq.1.and.ksmall.eq.0)Kbigst=1
          Kbiged=min(jsmall,Jbig)
c         write(6,*)Kbigst,Kbiged
          do Kbig=Kbigst,Kbiged
c ...       loop 2: the dumb if statement below is to avoid the situation
c ...       of K=0 and k<0, which should have been considered before coding
c ...       but not.
            if(.not.(ksmall.lt.0.and.Kbig.eq.0))then
            ib=ib+1
            icode(ib)=
     +      Kbig*100000000+jsmall*1000000+iabs(ksmall)*10000+ib*10
            if(ksmall.lt.0)icode(ib)=icode(ib)+1
c############## safety check for the code###############################
            Kbig2=mod(icode(ib)/100000000,100)
            jsmal2=mod(icode(ib)/1000000,100)
            ksmal2=mod(icode(ib)/10000,100)
            jb=mod(icode(ib)/10,1000)
            ksign=mod(icode(ib),10)
            ksmal2=((-1)**ksign)*ksmal2
c           write(6,*)ib,jsmall,ksmall,Kbig,jb,jsmal2,ksmal2,Kbig2,icode
            if(ib.ne.jb.or.jsmall.ne.jsmal2.or.ksmall.ne.ksmal2.or.
     +         Kbig.ne.Kbig2)then
              write(6,*)ib,jsmall,ksmall,Kbig,jb,jsmall2,ksmall2,Kbig2
              safe=.false.
            endif
c#################### end of safety check###############################
c ...       get the basis position for the eigen state analysis
            if(jsmall.eq.0.)ib000=ib
            if(jsmall.eq.1) then
              if(ksmall.eq.0.and.Kbig.eq.0)ib100=ib
              if(ksmall.eq.0.and.Kbig.eq.1)ib101=ib
              if(ksmall.eq.1.and.Kbig.eq.0)ib110=ib
              if(ksmall.eq.1.and.Kbig.eq.1)ib111=ib
              if(ksmall.eq.-1.and.Kbig.eq.1)jb111=ib
            endif
            endif
          enddo
        enddo
      enddo
      if(.not.safe)stop 'wrong in generating icode'
      if(ib.ne.numbas)then
        write(6,*)'wrong counting in generating icode',ib,numbas
        stop
      endif

      call bubble_sort(icode,numbas)

c ... group the Kbig
c ... just for printing
      write(6,9001)
      i=1
      jsmall=mod(icode(i)/1000000,100)
      ksmall=mod(icode(i)/10000,100)
      Kbig=mod(icode(i)/100000000,100)
      ib=mod(icode(i)/10,1000)
      ksign=mod(icode(i),10)
      ksmall=((-1)**ksign)*ksmall
      write(6,9002)i,icode(i),Kbig,jsmall,ksmall,ib
      nKbig=1
      Kbgind(1)=1
      do i=2,numbas
        jsmall=mod(icode(i)/1000000,100)
        ksmall=mod(icode(i)/10000,100)
        Kbig=mod(icode(i)/100000000,100)
        Kbg2=mod(icode(i-1)/100000000,100)
        ib=mod(icode(i)/10,1000)
        ksign=mod(icode(i),10)
        ksmall=((-1)**ksign)*ksmall
        if(Kbig.ne.Kbg2) then
c         write(6,*)'break'
          nKbig=nKbig+1
          Kbgind(nKbig)=i
        endif
        write(6,9002)i,icode(i),Kbig,jsmall,ksmall,ib
      enddo
c ... safety check
      if(nKbig-1.ne.Jbig) then
        write(6,9003)'nKbig-1 != Jbig',nKbig-1,Jbig
        stop
      endif
      Kbgind(nKbig+1)=numbas+1
      write(6,9003)nKbig,(Kbgind(i),i=1,nKbig)

c ... prepare Gauss quadrature grids
      call gaulegf(-one,one,glgrid,wgtgl,nlgrid,nlgrid)
c     call gaucheb(-one,one,gcgrid,wgtgc,ncgrid,ncgrid)
      call fftgrd(gcgrid,ncgrid,ncgrid)

c ... get the cosine and sine transformation matrix
      call transm(numbas,nKbig,icode,Kbgind,
     +            wgtgl,nlgrid,glgrid,ncgrid,gcgrid,Tsin,Tcos,
     +            maxfac,fact)

c ... test orthogonality of Tcos and Tsin
c     numlc=ncgrid*nlgrid
c     write(6,*)'punch in two integers between 1 and ',numbas
c     read(5,*)ib,jb
c     summ1=zero
c     summ2=zero
c     do lc=1,numlc
c       summ1=summ1+Tcos(lc,ib)*Tcos(lc,jb)
c       summ2=summ2+Tsin(lc,ib)*Tsin(lc,jb)
c     enddo
c     write(6,*)summ1,summ2

c ... read the potential energy surface on grid
      call potread(nchi,nthe,nrad,vpes,gchi,gthe,grad,potfil,lenfil)

c ... calculate step for radial grid
      rstep=(rlarge-rsmall)/dfloat(nrgrid)
c ... loop over radial grid to generate the grand potential energy matrix
      do ir=1,nrpont
        rinp=rsmall+rstep*dfloat(ir)
        write(6,*)ir,rinp
        call prev(rinp,vmat,nlgrid,ncgrid,glgrid,
     +            gcgrid,nchi,nthe,nrad,vpes,gchi,gthe,grad,ir,nrpont)
c       numlc=ncgrid*nlgrid
c       write(6,*)vmat((ir-1)*numlc+1)
      enddo

c ... calculate the r grids
      rgrid(1)=rsmall+rstep
      do ir=2,nrpont
        rgrid(ir)=rgrid(ir-1)+rstep
c       write(6,*)ir,rgrid(ir)
      enddo

c ... loop over to generate coriolis and H2O rotor matrix elements
c ... this is a stupid full loop.  Toby doesn't want to waste time
c ... on it.  Anyway, the loop is just called once.
c ... 100,110,120 form outer loop
      ncount=0
      ncorio=0
      nrotor=0
      irot=0
      ib=0
      do 100 jsmal1=0,jmax
      kmax1=kstcal(jsmal1,ik)
      do 110 ksmal1=-kmax1,kmax1,2
      Kst1=0
      if(mod(Jbig+ip,2).eq.1.and.ksmal1.eq.0)Kst1=1
      Ked1=min(jsmal1,Jbig)
      do 120 Kbig1=Kst1,Ked1
c ... loop 3: the dumb if statement below is to avoid the situation
c ... of K=0 and k<0, which should have been considered before coding
c ... but not. The if statement is closed before 120 continue
      if(.not.(ksmal1.lt.0.and.Kbig1.eq.0)) then
      ib=ib+1

c ... 105,115,125 form inner loop
      jb=0
      do 105 jsmal2=0,jmax
      kmax2=kstcal(jsmal2,ik)
      do 115 ksmal2=-kmax2,kmax2,2
      Kst2=0
      if(mod(Jbig+ip,2).eq.1.and.ksmal2.eq.0)Kst2=1
      Ked2=min(jsmal2,Jbig)
      do 125 Kbig2=Kst2,Ked2
c ... loop 4: the dumb if statement below is to avoid the situation
c ... of K=0 and k<0, which should have been considered before coding
c ... but not. The if statement is closed before 125 continue
      if(.not.(ksmal2.lt.0.and.Kbig2.eq.0))then
      jb=jb+1
c     write(6,*)jsmal1,ksmal1,Kbig1,ib,';',jsmal2,ksmal2,Kbig2,jb
      ncount=ncount+1
c ... calculate and store Coriolis matrix elements
      elemnt=corcal(jsmal1,ksmal1,Kbig1,jsmal2,ksmal2,Kbig2,Jbig,ip)
      if(abs(elemnt).gt.eps) then
        ncorio=ncorio+1
        coriol(ncorio)=elemnt*h2_2mu
        ijcori(ncorio,1)=jb
        ijcori(ncorio,2)=ib
        write(6,*)jb,jsmal2,ksmal2,Kbig2,ib,jsmal1,ksmal1,Kbig1,
     +            elemnt,ncorio,' Coriolis'
      endif
c ... calculate and store the H2O rotor matrix elements
      elemnt=rotcal(jsmal1,ksmal1,Kbig1,jsmal2,ksmal2,Kbig2,ip,Ah2o,
     +              Bh2o,Ch2o,Jbig)
      if(abs(elemnt).gt.eps)then
        nrotor=nrotor+1
        rotor(nrotor)=elemnt
        ijrot(nrotor,1)=jb
        ijrot(nrotor,2)=ib
        write(6,*)jb,jsmal2,ksmal2,Kbig2,ib,jsmal1,ksmal1,Kbig1,
     +            elemnt,nrotor,' Rotor'
      endif

      endif
  125 continue
  115 continue
  105 continue
      endif
  120 continue
  110 continue
  100 continue
c ... safety check
      if(ncount.ne.numbas*numbas) then
        write(6,*)'wrong count in Coriolis and rotor loop',ncount,
     +             numbas*numbas
        stop
      endif
      if(ncorio.gt.numbas*4) then
        write(6,*)'too many coriolis matrix elements',ncorio,numbas*4
        stop
      endif
      if(nrotor.gt.numbas*4) then
        write(6,*)'too many rotor matrix elements',nrotor,numbas*4
        stop
      endif

      write(6,9004)ncorio
      write(6,9005)nrotor

c ... loop over to generate all the radial kinetic energy matrix element
c ... and store in a row-wise left lower triangle matrix rkin
c ... for such a half matrix, the element ij is stored in the position
c ... n=i*(i-1)/2+j of the vector

      write(6,9006)

      numkin=nrpont*(nrpont+1)/2
c ... calculate some common factors
      cofac1=h2_2mu*pi*pi/(2.d0*(rlarge-rsmall)*(rlarge-rsmall))
      cofac2=dfloat(2*nrgrid*nrgrid+1)/3.d0
      pi_N=pi/dfloat(nrgrid)
      pi_2N=pi/dfloat(2*nrgrid)

      nelmnt=0
      do i=1,nrpont
        do j=1,i-1
          nelmnt=nelmnt+1
          phase=dfloat((-1)**(i-j))
          fac1=one/(sin(pi_2N*dfloat(i-j))*sin(pi_2N*dfloat(i-j)))
          fac2=one/(sin(pi_2N*dfloat(i+j))*sin(pi_2N*dfloat(i+j)))
          rkin(nelmnt)=cofac1*phase*(fac1-fac2)
c         write(6,*)cofac1*phase*(fac1-fac2),rkin(nelmnt)
        enddo
        nelmnt=nelmnt+1
        fac2=one/(sin(pi_N*i)*sin(pi_N*i))
        rkin(nelmnt)=cofac1*(cofac2-fac2)
      enddo
c     if(nelmnt.ne.numkin)then
c       write(6,*)'numkin != nelmnt',numkin,nelmnt
c       stop
c     endif

c ... printing the triangular radial kinetic matrix
c     call prtril(rkin,nrpont)

      write(6,9007)

 9000 format(/'PREPARING QUANTUM # CODE OF,K,j,k,ib,ksign')
 9001 format (/'ICOD     ICODE   K   j   k  IB')
 9002 format (i4,i10,4(i4))
 9003 format (/'NO. OF KBIG VALUES:',I4,3X,'STARTING INDICES:',11(I4))
 9004 format(/'NO. OF CORIOLIS MATRIX ELEMENTS',I5)
 9005 format(/'NO. OF H2O ROTOR MATRIX ELEMENTS',I5)
 9006 format (/'CALCULATING RADIAL KINETCI ENERGY MATRIX')
 9007 format(/'PREPARE FINISHED')
      return
      end
c-----------------------------------------------------------------------
c     Subroutine bubble_sort
c       this routine sorts the given data
c
      subroutine bubble_sort(data,count)
c     
c     argument:  count is a positive integer
      integer count
c     argument:  data is an array of size count
      integer data(count)
c
c     local variables:
      integer i
c       how many times we have passed through the array
      integer pass
c       flag variable: 1 if sorted; 0 if not sorted  
      integer sorted
c       temporary variable used for swapping       
      integer temp

      pass = 1
 1    continue
      sorted = 1
      do 2 i = 1,count-pass
        if(data(i) .gt. data(i+1)) then
          temp = data(i)
          data(i) = data(i+1)
          data(i+1) = temp
          sorted = 0
        endif
 2    continue
      pass = pass +1
      if(sorted .eq. 0) goto 1
      return
      end
c-----------------------------------------------------------------------
      double precision function corcal(jsmal1,ksmal1,Kbig1,jsmal2,
     +                                 ksmal2,Kbig2,Jbig,ip)
      implicit double precision(a-h,o-z)

      corcal=0.d0
      if(jsmal1.ne.jsmal2)return

      pre1=dfloat((1+kdel(Kbig2,0)*kdel(ksmal2,0))*
     +            (1+kdel(Kbig1,0)*kdel(ksmal1,0)))
      pre1=1.d0/sqrt(pre1)
      ifac1=(Jbig*(Jbig+1)+jsmal1*(jsmal1+1)-2*Kbig1*Kbig1)*
     +      kdel(Kbig1,Kbig2)
      ifac1=ifac1*(kdel(ksmal2,ksmal1)+((-1)**(Jbig+ksmal1+ip))*
     +      kdel(-ksmal2,ksmal1)*kdel(Kbig1,0))
      fac1=dfloat(ifac1)

      ifac2=(kdel(ksmal2,ksmal1)+((-1)**(Jbig+ksmal1+ip))*
     +      kdel(-ksmal2,ksmal1)*kdel(Kbig1,0))*kdel(Kbig2,Kbig1+1)
      fac2=dfloat(ifac2)
      fac2=fac2*cplus(Jbig,Kbig1)*cplus(jsmal1,Kbig1)

      ifac3=(kdel(ksmal2,ksmal1)+((-1)**(Jbig+ksmal1+ip))*
     +      kdel(-ksmal2,ksmal1)*kdel(Kbig1,1))*kdel(Kbig2,Kbig1-1)
      fac3=dfloat(ifac3)
      fac3=fac3*cminus(Jbig,Kbig1)*cminus(jsmal1,Kbig1)

      corcal=pre1*(fac1-fac2-fac3)

      return
      end
c------------------------------------------------------------------
      integer function kdel(i,j)
      implicit double precision(a-h,o-z)

c ... shift-up for the situation of i=j=0
c ... because of the shifting, the delta function
c ... is ILL-DEFINED for the case of i=j=-1000
      ii=i+1000
      jj=j+1000
      kdel=((ii+jj)-iabs(ii-jj))/((ii+jj)+iabs(ii-jj))

      return
      end
c------------------------------------------------------------------

      double precision function cplus(j,k)
      implicit double precision(a-h,o-z)

      parameter(one=1.0d+00,zero=0.0d+00)

c ... in case k runs out of the range
      if(k.ge.j.or.k.lt.(-j)) then
        cplus=zero
        return
      endif

      dj=dfloat(j)
      dk=dfloat(k)

      cplus=sqrt(dj*(dj+one)-dk*(dk+one))

      return
      end
c------------------------------------------------------------------

      double precision function cminus(j,k)
      implicit double precision(a-h,o-z)

      parameter(one=1.0d+00,zero=0.0d+00)

c ... in case k runs out of the range
      if(k.le.(-j).or.k.gt.j)then
        cminus=zero
        return
      endif

      dj=dfloat(j)
      dk=dfloat(k)

      cminus=sqrt(dj*(dj+one)-dk*(dk-one))

      return
      end
c-----------------------------------------------------------------------
      double precision function rotcal(jsmal1,ksmal1,Kbig1,jsmal2,
     +                        ksmal2,Kbig2,ip,Ah2o,Bh2o,Ch2o,Jbig)
      implicit double precision(a-h,o-z)

      rotcal=0.d0
      if(jsmal1.ne.jsmal2.or.Kbig1.ne.Kbig2)return

      pre1=dfloat((1+kdel(Kbig2,0)*kdel(ksmal2,0))*
     +            (1+kdel(Kbig1,0)*kdel(ksmal1,0)))
      pre1=1.d0/sqrt(pre1)

      fac1=(Ah2o+Ch2o)*0.5d0*dfloat(jsmal1*(jsmal1+1))+
     +     (Bh2o-0.5d0*(Ah2o+Ch2o))*dfloat(ksmal1*ksmal1)
      ifac1=kdel(ksmal2,ksmal1)+((-1)**(Jbig+ksmal1+ip))*kdel(Kbig1,0)*
     +      kdel(ksmal2,-ksmal1)
      fac1=fac1*dfloat(ifac1)

      ifac2=kdel(ksmal2,ksmal1+2)+((-1)**(Jbig+ksmal1+ip))*kdel(Kbig1,0)
     +      *kdel(-ksmal2,ksmal1+2)
      fac2=cplus(jsmal1,ksmal1)*cplus(jsmal1,ksmal1+1)*dfloat(ifac2)

      ifac3=kdel(ksmal2,ksmal1-2)+((-1)**(Jbig+ksmal1+ip))*kdel(Kbig1,0)
     +      *kdel(-ksmal2,ksmal1-2)
      fac3=cminus(jsmal1,ksmal1)*cminus(jsmal1,ksmal1-1)*dfloat(ifac3)

      fac2=(fac2+fac3)*0.25d0*(Ah2o-Ch2o)

      rotcal=(fac1+fac2)*pre1

      return
      end
c-----------------------------------------------------------------
      subroutine calfac(fact,maxfac)
      implicit double precision(a-h,o-z)
      dimension fact(0:maxfac)

      fact(0)=1.0d+00

      do i=1,maxfac
        fact(i)=fact(i-1)*dfloat(i)
      enddo
      return
      end
c-----------------------------------------------------------------
      SUBROUTINE GAULEGF(X1,X2,X,W,N,ng)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 X1,X2,X(Ng),W(Ng)
      PARAMETER (EPS=3.D-14)
      parameter(Pi=3.14159265358979323846d+00)
      M=(N+1)/2
      XM=0.5D0*(X2+X1)
      XL=0.5D0*(X2-X1)
      DO 12 I=1,M
        Z=COS(pi*(I-.25D0)/(N+.5D0))
1       CONTINUE
          P1=1.D0
          P2=0.D0
          DO 11 J=1,N
            P3=P2
            P2=P1
            P1=((2.D0*J-1.D0)*Z*P2-(J-1.D0)*P3)/J
11        CONTINUE
          PP=N*(Z*P1-P2)/(Z*Z-1.D0)
          Z1=Z
          Z=Z1-P1/PP
        IF(ABS(Z-Z1).GT.EPS)GO TO 1
        X(I)=XM-XL*Z
        X(N+1-I)=XM+XL*Z
        W(I)=2.D0*XL/((1.D0-Z*Z)*PP*PP)
        W(N+1-I)=W(I)
12    CONTINUE
      RETURN
      END
c-----------------------------------------------------------------
      SUBROUTINE GAUCHEB(X1,X2,X,W,N,ng)
      implicit double precision(a-h,o-z)
      dimension x(ng)
      parameter(Pi=3.14159265358979323846d+00)

      w=pi/dfloat(n)
      xl=(x2-x1)/2.d0
      do i=1,n
        phi=pi*dfloat(2*i-1)/dfloat(2*n)
        x(i)=x1+xl*(cos(phi)+1.d0)
      enddo

      return
      end

c-----------------------------------------------------------------
      subroutine transm(numbas,nKbig,icode,Kbgind,
     +                  wgtgl,nlgrid,glgrid,ncgrid,gcgrid,Tsin,Tcos,
     +                  maxfac,fact)
      implicit double precision(a-h,o-z)
c ... prepare cosine and sine transformation matrix for the potential
c ... matrix calculations

      dimension icode(numbas),Kbgind(nKbig+1),wgtgl(nlgrid),
     +          glgrid(nlgrid),gcgrid(ncgrid),
     +          Tsin(nlgrid*ncgrid,numbas),Tcos(nlgrid*ncgrid,numbas),
     +          fact(0:maxfac)

      write(6,*)'in transm'
      do i=1,numbas
c ... extract j,k,K
        jsmall=mod(icode(i)/1000000,100)
        ksmall=mod(icode(i)/10000,100)
        Kbig=mod(icode(i)/100000000,100)
        ib=mod(icode(i)/10,1000)
        ksign=mod(icode(i),10)
        ksmall=ksmall*((-1)**ksign)
c ... loop over chebychev and legendre grids
        do ic=1,ncgrid
          chi=acos(gcgrid(ic))
          do il=1,nlgrid
c ... the composite index for the grid points
            lc=(ic-1)*nlgrid+il
            theta=acos(glgrid(il))
            dsmall=wigd(jsmall,Kbig,ksmall,theta,maxfac,fact)
            coselm=cos(ksmall*chi)*dsmall
            sinelm=sin(ksmall*chi)*dsmall
            prefac=wgtgl(il)*dfloat(2*jsmall+1)/
     +             dfloat(ncgrid*(1+kdel(Kbig,0)*kdel(ksmall,0)))
            prefac=sqrt(prefac)
            Tcos(lc,ib)=prefac*coselm
            Tsin(lc,ib)=prefac*sinelm
c           write(6,*)ic,il,lc,ib,Tcos(lc,ib),Tsin(lc,ib)
          enddo
        enddo
      enddo

      return
      end
c------------------------------------------------------------------
      double precision function wigd(j,m,k,theta,maxfac,fact)
      implicit double precision(a-h,o-z)
      dimension fact(0:maxfac)
c ... this function calculates the wigner d-matrix element.
c ... It takes Eq. 3.57 of Zare, 1988.

      pre1=fact(j+k)*fact(j-k)*fact(j+m)*fact(j-m)
      pre1=sqrt(pre1)

c ... judge the upper bound of nu, the summing index
      nulow=max(0,k-m)
      nuup=min(j+k,j-m)
      thehlf=0.5d+00*theta
      wigd=0.0d+00

c ... summation over nu
      do nu=nulow,nuup
        denorm=fact(j-m-nu)*fact(j+k-nu)*fact(nu+m-k)*fact(nu)
     +          *(-1)**(nu)
        pre2=1.0d+00/denorm
        cosfac=cos(thehlf)
        sinfac=-sin(thehlf)
        cosfac=cosfac**(2*j+k-m-2*nu)
        sinfac=sinfac**(m-k+2*nu)
        wigd=wigd+pre2*cosfac*sinfac
      enddo

      wigd=wigd*pre1

      return
      end
c------------------------------------------------------------------
      subroutine potread(nchi,nthe,nrad,vpes,gchi,gthe,grad,potfil,
     +           lenfil)
      implicit double precision(a-h,o-z)

      dimension vpes(-4:nchi+5,-4:nthe+5,nrad),
     +          gchi(-4:nchi+5),gthe(-4:nthe+5),
     +          grad(nrad)

      character potfil*30

      write(6,'(a)')potfil(1:lenfil)

c     open(9,file='SO2H2_hind.pot',status='old')
      open(9,file=potfil(1:lenfil),status='old')

c ... initialize the pes minimum
      pesmin=100.d0
      do ichi=1,nchi
        do ithe=1,nthe
          do irad=1,nrad
            read(9,*)gchi(ichi),gthe(ithe),grad(irad),
     +               vpes(ichi,ithe,irad)
c           write(6,*)gchi(ichi),gthe(ithe),grad(irad),
c    +                vpes(ichi,ithe,irad)
            if(vpes(ichi,ithe,irad).lt.pesmin) then
              pesmin=vpes(ichi,ithe,irad)
              chimin=gchi(ichi)
              themin=gthe(ithe)
              radmin=grad(irad)
            endif
          enddo
        enddo
      enddo

      write(6,9001)pesmin,radmin,themin,chimin

      close(9,status='keep')

c ... extend to outer range of theta for spline
      do ithe=-4,0
        jthe=2+iabs(ithe)
        do ichi=1,nchi
          do irad=1,nrad
            vpes(ichi,ithe,irad)=vpes(ichi,jthe,irad)
            gthe(ithe)=-gthe(jthe)
          enddo
        enddo
      enddo

      do ithe=nthe+1,nthe+5
        jthe=2*nthe-ithe
        do ichi=1,nchi
          do irad=1,nrad
            vpes(ichi,ithe,irad)=vpes(ichi,jthe,irad)
            gthe(ithe)=180+(ithe-nthe)*10
          enddo
        enddo
      enddo

c ... extend to outer range of chi for spline
      do ichi=-4,0
        jchi=2+iabs(ichi)
        do ithe=-4,nthe+5
          do irad=1,nrad
            vpes(ichi,ithe,irad)=vpes(jchi,ithe,irad)
            gchi(ichi)=-gchi(jchi)
          enddo
        enddo
      enddo

      do ichi=nchi+1,nchi+5
        jchi=2*nchi-ichi
        do ithe=-4,nthe+5
          do irad=1,nrad
            vpes(ichi,ithe,irad)=vpes(jchi,ithe,irad)
            gchi(ichi)=90+(ichi-nchi)*10
          enddo
        enddo
      enddo

 9001 format(/'PES MINIMUM',f15.7,1x,'R=',f7.3,1x,'THETA=',
     +       f10.3,1x,'CHI=',f10.3)
      return
      end
C##############################################################################
      subroutine prev(rinp,vmat,nlgrid,ncgrid,glgrid,
     +                gcgrid,nchi,nthe,nrad,vpes,gchi,gthe,grad,ir,
     +                nrpont)
      implicit double precision(a-h,o-z)

      dimension vmat(nrpont*nlgrid*ncgrid),glgrid(nlgrid),
     +          gcgrid(ncgrid),vpes(-4:nchi+5,-4:nthe+5,nrad),
     +          gchi(-4:nchi+5),gthe(-4:nthe+5),grad(nrad)
      parameter(Pi=3.14159265358979323846d+00)

c ... calculate the initial position for the vector vmat to store matrix element
      numlc=ncgrid*nlgrid
      lc0=(ir-1)*numlc

      do ichi=1,ncgrid
          rchi=acos(gcgrid(ichi))
          dchi=180.d0*rchi/Pi
c ... map dchi into the range of 0:90 degrees
        if(dchi.gt.90.d0.and.dchi.le.180.d0) then
          dchi=180.d0-dchi
        elseif(dchi.gt.180.d0.and.dchi.le.270.d0) then
          dchi=dchi-180.d0
        elseif(dchi.gt.270.d0) then
          dchi=360.d0-dchi
        endif

        do ithe=1,nlgrid
          rthe=acos(glgrid(ithe))
          dthe=180.d0*rthe/Pi
c         write(6,*)ichi,rchi,dchi,ithe,rthe,dthe
          lc=(ichi-1)*nlgrid+ithe
          call h2oh2pes(rinp,dthe,dchi,vpot,nchi,nthe,nrad,vpes,gchi,
     +                  gthe,grad)
          vmat(lc0+lc)=vpot
          if (vmat(lc0+lc) .ge. 1000.d0) vmat(lc0+lc)=1000.d0
c         write(6,*)lc,dthe,dchi,vpot
        enddo
      enddo

      return
      end
c---------------------------------------------------------------------
      subroutine h2oh2pes(radinp,theinp,chiinp,vpot,nchi,nthe,nrad,vpes,
     +                    gchi,gthe,grad)
      implicit double precision(a-h,o-z)
      parameter(rmin=3.0d0,rmax=26.d0)
      dimension vpes(-4:nchi+5,-4:nthe+5,nrad),
     +          gchi(-4:nchi+5),gthe(-4:nthe+5),
     +          grad(nrad)
      parameter (m=100,n=100,l=100,ir2=2)
      dimension xt(m)
      dimension dty(l),ddty(l),s1(1),ds1(1),dds1(1),h1(l)
      dimension dny(n),ddny(n),s2(1),ds2(1),dds2(1),h2(n)
      dimension dhy(m),ddhy(m),s3(1),ds3(1),dds3(1),h3(m)
      dimension y(m),ss(m),sss(m),y2(l)

      if(radinp.gt.rmax)radinp=rmax
      if(radinp.lt.rmin)radinp=rmin

      do irad=1,nrad
        do ithe=-4,nthe+5
          nch=0
          do ichi=-4,nchi+5
            nch=nch+1
            xt(nch)=gchi(ichi)
            y(nch)=vpes(ichi,ithe,irad)
          enddo
          call spline(xt,y,nch,dy1,dyn,y2)
          call splint(xt,y,y2,nch,chiinp,y3)
          ss(ithe)=y3
        enddo
        nth=0
        do ithe=-4,nthe+5
          nth=nth+1
          xt(nth)=gthe(ithe)
          y(nth)=ss(ithe)
        enddo
        call spline(xt,y,nth,dy1,dyn,y2)
        call splint(xt,y,y2,nth,theinp,y3)
        sss(irad)=y3
      enddo
c ... linear interpolation for large r
c     if(radinp.ge.15.d0) then
c       do irad=1,nrad
c         if(grad(irad).ge.radinp) goto 300
c       enddo
c 300   slop=(sss(irad)-sss(irad-1))/(grad(irad)-grad(irad-1))
c       vpot=(radinp-grad(irad-1))*slop
c       return
c     endif
      do irad=1,nrad
        xt(irad)=grad(irad)
        y(irad)=sss(irad)
      enddo
      call spline(xt,y,nrad,dy1,dyn,y2)
      call splint(xt,y,y2,nrad,radinp,y3)
      vpot=y3

      return
      end
C##################################################################
C# SPLINE ROUTINES
C#            Numerical recipes in fortran
C#            Cambrige University Press
C#            York, 2nd edition, 1992.
C##################################################################
      SUBROUTINE splint(xa,ya,y2a,n,x,y)
      implicit double precision  (a-h,o-z)
      DIMENSION xa(n),y2a(n),ya(n)
      klo=1
      khi=n
 1    if (khi-klo.gt.1) then
         k=(khi+klo)/2
         if(xa(k).gt.x)then
            khi=k
          else
            klo=k
          endif
         goto 1
         endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.0d0) write(6,*) 'bad xa input in splint'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))
     1                                                   *(h**2)/6.0d0
      return
      END
C##############################################################################
      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      implicit double precision  (a-h,o-z)
      DIMENSION x(n),y(n),y2(n)
      PARAMETER (NMAX=100)
      DIMENSION u(NMAX)
      if (yp1.gt..99d30) then
         y2(1)=0.0d0
         u(1)=0.0d0
       else
         y2(1)=-0.5d0
         u(1)=(3.0d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
       endif
      do i=2,n-1
         sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
         p=sig*y2(i-1)+2.0d0
         y2(i)=(sig-1.0d0)/p
         u(i)=(6.0d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/
     1                    (x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
         enddo
      if (ypn.gt..99d30) then
          qn=0.0d0
          un=0.0d0
        else
          qn=0.5d0
          un=(3.0d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
        endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.0d0)
      do k=n-1,1,-1
          y2(k)=y2(k)*y2(k+1)+u(k)
          enddo
      return
      END
c------------------------------------------------------------------
      subroutine prtril(d,n)
      implicit double precision(a-h,o-z)

      dimension d(*)

      max=5
      mm1=max-1
      do 120 i0=1,n,max
        il=min(n,i0+mm1)
        write(6,9008)
        write(6,9028)(i,i=i0,il)
        write(6,9008)
        il=-1
        do 100 i=i0,n
          il=il+1
          j0=i0+(i*i-i)/2
          jl=j0+min(il,mm1)
          write(6,9048)i,(d(j),j=j0,jl)
  100   continue
  120 continue

 9008 format(1x)
 9028 FORMAT(15X,10(4X,I4,3X))
 9048 FORMAT(I5,10X,10F11.6)

      return
      end
c--------------------------------------------------------------------------
      subroutine Hv(vec,uec,vmat,ncgrid,nlgrid,rgrid,
     +              nrpont,numbas,Jbig,ip,Kbgind,icode,Tcos,
     +              Tsin,nKbig,coriol,ijcori,rotor,ijrot,
     +              ncorio,nrotor,rkin)
      implicit double precision(a-h,o-z)

      dimension vec(numbas*nrpont),uec(numbas*nrpont),
     +          vmat(nrpont*nlgrid*ncgrid),rgrid(nrpont),
     +          Kbgind(nKbig+1),icode(numbas),
     +          Tsin(nlgrid*ncgrid,numbas),Tcos(nlgrid*ncgrid,numbas),
     +          coriol(ncorio),rotor(nrotor),ijcori(numbas*4,2),
     +          ijrot(numbas*4,2),rkin(nrpont*(nrpont+1)/2)
c ... two work spaces
      dimension work1(ncgrid*nlgrid),work2(ncgrid*nlgrid)
      parameter(half=0.5d0,zero=0.0d0,one=1.d0)

c     do i=1,ncorio
c       write(6,*)coriol(i),(ijcori(i,j),j=1,2)
c     enddo
c     do i=1,nrotor
c       write(6,*)rotor(i),(ijrot(i,j),j=1,2)
c     enddo

c ... add V*v to u

      numlc=ncgrid*nlgrid
      do 100 ir=1,nrpont
c ... initial position for the angular basis
      ib0=(ir-1)*numbas
c ... initial position for the potential matrix element
      lc0=(ir-1)*numlc
c ... for K=0, we need an if statement for using cos or sin transformation
      iKbig=1
      if(mod(Jbig+ip,2).eq.0) then
c       write(6,*)'cos will be used for K=0'
        do icod=Kbgind(iKbig),Kbgind(iKbig+1)-1
          ib=mod(icode(icod)/10,1000)
          do lc=1,numlc
            work1(lc)=zero
            do jcod=Kbgind(iKbig),Kbgind(iKbig+1)-1
              jb=mod(icode(jcod)/10,1000)
              work1(lc)=work1(lc)+Tcos(lc,jb)*vec(jb+ib0)
            enddo
          enddo
          do lc=1,numlc
            uec(ib+ib0)=uec(ib+ib0)+Tcos(lc,ib)*vmat(lc+lc0)*work1(lc)
c           uec(ib+ib0)=uec(ib+ib0)+Tcos(lc,ib)*work1(lc)
          enddo
c         write(6,*)icod,ib,ib+ib0,uec(ib+ib0)
        enddo
      else
c       write(6,*)'sin will be used for K=0'
        do icod=Kbgind(iKbig),Kbgind(iKbig+1)-1
          ib=mod(icode(icod)/10,1000)
          do lc=1,numlc
            work1(lc)=zero
            do jcod=Kbgind(iKbig),Kbgind(iKbig+1)-1
              jb=mod(icode(jcod)/10,1000)
              work1(lc)=work1(lc)+Tsin(lc,jb)*vec(jb+ib0)
            enddo
          enddo
          do lc=1,numlc
            uec(ib+ib0)=uec(ib+ib0)+Tsin(lc,ib)*vmat(lc+lc0)*work1(lc)
c           uec(ib+ib0)=uec(ib+ib0)+Tsin(lc,ib)*work1(lc)
          enddo
c         write(6,*)icod,ib,ib+ib0,uec(ib+ib0)
        enddo
      endif

c     write(6,*)
      do iKbig=2,nKbig
        do icod=Kbgind(iKbig),Kbgind(iKbig+1)-1
          ib=mod(icode(icod)/10,1000)
          do lc=1,numlc
            work1(lc)=zero
            work2(lc)=zero
            do jcod=Kbgind(iKbig),Kbgind(iKbig+1)-1
              jb=mod(icode(jcod)/10,1000)
              work1(lc)=work1(lc)+Tcos(lc,jb)*vec(jb+ib0)
              work2(lc)=work2(lc)+Tsin(lc,jb)*vec(jb+ib0)
            enddo
          enddo
          do lc=1,numlc
c ... the half can be done outside of the loop, but because it require to run
c ... icode, this is not efficient.
            uec(ib+ib0)=uec(ib+ib0)+half*vmat(lc0+lc)*
     +             (Tcos(lc,ib)*work1(lc)+Tsin(lc,ib)*work2(lc))
c           uec(ib+ib0)=uec(ib+ib0)+half*
c    +             (Tcos(lc,ib)*work1(lc)+Tsin(lc,ib)*work2(lc))
          enddo
c         write(6,*)icod,ib,ib+ib0,uec(ib+ib0)
        enddo
c       write(6,*)
      enddo

c ... the following is for the Coriolis matrix multiplication
c ... get the radius minus square, r^{-2}
      rmsqr=one/(rgrid(ir)*rgrid(ir))
      do imatrx=1,ncorio
        elemnt=coriol(imatrx)*rmsqr
        iib=ijcori(imatrx,2)
        jjb=ijcori(imatrx,1)
        ib=iib+ib0
        jb=jjb+ib0
        uec(jb)=uec(jb)+elemnt*vec(ib)
      enddo

c ... the following is for the H2O rotor matrix multiplication
      do imatrx=1,nrotor
        elemnt=rotor(imatrx)
        iib=ijrot(imatrx,2)
        jjb=ijrot(imatrx,1)
        ib=iib+ib0
        jb=jjb+ib0
        uec(jb)=uec(jb)+elemnt*vec(ib)
      enddo
  100 continue

c ... the following loop is for the radial kinetic energy matrix element
      nelmnt=0
      do ir=1,nrpont
        ib0=(ir-1)*numbas
        do jr=1,ir-1
          jb0=(jr-1)*numbas
          nelmnt=nelmnt+1
          elemnt=rkin(nelmnt)
c ... loop over all the angular basis since this interaction is diagonal
c ... for the angular basis
          do iib=1,numbas
            ib=iib+ib0
            jb=iib+jb0
            uec(ib)=uec(ib)+elemnt*vec(jb)
c ...   Hermiticity is used in the next line
            uec(jb)=uec(jb)+elemnt*vec(ib)
          enddo
        enddo
c ...   for the diagonal multiplication
        nelmnt=nelmnt+1
        elemnt=rkin(nelmnt)
        do iib=1,numbas
          ib=iib+ib0
          uec(ib)=uec(ib)+elemnt*vec(ib)
        enddo
      enddo

      return
      end
c-----------------------------------------------------------------
      subroutine fftgrd(X,N,ng)
      implicit double precision(a-h,o-z)
      dimension x(ng)
      parameter(Pi=3.141592653589793d+00)

      do ibeta=1,n
        chi=Pi*(dfloat(ibeta-1)+0.0d0)/n
        x(ibeta)=cos(chi)
      enddo


      return
      end
