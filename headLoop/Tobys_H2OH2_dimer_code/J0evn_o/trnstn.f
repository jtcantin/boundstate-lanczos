      program trnstn
      implicit double precision(a-h,o-z)
c ... this program is to calculate the transition line strength for two sets of
c ... states from h2oh2.x calculations

      character file1*30,file2*30,argum*30
      parameter (mxrgrd=100,maxbas=1000,maxeig=10,maxfac=50)
      parameter(eps=1.d-14,eps2=1.d-4)
      dimension eigvc1(mxrgrd*maxbas*maxeig),
     +          eigvc2(mxrgrd*maxbas*maxeig),
c ... trndip stores the dipole matrix elements between the angular basis of the
c ... two states
     +          trndip(maxbas*maxbas),fact(0:maxfac),
     +          ijdip(maxbas*maxbas,2)

c ... prepare factorial
      call calfac(fact,maxfac)

c ... test 3j symbol
c     read(5,*)j1,j2,j3,m1,m2,m3
c     call cgc(j1,m1,j2,m2,j3,m3,w3j,1,fact,maxfac)
c     write(6,*)w3j

      call getarg(1,argum)
      read(argum,*)file1
      lenfl1=lastch(file1,30)
      call getarg(2,argum)
      read(argum,*)file2
      lenfl2=lastch(file2,30)

c     write(6,'(a)')file1(1:lenfl1)
c     write(6,'(a)')file2(1:lenfl2)
      ifil1=7
      ifil2=8

      open(ifil1,file=file1(1:lenfl1),status='old')
      read(ifil1,*)Jbig1,jmax1,ik1,ip1,numbs1,nrpt1,ntotb1,neign1
      write(6,*)Jbig1,jmax1,ik1,ip1,numbs1,nrpt1,ntotb1,neign1
      do istate=1,neign1
        ib0=(istate-1)*ntotb1
        read(ifil1,9001)(eigvc1(ib+ib0),ib=1,ntotb1)
      enddo

      close(ifil1,status='keep')

      open(ifil2,file=file2(1:lenfl2),status='old')
      read(ifil2,*)Jbig2,jmax2,ik2,ip2,numbs2,nrpt2,ntotb2,neign2
      write(6,*)Jbig2,jmax2,ik2,ip2,numbs2,nrpt2,ntotb2,neign2
      do istate=1,neign2
        ib0=(istate-1)*ntotb2
        read(ifil2,9001)(eigvc2(ib+ib0),ib=1,ntotb2)
      enddo

      close(ifil2,status='keep')

c ... safety check
      if(nrpt1.ne.nrpt2)stop'nrpt1 |= nrpt2'

c ... loop over all basis pairs to generate the dipole matrix for the
c ... angular basis

c ... loop over state 1, the bra state
      imat=0
      ib=0
      do jsmal1=0,jmax1
        ksmmx1=kstcal(jsmal1,ik1)
        do ksmal1=-ksmmx1,ksmmx1,2
          Kbgst1=0
          if(mod(Jbig1+ip1,2).eq.1.and.ksmal1.eq.0)Kbgst1=1
          Kbged1=min(jsmal1,Jbig1)
          do Kbig1=Kbgst1,Kbged1
            if(.not.(ksmal1.lt.0.and.Kbig1.eq.0)) then
              ib=ib+1

c ... loop over state 2, the ket state
              jb=0
              do jsmal2=0,jmax2
                ksmmx2=kstcal(jsmal2,ik2)
                do ksmal2=-ksmmx2,ksmmx2,2
                  Kbgst2=0
                  if(mod(Jbig2+ip2,2).eq.1.and.ksmal2.eq.0)Kbgst2=1
                  Kbged2=min(jsmal2,Jbig2)
                  do Kbig2=Kbgst2,Kbged2
                    if(.not.(ksmal2.lt.0.and.Kbig2.eq.0)) then
                      jb=jb+1
c ...                 calculate the matrix element
c                     write(6,*)jsmal1,ksmal1,Kbig1,jsmal2,ksmal2,Kbig2
                      call calmat(jsmal1,ksmal1,Kbig1,Jbig1,ip1,
     +                jsmal2,ksmal2,Kbig2,Jbig2,ip2,elemnt,fact,maxfac)
c                     call calovp(jsmal1,ksmal1,Kbig1,Jbig1,ip1,
c    +                jsmal2,ksmal2,Kbig2,Jbig2,ip2,elemnt)
                      if(abs(elemnt).gt.eps) then
                        write(6,*)jsmal1,ksmal1,Kbig1,Jbig1,jsmal2,
     +                            ksmal2,Kbig2,Jbig2,elemnt
                        imat=imat+1
                        trndip(imat)=elemnt
                        ijdip(imat,1)=ib
                        ijdip(imat,2)=jb
                      endif
                    endif
                  enddo
                enddo
              enddo

            endif
          enddo
        enddo
      enddo

      if(ib.ne.numbs1)stop 'ib != numbs1'
      if(jb.ne.numbs2)stop 'jb != numbs2'

      nmat=imat
      write(6,*)'# of nonzero matrix elements:',nmat

c ... get the Jmin for sum over M
      Jmin=min(Jbig1,Jbig2)
c     write(6,*)Jmin
c ... loop over pairs of states to calculate line strengths
      do ist1=1,neign1
        ib00=(ist1-1)*ntotb1
        do ist2=1,neign2
          jb00=(ist2-1)*ntotb2
          suma=0.d0
          do imat=1,nmat
            dipole=trndip(imat)
            ib=ijdip(imat,1)
            jb=ijdip(imat,2)
            do ir=1,nrpt1
              ib0=(ir-1)*numbs1
              coef1=eigvc1(ib00+ib0+ib)
              jb0=(ir-1)*numbs2
              coef2=eigvc2(jb00+jb0+jb)
              suma=suma+coef1*coef2*dipole
c             if(abs(coef1*coef2*dipole).gt.eps2)write(6,*)
c    +          ist1,ist2,ib,jb,ir1,ir2,coef1,coef2,dipole
            enddo
          enddo
c         write(6,*)ist1,ist2,'suma=',suma
          facM=0.d0
          streng=0.d0
          do Mbig=-Jmin,Jmin
            sgn=dfloat((-1)**Mbig)
            call cgc(Jbig1,-Mbig,1,0,Jbig2,Mbig,w3j,1,fact,maxfac)
c           write(6,*)w3j
            facM=facM+w3j*w3j
          enddo
          streng=facM*suma*suma*3.d0
c         streng=streng*streng*3.d0
          write(6,*)ist1,ist2,streng,suma,facM
c         write(6,*)ist1,ist2,streng
        enddo
      enddo

 9001 format(1P,5E15.8)

      end
c---------------------------------------------------------------------
      integer function lastch(line,len)
      implicit double precision(a-h,o-z)
      character*(*) line

      do i=len,1,-1
         if(line(i:i).ne.' ') then
           lastch=i
           return
         endif
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
c----------------------------------------------------------------------
      subroutine calmat(jsmal1,ksmal1,Kbig1,Jbig1,ip1,
     +           jsmal2,ksmal2,Kbig2,Jbig2,ip2,elemnt,fact,maxfac)
      implicit double precision(a-h,o-z)

      dimension fact(0:maxfac)

      elemnt=0.d0

      pre=2.d0*sqrt(dfloat((1+kdel(ksmal1,0)*kdel(Kbig1,0))*
     +    (1+kdel(ksmal2,0)*kdel(Kbig2,0))))
      pre=1.d0/pre
      pre=pre*dfloat((-1)**ksmal1)*sqrt(dfloat(2*jsmal2+1))*
     +    sqrt(dfloat(2*jsmal1+1))*sqrt(dfloat(2*Jbig2+1))*
     +    sqrt(dfloat(2*Jbig1+1))

      sgn1=dfloat((-1)**(Jbig2+ip2+ksmal2))
      sgn2=dfloat((-1)**(Jbig1+ip1+ksmal1))
      sgn3=sgn1*sgn2

c     write(6,*)pre,sgn1,sgn2,sgn3

      summa=0.d0
      do isigma=-1,1
c ...   calculate the first term in the summation
c       write(6,*)'isigma=',isigma
        call cgc(jsmal1,-ksmal1,1,0,jsmal2,ksmal2,w3j1,1,fact,maxfac)
        call cgc(jsmal1,-Kbig1,1,isigma,jsmal2,Kbig2,w3j2,1,fact,
     +           maxfac)
        call cgc(Jbig1,-Kbig1,1,isigma,Jbig2,Kbig2,w3j3,1,fact,
     +           maxfac)
        term1=w3j1*w3j2*w3j3
c       write(6,*)w3j1,w3j2,w3j3,term1
        call cgc(jsmal1,-ksmal1,1,0,jsmal2,-ksmal2,w3j1,1,fact,maxfac)
        call cgc(jsmal1,-Kbig1,1,isigma,jsmal2,-Kbig2,w3j2,1,fact,
     +           maxfac)
        call cgc(Jbig1,-Kbig1,1,isigma,Jbig2,-Kbig2,w3j3,1,fact,maxfac)
        term2=sgn1*w3j1*w3j2*w3j3
c       write(6,*)w3j1,w3j2,w3j3,term2
        call cgc(jsmal1,ksmal1,1,0,jsmal2,ksmal2,w3j1,1,fact,maxfac)
        call cgc(jsmal1,Kbig1,1,isigma,jsmal2,Kbig2,w3j2,1,fact,maxfac)
        call cgc(Jbig1,Kbig1,1,isigma,Jbig2,Kbig2,w3j3,1,fact,maxfac)
        term3=sgn2*w3j1*w3j2*w3j3
c       write(6,*)w3j1,w3j2,w3j3,term3
        call cgc(jsmal1,ksmal1,1,0,jsmal2,-ksmal2,w3j1,1,fact,maxfac)
        call cgc(jsmal1,Kbig1,1,isigma,jsmal2,-Kbig2,w3j2,1,fact,maxfac)
        call cgc(Jbig1,Kbig1,1,isigma,Jbig2,-Kbig2,w3j3,1,fact,maxfac)
        term4=sgn3*w3j1*w3j2*w3j3
c       write(6,*)w3j1,w3j2,w3j3,term4
        summa=summa+term1+term2+term3+term4
c       write(6,*)term1,terms,term3,term4,summa
      enddo

      elemnt=summa*pre

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

c Calculate the Clebsh-Gordan coefficient (or the 3-j symbol)
c The parameter ind3j.eq.1 indicates that the 3-J symbol is returned
c
        subroutine cgc(j1,m1,j2,m2,j,m3,value,ind3j,f,maxfac)
        implicit real*8 (a-h,o-z)
        dimension f(0:maxfac)
c
c ...   swabbing
        m=m3

        d3jfact = 1.d0
        if(ind3j.eq.1) then
         d3jfact = ((-1.d0)**(j1-j2-m))/dsqrt(dfloat(2*j+1))
         m = -m
        endif

        value=0.d0

c
c Check the triangle conditions
c
c       if(j.gt.(j1+j2)) write(6,*)'triangle violated'
        if(j.gt.(j1+j2)) then
          value=0.d0
          return
        endif
c       if(j.lt.abs(j1-j2)) write(6,*)'triangle violated'
        if(j.lt.abs(j1-j2)) then
          value=0.d0
          return
        endif
        if((m1+m2).ne.m) then
          value = 0.d0
          return
        endif

c
c Calculation proper... the pre-sum factor....
c
        facn = (2*j+1)*f(j1+j2-j)*f(j1-m1)*f(j2-m2)*f(j+m)*f(j-m)
        facd = f(j1+j2+j+1)*f(j+j1-j2)*f(j+j2-j1)*f(j1+m1)*f(j2+m2)
        fac = dsqrt(facn/facd)

c
c determine the limit of k summation...
c
        kmax = min(j2+j-m1,j-m,j1-m1)
        if(kmax.lt.0) kmax = 0
        kmin = max(-j1-m1,-j2+j-m1,0)

c
c perform the summation (at least one cycle must be completed...
c
        sum = 0.d0
        do k = kmin,kmax
         facn = f(j1+m1+k)*f(j2+j-m1-k)
         facd = f(k)*f(j-m-k)*f(j1-m1-k)*f(j2-j+m1+k)
         sum = sum + (facn/facd)*(-1)**k
        end do
        value = d3jfact*fac*sum*(-1)**(j1-m1)
c       write(10,*)j1,m1,j2,m2,j,m,value,ind3j
       return
       end
c------------------------------------------------------------------------
      subroutine calovp(jsmal1,ksmal1,Kbig1,Jbig1,ip1,
     +                jsmal2,ksmal2,Kbig2,Jbig2,ip2,elemnt)
      implicit double precision(a-h,o-z)

      elemnt=0.d0

      pre=2.d0*sqrt(dfloat((1+kdel(ksmal1,0)*kdel(Kbig1,0))*
     +    (1+kdel(ksmal2,0)*kdel(Kbig2,0))))
      pre=1.d0/pre

      pre=pre*dfloat(kdel(jsmal1,jsmal2)*kdel(Kbig1,Kbig2)*
     +               kdel(Jbig1,Jbig2))

      iterm1=1+(-1)**(Jbig1+Jbig2+ip1+ip2+ksmal1+ksmal2)
      iterm1=kdel(ksmal1,ksmal2)*iterm1

      iterm2=(-1)**(Jbig1+ip1+ksmal1)+(-1)**(Jbig2+ip2+ksmal2)
      iterm2=kdel(-ksmal1,ksmal2)*kdel(Kbig2,0)*iterm2
      iterm=iterm1+iterm2

      elemnt=dfloat(iterm)*pre


      return
      end
