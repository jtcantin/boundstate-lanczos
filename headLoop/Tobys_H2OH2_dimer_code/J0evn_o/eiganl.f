      subroutine eiganl(eigvec,ntotbs,neigen,nrpont,numbas,icode,Kbgind,
     +                  Jbig,ib000,ib100,ib101,ib110,ib111,jb111,ip,
     +                  eigval,niter,rgrid,h2_2mu,potfil,Tsin,Tcos,
     +                  ncgrid,nlgrid,wgtgl,gcgrid,glgrid,maxfac,fact,
     +                  jmax,ik,vmat,iplt3d)
      implicit double precision(a-h,o-z)

      dimension eigvec(ntotbs*neigen),icode(numbas),Kbgind(Jbig+2),
     +          eigval(niter),rgrid(nrpont),r2k(-3:2),
     +          Tsin(nlgrid*ncgrid,numbas),Tcos(nlgrid*ncgrid,numbas),
     +          wgtgl(nlgrid),glgrid(nlgrid),gcgrid(ncgrid),
     +          fact(0:maxfac),vmat(nrpont*nlgrid*ncgrid)

c ... rhogrd stores the density at the DVR and Gaussian quadrature points, and
c ... the converted x,y,z coordinates
      dimension rhogrd(nrpont*nlgrid*ncgrid)
      character potfil*30,cubfil*30,cstate*1
c ... define the small and large end of a box
      parameter(xmin=-15.d0,xmax=15.d0,ymin=-15.d0,ymax=15.d0,
     +          zmin=-15.d0,zmax=15.d0)
c ... define the density of plotting grids and the Gaussian width
      parameter(numx=101,numy=101,numz=101,gwidth=1./(2.*2.*1.),
     +          eps=1.d-5,pi=3.14159265358979323846d+00)
c ... conversion factor
      parameter(bo2ang=0.5291772,H2mass=2.015650642)
c ... rhoplt stores the density calculated by Gaussian distributions from
c ... the DVR and quadrature grid density
      dimension rhoplt(numx*numy*numz)
c ... rho2d stores the density at the angular quadrature grid point after
c ... radial normalization
      dimension rho2d(ncgrid*nlgrid),rhochi(2*ncgrid,2)
c ... boltzmann factor in the unit of cm-1/Kelvin
      parameter(boltz=0.6950356d0)
c ... logical flag for constructing and plotting 3d density
      logical plot3d

      write(6,9001)

      plot3d=.true.
      if(iplt3d.eq.0)plot3d=.false.
      if(plot3d) then
        write(6,'(''3D density is to be constructed and plotted'')')
      endif

c     write(6,*)'in eiganl',ntotbs,neigen,nrpont,numbas,eigvec(1),
c    +          eigvec(ntotbs*neigen),icode(1),icode(numbas),Jbig,
c    +          Kbgind(1),Kbgind(Jbig+2)
c     write(6,*)ib000,ib100,ib101,ib110,ib111,jb111,ip
c     write(6,*)(eigval(i),i=1,neigen)
c     write(6,*)(rgrid(i),i=1,nrpont)
c     write(6,*)Tsin(1,1),Tcos(1,1),Tsin(nlgrid*ncgrid,numbas),
c    +          Tcos(nlgrid*ncgrid,numbas)
c     write(6,*)wgtgl(nlgrid),glgrid(nlgrid),gcgrid(ncgrid)
c     write(6,*)maxfac
c     write(6,*)gwidth

c ... get the x y z step size
      xstep=(xmax-xmin)/dfloat(numx-1)
      ystep=(ymax-ymin)/dfloat(numy-1)
      zstep=(zmax-zmin)/dfloat(numz-1)
      dx=(xmax-xmin)/dfloat(numx)
      dy=(xmax-xmin)/dfloat(numy)
      dz=(xmax-xmin)/dfloat(numz)

c ... get the useful info from potfil
      do ich=3,30
        if(potfil(ich-2:ich).eq.'pot')lencha=ich-3
      enddo
      write(cubfil(1:lencha),'(a)')potfil(1:lencha)

      do 100 istate=1,neigen

c ... nullify the expectation values of r^2
      do ipow=-3,2
        r2k(ipow)=0.d0
      enddo
      write(6,9002)istate,eigval(istate),eigval(istate)/boltz
      ib00=(istate-1)*ntotbs
c ... get the K contributions
      do Kbig=0,Jbig
        iKbig=Kbig+1
        istart=Kbgind(iKbig)
        iend=Kbgind(iKbig+1)-1
c       write(6,*)istart,iend
        compnt=0.d0
        do icod=istart,iend
          iib=mod(icode(icod)/10,1000)
          do ir=1,nrpont
            ib0=(ir-1)*numbas
            ib=iib+ib0+ib00
            coef=eigvec(ib)
            compnt=compnt+coef*coef
          enddo
        enddo
        write(6,9003)Kbig,compnt
      enddo

c ... get the contributions from jKaKc states
c ... |000>
      if(ib000.eq.0.or.(mod(Jbig+ip,2).eq.1)) then
        cmp000=0.d0
      else
        cmp000=0.d0
        do ir=1,nrpont
          ib0=(ir-1)*numbas
          ib=ib0+ib00+ib000
          coef=eigvec(ib)
          cmp000=cmp000+coef*coef
        enddo
      endif
      write(6,9004)'000',cmp000

c ... |111>
      if(ib100.eq.0.or.(mod(Jbig+ip,2).eq.1)) then
        compo1=0.d0
      else
        compo1=0.d0
        do ir=1,nrpont
          ib0=(ir-1)*numbas
          ib=ib0+ib00+ib100
          coef=eigvec(ib)
          compo1=compo1+coef*coef
        enddo
      endif

      if(ib101.eq.0) then
        compo2=0.d0
      else
        compo2=0.d0
        do ir=1,nrpont
          ib0=(ir-1)*numbas
          ib=ib0+ib00+ib101
          coef=eigvec(ib)
          compo2=compo2+coef*coef
        enddo
      endif
      cmp111=compo1+compo2
      write(6,9004)'111',cmp111

c ...|101> and |110>
      if(ib110.eq.0) then
        compo1=0.d0
      else
        compo1=0.d0
        do ir=1,nrpont
          ib0=(ir-1)*numbas
          ib=ib0+ib00+ib110
          coef=eigvec(ib)
          compo1=compo1+coef*coef
        enddo
      endif

c ... use the fact that if there's basis 111, there must also be basis 1-11
      if(ib111.eq.0) then
        compo2=0.d0
        compo3=0.d0
      else
        compo2=0.d0
        compo3=0.d0
        do ir=1,nrpont
          ib0=(ir-1)*numbas
          ib=ib0+ib00+ib111
          jb=ib0+ib00+jb111
          coef1=eigvec(ib)
          coef2=eigvec(jb)
          compo2=compo2+(coef1-coef2)*(coef1-coef2)
          compo3=compo3+(coef1+coef2)*(coef1+coef2)
        enddo
      endif
      compo2=0.5d0*compo2
      compo3=0.5d0*compo3

      cmp101=dfloat((1+(-1)**(Jbig+ip))/2)*compo1+compo2
      cmp110=dfloat((1-(-1)**(Jbig+ip))/2)*compo1+compo3

      write(6,9004)'101',cmp101
      write(6,9004)'110',cmp110

c ... calculate the expectation values of R^n
      do ir=1,nrpont
        ib0=(ir-1)*numbas
        compo=0.d0
        do iib=1,numbas
          ib=ib00+ib0+iib
          coef=eigvec(ib)
          compo=compo+coef*coef
        enddo
        rval=rgrid(ir)
        do ipow=-3,2
          r2k(ipow)=r2k(ipow)+compo*(rval**dfloat(ipow))
        enddo
      enddo
      write(6,9005)(ipow,r2k(ipow),ipow=-3,2)
      rotcon=r2k(-2)*h2_2mu
      write(6,9006)rotcon

c ... not plotting 3d density
c     goto 500
c ... calculate and write the density file
c ... clean rhoplt
      do igrid=1,numx*numy*numz
        rhoplt(igrid)=0.d0
      enddo
c ... clean rho2d
      do lc=1,ncgrid*nlgrid
        rho2d(lc)=0.d0
      enddo
c ... clean rhogrd
      do lrc=1,nrpont*nlgrid*ncgrid
        rhogrd(lrc)=0.d0
      enddo

      numlc=nlgrid*ncgrid
      do ir=1,nrpont
        ib0=(ir-1)*numbas+ib00
        do il=1,nlgrid
          do ic=1,ncgrid
            call calrho(eigvec,ntotbs,neigen,nrpont,numbas,icode,Kbgind,
     +                  Jbig,ip,Tsin,Tcos,ncgrid,nlgrid,wgtgl,gcgrid,
     +                  glgrid,ib0,il,ic,rho)
c           if(rho.lt.0.d0) then
c             write(6,*)ir,il,ic,rho
c             stop 'density cannot be negative'
c           endif
            lcr=(ir-1)*numlc+(ic-1)*nlgrid+il
            rhogrd(lcr)=rho
c           theta=acos(glgrid(il))
c           chi=acos(gcgrid(ic))
c           call rhodir(icode,numbas,ib0,theta,chi,rho,maxfac,fact,
c    +                  Jbig,Kbgind,eigvec,neigen,ntotbs,ip)
c ...       check the coordinates of the DVR/quadrature grids
c           if(istate.eq.1.and.ir.eq.1) then
c ...       only punch the density at grids
c           theta=acos(glgrid(il))*180.d0/Pi
c           chi=acos(gcgrid(ic))*180.0/Pi
c           rpoint=rgrid(ir)
c           write(6,*)rpoint,theta,chi,rho
            if(rho.gt.eps) then
c ...         accumulate 2d density
              lc=(ic-1)*nlgrid+il
              rho2d(lc)=rho2d(lc)+rho
c ...         not plotting 3d rho
              if(.not.plot3d)goto 400
c             write(6,*)rho
              theta=acos(glgrid(il))
              chi=acos(gcgrid(ic))
              rpoint=rgrid(ir)
              x=rpoint*sin(theta)*cos(chi)
              y=rpoint*sin(theta)*sin(chi)
              z=rpoint*cos(theta)
c             igrid=0
              do ix=1,numx
                xpt=xmin+(ix-1)*xstep
                iz0=(ix-1)*numy*numz
                do iy=1,numy
                  ypt=ymin+(iy-1)*ystep
                  iz00=(iy-1)*numz
                  do iz=1,numz
                    zpt=zmin+(iz-1)*zstep
                    r1sq=(xpt-x)*(xpt-x)+(ypt-y)*(ypt-y)+(zpt-z)*(zpt-z)
                    r2sq=(xpt-x)*(xpt-x)+(ypt+y)*(ypt+y)+(zpt-z)*(zpt-z)
c                   write(6,*)x,y,z,xpt,ypt,zpt,
c    +                        'r1sq=',r1sq,'r2sq=',r2sq
                    igrid=iz+iz0+iz00
                    rhoplt(igrid)=rhoplt(igrid)+
     +                         rho*(exp(-gwidth*r1sq)+exp(-gwidth*r2sq))
c                   write(6,*)'rhoplt=',rhoplt(igrid)
c ...               check
c                   if(.not.(rhoplt(igrid).ge.0.0))
c    +                   write(6,*)ix,iy,iz,igrid,
c    +                   rhoplt(igrid)
                  enddo
                enddo
              enddo
  400         continue
            endif    
c             write(6,*)ir,il,ic,rpoint,theta,chi,x,y,z
c           endif
          enddo
        enddo
      enddo

      if(.not.plot3d)goto 600

      icube=10+istate
      write(cstate,'(i1)')istate
      open(icube,file=cubfil(1:lencha)//cstate//'.cube',
     +     status='unknown')

      write(icube,9008)
      write(icube,9009)3,xmin,ymin,zmin
      write(icube,9010)numx,xstep,0.d0,0.d0
      write(icube,9010)numy,0.d0,ystep,0.d0
      write(icube,9010)numz,0.d0,0.d0,zstep
      write(icube,9011)
c     write(icube,'(''6 6.0 0 0 0'')')

      do ix=1,numx
        iz0=(ix-1)*numy*numz
        do iy=1,numy
          iz00=(iy-1)*numz
          write(icube,9007)(rhoplt(iz+iz0+iz00),iz=1,numz)
c         do iz=1,numz
c           if(.not.(rhoplt(iz+iz0+iz00).ge.0.d0))
c    +          write(6,*)ix,iy,iz,rhoplt(iz+iz0+iz00)
c         enddo
        enddo
      enddo

      close(icube,status='keep')

  600 continue

      i2d=20+istate
      write(cstate,'(i1)')istate
      open(i2d,file=cubfil(1:lencha)//cstate//'.contour',
     +     status='unknown')

      do il=1,nlgrid
        theta=acos(glgrid(il))
        theta=180.d0*theta/Pi
        do ic=1,ncgrid
          chi=acos(gcgrid(ic))
          chi=180.d0*chi/Pi
          lc=(ic-1)*nlgrid+il
          if(ic.eq.1.) then
            write(i2d,9013)theta,0.d0,rho2d(lc)
            write(i2d,9013)theta,chi,rho2d(lc)
          elseif(ic.eq.ncgrid) then
            write(i2d,9013)theta,chi,rho2d(lc)
            write(i2d,9013)theta,180.d0,rho2d(lc)
          else
            write(i2d,9013)theta,chi,rho2d(lc)
          endif
        enddo
        write(i2d,*)
      enddo
      close(i2d,status='keep')

  500 continue
c ... check the normalization of rho2d
      sumt=0.d0
      do il=1,nlgrid
        sumc=0.d0
        do ic=1,ncgrid
          lc=(ic-1)*nlgrid+il
          sumc=sumc+rho2d(lc)*2.0
        enddo
        sumt=sumt+wgtgl(il)*sumc
      enddo
      sumt=sumt*pi/dfloat(ncgrid)
      write(6,*)'integral of rho2d',sumt

c ... check the normalization of rhogrd and calculate average potential.
c ... calculate the classical moment of inertia.
      sumrzz=0.d0
      sumryy=0.d0
      sumrxx=0.d0
      sumrxy=0.d0
      sumrxz=0.d0
      sumryz=0.d0
      sumr=0.d0
      do ir=1,nrpont
        sumt=0.d0
        sumtxx=0.d0
        sumtyy=0.d0
        sumtzz=0.d0
        sumtxy=0.d0
        sumtxz=0.d0
        sumtyz=0.d0
        do il=1,nlgrid
          sumc=0.d0
          sumcxx=0.d0
          sumcyy=0.d0
          sumczz=0.d0
          sumcxy=0.d0
          sumcxz=0.d0
          sumcyz=0.d0
          theta=acos(glgrid(il))
          do ic=1,ncgrid
            chi=acos(gcgrid(ic))
            lrc=(ir-1)*nlgrid*ncgrid+(ic-1)*nlgrid+il
            sumc=sumc+rhogrd(lrc)*2.d0*vmat(lrc)
            sumcxx=sumcxx+rhogrd(lrc)*2.d0*(sin(chi)**2)
            sumcyy=sumcyy+rhogrd(lrc)*2.d0*(cos(chi)**2)
            sumczz=sumczz+rhogrd(lrc)*2.d0
            sumcxy=sumcxy-rhogrd(lrc)*(sin(chi)*cos(chi)+sin(2*Pi-chi)*
     +             cos(2*Pi-chi))
            sumcxz=sumcxz-rhogrd(lrc)*(cos(chi)+cos(2*Pi-chi))
            sumcyz=sumcyz-rhogrd(lrc)*(sin(chi)+sin(2*Pi-chi))
          enddo
          sumt=sumt+sumc*wgtgl(il)
          sumtxx=sumtxx+(sumczz*(cos(theta)**2)+sumcxx*(sin(theta)**2))
     +           *wgtgl(il)
          sumtyy=sumtyy+(sumczz*(cos(theta)**2)+sumcyy*(sin(theta)**2))
     +           *wgtgl(il)
          sumtzz=sumtzz+sumczz*wgtgl(il)*(sin(theta)**2)
          sumtxy=sumtxy+sumcxy*wgtgl(il)*(sin(theta)**2)
          sumtxz=sumtxz+sumcxz*wgtgl(il)*sin(theta)*cos(theta)
          sumtyz=sumtyz+sumcyz*wgtgl(il)*sin(theta)*cos(theta)
        enddo
        sumr=sumr+sumt
        sumrxx=sumrxx+sumtxx*(rgrid(ir)**2)
        sumryy=sumryy+sumtyy*(rgrid(ir)**2)
        sumrzz=sumrzz+sumtzz*(rgrid(ir)**2)
        sumrxy=sumrxy+sumtxy*(rgrid(ir)**2)
        sumrxz=sumrxz+sumtxz*(rgrid(ir)**2)
        sumryz=sumryz+sumtyz*(rgrid(ir)**2)
      enddo
      sumr=sumr*pi/dfloat(ncgrid)
      sumrxx=sumrxx*pi/dfloat(ncgrid)
      sumryy=sumryy*pi/dfloat(ncgrid)
      sumrzz=sumrzz*pi/dfloat(ncgrid)
      sumrxy=sumrxy*pi/dfloat(ncgrid)
      sumrxz=sumrxz*pi/dfloat(ncgrid)
      sumryz=sumryz*pi/dfloat(ncgrid)
c ... conver the unit of sumrzz to amu*Angs**2
      sumrxx=sumrxx*bo2ang*bo2ang*H2mass
      sumryy=sumryy*bo2ang*bo2ang*H2mass
      sumrzz=sumrzz*bo2ang*bo2ang*H2mass
      sumrxy=sumrxy*bo2ang*bo2ang*H2mass
      sumrxz=sumrxz*bo2ang*bo2ang*H2mass
      sumryz=sumryz*bo2ang*bo2ang*H2mass
      write(6,*)'integral of rhogrd for average potential:',sumr,
     +           " cm-1; ",sumr/boltz,"K"
      write(6,*)'Ixx of H2 =',sumrxx,' Amu*Angs**2'
      write(6,*)'Iyy of H2 =',sumryy,' Amu*Angs**2'
      write(6,*)'Izz of H2 =',sumrzz,' Amu*Angs**2'
      write(6,*)'Ixy of H2 =',sumrxy,' Amu*Angs**2'
      write(6,*)'Ixz of H2 =',sumrxz,' Amu*Angs**2'
      write(6,*)'Iyz of H2 =',sumryz,' Amu*Angs**2'

c ... calculate the chi density
      i1dchi=30+istate
      open(i1dchi,file=cubfil(1:lencha)//cstate//'.chirho',
     +     status='unknown')
      sumc=0.d0
      do ic=1,ncgrid
        chi=acos(gcgrid(ic))
        sumt=0.d0
        do il=1,nlgrid
          lc=(ic-1)*nlgrid+il
          sumt=sumt+rho2d(lc)*wgtgl(il)
        enddo
        rhochi(ic,2)=sumt*Pi/180.d0
        rhochi(2*ncgrid-ic+1,2)=sumt*Pi/180.d0
        rhochi(ic,1)=chi*180d0/Pi
        rhochi(2*ncgrid-ic+1,1)=(2*Pi-chi)*180d0/Pi
        sumc=sumc+2.0*sumt
c       write(6,*)'sumc=',sumc
      enddo
      sumc=sumc*Pi/dfloat(ncgrid)
      write(6,*)'rhochi normalized to',sumc
      do ic=1,2*ncgrid
        write(i1dchi,9014)(rhochi(ic,i),i=1,2)
      enddo
      close(i1dchi,status='keep')

c ... calculate the theta density
      d180=180.0
      zero=0.0
      i1dthe=40+istate
      sumt=0.d0
      open(i1dthe,file=cubfil(1:lencha)//cstate//'.therho',
     +     status='unknown')
      write(i1dthe,9014)zero,zero
      do it=nlgrid,1,-1
        theta=acos(glgrid(it))
        sumc=0.d0
        do ic=1,ncgrid
          lc=(ic-1)*nlgrid+it
          sumc=sumc+rho2d(lc)*2.d0
        enddo
        sumc=sumc*Pi/dfloat(ncgrid)
        write(i1dthe,9014)theta*180.d0/Pi,sumc*sin(theta)*Pi/180.d0
c       write(i1dthe,9014)theta*180.d0/Pi,sumc*sin(theta)*Pi/180.d0*
c    +    wgtgl(it)/(Pi/(nlgrid-1))
        sumt=sumt+sumc*wgtgl(it)
      enddo
      write(i1dthe,9014)d180,zero
      close(i1dthe,status='keep')
      write(6,*)'rhothe normalized to ',sumt
      write(6,*)
  100 continue

c ... write eigenvector information for the calculation of transition line
c ... strength
      ieigen=21
      open(ieigen,file='eigvec',status='unknown')
      write(ieigen,'(9I6)')Jbig,jmax,ik,ip,numbas,nrpont,ntotbs,neigen
      do istate=1,neigen
        ib0=(istate-1)*ntotbs
        write(ieigen,9012)(eigvec(ib+ib0),ib=1,ntotbs)
      enddo
      close(ieigen,status='keep')

 9001 format(/'ENTERING EIGEN STATE ANALYSIS')
 9002 format(/'STATE:',I4,1x,'ENERGY:',f9.4,'CM-1',f9.4,'K')
 9003 format('   K=',I4,1x,'CONTRIBUTION:',f9.4)
 9004 format('jKaKc=',a,1x,'CONTRIBUTION:',f9.4)
 9005 format(6('R**',I2,'=',f9.4,2x))
 9006 format('EFFECTIVE ROTATIONAL CONSTANT=',f9.4,'cm-1')
 9007 format(6(1x,F15.8))
 9008 format('SO2H2 B CUBE FILE'/
     +       'OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z')
 9009 format(I1,3(1x,f6.2))
 9010 format(I3,3(1x,f6.2))
 9011 format('16 16.0 0 0 0.68179'/'8 8.0  2.33754 0 -0.68142'/
     +       '8 8.0 -2.33754 0 -0.68142')
 9012 format(1P,5E15.8)
 9013 format(3(1x,f12.6))
 9014 format(2(1x,f12.6))

      return
      end
c----------------------------------------------------------------------------
      subroutine calrho(eigvec,ntotbs,neigen,nrpont,numbas,icode,Kbgind,
     +                  Jbig,ip,Tsin,Tcos,ncgrid,nlgrid,wgtgl,gcgrid,
     +                  glgrid,ib0,il,ic,rho)
      implicit double precision(a-h,o-z)
      dimension eigvec(ntotbs*neigen),icode(numbas),Kbgind(Jbig+2),
     +          Tsin(nlgrid*ncgrid,numbas),Tcos(nlgrid*ncgrid,numbas),
     +          wgtgl(nlgrid),glgrid(nlgrid),gcgrid(ncgrid)

      parameter(Pi=3.14159265358979323846d+00)

      rho=0.d0
      lc=(ic-1)*nlgrid+il
c ... for K=0, we need an if statement for using cos or sin transformation
      iKbig=1
      if(mod(Jbig+ip,2).eq.0) then
c ... cos will be used
        do icod=Kbgind(iKbig),Kbgind(iKbig+1)-1
          ib=mod(icode(icod)/10,1000)
c ...     iib is to extract the coefficient.  The same as jjb below.
          iib=ib+ib0
          coef1=eigvec(iib)
          do jcod=Kbgind(iKbig),Kbgind(iKbig+1)-1
            jb=mod(icode(jcod)/10,1000)
            jjb=jb+ib0
            coef2=eigvec(jjb)
            rho=rho+coef1*coef2*Tcos(lc,ib)*Tcos(lc,jb)
          enddo
        enddo
      else
c ... sin will be used
        do icod=Kbgind(iKbig),Kbgind(iKbig+1)-1
          ib=mod(icode(icod)/10,1000)
c ...     iib is to extract the coefficient.  The same as jjb below.
          iib=ib+ib0
          coef1=eigvec(iib)
          do jcod=Kbgind(iKbig),Kbgind(iKbig+1)-1
            jb=mod(icode(jcod)/10,1000)
            jjb=jb+ib0
            coef2=eigvec(jjb)
            rho=rho+coef1*coef2*Tsin(lc,ib)*Tsin(lc,jb)
          enddo
        enddo
      endif

      nKbig=Jbig+1
      do iKbig=2,nKbig
        do icod=Kbgind(iKbig),Kbgind(iKbig+1)-1
          ib=mod(icode(icod)/10,1000)
          iib=ib+ib0
          coef1=eigvec(iib)
          do jcod=Kbgind(iKbig),Kbgind(iKbig+1)-1
            jb=mod(icode(jcod)/10,1000)
            jjb=jb+ib0
            coef2=eigvec(jjb)
            rho=rho+0.5d0*coef1*coef2*(Tsin(lc,ib)*Tsin(lc,jb)+
     +                                 Tcos(lc,ib)*Tcos(lc,jb))
          enddo
        enddo
      enddo

      rho=rho*dfloat(ncgrid)/(2.d0*Pi*wgtgl(il))

      return
      end
c----------------------------------------------------------------------------
      subroutine rhodir(icode,numbas,ib0,theta,chi,rho,maxfac,fact,
     +                  Jbig,Kbgind,eigvec,neigen,ntotbs,ip)
      implicit double precision(a-h,o-z)

      dimension icode(numbas),
     +          Kbgind(Jbig+2),eigvec(neigen*ntotbs),
     +          fact(0:maxfac)

      rho=0.d0
c     write(6,*)theta,chi
      do Kbig=0,Jbig
c       iKbig=Kbig+1
c       istart=Kbgind(iKbig)
c       iend=Kbgind(iKbig+1)-1
c       do icod=istart,iend
c         iib=mod(icode(icod)/10,1000)
c         ib=iib+ib0
c         coef1=eigvec(ib)
c         jsmal1=mod(icode(icod)/1000000,100)
c         ksmal1=mod(icode(icod)/10000,100)
c         do jcod=istart,iend
c           jjb=mod(icode(jcod)/10,1000)
c           jb=jjb+ib0
c           coef2=eigvec(jb)
c           jsmal2=mod(icode(jcod)/1000000,100)
c           ksmal2=mod(icode(jcod)/10000,100)
c           pre=dfloat((2*jsmal1+1)*(2*jsmal2+1))/
c    +          dfloat((1+kdel(ksmal1,0)*kdel(Kbig,0))*
c    +                 (1+kdel(ksmal2,0)*kdel(Kbig,0)))
c           pre=sqrt(pre)*coef1*coef2/(4.d0*Pi)
c           fac1=cos(ksmal2*chi)*cos(ksmal1*chi)
c         enddo
c       enddo
      enddo

      return
      end
