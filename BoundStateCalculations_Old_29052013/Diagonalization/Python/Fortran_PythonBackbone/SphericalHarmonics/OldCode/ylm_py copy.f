c-----------------------------------------------------------------------
c To get a nice numbering sequence for (l m), use the mapping 
c        f: (l m) -> n, where n = f(l,m) = l^2 + l + m + 1 
        subroutine ylm(L,M,th1,phi,ylmout,p)
        implicit none

        real*8 fact
        common/factorial/ fact(0:40)

        real*8 pifact
        data pifact/12.566370614359d0/

c Interface variables
        real*8 ylmout
        real*8 M
        real*8 L
        real*8 th1
        real*8 phi

c Internal variables
        real*8 MM
        real*8 p
        dimension p(0:50,0:50)
        real*8 costh1

        costh1=dcos(th1)
        call plmrb(p,costh1,L)

c        Normal coeffcient builed in p(L,M)
        if(M.eq.0)then
        ylmout=p(L,M)
        else
         if(M.gt.0)then
          ylmout=sqrt(2.0)*p(L,M)*dcos(M*phi)
         else
          MM=-M
          ylmout=sqrt(2.0)*(-1)**(MM)*p(L,MM)*dsin(-MM*phi)
         endif
        endif
        return
        end

c-----------------------------------------------------------------------
c
c Compute the set of associated Legendre polynomials P_lm
c for l=0,1,...,lmax, and m=0,1,...,l. First the standard
c polynomials
c
c   P^m_l(x) = (1/2^l l!)(1-x^2)^(m/2) (d^(l+m)/d x^(l+m))(x^2 -1)^l
c
c are computed, and then multiplied by
c
c  (-1)^m sqrt[(2l+1)(l-m)!/2(l+m)!]/sqrt(2Pi)
c
c to get the P_lm polynomials....
c
        subroutine plmrb(p,x,lmax)
        implicit none

        real*8 fact
        common/factorial/fact(0:40)
c inverse of dsqrt(2Pi)
        real*8 twopinv
        data twopinv /0.3989422804014d0/

c Interface variables
        real*8 p
        dimension p(0:50,0:50)
        real*8 x
        real*8 lmax
        
c Internal variables
        real*8 u
        real*8 m
        real*8 l
        real*8 pp
        real*8 mm
        real*8 dnorm
        
c
c starting value
c
        p(0,0) = 1.d0
        u = dsqrt(1-x*x)
c
c compute the diagonal elements
c
        do l=1,lmax
         p(l,l) = (2*l-1)*p(l-1,l-1)*u
        end do
c
c compute P_lm along the columns with fixed m
c
        do m = 0,lmax-1
        do l = m,lmax-1
         if((l-1).lt.m) then
           pp = 0
         else
           pp = p(l-1,m)
         endif
         p(l+1,m) = ((2*l+1)*x*p(l,m)-(l+m)*pp)/(l-m+1)
        end do
        end do
c
c Renormalize values...
c
        do l=0,lmax
        mm = 1
        do m=0,l
         dnorm = fact(l-m)*(2*l+1)/(2*fact(l+m))
         p(l,m) = mm*twopinv*dsqrt(dnorm)*p(l,m)
         mm = -mm
        end do
        end do
c
        return
        end

