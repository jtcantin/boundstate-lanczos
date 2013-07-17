c Module that contains the functions required for calculation of the 
c Real Spherical Harmonics, a.k.a. Tesseral Spherical Harmonics
c-----------------------------------------------------------------------
c Code originally written by Toby Zeng, but modified to have a Python
c interface by Joshua Cantin on 19 Apr 2013
c-----------------------------------------------------------------------
c To get a nice numbering sequence for (l m), use the mapping 
c        f: (l m) -> n, where n = f(l,m) = l^2 + l + m + 1 
c-----------------------------------------------------------------------
c NOTE: this code could be optimized further by having the p array
c as an input that was already created, instead of having Fortran
c recreate it every time; it may save up to half of the computational 
c time (see http://www.ucs.cam.ac.uk/docs/course-notes/unix-courses/
c pythonfortran/files/f2py.pdf, pg. 55 and pg. 48 to 55)

c This module can be recompiled by deleting "ylm_py2.pyf" and then running:
c f2py -h ylm_py2.pyf -m SphFort ylm_py2.f
c to remake the signature file (ylm_py2.pyf)
c
c Then, in ylm_py2.pyf, change:
c "real*8 :: ylmout" to "real*8, intent(out) :: ylmout"
c and
c "real*8 dimension(51,51) :: p" to "real*8 dimension(51,51), intent(inout) :: p"
c
c Then, run:
c f2py -c -m SphFort ylm_py2.pyf ylm_py2.f
c to compile and make the file SphFort.so
c 
c The module can be imported into Python via:
c "import SphFort"
c This is assuming the file SphFort.so is in Python's Path.

c-----------------------------------------------------------------------
c Subroutine that calculates the Tesseral Spherical Harmonics
c given:
c L : total angular momentum
c M : projection onto z-axis
c th1 : polar angle
c phi : azimuthal angle
c ylmout is the variable returned
c-----------------------------------------------------------------------

        subroutine ylm(L,M,th1,phi,ylmout)
        implicit none

        real*8 fact, pifact
c        common/factorial/ fact(0:40)
        data pifact/12.566370614359d0/

c Interface variables
        real*8 ylmout,M,L,th1,phi

c Internal variables
        real*8 MM,p,costh1
        dimension p(0:50,0:50)

        costh1=dcos(th1)
        call plmrb(p,costh1,L)

c Normal coeffcient builed in p(L,M)
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
c Subroutine that computes the set of associated Legendre polynomials P_lm
c for l=0,1,...,lmax, and m=0,1,...,l. First the standard
c polynomials
c
c   P^m_l(x) = (1/2^l l!)(1-x^2)^(m/2) (d^(l+m)/d x^(l+m))(x^2 -1)^l
c
c are computed, and then multiplied by
c
c  (-1)^m sqrt[(2l+1)(l-m)!/2(l+m)!]/sqrt(2Pi)
c
c to get the P_lm polynomials.
c-----------------------------------------------------------------------

        subroutine plmrb(p,x,lmax)
        implicit none

        real*8 fact
c        common/factorial/ fact(0:40)

c Inverse of dsqrt(2Pi)
        real*8 twopinv
        data twopinv /0.3989422804014d0/

c Interface variables
        real*8 p,x,lmax
        dimension p(0:50,0:50)
       
c Internal variables
        real*8 u,m,l,pp,mm,dnorm     

c Starting value
        p(0,0) = 1.d0
        u = dsqrt(1-x*x)

c Compute the diagonal elements
        do l=1,lmax
         p(l,l) = (2*l-1)*p(l-1,l-1)*u
        end do

c Compute P_lm along the columns with fixed m
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

c Renormalize values
        do l=0,lmax
        mm = 1
        do m=0,l
         dnorm = fact(l-m)*(2*l+1)/(2*fact(l+m))
         p(l,m) = mm*twopinv*dsqrt(dnorm)*p(l,m)
         mm = -mm
        end do
        end do

        return
        end

c-----------------------------------------------------------------------
c Factorial Function
c Algorithm from http://www.livephysics.com/computational-physics/
c fortran/fortran-subroutines-functions/
c-----------------------------------------------------------------------

        real*8 function fact(n)
         implicit none

         real*8 n,i
         
         fact = 1

         do i = 2, n
            fact = fact * i
         end do
         return
        end function fact
