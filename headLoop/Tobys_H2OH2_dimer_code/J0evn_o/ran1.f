      function ran1(idum)
      
      implicit double precision (a-h,o-z)
      
      dimension r(97)
      
      parameter (M1=259200, IA1=7141, IC1=54773, RM1=3.8580247d-6)
      parameter (M2=134456, IA2=8121, IC2=28411, RM2=7.4373773d-6)
      parameter (M3=243000, IA3=4561, IC3=51349)
      
      save
      
      data iff/0/
      
      if (idum.lt.0.or.iff.eq.0) then
         iff = 1
         IX1 = mod(IC1-idum,M1)
         IX1 = mod(IA1*IX1+IC1,M1)
         IX2 = mod(IX1,M2)
         IX1 = mod(IA1*IX1+IC1,M1)
         IX3 = mod(IX1,M3)
         do 10 j=1,97
            IX1 = mod(IA1*IX1+IC1,M1)
            IX2 = mod(IA2*IX2+IC2,M2)
            r(j) = (dfloat(IX1)+dfloat(IX2)*RM2)*RM1
 10      continue
         idum = 1
      endif
      IX1 = mod(IA1*IX1+IC1,M1)
      IX2 = mod(IA2*IX2+IC2,M2)
      IX3 = mod(IA3*IX3+IC3,M3)
      j = 1+(97*IX3)/M3
      if (j.gt.97.or.j.lt.1) PAUSE
      ran1 = r(j)
      r(j) = (dfloat(IX1)+dfloat(IX2)*RM2)*RM1
      
      return
      end


