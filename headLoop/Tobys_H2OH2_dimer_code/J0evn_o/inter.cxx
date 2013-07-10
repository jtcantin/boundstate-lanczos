#include "inter.h"
double Interp::interp(const double &x)
{
 double retval;
 //
 // before list - extrapolate
 //
 if (x<=xlist(0)) {
   retval = ylist(0) + (x-xlist(0)) * 
     (ylist(1)-ylist(0)) / (xlist(1)-xlist(0));

   return retval;
 }

 //
 // after list - extrapolate
 //
 if (x>=xlist(cnt-1)) {
   retval = ylist(cnt-1) + (x-xlist(cnt-1)) * 
     (ylist(cnt-2)-ylist(cnt-1)) / (xlist(cnt-2)-xlist(cnt-1));

   return retval;
 }

 //
 // maintain static position counter - this makes this routine most efficient
 // if nearby positions are asked for in succession
 //

 //
 // get x between xlist(pos) and xlist(pos+1)
 //
 while (pos<cnt-2 && x>xlist(pos+1)) pos++;
 while (pos>0     && x<xlist(pos  )) pos--;

 retval = ylist(pos) + (x-xlist(pos)) * 
          (ylist(pos+1)-ylist(pos)) / (xlist(pos+1)-xlist(pos));

 return retval;
} /* Interp::interp */
double Interp2d::interp2d(const double &x,const double &y)
{
    int i;
  // figure out where the point is
  double dx=xv(1)-xv(0);
  double dy=yv(1)-yv(0);
  int ixl=(int) floor((x-xv(0))/dx);
  int ixh=(int) ceil((x-xv(0))/dx);
  int iyl=(int) floor((y-yv(0))/dy);
  int iyh=(int) ceil((y-yv(0))/dy);

  // cout<<ixl<<endl;

  if (ixl == ixh) ixh+=1;
  if (iyl == iyh) iyh+=1;
  
  // extrapolate
  if ((ixh > size-1 || ixl <0) || (iyh > size-1 || iyl <0)) {
    //cerr<<"extrapolation not implemented yet!!"<<endl;
    //exit(0);
    //should be close to zmax....for my case
    return zm(size-1,size-1);
  }

  vector subx(2);
  vector suby(2);
  for (i=0;i<2;i++){ 
    subx(i)=xv(ixl+i);
    suby(i)=yv(iyl+i);
  }

  vector subz(2);
  vector ztemp(2);
  for (i=0;i<2;i++){ 
    for (int j=0;j<2;j++)
      subz(j)=zm(ixl+i,iyl+j);
    Interp inty(2,suby,subz);
    ztemp(i)=inty.interp(y);
  }
  Interp intx(2,subx,ztemp);

  double retval=intx.interp(x);
  return retval;

}
