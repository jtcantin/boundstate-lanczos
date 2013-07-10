/* perform cubic spline on input data 
   and spit out the requested number of points sampled uniformly in the 
   input range.
 */
#include "matvec.h"
class Interp {
  int cnt;
  vector xlist;
  vector ylist;
  int pos;

public:
  Interp(int cnt,
	 const vector &xlist,
	 const vector &ylist) : xlist(cnt), ylist(cnt)
    {this->cnt=cnt; pos=0; this->xlist=xlist, this->ylist=ylist;}
  
  double interp(const double &x);
};

class Interp2d {
  int size;
  vector xv;
  vector yv;
  matrix zm;
  int posx;
  int posy;

public:
  Interp2d(int size,const vector &xv, const vector &yv,const matrix &zm) 
    : xv(size), yv(size), zm(size,size)
    {
      this->size=size; 
      posx=0; posy=0;
      this->xv=xv, this->yv=yv, this->zm=zm;
    }
  double interp2d(const double &x,const double &y);
};

