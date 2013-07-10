class Rand {
 public:
    Rand() {}
    Rand(const unsigned seed) {srandom(seed);}
    void reseed(const unsigned seed) {srandom(seed);}
    
    int irange(const int brange,
	       const int erange);
    double frand();
    unsigned long irand();
};
static const unsigned long RANDMAX = 2147483647;
int Rand::irange(const int brange,
                 const int erange)
    //
    // generate a random integer in the range brange..erange inclusive
    //
{
    int ret = (int) ((erange-brange+1.0)*random()/(RANDMAX+1.0)) + (brange);
    //(random()%(erange-brange+1));
    return ret;
} /*Rand::irange*/
double Rand::frand()
    //return a float in the range [0..1)
{
    // static const unsigned long RAND_MAX = (unsigned long)2<<30;
    unsigned long rand = random();
    
    return (double)rand / (RANDMAX+1.0);
} /*Rand::frand*/

unsigned long Rand::irand()
    //return a random integer
{
    return random();
} /*Rand::irand*/
