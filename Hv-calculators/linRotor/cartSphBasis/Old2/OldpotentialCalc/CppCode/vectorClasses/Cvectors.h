#ifndef CVECTORS_H
#define CVECTORS_H

// ND Cartesian vectors #######################################################
struct CNvect 
{
	//NOTE: NO cross product!
	double *coor;
	int num_dim;
	
	//Constructors
	CNvect ();
	CNvect (int);
	CNvect (double*, int);
	
	//Destructor
	//~CNvect();
	
	//Functions
	double eNorm () const;
	void printVec () const;
	void setCoor(double*, int);
	void setDim(int);
	
	//Overloaded Operators	 
	CNvect operator + (const CNvect&) const;
	CNvect operator - (const CNvect&) const;
	
	//Scalar Multiplication	
	friend CNvect operator * (const CNvect&, int); 
	friend CNvect operator * (int, const CNvect&); 
	
	
	double operator * (const CNvect&) const; //Dot product
	CNvect& operator = (const CNvect&);
	
};


// 3D Cartesian vectors #######################################################
struct C3vect 
{
	//NOTE: NO cross product!
	double x, y, z;
	C3vect () 
	{
		x=0; y=0; z=0;
	};
	
	C3vect (double, double, double);
	double eNorm () const;
	void printVec () const;
	void setCor(double, double, double);
	
	//Overloaded Operators	
	C3vect operator + (const C3vect&) const;
	C3vect operator - (const C3vect&) const;
	
	//Scalar Multiplication	
	friend C3vect operator * (const C3vect&, int); 
	friend C3vect operator * (int, const C3vect&); 
	
	
	double operator * (const C3vect&) const; //Dot product
	C3vect& operator = (const C3vect&);
	
};

// 2D Cartesian vectors #######################################################

struct C2vect 
{
	//NOTE: NO cross product!
	double x, y;
	C2vect () 
	{
		x=0; y=0;
	};
	
	C2vect (double, double);
	double eNorm () const;
	void printVec () const;
	void setCor(double, double);
	
	//Overloaded Operators	
	C2vect operator + (const C2vect&) const;
	C2vect operator - (const C2vect&) const;
	
	//Scalar Multiplication	
	friend C2vect operator * (const C2vect&, int); 
	friend C2vect operator * (int, const C2vect&); 
	
	
	double operator * (const C2vect&) const; //Dot product
	C2vect& operator = (const C2vect&);
	
};

#endif