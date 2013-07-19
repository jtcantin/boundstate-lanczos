#ifndef INTERFACECONTAINERS_H
#define	INTERFACECONTAINERS_H

#include "vectClass.h" //This is to select which linear algebra class should be used.

struct lanczosStor {
	int total_basis_size;
	string sim_descr;
	string sim_descr_short;
};

struct generalStor {
	double dummy;
};

#endif
