#include "interface.h"


void HvInterfaceSetup(string HvCalculatorSwitch, generalStor **general_data, void (**HvPrepPtr)(int, char**, generalStor*, lanczosStor*), void (**HvPtr)(int, char**, generalStor*, lanczosStor*, double*, double*)){

	//Allocate transfer container and function pointers
	if (HvCalculatorSwitch == "linRotCartSph_Alavi_TIP4P") {
		interfaceStor *linRotCartSph_Alavi = new interfaceStor();
		
		//Potential function pointers 
		linRotCartSph_Alavi->fcnPointers = new fcnPointerStor();
		linRotCartSph_Alavi->fcnPointers->linearMoleculePotential = &Alavi_H2_Eng_Point;
		linRotCartSph_Alavi->fcnPointers->preCalcPotential = &preCalcPotential_Alavi_TIP4P;
		
		//Reinterpret cast for the sake of allowing C++ to pass it around without knowing what is inside
		(*general_data) = reinterpret_cast<generalStor*> (linRotCartSph_Alavi);
		
		// Lanczos Interface function pointers
		(*HvPrepPtr) = &Hv_Prep_linRotCartSph;
		(*HvPtr) = &Hv_linRotCartSph;
		
	}
	else {
		cerr << "HvCalculator '" << HvCalculatorSwitch << "' not recognized." << endl;
		exit(1);
	}
}

void Hv_Prep_linRotCartSph(int argc, char **argv, generalStor *general_data, lanczosStor *lanczos_data) {
	
	interfaceStor *Hv_data;
	Hv_data = reinterpret_cast<interfaceStor*> (general_data);
	
	HvPrep_Internal(argc, argv, Hv_data, lanczos_data);
	
};

void Hv_linRotCartSph(int argc, char **argv, generalStor *general_data, lanczosStor *lanczos_data, double *vec, double *uec) {
	int i;
	double *uec1;
	
	interfaceStor *Hv_data;
	Hv_data = reinterpret_cast<interfaceStor*> (general_data);
	
	uec1 = Hv_5D_oneCompositeIndex(Hv_data, vec);
	
	for (i=0; i<lanczos_data->total_basis_size; i++) {
		uec[i] += uec1[i];
	}
	
	delete [] uec1;
	
};