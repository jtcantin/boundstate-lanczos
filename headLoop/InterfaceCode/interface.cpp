#include "interface.h"

void HvInterfaceSetup(string HvCalculatorSwitch, generalStor **general_data, void (**HvPrepPtr)(int, char**, generalStor*, lanczosStor*), void (**HvPtr)(int, char**, generalStor*, lanczosStor*, double*, double*)){

	//Allocate transfer container and function pointers
	if (HvCalculatorSwitch == "linRotCartSph_Alavi_TIP4P") {
		interfaceStor *linRotCartSph_Alavi = new interfaceStor();
		
		//Potential function pointers 
		linRotCartSph_Alavi->fcnPointers = new fcnPointerStor();
		
		//Alavi Model-specific function pointers
		linRotCartSph_Alavi->fcnPointers->linearMoleculePotential = &Alavi_H2_Eng_Point;
		linRotCartSph_Alavi->fcnPointers->preCalcPotential = &preCalcPotential_Alavi;
		
		//TIP4P Model-specific function pointers
		linRotCartSph_Alavi->fcnPointers->SummedCoulombPotential = &Q_TIP4P_Eng;
		linRotCartSph_Alavi->fcnPointers->SummedLJPotential = &LJ_TIP4P_Eng;
		linRotCartSph_Alavi->fcnPointers->SummedLJPotentialFast = &LJ_TIP4P_Eng_Fast;
		linRotCartSph_Alavi->fcnPointers->getAtomGeometry = &getTIP4Patoms;
		
		//Reinterpret cast for the sake of allowing C++ to pass it around without knowing what is inside
		(*general_data) = reinterpret_cast<generalStor*> (linRotCartSph_Alavi);
		
		//Hv-specific function pointers
		(*HvPrepPtr) = &Hv_Prep_linRotCartSph;
		(*HvPtr) = &Hv_linRotCartSph;
		
	}
	else if (HvCalculatorSwitch == "linRotCartSph_Alavi_SPCE") {
		interfaceStor *linRotCartSph_Alavi = new interfaceStor();
		
		//Potential function pointers 
		linRotCartSph_Alavi->fcnPointers = new fcnPointerStor();
		
		//Alavi Model-specific function pointers
		linRotCartSph_Alavi->fcnPointers->linearMoleculePotential = &Alavi_H2_Eng_Point;
		linRotCartSph_Alavi->fcnPointers->preCalcPotential = &preCalcPotential_Alavi;
		
		//TIP4P Model-specific function pointers
		linRotCartSph_Alavi->fcnPointers->SummedCoulombPotential = &Q_SPCE_Eng;
		linRotCartSph_Alavi->fcnPointers->SummedLJPotential = &LJ_SPCE_Eng;
		linRotCartSph_Alavi->fcnPointers->SummedLJPotentialFast = &LJ_SPCE_Eng_Fast;
		linRotCartSph_Alavi->fcnPointers->getAtomGeometry = &getSPCEatoms;
		
		//Reinterpret cast for the sake of allowing C++ to pass it around without knowing what is inside
		(*general_data) = reinterpret_cast<generalStor*> (linRotCartSph_Alavi);
		
		//Hv-specific function pointers
		(*HvPrepPtr) = &Hv_Prep_linRotCartSph;
		(*HvPtr) = &Hv_linRotCartSph;
		
	}
	else if (HvCalculatorSwitch == "linRotCartSph_Coulomb") {
		interfaceStor *linRotCartSph = new interfaceStor();
		
		//Potential function pointers 
		linRotCartSph->fcnPointers = new fcnPointerStor();
		
		//Coulomb-specific function pointers
		linRotCartSph->fcnPointers->linearMoleculePotential = &CoulombPotential;
		linRotCartSph->fcnPointers->preCalcPotential = &preCalcPotential_Coulomb;
		
		//Reinterpret cast for the sake of allowing C++ to pass it around without knowing what is inside
		(*general_data) = reinterpret_cast<generalStor*> (linRotCartSph);
		
		//Hv-specific function pointers
		(*HvPrepPtr) = &Hv_Prep_linRotCartSph;
		(*HvPtr) = &Hv_linRotCartSph;
		
	}
	else if (HvCalculatorSwitch == "linRotCartSph_IsoHarmOsc") {
		interfaceStor *linRotCartSph = new interfaceStor();
		
		//Potential function pointers 
		linRotCartSph->fcnPointers = new fcnPointerStor();
		
		//Coulomb-specific function pointers
		linRotCartSph->fcnPointers->linearMoleculePotential = &IsoHarmOscEng;
		linRotCartSph->fcnPointers->preCalcPotential = &preCalcPotential_IsoHarmOsc;
		
		//Reinterpret cast for the sake of allowing C++ to pass it around without knowing what is inside
		(*general_data) = reinterpret_cast<generalStor*> (linRotCartSph);
		
		//Hv-specific function pointers
		(*HvPrepPtr) = &Hv_Prep_linRotCartSph;
		(*HvPtr) = &Hv_linRotCartSph;
		
	}
	else {
		cerr << "HvCalculator '" << HvCalculatorSwitch << "' not recognized." << endl;
		exit(1);
	}
};

