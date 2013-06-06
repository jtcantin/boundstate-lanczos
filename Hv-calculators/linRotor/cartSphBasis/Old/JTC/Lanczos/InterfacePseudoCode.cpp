struct AlaviContainer {
	int i;
	double *potential;
}

struct genType {
	double genDouble;
}

void AlaviPrepare(int argc, char **argv, genType *container){
	AlaviContainer *MyCont;
	
	//Copy pointer value over for use as the desired structure
	MyCont = reinterpret_cast <AlaviContainer> container;
	
	//Use contents of MyCont
	MyCont->potential = new double [MyCont->i];
	
	//Continue on with C++ code
	for (j=0; j<i; j++) {
		MyCont->potential[j] += CoulombPotential();
	}
	
}

void AlaviHv(int argc, char **argv, genType *container){
	AlaviContainer *MyCont;
	
	//Copy pointer value over for use as the desired structure
	MyCont = reinterpret_cast <AlaviContainer> container;
	
	//Use contents of MyCont
	MyCont->potential = new double [MyCont->i];
	
	//Continue on with Fortran code
	FORTRAN(HvProd) (MyCont->potential, &(MyCont->i));
	
}

//Start of Lanczos Code
int main(int argc, char** argc){
	
	//Declare generic container pointer
	
	genType *container;
	
	//Declare function pointer
	void (*Prepare)(int, char **, genType) = NULL;
	void (*Hv)(int, char **, genType) = NULL;
	
	//Determine which functions to use based on commandline argument:
	switch (argv[3]) {
		case "Alavi":
			
			//Assign function pointers
			Prepare = &AlaviPrepare;
			Hv = &AlaviHv;
			
			//Declare and assign specific container
			AlaviContainter *AlCont;

			AlCont = new AlaviContainter;
			container = reinterpret_cast <genType> AlCont;
			
			
			break;
		default:
			cerr << "Not recognized" < endl;
			break;
	}

	
	//Prepare for Lanczos diagonalization
	(*Prepare) (argc, argv, container);
	
	//Lanczos first loop
	(*Hv) (argc, argv, container);
	
	//Lanczos second loop
	(*Hv) (argc, argv, container);
	
}

