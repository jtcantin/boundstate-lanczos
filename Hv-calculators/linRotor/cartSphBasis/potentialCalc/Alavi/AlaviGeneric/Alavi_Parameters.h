#ifndef ALAVI_PARAMETERS_H
#define	ALAVI_PARAMETERS_H

//Alavi H2 parameters (NOTE: Lennard-Jones Parameters not included here as only mixed species are currently used and the LJ parameters are dependent on these other species)
#define Q_H2_H 0.4932 //e, from Alavi, S.; Ripmeester, J. A.; Klug, D. D. Molecular-dynamics study of structure II hydrogen clathrates. J. Chem. Phys. 2005, 123, 024507.
#define Q_H2_CM -0.9864 //e, from Alavi, S.; Ripmeester, J. A.; Klug, D. D. Molecular-dynamics study of structure II hydrogen clathrates. J. Chem. Phys. 2005, 123, 024507.

#define H2_H1_X 0.
#define H2_H1_Y 0.
#define H2_H1_Z 0.03707 //nm; Bondlength is 0.7414 Ang from Alavi_2005

#define H2_H2_X 0.
#define H2_H2_Y 0.
#define H2_H2_Z -0.03707 //nm; Bondlength is 0.7414 Ang from Alavi_2005 (Note: this is assumed to be equal to H2_H1_Z * -1.0 in the code)

#endif
