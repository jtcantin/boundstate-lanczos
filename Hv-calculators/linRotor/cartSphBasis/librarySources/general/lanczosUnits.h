#ifndef LANCZOSUNITS_H
#define	LANCZOSUNITS_H

//Uses the same unit system as MMTK (length - nm, time - ps, energy - kJ/mol)
#define PI 3.141592653589793 //From NumPy (16 sig figs)
#define EPS0 0.0005727656384448188 //From MMTK (16 sig figs)
//#define K_CONST 138.93548461110962 // 1/(4*PI*EPS0) (17 sig figs)
#define K_CONST 138.9354846111096 // 1/(4*PI*EPS0) (16 sig figs)

#define CHUNK_RATIO 100 //How large the chunks should be relative to the total data size; 10 or 100 seem to be reasonable values.
#define NM_PER_BOHR 0.052918
#define ANG_PER_BOHR 0.52918
#define NM_PER_ANG 0.1



#endif