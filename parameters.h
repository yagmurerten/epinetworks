#ifndef PARAMETERS_H
#define PARAMETERS_H


//PARAMETERS FOR GILLESPIE
//const double MUTATION_SD = 0.05;			//s.d. in mutation distribution
const double MUTATION_SD = 0.05;


//INPUT OPTIONS 
const bool OPTION_NETWORK_INPUT = 0;		// 0 = no input; 1 = input file to be read to construct network 
const bool OPTION_VIRULENCE_INPUT = 0;		// 0 = no input; 1 = parameters to be read

//PARAMETERS FOR PATHOGEN
//const double C = 0.0003;	// Fully Connected (for R0 of 2)
const double C = 0.003;
const double D = 0.75;		// Mean k = 4

const double E = 0.0001;		// SIR


const double RECOVERY = 1;				// Recovery


#endif //PARAMETERS_H
