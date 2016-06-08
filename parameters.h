#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "networkGenerator.h"
#include "Dynamics.h"

namespace epinetworks {

    //INPUT OPTIONS 
    static const bool OPTION_NETWORK_INPUT = 0;     // 0 = no input
    static const bool INPUT_PARAMETERS = 1;         // 0 = no input
    static const bool MUTATIONS = false;

	struct SIRparameters {
		static const double SIR_full_coef;   // SIR fully-connected
		static const double SIR_hom_coef;       // SIR k=4 (for 4-1=3)
	};

	struct SISparameters {
		static const double SIS_full_coef;
		static const double SIS_hom_coef;		// Mean k = 4
	};

	struct evoParameters {
		static const double MUTATION_SD;   //s.d. in mutation distribution
		static const int MUTATION_FRAC;      // pre-defined mutation rate
	};
    
	//DEFAULT NETWORK OPTIONS
	static const NetworkGenerator::NetworkType NETWORK_TYPE = NetworkGenerator::NetworkType::Homogeneous;
	static double VARIANCE_K = 4;
	static int NETWORK_SIZE = 10000;                 // # of individuals in the network

	//DEFAULT PATHOGEN PARAMETERS
	static double INITIAL_VIRULENCE = 0.3;          // initial value of virulence
	static double INITIAL_TRANSMISSION = 110;			// initial value of virulenc
	static const double RECOVERY = 1;				// Recovery
	const double evoParameters::MUTATION_SD = 0.01;
	const int evoParameters::MUTATION_FRAC = 0;

	//SIMULATION OPTIONS
	static const std::size_t NUMBER_OF_REPLICATES = 100;   // # of replicates on each network
	static const int ENDTIME=1000;                              // # of time steps

	static const Dynamics::DynamicsType DYNAMICS_TYPE = Dynamics::DynamicsType::SIR;
	const double SIRparameters::SIR_full_coef = 0.0001;
	const double SIRparameters::SIR_hom_coef = 0.5;
	const double SISparameters::SIS_full_coef = 0.003;
	const double SISparameters::SIS_hom_coef = 0.75;

	static void assignCoefficient(double &coef, NetworkGenerator::NetworkType networkType,
		Dynamics::DynamicsType dynamicsType) {
		if (dynamicsType == Dynamics::DynamicsType::SIS) {
			if (networkType == NetworkGenerator::NetworkType::FullyConnected)
				coef = SISparameters::SIS_full_coef;
			else
				coef = SISparameters::SIS_hom_coef;
		}
		else
		if (networkType == NetworkGenerator::NetworkType::FullyConnected)
			coef = SIRparameters::SIR_full_coef;
		else
			coef = SIRparameters::SIR_hom_coef;
	}
}

#endif //PARAMETERS_H
