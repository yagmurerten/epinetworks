#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "networkGenerator.h"
#include "Gillespie.h"

namespace epinetworks {

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

	//INPUT OPTIONS 
	static const bool OPTION_NETWORK_INPUT = 0;     // 0 = no input
	static const bool OPTION_VIRULENCE_INPUT = 0;   // 0 = no input

	//NETWORK OPTIONS
	static const NetworkGenerator::NetworkType NETWORK_TYPE = NetworkGenerator::NetworkType::Homogeneous;
	static double VARIANCE_K = 4;
	static int NETWORK_SIZE = 10000;                 // # of individuals in the network

	//PATHOGEN PARAMETERS
	static double INITIAL_VIRULENCE = 0.25;          // initial value of virulence
	static double INITIAL_TRANSMISSION = 1.1;			// initial value of virulenc
	static const double RECOVERY = 1;				// Recovery
	const double evoParameters::MUTATION_SD = 0.01;
	const int evoParameters::MUTATION_FRAC = 0;

	//SIMULARION OPTIONS
	static const std::size_t NUMBER_OF_REPLICATES = 100;   // # of replicates on each network
	static const int ENDTIME = 1500;                 // # of time steps

	//PARAMETERS FOR GILLESPIE
	static const Gillespie::DynamicsType DYNAMICS_TYPE = Gillespie::DynamicsType::SIR;
	const double SIRparameters::SIR_full_coef = 0.0001;
	const double SIRparameters::SIR_hom_coef = 0.5;
	const double SISparameters::SIS_full_coef = 0.003;
	const double SISparameters::SIS_hom_coef = 0.75;

	static void assignCoefficient(double &coef, NetworkGenerator::NetworkType networkType,
		Gillespie::DynamicsType dynamicsType) {
		if (dynamicsType == Gillespie::DynamicsType::SIS) {
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
