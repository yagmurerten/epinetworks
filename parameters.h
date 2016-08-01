#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "NetworkConstructor.h"
#include "Dynamics.h"

namespace epinetworks {

    //INPUT OPTIONS 
    static bool OPTION_NETWORK_INPUT = false;     
    static bool INPUT_PARAMETERS = true;         
    static bool MUTATIONS = true;
    static bool MORTALITY = true;
    static bool NETWORK_DEBUG = false;
    //SIMULATION OPTIONS
    static const std::size_t NUMBER_OF_REPLICATES = 100;
    static Dynamics::DynamicsType DYNAMICS_TYPE = Dynamics::DynamicsType::SIS;

	struct SIRparameters {
		static const double SIR_full_coef;      // SIR fully-connected
		static const double SIR_hom_coef;       // SIR k=4 (for 4-1=3)
        static const double INITIAL_VIRULENCE;
        static const double ENDTIME;
        static const double TRANSMISSION;
        static const double RECOVERY;
	};

	struct SISparameters {
		static const double SIS_full_coef;
		static const double SIS_hom_coef;
        static const double INITIAL_VIRULENCE;
        static const double ENDTIME;
        static const double RECOVERY;
	};

	struct evoParameters {
		static const double MUTATION_SD;     //s.d. in mutation distribution
		static const int MUTATION_FRAC;      // pre-defined mutation rate
	};
    
	//DEFAULT NETWORK OPTIONS
	static NetworkConstructor::NetworkType NETWORK_TYPE = NetworkConstructor::NetworkType::Gamma;
	static double VARIANCE_K = 4;
	static int NETWORK_SIZE = 10000;                 // # of individuals in the network

	//DEFAULT PATHOGEN PARAMETERS
	const double evoParameters::MUTATION_SD = 0.01;
	const int evoParameters::MUTATION_FRAC = 10;
    
	const double SIRparameters::SIR_full_coef = 0.0001;
	const double SIRparameters::SIR_hom_coef = 0.5;
	const double SISparameters::SIS_full_coef = 0.003;
	const double SISparameters::SIS_hom_coef = 0.75;

    const double SIRparameters::INITIAL_VIRULENCE = 0.;
    const double SISparameters::INITIAL_VIRULENCE = 0.35;

    const double SIRparameters::ENDTIME = 1000;
    const double SISparameters::ENDTIME = 25000;

    const double SIRparameters::TRANSMISSION = 300;

    const double SIRparameters::RECOVERY = 1.0;
    const double SISparameters::RECOVERY = 0.5;


	static void assignCoefficient(double &coef, NetworkConstructor::NetworkType networkType,
		Dynamics::DynamicsType dynamicsType) {
		if (dynamicsType == Dynamics::DynamicsType::SIS) {
			if (networkType == NetworkConstructor::NetworkType::FullyConnected)
				coef = SISparameters::SIS_full_coef;
			else
				coef = SISparameters::SIS_hom_coef;
		}
		else
		if (networkType == NetworkConstructor::NetworkType::FullyConnected)
			coef = SIRparameters::SIR_full_coef;
		else
			coef = SIRparameters::SIR_hom_coef;
	}
}

#endif //PARAMETERS_H
