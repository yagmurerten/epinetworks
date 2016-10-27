#ifndef PARAMETERS_SIR_H
#define PARAMETERS_SIR_H

#include "NetworkConstructor.h"
#include "Dynamics.h"
#include "print.h"

namespace epinetworks {

    static Dynamics::DynamicsType DYNAMICS_TYPE = Dynamics::DynamicsType::SIR;

    //DEFAULT NETWORK OPTIONS
    static NetworkConstructor::NetworkType NETWORK_TYPE = NetworkConstructor::NetworkType::Homogeneous;
    static double VARIANCE_K = 4.;  // # of NEIGHBOURS
    static int NETWORK_SIZE = 10000;                 // # of individuals in the network
   
    //SIMULATION OPTIONS
    static const std::size_t NUMBER_OF_REPLICATES = 10;
    static const bool NETWORK_DEBUG = false;
    const bool OPTION_NETWORK_INPUT = false;

    const double full_coef = 0.0001;      // SIR fully-connected
	const double hom_coef = 0.5;       // SIR k=4 (for 4-1=3)
    const double INITIAL_VIRULENCE = 0.;
    const double ENDTIME = 1000;
    const double NO_INPUT_TRANSMISSION = 2.5;
    const double RECOVERY = 1.0;

    const bool INPUT_PARAMETERS = false;
    const bool MUTATIONS = false;
    const bool MORTALITY = false;

	struct evoParameters {
		static const double MUTATION_SD;     //s.d. in mutation distribution
		static const int MUTATION_FRAC;      // pre-defined mutation rate
	};


	//DEFAULT PATHOGEN PARAMETERS
	const double evoParameters::MUTATION_SD = 0.01;
	const int evoParameters::MUTATION_FRAC = 10;
    

	static void assignCoefficient(double &coef, NetworkConstructor::NetworkType networkType) {
			if (networkType == NetworkConstructor::NetworkType::FullyConnected)
				coef = full_coef;
			else
				coef = hom_coef;
		}

    static void parametersNoInput(double &mutationFrac, double &var, double &endtime, double &recovery, 
        double &transmission, double &virulence, bool &mortality, bool &mutations, std::string parameterFileName) {
            recovery = RECOVERY;
            endtime = ENDTIME;
            var = VARIANCE_K;
            mortality = MORTALITY;
            mutations = MUTATIONS;
            transmission = NO_INPUT_TRANSMISSION;
            if (mutations)
                mutationFrac = evoParameters::MUTATION_FRAC;
            if (mortality)
                virulence = INITIAL_VIRULENCE;

        Print::printParameter(endtime, "endtime", parameterFileName);
        Print::printParameter(epinetworks::NETWORK_SIZE, "network size", parameterFileName);
        Print::printParameter(recovery, "recovery rate", parameterFileName);
        Print::printParameter(transmission, "transmission", parameterFileName);
        Print::printParameter(transmission / recovery, "R0", parameterFileName);
        if (mortality)
            Print::printParameter(virulence, "initial virulence", parameterFileName);
        if (mutations)
            Print::printParameter(epinetworks::evoParameters::MUTATION_SD, "mutation sd", parameterFileName);
    }

    static void parametersInputTrans(double &mutationFrac, double &var, double &endtime, double &recovery,
        double &transmission, double &virulence, bool &mortality, bool &mutations, 
        double argv1, double argv2, std::string parameterFileName) {
            double fracTransmission;
            recovery = RECOVERY;
            endtime = ENDTIME;
            var = VARIANCE_K;
            mortality = MORTALITY;
            mutations = MUTATIONS;
            if (mortality)
                virulence = INITIAL_VIRULENCE;
            transmission = argv1 / 100.0;

        Print::printParameter(endtime, "endtime", parameterFileName);
        Print::printParameter(epinetworks::NETWORK_SIZE, "network size", parameterFileName);
        Print::printParameter(recovery, "recovery rate", parameterFileName);
        Print::printParameter(transmission, "transmission", parameterFileName);
        Print::printParameter(transmission / recovery, "R0", parameterFileName);
        if (mortality)
            Print::printParameter(virulence, "initial virulence", parameterFileName);
        if (mutations)
            Print::printParameter(epinetworks::evoParameters::MUTATION_SD, "mutation sd", parameterFileName);
       
    }
}

#endif //PARAMETERS_SIR_H
