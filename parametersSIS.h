#ifndef PARAMETERS_SIS_H
#define PARAMETERS_SIS_H

#include "NetworkConstructor.h"
#include "Dynamics.h"
#include "print.h"

#define _EPITRANS

namespace epinetworks {

    static Dynamics::DynamicsType DYNAMICS_TYPE = Dynamics::DynamicsType::SIS;
    static int NETWORK_SIZE = 10000;                 // # of individuals in the network
   
    //SIMULATION OPTIONS
    static const bool NETWORK_DEBUG = false;
    const bool OPTION_NETWORK_INPUT = true;

    const double full_coef = 0.003;
    const double hom_coef = 0.75;

#ifdef _EPITRANS
        static NetworkConstructor::NetworkType NETWORK_TYPE = NetworkConstructor::NetworkType::Homogeneous;
        static const std::size_t NUMBER_OF_REPLICATES = 10;
        static double VARIANCE_K = 4.;
        const double INITIAL_VIRULENCE = 0.;
        const double ENDTIME = 100;
        const double RECOVERY = 1.0;
        const bool INPUT_PARAMETERS = true;
        const bool MUTATIONS = false;
        const bool MORTALITY = false;
        const double NO_INPUT_TRANSMISSION = 2.5;
#endif

#ifdef _EPIVIR
        static NetworkConstructor::NetworkType NETWORK_TYPE = NetworkConstructor::NetworkType::Gamma;
        static const std::size_t NUMBER_OF_REPLICATES = 10000;
        static double VARIANCE_K = 0.;
        const double INITIAL_VIRULENCE = 0.4;
        const double ENDTIME = 50000.;
        const double RECOVERY = 0.5;
        const bool INPUT_PARAMETERS = true;
        const bool MUTATIONS = true;
        const bool MORTALITY = true;
#endif

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
        double &transmission, double &virulence, bool &mortality, 
        bool &mutations, std::string parameterFileName) {

            recovery = RECOVERY;
            endtime = ENDTIME;
            mortality = MORTALITY;
            mutations = MUTATIONS;
            transmission = NO_INPUT_TRANSMISSION;
            var = VARIANCE_K;
            if (mutations)
                mutationFrac = evoParameters::MUTATION_FRAC;
            if (mortality)
                virulence = INITIAL_VIRULENCE;

        Print::printParameter(endtime, "endtime", parameterFileName);
        Print::printParameter(NETWORK_SIZE, "network size", parameterFileName);
        Print::printParameter(recovery, "recovery rate", parameterFileName);
        if (mortality)
            Print::printParameter(virulence, "initial virulence", parameterFileName);
        if (mutations)
            Print::printParameter(evoParameters::MUTATION_SD, "mutation sd", parameterFileName);
    }

    static void parametersInputVir(double &mutationFrac, double &var, double &endtime, double &recovery,
        double &transmission, double &virulence, bool &mortality, bool &mutations, 
        double argv1, double argv2, std::string parameterFileName) {

            recovery = RECOVERY;
            endtime = ENDTIME;
            mortality = MORTALITY;
            mutations = MUTATIONS;
            var = argv1;
            if (mutations)
                mutationFrac = argv2;
            if (mortality)
                virulence = INITIAL_VIRULENCE;

        Print::printParameter(endtime, "endtime", parameterFileName);
        Print::printParameter(NETWORK_SIZE, "network size", parameterFileName);
        Print::printParameter(recovery, "recovery rate", parameterFileName);
        if (mortality)
            Print::printParameter(virulence, "initial virulence", parameterFileName);
        if (mutations)
            Print::printParameter(evoParameters::MUTATION_SD, "mutation sd", parameterFileName);
       
    }

    static void parametersInputTrans(double &mutationFrac, double &var, double &endtime, double &recovery,
        double &transmission, double &virulence, bool &mortality, bool &mutations,
        double argv1, std::string parameterFileName) {

        recovery = RECOVERY;
        endtime = ENDTIME;
        mortality = MORTALITY;
        mutations = MUTATIONS;
        var = VARIANCE_K;
        transmission = argv1 / 100;
        Print::printParameter(endtime, "endtime", parameterFileName);
        Print::printParameter(NETWORK_SIZE, "network size", parameterFileName);
        Print::printParameter(recovery, "recovery rate", parameterFileName);
        if (mortality)
            Print::printParameter(virulence, "initial virulence", parameterFileName);
        if (mutations)
            Print::printParameter(evoParameters::MUTATION_SD, "mutation sd", parameterFileName);

    }
}

#endif //PARAMETERS_H
