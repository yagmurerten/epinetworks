#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "NetworkConstructor.h"
#include "Dynamics.h"
#include "print.h"

namespace epinetworks {


    //DEFAULT NETWORK OPTIONS
    static NetworkConstructor::NetworkType NETWORK_TYPE = NetworkConstructor::NetworkType::Homogeneous;
    static double VARIANCE_K = 4;
    static int NETWORK_SIZE = 10000;                 // # of individuals in the network
   
    //SIMULATION OPTIONS
    static const std::size_t NUMBER_OF_REPLICATES = 100;
    static const bool NETWORK_DEBUG = false;
    static Dynamics::DynamicsType DYNAMICS_TYPE = Dynamics::DynamicsType::SIR;
    const bool OPTION_NETWORK_INPUT = false;

	struct SIRparameters {
		static const double SIR_full_coef;      // SIR fully-connected
		static const double SIR_hom_coef;       // SIR k=4 (for 4-1=3)
        static const double INITIAL_VIRULENCE;
        static const double ENDTIME;
        static const double TRANSMISSION;
        static const double RECOVERY;
        static const bool INPUT_PARAMETERS;
        static const bool MUTATIONS;
        static const bool MORTALITY;
	};

	struct SISparameters {
		static const double SIS_full_coef;
		static const double SIS_hom_coef;
        static const double INITIAL_VIRULENCE;
        static const double ENDTIME;
        static const double RECOVERY;
        static const bool INPUT_PARAMETERS;
        static const bool MUTATIONS;
        static const bool MORTALITY;
	};

	struct evoParameters {
		static const double MUTATION_SD;     //s.d. in mutation distribution
		static const int MUTATION_FRAC;      // pre-defined mutation rate
	};

    const bool SIRparameters::INPUT_PARAMETERS = false;
    const bool SIRparameters::MUTATIONS = false;
    const bool SIRparameters::MORTALITY = false;

    const bool SISparameters::INPUT_PARAMETERS = true;
    const bool SISparameters::MUTATIONS = true;
    const bool SISparameters::MORTALITY = true;

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
    static void parametersNoInput(double &mutationFrac, double &var, double &endtime, double &recovery, 
        double &transmission, double &virulence, bool &mortality, 
        bool &mutations, Dynamics::DynamicsType dynamicsType, std::string parameterFileName) {

        if (dynamicsType == Dynamics::DynamicsType::SIR) {
            double fracTransmission;
            recovery = SIRparameters::RECOVERY;
            endtime = SIRparameters::ENDTIME;
            var = VARIANCE_K;
            mortality = SIRparameters::MORTALITY;
            mutations = SIRparameters::MUTATIONS;
            fracTransmission = SIRparameters::TRANSMISSION;
            if (mutations)
                mutationFrac = evoParameters::MUTATION_FRAC;
            if (mortality)
                virulence = SIRparameters::INITIAL_VIRULENCE;
            transmission = fracTransmission / 100.0;
        }

        if (dynamicsType == epinetworks::Dynamics::DynamicsType::SIS) {
            recovery = SISparameters::RECOVERY;
            endtime = SISparameters::ENDTIME;
            mortality = SISparameters::MORTALITY;
            mutations = SISparameters::MUTATIONS;
            var = epinetworks::VARIANCE_K;
            if (mutations)
                mutationFrac = epinetworks::evoParameters::MUTATION_FRAC;
            if (mortality)
                virulence = epinetworks::SISparameters::INITIAL_VIRULENCE;
        }

        Print::printParameter(endtime, "endtime", parameterFileName);
        Print::printParameter(epinetworks::NETWORK_SIZE, "network size", parameterFileName);
        Print::printParameter(recovery, "recovery rate", parameterFileName);
        if (epinetworks::DYNAMICS_TYPE == epinetworks::Dynamics::DynamicsType::SIR) {
            Print::printParameter(transmission, "transmission", parameterFileName);
            Print::printParameter(transmission / recovery, "R0", parameterFileName);
        }
        if (mortality)
            Print::printParameter(virulence, "initial virulence", parameterFileName);
        if (mutations)
            Print::printParameter(epinetworks::evoParameters::MUTATION_SD, "mutation sd", parameterFileName);
    }

    static void parametersInput(double &mutationFrac, double &var, double &endtime, double &recovery,
        double &transmission, double &virulence, bool &mortality, bool &mutations, 
        double argv1, double argv2, Dynamics::DynamicsType dynamicsType, std::string parameterFileName) {
        if (dynamicsType == Dynamics::DynamicsType::SIR) {
            double fracTransmission;
            recovery = SIRparameters::RECOVERY;
            endtime = SIRparameters::ENDTIME;
            var = VARIANCE_K;
            mortality = SIRparameters::MORTALITY;
            mutations = SIRparameters::MUTATIONS;
            fracTransmission = argv1;
            if (mutations)
                mutationFrac = argv2;
            if (mortality)
                virulence = SIRparameters::INITIAL_VIRULENCE;
            transmission = fracTransmission / 100.0;
        }

        if (dynamicsType == epinetworks::Dynamics::DynamicsType::SIS) {
            recovery = SISparameters::RECOVERY;
            endtime = SISparameters::ENDTIME;
            mortality = SISparameters::MORTALITY;
            mutations = SISparameters::MUTATIONS;
            var = argv1;
            if (mutations)
                mutationFrac = argv2;
            if (mortality)
                virulence = epinetworks::SISparameters::INITIAL_VIRULENCE;
        }

        Print::printParameter(endtime, "endtime", parameterFileName);
        Print::printParameter(epinetworks::NETWORK_SIZE, "network size", parameterFileName);
        Print::printParameter(recovery, "recovery rate", parameterFileName);
        if (epinetworks::DYNAMICS_TYPE == epinetworks::Dynamics::DynamicsType::SIR) {
            Print::printParameter(transmission, "transmission", parameterFileName);
            Print::printParameter(transmission / recovery, "R0", parameterFileName);
        }
        if (mortality)
            Print::printParameter(virulence, "initial virulence", parameterFileName);
        if (mutations)
            Print::printParameter(epinetworks::evoParameters::MUTATION_SD, "mutation sd", parameterFileName);
       
    }
}

#endif //PARAMETERS_H
