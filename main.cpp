//////////////////////////////////////////////////////////////////////////////////
/////////// EpiNetworks - E. Yagmur Erten ////////////////////////////////////////
/////////// 08.06.2016 ///////////////////////////////////////////////////////////
/////////// Code for simulating epidemics (SIR and SIS) //////////////////////////
/////////// on networks //////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

#include "Debug.h"
#include "Gillespie.h"
#include "Dynamics.h"
#include "individual.h"
#include "infecteds.h"
#include "network.h"
#include "networkGenerator.h"
#include "pathogen.h"
#include "print.h"

#include "random.h"

#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <memory>

// determine which type of epidemic dynamics we are going to look at
#define _SIS

#ifdef _SIR
#include "parametersSIR.h"
#endif

#ifdef _SIS
#include "parametersSIS.h"
#endif

// this is something I added for criticality & information dynamics analysis
// a vector that keeps track of states of an individual and its neighbours
typedef std::vector<std::vector<int>> States;
const bool EPITRANS = false;

int main(int, char *argv[]){
    double mutationFrac = 0.;
    double mutationRate = 0.;
    double var = 0.;
    double endtime = 0.;
    double recovery = 0.;
    double transmission = 0.;
    double virulence = 0.;
    bool mortality = 0;
    bool mutations = 0;
	// these are the parameters that will be assigned later depending on the type of epidemics
	// and the parasite & dynamics we want to check out

    std::string parameterFileName = "parameters.txt";
    std::ofstream fileLog("outputLog.txt", std::ios_base::app);
    
	// if parameters are not input of the simulation, then they are hard-coded in headers
	// parametersSIS.h and parametersSIR.h
	// they are included above in this file
	// I need to write a parametersInput version of this for the evolution model
	if (!epinetworks::INPUT_PARAMETERS) {
		epinetworks::parametersNoInput(mutationFrac, var, endtime, recovery, transmission, virulence,
			mortality, mutations, parameterFileName);
	}
	else if (epinetworks::INPUT_PARAMETERS && EPITRANS) {
		double argv1 = static_cast<double>(atoi(argv[1]));
		epinetworks::parametersInputTrans(mutationFrac, var, endtime, recovery, transmission, virulence,
			mortality, mutations, argv1, parameterFileName);
	}
	else {
		double argv1 = static_cast<double>(atoi(argv[1]));
		double argv2 = static_cast<double>(atoi(argv[2]));
		epinetworks::parametersInputVir(mutationFrac, var, endtime, recovery,
			transmission, virulence, mortality, mutations,
			argv1, argv2, parameterFileName);
	}

	// if we allow virulence to evolve, mutation rate is determined here
    if (mutations)
        mutationRate = 1.0 / mutationFrac;

    epinetworks::Print::printParameter(mutationRate, "mutation rate", parameterFileName);

    int succesfulRuns = 0;
    int extinctionCount = 0; 

	const std::vector<int> seed;
    epinetworks::RandomNumberGenerator rng = epinetworks::create_random_number_generator(seed);

	const std::string filenameNetwork = "networkNew";
    epinetworks::Network network(epinetworks::NETWORK_SIZE);

	// if we want to construct a network from scratch we use this code
    if (!epinetworks::OPTION_NETWORK_INPUT) {
        epinetworks::NetworkGenerator::networkWithoutInput(network, epinetworks::NETWORK_TYPE, var, rng, parameterFileName);
        if (epinetworks::NETWORK_TYPE != epinetworks::NetworkConstructor::NetworkType::FullyConnected){
            epinetworks::Print::networkOutput(filenameNetwork, network);
            epinetworks::Print::printNetwork(filenameNetwork, network);
            epinetworks::Print::printEdgeList(filenameNetwork, network);
        }
    }   

	// if we want to input the same network from a file we use this code
    if (epinetworks::OPTION_NETWORK_INPUT) {
        epinetworks::Print::printParameter("true", "option network ", parameterFileName);
        epinetworks::NetworkGenerator::inputNetwork(network, filenameNetwork);
		DEBUG_ASSERT(network.isValid());
		std::string filenameNetwork2 = filenameNetwork + "2";
        epinetworks::Print::networkOutput(filenameNetwork2, network);
        epinetworks::Print::printNetwork(filenameNetwork2, network);
        epinetworks::Print::printEdgeList(filenameNetwork2, network);
	}

	// setting neighbours to be susceptible, this is required for the epitrans version
	for (std::size_t t = 0; t < network.size(); ++t) {
		network[t].setSusceptibleNumber();
	}

    if (epinetworks::NETWORK_DEBUG)
        return 10;
   
    std::unique_ptr<epinetworks::Dynamics> currentDynamics = epinetworks::Gillespie::createDynamics(epinetworks::DYNAMICS_TYPE);

	// for SIR dynamics
    bool epiComplete = false;

    for (std::size_t replicate = 1u; replicate < epinetworks::NUMBER_OF_REPLICATES + 1; ++replicate) {
        const std::string filenameVirulence = "virulence" + std::to_string(replicate) + ".txt";
		const std::string filenameFinalSize = "finalsize.txt";
        const std::string filenameMaxIncidence = "maxIncidence.txt";
        const std::string filenameSnapshot = "degreeVir" + std::to_string(replicate) + ".csv";
     	std::ofstream finalSize(filenameFinalSize, std::ios_base::app);
        std::ofstream maxInc(filenameMaxIncidence, std::ios_base::app);
        
		double t = 0.;
		double coefficient;
        epinetworks::assignCoefficient(coefficient, epinetworks::NETWORK_TYPE);
        int i = epinetworks::getRandom(network.size(), rng);
        epinetworks::Individual &patientZero = network[i];
        if ((epinetworks::DYNAMICS_TYPE == epinetworks::Dynamics::DynamicsType::SIR) || 
            (epinetworks::DYNAMICS_TYPE == epinetworks::Dynamics::DynamicsType::SIS && virulence == 0.)) {
            // version with SIR dynamics or no virulence
			epinetworks::Pathogen initialPathogen(transmission*coefficient, recovery);
            patientZero.getInfected(initialPathogen);
        } else {
			// version with virulence-transmission trade-off & evolution
            epinetworks::Pathogen initialPathogen(virulence, coefficient, recovery);
            patientZero.getInfected(initialPathogen);
        }
       
        epinetworks::Infecteds infecteds(&patientZero, epinetworks::NETWORK_SIZE);
        patientZero.updateSusceptibleNeigbours(-1); 
		// one less susceptible neighbour for the patient zero's neighbours
		bool endEpidemics = false;
		//std::ofstream fileExtinctions("extinctions.txt", std::ios_base::app);
        bool outPutTaken =0;
        States states;
        if (mutations) {
            epinetworks::Print::virulenceOutput(filenameVirulence, 0, infecteds);
            epinetworks::Print::virulenceSnapShot(filenameSnapshot, 0, infecteds);
            outPutTaken = 1;
        }
        int maxIncidence = 0;
        double tau;
        double tprev = 0.;
        int prevTime = -1;
		do {
            double randDenom;
            do {
                randDenom = epinetworks::getRandomUniform(rng);
            } while (randDenom == 0); // to make sure that random is not divided by 0
            double random = -log(randDenom);
            double rateSum = epinetworks::Gillespie::rateSum(infecteds);
            tau = random / rateSum;
            epinetworks::Gillespie::selectEvent(infecteds, rng, mutationRate, epinetworks::evoParameters::MUTATION_SD, currentDynamics, mutations, mortality, states);
            tprev = t;
			t += tau;
				
                if (mutations) {
                    if (static_cast<int>(round(t)) % 5 == 0) {
                        if (outPutTaken == 0) {
                            epinetworks::Print::virulenceOutput(filenameVirulence, static_cast<int>(round(t)), infecteds);
                            epinetworks::Print::virulenceSnapShot(filenameSnapshot, static_cast<int>(round(t)), infecteds);
                            outPutTaken = 1;
                        }
                    }

                    if (static_cast<int>(round(t)) % 5 == 1) {
                        outPutTaken = 0;
                    }
                } else {
                    //epinetworks::Print::printStates(t, states, replicate);
					// implicitly, this is the 'epitrans' version, I need to make it clearer
                    int time = floor(t);
                    if (time != prevTime) {
                        epinetworks::Print::virulenceOutput(filenameVirulence, t, infecteds);
                        epinetworks::Print::printStatesAll(t, time, replicate, network);
                        prevTime = time;
                    }

                }
				
                int infectedSize = infecteds.getSizeInfected();
                if (infectedSize > maxIncidence) {
                    maxIncidence = infectedSize;
                }

				if (infecteds.getSizeInfected() == 0) {
                    if (epinetworks::DYNAMICS_TYPE == epinetworks::Dynamics::DynamicsType::SIR)
                        fileLog << replicate << "\t" << infecteds.getSizeInfected() << "\t" << t << std::endl;
                    else
                        fileLog << replicate <<  " ended with " << infecteds.getSizeInfected() << " infecteds at t: " << t << std::endl;
					//std::cout << "end of epidemics" << std::endl;
                    if (epinetworks::DYNAMICS_TYPE == epinetworks::Dynamics::DynamicsType::SIR) {
                        int finalRecovered = 0;
                        for (std::size_t i = 0u; i < network.size(); ++i) {
                            if ((network[i]).getStatus() == epinetworks::Individual::Status::recovered)
                                ++finalRecovered;
                        }
                        finalSize << replicate << "\t" << finalRecovered << std::endl;
                        if (finalRecovered > 2000) {
                            epiComplete = true;
                        }
                    }
                    maxInc << replicate << "\t" << maxIncidence << std::endl;
					endEpidemics = true;
                    break;
				}

		} while (t < endtime && endEpidemics==false);
        if (epinetworks::DYNAMICS_TYPE == epinetworks::Dynamics::DynamicsType::SIS) {
            epinetworks::Print::virulenceOutput(filenameVirulence, static_cast<int>(round(t)), infecteds);
            epinetworks::Print::virulenceSnapShot(filenameSnapshot, static_cast<int>(round(t)), infecteds);
        }
        if (t >= endtime) {
            fileLog << "replicate " << replicate << " ended at: " << t << " last tau is " << tau;
            fileLog << " prev t is: " << tprev << " infected size is " << infecteds.getSizeInfected() << std::endl;
            ++succesfulRuns;
        } else
            ++extinctionCount;
        for (std::size_t t = 0; t < network.size(); ++t) {
            network[t].getSusceptible();
            network[t].setSusceptibleNumber();
        }
        if (epinetworks::DYNAMICS_TYPE == epinetworks::Dynamics::DynamicsType::SIS && succesfulRuns >= 10 /*normally 10*/) {
            break;
        }
        if (epinetworks::DYNAMICS_TYPE == epinetworks::Dynamics::DynamicsType::SIR && epiComplete == true) {
            fileLog << "complete run: " << replicate << std::endl;
            ++succesfulRuns;
            break;
        }
	}
    fileLog << "succesful runs: " << succesfulRuns << std::endl;
    fileLog << "extinction count: " << extinctionCount << std::endl;
}


