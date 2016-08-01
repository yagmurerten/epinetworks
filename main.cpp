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
#include "NetworkGenerator.h"
#include "parameters.h"
#include "pathogen.h"
#include "print.h"

#include "random.h"

#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <memory>

//using namespace epinetworks;

typedef std::vector<std::vector<int>> States;

int main(int, char *argv[]){
    double mutationFrac;
    double mutationRate = 0.;
    double var;
    double endtime;
    double recovery;
    double transmission;
    double virulence;
    
    if (epinetworks::DYNAMICS_TYPE == epinetworks::Dynamics::DynamicsType::SIR) {
        double fracTransmission;
        recovery = epinetworks::SIRparameters::RECOVERY;
        endtime = epinetworks::SIRparameters::ENDTIME;
        var = epinetworks::VARIANCE_K;
        if (epinetworks::INPUT_PARAMETERS) {
            fracTransmission = static_cast<double>(atoi(argv[1]));
            if (epinetworks::MUTATIONS)
                mutationFrac = static_cast<double>(atoi(argv[2]));
        }
        else {
            fracTransmission = epinetworks::SIRparameters::TRANSMISSION;
            if (epinetworks::MUTATIONS)
                mutationFrac = epinetworks::evoParameters::MUTATION_FRAC;
        }
        if (epinetworks::MORTALITY)
            virulence = epinetworks::SIRparameters::INITIAL_VIRULENCE;
        transmission = fracTransmission / 100.0;
    }

    if (epinetworks::DYNAMICS_TYPE == epinetworks::Dynamics::DynamicsType::SIS) {
        recovery = epinetworks::SISparameters::RECOVERY;
        endtime = epinetworks::SISparameters::ENDTIME;
        if (epinetworks::INPUT_PARAMETERS) {
            var = static_cast<double>(atoi(argv[1]));
            if (epinetworks::MUTATIONS)
                mutationFrac = static_cast<double>(atoi(argv[2]));
        }
        else {
            var = epinetworks::VARIANCE_K;
            if (epinetworks::MUTATIONS)
                mutationFrac = epinetworks::evoParameters::MUTATION_FRAC;
        }
        if (epinetworks::MORTALITY)
            virulence = epinetworks::SISparameters::INITIAL_VIRULENCE;
    }
    
    if (epinetworks::MUTATIONS)
        mutationRate = 1.0 / mutationFrac;

    int succesfulRuns = 0;
    int extinctionCount = 0;

    std::string parameterFileName = "parameters.txt";
    std::ofstream fileLog("outputLog.txt", std::ios_base::app);
    std::ofstream logParameters(parameterFileName, std::ios_base::app);
    
    logParameters << "mutation rate: " << mutationRate << std::endl;
    logParameters << "mutation sd: " << epinetworks::evoParameters::MUTATION_SD << std::endl;
    logParameters << "endtime: " << endtime << std::endl;
    logParameters << "network size: " << epinetworks::NETWORK_SIZE << std::endl;
    logParameters << "recovery rate: " << recovery << std::endl;
    if (epinetworks::DYNAMICS_TYPE == epinetworks::Dynamics::DynamicsType::SIR) {
        logParameters << "transmission: " << transmission << std::endl;
        logParameters << "R0: " << transmission / recovery << std::endl;
    }
    if (epinetworks::MORTALITY)
        logParameters << "initial virulence: " << virulence << std::endl;

	const std::vector<int> seed;
    epinetworks::RandomNumberGenerator rng = epinetworks::create_random_number_generator(seed);

	const std::string filenameNetwork = "networkNew";
    epinetworks::Network network(epinetworks::NETWORK_SIZE);

    if (epinetworks::OPTION_NETWORK_INPUT == false) {
        epinetworks::NetworkGenerator::networkWithoutInput(network, epinetworks::NETWORK_TYPE, var, rng, parameterFileName);
        if (epinetworks::NETWORK_TYPE != epinetworks::NetworkConstructor::NetworkType::FullyConnected){
            epinetworks::Print::networkOutput(filenameNetwork, network);
            epinetworks::Print::printNetwork(filenameNetwork, network);
            epinetworks::Print::printEdgeList(filenameNetwork, network);
        }
    }   

    if (epinetworks::OPTION_NETWORK_INPUT == true) {
		logParameters << "option network: true" << std::endl;
        epinetworks::NetworkGenerator::inputNetwork(network, filenameNetwork);
		DEBUG_ASSERT(network.isValid());
		std::string filenameNetwork2 = filenameNetwork + "2";
        epinetworks::Print::networkOutput(filenameNetwork2, network);
        epinetworks::Print::printNetwork(filenameNetwork2, network);
        epinetworks::Print::printEdgeList(filenameNetwork2, network);
	}

	for (std::size_t t = 0; t < network.size(); ++t) {
		network[t].setSusceptibleNumber();
	}

    if (epinetworks::NETWORK_DEBUG)
        return 10;
   
    std::unique_ptr<epinetworks::Dynamics> currentDynamics = epinetworks::Gillespie::createDynamics(epinetworks::DYNAMICS_TYPE);

    bool epiComplete = false;

    for (std::size_t replicate = 1u; replicate < epinetworks::NUMBER_OF_REPLICATES + 1; ++replicate) {
        const std::string filenameVirulence = "virulence" + std::to_string(replicate) + ".txt";
		const std::string filenameFinalSize = "finalsize.txt";
        const std::string filenameMaxIncidence = "maxIncidence.txt";
        const std::string filenameSnapshot = "degreeVir" + std::to_string(replicate) + ".csv";
     	std::ofstream finalSize(filenameFinalSize, std::ios_base::app);
        std::ofstream maxInc(filenameMaxIncidence, std::ios_base::app);

        States states;
        states.reserve(10000);
        std::vector<int> initialStates;
        initialStates.reserve(5);
        for (size_t i = 0; i < 5; ++i)
            initialStates.push_back(0);
        for (size_t i = 0; i < network.size(); ++i)
            states.push_back(initialStates);
        
		double t = 0.;
		double coefficient;
        epinetworks::assignCoefficient(coefficient, epinetworks::NETWORK_TYPE, epinetworks::DYNAMICS_TYPE);
        int i = epinetworks::getRandom(network.size(), rng);
        epinetworks::Individual &patientZero = network[i];
        if (epinetworks::DYNAMICS_TYPE == epinetworks::Dynamics::DynamicsType::SIR) {
            epinetworks::Pathogen initialPathogen(transmission*coefficient, recovery);
            patientZero.getInfected(initialPathogen);
        }
        else {
            epinetworks::Pathogen initialPathogen(virulence, coefficient, recovery);
            patientZero.getInfected(initialPathogen);
        }
        epinetworks::Infecteds infecteds(&patientZero, epinetworks::NETWORK_SIZE);
        patientZero.updateSusceptibleNeigbours(-1);
		bool endEpidemics = false;
		//std::ofstream fileExtinctions("extinctions.txt", std::ios_base::app);
        bool outPutTaken =0;
        if (epinetworks::MUTATIONS) {
            epinetworks::Print::virulenceOutput(filenameVirulence, 0, infecteds);
            epinetworks::Print::virulenceSnapShot(filenameSnapshot, 0, infecteds);
            outPutTaken = 1;
        }
        int maxIncidence = 0;
        double tau;
        double tPrev = 0.;
        double ratePrev = 0.;
		do {
            double random = -log(epinetworks::getRandomUniform(rng));
            double rateSum = epinetworks::Gillespie::rateSum(infecteds);
            tau = random / rateSum;
            if ((rateSum == 0 || tau > endtime) && (t > 1000)) {
               fileLog << "ERROR in calculating tau with # infecteds: " << infecteds.getSizeInfected();
               fileLog << " random was " << random << " rateSum was " << rateSum << std::endl;
               std::ofstream infectedRates("errorRatesum.csv", std::ios_base::app);
               double newRateSum = 0.;
               for (size_t i = 0; i < infecteds.getSizeInfected(); ++i) {
                   newRateSum += (infecteds.returnIndividual(i)).getEventRate();
                   infectedRates << newRateSum << "," << i << "," << infecteds.returnIndividual(i).getEventRate() << std::endl;
                   return 1;
               }
            }
            epinetworks::Gillespie::selectEvent(infecteds, rng, mutationRate, epinetworks::evoParameters::MUTATION_SD, currentDynamics, epinetworks::MUTATIONS, epinetworks::MORTALITY, states);
            tPrev = t;
			t += tau;
				
                if (epinetworks::MUTATIONS) {
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
                }
                else {
                    epinetworks::Print::virulenceOutput(filenameVirulence, t, infecteds);
                    epinetworks::Print::printStates(states, replicate);
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
                        if (finalRecovered > 100) {
                            epiComplete = true;
                        }
                    }
                    maxInc << replicate << "\t" << maxIncidence << std::endl;
					endEpidemics = true;
                    break;
				}

		} while (t < endtime && endEpidemics==false); 
        epinetworks::Print::virulenceOutput(filenameVirulence, static_cast<int>(round(t)), infecteds);
        epinetworks::Print::virulenceSnapShot(filenameSnapshot, static_cast<int>(round(t)), infecteds);
		for (std::size_t t = 0; t < network.size(); ++t) {
			network[t].getSusceptible();
            network[t].setSusceptibleNumber();
		}
        if (t >= endtime) {
            fileLog << "ended at: " << t << "last tau is " << tau << "previous t is " << tPrev << std::endl;
            ++succesfulRuns;
        }
        else
            ++extinctionCount;
        if (epinetworks::DYNAMICS_TYPE == epinetworks::Dynamics::DynamicsType::SIS && succesfulRuns >= 10) {
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


