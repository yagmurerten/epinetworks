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

using namespace epinetworks;

int main(int, char *argv[]){
    double mutationFrac;
    double mutationRate = 0.;
    double var;
    double endtime;
    double recovery;
    double transmission;
    double virulence;
    
    if (DYNAMICS_TYPE == Dynamics::DynamicsType::SIR) {
        double fracTransmission;
        recovery = SIRparameters::RECOVERY;
        endtime = SIRparameters::ENDTIME;
        var = VARIANCE_K;
        if (INPUT_PARAMETERS) {
            fracTransmission = static_cast<double>(atoi(argv[1]));
            if (MUTATIONS)
                mutationFrac = static_cast<double>(atoi(argv[2]));
        }
        else {
            fracTransmission = SIRparameters::TRANSMISSION;
            if (MUTATIONS)
                mutationFrac = evoParameters::MUTATION_FRAC;
            if (MORTALITY)
                virulence = SIRparameters::INITIAL_VIRULENCE;
        }
        transmission = fracTransmission / 100.0;
    }

    if (DYNAMICS_TYPE == Dynamics::DynamicsType::SIS) {
        recovery = SISparameters::RECOVERY;
        endtime = SISparameters::ENDTIME;
        if (INPUT_PARAMETERS) {
            var = static_cast<double>(atoi(argv[1]));
            if (MUTATIONS)
                mutationFrac = static_cast<double>(atoi(argv[2]));
        }
        else {
            var = VARIANCE_K;
            if (MUTATIONS)
                mutationFrac = evoParameters::MUTATION_FRAC;
            if (MORTALITY)
                virulence = SISparameters::INITIAL_VIRULENCE;
        }
    }
    
    if (MUTATIONS)
        mutationRate = 1.0 / mutationFrac;

    int succesfulRuns = 0;
    int extinctionCount = 0;

    std::string parameterFileName = "parameters.txt";
    std::ofstream fileLog("outputLog.txt", std::ios_base::app);
    std::ofstream logParameters(parameterFileName, std::ios_base::app);
    
    logParameters << "mutation rate: " << mutationRate << std::endl;
    logParameters << "mutation sd: " << evoParameters::MUTATION_SD << std::endl;
    logParameters << "endtime: " << endtime << std::endl;
    logParameters << "network size: " << NETWORK_SIZE << std::endl;
    logParameters << "recovery rate: " << recovery << std::endl;
    if (DYNAMICS_TYPE == Dynamics::DynamicsType::SIR) {
        logParameters << "transmission: " << transmission << std::endl;
        logParameters << "R0: " << transmission / recovery << std::endl;
    }
    if (MORTALITY)
        logParameters << "initial virulence: " << virulence << std::endl;

	const std::vector<int> seed;
	RandomNumberGenerator rng = create_random_number_generator(seed);

	const std::string filenameNetwork = "networkNew";
    Network network(NETWORK_SIZE);

    if (OPTION_NETWORK_INPUT == false) {
        NetworkGenerator::networkWithoutInput(network, NETWORK_TYPE, var, rng, parameterFileName);
        if (NETWORK_TYPE != NetworkConstructor::NetworkType::FullyConnected){
            Print::networkOutput(filenameNetwork, network);
            Print::printNetwork(filenameNetwork, network);
            Print::printEdgeList(filenameNetwork, network);
        }
    }   

	if (OPTION_NETWORK_INPUT == true) {
		logParameters << "option network: true" << std::endl;
        NetworkGenerator::inputNetwork(network, filenameNetwork);
		DEBUG_ASSERT(network.isValid());
		std::string filenameNetwork2 = filenameNetwork + "2";
		Print::networkOutput(filenameNetwork2, network);
		Print::printNetwork(filenameNetwork2, network);
		Print::printEdgeList(filenameNetwork2, network);
	}

	for (std::size_t t = 0; t < network.size(); ++t) {
		Individual::setSusceptibleNumber(dynamic_cast<Individual &>(network[t]));
	}
   
	for (std::size_t replicate = 1u; replicate < NUMBER_OF_REPLICATES+1; ++replicate) {
        const std::string filenameVirulence = "virulence" + std::to_string(replicate) + ".txt";
		const std::string filenameFinalSize = "finalsize.txt";
        const std::string filenameMaxIncidence = "maxIncidence.txt";
        const std::string filenameSnapshot = "degreeVir" + std::to_string(replicate) + ".csv";
     	std::ofstream finalSize(filenameFinalSize, std::ios_base::app);
        std::ofstream maxInc(filenameMaxIncidence, std::ios_base::app);
		double t = 0.;
		double coefficient;
		assignCoefficient(coefficient, NETWORK_TYPE, DYNAMICS_TYPE);
        int i = getRandom(network.size(), rng);
        Individual &patientZero = dynamic_cast<Individual&>(network[i]);
        if (DYNAMICS_TYPE == Dynamics::DynamicsType::SIR) {
            Pathogen initialPathogen(transmission*coefficient, recovery);   
            patientZero.getInfected(initialPathogen);
        }
        else {
            Pathogen initialPathogen(virulence, coefficient, recovery);
            patientZero.getInfected(initialPathogen);
        }
        Infecteds infecteds(&patientZero,NETWORK_SIZE);
        Individual::updateSusceptibleNeigbours(patientZero, -1);
		bool endEpidemics = false;
		//std::ofstream fileExtinctions("extinctions.txt", std::ios_base::app);
        bool outPutTaken =0;
        if (MUTATIONS) {
            Print::virulenceOutput(filenameVirulence, 0, infecteds);
            Print::virulenceSnapShot(filenameSnapshot, 0, infecteds);
            outPutTaken = 1;
        }
        int maxIncidence = 0;
		do {
                double tau = -log(getRandomUniform(rng)) / Gillespie::rateSum(infecteds);

				Gillespie::selectEvent(infecteds, rng, mutationRate, evoParameters::MUTATION_SD, DYNAMICS_TYPE, MUTATIONS, MORTALITY);

				t += tau;
				
                if (MUTATIONS) {
                    if (static_cast<int>(round(t)) % 5 == 0) {
                        if (outPutTaken == 0) {
                            Print::virulenceOutput(filenameVirulence, static_cast<int>(round(t)), infecteds);
                            Print::virulenceSnapShot(filenameSnapshot, static_cast<int>(round(t)), infecteds);
                            outPutTaken = 1;
                        }
                    }

                    if (static_cast<int>(round(t)) % 5 == 1) {
                        outPutTaken = 0;
                    }
                }
				
                int infectedSize = infecteds.getSizeInfected();
                if (infectedSize > maxIncidence) {
                    maxIncidence = infectedSize;
                }
				if (infecteds.getSizeInfected() == 0) {
					//std::cout << "end of epidemics" << std::endl;
                    if (DYNAMICS_TYPE == Dynamics::DynamicsType::SIR) {
                        int finalRecovered = 0;
                        for (std::size_t i = 0u; i < network.size(); ++i) {
                            if (dynamic_cast<Individual&>(network[i]).getStatus() == Individual::Status::recovered)
                                ++finalRecovered;
                        }
                        finalSize << replicate << "\t" << finalRecovered << std::endl;
                    }
                    maxInc << replicate << "\t" << maxIncidence << std::endl;
					endEpidemics = true;
					break;
				}

		} while (t < endtime && endEpidemics==false); 
		for (std::size_t t = 0; t < network.size(); ++t) {
			dynamic_cast<Individual &>(network[t]).getSusceptible();
			Individual::setSusceptibleNumber(dynamic_cast<Individual &>(network[t]));
		}
        if (t >= endtime) {
            ++succesfulRuns;
        }
        else
            ++extinctionCount;
        if (DYNAMICS_TYPE == Dynamics::DynamicsType::SIS && succesfulRuns >= 10) {
            break;
        }
	}
    fileLog << "succesful runs: " << succesfulRuns << std::endl;
    fileLog << "extinction count: " << extinctionCount << std::endl;
}


