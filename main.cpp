//////////////////////////////////////////////////////////////////////////////////
/////////// EpiNetworks - E. Yagmur Erten ////////////////////////////////////////
/////////// 08.04.2016 ///////////////////////////////////////////////////////////
/////////// Code for simulating epidemics (SIR) on heterogeneous random networks /
//////////////////////////////////////////////////////////////////////////////////

#include "connectionPool.h"
#include "Debug.h"
#include "Gillespie.h"
#include "individual.h"
#include "network.h"
#include "networkGenerator.h"
#include "parameters.h"
#include "pathogen.h"
#include "print.h"

#include "random.h"

#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>

using namespace epinetworks;
// some parameters for main.cpp
std::size_t NUMBER_OF_REPLICATES = 100;	// # of replicates on each network
const int ENDTIME = 1500;				// # of time steps
double INITIAL_VIRULENCE = 0.25;			// initial value of virulence
int NETWORK_SIZE = 10000;				// # of individuals in the network
int MUTATION_FRAC = 0;
double VARIANCE_K = 4;
double INITIAL_TRANSMISSION = 1.1;			// initial value of virulenc
int NETWORK_TYPE = 2;

// Initializes a network of a given size and type
// using var as the variance in the degree distribution
// and 3.0 as mean of the degree distribution
void initializeNetwork(Network &network, NetworkGenerator::NetworkType type, double var, RandomNumberGenerator &rng) {
	if (type == NetworkGenerator::NetworkType::Gamma) {
		double beta = 3.0 / var;
		const double shapeParameter = 3.0*beta;			//alpha
		const double scaleParameter = 1.0 / beta;		//beta
		for (std::size_t i = 0; i < network.size(); ++i) {
				int numberOfContacts = static_cast<int>(round(
					getRandomGamma(shapeParameter, scaleParameter, rng))) + 1;
				network.setNode(std::unique_ptr<Individual>(new Individual(numberOfContacts)));
				network[i].setCoordinates(i);
		}
	}
	if (type == NetworkGenerator::NetworkType::FullyConnected) {
		int numberOfContacts = network.size() - 1;
		for (std::size_t i = 0; i < network.size(); ++i) {
			network.setNode(std::unique_ptr<Individual>(new Individual(numberOfContacts)));
			network[i].setCoordinates(i);
		}
	}
	if (type == NetworkGenerator::NetworkType::Homogeneous) {
		int NUMBER_OF_CONTACTS = static_cast<int>(var);
		for (std::size_t i = 0; i < network.size(); ++i) {
			network.setNode(std::unique_ptr<Individual>(new Individual(NUMBER_OF_CONTACTS)));
			network[i].setCoordinates(i);
		}

	}
}

//Initialize and generate network with input
//input=networkNew
void inputNetwork(Network &network, const std::string &input) {
	std::ifstream inputNumbers(input + ".csv");
	std::ifstream inputContacts(input + "Matrix.txt");
	for (std::size_t i = 0; i < network.size(); ++i) {
		int numberOfContacts;
		std::string line;
		std::getline(inputNumbers, line);
		std::stringstream ss(line);
		ss >> numberOfContacts;
		network.setNode(std::unique_ptr<Individual>(new Individual(numberOfContacts)));
		network[i].setCoordinates(i);
	}
	for (std::size_t i = 0; i < network.size(); ++i) {
		std::string line2;
		std::getline(inputContacts, line2);
		std::stringstream ss2(line2);
		if (i == 54) 
			int k = 120;
		int stubs = network[i].getStubs();
		for (std::size_t j = 0; j < stubs; ++j) {
			int stubs2 = network[i].getStubs();
			int coordContact;
			ss2 >> coordContact;
			NetworkNode &first = network[i];
			NetworkNode &second = network[coordContact];
			if (first.isNeighbour(second) == false){
				first.setNeighbour(second);
				second.setNeighbour(first);
				second.reduceStubs();
				first.reduceStubs();
			}
			else
				--j;
		}
	}
}

int main(int, char *argv[]){
	double fracTransmission = static_cast<double>(atoi(argv[1]));
	double transmission = fracTransmission / 100.0;
	std::ofstream logParameters("parameters.txt", std::ios_base::app);

	const std::vector<int> seed;

	RandomNumberGenerator rng = create_random_number_generator(seed);
	const std::string filenameNetwork = "networkNew";

	Network network(NETWORK_SIZE);
	double coefficient;
	double mutationFrac;
	if (OPTION_NETWORK_INPUT == false) {
		logParameters << "option network: false" << std::endl;
		double var = static_cast<double>(atoi(argv[2]));
		int networkType = NETWORK_TYPE; // atoi(argv[3]);
		switch (networkType) {
		case 0: {
					NetworkGenerator::NetworkType TYPE = NetworkGenerator::NetworkType::FullyConnected;
					initializeNetwork(network, TYPE, var, rng);
					DEBUG_ASSERT(network.isValid());
					NetworkGenerator::generate(network, TYPE, rng);
					DEBUG_ASSERT(network.isValid());
					coefficient = E;
					logParameters << "network type: fully connected" << std::endl; }
			break;
		case 1: {
					NetworkGenerator::NetworkType TYPE = NetworkGenerator::NetworkType::Gamma;
					initializeNetwork(network, TYPE, var, rng);
					DEBUG_ASSERT(network.isValid());
					NetworkGenerator::generate(network, TYPE, rng);
					DEBUG_ASSERT(network.isValid());
					Print::networkOutput(filenameNetwork, network);
					Print::printNetwork(filenameNetwork, network);
					Print::printEdgeList(filenameNetwork, network);
					coefficient = F; 
					logParameters << "network type: gamma distributed" << std::endl;
					logParameters << "variance: " << var << std::endl; }
			break;
		case 2: {
					NetworkGenerator::NetworkType TYPE = NetworkGenerator::NetworkType::Homogeneous;
					initializeNetwork(network, TYPE, var, rng);
					DEBUG_ASSERT(network.isValid());
					NetworkGenerator::generate(network, TYPE, rng);
					DEBUG_ASSERT(network.isValid());
					Print::networkOutput(filenameNetwork, network);
					Print::printNetwork(filenameNetwork, network);
					Print::printEdgeList(filenameNetwork, network);
					coefficient = F; 
					logParameters << "network type: homogeneous" << std::endl;
					logParameters << "number of contacts: " << var << std::endl; }
			break;
		}
	}

	//assuming we only want it for heterogeneous networks
	if (OPTION_NETWORK_INPUT == true) {
		logParameters << "option network: true" << std::endl;
		coefficient = D;
		inputNetwork(network, filenameNetwork);
		DEBUG_ASSERT(network.isValid());
		std::string filenameNetwork2 = filenameNetwork + "2";
		Print::networkOutput(filenameNetwork2, network);
		Print::printNetwork(filenameNetwork2, network);
		Print::printEdgeList(filenameNetwork2, network);
	}

	mutationFrac = static_cast<double>(MUTATION_FRAC); //static_cast<double>(atoi(argv[1]));

	double mutationRate = 1.0 / mutationFrac;
	logParameters << "mutation rate: " << mutationFrac << std::endl;
	logParameters << "mutation sd: " << MUTATION_SD << std::endl;
	logParameters << "endtime: " << ENDTIME << std::endl;
	logParameters << "transmission: " << transmission << std::endl;
	logParameters << "network size: " << NETWORK_SIZE << std::endl;
	logParameters << "recovery rate: " << RECOVERY << std::endl;
	logParameters << "R0: " << transmission / RECOVERY << std::endl;


	for (std::size_t t = 0; t < network.size(); ++t) {
		Individual::setSusceptibleNumber(dynamic_cast<Individual &>(network[t]));
	}

	for (std::size_t i = 0u; i < network.size(); ++i) {
		// Throws an exception if the node is not an Individual.
		Individual &individual = dynamic_cast<Individual &>(network[i]);
		// Return null if the node is not an Individual.
		Individual *individualPtr = dynamic_cast<Individual*>(&network[i]);
		DEBUG_ASSERT(individualPtr != nullptr);
	}
    
	for (std::size_t replicate = 1u; replicate < NUMBER_OF_REPLICATES+1; ++replicate) {
		const std::string filenameVirulence2 = filenameNetwork + "virulence" + std::to_string(replicate) + ".txt";
		const std::string filenameFinalSize = "finalsize.txt";
		std::ofstream finalSize(filenameFinalSize, std::ios_base::app);
		double t = 0.;
		Gillespie::Infecteds infecteds;
		Pathogen initialPathogen(transmission*coefficient);
		int i = getRandom(network.size(), rng);
		Individual &patientZero = dynamic_cast<Individual&>(network[i]);
		patientZero.getInfected(initialPathogen);
		infecteds.push_back(&patientZero);
		Individual::updateSusceptibleNeigbours(patientZero, Individual::UpdateRule::Down);
		Print::virulenceOutput(filenameVirulence2, t, infecteds, network);
		bool endEpidemics = false;
		//std::ofstream fileExtinctions("extinctions.txt", std::ios_base::app);
        bool outPutTaken =0;
		do {

				double rTotal = Gillespie::rateSum(infecteds);

				double tau = -log(getRandomUniform(rng)) / rTotal;

				Gillespie::selectEventSIR(network, infecteds, rTotal, rng);

				t += tau;
				
				Print::virulenceOutput(filenameVirulence2, t, infecteds, network);
							
				if (infecteds.size() == 0) {
					std::cout << "end of epidemics" << std::endl;
					int finalRecovered = 0;
					for (std::size_t i = 0u; i < network.size(); ++i) {
						if (dynamic_cast<Individual&>(network[i]).getStatus() == Individual::Status::recovered)
							++finalRecovered;
						}
					finalSize << replicate << "\t" << finalRecovered << std::endl;
					endEpidemics = true;
					break;
					
				}
		} while (t < ENDTIME && endEpidemics==false); 
		for (std::size_t t = 0; t < network.size(); ++t) {
			dynamic_cast<Individual &>(network[t]).getSusceptible();
			Individual::setSusceptibleNumber(dynamic_cast<Individual &>(network[t]));
		}
	}

}

