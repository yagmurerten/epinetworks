
#include "connectionPool.h"
#include "individual.h"
#include "network.h"
#include "networkConstructor.h"
#include "networkGenerator.h"
#include "random.h"

#include <algorithm>
#include <random>
#include <iostream>
#include <memory>
#include <fstream>
#include <sstream>

namespace epinetworks {

	
    // Initializes a network of a given size and type
    // using var as the variance in the degree distribution
    // and 3.0 as mean of the degree distribution
    void initializeNetwork(Network &network, NetworkConstructor::NetworkType type, double var, RandomNumberGenerator &rng) {
        if (type == NetworkConstructor::NetworkType::Gamma && var != 0) {
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
        if (type == NetworkConstructor::NetworkType::FullyConnected) {
            int numberOfContacts = network.size() - 1;
            for (std::size_t i = 0; i < network.size(); ++i) {
                network.setNode(std::unique_ptr<Individual>(new Individual(numberOfContacts)));
                network[i].setCoordinates(i);
            }
        }
        if (type == NetworkConstructor::NetworkType::Homogeneous) {
            int NUMBER_OF_CONTACTS = static_cast<int>(var);
            for (std::size_t i = 0; i < network.size(); ++i) {
                network.setNode(std::unique_ptr<Individual>(new Individual(NUMBER_OF_CONTACTS)));
                network[i].setCoordinates(i);
            }
           
        }
        if (type == NetworkConstructor::NetworkType::Gamma && var == 0) {
            int NUMBER_OF_CONTACTS = 4;
            for (std::size_t i = 0; i < network.size(); ++i) {
                network.setNode(std::unique_ptr<Individual>(new Individual(NUMBER_OF_CONTACTS)));
                network[i].setCoordinates(i);
            }
        }
    }

    void NetworkGenerator::networkWithoutInput(Network &network, NetworkConstructor::NetworkType networkType, double var, RandomNumberGenerator rng, std::string logFile) {
        std::ofstream logParameters(logFile, std::ios_base::app);
            logParameters << "option network: false" << std::endl;
            if (networkType == NetworkConstructor::NetworkType::FullyConnected) {
                initializeNetwork(network, networkType, var, rng);
                DEBUG_ASSERT(network.isValid());
                NetworkConstructor::generate(network, networkType, rng);
                DEBUG_ASSERT(network.isValid());
                logParameters << "network type: fully connected" << std::endl;
            }
            else if (networkType == NetworkConstructor::NetworkType::Gamma) {
                double varGamma = pow(var, 2);
                initializeNetwork(network, networkType, varGamma, rng);
                DEBUG_ASSERT(network.isValid());
                NetworkConstructor::generate(network, networkType, rng);
                DEBUG_ASSERT(network.isValid());
                logParameters << "network type: gamma distributed" << std::endl;
                logParameters << "variance: " << varGamma << std::endl;
            }
            else if (networkType == NetworkConstructor::NetworkType::Homogeneous) {
                initializeNetwork(network, networkType, var, rng);
                DEBUG_ASSERT(network.isValid());
                NetworkConstructor::generate(network, networkType, rng);
                DEBUG_ASSERT(network.isValid());
                logParameters << "network type: homogeneous" << std::endl;
                logParameters << "number of contacts: " << var << std::endl;
            }
            else {
                exit(5);
            }
    }

    //Initialize and generate network with input
    //input=networkNew
    void NetworkGenerator::inputNetwork(Network &network, const std::string &input) {
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
            int stubs = network[i].getStubs();
            for (std::size_t j = 0; j < stubs; ++j) {
                int coordContact;
                ss2 >> coordContact;
                Individual &first = network[i];
                Individual &second = network[coordContact];
                if (first.isNeighbour(second) == false){
                    first.setNeighbourCoord(second);
                    second.setNeighbourCoord(first);
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

}




