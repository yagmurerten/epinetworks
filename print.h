#ifndef PRINT_H_INCLUDED
#define PRINT_H_INCLUDED

#include "percolationCentrality.h"

#include <string>
#include <fstream>
#include <vector>

// Functions to print various output

namespace epinetworks {

	class Individual;
	class Network;

	class Print {
	public:
		// Outputs virulence in every time step it is called
		static void virulenceOutput(std::string filenameVirulence, double time, 
			std::vector<Individual*> &infecteds, Network &network){
			std::ofstream myOutput2(filenameVirulence, std::ios_base::app);
			myOutput2 << time;
			double totalVirulence = 0.;
			for (std::size_t i = 0; i < infecteds.size(); ++i) {
				Individual &ind = *infecteds[i];
				totalVirulence += ind.getPathogen().getVirulence();
			}
			double averageVirulence = totalVirulence / infecteds.size();
			double sd = 0.;
			for (std::size_t i = 0; i < infecteds.size(); ++i) {
				Individual &ind = *infecteds[i];
				sd += pow((ind.getPathogen().getVirulence() - averageVirulence), 2.0);
			}
			double totalNeighbours = 0;
			for (std::size_t i = 0; i < infecteds.size(); ++i) {
				Individual &ind = *infecteds[i];
				totalNeighbours += ind.sizeNeighbour();
			}
			double averageNeighbours = totalNeighbours / infecteds.size();
			sd = sqrt(sd / infecteds.size());
			double virSuperspread = 0.;
			int indexSuperspread = 0;
			for (std::size_t i = 0; i < infecteds.size(); ++i) {
				Individual &ind = *infecteds[i];
				if (ind.sizeNeighbour() >= 50) {
					virSuperspread += ind.getPathogen().getVirulence();
					++indexSuperspread;
				}
			}
			double averageSuperVirulence = 0;
			if (indexSuperspread != 0)
				averageSuperVirulence = virSuperspread / static_cast<double>(indexSuperspread);
			double virOther = 0.;
			int indexOther = 0;
			for (std::size_t i = 0; i < infecteds.size(); ++i) {
				Individual &ind = *infecteds[i];
				if (ind.sizeNeighbour() < 50) {
					virOther += ind.getPathogen().getVirulence();
					++indexOther;
				}
			}
			double averageOtherVirulence = 0;
			if (indexOther != 0)
				averageOtherVirulence = virOther/ indexOther;
			double ratio = PercolationCentrality::calculatePCBCratio(network);
			myOutput2 << "\t" << averageVirulence << "\t" << sd << "\t" << averageNeighbours << "\t" <<
				averageSuperVirulence << "\t" << averageOtherVirulence << "\t" << infecteds.size() << "\t" << 
				ratio;
			myOutput2 << std::endl;
			//std::cout << averageVirulence << " " << infecteds.size() << std::endl;

		}

		/*static void percolationOutput(std::string filenameVirulenceStatus, Network &network){
			std::ofstream myOutput(filenameVirulenceStatus, std::ios_base::app);
			double ratio = PercolationCentrality::calculatePCBCratio(network);
			myOutput << ratio << std::endl;
		}*/

		static void statusOutput(std::string filenameVirulenceStatus, Network &network){
			std::ofstream myOutput(filenameVirulenceStatus, std::ios_base::app);
			for (std::size_t i = 0; i < network.size(); ++i) {
				Individual &ind = dynamic_cast<Individual &>(network[i]);
				int status = static_cast<int>(ind.getStatus());
				myOutput << status << "\t";
			}
			myOutput << std::endl;
			//std::cout << averageVirulence << " " << infecteds.size() << std::endl;
		}

		static void virLevel(std::string filenameVirLevel, Network &network){
			std::ofstream myOutput(filenameVirLevel, std::ios_base::app);
			for (std::size_t i = 0; i < network.size(); ++i) {
				Individual &ind = dynamic_cast<Individual &>(network[i]);
				double level = ind.getPathogen().getVirulence();
				myOutput << level << "\t";
			}
			myOutput << std::endl;
			//std::cout << averageVirulence << " " << infecteds.size() << std::endl;
		}


		// Outputs # of contacts for each individual
		static void networkOutput(std::string filenameNetwork, Network &network) {
			std::ofstream myOutputNetwork(filenameNetwork + ".csv");
			for (std::size_t i = 0u; i < network.size(); ++i) {
				myOutputNetwork << network[i].getNumberOfContacts() << std::endl;
			}
		}

		// Outputs neighbours of all the individiuals
		static void printNetwork(std::string filenameNetwork, Network &network) {
			std::ofstream myOutputMatrix(filenameNetwork + "Matrix.txt");
			for (std::size_t i = 0u; i < network.size(); ++i) {
				for (std::size_t j = 0u; j < network[i].sizeNeighbour(); ++j) {
					myOutputMatrix << network[i].getNeighbour(j).getCoordinate() << "\t";
				}
				myOutputMatrix << std::endl;
			}
		}

		// Outputs neighbours of all the individuals as .csv
		static void printNetworkCsv(std::string filenameNetwork, Network &network) {
			std::ofstream myOutputMatrix(filenameNetwork + "Matrix.csv");
			for (std::size_t i = 0u; i < network.size(); ++i) {
				for (std::size_t j = 0u; j < network[i].sizeNeighbour(); ++j) {
					myOutputMatrix << network[i].getNeighbour(j).getCoordinate() << ",";
				}
				myOutputMatrix << std::endl;
			}
		}

		// Outputs edge list 
		static void printEdgeList(std::string filenameNetwork, Network &network) {
			std::ofstream myOutputMatrix(filenameNetwork + "EdgeList.csv");
			for (std::size_t i = 0u; i < network.size(); ++i) {
				for (std::size_t j = 0u; j < network[i].sizeNeighbour(); ++j) {
					std::size_t index = network[i].getNeighbour(j).getCoordinate();
					if (index>i) {
						myOutputMatrix << network[i].getCoordinate() << ",";
						myOutputMatrix << index << std::endl;
					}
				}
			}
		}

		// Outputs adjacency matrix
		static void printAdjMatrix(std::string filenameNetwork, Network &network){
			std::ofstream myOutputMatrix(filenameNetwork + "AdjMatrix.csv");
			myOutputMatrix << " ,";

			for (std::size_t i = 0u; i < network.size(); ++i) {
				myOutputMatrix << std::to_string(i) << ",";
			}

			myOutputMatrix << std::endl;
			for (std::size_t i = 0u; i < network.size(); ++i) {
				myOutputMatrix << std::to_string(i) << ",";
				for (std::size_t j = 0u; j < network.size(); ++j) {
					if (network[i].isNeighbour(network[j]) == true)
						myOutputMatrix << "1,";
					else if (network[i].isNeighbour(network[j]) == false)
						myOutputMatrix << "0,";
				}
				myOutputMatrix << std::endl;
			}
		}

		// Outputs all the levels of virulence at a given time step
		static void virulenceSnapShot(int time, std::size_t rep, std::vector<Individual*> &infecteds){
			std::ofstream myOutputVir(std::to_string(rep) + "pathogens" + std::to_string(time) + ".csv", std::ios_base::app);
			double totalVirulence = 0.;
			for (std::size_t i = 0; i < infecteds.size(); ++i) {
				Individual &ind = *infecteds[i];
				myOutputVir << ind.getPathogen().getVirulence() << std::endl;
			}
		}

		// Outputs all the infected nodes at a given time step
		static void infectedNodesOutput(int time, std::size_t rep, std::vector<Individual*> &infecteds) {
			std::ofstream outputInfecteds(std::to_string(rep) + "infecteds" + std::to_string(time) + ".csv", std::ios_base::app);
			for (std::size_t i = 0; i < infecteds.size(); ++i) {
				Individual &ind = *infecteds[i];
				outputInfecteds << ind.getCoordinate() << "," << ind.sizeNeighbour() << std::endl;
			}
		}
	};
}
#endif // PRINT_H