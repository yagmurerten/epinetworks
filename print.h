#ifndef PRINT_H_INCLUDED
#define PRINT_H_INCLUDED

#include "percolationCentrality.h"
#include "infecteds.h"

#include <string>
#include <fstream>
#include <vector>

// Functions to print various output

namespace epinetworks {

	class Individual;
	class Network;

	class Print {
	public:

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

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Outputs virulence in every time step it is called
		static void virulenceOutputWithTau(std::string filenameVirulence, double time, double tau,
			Infecteds &infecteds){
			std::ofstream myOutput(filenameVirulence, std::ios_base::app);
            if (time == 0.)
                myOutput << "time,avg_virulence, sd_virulence, avg_neighbours, number_infected, tau " << std::endl;
			myOutput << time;
            double averageVirulence = infecteds.getAverageVirulence();
            double sd = infecteds.getSDVirulence(averageVirulence);
            double averageNeighbours = infecteds.getAverageNumberNeighbour();
            myOutput << "\t" << averageVirulence << "\t" << sd << "\t" << averageNeighbours << "\t";
            myOutput<< infecteds.getSizeInfected() << "\t" << tau << std::endl;
		}

        static void virulenceOutput(std::string filenameVirulence, double time,
            Infecteds &infecteds){
            std::ofstream myOutput(filenameVirulence, std::ios_base::app);
            if (time == 0.)
                myOutput << "time,avg_virulence, sd_virulence, avg_neighbours, number_infected, tau " << std::endl;
            myOutput << time;
            double averageVirulence = infecteds.getAverageVirulence();
            double sd = infecteds.getSDVirulence(averageVirulence);
            double averageNeighbours = infecteds.getAverageNumberNeighbour();
            myOutput << "\t" << averageVirulence << "\t" << sd << "\t" << averageNeighbours << "\t";
            myOutput << infecteds.getSizeInfected() << std::endl;
        }

		static void statusOutput(std::string filenameVirulenceStatus, Network &network){
			std::ofstream myOutput(filenameVirulenceStatus, std::ios_base::app);
			for (std::size_t i = 0; i < network.size(); ++i) {
				Individual &ind = (network[i]);
				int status = static_cast<int>(ind.getStatus());
				myOutput << status << "\t";
			}
			myOutput << std::endl;
		}

		static void virLevel(std::string filenameVirLevel, Network &network){
			std::ofstream myOutput(filenameVirLevel, std::ios_base::app);
			for (std::size_t i = 0; i < network.size(); ++i) {
				Individual &ind = (network[i]);
				double level = ind.getPathogen().getVirulence();
				myOutput << level << "\t";
			}
			myOutput << std::endl;
		}


		// Outputs # of contacts for each individual
		static void networkOutput(std::string filenameNetwork, Network &network) {
			std::ofstream myOutputNetwork(filenameNetwork + ".csv");
			for (std::size_t i = 0u; i < network.size(); ++i) {
				myOutputNetwork << network[i].getNumberOfContacts() << std::endl;
			}
		}


		// Outputs all the levels of virulence at a given time step
        static void virulenceSnapShot(std::string filenameVirSnap, int time, Infecteds &infecteds){
            std::ofstream myOutputVir(filenameVirSnap, std::ios_base::app);
            if (time == 0)
                myOutputVir << "time, degree_4, avg_vir, degree_20, avg_vir, degree_50, avg_vir, degree_100, avg_vir, degree_sup, avg_vir," << std::endl;
            myOutputVir << time << ",";
            std::vector<int> values = { 0, 4, 20, 50, 100, 10000 };
            for (size_t index = 0; index < values.size()-1; ++index) {
                double averageVirulence = infecteds.getAverageVirulenceK(values[index], values[index+1]);
                int count = infecteds.getIndividualsWithKNeighbours(values[index], values[index+1]);
                myOutputVir << count << "," << averageVirulence << ",";
            }
            myOutputVir << std::endl;
        }

		// Outputs all the infected nodes at a given time step
		static void infectedNodesOutput(int time, std::size_t rep, std::vector<Individual*> &infecteds) {
			std::ofstream outputInfecteds(std::to_string(rep) + "infecteds" + std::to_string(time) + ".csv", std::ios_base::app);
			for (std::size_t i = 0; i < infecteds.size(); ++i) {
				Individual &ind = *infecteds[i];
				outputInfecteds << ind.getCoordinate() << "," << ind.sizeNeighbour() << std::endl;
			}
		}

        static void printStates(double time, const std::vector<std::vector<int>> &states, int replicate) {
            for (std::size_t i = 0u; i < states.size(); ++i) {
                int coordLoc = states[i].size() - 1;
                int index = states[i][coordLoc];
                std::ofstream outputState("STATES/state" + std::to_string(replicate) + "_" + std::to_string(index) + ".csv", std::ios_base::app);
                outputState << time;
                for (std::size_t j = 0u; j < states[i].size()-1; ++j) {                
                    outputState << "," << states[i][j];
                }
                outputState << std::endl;
            }
        }

        static void printStatesAll(double t, int time, int replicate, Network &network) {
            for (std::size_t i = 0u; i < network.size(); ++i) {
                Individual &ind = (network[i]);
                std::ofstream outputState("STATES/state" + std::to_string(replicate) + "_" + std::to_string(i) + ".csv", std::ios_base::app);
                outputState << t << "," << time;
                for (std::size_t j = 0u; j < ind.getStates().size() ; ++j) {
                    outputState << "," << ind.getStates()[j];
                }
                outputState << std::endl;
            }
        }

        template <class T>
        static void printParameter(T parameter, const char type[], std::string &fileName) {
            std::ofstream logParameters(fileName, std::ios_base::app);
            logParameters << type << ": " << parameter << std::endl;
        }
	};
}
#endif // PRINT_H