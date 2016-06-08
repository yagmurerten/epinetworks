
#include "connectionPool.h"
#include "network.h"
#include "networkGenerator.h"
#include "random.h"

#include <algorithm>
#include <random>
#include <iostream>
//

namespace epinetworks {

	// Makes a fully connected network (well-mixed model)
	static void makeFullLinks(Network &network) {
		for (std::size_t i = 0u; i< network.size(); ++i) {
			//std::cout << i << std::endl;
			for (std::size_t j = 0u; j < network.size(); ++j) {
						if (i != j) {
							network[i].setNeighbour(network[j]);
							network[i].reduceStubs();
				}
			}
		}
		std::cout << "done\n";
	}

	// Creates connection pool from a given network
	static void makeConnectionPool(ConnectionPool &connectionPool, Network &network) {
		for (std::size_t i = 0u; i < network.size(); ++i) {
				connectionPool.add(network[i]);
		}
	}

	// Returns total # of stubs in the network
	static const int totalStubs(Network &network) {
		int totalStubs = 0;
		for (std::size_t i = 0u; i < network.size(); ++i) {
				totalStubs += network[i].getStubs();
			}
		return totalStubs;
	}

	// Returns total # of stubs in the pool
	static const int totalStubsPool(ConnectionPool &connectionPool) {
		int totalStubs = 0;
		for (size_t i = 0u; i < connectionPool.size(); ++i) {
			totalStubs += connectionPool[i].getStubs();
		}
		return totalStubs;
	}

	// Factorial (X!)
	inline int factorial(int n){
		int factorial = 1;
		for (int i = n; i>1; --i)
			factorial *= i;
		return factorial;
	}

	// Checks if all the individuals left in the pool are already neighbours
	static bool checkIfAllNeighbour(ConnectionPool &connectionPool) {
		int count = 0;
		for (std::size_t i = 0u; i < connectionPool.size(); ++i) {
			for (std::size_t j = 0u; j < connectionPool.size(); ++j) {
				if (i != j) {
					if (connectionPool[i].isNeighbour(connectionPool[j]) == true)
						++count;
				}
			}
		}
		int combination = factorial(connectionPool.size()) / (2 * factorial(connectionPool.size() - 2));
		if (count == combination * 2)
			return true;
		else
			return false;
	}

	// Connects two individuals
	static void connect(NetworkNode &first, NetworkNode &second, Network &network)  {
		first.setNeighbour(second);
		second.setNeighbour(first);
		second.reduceStubs();
		first.reduceStubs();
	}

	// Makes a randomly linked network in which individuals' # of
	// contacts is predetermined
	static void makeRandomLinks(Network &network, RandomNumberGenerator &rng)  {
		int totalNumberOfStubs = totalStubs(network);
		ConnectionPool connectionPool;
		makeConnectionPool(connectionPool, network);
		int totalNumberOfStubsPool = totalStubsPool(connectionPool);
		while (connectionPool.size() > 1) {
			//int i = getRandom(connectionPool.size());
			int i = 0;
			for (std::size_t index = 0; index < connectionPool.size(); ++index) {
				if (connectionPool[index].getStubs() >= connectionPool[i].getStubs())
					i = index;
			}
			if (connectionPool.size() < 5 && connectionPool.size() > 1) {
				bool allNeighbour = checkIfAllNeighbour(connectionPool);
				if (allNeighbour == true) {
					std::cout << "all neighbour" << std::endl;
					while (connectionPool.size() > 0) {
						NetworkNode &node = connectionPool[connectionPool.size() - 1];
						node.updateConnectionZero();
						connectionPool.remove(connectionPool.size() - 1);
					}
					break;
				}
			}
			int j;
			NetworkNode &first = connectionPool[i];
			do {
				j = getRandom(connectionPool.size(), rng);
			} while ((i == j));
			NetworkNode &second = connectionPool[j];
			if (first.isNeighbour(second) == false) {
				if ((second.getStubs() != 0) && (first.getStubs() != 0))
					connect(first, second, network);
			}
			// some debugging code
			/*
			if (second.getStubs() < 0) {
				std::cout << "error4\n";
				system("pause");
			}
			
			
			totalNumberOfStubsPool = totalStubsPool(connectionPool);
			totalNumberOfStubs = totalStubs(network);
			if (totalNumberOfStubs != totalNumberOfStubsPool)
				std::cout << "problem" << std::endl;
			*/
			if (first.getStubs() == 0) {
				connectionPool.remove(i);
				if (network[i].getStubs() != 0)
					std::cout << "problem" << std::endl;
			}
			if (second.getStubs() == 0) {
				connectionPool.remove(j);
				if (network[j].getStubs() != 0)
					std::cout << "problem" << std::endl;
			}
			}
			if (connectionPool.size() == 1) {
				connectionPool[0].updateConnectionZero();
				connectionPool.remove(0);
			}
			int lastNumberOfStubs = totalStubs(network);

			if (lastNumberOfStubs != 0) {
				std::cout << lastNumberOfStubs << std::endl;
				for (std::size_t t = 0; t < network.size(); ++t) {
					if (network[t].getStubs() != 0) {
						network[t].updateConnectionZero();
					}
				}
				// some debugging code
				/*
				lastNumberOfStubs = totalStubs(network);
				if (lastNumberOfStubs != 0) {
					std::cout << "error5\n";
					system("pause");
				}
				*/
			}
	}

	// Other types of networks are not yet implemented
	static void makeSmallWorldLinks(Network &network) {};

	static void makePowerLawLinks(Network &network) {};

	static void makeErdosRenyiLinks(Network &network) {};

	// Calls the generator for a given type
	void NetworkGenerator::generate(Network &network, NetworkType type, RandomNumberGenerator &rng) {
		switch (type) {
		case NetworkType::FullyConnected:
			makeFullLinks(network);
			break;
		case NetworkType::Gamma:
			makeRandomLinks(network, rng);
			break;
		case NetworkType::PowerLaw:
			makePowerLawLinks(network);
			break;
		case NetworkType::SmallWorld:
			makeSmallWorldLinks(network);
			break;
		case NetworkType::ErdosRenyi:
			makeErdosRenyiLinks(network);
			break;
		case NetworkType::Homogeneous:
			makeRandomLinks(network, rng);
			break;
		}
	}	
}




