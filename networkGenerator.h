#ifndef NETWORKGENERATOR_H_INCLUDED
#define NETWORKGENERATOR_H_INCLUDED

#include "NetworkConstructor.h"
#include "random.h"

namespace epinetworks {

	class Network; 

	class NetworkGenerator {
	public:

        static void networkWithoutInput(Network &network, NetworkConstructor::NetworkType networkType, double var, RandomNumberGenerator rng, std::string logFile);
        static void inputNetwork(Network &network, const std::string &input);
	};

}



#endif // NETWORKGENERATOR