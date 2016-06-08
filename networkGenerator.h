#ifndef NETWORKGENERATOR_H_INCLUDED
#define NETWORKGENERATOR_H_INCLUDED

#include "random.h"

namespace epinetworks {

	class Network; 
	class NetworkNode;

	class NetworkGenerator {
	public:

		enum class NetworkType {
			FullyConnected,
			Gamma,
			Homogeneous,
			PowerLaw, 
			SmallWorld, 
			ErdosRenyi
		};

		static void generate(Network &network, NetworkType type, RandomNumberGenerator &rng);
		static const bool areNeighbours(NetworkNode &focal, NetworkNode &potentialNeighbour);
	};

}



#endif // NETWORKGENERATOR