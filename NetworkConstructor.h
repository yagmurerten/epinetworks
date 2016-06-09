#ifndef NETWORKCONSTRUCTOR_H_INCLUDED
#define NETWORKCONSTRUCTOR_H_INCLUDED

#include "random.h"

namespace epinetworks {

	class Network; 
	class NetworkNode;

	class NetworkConstructor {
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



#endif // NETWORKCONSTRUCTOR