#ifndef GILLESPIE_H_INCLUDED
#define GILLESPIE_H_INCLUDED

#include "random.h"
#include "Dynamics.h"

#include <memory>
#include <vector>

namespace epinetworks {

	class Network;
	class Individual; 

	class Gillespie {
	public:
        static double rateSum(Infecteds &infecteds);
        static void selectEvent(Infecteds &infecteds, RandomNumberGenerator &rng, double mutationRate, double mutationSD, std::unique_ptr<Dynamics> &dynamics, bool evolution, bool mortality, std::vector<std::vector<int>> &states);
        static std::unique_ptr<Dynamics> createDynamics(Dynamics::DynamicsType type);
	private:
		//static double getEventRate(Individual &individual);
	};
}

#endif // GILLESPIE_H_INCLUDED