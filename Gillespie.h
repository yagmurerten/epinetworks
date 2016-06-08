#ifndef GILLESPIE_H_INCLUDED
#define GILLESPIE_H_INCLUDED

#include "random.h"
#include "Dynamics.h"

#include <vector>

namespace epinetworks {

	class Network;
	class Individual; 

	class Gillespie {
	public:
        static double rateSum(Infecteds &infecteds);
        static void selectEvent(Infecteds &infecteds, double rateSum, RandomNumberGenerator &rng, double mutationRate, double mutationSD, Dynamics::DynamicsType type, bool evolution);
	private:
		static double getEventRate(Individual &individual);
	};
}

#endif // GILLESPIE_H_INCLUDED