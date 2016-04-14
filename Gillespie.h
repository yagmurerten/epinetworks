#ifndef GILLESPIE_H_INCLUDED
#define GILLESPIE_H_INCLUDED

#include "random.h"

#include <vector>

namespace epinetworks {

	class Network;
	class Individual; 

	class Gillespie {
	public:
		typedef std::vector<Individual*> Infecteds;

		static double rateSum(const Infecteds &infecteds);
		static void selectEvent(Network &population, Infecteds &infecteds, double rateSum, RandomNumberGenerator &rng, double mutationRate);
		static void selectEventSIR(Network &population, Infecteds &infecteds, double rateSum, RandomNumberGenerator &rng);

	private:
		static double getEventRate(Individual &individual);
	};
}

#endif // GILLESPIE_H_INCLUDED