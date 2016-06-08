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

		enum class DynamicsType {
			SIS,
			SIR
		};

		static double rateSum(const Infecteds &infecteds, double recoveryRate);
		static void selectEvent(Network &population, Infecteds &infecteds, double rateSum, RandomNumberGenerator &rng, double mutationRate, double mutationSD, double recoveryRate);
		static void selectEventSIR(Network &population, Infecteds &infecteds, double rateSum, RandomNumberGenerator &rng, double mutationRate, double mutationSD, double recoveryRate);

	private:
		static double getEventRate(Individual &individual, double recoveryRate);
	};
}

#endif // GILLESPIE_H_INCLUDED