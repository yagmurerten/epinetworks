#include "Debug.h"
#include "Gillespie.h"
#include "individual.h"
#include "network.h"
#include "parameters.h"
#include "pathogen.h"
#include "random.h"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>

namespace epinetworks {
	// Transmission event with mutation
	static void mutationTransmission(Individual &focal, Gillespie::Infecteds &infecteds, RandomNumberGenerator &rng) {
		if (focal.sizeSusceptibleNeighbours() > 0) {
			int newInfected;
			do {
				newInfected = getRandom(focal.getNumberOfContacts(), rng);
			} while (dynamic_cast<Individual &>(focal.getNeighbour(newInfected)).getStatus() == Individual::Status::susceptible);
			Individual &newInfectedNeighbour = dynamic_cast<Individual &>(focal.getNeighbour(newInfected));
			Pathogen mutatedPathogen = Pathogen::mutatePathogen(focal.getPathogen(), MUTATION_SD, rng);
			newInfectedNeighbour.getInfected(mutatedPathogen);
			Individual::updateSusceptibleNeigbours(newInfectedNeighbour, Individual::UpdateRule::Down);
			infecteds.push_back(&newInfectedNeighbour);
		}
	}

	// Transmission event without mutation
	static void noMutationTransmission(Individual &focal, Gillespie::Infecteds &infecteds, RandomNumberGenerator &rng) {
		if (focal.sizeSusceptibleNeighbours() > 0) {
			int newInfected;
			do {
				newInfected = getRandom(focal.getNumberOfContacts(), rng);
			} while (dynamic_cast<Individual &>(focal.getNeighbour(newInfected)).getStatus() != Individual::Status::susceptible);
			Individual &newInfectedNeighbour = dynamic_cast<Individual &>(focal.getNeighbour(newInfected));
			Pathogen oldPathogen = focal.getPathogen();
			newInfectedNeighbour.getInfected(oldPathogen);
			if (newInfectedNeighbour.sizeSusceptibleNeighbours()<0)
				std::cout << "stop 2" << std::endl;
			Individual::updateSusceptibleNeigbours(newInfectedNeighbour, Individual::UpdateRule::Down);
			infecteds.push_back(&newInfectedNeighbour);
		}
	}

	/*
	// Disease induced mortality event
	static void virulence(Individual &focal, Gillespie::Infecteds &infecteds, std::size_t index) {
		focal.getSusceptible();
		Individual::updateSusceptibleNeigbours(focal, Individual::UpdateRule::Up);
		std::swap(infecteds[index], infecteds.back());
		infecteds.pop_back();
	}
	*/

	/*
	// Recovery event SIS
	static void recovery(Individual &focal, Gillespie::Infecteds &infecteds, std::size_t index) {
		focal.getSusceptible();
		Individual::updateSusceptibleNeigbours(focal, Individual::UpdateRule::Up);
		std::swap(infecteds[index], infecteds.back());
		infecteds.pop_back();
	}
	*/

	static void recovery(Individual &focal, Gillespie::Infecteds &infecteds, std::size_t index) {
		focal.getRecovered();
		//Individual::updateSusceptibleNeigbours(focal, Individual::UpdateRule::Up);
		std::swap(infecteds[index], infecteds.back());
		infecteds.pop_back();
	}

	// Gets the total event rate of a single individual 
	double Gillespie::getEventRate(Individual &individual) {
		double eventRate = individual.getPathogen().getTransmission()*individual.sizeSusceptibleNeighbours()
			+ individual.getPathogen().getVirulence() + RECOVERY;
		return eventRate;
	}

	// Gets the total event rate in the population
	double Gillespie::rateSum(const Infecteds &infecteds) {
		double rateSum = 0.;
		for (size_t i = 0; i < infecteds.size(); ++i) {
			rateSum += getEventRate(*infecteds[i]);
		}
		return rateSum;
	}

	// Selects the next individual to have an event
	// and the next event to happen.
	void Gillespie::selectEventSIR(Network &population, Infecteds &infecteds, double rateSum, RandomNumberGenerator &rng) {
		double rand = getRandomUniform(rng);
		double threshold = rateSum*rand;
		double sp = 0.;
		int index = 0;
		double eventRate = 0.;
		while (sp <= threshold) {
			eventRate = getEventRate(*infecteds[index]);
			sp += eventRate;
			if (sp <= threshold) {
				++index;
				if (index > infecteds.size() - 1){
					std::cout << "no more susceptible neighbours" << std::endl;
					std::cout << "infected size" << infecteds.size() << std::endl;
					exit(1);
				}
			}
		}
		sp = 0.;
		rand = getRandomUniform(rng);
		threshold = rand*eventRate;
		Individual &focal = *infecteds[index];
		sp += focal.getPathogen().getTransmission()*focal.sizeSusceptibleNeighbours();
		if (threshold < sp){
			noMutationTransmission(focal, infecteds, rng);
		}
			else {
				sp += RECOVERY;
				if (threshold <= sp) {
					recovery(focal, infecteds, index);
					}
				}
			}
		}





