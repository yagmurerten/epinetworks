#include "Debug.h"
#include "Gillespie.h"
#include "Dynamics.h"
#include "DynamicsSIR.h"
#include "DynamicsSIS.h"

#include "infecteds.h"
#include "individual.h"
#include "network.h"
#include "NetworkConstructor.h"

#include "pathogen.h"
#include "random.h"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>

namespace epinetworks {

	// Gets the total event rate of a single individual 
	//double Gillespie::getEventRate(Individual &individual) {
	//	double eventRate = individual.getPathogen().getTransmission()*individual.sizeSusceptibleNeighbours()
 //           + individual.getPathogen().getVirulence() + individual.getPathogen().getRecoveryRate();
	//	return eventRate;
	//}

	// Gets the total event rate in the population
    double Gillespie::rateSum(Infecteds &infecteds) {
		double rateSum = 0.;
		for (size_t i = 0; i < infecteds.getSizeInfected(); ++i) {
			rateSum += (infecteds.returnIndividual(i)).getEventRate();
		}
		return rateSum;
	}

	// Selects the next individual to have an event
	// and the next event to happen.
   
    std::unique_ptr<Dynamics> Gillespie::createDynamics(Dynamics::DynamicsType type) {
        if (type == Dynamics::DynamicsType::SIS)
            return std::unique_ptr<Dynamics>(new DynamicsSIS);
        return std::unique_ptr<Dynamics>(new DynamicsSIR);
    }

    void Gillespie::selectEvent(Infecteds &infecteds, RandomNumberGenerator &rng, double mutationRate, double mutationSD, 
        std::unique_ptr<Dynamics> &dynamics, bool evolution, bool mortality, std::vector<std::vector<int>> &states) {
        double rand = getRandomUniform(rng);
        double rateSum = Gillespie::rateSum(infecteds);
        double threshold = rateSum*rand;
        double sp = 0.;
        int index = 0;
        double eventRate = 0.;
        while (sp <= threshold) {
            eventRate = (infecteds.returnIndividual(index)).getEventRate();
            sp += eventRate;
            if (sp <= threshold) {
                ++index;
                if (index > infecteds.getSizeInfected() - 1){
                    //double newRate = Gillespie::rateSum(infecteds);
                    std::cout << "no more susceptible neighbours" << std::endl;
                    std::cout << "infected size" << infecteds.getSizeInfected() << std::endl;
                    exit(7);
                }
            }
        }
        sp = 0.;
        rand = getRandomUniform(rng);
        threshold = rand*eventRate;
        Individual &focal = infecteds.returnIndividual(index);

        if (evolution) {
            if (mortality) {
            sp += mutationRate*focal.getPathogen().getTransmission()*focal.sizeSusceptibleNeighbours();
            if (threshold < sp){
                dynamics->transmission(focal, infecteds, rng, mutationSD, 1, states);
            }
            else {
                sp += (1 - mutationRate)*focal.getPathogen().getTransmission()*focal.sizeSusceptibleNeighbours();
                if (threshold < sp) {
                    dynamics->transmission(focal, infecteds, rng, mutationSD, 0, states);
                }
                else {
                    sp += focal.getPathogen().getRecoveryRate();
                    if (threshold < sp) {
                        dynamics->recovery(focal, infecteds, index, states);
                     }
                    else {
                        sp += focal.getPathogen().getVirulence();
                        if (threshold <= sp) {
                            dynamics->virulence(focal, infecteds, index);
                        }
                    }
                   }
               }
           }
            else {
                sp += mutationRate*focal.getPathogen().getTransmission()*focal.sizeSusceptibleNeighbours();
                if (threshold < sp){
                    dynamics->transmission(focal, infecteds, rng, mutationSD, 1, states);
                }
                else {
                    sp += (1 - mutationRate)*focal.getPathogen().getTransmission()*focal.sizeSusceptibleNeighbours();
                    if (threshold < sp) {
                        dynamics->transmission(focal, infecteds, rng, mutationSD, 0, states);
                    }
                    else {
                        sp += focal.getPathogen().getRecoveryRate();
                        if (threshold <= sp) {
                            dynamics->recovery(focal, infecteds, index, states);
                        }
                    }
                }
            }
        }
        // state updates are implemented dynamics without mutations or mortality!
        else {
            sp += focal.getPathogen().getTransmission()*focal.sizeSusceptibleNeighbours();
            if (threshold < sp){
                dynamics->transmission(focal, infecteds, rng, mutationSD, 0, states);
                // when transmission, states of the new infected individual and its neighbours' change
                //updateStates(infecteds.returnIndividual(infecteds.getSizeInfected() - 1), states);
            }
            else {
                sp += focal.getPathogen().getRecoveryRate();
                if (threshold <= sp) {
                    dynamics->recovery(focal, infecteds, index, states);
                    // when recovery, states of the focal individual and its neighbours' change
                    //updateStates(focal, states);
                }
            }
        }
    }
       
 }
	





