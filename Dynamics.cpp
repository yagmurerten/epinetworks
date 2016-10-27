#include "Dynamics.h"
#include "random.h"
#include "individual.h"
#include "infecteds.h"
#include "pathogen.h"

#include <iostream> 

namespace epinetworks {

    void Dynamics::updateStates(Individual &ind, std::vector<std::vector<int>> &states){
        states.clear();
        states.reserve(ind.sizeNeighbour() + 1);
        std::vector<int> coordStates;
        coordStates.reserve(ind.sizeNeighbour() + 2);
        coordStates = ind.getStates();
        coordStates.push_back(ind.getCoordinate());
        states.push_back(coordStates);
        for (size_t i = 0; i < ind.sizeNeighbour(); ++i) {
            std::vector<int> coordStates;
            coordStates.reserve(ind.getNeighbour(i).sizeNeighbour() + 2);
            coordStates = ind.getNeighbour(i).getStates();
            coordStates.push_back(ind.getNeighbour(i).getCoordinate());
            states.push_back(coordStates);
        }
    }

    void Dynamics::transmission(Individual &focal, Infecteds &infecteds, RandomNumberGenerator &rng, 
        double mutationSD, bool mutate, std::vector<std::vector<int>> &states)  {
        if (focal.sizeSusceptibleNeighbours() > 0) {
            int newInfected;
            do {
                newInfected = getRandom(focal.getNumberOfContacts(), rng);
            } while (focal.getNeighbour(newInfected).getStatus() != Individual::Status::susceptible);
            Individual &newInfectedNeighbour = focal.getNeighbour(newInfected);
            if (mutate) {
                Pathogen mutatedPathogen = Pathogen::mutatePathogen(focal.getPathogen(), mutationSD, rng);
                newInfectedNeighbour.getInfected(mutatedPathogen);
            }
            else {
                Pathogen oldPathogen = focal.getPathogen();
                newInfectedNeighbour.getInfected(oldPathogen);
                if (newInfectedNeighbour.sizeSusceptibleNeighbours() < 0)
                    std::cout << "stop 2" << std::endl;
            }
            newInfectedNeighbour.updateSusceptibleNeigbours(-1);
            infecteds.addInfected(&newInfectedNeighbour);
            //updateStates(newInfectedNeighbour, states);
        }
    }

    
}

