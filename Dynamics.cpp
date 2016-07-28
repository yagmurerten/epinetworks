#include "Dynamics.h"
#include "random.h"
#include "individual.h"
#include "infecteds.h"
#include "pathogen.h"

#include <iostream> 

namespace epinetworks {

    void Dynamics::transmission(Individual &focal, Infecteds &infecteds, RandomNumberGenerator &rng, double mutationSD, bool mutate)  {
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
        }
    }

    
}

