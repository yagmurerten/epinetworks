#include "DynamicsSIS.h"
#include "individual.h"
#include "infecteds.h"

#include <iostream>


namespace epinetworks {

    void DynamicsSIS::recovery(Individual &focal, Infecteds &infecteds, std::size_t index)  {
        focal.getRecovered(DynamicsType::SIS);
        Individual::updateSusceptibleNeigbours(focal, Individual::UpdateRule::Up);
        infecteds.remove(index);
    }

    // Disease induced mortality event with S replacement
    void DynamicsSIS::virulence(Individual &focal, Infecteds &infecteds, std::size_t index) {
        focal.getSusceptible();
        Individual::updateSusceptibleNeigbours(focal, Individual::UpdateRule::Up);
        infecteds.remove(index);
    }

    void DynamicsSIS::transmission(Individual &focal, Infecteds &infecteds, RandomNumberGenerator &rng, double mutationSD, bool mutate)  {
        if (focal.sizeSusceptibleNeighbours() > 0) {
            int newInfected;
            do {
                newInfected = getRandom(focal.getNumberOfContacts(), rng);
            } while (dynamic_cast<Individual &>(focal.getNeighbour(newInfected)).getStatus() != Individual::Status::susceptible);
            Individual &newInfectedNeighbour = dynamic_cast<Individual &>(focal.getNeighbour(newInfected));
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
            Individual::updateSusceptibleNeigbours(newInfectedNeighbour, Individual::UpdateRule::Down);
            infecteds.addInfected(&newInfectedNeighbour);
        }
    }
}





