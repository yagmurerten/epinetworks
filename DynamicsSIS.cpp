#include "DynamicsSIS.h"
#include "individual.h"
#include "infecteds.h"

#include <iostream>


namespace epinetworks {

    void DynamicsSIS::recovery(Individual &focal, Infecteds &infecteds, std::size_t index, std::vector<std::vector<int>> &states)  {
        focal.getSusceptible();
        focal.updateSusceptibleNeigbours(1);
        infecteds.remove(index);
        updateStates(focal, states);
    }

    // Disease induced mortality event with S replacement
    void DynamicsSIS::virulence(Individual &focal, Infecteds &infecteds, std::size_t index) {
        focal.getSusceptible();
        focal.updateSusceptibleNeigbours(1);
        infecteds.remove(index);
    }

}





