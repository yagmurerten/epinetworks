#include "DynamicsSIR.h"
#include "individual.h"
#include "infecteds.h"

#include <algorithm>
#include <iostream>



namespace epinetworks {

    class Individual;
    class Infecteds;

    void DynamicsSIR::recovery(Individual &focal, Infecteds &infecteds, std::size_t index, std::vector<std::vector<int>> &states)  {
        focal.getRecovered();
        focal.updateSusceptibleNeigbours(0);
        infecteds.remove(index);
        updateStates(focal, states);
    }

    // Disease induced mortality event with S replacement
    void DynamicsSIR::virulence(Individual &focal, Infecteds &infecteds, std::size_t index) {
        focal.getSusceptible();
        focal.updateSusceptibleNeigbours(1);
        infecteds.remove(index);
    }

    
}





