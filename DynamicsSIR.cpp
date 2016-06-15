#include "DynamicsSIR.h"
#include "individual.h"
#include "infecteds.h"

#include <algorithm>
#include <iostream>



namespace epinetworks {

    class Individual;
    class Infecteds;

    void DynamicsSIR::recovery(Individual &focal, Infecteds &infecteds, std::size_t index)  {
        focal.getRecovered();
        infecteds.remove(index);
    }

    // Disease induced mortality event with S replacement
    void DynamicsSIR::virulence(Individual &focal, Infecteds &infecteds, std::size_t index) {
        focal.getSusceptible();
        Individual::updateSusceptibleNeigbours(focal, 1);
        infecteds.remove(index);
    }

    
}





