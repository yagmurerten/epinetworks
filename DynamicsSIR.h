#ifndef DYNAMICSSIR_H_INCLUDED
#define DYNAMICSSIR_H_INCLUDED

#include "Dynamics.h"


namespace epinetworks {

    class Individual;
    class Infecteds;
    
    class DynamicsSIR : public Dynamics {
    public:
        void recovery(Individual &focal, Infecteds &infecteds, std::size_t index, std::vector<std::vector<int>> &states);
        void virulence(Individual &focal, Infecteds &infecteds, std::size_t index);
    };

}

#endif // DYNAMICSSIR_H_INCLUDED