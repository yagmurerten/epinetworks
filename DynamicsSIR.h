#ifndef DYNAMICSSIR_H_INCLUDED
#define DYNAMICSSIR_H_INCLUDED

#include "Dynamics.h"


namespace epinetworks {

    class Individual;
    class Infecteds;
    
    class DynamicsSIR : public Dynamics {
    public:
        void transmission(Individual &focal, Infecteds &infecteds, RandomNumberGenerator &rng, double mutationSD, bool mutate) ;
        void recovery(Individual &focal, Infecteds &infecteds, std::size_t index);
        void virulence(Individual &focal, Infecteds &infecteds, std::size_t index);
    };

}

#endif // DYNAMICSSIR_H_INCLUDED