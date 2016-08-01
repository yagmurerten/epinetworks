#ifndef DYNAMICS_H_INCLUDED
#define DYNAMICS_H_INCLUDED

#include "random.h"

#include <vector>

namespace epinetworks {

    class Individual;
    class Infecteds;
    
    class Dynamics {
    public:
        virtual ~Dynamics() {};

        enum class DynamicsType {
            SIS,
            SIR
        };

        virtual void recovery(Individual &focal, Infecteds &infecteds, std::size_t index, std::vector<std::vector<int>> &states) = 0;

        void transmission(Individual &focal, Infecteds &infecteds, RandomNumberGenerator &rng, 
            double mutationSD, bool mutate, std::vector<std::vector<int>> &states);

        virtual void virulence(Individual &focal, Infecteds &infecteds, std::size_t index) = 0;
        
        void updateStates(Individual &ind, std::vector<std::vector<int>> &states);

    };

	
}

#endif // DYNAMICS_H_INCLUDED