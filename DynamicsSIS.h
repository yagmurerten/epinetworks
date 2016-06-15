#ifndef DYNAMICSSIS_H_INCLUDED
#define DYNAMICSSIS_H_INCLUDED

#include "Dynamics.h"

namespace epinetworks {

    class Individual;
    class Infecteds;

    class DynamicsSIS : public Dynamics {
    public:
        void recovery(Individual &focal, Infecteds &infecteds, std::size_t index);
        void virulence(Individual &focal, Infecteds &infecteds, std::size_t index);
    };

}

#endif // DYNAMICSSIS_H_INCLUDED