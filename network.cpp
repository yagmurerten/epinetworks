#include "network.h"

#include "random.h"

#include <fstream>
#include <iostream>
#include <stddef.h>
#include <sstream>
#include <string>

namespace epinetworks {

	// Constructor for network
	Network::Network(size_t i) { 
		_vector.reserve(i);
		_size = i;
	}


	// Checks if the network is fully initialized
	bool Network::isValid() const{
		for (auto it = _vector.begin(); it != _vector.end(); ++it){
			if (*it == nullptr)
				return false;
		}
		return true;
	}

	// const [] operator
    const Individual &Network::operator[](size_t i) const {
		DEBUG_ASSERT(i < _vector.size());
        const Individual *ptr = _vector[i].get();
		DEBUG_ASSERT(ptr != nullptr);
		return *ptr;
	}

	// [] operator
    Individual &Network::operator[](size_t i){
        return const_cast<Individual &>(
			static_cast<const Network &>(*this)[i]);
	}
}







