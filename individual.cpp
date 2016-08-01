#include "Individual.h"
#include "pathogen.h"
#include "random.h"

#include <algorithm>
#include <iostream>
#include <vector>



namespace epinetworks {

	Pathogen nullPathogen(0.,0.);

	// Constructor for uninfected individuals
	Individual::Individual (int numberOfContacts) :
		_status(Status::susceptible), _pathogen(nullPathogen), 
		NetworkNode(numberOfContacts), _susceptibleNeighbours(numberOfContacts),
        _eventRate(0.) { 
        _states.reserve(5);
        for (size_t i = 0u; i < numberOfContacts; ++i)
            _states.push_back(0);
    };

	 //Sets # of susceptible neighbours for one individual
	void Individual::setSusceptibleNumber(){
		_susceptibleNeighbours = sizeNeighbour();
	}

	// Becomes susceptible
	void Individual::getSusceptible() {
		_status = Status::susceptible;
		_pathogen = nullPathogen;
        _eventRate = 0.;
	}
	

	// Becomes recovered
    void Individual::getRecovered() {
        _status = Status::recovered;
		_pathogen = nullPathogen;
        _eventRate = 0.;
	}


	// Becomes infected
	void Individual::getInfected(Pathogen &pathogen) {
		_status = Status::infected;
		_pathogen = pathogen;
        _eventRate = _pathogen.getTransmission()*_susceptibleNeighbours +
            _pathogen.getRecoveryRate() + _pathogen.getVirulence();
	}

    double Individual::getEventRate(){
        return _eventRate;
    }

    void Individual::updateEventRate(){
        _eventRate = _pathogen.getTransmission()*_susceptibleNeighbours +
            _pathogen.getRecoveryRate() + _pathogen.getVirulence();
    }

	// Returns the status of an individual
	Individual::Status Individual::getStatus() const {
		return _status;
	}

	// Return the pathogen that infects that individual
	Pathogen Individual::getPathogen() {
		return _pathogen;
	}

	// Returns the size of susceptible neighbours of one individual
	int Individual::sizeSusceptibleNeighbours(){
		return _susceptibleNeighbours;
	}
	
    void Individual::updateStates(){
        _states.clear();
        _states.push_back(static_cast<int>(getStatus()));
        for (std::size_t i = 0; i < sizeNeighbour(); ++i) {
            Individual &neighbour = getNeighbour(i);
            _states.push_back(static_cast<int>(neighbour.getStatus()));
        }
    }

	// Decreases or increases # of susceptible neighbours
	void Individual::updateSusceptibleNeigbours(int update){
        updateStates();
       for (std::size_t i = 0; i < sizeNeighbour(); ++i) {
            Individual &neighbour = getNeighbour(i);
            neighbour._susceptibleNeighbours += update;
            neighbour.updateEventRate();
            neighbour.updateStates();
        }
	}

    std::vector<int> Individual::getStates() {
        return _states;
    }

}
