#include "individual.h"
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
		NetworkNode(numberOfContacts) {};

	// Sets # of susceptible neighbours for one individual
	void Individual::setSusceptibleNumber(Individual &individual){
		individual._susceptibleNeighbours = individual.sizeNeighbour();
	}

	
	// Becomes susceptible
	void Individual::getSusceptible() {
		_status = Status::susceptible;
		_pathogen = nullPathogen;
	}
	

	// Becomes recovered
	void Individual::getRecovered() {
		_status = Status::recovered;
		_pathogen = nullPathogen;
	}


	// Becomes infected
	void Individual::getInfected(Pathogen &pathogen) {
		_status = Status::infected;
		_pathogen = pathogen;
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
	
	// Decreases or increases # of susceptible neighbours
	void Individual::updateSusceptibleNeigbours(Individual &individual,  UpdateRule rule){
		for (std::size_t i = 0; i < individual.sizeNeighbour(); ++i) {
			Individual &neighbour = dynamic_cast<Individual &>(individual.getNeighbour(i));
			if (rule == UpdateRule::Up)
				++neighbour._susceptibleNeighbours;
			if (rule == UpdateRule::Down) {
				--neighbour._susceptibleNeighbours;
			}
		}
	}
}
