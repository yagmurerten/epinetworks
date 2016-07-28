#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H

#include "pathogen.h"
#include "networkNode.h"

#include <vector>

namespace epinetworks {

	class Individual : public NetworkNode {
	public:
		explicit Individual(int numberOfContacts);
		enum class Status {
			susceptible,
			infected,
			recovered
		};
		void getInfected(Pathogen &pathogen);
		void getSusceptible();
		void getRecovered();
        double getEventRate();
        void updateEventRate();
		Individual::Status getStatus() const;
		Pathogen getPathogen(); 
		void updateSusceptibleNeigbours(int update);
		int sizeSusceptibleNeighbours();
		void setSusceptibleNumber();
        void updateStates();
        std::vector<int> getStates();
        Individual& getNeighbour(const int i){
            return *_neighbours[i];
        }
        void setNeighbour(Individual &node){
            _neighbours.push_back(&node);
        }
	private:
		Individual(const Individual &) = delete;
		Individual &operator=(const Individual &) = delete;
		int _susceptibleNeighbours;
		Pathogen _pathogen;
		Individual::Status _status;
        double _eventRate;
        std::vector<int> _states;
        std::vector<Individual *> _neighbours;
	};
}

#endif // INDIVIDUAL_H
