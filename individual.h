#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H

#include "Gillespie.h"
#include "pathogen.h"
#include "networkNode.h"

namespace epinetworks {

	class Individual : public NetworkNode {
	public:
		explicit Individual(int numberOfContacts);
		enum class Status {
			susceptible,
			infected,
			recovered
		};
		enum class UpdateRule {
			Up,
			Down
		};
		void getInfected(Pathogen &pathogen);
		void getSusceptible();
		void getRecovered(Gillespie::DynamicsType type);
		Individual::Status getStatus() const;
		Pathogen getPathogen(); 
		static void updateSusceptibleNeigbours(Individual &individual, UpdateRule rule);
		int sizeSusceptibleNeighbours();
		static void setSusceptibleNumber(Individual &individual);
	private:
		Individual(const Individual &) = delete;
		Individual &operator=(const Individual &) = delete;

		int _susceptibleNeighbours;
		Pathogen _pathogen;
		Individual::Status _status;
	};
}

#endif // INDIVIDUAL_H
