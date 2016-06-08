#include "pathogen.h"
#include "random.h"
#include <cmath> 

namespace epinetworks {

	// SIS with virulence 
	Pathogen::Pathogen(double virulence, double coefficient) :
		_virulence(virulence), _transmission(coefficient*sqrt(virulence)), 
		_coefficient(coefficient) {};

	//SIR without virulence
	Pathogen::Pathogen(double transmission) : 
		_virulence(0), _transmission(transmission),
		_coefficient(0) {};

	// Mutates pathogen with a given s.d. around 0 mean mutation effect
	Pathogen Pathogen::mutatePathogen(Pathogen const &initial, const double sdMutation,
		RandomNumberGenerator &rng) {
		double mutatedVirulence = initial._virulence;
		double mutationStep = getRandomNormal(0, sdMutation, rng);
		Pathogen mutatedPathogen(mutatedVirulence += mutationStep, initial._coefficient);
		Pathogen nullPathogen(0.0, 0.0);
		if (mutatedPathogen.getVirulence() <= 0.0)
			return nullPathogen;
		else
			return mutatedPathogen;
	}

	// Returns virulence level of a pathogen
	double Pathogen::getVirulence() const {
		return _virulence;
	}

	// Returns transmission level of a pathogen
	double Pathogen::getTransmission() const {
		return _transmission;
	}
}