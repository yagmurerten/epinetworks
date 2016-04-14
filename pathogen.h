#ifndef PATHOGEN_H
#define PATHOGEN_H

#include "random.h"

namespace epinetworks {

	class Pathogen {
	public:
		explicit Pathogen(double virulence, double coefficient);
		explicit Pathogen(double transmission);
		static Pathogen mutatePathogen(Pathogen const &initial, 
			const double sdMutation, RandomNumberGenerator &rng);
		double getVirulence() const;
		double getTransmission() const;
	private:
		double _virulence;
		double _transmission;
		double _coefficient;
	};
}


#endif // PATHOGEN_H
