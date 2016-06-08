#ifndef PATHOGEN_H
#define PATHOGEN_H

#include "random.h"

namespace epinetworks {

	class Pathogen {
	public:
        explicit Pathogen(double virulence, double coefficient, double recoveryRate);
        explicit Pathogen(double transmission, double recoveryRate);
		static Pathogen mutatePathogen(Pathogen const &initial, 
			const double sdMutation, RandomNumberGenerator &rng);
		double getVirulence() const;
        double getRecoveryRate() const;
		double getTransmission() const;
	private:
		double _virulence;
		double _transmission;
		double _coefficient;
        double _recoveryRate;
	};
}


#endif // PATHOGEN_H
