#ifndef INFECTEDS_H
#define INFECTEDS_H

#include <vector>
#include <memory>


namespace epinetworks {

    class Individual;

    class Infecteds {
    public:
        Infecteds(Individual *patientZero, int sizeNetwork);
        double getAverageVirulence();
        double getAverageVirulenceK(int minSizeNeighbour, int maxSizeNeighbour);
        int getIndividualsWithKNeighbours(int minSizeNeighbour, int maxSizeNeighbour);
        double getSDVirulence();
        double getSDVirulence(double averageVirulence);
        double getAverageNumberNeighbour();
        int getSizeInfected();
        void addInfected(Individual *newInfected);
        Individual& returnIndividual(int index);
        void remove(int index);
    private:
        std::vector<Individual*> _infecteds;
        Infecteds(const Infecteds &);
        Infecteds &operator=(const Infecteds &);
    };
}

#endif // INFECTEDS_H
