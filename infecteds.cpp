#include "infecteds.h"
#include "individual.h"
#include "pathogen.h"
#include "random.h"

#include <algorithm>
#include <iostream>
#include <vector>



namespace epinetworks {
    Infecteds::Infecteds(Individual *patientZero, int sizeNetwork) {
        _infecteds.reserve(sizeNetwork);
        _infecteds.push_back(patientZero);
    }

    Individual& Infecteds::returnIndividual(int index) {
        return *_infecteds[index];
    }
    int Infecteds::getSizeInfected(){
        return _infecteds.size();
    }
    double Infecteds::getAverageVirulence() {
        double totalVirulence = 0.;
        for (std::size_t i = 0; i < _infecteds.size(); ++i) {
            Individual &ind = *_infecteds[i];
            totalVirulence += ind.getPathogen().getVirulence();
        }
        return totalVirulence / _infecteds.size();
    }

    double Infecteds::getAverageVirulenceK(int minSizeNeighbour, int maxSizeNeighbour) {
        double totalVirulenceK = 0.;
        double count = 0.;
        for (std::size_t i = 0; i < _infecteds.size(); ++i) {
            Individual &ind = *_infecteds[i];
            if (ind.sizeNeighbour() > minSizeNeighbour && ind.sizeNeighbour() <= maxSizeNeighbour) {
                totalVirulenceK += ind.getPathogen().getVirulence();
                ++count;
            }
        }
        if (count != 0)
            return totalVirulenceK / count;
        return 0.;
    }

    int Infecteds::getIndividualsWithKNeighbours(int minSizeNeighbour, int maxSizeNeighbour) {
        double count = 0.;
        for (std::size_t i = 0; i < _infecteds.size(); ++i) {
            Individual &ind = *_infecteds[i];
            if (ind.sizeNeighbour() > minSizeNeighbour && ind.sizeNeighbour() <= maxSizeNeighbour) {
                ++count;
            }
        }
        return count;
    }

    double Infecteds::getSDVirulence(double averageVirulence) {
        double sd = 0.;
        for (std::size_t i = 0; i < _infecteds.size(); ++i) {
            Individual &ind = *_infecteds[i];
            sd += pow((ind.getPathogen().getVirulence() - averageVirulence), 2.0);
        }
        sd = sqrt(sd / _infecteds.size());
        return sd;
    }

    double Infecteds::getAverageNumberNeighbour(){
        double totalNeighbours = 0;
        for (std::size_t i = 0; i < _infecteds.size(); ++i) {
            Individual &ind = *_infecteds[i];
            totalNeighbours += ind.sizeNeighbour();
        }
        double averageNeighbours = totalNeighbours / _infecteds.size();
        return averageNeighbours;
    }


    void Infecteds::addInfected(Individual *newInfected) {
        _infecteds.push_back(newInfected);
    }

    void Infecteds::remove(int index) {
        std::swap(_infecteds[index], _infecteds.back());
        _infecteds.pop_back();
    }

    double Infecteds::getSDVirulence(){
        double totalVirulence = 0.;
        for (std::size_t i = 0; i < _infecteds.size(); ++i) {
            Individual &ind = *_infecteds[i];
            totalVirulence += ind.getPathogen().getVirulence();
        }
        double averageVirulence = totalVirulence / _infecteds.size();
        double sd = 0.;
        for (std::size_t i = 0; i < _infecteds.size(); ++i) {
            Individual &ind = *_infecteds[i];
            sd += pow((ind.getPathogen().getVirulence() - averageVirulence), 2.0);
        }
        return sd;
    }
}
