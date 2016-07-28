#ifndef NETWORKNODE_H_INCLUDED
#define NETWORKNODE_H_INCLUDED

#include "Debug.h"

#include <algorithm>
#include <vector>

namespace epinetworks {

	class NetworkNode {
	public:

		virtual ~NetworkNode() = 0;

		// Returns coordinates of a given node
		// Not sure if necessary, to be checked
		const size_t getCoordinate() {
			return _coordinate;
		}

		// Sets stubs to the number of contacts
		// In order to be used to see available connections
		void setStubs() {
			_stubs = _numberOfContacts;
		}

		// Reduces the number of stubs when a connection is established
		void reduceStubs() {
			--_stubs;
		}

		// Sets the number of stubs to zero when no available connections left
		// Updates connections to zero while setting number of contacts accordingly
		void updateConnectionZero() {
            _numberOfContacts = _neighbourCoordinates.size();
			_stubs = 0;
		}

		// Returns number of stubs
		int getStubs() {
			return _stubs;
		}

		// Returns the size of neighbours
		std::size_t sizeNeighbour(){
            return _neighbourCoordinates.size();
		}

		// Gets the neighbour i
		int& getNeighbourCoord(const int i){
            return _neighbourCoordinates[i];
		}

		// Sets node to be neighbour
		void setNeighbourCoord(NetworkNode &node){
            _neighbourCoordinates.push_back(node._coordinate);
		}

		// Gets # of contacts
		std::size_t getNumberOfContacts() {
			return _numberOfContacts;
		}

		// Sets coordinates 
		void setCoordinates(std::size_t coordinate){
			_coordinate = coordinate;
		}

		// Checks if two nodes are already neighbours
		bool isNeighbour(NetworkNode &potentialNeighbour) {
            if (_neighbourCoordinates.size() == 0)
				return false;
			else {
				std::vector<int>::iterator it;
                it = find(_neighbourCoordinates.begin(), _neighbourCoordinates.end(), potentialNeighbour._coordinate);
                if (it != _neighbourCoordinates.end())
					return true;
				else
					return false;
			}
		}
	protected:
		explicit NetworkNode(int numberOfContacts) {
			_numberOfContacts = numberOfContacts;
			_stubs = numberOfContacts;
            _neighbourCoordinates.reserve(numberOfContacts);
		}
		std::size_t _numberOfContacts;
		std::size_t _coordinate;
		//
	private:
		NetworkNode(const NetworkNode &) = delete;
		NetworkNode &operator=(const NetworkNode &) = delete;
		int _stubs;
        std::vector<int> _neighbourCoordinates;
	};

}

#endif // NETWORKNODE__H_INCLUDED
