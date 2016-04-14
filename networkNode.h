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
			_numberOfContacts = _neighbours.size();
			_stubs = 0;
		}

		// Returns number of stubs
		int getStubs() {
			return _stubs;
		}

		// Returns the size of neighbours
		std::size_t sizeNeighbour(){
			return _neighbours.size();
		}

		// Gets the neighbour i
		NetworkNode& getNeighbour(const int i){
			return *_neighbours[i];
		}

		// Sets node to be neighbour
		void setNeighbour(NetworkNode &node){
			_neighbours.push_back(&node);
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
			if (_neighbours.size() == 0)
				return false;
			else {
				std::vector<int> neighbourCoordinates;
				neighbourCoordinates.reserve(_neighbours.size());
				for (std::size_t i = 0u; i < _neighbours.size(); ++i) {
					NetworkNode& neighbour = *_neighbours[i];
					neighbourCoordinates.push_back(neighbour._coordinate);
				}
				std::vector<int>::iterator it;
				it = find(neighbourCoordinates.begin(), neighbourCoordinates.end(), potentialNeighbour._coordinate);
				if (it != neighbourCoordinates.end())
					return true;
				else
					return false;
			}
		}
	protected:
		explicit NetworkNode(int numberOfContacts) {
			_numberOfContacts = numberOfContacts;
			_stubs = numberOfContacts;
			_neighbours.reserve(numberOfContacts);
		}
		std::size_t _numberOfContacts;
		std::size_t _coordinate;
		//
	private:
		NetworkNode(const NetworkNode &) = delete;
		NetworkNode &operator=(const NetworkNode &) = delete;
		int _stubs;
		std::vector<NetworkNode*> _neighbours;
	};

}

#endif // NETWORKNODE__H_INCLUDED
