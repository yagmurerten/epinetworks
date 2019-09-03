#ifndef NETWORK_H_INCLUDED
#define NETWORK_H_INCLUDED

// #include "Array2D.h"
#include "networkNode.h"

#include <memory>
#include <string>
#include <vector>

#include "individual.h"

namespace epinetworks {

	class Network {
	private:
		typedef std::unique_ptr<Individual> UpNetworkNode;
		typedef std::vector<UpNetworkNode> NodeVector;

	public:
		typedef NodeVector::iterator iterator;

		explicit Network(size_t i);

		void setNode(UpNetworkNode node) {
			_vector.push_back(std::move(node));
		}

		NodeVector::iterator begin(){
			return _vector.begin();
		}

		NodeVector::iterator end(){
			return _vector.end();
		}
		size_t size() const {
			return _size;
		}

		bool isValid() const;

        const Individual &operator[](size_t i) const;
        Individual &operator[](size_t i);


	private:
		Network(const Network &);
		Network &operator=(const Network &);
		size_t _size;
		NodeVector _vector;
	};

}

#endif // NETWORK