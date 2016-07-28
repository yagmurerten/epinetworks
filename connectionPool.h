#ifndef CONNECTIONPOOL_H_INCLUDED
#define CONNECTIONPOOL_H_INCLUDED

#include "Debug.h"

#include <vector>
#include <algorithm>
#include <cassert>

namespace epinetworks {
	
	class Individual;

	class ConnectionPool {
	public:

		bool empty() const{
			return _pool.empty();
		}

        void add(Individual &node){
			_pool.push_back(&node);
		}

		void pop_back(){
			_pool.pop_back();
		}

		void remove(std::size_t i) {
			assert(i <= _pool.size());
			if ((i + 1u) < _pool.size())
				std::swap(_pool[i], _pool.back());
			_pool.pop_back();
		}

		std::size_t size() {
			return _pool.size();
		}

        const Individual& operator[](int i) const {
            const Individual *ptr = _pool[i];
			return *ptr;
		}

        Individual& operator[](int i){
            return const_cast<Individual &>(
				static_cast<const ConnectionPool &>(*this)[i]);
		}

		void reserve(std::size_t size){
			_pool.reserve(size);
		}

        std::vector<Individual*>::iterator begin(){
			return _pool.begin();
		}

        std::vector<Individual*>::iterator end(){
			return _pool.end();
		}

	private:

		std::vector<Individual*> _pool;

	};

}

#endif // CONNECTIONPOOL_H_INCLUDED